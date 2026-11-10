"""Retrieve and curate nuclear/eukaryotic SSU (18S) DNA by feature annotation.

Length is a tier, not a synonym or proof of completeness. Preserve genomic
introns and circular-origin continuity; never construct spliced DNA targets.
"""
from __future__ import annotations

from runtime import storage_path, settings
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import hashlib
import html
import json
from pathlib import Path
import re
from urllib.parse import urlencode
import xml.etree.ElementTree as ET

import amoebacrisp as app
from ssu18s_policy import filter_records,MIN_18S_BP,MAX_18S_BP

EUTILS='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
ALIASES=['18S','18S ribosomal RNA','18S rRNA','18S ribosomal DNA','18S rDNA',
         'small subunit ribosomal RNA','small-subunit ribosomal RNA','small subunit rRNA',
         'small subunit ribosomal DNA','SSU rRNA','SSU rDNA','SSUrDNA','SSUrRNA',
         'SSU ribosomal RNA','small ribosomal subunit RNA','small subunit RNA','small subunit rDNA']


def qualifier_map(feature):
    result={}
    for q in feature.findall('GBFeature_quals/GBQualifier'):
        result.setdefault(q.findtext('GBQualifier_name',''),[]).append(q.findtext('GBQualifier_value',''))
    return result


def identifies_18s(text):
    normalized=re.sub(r'[-_\s]+',' ',re.sub(r'\([^)]*\)',' ',text.casefold()))
    return bool(re.search(r'\b18\s*s(?=\b|rrna|rdna)',normalized) or
                re.search(r'\bsmall subunit (?:ribosomal (?:rna|dna)|rrna|rdna)\b',normalized) or
                re.search(r'\bsmall ribosomal subunit (?:rna|dna)\b',normalized) or
                re.search(r'\bssu ribosomal (?:rna|dna)\b',normalized) or
                re.search(r'\bssu\s*(?:rrna|rdna)\b',normalized))


def contiguous_feature_dna(sequence,feature,circular=False,accession=None):
    intervals=[]
    strands=set()
    for interval in feature.findall('GBFeature_intervals/GBInterval'):
        remote=interval.findtext('GBInterval_accession')
        if accession and remote and remote.split('.')[0]!=accession.split('.')[0]:
            raise ValueError('remote_feature_requires_review')
        begin=interval.findtext('GBInterval_from') or interval.findtext('GBInterval_point')
        end=interval.findtext('GBInterval_to') or begin
        if not begin:
            raise ValueError('feature_has_no_local_intervals')
        a,b=sorted([int(begin)-1,int(end)-1])
        if a<0 or b>=len(sequence):
            raise ValueError('feature_out_of_bounds')
        intervals.append((a,b+1))
        flag=interval.find('GBInterval_iscomp')
        strands.add(flag is not None and flag.get('value','true')=='true')
    if not intervals or len(strands)!=1:
        raise ValueError('missing_or_mixed_strand_feature')
    intervals.sort()
    reverse=next(iter(strands))
    start,end=intervals[0][0],max(b for a,b in intervals)
    wraps=False
    if circular and len(intervals)>1:
        # Exclude the largest gap around a circular molecule. Retain internal
        # gaps as genomic DNA (e.g. group I introns), not exon junctions.
        gaps=[(intervals[i+1][0]-intervals[i][1],i) for i in range(len(intervals)-1)]
        gaps.append((intervals[0][0]+len(sequence)-intervals[-1][1],len(intervals)-1))
        _,cut=max(gaps)
        if cut<len(intervals)-1:
            start,end=intervals[cut+1][0],intervals[cut][1]
            wraps=True
    dna=sequence[start:]+sequence[:end] if wraps else sequence[start:end]
    if reverse:dna=app.rc(dna)
    return dna,start+1,end,'-' if reverse else '+',wraps,len(intervals)


def source_info(record):
    organism=record.findtext('GBSeq_organism','')
    source=next((f for f in record.findall('GBSeq_feature-table/GBFeature')
                 if f.findtext('GBFeature_key')=='source'),None)
    q=qualifier_map(source) if source is not None else {}
    taxid=next((v.split(':',1)[1] for v in q.get('db_xref',[]) if v.startswith('taxon:')),'')
    organelle=';'.join(q.get('organelle',[]))
    return {'organism':organism,'taxid':taxid,'strain':';'.join(q.get('strain',[])),
            'isolate':';'.join(q.get('isolate',[])),'organelle':organelle,
            'mol_type':record.findtext('GBSeq_moltype','')}


def classify(organism):
    name=organism.casefold()
    if name=='naegleria gruberi' or name.startswith('naegleria gruberi '):return 'target_N_gruberi'
    if name.startswith('entamoeba '):return 'background_Entamoeba'
    if name.startswith('naegleria '):
        if re.match(r'naegleria\s+(?:sp\.?\b|cf\.?\b|aff\.?\b|uncultured\b)',name):
            return 'unresolved_Naegleria'
        return 'background_Naegleria_named'
    return 'taxonomy_review'


def curate_record(record):
    accession=record.findtext('GBSeq_accession-version','')
    title=record.findtext('GBSeq_definition','')
    sequence=record.findtext('GBSeq_sequence','').upper()
    info=source_info(record)
    base={'accession':accession,'title':title,**info,'source_record_length':len(sequence),
          'role':classify(info['organism'])}
    rejects=[]
    if not sequence:return [],[{**base,'reason':'no_sequence'}]
    if re.search('mitochond|plastid|chloroplast',info['organelle'],re.I) or re.search(r'\bmitochondri(?:al|on)\b',title,re.I):
        return [],[{**base,'reason':'organelle_SSU_not_nuclear_18S'}]
    features=[]
    for feature in record.findall('GBSeq_feature-table/GBFeature'):
        if feature.findtext('GBFeature_key') not in ['rRNA','misc_RNA']:continue
        q=qualifier_map(feature)
        # Do not infer an 18S rRNA feature from a surrounding operon/ITS note.
        label='; '.join(q.get('product',[])+q.get('gene',[])+q.get('standard_name',[]))
        mixed=re.search(r'internal transcribed|5[.]8s|large subunit|28s',label,re.I)
        if identifies_18s(label) and not mixed:features.append((feature,q,label))
    if not features:
        return [],[{**base,'reason':'no_explicit_18S_rRNA_feature_requires_review'}]
    rows=[]
    for number,(feature,q,label) in enumerate(features,1):
        location=feature.findtext('GBFeature_location','')
        try:
            dna,start,end,strand,wraps,parts=contiguous_feature_dna(sequence,feature,record.findtext('GBSeq_topology')=='circular',accession)
        except ValueError as exc:
            rejects.append({**base,'reason':str(exc)})
            continue
        partial=('<' in location or '>' in location or 'partial' in q or
                 ('partial' in title.casefold() and len(features)==1 and len(sequence)==len(dna)))
        length=len(dna)
        tier='eligible_800_4000' if MIN_18S_BP<=length<=MAX_18S_BP else ('below_800_excluded' if length<MIN_18S_BP else 'above_4000_excluded')
        row={**base,'id':app.slug(accession)+f'__18S_{number}','product':label,
             'feature_key':feature.findtext('GBFeature_key'),
             'feature_location':location,'record_start_1based':start,'record_end_1based':end,'strand':strand,
             'circular_origin_wrap':wraps,'annotation_intervals':parts,'genomic_introns_preserved':True,
             'length':length,'length_tier':tier,'annotated_partial':partial,
             'ambiguous_bases':sum(c not in 'ACGT' for c in dna),'sha256':hashlib.sha256(dna.encode()).hexdigest(),
             'DNA_supported':info['mol_type'].casefold() in ('dna','genomic dna','other dna'),
             'sequence':dna}
        if length<20:
            rejects.append({**base,'reason':'18S_fragment_shorter_than_20nt'})
        else:rows.append(row)
    return rows,rejects


def retrieve(root):
    root.mkdir(parents=True,exist_ok=True)
    queries={}
    for genus in ['Naegleria','Entamoeba']:
        # Broad discovery; length is intentionally not part of this query.
        aliases=' OR '.join('"'+alias+'"[All Fields]' for alias in ALIASES)
        term=f'{genus}[Organism:exp] AND ({aliases})'
        result=json.loads(app.fetch(EUTILS+'esearch.fcgi?'+urlencode({'db':'nuccore','term':term,'retmode':'json','retmax':10000})))
        total=int(result['esearchresult']['count'])
        ids=result['esearchresult']['idlist']
        if total>len(ids):raise RuntimeError('检索超过当前批量上限，不能静默截断')
        app.write_json(root/(genus+'_search.json'),{'term':term,'result':result})
        queries[genus]=ids
        app.log(f'{genus}：宽检索 {total} 条记录，待逐条核对 rRNA 注释。')
    ids=sorted(set(v for values in queries.values() for v in values),key=int)
    jobs=[]
    raw=root/'raw_genbank';raw.mkdir(exist_ok=True)
    for i in range(0,len(ids),100):
        block=ids[i:i+100]
        identity=hashlib.sha256(','.join(block).encode()).hexdigest()[:14]
        url=EUTILS+'efetch.fcgi?'+urlencode({'db':'nuccore','id':','.join(block),'rettype':'gb','retmode':'xml'})
        jobs.append((url,raw/f'{i//100:04d}_{identity}.xml'))
    def download(job):
        url,path=job
        app.fetch(url,path)
        tree=ET.parse(path)
        if tree.getroot().tag!='GBSet':raise RuntimeError(f'非 GenBank XML 返回：{path}')
        return path
    files=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending=[pool.submit(download,job) for job in jobs]
        for index,future in enumerate(as_completed(pending),1):
            files.append(future.result())
            app.log(f'下载与验证 GenBank XML：{index}/{len(jobs)} 批')
    all_rows,rejected,accessions=[],[],set()
    for path in sorted(files):
        for record in ET.parse(path).getroot().findall('GBSeq'):
            accession=record.findtext('GBSeq_accession-version','')
            if accession in accessions:continue
            accessions.add(accession)
            rows,skip=curate_record(record)
            all_rows.extend(rows);rejected.extend(skip)
    # Keep every accession in the manifest. FASTA representatives are deduplicated
    # per organism and DNA sequence, never across named species.
    groups={}
    for row in all_rows:
        key=(row['organism'],row['sequence'])
        groups.setdefault(key,[]).append(row)
    unique=[]
    for members in groups.values():
        representative=min(members,key=lambda r:(not r['DNA_supported'],r['annotated_partial'],r['ambiguous_bases'],r['id']))
        representative={**representative,'source_accessions':';'.join(sorted({m['accession'] for m in members})),
                        'source_record_count':len(members)}
        unique.append(representative)
    unique.sort(key=lambda r:(r['role'],r['organism'],r['id']))
    app.write_json(root/'curated_records.json',unique)
    columns=[k for k in unique[0] if k!='sequence'] if unique else ['id']
    app.write_tsv(root/'curated_18S_manifest.tsv',unique,columns)
    app.write_tsv(root/'excluded_or_review.tsv',rejected,
                  ['accession','organism','taxid','title','source_record_length','reason'])
    app.write_fasta(root/'all_curated_18S.fasta',((r['id'],r['sequence']) for r in unique))
    eligible=filter_records(unique)
    panels={
        'target_priority': [r for r in eligible if r['role']=='target_N_gruberi' and not r['annotated_partial'] and r['ambiguous_bases']==0 and r['DNA_supported']],
        'target_partial_or_review': [r for r in eligible if r['role']=='target_N_gruberi' and not (not r['annotated_partial'] and r['ambiguous_bases']==0 and r['DNA_supported'])],
        'known_background': [r for r in eligible if r['role'].startswith('background_')],
        'unresolved_Naegleria': [r for r in eligible if r['role']=='unresolved_Naegleria'],
    }
    for name,rows in panels.items():
        app.write_fasta(root/(name+'.fasta'),((r['id'],r['sequence']) for r in rows))
    species={}
    for row in unique:
        stats=species.setdefault(row['organism'],{'organism':row['organism'],'role':row['role'],'unique_18S_sequences':0,
                                                 'priority_length':0,'short_fragments':0,'long_annotated':0,'partial_annotation':0})
        stats['unique_18S_sequences']+=1
        stats['priority_length']+=MIN_18S_BP<=row['length']<=MAX_18S_BP
        stats['short_fragments']+=row['length']<MIN_18S_BP
        stats['long_annotated']+=row['length']>MAX_18S_BP
        stats['partial_annotation']+=row['annotated_partial']
    app.write_tsv(root/'species_coverage.tsv',list(species.values()),
                  ['organism','role','unique_18S_sequences','priority_length','short_fragments','long_annotated','partial_annotation'])
    summary={'retrieved_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'aliases':ALIASES,
             'min_input_bp':MIN_18S_BP,'max_input_bp':MAX_18S_BP,'eligible_by_length':len(eligible),'excluded_by_length':len(unique)-len(eligible),
             'raw_search_counts':{k:len(v) for k,v in queries.items()},'downloaded_accessions':len(accessions),
             'annotated_18S_features':len(all_rows),'unique_sequences_within_organism':len(unique),
             'organism_labels':len(species),'panel_counts':{k:len(v) for k,v in panels.items()},
             'partial_background_is_not_evidence_for_uncovered_sites':True,
             'full_length_not_proven_by_length_alone':True,'introns_preserved_for_DNA_design':True}
    app.write_json(root/'curation_summary.json',summary)
    app.log(json.dumps(summary,ensure_ascii=False,indent=2))
    write_report(root,summary,species)
    return root


def write_report(root,summary,species):
    rows=''.join('<tr>'+''.join('<td>'+html.escape(str(value))+'</td>' for value in
                 [r['organism'],r['role'],r['unique_18S_sequences'],r['priority_length'],r['short_fragments'],r['long_annotated']])+'</tr>'
                 for r in sorted(species.values(),key=lambda r:r['organism']))
    page=f'''<!doctype html><meta charset="utf-8"><title>18S rDNA 数据核对</title>
    <style>body{{font:16px system-ui;margin:36px;line-height:1.7}}td,th{{border:1px solid #ccc;padding:6px}}table{{border-collapse:collapse}}</style>
    <h1>奈格里属与内阿米巴属：18S rDNA 数据核对</h1>
    <p>按 18S/SSU 多别称检索，再按明确的 rRNA 或单一SSU产品的 misc_RNA feature 提取；剔除线粒体/质体 SSU。无法明确定位18S、或仅注释整段混合rDNA单元的记录放入待复核表。</p>
    <p>800–4000 bp（含边界）是比较、脱靶和设计输入的硬范围。范围外记录仅在原始数据中留档，不进入任何筛选。DNA内含子保留后计算长度；不构造剪接靶点。</p>
    <p>去重后明确注释序列：{summary['unique_sequences_within_organism']}；目标优先序列：{summary['panel_counts']['target_priority']}；已知背景：{summary['panel_counts']['known_background']}；未定种奈格里序列：{summary['panel_counts']['unresolved_Naegleria']}。</p>
    <p>同一物种的相同序列合并为一个代表，原 accession 保留在清单。不同物种即使序列相同也不合并。记录数不等于独立株系数。</p>
    <p>局部背景在某位置没有序列覆盖时，不能当作该位置阴性。未定种奈格里序列可能含目标种，单独保留。</p>
    <p><a href="curated_18S_manifest.tsv">序列清单</a> · <a href="species_coverage.tsv">物种数据覆盖</a> · <a href="all_curated_18S.fasta">18S FASTA</a> · <a href="excluded_or_review.tsv">排除与待复核</a></p>
    <table><tr><th>物种标签</th><th>分组</th><th>原始不同序列</th><th>800–4000 bp</th><th>低于800，排除</th><th>高于4000，排除</th></tr>{rows}</table>'''
    (root/'18S数据报告.html').write_text(page,encoding='utf-8')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default=str(storage_path('18S_analysis')))
    args=parser.parse_args()
    retrieve(Path(args.output))


if __name__=='__main__':main()

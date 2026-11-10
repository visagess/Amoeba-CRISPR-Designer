"""18S-only comparative profiles and actual CaSilico designs for benign N. gruberi."""
from __future__ import annotations

from runtime import storage_path, settings
import argparse
from array import array
import csv
import html
import json
import os
from pathlib import Path
import re
import subprocess

import amoebacrisp as app
import region_pipeline as regions
from ssu18s_policy import filter_records,MIN_18S_BP,MAX_18S_BP

BASE_BITS={'A':1,'C':2,'G':4,'T':8,'-':32}


def runs(mask,minimum=20):
    start=None
    for i,good in enumerate(list(mask)+[False]):
        if good and start is None:start=i
        elif not good and start is not None:
            if i-start>=minimum:yield start,i
            start=None


def conserved_profile(reference,cohort,observations,aligned_positions,aligned_steps,destination,prefix):
    conserved=[]
    for position,base in enumerate(reference):
        conserved.append(bool(cohort) and all(observations[r['id']][position]==BASE_BITS.get(base,0)
                         and aligned_positions[r['id']][position]>=0 for r in cohort))
    segments=[];start=None
    for i,good in enumerate(conserved+[False]):
        if not good:
            if start is not None and i-start>=20:segments.append((start,i))
            start=None
        elif start is None:start=i
        elif any(aligned_steps[r['id']][i]!=aligned_steps[r['id']][i-1]
                 or aligned_positions[r['id']][i]-aligned_positions[r['id']][i-1]!=aligned_steps[r['id']][i]
                 for r in cohort):
            if i-start>=20:segments.append((start,i))
            start=i
    rows=[]
    for i,(start,end) in enumerate(segments,1):
        fragment=reference[start:end]
        if any(fragment not in r['sequence'] and app.rc(fragment) not in r['sequence'] for r in cohort):
            raise RuntimeError('保守区未在每条DNA序列中连续存在')
        rows.append({'region':f'{prefix}_{i:03d}','reference_start_1based':start+1,'reference_end_1based_inclusive':end,
                     'length':end-start,'sequences_in_cohort':len(cohort),'sequence':reference[start:end]})
    app.write_tsv(destination/(prefix+'.tsv'),rows,
                  ['region','reference_start_1based','reference_end_1based_inclusive','length','sequences_in_cohort','sequence'])
    app.write_fasta(destination/(prefix+'.fasta'),((r['region'],r['sequence']) for r in rows))
    return rows


def sequence_covered(row,span,hsp_intervals):
    # Ambiguity or RNA-only evidence cannot establish a negative DNA result.
    return (row.get('DNA_supported',row.get('mol_type')=='DNA') and row['ambiguous_bases']==0 and
            any(left<=span[0] and right>=span[1] for left,right in hsp_intervals.get(row['id'],[])))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',default=str(storage_path('18S_analysis')))
    parser.add_argument('--threads',type=int,default=4)
    args=parser.parse_args()
    root=Path(args.data);out=root/'comparison';out.mkdir(exist_ok=True)
    temporary=Path(str(storage_path('temp')));temporary.mkdir(parents=True,exist_ok=True)
    os.environ.update(TEMP=str(temporary),TMP=str(temporary))
    config=app.default_config()
    all_records=json.loads((root/'curated_records.json').read_text())
    records=filter_records(all_records)
    app.write_json(out/'analysis_input_records.json',records)
    analysis_fasta=out/'analysis_input_18S.fasta'
    app.write_fasta(analysis_fasta,((r['id'],r['sequence']) for r in records))
    by_id={r['id']:r for r in records}
    targets=[r for r in records if r['role']=='target_N_gruberi' and
             not r['annotated_partial'] and r['ambiguous_bases']==0 and r.get('DNA_supported',False)]
    if not targets:raise RuntimeError('没有适合DNA设计的目标近全长记录')
    reference=next((r for r in targets if r['accession']=='AB298288.1'),targets[0])
    refdna=reference['sequence']
    app.write_fasta(out/'reference_18S.fasta',[(reference['id'],refdna)])
    app.log(f"参考18S：{reference['accession']}，{len(refdna)} bp；不同优先目标序列 {len(targets)} 条。")
    db=out/'sequence_db'/'database'
    regions.build_database(analysis_fasta,db,config,out/'makeblastdb.log')
    alignment=out/'reference_alignments.tsv'
    regions.run_tool('blastn',['-query',app.linux_path(out/'reference_18S.fasta'),'-db',app.linux_path(db),
                     '-task','blastn','-word_size','11','-strand','plus','-dust','no','-evalue','10',
                     '-max_target_seqs',str(max(1000,len(records))),'-num_threads',str(args.threads),
                     '-outfmt','6 qseqid sseqid qstart qend sstart send pident qseq sseq bitscore',
                     '-out',app.linux_path(alignment)],out/'alignment.log',config)
    observations={r['id']:bytearray(len(refdna)) for r in records}
    best_scores={r['id']:[-1.0]*len(refdna) for r in records}
    aligned_positions={r['id']:array('i',[-1])*len(refdna) for r in records}
    aligned_steps={r['id']:array('b',[0])*len(refdna) for r in records}
    hsp_intervals={}
    with alignment.open() as handle:
        for line in handle:
            _,subject,start,end,ss,se,identity,qseq,sseq,score=line.rstrip().split('\t')
            if subject not in by_id:raise RuntimeError('比对序列ID不能映射回清单')
            if float(identity)<70:continue
            score=float(score)
            start,end=int(start)-1,int(end)
            if float(identity)>=80 and end-start>=50:hsp_intervals.setdefault(subject,[]).append((start,end))
            position=start
            subject_position=int(ss)-1
            subject_step=1 if int(ss)<=int(se) else -1
            for qbase,sbase in zip(qseq.upper(),sseq.upper()):
                if qbase=='-':
                    if sbase!='-':subject_position+=subject_step
                    continue
                if not 0<=position<len(refdna):raise ValueError('参考比对坐标越界')
                bit=BASE_BITS.get(sbase,16)
                if score>best_scores[subject][position]:
                    observations[subject][position]=bit
                    best_scores[subject][position]=score
                    aligned_positions[subject][position]=subject_position if sbase!='-' else -1
                    aligned_steps[subject][position]=subject_step
                elif score==best_scores[subject][position]:
                    observations[subject][position] |= bit
                    if (aligned_positions[subject][position]!=subject_position or
                            aligned_steps[subject][position]!=subject_step):
                        aligned_positions[subject][position]=-2
                        aligned_steps[subject][position]=0
                if sbase!='-':subject_position+=subject_step
                position+=1
            if position!=end:raise ValueError('比对序列与参考坐标长度不一致')
    full_genus=[r for r in records if r['role'] in ['target_N_gruberi','background_Naegleria_named']
                and r['ambiguous_bases']==0 and r.get('DNA_supported',False)]
    genus_regions=conserved_profile(refdna,full_genus,observations,aligned_positions,aligned_steps,out,'Naegleria_observed_conserved_18S')
    target_regions=conserved_profile(refdna,targets,observations,aligned_positions,aligned_steps,out,'Ngruberi_observed_conserved_18S')
    app.write_json(out/'comparison_cohorts.json',{'reference':reference['id'],'target_priority':[r['id'] for r in targets],
                   'genus_priority':[r['id'] for r in full_genus],
                   'genus_organism_labels':sorted({r['organism'] for r in full_genus}),
                   'conservation_rule':'Identical unambiguous aligned bases AND consecutive physical DNA positions in every cohort sequence; no intron/junction bridging; missing/unaligned bases fail.',
                   'scope':'observed_sequences_not_all_described_species_or_all_strains',
                   'min_input_bp':MIN_18S_BP,'max_input_bp':MAX_18S_BP})
    design=out/'CaSilico';design.mkdir(exist_ok=True)
    target_samples=[]
    for r in targets:
        path=out/'target_records'/(r['id']+'.fasta')
        app.write_fasta(path,[(r['id'],r['sequence'])])
        target_samples.append({'id':r['id'],'name':r['organism']+' '+r['accession'],'genome':str(path)})
    group='Ngruberi_18S' if reference['organism']=='Naegleria gruberi' else app.slug(reference['organism'])+'_18S'
    group_dir=design/'casilico'/group
    previous=group_dir/'targets.fasta'
    tables=list(group_dir.rglob('Spacers_Information.Html')) if group_dir.exists() else []
    expected=[(r['id'],r['sequence']) for r in targets]
    if (group_dir/'completed.txt').exists() and previous.exists() and list(app.read_fasta(previous))==expected and len(tables)==1:
        raw=app.parse_casilico(tables[0])
    else:
        app.log('最终设计工具：CaSilico casV_A；输入为 N. gruberi 优先近全长18S DNA序列。')
        raw=app.design(group,[r['id'] for r in targets],{r['id']:{'dna':r['sequence']} for r in targets},design,config,args.threads)
    candidates=[{**r,'group':group} for r in raw]
    app.write_json(design/'casilico_candidates.json',candidates)
    target_stats=app.summarize(candidates,target_samples,design,config,[])
    with (design/'crRNA_designs.tsv').open(encoding='utf-8-sig') as handle:guide_rows=list(csv.DictReader(handle,delimiter='\t'))
    query=out/'guide_queries.fasta'
    app.write_fasta(query,((r['id'],r['spacer_dna_5to3']) for r in guide_rows))
    known=[r for r in records if r['role'].startswith('background_')]
    unresolved=[r for r in records if r['role']=='unresolved_Naegleria']
    counters={}
    for name,cohort in [('known_background',known),('unresolved_Naegleria',unresolved)]:
        work=out/name;work.mkdir(exist_ok=True)
        fasta=work/'18S.fasta'
        app.write_fasta(fasta,((r['id'],r['sequence']) for r in cohort))
        if guide_rows and cohort:
            cmd=app.linux_command(config['linux_env_bin']+'/python',[app.linux_path(Path(__file__).with_name('pam_screen.py')),
                        '--background',app.linux_path(fasta),'--queries',app.linux_path(query),'--index',app.linux_path(work/'pam.bin'),
                        '--output',app.linux_path(work/'matches.tsv'),'--summary',app.linux_path(work/'summary.json'),'--mismatches','3'],config)
            app.external(cmd,work/'screen.log')
            counters[name]=json.loads((work/'summary.json').read_text())['queries']
        else:counters[name]={}
    named_naegleria=sorted({r['organism'] for r in known if r['role']=='background_Naegleria_named'})
    partial_targets=[r for r in records if r['role']=='target_N_gruberi' and r['id'] not in {v['id'] for v in targets}]
    final=[]
    coverage_rows=[]
    for row in guide_rows:
        guide=row['spacer_dna_5to3']
        sites=app.exact_sites(refdna,guide)
        if not sites:continue
        site=sites[0]
        span=(site['start']-5,site['end']) if site['strand']=='+' else (site['start']-1,site['end']+4)
        supported=0;discordant=0
        for partial in partial_targets:
            if sequence_covered(partial,span,hsp_intervals):
                if app.exact_sites(partial['sequence'],guide):supported+=1
                else:discordant+=1
        covered_species={r['organism'] for r in known if r['role']=='background_Naegleria_named' and sequence_covered(r,span,hsp_intervals)}
        known_hits=counters['known_background'].get(row['id'],{}).get('near_hits',0)
        unknown_hits=counters['unresolved_Naegleria'].get(row['id'],{}).get('near_hits',0)
        row.update({'reference_accession':reference['accession'],'reference_start_1based':site['start'],
                    'reference_end_1based_inclusive':site['end'],'reference_strand':site['strand'],
                    'reference_PAM':site['pam'],'known_background_hits_le3':known_hits,'unresolved_Naegleria_hits_le3':unknown_hits,
                    'target_partial_support':supported,'target_partial_discordance':discordant,
                    'named_Naegleria_labels_with_site_coverage':len(covered_species),
                    'named_Naegleria_background_labels':len(named_naegleria),
                    'named_labels_missing_site_coverage':';'.join(s for s in named_naegleria if s not in covered_species),
                    'passes_observed_data':row['passes_target_filters']=='True' and known_hits==0 and discordant==0,
                    'whole_genus_specificity_established':False})
        final.append(row)
        for species in named_naegleria:
            coverage_rows.append({'guide':row['id'],'organism_label':species,'has_aligned_DNA_site_coverage':species in covered_species})
    app.write_tsv(out/'18S_crRNA_results.tsv',final,list(final[0]) if final else ['id'])
    app.write_tsv(out/'per_guide_species_coverage.tsv',coverage_rows,['guide','organism_label','has_aligned_DNA_site_coverage'])
    selected=[r for r in final if r['passes_observed_data']]
    strict=[r for r in selected if r['unresolved_Naegleria_hits_le3']==0]
    app.write_fasta(out/'Ngruberi_18S_spacers_candidates.fasta',((r['id'],r['spacer_rna_5to3']) for r in selected))
    app.write_fasta(out/'Ngruberi_18S_full_RNA_CaSilico_default.fasta',((r['id'],r['full_crrna_rna_5to3']) for r in selected))
    app.write_fasta(out/'Ngruberi_18S_spacers_no_unresolved_hits.fasta',((r['id'],r['spacer_rna_5to3']) for r in strict))
    summary={'target':reference['organism'],'reference':reference['accession'],'reference_18S_length':len(refdna),'target_priority_sequences':len(targets),
             'min_input_bp':MIN_18S_BP,'max_input_bp':MAX_18S_BP,'input_records':len(records),'excluded_by_length':len(all_records)-len(records),
             'other_target_sequences':len(partial_targets),'known_background_sequences':len(known),
             'unresolved_Naegleria_sequences':len(unresolved),'genus_conservation_cohort_sequences':len(full_genus),
             'genus_conservation_cohort_labels':len({r['organism'] for r in full_genus}),
             'genus_conserved_runs_ge20':len(genus_regions),'target_conserved_runs_ge20':len(target_regions),
             'CaSilico_candidates':target_stats['candidates'],'observed_data_RNA_candidates':len(selected),
             'candidates_without_unresolved_hits':len(strict),'all_described_species_specificity_proven':False,
             'background_conditions':'TTTV+20nt; exhaustive <=3 substitutions; no indels; partial records do not establish uncovered-site negatives'}
    app.write_json(out/'status.json',{'status':'completed',**summary})
    body=''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in
                 [r['id'],r['spacer_rna_5to3'],r['reference_start_1based'],r['reference_end_1based_inclusive'],
                  r['reference_strand'],r['known_background_hits_le3'],r['unresolved_Naegleria_hits_le3'],
                  f"{r['named_Naegleria_labels_with_site_coverage']}/{r['named_Naegleria_background_labels']}"])+'</tr>' for r in selected)
    page=f'''<!doctype html><meta charset="utf-8"><title>N. gruberi 18S 比较</title>
    <style>body{{font:16px system-ui;margin:36px;line-height:1.7}}table{{border-collapse:collapse}}td,th{{border:1px solid #ccc;padding:7px}}.note{{background:#fff2cf;padding:14px}}</style>
    <h1>仅限18S rDNA：全属比较，落脚 N. gruberi</h1>
    <p>参考 {reference['accession']}，提取18S DNA {len(refdna)} bp。优先目标不同序列 {len(targets)} 条，其他目标片段 {len(partial_targets)} 条；已知背景序列 {len(known)} 条，未定种奈格里序列 {len(unresolved)} 条。</p>
    <p>800–4000 bp（含边界）、明确DNA且无模糊碱基的奈格里比较集：{len(full_genus)} 条序列、{len({r['organism'] for r in full_genus})} 个物种标签。该比较集包含范围内partial记录，仅在每条序列均有比对支撑的位置评价严格保守。严格相同的 ≥20 bp 连续区：全比较集 {len(genus_regions)} 段，目标优先序列 {len(target_regions)} 段。范围外记录不参与比较、脱靶或设计。</p>
    <p>CaSilico 生成 {target_stats['candidates']} 条不同候选；通过当前已观察目标变异和已知背景筛查：{len(selected)} 条；其中未定种奈格里也无阈值内匹配：{len(strict)} 条。</p>
    <p class="note">这是当前公开序列范围内的计算结果，不代表覆盖所有已认可种或全部 N. gruberi 株系。片段没有覆盖位点时记为未评估，绝不当阴性。未定种奈格里可能包含目标种，匹配结果单列。每条候选均给出同源位点实际覆盖的背景物种标签数。</p>
    <p>背景筛查条件：明确碱基、两链、TTTV+20 nt、至多3个替换；不包含插入缺失或非典型PAM。序列比对仅用于位置覆盖与保守性描述，不替代实际DNA匹配检查。DNA内含子保留；完整RNA使用CaSilico默认DR，须核对所用Cas12a兼容性。</p>
    <p><a href="18S_crRNA_results.tsv">全部RNA与位点覆盖</a> · <a href="Ngruberi_18S_full_RNA_CaSilico_default.fasta">完整候选RNA</a> · <a href="Ngruberi_observed_conserved_18S.tsv">目标观察到的保守区</a> · <a href="Naegleria_observed_conserved_18S.tsv">奈格里比较集保守区</a> · <a href="../18S数据报告.html">18S数据核对</a></p>
    <table><tr><th>编号</th><th>spacer RNA</th><th>参考起点</th><th>终点</th><th>链</th><th>已知背景匹配</th><th>未定种匹配</th><th>具位点覆盖的奈格里标签数</th></tr>{body}</table>'''
    (out/'18S结果.html').write_text(page,encoding='utf-8')
    app.log(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()

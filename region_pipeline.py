"""Target/background nucleotide-region discovery, followed by actual CaSilico.

Region specificity is always conditional on the retrieved assemblies and stated
alignment thresholds. One target assembly cannot establish species conservation.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import html
import json
import os
from pathlib import Path
import re
import shutil
import sys
from urllib.parse import urlencode
import zipfile

import amoebacrisp as app


def unique_records(records):
    by_pair = {}
    for record in records:
        ids = tuple(sorted([record['accession'], record.get('paired_accession', record['accession'])]))
        if ids not in by_pair or record['accession'].startswith('GCF_'):
            by_pair[ids] = record
    return sorted(by_pair.values(), key=lambda r: r['accession'])


def same_species(record, target_name, target_taxids, target_accessions):
    organism = record['organism']
    name = organism['organism_name'].casefold()
    target = target_name.casefold().strip()
    return (organism['tax_id'] in target_taxids or record['accession'] in target_accessions
            or record.get('paired_accession') in target_accessions
            or name == target or name.startswith(target + ' '))


def genome_only(record, cache):
    accession = record['accession']
    folder = cache / 'ncbi' / accession
    genome = folder / 'genome.fna'
    full_package_stamp = folder / 'dataset.zip.source.json'
    if not genome.exists() or not ((folder / 'genome_only.zip.source.json').exists() or full_package_stamp.exists()):
        url = f"{app.NCBI}/genome/accession/{accession}/download?" + urlencode({'include_annotation_type': 'GENOME_FASTA'})
        bundle = app.fetch(url, folder / 'genome_only.zip')
        with zipfile.ZipFile(bundle) as archive:
            names = [n for n in archive.namelist() if n.endswith('_genomic.fna') and 'cds_from' not in n and 'rna_from' not in n]
            if len(names) != 1:
                raise RuntimeError(f'{accession} 下载包中没有唯一的基因组 FASTA')
            with archive.open(names[0]) as source, genome.open('wb') as output:
                shutil.copyfileobj(source, output)
    if genome.stat().st_size == 0:
        raise RuntimeError(f'{accession} 的基因组为空')
    return {'id': accession, 'name': record['organism']['organism_name'],
            'taxid': record['organism']['tax_id'], 'genome': str(genome),
            'source': 'ncbi', 'metadata': record, 'sha256': app.checksum(genome)}


def download_ncbi(target_name, genera, cache, run, inventory=None):
    names = [target_name] + genera
    if inventory:
        query_records = json.loads(Path(inventory).read_text(encoding='utf-8'))['queries']
        if any(name not in query_records for name in names):
            raise ValueError('已有清单缺少所请求的目标或背景属')
    else:
        with ThreadPoolExecutor(max_workers=3) as executor:
            query_records = dict(zip(names, executor.map(lambda name: app.ncbi_records(name, 0), names)))
    app.write_json(run / 'database_inventory.json', query_records)
    targets = unique_records(query_records[target_name])
    taxids = {r['organism']['tax_id'] for r in targets}
    accessions = {v for r in targets for v in (r['accession'], r.get('paired_accession', r['accession']))}
    background = unique_records([r for genus in genera for r in query_records[genus]
                                  if not same_species(r, target_name, taxids, accessions)])
    if not targets or not background:
        raise ValueError('目标或排除目标种后的背景组装为空')
    app.log(f'目标 {len(targets)} 个组装；背景 {len(background)} 个组装；全部下载，不截断样本数。')
    downloaded = {}
    with ThreadPoolExecutor(max_workers=3) as executor:
        pending = {executor.submit(genome_only, r, cache): r for r in targets + background}
        for index, future in enumerate(as_completed(pending), 1):
            record = pending[future]
            sample = future.result()  # Fail explicitly if any declared assembly is unavailable.
            downloaded[record['accession']] = sample
            app.log(f"基因组 {index}/{len(pending)}：{sample['id']} {sample['name']}")
    return [downloaded[r['accession']] for r in targets], [downloaded[r['accession']] for r in background]


def merge_intervals(intervals, length=None):
    """Merge zero-based half-open intervals."""
    merged = []
    for start, end in sorted((int(a),int(b)) for a,b in intervals):
        if length is not None:
            start, end = max(0, start), min(length, end)
        if start >= end:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def complement(intervals, length):
    cursor = 0
    result = []
    for start, end in merge_intervals(intervals, length):
        if cursor < start:
            result.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < length:
        result.append((cursor, length))
    return result


def intersect_intervals(left, right):
    result, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        if max(a, c) < min(b, d):
            result.append((max(a, c), min(b, d)))
        if b < d:
            i += 1
        else:
            j += 1
    return merge_intervals(result)


def write_background(samples, destination, run):
    contigs, bases, mapping = 0, 0, []
    with destination.open('w', encoding='ascii', newline='\n') as output:
        for sample in samples:
            for original, sequence in app.read_fasta(sample['genome']):
                contigs += 1
                identifier = f'BG{contigs:08d}'
                output.write('>' + identifier + '\n')
                for offset in range(0, len(sequence), 80):
                    output.write(sequence[offset:offset+80] + '\n')
                bases += len(sequence)
                mapping.append({'subject': identifier, 'sample': sample['id'], 'organism': sample['name'],
                                'original_contig': original, 'length': len(sequence)})
    app.write_tsv(run / 'background_contigs.tsv', mapping,
                  ['subject', 'sample', 'organism', 'original_contig', 'length'])
    return contigs, bases


def blast_tool(name, config):
    return config['linux_env_bin'].rstrip('/') + '/' + name


def run_tool(name, args, logfile, config):
    app.external(app.linux_command(blast_tool(name, config), args, config), logfile)


def build_database(fasta, prefix, config, log):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    run_tool('makeblastdb', ['-in', app.linux_path(fasta), '-dbtype', 'nucl',
                           '-parse_seqids', '-out', app.linux_path(prefix)], log, config)


def background_mask(reference, database, contig_count, run, args, config):
    output = run / 'background_alignments.tsv'
    queries, query_map = [], {}
    overlap = args.background_min_alignment-1
    if args.query_chunk <= overlap:
        raise ValueError('参考分块长度必须大于最小对齐长度')
    for contig,sequence in app.read_fasta(reference):
        for offset in range(0,len(sequence),args.query_chunk-overlap):
            fragment=sequence[offset:offset+args.query_chunk]
            if len(fragment)<args.background_min_alignment:
                continue
            key=f'REFQ{len(queries)+1:07d}'
            queries.append((key,fragment))
            query_map[key]={'contig':contig,'offset_0based':offset,'length':len(fragment)}
            if offset+len(fragment)==len(sequence):
                break
    query_file=run/'reference_queries.fasta'
    app.write_fasta(query_file,queries)
    app.write_json(run/'reference_query_map.json',query_map)
    # No cap on HSPs. Target-sequence limit covers every background contig.
    command = ['-query', app.linux_path(query_file), '-db', app.linux_path(database),
               '-task', args.blast_task, '-word_size', str(args.word_size), '-perc_identity', str(args.background_identity),
               '-evalue', '10', '-dust', 'yes', '-soft_masking', 'true',
               '-max_target_seqs', str(max(1000, contig_count)), '-num_threads', str(args.threads), '-mt_mode', '1',
               '-outfmt', '6 qseqid qstart qend pident length sseqid sstart send evalue bitscore',
               '-out', app.linux_path(output)]
    run_tool('blastn', command, run / 'background_blast.log', config)
    masks, accepted = {}, 0
    with output.open(encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip().split('\t')
            qid, a, b, identity, length = fields[:5]
            if int(length) >= args.background_min_alignment and float(identity) >= args.background_identity:
                start, end = sorted([int(a), int(b)])
                mapping=query_map[qid]
                masks.setdefault(mapping['contig'], []).append((mapping['offset_0based']+start-1,
                                                                mapping['offset_0based']+end))
                accepted += 1
    app.log(f'背景比对完成：{accepted} 个达到阈值的相似片段，按其目标坐标排除。')
    return masks, accepted


def dust_masks(reference, run, config):
    output = run / 'low_complexity.intervals'
    run_tool('dustmasker', ['-in', app.linux_path(reference), '-outfmt', 'interval',
                          '-out', app.linux_path(output)], run / 'dustmasker.log', config)
    result, current = {}, None
    for line in output.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            current = line[1:].split()[0]
        elif line.strip() and current:
            match = re.fullmatch(r'\s*(\d+)\s*-\s*(\d+)\s*', line)
            if not match:
                raise ValueError(f'不能解析 dustmasker 输出：{line}')
            a, b = map(int, match.groups())
            result.setdefault(current, []).append((a, b+1))
    return result


def exclude_unassessable_queries(masks,run):
    """A failed query is unassessed, never evidence of background absence."""
    log=run/'background_blast.log'
    mapping_file=run/'reference_query_map.json'
    if not log.exists() or not mapping_file.exists():
        return masks
    mapping=json.loads(mapping_file.read_text())
    bad=set()
    for line in log.read_text(encoding='utf-8',errors='replace').splitlines():
        if 'Karlin-Altschul parameters' in line or 'invalid query sequence' in line:
            bad.update(key for key in re.findall(r'REFQ\d+',line) if key in mapping)
    rows=[]
    for key in sorted(bad):
        m=mapping[key]
        a,b=m['offset_0based'],m['offset_0based']+m['length']
        masks.setdefault(m['contig'],[]).append((a,b))
        rows.append({'query':key,'contig':m['contig'],'start_1based':a+1,'end_1based_inclusive':b,
                     'reason':'BLAST_could_not_evaluate_query'})
    app.write_tsv(run/'unassessed_reference_chunks.tsv',rows,
                  ['query','contig','start_1based','end_1based_inclusive','reason'])
    if rows:app.log(f'另外排除 {len(rows)} 个无法完成有效背景比对的参考分块，避免将未评估误判为特异。')
    return masks


def conservation_intervals(reference, targets, sequences, run, args, config):
    intervals = {name: [(0, len(seq))] for name, seq in sequences.items()}
    if len(targets) == 1:
        app.log('目标种仅有 1 个独立组装：不将任何区域标为已证实的物种内保守区。')
        return intervals, 'single_reference_not_validated'
    for sample in targets[1:]:
        work = run / 'target_comparisons' / sample['id']
        work.mkdir(parents=True,exist_ok=True)
        db = work / 'database'
        build_database(sample['genome'], db, config, work / 'makeblastdb.log')
        hits = work / 'hits.tsv'
        run_tool('blastn', ['-query', app.linux_path(reference), '-db', app.linux_path(db), '-task', 'blastn',
                           '-word_size', '11', '-perc_identity', str(args.target_identity), '-dust', 'yes',
                           '-num_threads', str(args.threads), '-max_target_seqs', '1000000',
                           '-outfmt', '6 qseqid qstart qend pident length', '-out', app.linux_path(hits)],
                 work / 'blast.log', config)
        supported = {}
        for line in hits.read_text(encoding='utf-8').splitlines():
            name, start, end, identity, length = line.split('\t')
            if float(identity) >= args.target_identity and int(length) >= args.min_region:
                a, b = sorted([int(start), int(end)])
                supported.setdefault(name, []).append((a-1, b))
        intervals = {name: intersect_intervals(intervals[name], merge_intervals(supported.get(name, [])))
                     for name in intervals}
    return intervals, f'alignments_in_all_{len(targets)}_assemblies_identity>={args.target_identity}'


def discover_regions(sequences, masks, low_complexity, conserved, run, args, conservation_status):
    rows, fasta, number, blocks, block_fasta = [], [], 0, [], []
    for contig, seq in sequences.items():
        invalid = [(m.start(), m.end()) for m in re.finditer('[^ACGT]+', seq)]
        excluded = masks.get(contig, []) + low_complexity.get(contig, []) + invalid
        available = intersect_intervals(complement(excluded, len(seq)), conserved[contig])
        for start, end in available:
            if end - start < args.min_region:
                continue
            block_id=f'block_{len(blocks)+1:06d}'
            block_dna=seq[start:end]
            blocks.append({'block':block_id,'reference_contig':contig,'start_1based':start+1,
                           'end_1based_inclusive':end,'length':end-start,
                           'conservation':conservation_status if conservation_status.startswith('single_reference')
                           else 'preliminary_target_HSP_overlap_not_whole_block_validated'})
            block_fasta.append((block_id,block_dna))
            # Partition a maximal interval into windows without dropping a short tail.
            count = (end-start + args.max_region-1) // args.max_region
            for index in range(count):
                if (end-start)//count >= args.min_region:
                    a = start + index*(end-start)//count
                    b = start + (index+1)*(end-start)//count
                else:
                    # A narrow permitted length band needs overlapping windows.
                    a = start + index*((end-start)-args.max_region)//max(1,count-1)
                    b = min(end,a+args.max_region)
                if b-a < args.min_region:
                    continue
                dna = seq[a:b]
                number += 1
                identifier = f'region_{number:06d}'
                gc = round(100*(dna.count('G')+dna.count('C'))/len(dna), 2)
                rows.append({'region': identifier, 'reference_contig': contig, 'start_1based': a+1,
                             'end_1based_inclusive': b, 'length': len(dna), 'gc_percent': gc,
                             'parent_block':block_id,
                             'conservation': conservation_status,
                             'background_result': 'no_detected_alignment_at_declared_thresholds'})
                fasta.append((identifier, dna))
    app.write_tsv(run / 'candidate_regions.tsv', rows,
                  ['region', 'reference_contig', 'start_1based', 'end_1based_inclusive', 'length',
                   'gc_percent', 'parent_block', 'conservation', 'background_result'])
    app.write_fasta(run / 'candidate_regions.fasta', fasta)
    app.write_tsv(run/'maximal_regions.tsv',blocks,
                  ['block','reference_contig','start_1based','end_1based_inclusive','length','conservation'])
    app.write_fasta(run/'maximal_regions.fasta',block_fasta)
    return rows, dict(fasta)


def final_conservation(rows, dna, targets, run, args, config):
    """For >1 targets require one HSP covering a whole candidate, not a mosaic."""
    if len(targets) == 1 or not rows:
        return rows, dna
    valid = set(dna)
    for sample in targets[1:]:
        work = run / 'target_comparisons' / sample['id']
        hits = work / 'candidate_coverage.tsv'
        run_tool('blastn', ['-query', app.linux_path(run / 'candidate_regions.fasta'),
                           '-db', app.linux_path(work / 'database'), '-task', 'blastn', '-word_size', '11',
                           '-perc_identity', str(args.target_identity), '-qcov_hsp_perc', str(args.target_coverage),
                           '-num_threads', str(args.threads), '-max_target_seqs', '1000000',
                           '-outfmt', '6 qseqid pident qcovhsp', '-out', app.linux_path(hits)], work / 'candidate_coverage.log', config)
        supported = set()
        for line in hits.read_text().splitlines():
            name, identity, coverage = line.split('\t')
            if float(identity) >= args.target_identity and float(coverage) >= args.target_coverage:
                supported.add(name)
        valid &= supported
    rows = [row for row in rows if row['region'] in valid]
    dna = {key:value for key,value in dna.items() if key in valid}
    app.write_fasta(run / 'candidate_regions.fasta', dna.items())
    app.write_tsv(run / 'candidate_regions.tsv', rows, list(rows[0]) if rows else ['region'])
    return rows, dna


def guide_background_screen(candidates, background_fasta, contig_count, run, config, threads, mismatches):
    """Exact Hamming screening of every valid PAM-bearing 20 nt site."""
    guides = sorted({str(c['Guide_RNA']).upper().replace('U', 'T') for c in candidates})
    query = run / 'guide_queries.fasta'
    id_to_guide = {f'guide_{i:06d}': g for i,g in enumerate(guides, 1)}
    app.write_fasta(query, id_to_guide.items())
    if not guides:
        return {}, []
    output=run/'guide_background_sites.tsv'
    summary=run/'guide_background_summary.json'
    if summary.exists() and output.exists():
        previous=json.loads(summary.read_text())
        previous_queries={key:value['sequence'] for key,value in previous.get('queries',{}).items()}
        if (previous.get('method')=='exhaustive_2bit_Hamming' and previous.get('max_substitutions')==mismatches
                and previous_queries==id_to_guide and previous.get('background_sha256')==app.checksum(background_fasta)):
            app.log('复用背景校验值、RNA 序列及错配条件均一致的完整筛查结果。')
            return {id_to_guide[key]:value['near_hits'] for key,value in previous['queries'].items()}, []
    command=app.linux_command(config['linux_env_bin']+'/python',
              [app.linux_path(Path(__file__).with_name('pam_screen.py')),
               '--background',app.linux_path(background_fasta),'--queries',app.linux_path(query),
               '--index',app.linux_path(run/'background_pam_index.bin'),
               '--output',app.linux_path(output),'--summary',app.linux_path(summary),
               '--mismatches',str(mismatches)],config)
    app.external(command,run/'guide_background_screen.log')
    result=json.loads(summary.read_text())
    return {id_to_guide[key]:value['near_hits'] for key,value in result['queries'].items()}, []


def design_regions(rows, dna, targets, background, background_fasta, contig_count, run, args, config):
    if not args.design_regions or not rows:
        return {'status':'not_requested','regions':0,'selected':0}
    # Rank by length then GC proximity; explicitly not a claim of best assay activity.
    usable=[r for r in rows if re.search('TTT[ACG][ACGT]{20}',dna[r['region']])
            or re.search('TTT[ACG][ACGT]{20}',app.rc(dna[r['region']]))]
    selected = sorted(usable,key=lambda r:(-r['length'],abs(r['gc_percent']-50),r['region']))[:args.design_regions]
    design_dir = run / 'casilico_designs'
    design_dir.mkdir(exist_ok=True)
    app.write_json(design_dir / 'groups.json',{'eligible':len(rows),'processed':len(selected)})
    candidates=[]
    for i,row in enumerate(selected,1):
        key=row['region']
        app.log(f'候选区域 → CaSilico {i}/{len(selected)}：{key}')
        folder=design_dir/'casilico'/key
        previous=folder/'targets.fasta'
        tables=list(folder.rglob('Spacers_Information.Html')) if folder.exists() else []
        if args.resume and (folder/'completed.txt').exists() and previous.exists() and len(tables)==1 and list(app.read_fasta(previous))==[(key,dna[key])]:
            raw=app.parse_casilico(tables[0])
            app.log('复用已完成且输入序列一致的 CaSilico 结果。')
        else:
            raw=app.design(key,[key],{key:{'dna':dna[key]}},design_dir,config,args.threads)
        candidates.extend({**r,'group':key} for r in raw)
    app.write_json(design_dir / 'casilico_candidates.json',candidates)
    stats=app.summarize(candidates,targets,design_dir,config,[])
    counts,_=guide_background_screen(candidates,background_fasta,contig_count,run,config,args.threads,args.guide_mismatches)
    # Preserve the target-side output and create a separate, explicitly provisional background-screened export.
    import csv
    with (design_dir/'crRNA_designs.tsv').open(encoding='utf-8-sig') as handle:
        guide_rows=list(csv.DictReader(handle,delimiter='\t'))
    for row in guide_rows:
        row['background_pam_hits']=counts.get(row['spacer_dna_5to3'],0)
        row['background_screen_method']=f'exhaustive_Hamming<={args.guide_mismatches};20nt;TTTV;no_indels'
        row['background_genomes']=len(background)
    app.write_tsv(run/'crRNA_designs_with_background.tsv',guide_rows,list(guide_rows[0]) if guide_rows else ['id'])
    passed=[r for r in guide_rows if r['passes_target_filters']=='True' and r['background_pam_hits']==0]
    app.write_fasta(run/'crRNA_spacers_provisional.fasta',((r['id'],r['spacer_rna_5to3']) for r in passed))
    app.write_fasta(run/'crRNA_full_CaSilico_default_provisional.fasta',((r['id'],r['full_crrna_rna_5to3']) for r in passed))
    return {'status':'completed','regions':len(selected),'casilico_candidates':stats['candidates'],
            'selected':len(passed),'background_method':f'exhaustive_Hamming<={args.guide_mismatches};20nt;TTTV;no_indels'}


def report(run, args, targets, background, rows, stats, bg_stats, conservation_status):
    conserved = '尚未验证：只有一个目标参考组装' if len(targets)==1 else f'在 {len(targets)} 个目标组装中通过核查'
    listing=''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [r['region'],r['reference_contig'],r['start_1based'],r['end_1based_inclusive'],r['length'],r['gc_percent']])+'</tr>' for r in sorted(rows,key=lambda r:(-r['length'],r['region']))[:100])
    block_count=len({r.get('parent_block',r['region']) for r in rows})
    bad_file=run/'unassessed_reference_chunks.tsv'
    if bad_file.exists():
        with bad_file.open(encoding='utf-8-sig') as handle:unassessed=max(0,sum(1 for _ in handle)-1)
    else:unassessed=0
    text=f'''<!doctype html><meta charset="utf-8"><title>目标种与背景属比较</title>
    <style>body{{font:16px system-ui;margin:36px;line-height:1.7;color:#173247}}table{{border-collapse:collapse}}td,th{{border:1px solid #ccd6de;padding:7px}}.note{{background:#fff4d8;padding:16px}}</style>
    <h1>{html.escape(args.target)}：候选特异区域</h1>
    <p>背景属：{html.escape(', '.join(args.background_genus))}；从背景排除目标种。</p>
    <p>目标 {len(targets)} 个独立组装；背景 {len(background)} 个组装、{bg_stats['contigs']:,} 个 contig、{bg_stats['bases']:,} bp。</p>
    <p class="note">目标种保守性：{conserved}。特异性仅相对于本次下载的背景，并依赖所列参数；不代表实验检测性能。</p>
    <p>排除背景 BLAST 相似片段：一致性 ≥{args.background_identity}%、对齐长度 ≥{args.background_min_alignment} bp；task {args.blast_task}、word size {args.word_size}、E-value ≤10；参考分块 {args.query_chunk} bp，相邻块重叠 {args.background_min_alignment-1} bp。另排除 DUST 低复杂度区和非 A/C/G/T 碱基。</p>
    <p>无法完成有效比对的参考分块：{unassessed} 个，均不作为特异性证据，已从候选范围排除。</p>
    <p>候选窗口：{len(rows):,} 个，长度 {args.min_region}–{args.max_region} bp，来自 {block_count:,} 个连续候选片段；坐标为参考基因组正链、1-based、两端包含。比对为启发式检索，“未检出”不等于证明不存在相似序列。</p>
    <p>CaSilico：设计 {stats.get('regions',0)} 个区域，候选 {stats.get('casilico_candidates',0)} 条；背景筛查后暂留 {stats.get('selected',0)} 条。RNA 背景检查穷举两条链所有明确的 TTTV+20 nt 位点，允许至多 {args.guide_mismatches} 个替换，不包含插入缺失、非典型 PAM 或模糊碱基。RNA 文件为未经实验验证的计算候选。</p>
    <p><a href="candidate_regions.tsv">全部窗口与坐标</a> · <a href="candidate_regions.fasta">窗口 FASTA</a> · <a href="maximal_regions.tsv">连续片段坐标</a> · <a href="maximal_regions.fasta">连续片段 FASTA</a> · <a href="samples.tsv">样本清单</a> · <a href="crRNA_designs_with_background.tsv">RNA 结果表</a> · <a href="crRNA_spacers_provisional.fasta">RNA spacer</a> · <a href="crRNA_full_CaSilico_default_provisional.fasta">含默认 DR 的完整 RNA</a></p>
    <h2>前 100 个区域（按长度排列）</h2><table><tr><th>区域</th><th>参考 contig</th><th>起点</th><th>终点</th><th>长度</th><th>GC%</th></tr>{listing}</table>'''
    (run/'结果.html').write_text(text,encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description='目标种共有区域减去背景属相似区域，最后调用 CaSilico')
    parser.add_argument('--target',default='Naegleria gruberi')
    parser.add_argument('--background-genus',action='append',required=True)
    parser.add_argument('--database',choices=['ncbi','amoebadb'],default='ncbi')
    parser.add_argument('--background-fasta',action='append',default=[],help='额外本地背景基因组')
    parser.add_argument('--inventory',help='复用本次已检索的 NCBI 清单')
    parser.add_argument('--config',default=str(Path(__file__).resolve().parent/'config.json'))
    parser.add_argument('--threads',type=int,default=8)
    parser.add_argument('--background-identity',type=float,default=80)
    parser.add_argument('--background-min-alignment',type=int,default=50)
    parser.add_argument('--blast-task',choices=['blastn','megablast','dc-megablast'],default='megablast')
    parser.add_argument('--word-size',type=int,default=16)
    parser.add_argument('--query-chunk',type=int,default=10000)
    parser.add_argument('--target-identity',type=float,default=95)
    parser.add_argument('--target-coverage',type=float,default=95)
    parser.add_argument('--min-region',type=int,default=200)
    parser.add_argument('--max-region',type=int,default=1000)
    parser.add_argument('--design-regions',type=int,default=10,help='交给 CaSilico 的区域数；0=仅输出区域')
    parser.add_argument('--guide-mismatches',type=int,default=3,choices=range(5))
    parser.add_argument('--resume',help='恢复已有本次任务目录，复用下载和已经完成的比对')
    args=parser.parse_args()
    if not 50<=args.min_region<=args.max_region or args.threads<1 or args.design_regions<0:
        parser.error('区域长度、线程或设计区域数无效')
    config=app.default_config()
    if Path(args.config).exists():config.update(json.loads(Path(args.config).read_text(encoding='utf-8-sig')))
    storage=Path(config['storage'])
    run=Path(args.resume) if args.resume else storage/'region_analysis'/(dt.datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+app.slug(args.target))
    run.mkdir(parents=True,exist_ok=True)
    app.log(f'输出目录：{run}')
    request={**vars(args),'config':config}
    if args.resume and (run/'request.json').exists():
        old=json.loads((run/'request.json').read_text(encoding='utf-8'))
        critical=['target','database','background_genus','background_fasta','background_identity','background_min_alignment','target_identity','target_coverage','min_region','max_region']
        if (run/'background_mask.json').exists():critical+=['blast_task','word_size','query_chunk']
        if any(old.get(k, [] if k=='background_fasta' else None)!=request.get(k) for k in critical):
            raise ValueError('恢复任务时不能更改筛选条件，请新建任务')
    app.write_json(run/'request.json',request)
    app.write_json(run/'status.json',{'status':'running','stage':'download'})
    try:
        manifests=run/'samples.json'
        if args.resume and manifests.exists():
            data=json.loads(manifests.read_text());targets,background=data['targets'],data['background']
        elif args.database=='ncbi':
            targets,background=download_ncbi(args.target,args.background_genus,storage/'cache',run,args.inventory)
        else:
            targets=app.amoebadb_download(args.target,storage/'cache',0,config['amoebadb_base'])
            target_ids={s['id'] for s in targets}
            downloaded=[s for genus in args.background_genus for s in app.amoebadb_download(genus,storage/'cache',0,config['amoebadb_base'])]
            background=list({s['id']:s for s in downloaded if s['id'] not in target_ids}.values())
            if not background:raise ValueError('背景为空')
        if not (args.resume and manifests.exists()):
            target_hashes={app.checksum(s['genome']) for s in targets}
            for index,filename in enumerate(args.background_fasta,1):
                digest=app.checksum(filename)
                if digest in target_hashes:raise ValueError('额外背景 FASTA 与目标基因组相同，不能作为排除背景')
                background.append({'id':f'LOCAL_{index:03d}_{app.slug(Path(filename).stem)}',
                                   'name':Path(filename).name,'genome':str(Path(filename).resolve()),
                                   'sha256':digest,'source':'local_background'})
        app.write_json(manifests,{'targets':targets,'background':background})
        sample_rows=[{'role':role,'id':s['id'],'organism':s['name'],'genome':s['genome'],'sha256':s.get('sha256') or app.checksum(s['genome'])}
                     for role,samples in [('target',targets),('background',background)] for s in samples]
        app.write_tsv(run/'samples.tsv',sample_rows,['role','id','organism','genome','sha256'])
        reference=targets[0]['genome']
        sequences=dict(app.read_fasta(reference))
        background_fasta=run/'background.fasta'
        bg_stats_file=run/'background_stats.json'
        if args.resume and background_fasta.exists() and bg_stats_file.exists():
            bg_stats=json.loads(bg_stats_file.read_text())
        else:
            count,bases=write_background(background,background_fasta,run)
            bg_stats={'contigs':count,'bases':bases}
            app.write_json(bg_stats_file,bg_stats)
        db=run/'background_db'/'database'
        db_stamp=run/'background_db'/'completed.json'
        if not (args.resume and db_stamp.exists()):
            build_database(background_fasta,db,config,run/'makeblastdb.log')
            app.write_json(db_stamp,{'completed':True})
        app.write_json(run/'status.json',{'status':'running','stage':'background_alignment','targets':len(targets),'background':len(background)})
        mask_file=run/'background_mask.json'
        if args.resume and mask_file.exists():
            data=json.loads(mask_file.read_text());masks,hit_count=data['masks'],data['hits']
        else:
            masks,hit_count=background_mask(reference,db,bg_stats['contigs'],run,args,config)
            app.write_json(mask_file,{'masks':masks,'hits':hit_count})
        masks=exclude_unassessable_queries(masks,run)
        low_complexity=dust_masks(reference,run,config)
        conserved,conservation_status=conservation_intervals(reference,targets,sequences,run,args,config)
        rows,dna=discover_regions(sequences,masks,low_complexity,conserved,run,args,conservation_status)
        rows,dna=final_conservation(rows,dna,targets,run,args,config)
        app.log(f'得到 {len(rows)} 个候选区域；保守性状态：{conservation_status}')
        app.write_json(run/'status.json',{'status':'running','stage':'casilico','regions':len(rows)})
        stats=design_regions(rows,dna,targets,background,background_fasta,bg_stats['contigs'],run,args,config)
        report(run,args,targets,background,rows,stats,bg_stats,conservation_status)
        app.write_json(run/'status.json',{'status':'completed','targets':len(targets),'background':len(background),
                       'background_bases':bg_stats['bases'],'candidate_regions':len(rows),'conservation':conservation_status,
                       'maximal_blocks':len({r.get('parent_block',r['region']) for r in rows}),
                       'background_alignment_hits':hit_count,'casilico':stats})
        app.log(f"完成：{run/'结果.html'}")
    except Exception as exc:
        app.write_json(run/'status.json',{'status':'failed','error':str(exc)})
        raise


if __name__=='__main__':
    main()

"""Four clean 18S outputs. All-species claims require site coverage."""
from runtime import storage_path, settings
from array import array
from collections import Counter
import json
from pathlib import Path
import re
import amoebacrisp as app
from ssu18s_policy import filter_records,MIN_18S_BP,MAX_18S_BP,gc_percent

GROUPS=('保守序列_有PAM','保守序列_无PAM','特异序列_有PAM','特异序列_无PAM')


def mismatch_count(a,b,limit=None):
    count=0
    for left,right in zip(a,b):
        if left!=right:
            count+=1
            if limit is not None and count>limit:break
    return count


def near_hits(queries,cohort,radius=2):
    """Pigeonhole seeds find all substitution-only matches, with or without PAM."""
    if not queries:return Counter()
    k=len(next(iter(queries)))
    if not 0<=radius<k:raise ValueError('Mismatch radius must be smaller than window length')
    boundaries=[i*k//(radius+1) for i in range(radius+2)]
    indexes=[]
    for a,b in zip(boundaries,boundaries[1:]):
        index={}
        for query in queries:index.setdefault(query[a:b],set()).add(query)
        indexes.append((a,b,index))
    result=Counter()
    for row in cohort:
        found=set()
        for sequence in (row['sequence'],app.rc(row['sequence'])):
            for i in range(len(sequence)-k+1):
                word=sequence[i:i+k]
                if not re.fullmatch('[ACGT]+',word):continue
                candidates=set()
                for a,b,index in indexes:candidates.update(index.get(word[a:b],()))
                for query in candidates-found:
                    if mismatch_count(query,word,radius)<=radius:found.add(query)
        result.update(found)
    return result


def species_matches(queries,cohort,radius):
    """Return per-species record counts and the nearest observed match."""
    if not queries:return {}
    k=len(next(iter(queries)))
    boundaries=[i*k//(radius+1) for i in range(radius+2)]
    indexes=[]
    for left,right in zip(boundaries,boundaries[1:]):
        seeds={}
        for query in queries:seeds.setdefault(query[left:right],set()).add(query)
        indexes.append((left,right,seeds))
    result={}
    for row in cohort:
        record_matches={}
        for strand,sequence in (('+',row['sequence']),('-',app.rc(row['sequence']))):
            for start in range(len(sequence)-k+1):
                word=sequence[start:start+k]
                if not re.fullmatch('[ACGT]+',word):continue
                candidates=set()
                for left,right,seeds in indexes:candidates.update(seeds.get(word[left:right],()))
                for query in candidates:
                    previous=record_matches.get(query)
                    if previous and previous['min_mismatches']==0:continue
                    distance=mismatch_count(query,word,radius)
                    if distance>radius or (previous and distance>=previous['min_mismatches']):continue
                    begin,end=(start+1,start+k) if strand=='+' else (len(sequence)-start-k+1,len(sequence)-start)
                    record_matches[query]={'min_mismatches':distance,'example_accession':row.get('accession',row['id']),
                        'example_start':begin,'example_end':end,'example_strand':strand,'example_sequence':word}
        for query,match in record_matches.items():
            species=result.setdefault(query,{})
            previous=species.get(row['organism'])
            count=previous['records_matching']+1 if previous else 1
            if previous and previous['min_mismatches']<=match['min_mismatches']:match=previous.copy()
            species[row['organism']]={**match,'records_matching':count}
    return result


def unresolved_name(name):
    return bool(re.search(r'\b(?:sp|cf|aff)\.?\b|\buncultured\b',name,re.I))


def pam_options(reference,start,k=20):
    """Return oriented upstream TTTVs; coordinates still refer to plus DNA."""
    found=[]
    if start>=4 and re.fullmatch('TTT[ACG]',reference[start-4:start]):
        found.append(('+',reference[start-4:start],reference[start:start+k]))
    if start+k+4<=len(reference):
        pam=app.rc(reference[start+k:start+k+4])
        if re.fullmatch('TTT[ACG]',pam):found.append(('-',pam,app.rc(reference[start:start+k])))
    return found


def nearest_background(queries,cohort):
    """Exact minimum substitutions over all unambiguous windows on both strands."""
    import numpy as np
    if not queries:return {}
    words=sorted(queries);k=len(words[0])
    encoded=np.array([list(word.encode('ascii')) for word in words],dtype=np.uint8)
    result={word:{'mismatches':k+1,'accession':'','organism':''} for word in words}
    for row in cohort:
        for sequence in (row['sequence'],app.rc(row['sequence'])):
            if len(sequence)<k:continue
            windows=np.lib.stride_tricks.sliding_window_view(np.frombuffer(sequence.encode('ascii'),dtype=np.uint8),k)
            windows=windows[np.all(np.isin(windows,[65,67,71,84]),axis=1)]
            if not len(windows):continue
            for i,word in enumerate(words):
                if result[word]['mismatches']==0:continue
                minimum=int(np.count_nonzero(windows!=encoded[i],axis=1).min())
                if minimum<result[word]['mismatches']:
                    result[word]={'mismatches':minimum,'accession':row.get('accession',row['id']),'organism':row['organism']}
    return result


def projections(records,alignment,reference):
    """Highest scoring reliable HSP, preserving physical DNA adjacency."""
    n=len(reference['sequence']); byid={r['id']:r for r in records}
    bases={i:['N']*n for i in byid}
    positions={i:array('i',[-1])*n for i in byid}
    scores={i:array('d',[-1])*n for i in byid}
    with Path(alignment).open() as handle:
        for line in handle:
            _,sid,a,b,ss,se,identity,qseq,sseq,score=line.rstrip().split('\t')
            if sid not in byid:continue
            a,b=int(a)-1,int(b)
            if float(identity)<80 or b-a<50:continue
            score=float(score); q=a; p=int(ss)-1; step=1 if int(ss)<=int(se) else -1
            for qb,sb in zip(qseq.upper(),sseq.upper()):
                if qb=='-':
                    if sb!='-':p+=step
                    continue
                if score>scores[sid][q]:
                    bases[sid][q]=sb if sb in 'ACGT' else 'N'
                    positions[sid][q]=p if sb in 'ACGT' else -1
                    scores[sid][q]=score
                elif score==scores[sid][q] and (positions[sid][q]!=p or bases[sid][q]!=sb):
                    bases[sid][q]='N';positions[sid][q]=-1
                if sb!='-':p+=step
                q+=1
    # The reference is directly observed DNA, not inferred from a local HSP.
    bases[reference['id']]=list(reference['sequence'])
    positions[reference['id']]=array('i',range(n))
    # Very short genuine 18S fragments cannot meet a 50bp HSP threshold.
    # Use only a unique, exact full-fragment anchor, never an arbitrary location.
    for row in records:
        seq=row['sequence']
        if not 20<=len(seq)<50 or not re.fullmatch('[ACGT]+',seq):continue
        anchors=[]
        for step,word in ((1,seq),(-1,app.rc(seq))):
            at=reference['sequence'].find(word)
            while at>=0:
                anchors.append((at,step));at=reference['sequence'].find(word,at+1)
        if len(anchors)==1:
            at,step=anchors[0];word=seq if step==1 else app.rc(seq)
            for offset,base in enumerate(word):
                bases[row['id']][at+offset]=base
                positions[row['id']][at+offset]=offset if step==1 else len(seq)-1-offset
    return bases,positions


def contiguous_windows(bases,positions,k=20):
    """Never bridge an insertion, deletion, ambiguous base or alignment gap."""
    start=0; previous_step=None
    for i,base in enumerate(bases):
        if base not in 'ACGT' or positions[i]<0:
            start=i+1;previous_step=None;continue
        if i>start:
            step=positions[i]-positions[i-1]
            if abs(step)!=1 or (previous_step is not None and previous_step!=step):
                start=i;previous_step=None
            else:previous_step=step
        if i-start+1>=k:yield i-k+1,''.join(bases[i-k+1:i+1])


def analyze(root,conserved_mismatches=2,target_mismatches=2,background_mismatches=1,window_length=20,write_outputs=True,
            coverage_min=None,coverage_max=None,alignment_method='mafft',excluded_species=None,max_exception_species=0):
    if any(not 0<=value<=4 for value in (conserved_mismatches,target_mismatches)) or not 0<=background_mismatches<window_length:
        raise ValueError('目标错配需为0–4，背景排除错配数需小于窗口长度')
    if not isinstance(max_exception_species,int) or max_exception_species<0:raise ValueError('允许不特异物种数需为非负整数')
    root=Path(root);out=root/'comparison';k=window_length
    all_records=json.loads((root/'curated_records.json').read_text())
    length_records=filter_records(all_records)
    if excluded_species:
        raise ValueError('请保留全部背景，并使用每段允许不特异物种数设置例外上限')
    excluded=[]
    records=length_records
    cohorts=json.loads((out/'comparison_cohorts.json').read_text())
    byid={r['id']:r for r in records}; ref=byid[cohorts['reference']];dna=ref['sequence'];n=len(dna)-k+1
    valid=[r for r in records if re.search('[ACGT]{20}',r['sequence'])]
    priority=set(cohorts['target_priority']) & set(byid)
    if not priority:raise ValueError('800–4000bp范围内没有优先目标序列')
    genus=[r for r in valid if r['role'] in ('target_N_gruberi','background_Naegleria_named')]
    genus_labels={r['organism'] for r in genus}
    # Near-full 18S supplies a common comparison backbone. Short fragments add
    # evidence wherever they cover, but never veto an unrelated missing site.
    backbone_labels=genus_labels
    minimum=len(genus_labels) if coverage_min is None else coverage_min
    maximum=len(genus_labels) if coverage_max is None else coverage_max
    if not 1<=minimum<=maximum<=len(genus_labels):raise ValueError(f'保守物种覆盖区间需在1–{len(genus_labels)}内')
    background=[r for r in records if r['role'].startswith('background_') and not unresolved_name(r['organism'])]
    background_labels={r['organism'] for r in background}
    unresolved=[r for r in records if r['role']=='unresolved_Naegleria' or (r['role'].startswith('background_') and unresolved_name(r['organism']))]
    known_labels_by_start=[set() for _ in range(n)]
    genus_labels_by_start=[set() for _ in range(n)]
    genus_conflict=[False]*n;target_conflict=[False]*n;target_support=[set() for _ in range(n)]
    genus_max_mismatches=[0]*n;target_max_mismatches=[0]*n
    species_mismatches=[{} for _ in range(n)]
    if alignment_method=='mafft':
        from ssu18s_alignment import prepare,project
        bases,positions=project(valid,prepare(root),ref)
    elif alignment_method=='blast':bases,positions=projections(valid,out/'reference_alignments.tsv',ref)
    else:raise ValueError('Unknown alignment method')
    for r in valid:
        for start,word in contiguous_windows(bases[r['id']],positions[r['id']],k):
            differences=mismatch_count(word,dna[start:start+k])
            if r['role']=='target_N_gruberi':
                target_max_mismatches[start]=max(target_max_mismatches[start],differences)
                if differences>target_mismatches:target_conflict[start]=True
                if r['id'] in priority and differences<=target_mismatches:target_support[start].add(r['id'])
            if r['role'] in ('target_N_gruberi','background_Naegleria_named'):
                genus_labels_by_start[start].add(r['organism'])
                genus_max_mismatches[start]=max(genus_max_mismatches[start],differences)
                species_mismatches[start][r['organism']]=max(species_mismatches[start].get(r['organism'],0),differences)
                if differences>conserved_mismatches:genus_conflict[start]=True
            if r['role'].startswith('background_'):known_labels_by_start[start].add(r['organism'])
    target_ok=[len(target_support[i])==len(priority) and not target_conflict[i] for i in range(n)]
    priority_ok=[len(target_support[i])==len(priority) for i in range(n)]
    wanted={dna[i:i+k] for i in range(n) if priority_ok[i]}
    def exact_hits(cohort):
        hits=Counter()
        for r in cohort:
            seq=r['sequence']; found=set()
            for strand in (seq,app.rc(seq)):
                found.update(strand[i:i+k] for i in range(len(strand)-k+1) if strand[i:i+k] in wanted)
            hits.update(found)
        return hits
    known=exact_hits(background);unknown=exact_hits(unresolved)
    diagnostic_queries={dna[i:i+k] for i in range(n) if target_ok[i] and not known[dna[i:i+k]]}
    if write_outputs:
        closest_known=nearest_background(diagnostic_queries,background)
        closest_unknown=nearest_background(diagnostic_queries,unresolved)
        diagnostics=[{'sequence':word,'known_min_mismatches':closest_known[word]['mismatches'],
                      'known_accession':closest_known[word]['accession'],'known_species':closest_known[word]['organism'],
                      'unresolved_min_mismatches':closest_unknown[word]['mismatches'] if closest_unknown else None,
                      'unresolved_accession':closest_unknown[word]['accession'] if closest_unknown else ''} for word in sorted(diagnostic_queries)]
        app.write_json(out/'specificity_diagnostics.json',{'minimum_required':background_mismatches+1,'windows':diagnostics,
                       'scope':'All unambiguous contiguous windows, both strands, no PAM requirement. Unknown bases and indels not counted as substitutions.'})
    eligible={dna[i:i+k] for i in range(n) if target_ok[i]}
    print(f'无PAM限制的背景近似匹配筛查：{len(eligible)}个窗口，最多{background_mismatches}个替换',flush=True)
    known_matches=species_matches(eligible,background,background_mismatches)
    unresolved_matches=species_matches(eligible,unresolved,background_mismatches)
    known_near={query:sum(r['records_matching'] for r in matches.values()) for query,matches in known_matches.items()}
    unknown_near={query:sum(r['records_matching'] for r in matches.values()) for query,matches in unresolved_matches.items()}
    rows={name:[] for name in GROUPS};audit=Counter()
    for i in range(n):
        if priority_ok[i]:
            audit['priority_conserved_windows']+=1
            if not known[dna[i:i+k]]:
                audit['priority_no_known_exact_match_windows']+=1
                if target_conflict[i]:audit['priority_unique_excluded_by_target_partial_conflict']+=1
        if not target_ok[i]:continue
        word=dna[i:i+k];audit['target_conserved_windows']+=1
        matching_species={label for label,distance in species_mismatches[i].items() if distance<=conserved_mismatches}
        conserved=ref['organism'] in matching_species and minimum<=len(matching_species)<=maximum
        # DNA specificity is screened over the actual eligible DNA on both
        # strands. It does not require BLAST/MAFFT to align a homologous site.
        exceptions=known_matches.get(word,{})
        unresolved_hits=unresolved_matches.get(word,{})
        assessed_background=known_labels_by_start[i] & background_labels
        distinguishable=assessed_background-set(exceptions)
        unassessed=background_labels-assessed_background-set(exceptions)
        specific=bool(distinguishable) and len(exceptions)<=max_exception_species
        audit[f'windows_with_{len(exceptions)}_exception_species']+=1
        if max_exception_species==0 and unresolved_hits:specific=False
        if not known[word]:audit['no_known_exact_match_windows']+=1
        if not known[word] and unknown[word]:audit['unresolved_match_windows']+=1
        if known_near.get(word,0) or unknown_near.get(word,0):audit['near_background_match_windows']+=1
        if not known[word] and known_labels_by_start[i]!=background_labels:audit['background_coverage_incomplete_windows']+=1
        if not conserved:audit['not_all_genus_labels_conserved_or_covered']+=1
        if genus_labels_by_start[i]!=genus_labels:audit['genus_label_coverage_incomplete_windows']+=1
        if genus_conflict[i]:audit['genus_discordance_windows']+=1
        options=pam_options(dna,i,k)
        for kind,passed in [('保守序列',conserved),('特异序列',specific)]:
            if not passed:continue
            name=kind+('_有PAM' if options else '_无PAM')
            for strand,pam,oriented in options or [('+','',word)]:
                rows[name].append({'id':f'{name}_{len(rows[name])+1:04d}','reference':ref['accession'],
                    'start':i+1,'end':i+k,'strand':strand,'sequence':oriented,'pam':pam,
                    'gc_percent':gc_percent(oriented),
                    'target_support':len(priority),'compared_species':len(matching_species) if kind=='保守序列' else len(background_labels),
                    'assessed_species':len(genus_labels_by_start[i]) if kind=='保守序列' else len(known_labels_by_start[i]),
                    'matching_species':'; '.join(sorted(matching_species)) if kind=='保守序列' else ref['organism'],
                    'total_species':len(genus_labels) if kind=='保守序列' else len(background_labels),
                    'max_mismatches':max(species_mismatches[i][s] for s in matching_species) if kind=='保守序列' else target_max_mismatches[i],
                    'background_min_mismatches_lower_bound':background_mismatches+1 if kind=='特异序列' else None,
                    'exception_species_count':len(exceptions) if kind=='特异序列' else None,
                    'exception_species':'; '.join(sorted(exceptions)) if kind=='特异序列' else '',
                    'exception_matches':exceptions if kind=='特异序列' else {},
                    'distinguishable_species_count':len(distinguishable) if kind=='特异序列' else None,
                    'distinguishable_species':'; '.join(sorted(distinguishable)) if kind=='特异序列' else '',
                    'unresolved_matching_labels':'; '.join(sorted(unresolved_hits)) if kind=='特异序列' else '',
                    'unresolved_matching_records':sum(h['records_matching'] for h in unresolved_hits.values()) if kind=='特异序列' else None,
                    'specificity_status':('条件性候选；未定种命中待核查' if unresolved_hits else '允许例外的条件性候选' if exceptions else '当前已评估背景无近似匹配') if kind=='特异序列' else '',
                    'unassessed_species':'; '.join(sorted((genus_labels-genus_labels_by_start[i]) if kind=='保守序列' else unassessed)),
                    'source':ref.get('source_url','https://www.ncbi.nlm.nih.gov/nuccore/'+ref['accession'])})
    headers=['id','reference','start','end','strand','sequence','pam','gc_percent','target_support','compared_species','total_species','max_mismatches','unassessed_species','source']
    headers+=['exception_species_count','exception_species','distinguishable_species_count','distinguishable_species','unresolved_matching_labels','unresolved_matching_records','specificity_status']
    for name,items in (rows.items() if write_outputs else []):
        app.write_tsv(out/(name+'.tsv'),items,headers)
        app.write_fasta(out/(name+'.fasta'),((r['id'],r['sequence']) for r in items))
    result={'groups':rows,'counts':{name:len(items) for name,items in rows.items()},'reference':ref['accession'],
        'target':ref['organism'],'window_length':k,'genus_species':sorted(genus_labels),
        'DNA_source_records':sum(r.get('DNA_supported',False) for r in valid),
        'rRNA_or_other_source_records':sum(not r.get('DNA_supported',False) for r in valid),
        'min_input_bp':MIN_18S_BP,'max_input_bp':MAX_18S_BP,'input_records':len(records),'excluded_by_length':len(all_records)-len(length_records),
        'excluded_species':excluded,'excluded_by_species_records':len(length_records)-len(records),
        'background_species':sorted(background_labels),'valid_DNA_records':sum(r.get('DNA_supported',False) for r in valid),
        'valid_sequence_records':len(valid),'audit':dict(audit),
        'backbone_species':sorted(backbone_labels),'conserved_mismatches':conserved_mismatches,'target_mismatches':target_mismatches,'background_mismatches':background_mismatches,
        'coverage_min':minimum,'coverage_max':maximum,'alignment_method':alignment_method,
        'max_exception_species':max_exception_species,'unresolved_background_labels':sorted({r['organism'] for r in unresolved}),
        'rules':{'conserved':'Only extracted 18S 800–4000bp inclusive enter comparison. A species counts as conserved when every assessed record at this site is within tolerance to the reference. Conserved species count must fall in the inclusive user interval and include the target. Unassessed species remain missing.',
                 'specific':'Keep all background records. Each named species with any <=radius match on either strand is one per-window exception. Accept up to the configured count and require at least one covered, non-exception species. Unresolved labels are screened separately and marked for review, never counted as named species. Strict zero-exception mode also rejects unresolved hits. Missing homologous sites remain unassessed. No indels modeled.',
                 'pam':'Reference upstream TTTV on either strand; no-PAM means neither strand has TTTV for that 20nt interval. PAM conservation across other species is not required.',
                 'scope':'Downloaded annotated 18S only, not all public species. Exact DNA specificity is not a near-match crRNA safety assessment. CaSilico guide filtering remains separate.'}}
    if write_outputs:app.write_json(out/'four_groups.json',result)
    print(json.dumps({'counts':result['counts'],'genus_species':len(genus_labels),'background_labels':len(background_labels),'audit':dict(audit)},ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',default=str(storage_path('18S_analysis')))
    parser.add_argument('--mismatches',type=int,default=2)
    parser.add_argument('--background-mismatches',type=int,default=1)
    parser.add_argument('--coverage-min',type=int)
    parser.add_argument('--coverage-max',type=int)
    parser.add_argument('--exclude-species',action='append',default=[])
    parser.add_argument('--max-exception-species',type=int,default=3)
    args=parser.parse_args()
    analyze(args.data,args.mismatches,args.mismatches,args.background_mismatches,coverage_min=args.coverage_min,coverage_max=args.coverage_max,excluded_species=args.exclude_species,max_exception_species=args.max_exception_species)

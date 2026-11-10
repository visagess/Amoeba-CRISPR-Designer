"""User-selected target/background FASTA inputs, isolated from retrieved data."""
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import amoebacrisp as app
from ssu18s_policy import in_length_range,length_limits,MIN_18S_BP,MAX_18S_BP


def infer_species(text):
    text=text.replace('_',' ')
    match=re.search(r'\b(Naegleria|Entamoeba)[ .]+([a-z][a-z-]+)\b',text,re.I)
    if match and match[2].lower() not in ('sp','cf','aff'):
        return match[1].capitalize()+' '+match[2].lower()
    match=re.search(r'\bN[. ]+([a-z][a-z-]+)\b',text)
    if match:return 'Naegleria '+match[1].lower()
    match=re.search(r'\b([A-Z][a-z]+)\s+([a-z][a-z-]+)\b',text)
    return match[1]+' '+match[2] if match and match[2] not in ('sp','cf','aff') else ''


def read_records(path):
    opener=gzip.open if str(path).lower().endswith('.gz') else open
    header=None;parts=[]
    with opener(path,'rt',encoding='utf-8-sig') as handle:
        for line in handle:
            line=line.strip()
            if not line:continue
            if line.startswith('>'):
                if header is not None:yield header,''.join(parts).upper()
                header=line[1:].strip();parts=[]
                if not header:raise ValueError('FASTA序列缺少名称')
            else:
                if header is None:raise ValueError('文件不是FASTA格式：'+str(path))
                parts.append(re.sub(r'\s+','',line))
    if header is not None:yield header,''.join(parts).upper()


def import_inputs(manifest,root,min_bp=MIN_18S_BP,max_bp=MAX_18S_BP):
    from ssu18s import classify
    min_bp,max_bp=length_limits(min_bp,max_bp)
    manifest=Path(manifest);root=Path(root)
    config=json.loads(manifest.read_text(encoding='utf-8'))
    files=config['files'];targets=[f for f in files if f.get('target')]
    if len(targets)!=1:raise ValueError('请选择恰好一个目标FASTA文件')
    target_name=targets[0]['species'].strip()
    if len(target_name.split())<2:raise ValueError('请填写目标种完整学名')
    target_genus=target_name.split()[0].casefold()
    rows=[];rejected=[];seen=set()
    for item in files:
        path=Path(item['path']).resolve();species=item['species'].strip()
        if not species:raise ValueError('每个文件需要物种名称：'+path.name)
        if not item.get('target') and species.casefold()==target_name.casefold():
            raise ValueError('目标种的其他序列请合并进目标文件，不能将目标种设为背景')
        count=0
        for number,(header,sequence) in enumerate(read_records(path),1):
            count+=1
            inferred=infer_species(header)
            if inferred and inferred.casefold()!=species.casefold():
                raise ValueError(f'{path.name}包含{inferred}，与文件物种标签{species}不一致，请按物种拆分')
            if set(sequence)-set('ACGTRYSWKMBDHVNU'):
                raise ValueError(f'{path.name}的{header}含空位、蛋白字符或其他非核酸字符，请导入未比对的18S FASTA')
            sequence=sequence.replace('U','T')
            source_id=header.split()[0]
            key=(species,sequence)
            if key in seen:continue
            seen.add(key)
            role=classify(species,target_name,[item['species'].split()[0] for item in files])
            digest=hashlib.sha256((str(path)+'\0'+header+'\0'+sequence).encode()).hexdigest()[:16]
            row={'id':'local_'+digest,'accession':'local_'+digest,'organism':species,'role':role,'sequence':sequence,'length':len(sequence),
                 'title':header,'strain':'','isolate':'','taxid':'','source_file':str(path),'source_record_id':source_id,
                 'source_url':str(path),'source_accessions':source_id,'source_record_count':1,
                 'DNA_supported':True,'mol_type':'user-provided nucleotide FASTA','annotated_partial':False,
                 'ambiguous_bases':sum(b not in 'ACGT' for b in sequence),'sha256':hashlib.sha256(sequence.encode()).hexdigest(),
                 'annotation_method':'user_declared_18S_input','genomic_introns_preserved':True,
                 'record_start_1based':1,'record_end_1based':len(sequence),'strand':'+','circular_origin_wrap':False,'annotation_intervals':1}
            if in_length_range(row,min_bp,max_bp):rows.append(row)
            else:rejected.append({'file':str(path),'record':source_id,'length':len(sequence),'reason':f'outside_{min_bp}_{max_bp}'})
        if count==0:raise ValueError('FASTA为空：'+str(path))
    if not any(r['role']=='target' and not r['ambiguous_bases'] for r in rows):
        raise ValueError(f'目标文件没有{min_bp}–{max_bp}bp且无模糊碱基的18S序列')
    if not any(r['role'].startswith(('background_','unresolved_')) for r in rows):raise ValueError(f'背景文件没有{min_bp}–{max_bp}bp的有效18S序列')
    root.mkdir(parents=True,exist_ok=True)
    app.write_json(root/'curated_records.json',rows)
    app.write_json(root/'local_import.json',{'target':target_name,'files':files,'excluded_by_length':rejected,
                   'min_input_bp':min_bp,'max_input_bp':max_bp,
                   'database':'local','background_genera':sorted({item['species'].split()[0] for item in files}),
                   'input_assumption':'User-supplied sequences are nuclear 18S nucleotide sequences, one species per file. No GenBank feature annotation is inferred.'})
    app.write_fasta(root/'all_curated_18S.fasta',((r['id'],r['sequence']) for r in rows))
    inputs=root/'imported_files';inputs.mkdir(exist_ok=True)
    for i,item in enumerate(files):shutil.copy2(item['path'],inputs/(str(i+1)+'_'+Path(item['path']).name))
    print(f'本地FASTA导入：目标{target_name}，长度合格序列{len(rows)}条，范围外{len(rejected)}条',flush=True)
    return root

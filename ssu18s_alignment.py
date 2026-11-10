"""MAFFT nucleotide MSA, with reversible reference and physical DNA coordinates."""
from runtime import storage_path, settings
from array import array
import json
from pathlib import Path
import subprocess
import amoebacrisp as app
from ssu18s_policy import filter_records,run_limits


def prepare(root,threads=4):
    root=Path(root);out=root/'comparison';work=out/'MAFFT';work.mkdir(exist_ok=True)
    records=filter_records(json.loads((root/'curated_records.json').read_text()),root)
    fasta=work/'input.fasta';alignment=work/'aligned.fasta';meta=work/'alignment.json'
    app.write_fasta(fasta,((r['id'],r['sequence']) for r in records))
    digest=app.checksum(fasta)
    if alignment.exists() and meta.exists():
        previous=json.loads(meta.read_text())
        if previous.get('input_sha256')==digest and previous.get('alignment_sha256')==app.checksum(alignment):return alignment
    config=app.default_config()
    executable=config['linux_env_bin']+'/mafft'
    params=['--auto','--adjustdirectionaccurately','--thread',str(threads),app.linux_path(fasta)]
    cmd=app.linux_command(executable,params,config)
    temporary=alignment.with_suffix('.fasta.part')
    min_bp,max_bp=run_limits(root)
    print(f'MAFFT：{min_bp}–{max_bp}bp核18S序列（含rRNA分子类型记录），多序列比对及方向校正',flush=True)
    # MAFFT reopens /dev/stderr. Open its descriptors inside WSL so that
    # inherited Windows handles cannot cause /dev/stderr permission errors.
    wrapper="import subprocess,sys; out=open(sys.argv[1],'w'); log=open(sys.argv[2],'w'); subprocess.run(sys.argv[3:],stdout=out,stderr=log,check=True)"
    launch=app.linux_command(config['linux_env_bin']+'/python',['-c',wrapper,app.linux_path(temporary),app.linux_path(work/'mafft.log'),executable]+params,config)
    subprocess.run(launch,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
    aligned=list(app.read_fasta(temporary))
    if len(aligned)!=len(records) or len({len(seq) for _,seq in aligned})!=1:raise ValueError('MAFFT output is incomplete')
    temporary.replace(alignment)
    app.write_json(meta,{'input_sha256':digest,'alignment_sha256':app.checksum(alignment),'records':len(records),'command':cmd,
                        'method':'MAFFT --auto --adjustdirectionaccurately','scope':'eligible nuclear 18S nucleotide sequences, including rRNA source records'})
    return alignment


def project(records,alignment,reference):
    mapping={}
    for label,sequence in app.read_fasta(alignment):
        reverse=label.startswith('_R_');key=label[3:] if reverse else label
        if key in mapping:raise ValueError('Duplicate MAFFT sequence ID')
        mapping[key]=(sequence.upper(),reverse)
    ref_aligned,ref_reverse=mapping[reference['id']]
    if ref_reverse:
        mapping={key:(app.rc(seq),not reverse) for key,(seq,reverse) in mapping.items()}
        ref_aligned,_=mapping[reference['id']]
    if ref_aligned.replace('-','')!=reference['sequence']:raise ValueError('Reference DNA changed by alignment')
    columns=[i for i,b in enumerate(ref_aligned) if b!='-']
    bases={};positions={}
    for row in records:
        aligned,reverse=mapping[row['id']]
        expected=app.rc(row['sequence']) if reverse else row['sequence']
        if aligned.replace('-','')!=expected:raise ValueError('DNA changed by alignment: '+row['id'])
        physical=[];p=0
        for base in aligned:
            physical.append((len(row['sequence'])-1-p if reverse else p) if base!='-' else -1)
            if base!='-':p+=1
        bases[row['id']]=[aligned[c] if aligned[c] in 'ACGT' else 'N' for c in columns]
        positions[row['id']]=array('i',(physical[c] for c in columns))
    return bases,positions

if __name__=='__main__':prepare(str(storage_path('18S_analysis')))

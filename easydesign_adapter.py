"""Local Cas12a official design + official activity prediction, both strands."""
from runtime import storage_path, settings
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import amoebacrisp as app
from ssu18s_four_groups import analyze
from ssu18s_policy import filter_records

SOURCE=Path(settings().get('easydesign_source', storage_path('software/EasyDesign')))
PYTHON=Path(settings().get('easydesign_python', storage_path('software/easydesign-env/python.exe')))


def model_context(reference,row,context=10):
    spacer=row['sequence'];pam=row['pam']
    dna=reference if row['strand']=='+' else app.rc(reference)
    start=row['start']-1 if row['strand']=='+' else len(reference)-row['end']
    model_start=start-4;model_end=start+len(spacer)
    if model_start-context<0 or model_end+context>len(dna):return None
    model=dna[model_start:model_end]
    if model!=pam+spacer:raise ValueError('PAM/spacer orientation inconsistent with reference DNA')
    return model,dna[model_start-context:model_end+context]


def run(root,mismatches=2,coverage_min=None,coverage_max=None,background_min_mismatches=2,excluded_species=None,max_exception_species=3):
    root=Path(root);out=root/'comparison';work=out/'EasyDesign';work.mkdir(exist_ok=True)
    if not PYTHON.exists() or not (SOURCE/'easyDesign/models/regress/cas12a/saved_model.pb').exists():
        raise RuntimeError('EasyDesign Cas12a环境或模型缺失，未生成虚假活性分数')
    candidate=analyze(root,mismatches,mismatches,background_min_mismatches-1,window_length=21,write_outputs=False,
                      coverage_min=coverage_min,coverage_max=coverage_max,excluded_species=excluded_species,max_exception_species=max_exception_species)
    app.write_json(work/'candidate_regions_21nt.json',candidate)
    records=filter_records(json.loads((root/'curated_records.json').read_text()),root)
    records=[r for r in records if r['organism'] not in candidate['excluded_species']]
    reference=next(r for r in records if r['accession']==candidate['reference'])
    rows=[];skipped=0
    for group,items in candidate['groups'].items():
        if not group.endswith('_有PAM'):continue
        for item in items:
            context=model_context(reference['sequence'],item)
            if context is None:skipped+=1;continue
            model,target=context
            if len(model)!=25 or len(target)!=45:raise ValueError('Cas12a trained model requires 25+20 nt')
            rows.append({**item,'group':group,'model_sequence':model,'target_with_context':target,
                         'spacer_rna':item['sequence'].replace('T','U'),
                         'full_crrna_rna':app.default_config()['direct_repeat_rna']+item['sequence'].replace('T','U')})
    env=dict(os.environ,PYTHONPATH=str(SOURCE),TF_CPP_MIN_LOG_LEVEL='2',CUDA_VISIBLE_DEVICES='-1',
             TF_NUM_INTRAOP_THREADS='2',TF_NUM_INTEROP_THREADS='2',PYTHONUTF8='1',
             TEMP=str(storage_path('temp')),TMP=str(storage_path('temp')),PIP_CACHE_DIR=str(storage_path('software/pip-cache')))
    env['PATH']=str(PYTHON.parent/'Library/bin')+os.pathsep+str(PYTHON.parent)+os.pathsep+env.get('PATH','')
    cli_guides={}
    for strand,sequence in [('plus',reference['sequence']),('minus',app.rc(reference['sequence']))]:
        fasta=work/(strand+'.fasta');app.write_fasta(fasta,[(reference['id']+'_'+strand,sequence)])
        prefix=work/('official_design_'+strand)
        cmd=[str(PYTHON),str(SOURCE/'bin/design.py'),'sliding-window','fasta',str(fasta),
             '--window-size','25','--window-step','1','--out-tsv',str(prefix),
             '--obj','maximize-activity','--ai-model','cas12a','--guide-length','25',
             '--maximization-algorithm','greedy','--seed','1','--predict-activity-model-path',
             str(SOURCE/'easyDesign/models/classify/cas12a'),str(SOURCE/'easyDesign/models/regress/cas12a')]
        print('EasyDesign官方Cas12a设计：'+strand,flush=True)
        with (work/('design_'+strand+'.log')).open('w',encoding='utf-8') as log:
            subprocess.run(cmd,cwd=SOURCE,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        with Path(str(prefix)+'.tsv').open(encoding='utf-8-sig') as f:
            for item in csv.DictReader(f,delimiter='\t'):
                for guide,score in zip(item['target-sequences'].split(),item['guide-expected-activities'].split()):
                    cli_guides[(strand,guide)]=float(score)
    designed=[]
    for row in rows:
        key=('plus' if row['strand']=='+' else 'minus',row['model_sequence'])
        if key not in cli_guides:continue
        row['official_cli_activity']=cli_guides[key];designed.append(row)
    request=work/'prediction_inputs.json';response=work/'predictions.json'
    app.write_json(request,{'rows':designed})
    with (work/'prediction.log').open('w',encoding='utf-8') as log:
        subprocess.run([str(PYTHON),str(Path(__file__).with_name('easydesign_worker.py')),
                        '--source',str(SOURCE),'--input',str(request),'--output',str(response)],
                       env=env,cwd=SOURCE,stdout=log,stderr=subprocess.STDOUT,check=True)
    result=json.loads(response.read_text())
    for row in result['rows']:
        if not math.isclose(row['activity'],row['official_cli_activity'],rel_tol=1e-5,abs_tol=1e-5):
            raise RuntimeError('Official CLI and Predictor activity scores disagree')
    result['rows'].sort(key=lambda r:(-r['activity'],r['start'],r['strand']))
    for i,row in enumerate(result['rows'],1):row['design_id']=f'EasyDesign_Cas12a_{i:03d}'
    result.update(min_input_bp=candidate['min_input_bp'],max_input_bp=candidate['max_input_bp'],spacer_length=21,not_scored_20nt=True,
                  max_exception_species=max_exception_species,
                  excluded_species=candidate['excluded_species'],
                  coverage_min=candidate['coverage_min'],coverage_max=candidate['coverage_max'],background_min_mismatches=background_min_mismatches,
                  alignment_method=candidate['alignment_method'],
                  input_fasta_sha256=app.checksum(out/'analysis_input_18S.fasta'),
                  context_insufficient=skipped,input_records=len(records),source_url='https://github.com/scRNA-Compt/EasyDesign',
                  model_hashes={str(p.relative_to(SOURCE)):app.checksum(p) for p in (SOURCE/'easyDesign/models').rglob('*') if p.is_file()},
                  source_archive_sha256=app.checksum(SOURCE.parent/'EasyDesign-main.zip') if (SOURCE.parent/'EasyDesign-main.zip').exists() else None,
                  scope='Reference-target activity prediction; conserved genus candidates are not necessarily target-species-specific. No claim of experimentally measured activity or uniform activity across species. DR comes from CaSilico default, not a CaSilico validation of 21nt spacers.')
    app.write_json(out/'easydesign_results.json',result)
    fields=['design_id','group','reference','start','end','strand','pam','sequence','spacer_rna','full_crrna_rna',
            'gc_percent','activity','classification_score','regression_raw','highly_active','compared_species','total_species','max_mismatches','official_cli_activity']
    app.write_tsv(work/'Cas12a_activity_designs.tsv',result['rows'],fields)
    app.write_fasta(work/'Cas12a_crRNA_21nt.fasta',((r['design_id'],r['spacer_rna']) for r in result['rows']))
    app.write_fasta(work/'Cas12a_full_crRNA_default_DR.fasta',((r['design_id'],r['full_crrna_rna']) for r in result['rows']))
    print('EasyDesign Cas12a：实际设计并预测 '+str(len(result['rows']))+' 条21nt crRNA',flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',default=str(storage_path('18S_analysis')))
    parser.add_argument('--mismatches',type=int,default=2)
    parser.add_argument('--coverage-min',type=int)
    parser.add_argument('--coverage-max',type=int)
    parser.add_argument('--background-min-mismatches',type=int,default=2)
    parser.add_argument('--exclude-species',action='append',default=[])
    parser.add_argument('--max-exception-species',type=int,default=3)
    args=parser.parse_args();run(args.data,args.mismatches,args.coverage_min,args.coverage_max,args.background_min_mismatches,args.exclude_species,args.max_exception_species)

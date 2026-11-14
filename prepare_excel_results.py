"""Prepare only the four requested valid result groups."""
from runtime import storage_path, settings
import argparse
import json
import hashlib
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--data',default=str(storage_path('18S_analysis')))
parser.add_argument('--nopam-json')
parser.add_argument('--output-dir')
args=parser.parse_args()
out=Path(args.data)/'comparison'
dest=Path(args.output_dir) if args.output_dir else out/'excel_export';dest.mkdir(parents=True,exist_ok=True)
result=json.loads((out/'four_groups.json').read_text(encoding='utf-8'))
headers=['序列ID','参考Accession','起点（bp）','终点（bp）','链','DNA（5′→3′，20nt）','PAM','GC含量','目标优先序列支持数','保守覆盖物种数或背景已评估数','有效物种标签总数','覆盖物种内最大错配数','本窗口未评估物种','来源','保守覆盖物种名单','位点已评估物种数','非例外背景错配下界']
def exception_text(row):
    return '；'.join(f"{name}，最少{hit['min_mismatches']}错配，{hit['records_matching']}条序列，例{hit['example_accession']}:{hit['example_start']}-{hit['example_end']}({hit['example_strand']}) {hit['example_sequence']}" for name,hit in sorted(row.get('exception_matches',{}).items()))

tabs=[]
for name,items in result['groups'].items():
    specific=name.startswith('特异')
    extra_headers=['可区分物种数（有位点覆盖）','不特异物种数','不特异物种名单','匹配证据','未定种命中标签','未定种命中序列数','条件性结论'] if specific else []
    extra_widths=[32,22,78,110,100,28,48] if specific else []
    tabs.append({'name':name,'headers':headers+extra_headers,
        'rows':[[r['id'],r['reference'],r['start'],r['end'],r['strand'],r['sequence'],r['pam'],r['gc_percent']/100,
                 r['target_support'],r['compared_species'],r['total_species'],r['max_mismatches'],r['unassessed_species'],r['source'],
                 r['matching_species'],r['assessed_species'],r['background_min_mismatches_lower_bound']]+([r['distinguishable_species_count'],r['exception_species_count'],r['exception_species'],exception_text(r),r['unresolved_matching_labels'],r['unresolved_matching_records'],r['specificity_status']] if specific else []) for r in items],
        'widths':[27,22,14,14,8,35,12,14,26,32,24,28,72,72,90,26,28]+extra_widths,'formats':{'H':'0.0%'},'wrap_columns':['M','O']+(['T','U','V','X'] if specific else [])})
activity=json.loads((out/'easydesign_results.json').read_text(encoding='utf-8'))
if activity['input_records']!=result['input_records']:raise ValueError('EasyDesign与比较分析输入记录不一致')
if activity.get('excluded_species',[])!=result.get('excluded_species',[]):raise ValueError('例外物种已变化，请重新运行预测')
if activity.get('max_exception_species',0)!=result.get('max_exception_species',0):raise ValueError('逐段例外物种上限已变化，请重新预测')
for key in ('coverage_min','coverage_max','alignment_method','min_input_bp','max_input_bp'):
    if activity.get(key)!=result.get(key):raise ValueError('EasyDesign筛选设置已变化，请重新运行预测')
if activity.get('background_min_mismatches')!=result['background_mismatches']+1:raise ValueError('背景错配阈值已变化，请重新预测')
if activity['input_fasta_sha256']!=hashlib.sha256((out/'analysis_input_18S.fasta').read_bytes()).hexdigest():
    raise ValueError('EasyDesign输入已变化，请重新运行预测')
tabs.append({'name':'EasyDesign_21nt设计',
    'headers':['设计ID','类别','参考Accession','起点（bp）','终点（bp）','链','PAM','spacer DNA（21nt）','spacer RNA（21nt）',
               '完整crRNA（默认DR＋21nt）','GC含量','EasyDesign活性分数','分类原始分数','回归原始分数','官方阈值下高活性',
               '已覆盖物种数','有效物种标签总数','相对参考最大错配数','PAM＋spacer（模型25nt）','官方CLI核对分数','可区分物种数（有位点覆盖）','不特异物种数','不特异物种名单','例外匹配证据','未定种命中标签','未定种命中序列数','条件性结论'],
    'rows':[[r['design_id'],r['group'],r['reference'],r['start'],r['end'],r['strand'],r['pam'],r['sequence'],r['spacer_rna'],r['full_crrna_rna'],
             r['gc_percent']/100,r['activity'],r['classification_score'],r['regression_raw'],'是' if r['highly_active'] else '否',
             r['compared_species'],r['total_species'],r['max_mismatches'],r['model_sequence'],r['official_cli_activity'],r.get('distinguishable_species_count'),r.get('exception_species_count'),r.get('exception_species',''),exception_text(r),r.get('unresolved_matching_labels',''),r.get('unresolved_matching_records'),r.get('specificity_status','')] for r in activity['rows']],
    'widths':[32,26,22,14,14,8,12,34,34,60,14,24,23,23,26,22,24,28,42,26,32,22,78,110,100,28,48],
    'formats':{'K':'0.0%','L':'0.0000','M':'0.0000','N':'0.0000','T':'0.0000'},'wrap_columns':['W','X','Y','AA']})
payload={'four_groups':result,'activity':activity,'tabs':tabs}
if args.nopam_json:
    from nopam_design import result_tabs
    nopam=json.loads(Path(args.nopam_json).read_text(encoding='utf-8'))
    if nopam['source_sha256']!=hashlib.sha256((out/'four_groups.json').read_bytes()).hexdigest():
        raise ValueError('无PAM二次设计与原比较结果不一致，请重新运行')
    tabs.extend(result_tabs(nopam))
    payload['nopam_design']={key:value for key,value in nopam.items() if key!='rows'}
(dest/'inputs.json').write_text(json.dumps(payload,ensure_ascii=False),encoding='utf-8')
print(json.dumps(result['counts'],ensure_ascii=False))

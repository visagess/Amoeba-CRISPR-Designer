"""Configurable online and local nuclear 18S comparison/design workflow."""
import argparse
import datetime as dt
import os
from pathlib import Path
import subprocess
import sys
from runtime import storage_path
import amoebacrisp as app
from ssu18s import retrieve
from ssu18s_policy import length_limits,run_settings,run_limits


def parse_args(argv=None):
    parser=argparse.ArgumentParser()
    parser.add_argument('--data',help='任务目录；新任务默认使用独立目录')
    parser.add_argument('--target',help='目标种完整学名，默认Naegleria gruberi')
    parser.add_argument('--background-genus',action='append',help='可重复指定背景属，目标属自动纳入')
    parser.add_argument('--database',choices=['ncbi','amoebadb'])
    parser.add_argument('--min-bp',type=int,help='提取的18S最短长度，默认800')
    parser.add_argument('--max-bp',type=int,help='提取的18S最长长度，默认4000')
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--mismatches',type=int,default=2,choices=range(5))
    parser.add_argument('--coverage-min',type=int)
    parser.add_argument('--coverage-max',type=int)
    parser.add_argument('--background-min-mismatches',type=int,default=2,choices=range(1,21))
    parser.add_argument('--local-manifest')
    parser.add_argument('--reuse-data',action='store_true')
    parser.add_argument('--max-exception-species',type=int,default=3)
    return parser.parse_args(argv)


def prepare_data(args):
    if args.threads<1:raise ValueError('线程数至少为1')
    if args.local_manifest and args.reuse_data:raise ValueError('本地导入和复用数据不能同时选择')
    target=args.target or 'Naegleria gruberi'
    stamp=dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    root=Path(args.data) if args.data else storage_path('18S_analysis','runs',stamp+'_'+app.slug(target))
    if args.reuse_data:
        if not args.data:raise ValueError('复用数据时请指定--data任务目录')
        if not (root/'curated_records.json').exists():raise ValueError('复用目录没有curated_records.json')
        saved=run_settings(root);minimum,maximum=run_limits(root)
        expected={'target':args.target,'database':args.database,'min_input_bp':args.min_bp,'max_input_bp':args.max_bp}
        saved={**saved,'min_input_bp':minimum,'max_input_bp':maximum}
        for key,value in expected.items():
            if value is not None and str(value).casefold()!=str(saved.get(key,'')).casefold():
                raise ValueError(f'复用数据的{key}与原任务不同，请创建新任务')
        if args.background_genus:
            requested={g.casefold() for g in args.background_genus}
            requested.add(saved.get('target',target).split()[0].casefold())
            if requested!={g.casefold() for g in saved.get('background_genera',[])}:
                raise ValueError('复用数据的背景范围与原任务不同，请创建新任务')
        return root
    if (root/'curated_records.json').exists() or (root/'analysis_request.json').exists():
        raise ValueError('任务目录已存在，请使用新的--data目录，或--reuse-data复用原任务')
    minimum,maximum=length_limits(args.min_bp if args.min_bp is not None else 800,args.max_bp if args.max_bp is not None else 4000)
    if args.local_manifest:
        from ssu18s_import import import_inputs
        import_inputs(args.local_manifest,root,minimum,maximum)
    else:
        retrieve(root,target,args.background_genus,minimum,maximum,args.database or 'ncbi')
    return root


def main(argv=None):
    args=parse_args(argv)
    root=prepare_data(args)
    excel_directory=storage_path('outputs',root.name)
    print(f'输出目录：{excel_directory}',flush=True)
    temporary=storage_path('temp');temporary.mkdir(parents=True,exist_ok=True)
    os.environ.update(TEMP=str(temporary),TMP=str(temporary))
    import ssu18s_compare
    old_argv=sys.argv
    try:
        sys.argv=[old_argv[0],'--data',str(root),'--threads',str(args.threads)]
        ssu18s_compare.main()
    finally:
        sys.argv=old_argv
    from ssu18s_four_groups import analyze
    result=analyze(root,args.mismatches,args.mismatches,args.background_min_mismatches-1,
                   coverage_min=args.coverage_min,coverage_max=args.coverage_max,max_exception_species=args.max_exception_species)
    from easydesign_adapter import run
    run(root,args.mismatches,args.coverage_min,args.coverage_max,args.background_min_mismatches,max_exception_species=args.max_exception_species)
    here=Path(__file__).resolve().parent
    output=excel_directory/f"18S结果_{app.slug(result['target'])}_{result['min_input_bp']}-{result['max_input_bp']}bp_例外最多{args.max_exception_species}种.xlsx"
    subprocess.run([sys.executable,str(here/'prepare_excel_results.py'),'--data',str(root)],check=True)
    subprocess.run([sys.executable,str(here/'export_18s_excel.py'),str(root/'comparison'/'excel_export'/'inputs.json'),str(output)],check=True)


if __name__=='__main__':
    main()

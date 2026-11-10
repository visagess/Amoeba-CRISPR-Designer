"""Verified 18S-specific workflow for the N. gruberi pilot."""
from runtime import storage_path, settings
import argparse
import json
import os
from pathlib import Path
import sys
import subprocess
import datetime as dt
from ssu18s import retrieve
import ssu18s_compare

parser=argparse.ArgumentParser()
parser.add_argument('--data',default=str(storage_path('18S_analysis')))
parser.add_argument('--threads',type=int,default=4)
parser.add_argument('--mismatches',type=int,default=2,choices=range(5))
parser.add_argument('--coverage-min',type=int)
parser.add_argument('--coverage-max',type=int)
parser.add_argument('--background-min-mismatches',type=int,default=2,choices=range(1,21))
parser.add_argument('--local-manifest')
parser.add_argument('--reuse-data',action='store_true')
parser.add_argument('--exclude-species',action='append',default=[])
parser.add_argument('--max-exception-species',type=int,default=3)
args=parser.parse_args()
root=Path(args.data)
if args.local_manifest and args.data==str(storage_path('18S_analysis')):
    root=root/'local_runs'/dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
excel_directory=Path(str(storage_path('outputs')))
if args.local_manifest:excel_directory=excel_directory/root.name
print(f'输出目录：{excel_directory}',flush=True)
temporary=Path(str(storage_path('temp')));temporary.mkdir(parents=True,exist_ok=True)
os.environ.update(TEMP=str(temporary),TMP=str(temporary))
if args.local_manifest:
    from ssu18s_import import import_inputs
    import_inputs(args.local_manifest,root)
elif not args.reuse_data:retrieve(root)
sys.argv=[sys.argv[0],'--data',str(root),'--threads',str(args.threads)]
ssu18s_compare.main()
from ssu18s_four_groups import analyze
analyze(root,args.mismatches,args.mismatches,args.background_min_mismatches-1,coverage_min=args.coverage_min,coverage_max=args.coverage_max,excluded_species=args.exclude_species,max_exception_species=args.max_exception_species)
from easydesign_adapter import run as easydesign_run
easydesign_run(root,args.mismatches,args.coverage_min,args.coverage_max,args.background_min_mismatches,args.exclude_species,args.max_exception_species)
here=Path(__file__).resolve().parent
result=json.loads((root/'comparison/four_groups.json').read_text())
output=excel_directory/f"18S结果_逐段例外最多{args.max_exception_species}种_覆盖{result['coverage_min']}-{result['coverage_max']}_MAFFT.xlsx"
subprocess.run([sys.executable,str(here/'prepare_excel_results.py'),'--data',str(root)],check=True)
subprocess.run([sys.executable,str(here/'export_18s_excel.py'),str(root/'comparison'/'excel_export'/'inputs.json'),str(output)],check=True)

"""Run with EasyDesign's isolated Python 3.8; official predictor only."""
import argparse
import json
from pathlib import Path
import sys

parser=argparse.ArgumentParser()
parser.add_argument('--source',required=True)
parser.add_argument('--input',required=True)
parser.add_argument('--output',required=True)
args=parser.parse_args()
sys.path.insert(0,args.source)
import numpy as np
import tensorflow as tf
from easyDesign.utils.predict_activity import Predictor

models=Path(args.source)/'easyDesign'/'models'
predictor=Predictor(str(models/'classify'/'cas12a'),str(models/'regress'/'cas12a'))
payload=json.loads(Path(args.input).read_text(encoding='utf-8'))
rows=payload['rows']
pairs=[(r['target_with_context'],r['model_sequence']) for r in rows]
if any(len(guide)!=predictor.guide_length or len(target)!=predictor.guide_length+2*predictor.context_nt for target,guide in pairs):
    raise ValueError('Model dimensions do not match actual target context; no padding/truncation permitted')
if pairs:
    encoded=predictor._model_input_from_nt(pairs)
    classification=predictor._predict_from_onehot(predictor.classification_model,encoded)
    regression=predictor._predict_from_onehot(predictor.regression_model,encoded)
    activities=predictor.compute_activity(0,pairs)
    highly_active=predictor.determine_highly_active(0,pairs)
    for row,cla,reg,activity,high in zip(rows,classification,regression,activities,highly_active):
        if not all(np.isfinite(x) for x in (cla,reg,activity)):raise ValueError('Non-finite model result')
        row.update(classification_score=float(cla),regression_raw=float(reg),activity=float(activity),highly_active=bool(high))
result={'rows':rows,'tensorflow':tf.__version__,'guide_length_with_PAM':predictor.guide_length,
        'context_nt':predictor.context_nt,'classification_threshold':predictor.classification_threshold,
        'regression_threshold_shifted':predictor.regression_threshold,
        'score_definition':'Official Predictor.compute_activity: inactive=0; otherwise max(0, raw regression + 4). This is a model score, not an experimental efficiency or a percentage.'}
Path(args.output).write_text(json.dumps(result,ensure_ascii=False),encoding='utf-8')
print('EasyDesign Cas12a model predictions:',len(rows),flush=True)

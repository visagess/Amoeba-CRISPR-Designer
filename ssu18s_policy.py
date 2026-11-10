"""Run-scoped inclusive length limits and compatibility with earlier datasets."""
import json
from pathlib import Path
MIN_18S_BP=800
MAX_18S_BP=4000


def length_limits(min_bp=MIN_18S_BP, max_bp=MAX_18S_BP):
    if type(min_bp) is not int or type(max_bp) is not int or not 20 <= min_bp <= max_bp:
        raise ValueError('18S长度需为整数，且20 ≤ 最短bp ≤ 最长bp')
    return min_bp, max_bp


def run_settings(root):
    for name in ('analysis_request.json', 'local_import.json', 'curation_summary.json'):
        path = Path(root) / name
        if path.exists():
            return json.loads(path.read_text(encoding='utf-8-sig'))
    return {}


def run_limits(root):
    config = run_settings(root)
    return length_limits(config.get('min_input_bp', MIN_18S_BP), config.get('max_input_bp', MAX_18S_BP))


def in_length_range(row, min_bp=MIN_18S_BP, max_bp=MAX_18S_BP):
    return min_bp <= len(row['sequence']) <= max_bp


def filter_records(records, root=None):
    minimum, maximum = run_limits(root) if root is not None else (MIN_18S_BP, MAX_18S_BP)
    legacy_roles = {'target_N_gruberi': 'target', 'background_Naegleria_named': 'background_genus',
                    'background_Entamoeba': 'background_other', 'unresolved_Naegleria': 'unresolved_genus'}
    return [{**r, 'role': legacy_roles.get(r['role'], r['role'])}
            for r in records if in_length_range(r, minimum, maximum)]


def gc_percent(sequence):
    if not sequence or set(sequence)-set('ACGT'):
        raise ValueError('GC calculation requires nonempty unambiguous DNA')
    return 100*(sequence.count('G')+sequence.count('C'))/len(sequence)

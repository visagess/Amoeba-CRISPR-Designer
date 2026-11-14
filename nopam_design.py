"""Design complementary spacer RNAs from existing no-PAM comparison windows."""
import argparse
import datetime as dt
import json
import math
from pathlib import Path
import subprocess
import sys
import amoebacrisp as app
from runtime import settings, storage_path


def candidates(result, five_prime_window=5):
    rows = []
    rejected_gc = 0
    for group, source_rows in result['groups'].items():
        if not group.endswith('_无PAM'):
            continue
        for source in source_rows:
            dna = source['sequence'].upper()
            if set(dna) - set('ACGT') or not 1 <= five_prime_window <= len(dna):
                raise ValueError('spacer长度或碱基无法用于5′端GC计算')
            if 100 * (dna.count('G') + dna.count('C')) >= 60 * len(dna):
                rejected_gc += 1
                continue
            original_strand = source['strand']
            orientations = [(original_strand, dna), ('-' if original_strand == '+' else '+', app.rc(dna))]
            for strand, target in orientations:
                spacer = app.rc(target).replace('T', 'U')
                prefix = spacer[:five_prime_window]
                rows.append({
                    'source_id': source['id'], 'group': group,
                    'reference': source['reference'], 'start': source['start'], 'end': source['end'],
                    'target_strand': strand, 'target_dna': target, 'spacer_rna': spacer,
                    'length': len(spacer), 'gc_percent': 100 * (spacer.count('G') + spacer.count('C')) / len(spacer),
                    'five_prime_gc_percent': 100 * (prefix.count('G') + prefix.count('C')) / len(prefix),
                    'target_support': source['target_support'], 'compared_species': source['compared_species'],
                    'matching_species': source.get('matching_species', ''),
                    'exception_species': source.get('exception_species', ''),
                    'unassessed_species': source.get('unassessed_species', ''),
                    'unresolved_matching_labels': source.get('unresolved_matching_labels', ''),
                    'specificity_status': source.get('specificity_status', ''),
                    'targeting_status': '无PAM候选，适用的Cas12a检测体系需确认',
                })
    return rows, rejected_gc


def fold_sequences(sequences, temperature, five_prime_window):
    import RNA
    RNA.params_load_RNA_Turner2004()
    model = RNA.md()
    model.temperature = temperature
    folds = {}
    for sequence in sorted(set(sequences)):
        if not sequence or set(sequence) - set('ACGU'):
            raise ValueError('结构预测需要明确的RNA碱基')
        compound = RNA.fold_compound(sequence, model)
        structure, mfe = compound.mfe()
        compound.exp_params_rescale(mfe)
        _, ensemble_energy = compound.pf()
        probabilities = compound.bpp()
        unpaired = []
        for position in range(1, five_prime_window + 1):
            paired_probability = sum(probabilities[min(position, other)][max(position, other)]
                                     for other in range(1, len(sequence) + 1) if other != position)
            unpaired.append(max(0.0, min(1.0, 1.0 - paired_probability)))
        values = (float(mfe), float(ensemble_energy), sum(unpaired) / len(unpaired))
        if len(structure) != len(sequence) or not all(math.isfinite(value) for value in values):
            raise ValueError('RNAfold返回无效结构或自由能')
        folds[sequence] = {'mfe_kcal_mol': values[0], 'ensemble_energy_kcal_mol': values[1],
                           'structure': structure, 'paired_nt': sum(base != '.' for base in structure),
                           'five_prime_mean_unpaired_probability': values[2]}
    return {'engine': 'ViennaRNA/RNAfold', 'version': RNA.__version__, 'material': 'RNA',
            'parameters': 'Turner2004', 'temperature_celsius': temperature, 'folds': folds}


def select_designs(rows, folding, minimum_mfe):
    selected = []
    for row in rows:
        if row['spacer_rna'] not in folding['folds']:
            raise ValueError('存在未完成结构计算的spacer')
        folded = folding['folds'][row['spacer_rna']]
        if not math.isfinite(folded['mfe_kcal_mol']):
            raise ValueError('结构自由能无效')
        if folded['mfe_kcal_mol'] >= minimum_mfe:
            selected.append({**row, **folded})
    selected.sort(key=lambda row: (row['five_prime_gc_percent'], row['paired_nt'],
                                   -row['mfe_kcal_mol'], -row['five_prime_mean_unpaired_probability'],
                                   row['gc_percent'], row['reference'], row['start'], row['target_strand'], row['group']))
    for index, row in enumerate(selected, 1):
        row['design_id'] = f'NP_{index:05d}'
    return selected


def result_tabs(result):
    headers = ['设计ID', '来源序列ID', '参考Accession', '起点（bp）', '终点（bp）', '靶DNA链',
               '靶DNA（5′→3′）', 'spacer RNA（5′→3′）', '长度（nt）', '总体GC', 'spacer 5′窗口GC',
               'MFE（kcal/mol）', '集合自由能（kcal/mol）', 'spacer二级结构', 'MFE配对碱基数',
               '5′窗口平均未配对概率', '目标支持数', '比较物种数', '覆盖物种', '例外物种',
               '未评估物种', '未定种命中', '特异性状态', '检测体系说明']
    fields = ['design_id', 'source_id', 'reference', 'start', 'end', 'target_strand', 'target_dna',
              'spacer_rna', 'length', 'gc_percent', 'five_prime_gc_percent', 'mfe_kcal_mol',
              'ensemble_energy_kcal_mol', 'structure', 'paired_nt', 'five_prime_mean_unpaired_probability',
              'target_support', 'compared_species', 'matching_species', 'exception_species',
              'unassessed_species', 'unresolved_matching_labels', 'specificity_status', 'targeting_status']
    tabs = []
    for label in ('保守序列', '特异序列'):
        table_rows = []
        for row in result['rows']:
            if row['group'] != label + '_无PAM':
                continue
            values = [row.get(field, '') for field in fields]
            values[9] /= 100
            values[10] /= 100
            table_rows.append(values)
        tabs.append({'name': 'RNAfold_无PAM' + label[:2], 'headers': headers, 'rows': table_rows,
                     'widths': [18, 30, 24, 13, 13, 12, 32, 32, 14, 14, 22, 21, 26, 32, 22, 28,
                                20, 22, 80, 80, 80, 80, 48, 52],
                     'formats': {'J': '0.0%', 'K': '0.0%', 'L': '0.00', 'M': '0.00', 'P': '0.0%'},
                     'wrap_columns': ['S', 'T', 'U', 'V', 'W', 'X']})
    return tabs


def run(root, five_prime_window=5, minimum_mfe=-2.0, temperature=37.0):
    if type(five_prime_window) is not int or five_prime_window < 1:
        raise ValueError('5′端窗口长度需为正整数')
    if not math.isfinite(minimum_mfe) or not math.isfinite(temperature) or not 0 <= temperature <= 100:
        raise ValueError('请输入有效自由能阈值和0–100°C温度')
    root = Path(root).resolve()
    source = root / 'comparison' / 'four_groups.json'
    source_digest = app.checksum(source)
    source_result = json.loads(source.read_text(encoding='utf-8'))
    python = Path(settings().get('structure_python', storage_path('software/structure-env/Scripts/python.exe')))
    if not python.exists():
        raise RuntimeError('请安装ViennaRNA并在config.json设置structure_python')
    stamp = dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    work = storage_path('outputs', root.name, 'noPAM_' + stamp)
    work.mkdir(parents=True, exist_ok=False)
    print(f'输出目录：{work}', flush=True)
    rows, rejected_gc = candidates(source_result, five_prime_window)
    payload = {'sequences': sorted({row['spacer_rna'] for row in rows}),
               'temperature': temperature, 'five_prime_window': five_prime_window}
    app.write_json(work / 'fold_input.json', payload)
    subprocess.run([str(python), str(Path(__file__).resolve()), '--worker', str(work / 'fold_input.json'),
                    str(work / 'folds.json')], check=True)
    folding = json.loads((work / 'folds.json').read_text(encoding='utf-8'))
    selected = select_designs(rows, folding, minimum_mfe)
    if app.checksum(source) != source_digest:
        raise RuntimeError('原比较结果已变化，请重新运行二次设计')
    result = {'target': source_result['target'], 'source_sha256': source_digest,
              'source_path': str(source), 'five_prime_window': five_prime_window,
              'minimum_mfe_kcal_mol': minimum_mfe, 'gc_upper_bound_exclusive': 60,
              'temperature_celsius': temperature, 'engine': folding['engine'], 'engine_version': folding['version'],
              'parameters': folding['parameters'], 'rejected_source_windows_by_gc': rejected_gc,
              'gc_passed_oriented_candidates': len(rows), 'rejected_by_mfe': len(rows) - len(selected),
              'selected_designs': len(selected), 'rows': selected,
              'ranking': 'spacer 5prime GC ascending; paired nt ascending; MFE descending; 5prime unpaired probability descending',
              'threshold_basis': 'Adjustable initial engineering screen requiring assay-specific experimental validation.',
              'scope': 'Isolated spacer RNA intramolecular structure; both complementary targeting orientations; source comparison scope and exception species retained.'}
    app.write_json(work / 'results.json', result)
    app.write_fasta(work / 'spacer_RNA.fasta', ((row['design_id'], row['spacer_rna']) for row in selected))
    here = Path(__file__).resolve().parent
    subprocess.run([sys.executable, str(here / 'prepare_excel_results.py'), '--data', str(root),
                    '--nopam-json', str(work / 'results.json'), '--output-dir', str(work / 'excel_export')], check=True)
    output = work / ('18S_无PAM_spacer设计_' + app.slug(result['target']) + '.xlsx')
    subprocess.run([sys.executable, str(here / 'export_18s_excel.py'),
                    str(work / 'excel_export/inputs.json'), str(output)], check=True)
    print(f'无PAM二次设计完成：总体GC<60%，5′前{five_prime_window}nt；结构筛选后{len(selected)}条', flush=True)
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data')
    parser.add_argument('--five-prime-window', type=int, default=5)
    parser.add_argument('--min-mfe', type=float, default=-2.0)
    parser.add_argument('--temperature', type=float, default=37.0)
    parser.add_argument('--worker', nargs=2, metavar=('INPUT', 'OUTPUT'))
    args = parser.parse_args()
    if args.worker:
        payload = json.loads(Path(args.worker[0]).read_text(encoding='utf-8'))
        output = fold_sequences(payload['sequences'], payload['temperature'], payload['five_prime_window'])
        app.write_json(args.worker[1], output)
    elif args.data:
        run(args.data, args.five_prime_window, args.min_mfe, args.temperature)
    else:
        parser.error('请指定--data分析任务目录')


if __name__ == '__main__':
    main()

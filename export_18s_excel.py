"""Portable Excel export using the public XlsxWriter dependency."""
import argparse
import json
import math
from pathlib import Path
import xlsxwriter


def export(source, destination):
    data = json.loads(Path(source).read_text(encoding='utf-8'))
    result = data['four_groups']
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with xlsxwriter.Workbook(str(destination), {'strings_to_formulas': False, 'strings_to_urls': False}) as book:
        header = book.add_format({'bold': True, 'bg_color': '#243B53', 'font_color': 'white', 'text_wrap': True, 'valign': 'vcenter'})
        wrap = book.add_format({'text_wrap': True, 'valign': 'top'})
        overview = book.add_worksheet('结果概览')
        overview.set_column(0, 0, 34)
        overview.set_column(1, 1, 28)
        overview.set_column(2, 2, 105, wrap)
        overview.write_row(0, 0, ['结果分组', '有效窗口数', '说明'], header)
        for index, tab in enumerate(data['tabs'], 1):
            overview.write_row(index, 0, [tab['name'], len(tab['rows']), '窗口可相互重叠，详细条件见下方及各表证据'])
            sheet = book.add_worksheet(tab['name'])
            sheet.freeze_panes(1, 1)
            sheet.hide_gridlines(2)
            for column, width in enumerate(tab['widths']):
                letter = xlsxwriter.utility.xl_col_to_name(column)
                properties = {'num_format': tab.get('formats', {}).get(letter, 'General'), 'valign': 'top'}
                if letter in tab.get('wrap_columns', []):
                    properties['text_wrap'] = True
                sheet.set_column(column, column, width, book.add_format(properties))
            sheet.set_row(0, 48)
            if tab['rows']:
                sheet.add_table(0, 0, len(tab['rows']), len(tab['headers']) - 1, {
                    'columns': [{'header': value, 'header_format': header} for value in tab['headers']],
                    'data': tab['rows'], 'style': 'Table Style Medium 2'})
                for row_number, values in enumerate(tab['rows'], 1):
                    lines = max((math.ceil(len(str(value or '')) / max(1, tab['widths'][column]))
                                 for column, value in enumerate(values)
                                 if xlsxwriter.utility.xl_col_to_name(column) in tab.get('wrap_columns', [])), default=1)
                    sheet.set_row(row_number, min(409, max(25, 16 * lines)))
            else:
                sheet.write_row(0, 0, tab['headers'], header)
        conditions = [
            ['目标种', result['target'], '参考 ' + result['reference']],
            ['输入长度', '800–4000 bp', '包含边界，范围外不参与分析'],
            ['保守覆盖物种区间', f"{result['coverage_min']}–{result['coverage_max']}", f"目标/保守窗口允许最多{result['conserved_mismatches']}个替换"],
            ['背景最少错配', result['background_mismatches'] + 1, '双链全部窗口筛查，任一记录近似命中使该物种成为例外'],
            ['每段例外物种上限', result['max_exception_species'], '各段分别统计，全部背景保留，同种多记录计一次'],
            ['未定种处理', '单列待核查', 'sp./cf./aff.命中不计入已定种数，命中候选为条件性结果；上限0时排除'],
            ['缺少位点覆盖', '未评估', '仅有位点覆盖且没有近似命中的已定种计入可区分种数'],
            ['比对方法', result['alignment_method'], 'MAFFT --auto --adjustdirectionaccurately，保留连续DNA坐标'],
            ['PAM', 'TTTV', '无PAM表是比较候选；EasyDesign仅预测有PAM的21nt spacer'],
            ['EasyDesign', '官方Cas12a模型', '4nt PAM＋21nt spacer，两侧各10nt实际上下文；独立重算21nt比较'],
            ['活性分数', '模型分数', '不是效率百分比。无活性分类为0，否则max(0, 原始回归值＋4)'],
            ['高活性阈值', data['activity']['regression_threshold_shifted'], '采用官方模型阈值'],
            ['CaSilico', '独立20nt设计', '完整21nt RNA为默认DR＋spacer组装，非CaSilico对21nt的验证'],
            ['数据边界', result['input_records'], '仅针对本轮输入记录，缺失数据不能作为特异证据'],
        ]
        start = len(data['tabs']) + 3
        overview.write_row(start, 0, ['分析条件', '设置', '说明'], header)
        for row, values in enumerate(conditions, start + 1):
            overview.write_row(row, 0, values)
            overview.set_row(row, 44)
        overview.freeze_panes(1, 0)
    return destination


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    args = parser.parse_args()
    print(export(args.input, args.output))

import json
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from export_18s_excel import export


class ExcelExportTests(unittest.TestCase):
    def test_empty_group_and_formula_like_text_are_preserved(self):
        result = {'target': 'Example species', 'reference': 'example', 'coverage_min': 1,
                  'coverage_max': 1, 'conserved_mismatches': 2, 'background_mismatches': 1,
                  'max_exception_species': 3, 'alignment_method': 'mafft', 'input_records': 2}
        tabs = [{'name': 'empty', 'headers': ['id'], 'rows': [], 'widths': [20]},
                {'name': 'text', 'headers': ['id'], 'rows': [['=1+1']], 'widths': [20]}]
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.json'
            source.write_text(json.dumps({'four_groups': result, 'activity': {'regression_threshold_shifted': 4}, 'tabs': tabs}), encoding='utf-8')
            destination = export(source, Path(directory) / 'result.xlsx')
            ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
            with zipfile.ZipFile(destination) as archive:
                empty = ET.fromstring(archive.read('xl/worksheets/sheet2.xml'))
                text = ET.fromstring(archive.read('xl/worksheets/sheet3.xml'))
                self.assertEqual(len(empty.findall('.//s:sheetData/s:row', ns)), 1)
                self.assertEqual(text.findall('.//s:f', ns), [])
                self.assertIn('=1+1', archive.read('xl/sharedStrings.xml').decode())

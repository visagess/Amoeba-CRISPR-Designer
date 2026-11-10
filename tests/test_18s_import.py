import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ssu18s_import as imp
import amoebacrisp as app


class ImportTests(unittest.TestCase):
    def test_filename_and_header_species(self):
        self.assertEqual(imp.infer_species('M18732.1_N.gruberi_18S'),'Naegleria gruberi')
        self.assertEqual(imp.infer_species('U80060.1_Naegleria_italica'),'Naegleria italica')
        self.assertEqual(imp.infer_species('unknown_accession'),'')

    def make_manifest(self,folder):
        target=folder/'target.fasta';background=folder/'background.fasta'
        app.write_fasta(target,[('target1','ACGT'*375),('short','A'*799)])
        app.write_fasta(background,[('background1','TGCA'*375),('long','A'*4001)])
        manifest=folder/'manifest.json'
        manifest.write_text(json.dumps({'files':[{'path':str(target),'species':'Naegleria gruberi','target':True},
                                                {'path':str(background),'species':'Naegleria italica','target':False}]}))
        return manifest

    def test_local_target_choice_length_filter_and_source_trace(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);manifest=self.make_manifest(folder);out=folder/'run'
            imp.import_inputs(manifest,out)
            rows=json.loads((out/'curated_records.json').read_text())
            self.assertEqual(len(rows),2)
            self.assertEqual([r['role'] for r in rows],['target','background_genus'])
            self.assertEqual(len(json.loads((out/'local_import.json').read_text())['excluded_by_length']),2)
            self.assertTrue(all(Path(r['source_file']).exists() for r in rows))
            self.assertEqual(len(list((out/'imported_files').glob('*'))),2)

    def test_target_file_mixed_species_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);manifest=self.make_manifest(folder)
            app.write_fasta(folder/'target.fasta',[('Naegleria fowleri','A'*1500)])
            with self.assertRaisesRegex(ValueError,'不一致'):imp.import_inputs(manifest,folder/'run')

    def test_target_must_be_selected_once(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);manifest=self.make_manifest(folder)
            data=json.loads(manifest.read_text());data['files'][1]['target']=True
            manifest.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'恰好一个'):imp.import_inputs(manifest,folder/'run')

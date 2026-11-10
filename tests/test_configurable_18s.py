import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import parse_qs,urlparse
import xml.etree.ElementTree as ET

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ssu18s as s
import ssu18s_import as importer
import ssu18s_run as workflow
from ssu18s_policy import filter_records,run_limits,length_limits


def genbank(organism,accession,length):
    record=ET.Element('GBSeq')
    for name,value in {'organism':organism,'accession-version':accession,'definition':'18S ribosomal RNA gene',
                       'sequence':'a'*length,'moltype':'DNA','taxonomy':'Eukaryota; Vahlkampfia'}.items():
        ET.SubElement(record,'GBSeq_'+name).text=value
    table=ET.SubElement(record,'GBSeq_feature-table');feature=ET.SubElement(table,'GBFeature')
    ET.SubElement(feature,'GBFeature_key').text='rRNA'
    ET.SubElement(feature,'GBFeature_location').text=f'1..{length}'
    intervals=ET.SubElement(feature,'GBFeature_intervals');interval=ET.SubElement(intervals,'GBInterval')
    ET.SubElement(interval,'GBInterval_from').text='1';ET.SubElement(interval,'GBInterval_to').text=str(length)
    quals=ET.SubElement(feature,'GBFeature_quals');qual=ET.SubElement(quals,'GBQualifier')
    ET.SubElement(qual,'GBQualifier_name').text='product';ET.SubElement(qual,'GBQualifier_value').text='18S ribosomal RNA'
    return record


class ConfigurableTests(unittest.TestCase):
    def test_arbitrary_target_and_unresolved_background(self):
        self.assertEqual(s.classify('Vahlkampfia avara','Vahlkampfia avara',['Hartmannella']),'target')
        self.assertEqual(s.classify('Vahlkampfia inornata','Vahlkampfia avara',['Hartmannella']),'background_genus')
        self.assertEqual(s.classify('Hartmannella cf. vermiformis','Vahlkampfia avara',['Hartmannella']),'unresolved_other')
        self.assertEqual(s.classify('Naegleria gruberi','Vahlkampfia avara',['Hartmannella']),'taxonomy_review')

    def test_pagination_keeps_all_records(self):
        def fetch(url):
            start=int(parse_qs(urlparse(url).query)['retstart'][0])
            return json.dumps({'esearchresult':{'count':'1003','idlist':[str(i) for i in range(start,min(start+1000,1003))]}}).encode()
        with patch.object(s.app,'fetch',side_effect=fetch):
            self.assertEqual(len(s.search_ids('nuccore','query')),1003)

    def test_incomplete_search_is_not_a_negative(self):
        with patch.object(s.app,'fetch',return_value=b'{"esearchresult":{"count":"1","idlist":[]}}'):
            with self.assertRaisesRegex(RuntimeError,'不完整'):s.search_ids('nuccore','query')

    def test_online_retrieval_custom_species_and_bp(self):
        taxa=[{'name':'Vahlkampfia','taxid':'1'},{'name':'Hartmannella','taxid':'2'}]
        def fetch(url,destination=None):
            query=parse_qs(urlparse(url).query)
            if 'term' in query:
                ids=['10','11'] if 'txid1[' in query['term'][0] else ['12']
                return json.dumps({'esearchresult':{'count':str(len(ids)),'idlist':ids}}).encode()
            tree=ET.Element('GBSet')
            for identifier in query['id'][0].split(','):
                organism,length={'10':('Vahlkampfia avara',600),'11':('Vahlkampfia inornata',650),'12':('Hartmannella vermiformis',1000)}[identifier]
                tree.append(genbank(organism,identifier+'.1',length))
            Path(destination).write_bytes(ET.tostring(tree))
            return destination
        with tempfile.TemporaryDirectory() as directory, patch.object(s,'search_scope',return_value=({'name':'Vahlkampfia avara','taxid':'10'},taxa)),patch.object(s.app,'fetch',side_effect=fetch):
            root=Path(directory)
            s.retrieve(root,'Vahlkampfia avara',['Hartmannella'],500,700)
            rows=filter_records(json.loads((root/'curated_records.json').read_text()),root)
            self.assertEqual(run_limits(root),(500,700))
            self.assertEqual({r['organism'] for r in rows},{'Vahlkampfia avara','Vahlkampfia inornata'})
            self.assertIn('500–700',(root/'18S数据报告.html').read_text(encoding='utf-8'))
            self.assertEqual(json.loads((root/'curation_summary.json').read_text())['target'],'Vahlkampfia avara')

    def test_local_import_bp_and_legacy_compatibility(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);target=root/'target.fa';background=root/'background.fa'
            target.write_text('>Vahlkampfia avara\n'+'A'*500+'\n')
            background.write_text('>Vahlkampfia inornata\n'+'C'*700+'\n')
            manifest=root/'inputs.json';manifest.write_text(json.dumps({'files':[
                {'path':str(target),'species':'Vahlkampfia avara','target':True},
                {'path':str(background),'species':'Vahlkampfia inornata','target':False}]}))
            importer.import_inputs(manifest,root/'run',500,700)
            self.assertEqual(run_limits(root/'run'),(500,700))
            rows=json.loads((root/'run/curated_records.json').read_text())
            self.assertEqual(len(filter_records(rows,root/'run')),2)
            legacy=[{'sequence':'A'*800,'role':'target_N_gruberi'}]
            self.assertEqual(filter_records(legacy)[0]['role'],'target')

    def test_reuse_rejects_changed_length(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'curated_records.json').write_text('[]')
            (root/'analysis_request.json').write_text(json.dumps({'target':'Vahlkampfia avara','min_input_bp':500,'max_input_bp':700}))
            args=workflow.parse_args(['--data',str(root),'--reuse-data','--min-bp','800'])
            with self.assertRaisesRegex(ValueError,'原任务不同'):workflow.prepare_data(args)

    def test_invalid_length_ranges(self):
        for values in [(0,100),(800,799),(20.5,100)]:
            with self.assertRaises(ValueError):length_limits(*values)

    def test_amoebadb_gff_locus_preserves_introns(self):
        from ssu18s_amoebadb import retrieve_features
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);genome=root/'genome.fa';gff=root/'annotation.gff'
            genome.write_text('>contig1\n'+'A'*100+'C'*100+'G'*100+'\n')
            gff.write_text('contig1\tDB\trRNA\t1\t100\t.\t+\t.\tID=ssu;product=18S ribosomal RNA\n'
                           'contig1\tDB\trRNA\t201\t300\t.\t+\t.\tID=ssu;product=18S ribosomal RNA\n')
            sample={'id':'AmoebaDB-68_VavaraStrain','genome':str(genome),'gff':str(gff),'release':68,'url':'https://example.org/genome.fa'}
            xml=b'<TaxaSet><Taxon><ScientificName>Vahlkampfia avara</ScientificName></Taxon></TaxaSet>'
            with patch.object(s,'search_ids',return_value=['1']),patch.object(s.app,'fetch',return_value=xml),patch.object(s.app,'amoebadb_download',return_value=[sample]):
                rows,rejected,queries,accessions=retrieve_features(root,{'name':'Vahlkampfia avara'},[{'name':'Vahlkampfia','taxid':'1'}],200,400)
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['sequence'],'A'*100+'C'*100+'G'*100)
            self.assertEqual(rows[0]['role'],'target')

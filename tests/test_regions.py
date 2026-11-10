import csv
import json
from pathlib import Path
import random
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import amoebacrisp as a
import pam_screen as p
import region_pipeline as r


class RegionTests(unittest.TestCase):
    def test_masks_are_half_open_and_keep_reference_endpoints(self):
        merged=r.merge_intervals([(8,12),(0,2),(2,4),(-3,1),(10,30)],20)
        self.assertEqual(merged,[(0,4),(8,20)])
        self.assertEqual(r.complement(merged,20),[(4,8)])
        self.assertEqual(r.complement([],20),[(0,20)])
        self.assertEqual(r.intersect_intervals([(0,10),(20,30)],[(5,25)]),[(5,10),(20,25)])
        self.assertEqual(r.merge_intervals([[0,4],(2,8)],10),[(0,8)])

    def test_target_species_never_enters_background(self):
        record={'accession':'GCA_x','paired_accession':'GCF_x','organism':{'tax_id':999,'organism_name':'Naegleria gruberi strain XX'}}
        self.assertTrue(r.same_species(record,'Naegleria gruberi',{5762},{'GCF_x'}))
        record['organism']['organism_name']='Naegleria lovaniensis'
        record.pop('paired_accession')
        self.assertFalse(r.same_species(record,'Naegleria gruberi',{5762},{'GCF_x'}))

    def test_region_coordinates_exclude_masks_and_do_not_invent_conservation(self):
        seq='ACGT'*200
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as d:
            root=Path(d)
            rows,dna=r.discover_regions({'c':seq},{'c':[(100,200)]},{'c':[(400,500)]},
                                       {'c':[(0,len(seq))]},root,SimpleNamespace(min_region=50,max_region=200),
                                       'single_reference_not_validated')
            self.assertTrue(rows)
            for row in rows:
                start,end=row['start_1based']-1,row['end_1based_inclusive']
                self.assertEqual(dna[row['region']],seq[start:end])
                self.assertFalse(start<200 and end>100)
                self.assertFalse(start<500 and end>400)
                self.assertEqual(row['conservation'],'single_reference_not_validated')
            rows,dna=r.discover_regions({'c':seq[:301]},{},{},{'c':[(0,301)]},root,
                                       SimpleNamespace(min_region=200,max_region=300),'single_reference_not_validated')
            self.assertEqual(len(rows),2)
            self.assertEqual(rows[0]['start_1based'],1)
            self.assertEqual(rows[-1]['end_1based_inclusive'],301)

    def test_failed_background_query_is_excluded_instead_of_called_unique(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as d:
            root=Path(d)
            a.write_json(root/'reference_query_map.json',{'REFQ0000001':{'contig':'c','offset_0based':40,'length':100}})
            (root/'background_blast.log').write_text('Warning: REFQ0000001: invalid query sequence; Karlin-Altschul parameters\n')
            masks=r.exclude_unassessable_queries({},root)
            self.assertEqual(masks,{'c':[(40,140)]})

    def test_packed_exhaustive_scan_equals_bruteforce_both_strands(self):
        random.seed(7398)
        guide='ACGTCAGTACGATCGTACGT'
        opposite='TGCATGACTGCATGCATGCA'
        changed=''.join(('A' if c!='A' else 'C') if i in (0,7,19) else c for i,c in enumerate(guide))
        # Deliberately dispersed substitutions, including both ends; no reliance on BLAST extension.
        seq=''.join(random.choices('ACGT',k=5000))+'TTTA'+changed
        seq+=a.rc('TTTC'+guide)+'TTTG'+guide+'N'*8+'TTTT'+guide
        seq+='TTTA'+guide[:7]+'N'+guide[8:]
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as d:
            root=Path(d)
            a.write_fasta(root/'b.fa',[('bg',seq)])
            a.write_fasta(root/'q.fa',[('q1',guide),('q2',opposite)])
            metadata=p.build_index(root/'b.fa',root/'index.bin')
            result=p.scan_queries(root/'q.fa',root/'index.bin',metadata,root/'hits.tsv',root/'summary.json',3,chunk_size=7)
            with (root/'hits.tsv').open() as handle:
                actual=list(csv.DictReader(handle,delimiter='\t'))
            expected=[]
            for name,g in [('q1',guide),('q2',opposite)]:
                for strand,dna in [('+',seq),('-',a.rc(seq))]:
                    for i in range(4,len(dna)-19):
                        pam=dna[i-4:i]
                        target=dna[i:i+20]
                        if pam[:3]!='TTT' or pam[3] not in 'ACG' or not set(target)<=set('ACGT'):
                            continue
                        distance=sum(x!=y for x,y in zip(g,target))
                        if distance<=3:
                            start=i+1 if strand=='+' else len(seq)-i-19
                            expected.append((name,start,strand,pam,distance))
            canonical=lambda hit:(hit['query'],int(hit['start']),hit['strand'],hit['pam'],int(hit['mismatches']))
            self.assertEqual(sorted(map(canonical,actual)),sorted(expected))
            self.assertTrue(any(hit[-1]==3 for hit in expected))
            self.assertEqual(result['queries']['q1']['near_hits'],sum(x[0]=='q1' for x in expected))


if __name__=='__main__':
    unittest.main()

import json
from pathlib import Path
import random
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ssu18s_four_groups as four


class FourGroupTests(unittest.TestCase):
    def test_pam_orientation_and_sequence_edges(self):
        self.assertEqual(four.pam_options('TTTA'+'C'*20,4),[('+','TTTA','C'*20)])
        self.assertEqual(four.pam_options('C'*20+'CAAA',0),[('-','TTTG','G'*20)])
        self.assertEqual(four.pam_options('C'*20,0),[])
        self.assertEqual(four.pam_options('TTTT'+'C'*20,4),[])

    def test_six_mismatch_boundary_and_nearest_background(self):
        query='A'*20
        background=[{'id':'five','sequence':'C'*5+'A'*15,'organism':'background','accession':'five'}]
        self.assertEqual(four.near_hits({query},background,5)[query],1)
        six=[{'id':'six','sequence':'C'*6+'A'*14,'organism':'background','accession':'six'}]
        self.assertEqual(four.near_hits({query},six,5)[query],0)
        self.assertEqual(four.nearest_background({query},six)[query]['mismatches'],6)

    def test_mafft_projection_preserves_insertions_and_reversed_reference(self):
        import ssu18s_alignment as msa
        reference={'id':'ref','sequence':'A'*20+'C'*20}
        other={'id':'other','sequence':'A'*20+'GG'+'C'*20}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'aligned.fasta'
            four.app.write_fasta(p,[('ref','A'*20+'--'+'C'*20),('other',other['sequence'])])
            bases,positions=msa.project([reference,other],p,reference)
            self.assertEqual([i for i,_ in four.contiguous_windows(bases['other'],positions['other'])],[0,20])
            four.app.write_fasta(p,[('_R_ref',four.app.rc('A'*20+'--'+'C'*20)),('_R_other',four.app.rc(other['sequence']))])
            reversed_bases,reversed_positions=msa.project([reference,other],p,reference)
            self.assertEqual(bases,reversed_bases)
            self.assertEqual(positions,reversed_positions)

    def test_windows_cannot_bridge_insertions_or_uncertain_bases(self):
        self.assertEqual(len(list(four.contiguous_windows('A'*40,list(range(40))))),21)
        positions=list(range(20))+list(range(24,44))
        self.assertEqual([i for i,_ in four.contiguous_windows('A'*40,positions)],[0,20])
        self.assertEqual(list(four.contiguous_windows('A'*19+'N'+'A'*19,list(range(39)))),[])

    def fixture(self,root,extra=None,tolerance=0,coverage_min=None,coverage_max=None,excluded_species=None):
        rng=random.Random(7)
        ref=list(''.join(rng.choice('ACGT') for _ in range(1500)))
        ref[:4]='TTTA';ref[46:50]='TTTC';ref=''.join(ref)
        other=list(ref)
        for i in (50,80):other[i]='A' if ref[i]!='A' else 'C'
        other=''.join(other)
        def record(sid,seq,role,organism):
            return dict(id=sid,sequence=seq,role=role,organism=organism,DNA_supported=True,accession=sid)
        rows=[record('ref',ref,'target_N_gruberi','Naegleria gruberi'),
              record('target2',ref,'target_N_gruberi','Naegleria gruberi'),
              record('other',other,'background_Naegleria_named','Naegleria other')]
        if extra=='unknown':rows.append(record('unknown',ref,'unresolved_Naegleria','Naegleria sp.'))
        if extra=='missing':rows.append(record('missing','A'*20,'background_Naegleria_named','Naegleria missing'))
        if extra=='target_conflict':rows.append(record('partial',other,'target_N_gruberi','Naegleria gruberi'))
        root=Path(root);out=root/'comparison';out.mkdir()
        (root/'curated_records.json').write_text(json.dumps(rows))
        (out/'comparison_cohorts.json').write_text(json.dumps({'reference':'ref','target_priority':['ref','target2']}))
        (out/'reference_alignments.tsv').write_text('\n'.join('\t'.join(map(str,['ref',r['id'],1,1500,1,1500,98,ref,r['sequence'],200]))
            for r in rows if len(r['sequence'])==1500))
        return four.analyze(root,tolerance,tolerance,tolerance,alignment_method='blast',coverage_min=coverage_min,coverage_max=coverage_max,excluded_species=excluded_species)

    def test_global_species_exclusions_are_rejected(self):
        with tempfile.TemporaryDirectory() as d,self.assertRaisesRegex(ValueError,'保留全部背景'):
            self.fixture(d,excluded_species=['Naegleria other'])

    def test_exception_count_is_per_species_and_records_are_aggregated(self):
        query='A'*20
        rows=[{'id':'a','accession':'a','organism':'Species alpha','sequence':'C'+'A'*19},
              {'id':'b','accession':'b','organism':'Species alpha','sequence':'A'*20},
              {'id':'c','accession':'c','organism':'Species beta','sequence':'CC'+'A'*18}]
        hits=four.species_matches({query},rows,1)[query]
        self.assertEqual(set(hits),{'Species alpha'})
        self.assertEqual(hits['Species alpha']['records_matching'],2)
        self.assertEqual(hits['Species alpha']['min_mismatches'],0)
        self.assertEqual(hits['Species alpha']['example_accession'],'b')
        self.assertTrue(four.unresolved_name('Naegleria sp. 1'))
        self.assertTrue(four.unresolved_name('Entamoeba cf. bovis'))
        self.assertFalse(four.unresolved_name('Naegleria pagei'))

    def test_per_window_exception_budget_retains_background_and_changes_results(self):
        rng=random.Random(61)
        reference=''.join(rng.choice('ACGT') for _ in range(1500))
        divergent=list(reference)
        for i in range(50,80):divergent[i]='C' if reference[i]!='C' else 'G'
        divergent=''.join(divergent)
        rows=[dict(id='ref',sequence=reference,role='target_N_gruberi',organism='Naegleria gruberi',DNA_supported=True,accession='ref'),
              dict(id='a1',sequence=reference,role='background_Naegleria_named',organism='Naegleria alpha',DNA_supported=True,accession='a1'),
              dict(id='a2',sequence=reference,role='background_Naegleria_named',organism='Naegleria alpha',DNA_supported=True,accession='a2'),
              dict(id='b',sequence=divergent,role='background_Naegleria_named',organism='Naegleria beta',DNA_supported=True,accession='b')]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'comparison';out.mkdir()
            (root/'curated_records.json').write_text(json.dumps(rows))
            (out/'comparison_cohorts.json').write_text(json.dumps({'reference':'ref','target_priority':['ref']}))
            (out/'reference_alignments.tsv').write_text('\n'.join('\t'.join(map(str,['ref',r['id'],1,1500,1,1500,98,reference,r['sequence'],2000])) for r in rows))
            strict=four.analyze(root,alignment_method='blast',max_exception_species=0,write_outputs=False)
            relaxed=four.analyze(root,alignment_method='blast',max_exception_species=1,write_outputs=False)
        self.assertEqual(strict['counts']['特异序列_有PAM']+strict['counts']['特异序列_无PAM'],0)
        kept=relaxed['groups']['特异序列_有PAM']+relaxed['groups']['特异序列_无PAM']
        self.assertTrue(kept)
        self.assertEqual(relaxed['input_records'],4)
        self.assertTrue(all(r['exception_species_count']==1 and r['exception_species']=='Naegleria alpha' for r in kept))
        self.assertTrue(all(r['exception_matches']['Naegleria alpha']['records_matching']==2 for r in kept))
        self.assertTrue(all(r['distinguishable_species_count']==1 and r['distinguishable_species']=='Naegleria beta' for r in kept))

    def test_species_coverage_interval_counts_matching_species_not_just_data(self):
        with tempfile.TemporaryDirectory() as d:full=self.fixture(d)
        with tempfile.TemporaryDirectory() as d:subset=self.fixture(d,coverage_min=1,coverage_max=1)
        subset_rows=subset['groups']['保守序列_有PAM']+subset['groups']['保守序列_无PAM']
        full_rows=full['groups']['保守序列_有PAM']+full['groups']['保守序列_无PAM']
        self.assertTrue(subset_rows)
        self.assertTrue(all(r['compared_species']==1 and r['assessed_species']==2 for r in subset_rows))
        self.assertTrue(all(r['matching_species']=='Naegleria gruberi' for r in subset_rows))
        self.assertFalse({r['start'] for r in subset_rows}&{r['start'] for r in full_rows})

    def test_substitution_tolerance_increases_conserved_windows(self):
        with tempfile.TemporaryDirectory() as d:strict=self.fixture(d)
        with tempfile.TemporaryDirectory() as d:relaxed=self.fixture(d,tolerance=2)
        self.assertGreater(relaxed['counts']['保守序列_无PAM'],strict['counts']['保守序列_无PAM'])
        self.assertEqual(relaxed['counts']['特异序列_有PAM']+relaxed['counts']['特异序列_无PAM'],0)
        self.assertTrue(all(r['max_mismatches']<=2 for group in relaxed['groups'].values() for r in group))

    def test_near_matches_agree_with_exhaustive_both_strand_comparison(self):
        rng=random.Random(29)
        seed=''.join(rng.choice('ACGT') for _ in range(20))
        one=('A' if seed[0]!='A' else 'C')+seed[1:]
        two=one[:11]+('A' if one[11]!='A' else 'C')+one[12:]
        records=[{'sequence':'GGG'+four.app.rc(two)+'CCC'},{'sequence':seed},{'sequence':'N'*22}]
        queries={seed,one,two,'T'*20}
        for radius in (0,1,2,3,5):
            actual=four.near_hits(queries,records,radius)
            for query in queries:
                expected=sum(any(set(word)<=set('ACGT') and four.mismatch_count(word,query)<=radius
                    for seq in (r['sequence'],four.app.rc(r['sequence']))
                    for word in (seq[i:i+20] for i in range(len(seq)-19))) for r in records)
                self.assertEqual(actual[query],expected,(query,radius))

    def test_missing_short_fragment_cannot_veto_near_full_backbone(self):
        rng=random.Random(39)
        seq=''.join(rng.choice('ACGT') for _ in range(1500))
        rows=[dict(id='ref',sequence=seq,role='target_N_gruberi',organism='Naegleria gruberi',DNA_supported=True,accession='ref',length=1500),
              dict(id='other',sequence=seq,role='background_Naegleria_named',organism='Naegleria other',DNA_supported=True,accession='other',length=1500),
              dict(id='short',sequence='A'*21,role='background_Naegleria_named',organism='Naegleria short',DNA_supported=True,accession='short',length=21)]
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'comparison';out.mkdir()
            (root/'curated_records.json').write_text(json.dumps(rows))
            (out/'comparison_cohorts.json').write_text(json.dumps({'reference':'ref','target_priority':['ref']}))
            (out/'reference_alignments.tsv').write_text('\n'.join('\t'.join(map(str,['ref',r['id'],1,1500,1,1500,100,seq,seq,2000])) for r in rows[:2]))
            result=four.analyze(root,alignment_method='blast')
        conserved=result['groups']['保守序列_有PAM']+result['groups']['保守序列_无PAM']
        self.assertTrue(conserved)
        self.assertTrue(all(r['compared_species']==2 and r['total_species']==2 and r['unassessed_species']=='' for r in conserved))
        self.assertEqual(result['excluded_by_length'],1)

    def test_four_groups_include_valid_windows_and_no_fake_negatives(self):
        with tempfile.TemporaryDirectory() as d:
            baseline=self.fixture(d)
        self.assertTrue(all(n>0 for n in baseline['counts'].values()),baseline['counts'])
        for extra in ('unknown','target_conflict'):
            with self.subTest(extra=extra),tempfile.TemporaryDirectory() as d:
                result=self.fixture(d,extra)
            self.assertEqual(result['counts']['特异序列_有PAM']+result['counts']['特异序列_无PAM'],0)
        with tempfile.TemporaryDirectory() as d:short=self.fixture(d,'missing')
        self.assertEqual(short['counts'],baseline['counts'])

    def test_hard_length_boundaries_and_gc(self):
        import ssu18s_policy as policy
        for length,expected in [(799,False),(800,True),(4000,True),(4001,False)]:
            self.assertEqual(policy.in_length_range({'sequence':'A'*length}),expected)
        self.assertEqual(policy.gc_percent('ACGT'),50)
        self.assertEqual(policy.gc_percent('GGCC'),100)
        self.assertEqual(policy.gc_percent('AAAA'),0)

if __name__=='__main__':unittest.main()

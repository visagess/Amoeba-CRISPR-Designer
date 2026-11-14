from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import amoebacrisp as app
from nopam_design import candidates, select_designs, result_tabs


def source(sequence):
    return {'id':'source1','sequence':sequence,'strand':'+','reference':'reference1',
            'start':21,'end':20+len(sequence),'target_support':2,'compared_species':5,
            'exception_species':'Example exception','specificity_status':'conditional'}


class NoPamTests(unittest.TestCase):
    def test_strict_gc_and_reverse_complement_orientation(self):
        result={'groups':{'特异序列_无PAM':[source('G'*12+'A'*8),source('G'*11+'A'*9)],
                          '保守序列_有PAM':[source('A'*20)]}}
        rows,rejected=candidates(result,5)
        self.assertEqual(rejected,1)
        self.assertEqual(len(rows),2)
        self.assertEqual({row['target_strand'] for row in rows},{'+','-'})
        self.assertEqual({row['five_prime_gc_percent'] for row in rows},{0,100})
        for row in rows:
            self.assertEqual(row['spacer_rna'],app.rc(row['target_dna']).replace('T','U'))
            self.assertLess(row['gc_percent'],60)
            self.assertEqual(row['exception_species'],'Example exception')

    def test_energy_boundary_and_low_five_prime_gc_ranking(self):
        rows,_=candidates({'groups':{'保守序列_无PAM':[source('GGGGG'+'A'*15)]}},5)
        folds={row['spacer_rna']:{'mfe_kcal_mol':-2.0,'paired_nt':4,
                                'five_prime_mean_unpaired_probability':0.8} for row in rows}
        selected=select_designs(rows,{'folds':folds},-2.0)
        self.assertEqual(len(selected),2)
        self.assertEqual(selected[0]['five_prime_gc_percent'],0)
        folds[rows[0]['spacer_rna']]['mfe_kcal_mol']=-2.01
        self.assertEqual(len(select_designs(rows,{'folds':folds},-2.0)),1)

    def test_missing_structure_cannot_pass(self):
        rows,_=candidates({'groups':{'保守序列_无PAM':[source('A'*20)]}},5)
        with self.assertRaisesRegex(ValueError,'未完成结构计算'):
            select_designs(rows,{'folds':{}},-2.0)

    def test_window_and_empty_output(self):
        with self.assertRaises(ValueError):candidates({'groups':{'保守序列_无PAM':[source('A'*20)]}},21)
        self.assertEqual([len(tab['rows']) for tab in result_tabs({'rows':[]})],[0,0])

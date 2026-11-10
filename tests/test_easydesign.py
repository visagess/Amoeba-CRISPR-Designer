import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import easydesign_adapter as ed


class EasyDesignInputTests(unittest.TestCase):
    def test_real_context_and_pam_are_oriented_on_both_strands(self):
        spacer='ACGTCGATCGATCGTAGCTAG'
        plus='C'*10+'TTTA'+spacer+'G'*10
        row={'start':15,'end':35,'strand':'+','sequence':spacer,'pam':'TTTA'}
        self.assertEqual(ed.model_context(plus,row),('TTTA'+spacer,plus))
        minus=ed.app.rc(plus)
        row.update(start=11,end=31,strand='-')
        self.assertEqual(ed.model_context(minus,row),('TTTA'+spacer,plus))
        self.assertEqual(len(ed.model_context(minus,row)[0]),25)
        self.assertEqual(len(ed.model_context(minus,row)[1]),45)

    def test_insufficient_context_is_not_padded_or_fabricated(self):
        spacer='ACGTCGATCGATCGTAGCTAG'
        row={'start':5,'end':25,'strand':'+','sequence':spacer,'pam':'TTTA'}
        self.assertIsNone(ed.model_context('TTTA'+spacer,row))

if __name__=='__main__':unittest.main()

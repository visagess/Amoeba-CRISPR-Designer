import sys
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ssu18s as s
import ssu18s_compare as compare

def feature(intervals):
    f=ET.Element('GBFeature');parts=ET.SubElement(f,'GBFeature_intervals')
    for start,end in intervals:
        i=ET.SubElement(parts,'GBInterval')
        ET.SubElement(i,'GBInterval_from').text=str(start)
        ET.SubElement(i,'GBInterval_to').text=str(end)
    return f

class SSUTests(unittest.TestCase):
    def test_aliases_and_other_rrna_are_distinct(self):
        for name in ['18S ribosomal RNA','18SrRNA','SSU rRNA','SSU ribosomal RNA','SSUrDNA','small-subunit ribosomal RNA','small subunit (SSU) ribosomal RNA']:
            self.assertTrue(s.identifies_18s(name),name)
        for name in ['28S ribosomal RNA','5.8S ribosomal RNA','internal transcribed spacer 1']:
            self.assertFalse(s.identifies_18s(name),name)
    def test_genomic_introns_are_not_spliced_into_artificial_targets(self):
        dna,start,end,_,wraps,_=s.contiguous_feature_dna('AAAACCCCGGGGTTTT',feature([(1,4),(9,12)]))
        self.assertEqual(dna,'AAAACCCCGGGG')
        self.assertEqual((start,end,wraps),(1,12,False))
    def test_circular_origin_is_contiguous_on_the_molecule(self):
        dna,start,end,_,wraps,_=s.contiguous_feature_dna('AAAACCCCGGGGTTTT',feature([(13,16),(1,4)]),True)
        self.assertEqual(dna,'TTTTAAAA')
        self.assertEqual((start,end,wraps),(13,4,True))
    def test_unknown_species_not_automatically_a_negative(self):
        self.assertEqual(s.classify('Naegleria sp. PL0398'),'unresolved_Naegleria')
        self.assertEqual(s.classify('Naegleria cf. gruberi'),'unresolved_Naegleria')
        self.assertEqual(s.classify('Naegleria spitzbergeniensis'),'background_Naegleria_named')
    def test_conserved_dna_regions_cannot_bridge_an_insertion(self):
        reference='A'*25+'C'*25
        cohort=[{'id':'a','sequence':reference},{'id':'b','sequence':'A'*25+'GGGG'+'C'*25}]
        observations={r['id']:bytearray([1]*25+[2]*25) for r in cohort}
        positions={'a':list(range(50)),'b':list(range(25))+list(range(29,54))}
        steps={r['id']:[1]*50 for r in cohort}
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as d:
            rows=compare.conserved_profile(reference,cohort,observations,positions,steps,Path(d),'test')
            self.assertEqual([r['length'] for r in rows],[25,25])

if __name__=='__main__':unittest.main()

import itertools
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import amoebacrisp as app


class PipelineTests(unittest.TestCase):
    def test_pam_orientation_and_coordinates(self):
        guide = "ACGTACGTACGTACGTACGA"
        seq = "GG" + "TTTA" + guide + "CC"
        self.assertEqual(app.exact_sites(seq, guide), [{"start": 7, "end": 26, "strand": "+", "pam": "TTTA"}])
        self.assertEqual(app.exact_sites(app.rc(seq), guide), [{"start": 3, "end": 22, "strand": "-", "pam": "TTTA"}])
        self.assertEqual(app.exact_sites("TTTT" + guide, guide), [])
        self.assertEqual(app.exact_sites(guide, guide), [])

    def test_seeded_mismatch_search_against_bruteforce(self):
        random.seed(902)
        guide = "ACGTCAGTACGATCGTACGT"
        mutated = guide[:3] + "A" + guide[4:10] + "T" + guide[11:]
        seq = "".join(random.choices("ACGT", k=1000)) + "TTTC" + mutated
        seq += app.rc("TTTG" + guide) + "TTTT" + guide
        for k in range(5):
            expected = []
            for strand, dna in [("+", seq), ("-", app.rc(seq))]:
                for i in range(4, len(dna)-19):
                    pam = dna[i-4:i]
                    mismatches = sum(a != b for a,b in zip(guide, dna[i:i+20]))
                    if pam[:3] == "TTT" and pam[3] in "ACG" and mismatches <= k:
                        start, end = (i+1, i+20) if strand == "+" else (len(seq)-i-19, len(seq)-i)
                        expected.append({"start": start, "end": end, "strand": strand, "pam": pam, "mismatches": mismatches})
            key = lambda h: (h["strand"], h["start"])
            self.assertEqual(sorted(app.mismatch_sites(seq, guide, k), key=key), sorted(expected, key=key))

    def test_spliced_proteins_and_contiguous_dna_isoform_selection(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            coding = "ATG" + "GCT"*35 + "TAA"
            genome = coding[:31] + "CCCCCCCCCC" + coding[31:]
            app.write_fasta(root / "genome.fa", [("chr1", genome)])
            gff = ["##gff-version 3", f"chr1\tx\tgene\t1\t{len(genome)}\t.\t+\t.\tID=g",
                   f"chr1\tx\tmRNA\t1\t{len(genome)}\t.\t+\t.\tID=t;Parent=g",
                   "chr1\tx\tCDS\t1\t31\t.\t+\t0\tID=c;Parent=t",
                   f"chr1\tx\tCDS\t42\t{len(genome)}\t.\t+\t2\tID=c;Parent=t"]
            (root / "ann.gff").write_text("\n".join(gff))
            models = app.extract_models({"id": "a", "genome": str(root / "genome.fa"), "gff": str(root / "ann.gff")}, root / "proteins")
            self.assertEqual(len(models), 1)
            model = next(iter(models.values()))
            self.assertEqual(model["protein"], "M" + "A"*35)
            self.assertEqual(model["dna"], genome)
            self.assertNotEqual(model["dna"], coding)
            self.assertEqual(model["gene"], "g")

    def test_html_widget_full_data_not_visible_page(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "result.html"
            widget = {"x": {"container": "<table><tr><th></th><th>Guide_RNA</th><th>Guide_RNA_ID</th></tr></table>",
                             "data": [[str(i) for i in range(100)], ["ACGT"*5]*100, [[f"id{i}"] for i in range(100)]]}}
            path.write_text('<script type="application/json">'+json.dumps(widget)+'</script>')
            rows = app.parse_casilico(path)
            self.assertEqual(len(rows), 100)
            self.assertEqual(rows[-1]["Guide_RNA_ID"], "id99")

    def test_minus_strand_cds_and_longest_isoform(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            coding = "ATG" + "GCT"*40 + "TAA"
            oriented = coding[:34] + "CCCCCCC" + coding[34:]
            genome = app.rc(oriented)
            n = len(genome)
            app.write_fasta(root / "genome.fa", [("c", genome)])
            features = [("gene", 1, n, ".", "ID=g"), ("mRNA", 1, n, ".", "ID=t;Parent=g"),
                        ("CDS", n-33, n, "0", "ID=x;Parent=t"),
                        ("CDS", 1, n-41, "2", "ID=x;Parent=t"),
                        ("mRNA", n-92, n, ".", "ID=short;Parent=g"),
                        ("CDS", n-92, n, "0", "ID=y;Parent=short")]
            (root / "g.gff").write_text("\n".join(f"c\tx\t{k}\t{s}\t{e}\t.\t-\t{p}\t{a}" for k,s,e,p,a in features))
            models = app.extract_models({"id": "negative", "genome": str(root / "genome.fa"), "gff": str(root / "g.gff")}, root / "p")
            self.assertEqual(len(models), 1)
            model = next(iter(models.values()))
            self.assertEqual(model["protein"], "M"+"A"*40)
            self.assertEqual(model["dna"], oriented)
            self.assertEqual(model["transcript"], "t")

    def test_paired_assemblies_not_double_counted(self):
        old = app.fetch
        reports = [{"accession": "GCA_1.1", "paired_accession": "GCF_1.1"},
                   {"accession": "GCF_1.1", "paired_accession": "GCA_1.1"}]
        app.fetch = lambda *a, **kw: json.dumps({"reports": reports}).encode()
        try:
            self.assertEqual([r["accession"] for r in app.ncbi_records("test", 0)], ["GCF_1.1"])
        finally:
            app.fetch = old


if __name__ == "__main__":
    unittest.main()

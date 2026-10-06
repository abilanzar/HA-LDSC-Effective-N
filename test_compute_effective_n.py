import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from compute_effective_n import convert, effective_n


class EffectiveNTests(unittest.TestCase):
    def test_formula(self):
        self.assertAlmostEqual(effective_n(5, 15), 15.0)

    def test_preserves_snp_alleles_z_and_sums_contributors(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "counts.tsv").write_text("cohort\tcases\tcontrols\nA\t10\t90\nB\t20\t80\n")
            with gzip.open(root / "source.gz", "wt") as f:
                f.write("SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t200\t1.25\nrs2\tC\tT\t100\t-0.5\n")
            with gzip.open(root / "membership.gz", "wt") as f:
                f.write("SNP\tA\tB\nrs1\t1\t1\nrs2\t0\t1\n")
            out = root / "neff.gz"
            self.assertEqual(convert(root / "source.gz", root / "membership.gz", root / "counts.tsv", out), 2)
            with gzip.open(out, "rt") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual([(r["SNP"], r["A1"], r["A2"], r["Z"]) for r in rows],
                             [("rs1", "A", "G", "1.25"), ("rs2", "C", "T", "-0.5")])
            self.assertAlmostEqual(float(rows[0]["N"]), 100.0)
            self.assertAlmostEqual(float(rows[1]["N"]), 64.0)
            with self.assertRaises(FileExistsError):
                convert(root / "source.gz", root / "membership.gz", root / "counts.tsv", out)

    def test_mismatched_snp_fails_without_output(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "counts.tsv").write_text("cohort\tcases\tcontrols\nA\t10\t90\n")
            with gzip.open(root / "source.gz", "wt") as f:
                f.write("SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t100\t1\n")
            with gzip.open(root / "membership.gz", "wt") as f:
                f.write("SNP\tA\nrs2\t1\n")
            out = root / "neff.gz"
            with self.assertRaises(ValueError):
                convert(root / "source.gz", root / "membership.gz", root / "counts.tsv", out)
            self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()

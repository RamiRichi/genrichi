"""Tests for workflow/scripts/build_cds_bed.py and the panel/CDS pairing in the configs.

Synthetic GTF cases run anywhere; the pairing checks read the committed panel
BEDs and the local MANE GTF (skipped if the GTF is absent).
"""
import importlib.util
import re
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_cds_bed", REPO / "workflow/scripts/build_cds_bed.py")
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)

PHASE1_BED = REPO / "resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1.bed"
PHASE1_CDS = REPO / "resources/panel/phase1_solid_tumor/solid_tumor_phase1_v1_cds.bed"
C62_BED = REPO / "resources/panel/comprehensive_genes.bed"
C62_CDS = REPO / "resources/panel/comprehensive_genes_cds.bed"
MANE = REPO / "resources/reference/mane/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz"
# The 62-gene pair is separate from the deployed Phase 1 default; its checks run
# only in trees that carry config/comprehensive_62gene_config.yaml.
HAVE_62_CONFIG = (REPO / "config/comprehensive_62gene_config.yaml").exists()


def gtf_line(feature, chrom, s, e, gene, tx, mane=True):
    tag = ' tag "MANE_Select";' if mane else ""
    return (f'{chrom}\tt\t{feature}\t{s}\t{e}\t.\t+\t0\t'
            f'gene_id "G_{gene}"; transcript_id "{tx}"; gene_name "{gene}";{tag}\n')


class SyntheticBuild(unittest.TestCase):
    def build(self, gtf_lines, bed_lines):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "m.gtf").write_text("".join(gtf_lines))
            (d / "p.bed").write_text("".join(bed_lines))
            return b.run(str(d / "m.gtf"), str(d / "p.bed"), str(d / "o.bed"))

    def test_coordinates_clip_merge_and_sort(self):
        gtf = [
            gtf_line("transcript", "chrX", 1, 1000, "GB", "TB"),
            gtf_line("CDS", "chrX", 101, 200, "GB", "TB"),      # 1-based -> [100,200)
            gtf_line("transcript", "chr2", 1, 1000, "GA", "TA"),
            gtf_line("CDS", "chr2", 101, 200, "GA", "TA"),
            gtf_line("CDS", "chr2", 201, 300, "GA", "TA"),      # adjacent -> merged with previous
            gtf_line("CDS", "chr2", 501, 600, "GA", "TA"),
            gtf_line("transcript", "chr2", 1, 1000, "GA", "TZ", mane=False),
            gtf_line("CDS", "chr2", 701, 800, "GA", "TZ", mane=False),   # non-MANE transcript ignored
        ]
        bed = ["chrX\t0\t1000\tGB\n", "chr2\t0\t550\tGA\n"]
        recs, bp = self.build(gtf, bed)
        self.assertEqual(recs, [("chr2", 100, 300, "GA"), ("chr2", 500, 550, "GA"), ("chrX", 100, 200, "GB")])
        self.assertEqual(bp, 200 + 50 + 100)

    def test_genes_never_merged_together(self):
        gtf = [gtf_line("transcript", "chr1", 1, 900, "G1", "T1"), gtf_line("CDS", "chr1", 101, 200, "G1", "T1"),
               gtf_line("transcript", "chr1", 1, 900, "G2", "T2"), gtf_line("CDS", "chr1", 201, 300, "G2", "T2")]
        recs, _ = self.build(gtf, ["chr1\t0\t900\tG1\n", "chr1\t0\t900\tG2\n"])
        self.assertEqual(recs, [("chr1", 100, 200, "G1"), ("chr1", 200, 300, "G2")])

    def test_clipping_uses_only_same_gene_regions(self):
        gtf = [gtf_line("transcript", "chr1", 1, 900, "G1", "T1"), gtf_line("CDS", "chr1", 101, 400, "G1", "T1"),
               gtf_line("transcript", "chr1", 1, 900, "G2", "T2"), gtf_line("CDS", "chr1", 601, 700, "G2", "T2")]
        recs, _ = self.build(gtf, ["chr1\t0\t250\tG1\n", "chr1\t500\t800\tG2\n"])
        self.assertEqual(recs, [("chr1", 100, 250, "G1"), ("chr1", 600, 700, "G2")])

    def test_missing_mane_transcript_fails(self):
        gtf = [gtf_line("transcript", "chr1", 1, 900, "G1", "T1", mane=False),
               gtf_line("CDS", "chr1", 101, 200, "G1", "T1", mane=False)]
        with self.assertRaisesRegex(b.CdsBuildError, "no MANE_Select transcript.*G1"):
            self.build(gtf, ["chr1\t0\t900\tG1\n"])

    def test_duplicate_mane_transcript_fails(self):
        gtf = [gtf_line("transcript", "chr1", 1, 900, "G1", "T1"), gtf_line("transcript", "chr1", 1, 900, "G1", "T2")]
        with self.assertRaisesRegex(b.CdsBuildError, "more than one MANE_Select.*G1"):
            self.build(gtf, ["chr1\t0\t900\tG1\n"])

    def test_gene_without_cds_in_panel_fails(self):
        gtf = [gtf_line("transcript", "chr1", 1, 900, "G1", "T1"), gtf_line("CDS", "chr1", 101, 200, "G1", "T1")]
        with self.assertRaisesRegex(b.CdsBuildError, "no CDS overlaps.*G1"):
            self.build(gtf, ["chr1\t500\t900\tG1\n"])

    def test_bed3_and_empty_interval_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "p.bed"
            p.write_text("chr1\t0\t100\n")
            with self.assertRaisesRegex(b.CdsBuildError, "gene name"):
                b.read_panel_bed(str(p))
            p.write_text("chr1\t100\t100\tG\n")
            with self.assertRaisesRegex(b.CdsBuildError, "empty or invalid"):
                b.read_panel_bed(str(p))

    def test_chromosome_order(self):
        order = sorted(["chrM", "chrY", "chrX", "chr10", "chr2", "chr22", "chr1"], key=b.chrom_key)
        self.assertEqual(order, ["chr1", "chr2", "chr10", "chr22", "chrX", "chrY", "chrM"])


def read_bed(path):
    return [tuple(line.rstrip("\r\n").split("\t")) for line in open(path) if line.strip() and not line.startswith("#")]


class CommittedPanelPairs(unittest.TestCase):
    def check_pair(self, bed, cds, n_genes):
        panel = {}
        for c, s, e, g in read_bed(bed):
            panel.setdefault(g, []).append((c, int(s), int(e)))
        rows = read_bed(cds)
        self.assertEqual({r[3] for r in rows}, set(panel))
        self.assertEqual(len(panel), n_genes)
        for c, s, e, g in rows:
            s, e = int(s), int(e)
            self.assertLess(s, e)
            self.assertTrue(any(c == pc and s >= ps and e <= pe for pc, ps, pe in panel[g]), (c, s, e, g))

    def test_phase1_pair(self):
        self.check_pair(PHASE1_BED, PHASE1_CDS, 55)
        rows = read_bed(PHASE1_CDS)
        self.assertEqual(len(rows), 901)
        self.assertEqual(sum(int(r[2]) - int(r[1]) for r in rows), 164745)

    @unittest.skipUnless(HAVE_62_CONFIG, "62-gene config not present in this tree")
    def test_62_gene_pair(self):
        self.check_pair(C62_BED, C62_CDS, 62)
        rows = read_bed(C62_CDS)
        self.assertEqual(len(rows), 1107)
        self.assertEqual(sum(int(r[2]) - int(r[1]) for r in rows), 209109)

    @unittest.skipUnless(MANE.exists(), "local MANE GTF not present")
    def test_committed_cds_reproducible_from_mane(self):
        pairs = [(PHASE1_BED, PHASE1_CDS)] + ([(C62_BED, C62_CDS)] if HAVE_62_CONFIG else [])
        for bed, cds in pairs:
            with tempfile.TemporaryDirectory() as d:
                out = Path(d) / "o.bed"
                b.run(str(MANE), str(bed), str(out))
                self.assertEqual(out.read_text(), Path(cds).read_text().replace("\r\n", "\n"))

    def test_configs_reference_matching_pairs(self):
        def cfg(name):
            t = (REPO / "config" / name).read_text(encoding="utf-8")
            bed = re.search(r"^  bed: (\S+)", t, re.M).group(1)
            cds = re.search(r"^  coding_bed: (\S+)", t, re.M).group(1)
            return bed, cds

        def rel(p):
            return str(p.relative_to(REPO)).replace("\\", "/")

        self.assertEqual(cfg("comprehensive_config.yaml"), (rel(PHASE1_BED), rel(PHASE1_CDS)))
        if HAVE_62_CONFIG:
            self.assertEqual(cfg("comprehensive_62gene_config.yaml"), (rel(C62_BED), rel(C62_CDS)))


if __name__ == "__main__":
    unittest.main()

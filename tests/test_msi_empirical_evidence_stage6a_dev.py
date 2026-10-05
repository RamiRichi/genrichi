"""
SYNTHETIC, DEVELOPMENT-ONLY tests for dev/msi_redesign/empirical_evidence_stage6a.py
and its independent verifier verify_empirical_evidence_stage6a.py.

*** Fabricated reference, BAMs, loci and manifests only. Real results come ***
*** from a separate run over the hash-verified Stage-5 chain.              ***

Tests that need `samtools` (synthetic BAM / FASTA) are skipped when it is
not installed (e.g. on the Windows host); they run in the WSL suite.

Run with:
    python -m unittest tests.test_msi_empirical_evidence_stage6a_dev -v
"""

import copy
import csv
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import empirical_evidence_stage6a as e6  # noqa: E402
import verify_empirical_evidence_stage6a as v6  # noqa: E402
import bam_extract  # noqa: E402  (prototype; used ONLY to cross-check the CIGAR logic)
import msi_locus_model  # noqa: E402

SAMTOOLS = shutil.which("samtools")
needs_samtools = unittest.skipUnless(SAMTOOLS, "samtools not installed")


# --------------------------------------------------------------------------- pure functions
def read(mapq=60, fb=20, fa=20, obs=12, rev=False, name="r", off=20, rl=40):
    return {"qname": name, "excluded_flag": False, "mapq": mapq, "reverse": rev, "flank_before": fb,
            "flank_after": fa, "observed_len": obs, "offset": off, "read_len": rl}


class TestPureAnalysis(unittest.TestCase):
    def test_parameters_are_the_provisional_prototype_values_and_unchanged(self):
        self.assertEqual((e6.MIN_MAPQ, e6.MIN_FLANK, e6.MIN_LOCUS_DEPTH), (20, 5, 20))
        self.assertEqual((bam_extract.MIN_MAPQ, bam_extract.MIN_FLANK, msi_locus_model.MIN_LOCUS_DEPTH), (20, 5, 20))
        self.assertEqual((v6.source_const("bam_extract.py", "MIN_MAPQ"), v6.source_const("bam_extract.py", "MIN_FLANK"),
                          v6.source_const("msi_locus_model.py", "MIN_LOCUS_DEPTH")), (20, 5, 20))

    def test_span_read(self):
        r = e6.analyse_alignment(0, 81, 60, "40M", 100, 112)      # covers ref 80..120
        self.assertEqual((r["observed_len"], r["flank_before"], r["flank_after"], r["offset"], r["read_len"]), (12, 20, 8, 20, 40))

    def test_deletion_and_insertion(self):
        self.assertEqual(e6.analyse_alignment(0, 81, 60, "20M2D18M", 100, 112)["observed_len"], 10)
        self.assertEqual(e6.analyse_alignment(0, 81, 60, "20M2I20M", 100, 112)["observed_len"], 14)

    def test_unmapped_and_non_overlapping_are_skipped(self):
        self.assertIsNone(e6.analyse_alignment(4, 81, 60, "*", 100, 112))
        self.assertIsNone(e6.analyse_alignment(0, 1, 60, "40M", 100, 112))
        self.assertTrue(e6.analyse_alignment(1024, 81, 60, "40M", 100, 112)["excluded_flag"])

    def test_agrees_with_the_prototype_on_random_cigars(self):
        rng = random.Random(6)
        checked = 0
        for _ in range(4000):
            ops = []
            for _ in range(rng.randint(1, 6)):
                ops.append((rng.randint(1, 30), rng.choice("MMMMIDS")))
            cigar = "".join(f"{n}{o}" for n, o in ops)
            pos = rng.randint(1, 120)
            start = rng.randint(60, 110)
            end = start + rng.randint(8, 30)
            mine = e6.analyse_alignment(0, pos, 60, cigar, start, end)
            proto = bam_extract.observed_repeat_length(bam_extract.SamRecord(0, pos, 60, cigar), start, end, min_flank=5)
            ok = mine is not None and mine["flank_before"] >= 5 and mine["flank_after"] >= 5
            if ok:
                self.assertEqual(proto, mine["observed_len"], (cigar, pos, start, end))
                checked += 1
            else:
                self.assertIsNone(proto, (cigar, pos, start, end))
        self.assertGreater(checked, 100)

    def test_independent_verifier_cigar_agrees_with_implementation(self):
        rng = random.Random(7)
        for _ in range(3000):
            ops = [(rng.randint(1, 30), rng.choice("MMMMIDS")) for _ in range(rng.randint(1, 6))]
            cigar = "".join(f"{n}{o}" for n, o in ops)
            pos, start = rng.randint(1, 120), rng.randint(60, 110)
            end = start + rng.randint(8, 30)
            flag = rng.choice([0, 16, 1024, 256, 4])
            a = e6.analyse_alignment(flag, pos, 42, cigar, start, end)
            b = v6.analyse(["q", str(flag), "chr1", str(pos), "42", cigar], start, end)
            self.assertEqual(a is None, b is None, cigar)
            if a:
                self.assertEqual((a["observed_len"], a["flank_before"], a["flank_after"], a["offset"], a["read_len"],
                                  a["excluded_flag"], a["reverse"]),
                                 (b["obs"], b["fb"], b["fa"], b["off"], b["rl"], b["excl"], b["rev"]))

    def test_summary_counts_and_statuses(self):
        reads = ([read(name=f"a{i}", rev=i % 3 == 0) for i in range(30)] + [read(obs=10, name=f"d{i}") for i in range(6)]
                 + [read(fb=3, name=f"f{i}") for i in range(3)] + [read(mapq=5, name=f"m{i}") for i in range(2)])
        s = e6.summarise_locus_sample(reads, 4, 12)
        self.assertEqual((s["n_overlap_alignments"], s["n_excluded_flag"], s["n_excluded_mapq"], s["n_pass_filters"],
                          s["n_insufficient_flank"], s["n_informative"]), (45, 4, 2, 39, 3, 36))
        self.assertEqual(s["status"], e6.FS_ASSESSED)
        self.assertEqual(s["length_hist"], "10:6;12:30")
        self.assertEqual((s["modal_length"], round(s["modal_fraction"], 6)), (12, round(30 / 36, 6)))
        self.assertEqual(s["n_informative_fwd"] + s["n_informative_rev"], 36)

    def test_insufficient_missing_and_never_zero_stability(self):
        few = e6.summarise_locus_sample([read(name=f"x{i}") for i in range(5)], 0, 12)
        self.assertEqual(few["status"], e6.FS_INSUFFICIENT)
        self.assertIsNone(few["modal_length"])            # not 0, not "stable"
        self.assertIsNone(few["modal_fraction"])
        none = e6.summarise_locus_sample([], 0, 12)
        self.assertEqual(none["status"], e6.FS_MISSING)
        self.assertIsNone(none["callable_fraction"])       # undefined, not 0.0
        self.assertIsNone(none["length_hist"])

    def test_tie_break_is_the_shorter_length_and_tvd(self):
        st = e6.allele_stats({11: 5, 12: 5}, 12)
        self.assertEqual(st["modal_length"], 11)
        self.assertAlmostEqual(e6.tvd({12: 10}, {12: 10}), 0.0)
        self.assertAlmostEqual(e6.tvd({12: 10}, {11: 10}), 1.0)
        self.assertAlmostEqual(e6.tvd({12: 5, 11: 5}, {12: 10}), 0.5)

    def test_sensitivity_is_descriptive_counts_only(self):
        loci = {"a": [read(name=f"r{i}") for i in range(25)], "b": [read(name=f"s{i}", fb=2) for i in range(25)]}
        out = e6.sensitivity({"f1": loci}, {"a": 12, "b": 12})
        d = out["f1"]["n_loci_with_at_least_depth_at_default_mapq_flank"]
        self.assertEqual(d["20"], 1)                         # only locus a at flank 5
        self.assertEqual(out["f1"]["n_loci_meeting_provisional_depth_by_mapq_flank"]["mapq20_flank1"], 2)
        self.assertNotIn("best", json.dumps(out).lower())


class TestRegistry(unittest.TestCase):
    def base(self):
        return {"off_scope_policy": "NOT_ASSESSED_residual_reads_ignored",
                "biological_units": {"U1": {"label": "MSS", "label_source": "external", "assay": "WES", "dataset_origin": "X"}},
                "files": [{"file_id": "f1", "biological_unit_id": "U1", "role": "tumor", "order_id": "o1", "bam": "/a.bam"},
                          {"file_id": "f2", "biological_unit_id": "U1", "role": "normal", "order_id": "o1", "bam": "/b.bam"}]}

    def load(self, reg):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "r.json"
            p.write_text(json.dumps(reg))
            return e6.load_registry(str(p))

    def test_valid(self):
        self.assertEqual(len(self.load(self.base())["files"]), 2)

    def test_rejections(self):
        for mutate in (lambda r: r["files"][1].update(file_id="f1"),
                       lambda r: r["files"][0].update(biological_unit_id="U9"),
                       lambda r: r["files"][0].update(role="tumour"),
                       lambda r: r["files"][1].update(bam="/a.bam"),
                       lambda r: r["biological_units"]["U1"].update(label="maybe"),
                       lambda r: r["biological_units"]["U1"].pop("label_source"),
                       lambda r: r["files"][1].update(role="tumor"),
                       lambda r: r.pop("off_scope_policy")):
            reg = copy.deepcopy(self.base())
            mutate(reg)
            with self.assertRaises(SystemExit):
                self.load(reg)

    def test_real_registry_is_valid_and_one_unit(self):
        reg = e6.load_registry(str(DEV_DIR / "stage6a_dataset_registry.v1.json"))
        self.assertEqual({f["biological_unit_id"] for f in reg["files"]}, {"HCC1395"})
        self.assertEqual(reg["analysis_scope_contigs"], ["chr17"])
        self.assertEqual(reg["biological_units"]["HCC1395"]["label"], "MSS")


class TestUpstreamPin(unittest.TestCase):
    def make(self, d, **over):
        out = Path(d) / "s5_out.tsv"
        out.write_text("x\n")
        m = {"final_panel_size_fixed": False, "panel_size_selected": False, "membership_files_written": False,
             "hcc1395_results_used": False, "classification_performed": False, "production_changes": False,
             "upstream_outputs_modified": False,
             "outputs": {"o": {"path": str(out), "sha256": hashlib.sha256(b"x\n").hexdigest()}},
             "stage4_manifest_path": "nowhere", "verified_upstream_hashes": {"stage4_manifest": "0" * 64}}
        m.update(over)
        p = Path(d) / "m5.json"
        p.write_text(json.dumps(m))
        return p, out

    def test_pin_mismatch_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p, _ = self.make(d)
            with self.assertRaises(SystemExit) as cm:
                e6.verify_stage5(str(p), "f" * 64)
            self.assertIn("pinned", str(cm.exception))

    def test_safety_flag_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p, _ = self.make(d, panel_size_selected=True)
            with self.assertRaises(SystemExit):
                e6.verify_stage5(str(p), hashlib.sha256(p.read_bytes()).hexdigest())

    def test_tampered_stage5_output_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self.make(d)
            pin = hashlib.sha256(p.read_bytes()).hexdigest()
            out.write_text("tampered\n")
            with self.assertRaises(SystemExit) as cm:
                e6.verify_stage5(str(p), pin)
            self.assertIn("recorded hash", str(cm.exception))

    def test_default_pin_is_the_approved_stage5_manifest(self):
        self.assertEqual(e6.EXPECTED_STAGE5_MANIFEST_SHA256, v6.PINNED_STAGE5)
        self.assertTrue(v6.PINNED_STAGE5.startswith("9a5ecebf"))


class TestHashChain(unittest.TestCase):
    """A synthetic 5-level chain (stage 5 -> 4 -> 3 -> 2 -> 250 bp) with Windows-style backslash paths,
    exactly as the immutable real manifests record them."""

    F5 = ("final_panel_size_fixed", "panel_size_selected", "membership_files_written", "hcc1395_results_used",
          "classification_performed", "production_changes", "upstream_outputs_modified")
    F4 = ("final_panel_size_fixed", "panel_size_selected", "hcc1395_results_used", "classification_performed",
          "production_changes", "ranking_algorithm_rewritten", "stage3_outputs_modified", "chromosome_is_a_ranking_key",
          "chromosome_quotas", "chrY_removed")
    F3 = ("hcc1395_results_used", "classification_performed", "production_changes", "final_panel_size_fixed",
          "panel_size_selected", "stage2_outputs_modified")

    def build(self, d):
        d = Path(d)
        bs = chr(92)                       # Windows-style separator, as the real manifests record

        def data(name, text):
            (d / name).write_text(text)
            return {"path": d.as_posix() + bs + name, "sha256": hashlib.sha256(text.encode()).hexdigest()}

        def manifest(name, obj):
            (d / name).write_text(json.dumps(obj))
            return d.as_posix() + bs + name, hashlib.sha256((d / name).read_bytes()).hexdigest()
        prior_path, prior_sha = manifest("m250.json", {"outputs": {"viewA": {"panel": data("p250.tsv", "p")}}})
        p2, h2 = manifest("m2.json", {"outputs": {"trace": data("t2.tsv", "t")},
                                      "prior_manifest": {"path": prior_path, "sha256": prior_sha}})
        p3, h3 = manifest("m3.json", {**{f: False for f in self.F3}, "outputs": {"ranked_union": data("u3.tsv", "u")},
                                      "stage2_manifest_path": p2, "inputs_verified_hashes": {"stage2_manifest": h2}})
        p4, h4 = manifest("m4.json", {**{f: False for f in self.F4}, "outputs": {"cmp": data("c4.tsv", "c")},
                                      "stage3_manifest_path": p3, "verified_upstream_hashes": {"stage3_manifest": h3}})
        p5, h5 = manifest("m5.json", {**{f: False for f in self.F5}, "outputs": {"hybrid": data("h5.tsv", "h")},
                                      "stage4_manifest_path": p4, "verified_upstream_hashes": {"stage4_manifest": h4}})
        return p5, h5

    def test_chain_verifies_with_backslash_paths_and_records_every_hash(self):
        with tempfile.TemporaryDirectory() as d:
            p5, h5 = self.build(d)
            m5, m4, m3, hashes = e6.verify_stage5(p5, h5)
            for key in ("stage5_manifest", "stage5_hybrid", "stage4_manifest", "stage4_cmp", "stage3_manifest",
                        "stage3_ranked_union", "stage2_manifest", "stage2_trace", "spacing250_manifest", "spacing250_viewA_panel"):
                self.assertIn(key, hashes)

    def test_tampering_at_any_level_is_refused(self):
        for victim in ("h5.tsv", "c4.tsv", "u3.tsv", "t2.tsv", "p250.tsv", "m4.json", "m3.json", "m2.json", "m250.json"):
            with tempfile.TemporaryDirectory() as d:
                p5, h5 = self.build(d)
                f = Path(d) / victim
                f.write_text(f.read_text() + " ")
                with self.assertRaises(SystemExit, msg=victim):
                    e6.verify_stage5(p5, h5)

    def test_safety_flag_at_a_lower_level_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p5, h5 = self.build(d)
            m4 = json.loads((Path(d) / "m4.json").read_text())
            m4["panel_size_selected"] = True
            (Path(d) / "m4.json").write_text(json.dumps(m4))
            with self.assertRaises(SystemExit):
                e6.verify_stage5(p5, h5)          # also a hash mismatch: refused either way


# --------------------------------------------------------------------------- synthetic BAM fixtures
REF_LEN = {"chr1": 1000, "chr2": 800}


def build_reference(d):
    seqs = {}
    for name, n in REF_LEN.items():
        seq = list(("ACGTGCATCG" * (n // 10 + 1))[:n])
        seqs[name] = seq
    def run_(chrom, s, e, base):
        seqs[chrom][s - 1] = "C" if base != "C" else "G"
        seqs[chrom][e] = "C" if base != "C" else "G"
        for i in range(s, e):
            seqs[chrom][i] = base
    run_("chr1", 100, 112, "A")
    run_("chr1", 300, 315, "T")
    run_("chr1", 500, 512, "A")
    run_("chr2", 200, 214, "A")
    fa = Path(d) / "ref.fa"
    with open(fa, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            s = "".join(seq)
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    subprocess.run([SAMTOOLS, "faidx", str(fa)], check=True)
    return str(fa)


def sam_line(name, flag, chrom, pos, mapq, cigar):
    qlen = sum(int(n) for n, o in __import__("re").findall(r"(\d+)([MIDNSHP=X])", cigar) if o in "MIS=X") if cigar != "*" else 40
    return "\t".join([name, str(flag), chrom, str(pos), str(mapq), cigar, "*", "0", "0", "A" * qlen, "I" * qlen, "RG:Z:rg1"])


def locus1_reads(tag):
    """chr1 [100,112): the designed mixture described in the test docstring."""
    out, i = [], 0
    def add(n, flag, pos, mapq, cigar):
        nonlocal i
        for _ in range(n):
            i += 1
            out.append(sam_line(f"{tag}_{i:04d}", flag, "chr1", pos, mapq, cigar))
    add(20, 0, 81, 60, "40M")            # spanning, forward
    add(10, 16, 81, 60, "40M")           # spanning, reverse
    add(6, 0, 81, 60, "20M2D18M")        # deletion -> 10
    add(4, 0, 81, 60, "20M2I20M")        # insertion -> 14
    add(3, 0, 98, 60, "40M")             # only 3 bp of left flank -> insufficient flank
    add(2, 0, 81, 5, "40M")              # low mapq
    add(2, 1024, 81, 60, "40M")          # duplicate
    add(1, 256, 81, 60, "40M")           # secondary
    add(1, 2048, 81, 60, "40M")          # supplementary
    add(1, 4, 101, 0, "*")               # unmapped, placed
    return out


def locus2_reads(tag, n=8):
    return [sam_line(f"{tag}_L2_{i:03d}", 0, "chr1", 281, 60, "40M") for i in range(n)]   # 8 reads -> INSUFFICIENT


def write_bam(d, name, records, contigs=None, ref="/synthetic/ref.fa"):
    contigs = contigs or REF_LEN
    sam = Path(d) / f"{name}.sam"
    hdr = ["@HD\tVN:1.6\tSO:coordinate"] + [f"@SQ\tSN:{c}\tLN:{n}" for c, n in contigs.items()]
    hdr.append("@RG\tID:rg1\tSM:s1")
    hdr.append(f"@PG\tID:bwa\tPN:bwa\tCL:bwa mem -M {ref} r1.fq r2.fq")
    key = {c: i for i, c in enumerate(contigs)}
    records = sorted(records, key=lambda l: (key[l.split("\t")[2]], int(l.split("\t")[3])))
    sam.write_text("\n".join(hdr + records) + "\n")
    bam = Path(d) / f"{name}.bam"
    subprocess.run([SAMTOOLS, "view", "-b", "-o", str(bam), str(sam)], check=True, capture_output=True)
    subprocess.run([SAMTOOLS, "index", str(bam)], check=True, capture_output=True)
    return str(bam)


UNION_HEADER = ["marker_id", "source_view_class", "in_viewA_strong", "in_viewB_strong", "chrom", "start", "end", "span",
                "span_subband", "motif", "ref_repeat_length", "umap_status", "rmsk_self_category",
                "read_measurable_fraction", "K2_promega_span_distance", "rank_viewA", "rank_viewB",
                "rank_viewA_no_shared_key", "rank_viewB_no_shared_key"]


def union_rows(perturb=False):
    spec = [("L1", "A_and_B", "chr1", 100, 112, "A", 1, 1), ("L2", "A_and_B", "chr1", 300, 315, "T", 2, 2),
            ("L3", "B_only", "chr1", 500, 512, "A", None, 3), ("L4", "B_only", "chr2", 200, 214, "A", None, 4)]
    rows = []
    for mid, cls, chrom, s, e, motif, ra, rb in spec:
        rows.append({"marker_id": mid, "source_view_class": cls, "in_viewA_strong": str(ra is not None),
                     "in_viewB_strong": "True", "chrom": chrom, "start": str(s), "end": str(e), "span": str(e - s),
                     "span_subband": "31_40" if perturb else "21_30", "motif": motif, "ref_repeat_length": str(e - s),
                     "umap_status": "PERTURBED" if perturb else "TIER0_MEAN_1.0",
                     "rmsk_self_category": "no_rmsk_hit" if perturb else "microsatellite_class_only",
                     "read_measurable_fraction": "0.01" if perturb else "0.5323", "K2_promega_span_distance": "9" if perturb else "0",
                     "rank_viewA": "" if ra is None else str(ra), "rank_viewB": str(rb),
                     "rank_viewA_no_shared_key": "" if ra is None else str(ra), "rank_viewB_no_shared_key": str(rb)})
    return rows


def as_loaded(rows):
    out = []
    for r in copy.deepcopy(rows):
        for c in ("rank_viewA", "rank_viewB", "rank_viewA_no_shared_key", "rank_viewB_no_shared_key"):
            r[c] = int(r[c]) if r[c] != "" else None
        r["_shared"] = r["source_view_class"] == "A_and_B"
        r["_inA"], r["_inB"] = r["in_viewA_strong"] == "True", r["in_viewB_strong"] == "True"
        out.append(r)
    return out


class Fixture:
    """Builds reference, BAMs, registry, union table in a temp dir."""

    def __init__(self, d, dup_other_unit=False, second_unit=False, bad_contig=False, reads_out_of_scope=False):
        self.d = d
        self.fasta = build_reference(d)
        a = locus1_reads("A") + locus2_reads("A")
        if reads_out_of_scope:
            a = a + [sam_line("oos", 0, "chr2", 181, 60, "40M")]
        self.bam_t = write_bam(d, "tumor", a, ref=self.fasta)
        # normal: same designed locus-1 mixture, different read names (not a duplicate of the tumor)
        self.bam_n = write_bam(d, "normal", locus1_reads("N") + locus2_reads("N", 30), ref=self.fasta)
        self.bam_dup = write_bam(d, "tumor_dup", a, ref=self.fasta)      # identical alignments, another file
        contigs = dict(REF_LEN)
        if bad_contig:
            contigs["chr2"] = 801
        self.bam_bad = write_bam(d, "bad", locus1_reads("B"), contigs=contigs, ref=self.fasta)
        files = [{"file_id": "T1", "biological_unit_id": "U1", "role": "tumor", "order_id": "o1", "bam": self.bam_t},
                 {"file_id": "N1", "biological_unit_id": "U1", "role": "normal", "order_id": "o1", "bam": self.bam_n},
                 {"file_id": "T1dup", "biological_unit_id": "U2" if dup_other_unit else "U1", "role": "tumor",
                  "order_id": "o2", "bam": self.bam_dup}]
        if bad_contig:
            files = [{"file_id": "Bad", "biological_unit_id": "U1", "role": "tumor", "order_id": "o3", "bam": self.bam_bad}]
        units = {"U1": {"label": "MSS", "label_source": "synthetic external label", "label_method": "synthetic",
                        "assay": "synthetic", "dataset_origin": "SYN", "use": "development_only"}}
        if dup_other_unit or second_unit:
            units["U2"] = dict(units["U1"])
        self.registry = {"registry_version": "t", "reference": {"build": "GRCh38", "fasta": self.fasta},
                         "analysis_scope_contigs": ["chr1"], "biological_units": units, "files": files,
                         "off_scope_policy": "NOT_ASSESSED_residual_reads_ignored"}
        self.registry_path = Path(d) / "registry.json"
        self.registry_path.write_text(json.dumps(self.registry))
        self.fai = e6.read_fai(self.fasta + ".fai")

    def union(self, perturb=False):
        return as_loaded(union_rows(perturb))

    def run(self, out, perturb=False):
        return e6.run_analysis(self.union(perturb), self.registry, self.fai, SAMTOOLS, str(out), [2, 3], [1], [3])


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@needs_samtools
class TestSyntheticEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.fx = Fixture(cls.tmp.name)
        cls.out = Path(cls.tmp.name) / "out"
        cls.frag, cls.paths = cls.fx.run(cls.out)
        cls.ls = {(r["file_id"], r["marker_id"]): r for r in read_tsv(cls.paths["locus_sample_technical"])}
        cls.em = {r["candidate"]: r for r in read_tsv(cls.paths["evidence_matrix"])}
        cls.lts = {r["marker_id"]: r for r in read_tsv(cls.paths["locus_technical_summary"])}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_designed_locus1_counts(self):
        r = self.ls[("T1", "L1")]
        self.assertEqual((r["n_overlap_alignments"], r["n_excluded_flag"], r["n_excluded_mapq"], r["n_pass_filters"],
                         r["n_insufficient_flank"], r["n_informative"]), ("49", "4", "2", "43", "3", "40"))
        self.assertEqual(r["length_hist"], "10:6;12:30;14:4")
        self.assertEqual((r["modal_length"], r["modal_fraction"], r["n_informative_fwd"], r["n_informative_rev"]),
                         ("12", "0.750000", "30", "10"))
        self.assertEqual(r["status"], "ASSESSED")
        self.assertEqual(r["callable_fraction"], format(40 / 49, ".6f"))

    def test_duplicate_file_in_same_unit_is_flagged_and_not_double_counted(self):
        self.assertEqual(self.frag["duplicate_checks"]["T1dup"]["exact_duplicate_of"], "T1")
        self.assertEqual(self.frag["analysis_file_ids"], ["T1", "N1"])
        self.assertNotIn(("T1dup", "L1"), self.ls)

    def test_missing_and_insufficient_are_distinct_states(self):
        self.assertEqual(self.ls[("T1", "L2")]["status"], e6.FS_INSUFFICIENT)     # 8 reads
        self.assertEqual(self.ls[("N1", "L2")]["status"], e6.FS_ASSESSED)         # 30 reads
        self.assertEqual(self.ls[("T1", "L3")]["status"], e6.FS_MISSING)          # no reads
        self.assertEqual(self.ls[("T1", "L3")]["callable_fraction"], "NA")
        self.assertEqual(self.ls[("T1", "L2")]["modal_length"], "NA")             # never 0

    def test_evidence_matrix_statuses(self):
        l1, l2, l3, l4 = (self.em[k] for k in ("L1", "L2", "L3", "L4"))
        self.assertEqual((l1["coverage_evidence"], l1["stability_evidence"], l1["technical_evidence"]),
                         ("ASSESSED", "ASSESSED", "ASSESSED"))
        self.assertEqual((l2["coverage_evidence"], l2["stability_evidence"]), ("INSUFFICIENT", "INSUFFICIENT"))
        self.assertEqual(l2["technical_evidence"], "ASSESSED")                    # normal reaches depth
        self.assertEqual((l3["coverage_evidence"], l3["technical_evidence"]), ("INSUFFICIENT", "INSUFFICIENT"))
        self.assertEqual(l3["missingness"], "1.000000")
        self.assertEqual(l2["cov_tumor_informative_median"], "NA")                # not 0
        self.assertEqual(l2["stab_tn_tvd_median"], "NA")

    def test_out_of_scope_is_not_assessed_everywhere_and_never_zero(self):
        l4 = self.em["L4"]
        self.assertEqual(l4["data_scope"], e6.SCOPE_OUT)
        for c in e6.EM_EVIDENCE_COLS:
            self.assertEqual(l4[c], "NOT_ASSESSED", c)
        for c in e6.EM_SUPPORT_COLS:
            self.assertEqual(l4[c], "NA", c)
        self.assertNotIn("L4", self.lts)

    def test_msi_separation_never_assessed_and_reproducibility_insufficient(self):
        self.assertTrue(all(r["MSI_separation_evidence"] == "NOT_ASSESSED" for r in self.em.values()))
        for k in ("L1", "L2", "L3"):
            self.assertEqual(self.em[k]["reproducibility_evidence"], "INSUFFICIENT")
        self.assertEqual((self.em["L1"]["sample_count"], self.em["L1"]["ground_truth_count"],
                          self.em["L1"]["ground_truth_MSS_units"], self.em["L1"]["ground_truth_MSI_H_units"]),
                         ("1", "1", "1", "0"))

    def test_tumor_normal_tvd_and_summary(self):
        s = self.lts["L1"]
        self.assertEqual((s["tn_pairs_assessed"], s["tn_tvd_median"]), ("1", "0.000000"))   # same designed mixture
        self.assertEqual(s["tumor_modal_length_values"], "12")

    def test_scenarios_report_counts_only(self):
        sc = {r["scenario_id"]: r for r in read_tsv(self.paths["scenario_comparison"])}
        row = sc["viewB_quality_first_3"]
        self.assertEqual((row["n_loci"], row["n_in_scope"], row["n_out_of_scope_NOT_ASSESSED"]), ("3", "3", "0"))
        self.assertEqual(sc["viewB_quality_first_2"]["coverage_ASSESSED"], "1")
        self.assertIn("viewB_hybrid_core1_total3", sc)
        self.assertEqual(sc["viewB_hybrid_core1_total3"]["core_in_scope"], "1")
        header = list(sc["viewB_quality_first_3"].keys())
        self.assertFalse([c for c in header if "score" in c.lower() or "winner" in c.lower() or "best" in c.lower()])
        self.assertTrue(all(r["evidence_scope_caveat"].startswith("CHR17_ONLY") for r in sc.values()))

    def test_undefined_stage3_rank_is_written_as_NA_in_the_strat_column(self):
        rows = union_rows()
        rows[3]["rank_viewB_no_shared_key"] = ""        # a locus without a View-B rank (like the A-only loci)
        rows[3]["in_viewB_strong"] = "False"
        rows[3]["in_viewA_strong"] = "True"
        rows[3]["rank_viewA"] = rows[3]["rank_viewA_no_shared_key"] = "9"
        with tempfile.TemporaryDirectory() as d2:
            frag, paths = e6.run_analysis(as_loaded(rows), self.fx.registry, self.fx.fai, SAMTOOLS, d2, [2], [1], [3])
            em = {r["candidate"]: r for r in read_tsv(paths["evidence_matrix"])}
            self.assertEqual(em["L4"]["strat_rank_viewB_no_shared_key"], "NA")

    def test_no_membership_files_or_composite_score(self):
        names = {p.name for p in Path(self.out).iterdir()}
        self.assertFalse([n for n in names if "membership" in n.lower() or "panel" in n.lower()])
        for p in Path(self.out).iterdir():
            head = open(p, encoding="utf-8").readline().lower()
            self.assertNotIn("composite", head)

    def test_circularity_stage3_attributes_do_not_change_any_evidence_value(self):
        with tempfile.TemporaryDirectory() as d2:
            _, paths2 = self.fx.run(Path(d2), perturb=True)
            for name, cols in (("locus_sample_technical", None), ("locus_technical_summary", None), ("evidence_matrix", None)):
                a = read_tsv(self.paths[name])
                b = read_tsv(paths2[name])
                self.assertEqual(len(a), len(b))
                for ra, rb in zip(a, b):
                    for c in ra:
                        if c.startswith("strat_"):
                            continue
                        self.assertEqual(ra[c], rb[c], (name, c))
            a = read_tsv(self.paths["evidence_matrix"])
            b = read_tsv(paths2["evidence_matrix"])
            self.assertNotEqual([r["strat_umap_status"] for r in a], [r["strat_umap_status"] for r in b])   # they DO differ

    def test_evidence_columns_contain_no_stage3_attribute_names(self):
        forbidden = ("umap", "rmsk", "repeatmasker", "read_measurable", "span_distance", "K2", "K3", "K4", "K5", "K6", "K7",
                     "rank", "hash", "uniq")
        for c in e6.EM_EVIDENCE_COLS + e6.EM_SUPPORT_COLS + e6.LS_COLS:
            self.assertFalse(any(f.lower() in c.lower() for f in forbidden), c)
        self.assertTrue(all(c.startswith("strat_") for c in e6.STRAT_COLS))

    def test_verifier_passes_on_the_synthetic_run_and_detects_tampering(self):
        with tempfile.TemporaryDirectory() as d2:
            out = Path(d2) / "out"
            frag, paths = self.fx.run(out)
            union_path = Path(d2) / "union.tsv"
            with open(union_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=UNION_HEADER, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(union_rows())
            stage1 = Path(d2) / "scan_manifest.json"
            stage1.write_text("{}")
            mp = e6.write_manifest(out_dir=str(out), samtools=SAMTOOLS, registry_path=str(self.fx.registry_path),
                                   registry=self.fx.registry, fai_path=self.fx.fasta + ".fai", stage1_path=str(stage1),
                                   union_path=str(union_path), upstream={"stage5_manifest": "synthetic"},
                                   stage5_manifest_path="synthetic", expected_pin="synthetic", sizes=[2, 3],
                                   hybrid_cores=[1], hybrid_totals=[3], frag=frag, paths=paths)
            self.assertEqual(v6.verify(mp, skip_upstream=True), [])
            # tamper one cell AND re-hash so only the independent recomputation can catch it
            p = Path(paths["locus_sample_technical"])
            text = p.read_text().replace("\t10:6;12:30;14:4\t", "\t10:6;12:31;14:4\t", 1)
            p.write_text(text)
            m = json.loads(Path(mp).read_text())
            m["outputs"]["locus_sample_technical"]["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
            Path(mp).write_text(json.dumps(m))
            fails = v6.verify(mp, skip_upstream=True)
            self.assertTrue(any("length_hist" in f for f in fails), fails[:3])

    def test_verifier_detects_a_zero_substituted_not_assessed_cell(self):
        with tempfile.TemporaryDirectory() as d2:
            out = Path(d2) / "out"
            frag, paths = self.fx.run(out)
            union_path = Path(d2) / "union.tsv"
            with open(union_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=UNION_HEADER, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(union_rows())
            stage1 = Path(d2) / "s.json"
            stage1.write_text("{}")
            mp = e6.write_manifest(out_dir=str(out), samtools=SAMTOOLS, registry_path=str(self.fx.registry_path),
                                   registry=self.fx.registry, fai_path=self.fx.fasta + ".fai", stage1_path=str(stage1),
                                   union_path=str(union_path), upstream={}, stage5_manifest_path="s", expected_pin="s",
                                   sizes=[2, 3], hybrid_cores=[1], hybrid_totals=[3], frag=frag, paths=paths)
            p = Path(paths["evidence_matrix"])
            lines = p.read_text().splitlines()
            head = lines[0].split("\t")
            i = next(k for k, l in enumerate(lines) if l.startswith("L4\t"))
            cells = lines[i].split("\t")
            cells[head.index("missingness")] = "0.000000"
            lines[i] = "\t".join(cells)
            p.write_text("\n".join(lines) + "\n")
            m = json.loads(Path(mp).read_text())
            m["outputs"]["evidence_matrix"]["sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
            Path(mp).write_text(json.dumps(m))
            self.assertTrue(any("must be NOT_ASSESSED" in f for f in v6.verify(mp, skip_upstream=True)))


@needs_samtools
class TestIntegrityChecks(unittest.TestCase):
    def test_contig_length_mismatch_fails_loudly(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d, bad_contig=True)
            with self.assertRaises(SystemExit) as cm:
                fx.run(Path(d) / "o")
            self.assertIn("@SQ", str(cm.exception))

    def test_sparse_residual_reads_outside_scope_are_reported_not_assessed(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d, reads_out_of_scope=True)
            frag, paths = fx.run(Path(d) / "o")
            info = frag["reference_checks"]["T1"]
            self.assertEqual(info["mapped_reads_off_scope_primary_contigs"], 1)
            self.assertGreater(info["off_scope_read_fraction_of_primary"], 0)
            em = {r["candidate"]: r for r in read_tsv(paths["evidence_matrix"])}
            self.assertEqual(em["L4"]["technical_evidence"], "NOT_ASSESSED")   # a chr2 locus with a residual read is still not assessed

    def test_scope_not_holding_the_reads_fails_loudly(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d)
            fx.registry["analysis_scope_contigs"] = ["chr2"]
            with self.assertRaises(SystemExit) as cm:
                e6.check_reference_compatibility(SAMTOOLS, fx.registry, fx.fai, fx.registry["files"])
            self.assertIn("hold no reads", str(cm.exception))
            fx.registry["analysis_scope_contigs"] = ["chr1", "chr2"]     # chr2 has no reads in the tumor BAM
            with self.assertRaises(SystemExit) as cm:
                e6.check_reference_compatibility(SAMTOOLS, fx.registry, fx.fai, fx.registry["files"])
            self.assertIn("hold no reads", str(cm.exception))

    def test_contig_with_most_reads_outside_scope_fails_loudly(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d)
            recs = ([sam_line(f"c1_{i}", 0, "chr1", 81, 60, "40M") for i in range(5)]
                    + [sam_line(f"c2_{i}", 0, "chr2", 181, 60, "40M") for i in range(50)])
            bam = write_bam(d, "skewed", recs, ref=fx.fasta)
            files = [{"file_id": "S", "biological_unit_id": "U1", "role": "tumor", "order_id": "os", "bam": bam}]
            with self.assertRaises(SystemExit) as cm:
                e6.check_reference_compatibility(SAMTOOLS, fx.registry, fx.fai, files)
            self.assertIn("holding most reads", str(cm.exception).replace("contig holding most reads", "holding most reads"))

    def test_reference_file_not_named_in_bam_fails(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d)
            fx.registry["reference"]["fasta"] = "/somewhere/else.fa"
            with self.assertRaises(SystemExit) as cm:
                e6.check_reference_compatibility(SAMTOOLS, fx.registry, fx.fai, fx.registry["files"])
            self.assertIn("@PG", str(cm.exception))

    def test_locus_coordinate_and_reference_allele_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d)
            rows = fx.union()
            rows[0]["motif"] = "C"                       # the reference is a poly-A here
            with self.assertRaises(SystemExit) as cm:
                e6.check_locus_coordinates(SAMTOOLS, fx.fasta, fx.fai, rows, d)
            self.assertIn("mismatches", str(cm.exception))
            rows = fx.union()
            rows[1]["ref_repeat_length"] = "14"           # end - start = 15
            with self.assertRaises(SystemExit):
                e6.check_locus_coordinates(SAMTOOLS, fx.fasta, fx.fai, rows, d)
            rows = fx.union()
            rows[1]["chrom"] = "chr9"
            with self.assertRaises(SystemExit):
                e6.check_locus_coordinates(SAMTOOLS, fx.fasta, fx.fai, rows, d)

    def test_identical_content_under_different_biological_units_fails(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d, dup_other_unit=True)
            with self.assertRaises(SystemExit) as cm:
                fx.run(Path(d) / "o")
            self.assertIn("different biological units", str(cm.exception))

    def test_shared_read_names_under_different_units_fail(self):
        with tempfile.TemporaryDirectory() as d:
            fx = Fixture(d, second_unit=True)
            # same read names (a re-alignment of the same FASTQ) but different content, filed as another unit
            recs = [sam_line(f"A_{i:04d}", 0, "chr1", 85, 60, "40M") for i in range(1, 30)]
            other = write_bam(d, "other", recs, ref=fx.fasta)
            files = [fx.registry["files"][0], {"file_id": "X", "biological_unit_id": "U2", "role": "tumor",
                                               "order_id": "ox", "bam": other}]
            with self.assertRaises(SystemExit) as cm:
                e6.check_duplicates(files, SAMTOOLS, "chr1")
            self.assertIn("shared read names", str(cm.exception))


if __name__ == "__main__":
    unittest.main()

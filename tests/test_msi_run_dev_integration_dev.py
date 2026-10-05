"""Synthetic-only end-to-end wiring checks for the development MSI runner.

Marker coordinates and SAM records are fabricated. The samtools process
boundary is simulated, while SAM parsing, read filtering, CIGAR repeat-length
extraction, observation construction, classification, and report writing are
the real development code paths. No BAM or biological data are used.
"""

import contextlib
import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

DEV_DIR = Path(__file__).resolve().parents[1] / "dev" / "msi_redesign"
sys.path.insert(0, str(DEV_DIR))

import run_msi_dev  # noqa: E402


class TestSyntheticRunnerIntegration(unittest.TestCase):
    def _marker_file(self, directory):
        marker_path = Path(directory) / "synthetic_markers.tsv"
        with marker_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["marker_id", "chrom", "start", "end", "motif",
                             "ref_repeat_length", "reference_build", "gene", "kind"])
            for i in range(5):
                start = 100 + i * 100
                writer.writerow([f"M{i}", "synthetic_contig", start, start + 10,
                                 "A", 10, "GRCh38", "SYNTHETIC", "mononucleotide"])
        return marker_path

    def _run_with_tumor_lengths(self, tumor_lengths_by_locus):
        temp = tempfile.TemporaryDirectory(prefix="msi-runner-synthetic-")
        self.addCleanup(temp.cleanup)
        marker_path = self._marker_file(temp.name)

        def sam_lines(length, marker_index):
            start = 100 + marker_index * 100
            pos = start - 4  # 0-based start-5 flank -> SAM 1-based POS
            cigar = "5M10M5M" if length == 10 else "5M1I10M5M"
            return "".join(
                f"read{i}\t0\tsynthetic_contig\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*\n"
                for i in range(20)
            )

        def simulated_samtools(command, **_kwargs):
            # Expected invocation: samtools view <bam> <chrom:start-end>
            self.assertEqual(command[1], "view")
            bam_path = command[2]
            region = command[3]
            start_1based = int(region.rsplit(":", 1)[1].split("-", 1)[0])
            marker_start = start_1based - 1
            marker_index = (marker_start - 100) // 100
            if bam_path == "tumor.synthetic.bam":
                read_lengths = tumor_lengths_by_locus[marker_index]
            elif bam_path == "normal.synthetic.bam":
                read_lengths = [10] * 20
            else:
                self.fail(f"unexpected fabricated BAM path: {bam_path}")
            self.assertEqual(len(read_lengths), 20)
            self.assertTrue(set(read_lengths).issubset({10, 11}))
            stdout = "".join(sam_lines(length, marker_index) for length in read_lengths)
            return type("Completed", (), {"returncode": 0, "stdout": stdout, "stderr": ""})()

        output = io.StringIO()
        with patch.object(run_msi_dev.be.subprocess, "run", side_effect=simulated_samtools):
            with contextlib.redirect_stdout(output):
                run_msi_dev.main([
                    "--markers", str(marker_path),
                    "--tumor-bam", "tumor.synthetic.bam",
                    "--normal-bam", "normal.synthetic.bam",
                    "--out-dir", temp.name,
                    "--sample-id", "SYNTHETIC",
                ])

        summary_path = Path(temp.name) / "SYNTHETIC.msi_dev.summary.json"
        per_locus_path = Path(temp.name) / "SYNTHETIC.msi_dev.per_locus.tsv"
        return json.loads(summary_path.read_text(encoding="utf-8")), per_locus_path

    def test_marker_tsv_to_classifier_reports_synthetic_msi_h(self):
        tumor = [[11] * 20, [11] * 20, [10] * 20, [10] * 20, [10] * 20]
        summary, per_locus_path = self._run_with_tumor_lengths(tumor)
        self.assertEqual(summary["sample_status"], "VALID_MSI_H")
        self.assertEqual((summary["n_markers_in_panel"], summary["n_callable"],
                          summary["n_unstable"]), (5, 5, 2))
        with per_locus_path.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 5)
        self.assertEqual([r["status"] for r in rows],
                         ["UNSTABLE", "UNSTABLE", "STABLE", "STABLE", "STABLE"])

    def test_marker_tsv_to_classifier_reports_synthetic_mss(self):
        tumor = [[10] * 20 for _ in range(5)]
        summary, _ = self._run_with_tumor_lengths(tumor)
        self.assertEqual(summary["sample_status"], "VALID_MSS")
        self.assertEqual((summary["n_callable"], summary["n_unstable"]), (5, 0))

    def test_mixed_repeat_lengths_are_independent_of_read_order(self):
        mixed_forward = [10] * 10 + [11] * 10
        mixed_reverse = [11] * 10 + [10] * 10
        tumor_forward = [mixed_forward] + [[10] * 20 for _ in range(4)]
        tumor_reverse = [mixed_reverse] + [[10] * 20 for _ in range(4)]

        forward_summary, forward_path = self._run_with_tumor_lengths(tumor_forward)
        reverse_summary, reverse_path = self._run_with_tumor_lengths(tumor_reverse)

        self.assertEqual(forward_summary["sample_status"], reverse_summary["sample_status"])
        self.assertEqual(forward_summary["n_unstable"], reverse_summary["n_unstable"])
        with forward_path.open(encoding="utf-8") as handle:
            forward_rows = list(csv.DictReader(handle, delimiter="\t"))
        with reverse_path.open(encoding="utf-8") as handle:
            reverse_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(forward_rows[0]["status"], reverse_rows[0]["status"])


if __name__ == "__main__":
    unittest.main()

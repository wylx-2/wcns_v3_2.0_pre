#!/usr/bin/env python3
"""Regression tests for the bounded local-probe checks in analyze_rans_run.py."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class RansAnalyzerTests(unittest.TestCase):
    analyzer: Path

    def analyze(self, maximum_peak_factor: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            history = output / "probe.history.r1.txt"
            history.write_text(
                "# step wall_time total_l2 accepted_dt k_floor_repairs "
                "omega_floor_repairs stop_reason\n"
                "1 1 10 1 0 0 running\n"
                "2 2 14 1 0 0 running\n"
                "3 3 13 1 0 0 maximum_steps\n",
                encoding="utf-8",
            )
            return subprocess.run(
                [
                    sys.executable,
                    str(self.analyzer),
                    "--output",
                    str(output),
                    "--cell-count",
                    "8",
                    "--minimum-residual-decades",
                    "-1",
                    "--last-steps-without-projection",
                    "2",
                    "--maximum-residual-peak-factor",
                    maximum_peak_factor,
                    "--maximum-final-residual-factor",
                    "1.5",
                    "--maximum-wall-time",
                    "10",
                    "--maximum-accepted-steps",
                    "3",
                ],
                capture_output=True,
                text=True,
                check=False,
            )

    def test_accepts_bounded_transient(self) -> None:
        result = self.analyze("1.5")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('"residual_peak_factor": 1.4', result.stdout)

    def test_rejects_excessive_peak(self) -> None:
        result = self.analyze("1.3")
        self.assertEqual(result.returncode, 1)
        self.assertIn('"residual_peak_factor": false', result.stdout)


def main() -> int:
    arguments = argparse.ArgumentParser()
    arguments.add_argument("--analyzer", type=Path, required=True)
    args, unused = arguments.parse_known_args()
    RansAnalyzerTests.analyzer = args.analyzer.resolve()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(RansAnalyzerTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())

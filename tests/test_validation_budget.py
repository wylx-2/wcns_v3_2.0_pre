#!/usr/bin/env python3
"""CLI regression tests for tools/check_validation_budget.py."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class ValidationBudgetTests(unittest.TestCase):
    checker: Path

    def run_checker(self, config_text: str, *extra: str) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / "case.wcns"
            config.write_text(config_text, encoding="utf-8")
            command = [
                sys.executable,
                str(self.checker),
                "--config",
                str(config),
                "--kind",
                "probe",
                "--dimension",
                "2",
                "--cells",
                "14000",
                *extra,
            ]
            return subprocess.run(command, capture_output=True, text=True, check=False)

    def test_accepts_bounded_probe(self) -> None:
        result = self.run_checker(
            "run.max_steps = 500\nrun.max_wall_time = 900\n"
            "output.checkpoint.enabled = true\n",
            "--phase-wall-seconds",
            "1800",
            "--estimated-memory-gib",
            "1.5",
            "--planned-artifacts-gib",
            "0.25",
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('"status": "PASS"', result.stdout)

    def test_rejects_unbounded_config_and_resource_request(self) -> None:
        result = self.run_checker(
            "run.max_steps = 5000\nrun.max_wall_time = 0\n"
            "output.checkpoint.enabled = false\n",
            "--phase-wall-seconds",
            "7201",
            "--concurrent-solvers",
            "3",
            "--estimated-memory-gib",
            "4.1",
            "--planned-artifacts-gib",
            "1.1",
        )
        self.assertEqual(result.returncode, 2)
        for fragment in (
            "run.max_steps",
            "run.max_wall_time",
            "output.checkpoint.enabled",
            "phase_wall_seconds",
            "concurrent_solvers",
            "estimated_memory_gib",
            "planned_artifacts_gib",
        ):
            self.assertIn(fragment, result.stdout)

    def test_rejects_oversized_grid(self) -> None:
        result = self.run_checker(
            "run.max_steps = 20\nrun.max_wall_time = 60\n"
            "output.checkpoint.enabled = true\n",
            "--cells",
            "32769",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("cells=32769", result.stdout)

    def test_rejects_duplicate_config_key(self) -> None:
        result = self.run_checker(
            "run.max_steps = 20\nrun.max_steps = 10\nrun.max_wall_time = 60\n"
            "output.checkpoint.enabled = true\n"
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("duplicate key", result.stdout)


def main() -> int:
    arguments = argparse.ArgumentParser()
    arguments.add_argument("--checker", type=Path, required=True)
    args, unittest_args = arguments.parse_known_args()
    ValidationBudgetTests.checker = args.checker.resolve()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ValidationBudgetTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())

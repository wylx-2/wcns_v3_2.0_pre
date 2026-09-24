#!/usr/bin/env python3
"""Regression tests for Case06 parameter conversion and config materialization."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class Case06PreparationTests(unittest.TestCase):
    generator: Path
    manifest_path: Path

    def run_generator(self, *extra: str) -> tuple[subprocess.CompletedProcess[str], str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            grid = root / "grid.cgns"
            grid.write_bytes(b"test fixture: generator checks identity later in dry-run")
            config = root / "case.wcns"
            command = [
                sys.executable,
                str(self.generator),
                "--grid", str(grid),
                "--output-directory", str(root / "output"),
                "--config", str(config),
                "--case-name", "case06-test",
                *extra,
            ]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            text = config.read_text(encoding="utf-8") if config.is_file() else ""
            return result, text

    def test_manifest_conversion_identities(self) -> None:
        manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        physics = manifest["physics"]
        scales = manifest["normalization"]
        model = manifest["models"]["sst_2003m"]
        mach = physics["mach"]
        velocity = mach * math.sqrt(
            physics["gamma"] * physics["specific_gas_constant"]
            * physics["reference_temperature_kelvin"]
        )
        self.assertAlmostEqual(scales["reference_velocity"], velocity, places=13)
        reynolds = scales["reference_density"] * velocity * scales["reference_length"] \
            / scales["reference_viscosity"]
        self.assertAlmostEqual(reynolds, physics["reynolds_chord"], places=8)

        intensity = model["freestream_turbulence_intensity"]
        k = 1.5 * intensity * intensity
        omega = physics["reynolds_chord"] * k \
            / model["freestream_eddy_viscosity_ratio"]
        length_scale = math.sqrt(k) / (0.09 ** 0.25 * omega)
        self.assertAlmostEqual(model["freestream_k_over_u_squared"], k, places=20)
        self.assertAlmostEqual(model["freestream_omega_length_over_u"], omega, places=13)
        self.assertAlmostEqual(
            model["freestream_length_scale_over_chord"], length_scale, places=20
        )
        self.assertAlmostEqual(k * mach * mach, 9.126e-9, places=20)
        self.assertAlmostEqual(omega / physics["reynolds_chord"] * mach * mach,
                               1.014e-6, places=18)

    def test_materializes_sst_angle_and_local_limits(self) -> None:
        result, text = self.run_generator(
            "--model", "sst_2003m", "--mode", "local", "--angle-deg", "0",
            "--max-steps", "20", "--max-wall-time", "300",
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("initial.u = 1", text)
        self.assertIn("initial.v = 0", text)
        self.assertIn("turbulence.freestream.intensity = 0.00052", text)
        self.assertIn("turbulence.freestream.length_scale = 4.300130725961134e-6", text)
        self.assertIn("run.max_steps = 20", text)
        self.assertIn("preconditioner.type = weiss_smith", text)

    def test_materializes_server_sa_reference_point_vortex(self) -> None:
        result, text = self.run_generator(
            "--model", "sa_neg", "--mode", "server", "--angle-deg", "10",
            "--max-steps", "1000", "--max-wall-time", "3600",
            "--reference-point-vortex", "--preconditioner", "none",
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("boundary.point_vortex.enabled = true", text)
        self.assertIn("boundary.point_vortex.lift_coefficient = 1.0909146672", text)
        self.assertIn("preconditioner.type = none", text)
        self.assertNotIn("preconditioner.mach_cutoff", text)
        self.assertNotIn("preconditioner.viscous_cutoff", text)

    def test_rejects_unbounded_or_misclassified_runs(self) -> None:
        local, _ = self.run_generator(
            "--model", "sa_neg", "--mode", "local", "--max-steps", "21",
        )
        self.assertNotEqual(local.returncode, 0)
        self.assertIn("1..20", local.stderr)
        server, _ = self.run_generator(
            "--model", "sst_2003m", "--mode", "server", "--max-steps", "1000",
        )
        self.assertNotEqual(server.returncode, 0)
        self.assertIn("--max-wall-time", server.stderr)
        point_vortex, _ = self.run_generator(
            "--model", "sa_neg", "--mode", "local", "--reference-point-vortex",
        )
        self.assertNotEqual(point_vortex.returncode, 0)
        self.assertIn("server-only", point_vortex.stderr)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args, unittest_args = parser.parse_known_args()
    Case06PreparationTests.generator = args.generator.resolve()
    Case06PreparationTests.manifest_path = args.manifest.resolve()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Case06PreparationTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())

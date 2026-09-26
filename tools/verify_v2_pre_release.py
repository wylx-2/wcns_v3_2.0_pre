#!/usr/bin/env python3
"""Check the source/document contracts required by the v2.0_pre preview."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def require_text(root: Path, relative: str, fragments: tuple[str, ...]) -> int:
    path = root / relative
    if not path.is_file():
        raise RuntimeError(f"missing required file: {relative}")
    text = path.read_text(encoding="utf-8")
    missing = [fragment for fragment in fragments if fragment not in text]
    if missing:
        raise RuntimeError(f"{relative} lacks required contracts: {missing}")
    return len(fragments)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.repository.resolve()
    checks = 0

    cmake = (root / "CMakeLists.txt").read_text(encoding="utf-8")
    contracts = (
        r"project\s*\(\s*wcns\s+VERSION\s+2\.0\.0",
        r'set\s*\(\s*WCNS_PROGRAM_VERSION\s+"2\.0_pre"',
        r'option\s*\(\s*WCNS_BUILD_TESTS\s+"[^"]+"\s+OFF\s*\)',
        r'option\s*\(\s*WCNS_INSTALL_EXAMPLES\s+"[^"]+"\s+OFF\s*\)',
    )
    for contract in contracts:
        if re.search(contract, cmake, flags=re.IGNORECASE) is None:
            raise RuntimeError(f"CMake release contract is missing: {contract}")
        checks += 1

    checks += require_text(
        root,
        "src/solver/turbulence_model.cpp",
        (
            'name == "sa_neg"',
            'name == "k_omega_sst"',
            'name == "k_epsilon"',
            'name == "smagorinsky"',
            'name == "scale_similarity"',
            'name == "mixed_smagorinsky_similarity"',
            'name == "dynamic_smagorinsky"',
            'name == "wale"',
            '"sgs_energy_transfer"',
            '"les_dynamic_coefficient"',
            '"les_grid_anisotropy"',
        ),
    )
    checks += require_text(
        root,
        "src/runtime/case_config.cpp",
        (
            '"time.integrator"',
            '"preconditioner.type"',
            '"statistics.time.weight"',
            '"output.boundary.span_bin_edges"',
            '"output.statistics.xz_planes.enabled"',
            '"output.statistics.yz_planes.enabled"',
            '"output.statistics.channel_walls.enabled"',
        ),
    )
    checks += require_text(
        root,
        "src/runtime/quantity_registry.cpp",
        (
            '"total_mass"',
            '"total_momentum_x"',
            '"total_momentum_y"',
            '"total_momentum_z"',
            '"total_energy"',
            '"xz_mean_u_j"',
            '"yz_mean_u"',
            '"channel_wall_shear_lower"',
            '"channel_friction_velocity"',
            '"channel_re_tau"',
        ),
    )
    checks += require_text(
        root,
        "src/runtime/weighted_statistics.cpp",
        (
            "WeightedMomentState::add",
            "WeightedCovarianceState::add",
            "AcceptedTimeStatistics::sample",
            '"wcns_weighted_statistics_v2"',
        ),
    )
    checks += require_text(
        root,
        "src/runtime/boundary_output.cpp",
        (
            '"Cp"',
            '"Cf"',
            '"q_wall"',
            '"friction_velocity"',
            '"wall_y_plus"',
            '"pressure_traction_"',
            '"viscous_traction_"',
            '"traction_"',
            '".spanwise_loads.r"',
        ),
    )
    checks += require_text(
        root,
        "算法补充.md",
        (
            "### 12.1 Favre 平均 RANS 系统",
            "### 12.7 Weiss--Smith 型低 Mach 预处理",
            "### 12.9 壁面量与时间统计",
            "### 12.12 可压缩 Favre 滤波 LES 系统",
            "### 12.13 Smagorinsky、尺度相似、动态与 WALE 模型",
            "### 12.15 LU-SGS 定常与非定常隐式推进",
            "### 12.16 LES 与 LU-SGS 主要规格来源",
        ),
    )
    checks += require_text(
        root,
        "docs/v2.0-pre-capability-matrix.md",
        (
            "物理待验证",
            "accepted_dt",
            "Reynolds/Favre covariance",
            "展向载荷",
            "明确未内建的统计功能",
            "目标翼型未验证",
        ),
    )
    checks += require_text(
        root,
        "tools/package_v2_pre.py",
        (
            'PACKAGE_NAME = "WCNS_v2.0_pre"',
            '"cases"',
            '"examples"',
            '"tests"',
            '".cgns"',
            "require_clean_payload(entries)",
            "PACKAGE_CONTENTS.sha256",
            "verify_markdown_links(directory)",
        ),
    )

    print(f"v2.0_pre release contracts verified: {checks} checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare the current MDCD branches with ERCOFTAC Case 032 / MKM DNS."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_mdcd_dissipation_branches import (
    concatenate_tables,
    pooled_profiles,
    read_table,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
PROFILE_DIRECTORY = RESULTS / "branch-mdcd-dissipation-comparison" / "profiles"
OUTPUT = RESULTS / "branch-mdcd-dissipation-comparison" / "dns-case032-comparison"
REFERENCE = ROOT / "reference" / "ercoftac_case032_mkm"
REFERENCE_RE_TAU = 180.0
PROFILE_TIMES = (250, 255, 260, 265, 270, 275)
CASES = {
    "diss=0.001": (
        "branch_scmm6_mdcd0p001_roe_t223",
        [
            RESULTS / "branch-scmm6_mdcd0p001_roe-t223",
            RESULTS / "branch-scmm6_mdcd0p001_roe-t223_segment02",
        ],
    ),
    "diss=0.01": (
        "branch_scmm6_mdcd0p01_roe_t223",
        [
            RESULTS / "branch-scmm6_mdcd0p01_roe-t223",
            RESULTS / "branch-scmm6_mdcd0p01_roe-t223_segment02",
        ],
    ),
}


def load_dns() -> dict[str, np.ndarray | float]:
    lines = (REFERENCE / "simul1.dat").read_text(encoding="utf-8").splitlines()
    marker = lines.index("Mean and mean-square fluctuations:")
    rows = []
    for line in lines[marker + 1 :]:
        fields = line.split()
        if len(fields) == 8 and fields[0].isdigit():
            rows.append([float(value) for value in fields])
            if len(rows) == 65:
                break
    if len(rows) != 65:
        raise ValueError(f"expected 65 primary DNS profile rows, found {len(rows)}")
    data = np.asarray(rows)
    y = data[:, 1]
    u_plus = data[:, 3]
    bulk_plus = float(np.trapezoid(u_plus, y))
    return {
        "eta": y,
        "y_plus": data[:, 2],
        "mean_u_plus": u_plus,
        "rms_u_plus": np.sqrt(np.maximum(0.0, data[:, 4])),
        "rms_v_plus": np.sqrt(np.maximum(0.0, data[:, 5])),
        "rms_w_plus": np.sqrt(np.maximum(0.0, data[:, 6])),
        "minus_uv_plus": -data[:, 7],
        "bulk_velocity_plus": bulk_plus,
        "skin_friction_coefficient": 2.0 / bulk_plus**2,
    }


def find_profiles(stem: str) -> list[Path]:
    paths = []
    for time in PROFILE_TIMES:
        token = f"time2p{time - 200:02d}0000000eP02"
        matches = sorted(PROFILE_DIRECTORY.glob(f"*{stem}*{token}*.profile.txt"))
        if len(matches) != 1:
            raise RuntimeError(
                f"expected one profile for {stem} at t={time}, found {len(matches)}"
            )
        paths.append(matches[0])
    return paths


def fold_channel(profile: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    count = profile["y"].size
    if count % 2 != 0:
        raise ValueError("full-channel profile must have an even cell count")
    lower = np.arange(count // 2)
    upper = np.arange(count - 1, count // 2 - 1, -1)
    eta_lower = 1.0 + profile["y"][lower]
    eta_upper = 1.0 - profile["y"][upper]
    if not np.allclose(eta_lower, eta_upper):
        raise ValueError("wall-normal grid is not mirror symmetric")
    result = {
        "eta": 0.5 * (eta_lower + eta_upper),
        "mean_u": 0.5 * (profile["mean_u"][lower] + profile["mean_u"][upper]),
        "minus_uv": -0.5
        * (profile["reynolds_uv"][lower] - profile["reynolds_uv"][upper]),
    }
    for name in ("rms_u", "rms_v", "rms_w"):
        result[name] = np.sqrt(
            0.5 * (profile[name][lower] ** 2 + profile[name][upper] ** 2)
        )
    return result


def relative_l2(current: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(current - reference) / np.linalg.norm(reference))


def cell_widths_from_centers(eta: np.ndarray) -> np.ndarray:
    """Recover symmetric half-channel cell widths from cell-centre locations."""
    faces = np.empty(eta.size + 1)
    faces[0] = 0.0
    for index, center in enumerate(eta):
        faces[index + 1] = 2.0 * center - faces[index]
    if not np.isclose(faces[-1], 1.0, atol=1.0e-12):
        raise ValueError(f"reconstructed centreline face is {faces[-1]}, expected 1")
    return np.diff(faces)


def comparison_metrics(
    current: dict[str, np.ndarray | float], dns: dict[str, np.ndarray | float]
) -> dict[str, object]:
    eta = np.asarray(current["eta"])
    y_plus = np.asarray(current["y_plus"])
    overlap = y_plus <= float(np.asarray(dns["y_plus"])[-1]) + 1.0e-12
    result: dict[str, object] = {
        "friction_velocity": float(current["friction_velocity"]),
        "re_tau": float(current["re_tau"]),
        "re_tau_relative_error": float(current["re_tau"] / REFERENCE_RE_TAU - 1.0),
        "bulk_velocity": float(current["bulk_velocity"]),
        "mass_weighted_bulk_velocity": float(current["mass_weighted_bulk_velocity"]),
        "bulk_velocity_plus": float(current["bulk_velocity_plus"]),
        "bulk_velocity_plus_relative_error": float(
            current["bulk_velocity_plus"] / dns["bulk_velocity_plus"] - 1.0
        ),
        "skin_friction_coefficient": float(current["skin_friction_coefficient"]),
        "skin_friction_coefficient_relative_error": float(
            current["skin_friction_coefficient"] / dns["skin_friction_coefficient"] - 1.0
        ),
        "inner_scaled_profile_errors_on_common_y_plus": {},
        "outer_scaled_mean_velocity_error": {},
    }
    for name in (
        "mean_u_plus",
        "rms_u_plus",
        "rms_v_plus",
        "rms_w_plus",
        "minus_uv_plus",
    ):
        reference = np.interp(
            y_plus[overlap], np.asarray(dns["y_plus"]), np.asarray(dns[name])
        )
        values = np.asarray(current[name])[overlap]
        result["inner_scaled_profile_errors_on_common_y_plus"][name] = {
            "relative_l2": relative_l2(values, reference),
            "mean_absolute": float(np.mean(np.abs(values - reference))),
            "max_absolute": float(np.max(np.abs(values - reference))),
        }
    dns_outer = np.interp(
        eta, np.asarray(dns["eta"]), np.asarray(dns["mean_u_plus"])
    ) / float(dns["bulk_velocity_plus"])
    current_outer = np.asarray(current["mean_u_over_bulk"])
    result["outer_scaled_mean_velocity_error"] = {
        "relative_l2": relative_l2(current_outer, dns_outer),
        "mean_absolute": float(np.mean(np.abs(current_outer - dns_outer))),
        "max_absolute": float(np.max(np.abs(current_outer - dns_outer))),
    }
    result["profile_peaks"] = {}
    for name in ("rms_u_plus", "rms_v_plus", "rms_w_plus", "minus_uv_plus"):
        dns_index = int(np.argmax(np.asarray(dns[name])))
        current_index = int(np.argmax(np.asarray(current[name])))
        dns_peak = float(np.asarray(dns[name])[dns_index])
        current_peak = float(np.asarray(current[name])[current_index])
        result["profile_peaks"][name] = {
            "current": current_peak,
            "dns": dns_peak,
            "relative_error": current_peak / dns_peak - 1.0,
            "current_y_plus": float(y_plus[current_index]),
            "dns_y_plus": float(np.asarray(dns["y_plus"])[dns_index]),
        }
    return result


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    dns = load_dns()
    current_cases = {}
    metrics = {
        "reference": {
            "name": "ERCOFTAC Case 032 simul1.dat / Kim-Moin-Moser Re_tau=180",
            "actual_re_tau": REFERENCE_RE_TAU,
            "normalization": "u_tau and channel half-height h",
            "bulk_velocity_plus": dns["bulk_velocity_plus"],
            "skin_friction_coefficient": dns["skin_friction_coefficient"],
        },
        "current_profile_snapshot_times": list(PROFILE_TIMES),
        "current_statistics_scaling": (
            "equal-weight mean of u_tau and Re_tau interpolated at field times; "
            "bulk velocity integrated from the arithmetic mean profile"
        ),
        "cases": {},
    }
    for label, (stem, directories) in CASES.items():
        pooled = pooled_profiles(find_profiles(stem))
        folded = fold_channel(pooled)
        statistics = concatenate_tables(
            [read_table(next(directory.glob("*.statistics.r64.txt"))) for directory in directories]
        )
        times = np.asarray(PROFILE_TIMES, dtype=float)
        friction_velocity = float(
            np.mean(np.interp(times, statistics["time"], statistics["channel_friction_velocity"]))
        )
        re_tau = float(
            np.mean(np.interp(times, statistics["time"], statistics["channel_re_tau"]))
        )
        mass_weighted_bulk_velocity = float(
            np.mean(
                np.interp(
                    times,
                    statistics["time"],
                    statistics["total_momentum_x"] / statistics["total_mass"],
                )
            )
        )
        bulk_velocity = float(
            np.sum(folded["mean_u"] * cell_widths_from_centers(folded["eta"]))
        )
        current = {
            "eta": folded["eta"],
            "y_plus": folded["eta"] * re_tau,
            "mean_u_plus": folded["mean_u"] / friction_velocity,
            "mean_u_over_bulk": folded["mean_u"] / bulk_velocity,
            "rms_u_plus": folded["rms_u"] / friction_velocity,
            "rms_v_plus": folded["rms_v"] / friction_velocity,
            "rms_w_plus": folded["rms_w"] / friction_velocity,
            "minus_uv_plus": folded["minus_uv"] / friction_velocity**2,
            "friction_velocity": friction_velocity,
            "re_tau": re_tau,
            "bulk_velocity": bulk_velocity,
            "mass_weighted_bulk_velocity": mass_weighted_bulk_velocity,
            "bulk_velocity_plus": bulk_velocity / friction_velocity,
            "skin_friction_coefficient": 2.0 * (friction_velocity / bulk_velocity) ** 2,
        }
        current_cases[label] = current
        metrics["cases"][label] = comparison_metrics(current, dns)

    (OUTPUT / "dns_comparison_metrics.json").write_text(
        json.dumps(metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    colors = {"diss=0.001": "#d55e00", "diss=0.01": "#0072b2"}
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    axes[0].semilogx(
        dns["y_plus"], dns["mean_u_plus"], "ko", markersize=3.2, label="DNS, Re_tau=180"
    )
    axes[1].plot(
        np.asarray(dns["mean_u_plus"]) / float(dns["bulk_velocity_plus"]),
        dns["eta"],
        "ko",
        markersize=3.2,
        label="DNS",
    )
    for label, current in current_cases.items():
        axes[0].semilogx(
            current["y_plus"], current["mean_u_plus"], color=colors[label], linewidth=1.8, label=label
        )
        axes[1].plot(
            current["mean_u_over_bulk"], current["eta"], color=colors[label], linewidth=1.8, label=label
        )
    axes[0].set_xlabel(r"$y^+$")
    axes[0].set_ylabel(r"$U^+$")
    axes[0].set_xlim(0.8, 220.0)
    axes[1].set_xlabel(r"$U/U_b$")
    axes[1].set_ylabel(r"wall distance $y/h$")
    for axis in axes:
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(frameon=False)
    fig.suptitle("Case05 late-window mean velocity vs ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_mean_velocity_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.8), sharex=True)
    fields = (
        ("rms_u_plus", r"$u'_{rms}/u_\tau$"),
        ("rms_v_plus", r"$v'_{rms}/u_\tau$"),
        ("rms_w_plus", r"$w'_{rms}/u_\tau$"),
        ("minus_uv_plus", r"$-\langle u'v'\rangle/u_\tau^2$"),
    )
    for axis, (name, ylabel) in zip(axes.flat, fields):
        axis.semilogx(
            dns["y_plus"], dns[name], "ko", markersize=3.0, label="DNS, Re_tau=180"
        )
        for label, current in current_cases.items():
            axis.semilogx(
                current["y_plus"], current[name], color=colors[label], linewidth=1.7, label=label
            )
        axis.set_ylabel(ylabel)
        axis.set_xlim(0.8, 220.0)
        axis.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].set_xlabel(r"$y^+$")
    axes[1, 1].set_xlabel(r"$y^+$")
    fig.suptitle("Case05 late-window turbulence statistics vs ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_turbulence_statistics_comparison.png", dpi=180)
    plt.close(fig)

    profile_names = (
        ("mean_u_plus", r"$U^+$"),
        ("rms_u_plus", r"$u'_{rms}/u_\tau$"),
        ("rms_v_plus", r"$v'_{rms}/u_\tau$"),
        ("rms_w_plus", r"$w'_{rms}/u_\tau$"),
        ("minus_uv_plus", r"$-\langle u'v'\rangle/u_\tau^2$"),
    )
    scalar_names = (
        ("re_tau_relative_error", r"$Re_\tau$"),
        ("bulk_velocity_plus_relative_error", r"$U_b^+$"),
        ("skin_friction_coefficient_relative_error", r"$C_f$"),
    )
    x_profile = np.arange(len(profile_names))
    x_scalar = np.arange(len(scalar_names))
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    for offset, label in zip((-0.5, 0.5), CASES):
        case_metrics = metrics["cases"][label]
        profile_errors = [
            100.0
            * case_metrics["inner_scaled_profile_errors_on_common_y_plus"][name][
                "relative_l2"
            ]
            for name, _ in profile_names
        ]
        scalar_errors = [
            100.0 * case_metrics[name] for name, _ in scalar_names
        ]
        axes[0].bar(
            x_profile + offset * width,
            profile_errors,
            width,
            color=colors[label],
            label=label,
        )
        axes[1].bar(
            x_scalar + offset * width,
            scalar_errors,
            width,
            color=colors[label],
            label=label,
        )
    axes[0].set_xticks(x_profile, [label for _, label in profile_names])
    axes[0].set_ylabel("relative L2 error (%)")
    axes[0].set_title(r"Profiles on common $y^+$ range")
    axes[1].set_xticks(x_scalar, [label for _, label in scalar_names])
    axes[1].set_ylabel("signed relative error (%)")
    axes[1].set_title("Integral / friction quantities")
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    for axis in axes:
        axis.grid(True, axis="y", alpha=0.25)
        axis.legend(frameon=False)
    fig.suptitle("Error summary relative to ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_error_summary.png", dpi=180)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

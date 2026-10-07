#!/usr/bin/env python3
"""Analyze the completed extended Case05 branches and compare them with DNS."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_mdcd_dissipation_branches import pooled_profiles, profile_summary
from analyze_segment03_dns import (
    find_profile,
    load_concatenated,
    make_current_profile,
    scalar_summary,
)
from compare_branches_to_dns import comparison_metrics, load_dns


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = RESULTS / "branch-mdcd-dissipation-comparison" / "extended-final-analysis"
BODY_ACCELERATION = 0.0041720265499730312

DIRECTORIES = {
    "diss=0.001": [
        RESULTS / "branch-scmm6_mdcd0p001_roe-t223",
        RESULTS / "branch-scmm6_mdcd0p001_roe-t223_segment02",
        RESULTS / "branch-scmm6_mdcd0p001_roe-t223_segment03",
        RESULTS / "branch-scmm6_mdcd0p001_roe-t223_segment04",
    ],
    "diss=0.01": [
        RESULTS / "branch-scmm6_mdcd0p01_roe-t223",
        RESULTS / "branch-scmm6_mdcd0p01_roe-t223_segment02",
        RESULTS / "branch-scmm6_mdcd0p01_roe-t223_segment03",
    ],
}
STEMS = {
    "diss=0.001": "branch_scmm6_mdcd0p001_roe_t223",
    "diss=0.01": "branch_scmm6_mdcd0p01_roe_t223",
}
COMMON_PROFILE_TIMES = (300, 305, 310, 315, 320, 325)
LATEST_PROFILE_TIMES = {
    "diss=0.001": (410, 415, 420, 425, 430, 435),
    "diss=0.01": COMMON_PROFILE_TIMES,
}
PREVIOUS_PROFILE_TIMES = {
    "diss=0.001": (365, 370, 375, 380, 385, 390),
    "diss=0.01": (250, 255, 260, 265, 270, 275),
}
COLORS = {"diss=0.001": "#d55e00", "diss=0.01": "#0072b2"}


def continuity(previous: dict[str, np.ndarray], current: dict[str, np.ndarray]) -> dict[str, float]:
    return {
        name: float(current[name][0] - previous[name][-1])
        for name in (
            "step",
            "time",
            "total_mass",
            "total_momentum_x",
            "total_energy",
            "channel_re_tau",
        )
    }


def health_summary(history: dict[str, np.ndarray]) -> dict[str, float | int]:
    positive_dt = history["dt"][history["dt"] > 0.0]
    result: dict[str, float | int] = {
        name: int(np.nanmax(history[name]))
        for name in (
            "reconstruction_fallbacks",
            "riemann_fallbacks",
            "troubled_cells",
            "local_recomputations",
            "step_retries",
        )
    }
    result.update(
        {
            "last_time": float(history["time"][-1]),
            "last_step": int(history["step"][-1]),
            "wall_time": float(history["wall_time"][-1]),
            "dt_median_positive": float(np.median(positive_dt)),
            "dt_min_positive": float(np.min(positive_dt)),
            "dt_max": float(np.max(history["dt"])),
        }
    )
    return result


def aligned_blocks(
    table: dict[str, np.ndarray], start: float, end: float, width: float = 10.0
) -> list[dict]:
    blocks = []
    block_start = float(np.ceil(start / width) * width)
    while block_start + width <= end + 1.0e-10:
        blocks.append(scalar_summary(table, block_start, block_start + width))
        block_start += width
    if end - block_start >= 2.0:
        blocks.append(scalar_summary(table, block_start, end))
    return blocks


def momentum_balance(window: dict) -> dict[str, float]:
    wall_shear = float(window["wall_shear_mean"]["mean"])
    measured = float(window["bulk_velocity"]["slope_per_time"])
    predicted = BODY_ACCELERATION - wall_shear
    return {
        "body_acceleration": BODY_ACCELERATION,
        "mean_wall_shear": wall_shear,
        "wall_shear_relative_to_body_force": wall_shear / BODY_ACCELERATION - 1.0,
        "predicted_bulk_velocity_slope": predicted,
        "measured_bulk_velocity_slope": measured,
        "prediction_relative_difference": measured / predicted - 1.0
        if abs(predicted) > 1.0e-15
        else float("nan"),
    }


def pooled_quality(label: str, times: tuple[int, ...]) -> dict[str, float]:
    return profile_summary(
        pooled_profiles([find_profile(STEMS[label], time) for time in times])
    )


def profile_window_change(previous: dict, current: dict) -> dict[str, float]:
    result = {}
    for name in (
        "mean_u_over_bulk",
        "mean_u_plus",
        "rms_u_plus",
        "rms_v_plus",
        "rms_w_plus",
        "minus_uv_plus",
    ):
        old = np.asarray(previous[name])
        new = np.asarray(current[name])
        result[name] = float(np.linalg.norm(new - old) / np.linalg.norm(old))
    return result


def plot_histories(statistics: dict[str, dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(11.0, 9.3), sharex=True)
    for label, table in statistics.items():
        bulk = table["total_momentum_x"] / table["total_mass"]
        asymmetry = (
            np.abs(table["channel_wall_shear_lower"] - table["channel_wall_shear_upper"])
            / table["channel_wall_shear_mean"]
        )
        axes[0].plot(table["time"], bulk, color=COLORS[label], linewidth=1.2, label=label)
        axes[1].plot(
            table["time"], table["channel_re_tau"], color=COLORS[label], linewidth=1.2, label=label
        )
        axes[2].plot(table["time"], 100.0 * asymmetry, color=COLORS[label], linewidth=1.0, label=label)
    axes[1].axhline(180.0, color="black", linestyle=":", linewidth=1.1, label="DNS 180")
    axes[0].set_ylabel(r"$P_x/M$")
    axes[1].set_ylabel(r"$Re_\tau$")
    axes[2].set_ylabel("wall-shear asymmetry (%)")
    axes[2].set_xlabel("time")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False, ncol=3)
    fig.suptitle("Case05 complete branch histories")
    fig.tight_layout()
    fig.savefig(OUTPUT / "extended_branch_histories.png", dpi=180)
    plt.close(fig)


def plot_dns_profiles(
    dns: dict[str, np.ndarray | float],
    cases: dict[str, dict],
    prefix: str,
    title: str,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    axes[0].semilogx(
        dns["y_plus"], dns["mean_u_plus"], "ko", markersize=3.2, label=r"DNS, $Re_\tau=180$"
    )
    axes[1].plot(
        np.asarray(dns["mean_u_plus"]) / float(dns["bulk_velocity_plus"]),
        dns["eta"],
        "ko",
        markersize=3.2,
        label="DNS",
    )
    for label, current in cases.items():
        base_label = label.split(",")[0]
        axes[0].semilogx(
            current["y_plus"], current["mean_u_plus"], color=COLORS[base_label], linewidth=1.8, label=label
        )
        axes[1].plot(
            current["mean_u_over_bulk"], current["eta"], color=COLORS[base_label], linewidth=1.8, label=label
        )
    axes[0].set_xlabel(r"$y^+$")
    axes[0].set_ylabel(r"$U^+$")
    axes[0].set_xlim(0.8, 220.0)
    axes[1].set_xlabel(r"$U/U_b$")
    axes[1].set_ylabel(r"wall distance $y/h$")
    for axis in axes:
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(OUTPUT / f"{prefix}_mean_velocity.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.8), sharex=True)
    fields = (
        ("rms_u_plus", r"$u'_{rms}/u_\tau$"),
        ("rms_v_plus", r"$v'_{rms}/u_\tau$"),
        ("rms_w_plus", r"$w'_{rms}/u_\tau$"),
        ("minus_uv_plus", r"$-\langle u'v'\rangle/u_\tau^2$"),
    )
    for axis, (name, ylabel) in zip(axes.flat, fields):
        axis.semilogx(dns["y_plus"], dns[name], "ko", markersize=3.0, label=r"DNS, $Re_\tau=180$")
        for label, current in cases.items():
            base_label = label.split(",")[0]
            axis.semilogx(
                current["y_plus"], current[name], color=COLORS[base_label], linewidth=1.7, label=label
            )
        axis.set_ylabel(ylabel)
        axis.set_xlim(0.8, 220.0)
        axis.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].set_xlabel(r"$y^+$")
    axes[1, 1].set_xlabel(r"$y^+$")
    fig.suptitle(title.replace("mean velocity", "turbulence statistics"))
    fig.tight_layout()
    fig.savefig(OUTPUT / f"{prefix}_turbulence_statistics.png", dpi=180)
    plt.close(fig)


def plot_error_summary(metrics: dict[str, dict], filename: str, title: str) -> None:
    profile_fields = (
        ("mean_u_plus", r"$U^+$"),
        ("rms_u_plus", r"$u'_{rms}/u_\tau$"),
        ("rms_v_plus", r"$v'_{rms}/u_\tau$"),
        ("rms_w_plus", r"$w'_{rms}/u_\tau$"),
        ("minus_uv_plus", r"$-\langle u'v'\rangle/u_\tau^2$"),
    )
    scalar_fields = (
        ("re_tau_relative_error", r"$Re_\tau$"),
        ("bulk_velocity_plus_relative_error", r"$U_b^+$"),
        ("skin_friction_coefficient_relative_error", r"$C_f$"),
    )
    x_profile = np.arange(len(profile_fields))
    x_scalar = np.arange(len(scalar_fields))
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    for index, (label, case) in enumerate(metrics.items()):
        base_label = label.split(",")[0]
        profile_values = [
            100.0 * case["inner_scaled_profile_errors_on_common_y_plus"][name]["relative_l2"]
            for name, _ in profile_fields
        ]
        scalar_values = [100.0 * case[name] for name, _ in scalar_fields]
        offset = (index - 0.5) * width
        axes[0].bar(x_profile + offset, profile_values, width, color=COLORS[base_label], label=label)
        axes[1].bar(x_scalar + offset, scalar_values, width, color=COLORS[base_label], label=label)
    axes[0].set_xticks(x_profile, [label for _, label in profile_fields])
    axes[0].set_ylabel("relative L2 error (%)")
    axes[0].set_title(r"Profiles on common $y^+$ range")
    axes[1].set_xticks(x_scalar, [label for _, label in scalar_fields])
    axes[1].set_ylabel("signed relative error (%)")
    axes[1].set_title("Integral / friction quantities")
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    for axis in axes:
        axis.grid(True, axis="y", alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(OUTPUT / filename, dpi=180)
    plt.close(fig)


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    dns = load_dns()
    statistics: dict[str, dict[str, np.ndarray]] = {}
    metadata = {}
    histories = {}
    history_metadata = {}
    for label, directories in DIRECTORIES.items():
        statistics[label], metadata[label] = load_concatenated(directories, "statistics")
        histories[label], history_metadata[label] = load_concatenated(directories, "history")

    common_profiles = {
        label: make_current_profile(STEMS[label], COMMON_PROFILE_TIMES, statistics[label])
        for label in DIRECTORIES
    }
    latest_profiles = {
        label: make_current_profile(STEMS[label], LATEST_PROFILE_TIMES[label], statistics[label])
        for label in DIRECTORIES
    }
    previous_profiles = {
        label: make_current_profile(STEMS[label], PREVIOUS_PROFILE_TIMES[label], statistics[label])
        for label in DIRECTORIES
    }
    common_metrics = {
        f"{label}, t={COMMON_PROFILE_TIMES[0]}-{COMMON_PROFILE_TIMES[-1]}": comparison_metrics(profile, dns)
        for label, profile in common_profiles.items()
    }
    latest_metrics = {
        f"{label}, t={LATEST_PROFILE_TIMES[label][0]}-{LATEST_PROFILE_TIMES[label][-1]}": comparison_metrics(
            latest_profiles[label], dns
        )
        for label in DIRECTORIES
    }
    previous_metrics = {
        f"{label}, t={PREVIOUS_PROFILE_TIMES[label][0]}-{PREVIOUS_PROFILE_TIMES[label][-1]}": comparison_metrics(
            previous_profiles[label], dns
        )
        for label in DIRECTORIES
    }

    common_start = 280.0
    common_end = min(float(table["time"][-1]) for table in statistics.values())
    branch_report = {}
    for label in DIRECTORIES:
        segment = metadata[label]["segments"][-1]
        previous_segment = metadata[label]["segments"][-2]
        history_segment = history_metadata[label]["segments"][-1]
        last_time = float(segment["time"][-1])
        last_20 = scalar_summary(statistics[label], last_time - 20.0, last_time)
        branch_report[label] = {
            "new_segment": {
                "start_time": float(segment["time"][0]),
                "end_time": last_time,
                "start_step": int(segment["step"][0]),
                "end_step": int(segment["step"][-1]),
                "duration": last_time - float(segment["time"][0]),
                "stop_reason": "wall_time_checkpoint",
                "restart_continuity_delta": continuity(previous_segment, segment),
                "mass_relative_drift": float(
                    (segment["total_mass"][-1] - segment["total_mass"][0])
                    / segment["total_mass"][0]
                ),
                "health": health_summary(history_segment),
                "last_20_time_units": last_20,
                "last_20_momentum_balance": momentum_balance(last_20),
                "ten_time_unit_blocks": aligned_blocks(
                    statistics[label], float(segment["time"][0]), last_time
                ),
            },
            "common_integral_window": scalar_summary(
                statistics[label], common_start, common_end
            ),
            "common_profile_times": list(COMMON_PROFILE_TIMES),
            "latest_profile_times": list(LATEST_PROFILE_TIMES[label]),
            "previous_profile_times": list(PREVIOUS_PROFILE_TIMES[label]),
            "common_profile_quality": pooled_quality(label, COMMON_PROFILE_TIMES),
            "latest_profile_quality": pooled_quality(label, LATEST_PROFILE_TIMES[label]),
            "latest_vs_previous_profile_change": profile_window_change(
                previous_profiles[label], latest_profiles[label]
            ),
        }

    report = {
        "reference": {
            "name": "ERCOFTAC Case 032 simul1.dat / Kim-Moin-Moser Re_tau=180",
            "re_tau": 180.0,
            "bulk_velocity_plus": float(dns["bulk_velocity_plus"]),
            "skin_friction_coefficient": float(dns["skin_friction_coefficient"]),
        },
        "common_integral_window": [common_start, common_end],
        "common_profile_times": list(COMMON_PROFILE_TIMES),
        "branches": branch_report,
        "dns_comparison": {
            "strict_common_window_300_325": common_metrics,
            "latest_available_windows": latest_metrics,
            "previous_windows": previous_metrics,
        },
    }
    (OUTPUT / "extended_dns_metrics.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    plot_histories(statistics)
    common_plot_cases = {
        f"{label}, t=300-325": common_profiles[label] for label in DIRECTORIES
    }
    latest_plot_cases = {
        f"{label}, t={LATEST_PROFILE_TIMES[label][0]}-{LATEST_PROFILE_TIMES[label][-1]}": latest_profiles[label]
        for label in DIRECTORIES
    }
    plot_dns_profiles(
        dns,
        common_plot_cases,
        "dns_common_300_325",
        "Common-window mean velocity vs ERCOFTAC Case 032 DNS",
    )
    plot_dns_profiles(
        dns,
        latest_plot_cases,
        "dns_latest_available",
        "Latest-available mean velocity vs ERCOFTAC Case 032 DNS",
    )
    plot_error_summary(
        common_metrics,
        "dns_common_300_325_error_summary.png",
        "Common-window errors relative to ERCOFTAC Case 032 DNS",
    )
    plot_error_summary(
        latest_metrics,
        "dns_latest_available_error_summary.png",
        "Latest-available errors relative to ERCOFTAC Case 032 DNS",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Analyze the completed diss=0.001 segment03 and compare both branches with DNS."""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_mdcd_dissipation_branches import (
    concatenate_tables,
    pooled_profiles,
    profile_summary,
)
from compare_branches_to_dns import (
    REFERENCE_RE_TAU,
    cell_widths_from_centers,
    comparison_metrics,
    fold_channel,
    load_dns,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
PROFILE_DIRECTORY = RESULTS / "branch-mdcd-dissipation-comparison" / "profiles"
OUTPUT = RESULTS / "branch-mdcd-dissipation-comparison" / "segment03-final-analysis"

BRANCH_001_DIRECTORIES = [
    RESULTS / "branch-scmm6_mdcd0p001_roe-t223",
    RESULTS / "branch-scmm6_mdcd0p001_roe-t223_segment02",
    RESULTS / "branch-scmm6_mdcd0p001_roe-t223_segment03",
]
BRANCH_01_DIRECTORIES = [
    RESULTS / "branch-scmm6_mdcd0p01_roe-t223",
    RESULTS / "branch-scmm6_mdcd0p01_roe-t223_segment02",
]
COMMON_TIMES = (250, 255, 260, 265, 270, 275)


def find_table(directory: Path, kind: str) -> Path:
    matches = sorted(directory.glob(f"*.{kind}.r64.txt"))
    if not matches:
        matches = sorted(directory.glob(f"*.{kind}.r64.txt.tmp"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {kind} table in {directory}, found {len(matches)}")
    return matches[0]


def read_live_table(path: Path) -> tuple[dict[str, np.ndarray], int]:
    """Read only complete numeric rows, allowing a solver to append concurrently."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].lstrip("# ").split()
    numeric_count = len(header) - (1 if header[-1] == "stop_reason" else 0)
    names = header[:numeric_count]
    rows: list[list[float]] = []
    discarded = 0
    for line in lines[1:]:
        fields = line.split()
        if not fields or line.lstrip().startswith("#"):
            continue
        if len(fields) < numeric_count:
            discarded += 1
            continue
        try:
            rows.append([float(value) for value in fields[:numeric_count]])
        except ValueError:
            discarded += 1
    if not rows:
        raise ValueError(f"no complete data rows in {path}")
    values = np.asarray(rows)
    return ({name: values[:, index] for index, name in enumerate(names)}, discarded)


def load_concatenated(directories: list[Path], kind: str) -> tuple[dict[str, np.ndarray], dict]:
    tables = []
    sources = []
    for directory in directories:
        path = find_table(directory, kind)
        table, discarded = read_live_table(path)
        tables.append(table)
        sources.append(
            {
                "path": str(path.relative_to(ROOT)),
                "first_time": float(table["time"][0]),
                "last_complete_time": float(table["time"][-1]),
                "complete_rows": int(table["time"].size),
                "discarded_incomplete_rows": discarded,
            }
        )
    return concatenate_tables(tables), {"sources": sources, "segments": tables}


def profile_time(path: Path) -> float | None:
    match = re.search(r"time([0-9]+)p([0-9]+)eP([0-9]+)", path.name)
    if match is None:
        return None
    mantissa = float(f"{match.group(1)}.{match.group(2)}")
    return mantissa * 10.0 ** int(match.group(3))


def find_profile(stem: str, time: int) -> Path:
    candidates = []
    for path in PROFILE_DIRECTORY.glob(f"*{stem}*.profile.txt"):
        parsed = profile_time(path)
        if parsed is not None and abs(parsed - time) < 1.0e-7:
            candidates.append(path)
    if len(candidates) != 1:
        raise RuntimeError(f"expected one profile for {stem} at t={time}, found {len(candidates)}")
    return candidates[0]


def latest_segment03_times() -> tuple[int, ...]:
    times = []
    for path in PROFILE_DIRECTORY.glob("*mdcd0p001*segment03*.profile.txt"):
        parsed = profile_time(path)
        if parsed is not None and abs(parsed / 5.0 - round(parsed / 5.0)) < 1.0e-8:
            times.append(int(round(parsed)))
    times = sorted(set(times))
    if len(times) < 6:
        raise RuntimeError(f"need at least six complete segment03 profiles, found {times}")
    return tuple(times[-6:])


def make_current_profile(
    stem: str,
    times: tuple[int, ...],
    statistics: dict[str, np.ndarray],
) -> dict[str, np.ndarray | float]:
    pooled = pooled_profiles([find_profile(stem, time) for time in times])
    folded = fold_channel(pooled)
    sample_times = np.asarray(times, dtype=float)
    friction_velocity = float(
        np.mean(
            np.interp(
                sample_times,
                statistics["time"],
                statistics["channel_friction_velocity"],
            )
        )
    )
    re_tau = float(
        np.mean(
            np.interp(sample_times, statistics["time"], statistics["channel_re_tau"])
        )
    )
    bulk_velocity = float(
        np.sum(folded["mean_u"] * cell_widths_from_centers(folded["eta"]))
    )
    mass_weighted_bulk_velocity = float(
        np.mean(
            np.interp(
                sample_times,
                statistics["time"],
                statistics["total_momentum_x"] / statistics["total_mass"],
            )
        )
    )
    return {
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


def scalar_summary(table: dict[str, np.ndarray], start: float, end: float) -> dict:
    mask = (table["time"] >= start - 1.0e-10) & (table["time"] <= end + 1.0e-10)
    time = table["time"][mask]
    if not np.any(mask):
        raise ValueError(f"empty statistics window {start}--{end}")
    bulk = table["total_momentum_x"][mask] / table["total_mass"][mask]
    lower = table["channel_wall_shear_lower"][mask]
    upper = table["channel_wall_shear_upper"][mask]
    mean_shear = table["channel_wall_shear_mean"][mask]

    def summarize(values: np.ndarray) -> dict[str, float]:
        return {
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
            "slope_per_time": float(np.polyfit(time, values, 1)[0])
            if values.size > 1
            else float("nan"),
        }

    return {
        "start": float(time[0]),
        "end": float(time[-1]),
        "samples": int(time.size),
        "bulk_velocity": summarize(bulk),
        "re_tau": summarize(table["channel_re_tau"][mask]),
        "friction_velocity": summarize(table["channel_friction_velocity"][mask]),
        "wall_shear_mean": summarize(mean_shear),
        "wall_shear_relative_asymmetry": {
            "mean": float(np.mean(np.abs(lower - upper) / mean_shear)),
            "max": float(np.max(np.abs(lower - upper) / mean_shear)),
        },
    }


def plot_global_history(stats_001: dict[str, np.ndarray], stats_01: dict[str, np.ndarray]) -> None:
    colors = {"diss=0.001": "#d55e00", "diss=0.01": "#0072b2"}
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True)
    for label, table in (("diss=0.001", stats_001), ("diss=0.01", stats_01)):
        bulk = table["total_momentum_x"] / table["total_mass"]
        axes[0].plot(table["time"], bulk, color=colors[label], linewidth=1.2, label=label)
        axes[1].plot(
            table["time"], table["channel_re_tau"], color=colors[label], linewidth=1.2, label=label
        )
    for axis in axes:
        axis.axvline(279.3638236005483, color="0.35", linestyle="--", linewidth=1.0)
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False)
    axes[1].axhline(REFERENCE_RE_TAU, color="black", linestyle=":", linewidth=1.1, label="DNS 180")
    axes[1].legend(frameon=False)
    axes[0].set_ylabel(r"$P_x/M$")
    axes[1].set_ylabel(r"$Re_\tau$")
    axes[1].set_xlabel("time")
    fig.suptitle("Case05 branch histories including completed segment03")
    fig.tight_layout()
    fig.savefig(OUTPUT / "segment03_global_history.png", dpi=180)
    plt.close(fig)


def plot_dns_profiles(
    dns: dict[str, np.ndarray | float], cases: dict[str, dict[str, np.ndarray | float]]
) -> None:
    colors = dict(zip(cases, ("#d55e00", "#0072b2")))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    axes[0].semilogx(dns["y_plus"], dns["mean_u_plus"], "ko", markersize=3.2, label=r"DNS, $Re_\tau=180$")
    axes[1].plot(
        np.asarray(dns["mean_u_plus"]) / float(dns["bulk_velocity_plus"]),
        dns["eta"],
        "ko",
        markersize=3.2,
        label="DNS",
    )
    for label, current in cases.items():
        axes[0].semilogx(current["y_plus"], current["mean_u_plus"], color=colors[label], linewidth=1.8, label=label)
        axes[1].plot(current["mean_u_over_bulk"], current["eta"], color=colors[label], linewidth=1.8, label=label)
    axes[0].set_xlabel(r"$y^+$")
    axes[0].set_ylabel(r"$U^+$")
    axes[0].set_xlim(0.8, 220.0)
    axes[1].set_xlabel(r"$U/U_b$")
    axes[1].set_ylabel(r"wall distance $y/h$")
    for axis in axes:
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle("Latest available branch statistics vs ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_latest_mean_velocity.png", dpi=180)
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
            axis.semilogx(current["y_plus"], current[name], color=colors[label], linewidth=1.7, label=label)
        axis.set_ylabel(ylabel)
        axis.set_xlim(0.8, 220.0)
        axis.grid(True, which="both", alpha=0.25)
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].set_xlabel(r"$y^+$")
    axes[1, 1].set_xlabel(r"$y^+$")
    fig.suptitle("Latest available turbulence statistics vs ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_latest_turbulence_statistics.png", dpi=180)
    plt.close(fig)


def plot_evolution(
    dns: dict[str, np.ndarray | float],
    old: dict[str, np.ndarray | float],
    new: dict[str, np.ndarray | float],
    latest_label: str,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4))
    fields = (
        ("mean_u_plus", r"$U^+$"),
        ("rms_u_plus", r"$u'_{rms}/u_\tau$"),
        ("minus_uv_plus", r"$-\langle u'v'\rangle/u_\tau^2$"),
    )
    for axis, (name, ylabel) in zip(axes, fields):
        axis.semilogx(dns["y_plus"], dns[name], "ko", markersize=3.0, label="DNS")
        axis.semilogx(old["y_plus"], old[name], color="#e69f00", linewidth=1.7, label="0.001, t=250-275")
        axis.semilogx(new["y_plus"], new[name], color="#d55e00", linewidth=1.8, label=latest_label)
        axis.set_xlabel(r"$y^+$")
        axis.set_ylabel(ylabel)
        axis.set_xlim(0.8, 220.0)
        axis.grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Evolution of the diss=0.001 branch relative to DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "diss0p001_profile_evolution_dns.png", dpi=180)
    plt.close(fig)


def plot_error_summary(metrics: dict[str, dict]) -> None:
    labels = tuple(metrics)
    colors = ("#d55e00", "#0072b2")
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
    for index, (label, color) in enumerate(zip(labels, colors)):
        case = metrics[label]
        profile_values = [
            100.0
            * case["inner_scaled_profile_errors_on_common_y_plus"][name]["relative_l2"]
            for name, _ in profile_fields
        ]
        scalar_values = [100.0 * case[name] for name, _ in scalar_fields]
        offset = (index - 0.5) * width
        axes[0].bar(x_profile + offset, profile_values, width, color=color, label=label)
        axes[1].bar(x_scalar + offset, scalar_values, width, color=color, label=label)
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
    fig.suptitle("Latest available errors relative to ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_latest_error_summary.png", dpi=180)
    plt.close(fig)


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    dns = load_dns()
    statistics_001, metadata_001 = load_concatenated(BRANCH_001_DIRECTORIES, "statistics")
    statistics_01, metadata_01 = load_concatenated(BRANCH_01_DIRECTORIES, "statistics")
    history_001, history_metadata_001 = load_concatenated(BRANCH_001_DIRECTORIES, "history")

    latest_times = latest_segment03_times()
    latest_001_label = f"diss=0.001, t={latest_times[0]}-{latest_times[-1]}"
    latest_01_label = f"diss=0.01, t={COMMON_TIMES[0]}-{COMMON_TIMES[-1]}"
    common_001 = make_current_profile(
        "branch_scmm6_mdcd0p001_roe_t223", COMMON_TIMES, statistics_001
    )
    common_01 = make_current_profile(
        "branch_scmm6_mdcd0p01_roe_t223", COMMON_TIMES, statistics_01
    )
    latest_001 = make_current_profile(
        "branch_scmm6_mdcd0p001_roe_t223_segment03", latest_times, statistics_001
    )
    latest_profile_quality = profile_summary(
        pooled_profiles(
            [
                find_profile("branch_scmm6_mdcd0p001_roe_t223_segment03", time)
                for time in latest_times
            ]
        )
    )

    segment03 = metadata_001["segments"][-1]
    segment02 = metadata_001["segments"][-2]
    join_fields = (
        "step",
        "time",
        "total_mass",
        "total_momentum_x",
        "channel_re_tau",
    )
    continuity = {
        name: float(segment03[name][0] - segment02[name][-1]) for name in join_fields
    }
    block_windows = []
    block_start = 280.0
    last_complete = float(segment03["time"][-1])
    while block_start < last_complete - 1.0e-10:
        block_end = min(block_start + 10.0, last_complete)
        block_windows.append(scalar_summary(statistics_001, block_start, block_end))
        block_start += 10.0

    segment03_history = history_metadata_001["segments"][-1]
    health = {
        name: int(np.nanmax(segment03_history[name]))
        for name in (
            "reconstruction_fallbacks",
            "riemann_fallbacks",
            "troubled_cells",
            "local_recomputations",
            "step_retries",
        )
    }
    positive_dt = segment03_history["dt"][segment03_history["dt"] > 0.0]
    health.update(
        {
            "history_last_complete_time": float(history_001["time"][-1]),
            "history_last_complete_step": int(history_001["step"][-1]),
            "dt_median_positive": float(np.median(positive_dt)),
            "dt_min_positive": float(np.min(positive_dt)),
            "dt_max": float(np.max(segment03_history["dt"])),
        }
    )

    fair_metrics = {
        "diss=0.001": comparison_metrics(common_001, dns),
        "diss=0.01": comparison_metrics(common_01, dns),
    }
    latest_metrics = {
        latest_001_label: comparison_metrics(latest_001, dns),
        latest_01_label: comparison_metrics(common_01, dns),
    }
    report = {
        "analysis_is_live_snapshot": False,
        "reference": {
            "name": "ERCOFTAC Case 032 simul1.dat / Kim-Moin-Moser Re_tau=180",
            "re_tau": REFERENCE_RE_TAU,
            "bulk_velocity_plus": float(dns["bulk_velocity_plus"]),
            "skin_friction_coefficient": float(dns["skin_friction_coefficient"]),
        },
        "input_tables": {
            "diss_0p001_statistics": metadata_001["sources"],
            "diss_0p01_statistics": metadata_01["sources"],
            "diss_0p001_history": history_metadata_001["sources"],
        },
        "segment03": {
            "statistics_first_time": float(segment03["time"][0]),
            "statistics_last_complete_time": last_complete,
            "statistics_last_complete_step": int(segment03["step"][-1]),
            "stop_reason": "wall_time_checkpoint",
            "restart_continuity_delta": continuity,
            "mass_relative_drift": float(
                (segment03["total_mass"][-1] - segment03["total_mass"][0])
                / segment03["total_mass"][0]
            ),
            "health": health,
            "last_20_time_units": scalar_summary(
                statistics_001, max(segment03["time"][0], last_complete - 20.0), last_complete
            ),
            "ten_time_unit_blocks": block_windows,
            "latest_complete_profile_times": list(latest_times),
            "latest_pooled_profile_quality": latest_profile_quality,
        },
        "dns_comparison": {
            "fair_common_window_250_275": fair_metrics,
            "latest_available_windows": latest_metrics,
        },
    }
    (OUTPUT / "segment03_dns_metrics.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    plot_global_history(statistics_001, statistics_01)
    plot_dns_profiles(
        dns,
        {
            latest_001_label: latest_001,
            latest_01_label: common_01,
        },
    )
    plot_evolution(dns, common_001, latest_001, latest_001_label)
    plot_error_summary(latest_metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare the Case05 SCMM6/MDCD-linear/Roe dissipation branches."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = RESULTS / "branch-mdcd-dissipation-comparison"
CASES = {
    "diss=0.001": [
        RESULTS / "branch-scmm6_mdcd0p001_roe-t223",
        RESULTS / "branch-scmm6_mdcd0p001_roe-t223_segment02",
    ],
    "diss=0.01": [
        RESULTS / "branch-scmm6_mdcd0p01_roe-t223",
        RESULTS / "branch-scmm6_mdcd0p01_roe-t223_segment02",
    ],
}
PROFILE_WINDOWS = {
    "225_245": (225, 230, 235, 240, 245),
    "250_275": (250, 255, 260, 265, 270, 275),
}


def read_table(path: Path) -> dict[str, np.ndarray]:
    header = path.read_text(encoding="utf-8").splitlines()[0].lstrip("# ").split()
    numeric_count = len(header) - (1 if header[-1] == "stop_reason" else 0)
    values = np.loadtxt(path, comments="#", usecols=range(numeric_count))
    header = header[:numeric_count]
    return {name: values[:, index] for index, name in enumerate(header)}


def fit_slope(time: np.ndarray, value: np.ndarray) -> float:
    if time.size < 2:
        return float("nan")
    return float(np.polyfit(time, value, 1)[0])


def window_summary(table: dict[str, np.ndarray], start: float, end: float) -> dict:
    time = table["time"]
    mask = (time >= start - 1.0e-10) & (time <= end + 1.0e-10)
    result: dict[str, object] = {
        "start": start,
        "end": end,
        "samples": int(mask.sum()),
    }
    for name in (
        "total_mass",
        "total_momentum_x",
        "total_energy",
        "yz_mean_u_plane0",
        "yz_mean_u_plane1",
        "yz_mass_flow_x_plane0",
        "yz_mass_flow_x_plane1",
        "channel_wall_shear_lower",
        "channel_wall_shear_upper",
        "channel_wall_shear_mean",
        "channel_friction_velocity",
        "channel_re_tau",
    ):
        selected = table[name][mask]
        result[name] = {
            "mean": float(np.mean(selected)),
            "std": float(np.std(selected, ddof=1)) if selected.size > 1 else 0.0,
            "min": float(np.min(selected)),
            "max": float(np.max(selected)),
            "slope_per_time": fit_slope(time[mask], selected),
        }
    bulk = table["total_momentum_x"][mask] / table["total_mass"][mask]
    result["mass_weighted_bulk_velocity"] = {
        "mean": float(np.mean(bulk)),
        "std": float(np.std(bulk, ddof=1)) if bulk.size > 1 else 0.0,
        "slope_per_time": fit_slope(time[mask], bulk),
    }
    lower = table["channel_wall_shear_lower"][mask]
    upper = table["channel_wall_shear_upper"][mask]
    mean_shear = table["channel_wall_shear_mean"][mask]
    asymmetry = np.abs(lower - upper) / mean_shear
    result["wall_shear_relative_asymmetry"] = {
        "mean": float(np.mean(asymmetry)),
        "max": float(np.max(asymmetry)),
    }
    plane_u_difference = np.abs(
        table["yz_mean_u_plane0"][mask] - table["yz_mean_u_plane1"][mask]
    )
    plane_flow_difference = np.abs(
        table["yz_mass_flow_x_plane0"][mask]
        - table["yz_mass_flow_x_plane1"][mask]
    )
    result["yz_plane_u_max_abs_difference"] = float(np.max(plane_u_difference))
    result["yz_plane_flow_max_abs_difference"] = float(
        np.max(plane_flow_difference)
    )
    return result


def concatenate_tables(tables: list[dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    names = tuple(tables[0])
    if any(tuple(table) != names for table in tables[1:]):
        raise ValueError("table columns do not match")
    pieces: dict[str, list[np.ndarray]] = {name: [] for name in names}
    last_time = -np.inf
    for table in tables:
        mask = table["time"] > last_time + 1.0e-10
        for name in names:
            pieces[name].append(table[name][mask])
        if np.any(mask):
            last_time = float(table["time"][mask][-1])
    return {name: np.concatenate(values) for name, values in pieces.items()}


def history_summary(tables: list[dict[str, np.ndarray]]) -> dict:
    table = concatenate_tables(tables)
    positive_dt = table["dt"][table["dt"] > 0.0]
    advanced_steps = table["step"][-1] - table["step"][0]
    elapsed_wall = sum(
        float(segment["wall_time"][-1] - segment["wall_time"][0])
        for segment in tables
    )
    return {
        "rows": int(table["step"].size),
        "segments": len(tables),
        "dt_min_positive": float(np.min(positive_dt)),
        "dt_median_positive": float(np.median(positive_dt)),
        "dt_mean_positive": float(np.mean(positive_dt)),
        "dt_max": float(np.max(table["dt"])),
        "wall_seconds": float(elapsed_wall),
        "steps_per_wall_second": float(advanced_steps / elapsed_wall),
        "reconstruction_fallbacks_max": int(np.max(table["reconstruction_fallbacks"])),
        "riemann_fallbacks_max": int(np.max(table["riemann_fallbacks"])),
        "troubled_cells_max": int(np.max(table["troubled_cells"])),
        "local_recomputations_max": int(np.max(table["local_recomputations"])),
        "step_retries_max": int(np.max(table["step_retries"])),
    }


def read_profile(path: Path) -> dict[str, np.ndarray]:
    header = path.read_text(encoding="utf-8").splitlines()[0].lstrip("# ").split()
    values = np.loadtxt(path, comments="#")
    return {name: values[:, index] for index, name in enumerate(header)}


def pooled_profiles(paths: list[Path]) -> dict[str, np.ndarray]:
    """Pool x-z moments over equally spaced field snapshots."""
    profiles = [read_profile(path) for path in paths]
    y = profiles[0]["y"]
    if any(not np.allclose(profile["y"], y) for profile in profiles[1:]):
        raise ValueError("profile y coordinates do not match")

    result = {"y": y}
    for name in ("mean_u", "mean_v", "mean_w", "mean_rho", "mean_temperature", "mean_mach"):
        result[name] = np.mean([profile[name] for profile in profiles], axis=0)
    for velocity, rms in (("mean_u", "rms_u"), ("mean_v", "rms_v"), ("mean_w", "rms_w")):
        second_moment = np.mean(
            [profile[rms] ** 2 + profile[velocity] ** 2 for profile in profiles],
            axis=0,
        )
        result[rms] = np.sqrt(np.maximum(0.0, second_moment - result[velocity] ** 2))
    mean_uv = np.mean(
        [profile["reynolds_uv"] + profile["mean_u"] * profile["mean_v"] for profile in profiles],
        axis=0,
    )
    result["reynolds_uv"] = mean_uv - result["mean_u"] * result["mean_v"]
    return result


def profile_summary(profile: dict[str, np.ndarray]) -> dict[str, float]:
    center = np.argsort(np.abs(profile["y"]))[:2]
    result = {
        "centerline_mean_u": float(np.mean(profile["mean_u"][center])),
        "peak_rms_u": float(np.max(profile["rms_u"])),
        "peak_rms_v": float(np.max(profile["rms_v"])),
        "peak_rms_w": float(np.max(profile["rms_w"])),
        "peak_abs_reynolds_uv": float(np.max(np.abs(profile["reynolds_uv"]))),
        "max_mean_mach": float(np.max(profile["mean_mach"])),
    }
    for name in ("mean_u", "rms_u", "rms_v", "rms_w"):
        result[f"{name}_mirror_l2_relative"] = float(
            np.linalg.norm(profile[name] - profile[name][::-1])
            / np.linalg.norm(profile[name])
        )
    result["reynolds_uv_antisymmetry_l2_relative"] = float(
        np.linalg.norm(profile["reynolds_uv"] + profile["reynolds_uv"][::-1])
        / np.linalg.norm(profile["reynolds_uv"])
    )
    return result


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    statistics = {}
    histories = {}
    segment_statistics = {}
    for label, directories in CASES.items():
        segment_statistics[label] = [
            read_table(next(directory.glob("*.statistics.r64.txt")))
            for directory in directories
        ]
        histories[label] = [
            read_table(next(directory.glob("*.history.r64.txt")))
            for directory in directories
        ]
        statistics[label] = concatenate_tables(segment_statistics[label])

    source_directory = RESULTS / "lowmach-longrun-segment02"
    source_statistics = read_table(
        next(source_directory.glob("*.statistics.r64.txt"))
    )

    common_start = max(table["time"][0] for table in statistics.values())
    common_end = min(table["time"][-1] for table in statistics.values())
    report: dict[str, object] = {
        "common_start": float(common_start),
        "common_end": float(common_end),
        "common_duration": float(common_end - common_start),
        "source_last_20": window_summary(source_statistics, 203.0, common_start),
        "segment02_common_start": float(
            max(tables[1]["time"][0] for tables in segment_statistics.values())
        ),
        "cases": {},
    }
    for label, table in statistics.items():
        initial = {name: float(values[0]) for name, values in table.items()}
        final = {name: float(values[-1]) for name, values in table.items()}
        report["cases"][label] = {
            "initial": initial,
            "final": final,
            "duration": float(table["time"][-1] - table["time"][0]),
            "advanced_steps": int(table["step"][-1] - table["step"][0]),
            "mass_relative_drift": float(
                (table["total_mass"][-1] - table["total_mass"][0])
                / table["total_mass"][0]
            ),
            "momentum_x_relative_change": float(
                (table["total_momentum_x"][-1] - table["total_momentum_x"][0])
                / table["total_momentum_x"][0]
            ),
            "history": history_summary(histories[label]),
            "segments": [
                {
                    "start_step": int(segment["step"][0]),
                    "end_step": int(segment["step"][-1]),
                    "start_time": float(segment["time"][0]),
                    "end_time": float(segment["time"][-1]),
                    "rows": int(segment["time"].size),
                }
                for segment in segment_statistics[label]
            ],
            "restart_boundary": {
                "step_delta": int(
                    segment_statistics[label][1]["step"][0]
                    - segment_statistics[label][0]["step"][-1]
                ),
                "time_delta": float(
                    segment_statistics[label][1]["time"][0]
                    - segment_statistics[label][0]["time"][-1]
                ),
                "value_deltas": {
                    name: float(
                        segment_statistics[label][1][name][0]
                        - segment_statistics[label][0][name][-1]
                    )
                    for name in (
                        "total_mass",
                        "total_momentum_x",
                        "total_energy",
                        "channel_wall_shear_mean",
                        "channel_re_tau",
                    )
                },
            },
            "windows": {
                "full_common": window_summary(table, common_start, common_end),
                "223_245": window_summary(table, common_start, 245.0),
                "segment02_common": window_summary(
                    table, report["segment02_common_start"], common_end
                ),
                "250_260": window_summary(table, 250.0, 260.0),
                "260_270": window_summary(table, 260.0, 270.0),
                "270_common_end": window_summary(table, 270.0, common_end),
                "275_common_end": window_summary(table, 275.0, common_end),
            },
            "five_time_unit_blocks": {
                f"{start}_{start + 5}": window_summary(
                    table, float(start), float(start + 5)
                )
                for start in range(250, 275, 5)
            },
        }

    shared_times = statistics["diss=0.01"]["time"]
    shared_mask = shared_times <= common_end + 1.0e-10
    comparison: dict[str, object] = {}
    for name in (
        "total_momentum_x",
        "total_energy",
        "channel_wall_shear_mean",
        "channel_re_tau",
    ):
        control = statistics["diss=0.01"][name][shared_mask]
        changed = np.interp(
            shared_times[shared_mask],
            statistics["diss=0.001"]["time"],
            statistics["diss=0.001"][name],
        )
        delta = changed - control
        comparison[name] = {
            "delta_at_common_end": float(delta[-1]),
            "relative_delta_at_common_end": float(delta[-1] / control[-1]),
            "max_abs_delta": float(np.max(np.abs(delta))),
        }
    report["pointwise_diss0p001_minus_diss0p01"] = comparison

    profile_directory = OUTPUT / "profiles"
    pooled_by_window = {}
    for window_name, profile_times in PROFILE_WINDOWS.items():
        pooled = {}
        for label, directories in CASES.items():
            stem_fragment = directories[0].name.replace("-", "_")
            paths = []
            for time in profile_times:
                token = f"time2p{time - 200:02d}0000000eP02"
                matches = sorted(profile_directory.glob(f"*{stem_fragment}*{token}*.txt"))
                if len(matches) != 1:
                    raise RuntimeError(
                        f"expected one {label} profile at t={time}, found {len(matches)}"
                    )
                paths.append(matches[0])
            pooled[label] = pooled_profiles(paths)
        pooled_by_window[window_name] = pooled
        profile_report = {
            "snapshot_times": list(profile_times),
            "pooling": (
                f"equal-weight pooled x-z moments over {len(profile_times)} field snapshots"
            ),
            "cases": {
                label: profile_summary(profile) for label, profile in pooled.items()
            },
        }
        changed = pooled["diss=0.001"]
        control = pooled["diss=0.01"]
        profile_report["diss0p001_minus_diss0p01"] = {
            name: {
                "l2_absolute": float(np.linalg.norm(changed[name] - control[name])),
                "l2_relative_to_diss0p01": float(
                    np.linalg.norm(changed[name] - control[name])
                    / np.linalg.norm(control[name])
                ),
                "max_abs_difference": float(
                    np.max(np.abs(changed[name] - control[name]))
                ),
            }
            for name in ("mean_u", "rms_u", "rms_v", "rms_w", "reynolds_uv")
        }
        report[f"profiles_{window_name}"] = profile_report

    report["profile_evolution_225_245_to_250_275"] = {}
    for label in CASES:
        early = pooled_by_window["225_245"][label]
        late = pooled_by_window["250_275"][label]
        report["profile_evolution_225_245_to_250_275"][label] = {
            name: {
                "l2_absolute": float(np.linalg.norm(late[name] - early[name])),
                "l2_relative_to_early": float(
                    np.linalg.norm(late[name] - early[name]) / np.linalg.norm(early[name])
                ),
                "max_abs_difference": float(np.max(np.abs(late[name] - early[name]))),
            }
            for name in ("mean_u", "rms_u", "rms_v", "rms_w", "reynolds_uv")
        }

    (OUTPUT / "comparison_metrics.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    colors = {"diss=0.001": "#d55e00", "diss=0.01": "#0072b2"}
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 9.5), sharex=True)
    source_mask = source_statistics["time"] >= 200.0
    axes[0].plot(
        source_statistics["time"][source_mask],
        (
            source_statistics["total_momentum_x"][source_mask]
            / source_statistics["total_mass"][source_mask]
        ),
        color="0.55",
        linewidth=1.1,
        label="source diss=0.01",
    )
    axes[1].plot(
        source_statistics["time"][source_mask],
        source_statistics["channel_re_tau"][source_mask],
        color="0.55",
        linewidth=1.1,
        label="source diss=0.01",
    )
    source_asymmetry = np.abs(
        source_statistics["channel_wall_shear_lower"][source_mask]
        - source_statistics["channel_wall_shear_upper"][source_mask]
    ) / source_statistics["channel_wall_shear_mean"][source_mask]
    axes[2].plot(
        source_statistics["time"][source_mask],
        100.0 * source_asymmetry,
        color="0.55",
        linewidth=1.1,
        label="source diss=0.01",
    )
    for label, table in statistics.items():
        time = table["time"]
        axes[0].plot(
            time,
            table["total_momentum_x"] / table["total_mass"],
            label=label,
            color=colors[label],
            linewidth=1.6,
        )
        axes[1].plot(
            time,
            table["channel_re_tau"],
            label=label,
            color=colors[label],
            linewidth=1.6,
        )
        asymmetry = np.abs(
            table["channel_wall_shear_lower"]
            - table["channel_wall_shear_upper"]
        ) / table["channel_wall_shear_mean"]
        axes[2].plot(
            time,
            100.0 * asymmetry,
            label=label,
            color=colors[label],
            linewidth=1.6,
        )
    axes[0].set_ylabel(r"$P_x/M$ (bulk proxy)")
    axes[1].set_ylabel(r"instantaneous $Re_\tau$")
    axes[2].set_ylabel("wall-shear asymmetry [%]")
    axes[2].set_xlabel(r"time $tU_{b,0}/h$")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.axvline(common_start, color="0.15", linestyle=":", linewidth=1.0)
        axis.axvline(245.5, color="0.4", linestyle="--", linewidth=0.8)
        axis.axvline(270.0, color="0.65", linestyle=":", linewidth=0.8)
    axes[0].legend(frameon=False, ncol=3)
    fig.suptitle("Case05 SCMM6/MDCD-linear/Roe dissipation branches")
    fig.tight_layout()
    fig.savefig(OUTPUT / "branch_global_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 6.8), sharex=True)
    for label, table in statistics.items():
        time = table["time"]
        axes[0].plot(
            time,
            table["channel_wall_shear_lower"],
            color=colors[label],
            linewidth=1.25,
            label=f"{label}, lower",
        )
        axes[0].plot(
            time,
            table["channel_wall_shear_upper"],
            color=colors[label],
            linewidth=1.25,
            linestyle="--",
            label=f"{label}, upper",
        )
        plane_average = 0.5 * (
            table["yz_mass_flow_x_plane0"]
            + table["yz_mass_flow_x_plane1"]
        )
        axes[1].plot(
            time,
            plane_average,
            color=colors[label],
            linewidth=1.5,
            label=label,
        )
    axes[0].set_ylabel("wall shear")
    axes[1].set_ylabel("mean of two x-plane mass flows")
    axes[1].set_xlabel(r"time $tU_{b,0}/h$")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(OUTPUT / "branch_wall_and_flow_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), sharex=True)
    block_starts = list(range(250, 275, 5))
    block_centers = np.asarray(block_starts, dtype=float) + 2.5
    block_fields = (
        (axes[0, 0], "mass_weighted_bulk_velocity", r"block mean $P_x/M$"),
        (axes[0, 1], "channel_re_tau", r"block mean $Re_\tau$"),
        (axes[1, 0], "channel_wall_shear_mean", "block mean wall shear"),
    )
    for label, table in statistics.items():
        summaries = [
            window_summary(table, float(start), float(start + 5))
            for start in block_starts
        ]
        for axis, field, ylabel in block_fields:
            means = [summary[field]["mean"] for summary in summaries]
            stds = [summary[field]["std"] for summary in summaries]
            axis.errorbar(
                block_centers,
                means,
                yerr=stds,
                color=colors[label],
                marker="o",
                capsize=2,
                linewidth=1.4,
                label=label,
            )
            axis.set_ylabel(ylabel)
        asymmetries = [
            100.0 * summary["wall_shear_relative_asymmetry"]["mean"]
            for summary in summaries
        ]
        axes[1, 1].plot(
            block_centers,
            asymmetries,
            color=colors[label],
            marker="o",
            linewidth=1.4,
            label=label,
        )
    axes[1, 1].set_ylabel("block mean wall asymmetry [%]")
    axes[0, 1].axhline(180.0, color="0.5", linestyle=":", linewidth=1.0)
    for axis in axes.flat:
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False)
    axes[1, 0].set_xlabel(r"block center $tU_{b,0}/h$")
    axes[1, 1].set_xlabel(r"block center $tU_{b,0}/h$")
    fig.suptitle("Five-time-unit block statistics in segment02")
    fig.tight_layout()
    fig.savefig(OUTPUT / "branch_segment02_block_means.png", dpi=180)
    plt.close(fig)

    profile_fields = (
        ("mean_u", r"$\langle u\rangle_{xzt}$"),
        ("rms_u", r"$u'_{\mathrm{rms},xzt}$"),
        ("rms_v", r"$v'_{\mathrm{rms},xzt}$"),
        ("rms_w", r"$w'_{\mathrm{rms},xzt}$"),
    )
    for window_name, pooled in pooled_by_window.items():
        fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.0), sharey=True)
        for axis, (name, xlabel) in zip(axes.flat, profile_fields):
            for label, profile in pooled.items():
                axis.plot(
                    profile[name],
                    profile["y"],
                    color=colors[label],
                    linewidth=1.7,
                    label=label,
                )
            axis.set_xlabel(xlabel)
            axis.grid(True, alpha=0.25)
        axes[0, 0].set_ylabel(r"wall-normal coordinate $y/h$")
        axes[1, 0].set_ylabel(r"wall-normal coordinate $y/h$")
        axes[0, 0].legend(frameon=False)
        times_text = ", ".join(str(value) for value in PROFILE_WINDOWS[window_name])
        fig.suptitle(f"Pooled x-z profiles at t = {times_text}")
        fig.tight_layout()
        fig.savefig(OUTPUT / f"branch_profiles_{window_name}.png", dpi=180)
        plt.close(fig)

        fig, axis = plt.subplots(figsize=(6.2, 5.4))
        for label, profile in pooled.items():
            axis.plot(
                -profile["reynolds_uv"],
                profile["y"],
                color=colors[label],
                linewidth=1.7,
                label=label,
            )
        axis.set_xlabel(r"$-\langle u'v'\rangle_{xzt}$")
        axis.set_ylabel(r"wall-normal coordinate $y/h$")
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False)
        axis.set_title(f"Pooled snapshots: t = {times_text}")
        fig.tight_layout()
        fig.savefig(OUTPUT / f"branch_reynolds_shear_{window_name}.png", dpi=180)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 8.0), sharey=True)
    line_styles = {"225_245": "--", "250_275": "-"}
    for axis, (name, xlabel) in zip(axes.flat, profile_fields):
        for window_name, pooled in pooled_by_window.items():
            for label, profile in pooled.items():
                axis.plot(
                    profile[name],
                    profile["y"],
                    color=colors[label],
                    linestyle=line_styles[window_name],
                    linewidth=1.5,
                    label=f"{label}, t={window_name.replace('_', '--')}",
                )
        axis.set_xlabel(xlabel)
        axis.grid(True, alpha=0.25)
    axes[0, 0].set_ylabel(r"wall-normal coordinate $y/h$")
    axes[1, 0].set_ylabel(r"wall-normal coordinate $y/h$")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Profile evolution from segment01 to segment02")
    fig.tight_layout()
    fig.savefig(OUTPUT / "branch_profile_evolution.png", dpi=180)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build the data products used by the final Case05 LaTeX report.

The report treats the executed manifests and text histories as the source of
truth.  It summarizes every completed segment, verifies restart continuity,
plots the complete campaign, and evaluates selected field windows against the
ERCOFTAC Case 032 Re_tau=180 DNS data.
"""

from __future__ import annotations

import json
import re
from collections import OrderedDict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from analyze_mdcd_dissipation_branches import (
    concatenate_tables,
    pooled_profiles,
    profile_summary,
    read_table,
)
from analyze_segment03_dns import find_profile, make_current_profile, scalar_summary
from compare_branches_to_dns import (
    cell_widths_from_centers,
    comparison_metrics,
    fold_channel,
    load_dns,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = RESULTS / "final-campaign-report"
PROFILE_OUTPUT = OUTPUT / "profiles"
BODY_ACCELERATION = 0.0041720265499730312

SEGMENTS = (
    {
        "id": "F0",
        "group": "feasibility",
        "directory": "lowmach-feasibility-r4",
        "description": "4-rank feasibility gate",
    },
    {
        "id": "B1",
        "group": "baseline",
        "directory": "lowmach-longrun-segment01",
        "description": "baseline segment 01",
    },
    {
        "id": "B2",
        "group": "baseline",
        "directory": "lowmach-longrun-segment02",
        "description": "baseline segment 02",
    },
    {
        "id": "L1",
        "group": "diss=0.001",
        "directory": "branch-scmm6_mdcd0p001_roe-t223",
        "description": "low-dissipation branch segment 01",
    },
    {
        "id": "L2",
        "group": "diss=0.001",
        "directory": "branch-scmm6_mdcd0p001_roe-t223_segment02",
        "description": "low-dissipation branch segment 02",
    },
    {
        "id": "L3",
        "group": "diss=0.001",
        "directory": "branch-scmm6_mdcd0p001_roe-t223_segment03",
        "description": "low-dissipation branch segment 03",
    },
    {
        "id": "L4",
        "group": "diss=0.001",
        "directory": "branch-scmm6_mdcd0p001_roe-t223_segment04",
        "description": "low-dissipation branch segment 04",
    },
    {
        "id": "H1",
        "group": "diss=0.01",
        "directory": "branch-scmm6_mdcd0p01_roe-t223",
        "description": "reference-dissipation branch segment 01",
    },
    {
        "id": "H2",
        "group": "diss=0.01",
        "directory": "branch-scmm6_mdcd0p01_roe-t223_segment02",
        "description": "reference-dissipation branch segment 02",
    },
    {
        "id": "H3",
        "group": "diss=0.01",
        "directory": "branch-scmm6_mdcd0p01_roe-t223_segment03",
        "description": "reference-dissipation branch segment 03",
    },
)

GROUP_SEGMENTS = {
    "baseline": ("B1", "B2"),
    "diss=0.001": ("L1", "L2", "L3", "L4"),
    "diss=0.01": ("H1", "H2", "H3"),
}

COLORS = {
    "baseline": "#5b5b5b",
    "diss=0.001": "#d55e00",
    "diss=0.01": "#0072b2",
}


def parse_manifest(path: Path) -> tuple[dict[str, str], str]:
    text = path.read_text(encoding="utf-8")
    direct: dict[str, str] = {}
    for line in text.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key not in direct:
                direct[key] = value
    return direct, text


def last_match(pattern: str, text: str, default: str = "") -> str:
    matches = re.findall(pattern, text)
    return matches[-1] if matches else default


def segment_files(directory: Path) -> tuple[Path, Path, Path]:
    manifest = next(directory.glob("*.manifest.*.txt"))
    history = next(directory.glob("*.history.*.txt"))
    statistics = next(directory.glob("*.statistics.*.txt"))
    return manifest, history, statistics


def endpoint(table: dict[str, np.ndarray], index: int) -> dict[str, float]:
    return {name: float(values[index]) for name, values in table.items()}


def finite_max(table: dict[str, np.ndarray], name: str) -> int:
    values = np.asarray(table[name])
    return int(np.nanmax(values)) if values.size else 0


def load_segment(spec: dict[str, str]) -> dict[str, object]:
    directory = RESULTS / spec["directory"]
    manifest_path, history_path, statistics_path = segment_files(directory)
    manifest, manifest_text = parse_manifest(manifest_path)
    history = read_table(history_path)
    statistics = read_table(statistics_path)
    positive_dt = history["dt"][history["dt"] > 0.0]
    wall_seconds = float(history["wall_time"][-1] - history["wall_time"][0])
    advanced_steps = int(round(history["step"][-1] - history["step"][0]))
    field_files = sorted(directory.glob("*.field.*.cgns"))
    checkpoint_files = sorted(directory.glob("*.checkpoint.step*.cgns"))
    latest_checkpoints = sorted(directory.glob("*.checkpoint.latest.cgns"))
    all_files = [path for path in directory.iterdir() if path.is_file()]
    signature = manifest.get("restart_signature", "")
    config_summary = manifest.get("config_summary", "")
    reconstruction_matches = re.findall(r"reconstruction=([^;]+)", signature)
    reconstruction = reconstruction_matches[-1] if reconstruction_matches else "unknown"
    compatibility = last_match(
        r"restart\.compatibility=([^,;\)]+)", manifest_text, "strict/default"
    )
    restart_path = last_match(
        r"restart\.path=([^,]+?)(?:,restart\.compatibility|,digest)",
        manifest_text,
        "<none>",
    )
    start = endpoint(statistics, 0)
    end = endpoint(statistics, -1)
    result: dict[str, object] = {
        **spec,
        "path": str(directory.relative_to(ROOT)).replace("\\", "/"),
        "program_version": manifest.get("program_version", "unknown"),
        "git_commit": manifest.get("git_commit", "unknown"),
        "compiler": manifest.get("compiler", "unknown"),
        "mpi_ranks": int(manifest.get("mpi_ranks", 0)),
        "case": manifest.get("case", "unknown"),
        "profile": last_match(r"name=([^;]+)", signature, "unknown"),
        "reconstruction": reconstruction,
        "variables": last_match(r"variables=([^;]+)", signature, "unknown"),
        "mdcd_dispersion": float(last_match(r"mdcd_dispersion=([^;]+)", signature, "nan")),
        "mdcd_dissipation": float(last_match(r"mdcd_dissipation=([^;]+)", signature, "nan")),
        "riemann": last_match(r"riemann_solver=([^;]+)", signature, "unknown"),
        "cfl": float(np.nanmedian(history["cfl"])),
        "restart_path": restart_path,
        "restart_compatibility": compatibility,
        "start": start,
        "end": end,
        "duration": float(end["time"] - start["time"]),
        "advanced_steps": advanced_steps,
        "wall_seconds": wall_seconds,
        "wall_hours": wall_seconds / 3600.0,
        "core_hours": wall_seconds * int(manifest.get("mpi_ranks", 0)) / 3600.0,
        "steps_per_wall_second": advanced_steps / wall_seconds if wall_seconds > 0 else float("nan"),
        "stop_reason": manifest.get("stop_reason", "unknown"),
        "dt_min_positive": float(np.min(positive_dt)),
        "dt_median_positive": float(np.median(positive_dt)),
        "dt_max": float(np.max(history["dt"])),
        "mass_relative_drift": (end["total_mass"] - start["total_mass"]) / start["total_mass"],
        "energy_relative_change": (end["total_energy"] - start["total_energy"]) / start["total_energy"],
        "reconstruction_fallbacks": finite_max(history, "reconstruction_fallbacks"),
        "riemann_fallbacks": finite_max(history, "riemann_fallbacks"),
        "troubled_cells": finite_max(history, "troubled_cells"),
        "local_recomputations": finite_max(history, "local_recomputations"),
        "step_retries": finite_max(history, "step_retries"),
        "field_count": len(field_files),
        "checkpoint_count": len(checkpoint_files),
        "latest_checkpoint_count": len(latest_checkpoints),
        "file_count": len(all_files),
        "bytes": sum(path.stat().st_size for path in all_files),
        "manifest_path": str(manifest_path.relative_to(ROOT)).replace("\\", "/"),
        "history_path": str(history_path.relative_to(ROOT)).replace("\\", "/"),
        "statistics_path": str(statistics_path.relative_to(ROOT)).replace("\\", "/"),
        "_history": history,
        "_statistics": statistics,
    }
    return result


def clean_segment(record: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in record.items() if not key.startswith("_")}


def continuity(previous: dict[str, object], current: dict[str, object]) -> dict[str, object]:
    previous_end = previous["end"]
    current_start = current["start"]
    names = (
        "step",
        "time",
        "total_mass",
        "total_momentum_x",
        "total_momentum_y",
        "total_momentum_z",
        "total_energy",
        "channel_wall_shear_lower",
        "channel_wall_shear_upper",
        "channel_re_tau",
    )
    differences = {
        name: float(current_start[name] - previous_end[name]) for name in names
    }
    return {
        "from": previous["id"],
        "to": current["id"],
        "restart_compatibility": current["restart_compatibility"],
        "differences": differences,
        "max_absolute_difference": max(abs(value) for value in differences.values()),
    }


def family_table(records: dict[str, dict[str, object]], ids: tuple[str, ...], kind: str) -> dict[str, np.ndarray]:
    key = "_statistics" if kind == "statistics" else "_history"
    return concatenate_tables([records[segment_id][key] for segment_id in ids])


def generic_current_profile(
    paths: list[Path], times: tuple[int, ...], statistics: dict[str, np.ndarray]
) -> tuple[dict[str, np.ndarray | float], dict[str, float]]:
    pooled = pooled_profiles(paths)
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
        np.mean(np.interp(sample_times, statistics["time"], statistics["channel_re_tau"]))
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
    current: dict[str, np.ndarray | float] = {
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
    return current, profile_summary(pooled)


def baseline_profile_paths(times: tuple[int, ...]) -> list[Path]:
    paths = []
    for time in times:
        token_map = {
            195: "time1p950000000eP02",
            200: "time2p000000000eP02",
            205: "time2p050000000eP02",
            210: "time2p100000000eP02",
            215: "time2p150000000eP02",
            220: "time2p200000000eP02",
        }
        matches = list(PROFILE_OUTPUT.glob(f"*{token_map[time]}*.profile.txt"))
        if len(matches) != 1:
            raise RuntimeError(f"expected one baseline profile at t={time}, found {len(matches)}")
        paths.append(matches[0])
    return paths


def profile_to_json(profile: dict[str, np.ndarray | float]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in profile.items():
        if isinstance(value, np.ndarray):
            result[key] = value.tolist()
        else:
            result[key] = float(value)
    return result


def plot_timeline(records: dict[str, dict[str, object]]) -> None:
    fig, axis = plt.subplots(figsize=(11.2, 4.0))
    lanes = {"baseline": 2, "diss=0.001": 1, "diss=0.01": 0}
    for group, ids in GROUP_SEGMENTS.items():
        for segment_id in ids:
            record = records[segment_id]
            start = float(record["start"]["time"])
            end = float(record["end"]["time"])
            axis.barh(
                lanes[group],
                end - start,
                left=start,
                height=0.52,
                color=COLORS[group],
                alpha=0.88,
                edgecolor="white",
            )
            axis.text(
                0.5 * (start + end),
                lanes[group],
                segment_id,
                ha="center",
                va="center",
                color="white",
                fontsize=9,
                fontweight="bold",
            )
    axis.axvline(223.00310784903166, color="black", linestyle="--", linewidth=1.2)
    axis.text(225.0, 2.20, "fork at t=223.0031", fontsize=8, va="top")
    axis.set_yticks([0, 1, 2], ["diss=0.01", "diss=0.001", "baseline"])
    axis.set_xlabel(r"nondimensional time $tU_{b,0}/h$")
    axis.set_xlim(-3.0, 447.0)
    axis.grid(True, axis="x", alpha=0.25)
    axis.set_title("Case05 computation campaign and completed segments")
    fig.tight_layout()
    fig.savefig(OUTPUT / "campaign_timeline.png", dpi=220)
    plt.close(fig)


def plot_campaign_histories(families: dict[str, dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(4, 1, figsize=(11.2, 11.0), sharex=True)
    for label, table in families.items():
        bulk = table["total_momentum_x"] / table["total_mass"]
        shear_mismatch = 100.0 * (table["channel_wall_shear_mean"] / BODY_ACCELERATION - 1.0)
        asymmetry = 100.0 * np.abs(
            table["channel_wall_shear_lower"] - table["channel_wall_shear_upper"]
        ) / table["channel_wall_shear_mean"]
        axes[0].plot(table["time"], bulk, color=COLORS[label], linewidth=1.15, label=label)
        axes[1].plot(
            table["time"], table["channel_re_tau"], color=COLORS[label], linewidth=1.05, label=label
        )
        axes[2].plot(table["time"], shear_mismatch, color=COLORS[label], linewidth=1.0, label=label)
        axes[3].plot(table["time"], asymmetry, color=COLORS[label], linewidth=0.95, label=label)
    for axis in axes:
        axis.axvline(223.00310784903166, color="0.25", linestyle="--", linewidth=0.9)
        axis.grid(True, alpha=0.24)
        axis.legend(frameon=False, ncol=3, fontsize=8)
    axes[1].axhline(180.0, color="black", linestyle=":", linewidth=1.0)
    axes[2].axhline(0.0, color="black", linestyle=":", linewidth=1.0)
    axes[0].set_ylabel(r"$P_x/M$")
    axes[1].set_ylabel(r"$Re_\tau$")
    axes[2].set_ylabel(r"$(\bar\tau_w/a_x-1)\times100\%$")
    axes[3].set_ylabel("wall asymmetry (%)")
    axes[3].set_xlabel(r"nondimensional time $tU_{b,0}/h$")
    axes[3].set_ylim(bottom=0.0)
    fig.suptitle("Complete integral evolution: baseline and both branches")
    fig.tight_layout()
    fig.savefig(OUTPUT / "campaign_integral_histories.png", dpi=220)
    plt.close(fig)


def plot_conservation(families: dict[str, dict[str, np.ndarray]]) -> None:
    reference_mass = float(families["baseline"]["total_mass"][0])
    reference_energy = float(families["baseline"]["total_energy"][0])
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.0), sharex=True)
    for label, table in families.items():
        mass_scaled = (table["total_mass"] / reference_mass - 1.0) * 1.0e12
        energy_percent = 100.0 * (table["total_energy"] / reference_energy - 1.0)
        py = table["total_momentum_y"] / table["total_mass"]
        pz = table["total_momentum_z"] / table["total_mass"]
        u_plane_error = 100.0 * np.abs(
            table["yz_mean_u_plane0"] - table["yz_mean_u_plane1"]
        ) / np.maximum(
            0.5 * np.abs(table["yz_mean_u_plane0"] + table["yz_mean_u_plane1"]), 1.0e-30
        )
        flow_plane_error = 100.0 * np.abs(
            table["yz_mass_flow_x_plane0"] - table["yz_mass_flow_x_plane1"]
        ) / np.maximum(
            0.5
            * np.abs(table["yz_mass_flow_x_plane0"] + table["yz_mass_flow_x_plane1"]),
            1.0e-30,
        )
        axes[0, 0].plot(table["time"], mass_scaled, color=COLORS[label], linewidth=1.0, label=label)
        axes[0, 1].plot(table["time"], energy_percent, color=COLORS[label], linewidth=1.0, label=label)
        axes[1, 0].plot(table["time"], py, color=COLORS[label], linewidth=0.9, label=f"{label}: Py/M")
        axes[1, 0].plot(table["time"], pz, color=COLORS[label], linewidth=0.9, linestyle="--", label=f"{label}: Pz/M")
        axes[1, 1].plot(table["time"], u_plane_error, color=COLORS[label], linewidth=0.9, label=f"{label}: mean U")
        axes[1, 1].plot(table["time"], flow_plane_error, color=COLORS[label], linewidth=0.9, linestyle="--", label=f"{label}: mass flow")
    axes[0, 0].set_ylabel(r"$(M/M_0-1)\times10^{12}$")
    axes[0, 1].set_ylabel(r"$(E/E_0-1)\times100\%$")
    axes[1, 0].set_ylabel("transverse momentum / mass")
    axes[1, 1].set_ylabel("two-plane mismatch (%)")
    for axis in axes.flat:
        axis.axvline(223.00310784903166, color="0.25", linestyle=":", linewidth=0.8)
        axis.grid(True, alpha=0.24)
        axis.set_xlabel(r"$tU_{b,0}/h$")
        axis.legend(frameon=False, fontsize=6.8, ncol=2)
    axes[1, 0].set_yscale("symlog", linthresh=1.0e-5)
    axes[1, 1].set_yscale("log")
    fig.suptitle("Conservation, transverse momentum, and streamwise-plane consistency")
    fig.tight_layout()
    fig.savefig(OUTPUT / "campaign_conservation.png", dpi=220)
    plt.close(fig)


def plot_numerical_health(
    histories: dict[str, dict[str, np.ndarray]], records: dict[str, dict[str, object]]
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.0))
    for label, table in histories.items():
        mask = table["dt"] > 0.0
        axes[0, 0].plot(table["time"][mask], table["dt"][mask], color=COLORS[label], linewidth=0.9, label=label)
        axes[0, 1].plot(table["time"], table["total_l2"], color=COLORS[label], linewidth=0.85, label=label)
    axes[0, 0].set_ylabel(r"accepted $\Delta t$")
    axes[0, 1].set_ylabel("unsteady update total L2")
    axes[0, 1].set_yscale("log")
    segment_ids = [spec["id"] for spec in SEGMENTS if spec["id"] != "F0"]
    x = np.arange(len(segment_ids))
    rates = [float(records[segment_id]["steps_per_wall_second"]) for segment_id in segment_ids]
    wall_hours = [float(records[segment_id]["wall_hours"]) for segment_id in segment_ids]
    bar_colors = [COLORS[records[segment_id]["group"]] for segment_id in segment_ids]
    axes[1, 0].bar(x, rates, color=bar_colors)
    axes[1, 1].bar(x, wall_hours, color=bar_colors)
    axes[1, 0].set_xticks(x, segment_ids)
    axes[1, 1].set_xticks(x, segment_ids)
    axes[1, 0].set_ylabel("advanced steps / wall second")
    axes[1, 1].set_ylabel("wall time per segment (h)")
    for axis in axes.flat:
        axis.grid(True, alpha=0.24)
        axis.legend(frameon=False, fontsize=8) if axis in axes[0, :] else None
    axes[0, 0].set_xlabel(r"$tU_{b,0}/h$")
    axes[0, 1].set_xlabel(r"$tU_{b,0}/h$")
    fig.suptitle("Time-step behavior and computational throughput")
    fig.tight_layout()
    fig.savefig(OUTPUT / "campaign_numerical_health.png", dpi=220)
    plt.close(fig)


def plot_dns_profiles(
    dns: dict[str, np.ndarray | float], cases: OrderedDict[str, dict[str, np.ndarray | float]]
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 8.0))
    axes = axes.flat
    fields = (
        ("mean_u_plus", r"$U^+$", True),
        ("rms_u_plus", r"$u_{rms}^+$", False),
        ("rms_v_plus", r"$v_{rms}^+$", False),
        ("rms_w_plus", r"$w_{rms}^+$", False),
        ("minus_uv_plus", r"$-\langle u'v'\rangle^+$", False),
    )
    colors = (COLORS["baseline"], COLORS["diss=0.001"], COLORS["diss=0.01"])
    for axis, (name, ylabel, logarithmic) in zip(axes[:5], fields):
        if logarithmic:
            axis.semilogx(dns["y_plus"], dns[name], "ko", markersize=2.7, label="DNS")
        else:
            axis.plot(dns["y_plus"], dns[name], "ko", markersize=2.7, label="DNS")
        for (label, current), color in zip(cases.items(), colors):
            if logarithmic:
                axis.semilogx(current["y_plus"], current[name], color=color, linewidth=1.6, label=label)
            else:
                axis.plot(current["y_plus"], current[name], color=color, linewidth=1.6, label=label)
        axis.set_xlabel(r"$y^+$")
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", alpha=0.24)
        axis.legend(frameon=False, fontsize=7)
    outer = axes[5]
    outer.plot(
        np.asarray(dns["mean_u_plus"]) / float(dns["bulk_velocity_plus"]),
        dns["eta"],
        "ko",
        markersize=2.7,
        label="DNS",
    )
    for (label, current), color in zip(cases.items(), colors):
        outer.plot(current["mean_u_over_bulk"], current["eta"], color=color, linewidth=1.6, label=label)
    outer.set_xlabel(r"$U/U_b$")
    outer.set_ylabel(r"wall distance $y/h$")
    outer.grid(True, alpha=0.24)
    outer.legend(frameon=False, fontsize=7)
    fig.suptitle("Final available profile windows compared with ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_final_profile_overview.png", dpi=220)
    plt.close(fig)


def plot_dns_error_matrix(metrics: OrderedDict[str, dict[str, object]]) -> None:
    columns = (
        ("Retau", lambda item: abs(float(item["re_tau_relative_error"]))),
        ("Ub+", lambda item: abs(float(item["bulk_velocity_plus_relative_error"]))),
        ("Cf", lambda item: abs(float(item["skin_friction_coefficient_relative_error"]))),
        ("U+", lambda item: item["inner_scaled_profile_errors_on_common_y_plus"]["mean_u_plus"]["relative_l2"]),
        ("u_rms+", lambda item: item["inner_scaled_profile_errors_on_common_y_plus"]["rms_u_plus"]["relative_l2"]),
        ("v_rms+", lambda item: item["inner_scaled_profile_errors_on_common_y_plus"]["rms_v_plus"]["relative_l2"]),
        ("w_rms+", lambda item: item["inner_scaled_profile_errors_on_common_y_plus"]["rms_w_plus"]["relative_l2"]),
        ("-uv+", lambda item: item["inner_scaled_profile_errors_on_common_y_plus"]["minus_uv_plus"]["relative_l2"]),
        ("U/Ub", lambda item: item["outer_scaled_mean_velocity_error"]["relative_l2"]),
    )
    values = 100.0 * np.asarray([[extractor(item) for _, extractor in columns] for item in metrics.values()])
    display = np.maximum(values, 0.08)
    fig, axis = plt.subplots(figsize=(11.4, 6.5))
    image = axis.imshow(
        display,
        aspect="auto",
        cmap="viridis",
        norm=mpl.colors.LogNorm(vmin=0.08, vmax=max(40.0, float(np.max(display)))),
    )
    axis.set_xticks(np.arange(len(columns)), [name for name, _ in columns])
    axis.set_yticks(np.arange(len(metrics)), list(metrics.keys()))
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            color = "white" if display[row, column] > 5.0 else "black"
            axis.text(column, row, f"{values[row, column]:.2f}", ha="center", va="center", fontsize=7.3, color=color)
    colorbar = fig.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label("absolute relative error (%) - logarithmic color scale")
    axis.set_title("DNS error matrix for all analyzed six-snapshot windows")
    axis.set_xlabel("quantity")
    axis.set_ylabel("case and field-snapshot window")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_error_matrix.png", dpi=220)
    plt.close(fig)


def f(value: float, digits: int = 5) -> str:
    return f"{value:.{digits}g}"


def sci(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}e}"


def tex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "_": r"\_",
        "%": r"\%",
        "&": r"\&",
        "#": r"\#",
    }
    result = text
    for source, target in replacements.items():
        result = result.replace(source, target)
    return result


def write_generated_tex(
    records: dict[str, dict[str, object]],
    continuities: list[dict[str, object]],
    final_windows: dict[str, dict[str, object]],
    dns_metrics: OrderedDict[str, dict[str, object]],
    totals: dict[str, float | int],
) -> None:
    lines: list[str] = []
    lines.extend(
        [
            "% This file is generated by generate_case05_final_report_data.py.",
            rf"\providecommand{{\CampaignAdvancedSteps}}{{{int(totals['advanced_steps']):,}}}",
            rf"\providecommand{{\CampaignWallHours}}{{{totals['wall_hours']:.2f}}}",
            rf"\providecommand{{\CampaignCoreHours}}{{{totals['core_hours']:.1f}}}",
            rf"\providecommand{{\CampaignStorageGiB}}{{{totals['storage_gib']:.2f}}}",
            rf"\providecommand{{\CampaignFieldCount}}{{{int(totals['field_count'])}}}",
            rf"\providecommand{{\CampaignCheckpointCount}}{{{int(totals['checkpoint_count'])}}}",
            "",
            r"\begin{landscape}",
            r"\small",
            r"\begin{longtable}{@{}llrrrrrl@{}}",
            r"\caption{全部计算段的时间范围、算法与终止信息。}\label{tab:segments}\\",
            r"\toprule",
            r"编号 & 算法简写 & $t_0$ & $t_1$ & $\Delta t$ & 推进步数 & MPI & 终止原因 \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"编号 & 算法简写 & $t_0$ & $t_1$ & $\Delta t$ & 推进步数 & MPI & 终止原因 \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for spec in SEGMENTS:
        record = records[spec["id"]]
        algorithm = f"{record['profile']}/{record['reconstruction']}/{record['riemann']}, d={record['mdcd_dissipation']:g}"
        lines.append(
            f"{record['id']} & {tex_escape(algorithm)} & {record['start']['time']:.6f} & "
            f"{record['end']['time']:.6f} & {record['duration']:.6f} & "
            f"{record['advanced_steps']:,} & {record['mpi_ranks']} & "
            f"{tex_escape(str(record['stop_reason']))} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\end{landscape}", ""])

    lines.extend(
        [
            r"\begin{landscape}",
            r"\begin{longtable}{@{}lrrrrrrrr@{}}",
            r"\caption{各计算段的计算代价和输出文件清单。检查点数不含 latest 别名文件。}\label{tab:performance}\\",
            r"\toprule",
            r"编号 & 墙钟/h & 核时/h & step/s & $\Delta t_{min}$ & $\Delta t_{50}$ & 场文件 & 检查点 & 容量/GiB \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"编号 & 墙钟/h & 核时/h & step/s & $\Delta t_{min}$ & $\Delta t_{50}$ & 场文件 & 检查点 & 容量/GiB \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for spec in SEGMENTS:
        record = records[spec["id"]]
        lines.append(
            f"{record['id']} & {record['wall_hours']:.3f} & {record['core_hours']:.1f} & "
            f"{record['steps_per_wall_second']:.3f} & {sci(record['dt_min_positive'])} & "
            f"{sci(record['dt_median_positive'])} & {record['field_count']} & "
            f"{record['checkpoint_count']} & {record['bytes'] / 2**30:.3f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\end{landscape}", ""])

    lines.extend(
        [
            r"\begin{longtable}{@{}llrrrrr@{}}",
            r"\caption{续算与分叉连接点连续性。状态差由当前段首行减前一段末行得到。}\label{tab:continuity}\\",
            r"\toprule",
            r"连接 & 模式 & $\Delta t$ & $\Delta M$ & $\Delta P_x$ & $\Delta E$ & $\Delta Re_\tau$ \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"连接 & 模式 & $\Delta t$ & $\Delta M$ & $\Delta P_x$ & $\Delta E$ & $\Delta Re_\tau$ \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for item in continuities:
        delta = item["differences"]
        lines.append(
            f"{item['from']}$\\to${item['to']} & {tex_escape(str(item['restart_compatibility']))} & "
            f"{sci(delta['time'])} & {sci(delta['total_mass'])} & {sci(delta['total_momentum_x'])} & "
            f"{sci(delta['total_energy'])} & {sci(delta['channel_re_tau'])} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", ""])

    lines.extend(
        [
            r"\begin{table}[H]",
            r"\centering\small",
            r"\caption{各主轨迹最后20个无量纲时间单位的判稳指标。}\label{tab:last20}",
            r"\begin{tabular}{@{}lrrrrrr@{}}",
            r"\toprule",
            r"轨迹 & $\overline U_b$ & $dU_b/dt$ & $\overline{Re_\tau}$ & $dRe_\tau/dt$ & 壁剪不对称 & 壁剪失配 \\",
            r"\midrule",
        ]
    )
    for label, window in final_windows.items():
        mismatch = float(window["wall_shear_mean"]["mean"]) / BODY_ACCELERATION - 1.0
        lines.append(
            f"{tex_escape(label)} & {window['bulk_velocity']['mean']:.6f} & "
            f"{sci(window['bulk_velocity']['slope_per_time'])} & {window['re_tau']['mean']:.3f} & "
            f"{sci(window['re_tau']['slope_per_time'])} & "
            f"{100.0 * window['wall_shear_relative_asymmetry']['mean']:.2f}\\% & "
            f"{100.0 * mismatch:+.2f}\\% \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])

    lines.extend(
        [
            r"\begin{landscape}",
            r"\begin{longtable}{@{}lrrrrrrrrr@{}}",
            r"\caption{各六场窗口相对 ERCOFTAC Case 032 DNS 的误差。前三项为带符号相对误差，其余为相对 $L_2$ 误差。}\label{tab:dns-errors}\\",
            r"\toprule",
            r"窗口 & $Re_\tau$ & $U_b^+$ & $C_f$ & $U^+$ & $u_{rms}^+$ & $v_{rms}^+$ & $w_{rms}^+$ & $-\langle u'v'\rangle^+$ & $U/U_b$ \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"窗口 & $Re_\tau$ & $U_b^+$ & $C_f$ & $U^+$ & $u_{rms}^+$ & $v_{rms}^+$ & $w_{rms}^+$ & $-\langle u'v'\rangle^+$ & $U/U_b$ \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for label, metric in dns_metrics.items():
        profile_errors = metric["inner_scaled_profile_errors_on_common_y_plus"]
        lines.append(
            f"{tex_escape(label)} & {100.0 * metric['re_tau_relative_error']:+.2f}\\% & "
            f"{100.0 * metric['bulk_velocity_plus_relative_error']:+.2f}\\% & "
            f"{100.0 * metric['skin_friction_coefficient_relative_error']:+.2f}\\% & "
            f"{100.0 * profile_errors['mean_u_plus']['relative_l2']:.2f}\\% & "
            f"{100.0 * profile_errors['rms_u_plus']['relative_l2']:.2f}\\% & "
            f"{100.0 * profile_errors['rms_v_plus']['relative_l2']:.2f}\\% & "
            f"{100.0 * profile_errors['rms_w_plus']['relative_l2']:.2f}\\% & "
            f"{100.0 * profile_errors['minus_uv_plus']['relative_l2']:.2f}\\% & "
            f"{100.0 * metric['outer_scaled_mean_velocity_error']['relative_l2']:.2f}\\% \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\end{landscape}", ""])

    lines.extend(
        [
            r"\begin{landscape}",
            r"\begin{longtable}{@{}lrrrrrr@{}}",
            r"\caption{所有计算段终点的整体守恒量与主要壁面量。$U_b=P_x/M$。}\label{tab:endpoints-a}\\",
            r"\toprule",
            r"编号 & $M$ & $U_b$ & $P_y/M$ & $P_z/M$ & $E$ & $Re_\tau$ \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"编号 & $M$ & $U_b$ & $P_y/M$ & $P_z/M$ & $E$ & $Re_\tau$ \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for spec in SEGMENTS:
        record = records[spec["id"]]
        end = record["end"]
        lines.append(
            f"{record['id']} & {end['total_mass']:.10f} & {end['total_momentum_x']/end['total_mass']:.8f} & "
            f"{end['total_momentum_y']/end['total_mass']:.3e} & {end['total_momentum_z']/end['total_mass']:.3e} & "
            f"{end['total_energy']:.8f} & {end['channel_re_tau']:.5f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\end{landscape}", ""])

    lines.extend(
        [
            r"\begin{landscape}",
            r"\begin{longtable}{@{}lrrrrrrr@{}}",
            r"\caption{所有计算段终点的截面与壁面统计量。}\label{tab:endpoints-b}\\",
            r"\toprule",
            r"编号 & $U_{yz,0}$ & $U_{yz,1}$ & $\dot m_{x,0}$ & $\dot m_{x,1}$ & $\tau_{w,l}$ & $\tau_{w,u}$ & $u_\tau$ \\",
            r"\midrule",
            r"\endfirsthead",
            r"\toprule",
            r"编号 & $U_{yz,0}$ & $U_{yz,1}$ & $\dot m_{x,0}$ & $\dot m_{x,1}$ & $\tau_{w,l}$ & $\tau_{w,u}$ & $u_\tau$ \\",
            r"\midrule",
            r"\endhead",
        ]
    )
    for spec in SEGMENTS:
        record = records[spec["id"]]
        end = record["end"]
        lines.append(
            f"{record['id']} & {end['yz_mean_u_plane0']:.8f} & {end['yz_mean_u_plane1']:.8f} & "
            f"{end['yz_mass_flow_x_plane0']:.8f} & {end['yz_mass_flow_x_plane1']:.8f} & "
            f"{end['channel_wall_shear_lower']:.7f} & {end['channel_wall_shear_upper']:.7f} & "
            f"{end['channel_friction_velocity']:.7f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\end{landscape}", ""])

    lines.extend(
        [
            r"\begin{table}[H]",
            r"\centering\small",
            r"\caption{逐段数值健康检查。全部 fallback、troubled cell、局部重算和 step retry 均为零。}\label{tab:health}",
            r"\begin{tabular}{@{}lrrrrr@{}}",
            r"\toprule",
            r"编号 & $\Delta M/M_0$ & $\Delta E/E_0$ & CFL & 程序版本 & Git提交 \\",
            r"\midrule",
        ]
    )
    for spec in SEGMENTS:
        record = records[spec["id"]]
        lines.append(
            f"{record['id']} & {sci(record['mass_relative_drift'])} & {sci(record['energy_relative_change'])} & "
            f"{record['cfl']:.2f} & {tex_escape(str(record['program_version']))} & "
            f"\\texttt{{{tex_escape(str(record['git_commit']))}}} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])

    (OUTPUT / "generated_report_data.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    records = {spec["id"]: load_segment(spec) for spec in SEGMENTS}
    statistics = {
        label: family_table(records, ids, "statistics") for label, ids in GROUP_SEGMENTS.items()
    }
    histories = {
        label: family_table(records, ids, "history") for label, ids in GROUP_SEGMENTS.items()
    }

    continuity_pairs = (
        ("B1", "B2"),
        ("B2", "L1"),
        ("L1", "L2"),
        ("L2", "L3"),
        ("L3", "L4"),
        ("B2", "H1"),
        ("H1", "H2"),
        ("H2", "H3"),
    )
    continuities = [continuity(records[left], records[right]) for left, right in continuity_pairs]

    final_windows: dict[str, dict[str, object]] = {}
    for label, table in statistics.items():
        end = float(table["time"][-1])
        final_windows[label] = scalar_summary(table, end - 20.0, end)

    dns = load_dns()
    baseline_times = (195, 200, 205, 210, 215, 220)
    baseline_current, baseline_quality = generic_current_profile(
        baseline_profile_paths(baseline_times), baseline_times, statistics["baseline"]
    )

    window_specs = OrderedDict(
        [
            ("baseline 195-220", ("baseline", baseline_times)),
            ("d=0.001 225-245", ("diss=0.001", (225, 230, 235, 240, 245))),
            ("d=0.001 250-275", ("diss=0.001", (250, 255, 260, 265, 270, 275))),
            ("d=0.001 300-325", ("diss=0.001", (300, 305, 310, 315, 320, 325))),
            ("d=0.001 365-390", ("diss=0.001", (365, 370, 375, 380, 385, 390))),
            ("d=0.001 410-435", ("diss=0.001", (410, 415, 420, 425, 430, 435))),
            ("d=0.01 225-245", ("diss=0.01", (225, 230, 235, 240, 245))),
            ("d=0.01 250-275", ("diss=0.01", (250, 255, 260, 265, 270, 275))),
            ("d=0.01 300-325", ("diss=0.01", (300, 305, 310, 315, 320, 325))),
        ]
    )
    stems = {
        "diss=0.001": "branch_scmm6_mdcd0p001_roe_t223",
        "diss=0.01": "branch_scmm6_mdcd0p01_roe_t223",
    }
    profiles: OrderedDict[str, dict[str, np.ndarray | float]] = OrderedDict()
    qualities: OrderedDict[str, dict[str, float]] = OrderedDict()
    profiles["baseline 195-220"] = baseline_current
    qualities["baseline 195-220"] = baseline_quality
    for label, (family, times) in list(window_specs.items())[1:]:
        current = make_current_profile(stems[family], times, statistics[family])
        profiles[label] = current
        quality_paths = [find_profile(stems[family], time) for time in times]
        qualities[label] = profile_summary(pooled_profiles(quality_paths))

    dns_metrics: OrderedDict[str, dict[str, object]] = OrderedDict(
        (label, comparison_metrics(profile, dns)) for label, profile in profiles.items()
    )

    totals: dict[str, float | int] = {
        "advanced_steps": sum(int(record["advanced_steps"]) for record in records.values()),
        "wall_hours": sum(float(record["wall_hours"]) for record in records.values()),
        "core_hours": sum(float(record["core_hours"]) for record in records.values()),
        "field_count": sum(int(record["field_count"]) for record in records.values()),
        "checkpoint_count": sum(int(record["checkpoint_count"]) for record in records.values()),
        "storage_gib": sum(int(record["bytes"]) for record in records.values()) / 2**30,
    }

    plot_timeline(records)
    plot_campaign_histories(statistics)
    plot_conservation(statistics)
    plot_numerical_health(histories, records)
    final_profile_cases = OrderedDict(
        [
            ("baseline, t=195-220", profiles["baseline 195-220"]),
            ("diss=0.001, t=410-435", profiles["d=0.001 410-435"]),
            ("diss=0.01, t=300-325", profiles["d=0.01 300-325"]),
        ]
    )
    plot_dns_profiles(dns, final_profile_cases)
    plot_dns_error_matrix(dns_metrics)
    write_generated_tex(records, continuities, final_windows, dns_metrics, totals)

    output = {
        "reference": {
            "name": "ERCOFTAC Case 032 simul1.dat / Kim-Moin-Moser Re_tau=180",
            "re_tau": 180.0,
            "bulk_velocity_plus": float(dns["bulk_velocity_plus"]),
            "skin_friction_coefficient": float(dns["skin_friction_coefficient"]),
        },
        "campaign_totals": totals,
        "segments": [clean_segment(records[spec["id"]]) for spec in SEGMENTS],
        "continuity": continuities,
        "final_twenty_time_unit_windows": final_windows,
        "dns_windows": {
            label: {
                "profile_times": list(window_specs[label][1]),
                "quality": qualities[label],
                "metrics": dns_metrics[label],
                "profile": profile_to_json(profiles[label]),
            }
            for label in window_specs
        },
    }
    (OUTPUT / "case05_final_campaign_metrics.json").write_text(
        json.dumps(output, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

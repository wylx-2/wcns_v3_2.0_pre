#!/usr/bin/env python3
"""Compare the new HLLC branch with the completed Roe branches and DNS.

The strict algorithm comparison uses the common six-snapshot window
t=275,280,...,300.  A second, explicitly non-common comparison places the
new HLLC result beside the latest windows from the final campaign report.
"""

from __future__ import annotations

import json
import re
import subprocess
from collections import OrderedDict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_mdcd_dissipation_branches import (
    concatenate_tables,
    pooled_profiles,
    profile_summary,
    read_table,
)
from analyze_segment03_dns import scalar_summary
from compare_branches_to_dns import (
    cell_widths_from_centers,
    comparison_metrics,
    fold_channel,
    load_dns,
)


ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[2]
RESULTS = ROOT / "results"
OUTPUT = RESULTS / "branch-riemann-comparison-hllc-vs-roe"
PROFILE_OUTPUT = OUTPUT / "profiles"
EXISTING_PROFILE_DIRECTORY = RESULTS / "branch-mdcd-dissipation-comparison" / "profiles"
FINAL_METRICS = RESULTS / "final-campaign-report" / "case05_final_campaign_metrics.json"

PROFILE_TIMES = (275, 280, 285, 290, 295, 300)
FINAL_SCALAR_WINDOW = (280.0, 300.0)
DNS_RE_TAU = 180.0

CASES = OrderedDict(
    [
        (
            "HLLC, d=0.001, CFL=0.6",
            {
                "short": "HLLC d=0.001",
                "directories": ("branch-scmm6_mdcd0p001_hllc-t223",),
                "profile_stem": "branch_scmm6_mdcd0p001_hllc_t223",
                "color": "#009e73",
                "linestyle": "-",
                "hllc": True,
            },
        ),
        (
            "Roe, d=0.001, CFL=0.3",
            {
                "short": "Roe d=0.001",
                "directories": (
                    "branch-scmm6_mdcd0p001_roe-t223",
                    "branch-scmm6_mdcd0p001_roe-t223_segment02",
                    "branch-scmm6_mdcd0p001_roe-t223_segment03",
                    "branch-scmm6_mdcd0p001_roe-t223_segment04",
                ),
                "profile_stem": "branch_scmm6_mdcd0p001_roe_t223",
                "color": "#d55e00",
                "linestyle": "--",
                "hllc": False,
            },
        ),
        (
            "Roe, d=0.01, CFL=0.3",
            {
                "short": "Roe d=0.01",
                "directories": (
                    "branch-scmm6_mdcd0p01_roe-t223",
                    "branch-scmm6_mdcd0p01_roe-t223_segment02",
                    "branch-scmm6_mdcd0p01_roe-t223_segment03",
                ),
                "profile_stem": "branch_scmm6_mdcd0p01_roe_t223",
                "color": "#0072b2",
                "linestyle": "-.",
                "hllc": False,
            },
        ),
    ]
)


def find_single(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {pattern} in {directory}, found {len(matches)}")
    return matches[0]


def parse_time(path: Path) -> float | None:
    match = re.search(r"time([0-9]+)p([0-9]+)eP([0-9]+)", path.name)
    if match is None:
        return None
    return float(f"{match.group(1)}.{match.group(2)}") * 10.0 ** int(match.group(3))


def find_time_file(directory: Path, pattern: str, target: float) -> Path:
    matches = [
        path
        for path in directory.glob(pattern)
        if (value := parse_time(path)) is not None and abs(value - target) < 1.0e-7
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one {pattern} at t={target:g} in {directory}, found {len(matches)}"
        )
    return matches[0]


def extractor_path() -> Path:
    candidates = (
        REPOSITORY / "build-rc-serial" / "extract_channel_profile.exe",
        REPOSITORY / "build-rc-serial" / "extract_channel_profile",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("extract_channel_profile executable was not found")


def ensure_hllc_profiles() -> list[Path]:
    PROFILE_OUTPUT.mkdir(parents=True, exist_ok=True)
    directory = RESULTS / CASES["HLLC, d=0.001, CFL=0.6"]["directories"][0]
    executable = extractor_path()
    paths: list[Path] = []
    for time in PROFILE_TIMES:
        field = find_time_file(directory, "*.field.*.cgns", float(time))
        profile = PROFILE_OUTPUT / field.name.replace(".cgns", ".profile.txt")
        if not profile.is_file() or profile.stat().st_mtime < field.stat().st_mtime:
            subprocess.run([str(executable), str(field), str(profile)], check=True)
        paths.append(profile)
    return paths


def existing_profiles(stem: str) -> list[Path]:
    return [
        find_time_file(EXISTING_PROFILE_DIRECTORY, f"*{stem}*.profile.txt", float(time))
        for time in PROFILE_TIMES
    ]


def load_case_tables(spec: dict[str, object]) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    statistics = []
    histories = []
    for name in spec["directories"]:
        directory = RESULTS / str(name)
        statistics.append(read_table(find_single(directory, "*.statistics.r64.txt")))
        histories.append(read_table(find_single(directory, "*.history.r64.txt")))
    return concatenate_tables(statistics), concatenate_tables(histories)


def make_profile(
    paths: list[Path], statistics: dict[str, np.ndarray]
) -> dict[str, np.ndarray | float]:
    pooled = pooled_profiles(paths)
    folded = fold_channel(pooled)
    sample_times = np.asarray(PROFILE_TIMES, dtype=float)
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


def parse_manifest(directory_name: str) -> dict[str, object]:
    directory = RESULTS / directory_name
    path = find_single(directory, "*.manifest.r64.txt")
    text = path.read_text(encoding="utf-8")

    def first(pattern: str, default: str = "unknown") -> str:
        match = re.search(pattern, text, flags=re.MULTILINE)
        return match.group(1) if match else default

    history = read_table(find_single(directory, "*.history.r64.txt"))
    files = [path for path in directory.iterdir() if path.is_file()]
    return {
        "directory": directory_name,
        "program_version": first(r"^program_version=(.+)$", "unknown"),
        "git_commit": first(r"^git_commit=(.+)$", "unknown"),
        "mpi_ranks": int(first(r"^mpi_ranks=(\d+)$", "0")),
        "riemann_solver": first(r"riemann_solver=([^;]+)"),
        "mdcd_dissipation": float(first(r"mdcd_dissipation=([^;]+)", "nan")),
        "cfl": float(np.nanmedian(history["cfl"])),
        "restart_compatibility": first(r"restart\.compatibility=([^,;\)]+)", "strict/default"),
        "stop_reason": first(r"^stop_reason=(.+)$", "unknown"),
        "field_count": len(list(directory.glob("*.field.*.cgns"))),
        "checkpoint_count": len(list(directory.glob("*.checkpoint.step*.cgns"))),
        "bytes": sum(item.stat().st_size for item in files),
    }


def common_run_summary(
    history: dict[str, np.ndarray], statistics: dict[str, np.ndarray], end: float = 300.0
) -> dict[str, object]:
    mask_h = (history["time"] >= history["time"][0] - 1.0e-10) & (
        history["time"] <= end + 1.0e-10
    )
    mask_s = (statistics["time"] >= statistics["time"][0] - 1.0e-10) & (
        statistics["time"] <= end + 1.0e-10
    )
    positive_dt = history["dt"][mask_h & (history["dt"] > 0.0)]
    start_index = int(np.flatnonzero(mask_s)[0])
    end_index = int(np.flatnonzero(mask_s)[-1])
    start_step = float(history["step"][np.flatnonzero(mask_h)[0]])
    end_step = float(history["step"][np.flatnonzero(mask_h)[-1]])

    # Wall time resets at segment boundaries.  Sum each monotone block.
    selected_wall = history["wall_time"][mask_h]
    selected_time = history["time"][mask_h]
    breaks = np.flatnonzero(np.diff(selected_wall) < 0.0) + 1
    blocks = np.split(np.arange(selected_wall.size), breaks)
    wall_seconds = float(
        sum(selected_wall[block[-1]] - selected_wall[block[0]] for block in blocks if block.size)
    )
    advanced_steps = int(round(end_step - start_step))
    first_mass = float(statistics["total_mass"][start_index])
    last_mass = float(statistics["total_mass"][end_index])
    first_energy = float(statistics["total_energy"][start_index])
    last_energy = float(statistics["total_energy"][end_index])
    return {
        "start_time": float(statistics["time"][start_index]),
        "end_time": float(statistics["time"][end_index]),
        "advanced_steps": advanced_steps,
        "wall_seconds": wall_seconds,
        "wall_hours": wall_seconds / 3600.0,
        "steps_per_wall_second": advanced_steps / wall_seconds,
        "simulated_time_per_wall_hour": (
            float(statistics["time"][end_index] - statistics["time"][start_index])
            / (wall_seconds / 3600.0)
        ),
        "dt_min_positive": float(np.min(positive_dt)),
        "dt_median_positive": float(np.median(positive_dt)),
        "dt_max": float(np.max(positive_dt)),
        "mass_relative_drift": (last_mass - first_mass) / first_mass,
        "energy_relative_change": (last_energy - first_energy) / first_energy,
        "reconstruction_fallbacks_max": int(np.nanmax(history["reconstruction_fallbacks"][mask_h])),
        "riemann_fallbacks_max": int(np.nanmax(history["riemann_fallbacks"][mask_h])),
        "troubled_cells_max": int(np.nanmax(history["troubled_cells"][mask_h])),
        "local_recomputations_max": int(np.nanmax(history["local_recomputations"][mask_h])),
        "step_retries_max": int(np.nanmax(history["step_retries"][mask_h])),
    }


def restart_continuity(hllc_statistics: dict[str, np.ndarray]) -> dict[str, object]:
    source_directory = RESULTS / "lowmach-longrun-segment02"
    source = read_table(find_single(source_directory, "*.statistics.r64.txt"))
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
        name: float(hllc_statistics[name][0] - source[name][-1]) for name in names
    }
    return {
        "source": "lowmach-longrun-segment02",
        "target": "branch-scmm6_mdcd0p001_hllc-t223",
        "differences": differences,
        "max_absolute_difference": max(abs(value) for value in differences.values()),
    }


def profile_to_json(profile: dict[str, np.ndarray | float]) -> dict[str, object]:
    result: dict[str, object] = {}
    for name, value in profile.items():
        result[name] = value.tolist() if isinstance(value, np.ndarray) else float(value)
    return result


def metric_row(metric: dict[str, object]) -> list[float]:
    errors = metric["inner_scaled_profile_errors_on_common_y_plus"]
    return [
        100.0 * float(metric["re_tau_relative_error"]),
        100.0 * float(metric["bulk_velocity_plus_relative_error"]),
        100.0 * float(metric["skin_friction_coefficient_relative_error"]),
        100.0 * float(errors["mean_u_plus"]["relative_l2"]),
        100.0 * float(errors["rms_u_plus"]["relative_l2"]),
        100.0 * float(errors["rms_v_plus"]["relative_l2"]),
        100.0 * float(errors["rms_w_plus"]["relative_l2"]),
        100.0 * float(errors["minus_uv_plus"]["relative_l2"]),
        100.0 * float(metric["outer_scaled_mean_velocity_error"]["relative_l2"]),
    ]


def mean_profile_error(metric: dict[str, object]) -> float:
    errors = metric["inner_scaled_profile_errors_on_common_y_plus"]
    return float(
        np.mean(
            [
                errors[name]["relative_l2"]
                for name in (
                    "mean_u_plus",
                    "rms_u_plus",
                    "rms_v_plus",
                    "rms_w_plus",
                    "minus_uv_plus",
                )
            ]
        )
    )


def pairwise_profile_difference(
    first: dict[str, np.ndarray | float], second: dict[str, np.ndarray | float]
) -> dict[str, float]:
    names = ("mean_u_plus", "rms_u_plus", "rms_v_plus", "rms_w_plus", "minus_uv_plus")
    result = {}
    for name in names:
        first_y = np.asarray(first["y_plus"])
        second_y = np.asarray(second["y_plus"])
        limit = min(float(first_y[-1]), float(second_y[-1]))
        mask = first_y <= limit + 1.0e-12
        interpolated = np.interp(first_y[mask], second_y, np.asarray(second[name]))
        values = np.asarray(first[name])[mask]
        result[f"{name}_relative_l2"] = float(
            np.linalg.norm(values - interpolated) / np.linalg.norm(interpolated)
        )
    return result


def plot_histories(statistics: dict[str, dict[str, np.ndarray]], dns: dict[str, object]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), sharex=True)
    for label, spec in CASES.items():
        table = statistics[label]
        mask = (table["time"] >= 223.0) & (table["time"] <= 300.0 + 1.0e-10)
        time = table["time"][mask]
        bulk = table["total_momentum_x"][mask] / table["total_mass"][mask]
        friction = table["channel_friction_velocity"][mask]
        asymmetry = (
            np.abs(table["channel_wall_shear_lower"][mask] - table["channel_wall_shear_upper"][mask])
            / table["channel_wall_shear_mean"][mask]
        )
        cf = 2.0 * (friction / bulk) ** 2
        values = (bulk / friction, table["channel_re_tau"][mask], 100.0 * asymmetry, cf)
        for axis, value in zip(axes.flat, values):
            axis.plot(
                time,
                value,
                color=spec["color"],
                linestyle=spec["linestyle"],
                linewidth=1.3,
                label=spec["short"],
            )
    axes[0, 0].axhline(float(dns["bulk_velocity_plus"]), color="black", linestyle=":", linewidth=1.1)
    axes[0, 1].axhline(DNS_RE_TAU, color="black", linestyle=":", linewidth=1.1)
    axes[1, 1].axhline(float(dns["skin_friction_coefficient"]), color="black", linestyle=":", linewidth=1.1)
    axes[0, 0].set_ylabel(r"$U_b^+$")
    axes[0, 1].set_ylabel(r"$Re_\tau$")
    axes[1, 0].set_ylabel("wall-shear asymmetry [%]")
    axes[1, 1].set_ylabel(r"$C_f$")
    axes[1, 0].set_xlabel("physical time")
    axes[1, 1].set_xlabel("physical time")
    for axis in axes.flat:
        axis.grid(True, alpha=0.25)
    axes[0, 0].legend(ncol=1, fontsize=8.5)
    fig.suptitle("Common-time evolution after the t=223.0031 branch point")
    fig.tight_layout()
    fig.savefig(OUTPUT / "riemann_common_history.png", dpi=220)
    plt.close(fig)


def plot_profiles(
    dns: dict[str, object], profiles: OrderedDict[str, dict[str, np.ndarray | float]]
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    axes[0].semilogx(dns["y_plus"], dns["mean_u_plus"], "ko", markersize=3.0, label="DNS")
    axes[1].plot(dns["eta"], np.asarray(dns["mean_u_plus"]) / float(dns["bulk_velocity_plus"]), "ko", markersize=3.0, label="DNS")
    for label, profile in profiles.items():
        spec = CASES[label]
        axes[0].semilogx(profile["y_plus"], profile["mean_u_plus"], color=spec["color"], linestyle=spec["linestyle"], linewidth=1.7, label=spec["short"])
        axes[1].plot(profile["eta"], profile["mean_u_over_bulk"], color=spec["color"], linestyle=spec["linestyle"], linewidth=1.7, label=spec["short"])
    axes[0].set_xlabel(r"$y^+$")
    axes[0].set_ylabel(r"$U^+$")
    axes[1].set_xlabel(r"$y/h$")
    axes[1].set_ylabel(r"$U/U_b$")
    axes[0].set_xlim(1.0, 220.0)
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.legend(fontsize=8.5)
    fig.suptitle("Common t=275-300 mean-velocity window vs ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_common_275_300_mean_velocity.png", dpi=220)
    plt.close(fig)

    quantities = (
        ("rms_u_plus", r"$u'_{rms}{}^+$"),
        ("rms_v_plus", r"$v'_{rms}{}^+$"),
        ("rms_w_plus", r"$w'_{rms}{}^+$"),
        ("minus_uv_plus", r"$-\overline{u'v'}{}^+$"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0))
    for axis, (name, ylabel) in zip(axes.flat, quantities):
        axis.semilogx(dns["y_plus"], dns[name], "ko", markersize=3.0, label="DNS")
        for label, profile in profiles.items():
            spec = CASES[label]
            axis.semilogx(profile["y_plus"], profile[name], color=spec["color"], linestyle=spec["linestyle"], linewidth=1.6, label=spec["short"])
        axis.set_xlim(1.0, 220.0)
        axis.set_xlabel(r"$y^+$")
        axis.set_ylabel(ylabel)
        axis.grid(True, alpha=0.25)
    axes[0, 0].legend(fontsize=8.3)
    fig.suptitle("Common t=275-300 turbulence statistics vs ERCOFTAC Case 032 DNS")
    fig.tight_layout()
    fig.savefig(OUTPUT / "dns_common_275_300_turbulence_statistics.png", dpi=220)
    plt.close(fig)


def plot_error_summary(
    metrics: OrderedDict[str, dict[str, object]], filename: str, title: str
) -> None:
    names = list(metrics)
    columns = (r"$Re_\tau$", r"$U_b^+$", r"$C_f$", r"$U^+$", r"$u'_{rms}{}^+$", r"$v'_{rms}{}^+$", r"$w'_{rms}{}^+$", r"$-uv^+$", r"$U/U_b$")
    values = np.asarray([metric_row(metrics[name]) for name in names])
    fig, axis = plt.subplots(figsize=(13.0, 5.4))
    x = np.arange(len(columns))
    width = 0.24 if len(names) <= 3 else 0.19
    offsets = (np.arange(len(names)) - 0.5 * (len(names) - 1)) * width
    def color(name: str) -> str:
        if name in CASES:
            return str(CASES[name]["color"])
        if name.startswith("HLLC"):
            return "#009e73"
        if name.startswith("baseline"):
            return "#777777"
        if name.startswith("Roe d=0.001"):
            return "#d55e00"
        if name.startswith("Roe d=0.01"):
            return "#0072b2"
        return "#777777"

    palette = [color(name) for name in names]
    for index, name in enumerate(names):
        axis.bar(x + offsets[index], values[index], width, color=palette[index], label=name)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xticks(x, columns)
    axis.set_ylabel("error relative to DNS [%]")
    axis.grid(True, axis="y", alpha=0.25)
    axis.legend(fontsize=8.0, ncol=2)
    axis.set_title(title)
    fig.tight_layout()
    fig.savefig(OUTPUT / filename, dpi=220)
    plt.close(fig)


def write_markdown(report: dict[str, object]) -> None:
    dns = report["reference"]
    common = report["common_window"]
    final = report["final_20_time_units"]
    performance = report["common_run_performance"]
    metrics = report["dns_common_window_metrics"]
    latest = report["dns_latest_available_metrics"]
    hllc_run = performance["HLLC, d=0.001, CFL=0.6"]
    roe_low_run = performance["Roe, d=0.001, CFL=0.3"]
    physical_throughput_ratio = (
        hllc_run["simulated_time_per_wall_hour"]
        / roe_low_run["simulated_time_per_wall_hour"]
    )
    step_throughput_ratio = (
        hllc_run["steps_per_wall_second"] / roe_low_run["steps_per_wall_second"]
    )
    lines = [
        "# Case05 HLLC 分支分析及 Roe/DNS 对比",
        "",
        "## 结论摘要",
        "",
        "本次新增分支从与两条 Roe 分支完全相同的 `t=223.003107849` 状态启动。"
        "restart 接口所比较状态量的最大绝对差为 "
        f"`{report['restart_continuity']['max_absolute_difference']:.3e}`，因此分叉初值一致。",
        "",
        f"严格公共窗口采用 `t={PROFILE_TIMES[0]}-{PROFILE_TIMES[-1]}` 的六个三维场。"
        "HLLC 同时使用了更高的 CFL=0.6，而 Roe 为 CFL=0.3，所以本次数据反映的是"
        "“Riemann 求解器 + 时间步设置”的联合差异，不能完全归因于 HLLC/Roe 本身。",
        "",
    ]
    ranking = sorted(metrics, key=lambda name: mean_profile_error(metrics[name]))
    lines.extend(
        [
            f"按五个内尺度剖面相对 L2 误差的算术平均，公共窗口排序为："
            + " < ".join(
                f"{CASES[name]['short']} ({100.0 * mean_profile_error(metrics[name]):.2f}%)"
                for name in ranking
            )
            + "。",
            "",
            "## 运行配置与计算量",
            "",
            "| 算例 | Riemann | d | CFL | 步数至 t=300 | 墙钟小时 | 物理时长/墙钟小时 | 中位 dt |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for name, config in report["configurations"].items():
        run = performance[name]
        lines.append(
            f"| {CASES[name]['short']} | {config['riemann_solver']} | {config['mdcd_dissipation']:.3g} | "
            f"{config['cfl']:.1f} | {run['advanced_steps']:,} | {run['wall_hours']:.3f} | "
            f"{run['simulated_time_per_wall_hour']:.3f} | {run['dt_median_positive']:.6g} |"
        )
    hllc_config = report["configurations"]["HLLC, d=0.001, CFL=0.6"]
    lines.extend(
        [
            "",
            f"HLLC 新增段包含 {hllc_config['field_count']} 个场文件、"
            f"{hllc_config['checkpoint_count']} 个编号 checkpoint，占用 "
            f"{hllc_config['bytes'] / 2**20:.1f} MiB；停止原因为 `{hllc_config['stop_reason']}`。",
            "",
            f"在本批实际墙钟记录中，HLLC 的物理时间推进率是同耗散 Roe 的 "
            f"{physical_throughput_ratio:.2f} 倍，每秒步数为 {step_throughput_ratio:.2f} 倍。"
            "这不是严格的核函数基准测试：CFL 同时从 0.3 改为 0.6，而且两次作业的节点负载可能不同。",
            "",
            "三种分支在公共时段内均未触发 reconstruction/Riemann fallback、troubled cell、"
            "local recomputation 或 step retry。HLLC 的质量相对漂移为 "
            f"`{performance['HLLC, d=0.001, CFL=0.6']['mass_relative_drift']:.3e}`。",
            "",
            "## 公共窗口积分量与稳态性",
            "",
            f"下表为最后 20 个时间单位 `t={FINAL_SCALAR_WINDOW[0]:g}-{FINAL_SCALAR_WINDOW[1]:g}`；"
            "斜率仍显著非零表示有限时间统计尚未严格收敛。",
            "",
            "| 算例 | Ub | dUb/dt | Re_tau | dRe_tau/dt | 壁面剪切不对称均值 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for name in CASES:
        item = final[name]
        lines.append(
            f"| {CASES[name]['short']} | {item['bulk_velocity']['mean']:.6f} | "
            f"{item['bulk_velocity']['slope_per_time']:+.3e} | {item['re_tau']['mean']:.3f} | "
            f"{item['re_tau']['slope_per_time']:+.3e} | "
            f"{100.0 * item['wall_shear_relative_asymmetry']['mean']:.2f}% |"
        )
    lines.extend(
        [
            "",
            "HLLC 在该时间段内的 Re_tau 正在明显回落，但平均值仍高于 DNS 目标 180；"
            "因此 `t=275-300` 应视为当前有限时窗，而非 HLLC 的最终稳态统计。",
            "",
            "## 与 ERCOFTAC Case 032 DNS 的严格公共窗口对比",
            "",
            f"DNS 参考量：`Re_tau={dns['re_tau']:.0f}`，`Ub+={dns['bulk_velocity_plus']:.6f}`，"
            f"`Cf={dns['skin_friction_coefficient']:.8f}`。百分数为相对 DNS 的有符号误差；"
            "剖面列为公共 y+ 区间的相对 L2 误差。",
            "",
            "| 算例 | Re_tau | Ub+ | Cf | U+ L2 | urms+ L2 | vrms+ L2 | wrms+ L2 | -uv+ L2 | U/Ub L2 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for name in CASES:
        row = metric_row(metrics[name])
        lines.append(
            f"| {CASES[name]['short']} | " + " | ".join(f"{value:+.2f}%" if index < 3 else f"{value:.2f}%" for index, value in enumerate(row)) + " |"
        )
    pair = report["hllc_vs_roe_d001_profile_difference"]
    lines.extend(
        [
            "",
            "在相同 d=0.001 的两条分支之间，HLLC 相对 Roe 的公共窗口剖面差异为："
            f"U+ {100.0 * pair['mean_u_plus_relative_l2']:.2f}%、"
            f"urms+ {100.0 * pair['rms_u_plus_relative_l2']:.2f}%、"
            f"vrms+ {100.0 * pair['rms_v_plus_relative_l2']:.2f}%、"
            f"wrms+ {100.0 * pair['rms_w_plus_relative_l2']:.2f}%、"
            f"-uv+ {100.0 * pair['minus_uv_plus_relative_l2']:.2f}%。",
            "",
            "HLLC 与同耗散 Roe 的总体结果非常接近：HLLC 的 U+、urms+、vrms+ 和外尺度"
            "平均速度误差略小，但 wrms+ 和 -uv+ 误差更大。Roe d=0.01 虽然 Re_tau 最接近"
            "目标值，其速度及雷诺应力剖面误差显著更大，且最后 20 时间单位的壁面剪切"
            "不对称均值达到 21.76%。",
            "",
            "## 与上一份报告的最新可用结果衔接",
            "",
            "下表的时间窗不同，只用于说明当前各算例最新状态，不能作为严格的同时间算法排序。",
            "",
            "| 结果窗口 | Re_tau误差 | Ub+误差 | Cf误差 | U+ L2 | 五剖面平均L2 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for name, item in latest.items():
        row = metric_row(item)
        lines.append(
            f"| {name} | {row[0]:+.2f}% | {row[1]:+.2f}% | {row[2]:+.2f}% | "
            f"{row[3]:.2f}% | {100.0 * mean_profile_error(item):.2f}% |"
        )
    lines.extend(
        [
            "",
            "## 图件",
            "",
            "- `riemann_common_history.png`：三个分支从共同分叉点到 t=300 的积分量历程。",
            "- `dns_common_275_300_mean_velocity.png`：公共窗口平均速度。",
            "- `dns_common_275_300_turbulence_statistics.png`：公共窗口雷诺应力统计。",
            "- `dns_common_275_300_error_summary.png`：严格公共窗口误差汇总。",
            "- `dns_latest_available_error_summary.png`：与上一份最终报告最新窗口的衔接比较。",
            "",
            "## 判断",
            "",
            "1. HLLC 分支数值运行健康，restart 连续且在 CFL=0.6 下未触发任何鲁棒性回退。",
            "2. HLLC 约用 Roe 一半的步数推进到 t=300；实测物理时间推进率更高，但性能结论同时包含 CFL 与运行环境差异。",
            "3. HLLC 当前窗口尚处于明显演化阶段，尤其 Re_tau 高于 180 且仍有趋势，不能认定稳态。",
            "4. 对算法本身作最终排序前，建议把 HLLC 至少续算到趋势减弱，并使用与 Roe 相同 CFL 做一段控制试验。",
        ]
    )
    (OUTPUT / "HLLC_ROE_DNS_ANALYSIS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    dns = load_dns()
    statistics: dict[str, dict[str, np.ndarray]] = {}
    histories: dict[str, dict[str, np.ndarray]] = {}
    configurations: dict[str, dict[str, object]] = {}
    profiles: OrderedDict[str, dict[str, np.ndarray | float]] = OrderedDict()
    qualities: dict[str, dict[str, float]] = {}

    hllc_paths = ensure_hllc_profiles()
    for label, spec in CASES.items():
        statistics[label], histories[label] = load_case_tables(spec)
        configurations[label] = parse_manifest(spec["directories"][0])
        paths = hllc_paths if spec["hllc"] else existing_profiles(str(spec["profile_stem"]))
        profiles[label] = make_profile(paths, statistics[label])
        qualities[label] = profile_summary(pooled_profiles(paths))

    common_metrics: OrderedDict[str, dict[str, object]] = OrderedDict(
        (label, comparison_metrics(profile, dns)) for label, profile in profiles.items()
    )
    common_scalars = {
        label: scalar_summary(statistics[label], float(PROFILE_TIMES[0]), float(PROFILE_TIMES[-1]))
        for label in CASES
    }
    final_scalars = {
        label: scalar_summary(statistics[label], *FINAL_SCALAR_WINDOW) for label in CASES
    }
    performance = {
        label: common_run_summary(histories[label], statistics[label]) for label in CASES
    }

    previous = json.loads(FINAL_METRICS.read_text(encoding="utf-8"))
    latest_metrics: OrderedDict[str, dict[str, object]] = OrderedDict(
        [
            ("HLLC d=0.001, t=275-300", common_metrics["HLLC, d=0.001, CFL=0.6"]),
            ("baseline Roe, t=195-220", previous["dns_windows"]["baseline 195-220"]["metrics"]),
            ("Roe d=0.001, t=410-435", previous["dns_windows"]["d=0.001 410-435"]["metrics"]),
            ("Roe d=0.01, t=300-325", previous["dns_windows"]["d=0.01 300-325"]["metrics"]),
        ]
    )

    report: dict[str, object] = {
        "reference": {
            "name": "ERCOFTAC Case 032 simul1.dat / Kim-Moin-Moser Re_tau=180",
            "re_tau": DNS_RE_TAU,
            "bulk_velocity_plus": float(dns["bulk_velocity_plus"]),
            "skin_friction_coefficient": float(dns["skin_friction_coefficient"]),
        },
        "profile_times": list(PROFILE_TIMES),
        "final_scalar_window": list(FINAL_SCALAR_WINDOW),
        "restart_continuity": restart_continuity(statistics["HLLC, d=0.001, CFL=0.6"]),
        "configurations": configurations,
        "common_run_performance": performance,
        "common_window": common_scalars,
        "final_20_time_units": final_scalars,
        "profile_quality": qualities,
        "dns_common_window_metrics": common_metrics,
        "dns_latest_available_metrics": latest_metrics,
        "hllc_vs_roe_d001_profile_difference": pairwise_profile_difference(
            profiles["HLLC, d=0.001, CFL=0.6"],
            profiles["Roe, d=0.001, CFL=0.3"],
        ),
        "profiles": {label: profile_to_json(profile) for label, profile in profiles.items()},
    }

    plot_histories(statistics, dns)
    plot_profiles(dns, profiles)
    plot_error_summary(
        common_metrics,
        "dns_common_275_300_error_summary.png",
        "Strict common-window errors: t=275-300",
    )
    plot_error_summary(
        latest_metrics,
        "dns_latest_available_error_summary.png",
        "Latest available windows (different sampling times)",
    )
    (OUTPUT / "hllc_roe_dns_metrics.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_markdown(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

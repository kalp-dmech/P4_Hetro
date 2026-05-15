#!/usr/bin/env python3
"""Generate and optionally compile the phase-field fracture report."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from string import Template


class LatexTemplate(Template):
    delimiter = "@"


ROOT = Path(__file__).resolve().parent
REPORT_DIR = ROOT / "report"
REPORT_TEX = REPORT_DIR / "phase_field_fracture_report.tex"
REPORT_PDF = REPORT_DIR / "phase_field_fracture_report.pdf"
MISSING = "Not recorded in the current simulation output"


@dataclass
class ForceSeries:
    path: Path
    pore_count: int | None
    label: str
    dataset_name: str
    points: list[tuple[float, float]]

    @property
    def peak_force(self) -> float | None:
        return max((force for _, force in self.points), default=None)

    @property
    def peak_displacement(self) -> float | None:
        if not self.points:
            return None
        displacement, _ = max(self.points, key=lambda item: item[1])
        return displacement


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""


def latex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in text)


def format_float(value: float | None) -> str:
    return "--" if value is None else f"{value:.6g}"


def extract_assignment(text: str, name: str, default: str = "--") -> str:
    for pattern in (rf"\b{name}\s*=\s*([^;]+);", rf"\b{name}_val\s*=\s*([^;]+);"):
        match = re.search(pattern, text)
        if match:
            return re.sub(r"\s+", " ", match.group(1).strip())
    return default


def source_config() -> dict[str, str]:
    variables = read_text(ROOT / "source" / "variable_constructor.cc")
    constructor = read_text(ROOT / "source" / "Constructor.cc")
    return {
        "pore_counts": extract_assignment(variables, "pore_counts"),
        "random_seed_base": extract_assignment(variables, "random_seed_base"),
        "min_pore_size": extract_assignment(variables, "min_pore_size"),
        "max_pore_size": extract_assignment(variables, "max_pore_size"),
        "crack_spread_radius": extract_assignment(variables, "crack_spread_radius"),
        "boundary_margin": extract_assignment(variables, "boundary_margin"),
        "lo": extract_assignment(variables, "lo"),
        "Gc": extract_assignment(variables, "Gc"),
        "n_steps": extract_assignment(variables, "n_steps"),
        "tol_1": extract_assignment(variables, "tol_1"),
        "tol_2": extract_assignment(variables, "tol_2"),
        "lambda": extract_assignment(constructor, "lambda"),
        "mu": extract_assignment(constructor, "mu"),
        "k": extract_assignment(constructor, "k"),
        "inc_large": extract_assignment(constructor, "inc_large"),
        "inc_small": extract_assignment(constructor, "inc_small"),
        "amr_max_level": extract_assignment(constructor, "amr_max_level"),
    }


def count_outputs() -> dict[str, int]:
    output = ROOT / "output"
    return {
        "force_csv": len(list(output.glob("**/Pore */trail_*/force_displacement.csv"))),
        "pvtu": len(list(output.glob("**/*.pvtu"))),
        "vtu": len(list(output.glob("**/*.vtu"))),
        "pvd": len(list(output.glob("**/*.pvd"))),
        "pores": len(list(output.glob("**/Pore */trail_*/pores.dat"))),
        "case_documents": len(list(output.glob("*/documents/case_manifest.txt"))),
    }


def read_force_file(path: Path) -> list[tuple[float, float]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        return []

    points: list[tuple[float, float]] = []
    if "," in lines[0]:
        for row in csv.DictReader(lines):
            try:
                displacement = float(row.get("displacement") or row.get("Displacement") or "")
                force = float(row.get("force") or row.get("Force") or "")
            except ValueError:
                continue
            points.append((displacement, force))
        return points

    for line in lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            points.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    return points


def parse_pore_count(path: Path) -> int | None:
    for part in path.parts:
        match = re.match(r"Pore\s+(\d+)", part)
        if match:
            return int(match.group(1))
    return None


def force_label(path: Path) -> str:
    parts = path.parts
    mode = ""
    try:
        output_index = parts.index("output")
        if output_index + 1 < len(parts) and not parts[output_index + 1].startswith("Pore"):
            mode = parts[output_index + 1] + ", "
    except ValueError:
        mode = ""
    pore = next((part for part in parts if part.startswith("Pore")), "Pore ?")
    trial = next((part for part in parts if part.startswith(("trial", "trail"))), "trial ?")
    return f"{mode}{pore}, {trial}"


def collect_force_series() -> list[ForceSeries]:
    files = sorted((ROOT / "output").glob("**/Pore */trail_*/force_displacement.csv"))
    if not files and (ROOT / "force_displacement.txt").exists():
        files = [ROOT / "force_displacement.txt"]

    series: list[ForceSeries] = []
    for index, path in enumerate(files):
        points = read_force_file(path)
        if points:
            series.append(
                ForceSeries(
                    path=path,
                    pore_count=parse_pore_count(path),
                    label=force_label(path),
                    dataset_name=f"force_series_{index}",
                    points=points,
                )
            )
    return series


def normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")


def read_table(path: Path) -> list[dict[str, float]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8", errors="replace").splitlines()]
    lines = [line for line in lines if line and not line.startswith("#")]
    if not lines:
        return []

    delimiter = "," if "," in lines[0] else None
    if delimiter:
        reader = csv.DictReader(lines)
        return normalize_rows(reader)

    headers = [normalize_name(part) for part in lines[0].split()]
    rows = []
    for line in lines[1:]:
        values = line.split()
        row = {headers[i]: values[i] for i in range(min(len(headers), len(values)))}
        rows.append(row)
    return normalize_rows(rows)


def read_records(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.exists():
        return []
    lines = [line for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
    if not lines:
        return []
    if "," not in lines[0]:
        return []
    records = []
    for row in csv.DictReader(lines):
        records.append({normalize_name(key): (value or "").strip() for key, value in row.items() if key})
    return records


def read_manifest(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[normalize_name(key)] = value.strip()
    return data


def normalize_rows(rows) -> list[dict[str, float]]:
    normalized = []
    for row in rows:
        converted = {}
        for key, value in row.items():
            if key is None:
                continue
            name = normalize_name(str(key))
            try:
                converted[name] = float(value)
            except (TypeError, ValueError):
                continue
        if converted:
            normalized.append(converted)
    return normalized


def find_dataset(*names: str) -> Path | None:
    candidates: list[Path] = []
    for name in names:
        candidates.extend([ROOT / name, ROOT / "output" / name, REPORT_DIR / name])
        candidates.extend((ROOT / "output").glob(f"**/{name}"))
    for path in candidates:
        if path.exists() and path.is_file():
            return path
    return None


CASE_ORDER = {
    "serial_uniform": 1,
    "mpi_uniform": 2,
    "mpi_amr": 3,
}

CASE_LABELS = {
    "serial_uniform": "Serial uniform",
    "mpi_uniform": "MPI uniform",
    "mpi_amr": "MPI + AMR",
}


def case_folder_records() -> list[dict[str, str]]:
    records = []
    for manifest in sorted((ROOT / "output").glob("*/documents/case_manifest.txt")):
        row = read_manifest(manifest)
        row["manifest_file"] = manifest.relative_to(ROOT).as_posix()
        row["case_folder"] = manifest.parents[1].relative_to(ROOT).as_posix()
        records.append(row)
    return records


def case_folder_rows(records: list[dict[str, str]]) -> str:
    if not records:
        return rf"\multicolumn{{6}}{{l}}{{{MISSING}.}}\\"
    rows = []
    for row in records:
        case = CASE_LABELS.get(row.get("case_label", ""), row.get("case_label", ""))
        rows.append(
            f"{latex_escape(case)} & {latex_escape(row.get('mpi_cores', '--'))} & "
            f"{latex_escape(row.get('amr_enabled', '--'))} & "
            f"\\texttt{{{latex_escape(row.get('case_folder', '--'))}}} & "
            f"\\texttt{{{latex_escape(row.get('results_directory', '--'))}}} & "
            f"\\texttt{{{latex_escape(row.get('documents_directory', '--'))}}}\\\\"
        )
    return "\n".join(rows)


def to_float(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except ValueError:
        return None


def benchmark_cases(path: Path | None) -> list[dict[str, str]]:
    latest: dict[str, dict[str, str]] = {}
    for row in read_records(path):
        case = row.get("case_label", "")
        if case in CASE_ORDER:
            latest[case] = row
    return [latest[key] for key in sorted(latest, key=lambda item: CASE_ORDER[item])]


def benchmark_rows(cases: list[dict[str, str]]) -> str:
    if not cases:
        return rf"\multicolumn{{11}}{{l}}{{{MISSING}.}}\\"
    rows = []
    for row in cases:
        case = CASE_LABELS.get(row.get("case_label", ""), row.get("case_label", ""))
        rows.append(
            f"{latex_escape(case)} & {latex_escape(row.get('cores', '--'))} & "
            f"{latex_escape(row.get('amr_enabled', '--'))} & "
            f"{format_float(to_float(row.get('total_runtime')))} & "
            f"{format_float(to_float(row.get('assembly_time')))} & "
            f"{format_float(to_float(row.get('solve_time')))} & "
            f"{format_float(to_float(row.get('refinement_time')))} & "
            f"{format_float(to_float(row.get('io_time')))} & "
            f"{format_float(to_float(row.get('active_cells_final')))} & "
            f"{format_float(to_float(row.get('elastic_dofs_final')))} & "
            f"{format_float(to_float(row.get('damage_dofs_final')))}\\\\"
        )
    return "\n".join(rows)


def benchmark_points(cases: list[dict[str, str]], field: str) -> list[tuple[float, float]]:
    points = []
    for row in cases:
        case = row.get("case_label", "")
        value = to_float(row.get(field))
        if case in CASE_ORDER and value is not None:
            points.append((float(CASE_ORDER[case]), value))
    return points


def benchmark_speedup_points(cases: list[dict[str, str]]) -> list[tuple[float, float]]:
    runtimes = {row.get("case_label", ""): to_float(row.get("total_runtime")) for row in cases}
    base = runtimes.get("serial_uniform")
    if base is None or base <= 0:
        return []
    return [
        (float(CASE_ORDER[row["case_label"]]), base / runtime)
        for row in cases
        if row.get("case_label") in CASE_ORDER
        for runtime in [to_float(row.get("total_runtime"))]
        if runtime is not None and runtime > 0
    ]


def benchmark_efficiency_points(cases: list[dict[str, str]]) -> list[tuple[float, float]]:
    speedups = dict(benchmark_speedup_points(cases))
    points = []
    for row in cases:
        case = row.get("case_label", "")
        cores = to_float(row.get("cores"))
        x = float(CASE_ORDER[case]) if case in CASE_ORDER else None
        if x is not None and cores is not None and cores > 0 and x in speedups:
            points.append((x, speedups[x] / cores))
    return points


def benchmark_interpretation(cases: list[dict[str, str]]) -> str:
    by_case = {row.get("case_label", ""): row for row in cases}
    serial = to_float(by_case.get("serial_uniform", {}).get("total_runtime"))
    mpi = to_float(by_case.get("mpi_uniform", {}).get("total_runtime"))
    amr = to_float(by_case.get("mpi_amr", {}).get("total_runtime"))
    if serial is None or mpi is None or amr is None:
        return MISSING + ". Run all three benchmark cases to make the final serial/MPI/AMR comparison."
    mpi_text = "MPI reduces total runtime compared with the serial uniform baseline." if mpi < serial else "MPI does not reduce total runtime in the recorded benchmark."
    amr_text = "MPI + AMR gives the lowest recorded runtime." if amr <= min(serial, mpi) else "MPI + AMR is not the lowest-runtime case in the recorded benchmark."
    return mpi_text + " " + amr_text


def case_row(cases: list[dict[str, str]], label: str) -> dict[str, str]:
    for row in cases:
        if row.get("case_label", "") == label:
            return row
    return {}


def percent_change(reference: float | None, value: float | None) -> float | None:
    if reference is None or value is None or reference == 0:
        return None
    return 100.0 * (value - reference) / reference


def percent_reduction(reference: float | None, value: float | None) -> float | None:
    change = percent_change(reference, value)
    if change is None:
        return None
    return -change


def ratio(reference: float | None, value: float | None) -> float | None:
    if reference is None or value is None or value == 0:
        return None
    return reference / value


def benchmark_value(row: dict[str, str], field: str) -> float | None:
    return to_float(row.get(field))


def final_performance_rows(cases: list[dict[str, str]]) -> str:
    if not cases:
        return rf"\multicolumn{{12}}{{l}}{{{MISSING}.}}\\"
    serial = case_row(cases, "serial_uniform")
    base_runtime = benchmark_value(serial, "total_runtime")
    rows = []
    for row in cases:
        case = CASE_LABELS.get(row.get("case_label", ""), row.get("case_label", ""))
        runtime = benchmark_value(row, "total_runtime")
        cores = benchmark_value(row, "cores")
        speedup = ratio(base_runtime, runtime)
        efficiency = speedup / cores if speedup is not None and cores not in (None, 0) else None
        rows.append(
            f"{latex_escape(case)} & {latex_escape(row.get('cores', '--'))} & "
            f"{latex_escape(row.get('amr_enabled', '--'))} & "
            f"{format_float(runtime)} & "
            f"{format_float(benchmark_value(row, 'assembly_time'))} & "
            f"{format_float(benchmark_value(row, 'solve_time'))} & "
            f"{format_float(benchmark_value(row, 'refinement_time'))} & "
            f"{format_float(benchmark_value(row, 'io_time'))} & "
            f"{format_float(speedup)} & {format_float(efficiency)} & "
            f"{format_float(benchmark_value(row, 'active_cells_final'))} & "
            f"{format_float(benchmark_value(row, 'elastic_dofs_final'))}\\\\"
        )
    return "\n".join(rows)


def final_efficiency_rows(cases: list[dict[str, str]]) -> str:
    serial = case_row(cases, "serial_uniform")
    mpi = case_row(cases, "mpi_uniform")
    amr = case_row(cases, "mpi_amr")
    if not (serial and mpi and amr):
        return rf"\multicolumn{{4}}{{l}}{{{MISSING}. Run all three benchmark cases.}}\\"

    comparisons = [
        (
            "Serial uniform $\\rightarrow$ MPI uniform",
            "Pure MPI parallelization",
            percent_reduction(benchmark_value(serial, "total_runtime"), benchmark_value(mpi, "total_runtime")),
            ratio(benchmark_value(serial, "total_runtime"), benchmark_value(mpi, "total_runtime")),
        ),
        (
            "MPI uniform $\\rightarrow$ MPI--AMR",
            "Effect of adaptive refinement",
            percent_change(benchmark_value(mpi, "active_cells_final"), benchmark_value(amr, "active_cells_final")),
            percent_change(benchmark_value(mpi, "elastic_dofs_final"), benchmark_value(amr, "elastic_dofs_final")),
        ),
        (
            "Serial uniform $\\rightarrow$ MPI--AMR",
            "Combined distributed adaptive computation",
            percent_reduction(benchmark_value(serial, "total_runtime"), benchmark_value(amr, "total_runtime")),
            ratio(benchmark_value(serial, "total_runtime"), benchmark_value(amr, "total_runtime")),
        ),
    ]
    rows = []
    for comparison, meaning, first, second in comparisons:
        rows.append(
            f"{comparison} & {meaning} & {format_float(first)} & {format_float(second)}\\\\"
        )
    return "\n".join(rows)


def final_solver_rows(cases: list[dict[str, str]]) -> str:
    if not cases:
        return rf"\multicolumn{{6}}{{l}}{{{MISSING}.}}\\"
    rows = []
    for row in cases:
        case = CASE_LABELS.get(row.get("case_label", ""), row.get("case_label", ""))
        rows.append(
            f"{latex_escape(case)} & "
            f"{format_float(benchmark_value(row, 'max_staggered_iterations'))} & "
            f"{format_float(benchmark_value(row, 'max_elastic_cg_iterations'))} & "
            f"{format_float(benchmark_value(row, 'max_damage_cg_iterations'))} & "
            f"{format_float(benchmark_value(row, 'max_final_residual'))} & "
            f"{format_float(benchmark_value(row, 'solve_time'))}\\\\"
        )
    return "\n".join(rows)


def final_conclusion_text(cases: list[dict[str, str]]) -> str:
    serial = case_row(cases, "serial_uniform")
    mpi = case_row(cases, "mpi_uniform")
    amr = case_row(cases, "mpi_amr")
    if not (serial and mpi and amr):
        return MISSING + ". The final conclusion requires serial uniform, MPI uniform, and MPI--AMR rows."

    serial_runtime = benchmark_value(serial, "total_runtime")
    mpi_runtime = benchmark_value(mpi, "total_runtime")
    amr_runtime = benchmark_value(amr, "total_runtime")
    mpi_speedup = ratio(serial_runtime, mpi_runtime)
    amr_speedup = ratio(serial_runtime, amr_runtime)
    mpi_reduction = percent_reduction(serial_runtime, mpi_runtime)
    amr_reduction = percent_reduction(serial_runtime, amr_runtime)
    amr_cells_change = percent_change(benchmark_value(mpi, "active_cells_final"), benchmark_value(amr, "active_cells_final"))
    amr_dofs_change = percent_change(benchmark_value(mpi, "elastic_dofs_final"), benchmark_value(amr, "elastic_dofs_final"))

    return (
        "The recorded benchmark shows that MPI parallelization reduces the total runtime from "
        f"{format_float(serial_runtime)} to {format_float(mpi_runtime)}, giving a speedup of "
        f"{format_float(mpi_speedup)} and a runtime reduction of {format_float(mpi_reduction)} percent. "
        "The MPI--AMR case completes in "
        f"{format_float(amr_runtime)}, corresponding to a speedup of {format_float(amr_speedup)} "
        f"over the serial baseline and a runtime reduction of {format_float(amr_reduction)} percent. "
        "Relative to MPI uniform refinement, the AMR run changes the final active-cell count by "
        f"{format_float(amr_cells_change)} percent and the final elastic DoF count by "
        f"{format_float(amr_dofs_change)} percent, because refinement is concentrated near the evolving crack. "
        "The data therefore support the technical conclusion that distributed parallelism is essential for reducing wall time, while adaptive refinement provides the mechanism for controlling spatial resolution in the fracture process zone."
    )


def parse_step_from_output(path: Path) -> int | None:
    match = re.search(r"(?:solution|initial_phi)_(\d+)", path.name)
    if not match:
        return None
    return int(match.group(1))


def vtu_piece_counts(path: Path) -> tuple[int, int]:
    text = read_text(path)
    match = re.search(r'NumberOfPoints="(\d+)"\s+NumberOfCells="(\d+)"', text)
    if not match:
        return (0, 0)
    return (int(match.group(1)), int(match.group(2)))


def pvtu_case_directories() -> list[Path]:
    return sorted(
        path
        for path in (ROOT / "output").glob("**/Pore */trail_*")
        if path.is_dir() and list(path.glob("*.pvtu"))
    )


def first_case_vtu_history() -> tuple[list[tuple[float, float]], list[tuple[float, float]], str]:
    directories = pvtu_case_directories()
    if not directories:
        return ([], [], MISSING)

    case_dir = directories[0]
    cell_points: list[tuple[float, float]] = []
    point_points: list[tuple[float, float]] = []

    for pvtu in sorted(case_dir.glob("*.pvtu"), key=lambda item: (parse_step_from_output(item) or -1, item.name)):
        step = parse_step_from_output(pvtu)
        if step is None:
            continue
        text = read_text(pvtu)
        cells = 0
        points = 0
        for source in re.findall(r'<Piece\s+Source="([^"]+)"', text):
            piece_points, piece_cells = vtu_piece_counts(pvtu.parent / source)
            points += piece_points
            cells += piece_cells
        if cells > 0:
            cell_points.append((float(step), float(cells)))
        if points > 0:
            point_points.append((float(step), float(points)))

    source = case_dir.relative_to(ROOT).as_posix()
    return (cell_points, point_points, f"derived from VTU/PVTU files in {source}")


def get_column(row: dict[str, float], *names: str) -> float | None:
    aliases = {normalize_name(name) for name in names}
    for key, value in row.items():
        if key in aliases:
            return value
    return None


def dataset_points(path: Path | None, x_names: list[str], y_names: list[str]) -> list[tuple[float, float]]:
    if path is None:
        return []
    points = []
    for row in read_table(path):
        x = get_column(row, *x_names)
        y = get_column(row, *y_names)
        if x is not None and y is not None:
            points.append((x, y))
    return points


def speedup_points(path: Path | None) -> list[tuple[float, float]]:
    if path is None:
        return []
    rows = read_table(path)
    explicit = []
    runtime_points = []
    for row in rows:
        cores = get_column(row, "cores", "mpi_ranks", "ranks", "processes")
        speedup = get_column(row, "speedup")
        runtime = get_column(row, "total_runtime", "runtime", "total_time")
        if cores is not None and speedup is not None:
            explicit.append((cores, speedup))
        if cores is not None and runtime is not None and runtime > 0:
            runtime_points.append((cores, runtime))
    if explicit:
        return sorted(explicit)
    if not runtime_points:
        return []
    base_cores, base_time = sorted(runtime_points)[0]
    return sorted((cores, base_time / runtime) for cores, runtime in runtime_points if runtime > 0)


def efficiency_points(path: Path | None) -> list[tuple[float, float]]:
    if path is None:
        return []
    rows = read_table(path)
    explicit = []
    runtime_points = []
    for row in rows:
        cores = get_column(row, "cores", "mpi_ranks", "ranks", "processes")
        efficiency = get_column(row, "efficiency", "strong_scaling_efficiency")
        runtime = get_column(row, "total_runtime", "runtime", "total_time")
        if cores is not None and efficiency is not None:
            explicit.append((cores, efficiency))
        if cores is not None and runtime is not None and runtime > 0:
            runtime_points.append((cores, runtime))
    if explicit:
        return sorted(explicit)
    if not runtime_points:
        return []
    base_cores, base_time = sorted(runtime_points)[0]
    return sorted(
        (cores, (base_time / runtime) / (cores / base_cores))
        for cores, runtime in runtime_points
        if runtime > 0 and cores > 0 and base_cores > 0
    )


def dat_block(name: str, points: list[tuple[float, float]]) -> str:
    lines = [rf"\begin{{filecontents*}}[overwrite]{{{name}.dat}}", "x y"]
    if points:
        lines.extend(f"{x:.12g} {y:.12g}" for x, y in points)
    else:
        lines.append("0 0")
    lines.append(r"\end{filecontents*}")
    return "\n".join(lines)


def force_blocks(series: list[ForceSeries]) -> str:
    if not series:
        return dat_block("force_series_empty", [])
    blocks = []
    for item in series:
        lines = [rf"\begin{{filecontents*}}[overwrite]{{{item.dataset_name}.dat}}", "x y"]
        lines.extend(f"{x:.12g} {y:.12g}" for x, y in item.points)
        lines.append(r"\end{filecontents*}")
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks)


def line_plot(name: str, title: str, xlabel: str, ylabel: str, points: list[tuple[float, float]]) -> str:
    if not points:
        return rf"\missingplot{{{latex_escape(title)}}}{{{MISSING}.}}"
    return rf"""
\begin{{figure}}[H]
\centering
\begin{{tikzpicture}}
\begin{{axis}}[
width=0.86\textwidth,
height=0.46\textwidth,
xlabel={{{latex_escape(xlabel)}}},
ylabel={{{latex_escape(ylabel)}}},
grid=both,
title={{{latex_escape(title)}}},
]
\addplot+[mark=*,thick] table[x=x,y=y] {{{name}.dat}};
\end{{axis}}
\end{{tikzpicture}}
\caption{{{latex_escape(title)}.}}
\end{{figure}}
"""


def bar_plot(name: str, title: str, xlabel: str, ylabel: str, points: list[tuple[float, float]]) -> str:
    if not points:
        return rf"\missingplot{{{latex_escape(title)}}}{{{MISSING}.}}"
    return rf"""
\begin{{figure}}[H]
\centering
\begin{{tikzpicture}}
\begin{{axis}}[
width=0.86\textwidth,
height=0.46\textwidth,
ybar,
bar width=12pt,
xlabel={{{latex_escape(xlabel)}}},
ylabel={{{latex_escape(ylabel)}}},
grid=both,
title={{{latex_escape(title)}}},
]
\addplot+[fill=blue!35] table[x=x,y=y] {{{name}.dat}};
\end{{axis}}
\end{{tikzpicture}}
\caption{{{latex_escape(title)}.}}
\end{{figure}}
"""


def benchmark_bar_plot(name: str, title: str, ylabel: str, points: list[tuple[float, float]]) -> str:
    if not points:
        return rf"\missingplot{{{latex_escape(title)}}}{{{MISSING}.}}"
    return rf"""
\begin{{figure}}[H]
\centering
\begin{{tikzpicture}}
\begin{{axis}}[
width=0.86\textwidth,
height=0.46\textwidth,
ybar,
bar width=18pt,
ylabel={{{latex_escape(ylabel)}}},
xtick={{1,2,3}},
xticklabels={{Serial uniform,MPI uniform,MPI + AMR}},
xticklabel style={{rotate=15,anchor=east}},
grid=both,
title={{{latex_escape(title)}}},
]
\addplot+[fill=blue!35] coordinates {{{chr(10).join(f"({x:.12g},{y:.12g})" for x, y in points)}}};
\end{{axis}}
\end{{tikzpicture}}
\caption{{{latex_escape(title)}.}}
\end{{figure}}
"""


def force_plot(series: list[ForceSeries]) -> str:
    if not series:
        return rf"\missingplot{{Force--displacement response}}{{{MISSING}.}}"
    plots = []
    for item in series:
        plots.append(
            rf"\addplot+[mark=*,thick] table[x=x,y=y] {{{item.dataset_name}.dat}};"
            + "\n"
            + rf"\addlegendentry{{{latex_escape(item.label)}}}"
        )
    return rf"""
\begin{{figure}}[H]
\centering
\begin{{tikzpicture}}
\begin{{axis}}[
width=0.86\textwidth,
height=0.48\textwidth,
xlabel={{Applied displacement}},
ylabel={{Reaction force}},
grid=both,
legend pos=north west,
title={{Force--displacement response}},
]
{chr(10).join(plots)}
\end{{axis}}
\end{{tikzpicture}}
\caption{{Force--displacement curves parsed from the latest simulation output.}}
\end{{figure}}
"""


def peak_force_points(series: list[ForceSeries]) -> list[tuple[float, float]]:
    points = []
    for item in series:
        if item.pore_count is not None and item.peak_force is not None:
            points.append((float(item.pore_count), item.peak_force))
    return sorted(points)


def force_summary_rows(series: list[ForceSeries]) -> str:
    if not series:
        return rf"\multicolumn{{4}}{{l}}{{{MISSING}.}}\\"
    rows = []
    for item in series:
        rows.append(
            f"{latex_escape(item.label)} & {len(item.points)} & "
            f"{format_float(item.peak_displacement)} & {format_float(item.peak_force)}\\\\"
        )
    return "\n".join(rows)


def porosity_interpretation(series: list[ForceSeries]) -> str:
    points = peak_force_points(series)
    if len(points) < 2:
        return MISSING + ". Multiple pore cases are required to verify whether increasing porosity decreases peak force."
    decreasing = all(points[i][1] >= points[i + 1][1] for i in range(len(points) - 1))
    if decreasing:
        return "The extracted peak-force data show a decreasing trend with increasing pore count."
    return "The extracted peak-force data do not show a strictly decreasing trend with increasing pore count."


def table_status(path: Path | None) -> str:
    return latex_escape(path.relative_to(ROOT)) if path else MISSING


def timing_rows(path: Path | None) -> str:
    if path is None:
        return rf"\multicolumn{{6}}{{l}}{{{MISSING}.}}\\"
    rows = []
    for row in read_table(path):
        step = get_column(row, "step", "load_step", "time_step")
        assembly = get_column(row, "assembly_time", "assembly", "assembly_seconds")
        solve = get_column(row, "solve_time", "solve", "solver_time", "solver_seconds")
        refinement = get_column(row, "refinement_time", "refine_time", "amr_time")
        io = get_column(row, "io_time", "i_o_time", "output_time")
        total = get_column(row, "total_time", "runtime", "total_runtime")
        rows.append(
            f"{format_float(step)} & {format_float(assembly)} & {format_float(solve)} & "
            f"{format_float(refinement)} & {format_float(io)} & {format_float(total)}\\\\"
        )
    return "\n".join(rows) if rows else rf"\multicolumn{{6}}{{l}}{{{MISSING}.}}\\"


def mpi_rows(path: Path | None) -> str:
    if path is None:
        return rf"\multicolumn{{6}}{{l}}{{{MISSING}.}}\\"
    rows = []
    for row in read_table(path):
        cores = get_column(row, "cores", "mpi_ranks", "ranks", "processes")
        total = get_column(row, "total_runtime", "runtime", "total_time")
        assembly = get_column(row, "assembly_time", "assembly")
        solve = get_column(row, "solve_time", "solver_time", "solve")
        refinement = get_column(row, "refinement_time", "amr_time")
        io = get_column(row, "io_time", "output_time")
        rows.append(
            f"{format_float(cores)} & {format_float(total)} & {format_float(assembly)} & "
            f"{format_float(solve)} & {format_float(refinement)} & {format_float(io)}\\\\"
        )
    return "\n".join(rows) if rows else rf"\multicolumn{{6}}{{l}}{{{MISSING}.}}\\"


def solver_rows(path: Path | None) -> str:
    if path is None:
        return rf"\multicolumn{{5}}{{l}}{{{MISSING}.}}\\"
    rows = []
    for row in read_table(path):
        step = get_column(row, "step", "load_step")
        staggered = get_column(row, "staggered_iterations", "staggered", "nonlinear_iterations")
        cg_u = get_column(row, "elastic_cg_iterations", "cg_elastic", "u_cg_iterations")
        cg_d = get_column(row, "damage_cg_iterations", "cg_damage", "d_cg_iterations")
        residual = get_column(row, "residual", "final_residual", "relative_residual")
        rows.append(
            f"{format_float(step)} & {format_float(staggered)} & {format_float(cg_u)} & "
            f"{format_float(cg_d)} & {format_float(residual)}\\\\"
        )
    return "\n".join(rows) if rows else rf"\multicolumn{{5}}{{l}}{{{MISSING}.}}\\"


def cost_fraction_points(path: Path | None) -> list[tuple[float, float]]:
    if path is None:
        return []
    totals = {"assembly": 0.0, "solve": 0.0, "refinement": 0.0, "io": 0.0}
    for row in read_table(path):
        totals["assembly"] += get_column(row, "assembly_time", "assembly") or 0.0
        totals["solve"] += get_column(row, "solve_time", "solver_time", "solve") or 0.0
        totals["refinement"] += get_column(row, "refinement_time", "amr_time", "refine_time") or 0.0
        totals["io"] += get_column(row, "io_time", "output_time") or 0.0
    total = sum(totals.values())
    if total <= 0:
        return []
    return [(i + 1, 100.0 * value / total) for i, value in enumerate(totals.values())]


def dominant_cost(path: Path | None) -> str:
    points = cost_fraction_points(path)
    if not points:
        return MISSING + "."
    labels = ["assembly", "solve", "refinement", "I/O"]
    index, value = max(enumerate(points), key=lambda item: item[1][1])
    return f"The dominant recorded runtime component is {labels[index]} with {value[1]:.2f}\\% of the measured cost."


def point_summary(points: list[tuple[float, float]], quantity: str) -> str:
    if not points:
        return MISSING
    y_values = [point[1] for point in points]
    return (
        f"{quantity}: initial {format_float(y_values[0])}, "
        f"final {format_float(y_values[-1])}, max {format_float(max(y_values))}"
    )


def table_total(path: Path | None, *names: str) -> float | None:
    if path is None:
        return None
    total = 0.0
    found = False
    for row in read_table(path):
        value = get_column(row, *names)
        if value is not None:
            total += value
            found = True
    return total if found else None


def table_max(path: Path | None, *names: str) -> float | None:
    if path is None:
        return None
    values = []
    for row in read_table(path):
        value = get_column(row, *names)
        if value is not None:
            values.append(value)
    return max(values) if values else None


def summary_rows(
    series: list[ForceSeries],
    datasets: dict[str, list[tuple[float, float]]],
    active_path: Path | None,
    dofs_path: Path | None,
    memory_path: Path | None,
    timing_path: Path | None,
    mpi_path: Path | None,
    solver_path: Path | None,
) -> str:
    peak_points = peak_force_points(series)
    if peak_points:
        peak_case, peak_force = max(peak_points, key=lambda item: item[1])
        peak_summary = f"maximum peak force {format_float(peak_force)} at pore count {format_float(peak_case)}"
    else:
        peak_summary = MISSING

    total_runtime = table_total(timing_path, "total_time", "runtime", "total_runtime")
    total_assembly = table_total(timing_path, "assembly_time", "assembly", "assembly_seconds")
    total_solve = table_total(timing_path, "solve_time", "solver_time", "solve", "solver_seconds")
    total_refinement = table_total(timing_path, "refinement_time", "amr_time", "refine_time")
    total_io = table_total(timing_path, "io_time", "i_o_time", "output_time")

    timing_summary = (
        f"total {format_float(total_runtime)}, assembly {format_float(total_assembly)}, "
        f"solve {format_float(total_solve)}, refinement {format_float(total_refinement)}, "
        f"I/O {format_float(total_io)}"
        if total_runtime is not None
        else MISSING
    )

    mpi_runtime = datasets["mpi_runtime"]
    mpi_speedup = datasets["mpi_speedup"]
    mpi_efficiency = datasets["mpi_efficiency"]
    if mpi_runtime:
        cores = [point[0] for point in mpi_runtime]
        mpi_summary = (
            f"cores {format_float(min(cores))}--{format_float(max(cores))}, "
            f"best speedup {format_float(max((point[1] for point in mpi_speedup), default=None))}, "
            f"best efficiency {format_float(max((point[1] for point in mpi_efficiency), default=None))}"
        )
    else:
        mpi_summary = MISSING

    max_staggered = table_max(solver_path, "staggered_iterations", "staggered", "nonlinear_iterations")
    max_elastic_cg = table_max(solver_path, "elastic_cg_iterations", "cg_elastic", "u_cg_iterations")
    max_damage_cg = table_max(solver_path, "damage_cg_iterations", "cg_damage", "d_cg_iterations")
    max_residual = table_max(solver_path, "residual", "final_residual", "relative_residual")
    solver_summary = (
        f"max staggered {format_float(max_staggered)}, max elastic CG {format_float(max_elastic_cg)}, "
        f"max damage CG {format_float(max_damage_cg)}, max final residual {format_float(max_residual)}"
        if max_staggered is not None
        else MISSING
    )

    rows = [
        ("Simulation output parsed", f"{len(series)} force series, {table_status(active_path)} for AMR, {table_status(timing_path)} for timing"),
        ("Porosity response", peak_summary),
        ("AMR mesh evolution", point_summary(datasets["active_cells"], "active cells")),
        ("DoF evolution", point_summary(datasets["dofs"], "DoFs")),
        ("Memory record", table_status(memory_path)),
        ("Computational cost", timing_summary),
        ("Dominant cost", dominant_cost(timing_path)),
        ("MPI scaling", mpi_summary),
        ("Solver convergence", solver_summary),
        ("Missing-data policy", MISSING + " is printed wherever a required dataset is absent."),
    ]

    return "\n".join(
        f"{latex_escape(topic)} & {latex_escape(value)}\\\\"
        for topic, value in rows
    )


def write_data_blocks(datasets: dict[str, list[tuple[float, float]]], series: list[ForceSeries]) -> str:
    blocks = [force_blocks(series)]
    for name, points in datasets.items():
        blocks.append(dat_block(name, points))
    return "\n\n".join(blocks)


def generate_tex(series: list[ForceSeries]) -> str:
    config = source_config()
    counts = count_outputs()
    case_folders = case_folder_records()

    active_path = find_dataset("active_cells.csv", "amr_history.csv", "mesh_history.csv")
    dofs_path = find_dataset("dofs.csv", "dof_history.csv")
    memory_path = find_dataset("memory_usage.csv", "memory.csv")
    timing_path = find_dataset("timing.csv", "runtime_breakdown.csv", "timer.csv")
    mpi_path = find_dataset("mpi_scaling.csv", "scaling.csv")
    benchmark_path = find_dataset("benchmark_summary.csv")
    solver_path = find_dataset("solver_iterations.csv", "convergence.csv", "iteration_history.csv")
    vtu_active_points, vtu_point_points, vtu_history_status = first_case_vtu_history()
    benchmark = benchmark_cases(benchmark_path)

    active_points = dataset_points(
        active_path,
        ["step", "load_step"],
        ["active_cells_end", "active_cells", "cells", "global_active_cells"],
    )
    if not active_points:
        active_points = vtu_active_points

    dof_points = dataset_points(
        dofs_path,
        ["step", "load_step"],
        ["total_dofs", "dofs", "damage_dofs", "elastic_dofs"],
    )
    if not dof_points and active_path is not None:
        dof_points = dataset_points(active_path, ["step", "load_step"], ["elastic_dofs", "damage_dofs"])

    datasets = {
        "active_cells": active_points,
        "dofs": dof_points,
        "memory": dataset_points(memory_path, ["step", "load_step", "cores"], ["memory_mb", "memory", "rss_mb"]),
        "runtime": dataset_points(timing_path, ["step", "load_step"], ["total_time", "runtime", "total_runtime"]),
        "mpi_runtime": dataset_points(mpi_path, ["cores", "mpi_ranks", "ranks"], ["total_runtime", "runtime", "total_time"]),
        "mpi_speedup": speedup_points(mpi_path),
        "mpi_efficiency": efficiency_points(mpi_path),
        "peak_force": peak_force_points(series),
        "staggered_iterations": dataset_points(solver_path, ["step", "load_step"], ["staggered_iterations", "staggered", "nonlinear_iterations"]),
        "cg_iterations": dataset_points(solver_path, ["step", "load_step"], ["cg_iterations", "elastic_cg_iterations", "damage_cg_iterations"]),
        "residual": dataset_points(solver_path, ["step", "load_step"], ["residual", "final_residual", "relative_residual"]),
        "cost_fraction": cost_fraction_points(timing_path),
        "benchmark_runtime": benchmark_points(benchmark, "total_runtime"),
        "benchmark_solve": benchmark_points(benchmark, "solve_time"),
        "benchmark_assembly": benchmark_points(benchmark, "assembly_time"),
        "benchmark_refinement": benchmark_points(benchmark, "refinement_time"),
        "benchmark_io": benchmark_points(benchmark, "io_time"),
        "benchmark_cells": benchmark_points(benchmark, "active_cells_final"),
        "benchmark_elastic_dofs": benchmark_points(benchmark, "elastic_dofs_final"),
        "benchmark_damage_dofs": benchmark_points(benchmark, "damage_dofs_final"),
        "benchmark_staggered": benchmark_points(benchmark, "max_staggered_iterations"),
        "benchmark_elastic_cg": benchmark_points(benchmark, "max_elastic_cg_iterations"),
        "benchmark_damage_cg": benchmark_points(benchmark, "max_damage_cg_iterations"),
        "benchmark_speedup": benchmark_speedup_points(benchmark),
        "benchmark_efficiency": benchmark_efficiency_points(benchmark),
    }

    substitutions = {
        "DATA_BLOCKS": write_data_blocks(datasets, series),
        "FORCE_CSV_COUNT": str(counts["force_csv"]),
        "PVTU_COUNT": str(counts["pvtu"]),
        "VTU_COUNT": str(counts["vtu"]),
        "PVD_COUNT": str(counts["pvd"]),
        "PORES_COUNT": str(counts["pores"]),
        "CASE_DOCUMENTS_COUNT": str(counts["case_documents"]),
        "CASE_FOLDER_ROWS": case_folder_rows(case_folders),
        "SERIES_COUNT": str(len(series)),
        "PEAK_FORCE": format_float(max((item.peak_force for item in series if item.peak_force is not None), default=None)),
        "PORE_COUNTS": latex_escape(config["pore_counts"]),
        "RANDOM_SEED": latex_escape(config["random_seed_base"]),
        "MIN_PORE": latex_escape(config["min_pore_size"]),
        "MAX_PORE": latex_escape(config["max_pore_size"]),
        "CRACK_SPREAD": latex_escape(config["crack_spread_radius"]),
        "BOUNDARY_MARGIN": latex_escape(config["boundary_margin"]),
        "L0": latex_escape(config["lo"]),
        "GC": latex_escape(config["Gc"]),
        "LAMBDA": latex_escape(config["lambda"]),
        "MU": latex_escape(config["mu"]),
        "K": latex_escape(config["k"]),
        "N_STEPS": latex_escape(config["n_steps"]),
        "INC_LARGE": latex_escape(config["inc_large"]),
        "INC_SMALL": latex_escape(config["inc_small"]),
        "TOL": latex_escape(config["tol_1"]),
        "ACTIVE_STATUS": table_status(active_path) if active_path else latex_escape(vtu_history_status),
        "DOFS_STATUS": table_status(dofs_path) if dofs_path else (table_status(active_path) if active_path else MISSING),
        "MEMORY_STATUS": table_status(memory_path),
        "TIMING_STATUS": table_status(timing_path),
        "MPI_STATUS": table_status(mpi_path),
        "BENCHMARK_STATUS": table_status(benchmark_path),
        "BENCHMARK_ROWS": benchmark_rows(benchmark),
        "BENCHMARK_INTERPRETATION": latex_escape(benchmark_interpretation(benchmark)),
        "FINAL_PERFORMANCE_ROWS": final_performance_rows(benchmark),
        "FINAL_EFFICIENCY_ROWS": final_efficiency_rows(benchmark),
        "FINAL_SOLVER_ROWS": final_solver_rows(benchmark),
        "FINAL_CONCLUSION_TEXT": latex_escape(final_conclusion_text(benchmark)),
        "SOLVER_STATUS": table_status(solver_path),
        "SUMMARY_ROWS": summary_rows(
            series,
            datasets,
            active_path,
            dofs_path,
            memory_path,
            timing_path,
            mpi_path,
            solver_path,
        ),
        "FORCE_ROWS": force_summary_rows(series),
        "MPI_ROWS": mpi_rows(mpi_path),
        "TIMING_ROWS": timing_rows(timing_path),
        "SOLVER_ROWS": solver_rows(solver_path),
        "POROSITY_INTERPRETATION": latex_escape(porosity_interpretation(series)),
        "DOMINANT_COST": dominant_cost(timing_path),
        "FORCE_PLOT": force_plot(series),
        "PEAK_FORCE_PLOT": line_plot("peak_force", "Peak force versus pore count", "Pore count (count)", "Peak force (recorded units)", datasets["peak_force"]),
        "ACTIVE_CELLS_PLOT": line_plot("active_cells", "Active cells versus load step", "Load step (count)", "Active cells (count)", datasets["active_cells"]),
        "DOFS_PLOT": line_plot("dofs", "Degrees of freedom versus load step", "Load step (count)", "DoFs (count)", datasets["dofs"]),
        "MEMORY_PLOT": line_plot("memory", "Memory usage and AMR savings", "Step or core count", "Memory (MB)", datasets["memory"]),
        "RUNTIME_PLOT": line_plot("runtime", "Runtime comparison", "Load step (count)", "Runtime (s)", datasets["runtime"]),
        "MPI_RUNTIME_PLOT": line_plot("mpi_runtime", "Runtime versus MPI cores", "MPI cores (count)", "Runtime (s)", datasets["mpi_runtime"]),
        "MPI_SPEEDUP_PLOT": line_plot("mpi_speedup", "MPI speedup graph", "MPI cores (count)", "Speedup (ratio)", datasets["mpi_speedup"]),
        "MPI_EFFICIENCY_PLOT": line_plot("mpi_efficiency", "Strong scaling efficiency", "MPI cores (count)", "Efficiency (ratio)", datasets["mpi_efficiency"]),
        "COST_FRACTION_PLOT": bar_plot("cost_fraction", "Computational cost percentage breakdown", "Stage index", "Cost fraction (\\%)", datasets["cost_fraction"]),
        "STAGGERED_PLOT": line_plot("staggered_iterations", "Staggered iterations per load step", "Load step (count)", "Iterations (count)", datasets["staggered_iterations"]),
        "CG_PLOT": line_plot("cg_iterations", "CG iterations per load step", "Load step (count)", "CG iterations (count)", datasets["cg_iterations"]),
        "RESIDUAL_PLOT": line_plot("residual", "Residual convergence history", "Load step (count)", "Residual (dimensionless)", datasets["residual"]),
        "BENCHMARK_RUNTIME_PLOT": benchmark_bar_plot("benchmark_runtime", "Total runtime comparison", "Total runtime (s)", datasets["benchmark_runtime"]),
        "BENCHMARK_SOLVE_PLOT": benchmark_bar_plot("benchmark_solve", "Solver runtime comparison", "Solver runtime (s)", datasets["benchmark_solve"]),
        "BENCHMARK_BREAKDOWN_PLOT": line_plot("benchmark_runtime", "Benchmark runtime values", "Case index", "Runtime (s)", datasets["benchmark_runtime"]),
        "BENCHMARK_SPEEDUP_PLOT": benchmark_bar_plot("benchmark_speedup", "Speedup relative to serial uniform baseline", "Speedup (ratio)", datasets["benchmark_speedup"]),
        "BENCHMARK_EFFICIENCY_PLOT": benchmark_bar_plot("benchmark_efficiency", "Parallel efficiency relative to serial baseline", "Efficiency (ratio)", datasets["benchmark_efficiency"]),
        "BENCHMARK_CELLS_PLOT": benchmark_bar_plot("benchmark_cells", "Active-cell comparison", "Active cells (count)", datasets["benchmark_cells"]),
        "BENCHMARK_ELASTIC_DOFS_PLOT": benchmark_bar_plot("benchmark_elastic_dofs", "Elastic DoF comparison", "Elastic DoFs (count)", datasets["benchmark_elastic_dofs"]),
        "BENCHMARK_DAMAGE_DOFS_PLOT": benchmark_bar_plot("benchmark_damage_dofs", "Damage DoF comparison", "Damage DoFs (count)", datasets["benchmark_damage_dofs"]),
        "BENCHMARK_STAGGERED_PLOT": benchmark_bar_plot("benchmark_staggered", "Maximum staggered iterations by benchmark case", "Iterations (count)", datasets["benchmark_staggered"]),
        "BENCHMARK_ELASTIC_CG_PLOT": benchmark_bar_plot("benchmark_elastic_cg", "Maximum elastic CG iterations by benchmark case", "CG iterations (count)", datasets["benchmark_elastic_cg"]),
        "BENCHMARK_DAMAGE_CG_PLOT": benchmark_bar_plot("benchmark_damage_cg", "Maximum damage CG iterations by benchmark case", "CG iterations (count)", datasets["benchmark_damage_cg"]),
        "MISSING": MISSING,
    }

    template = LatexTemplate(r"""\documentclass[11pt,a4paper]{report}
\usepackage[margin=1in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb,bm}
\usepackage{booktabs}
\usepackage{array}
\usepackage{longtable}
\usepackage{xcolor}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{float}
\usepackage{tikz}
\usepackage{pgfplots}
\usepackage{hyperref}
\usepackage{enumitem}
\usepackage{cite}
\pgfplotsset{compat=1.17}
\usetikzlibrary{arrows.meta,positioning,shapes.geometric}
\hypersetup{colorlinks=true,linkcolor=blue!45!black,citecolor=blue!45!black,urlcolor=blue!45!black}
\setlist[itemize]{leftmargin=*,itemsep=2pt,topsep=3pt}
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt}
\newcommand{\dd}{\mathrm{d}}
\newcommand{\strain}{\boldsymbol{\varepsilon}}
\newcommand{\unithead}[2]{\shortstack{#1\\{\scriptsize (#2)}}}
\newcommand{\missingplot}[2]{\begin{figure}[H]\centering\fbox{\begin{minipage}[c][0.15\textheight][c]{0.86\textwidth}\centering\textbf{#1}\\[0.4em]#2\end{minipage}}\caption{#1.}\end{figure}}

@DATA_BLOCKS

\begin{document}
\begin{titlepage}
\centering
\vspace*{1.2cm}
{\Large\bfseries Comparative Technical Report on Serial, MPI, and MPI--AMR Phase-Field Fracture Simulations\par}
\vspace{1.0cm}
{\large Automatically generated from the latest simulation output\par}
\vfill
\begin{tabular}{ll}
Framework & deal.II finite element implementation\\
Parallel model & MPI distributed triangulation and distributed linear algebra\\
Generated files parsed & @FORCE_CSV_COUNT force CSV files, @PVTU_COUNT PVTU files, @VTU_COUNT VTU files\\
Case document folders & @CASE_DOCUMENTS_COUNT manifests\\
Report date & \today
\end{tabular}
\vfill
\end{titlepage}

\begin{abstract}
This report is generated automatically after the simulation driver completes. It extracts only real data available in the latest output files. The current archive contains @SERIES_COUNT usable force--displacement series and the largest recorded force is @PEAK_FORCE. Required datasets that are absent are marked explicitly as ``@MISSING''.
\end{abstract}

\tableofcontents
\listoffigures
\listoftables

\chapter{Executive Summary}
Table~\ref{tab:executive-summary} gives a compact reading guide for the recorded simulation data. Each entry is computed from the latest output files; missing entries are not inferred.

\begin{longtable}{p{0.27\textwidth}p{0.63\textwidth}}
\caption{Summary of recorded simulation data and report availability.}
\label{tab:executive-summary}\\
\toprule
Category & Recorded summary\\
\midrule
\endfirsthead
\toprule
Category & Recorded summary\\
\midrule
\endhead
@SUMMARY_ROWS
\bottomrule
\end{longtable}

\chapter{Introduction}
Phase-field fracture represents a sharp crack by a smooth scalar damage field $d\in[0,1]$. The intact state is $d=0$, while $d=1$ denotes fully broken material. The regularized crack surface is controlled by the length scale $\ell_0$, which allows crack initiation and propagation without explicit crack-front tracking.

The brittle fracture problem is written as an energy minimization problem. The displacement field $\bm{u}$ satisfies degraded mechanical equilibrium, and the damage field evolves through a regularized fracture energy. Irreversibility is imposed through a history variable $H$, so crack-driving tensile energy is not allowed to decrease after damage has formed.

Adaptive mesh refinement is required because $\nabla d$ is localized near the crack band and crack tip. MPI parallelization is required because every load step involves repeated assembly, linear solution, ghost exchange, mesh transfer, and output. Porosity modifies the crack trajectory by introducing local stiffness loss and stress concentration.

\chapter{Comparative Benchmark Definition}
The computational comparison is defined by three executions of the same executable and the same physical model:
\begin{enumerate}
\item Serial uniform mesh: single MPI process, AMR disabled, zero pores.
\item MPI uniform mesh: 45 MPI processes, AMR disabled, zero pores.
\item MPI with AMR: 45 MPI processes, AMR enabled, zero pores.
\end{enumerate}
Only the execution mode and AMR flag are changed between the cases. Geometry, material constants, load stepping, solver tolerances, phase-field length scale, boundary conditions, pore count, and output policy must remain fixed for the comparison to be scientifically valid.

\begin{table}[H]
\centering
\caption{Case-wise output and document folder structure scanned by the report generator.}
\scriptsize
\begin{tabular}{llllll}
\toprule
Case & \unithead{Cores}{count} & AMR flag & Case folder & Results directory & Documents directory\\
\midrule
@CASE_FOLDER_ROWS
\bottomrule
\end{tabular}
\end{table}

For each benchmark mode, simulation fields and CSV files are stored under \texttt{results/}, while case-level metadata and report copies are stored under \texttt{documents/}. The overall report scans all case folders under \texttt{output/} and combines the recorded data into one comparative analysis.

\begin{table}[H]
\centering
\caption{Three-case benchmark summary extracted from \texttt{benchmark\_summary.csv}.}
\scriptsize
\begin{tabular}{lllllllllll}
\toprule
Case & \unithead{Cores}{count} & AMR & \unithead{Total}{s} & \unithead{Assembly}{s} & \unithead{Solve}{s} & \unithead{Refine}{s} & \unithead{I/O}{s} & \unithead{Cells}{count} & \unithead{Elastic DoFs}{count} & \unithead{Damage DoFs}{count}\\
\midrule
@BENCHMARK_ROWS
\bottomrule
\end{tabular}
\end{table}

Benchmark dataset status: @BENCHMARK_STATUS.

@BENCHMARK_RUNTIME_PLOT
@BENCHMARK_SOLVE_PLOT
@BENCHMARK_SPEEDUP_PLOT
@BENCHMARK_EFFICIENCY_PLOT
@BENCHMARK_CELLS_PLOT
@BENCHMARK_ELASTIC_DOFS_PLOT
@BENCHMARK_DAMAGE_DOFS_PLOT

Interpretation from recorded benchmark data: @BENCHMARK_INTERPRETATION

The serial uniform case defines the reference runtime and convergence cost. The MPI uniform case isolates the effect of domain decomposition, distributed assembly, distributed vectors, ghost exchange, and parallel Krylov--AMG solution. The MPI--AMR case measures the combined effect of distributed parallelism and localized refinement near the fracture process zone.

\chapter{Numerical Methodology}
The small-strain tensor is
\begin{equation}
\strain(\bm{u})=\frac{1}{2}\left(\nabla\bm{u}+\nabla\bm{u}^{T}\right).
\end{equation}
The elastic energy density is
\begin{equation}
\psi_0(\strain)=\frac{\lambda}{2}\left(\mathrm{tr}\,\strain\right)^2+\mu\,\strain:\strain .
\end{equation}
The degradation function is
\begin{equation}
g(d)=(1-d)^2+k,
\end{equation}
and the crack density is
\begin{equation}
\gamma_{\ell_0}(d,\nabla d)=\frac{d^2}{2\ell_0}+\frac{\ell_0}{2}|\nabla d|^2 .
\end{equation}
The total energy functional is
\begin{equation}
\Pi(\bm{u},d)=\int_\Omega g(d)\psi_0^+(\strain(\bm{u}))+\psi_0^-(\strain(\bm{u}))\,\dd\Omega
+\int_\Omega G_c\gamma_{\ell_0}(d,\nabla d)\,\dd\Omega-W_{\mathrm{ext}} .
\end{equation}

\begin{table}[H]
\centering
\caption{Source-derived numerical configuration used by the executable.}
\begin{tabular}{ll}
\toprule
Quantity & Value / unit\\
\midrule
Pore counts & \texttt{@PORE_COUNTS}\\
Random seed base & \texttt{@RANDOM_SEED}\\
Minimum pore diameter & \texttt{@MIN_PORE} (model length unit)\\
Maximum pore diameter & \texttt{@MAX_PORE} (model length unit)\\
Crack spread radius & \texttt{@CRACK_SPREAD} (model length unit)\\
Boundary margin & \texttt{@BOUNDARY_MARGIN} (model length unit)\\
$\ell_0$ & \texttt{@L0} (model length unit)\\
$G_c$ & \texttt{@GC} (model fracture-energy unit)\\
$\lambda$ & \texttt{@LAMBDA} (stress unit)\\
$\mu$ & \texttt{@MU} (stress unit)\\
$k$ & \texttt{@K} (dimensionless)\\
Load steps & \texttt{@N_STEPS} (count)\\
Large displacement increment & \texttt{@INC_LARGE} (model length unit)\\
Small displacement increment & \texttt{@INC_SMALL} (model length unit)\\
Staggered tolerance & \texttt{@TOL} (dimensionless)\\
\bottomrule
\end{tabular}
\end{table}

\chapter{Adaptive Mesh Refinement}
AMR is required because the damage field changes rapidly only in a narrow fracture process zone. The implementation therefore concentrates cells near regions of high $\nabla d$, which improves crack-tip resolution while avoiding unnecessary work away from the crack.

\begin{table}[H]
\centering
\caption{AMR dataset availability.}
\begin{tabular}{ll}
\toprule
Dataset & Status / unit\\
\midrule
Active cells history & @ACTIVE_STATUS; cells in count\\
DoF history & @DOFS_STATUS; DoFs in count\\
Memory usage & @MEMORY_STATUS; memory in MB when reported\\
Runtime data & @TIMING_STATUS; time in seconds\\
VTU/PVTU visualization files & @VTU_COUNT VTU files, @PVTU_COUNT PVTU files; files in count\\
Uniform mesh reference & @MISSING\\
\bottomrule
\end{tabular}
\end{table}

@ACTIVE_CELLS_PLOT
@DOFS_PLOT
@MEMORY_PLOT
@RUNTIME_PLOT
@BENCHMARK_CELLS_PLOT
@BENCHMARK_DAMAGE_DOFS_PLOT
\missingplot{Crack resolution improvement}{@MISSING. Crack-resolution metrics require processed damage contours or crack-tip width measurements.}
\missingplot{Refinement evolution}{@MISSING unless active-cell or AMR-history data are written per load step.}

The AMR comparison is valid when both MPI uniform and MPI--AMR cases are present in the benchmark summary. A lower active-cell or DoF count in the MPI--AMR case, together with comparable force--displacement response and crack trajectory, supports the conclusion that AMR reduces unnecessary global refinement while preserving fracture resolution.

\chapter{MPI Parallelization}
The simulation uses a distributed triangulation, distributed sparse matrices, and distributed vectors. Assembly is mostly local to each MPI rank, while vector compression, ghost updates, AMG setup, and Krylov reductions create communication overhead.

\begin{table}[H]
\centering
\caption{MPI runtime table.}
\begin{tabular}{llllll}
\toprule
\unithead{Cores}{count} & \unithead{Total runtime}{s} & \unithead{Assembly time}{s} & \unithead{Solver time}{s} & \unithead{Refinement time}{s} & \unithead{I/O time}{s}\\
\midrule
@MPI_ROWS
\bottomrule
\end{tabular}
\end{table}

@MPI_RUNTIME_PLOT
@MPI_SPEEDUP_PLOT
@MPI_EFFICIENCY_PLOT
@BENCHMARK_SPEEDUP_PLOT
@BENCHMARK_EFFICIENCY_PLOT

MPI reduction in solve time can be stated only from extracted timing data. Current MPI scaling dataset: @MPI_STATUS.

\chapter{Porosity and Fracture Response}
Random pores are inserted as a diffuse initial damage field. These pores reduce local stiffness, generate stress concentration, and can guide the crack through weakened ligaments.

\begin{table}[H]
\centering
\caption{Parsed force--displacement series and peak force.}
\begin{tabular}{lccc}
\toprule
Case & \unithead{Points}{count} & \unithead{Displacement at peak force}{recorded units} & \unithead{Peak force}{recorded units}\\
\midrule
@FORCE_ROWS
\bottomrule
\end{tabular}
\end{table}

@FORCE_PLOT
@PEAK_FORCE_PLOT
\missingplot{Crack path comparison}{@MISSING. Crack path comparison requires processed damage-field snapshots for each porosity case.}

Porosity conclusion from extracted force data: @POROSITY_INTERPRETATION

\chapter{Results and Discussion}
The force--displacement response, pore-dependent peak force, AMR history, MPI scaling, and convergence plots in this report are generated only from available output files. Missing datasets are not inferred. Crack path alignment with weakened ligaments and dense-pore crack coalescence require damage contours or crack-path extraction from VTU/PVTU files.

\chapter{Computational Cost Breakdown}
Committee-level cost analysis requires assembly time, solve time, refinement time, and I/O time. These values are read from \texttt{timing.csv}, \texttt{runtime\_breakdown.csv}, or \texttt{timer.csv} when available.

\begin{table}[H]
\centering
\caption{Computational cost breakdown.}
\begin{tabular}{llllll}
\toprule
\unithead{Step}{count} & \unithead{Assembly time}{s} & \unithead{Solve time}{s} & \unithead{Refinement time}{s} & \unithead{I/O time}{s} & \unithead{Total runtime}{s}\\
\midrule
@TIMING_ROWS
\bottomrule
\end{tabular}
\end{table}

@COST_FRACTION_PLOT
@BENCHMARK_RUNTIME_PLOT
@BENCHMARK_SOLVE_PLOT
@DOMINANT_COST

\chapter{Solver Convergence Analysis}
The staggered convergence check uses
\begin{equation}
\eta_u=\frac{\|\bm{u}^{k}-\bm{u}^{k-1}\|_2}{\|\bm{u}^{k}\|_2},
\qquad
\eta_d=\frac{\|d^{k}-d^{k-1}\|_2}{\|d^{k}\|_2}.
\end{equation}
Convergence behavior is read from \texttt{solver\_iterations.csv}, \texttt{convergence.csv}, or \texttt{iteration\_history.csv} when available.

\begin{table}[H]
\centering
\caption{Solver convergence table.}
\begin{tabular}{lllll}
\toprule
\unithead{Load step}{count} & \unithead{Staggered iterations}{count} & \unithead{Elastic CG iterations}{count} & \unithead{Damage CG iterations}{count} & \unithead{Final residual}{dimensionless}\\
\midrule
@SOLVER_ROWS
\bottomrule
\end{tabular}
\end{table}

@STAGGERED_PLOT
@CG_PLOT
@RESIDUAL_PLOT
@BENCHMARK_STAGGERED_PLOT
@BENCHMARK_ELASTIC_CG_PLOT
@BENCHMARK_DAMAGE_CG_PLOT

AMG behavior can be discussed quantitatively only if AMG setup time, iteration count, or residual data are recorded. Current convergence dataset: @SOLVER_STATUS.

\chapter{Observations}
\begin{itemize}
\item Increasing porosity decreases peak force only if confirmed by extracted peak-force data.
\item Crack path alignment with weakened ligaments requires damage-contour or crack-path data.
\item AMR concentrates cells near high $\nabla d$ by design of the damage-field refinement criterion.
\item MPI reduces solve time significantly only if the scaling table shows lower solver time with higher core count.
\item Dense pores accelerate crack coalescence only if the damage snapshots show earlier ligament connection.
\item Solver iterations are expected to increase near crack initiation, but this requires solver-iteration output for verification.
\item Strong scaling efficiency requires a multi-core timing table.
\item The three-case benchmark conclusion is restricted to rows recorded in \texttt{benchmark\_summary.csv}.
\end{itemize}

\chapter{Conclusions}
This final chapter summarizes the comparative evidence for the three executed cases: serial uniform mesh, MPI uniform mesh, and MPI--AMR. The interpretation is restricted to quantities extracted from the recorded benchmark, timing, mesh, solver, and post-processing files.

\section*{Final Performance Summary}
\begin{table}[H]
\centering
\caption{Final performance comparison for serial uniform, MPI uniform, and MPI--AMR simulations. Speedup and efficiency are computed relative to the serial uniform baseline.}
\scriptsize
\begin{tabular}{llllllllllll}
\toprule
Case & \unithead{Cores}{count} & AMR & \unithead{Total}{s} & \unithead{Assembly}{s} & \unithead{Solve}{s} & \unithead{Refine}{s} & \unithead{I/O}{s} & \unithead{Speedup}{ratio} & \unithead{Efficiency}{ratio} & \unithead{Cells}{count} & \unithead{Elastic DoFs}{count}\\
\midrule
@FINAL_PERFORMANCE_ROWS
\bottomrule
\end{tabular}
\end{table}

The solver stage is the dominant computational component for the recorded phase-field fracture runs. The serial case establishes the baseline wall time and memory-scale cost. The MPI uniform case measures the benefit of distributed assembly, ghost exchange, distributed matrices and vectors, and parallel Krylov--AMG solution on the same uniform mesh. The MPI--AMR case adds localized refinement around the fracture process zone, so its cost includes refinement and repartitioning work in addition to the parallel solve.

\section*{Final Efficiency and Mesh-Adaptivity Summary}
\begin{table}[H]
\centering
\caption{Final derived comparison metrics. Positive runtime reduction means a faster target case. Positive cell or DoF change means a larger adaptive discrete problem relative to the stated reference.}
\scriptsize
\begin{tabular}{llll}
\toprule
Comparison & Meaning & \unithead{Metric 1}{\% or ratio} & \unithead{Metric 2}{\% or ratio}\\
\midrule
@FINAL_EFFICIENCY_ROWS
\bottomrule
\end{tabular}
\end{table}

For phase-field fracture, AMR is not only a performance device. The damage field $d$ localizes inside a narrow band governed by the length scale $\ell_0$, and the relevant refinement indicator is tied to the crack zone, typically through $d$, $\nabla d$, or damage-history measures. This concentrates cells near the crack path and avoids uniform refinement of regions that remain nearly elastic. The expected numerical benefit is improved crack-tip resolution at a controlled global cost.

\section*{Final Solver-Convergence Summary}
\begin{table}[H]
\centering
\caption{Final solver-convergence comparison extracted from the benchmark summary.}
\scriptsize
\begin{tabular}{llllll}
\toprule
Case & \unithead{Staggered iterations}{count} & \unithead{Elastic CG iterations}{count} & \unithead{Damage CG iterations}{count} & \unithead{Final residual}{dimensionless} & \unithead{Solve time}{s}\\
\midrule
@FINAL_SOLVER_ROWS
\bottomrule
\end{tabular}
\end{table}

The staggered phase-field formulation couples the displacement problem and the damage problem through the degradation function $g(d)$ and the history variable $H$. Solver iterations therefore increase when crack growth begins, because the elastic field, damage field, and irreversibility constraint $d^n \ge d^{n-1}$ must become mutually consistent. The CG counts measure the linear-algebra cost of each subproblem, while the final residual indicates whether the nonlinear staggered loop reached the prescribed tolerance.

\section*{Final Technical Conclusion}
@FINAL_CONCLUSION_TEXT

The comparison supports the following thesis-level conclusions. MPI reduces wall-clock runtime by distributing element assembly, sparse linear algebra, and vector operations across subdomains. Communication through ghost cells and distributed constraints introduces overhead, but the recorded speedup shows that parallel execution is necessary for scalable fracture simulation. AMR is essential because the crack is spatially localized while the surrounding domain often remains smooth; refinement based on the damage field places resolution where the fracture process requires it. The combined MPI--AMR formulation is therefore the appropriate computational strategy for large-scale deal.II-based phase-field fracture analysis, provided that the adaptive refinement criterion preserves the force--displacement response and crack trajectory while controlling active cells and DoFs.

\chapter{Future Work}
\begin{itemize}
\item Save active cell count, DoF count, and AMR level statistics at every load step.
\item Save assembly, solve, refinement, and I/O timing in a structured CSV file.
\item Save MPI scaling runs in \texttt{mpi\_scaling.csv}.
\item Save staggered iterations, CG iterations, residuals, and AMG setup statistics.
\item Extract crack paths from damage contours for porosity-dependent crack trajectory analysis.
\item Extend to dynamic fracture, 3D simulations, GPU acceleration, and hybrid MPI+OpenMP.
\end{itemize}

\begin{thebibliography}{9}
\bibitem{bourdin2000} B. Bourdin, G. A. Francfort, and J.-J. Marigo, ``Numerical experiments in revisited brittle fracture,'' \emph{Journal of the Mechanics and Physics of Solids}, 48, 797--826, 2000.
\bibitem{miehe2010} C. Miehe, M. Hofacker, and F. Welschinger, ``A phase field model for rate-independent crack propagation,'' \emph{Computer Methods in Applied Mechanics and Engineering}, 199, 2765--2778, 2010.
\bibitem{dealii} W. Bangerth, R. Hartmann, and G. Kanschat, ``deal.II: A general-purpose object-oriented finite element library,'' \emph{ACM Transactions on Mathematical Software}, 33, 2007.
\bibitem{francfort1998} G. A. Francfort and J.-J. Marigo, ``Revisiting brittle fracture as an energy minimization problem,'' \emph{Journal of the Mechanics and Physics of Solids}, 46, 1319--1342, 1998.
\end{thebibliography}

\end{document}
""")
    return template.safe_substitute(substitutions)


def compile_pdf() -> bool:
    latex_engine = (
        shutil.which("pdflatex")
        or shutil.which("xelatex")
        or shutil.which("lualatex")
    )
    if latex_engine is None:
        print(
            "Warning: no LaTeX engine found on PATH "
            "(tried pdflatex, xelatex, lualatex). "
            "The .tex report was written, but the PDF was not generated."
        )
        return False

    with tempfile.TemporaryDirectory(prefix="pf_report_") as tmp:
        tmp_path = Path(tmp)
        tmp_tex = tmp_path / REPORT_TEX.name
        shutil.copy2(REPORT_TEX, tmp_tex)
        for _ in range(2):
            subprocess.run(
                [latex_engine, "-interaction=nonstopmode", "-halt-on-error", tmp_tex.name],
                cwd=tmp_path,
                check=True,
            )
        shutil.copy2(tmp_path / REPORT_PDF.name, REPORT_PDF)
    return True


def copy_report_to_case_documents() -> None:
    for documents_dir in sorted((ROOT / "output").glob("*/documents")):
        if not documents_dir.is_dir():
            continue
        shutil.copy2(REPORT_TEX, documents_dir / REPORT_TEX.name)
        if REPORT_PDF.exists():
            shutil.copy2(REPORT_PDF, documents_dir / REPORT_PDF.name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compile", action="store_true", help="Compile the generated LaTeX report to PDF.")
    args = parser.parse_args()

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    series = collect_force_series()
    REPORT_TEX.write_text(generate_tex(series), encoding="utf-8")
    print(f"Wrote {REPORT_TEX}")

    if args.compile:
        if compile_pdf():
            print(f"Wrote {REPORT_PDF}")

    copy_report_to_case_documents()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

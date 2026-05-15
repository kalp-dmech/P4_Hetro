#!/usr/bin/env python3
"""Create a PNG comparison of validation and generated force-displacement data."""

import argparse
import csv
import math
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
import matplotlib.pyplot as plt

OFFICE_NS = "urn:oasis:names:tc:opendocument:xmlns:office:1.0"
TABLE_NS = "urn:oasis:names:tc:opendocument:xmlns:table:1.0"
SERIES_COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#000000",
]
ROOT = Path(__file__).resolve().parent

def _attr(ns, name):
    return f"{{{ns}}}{name}"

def _to_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None

def read_ods_columns(path, sheet_name):
    path = Path(path) if Path(path).is_absolute() else ROOT / path
    if not path.exists():
        print(f"Warning: Validation file {path} not found. Skipping validation data.")
        return {}
    
    try:
        with zipfile.ZipFile(path) as ods:
            root = ET.fromstring(ods.read("content.xml"))
    except Exception as e:
        print(f"Warning: Could not read {path}: {e}")
        return {}

    table = None
    for candidate in root.iter(_attr(TABLE_NS, "table")):
        if candidate.attrib.get(_attr(TABLE_NS, "name")) == sheet_name:
            table = candidate
            break
    if table is None:
        print(f"Warning: Sheet {sheet_name!r} not found in {path}")
        return {}

    rows = []
    for row in table.findall(_attr(TABLE_NS, "table-row")):
        repeated_rows = int(row.attrib.get(_attr(TABLE_NS, "number-rows-repeated"), "1"))
        values = []
        for cell in row.findall(_attr(TABLE_NS, "table-cell")):
            repeated_cols = int(cell.attrib.get(_attr(TABLE_NS, "number-columns-repeated"), "1"))
            value = cell.attrib.get(_attr(OFFICE_NS, "value"))
            if value is None:
                value = "".join(cell.itertext()).strip()
            values.extend([value] * min(repeated_cols, 20))
        if any(values):
            rows.extend([values] * min(repeated_rows, 20))

    if not rows:
        return {}

    headers = [cell.strip() for cell in rows[0]]
    columns = {header: [] for header in headers if header}
    for row in rows[1:]:
        for index, header in enumerate(headers):
            if not header or index >= len(row):
                continue
            number = _to_float(row[index])
            if number is not None:
                columns[header].append(number)
    return columns

def read_generated(path):
    path = Path(path) if Path(path).is_absolute() else ROOT / path
    if not path.exists():
        raise FileNotFoundError(path)

    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not text:
        raise ValueError(f"{path} is empty")

    delimiter = "," if "," in text[0] else None
    x, y = [], []
    if delimiter == ",":
        reader = csv.DictReader(text)
        for row in reader:
            disp = _to_float(row.get("displacement", "") or row.get("Displacement", ""))
            force = _to_float(row.get("force", "") or row.get("Force", ""))
            if disp is not None and force is not None:
                x.append(disp)
                y.append(force)
    else:
        for line in text:
            parts = line.split()
            if len(parts) < 2: continue
            disp, force = _to_float(parts[0]), _to_float(parts[1])
            if disp is not None and force is not None:
                x.append(disp)
                y.append(force)
    return x, y

def make_plot(series, out, title):
    plt.figure(figsize=(10, 6))
    for name, xs, ys, color in series:
        plt.plot(xs, ys, label=name, color=color, linewidth=2, alpha=0.85)
    
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel("Displacement (mm)", fontsize=12)
    plt.ylabel("Force (kN)", fontsize=12)
    plt.legend(loc='best')
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Set axis to start from 0
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    
    # Add minor grid for better readability
    plt.minorticks_on()
    plt.grid(True, which='minor', linestyle=':', alpha=0.3)
    
    plt.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=300)
    plt.close()

def plot_one(validation, generated, output, title):
    generated_x, generated_y = read_generated(generated)
    series = []
    
    if validation:
        validation_x = validation.get("Dis(Paper)", [])
        validation_y = validation.get("Force (paper)", [])
        if validation_x and validation_y:
            count = min(len(validation_x), len(validation_y))
            series.append(("Reference (Paper)", validation_x[:count], validation_y[:count], "#0072B2"))
    
    series.append(("Simulation Results", generated_x, generated_y, "#D55E00"))
    
    make_plot(series, output, title)
    print(f"Plot saved to: {output}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation", type=Path, default=Path("force_displacement_VALIDATION.ods"))
    parser.add_argument("--validation-sheet", default="Sheet2")
    parser.add_argument("--generated", type=Path, default=Path("force_displacement.txt"))
    parser.add_argument("--output", type=Path, default=Path("force_displacement_comparison.png"))
    args = parser.parse_args()

    validation = read_ods_columns(args.validation, args.validation_sheet)
    output = args.output if args.output.is_absolute() else ROOT / args.output
    
    plot_one(validation, args.generated, output, "Force vs Displacement Comparison")

if __name__ == "__main__":
    main()

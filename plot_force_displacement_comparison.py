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
    """Read columns from an ODS file using pure python to avoid dependencies."""
    path = Path(path)
    if not path.is_absolute():
        path = ROOT / path
        
    if not path.exists():
        print(f"Warning: Validation file {path} not found.")
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
            values.extend([value] * min(repeated_cols, 100))
        if any(values):
            rows.extend([values] * min(repeated_rows, 100))

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
    """Read simulation output data (supports CSV or space/tab separated)."""
    path = Path(path)
    if not path.is_absolute():
        path = ROOT / path
        
    if not path.exists():
        print(f"Error: Simulation data {path} not found.")
        return [], []

    content = path.read_text(encoding="utf-8", errors="replace")
    lines = content.splitlines()
    if not lines:
        return [], []

    # Detect delimiter
    header = lines[0].lower()
    delimiter = "," if "," in header else None
    
    x, y = [], []
    if delimiter == ",":
        reader = csv.DictReader(lines)
        for row in reader:
            # Handle various header capitalizations
            d_val = row.get("displacement") or row.get("Displacement")
            f_val = row.get("force") or row.get("Force")
            dv, fv = _to_float(d_val), _to_float(f_val)
            if dv is not None and fv is not None:
                x.append(dv)
                y.append(fv)
    else:
        for line in lines:
            parts = line.split()
            if len(parts) < 2:
                continue
            dv, fv = _to_float(parts[0]), _to_float(parts[1])
            if dv is not None and fv is not None:
                x.append(dv)
                y.append(fv)
    
    return x, y

def make_plot(series_list, output_path, title):
    """Generate the Matplotlib plot."""
    plt.figure(figsize=(10, 7), dpi=150)
    
    for label, x, y, color, style in series_list:
        if not x or not y:
            continue
        plt.plot(x, y, label=label, color=color, linestyle=style, 
                 linewidth=2.5, marker='o', markersize=3, alpha=0.8)

    plt.title(title, fontsize=14, fontweight='bold', pad=15)
    plt.xlabel("Displacement (mm)", fontsize=12, labelpad=10)
    plt.ylabel("Force (kN)", fontsize=12, labelpad=10)
    
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.minorticks_on()
    plt.grid(True, which='minor', linestyle=':', alpha=0.2)
    
    # Force axes to start at 0
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    
    plt.legend(loc='best', frameon=True, shadow=True)
    plt.tight_layout()
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()
    print(f"Success: Comparison plot saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation", type=str, default="force_displacement_VALIDATION.ods")
    parser.add_argument("--sheet", default="Sheet2")
    parser.add_argument("--generated", type=str, default="force_displacement.txt")
    parser.add_argument("--output", type=str, default="images/force_displacement_comparison_hetero.png")
    args = parser.parse_args()

    # 1. Read Validation Data
    val_data = read_ods_columns(args.validation, args.sheet)
    val_x = val_data.get("Dis(Paper)", [])
    val_y = val_data.get("Force (paper)", [])

    # 2. Read Simulation Data
    sim_x, sim_y = read_generated(args.generated)

    # 3. Assemble Plot Series
    series = []
    if val_x and val_y:
        series.append(("Reference (Paper)", val_x, val_y, "#0072B2", "--"))
    else:
        print("Warning: No validation data points found to plot.")

    if sim_x and sim_y:
        series.append(("Simulation Results", sim_x, sim_y, "#D55E00", "-"))
    else:
        print("Warning: No simulation data points found to plot.")

    # 4. Create Plot
    if not series:
        print("Error: No data available to plot.")
        return

    make_plot(series, args.output, "Heterogeneous Phase-Field Fracture: Force vs Displacement")

if __name__ == "__main__":
    main()

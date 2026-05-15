#!/usr/bin/env python3
"""Aggregate per-trial force-displacement CSV files into a compact XLSX file."""

from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "output"


def collect_rows() -> list[list[str]]:
    rows = [["pore_count", "trial", "displacement", "force"]]
    for csv_file in sorted(OUTPUT_DIR.glob("**/Pore */trail_*/force_displacement.csv")):
        pore_match = re.search(r"Pore (\d+)", str(csv_file))
        trial_match = re.search(r"trail_(\d+)", str(csv_file))
        pore_count = pore_match.group(1) if pore_match else ""
        trial = trial_match.group(1) if trial_match else ""
        with csv_file.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rows.append(
                    [
                        pore_count,
                        trial,
                        row.get("displacement", ""),
                        row.get("force", ""),
                    ]
                )
    return rows


def cell_ref(row: int, column: int) -> str:
    name = ""
    while column:
        column, remainder = divmod(column - 1, 26)
        name = chr(ord("A") + remainder) + name
    return f"{name}{row}"


def sheet_xml(rows: list[list[str]]) -> str:
    xml_rows = []
    for r_index, row in enumerate(rows, start=1):
        cells = []
        for c_index, value in enumerate(row, start=1):
            ref = cell_ref(r_index, c_index)
            try:
                number = float(value)
                cells.append(f'<c r="{ref}"><v>{number:.16e}</v></c>')
            except ValueError:
                cells.append(
                    f'<c r="{ref}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'
                )
        xml_rows.append(f'<row r="{r_index}">{"".join(cells)}</row>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(xml_rows)}</sheetData></worksheet>'
    )


def write_xlsx(rows: list[list[str]], output: Path) -> None:
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(
            "[Content_Types].xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            "</Types>",
        )
        archive.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            "</Relationships>",
        )
        archive.writestr(
            "xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            '<sheets><sheet name="force_displacement" sheetId="1" r:id="rId1"/></sheets>'
            "</workbook>",
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
            "</Relationships>",
        )
        archive.writestr("xl/worksheets/sheet1.xml", sheet_xml(rows))


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_xlsx(collect_rows(), OUTPUT_DIR / "overall_force.xlsx")

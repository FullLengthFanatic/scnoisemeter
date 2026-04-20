"""
tests/test_sample_sheet.py
==========================
Unit tests for the Smart-seq / FLASH-seq sample-sheet parser.

Covers 96-well (A1..H12) and 384-well (A1..P24) plate geometries across
both input formats: Sample_Name (plate_well combined) and separate
PlateID / WellID columns.
"""

from __future__ import annotations

import csv
import string
from itertools import product
from pathlib import Path

from scnoisemeter.utils.sample_sheet import parse_sample_sheet


def _wells(rows: int, cols: int) -> list[str]:
    letters = list(string.ascii_uppercase[:rows])
    return [f"{r}{c}" for r, c in product(letters, range(1, cols + 1))]


def _write_sample_name_sheet(path: Path, plate_id: str, wells: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Sample_Name", "i7 sequence", "i5 sequence"])
        for i, well in enumerate(wells):
            w.writerow([f"{plate_id}_{well}", f"I7_{i:04d}", f"I5_{i:04d}"])


def _write_plate_well_sheet(path: Path, plate_id: str, wells: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["PlateID", "WellID", "i7_sequence", "i5_sequence"])
        for i, well in enumerate(wells):
            w.writerow([plate_id, well, f"I7_{i:04d}", f"I5_{i:04d}"])


def test_96_well_sample_name_format(tmp_path):
    wells = _wells(8, 12)
    assert len(wells) == 96
    sheet_path = tmp_path / "plate_96.csv"
    _write_sample_name_sheet(sheet_path, plate_id="P96", wells=wells)

    sheet = parse_sample_sheet(sheet_path)

    assert set(sheet.keys()) == {f"P96/{w}" for w in wells}
    # Spot-check 96-well corners.
    for corner in ("A1", "A12", "H1", "H12"):
        assert f"P96/{corner}" in sheet


def test_384_well_sample_name_format(tmp_path):
    wells = _wells(16, 24)
    assert len(wells) == 384
    sheet_path = tmp_path / "plate_384.csv"
    _write_sample_name_sheet(sheet_path, plate_id="P384", wells=wells)

    sheet = parse_sample_sheet(sheet_path)

    assert set(sheet.keys()) == {f"P384/{w}" for w in wells}
    # Corners and positions reachable only on 384-well geometry.
    for corner in ("A1", "A24", "P1", "P24"):
        assert f"P384/{corner}" in sheet
    for well_384_only in ("I13", "I24", "P13"):
        assert f"P384/{well_384_only}" in sheet


def test_384_well_plateid_wellid_format(tmp_path):
    wells = _wells(16, 24)
    sheet_path = tmp_path / "plate_384_format_a.csv"
    _write_plate_well_sheet(sheet_path, plate_id="P384", wells=wells)

    sheet = parse_sample_sheet(sheet_path)

    assert len(sheet) == 384
    assert "P384/A1" in sheet
    assert "P384/P24" in sheet
    entry = sheet["P384/A1"]
    assert entry["well_id"] == "A1"
    assert entry["plate_id"].upper() == "P384"
    assert entry["i7_seq"].startswith("I7_")
    assert entry["i5_seq"].startswith("I5_")


def test_wellid_only_384_has_wellid_keys(tmp_path):
    """Without a PlateID column, canonical keys are bare well IDs."""
    sheet_path = tmp_path / "plate_no_plate_id.csv"
    with open(sheet_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["WellID", "i7_sequence", "i5_sequence"])
        for i, well in enumerate(_wells(16, 24)):
            w.writerow([well, f"I7_{i:04d}", f"I5_{i:04d}"])

    sheet = parse_sample_sheet(sheet_path)

    assert len(sheet) == 384
    assert "A1" in sheet
    assert "P24" in sheet

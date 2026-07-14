"""Backfill Russell (2004) Table 3.4/3.9-3.11 per-row a_au/e/tof_days into catalogue.yaml.

Method validated in docs/notes/2026-07-15-596-russell-backfill-method-validated.md:
tof_days is cited DIRECTLY from the table (matches the existing Aldrin catalogue
precedent); a_au/e are DERIVED (kind: derive) via an AR + V_inf_Earth inversion
against cyclerfinder.search.free_return.free_return_geometry, doubly validated
against Aldrin's independently-sourced elements and cycler 2.5.1.+0's state-vector
ground truth (both < 0.3% error).

Table data transcribed verbatim in
docs/notes/2026-06-07-russell-2004-member-tables-transcription.md (task #142).
Row tuples: (designator, AR, TR, tof_days, vinf_e, vinf_m, table_source).
"""

from __future__ import annotations

import copy
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import yaml
from ruamel.yaml import YAML
from scipy.optimize import least_squares

from cyclerfinder.search.free_return import free_return_geometry

CATALOGUE_PATH = "data/catalogue.yaml"
A_MARS_RUSSELL = 1.52  # Russell's own rounded Mars sma, used to define AR (Ch.3.7)


@dataclass(frozen=True)
class RussellRow:
    designator: str
    ar: float
    tr: float
    tof_days: float
    vinf_e: float
    vinf_m: float
    table: str


# Table 3.4 (p.83) -- 44 rows, 2/3/4-synodic period cyclers
TABLE_3_4 = [
    ("1.0.1.-1", 1.47, 0.86, 146, 6.5, 9.7),
    ("2.1.1.+2", 0.95, 1.11, 207, 4.1, 2.0),
    ("2.3.1.+1", 1.08, 0.92, 143, 5.4, 5.3),
    ("2.5.1.+0", 1.44, 1.12, 94, 7.8, 9.9),
    ("3.1.1.+3", 1.07, 1.19, 174, 3.6, 4.6),
    ("3.1.1.+2", 1.43, 0.89, 115, 5.4, 9.2),
    ("3.1.2.+1", 1.07, 1.23, 181, 3.4, 4.6),
    ("3.1.3.+0", 1.43, 0.93, 123, 5.1, 9.1),
    ("3.3.1.+2", 1.19, 1.06, 141, 4.3, 6.8),
    ("3.5.1.+2", 0.94, 1.80, 231, 2.7, 1.5),
    ("3.5.1.+1", 1.43, 1.15, 115, 5.4, 9.2),
    ("3.5.2.+0", 1.43, 1.06, 121, 5.2, 9.2),
    ("3.7.1.+1", 1.07, 1.56, 175, 3.6, 4.6),
    ("3.9.1.+0", 1.43, 1.17, 116, 5.4, 9.2),
    ("4.0.3.+1", 1.07, 1.18, 160, 4.3, 4.9),
    ("4.1.1.-6", 0.94, 1.37, 256, 2.7, 1.6),
    ("4.1.1.-5", 1.15, 1.11, 173, 4.1, 6.1),
    ("4.1.1.-4", 1.44, 0.89, 137, 5.5, 9.3),
    ("4.1.2.-3", 0.94, 1.40, 250, 2.6, 1.5),
    ("4.1.2.-2", 1.43, 0.93, 132, 5.2, 9.2),
    ("4.1.4.-1", 1.43, 0.93, 129, 5.1, 9.2),
    ("4.3.1.-5", 0.99, 1.29, 268, 3.1, 2.5),
    ("4.3.1.-4", 1.26, 1.01, 154, 4.7, 7.6),
    ("4.5.1.-4", 1.07, 1.55, 196, 3.6, 4.7),
    ("4.5.1.-3", 1.44, 1.15, 137, 5.5, 9.3),
    ("4.5.2.-2", 1.07, 1.40, 191, 3.4, 4.6),
    ("4.5.3.-1", 1.43, 1.02, 130, 5.1, 9.2),
    ("4.6.1.-4", 0.91, 1.50, 154, 6.8, 2.1),
    ("4.6.3.+0", 1.43, 0.88, 105, 6.4, 9.5),
    ("4.7.1.-3", 1.20, 1.38, 163, 4.3, 6.8),
    ("4.7.1.-2", 1.77, 0.96, 120, 6.6, 11.4),
    ("4.8.1.+3", 0.96, 1.64, 164, 7.7, 3.1),
    ("4.8.1.+2", 1.31, 0.86, 76, 12.5, 10.7),
    ("4.9.1.-3", 0.94, 1.83, 256, 2.7, 1.6),
    ("4.9.1.-2", 1.44, 1.16, 137, 5.5, 9.3),
    ("4.9.2.-1", 1.44, 1.05, 132, 5.2, 9.2),
    ("4.10.1.-3", 0.92, 1.46, 263, 10.2, 3.6),
    ("4.10.1.+2", 1.03, 1.65, 131, 8.9, 5.0),
    ("4.11.1.-2", 1.07, 1.58, 195, 3.6, 4.7),
    ("4.12.1.-2", 0.97, 1.43, 268, 11.6, 4.8),
    ("4.12.1.+1", 1.16, 1.48, 93, 10.8, 8.2),
    ("4.13.1.-1", 1.44, 1.16, 137, 5.5, 9.3),
    ("4.14.1.-1", 1.12, 1.13, 199, 14.7, 9.4),
    ("4.14.1.+0", 1.49, 1.09, 66, 14.1, 12.7),
]

# Table 3.9 (p.90) -- 58 rows, 5-synodic period cyclers
TABLE_3_9 = [
    ("5.1.1.-7", 1.04, 0.97, 229, 5.0, 4.3),
    ("5.1.2.-3", 1.20, 1.00, 168, 4.7, 7.0),
    ("5.1.5.-1", 1.44, 0.92, 133, 5.2, 9.2),
    ("5.2.1.-7", 0.90, 1.07, 182, 4.5, 1.3),
    ("5.2.2.+2", 1.20, 0.94, 128, 5.2, 7.1),
    ("5.2.5.+0", 1.43, 0.91, 118, 5.3, 9.2),
    ("5.3.1.-7", 0.92, 1.17, 270, 3.8, 1.4),
    ("5.3.1.-6", 1.10, 0.90, 205, 5.5, 5.7),
    ("5.3.3.-2", 1.07, 1.19, 195, 3.6, 4.7),
    ("5.4.1.+6", 0.94, 1.45, 189, 4.9, 1.9),
    ("5.4.1.+5", 1.12, 1.06, 122, 7.0, 6.3),
    ("5.4.3.+1", 1.07, 1.45, 170, 3.8, 4.7),
    ("5.5.1.-6", 0.96, 1.95, 279, 4.3, 2.2),
    ("5.5.1.-5", 1.18, 1.44, 186, 6.2, 7.0),
    ("5.5.1.-4", 1.48, 1.08, 154, 8.0, 10.3),
    ("5.5.2.-3", 0.94, 1.79, 262, 3.0, 1.7),
    ("5.5.2.-2", 1.45, 1.19, 142, 5.9, 9.5),
    ("5.5.4.-1", 1.44, 1.10, 134, 5.3, 9.3),
    ("5.6.1.+5", 0.98, 1.74, 198, 5.4, 2.7),
    ("5.6.1.+4", 1.20, 1.23, 107, 7.7, 7.6),
    ("5.6.1.+3", 1.49, 0.90, 82, 9.8, 11.0),
    ("5.6.2.+2", 0.94, 1.36, 219, 3.3, 1.7),
    ("5.6.2.+1", 1.44, 0.87, 104, 6.4, 9.5),
    ("5.6.4.+0", 1.43, 1.16, 116, 5.4, 9.2),
    ("5.7.1.-5", 1.02, 1.76, 245, 4.8, 3.6),
    ("5.7.1.-4", 1.30, 1.27, 169, 7.0, 8.5),
    ("5.7.1.-3", 1.71, 0.93, 142, 9.1, 11.9),
    ("5.8.1.+4", 1.03, 1.91, 154, 6.1, 4.3),
    ("5.8.1.+3", 1.31, 1.30, 94, 8.6, 9.1),
    ("5.8.1.+2", 1.72, 0.92, 73, 11.0, 12.6),
    ("5.9.1.-4", 1.10, 2.09, 204, 5.6, 5.7),
    ("5.9.1.-3", 1.48, 1.38, 154, 8.0, 10.3),
    ("5.9.1.-2", 2.15, 0.93, 130, 10.5, 13.8),
    ("5.9.2.-2", 1.08, 1.57, 198, 4.0, 4.9),
    ("5.9.2.+1", 2.10, 0.90, 117, 7.9, 12.9),
    ("5.9.3.-1", 1.44, 1.16, 137, 5.5, 9.3),
    ("5.10.1.+3", 1.11, 1.94, 123, 6.9, 6.2),
    ("5.10.1.+2", 1.49, 1.25, 82, 9.8, 10.9),
    ("5.10.2.+1", 1.07, 1.71, 160, 4.3, 4.9),
    ("5.10.2.+0", 2.08, 0.87, 81, 8.6, 13.0),
    ("5.10.3.+0", 1.43, 1.24, 112, 5.7, 9.3),
    ("5.11.1.-3", 1.24, 1.76, 177, 6.6, 7.8),
    ("5.11.1.-2", 1.83, 1.09, 138, 9.5, 12.3),
    ("5.12.1.+2", 1.14, 1.00, 101, 9.6, 7.5),
    ("5.12.1.+1a", 1.24, 1.81, 101, 8.1, 8.3),  # duplicated index in source PDF (verbatim)
    ("5.12.1.+1b", 1.82, 1.07, 70, 11.5, 13.2),  # duplicated index in source PDF (verbatim)
    ("5.13.1.-3", 0.97, 2.69, 280, 4.3, 2.3),
    ("5.13.1.-2", 1.49, 1.39, 153, 8.1, 10.3),
    ("5.13.2.-1", 1.45, 1.20, 141, 5.9, 9.5),
    ("5.14.1.+2", 0.97, 3.06, 196, 5.3, 2.6),
    ("5.14.1.+1", 1.48, 1.41, 82, 9.8, 10.9),
    ("5.14.2.+0", 1.43, 1.22, 105, 6.4, 9.7),
    ("5.15.1.-2", 1.11, 2.13, 202, 5.6, 5.8),
    ("5.15.1.-1", 2.16, 0.93, 130, 10.5, 13.9),
    ("5.16.1.+1", 1.10, 2.37, 126, 6.8, 6.0),
    ("5.16.1.+0", 2.12, 0.91, 64, 12.6, 14.6),
    ("5.17.1.-1", 1.50, 1.37, 152, 8.1, 10.4),
    ("5.18.1.+0", 1.46, 1.45, 84, 9.6, 10.6),
]

# Table 3.10 (p.91) -- 51 rows, 6-synodic period cyclers, Part I
TABLE_3_10 = [
    ("6.0.1.+9", 0.92, 1.40, 213, 3.0, 1.2),
    ("6.0.1.+8", 1.03, 1.22, 179, 4.0, 3.9),
    ("6.0.1.+7", 1.17, 1.07, 133, 5.0, 6.7),
    ("6.0.1.+6", 1.34, 0.93, 111, 6.0, 8.7),
    ("6.1.2.-4", 1.09, 0.91, 203, 4.9, 5.4),
    ("6.1.3.-3", 0.95, 1.30, 264, 3.1, 1.7),
    ("6.1.4.-2", 1.07, 1.16, 197, 3.8, 4.8),
    ("6.1.6.-1", 1.44, 0.90, 135, 5.4, 9.3),
    ("6.2.1.+8", 0.94, 1.26, 220, 3.3, 1.7),
    ("6.2.1.+7", 1.08, 1.07, 158, 4.3, 5.0),
    ("6.2.1.+6", 1.24, 0.91, 123, 5.4, 7.6),
    ("6.2.2.+3", 1.07, 1.19, 174, 3.6, 4.6),
    ("6.2.2.+2", 1.43, 0.89, 115, 5.4, 9.2),
    ("6.2.3.+2", 0.94, 1.39, 235, 2.6, 1.5),
    ("6.2.3.+1", 1.43, 0.92, 119, 5.2, 9.2),
    ("6.2.4.+1", 1.07, 1.23, 181, 3.4, 4.6),
    ("6.2.6.+0", 1.43, 0.93, 123, 5.1, 9.1),
    ("6.3.1.-9", 0.92, 0.89, 279, 5.7, 1.8),
    ("6.3.4.+1", 1.07, 1.04, 156, 4.5, 5.0),
    ("6.4.1.+7", 0.98, 1.59, 227, 3.6, 2.4),
    ("6.4.1.+6", 1.13, 1.33, 142, 4.7, 6.0),
    ("6.4.1.+5", 1.33, 1.11, 113, 5.9, 8.5),
    ("6.4.1.+4", 1.58, 0.93, 96, 7.0, 10.6),
    ("6.5.1.-8", 0.95, 1.62, 283, 6.2, 2.4),
    ("6.5.1.-7", 1.11, 1.16, 213, 8.4, 6.6),
    ("6.5.1.-6", 1.31, 0.88, 180, 10.4, 9.8),
    ("6.5.5.-1", 1.44, 1.15, 137, 5.5, 9.3),
    ("6.6.1.-6", 1.02, 1.82, 189, 3.9, 3.5),
    ("6.6.1.+5", 1.20, 1.48, 128, 5.2, 7.1),
    ("6.6.1.+4", 1.45, 1.21, 104, 6.4, 9.6),
    ("6.6.1.+3", 1.78, 0.99, 89, 7.7, 11.7),
    ("6.6.2.+2", 1.19, 1.38, 141, 4.3, 6.8),
    ("6.6.2.+1", 1.77, 0.96, 99, 6.6, 11.3),
    ("6.6.5.+0", 1.43, 1.03, 122, 5.1, 9.2),
    ("6.7.1.-7", 0.98, 1.49, 289, 6.7, 3.1),
    ("6.7.1.-6", 1.17, 1.05, 199, 9.1, 7.8),
    ("6.7.2.+3", 0.91, 0.98, 176, 5.1, 1.5),
    ("6.7.3.-2", 1.08, 1.40, 199, 4.1, 5.0),
    ("6.7.5.+0", 1.43, 0.99, 107, 6.1, 9.4),
    ("6.8.1.+6", 0.91, 2.39, 211, 2.9, 1.0),
    ("6.8.1.+5", 1.08, 1.89, 158, 4.3, 5.0),
    ("6.8.1.+4", 1.30, 1.50, 116, 5.7, 8.3),
    ("6.8.1.+3", 1.62, 1.18, 95, 7.2, 10.8),
    ("6.8.1.+2", 2.09, 0.94, 81, 8.6, 13.1),
    ("6.8.3.+1", 1.07, 1.47, 179, 3.4, 4.6),
    ("6.9.1.-6", 1.03, 1.97, 248, 7.3, 4.5),
    ("6.9.1.-5", 1.26, 1.35, 186, 10.0, 9.1),
    ("6.9.1.-4", 1.57, 0.98, 161, 12.4, 12.4),
    ("6.9.2.-3", 0.96, 1.61, 274, 3.8, 2.0),
    ("6.9.2.-2", 1.47, 0.94, 150, 7.2, 10.0),
    ("6.9.4.-1", 1.45, 1.21, 139, 5.7, 9.4),
]

# Table 3.11 (p.92) -- 48 rows, 6-synodic period cyclers, Part II
TABLE_3_11 = [
    ("6.10.1.+5", 0.94, 2.31, 219, 3.3, 1.7),
    ("6.10.1.+4", 1.15, 1.76, 137, 4.9, 6.4),
    ("6.10.1.+3", 1.44, 1.34, 104, 6.4, 9.6),
    ("6.10.1.+2", 1.89, 1.02, 86, 8.1, 12.2),
    ("6.10.2.+2", 0.94, 1.84, 231, 2.7, 1.5),
    ("6.10.2.+1", 1.43, 1.17, 115, 5.4, 9.2),
    ("6.10.4.+0", 1.43, 1.06, 121, 5.2, 9.2),
    ("6.11.1.-5", 1.09, 1.75, 219, 8.2, 6.2),
    ("6.11.1.-4", 1.38, 1.17, 173, 11.0, 10.7),
    ("6.11.2.+2", 0.98, 1.61, 191, 6.1, 2.9),
    ("6.12.1.+4", 1.00, 2.14, 232, 3.7, 2.7),
    ("6.12.1.+3", 1.26, 1.57, 120, 5.5, 7.8),
    ("6.12.1.+2", 1.67, 1.15, 92, 7.4, 11.1),
    ("6.13.1.+5", 1.08, 0.94, 79, 16.7, 9.7),
    ("6.13.1.-5", 0.90, 3.44, 276, 5.4, 1.5),
    ("6.13.1.-4", 1.18, 1.92, 198, 9.2, 7.9),
    ("6.13.1.-3", 1.58, 1.21, 160, 12.4, 12.5),
    ("6.13.2.-2", 1.10, 1.65, 202, 5.0, 5.4),
    ("6.13.3.-1", 1.46, 1.06, 143, 6.1, 9.6),
    ("6.14.1.+3", 1.08, 1.93, 158, 4.3, 4.9),
    ("6.14.1.+2", 1.44, 1.34, 104, 6.4, 9.6),
    ("6.14.1.+1", 2.09, 0.94, 81, 8.6, 13.1),
    ("6.14.2.+1", 1.07, 1.59, 175, 3.6, 4.6),
    ("6.14.2.+0", 2.08, 0.85, 91, 7.4, 12.7),
    ("6.14.3.+0", 1.43, 1.10, 119, 5.2, 9.2),
    ("6.15.1.+4", 1.11, 1.90, 74, 17.2, 10.4),
    ("6.15.1.-4", 0.95, 3.00, 285, 6.3, 2.6),
    ("6.15.1.-3", 1.32, 1.57, 179, 10.5, 10.0),
    ("6.15.1.-2", 1.94, 0.94, 147, 14.3, 14.8),
    ("6.15.2.+1", 1.11, 1.27, 117, 7.7, 6.4),
    ("6.16.1.+2", 1.20, 1.67, 128, 5.2, 7.1),
    ("6.16.1.+1", 1.78, 1.08, 89, 7.7, 11.7),
    ("6.17.1.+3", 1.14, 1.06, 69, 17.9, 11.3),
    ("6.17.1.-3", 1.04, 2.91, 241, 7.5, 4.8),
    ("6.17.1.-2", 1.59, 1.24, 160, 12.5, 12.6),
    ("6.17.2.-1", 1.48, 1.30, 149, 7.3, 10.1),
    ("6.18.1.+2", 0.94, 2.32, 219, 3.3, 1.7),
    ("6.18.1.+1", 1.44, 1.34, 104, 6.4, 9.5),
    ("6.18.2.+0", 1.43, 1.17, 116, 5.4, 9.2),
    ("6.19.1.+2", 1.19, 0.97, 64, 18.8, 12.5),
    ("6.19.1.-2", 1.20, 2.09, 195, 9.4, 8.2),
    ("6.19.2.+0", 1.47, 1.18, 79, 10.6, 11.1),
    ("6.20.1.-4", 0.93, 1.09, 183, 12.8, 5.0),
    ("6.20.1.+1", 1.07, 1.94, 160, 4.3, 4.9),
    ("6.20.1.+0", 2.08, 0.94, 81, 8.6, 13.0),
    ("6.21.1.+1", 1.29, 1.01, 57, 20.3, 14.4),
    ("6.21.1.-1", 1.63, 1.19, 158, 12.7, 12.9),
    ("6.22.1.+0", 1.43, 1.35, 105, 6.4, 9.5),
]


def build_rows() -> list[RussellRow]:
    rows = []
    for table_name, table in (
        ("3.4", TABLE_3_4),
        ("3.9", TABLE_3_9),
        ("3.10", TABLE_3_10),
        ("3.11", TABLE_3_11),
    ):
        for designator, ar, tr, tof, vinf_e, vinf_m in table:
            rows.append(RussellRow(designator, ar, tr, float(tof), vinf_e, vinf_m, table_name))
    return rows


def normalize_designator(raw: str) -> str | None:
    m = re.search(r"(\d+)\.(\d+)\.(\d+)\.?([+-]\d+)", raw)
    if m:
        return f"{m.group(1)}.{m.group(2)}.{m.group(3)}.{m.group(4)}"
    return None


def invert_a_e(ar: float, vinf_e_target: float) -> tuple[float, float, float]:
    """AR + V_inf_Earth -> (a_au, e), validated method (see module docstring).

    Returns (a_au, e, cost) from the best of a coarse multi-start global search.
    """
    aphelion_target = ar * A_MARS_RUSSELL

    def resid(x):
        a_au, e = x
        if not (0.0 < e < 0.95) or a_au <= 0.0:
            return [1e3, 1e3]
        aphelion = a_au * (1.0 + e)
        try:
            g = free_return_geometry(a_au, e)
        except ValueError:
            return [1e3, 1e3]
        return [aphelion - aphelion_target, g.vinf["E"] - vinf_e_target]

    best = None
    for a0 in (1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.4, 2.7):
        for e0 in (0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75):
            sol = least_squares(resid, x0=[a0, e0], bounds=([1.0, 0.01], [3.5, 0.94]))
            if best is None or sol.cost < best[0]:
                best = (sol.cost, sol.x)
    cost, x = best
    return float(x[0]), float(x[1]), float(cost)


def h_value(designator: str) -> int:
    m = re.match(r"(\d+)\.(\d+)\.(\d+)\.?([+-]\d+)", designator)
    return int(m.group(2))


def analyze(write: bool = False) -> None:
    rows = build_rows()
    print(f"Parsed {len(rows)} Russell Table 3.4/3.9-3.11 rows.")

    with open(CATALOGUE_PATH) as f:
        catalogue = yaml.safe_load(f)

    cat_by_designator: dict[str, dict] = {}
    for entry in catalogue:
        eid = entry.get("id", "")
        if "russell" not in eid.lower() and "ocampo" not in eid.lower():
            continue
        norm = normalize_designator(eid)
        if norm:
            cat_by_designator.setdefault(norm, entry)

    results = []
    excluded = []
    flagged = []
    cost_gate = 1e-4
    for row in rows:
        norm = normalize_designator(row.designator)
        entry = cat_by_designator.get(norm)
        if entry is None:
            continue
        a_au, e, cost = invert_a_e(row.ar, row.vinf_e)
        h = h_value(row.designator)
        if cost >= cost_gate:
            # AR<1 near-ballistic rows: the ellipse's aphelion doesn't reach Mars,
            # so this crossing-based inversion doesn't apply -- a genuine model
            # boundary (confirmed: this correlates 1:1 with row.ar < 1.0 in the
            # full Table 3.4/3.9-3.11 set), NOT a solver bug. Leave as an open gap.
            excluded.append(
                {
                    "designator": row.designator,
                    "entry_id": entry.get("id"),
                    "ar": row.ar,
                    "cost": cost,
                }
            )
            continue
        try:
            g = free_return_geometry(a_au, e)
        except ValueError:
            excluded.append(
                {
                    "designator": row.designator,
                    "entry_id": entry.get("id"),
                    "ar": row.ar,
                    "cost": cost,
                }
            )
            continue
        vinf_m_err_pct = 100.0 * abs(g.vinf["M"] - row.vinf_m) / row.vinf_m if row.vinf_m else None
        tof_err_pct = 100.0 * abs(g.tof_em_days - row.tof_days) / row.tof_days
        rec = {
            "designator": row.designator,
            "entry_id": entry.get("id"),
            "table": row.table,
            "a_au": a_au,
            "e": e,
            "cost": cost,
            "vinf_m_emerged": g.vinf["M"],
            "vinf_m_target": row.vinf_m,
            "vinf_m_err_pct": vinf_m_err_pct,
            "tof_em_emerged": g.tof_em_days,
            "tof_target": row.tof_days,
            "tof_err_pct": tof_err_pct,
            "h": h,
        }
        results.append(rec)
        # Flag anything where the free (non-imposed) V_inf-at-Mars cross-check
        # is not tight -- these need eyeballing before any writeback, not blind trust.
        if vinf_m_err_pct is not None and vinf_m_err_pct > 5.0:
            flagged.append(rec)

    print(f"Rows EXCLUDED (cost >= {cost_gate}, model doesn't apply -- left open): {len(excluded)}")
    ar_of_excluded = [e["ar"] for e in excluded]
    if ar_of_excluded:
        all_below_1 = all(a < 1.0 for a in ar_of_excluded)
        print(f"  all AR < 1.0? {all_below_1} (max AR excluded: {max(ar_of_excluded):.2f})")
    print(f"Rows with a computed (a,e), accepted for writeback: {len(results)}")
    print(f"Flagged (V_inf_Mars cross-check error > 5%): {len(flagged)}")
    for rec in flagged:
        print(
            f"  {rec['designator']} ({rec['entry_id']}): "
            f"a={rec['a_au']:.4f} e={rec['e']:.4f} "
            f"vinf_M emerged={rec['vinf_m_emerged']:.3f} target={rec['vinf_m_target']:.3f} "
            f"err={rec['vinf_m_err_pct']:.1f}%"
        )

    h_zero = [r for r in results if r["h"] == 0]
    h_nonzero = [r for r in results if r["h"] != 0]
    print(f"h=0 rows (safe to backfill out-em + ret-me): {len(h_zero)}")
    print(f"h>0 rows (safe to backfill out-em + orbit_elements ONLY): {len(h_nonzero)}")

    import statistics

    errs = [r["vinf_m_err_pct"] for r in results if r["vinf_m_err_pct"] is not None]
    print(
        f"V_inf_Mars cross-check error: median={statistics.median(errs):.2f}% "
        f"mean={statistics.mean(errs):.2f}% max={max(errs):.2f}%"
    )

    if write:
        write_via_patch(results)


def apply_writeback(catalogue: list[dict], cat_by_designator: dict, results: list[dict]) -> None:
    by_id = {e.get("id"): e for e in catalogue}
    for rec in results:
        entry = by_id.get(rec["entry_id"])
        if entry is None:
            continue
        a_au, e = round(rec["a_au"], 4), round(rec["e"], 4)
        citation = (
            f"DERIVED (kind: derive) from Russell 2004 Table {rec.get('table', '3.4/3.9-3.11')} "
            f"row {rec['designator']}'s sourced Aphelion Ratio + V_inf_Earth, via the "
            "AR+V_inf_Earth -> (a_au, e) inversion validated in "
            "docs/notes/2026-07-15-596-russell-backfill-method-validated.md "
            f"(V_inf_Mars free cross-check: emerged {rec['vinf_m_emerged']:.3f} vs. "
            f"tabulated {rec['vinf_m_target']:.3f} km/s, {rec['vinf_m_err_pct']:.1f}% error)."
        )

        # NOTE: top-level orbit_elements.a_au/e is schema-restricted to
        # cycler_class == 'single-ellipse' ("one ellipse for the whole
        # cycler"). Every Russell-family entry here is 'multi-arc' (out-em/
        # ret-me/loop-ee-* are potentially DIFFERENT ellipses), so writing
        # orbit_elements.a_au/e would be a genuine schema violation, not a
        # style choice -- confirmed by test_validate_catalogue.py's
        # combined-gate failure the first time this was tried. Only
        # per-segment a_au/e are written below.
        wrote_out_em = False

        segs = (entry.get("trajectory") or {}).get("segments") or []
        seg_by_id = {s.get("id"): s for s in segs}

        out_em = seg_by_id.get("out-em")
        if out_em is not None and out_em.get("a_au") is None:
            out_em["a_au"] = a_au
            out_em["e"] = e
            out_em["note"] = (out_em.get("note") or "").rstrip() + " " + citation
            wrote_out_em = True

        wrote_ret_me = False
        if rec["h"] == 0:
            ret_me = seg_by_id.get("ret-me")
            if ret_me is not None and ret_me.get("a_au") is None:
                ret_me["a_au"] = a_au
                ret_me["e"] = e
                ret_me["note"] = (ret_me.get("note") or "").rstrip() + " " + citation
                wrote_ret_me = True

        if wrote_out_em or wrote_ret_me:
            sq = entry.setdefault("source_quotes", {})
            if wrote_out_em:
                sq["trajectory.segments[out-em].a_au"] = citation
                sq["trajectory.segments[out-em].e"] = citation
            if wrote_ret_me:
                sq["trajectory.segments[ret-me].a_au"] = citation
                sq["trajectory.segments[ret-me].e"] = citation

        # Narrow/remove data_gaps entries this writeback actually closes.
        # "orbit_elements.a_au" is NEVER touched here -- it stays null (that's
        # schema-correct for cycler_class=multi-arc, see the note above) and
        # any pre-existing gap entry for it is left exactly as-is.
        # "trajectory.segments[out-em].a_au" (scoped to out-em ONLY) is fully
        # closed by this writeback regardless of h. "trajectory.segments[*]
        # .a_au" (wildcard covering every segment) is only fully closed when
        # h==0 (both out-em and ret-me resolved); otherwise narrowed to note
        # that out-em specifically is now done.
        gaps = entry.get("data_gaps") or []
        new_gaps = []
        for gap in gaps:
            path = gap.get("path", "")
            if path == "trajectory.segments[out-em].a_au" and wrote_out_em:
                continue  # fully closed by this writeback
            if path == "trajectory.segments[*].a_au" and wrote_out_em:
                if rec["h"] == 0 and wrote_ret_me:
                    continue  # both out-em and ret-me now resolved -- fully closed
                gap = dict(gap)
                gap["note"] = (
                    (gap.get("note") or "")
                    + " [2026-07-15: out-em a_au/e now DERIVED via the validated "
                    "AR+V_inf_Earth inversion (docs/notes/2026-07-15-596-russell-"
                    "backfill-method-validated.md); remaining segments (ret-me / loop-ee-*) "
                    "still need the multi-rev Lambert rev-count resolution (task #54) before "
                    "their own a_au/e can be safely assigned -- NOT simply copied from "
                    "out-em, since they may involve additional revolutions.]"
                )
            new_gaps.append(gap)
        entry["data_gaps"] = new_gaps


def _make_ruamel() -> YAML:
    y = YAML(typ="rt")
    y.preserve_quotes = True
    y.width = 4096
    return y


def write_via_patch(results: list[dict]) -> None:
    """Apply the backfill to data/catalogue.yaml with minimal diff noise.

    A plain PyYAML load+dump (or even a ruamel round-trip of the WHOLE file)
    reformats/drifts far more of the file than this backfill actually touches
    (see docs/notes/2026-07-15-596-russell-backfill-method-validated.md's
    companion investigation -- confirmed empirically before writing this
    function). Strategy: ruamel-round-trip the file twice (once untouched,
    once with the backfill applied), diff those two round-trips against each
    other (isolating ONLY the actual backfill's changes from any round-trip
    reformatting noise), then apply that diff as a patch onto the pristine
    original file. Round-trip reformatting noise never touches the original.
    """
    yaml_rt = _make_ruamel()
    path = Path(CATALOGUE_PATH)
    original_text = path.read_text(encoding="utf-8")

    with open(CATALOGUE_PATH) as f:
        baseline = yaml_rt.load(f)
    modified = copy.deepcopy(baseline)

    cat_by_designator = {}
    for entry in modified:
        eid = entry.get("id", "")
        if "russell" not in eid.lower() and "ocampo" not in eid.lower():
            continue
        norm = normalize_designator(eid)
        if norm:
            cat_by_designator[norm] = entry

    apply_writeback(modified, cat_by_designator, results)

    with tempfile.TemporaryDirectory() as tmpdir:
        baseline_path = Path(tmpdir) / "baseline.yaml"
        modified_path = Path(tmpdir) / "modified.yaml"
        with open(baseline_path, "w") as f:
            yaml_rt.dump(baseline, f)
        with open(modified_path, "w") as f:
            yaml_rt.dump(modified, f)

        diff = subprocess.run(
            ["diff", "-u", str(baseline_path), str(modified_path)],
            capture_output=True,
            text=True,
        )
        patch_text = diff.stdout
        n_diff_lines = len(patch_text.splitlines())
        print(f"Isolated backfill diff: {n_diff_lines} lines (vs. round-trip baseline)")
        if not patch_text.strip():
            print("Nothing to write (empty diff).")
            return

        # Rewrite the patch's file headers to target the real file, then apply.
        patch_lines = patch_text.splitlines(keepends=True)
        patch_lines[0] = f"--- a/{CATALOGUE_PATH}\n"
        patch_lines[1] = f"+++ b/{CATALOGUE_PATH}\n"
        patch_path = Path(tmpdir) / "backfill.patch"
        patch_path.write_text("".join(patch_lines), encoding="utf-8")

        result = subprocess.run(
            ["patch", "--fuzz=3", str(path), str(patch_path)],
            capture_output=True,
            text=True,
        )
        print(result.stdout)
        if result.returncode != 0:
            print(result.stderr)
            raise SystemExit(f"patch failed (exit {result.returncode}) -- original untouched")

    new_text = path.read_text(encoding="utf-8")
    n_before, n_after = len(original_text.splitlines()), len(new_text.splitlines())
    print(f"Original file: {n_before} lines -> new file: {n_after} lines")


if __name__ == "__main__":
    import sys

    analyze(write="--write" in sys.argv)

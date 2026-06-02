"""Data-integrity (V0) gate for the Hollister-Menning 1970 15-orbit E-V family.

These entries are NOT yet V1-reproducible (they are 26-encounter E-V orbits,
the same multi-encounter reconstruction blocker as the SnLm rows), so this is
a *traceability* test, not an optimiser rediscovery: it asserts every
catalogue entry's V_inf, period, and family linkage trace exactly to the
auditable transcription of the primary Table 3
(``data/sources/hollister-menning-1970-table3.yaml``). It guards against the
catalogue drifting from its sourced numbers.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

_REPO = Path(__file__).resolve().parent.parent.parent
_CATALOGUE = _REPO / "data" / "catalogue.yaml"
_TABLE3 = _REPO / "data" / "sources" / "hollister-menning-1970-table3.yaml"

# Earth Mean Orbital Speed (km/s) — the EMOS unit Table 3 tabulates V_r in.
_EMOS_KMS = 29.785


def _orbits() -> list[dict[str, Any]]:
    d = yaml.safe_load(_TABLE3.read_text())
    return d if isinstance(d, list) else d.get("orbits", d)


def _entries() -> dict[str, dict[str, Any]]:
    return {r["id"]: r for r in yaml.safe_load(_CATALOGUE.read_text())}


def test_family_has_15_members() -> None:
    rows = _entries()
    ids = [f"hollister-menning-1970-ev-orbit-{n:02d}" for n in range(1, 16)]
    for cid in ids:
        assert cid in rows, f"missing family entry {cid}"
        fam = rows[cid].get("family") or {}
        assert fam.get("id") == "hollister-menning-1970-ev"
        assert fam.get("members") == 15


@pytest.mark.parametrize("orbit_n", range(1, 16))
def test_entry_vinf_traces_to_table3(orbit_n: int) -> None:
    """Each entry's V_inf == round(max per-body V_r * EMOS, 2) from Table 3."""
    orbit = next(o for o in _orbits() if o["orbit"] == orbit_n)
    encs = orbit["encounters"]
    max_e = max(e["vr_emos"] for e in encs if e["planet"] == "E")
    max_v = max(e["vr_emos"] for e in encs if e["planet"] == "V")
    expected = {"E": round(max_e * _EMOS_KMS, 2), "V": round(max_v * _EMOS_KMS, 2)}

    entry = _entries()[f"hollister-menning-1970-ev-orbit-{orbit_n:02d}"]
    got = {v["body"]: v["vinf_kms"] for v in entry["vinf_kms_at_encounters"]}
    assert got["E"] == pytest.approx(expected["E"], abs=0.01), f"orbit {orbit_n} E V_inf"
    assert got["V"] == pytest.approx(expected["V"], abs=0.01), f"orbit {orbit_n} V V_inf"

    period = entry["period"]
    assert period["k"] == 10
    assert period["years"] == pytest.approx(16.0)


def test_table3_date_checksums() -> None:
    """Earth-row date span repeats after 5844 d (16 yr) for orbits 1-14; orbit
    15 prints 5843 — a documented 1-day inconsistency in the source itself."""
    for o in _orbits():
        earth_dates = [e["date_j2440000"] for e in o["encounters"] if e["planet"] == "E"]
        span = earth_dates[-1] - earth_dates[0]
        expected = 5843 if o["orbit"] == 15 else 5844
        assert span == expected, f"orbit {o['orbit']} date span {span} != {expected}"

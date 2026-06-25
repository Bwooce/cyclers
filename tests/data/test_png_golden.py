"""#450 Task 1: sourced Png' hybrid-family golden (arXiv:2509.12671).

The golden EXPECTED side traces ONLY to arXiv:2509.12671's printed tables (Table 3
first planar case, Table 5 second planar case), transcribed in
``docs/notes/2026-06-13-high-order-transfer-map-2509.12671-mining.md``. Never a
value our own enumerator/corrector computed (feedback_golden_tests_sourced_only).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

_GOLDEN = Path(__file__).resolve().parents[2] / "data" / "golden" / "png_hybrid_family.yaml"

# Transcribed verbatim from arXiv:2509.12671 Table 3 (first planar case) and
# Table 5 (P3g', second planar case). order: (x0, xdot0, ydot0, period, n).
_EXPECTED_MEMBERS: dict[str, tuple[float, float, float, float, int]] = {
    "P5g'": (0.807357887647950, -0.0956081545978604, 0.433518861583397, 11.1751086919436, 5),
    "P7g'-I": (0.916929181700578, -0.175717632213615, 0.528042521610021, 11.6960585356669, 7),
    "P7g'-II": (0.843301779064331, 0.0264518382948230, 0.433467279267690, 11.6716427018773, 7),
    "P7g'": (0.800533327625822, -0.0783668113509505, 0.441933007711387, 15.0692523706070, 7),
    "P9g'": (0.807337935300132, -0.0956506138795539, 0.433522872600990, 20.9914771396290, 9),
    "P3g'": (0.852098052983502, -0.187721536949396, 0.368541113320107, 9.25206400352203, 3),
}
_FIELDS = ("x0", "xdot0", "ydot0", "period", "n")


def _load() -> dict[str, Any]:
    data: dict[str, Any] = yaml.safe_load(_GOLDEN.read_text())
    return data


def test_golden_exists_and_has_provenance() -> None:
    data = _load()
    assert data["source"] == "arXiv:2509.12671"
    assert "members" in data


def test_golden_has_six_named_members() -> None:
    data = _load()
    names = {m["label"] for m in data["members"]}
    assert names == set(_EXPECTED_MEMBERS), names


def test_each_member_carries_required_fields() -> None:
    data = _load()
    for m in data["members"]:
        for field in ("label", "x0", "xdot0", "ydot0", "period", "jacobi", "n"):
            assert field in m, (m.get("label"), field)


def test_member_values_match_paper_transcription() -> None:
    data = _load()
    by_label = {m["label"]: m for m in data["members"]}
    for label, exp in _EXPECTED_MEMBERS.items():
        m = by_label[label]
        for k, v in zip(_FIELDS, exp, strict=True):
            assert float(m[k]) == float(v), (label, k, m[k], v)


def test_jacobi_constants_match_paper_planar_cases() -> None:
    data = _load()
    by_label = {m["label"]: m for m in data["members"]}
    # First planar case C_J = 3.00022 (Table 3); P3g' from second case 3.020052 (Table 5).
    for label in ("P5g'", "P7g'-I", "P7g'-II", "P7g'", "P9g'"):
        assert float(by_label[label]["jacobi"]) == 3.00022
    assert float(by_label["P3g'"]["jacobi"]) == 3.020052

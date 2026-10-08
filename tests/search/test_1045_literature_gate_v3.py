"""#1045: literature_check gate v3 fixes G1-G5 -- pinned probes (fast, offline).

Pre-registration: docs/notes/2026-10-08-1045-literature-gate-v3-preregistration.md sec. 1. Each
probe is the case the #972 v2 re-review found (v2 note sec. 9.1); the expected side is the
pre-registered rule.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

from cyclerfinder.search import literature_check as lc
from cyclerfinder.search.literature_check import (
    TOUR_ONLY_CAP,
    CandidateSignature,
    check_literature,
    offline_corpus_search,
    signature_from_review_entry,
)

RM = frozenset({"repeated-moon"})


def _ratchet() -> Any:
    """The F13 ratchet's scanner from the #972 pinned tests (one implementation)."""
    path = Path(__file__).with_name("test_972_literature_gate_v2.py")
    spec = importlib.util.spec_from_file_location("_t972", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod._status_compare_violations_in


def _sig(primary: str, seq: tuple[str, ...], **kw: Any) -> CandidateSignature:
    base: dict[str, Any] = {
        "primary": primary,
        "sequence": seq,
        "period_k": None,
        "vinf_per_encounter_kms": (),
    }
    base.update(kw)
    return CandidateSignature(**base)


def _check(sig: CandidateSignature) -> lc.LiteratureCheckResult:
    return check_literature(sig, search=offline_corpus_search)


# --- G1: a tour-only anchor cannot make an untopologied signature "published" -----------------


def test_g1_untopologied_ev_not_published_by_hughes() -> None:
    for seq in (("E", "V"), ("E", "V", "E"), ("V", "V", "V", "E")):
        r = _check(_sig("Sun", seq))
        assert r.status != "published", (seq, r)
        assert r.confidence <= TOUR_ONLY_CAP, (seq, r)


def test_g1_topology_declared_tour_still_matches_tour_anchor() -> None:
    r = _check(
        _sig(
            "Uranus",
            ("Umbriel", "Oberon", "Umbriel"),
            topology_label=frozenset({"mga-tour"}),
        )
    )
    assert r.status == "published", r


# --- G2: F7 on body-set inclusion ----------------------------------------------------------------


def test_g2_single_moon_subset_goes_to_a_human() -> None:
    titan = _check(_sig("Saturn", ("Titan",), topology_label=RM))
    assert titan.status == "inconclusive", titan
    assert "Russell-Strange 2009 Titan-Enceladus" in titan.notes, titan
    gan = _check(_sig("Jupiter", ("Ganymede", "Ganymede"), topology_label=RM))
    assert gan.status == "inconclusive", gan


# --- G3: a missing review-entry primary is never guessed --------------------------------------


class _Entry:
    def __init__(self, audit: dict[str, Any]) -> None:
        self.sequence = ("Ganymede", "Callisto", "Ganymede")
        self.period_k = 3
        self.vinf_per_encounter_kms = (3.2, 3.0, 3.2)
        self.verdict_audit = audit


def test_g3_missing_primary_is_inconclusive() -> None:
    sig = signature_from_review_entry(_Entry({"n_rev": [0, 0]}))
    assert sig.primary == ""
    r = _check(sig)
    assert r.status == "inconclusive" and "Primary unknown" in r.notes, r


def test_g3_stamped_primary_unchanged() -> None:
    sig = signature_from_review_entry(_Entry({"primary": "Jupiter"}))
    assert sig.primary == "Jupiter"


# --- G4: HR inside the H&M scope ---------------------------------------------------------------


def test_g4_hr_return_matches_hm() -> None:
    r = _check(
        _sig(
            "Sun",
            ("E", "V"),
            topology_label=RM,
            working_bodies="two",
            return_types=frozenset({"FR", "HR"}),
        )
    )
    assert r.status == "published" and "Hollister" in str(r.citation), r


# --- G5: the widened F13 ratchet flags each evasion form ----------------------------------------


def test_g5_ratchet_flags_evasion_forms() -> None:
    flag = _ratchet()
    forms = {
        "name_in_set": "FRESH = {'x'}\nok = result_status in FRESH\n",
        "match": "match status:\n    case 'not-found':\n        pass\n",
        "startswith": "ok = s.startswith('not-found')\n",
        "attr_vs_name": "ok = r.status == WANT\n",
        "get_vs_name": "ok = d.get('status') != WANT\n",
        "subscript_vs_name": "ok = d['status'] == WANT\n",
        "literal": "ok = r.status == 'not-found'\n",
        # #1045 A1-4
        "alias": "st = r.status\nok = st in FRESH\n",
        "attr_set": "ok = r.status in mod.OTHER_SET\n",
        "attr_const": "ok = r.status == consts.NOT_FOUND\n",
    }
    for name, src in forms.items():
        assert flag(src), name
    assert not flag("ok = r.status in FRESH_STATUSES\n")
    assert not flag("ok = r.status in lc.FRESH_STATUSES\n")
    assert not flag("ok = is_literature_fresh(r.status)\n")


# --- #1045 A1: G1' (mixed-label tour anchors), A1-2, A1-6 -----------------------------------


def test_a1_jovian_untopologied_not_published_by_niehoff() -> None:
    r = _check(_sig("Jupiter", ("Ganymede", "Callisto", "Ganymede"), period_k=3))
    assert r.status != "published", r
    assert "Niehoff" not in str(r.citation), r


def test_a1_capped_best_hit_names_f7_anchor() -> None:
    r = _check(_sig("Sun", ("E", "V")))
    assert r.status == "inconclusive", r
    assert "Hollister / Hollister-Menning" in r.notes.split("Corpus anchors consulted")[0], r


def test_a1_empty_sequence_not_f7() -> None:
    r = _check(_sig("Jupiter", (), working_bodies="two"))
    assert "different declared architecture" not in r.notes, r

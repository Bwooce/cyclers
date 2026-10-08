"""#972: literature_check gate v2 fixes F2-F8, F13, F14 -- pinned probes (fast, offline).

Pre-registration: docs/notes/2026-10-07-972-literature-gate-v2-preregistration.md (sec. 1 the
fixes and their probes, sec. 6 F13/F14, sec. 7 the F7 reading and F8). Every probe below is the
case that exposed the defect; the expected side is the pre-registered rule, not a value read off
the code.
"""

from __future__ import annotations

import ast
import dataclasses
from pathlib import Path

import pytest

from cyclerfinder.search import literature_check as lc
from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    CandidateSignature,
    _architecture_anchors,
    _declared_scope_exclusions,
    _has_consecutive_same_body,
    check_literature,
    is_literature_fresh,
    is_novelty_claimable,
    offline_corpus_search,
)

REPO = Path(__file__).resolve().parents[2]
RM = frozenset({"repeated-moon"})


def _anchor(prefix: str) -> lc.CorpusAnchor:
    hits = [a for a in KNOWN_CORPUS if a.name.startswith(prefix)]
    assert len(hits) == 1, (prefix, [a.name for a in hits])
    return hits[0]


def _sig(primary: str, seq: tuple[str, ...], **kw: object) -> CandidateSignature:
    base: dict[str, object] = {
        "primary": primary,
        "sequence": seq,
        "period_k": None,
        "vinf_per_encounter_kms": (),
        "topology_label": RM,
    }
    base.update(kw)
    return CandidateSignature(**base)  # type: ignore[arg-type]


def _check(sig: CandidateSignature) -> lc.LiteratureCheckResult:
    return check_literature(sig, search=offline_corpus_search)


HM = "Hollister / Hollister-Menning"
RS_GC = "Russell-Strange 2009 Ganymede-Callisto"
GCGC = "Campagnola et al. Ganymede-Callisto GCGC"
JONES = "Jones et al. VEM triple cyclers"
LIANG = "Liang et al. Callisto-Ganymede-Europa"
HUGHES = "Hughes-Edelman-Longuski"


# --- each declared tag, positive and negative, with a mutated-tag check -------------------------


def test_n_bodies_tag() -> None:
    jones = _anchor(JONES)
    vem = _sig("Sun", ("E", "M", "E", "V", "V", "E"))
    assert "n-bodies" not in _declared_scope_exclusions(vem, jones)
    assert "n-bodies" in _declared_scope_exclusions(
        vem, dataclasses.replace(jones, n_bodies_scope=2)
    )


def test_working_bodies_tag() -> None:
    rs = _anchor(RS_GC)
    one = _sig("Jupiter", ("Ganymede", "Callisto"), working_bodies="one")
    two = _sig("Jupiter", ("Ganymede", "Callisto"), working_bodies="two")
    assert "working-bodies" not in _declared_scope_exclusions(one, rs)
    assert "working-bodies" in _declared_scope_exclusions(two, rs)
    assert "working-bodies" in _declared_scope_exclusions(
        one, dataclasses.replace(rs, working_bodies_scope="two")
    )


def test_return_types_tag() -> None:
    hm = _anchor(HM)
    ok = _sig("Sun", ("E", "V"), working_bodies="two", return_types=frozenset({"FR", "SY"}))
    gen = _sig("Sun", ("E", "V"), working_bodies="two", return_types=frozenset({"FR", "GEN"}))
    assert _declared_scope_exclusions(ok, hm) == []
    assert "return-types" in _declared_scope_exclusions(gen, hm)
    assert "return-types" in _declared_scope_exclusions(
        ok, dataclasses.replace(hm, return_types_scope=frozenset({"FR"}))
    )


def test_alternating_tag() -> None:
    gcgc = _anchor(GCGC)
    alt = _sig("Jupiter", ("Ganymede", "Callisto", "Ganymede", "Callisto"), working_bodies="two")
    rep = _sig("Jupiter", ("Ganymede", "Callisto", "Callisto", "Ganymede"), working_bodies="two")
    assert "alternating" not in _declared_scope_exclusions(alt, gcgc)
    assert "alternating" in _declared_scope_exclusions(rep, gcgc)
    assert "alternating" not in _declared_scope_exclusions(
        rep, dataclasses.replace(gcgc, alternating_scope=None)
    )


# --- F2: a treated system is never "never-treated" ---------------------------------------------


def test_f2_treated_system_is_not_known_architecture_new_system() -> None:
    sig = _sig(
        "Jupiter",
        ("Ganymede", "Callisto", "Ganymede"),
        working_bodies="one",
        topology_label=frozenset({"resonant-hopping"}),
    )
    assert _architecture_anchors(sig) == []
    assert _check(sig).status != "known-architecture-new-system"


def test_known_architecture_new_system_still_reachable() -> None:
    """The new status survives F2: a "one" architecture at a Jovian pair no "one" anchor holds."""
    sig = _sig("Jupiter", ("Io", "Callisto"), working_bodies="one")
    assert _architecture_anchors(sig)
    r = _check(sig)
    assert r.status == "known-architecture-new-system", r
    assert is_literature_fresh(r.status)
    assert is_novelty_claimable({"checked": True, "status": r.status})


# --- F3: closing repeat -------------------------------------------------------------------------


def test_f3_closing_repeat_dropped() -> None:
    assert not _has_consecutive_same_body(("G", "C", "G", "C", "G"))
    assert not _has_consecutive_same_body(("G", "C", "G", "C"))
    assert _has_consecutive_same_body(("G", "C", "C", "G"))
    r = _check(
        _sig(
            "Jupiter",
            ("Ganymede", "Callisto", "Ganymede", "Callisto", "Ganymede"),
            working_bodies="two",
        )
    )
    assert r.status == "published" and "Campagnola" in str(r.citation), r


# --- F4: an anchor-synthesised hit outside the footprint is not scored --------------------------


def test_f4_off_footprint_anchor_hit_not_scored() -> None:
    r = _check(_sig("Sun", ("E", "J")))
    assert "Koon" not in str(r.citation), r
    # #1045 G7: the exact outcome, not just "not the old one".
    assert r.status == "not-found", r
    assert "outside the candidate's primary/body-set footprint were not scored" in r.notes, r


# --- F7 (ruling (b)): different architecture at a treated system -> inconclusive ----------------


def test_f7_gc2_like_goes_to_a_human() -> None:
    r = _check(
        _sig(
            "Jupiter",
            ("Ganymede", "Callisto", "Callisto", "Ganymede"),
            working_bodies="two",
            return_types=frozenset({"SY"}),
        )
    )
    assert r.status == "inconclusive", r
    assert "Russell-Strange 2009 Ganymede-Callisto" in r.notes


def test_f7_reading_b_working_bodies_among_reasons() -> None:
    """ev-C-like: H&M excludes by working-bodies AND return-types; ruling (b) still escalates."""
    sig = _sig("Sun", ("V", "V", "E"), working_bodies="one", return_types=frozenset({"GEN"}))
    assert {"working-bodies", "return-types"} <= set(_declared_scope_exclusions(sig, _anchor(HM)))
    r = _check(sig)
    assert r.status == "inconclusive" and HM in r.notes, r


# --- F5 + F8: the Hughes anchor ------------------------------------------------------------------


def test_f5_f8_hughes_anchor_grounded() -> None:
    h = _anchor(HUGHES)
    assert h.n_bodies_scope is None
    assert h.doi == "10.2514/6.2014-4109"
    assert h.topology_label == frozenset({"mga-tour"})
    assert h.provenance == "verified-against-source"
    assert "AAS 14-822" not in h.citation


def test_f8_ev_a_like_not_claimed_by_hughes() -> None:
    sig = _sig(
        "Sun", ("V", "V", "V", "E"), working_bodies="two", return_types=frozenset({"FR", "GEN"})
    )
    r = _check(sig)
    assert "Hughes" not in str(r.citation)
    assert r.status == "not-found", r


# --- F14: architecture-scoped anchors need declared labels --------------------------------------


def test_f14_unlabelled_sun_ve_does_not_match_hm() -> None:
    r = _check(_sig("Sun", ("V", "E")))
    assert not (r.status == "published" and "Hollister" in str(r.citation)), r
    assert "working-bodies-unlabelled" in _declared_scope_exclusions(
        _sig("Sun", ("V", "E")), _anchor(HM)
    )


def test_f14_labelled_hm_member_matches_hm() -> None:
    r = _check(_sig("Sun", ("E", "V"), working_bodies="two", return_types=frozenset({"FR", "SY"})))
    assert r.status == "published" and "Hollister" in str(r.citation), r


# --- sec. 3 controls: triple cyclers -------------------------------------------------------------


def test_triple_controls_published() -> None:
    jones = _check(_sig("Sun", ("E", "M", "E", "V", "V", "E"), period_k=2))
    assert jones.status == "published" and "Jones" in str(jones.citation), jones
    liang = _check(_sig("Jupiter", ("Callisto", "Ganymede", "Callisto", "Europa", "Callisto")))
    assert liang.status == "published" and "Liang" in str(liang.citation), liang


def test_triple_control_mutated_tag_negative(monkeypatch: pytest.MonkeyPatch) -> None:
    mutated = tuple(
        dataclasses.replace(a, n_bodies_scope=2) if a.name.startswith(JONES) else a
        for a in KNOWN_CORPUS
    )
    monkeypatch.setattr(lc, "KNOWN_CORPUS", mutated)
    r = _check(_sig("Sun", ("E", "M", "E", "V", "V", "E"), period_k=2))
    assert "Jones" not in str(r.citation), r


# --- F13: one helper; no literal "not-found" comparison outside it --------------------------------


def test_f13_helper() -> None:
    assert is_literature_fresh("not-found")
    assert is_literature_fresh("known-architecture-new-system")
    assert not is_literature_fresh("published")
    assert not is_literature_fresh("inconclusive")
    assert not is_literature_fresh(None)


_NF = "not-found"


def _is_status_access(node: ast.expr) -> bool:
    """``x_status`` / ``obj.status`` / ``d["status"]`` / ``d.get("status")`` (#1045 G5)."""
    if isinstance(node, ast.Name):
        return node.id.endswith("status")
    if isinstance(node, ast.Attribute):
        return node.attr.endswith("status")
    if isinstance(node, ast.Subscript):
        return isinstance(node.slice, ast.Constant) and node.slice.value == "status"
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
        return (
            node.func.attr == "get"
            and bool(node.args)
            and isinstance(node.args[0], ast.Constant)
            and node.args[0].value == "status"
        )
    return False


def _status_compare_violations_in(text: str, filename: str = "<src>") -> list[int]:
    """Lines comparing a literature status outside ``is_literature_fresh`` (#972 F13,
    widened by #1045 G5): a literal "not-found" in a comparison, a ``match`` case or a
    ``startswith``/``endswith`` call; or a status access compared with a NAME other
    than ``FRESH_STATUSES``."""
    if _NF not in text and "status" not in text:
        return []
    tree = ast.parse(text, filename=filename)
    lines: list[int] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Compare):
            sides = (node.left, *node.comparators)
            for side in sides:
                consts = (
                    [side]
                    if isinstance(side, ast.Constant)
                    else list(side.elts)
                    if isinstance(side, (ast.Tuple, ast.List, ast.Set))
                    else []
                )
                if any(isinstance(c, ast.Constant) and c.value == _NF for c in consts):
                    lines.append(node.lineno)
                    break
            else:
                if any(_is_status_access(x) for x in sides) and any(
                    isinstance(x, ast.Name)
                    and not _is_status_access(x)
                    and x.id != "FRESH_STATUSES"
                    for x in sides
                ):
                    lines.append(node.lineno)
        elif isinstance(node, ast.MatchValue):
            if isinstance(node.value, ast.Constant) and node.value.value == _NF:
                lines.append(node.value.lineno)
        elif (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in ("startswith", "endswith")
            and any(isinstance(a, ast.Constant) and a.value == _NF for a in node.args)
        ):
            lines.append(node.lineno)
    return lines


def _literal_not_found_compares(path: Path) -> list[int]:
    return _status_compare_violations_in(path.read_text(), str(path))


def test_f13_no_literal_not_found_comparisons() -> None:
    allowed = REPO / "src" / "cyclerfinder" / "search" / "literature_check.py"
    bad = []
    for root in (REPO / "src", REPO / "scripts"):
        for path in sorted(root.rglob("*.py")):
            if path == allowed:
                continue
            bad += [f"{path.relative_to(REPO)}:{n}" for n in _literal_not_found_compares(path)]
    assert bad == [], f"compare a literature status via is_literature_fresh(): {bad}"

#!/usr/bin/env python3
"""Schema v4.5 (spec §16.7.12 / §14): back-fill the ``validation_level`` tag.

Applies spec §14's V0-V5 gauntlet definitions MECHANICALLY to ``data/catalogue.yaml``
from RECORDED in-repo test evidence — never aspirationally. The level a row earns
is the highest §14 gate for which a green, teeth-bearing test exists; when the
recorded evidence does not mechanically justify a higher gate, the row stays V0
(the internal-consistency floor). NO new physics value is introduced and NO
external source is consulted; this is a derived-metadata stamp keyed to the
evidence registry below.

Evidence applied (verified green against the live suite, 2026-06-06)
-------------------------------------------------------------------
* **V2** — ``aldrin-classic-em-k1-outbound`` (2026-06-07, the §14 V2 class-split
  amendment): the powered Aldrin outbound clears the amended **V2-powered** gate —
  >=3 consecutive in-family cycles, each achieving its encounters with the
  per-cycle maintenance applied (Mars-flyby V_inf continuity bounded;
  strictly-positive in-family maintenance dV) AND bounded intra-cycle drift vs the
  planned trajectory (per-leg Kepler forward-reprop residual). The original single
  V2 gate (cross-cycle rotating-frame-repeat drift) is structurally unsatisfiable
  for a per-cycle-retargeted cycler (~4.14e8 km / 3 laps, #134). Demonstrated with
  teeth by ``tests/verify/test_aldrin_v2_powered.py::
  test_aldrin_outbound_passes_v2_powered``.

* **V1** — the real-DE440 Aldrin INBOUND twin clears spec §14 V1 — every leg
  re-solved with ``lamberthub`` izzo2015 + gooding1990 agrees to <
  ``V1_TOLERANCE_MPS``, AND the Kepler forward re-propagation residual passes (the
  §14 "re-propagated with the Kepler propagator, planet positions met < tol"
  half). The outbound twin's V1 is subsumed by its V2-powered promotion above.
  Demonstrated with teeth by
  ``tests/verify/test_agreement_lamberthub.py::test_report_includes_lamberthub_path``
  and ``...::test_real_eph_paths_a_and_c_pass_b_flags_model_mismatch``.

* **V0** — the rows that have actually been exercised by the V0/real-closure
  gauntlet machinery (the M6b regression set + the two CONSTRUCTIBLE rows). Each
  carries internal-consistency evidence (constructed + propagated), but no
  recorded evidence reaches a higher §14 gate:
    - the Aldrin INBOUND row is powered and no test builds/cross-checks it on real
      ephemeris (its outbound twin is the binding gate) ⇒ V0;
    - the other regression rows are EXPECTED_SKIPS (incomplete leg data / wrong
      topology) and do not pass real-closure ⇒ V0;
    - the two CONSTRUCTIBLE rows are the Aldrin pair (already covered).

Everything NOT in the registry is left UNTAGGED — an absent ``validation_level``
is the explicit V0 floor for downstream views. Tagging only evidence-backed rows
keeps the stamp auditable and avoids asserting a gauntlet pass for the ~230 rows
the auto-pipeline has not run.

Idempotent: a row that already carries ``validation_level`` is skipped.

Usage::

    uv run python scripts/backfill_validation_level.py            # apply
    uv run python scripts/backfill_validation_level.py --dry-run  # report only
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import atomic_write_text
from cyclerfinder.data.validate import has_level_evidence

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"

# id -> validation_level, applied mechanically from recorded test evidence.
# V1 for the Aldrin outbound (real-DE440 lamberthub + Kepler reprop); V0 for the
# rows the gauntlet machinery has exercised at the internal-consistency floor.
_LEVEL_BY_ID: dict[str, str] = {
    # 2026-06-07: the §14 V2 class-split amendment promotes the powered Aldrin
    # outbound to V2-powered (>=3 consecutive in-family cycles, each achieving its
    # encounters with the per-cycle maintenance applied AND bounded intra-cycle
    # drift vs the planned trajectory). The inbound stays V1 (its real-window solve
    # lands on a ballistic dV~0 off-family neighbour; maintenance not demonstrated).
    "aldrin-classic-em-k1-outbound": "V2",
    "aldrin-classic-em-k1-inbound": "V1",  # #125 Part 1: real-DE440 inbound clears §14 V1
    # #137 Part 1: three Russell free-return rows whose single heliocentric ellipse
    # forms a CLOSED, V_inf-continuous E->M->E cycler clear §14 V1 mechanics
    # like-for-like on the circular ephemeris (lamberthub + Kepler reprop, Mars
    # V_inf continuity intact). The other matched rows are genuinely multi-arc
    # (forced return breaks continuity) and stay V0.
    "russell-ch4-5.30gGf3": "V1",
    "russell-ch4-9.94Gg3": "V1",
    "russell-ch4-5.75ggF3": "V1",
    # #137 Part 3: dense-phase-scan promotes the deep-aphelion 9.353Gg2 to a closed,
    # V_inf-continuous free-return arc that clears §14 V1 mechanics.
    "russell-ch4-9.353Gg2": "V1",
    # closer sweep 2026-06-08 (#142 continuation): six Russell 2004 Table 3.4 rows
    # whose single circular-coplanar ellipse closes to a V_inf-continuous E->M->E
    # cycler and clears §14 V1 mechanics like-for-like (emerged V_inf within
    # 0.5 km/s of the sourced anchor, lamberthub + Kepler reprop, Mars V_inf
    # continuity intact). See docs/notes/2026-06-08-closer-sweep-v1-candidates.md.
    "russell-ocampo-3.1.1+2": "V1",  # closer sweep 2026-06-08: closed free-return, §14 V1 pass
    "russell-ocampo-3.1.3+0": "V1",  # closer sweep 2026-06-08: closed free-return, §14 V1 pass
    "russell-ocampo-4.1.1-4": "V1",  # closer sweep 2026-06-08: closed free-return, §14 V1 pass
    "russell-ocampo-4.1.2-2": "V1",  # closer sweep 2026-06-08: closed free-return, §14 V1 pass
    "russell-ocampo-4.1.4-1": "V1",  # closer sweep 2026-06-08: closed free-return, §14 V1 pass
    "russell-ocampo-4.6.3+0": "V1",  # closer sweep 2026-06-08: closed free-return, §14 V1 pass
    # #181 ToF-fix (2026-06-10): four descriptor-bearing Russell Ch.4 rows close on
    # the real DE440 ephemeris via the joint (epoch, ToF) self-seed closer (the
    # Stage-B coplanar-branch-ToF artifact corrected) — emerged E/M v∞ within
    # 0.08 km/s of the sourced Russell 2004 Table 4.x anchor, lamberthub
    # izzo2015+gooding1990 + Kepler reprop pass (both §14 V1 halves). All four are
    # V3-CANDIDATES (single-leg REBOUND/IAS15 confirm in-band; the multi-lap
    # horizon-TCM is the named follow-up). The siblings 9.353Gg2 / 9.94Gg3 were
    # already V1 (circular like-for-like); the real-eph closure is ADDED EVIDENCE
    # for them, not a level change. See
    # docs/notes/2026-06-10-tof-fix-closure-results.md.
    "russell-ch4-3.78Gg3": "V1",
    "russell-ch4-6.44Gg3": "V1",
    "russell-ch4-3.64gGg3": "V1",
    "russell-ch4-5.30ggF3": "V1",
    # #167/#94: S1L1's FIRST V3 — corrected-topology closure CONFIRMED on the real
    # DE440 ephemeris (independent REBOUND/IAS15, all 7 Mars encounters in 3x Mars-SOI
    # at the published per-leg v∞). docs/notes/2026-06-08-s1l1-corrected-closure-results.md.
    "russell-ch4-4.991gG2": "V3",
    # #175: the §14 V3 class-split. russell-ch4-8.049gGf2 (App-C #188) is V3-POWERED
    # — encounters confirmed on real DE440 by an independent REBOUND/IAS15 integrator
    # (#170), and its continuous-from-one-seed TCM (163.6 m/s) is under its OWN
    # documented App-C ΔV (420 m/s, Russell 2004 Table 5.5; 0.39x budget), not the
    # 120 m/s ballistic bar. Encodes as the V3 enum value (the powered class lives in
    # the _LEVEL_EVIDENCE text, as V2-powered does under V2). The sibling #192
    # (russell-ch4-8.165Gfh-f2) is NOT promoted: its TCM 2040.6 m/s EXCEEDS its own
    # 1678 m/s budget (1.22x). See docs/notes/2026-06-08-v3-powered-classsplit.md.
    "russell-ch4-8.049gGf2": "V3",
    # #216 (2026-06-12, USER-approved): five Ross & Roberts-Tsoukkas 2025 (AAS
    # 25-621) stable prograde Earth-Moon (k1,k2)-cyclers clear §14 V1 like-for-like
    # in the planar CR3BP — same-model reproduction of a sourced same-model orbit
    # (corrector closes on the published C^stable/T^stable, Barden |nu|<1 STABLE
    # verdict reproduced, AND an independent Radau integrator closes the full period
    # dJ<1e-12). state_nd DERIVED (not a golden). Adopted #212b (commit 4be2375).
    # See docs/notes/2026-06-12-ross-adoption-results.md.
    # #229 (2026-06-13, USER-approved): promoted V1 -> V2 (V2-ballistic, spec §14):
    # 100-period inertial REBOUND/IAS15 bounded-band run (3x the measured-delta1
    # discrimination threshold) in the rows' DEFINING MODEL (CR3BP, like-for-like).
    # STABLE verdict from Barden nu + the long-span run, never the 5-period gate.
    # No real-ephemeris claim; V2-ballistic is this lane's ceiling.
    # docs/notes/2026-06-13-ross-v2-longspan-evidence.md.
    "ross-rt-em-cycler-11-2025": "V2",
    "ross-rt-em-cycler-21-2025": "V2",
    "ross-rt-em-cycler-31-2025": "V2",
    "ross-rt-em-cycler-32-2025": "V2",
    "ross-rt-em-cycler-33-2025": "V2",
    # #222 (2026-06-13, USER-approved): Liang et al. 2024 (JGCD, DOI
    # 10.2514/1.G008387) idealized CGE triple-cycler Members A/B/C clear §14 V1
    # like-for-like — same-model reconstruction (search/cge_scaffold.py) matches
    # the printed per-flyby V_inf (Tables 3/5/7) inside the Table 1
    # print-quantization tolerance; all legs 1-rev multi-rev Lambert as published.
    # Member D (ephemeris-2033) stays V0: no numeric per-flyby anchors published.
    # docs/notes/2026-06-13-liang-abc-reproduction.md.
    "liang-2024-cgcec-111-highperijove": "V1",
    "liang-2024-cgcec-110-highperijove": "V1",
    "liang-2024-cgcec-111-lowperijove": "V1",
    # #249 (2026-06-15, USER-approved writeback): three Braik & Ross 2026 (arXiv
    # 2605.31543) common-energy Earth-Moon CR3BP cycler reproductions at
    # C_J=3.1294. Clear spec §14 V1 like-for-like in the planar CR3BP -- same-
    # model fixed-Jacobi corrector closes (C enforced to machine eps, T matches
    # Braik-Ross Table-2 to 0.001%-0.06%), Barden nu reproduces the published
    # Floquet UNSTABLE verdict (sigma>0), independent Radau preserves Jacobi
    # (dJ<1e-12). NOT V2: unstable orbits cannot satisfy V2-ballistic's
    # bounded-drift gate. state_nd DERIVED (not a golden).
    # docs/notes/2026-06-14-249-unstable-member-recovery-plan.md (final
    # disposition section; commits a19eb24, 4a20243, f608c6b, 325c8a2).
    "braik-ross-c11a-cycler-2026": "V1",
    "braik-ross-c11b-cycler-2026": "V1",
    "braik-ross-c32-cycler-2026": "V1",
    # #494 (2026-06-30): Ross & Roberts-Tsoukkas 2026 mu-family Table-I
    # representatives + Pluto-Charon (3,2) instantiation (all V1, task #494).
    "ross-rt-mu001-cycler-11-2026": "V1",
    "ross-rt-mu01-cycler-32-2026": "V1",
    "ross-rt-mu03-cycler-31-2026": "V1",
    "ross-rt-mu05-cycler-11-2026": "V1",
    "ross-rt-pc-cycler-32-2026": "V2",
    "mcconaghy-2006-em-k2": "V0",
    "russell-ocampo-2.1.1+2-case2": "V0",
    "russell-ocampo-2.5.1+0": "V0",
}


def preflight_registry_drift(level_by_id: dict[str, str]) -> list[str]:
    """M1 guard (#196): refuse any (id, level>V0) absent from ``_LEVEL_EVIDENCE``.

    ``_LEVEL_BY_ID`` (this script) and ``_LEVEL_EVIDENCE``
    (``cyclerfinder.data.validate``) are maintained by hand in two files; if
    they drift, this script could stamp a level the over-claim guard would
    then reject — or worse, a level nothing ever checks. Every entry above
    the V0 floor must be backed by the SAME registry the validator reads
    (via :func:`has_level_evidence`). Returns a list of drift descriptions
    (empty when clean).
    """
    errors: list[str] = []
    for rid, lvl in sorted(level_by_id.items()):
        if not has_level_evidence(rid, lvl):
            errors.append(
                f"{rid}: level {lvl!r} is in _LEVEL_BY_ID but ({rid!r}, {lvl!r}) "
                f"has no entry in cyclerfinder.data.validate._LEVEL_EVIDENCE — "
                f"register the evidence pointer there first (registries drifted)"
            )
    return errors


def _render_line(level: str) -> str:
    return (
        f"  validation_level: {level}"
        "   # schema v4.5 (2026-06-06, spec §16.7.12 / §14); MECHANICAL gauntlet "
        "level from recorded test evidence (scripts/backfill_validation_level.py); "
        "see data/README.md\n"
    )


def main() -> None:
    dry_run = "--dry-run" in sys.argv

    # M1 preflight (#196): _LEVEL_BY_ID must agree with the validator's
    # evidence registry BEFORE any write — otherwise this script could stamp
    # an over-claim the suite would only catch after the file changed.
    drift = preflight_registry_drift(_LEVEL_BY_ID)
    if drift:
        print("ERROR: evidence-registry drift; refusing to write:", file=sys.stderr)
        for line in drift:
            print(f"  {line}", file=sys.stderr)
        sys.exit(1)

    raw_text = CATALOGUE_PATH.read_text(encoding="utf-8")
    rows = yaml.safe_load(raw_text)

    to_tag: dict[str, str] = {}
    already: set[str] = set()
    for row in rows:
        rid = row["id"]
        if rid not in _LEVEL_BY_ID:
            continue
        if "validation_level" in row:
            already.add(rid)
            continue
        to_tag[rid] = _LEVEL_BY_ID[rid]

    print(f"rows total: {len(rows)}")
    print(f"evidence-registry rows: {len(_LEVEL_BY_ID)}")
    print(f"already tagged (skipped): {sorted(already)}")
    print(f"rows to tag this run: {len(to_tag)}")
    for rid, lvl in sorted(to_tag.items()):
        print(f"  {rid}: {lvl}")

    if dry_run:
        print("\n[dry-run] no file written.")
        return
    if not to_tag:
        print("\nNothing to insert (all evidence-registry rows already tagged).")
        return

    # Line-by-line insertion after each row's cycler_class: line (matches the
    # provenance backfill's insertion site; preserves all comments/formatting).
    lines = raw_text.splitlines(keepends=True)
    out: list[str] = []
    current_id: str | None = None
    inserted = 0
    for line in lines:
        stripped = line.lstrip()
        if line.startswith("- id:"):
            slug = stripped.split(":", 1)[1].strip()
            if "#" in slug:
                slug = slug[: slug.index("#")].strip()
            current_id = slug.strip().strip('"')
        out.append(line)
        if line.startswith("  cycler_class:") and current_id is not None and current_id in to_tag:
            out.append(_render_line(to_tag[current_id]))
            inserted += 1

    if inserted != len(to_tag):
        print(
            f"ERROR: expected to insert {len(to_tag)} lines, inserted {inserted}. Aborting.",
            file=sys.stderr,
        )
        sys.exit(1)

    # N1 (#196): atomic replace — a writer dying mid-write must not truncate
    # the catalogue.
    atomic_write_text(CATALOGUE_PATH, "".join(out))
    print(f"\nInserted validation_level into {inserted} rows -> {CATALOGUE_PATH}")


if __name__ == "__main__":
    main()

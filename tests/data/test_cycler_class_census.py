"""Task 2.1: Census ratchet test for cycler_class tags on all catalogue rows.

Verifies:
1. Every row has a ``cycler_class`` key (no missing tags).
2. The class distribution is exactly {single-ellipse: 28, multi-arc: 240, non-keplerian: 9}.
3. The set of ids tagged multi-arc exactly equals the 240-id MULTI_ARC_ALLOWLIST from
   docs/notes/multi-arc-classification.md §9 (frozen ratchet).
4. The 9 non-keplerian ids match the §3 list.

Source of truth: docs/notes/multi-arc-classification.md §3 and §9.

2026-06-07 (#142 catalogue ingest, batch 1): +15 Rall 1970 (MIT TE-34) free-fall
periodic Earth-Mars orbits (rall-1970-*), all multi-arc. multi-arc 203->218,
total 237->252. See docs/notes/2026-06-07-catalogue-ingest-rall-russell.md.

2026-06-07 (#142 catalogue ingest, batch 2): +16 Russell 2004 Table 3.4
circular-coplanar cyclers (russell-ocampo-* not previously catalogued), all
multi-arc. multi-arc 218->234, total 252->268.

2026-06-07 (moon-tour Tier-1, task #76): the two Jovian patched-conic
family-seed rows re-tagged non-keplerian -> multi-arc (they are circular-coplanar
patched-conic moon tours about Jupiter, not CR3BP). multi-arc 234->236,
non-keplerian 6->4; total unchanged at 268.

2026-06-12 (#216, USER-approved writeback): +5 Ross & Roberts-Tsoukkas 2025
stable Earth-Moon CR3BP cyclers (non-keplerian 4->9) and +4 Liang et al. 2024
Callisto-Ganymede-Europa triple-cycler members (multi-arc 236->240); total
268->277.

2026-06-15 (#249, USER-approved writeback): +3 Braik & Ross 2026 common-energy
Earth-Moon CR3BP cycler reproductions (braik-ross-c11a/c11b/c32-cycler-2026,
non-keplerian 9->12) -- new (1,1)/(3,2) family members at the common Jacobi
C_J=3.1294, distinct from the per-family stable midpoints already catalogued
under ross-rt-em-cycler-*. The 4th Braik-Ross cycler (C21) already exists as
ross-rt-em-cycler-21-2025 (added there as corroborating_sources, no new row).
total 277->280.

2026-06-15 (#294, USER-approved writeback): +1 Tito 2018 Mars free-return
admitted as mga_tour under the v4.7 catalogue-scope expansion
(tito-2018-mars-free-return, cycler_class=multi-arc structural).
multi-arc 240->241; total 280->281.

2026-06-16 (#336): +1 Heaton-Longuski 2003 Uranian satellite tour U00-01
admitted as mga_tour (second mga_tour row after Tito 2018) -- 40-flyby
Galileo-style tour terminating at Ariel rendezvous V_inf=0.92 km/s
(heaton-longuski-2003-uranian-tour-u00-01, cycler_class=multi-arc structural).
multi-arc 241->242; total 281->282.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import yaml  # type: ignore[import-untyped]

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"

# ---------------------------------------------------------------------------
# Frozen ratchet: the 240-id MULTI_ARC_ALLOWLIST
# Source: docs/notes/multi-arc-classification.md §9
# (184 russell-ocampo-* + 14 russell-ch4-* + mcconaghy-2006-em-k2
#  + sanchez-net-2022-eem-cycler1 + sanchez-net-2022-em-cycler2
#  + jones-2017-vem-emevve-outbound + jones-2017-vem-meevem-inbound
#  + 15 rall-1970-* free-fall periodic orbits, #142 batch 1 = 218
#  + 16 russell-ocampo-* Table 3.4 circular-coplanar cyclers, #142 batch 2 = 234
#  + 2 Jovian patched-conic family-seed rows, moon-tour Tier-1 #76 = 236
#  + 4 Liang et al. 2024 CGE triple-cycler members, #216 = 240)
# ---------------------------------------------------------------------------
MULTI_ARC_ALLOWLIST: frozenset[str] = frozenset(
    [
        # russell-ocampo family (184 rows) - Russell 2004 Ch3, Tables 3.4-3.11
        "russell-ocampo-2.1.1+2-case2",
        "russell-ocampo-2.3.1+1-case3",
        "russell-ocampo-2.5.1+0",
        "russell-ocampo-3.1.1+3",
        "russell-ocampo-3.1.2+1",
        "russell-ocampo-3.3.1+2",
        "russell-ocampo-3.5.1+1",
        "russell-ocampo-3.5.2+0",
        "russell-ocampo-3.7.1+1",
        "russell-ocampo-3.9.1+0",
        "russell-ocampo-4.0.3+1",
        "russell-ocampo-4.1.1-5",
        "russell-ocampo-4.10.1+2",
        "russell-ocampo-4.11.1-2",
        "russell-ocampo-4.12.1+1",
        "russell-ocampo-4.13.1-1",
        "russell-ocampo-4.14.1+0",
        "russell-ocampo-4.14.1-1",
        "russell-ocampo-4.3.1-4",
        "russell-ocampo-4.3.1-5",
        "russell-ocampo-4.5.1-3",
        "russell-ocampo-4.5.1-4",
        "russell-ocampo-4.5.2-2",
        "russell-ocampo-4.5.3-1",
        "russell-ocampo-4.7.1-3",
        "russell-ocampo-4.9.1-2",
        "russell-ocampo-4.9.2-1",
        "russell-ocampo-5.1.1-7",
        "russell-ocampo-5.1.2-3",
        "russell-ocampo-5.1.5-1",
        "russell-ocampo-5.2.1+7",
        "russell-ocampo-5.2.2+2",
        "russell-ocampo-5.2.5+0",
        "russell-ocampo-5.3.1-6",
        "russell-ocampo-5.3.1-7",
        "russell-ocampo-5.3.3-2",
        "russell-ocampo-5.4.1+5",
        "russell-ocampo-5.4.1+6",
        "russell-ocampo-5.4.3+1",
        "russell-ocampo-5.5.1-4",
        "russell-ocampo-5.5.1-5",
        "russell-ocampo-5.5.1-6",
        "russell-ocampo-5.5.2-2",
        "russell-ocampo-5.5.2-3",
        "russell-ocampo-5.5.4-1",
        "russell-ocampo-5.6.1+3",
        "russell-ocampo-5.6.1+4",
        "russell-ocampo-5.6.1+5",
        "russell-ocampo-5.6.2+1",
        "russell-ocampo-5.6.2+2",
        "russell-ocampo-5.6.4+0",
        "russell-ocampo-5.7.1-3",
        "russell-ocampo-5.7.1-4",
        "russell-ocampo-5.7.1-5",
        "russell-ocampo-5.8.1+2",
        "russell-ocampo-5.8.1+3",
        "russell-ocampo-5.8.1+4",
        "russell-ocampo-5.9.1-2",
        "russell-ocampo-5.9.1-3",
        "russell-ocampo-5.9.1-4",
        "russell-ocampo-5.9.2-1",
        "russell-ocampo-5.9.2-2",
        "russell-ocampo-5.9.3-1",
        "russell-ocampo-5.10.1+2",
        "russell-ocampo-5.10.1+3",
        "russell-ocampo-5.10.2+0",
        "russell-ocampo-5.10.2+1",
        "russell-ocampo-5.10.3+0",
        "russell-ocampo-5.11.1-2",
        "russell-ocampo-5.11.1-3",
        "russell-ocampo-5.11.2+1",
        "russell-ocampo-5.12.1+1",
        "russell-ocampo-5.12.1+2",
        "russell-ocampo-5.13.1-2",
        "russell-ocampo-5.13.1-3",
        "russell-ocampo-5.13.2-1",
        "russell-ocampo-5.14.1+1",
        "russell-ocampo-5.14.1+2",
        "russell-ocampo-5.14.2+0",
        "russell-ocampo-5.15.1-1",
        "russell-ocampo-5.15.1-2",
        "russell-ocampo-5.16.1+0",
        "russell-ocampo-5.16.1+1",
        "russell-ocampo-5.17.1-1",
        "russell-ocampo-5.18.1+0",
        "russell-ocampo-6.0.1+6d",
        "russell-ocampo-6.0.1+7c",
        "russell-ocampo-6.0.1+8b",
        "russell-ocampo-6.0.1+9a",
        "russell-ocampo-6.1.2-4",
        "russell-ocampo-6.1.3-3",
        "russell-ocampo-6.1.4-2",
        "russell-ocampo-6.1.6-1",
        "russell-ocampo-6.2.1+6",
        "russell-ocampo-6.2.1+7",
        "russell-ocampo-6.2.1+8",
        "russell-ocampo-6.2.2+2",
        "russell-ocampo-6.2.2+3",
        "russell-ocampo-6.2.3+1",
        "russell-ocampo-6.2.3+2",
        "russell-ocampo-6.2.4+1",
        "russell-ocampo-6.2.6+0",
        "russell-ocampo-6.3.1-9",
        "russell-ocampo-6.3.4+1",
        "russell-ocampo-6.4.1+4",
        "russell-ocampo-6.4.1+5",
        "russell-ocampo-6.4.1+6",
        "russell-ocampo-6.4.1+7",
        "russell-ocampo-6.5.1-6",
        "russell-ocampo-6.5.1-7",
        "russell-ocampo-6.5.1-8",
        "russell-ocampo-6.5.5-1",
        "russell-ocampo-6.6.1+3",
        "russell-ocampo-6.6.1+4",
        "russell-ocampo-6.6.1+5",
        "russell-ocampo-6.6.1+6",
        "russell-ocampo-6.6.2+1",
        "russell-ocampo-6.6.2+2",
        "russell-ocampo-6.6.5+0",
        "russell-ocampo-6.7.1-6",
        "russell-ocampo-6.7.1-7",
        "russell-ocampo-6.7.2+3",
        "russell-ocampo-6.7.3-2",
        "russell-ocampo-6.7.5+0",
        "russell-ocampo-6.8.1+2",
        "russell-ocampo-6.8.1+3",
        "russell-ocampo-6.8.1+4",
        "russell-ocampo-6.8.1+5",
        "russell-ocampo-6.8.1+6",
        "russell-ocampo-6.8.3+1",
        "russell-ocampo-6.9.1-4",
        "russell-ocampo-6.9.1-5",
        "russell-ocampo-6.9.1-6",
        "russell-ocampo-6.9.2-2",
        "russell-ocampo-6.9.2-3",
        "russell-ocampo-6.9.4-1",
        "russell-ocampo-6.10.1+2",
        "russell-ocampo-6.10.1+3",
        "russell-ocampo-6.10.1+4",
        "russell-ocampo-6.10.1+5",
        "russell-ocampo-6.10.2+1",
        "russell-ocampo-6.10.2+2",
        "russell-ocampo-6.10.4+0",
        "russell-ocampo-6.11.1-4",
        "russell-ocampo-6.11.1-5",
        "russell-ocampo-6.11.2+2",
        "russell-ocampo-6.12.1+2",
        "russell-ocampo-6.12.1+3",
        "russell-ocampo-6.12.1+4",
        "russell-ocampo-6.13.1-3",
        "russell-ocampo-6.13.1-4",
        "russell-ocampo-6.13.1-5",
        "russell-ocampo-6.13.1+5",
        "russell-ocampo-6.13.2-2",
        "russell-ocampo-6.13.3-1",
        "russell-ocampo-6.14.1+1",
        "russell-ocampo-6.14.1+2",
        "russell-ocampo-6.14.1+3",
        "russell-ocampo-6.14.2+0",
        "russell-ocampo-6.14.2+1",
        "russell-ocampo-6.14.3+0",
        "russell-ocampo-6.15.1-2",
        "russell-ocampo-6.15.1-3",
        "russell-ocampo-6.15.1-4",
        "russell-ocampo-6.15.1+4",
        "russell-ocampo-6.15.2+1",
        "russell-ocampo-6.16.1+1",
        "russell-ocampo-6.16.1+2",
        "russell-ocampo-6.17.1-2",
        "russell-ocampo-6.17.1-3",
        "russell-ocampo-6.17.1+3",
        "russell-ocampo-6.17.2-1",
        "russell-ocampo-6.18.1+1",
        "russell-ocampo-6.18.1+2",
        "russell-ocampo-6.18.2+0",
        "russell-ocampo-6.19.1-2",
        "russell-ocampo-6.19.1+2",
        "russell-ocampo-6.19.2+0",
        "russell-ocampo-6.20.1-4",
        "russell-ocampo-6.20.1+0",
        "russell-ocampo-6.20.1+1",
        "russell-ocampo-6.21.1-1",
        "russell-ocampo-6.21.1+1",
        "russell-ocampo-6.22.1+0",
        # russell-ch4 family (14 rows) - Russell 2004 Ch4, Tables 4.9-4.13
        "russell-ch4-4.991gG2",
        "russell-ch4-8.049gGf2",
        "russell-ch4-8.165Gfh-f2",
        "russell-ch4-9.353Gg2",
        "russell-ch4-3.64gGg3",
        "russell-ch4-3.77Gh3",
        "russell-ch4-3.78Gg3",
        "russell-ch4-5.30gGf3",
        "russell-ch4-5.66Gfh3",
        "russell-ch4-9.94Gg3",
        "russell-ch4-3.66gfF3",
        "russell-ch4-5.30ggF3",
        "russell-ch4-5.75ggF3",
        "russell-ch4-6.44Gg3",
        # resolved from §7 — same physical cycler as russell-ch4-4.991gG2 (Russell Table 4.9)
        "mcconaghy-2006-em-k2",
        # Sanchez Net 2022 EEM near-ballistic real-date patched-conic cycler (Fig. 2a)
        "sanchez-net-2022-eem-cycler1",
        # Sanchez Net 2022 EM near-ballistic real-date patched-conic cycler (Fig. 2b)
        "sanchez-net-2022-em-cycler2",
        # Jones 2017 VEM triple cyclers (AAS 17-577 Tables 2 & 3) — real-ephemeris
        # patched-conic, 6 flybys per cycle, no single repeating (a,e) conic.
        "jones-2017-vem-emevve-outbound",
        "jones-2017-vem-meevem-inbound",
        # Rall 1970 (MIT TE-34) free-fall periodic Earth-Mars orbits (15 rows),
        # #142 catalogue ingest batch 1. App E circular-coplanar (M4-1, M6-1..3)
        # + App F eccentric-inclined (M4-1a, M5-1a..e, M5-2a..e). Multi-arc:
        # each is multiple E-M-E round trips + Earth direct returns (note §6).
        "rall-1970-m4-1",
        "rall-1970-m4-1a",
        "rall-1970-m5-1a",
        "rall-1970-m5-1b",
        "rall-1970-m5-1c",
        "rall-1970-m5-1d",
        "rall-1970-m5-1e",
        "rall-1970-m5-2a",
        "rall-1970-m5-2b",
        "rall-1970-m5-2c",
        "rall-1970-m5-2d",
        "rall-1970-m5-2e",
        "rall-1970-m6-1",
        "rall-1970-m6-2",
        "rall-1970-m6-3",
        # Russell 2004 Table 3.4 circular-coplanar cyclers (16 rows), #142 batch
        # 2. Eligible-but-absent rows from Table 3.4 (p.83) not previously
        # catalogued (the 27 already present + Aldrin 1.0.1.-1 are excluded by
        # descriptor). Multi-arc: free-return E-M cyclers with intermediate
        # Earth flybys (n_flybys-1 loops per cycle).
        "russell-ocampo-3.1.1+2",
        "russell-ocampo-3.1.3+0",
        "russell-ocampo-3.5.1+2",
        "russell-ocampo-4.1.1-6",
        "russell-ocampo-4.1.1-4",
        "russell-ocampo-4.1.2-3",
        "russell-ocampo-4.1.2-2",
        "russell-ocampo-4.1.4-1",
        "russell-ocampo-4.6.1-4",
        "russell-ocampo-4.6.3+0",
        "russell-ocampo-4.7.1-2",
        "russell-ocampo-4.8.1+3",
        "russell-ocampo-4.8.1+2",
        "russell-ocampo-4.9.1-3",
        "russell-ocampo-4.10.1-3",
        "russell-ocampo-4.12.1-2",
        # Re-tagged non-keplerian -> multi-arc by moon-tour Tier-1 (task #76):
        # circular-coplanar patched-conic multi-arc moon tours about Jupiter.
        # Family-seed null-numeric rows (exempt from the multi-arc completeness
        # invariants — see tests/data/test_multi_arc_invariants.py).
        "hernandez-2017-jovian-ieg-triple-family",
        "russell-strange-2009-jovian-multimoon-family",
        # #216 (2026-06-12): four Liang et al. 2024 (JGCD doi 10.2514/1.G008387)
        # Callisto-Ganymede-Europa triple-cycler members — patched-conic multi-rev
        # Lambert legs, repeated-moon CGCEC sequence about Jupiter (multi-arc).
        # Members A-C carry per-leg ToF (idealized circular-coplanar); member D is
        # the SPICE-ephemeris member (figures-only, exempt from the transit-times
        # completeness invariant — see _FAMILY_SEED_NULL_NUMERIC_MULTI_ARC).
        "liang-2024-cgcec-111-highperijove",
        "liang-2024-cgcec-110-highperijove",
        "liang-2024-cgcec-111-lowperijove",
        "liang-2024-cgcec-ephemeris-2033",
        # #294 (2026-06-15): Tito 2018 Mars free-return admitted as mga_tour
        # under the schema v4.7 scope expansion. cycler_class=multi-arc is
        # structural (two distinct heliocentric legs E->M and M->E), not a
        # cycler claim — see docs/notes/2026-06-16-catalogue-scope-taxonomy.md.
        # The cycler-specific invariants tests in test_multi_arc_invariants.py
        # filter orbit_class != cycler so this row is exempt from
        # invariants{}/transit_times_days completeness.
        "tito-2018-mars-free-return",
        # #336 (2026-06-16): Heaton-Longuski 2003 Uranian satellite tour U00-01
        # admitted as mga_tour under the schema v4.7 scope expansion (second
        # mga_tour row). cycler_class=multi-arc is structural (one heliocentric
        # arc E->J->U plus a 40-flyby Uranian-system tour with multiple distinct
        # arcs between satellite encounters), not a cycler claim. As with the
        # Tito row, the cycler-specific invariants tests filter orbit_class !=
        # cycler so this row is exempt from invariants{}/transit_times_days
        # completeness.
        "heaton-longuski-2003-uranian-tour-u00-01",
        # #339 (2026-06-17): umbriel-oberon-1-1-uranian-quasi-cycler-2026 admitted
        # as catalogue's first computed quasi_cycler row (Umbriel-Oberon-Umbriel
        # (1,1) Uranian near-5:1 synodic resonance, 84-yr validity window
        # 2000-2083 per #338 EFFECTIVELY_CYCLIC verdict). cycler_class=multi-arc
        # is structural (3 distinct Lambert legs Umbriel->Oberon->Umbriel), not
        # a strict-periodic cycler claim — the v4.7 quasi_cycler class is
        # closes-up-to-rotation. Same exemption pattern as Heaton-Longuski 2003
        # and Tito 2018 from the invariants{}/transit_times_days completeness
        # tests.
        "umbriel-oberon-1-1-uranian-quasi-cycler-2026",
        # #356 (2026-06-17): damario-1992-galileo-veega admitted as the third
        # computed mga_tour row (after Tito 2018 and Heaton-Longuski 2003).
        # Galileo VEEGA flown trajectory (D'Amario-Bright-Wolf 1992 SSR
        # 60(1-4):23-78, DOI 10.1007/bf00216849). cycler_class=multi-arc is
        # structural (three distinct heliocentric arcs E->V, V->E1, E1->E2,
        # E2->J), not a cycler claim. Same exemption pattern as Tito 2018 /
        # Heaton-Longuski 2003 from the invariants{}/transit_times_days
        # completeness tests.
        "damario-1992-galileo-veega",
        # #390 (2026-06-19): the catalogue's first SPK-derived mga_tour rows --
        # Voyager 1 (E-J-S) and Voyager 2 (E-J-S-U-N, the only four-giant-planet
        # Grand Tour). cycler_class=multi-arc is structural (two/four distinct
        # heliocentric gravity-assist arcs), not a cycler claim. Per-encounter
        # V_inf DERIVED from the NAIF reconstructed spacecraft SPK at the flown
        # flyby epochs (cyclerfinder.verify.mission_spk; data/390_mission_vinf.
        # jsonl). Same exemption pattern as Tito 2018 / Heaton-Longuski 2003 /
        # D'Amario 1992 from the invariants{}/transit_times_days completeness
        # tests.
        "voyager-1-jupiter-saturn-grand-tour",
        "voyager-2-grand-tour",
        # #399 (2026-06-19): 5 more SPK-derived mga_tour rows admitted via the
        # #390 extractor (Pioneer 10/11, Cassini, Juno, Mariner-10). Each is a
        # one-shot gravity-assist tour -> cycler_class=multi-arc structural, same
        # exemption pattern as the Voyager rows above.
        "pioneer-10-jupiter-flyby",
        "pioneer-11-jupiter-saturn-flyby",
        "cassini-huygens-vvejga",
        "juno-earth-flyby-jupiter",
        "mariner-10-venus-mercury",
        # BepiColombo (#399, 2026-06-19): 1 Earth + 2 Venus + 6 Mercury gravity
        # assists, ESA reconstructed MPO SPK (NAIF -121); closes the #345 backlog.
        "bepicolombo-earth-venus-mercury",
    ]
)

assert len(MULTI_ARC_ALLOWLIST) == 252, (
    f"Allowlist must have 252 entries, got {len(MULTI_ARC_ALLOWLIST)}"
)

# ---------------------------------------------------------------------------
# Frozen ratchet: the 4 non-keplerian ids
# Source: docs/notes/multi-arc-classification.md §3
# ---------------------------------------------------------------------------
NON_KEPLERIAN_IDS: frozenset[str] = frozenset(
    [
        "arenstorf-em-figure8-1963",
        "genova-aldrin-2015-em-3petal-cycler",
        "wittal-2022-em-cycler-family",
        # The two Jovian patched-conic family-seed rows were re-tagged
        # non-keplerian -> multi-arc by moon-tour Tier-1 (task #76): they are
        # circular-coplanar patched-conic multi-arc moon tours, not CR3BP. See
        # MULTI_ARC_ALLOWLIST below. The Saturnian seed stays non-keplerian (its
        # midsize members are genuinely CR3BP — Tier-2; Titan split deferred).
        "russell-strange-2009-saturnian-multimoon-family",
        # #216 (2026-06-12): five Ross & Roberts-Tsoukkas 2025 (AAS 25-621) stable
        # prograde Earth-Moon (k1,k2)-cyclers — genuine planar CR3BP periodic
        # orbits (non-keplerian; rotating-frame, Jacobi-constant identity).
        "ross-rt-em-cycler-11-2025",
        "ross-rt-em-cycler-21-2025",
        "ross-rt-em-cycler-31-2025",
        "ross-rt-em-cycler-32-2025",
        "ross-rt-em-cycler-33-2025",
        # #249 (2026-06-15): three Braik & Ross 2026 (arXiv 2605.31543) common-
        # energy Earth-Moon (1,1)/(3,2) cycler reproductions at C_J=3.1294 -- new
        # family members at a different energy from the Ross-RT 2025 anchors
        # (rotating-frame periodic, Jacobi-constant identity). The 4th Braik-Ross
        # cycler (C21) is already catalogued as ross-rt-em-cycler-21-2025.
        "braik-ross-c11a-cycler-2026",
        "braik-ross-c11b-cycler-2026",
        "braik-ross-c32-cycler-2026",
    ]
)

assert len(NON_KEPLERIAN_IDS) == 12


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _load_rows() -> list[dict]:  # type: ignore[type-arg]
    return yaml.safe_load(CATALOGUE_PATH.read_text())  # type: ignore[no-any-return]


def test_all_rows_have_cycler_class() -> None:
    """Every row in catalogue.yaml must have a cycler_class key."""
    rows = _load_rows()
    missing = [r["id"] for r in rows if "cycler_class" not in r]
    assert missing == [], f"Rows missing cycler_class ({len(missing)} rows): {missing[:10]}"


def test_census_distribution() -> None:
    """Exact class distribution: single-ellipse=46, multi-arc=246, non-keplerian=12.

    #390 (2026-06-19) admitted the catalogue's first two SPK-derived mga_tour
    rows -- voyager-1-jupiter-saturn-grand-tour (E-J-S) and voyager-2-grand-tour
    (E-J-S-U-N, the only four-giant-planet Grand Tour), both cycler_class=
    multi-arc (structural gravity-assist arcs): multi-arc 244->246.

    (moon-tour Tier-1 task #76 re-tagged the two Jovian patched-conic family-seed
    rows non-keplerian -> multi-arc: multi-arc 234->236, non-keplerian 6->4. #216
    added 5 Ross CR3BP cyclers (non-keplerian 4->9) and 4 Liang CGE triple-cycler
    members (multi-arc 236->240). #249 (2026-06-15) added 3 Braik-Ross 2026 CR3BP
    cycler reproductions (non-keplerian 9->12). #294 (2026-06-15) admitted Tito
    2018 Mars free-return as mga_tour with cycler_class=multi-arc (structural)
    under the v4.7 scope expansion: multi-arc 240->241. #336 (2026-06-16)
    admitted Heaton-Longuski 2003 Uranian satellite tour U00-01 as mga_tour
    with cycler_class=multi-arc (one heliocentric arc + 40 Uranian-system
    flybys): multi-arc 241->242. #339 (2026-06-17) admitted Umbriel-Oberon-
    Umbriel (1,1) Uranian quasi_cycler — catalogue's first computed
    quasi_cycler row — with cycler_class=multi-arc (3 Lambert legs structural,
    not strict-periodic): multi-arc 242->243. #356 (2026-06-17) admitted
    Galileo VEEGA (D'Amario-Bright-Wolf 1992) as the third computed mga_tour
    row, cycler_class=multi-arc (E->V, V->E1, E1->E2, E2->J heliocentric arcs):
    multi-arc 243->244. #367 (2026-06-17) admitted 7 Rogers 2015 Table 4
    precursor_mga insertion-trajectory rows (VISIT-1/2, Case 1/2/3, S1L1, U0L1),
    all cycler_class=single-ellipse: single-ellipse 28->35. #367 wave 2
    (2026-06-17) admitted the VISIT-1 5:4(3)- circular-coplanar Table-3
    precursor_mga sub-variant: single-ellipse 35->36; then VISIT-2 5:4(3)-
    and 3:2(2)- circular-coplanar precursor_mga: single-ellipse 36->38;
    then Case 1 5:4(3)- and 3:2(2)- circular-coplanar precursor_mga:
    single-ellipse 38->40; then Case 2 5:4(3)- and Case 3 3:2(2)-
    circular-coplanar precursor_mga: single-ellipse 40->42; then S1L1
    5:4(3)- and 3:2(2)- circular-coplanar precursor_mga: single-ellipse
    42->44; then U0L1 4:3(3)- and 2:1(1)- circular-coplanar precursor_mga:
    single-ellipse 44->46.)
    """
    rows = _load_rows()
    counts = Counter(r.get("cycler_class", "single-ellipse") for r in rows)
    expected = {"single-ellipse": 46, "multi-arc": 252, "non-keplerian": 12}
    assert dict(counts) == expected, (
        f"Census mismatch.\n  Expected: {expected}\n  Got:      {dict(counts)}"
    )


def test_multi_arc_ids_match_allowlist() -> None:
    """The exact set of multi-arc ids matches the 240-id MULTI_ARC_ALLOWLIST ratchet."""
    rows = _load_rows()
    actual = frozenset(r["id"] for r in rows if r.get("cycler_class") == "multi-arc")
    extra = actual - MULTI_ARC_ALLOWLIST
    missing_from_actual = MULTI_ARC_ALLOWLIST - actual
    assert actual == MULTI_ARC_ALLOWLIST, (
        f"multi-arc id mismatch.\n"
        f"  In catalogue but NOT in allowlist ({len(extra)}): {sorted(extra)}\n"
        f"  In allowlist but NOT in catalogue ({len(missing_from_actual)}):"
        f" {sorted(missing_from_actual)}"
    )


def test_non_keplerian_ids_match_ratchet() -> None:
    """The exact set of non-keplerian ids matches the 9-id NON_KEPLERIAN_IDS ratchet."""
    rows = _load_rows()
    actual = frozenset(r["id"] for r in rows if r.get("cycler_class") == "non-keplerian")
    extra = actual - NON_KEPLERIAN_IDS
    missing_from_actual = NON_KEPLERIAN_IDS - actual
    assert actual == NON_KEPLERIAN_IDS, (
        f"non-keplerian id mismatch.\n"
        f"  In catalogue but NOT in ratchet ({len(extra)}): {sorted(extra)}\n"
        f"  In ratchet but NOT in catalogue ({len(missing_from_actual)}):"
        f" {sorted(missing_from_actual)}"
    )

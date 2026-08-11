# `#810`: Pluto-Charon (5,1) fixed-hc re-sweep — STABLE MEMBER FOUND

**Task:** `#810`, registered 2026-08-09 (found during `#808`, `#656`'s own recommended follow-up,
never given a number until then). Re-run the Pluto-Charon (5,1) topology from `#656`'s genuine
prograde seed with the half-crossing count held FIXED through the C-sweep (`#504`'s "sweep upward
only, hc fixed" (3,2)-positive-control convention) instead of `sweep_family_grid`'s `hc=None`
auto-redetection — the measured branch-loss mechanism (`#656`, deterministically replayed by
`#808`). This was the ONE `k2<=k1<=5` PC topology whose negative rested on a known search-method
gap; the registration's own odds read "honest odds low". **The low-odds branch happened: the
family has a genuine, gate-passing stable member.** The last asterisk on the `#504`/`#549`/`#656`
15-topology census closes with a POSITIVE, not the expected clean negative.

**Script:** `scripts/run_810_pc_51_fixed_hc_sweep.py` (phases: control / seed / up / down /
window / gate / verdict; foreground, per-step JSONL checkpointing, resumable).
**Data:** `data/found/810_pc51_fixed_hc_sweep/{sweep_up,sweep_down,sweep_window,gates}.jsonl`
(full per-step record), `docs/notes/scratch/810_pc51_fixedhc_raw.txt` (raw log).
**Evidence tests:** `tests/scripts/test_run_810_pc_51_fixed_hc.py` (3 tests, NOT slow —
standalone independent reproduction of seed + member through fresh corrector / winding / Barden /
Radau / clearance calls).

---

## The discovered member (nu=0 midpoint of the family's one stable window)

| quantity | value |
|---|---|
| Jacobi constant C | **3.167935964707279** |
| x0 (ND, rotating frame) | **-0.7058054139668293** |
| ydot0 (ND) | -0.6722667009872901 |
| period | **24.305715846921732 TU = 24.7153 days** (t_s = 87855.81 s) |
| Barden nu | 5.98e-10 (~0; window midpoint) |
| winding topology | (k1,k2) = **(5,1)**, prograde, reaches_secondary=True |
| independent Radau crosscheck | PASS, dJ = 7.8e-13 |
| min distance to Pluto centre | 5128.5 km (radius 1188.3 km) — clears |
| min distance to Charon centre | **745.5 km (radius 606.0 km) — clears by only ~139.5 km** |
| stable window (|nu|<1) | C in [3.167773862, 3.168099844], width **3.26e-4** |
| mu | 0.10876473603280369 (`cr3bp_system("Pluto","Charon")`) |

The Charon clearance is narrow (~1.23 Charon radii from the centre, ~140 km above the surface)
but PASSING under the `#660` gate (sourced Nimmo 2017 radii, margin 0) — recorded explicitly per
the `#659` Antiope lesson, for the adjudicator to weigh (an ideal-point-mass CR3BP number; a real
mission would care about this margin).

## Run stages (all foreground, chunked; total ~35 min wall)

1. **Mandatory positive control PASSED** (18.7 s): `#504`'s own `sweep_32_positive_control()`,
   UNMODIFIED — C=3.5795150, x0=-0.693198287, T=11.83346 TU, nu=-1.2e-07, topo_ok, xcheck —
   matches the committed `ross-rt-pc-cycler-32-2026` row exactly.
2. **Seed independently re-verified** (not trusted from the bullet): `#656`'s recorded seed
   (x0=-0.6685146994, C=3.05, T=28.1427 TU) reconverges at residual 5.4e-13 with winding
   EXACTLY (w1,w2)=(+5.0000,+1.0000), prograde, reaches_secondary. Its own measured
   perpendicular-crossing index (nearest T/2) is **hc=5** — the `#656` bullet's prose "hc=4"
   was the LOST (4,0) branch's index, not the seed's own. hc=5 is what the sweep holds fixed.
3. **Fixed-dC control experiment**: a plain dC=0.005 fixed-grid sweep (the literal `#504`
   convention) jumps off this branch at its FIRST step (C=3.055 → a non-prograde (5,0);
   dx0/dC ~ 2.5, seed nu=+7.7e+03 — a far more fragile family than (3,2)). The production sweep
   is therefore ADAPTIVE: dc starts at 0.005, halves on any failed/topology-losing correction
   (retrying from the last good point) down to a 1e-5 floor, doubling back on success; hc stays
   fixed at 5 throughout and winding topology is verified at EVERY step — the exact check the
   `hc=None` path could not make.
4. **Up-sweep** (C=3.05 → fold): 47 good steps. nu falls monotonically +7.7e+03 → +19.8
   (C=3.165), crosses zero in [3.165,3.170] (brentq midpoint = the member above), bottoms at
   -58.8 (C=3.190), rises again, and the branch terminates at a **measured fold at
   C=3.24603516** (nu → +3.4e+03; past it even dc=1e-5 only recaptures the unrelated (5,0)).
5. **Down-sweep** (C=3.05 → 2.90, beyond-convention diagnostic since the seed sat AT `#656`'s
   grid floor): 100+ good steps at dc~1.25e-3 (dc=2.5e-3 attempts repeatedly capture (4,0) —
   the same fragility `#656` hit, now caught and retried instead of silently followed). nu rises
   monotonically +7.7e+03 → +5.4e+04: the family continues but NO stable member below the seed.
6. **Window re-walk** (C in [3.2400, 3.24604] at dc=5e-5, 121 steps) + an ultra-fine dc=2e-6
   eigenvalue probe: the up-sweep's apparent second nu sign change near C=3.2440 is REAL-
   EIGENVALUE SIGN FLICKER at |lambda| ~ 4.4-5.4 (|nu| >= 2.33 throughout the zone), NOT a zero
   crossing — **no second stable window**. The one stable window in step 4 is the family's only
   one across the entire swept extent C in [2.90, 3.24604].
7. **Stable-window edges** (brentq on nu-+1/nu=-1 from the recorded bracket steps):
   nu=+1 at C=3.167773862 (x0=-0.705889227, T=24.312466 TU), nu=-1 at C=3.168099844
   (x0=-0.705720543, T=24.298891 TU).
8. **`#660` clearance gate**: figures above; PASS.

## Literature-novelty gate (mandatory; run live 2026-08-11)

`search/literature_check.py` signature: primary=Pluto, sequence=(Charon,),
topology_label={repeated-moon}, resonances=("5:1",). Query trail run via live web search:
"Charon cycler trajectory", "Pluto-Charon CR3BP periodic orbit cycler", "Ross Roberts-Tsoukkas
stable prograde cycler binary mass parameter (5,1)", ""Pluto-Charon" cycler stable periodic
orbit winding (5,1)". Overlapping KNOWN_CORPUS anchors surfaced and checked: Persephone
(Howard/Stern 2021, PSJ 2(2):56 — CR3BP science orbits, no cycler families), Stern
Game-Changer (DPS 2018 — flyby tour), Brozovic (satellite orbit determination). The governing
closest prior art, **grounded against the actual source** (arXiv:2606.29189 HTML fetched, per
[[feedback_ground_citations_against_content]]): Ross & Roberts-Tsoukkas 2026 tabulates ONLY
(1,1), (3,1), (3,2), (3,3) — **no (5,1) family appears anywhere in the paper** — and computes
no Pluto-Charon (mu~0.1087) members at all. Their universal-stable-subfamily CONJECTURE
predicts something like this should exist; this member is a confirming INSTANCE of that
conjecture for a family they never computed, at a mu they never addressed. JPL SSD
Three-Body Periodic Orbits catalog: Pluto-Charon is NOT one of its 7 indexed systems
(`jpl_family_check` verdict would be "not-covered" — that gate cannot adjudicate this).

**Status: not-found** — necessary-not-sufficient per the module's own discipline. NOT certified
novel; `#836`'s adjudication governs. `data/catalogue.yaml` grep confirms `our_status` context:
the only PC row is `ross-rt-pc-cycler-32-2026` ((3,2), V1, known-class member) — this (5,1) is
not a duplicate of anything we hold.

## What this changes, and what it does not

- `#656`'s (5,1) UNSETTLED adjudication is now RESOLVED — the suspected false negative was
  real: the `hc=None` auto-redetection was hiding a genuine stable member. The
  `[[feedback_bugfix_invalidates_past_searches]]` discipline (a buggy/limited search path is a
  false-negative generator) scores a full hit here.
- The 15-topology census verdict flips from "0 of 15 beyond (3,2)" to **"1 of 15: (5,1) has a
  stable member"** — PC (3,2) is no longer structurally unique at this mu; adjudication of the
  new member decides how the census prose gets rewritten.
- The `pluto-charon-kk-45-cycler-sweep-2026-07-19` empty-region stamp is NOT edited (append-only
  registry); its (5,1) row already says UNSETTLED-not-certified-empty, which remains literally
  true. NO new empty-region stamp (a positive is not an empty region). NO catalogue writeback
  from this task — `#836` (registered) owns independent re-verification, adjudication, and any
  writeback.
- Other 8 higher topologies' certified-empty negatives are untouched (their grids found no seed
  at all; nothing downstream to lose).

## Follow-ups registered

- **`#836`**: Opus+Fable adjudication + catalogue writeback decision for the PC (5,1) stable
  member (independent reproduction beyond this task's own Radau/standalone-test checks, the
  narrow ~140 km Charon clearance as an explicit adjudication axis, V-tier assignment, and the
  census-prose rewrite in `#504`/`#549`/`#656`'s records if admitted).

## Verification

`tests/scripts` full suite (incl. the 3 new evidence tests + the `preflight_search` AST ratchet
on the new script), `tests/data tests/search -q` full ratchet, `ruff check` / `ruff format
--check`, full `mypy src tests` — see the `#810` commits in `git log`.

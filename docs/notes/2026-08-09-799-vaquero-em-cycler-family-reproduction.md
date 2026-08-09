# #799: Vaquero 2013's Earth-Moon 2:1 / 3:1 "Periodic Cycler" families — reproduced by direct CR3BP continuation

**Task:** `#799`, registered 2026-08-08 (spawned by `#798`'s decisive negative). Direct CR3BP
family-continuation reproduction of the two planar resonant "Earth-Moon Periodic Cycler"
families of Vaquero (2013), Purdue Ph.D. dissertation, Sec. 4.4.7 (pp. 169-172, Figs.
4.43-4.44). Source digest: `docs/notes/2026-08-08-787-vaquero-2013-earth-moon-chapter-digest.md`
(§4). Neither Vaquero nor Casoliva print a digit-grade IC table for these families anywhere —
Fig. 4.44's `x0` axis carries no value labels — so continuation against her PROSE values
(Jacobi ranges, endpoint TOFs, four selection criteria) was the only reproduction route.

**Verdict: REPRODUCED — both families, across their FULL printed Jacobi ranges, endpoint to
endpoint, with all four of her printed endpoint TOFs recovered (0.10%-1.83% relative) and
every member passing the full convergence gauntlet including an independent-Radau
cross-check.** This is a reproduction of a published result, NOT a novel finding — nothing
here is claimed novel and the `literature_check.py` novelty gate is not in play (`#780`'s own
precedent for reproduction-scoped modules; the reproduction target IS the published record).
Digit-grade, writeback-eligible ICs now exist for both families (below + full 129-member
record in `data/found/799_vaquero_em_cycler_families/results.json`). **No catalogue writeback
performed here** (out of `#799`'s scope; registered as `#811`).

---

## Method

New module `src/cyclerfinder/search/vaquero_em_cyclers.py` (tests:
`tests/search/test_vaquero_em_cyclers.py`; driver:
`scripts/screen_799_vaquero_em_cycler_families.py`). Model: planar CR3BP at this project's
own canonical registry Earth-Moon system (`cr3bp.cr3bp_system("Earth", "Moon")`,
`mu=0.01215058439469525`, `l*=384,400` km, `t*=375,190.26` s) — Vaquero's Purdue/Howell
frame convention (Earth at `(-mu,0)`, Moon at `(1-mu,0)`) matches this project's own, so
`x0` values compare directly with no flip (unlike the Casoliva flip `#780` documents).

**Seeding (two-body, per Vaquero's own Ch. 3.5.1 method lineage).** Both families are `p:1`
INTERIOR resonances (her Eq. 3.1 convention: spacecraft `p` revs per 1 Moon rev), so the
two-body ellipse has `a = p^(-2/3)` and full nondim period `~2*pi`. The established naive
seed (`jrf.two_body_resonant_seed`, periapsis pinned AT the secondary's radius) is documented
4x project-wide as never landing on the intended topology (`TWO_BODY_SEED_LINEAGE_NOTE`), so
this module instead pins the ellipse's APOAPSIS perpendicular crossing on the `+x` (Moon)
axis and fixes the eccentricity from the target Jacobi constant via the Tisserand relation
`C ~= 1/a + 2*sqrt(a(1-e^2))`. That construction IS the geometry of Vaquero's own families
(perigee in the LEO-GEO band, apogee at the Moon's vicinity: 3:1 apogee ~0.89 LU cislunar/L1,
2:1 apogee ~1.19 LU circumlunar/L2 — exactly Fig. 4.44's plotted `x0 ∈ [~0.6, ~1.2+]` band).
Seeded at mid-range energies (2:1 at `C=2.30`, 3:1 at `C=2.80`) chosen so her printed
endpoints land exactly on the `dC=0.01` continuation grid. **Both seeds converged on the
target family on the FIRST corrector attempt** (residuals `1e-11`-`3.5e-11`; no fallback to
the R31-S-anchored strategy the dispatch note held in reserve was needed).

**Continuation.** The existing `cyclerfinder.search.cr3bp_continuation.continue_family`
(no new machinery), `half_crossings=3` (measured on the converged seeds — 6 x-axis
crossings per period — then asserted at seed time, not assumed), `ydot0_sign=-1`,
`rtol=atol=1e-13`, walked from each seed in both directions to Vaquero's own printed range
endpoints. Every kept member passes the module's full gauntlet: corrector convergence,
period bounds, non-equilibrium, dedup, Jacobi conservation `<=1e-10`, and the
independent-Radau cross-check (different integrator family than the DOP853 corrector — this
project's standing second-integrator verification pattern).

**Why 1e-13 (not the 1e-12 default):** the low-C end of both families dives to ~7,000-8,500
km perigees, where a 1e-12 DOP853 propagation leaks just past the gauntlet's 1e-10
Jacobi-conservation gate. Measured directly: the `C=2.54` 3:1 member converges at residual
`7e-13` but fails that ONE gate at 1e-12; at 1e-13 the whole range passes. (First-run
artifact worth recording per the sweep-singleton discipline: the initial 1e-12 run stopped
one step short of `C=2.54` with `gauntlet_reject` — a tolerance edge, not a family
termination.)

**Earth→Moon TOF metric.** Shortest periodic-wraparound time from a USABLE perigee (any
local `r1` minimum inside her LEO-GEO insertion band) to a genuine close lunar approach
(local `r2` minimum within 5% of the closest). This matters: at the 3:1 family's low-C end
the three perigees per period are UNEQUAL — at `C=2.54` the global perigee (7,170 km,
below the 6,558 km LEO floor) sits half a period from the Moon encounter, while the two
in-band 8,300-km-class perigees sit 4.90 d from it. Vaquero's own printed 4.90 d at `C=2.54`
is unambiguously the in-band-perigee leg — measuring from the global perigee gives 13.6 d
and would misread the family. (An endpoint artifact was also fixed: the Moon encounter for
these apoapsis-seeded orbits is at `t=0` itself, so extrema detection must be circular over
`[0, T)`, not interior-only.)

## Results

**All four printed endpoint TOFs reproduce** (the only digit-grade values Vaquero prints
for these families):

| family | endpoint | her TOF | ours | rel. err |
|---|---|---|---|---|
| 2:1 | C=1.98 | 6.39 d | 6.399 d | 0.14% |
| 2:1 | C=2.66 | 4.91 d | 4.967 d | 1.16% |
| 3:1 | C=2.54 | 4.90 d | 4.905 d | 0.10% |
| 3:1 | C=3.13 | 5.04 d | 5.132 d | 1.83% |

(The ~1-2% misses are at the high-C endpoints; her values are 3-sig-fig prints and her mu is
not stated for this section — same rounded-display caveat as everywhere else in this
lineage. The low-C endpoints match to 1e-3.)

**2:1 family: 69 members, `C ∈ [1.9800, 2.6600]` complete** (both branches stopped at
`jacobi_bound`, i.e. her printed range walked end to end). Perigee 8,409 → 71,116 km
(monotone in C), apogee ~433,500-485,000 km (beyond the Moon — circumlunar, her L2-connected
family), closest lunar approach 86,911 → 17,675 km (well inside the lunar SOI at the high-C
end). Two-body `a` 0.6418-0.6564 LU vs. the Keplerian 2:1 value 0.63 (+2-4%, the expected
CR3BP shift). Stability: **linearly stable (|nu|<1) for C <= 2.46, small unstable modes
above (|lambda| <= 5.73)** — flip between C=2.46 and 2.47.

**3:1 family: 60 members, `C ∈ [2.5400, 3.1300]` complete** (both branches `jacobi_bound`).
Perigee 7,170 → 57,910 km, apogee ~317,400-351,400 km (inside the Moon's orbit — cislunar,
her L1-connected family), closest lunar approach 33,258-66,995 km (0.57-1.15x the Moon-L1
distance — the L1-LPO-connection vicinity of her criterion 2). Two-body `a` 0.4666-0.4882
LU vs. Keplerian 0.4807. Stability: **no linearly stable member anywhere in her printed
range; |lambda| ∈ [11.3, 16.9] throughout** — small unstable modes (e.g. vs. |lambda|=2513.2
for the 4:3 family she herself excludes as a non-cycler), out-of-plane stable
(`nu_z ∈ [0.46, 0.66]`, |nu_z|<1) everywhere.

**Representative digit-grade ICs** (perpendicular `+x`-axis crossing, this project's frame,
`state = (x0, 0, 0, 0, ydot0, 0)`, registry `mu=0.01215058439469525`; full 129-member record
with per-member criteria in `data/found/799_vaquero_em_cycler_families/results.json`):

| family | C | x0 | ydot0 | T (nd) |
|---|---|---|---|---|
| 2:1 | 1.98 | 1.2139445950765162 | -1.1011440751334407 | 6.211141776252044 |
| 2:1 | 2.30 | 1.1592937721573084 | -0.9339460706615066 | 6.126282583530318 |
| 2:1 | 2.66 | 1.0338302047346954 | -0.9089334377051435 | 5.662843584779122 |
| 3:1 | 2.54 | 0.9013301668020125 | -0.8462249954358775 | 6.269604424886022 |
| 3:1 | 2.80 | 0.8768242486174663 | -0.6404271138199539 | 6.306045182360639 |
| 3:1 | 3.13 | 0.8135643069819515 | -0.25304820538336525 | 6.45496522207971 |

Worst crossing residual across all 129 members: `1.95e-12`; worst independent-Radau Jacobi
drift: `2.10e-12` (2:1 first-run figure `1.74e-11` residual at 1e-12 tolerances; final run
tighter).

## Family-identity verification (all-constraints discipline, mirroring `#798`)

1. **Jacobi range**: members span exactly her printed `[1.98, 2.66]` / `[2.54, 3.13]`,
   endpoint to endpoint, both branches terminating on `jacobi_bound` (not folds/rejects).
2. **Endpoint TOFs**: all four printed values recovered (table above) — the only digit-grade
   prose values she gives, and the strongest single identity check available.
3. **Resonance**: two-body semi-major axis within 2-4% of `p^(-2/3)` for every member (same
   bookkeeping `#798` used to place R31-S on the 3:1 resonance).
4. **Geometry**: perigee in/near her LEO-GEO insertion band (48/69 resp. 48/60 members
   strictly in-band; the high-C tails exceed the GEO ceiling — see honest caveats), lunar
   approach at SOI/L1/L2 scale, 2:1 circumlunar vs. 3:1 cislunar exactly as she describes.
5. **Fig. 4.44 x0 band** (corroborating, figure-read): our 2:1 `x0 ∈ [1.03, 1.21]`, 3:1
   `x0 ∈ [0.81, 0.90]` — inside her plotted `[~0.6, ~1.2+]` axis span, on the correct sides.
   (Contrast R31-S's `x0=-0.808`: `#798`'s non-member verdict stands — R31-S's perpendicular
   crossings are on the opposite side of the barycenter from this entire family.)
6. **Stability criterion** ("stable or possess small unstable modes", p. 170): 2:1 stable
   through most of its range, small unstable modes at the high-C end; 3:1 small unstable
   modes throughout. Her own free-transfer statement (p. 172: unstable-to-unstable at equal
   C in the `[2.54, 2.66]` overlap band) REQUIRES unstable members of both families exactly
   there — and that is exactly what we find (2:1 |lambda| 2.6-5.7 and 3:1 |lambda| ~11.3
   in the overlap band). Her Fig. 4.44 ν3D coloring spans ~0.5-2.5+; our full-period
   out-of-plane trace `2*nu_z` spans 0.93-1.32 (3:1) and 2.07-5.76 (2:1, saturating the
   colorbar at the high-C end where her figure shows the orange/red arrows).

## Honest caveats (reported, not absorbed)

- **The 3:1 family has NO linearly stable member in our reproduction** anywhere in her
  printed range. Her criterion admits "small unstable modes" so this does not break the
  reproduction (and |lambda| <= 16.9 is small in her own 4:3-contrast sense), but if her
  Fig. 4.44 ν2D arrows encode stable 3:1 members somewhere, our family disagrees there —
  unverifiable without her per-member data (the figure prints no values).
- **High-C perigees exceed her GEO insertion ceiling** (2:1: >42,164 km for C >~ 2.45; 3:1:
  for C >~ 3.05). Her printed C ranges nonetheless extend there, and her own printed
  endpoint TOFs match ours at those same energies — read as: the FAMILY spans her printed
  range; her insertion-band criterion selects the members usable from LEO, it is not a
  family-membership condition at every C. (Same criteria-vs-family-extent tension `#798`
  documented from the other side.)
- Her exact mu for this chapter is unstated; registry mu used throughout (project policy +
  `#780`'s evidence that her/Casoliva-lineage internals sit closer to the registry value
  than to displayed roundings).
- TOF at the high-C endpoints misses by 1.2-1.8% (3-sig-fig prints, unknown mu, and her
  TOF's exact departure/arrival definition is unstated — ours is in-band-perigee to closest
  lunar approach).

## Verification run

- `tests/search/test_vaquero_em_cyclers.py`: 15 tests (sourced-constant verbatim, seed
  geometry, both seeds' convergence + criteria, independent-Radau cross-check, short
  two-directional continuation smoke) — all pass.
- Full `tests/search tests/data tests/scripts -q`, `uv run ruff check .`,
  `uv run ruff format --check .`, full `uv run mypy src tests` — all clean.
- Full-range driver: `scripts/screen_799_vaquero_em_cycler_families.py` (foreground, ~8 min),
  archived at `data/found/799_vaquero_em_cycler_families/results.json`.

## Follow-up registered

- `#811`: catalogue writeback of these two reproduced families (rows, validation-level
  gates, `#797`-style class bookkeeping) — writeback-eligible now that digit-grade,
  Radau-cross-checked ICs exist inside Vaquero's own printed ranges; explicitly out of
  `#799`'s own scope.

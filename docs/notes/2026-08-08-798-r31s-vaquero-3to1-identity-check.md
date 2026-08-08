# #798: Is `braik-ross-planar-r31-s-corridor` the same orbit as Vaquero's 3:1 Earth-Moon Periodic Cycler?

**Task:** `#798`, registered 2026-08-08 during `#787`'s Vaquero 2013 Ch.4.4 digest (see
`docs/notes/2026-08-08-787-vaquero-2013-earth-moon-chapter-digest.md` §6, "The R31-S
near-coincidence"). That digest found `data/catalogue.yaml`'s already-catalogued
`braik-ross-planar-r31-s-corridor` row (Braik & Ross 2026 Table 2, "Resonant 3to1 Stable",
`C_J=3.1294`) sits within `0.001` of Vaquero 2013's own printed 3:1 "Earth-Moon Periodic
Cycler" family upper bound (`C=3.13`, dissertation Sec. 4.4.7, p.171) and asked whether this
is a genuine orbit-identity coincidence. A Keplerian (two-body ellipse) estimate in that digest
was indicative but not decisive.

**Verdict: RULED OUT, decisively.** R31-S is a genuine 3:1 interior resonant CR3BP orbit
(confirmed below), but it fails BOTH of Vaquero's own printed selection criteria at EVERY point
along its full period — not just the single apse crossing the earlier Keplerian estimate
checked — and its own perpendicular x-axis crossings sit entirely outside the plotted range of
her own Fig. 4.44. The near-coincidence in `C` is judged coincidental. No catalogue.yaml edit
made; this is a pure identity-resolution/documentation task per its own registered scope.

---

## Method

Direct propagation of R31-S's own already-catalogued IC (`mass_ratio=0.0121505842705716`,
`state_nd=[-0.8081272737956029, 0, 0, 0, 0.1389495551262821, 0]`,
`period_nd=6.267093371035164`) over one FULL period in this project's own CR3BP dynamics
(`cyclerfinder.core.cr3bp.cr3bp_eom`, DOP853, `rtol=atol=1e-13`), sampled at 400,000 points
across the period (~15,700 samples per Earth-orbit loop — far finer than needed to resolve the
three periapsis/apoapsis structure of a 3:1 resonant orbit). This directly answers the
cheap-first-step scope in `#798`'s own registration: check whether the orbit EVER (not just at
the one printed apse) comes within Vaquero's own stated bounds.

**Sanity checks on the propagation itself** (a "did I reproduce the sourced orbit correctly"
gate, per this project's orbit-closure discipline):
- Periodic closure: `|state(T) - state(0)| = 2.83e-12` (nondimensional) — confirms this IS the
  catalogued periodic orbit, correctly reproduced.
- Jacobi constant conservation over the full trajectory: range `[3.12939999999961,
  3.12940000000226]`, spread `2.65e-12` — matches the catalogued `jacobi_constant=3.1294` to all
  printed digits, confirms clean numerical integration (not comparing against a mis-propagated
  trajectory).
- Three local Earth-distance (`r1`) minima found across the full period (63,585 / 62,882 / 63,585
  km at physical scale, actually 64,461 / 62,882 / 64,461 km — see below), consistent with the
  "3" in a 3:1 spacecraft:moon resonance (three Earth loops per one lunar-period-length orbit,
  Vaquero's own convention, Eq. 3.1 p.83-84, cross-validated in `#787`'s digest).

## Findings

### 1. Close-Earth-approach criterion — fails at ALL THREE perigees, not just the checked apse

Vaquero's own printed criterion (Sec. 4.4.7, p.169-170, quoted in `#787`'s digest): insertion
from LEO with `180 km ≤ r-r_E ≤ 35,786 km` (GEO altitude), i.e. Earth-center distance in
`[6558, 42164]` km (`R_Earth=6378.137` km + the stated altitude band; 42,164 km is a derived
sum, not itself a printed number).

Full-period scan of R31-S's Earth-center distance (`r1`):

| Perigee # | `t` (days) | `r1` (km) | vs. Vaquero's 42,164 km ceiling |
|---|---|---|---|
| 1 | 4.611 | 64,461.49 | +52.9% over |
| 2 | 13.626 | 62,882.45 (global min) | +49.1% over |
| 3 | 22.640 | 64,461.49 | +52.9% over |

**All three perigees exceed the ceiling by ~50%.** This is stronger than the earlier Keplerian
estimate suggested (`~63,600` km at the one checked apse) — the full-orbit scan confirms it is
not an artifact of which apse was checked; every perigee in the orbit fails.

### 2. Close-Moon-approach criterion — fails everywhere, by a wide margin

Vaquero's criterion: surface contact or a small-amplitude L1/L2 LPO connection. Computed the
Earth-Moon L1/L2 positions with this project's own `cyclerfinder.core.crnbp.cr3bp_collinear_point`
(1D Newton-Raphson root of the pseudo-potential gradient) at R31-S's own `mu`:

- `L1` at `x=0.836915` (nd) → Moon-L1 distance `58,019` km
- `L2` at `x=1.155682` (nd) → Moon-L2 distance `64,515` km

Full-period scan of R31-S's Moon-center distance (`r2`) found four local minima, all in the
288,000-328,500 km range (global min `288,002.28` km at `t=5.247` d). That is:
- **4.96x** the Moon-L1 distance, **4.46x** the Moon-L2 distance — R31-S never even approaches
  the vicinity of either collinear libration point, let alone the Moon's surface (`R_Moon=1737.4`
  km — the closest approach is >165 lunar radii away).

### 3. Family-membership check (not just criterion-failure) — is it at least a continuation of the same 3:1 family, past her selection cutoff?

Two independent checks, both against the "coincidence in `C` only" reading:

**(a) Two-body semi-major-axis consistency.** `(r_peri,min + r_apo,max)/2 = (62,882.45 +
305,973.44)/2 / 384,400 = 0.479781` LU, matching Casoliva's own printed 3:1-resonance two-body
value `a_s = (q/p)^(2/3) = (1/3)^(2/3) = 0.480750` LU to 0.2%. **This confirms R31-S genuinely
IS a 3:1 interior resonant orbit** by the same two-body bookkeeping Vaquero/Casoliva use to seed
their own families — it is on the right resonance, just not admitted by her further filters.

**(b) Smoothness argument against "continuation past her stated cutoff."** Vaquero's own printed
3:1 range is `C ∈ [2.54, 3.13]`; R31-S sits at `C=3.1294`, only `ΔC=0.0006` past her stated
endpoint. If R31-S were simply her own family continued a hair further in `C`, the perigee should
change smoothly and by a small amount across that `ΔC`. Instead perigee would have to jump from
`≤42,164` km (satisfying her filter at `C=3.13`) to `62,882` km at `C=3.1294` — a ~20,700 km
jump over a `0.02%` change in Jacobi constant. That is not smooth family continuation; it argues
R31-S sits on a different branch/family member, not simply just past her cutoff.

**(c) `x0`-axis corroboration (figure-read, corroborating only, not primary).** `#787`'s own
Fig. 4.44 vision read records the plotted `x0` axis spanning `~0.6` to `~1.2+` (positive,
interior/near-Moon side) for BOTH the 2:1 and 3:1 families. A full-period scan of R31-S found
exactly two perpendicular (`vx=0`, `y=0`) x-axis crossings: `x=-0.808127` (the catalogued IC
itself) and `x=-0.175737` (the opposite perigee-symmetric crossing). **Both are negative** —
neither lands anywhere near Vaquero's plotted `[0.6, 1.2+]` range; R31-S's own perpendicular
crossings sit on the OPPOSITE side of the barycenter from every member she plots. This is a
figure-read corroboration (own vision extraction, not a printed number) so it is not weighted
as heavily as the printed prose criteria above, but it points the same direction.

## Verdict

**Ruled out, not merely "argues against."** R31-S:
1. Fails Vaquero's close-Earth-approach criterion at all three perigees in its full period
   (~50% over the GEO ceiling, not just at one checked apse).
2. Fails Vaquero's close-Moon-approach criterion everywhere in its full period (~4.5-5x the
   Moon-L1/L2 distance at closest approach).
3. Is genuinely on the 3:1 resonance (semi-major-axis check confirms this), but the perigee
   discontinuity implied by treating it as "her family, ΔC=0.0006 further" is not physically
   smooth — it looks like a different family member/branch, not a simple continuation past her
   stated cutoff.
4. Its own perpendicular x-axis crossings sit outside the plotted range of her own Fig. 4.44
   (corroborating, not primary).

The near-coincidence in Jacobi constant (`C=3.1294` vs. her printed `3.13` upper bound) is judged
**coincidental** — not evidence of orbit identity or even same-family-different-member status in
any useful sense for a writeback. This resolves `#798` without needing the more expensive escalation
step (direct CR3BP family continuation of Vaquero's own 3:1 family) — the direct full-period
propagation was decisive on its own, not merely "inconclusive."

**No catalogue.yaml edit made** — this is a pure identity-resolution task per its own registered
scope; nothing here changes any existing row's `our_status` or invalidates the existing
`braik-ross-planar-r31-s-corridor` row (which remains correctly sourced to Braik & Ross 2026,
unrelated to Vaquero's work).

## New task registered

`#799` (registered in `data/OUTSTANDING.md`): a dedicated direct CR3BP family-continuation
reproduction of Vaquero's own 2:1/3:1 resonant-cycler families (two-body-seeded per her Ch.3.5.1
method, continued in `C` across her own stated ranges) to produce digit-grade ICs for these
families — the concrete, now-unblocked path to closing `#797`'s catalogue-class gap for the
Vaquero-specific families (distinct from `#797`'s own Casoliva-Table-3-writeback scope, which
already has digit-grade source data and doesn't need continuation).

## Verification

Pure investigative task; no catalogue.yaml or source-code changes. Ran the OUTSTANDING.md
structural ratchets after editing that file:

```
uv run pytest tests/data/test_outstanding_structure.py tests/data/test_outstanding_header_body_consistency.py -q
```

# `#780`: Earth-Moon planar CR3BP resonant/Lyapunov family gate — Casoliva 2008/2010

**Task:** `#780`, Stage A of a new discovery lane for this project's own flagship Earth-Moon
system (322 of 383 catalogue rows), which — despite that scale — had never had the
reproduce-a-published-family-table-then-connections pipeline this project already built twice
(Saturn-Titan `#765`/`#767`, Neptune-Triton `#776`/`#777`) applied to it. Dispatched by the user
("start the two paths in parallel, assign new task numbers"), anchored on:

1. Casoliva, J., Mondelo, J. M., Villac, B. F., Mease, K. D., Barrabés, E. & Ollé, M., "Two
   Classes of Cycler Trajectories in the Earth-Moon System," *JGCD* 33(5), 2010, pp. 1623-1640,
   DOI `10.2514/1.46856` ("2010").
2. Same authors, "Families of Cycler Trajectories in the Earth-Moon System," AIAA 2008-6434
   ("2008", the JGCD paper's conference precursor).

Both filed at `cyclers_pdf/papers/casoliva-mondelo-villac-mease-barrabes-olle-*`. Digest:
`docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md`.

**Code delivered:** `src/cyclerfinder/search/earth_moon_resonant_families.py` (new module, a
thin sibling of `saturn_titan_resonant_families.py`/`neptune_triton_resonant_families.py`,
reusing `cyclerfinder.search.cr3bp_periodic` (`correct_periodic`, `correct_symmetric_fixed_jacobi`,
`barden_stability`, `crosscheck_periodic`), `cyclerfinder.search.cr3bp_continuation.continue_family`,
and `cyclerfinder.search.jovian_resonant_families` (`two_body_resonant_seed`,
`converge_candidate`) directly — no reimplementation of shared machinery) +
`tests/search/test_earth_moon_resonant_families.py` (67 tests, all passing). Both clean under
`ruff check`, `ruff format --check`, and `uv run mypy src tests`.

---

## What was vendored

**Class 1** (Sec. IV both papers, high-energy near-Keplerian p-q resonant orbits): Table 3
(2010, p.1630), **ALL 16 rows**, verbatim, at Casoliva's own printed precision (10-11
significant figures). Casoliva's own footnotes flag 7 of the 16 as either not satisfying their
own p-q resonance relation (`satisfies_resonance=False`: 1-2a, 1-2b, 2-1c, 3-2a, 7-3d) or
flying through the Earth's own radius (`exists_in_em_system=False`, physically unusable as a
real cycler even though it is a mathematically valid PCR3BP periodic orbit: 2-1c, 2-1d, 3-2d,
7-3d) — both flags carried through verbatim, not silently dropped. The 9 "clean" rows (satisfy
resonance AND exist in the EM system) split 5 stable / 4 unstable per Casoliva's own printed
stability index: **stable** 1-2c, 1-2e, 2-1a, 3-2c, 7-3a; **unstable** 1-2d, 2-1b, 7-3b, 7-3c —
Casoliva's own documented "at least one unstable and one stable cycler per resonance class"
finding, confirmed directly. This is the "separate, uncatalogued stable family" the `#725`
digest flagged.

**Class 2** (Sec. V both papers, low-energy L1-Lyapunov homoclinic cyclers): Table 4 (2010,
p.1635), the **He1-family golden connection anchor** at energy `h = -1.45016232260699`
(printed to 14-15 significant figures) — the L1 Lyapunov periodic orbit itself (the BASE
OBJECT the homoclinic connection is built around), period `T = 6.706878522271349` (nondim,
29.1640 days per the paper's own text), unstable eigenvalue `108.5966557497375`. The
homoclinic CONNECTION itself (He1-4, Hm1-2 discrete points) is explicitly out of scope this
stage — see "Registered follow-up" below.

## Coordinate convention: a real hazard, resolved and verified, not assumed

Casoliva's own convention places Earth at `(+mu, 0)`, Moon at `(mu-1, 0)` — the OPPOSITE side
from this project's own `cyclerfinder.core.cr3bp` convention (primary at `(-mu, 0)`, secondary
at `(1-mu, 0)`). The two are related by an exact rigid 180° rotation about the barycenter:
`(x, y, vx, vy) -> (-x, -y, -vx, -vy)`. This was NOT just asserted — it was verified two
independent ways, and the second one caught a real mistake before it shipped:

1. **Jacobi-constant round-trip** (necessary but not sufficient): the flip-transformed IC
   reproduces the row's own printed `C_J` to 3.2e-7 relative (1-2c). But `C_J` only depends on
   `vx^2+vy^2` (velocity magnitude), so this check is BLIND to a wrong velocity sign — a wrong
   sign can still match `C_J` almost exactly by coincidence (measured directly during this
   task's own scoping: a wrong-sign Class 2 candidate matched the target `C_J` to 1.8e-6
   relative yet was NOT the right orbit at all).
2. **Direct-propagation periodicity** (decisive): propagating the transformed IC for the row's
   own printed period and checking closure. A wrong sign produces dramatic non-closure
   (chaotic/close-encounter blowup, residual O(1) or worse) within one period; the correct
   transform closes all 17 rows (16 Table 3 + the Class 2 anchor) to 1e-4-1e-6 (position/
   velocity-scale units) at this project's own registry mu.

An earlier working draft of this task's own scoping mis-transcribed Table 4's 4th
canonical-momentum component (`p_y`) as `+1.43284467834384` when the correct reading is
`p_y = -1.43284467834384` (the printed digits render without the colon-substituted minus sign
Table 3's negative numbers use, an inconsistency traced to this module's own OCR tooling on
Table 4's specific typesetting, not a defect in Casoliva's own PDF). The wrong-sign reading
matched the target Jacobi constant to 1.8e-6 relative (via the degenerate-velocity-magnitude
coincidence above) but propagated into a close-encounter blowup within one period. The correct
reading — confirmed by both direct-propagation closure (2.5e-6) AND, decisively, by the
recovered orbit's own dominant monodromy eigenvalue matching Table 4's printed
`108.5966557497375` to **5.3e-8 relative** — resolved it. (This finding surfaced during an
`advisor()` review call mid-task; caught before any vendored constant was committed.)

## mu / l* / t*: registry, per explicit task instruction — and numerically the better choice

Unlike the Saturn-Titan/Neptune-Triton modules (which use the SOURCE paper's own displayed mu
for those systems, which this project has no other canonical anchor for), this module uses
`cyclerfinder.core.cr3bp.cr3bp_system("Earth", "Moon")` throughout — the SAME registry mu/l*/t*
already anchoring 322 catalogue rows — NOT Casoliva's own displayed `mu_EM = 0.0121529529`
(2010) / `0.01215` (2008), per explicit task instruction ("do NOT introduce a new one"). This
turns out to be numerically better, not merely policy-compliant: the Class 2 golden anchor's
own printed energy reproduces the registry-mu Jacobi constant to **1.3e-9 relative**, TIGHTER
than at Casoliva's own displayed 7-significant-figure mu (1.8e-6 relative) — evidence Casoliva's
internal computation used a mu closer to this project's own DE440-registry value than to her
own paper's rounded display.

## The Class 1 gate: 12 of 16 rows fully pass; 4 honest misses, all on stability-index only

Gate criteria per row: (a) `correct_periodic` (general full-state Newton, since Casoliva's own
printed IC point is a generic y=0 Poincaré crossing — her own text: "we have not used this
property [the perpendicular crossings] in this paper" — NOT the perpendicular symmetric
crossing the standard symmetric corrector needs) converges; (b) recovered `x0`/period/Jacobi
AND the stability index `k` reproduce Casoliva's own printed Table 3 values within 1e-2
relative; (c) internal trace-identity-vs-eigenpair-selection AND independent-Radau-integrator
cross-checks both agree.

**Stability-index finding (the genuinely new piece of work this task did NOT anticipate
needing).** A naive full-period in-plane eigenvalue (`k_par = trace(Phi4_xy) - 2`) reproduces
Casoliva's own printed `k` for only 5 of 16 rows. Re-reading Casoliva's own Eq. 6-8 text
(2010 p.1626) resolved it: her printed "k" is `max(|kappa_par|, |kappa_perp|)`, the LARGER of
the in-plane AND out-of-plane (vertical, `z`/`vz`) stability indices — even for these strictly
PLANAR orbits, out-of-plane linear stability is a distinct, generally-larger-magnitude
quantity. Computing `k_perp` from the DECOUPLED 2x2 `(z,vz)` monodromy block (for a planar
orbit `z=vz=0` is an invariant subspace, so `trace(Phi_z) = kappa_perp + 1/kappa_perp` exactly,
no eigenvalue-pairing needed) and taking the signed max of `{k_par, k_perp}` reproduces
Casoliva's own printed k to **8.6e-8 to 1.7e-3 relative for 12 of the 16 rows**:

| Row | k_par | k_perp | k (signed max) | Printed k | rel. err | Passed |
|---|---|---|---|---|---|---|
| 1-2a | 0.536 | 1.997 | 1.9971 | 1.9971178106 | 6.5e-7 | ✓ |
| 1-2b | 0.534 | 1.998 | 1.9976 | 1.9975766156 | 4.6e-7 | ✓ |
| 1-2c | -1.029 | 2.000 | 2.0000 | 1.9995914782 | 8.6e-8 | ✓ |
| 1-2d | 4.857 | 1.999 | 4.8573 | 4.8578794987 | 1.2e-4 | ✓ |
| **1-2e** | -4.191 | 2.000 | -4.191 | 1.9997758212 | 3.10 | **✗ (k only)** |
| 2-1a | 1.283 | 1.926 | 1.9256 | 1.9255637447 | 2.1e-5 | ✓ |
| 2-1b | 1.513 | 2.037 | 2.0374 | 2.0373978625 | 7.0e-6 | ✓ |
| 2-1c | 103.05 | -32.80 | 103.05 | 102.8780940137 | 1.7e-3 | ✓ |
| 2-1d | 119.00 | -31.88 | 119.00 | 118.9539595702 | 3.9e-4 | ✓ |
| **3-2a** | -3.222 | 1.511 | -3.222 | 1.5113469213 | 3.13 | **✗ (k only)** |
| 3-2c | -1.246 | 1.875 | 1.8753 | 1.8751541804 | 6.8e-5 | ✓ |
| 3-2d | -0.334 | 2.001 | 2.0008 | 2.0008005398 | 4.2e-8 | ✓ |
| **7-3a** | -4.965 | -1.299 | -4.965 | -1.2990228617 | 2.82 | **✗ (k only)** |
| 7-3b | 57.356 | -2.270 | 57.356 | 57.3519357466 | 7.4e-5 | ✓ |
| 7-3c | 57.043 | -2.253 | 57.043 | 57.3519357463 | 5.4e-3 | ✓ |
| **7-3d** | -12.17 | -1.633 | -12.17 | 4.4005195470 | 3.77 | **✗ (k + IC)** |

Three of the four misses (1-2e, 3-2a, 7-3a) still reproduce IC/period/Jacobi to 1e-5-2e-5
relative — genuinely the same orbit, the miss is on the stability index alone, unexplained
(considered and not confirmed: a return-map/sub-period stability convention, or a different
eigenvalue-pair selection specific to the asymmetric printed IC point). The fourth, 7-3d, also
misses on IC (`x0` off by a factor of ~3.8) — it separately carries BOTH footnotes (does not
satisfy its own resonance relation AND flies through the Earth), the single most degenerate row
in the table, so this broader miss is less surprising. This is reported plainly, not smoothed
over — an unresolved 3-row (4 counting the degenerate one) gap in an otherwise strong
reproduction.

Internal cross-checks: every converged row's trace-identity `k_par` agrees with the direct
eigenpair-selection `k_eig` to relative tolerance (after fixing an initially-too-tight absolute
threshold that flagged 2 rows' large-eigenvalue near-agreement as a "disagreement" — see module
`git blame`/history), and every row's independent-Radau-integrator cross-check
(`crosscheck_periodic`) passes.

**Convergence note.** `correct_periodic` is an undamped min-norm Newton iteration; the two
near-Earth-singular rows (2-1c, 2-1d, both footnote-`d`) need the full `max_iter=1000,
tol=1e-8` to close — all 16 rows converge at these settings (module default).

## The Class 2 gate: passes on its primary (eigenvalue) criterion, small honest secondary miss

| Quantity | Recovered | Casoliva (Table 4) | rel. err |
|---|---|---|---|
| x0 (nondim) | 0.6280674139348356 | 0.6280674149446867 (derived) | 1.6e-9 |
| ydot0 (nondim) | 0.8047772671762206 | 0.8047772633991533 (derived) | 4.7e-9 |
| period (nondim) | 6.706878643958127 | 6.706878522271349 | 1.8e-8 |
| **unstable eigenvalue** | **108.59666149770102** | **108.5966557497375** | **5.3e-8** ✓ |
| period (days) | 29.12448558936166 | 29.1640 | 1.35e-3 |

The eigenvalue (a genuinely non-degenerate quantity, unlike Jacobi-constant-only checks) is
the primary reproduction criterion and passes at 5.3e-8 relative — six orders of magnitude
inside the 1e-4 gate tolerance. Barden half-period stability and the independent
`_planar_floquet` full-period cross-check agree. The period-in-days figure is reported as a
small, honest, SECONDARY discrepancy (0.135% relative) — NOT attributable to the registry-vs-
Casoliva `t_s` difference alone (that is only ~9e-6 relative, checked directly: converting the
same recovered nondim period through Casoliva's own stated `omega_EM = 2.66529e-6` still gives
~29.1248 days, not 29.1640) — unresolved this task, does not gate `passed`.

## Gate item (d): two-body-seed lineage check — clean honest negative (4th confirmation)

`jrf.two_body_resonant_seed(p, q, x0_sign=-1)` converges cleanly at its own natural Jacobi
constant for all 4 distinct (p,q) pairs in Table 3 ((1,2), (2,1), (3,2), (7,3)) via
`jrf.converge_candidate` (residuals ~1e-11 to 1e-14) — but EVERY ONE lands on a
`period/2pi ~= 0.99` orbit, not its own labeled multi-period p-q resonance, and none lands
near its nearest same-(p,q) vendored row (measured: (1,2) → x0=-0.840 vs 1-2e's own
x0=+2.160; (2,1) → x0=-1.346 vs 2-1a's +0.736; (3,2) → x0=-1.169 vs 3-2a's +0.819;
(7,3) → x0=-1.484 vs 7-3b's +0.840, all project-frame). The naive two-body-resonant-ellipse
construction converges to a genuine periodic orbit every time but never identifies the correct
topology or magnitude of any Table 3 member — the same qualitative finding Anderson & Lo 2011
(Jovian), Vaquero 2013 (Saturn-Titan), and Miceli & Bosanac 2026 (Neptune-Triton) each document
for their own analogous naive attempts. This is the 4th confirmation of this project-wide
pattern.

## Gate item (e): explicitly narrow-scoped, per the dispatch's own instruction

The Barrabés-Mondelo-Ollé homoclinic-CONNECTION continuation algorithm (2010 Eq. 20) is
explicitly OUT OF SCOPE this task (dispatch instruction: "Do NOT build the homoclinic-
connection stage or the Barrabés-Mondelo-Ollé continuation-in-energy primitive in this
dispatch"). What WAS built and run: a plain natural-parameter (Jacobi) continuation of the
Class 2 Lyapunov ORBIT family from the He1 golden anchor, via the already-existing
`cyclerfinder.search.cr3bp_continuation.continue_family` — a categorically different,
much simpler operation than the excluded connection-continuator. A 5-step smoke test
(`d_jacobi=0.01`) produces 6 gauntlet-passing family members with monotonically increasing
Jacobi constant and no fold/rejection — the family persists smoothly near the anchor, as
expected. This is reported as smoke-test evidence only, NOT a reproduction (no printed IC/
eigenvalue exists at these intermediate energies to reproduce against). Class 1's 16 rows are
each isolated one-per-(p,q,label) points in Table 3 (not multi-member sequences) — no
continuation is attempted on them.

## Literature novelty gate

Not triggered — this task is pure reproduction of Casoliva's own published tables, explicitly
scoped away from anything novel (`search/literature_check.py` not invoked; nothing here is
framed as a discovery).

## Verification

1. `uv run pytest tests/search/test_earth_moon_resonant_families.py -v` — 67 passed.
2. `uv run pytest tests/data tests/search -q` — full ratchet suite (see commit log for result;
   run serialized against a sibling `#781`/`#782` background test process per standing
   discipline).
3. `uv run ruff check .` / `uv run ruff format --check .` — clean.
4. `uv run mypy src tests` — 832 source files, clean.

## Registered follow-up

**`#783`** — Earth-Moon homoclinic-connection stage (Class 2's actual He1-4/Hm1-2 discrete
connection points) + the Barrabés-Mondelo-Ollé (2009, Nonlinearity, already acquired) numerical-
continuation-of-homoclinic-connections algorithm (2010 Eq. 20), mirroring exactly how `#767`
followed `#765` for Saturn-Titan and `#777` followed `#776` for Neptune-Triton. Golden targets:
Table 4 (continuation variables, already vendored as `HE1_*` constants in this module),
Tables 5-6 (pericenter/apocenter flight-times and orbital elements for the He1 connection at
`h=-1.450162`), and the LEO-rendezvous ΔV figures (`HE1_LEO_DV_MPS_2008`=703 m/s,
`HE1_LEO_DV_MPS_2010`=717.5 m/s, both already vendored). This module's `recover_he1_lyapunov`
supplies the base Lyapunov p.o. the connection stage needs as its starting object. Registered
in `data/OUTSTANDING.md` as "registered, NOT dispatched" — not dispatched as part of this task.

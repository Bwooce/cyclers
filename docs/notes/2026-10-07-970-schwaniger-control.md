# #970: Schwaniger 1963 periodic free return re-closed with the project corrector; positive control

Task `#970` (earthmoon-opus, 2026-10-07). Source: NASA TN D-1833 sec. III.E; digest
`docs/notes/2026-10-07-digest-schwaniger-1963-earth-moon-symmetrical-free-return.md`.
Script `scripts/run_970_schwaniger_control.py`, output `data/970_schwaniger/control.json`,
test `tests/core/test_970_schwaniger_control.py`, draft row
`data/970_schwaniger/catalogue_row_draft.yaml` (NOT inserted into the catalogue).

## 1. Result

The orbit re-closes with the project corrector at both mass ratios. The defining condition is
Schwaniger's own: the perigee is at r = 6555 km ("100 nautical mile altitude or 6555 km radius",
sec. II.C) and lies on the Earth-Moon line. In the one-parameter family of x-axis-symmetric
orbits, that is the member whose half-period perpendicular crossing is at r1 = 6555 km.

| quantity | mu = 0.012150 (digest) | registry mu = 0.01215058439 | paper |
|---|---|---|---|
| x0 (periselene start) | 0.9821202024620871 | 0.9821193418196401 | - |
| ydot0 | -2.4712785397261854 | -2.4712777483864743 | - |
| C (project, no mu(1-mu)) | 1.0854156154 | 1.0854167287 | - |
| C with mu(1-mu) | 1.0974179929 | 1.0974196763 | - |
| period | 6.0015025963 TU = 625.474 h = 26.0614 d | 6.0015027496 TU = 625.474 h | "about 650 hours" |
| perigee | 6555.000 km (176.86 km alt.) | 6555.000 km | 6555 km (definition) |
| periselene | 2202.534 km (465.13 km alt.) | 2202.640 km (465.24 km alt.) | "about 2150 km" |
| apogee | 531,167 km (1.382 L) | 531,167 km | - |
| perigee speed | 10,959.95 m/s Earth-relative inertial; 10,977.42 m/s rotating | same | - |
| sense | retrograde about the Earth (inertial h < 0) | same | "counter-rotational" |
| y = 0 crossings per half period | 2 (apogee side, then perigee) | 2 | - |

Units for the digest-mu column: L = 384,400 km, TU = 375,190.39 s (GM_E = 398,600.4 km^3/s^2,
the digest's choice). The digest's independent script gave 2202.5 km, 625.47 h, C = 1.08541562
and 10,960.0 m/s; all agree.

**The brief's and OUTSTANDING's C = 1.09727 is not this orbit.** 1.09727 = 1.085263 + mu(1-mu),
and 1.085263 is the Jacobi constant of the 6-digit rounded initial state printed in the digest.
That rounded state corrects (at its own C) to a neighbouring member with perigee 6560.2 km. The
6555 km member has C = 1.0854156 (1.0974180 with the mu(1-mu) term). A 5 km perigee shift from
6-digit rounding shows how sensitive the orbit is; the draft row stores the full-precision state.

Against the paper: periselene +2.4%, period -3.8%. Both printed values are graph readings
("about"). That puts the row at V0, as the digest proposed (outside the 1e-2 V1 tolerance).

**Floors.** Perigee altitude 176.86 km is 23 km BELOW the Earth registry floor of 200 km
(`PLANETS["E"].safe_alt_km`, Russell 2004). The lunar pass (465 km altitude) clears the 100 km
Moon floor. So the Schwaniger member itself is not a both-floors cycler-class candidate; the
`#997` perigee continuation tabulates the members that are.

## 2. (a) Corrector, regularisation, monodromy

- Corrector: `cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi` (fixed C, unknown
  x0, ydot0 from C, Newton on xdot at the 2nd y = 0 crossing; tol 1e-12), with an outer secant
  on C for r1(T/2) = 6555 km. Converges in 2-3 secant steps of 3 Newton iterations each; seeded
  from the digest's 6-digit state. Crossing residual |xdot(T/2)| = 9e-13 (digest mu), 5e-13
  (registry).
- **Regularisation: none.** The corrector integrates the CR3BP and its variational equations
  with DOP853 (rtol = atol = 1e-12) in physical time, at both primaries. The 177 km perigee
  (r1 = 0.0171 L) and the 465 km perilune (r2 = 0.0057 L) are well within the plain integrator's
  range; the `#652` close-encounter guard (r2 < 1e-5) does not fire.
- What core/ has (checked by grep for "kustaanheimo", "levi", "regulari"):
  - `core/cr3bp_regularized.py`: Sundman time transform dt/ds = r1 r2 (or r1, or r2). It is a
    time transformation only; its own docstring says the force law is not regularised. No STM.
  - `core/cr3bp_ks.py` + `core/ks.py`: Kustaanheimo-Stiefel regularisation centred on ONE body,
    the Moon (`MoonCentredCR3BP`), with a fixed-physical-time 6 x 6 transition matrix
    (`#928`). Passes near the Earth are not regularised there ("one centre only").
  - No Levi-Civita, Birkhoff or Thiele-Burrau transform; nothing regularises both primaries at once.
- Monodromy: full-period variational equations, `stm_mode="fixed_path"` (Pellegrini-Russell
  2016 mitigation; lambda ~ 500). No finite differences.

| | mu = 0.012150 | registry |
|---|---|---|
| det M (6 x 6) | 1 - 1.7e-8 | 1 - 2.3e-8 |
| planar eigenvalues | 525.328, 1.00079, 0.99921, 0.0019036 | 525.302, 1.00086, 0.99914, 0.0019037 |
| b_h = tr(M4) - 2 | 525.330 | 525.304 |
| Barden half-period 2 nu | 525.330 (nu = 262.665) | 525.304 (nu = 262.652) |
| vertical eigenvalues | -182.354, -0.0054838 | -182.345, -0.0054841 |
| b_v = tr(Mz) | -182.359 | -182.351 |

Strongly unstable in the plane (b_h >> 2) and vertically unstable (|b_v| >> 2, negative: a
period-doubling-type vertical eigenvalue). The trivial pair is recovered to 8e-4; with lambda =
525 that is the conditioning floor of a 1e-12 integration. The digest's finite-difference
estimate (largest eigenvalue about 5e2) is confirmed: 525.3.

## 3. (b) Positive control for a both-primary regularised corrector

Independent re-propagation of the corrected state over one period (coordination rule: a
different path, constants from the reference model):

| integrator | closure dr | closure dv | other |
|---|---|---|---|
| plain DOP853 (`core.cr3bp.propagate`) | 7.6e-7 km | 3.5e-7 m/s | Jacobi drift < 1e-10 |
| Moon-centred KS (`core.cr3bp_ks.propagate_ks`, with STM) | 4.9e-6 km | 2.6e-6 m/s | b_h agrees to 2.0e-8 relative; regularised Hamiltonian max 1.3e-11; r_min = 2202.534 km |
| Sundman r1 r2 (`core.cr3bp_regularized`, t_stop = T) | 8.8e-6 km | 3.1e-6 m/s | |

(digest-mu column; registry: 1.1e-5 km KS, 8.7e-6 km Sundman, b_h 2.9e-8.)

`tests/core/test_970_schwaniger_control.py` (about 1.5 s) pins the control:
- sourced: the perigee crossing is perpendicular at r1 = 6555 km between the Earth and the Moon
  (Schwaniger's definition); periselene and period within 5% of "about 2150 km" and "about 650
  hours";
- structural: full-period closure by all three integrators, Jacobi conservation, det M = 1, the
  trivial pair, b_h > 2 and |b_v| > 2, and the Barden half-period form equal to the full-period
  trace;
- cross-code (labelled, not a golden): 2202.5 km and 625.5 h from the digest's separate script.

**How strong a control is it?** Weak for regularisation as such: the plain integrator handles
it. Its value for `#948` R4 is as an agreement case at the real Earth-Moon mass ratio with a
close pass at each primary in one orbit: a both-primary-regularised corrector must reproduce the
table above (C, T, perigee, periselene, b_h, b_v). The stress cases are elsewhere: the `#997`
perigee continuation continued below the floors toward the Earth's surface, and the Bruno-Varin
collision orbits (`#992`). I did not build a two-centre regularised corrector here; that is `#948`'s
deliverable. (A straightforward route, noted for `#948`: an Earth-centred `KSModel` subclass
mirroring `MoonCentredCR3BP`, with the switch distance of Aarseth's criterion, composing the
fixed-time 6 x 6 segment matrices; the existing KS transition matrix is already fixed-time, so
segments compose by multiplication.)

## 4. (c) Draft catalogue row

`data/970_schwaniger/catalogue_row_draft.yaml`: id
`schwaniger-1963-em-cislunar-retrograde-periodic-free-return`, source literature, our_status
known-reproduction, orbit_class cycler, model cr3bp, V0, registry mu, full-precision state,
first_published Schwaniger 1963 (NTRS 19630007117), priority_date 1963-06-01, verbatim quotes
checked on the OCR text. It validates against `data/catalogue.schema.json` (jsonschema). Data gaps:
unstated mass ratio; the unresolved Fig. 9 speed offset (digest); the Earth-floor flag on
orbit_class. NOT inserted into `data/catalogue.yaml`; the lead schedules that and the full ratchets.

## 5. Verified vs assumed

- Verified: closure (three integrators), perpendicularity, the perigee condition, the monodromy
  identities, agreement with the digest's independent script, verbatim source quotes.
- Assumed: the mass ratio (Schwaniger states none; two values reported); the reading that his
  "periodic trajectory" is the 6555 km member (his own sec. II.C condition).

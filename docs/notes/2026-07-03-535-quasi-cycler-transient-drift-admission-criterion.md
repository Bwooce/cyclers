# #535 admission criterion: transient-drift-phase `quasi_cycler` candidates

Settled in writing BEFORE any sweep code is built, per the #339-style-criterion-trap
discipline (`data/OUTSTANDING.md`'s own repeated warning: decide the admission rule first,
or a search's "found N" / "found none" verdict is meaningless because nobody agreed what
counts). This note is the single source of truth #535's search code must implement exactly.

## Why a new criterion is needed

Every prior encounter-counting criterion in this project (`#523`, `#527`) was written for
CERTIFIED PERIODIC ORBITS: the trajectory closes exactly, so "n returns" is just the orbit's
own period count over a span. #535's targets are explicitly NOT periodic-orbit family
members — they are genuinely aperiodic/chaotic trajectories (a captured minimoon's
drift-then-libration transition; separatrix-adjacent chaotic transport near a resonance).
"How many times did it return" has no periodic-orbit closure to lean on; it must be counted
directly from the propagated trajectory's own crossing pattern.

## 1. Encounter definition (reused, not reinvented)

A trajectory is IN ENCOUNTER with the target body at any instant its distance from that body
is less than the body's Hill-sphere radius `R_hill = a * (mu/3)^(1/3)` — the SAME measure of
"gravitationally relevant close approach" already adopted by `#523` (Earth) and `#527`
(Jupiter). No new distance threshold; this one is already justified project-wide and reusing
it keeps `quasi_cycler` and `cycler`/`resonant_po` results comparable.

## 2. Return definition: one continuous residency = one return

A single **maximal continuous time interval** during which the trajectory stays inside the
Hill sphere counts as ONE return, regardless of how many local periapsis wiggles occur inside
that interval (a close quasi-satellite episode can graze in and out of strict Hill-sphere
containment on a short sub-timescale without that being a physically distinct "different"
close approach — e.g. 2006 RH120's ~0.7-year captured episode is described in the source
literature as ONE episode, not a train of separate encounters). Concretely: find every
`(t_enter, t_exit)` pair where distance crosses below `R_hill` at `t_enter` and back above it
at `t_exit`; each pair is one candidate return.

## 3. Minimum separation between returns: 1 year

Two candidate returns (per §2) count as DISTINCT returns only if they are separated by at
least **1 year** (`t_enter` of the second minus `t_exit` of the first) of the trajectory
remaining outside the Hill sphere. Rationale: Earth's own orbital period is 1 year, so any
gap shorter than that reflects numerical/short-timescale noise in a grazing crossing, not a
genuine new approach; the real physical return timescale this session already measured for
the horseshoe-libration mechanism is ~4.6-9.2 years (#523's positive-control integration), an
order of magnitude longer than this floor — the 1-year gate cleanly separates signal from
noise without presupposing the exact multi-year period in advance (which is the whole point:
a chaotic/aperiodic trajectory's return spacing is not assumed, only measured).

## 4. Admission window: search for a 10-15 year sub-window with 3-15 returns

The catalogue's own `quasi_cycler` schema (`docs/notes/2026-06-16-catalogue-scope-taxonomy.md`)
requires `epoch_locked=true`, a `validity_window` of 10-15 years, and `n_returns` an integer in
[3, 15]. Given a long propagated trajectory (multi-decade, to give the search room to find a
qualifying window at all), slide a window of length in [10, 15] years across the full
propagated span and count DISTINCT returns (per §§2-3) falling inside it. The candidate is
ADMISSIBLE iff at least one such window contains between 3 and 15 returns inclusive. Report
the window bounds and the actual return count/epochs used — do NOT silently pick the first
window that happens to qualify without recording why (an honesty requirement, not a technical
one: a cherry-picked window is a fabricated result).

## 5. Bounded-geometry check within the admission window (adapted V2 discipline)

Within the admitted window, record the trajectory's Earth-relative rotating-frame position at
each return's closest-approach instant. The candidate passes the bounded-geometry check iff
these positions do not drift by more than **a factor of 3 in radial distance** from the
window's own minimum closest-approach distance across all returns in the window (i.e. the
LOOSEST individual return in the window is no more than 3x farther than the closest) — a
direct, scoped analogue of the existing V2 "long-span bounded-drift" gate
(`docs/notes/2026-06-16-catalogue-scope-taxonomy.md` §V2: "for `quasi_cycler`, 'long-span'
means the full `validity_window`, not infinite"), here applied to POSITION spread across
returns rather than a periodic orbit's own closure residual, since there is no closure
residual to measure for a non-periodic trajectory. This is a genuinely NEW numeric choice
(the factor of 3) rather than a value carried over from an existing gate; if the actual
propagated data shows a natural, tighter clustering, tighten it — record whichever choice is
made and why once real trajectories are in hand, rather than defending 3x in the abstract.

## 6. What this criterion does NOT decide

- It does not decide encounter VELOCITY / flyby-quality requirements (no `dv_band`
  classification is attempted here — #535 is a discovery screen, not a full V0-V5 gauntlet
  pass; a candidate clearing this criterion still needs the standard catalogue admission
  pipeline before being trusted as more than a candidate).
- It does not decide which SEED trajectories to search over (a separate, genuinely open
  design question — broad seeding across the co-orbital/Hilda phase space, not just at
  already-certified periodic orbits, per the #535 OUTSTANDING.md entry's own scoping).
- It does not extend to the Sun-Jupiter Hilda separatrix case's specific geometry (mu, Hill
  radius, and expected return timescale all differ from the Earth co-orbital case) — the
  ENCOUNTER/RETURN/WINDOW/GEOMETRY structure above is body-agnostic, but the 1-year minimum
  separation floor (§3) is Earth-orbital-period-motivated and should be reconsidered
  (re-derived from Jupiter's own orbital period, ~11.86 yr, or the Hilda resonance's own
  libration timescale) before reuse at Jupiter.

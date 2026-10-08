# #1050: do `#1000` families A (5:2 prograde) and D (7:3 prograde) continue to Liang, Xu & Xu 2017's polygonal-like orbits?

Task `#1050` (earthmoon-opus, 2026-10-08, lead request). Parent: `#1000`
(`docs/notes/2026-10-08-1000-earth-moon-complement-search.md` sec. 5.1b/5.1c). The eta variant is
not used. Nothing is called novel.

## 0. Pre-registration (committed before any `#1050` run)

- **Question.** Liang, Xu & Xu 2017 (digest `2026-10-06-digest-liang-xu-xu-2017-...`) print a 5:2
  polygonal-like orbit and a 7:3 one (PLPOs), in modified Kepler elements:
  - 5:2: a = 0.5620, e = 0.4660, omega = 1.3634, C = 3.0996;
  - 7:3: a = 0.5604, e = 0.3430, omega = 4.1493, C = 3.1858.
  Both have apoapsis about 0.82 / 0.75 L and no lunar pass. Are they the high-C ends of families
  A and D (C 2.41-2.59 and 2.41-2.56), which do have a lunar pass?
- **Method.**
  - The `#1000` continuation (`scripts/run_1000_continue.py`): fixed-C multiple shooting with N
    = p arcs, the secant predictor with a fall-back, C step 0.01 halved up to 3 times, and b from
    the arc STMs.
  - Start from the last converged member of the existing upward runs: A at C 2.5804, D at
    C 2.5601. D's upward run stopped on loss of convergence; it is retried once from that
    member with the same rules.
  - The 40-step cap is lifted to 400. New output files, `*_p_ext.jsonl`; the `#1000` files are
    not changed.
- **Stop rules** (per family, upward only):
  - impact at either floor;
  - C >= 3.2;
  - loss of convergence after 3 halvings (recorded as "possible fold or turning point"; a
    fixed-C continuation cannot pass a fold in C);
  - b crossing -2 (period doubling);
  - a topology change (the winding numbers about the Earth or the Moon change);
  - becoming symmetric.
  - Foreground calls under 8 minutes, checkpointed and resumed.
- **Comparison with Liang 2017.**
  - Liang's state is built from (a, e, omega) at periapsis (theta = 0), as an Earth-centred
    ellipse with omega measured from the +x axis of the rotating frame at t = 0. Liang's
    orbits are frozen tori under averaging, not necessarily exact CR3BP periodic orbits, so
    the state is propagated and its perigee set recorded.
  - Match rule, following sec. 5.1b but adapted because the PLPOs have no lunar pass:
    - the family has a member at |dC| < 0.01 from Liang's C;
    - that member's osculating Earth-centred (a, e) at one of its perigees is within 2% of
      Liang's;
    - and omega at that perigee is within 0.05 rad of Liang's omega, or of omega plus a multiple
      of the polygon's apse spacing, 2 pi / p.
  - The period is not compared, because Liang prints none.
- **Expected outcomes, stated now.**
  - (i) **Connects:** a family reaches Liang's C and matches. Then A or D is the lunar-pass
    continuation of a published Liang PLPO family, and the class-level prior art is
    strengthened (still `#972`'s call).
  - (ii) **Does not connect:** the family ends by a stop rule below Liang's C, or reaches it
    without matching. The end point and the reason are recorded.
  - My prior: (ii) is more likely, at least for D.
    - The family members pass the Moon at 11,000-16,000 km altitude with apogee about 1.0 L.
    - A PLPO's apoapsis is about 0.8 L.
    - Losing the lunar pass smoothly would need the apogee to shrink by about 0.2 L over a C
      change of 0.5-0.6.

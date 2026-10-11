# #1022: the Pluto small-moon registry entries against plu060 and the JPL SSD pages (2026-10-11)

Verdict: the registry needs no data change. Every stored value equals the JPL SSD page value; the
kernel-fitted periods agree with the pages' printed periods to the digits printed. One property of the
registry is worth knowing: its derived mean motions for the four small moons are 1.5-3.7 percent slower than
the kernel-fitted ones, because the page's printed a is not the Kepler radius of its own period.

Sources: JPL SSD satellite pages `ssd.jpl.nasa.gov/sats/phys_par/` and `/sats/elem/` (ephemeris PLU060; page
reference Brozovic and Jacobson 2024, AJ 167:256), accessed 2026-10-11; `plu060.bsp` fitted periods from the
#998 note sec. 1 (2030-01-01 plus 3 yr, in-plane angle fit about barycentre 9).

## 1. Per value

| Moon | Value | Registry | JPL page | Match |
|---|---|---|---|---|
| Charon | GM | 106.1 | 106.1 +- 0.3 | yes |
| Charon | mean radius | 606.0 | 606.0 +- 0.5 | yes |
| Charon | a | 19600 | 19,600 | yes |
| Styx | GM | 0.0003 | < 0.0003 (limit) | yes (limit stored) |
| Styx | radius, a | 5.2, 43200 | 5.2, 43,200 | yes |
| Nix | GM, radius, a | 0.0015, 18.0, 49300 | 0.0015 +- 0.0005, 18.0 +- 1.0, 49,300 | yes |
| Kerberos | GM | 0.0002 | < 0.0002 (limit) | yes (limit stored) |
| Kerberos | radius, a | 6.0, 58300 | 6.0, 58,300 | yes |
| Hydra | GM, radius, a | 0.0020, 18.5, 65200 | 0.0020 +- 0.0003, 18.5 +- 1.0, 65,200 | yes |

The registry stores no period or eccentricity; `mean_motion_deg_day` is derived from the stored a and the
system GM (975.5 km^3/s^2).

## 2. Kernel-fitted periods against the pages' printed periods

| Moon | Fitted (d) | Printed (d) | Difference (d) |
|---|---|---|---|
| Charon | 6.38722 | 6.387222 | 0.000002 |
| Styx | 20.16195 | 20.16 | 0.002 |
| Nix | 24.85472 | 24.85 | 0.005 |
| Kerberos | 32.16798 | 32.17 | -0.002 |
| Hydra | 38.20202 | 38.20 | 0.002 |

All agree to the digits the page prints (two decimals; Charon six).

## 3. The a-versus-period offset (diagnosis, partial)

Kepler a from the fitted period about the system GM, against the printed a: Charon 19596 against 19600
(+0.02 percent, consistent); Styx +2.45 percent, Nix +1.69, Kerberos +1.26, Hydra +0.98. Equivalently the
registry-derived mean motion is below the fitted one by 3.7, 2.5, 1.9 and 1.5 percent.

What the page says: "Mean elements ... are simply the elements of a precessing ellipse which has been fit in a
least squares sense to the numerically integrated orbit", with the warning that they are "not intended for
ephemeris computation". It does not state the origin (Pluto or the system barycentre) and it does not tie
the printed a to the printed period. So barycentric against Pluto-centric cannot be settled from the page's own
notes. The origin does not explain it either: Charon, in the same table, is consistent.

What the numbers say: the offset times a^2 is nearly constant over the four moons (6.9e7, 6.2e7, 6.4e7,
6.3e7 km^2 for Styx, Nix, Kerberos, Hydra), the form of a quadrupole term from the Pluto-Charon binary
(eta = 0.097, binary separation 19596 km). The plain quadrupole shift of the mean motion,
(3/8) eta (a_b/a)^2, is 0.75, 0.57, 0.41 and 0.33 percent, about a fifth of the observed offset; if the
page's a were defined from the radial (epicyclic) frequency the shift would double, to about two fifths.
The remainder is unexplained here. This is a diagnosis from the numbers, not from a source that states the
definition; the page's a should be read as a fit parameter of a precessing ellipse, not as the Kepler radius
of the mean motion.

## 4. Consequence and what is not changed

- #998 uses the fitted periods for its circles; that choice stands.
- Registry-based sweeps that take phase evolution from `mean_motion_deg_day` of the small moons (for
  example the #320 Pluto sweep, which used Nix and Hydra) advance them 1.5-2.5 percent slower than the real
  moons. That affects phasing in those sweeps, not the registry's correctness against its source. Whether to
  store the fitted mean motion (or the period) for the small moons, so that sweeps use it, is a design
  decision; not made here.

## 5. Addendum 2026-10-11: sourced explanation of the a-versus-period offset (replaces "unexplained" in sec. 3)

The Brozovic and Jacobson (2024) digest (`docs/notes/2026-10-11-digest-brozovic-jacobson-2024-post-new-horizons-orbits-masses-pluto-satellites-aj.md`, section "Registry period question") settles it. The JPL SSD page
pairs two different element sets:
- the printed a (43200, 49300, 58300, 65200 km) are the paper's Table 10, MEAN OSCULATING a, rounded to 100 km
  (43163, 49340, 58272, 65175 km);
- the printed periods (20.16, 24.85, 32.17, 38.20 d) are Table 11, the PRECESSING-ELLIPSE periods, whose a
  (42410, 48689, 57749, 64719 km) is not the Kepler radius of its mean motion (Kepler gives a period 0.87, 0.67,
  0.45 and 0.36 percent too long for those a).
Each table is consistent with itself: Table 10's Kepler period from its a, with the kernel GM_P + GM_C, matches its
own period to 0.05 percent. So the registry (a from Table 10) reproduces the OSCULATING-mean periods (20.91, 25.49,
32.78, 38.76 d against Table 10's 20.89, 25.51, 32.75, 38.74 d), and the 1.5-3.7 percent offset in sec. 3 is the
difference between osculating-mean and precessing-ellipse (true mean-motion) periods, forced by the binary
(Table 11's free eccentricities are below 1e-5 for Styx, Nix and Kerberos; the osculating 0.009-0.025 is forced).
My quadrupole estimate (a fifth to two fifths of the offset) was a partial account of the same effect. The fitted
periods used by #998 are the Table 11 periods to the printed digits.

Consequence for the registry: it is not wrong against its source; it stores osculating-mean a. A period-based
consumer must take the Table 11 period (or n), a sourced value, not a fit of ours. The 1022b/1022c notes carry
the same correction.

A second, separate observation from the same digest: the registry stores Styx 0.0003 and Kerberos 0.0002 km^3/s^2,
the rounded 1-sigma upper bounds (Table 8), whereas PLU060 was integrated with held values 0.04e-3 and 0.06e-3
(kernel 4.045e-5 and 6.038e-5). The registry values are 7.5 and 3.3 times the held ones. It does not matter for
#998 (massless targets); it overstates the bending of any lane that uses these two as flyby bodies. Proposal only:
a comment at the registry entries, or the held values with the bound as metadata. No registry change was made.

## 6. Owner decision (2026-10-11)

Question as put (OUTSTANDING #1022, waiting on owner): "(a) store fitted periods plus a derived Kepler radius for
the small moons, or (b) leave and document" (impact and options in
`docs/notes/2026-10-11-1022b-registry-period-impact.md`).

Owner decision, verbatim label: "Store the fitted period beside a". That is option (a). The sourced period is
Brozovic and Jacobson (2024) Table 11 (sec. 5 above), not a fit of ours. Pre-registration of the change:
`docs/notes/2026-10-11-1022d-fitted-periods-prereg.md`.

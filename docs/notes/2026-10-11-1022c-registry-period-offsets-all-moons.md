# #1022c: registry-period offsets for every registry moon (2026-10-11)

Read-only. No registry change. Question left open in `docs/notes/2026-10-11-1022b-registry-period-impact.md`:
does the Pluto small-moon offset (Kepler period from the registry a and the system GM longer than the real
period) occur elsewhere in `src/cyclerfinder/core/satellites.py`?

Answer: no. Of the 28 registry moons that can be compared, exactly 4 exceed 0.5 percent: Styx +3.69,
Nix +2.54, Kerberos +1.89, Hydra +1.47 percent. The largest other offset is Amalthea at +0.20 percent. The 9 small-body
moons (Sylvia, Elektra, Eugenia, Kleopatra systems) are not compared against a page: their system GM was
itself fitted to those moons' Kepler orbits, so the period is consistent by construction (the registry comment
quotes Elektra Delta's published period, 1.642112 d, against the registry Kepler value 1.64232 d, +0.01 percent).

## Method and sources

- Kepler period = 2 pi sqrt(a^3 / mu_primary), with the registry a and `PRIMARIES[primary]`.
- Reference period, in this order: the kernel-fitted sidereal period where a kernel is on disk
  (`~/dev/references/kernels`: `jup365.bsp`, `sat441.bsp`, `plu060.bsp`; fit of the in-plane angle about the
  mean orbit normal over 2 yr from 2000-01-01 12:00 for Jupiter and Saturn, 3 yr from 2030-01-01 for Pluto, relative
  to the planet centre, 4000 samples); otherwise the period printed on the JPL SSD mean-elements page
  (`ssd.jpl.nasa.gov/sats/elem/`, accessed 2026-10-11, retrieved through a page-summarising fetch tool).
- Registry source comments (`satellites.py`): a, GM and radius from the JPL SSD satellite physical and orbital
  tables, accessed 2026-06-07 (original eleven moons) and 2026-06-14 (Mars, Uranus, Neptune, Pluto moons and the
  Jupiter/Saturn additions); Styx and Kerberos from the same pages, accessed 2026-10-11 (#1037).
- Caveat on the page values: the fetched Jupiter rows (for example Io 1.762732 d, Europa 3.525463 d) disagree with
  the jup365 fit (1.769138, 3.551181) and with the known sidereal periods, so they were treated as a fetch error and
  not used; the kernel fits are used for Jupiter and Saturn. The Uranian rows match the known periods to six
  digits. The Mars, Earth (DE405) and Neptune rows come from the same fetch; they agree with known values to the
  digits printed, but they are the least checked entries here.
- The Earth entry uses the Earth+Moon system GM (403503 km^3/s^2) as the registry does.

## Table

| Moon | Primary | Registry a (km) | Kepler period (d) | Reference period (d) | Reference | Offset |
|---|---|---|---|---|---|---|
| Moon | Earth | 384400 | 27.28461 | 27.322 | page, DE405 | -0.14 % |
| Phobos | Mars | 9375 | 0.31898 | 0.3187 | page, MAR099 | +0.09 % |
| Deimos | Mars | 23457 | 1.26244 | 1.2625 | page, MAR099 | -0.01 % |
| Io | Jupiter | 421800 | 1.76995 | 1.76914 | jup365 fit | +0.05 % |
| Europa | Jupiter | 671100 | 3.55207 | 3.55118 | jup365 fit | +0.03 % |
| Ganymede | Jupiter | 1.0704e+06 | 7.15517 | 7.15456 | jup365 fit | +0.01 % |
| Callisto | Jupiter | 1.8827e+06 | 16.69062 | 16.6891 | jup365 fit | +0.01 % |
| Amalthea | Jupiter | 181400 | 0.49918 | 0.498179 | jup365 fit | +0.20 % |
| Mimas | Saturn | 185540 | 0.94368 | 0.942443 | sat441 fit | +0.13 % |
| Enceladus | Saturn | 238040 | 1.37133 | 1.37022 | sat441 fit | +0.08 % |
| Tethys | Saturn | 294670 | 1.88874 | 1.8878 | sat441 fit | +0.05 % |
| Dione | Saturn | 377420 | 2.73782 | 2.73692 | sat441 fit | +0.03 % |
| Rhea | Saturn | 527070 | 4.51825 | 4.5175 | sat441 fit | +0.02 % |
| Titan | Saturn | 1.22187e+06 | 15.94797 | 15.9454 | sat441 fit | +0.02 % |
| Hyperion | Saturn | 1.4815e+06 | 21.29219 | 21.2632 | sat441 fit | +0.14 % |
| Iapetus | Saturn | 3.5617e+06 | 79.36953 | 79.3256 | sat441 fit | +0.06 % |
| Miranda | Uranus | 129846 | 1.41351 | 1.41348 | page, URA182 | +0.00 % |
| Ariel | Uranus | 190929 | 2.52037 | 2.52038 | page, URA182 | -0.00 % |
| Umbriel | Uranus | 265986 | 4.14423 | 4.14418 | page, URA182 | +0.00 % |
| Titania | Uranus | 436298 | 8.70624 | 8.70587 | page, URA182 | +0.00 % |
| Oberon | Uranus | 583511 | 13.46572 | 13.4632 | page, URA182 | +0.02 % |
| Triton | Neptune | 354800 | 5.87792 | 5.87699 | page, NEP097 | +0.02 % |
| Proteus | Neptune | 117600 | 1.12166 | 1.12231 | page, NEP097 | -0.06 % |
| Charon | Pluto | 19600 | 6.38905 | 6.38722 | plu060 fit | +0.03 % |
| Styx | Pluto | 43200 | 20.90634 | 20.162 | plu060 fit | +3.69 % **>0.5** |
| Nix | Pluto | 49300 | 25.48724 | 24.8547 | plu060 fit | +2.54 % **>0.5** |
| Kerberos | Pluto | 58300 | 32.77595 | 32.168 | plu060 fit | +1.89 % **>0.5** |
| Hydra | Pluto | 65200 | 38.76358 | 38.202 | plu060 fit | +1.47 % **>0.5** |

## Offsets above 0.5 percent

| Moon | Registry source | Offset | Likely cause |
|---|---|---|---|
| Styx | JPL SSD pages, accessed 2026-10-11 (#1037) | +3.69 % | circumbinary orbit about Pluto and Charon: the printed a is the fit parameter of a precessing ellipse, not the Kepler radius of its period |
| Nix | JPL SSD pages, accessed 2026-06-14 | +2.54 % | same |
| Kerberos | JPL SSD pages, accessed 2026-10-11 (#1037) | +1.89 % | same |
| Hydra | JPL SSD pages, accessed 2026-06-14 | +1.47 % | same |

The cause is the binary, not the choice of origin: Charon in the same table agrees to 0.03 percent, and the
offsets scale as a^-2 (the form of the Pluto-Charon quadrupole). The plain quadrupole mean-motion shift explains
about a fifth of each; the rest stays unexplained (#1022 note sec. 3). A barycentric-versus-Pluto-centric origin
cannot be read from the page.

## Offsets below the threshold, with a cause where one is clear

Oblateness explains the next group. A planet's J2 raises the mean motion by (3/2) J2 (R/a)^2, so the Kepler period
from the registry a is longer by 0.75 J2 (R/a)^2 (J2 and R from the standard values: Jupiter 0.014696, 71492 km;
Saturn 0.016298, 60268 km): Amalthea predicted 0.17 percent against 0.20 observed, Mimas 0.13 against 0.13,
Enceladus 0.08 against 0.08. Hyperion (+0.14 percent) has e about 0.1 and is solar-perturbed, so a Kepler comparison
is loose. The Moon (-0.14 percent, the other sign) is the Sun's perturbation of the lunar orbit; the page's 27.322 d
is the mean sidereal period. None of these needs action.

## Count

4 of the 28 compared moons are affected (all Pluto small moons); 0 of the 24 others exceed 0.5 percent;
the 9 small-body moons are consistent by construction. Where there is no page and no kernel on disk (Mars,
Uranus, Neptune, Earth), the check rests on the page values noted above.

## Addendum 2026-10-11: the cause of the Pluto offset is sourced

Brozovic and Jacobson (2024), Tables 10 and 11 (digest `docs/notes/2026-10-11-digest-brozovic-jacobson-2024-post-new-horizons-orbits-masses-pluto-satellites-aj.md`): the registry's a is the MEAN OSCULATING a, and its
Kepler period reproduces the paper's osculating period (Table 10) to 0.05-0.1 percent; the printed periods are the
precessing-ellipse (true mean-motion) periods of Table 11. The offsets in this note are the osculating-versus-mean
difference, not a registry error. Where this note says the cause is "the binary quadrupole, remainder unexplained",
read: the binary forces the osculating elements away from the mean ones, and the paper's tables give both. See
`docs/notes/2026-10-11-1022-pluto-small-moon-registry-verification.md` sec. 5. The numbers and the options are
unchanged; for option (a) the sourced period is Table 11's.

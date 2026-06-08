# S1L1 / 4.991gG2 source dig — Russell App-C real-eph states + topology reconciliation (#166)

**Date:** 2026-06-08
**Task:** #166 — comb every held description of the S1L1 / `russell-ch4-4.991gG2`
cycler for the geometry/phasing insight that explains or resolves the #165 DRIFT
(our two-arc construction crosses Mars's radius ~110° from where real DE440 Mars
actually is).
**Deliverable:** this insights note. No production code; read-only on all src/data.
**Sources combed:** Russell 2004 dissertation (Appendix C + §4.8 + §5.4 + Table 4.9);
McConaghy/Longuski/Byrnes 2002 (AIAA 2002-4420) nomenclature; Patel 2019 (Florida
Institute of Technology MSc thesis, "Earth-Mars Cycler Vehicle Conceptual Design").
All sources cited by author/title/venue only.

---

## BOTTOM LINE (the single most useful finding)

**Russell's Appendix C gives the real-ephemeris Mars-encounter epoch AND longitude
exactly, and it reconstructs on DE440 to a perfect intercept.** I transcribed the
`4.991gG2` (#83) "DATA NECESSARY TO REPRODUCE" block (DE405, J2000-referenced),
rebuilt the Mars-transit leg from its published `(epoch, time-start, v∞-vector)`,
and propagated it two-body-Sun over DE440 planet positions:

> Earth departure **2026-12-15**, Mars flyby **2027-06-13** (transit 179.8 d,
> matching the published transit table), spacecraft–Mars miss **1.1e-8 AU (≈1.7 km)**,
> reconstructed v∞ at Mars **5.2480 km/s** vs the published table value **5.248**
> (4-decimal match), spacecraft and Mars at the **same heliocentric ecliptic
> longitude 201.0°** (difference 4e-7 deg).

So the longitude-rendezvous constraint our #164/#165 construction omitted is **not
missing from the sources** — it is fully specified by Russell's per-leg state block.
The #165 ~110° miss is a property of OUR construction, not an unavoidable data gap.

**And the topology IS wrong.** Russell's own data shows the lowercase `g` arc is a
pure **Earth-to-Earth free return that never reaches Mars** (aphelion 1.27 AU,
closest approach to Mars 1.05 AU, returns to Earth to 0.0000 AU). Only the uppercase
`G` arc is the Mars-transit leg. Our two-arc construction modeled **both** arcs as
crossing Mars's radius — that is the structural error.

**Recommendation: (b) the corrected sequence topology** (with (a) Russell's exact
real-eph seed now in hand). Details in the final section.

---

## 1. Russell App-C `4.991gG2` (#83) — the real-eph node geometry (PRIORITY, transcribed)

Source: Russell 2004 dissertation, Appendix C ("Cycler Trajectories Using an
Ephemeris Model"). Model: **DE405**, positions/velocities referenced to Earth mean
equator/equinox of J2000 then rotated to the ecliptic by the mean obliquity
**0.409092629205 rad**. Epoch convention: "EPOCH TIME (days after J2000)". The
prose states "The 'Epoch time' and the detailed leg information are the only values
required to reproduce the solution."

This is the block the #142 transcription pass recorded as PRESENT-but-not-transcribed.
It is now transcribed verbatim.

### Header (verbatim)

```
================ PARENT CYCLER 4.991gG2       =======================
Parent cycler number                                                83
Approximate search space (synodic periods after J2000)              11
Number of steps to walk eccentricity/inclination              1 /    1
Number of cycles                                                     7
Total delta v over 30.07 years (km/s)                        0.000000
```

Russell §5.4 / Table 5.5 cross-reference: `4.991gG2`(#83), launch **Jun-25**, total
∆v **0** m/s, E-M transit avg **165 d**, avg v∞E− 5.37 / v∞E+ 5.37 / v∞M− 5.48 /
v∞M+ 5.48 km/s. §5.4 prose: *"Most notable is cycler 4.991gG2(#83)… Also known as
the 'S1L1' cycler, this two-synodic period cycler is essentially ballistic for all
launch dates."*

### EARTH-TO-MARS TRANSIT LEG CHARACTERISTICS (verbatim)

| LEG | E-M transit time (d) | Earth vinfm | Earth vinfp | Mars vinfm | Mars vinfp |
|---|---|---|---|---|---|
| 2  | 179.8 | 5.818 | 5.818 | 5.248 | 5.248 |
| 5  | 125.5 | 5.764 | 5.764 | 7.693 | 7.693 |
| 8  | 138.4 | 3.752 | 3.752 | 4.657 | 4.657 |
| 11 | 210.9 | 5.532 | 5.532 | 3.198 | 3.198 |
| 14 | 154.6 | 6.946 | 6.946 | 6.263 | 6.263 |
| 17 | 112.2 | 5.121 | 5.121 | 8.046 | 8.046 |
| 20 | 230.2 | 4.646 | 4.646 | 3.219 | 3.219 |
| AVERAGE | 164.5 | 5.368 | 5.368 | 5.475 | 5.475 |

Note `vinfm == vinfp` at every Mars leg → **ballistic Mars flyby** (no ∆v; the v∞
magnitude is unchanged across the flyby, only its direction bends). The seven transit
legs are the seven Mars encounters across the 7 propagated cycles; **one Mars
encounter per cycle**, with the per-cycle Mars v∞ ranging 3.2–8.0 km/s as the real
geometry drifts cycle-to-cycle (only ~repeating every ~15 yr per §5.4).

### DATA NECESSARY TO REPRODUCE (verbatim)

```
EPOCH TIME (days after J2000)                             9325.742435
LEG E/M time start        vinfx            vinfy            vinfz           time dv          dvx              dvy              dxz
         (days)           (km/s)           (km/s)           (km/s)          (days)           (km/s)           (km/s)           (km/s)
 1 E -0.174074986E+02 -0.238703220E+00 0.579754982E+01 -0.176632468E-02   0.240347612E+03 -0.246767343E-11 0.902746863E-12 0.301369564E-10
 2 E 0.519582315E+03 -0.227826290E+01 0.532198991E+01 0.573813715E+00     0.564535086E+03 -0.727488314E-13 0.229820172E-12 -0.412405272E-11
 3 M 0.699393400E+03 -0.499623059E+01 0.707421771E+00 0.144166283E+01     0.111468222E+04 -0.968220811E-11 0.177011134E-10 0.106068623E-10
 4 E 0.152997105E+04 -0.575383361E+01 0.386418107E+00 0.628610087E-04     0.178763540E+04 0.276478627E-10 -0.757910553E-11 0.224923555E-11
 5 E 0.206677179E+04 -0.550681663E+01 -0.160954336E+01 -0.550229544E+00   0.208560007E+04 0.410319943E-10 -0.107006918E-10 0.369882340E-11
 6 M 0.219229367E+04 -0.355470987E+01 -0.672927795E+01 0.112658626E+01    0.244848716E+04 0.240567159E-11 -0.249330086E-11 0.698061980E-11
 7 E 0.310727040E+04 0.130066966E+01 -0.353985845E+01 0.127838602E-02     0.336108937E+04 -0.279214975E-11 -0.115736777E-12 -0.298533208E-10
 8 E 0.363605992E+04 0.372962888E+01 -0.150726402E-01 -0.411042515E+00    0.365682578E+04 -0.370357687E-12 -0.279917663E-11 0.124055358E-11
 9 M 0.377449899E+04 0.463453129E+01 -0.454933848E+00 -0.467039188E-01    0.426411381E+04 -0.104973257E-10 0.235441673E-11 -0.311497269E-11
10 E 0.468119309E+04 0.398651382E+01 0.380698910E+01 0.599174579E-03      0.493839453E+04 0.215667218E-10 0.108875240E-10 0.209391049E-12
11 E 0.521702941E+04 0.195647577E+01 0.500455176E+01 0.131545295E+01      0.524866457E+04 0.177176472E-10 0.241691459E-10 -0.128901836E-11
12 M 0.542793047E+04 -0.106091993E+01 0.286375761E+01 0.950243209E+00     0.557829938E+04 0.652424747E-11 0.728315006E-11 -0.625360942E-11
13 E 0.621934578E+04 -0.346922206E+01 0.601287549E+01 -0.294263108E-02    0.647390151E+04 0.694255325E-11 -0.189742179E-10 -0.347896187E-11
14 E 0.676095373E+04 -0.494762984E+01 0.487023018E+01 0.229687664E+00     0.678415083E+04 0.306801663E-10 -0.483515188E-10 -0.203554640E-11
15 M 0.691560104E+04 -0.580355995E+01 -0.130028398E+01 0.196236847E+01    0.750624237E+04 -0.180880049E-11 -0.783703320E-12 -0.101848364E-11
16 E 0.777160296E+04 -0.455750356E+01 -0.237566423E+01 -0.883416315E-03   0.802804869E+04 0.972188929E-11 0.213121008E-11 -0.530653406E-11
17 E 0.830586489E+04 -0.289917665E+01 -0.412851216E+01 -0.880976268E+00   0.832269007E+04 0.425084649E-11 -0.139214809E-11 0.707337715E-12
18 M 0.841803274E+04 -0.474724840E+00 -0.793519915E+01 0.124113571E+01    0.890636978E+04 0.642835129E-11 0.118613663E-10 -0.129874748E-10
19 E 0.933942338E+04 0.322773002E+01 -0.335738552E+01 0.125353647E-02     0.960027089E+04 0.866703124E-11 -0.120333181E-10 -0.195906858E-10
20 E 0.987176523E+04 0.254872658E+01 0.229406465E+01 0.313437961E+01      0.990629794E+04 0.155550229E-10 -0.473859434E-11 -0.870795245E-11
21 M 0.101019833E+05 0.212149789E+01 0.424849528E+00 -0.238347620E+01     0.105258514E+05 -0.207499508E-11 -0.348202361E-11 0.729177865E-11
22 E 0.109670202E+05
```

All `dv` components are ~1e-11 to 1e-10 km/s — numerically zero, confirming the
ballistic claim. Russell's reproduction recipe (§App-C, verbatim): at each leg start,
"planet velocity at this time from DE405; spacecraft velocity = planet velocity +
given v∞", then propagate by Kepler's equation.

### Per-cycle structure read directly off the E/M column

The encounter sequence (E/M) over 7 cycles is:
`E E M | E E M | E E M | E E M | E E M | E E M | E E M | E`
i.e. **three encounters per cycle: E → E → M**. The Mars flybys are encounters
3, 6, 9, 12, 15, 18, 21 — exactly **one Mars encounter per cycle**. The transit legs
in the characteristics table (legs 2, 5, 8, 11, 14, 17, 20) are the **E#k → M#(k+1)**
outbound Earth-to-Mars legs (e.g. leg 2: E#2 at 519.58 d → M#3 at 699.39 d = 179.8 d
transit, matching the table).

This is the literal embodiment of the [[s1l1-nomenclature]] memory: **S/L are
Earth-to-Earth resonant intervals; Mars is a flyby of opportunity on ONE leg.**

---

## 2. Independent DE440 reconstruction — Russell's geometry CLOSES (the decisive check)

I rebuilt the first Mars-transit leg (leg 2: E#2 → M#3) from the App-C block on
**DE440** (astropy, ecliptic frame via Russell's obliquity), two-body Sun propagation:

| quantity | value | source / check |
|---|---|---|
| Earth-depart epoch | 2026-12-15 (J2000+9845.32 d) | epoch 9325.74 + time-start 519.58 |
| Mars-flyby epoch | 2027-06-13 (J2000+10025.14 d) | epoch + 699.39 |
| transit time | 179.81 d | matches published 179.8 d |
| Earth helio dist at depart | 0.984 AU | sanity |
| Mars helio dist at arrive | 1.618 AU | near cycler aphelion 1.64 AU |
| **spacecraft–Mars miss** | **1.1e-8 AU (~1.7 km)** | exact intercept |
| **reconstructed v∞ at Mars** | **5.2480 km/s** | published table 5.248 |
| **Mars heliocentric ecliptic longitude at flyby** | **201.0°** | the rendezvous constraint |
| SC longitude − Mars longitude | 4e-7 deg | exact rendezvous |

The reconstruction is a 4-decimal match on v∞ and a ~1.7 km position match — Russell's
published per-leg state, on a DIFFERENT ephemeris (DE440 vs his DE405), lands the
spacecraft on real Mars. **The published cycler is a real, longitude-correct,
ballistic Mars intercept** — there is no model we "can't reconstruct from held data";
the held data reconstructs cleanly.

This directly answers the #166 question: **Russell's real-eph Mars encounter happens
at heliocentric ecliptic longitude 201.0° on 2027-06-13, and our #164/#165
construction missed it because it never targeted Mars's true position (and modeled the
wrong topology).** The "~110° off" of #165 is the signature of a radius+speed+ToF
closure with no longitude target, exactly as #165 diagnosed.

---

## 3. Topology — the `g` arc does NOT go to Mars (our two-arc construction is wrong)

Russell §4.8 (Table 4.9 setup), verbatim:

> "All legs reported are Earth-Earth free-returns. The transit times and Mars v∞ are
> calculated using the designated transit leg, **as indicated by an uppercase
> descriptor letter**… for these cyclers **only one Mars encounter is guaranteed
> during each cycler period**. Therefore, a given cycler of this class must be
> designated as an **inbound or outbound** cycler."

> "An **uppercase** letter indicates that the transit times and Mars v∞ were
> calculated using this leg… the first number in the parenthesis is the time of flight
> in years for that Earth-Earth leg."

Table 4.9 row for `4.991gG2` (circular-coplanar model, verbatim):

```
v∞E   v∞M  tout tin aphel.  TR    Leg 1                   Leg 2
4.99  5.10 150  150  1.64   2.65  g(1.4612,526.02,Ll)     G(2.8096,651.46,U)
```

- **Leg 1 `g(1.4612,…)`** — lowercase `g` = a generic Earth-to-Earth free return,
  ToF 1.4612 yr, that **does not encounter Mars**.
- **Leg 2 `G(2.8096,…)`** — UPPERCASE `G` = the designated **Mars-transit** leg,
  ToF 2.8096 yr; this is the leg the v∞M = 5.10 and aphelion 1.64 AU are computed from.

I verified this against the App-C real-eph data by rebuilding **Leg 1 (the g-arc,
E#1→E#2)** on DE440:

| g-arc (Leg 1) property | value | interpretation |
|---|---|---|
| returns to Earth? | miss **0.0000 AU** | yes — it IS Earth-to-Earth |
| arc aphelion | **1.27 AU** | well below Mars's 1.52 AU (and the 1.64 AU cycler aphelion) |
| closest approach to Mars over the arc | **1.05 AU** | never anywhere near Mars |

**This is the structural error in #164/#165.** Our `continuation_chain` two-arc
construction (note `2026-06-08-continuation-chain-s1l1-results.md`) gave BOTH arcs a
"Mars crossing epoch" (#165 per-encounter table, arc1 AND arc2 each crossing Mars's
radius). Russell's cycler has the Mars encounter on **only the G-arc**; the g-arc is a
sub-Mars-aphelion Earth-to-Earth resonant return. Forcing the g-arc to cross Mars's
radius is geometrically impossible for the real cycler (its aphelion is 1.27 AU) and
is precisely why the closure was a `(radius, V∞, ToF)` artifact with no real-Mars
longitude solution.

Reconciliation with the descriptor ToF debate (multi-arc-classification.md §7/§12,
blocker schema-v4.1): the descriptor's two arcs `g(1.4612 yr) + G(2.8096 yr)` sum to
4.27 yr = 2 synodic periods — that part is right. What was wrong was treating them as
**two symmetric Mars-crossing arcs**. They are **one Earth-to-Earth resonant return
(g) + one Mars-transit-and-return leg (G)**.

---

## 4. Anchor-framing reconciliation — which V∞ does the real-eph cycler land at?

| framing | source | E v∞ | M v∞ | model |
|---|---|---|---|---|
| spec §9 / CPOM | Rogers 2012 Table 1 (a=1.30/e=0.257) | 5.65 | 3.05 | circular-coplanar single-ellipse |
| Russell coplanar | Russell 2004 Table 4.9 | 4.99 | 5.10 | circular-coplanar two-arc descriptor |
| Russell **real-eph** | Russell 2004 App-C #83 + Table 5.5 | **5.37** (avg E±) | **5.48** (avg M±) | **DE405 ephemeris** |
| McConaghy real-eph | (cited in blocker) | 4.7 | 5.0 | real-eph |
| Patel / McConaghy | Patel 2019 Table 1 | launch 2.492; **flyby 3.657** | — | real-eph (Earth-flyby v∞) |
| Sanchez Net 2022 | near-ballistic survey | Mars 5.2–7.3 | Earth 3.6–5.7 | real-date patched-conic |

Key reconciliation point: the **per-leg Mars v∞ in App-C ranges 3.2–8.0 km/s across
the 7 cycles** (table in §1), averaging 5.48. The single published "5.10" (Table 4.9)
is the **circular-coplanar idealization**; the single "3.05" (Rogers/CPOM) is a
**different idealization of the same cycler**; on the real ephemeris the Mars v∞ is
**not a single number** — it breathes 3.2–8.0 cycle-to-cycle. So the question "which
framing does a longitude-correct real-eph closure land at?" has a concrete answer:
**Russell's own real-eph values 5.37 (Earth) / 5.48 (Mars) averaged, with per-cycle
Mars v∞ 3.2–8.0** — NOT a clean reproduction of any single coplanar anchor. The
3.05 and 5.10 are coplanar artifacts; the real cycler does not hold either.

Patel's "Earth flyby v∞ 3.657 km/s" is the **Earth-flyby** v∞ (the resonant-return
turning), a different quantity from the Earth-departure/Mars v∞ — consistent, not
contradictory.

This is why every prior attempt to "close S1L1 at the anchors" landed off-family: the
anchors are model-specific idealizations, and the real cycler's v∞ is epoch-dependent.
The #164 closure at 4.99/5.10 was at the COPLANAR anchor — but the coplanar anchor is
not what the real ephemeris produces.

---

## 5. RECOMMENDATION

**Recommendation (b): build the corrected sequence topology — and we now also have
(a) Russell's exact real-eph seed to do it.**

The sources fully disambiguate. Concretely:

### (b) Corrected topology — the per-cycle sequence to build
Per cycle (2 synodic periods, ~4.27 yr), the sequence is **E → [g: Earth-to-Earth
generic free return, 1.4612 yr, aphelion ~1.27 AU, NO Mars] → E (flyby) → [G:
Mars-transit leg, 2.8096 yr, aphelion 1.64 AU, ONE Mars flyby] → E**. The Mars
encounter is on the **G leg only**. Any construction must:
1. patch a g free-return (Earth→Earth, sub-Mars aphelion) and a G Mars-transit leg at
   the intermediate Earth flyby;
2. enforce a **longitude rendezvous with real DE440 Mars on the G leg only** (a
   Lambert/shooting node targeting Mars's TRUE position, not just its radius);
3. leave the g leg free of any Mars-crossing constraint.

The #164/#165 construction violated (1) and (3) (both arcs Mars-crossing) and lacked
(2) (no longitude target) — the three errors that produced the 2.6 AU drift.

### (a) The exact real-eph seed (from Russell App-C #83, now transcribed)
- **Epoch:** J2000 + 9325.742435 d (2025-07-14 TDB).
- **Launch (E#1):** time-start −17.41 d (2025-06-26), v∞ = (−0.2387, 5.7975,
  −0.001766) km/s ecliptic.
- **G-leg Mars-transit (leg 2):** depart Earth E#2 at +519.58 d (2026-12-15), v∞ =
  (−2.2783, 5.3220, 0.5738) km/s; arrive Mars M#3 at +699.39 d (2027-06-13), Mars at
  ecliptic longitude **201.0°**, v∞M = 5.248 km/s, transit 179.8 d.
- This seed reconstructs on DE440 to a 1.7 km Mars intercept (§2) — it is a
  ready-made V3 seed/golden for a longitude-correct closure or cross-validation.

The harness already has the machinery (`nbody/propagator.py` RestrictedNBody over
DE440; `core/ephemeris.py`). A clean V3 path: seed the G-leg from the App-C state
above, propagate on DE440, and gate the spacecraft–Mars miss + v∞ against Russell's
published per-leg numbers (a legitimate golden — the EXPECTED side traces entirely to
Russell's printed App-C block, not to anything our code computed).

### What is NOT needed
No acquisition (#116) is required for S1L1 itself: Russell's App-C supplies the full
real-eph reproduction data, and it reconstructs cleanly. The earlier worry that "the
published values trace to a model we can't reconstruct from held data" is **refuted**
for this cycler — the held data (Russell App-C #83) reconstructs to 1.7 km / 4-decimal
v∞ on an independent ephemeris.

### Caveat on "exact anchor" reproduction
A real-eph S1L1 will **not** reproduce 5.10 (Russell coplanar) or 3.05 (Rogers/CPOM)
Mars v∞ — those are coplanar idealizations. It reproduces Russell's real-eph per-leg
values (Mars v∞ breathing 3.2–8.0, avg 5.48). Any golden gate must target the
**App-C real-eph numbers**, not the coplanar anchors, or it will mis-judge a correct
closure as a miss (the inverse of the #164 trap).

---

## Provenance / honesty notes

- Russell App-C `4.991gG2` block transcribed verbatim from the held dissertation
  text layer (clean digital typeset, not OCR-of-image). DE440 reconstruction used the
  project's astropy backend; the EXPECTED comparison values (179.8 d, 5.248 km/s) are
  Russell's printed numbers — the reconstruction is a cross-check of two published
  quantities against an independent ephemeris, not a fit.
- The reconstruction code was a throwaway probe run in this session; **no production
  code was added or modified.** Only this note is committed.
- Memory updates implied (for a later pass, not done here): [[s1l1-nomenclature]] is
  CONFIRMED by Russell's own data (g = Earth-Earth, G = Mars transit, one Mars
  encounter/cycle); the [[s1l1-realeph-closure-blocker]] "can't reconstruct from held
  data" lean is now **refuted** for 4.991gG2 — App-C reconstructs to 1.7 km.

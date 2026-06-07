# Rall (Hollister adv.) 1970 — Appendices E & F transcription (M4-1 / M5-1 / M5-2 + siblings)

Transcribed 2026-06-07 (Task #142). VERBATIM data-transcription pass of the
primary-source encounter tables. Companion to the method-mining note
`2026-06-07-hollister-rall-1970-periodic-orbits-mining.md` (read that first for
algorithm/context; this note is the numeric member data the mining pass flagged as
NOT-yet-extracted, "§7 honest not-extractable list").

**Source (cite exactly, no file path):**
Rall, C. S., "Free-Fall Periodic Orbits Connecting Earth and Mars," Sc.D. thesis,
MIT Dept. of Aeronautics and Astronautics, report **TE-34** (also published under
NASA Grant NGR 22-009-010, Measurement Systems Laboratory), October 1969 (pub.
1970). Thesis supervisors: W. M. Hollister, W. E. VanderVelde, J. E. Potter,
W. McKay. 225 pp.

> Scanned typescript, clean. Vision-read of front-matter/TOC (pp.v-ix) for
> appendix locations + column legends, **Appendix E** (Circular Coplanar
> Earth-Mars Periodic Orbits, report pp.179-188) in full, and **Appendix F**
> (Eccentric Inclined Earth-Mars Periodic Orbits, report pp.189-203) in full. All
> table numerals read unambiguously unless flagged. Page numbers below are the
> REPORT's printed page numbers (the appendices are also referenced by the TOC at
> report p.vii: Appendix E p.179, Appendix F p.189, References p.205).

---

## 0. Where the orbit-member data lives

Two appendices hold the periodic-orbit member data; **no other appendix holds
encounter dates**. (Appendix A = computer program; B = mean-anomaly calc for a
half-revolution return; C = number of combinations; D = combinations of direct
returns at Earth — these are method/enumeration, not orbit members.)

- **Appendix E** "CIRCULAR COPLANAR EARTH-MARS PERIODIC ORBITS" (pp.179-188):
  the headline circular-coplanar member list. Encounters at Earth/Mars,
  collapsed (only the Earth encounters immediately adjacent to Mars are listed —
  see legend). Covers **M4-1, M5-1, M5-2, M5-3, M5-4, M5-5, M5-6, M6-1 … M6-11**
  (the "M5-3"..."M6-11" rows are siblings/near-misses, included verbatim below).
- **Appendix F** "ECCENTRIC INCLINED EARTH-MARS PERIODIC ORBITS" (pp.189-203):
  full per-encounter listings (every Earth + Mars encounter, not collapsed) for
  the **11 orbits that converged in the eccentric-inclined model**: M4-1a, the
  five versions M5-1a..e, and the five versions M5-2a..e. Solar-system Model
  III.B (eccentric inclined, constant mean elements).

---

## 1. Date convention (VERBATIM — do NOT convert)

**TWO DIFFERENT date encodings between the two appendices.** Both are
**Julian Date**, NOT calendar. Neither appendix states ET vs UT explicitly →
flag for any future conversion (1969-era patched-conic on a mean-element series;
ephemeris-time vs UT distinction is unstated and must be treated as unknown).

### Appendix E (p.180, VERBATIM):
> "The second column gives the dates of encounter rounded to the nearest day
> corresponding to those which would occur if the solar system did correspond to
> the circular coplanar Model I.B.; the dates are listed as **Julian Date minus
> 2440000**."

So an Appendix-E date `D` ⇒ JD = D + 2440000. (e.g. M4-1 Earth `546` ⇒
JD 2440546; `3666` ⇒ JD 2443666.) Reference epoch JD 2440000.0 = 1968-05-23 (do
NOT bake this in; stated only to orient — conversion is a derived value).

Also from p.180-181: the column-3 signed numbers (e.g. `-148`, `+16`, `+622`)
are NOT dates — they are "the dates (rounded to the nearest day) of the first
Earth-Mars-Earth round trip **relative to the date of opposition**." Second
round trip dates are "the **negatives** of the dates given for the first round
trip segment in each case" (the reciprocity/symmetry stated in the mining note).

### Appendix F (p.189, VERBATIM):
> "The second column in each list gives the **Julian Date** of each encounter
> rounded to the nearest day. **Only the last five digits of the Julian Date are
> listed, and the first digit listed is separated from the other four by a
> hyphen.**"

So an Appendix-F date `d1-dddd` ⇒ JD = 244 d1 dddd (the 5 printed digits are the
trailing 5 of a JD whose leading group is "244"). e.g. M4-1a `4-0538` ⇒
JD 2440538; `8-7334` ⇒ JD 2448**7334**? No — exactly five digits total: 4,0,5,3,8
⇒ JD 2440538; `8-7334` ⇒ digits 8,7,3,3,4 ⇒ JD 2448**7334** is six → read as
JD 244 8 7334 has six trailing; the printed token `8-7334` = digits 8 7 3 3 4
(five) ⇒ JD 2448 7334 → **JD 2448734**? FLAG: the "8-7334" tokens are
unambiguous as 5 glyphs but the JD reconstruction rule "first digit + four
digits" gives 5 trailing digits appended to "244" = JD 244_87334 = 8 digits,
which is internally consistent (JD ~2.4M is 7 digits, so the leading group is
"244" + 4 digits normally; a 5-digit tail means the leading group is "24" — i.e.
JD 24 4 87334?). **The exact JD reconstruction is therefore a DERIVED value with
a documented ambiguity; do NOT convert without re-reading p.189 against a known
anchor.** What IS unambiguous and transcribed verbatim below is the printed token
itself. (Cross-check that resolves it: M4-1a first Earth `4-0538` must be ≈ the
Appendix-E M4-1 first listed Earth `546` ⇒ both ≈ JD 2440538/2440546, so the
correct rule is JD = 2440000 + (token with hyphen removed), i.e. `4-0538`→40538
is wrong; `4-0538`→0538 with leading "4" giving 2440538. The consistent rule:
**JD = 244 0000 + the five printed digits read as one number is too big; the
working rule matching Appendix E is JD-2440000 = the printed digits read as a
number minus the leading-digit place.** Leave as-is — see eligibility note: use
Appendix E's clean "JD-2440000" dates as the canonical anchor, treat Appendix F
tokens as the same JD family rounded under the eccentric model.)

**Bottom line for downstream code: trust Appendix E's "JD minus 2440000"
integers as the primary date provenance. Appendix F tokens are the per-encounter
refinement under the eccentric model; their JD reconstruction carries a flagged
ambiguity and must be confirmed before any conversion.**

---

## 2. Units (VERBATIM column legends)

### Appendix E legend (header p.182; prose pp.180-181):
Columns: **(1) planet** ("E" = Earth, "M" = Mars); **(2) date** (JD − 2440000,
nearest day); **(3) direct returns & date** — direct-return scheme symbols (per
Appendix D: HR/FR/SiSR/LiSR etc.) AND the signed date-relative-to-opposition;
**(4) speed** = "the hyperbolic excess speed at each encounter in **EMOS**"
(Earth Mean Orbital Speed units; EMOS = 29.77 km/s per mining-note Table 3-1
caption). "Only three speeds need to be listed, because the speeds at both Mars
encounters are the same by symmetry and because a series of direct return
trajectories at a planet in a circular orbit will result in the same hyperbolic
excess speed"; **(5) passing distance** = "planetary passing distances in units
of the radius of the planet encountered ... all of the **different** planetary
passing distances are listed."

### Appendix F legend (prose pp.189-190):
Columns: **(1) planet** (E/M; ALL encounters listed, not collapsed); **(2) date**
(Julian Date, last-five-digits hyphen form — see §1); **(3) speed** = hyperbolic
excess speed in **EMOS**; **(4) turn angle** in **degrees**; **(5) passing
distance** = planetocentric passing distance in units of the **local planetary
radii**; **(6) trajectory type** = type of trajectory to FOLLOW each encounter
(symbols as before + "**IP**" = interplanetary trajectory; FR/HR/SiSR/LiSR =
direct-return types).

> Appendix F note (p.190): "The series of full revolution returns in each case is
> optimal in that the minimum passing distance for the series is maximized.
> However, the combinations of full revolution returns with a half revolution
> return in all of the versions of periodic orbit M5-2 are not optimally
> arranged."

---

## 3. APPENDIX E — Circular Coplanar member list (VERBATIM)

Column order as printed: **planet | date(JD−2440000) | (direct-returns)/[signed
date rel. opposition] | speed(EMOS) | passing-dist(planet radii)**. Where a row
shows a parenthesised return scheme (e.g. `(S1SR)`, `(3FR)`) that token is the
column-3 direct-return symbol; where it shows a signed integer (`-148`) that is
the date-relative-to-opposition. Trailing prose per orbit quoted verbatim.

### M4-1 (p.182)
| planet | date | col3 | speed | pass.dist |
|---|---|---|---|---|
| E | 546 | (S1SR) | | |
| E | 1030 | −148 | 0.257 | 1.54 |
| M | 1194 | +16 | 0.314 | 3.77 |
| E | 1800 | +622 | 0.181 | 1.30 |
| E | 2531 | (3FR) | | |
| M | 3502 | | | |
| E | 3666 | | | |
> "The 'a' version exists for the best approximations used."

### M5-1 (p.182)
| E | 448 | (4FR) | | |
| E | 1901 | −48 | 0.249 | 1.78 |
| M | 2008 | +49 | 0.316 | 7.3 |
| E | 2580 | +622 | 0.211 | 1.55 |
| E | 3676 | (3FR) | | |
| M | 4249 | | | |
| E | 4348 | | | |
> "All versions exist for the best approximations used."

### M5-2 (p.183)
| E | 539 | (FR)(HR)(2FR) | | 1.42 |
| E | 1818 | −141 | 0.245 | 2.06 |
| M | 1978 | +20 | 0.314 | 4.79 |
| E | 2580 | +622 | 0.183 | 1.37 |
| E | 3676 | (3FR) | | |
| M | 4278 | | | |
| E | 4439 | | | |
> "All versions exist for the best approximations used."

### M5-3 (p.183) — sibling / non-existent "a"
| E | 574 | (2FR)(S1SR) | | 2.21, ∞, 8.95 |
| E | 1783 | −175 | 0.302 | 1.54 |
| M | 1958 | −0 | 0.316 | 2.00 |
| E | 2580 | +622 | 0.175 | 0.94 |
| E | 3676 | (3FR) | | |
| M | 4298 | | | |
| E | 4473 | | | |
> "The 'a' version does not exist, at least in the region of the circular
> coplanar solution. Discussed in Chapter 2."

### M5-4 (p.183)
| E | 566 | (FR)(S2SR) | | 2.19, 18.0 |
| E | 1791 | −168 | 0.289 | 1.79 |
| M | 1963 | +5 | 0.316 | 2.33 |
| E | 2580 | +622 | 0.176 | 1.06 |
| E | 3676 | (3FR) | | |
| M | 4293 | | | |
| E | 4466 | | | |
> "Existence can neither be confirmed nor denied because of the lack of
> convergence near 180°."

### M5-5 (p.184)
| E | 563 | (S3SR) | | |
| E | 1794 | −164 | 0.284 | 1.89 |
| M | 1965 | +7 | 0.315 | 2.49 |
| E | 2580 | +622 | 0.177 | 1.10 |
| E | 3676 | (3FR) | | |
| M | 4292 | | | |
| E | 4462 | | | |
> "The 'a' version, at least, does not exist, because it intersects Mars twice
> in 32 years."

### M5-6 (p.184)
| E | 451 | (FR)(HR)(FR)(S1SR) | | 1.46, 2.03 |
| E | 1906 | −52 | 0.237 | 2.70 |
| M | 2007 | +49 | 0.316 | 11.0 |
| E | 2580 | +622 | 0.210 | 1.55 |
| E | 3676 | (3FR) | | |
| M | 4250 | | | |
| E | 4350 | | | |
> "The 'a' version does not exist, because it intersects Earth."

### M6-1 (p.184)
| E | 473 | (6FR) | | 1.27, 4.90 |
| E | 2664 | −74 | 0.180 | 1.27 |
| M | 2782 | +43 | 0.315 | 14. |
| E | 3360 | +622 | 0.204 | 1.56 |
| E | 4456 | (3FR) | | |
| M | 5034 | | | |
| E | 5152 | | | |
> "All versions exist for the best approximations used."

### M6-2 (p.185)
| E | 534 | (S1SR)(3FR)(S1SR) | | 1.52, 2.05 |
| E | 2603 | −136 | 0.237 | 1.52 |
| M | 2760 | +22 | 0.314 | 5.84 |
| E | 3360 | +622 | 0.185 | 1.41 |
| E | 4456 | (3FR) | | |
| M | 5056 | | | |
| E | 5214 | | | |
> "The 'a' version exists for the best approximations used."

### M6-3 (p.185)
| E | 479 | (2FR)(HR)(2FR)(L1SR) | | .998, 3.74, 6.55 |
| E | 2657 | −81 | 0.175 | 1.22 |
| M | 2780 | +41 | 0.315 | 12.5 |
| E | 3360 | +622 | 0.202 | 1.56 |
| E | 4456 | (3FR) | | |
| M | 5036 | | | |
| E | 5158 | | | |
> "The 'a' version exists for the best approximations used, although
> rearrangement of the direct returns next to the short transfer to Mars is
> occasionally necessary."

### M6-4 (p.185)
| E | 500 | (FR)(HR)(FR)(S1SR)(2FR) | | 1.42, 3.32 |
| E | 2636 | −102 | 0.186 | 24., ∞, 1.47 |
| M | 2774 | +36 | 0.315 | 33. |
| E | 3360 | +622 | 0.196 | 1.54 |
| E | 4456 | (3FR) | | |
| M | 5042 | | | |
| E | 5180 | | | |
> "The 'a' version exists for solar system Model II."

### M6-5 (p.186)
| E | 500 | (2FR)(HR)(2FR)(S1SR) | | 1.47, 3.33, 5.83 |
| E | 2636 | −102 | 0.186 | 1.13 |
| M | 2774 | +36 | 0.315 | 33. |
| E | 3360 | +622 | 0.196 | 1.54 |
| E | 4456 | (3FR) | | |
| M | 5042 | | | |
| E | 5180 | | | |
> "The 'a' version exists for solar system Model II."

### M6-6 (p.186)
| E | 499 | (S2SR)(L1SR)(2FR) | | 1.24, 15.7 |
| E | 2637 | −101 | 0.185 | 1.32 |
| M | 2774 | +36 | 0.314 | 30. |
| E | 3360 | +622 | 0.196 | 1.54 |
| E | 4456 | (3FR) | | |
| M | 5042 | | | |
| E | 5179 | | | |
> "The 'a' version exists for solar system Model II."

### M6-7 (p.186)
| E | 499 | (S1SR)(L2SR)(2FR) | | 1.11, 15.2 |
| E | 2637 | −101 | 0.185 | 1.32 |
| M | 2774 | +36 | 0.314 | 30. |
| E | 3360 | +622 | 0.196 | 1.54 |
| E | 4456 | (3FR) | | |
| M | 5042 | | | |
| E | 5179 | | | |
> "The 'a' version exists for solar system Model II. Convergence was obtained
> only by using most of the encounter dates from M6-6 above."

### M6-8 (p.187)
| E | 506 | (FR)(S1SR)(L1SR)(2FR) | | 1.63, 22., 10.8 |
| E | 2630 | −108 | 0.193 | 1.35 |
| M | 2772 | +34 | 0.314 | ∞ |
| E | 3360 | +622 | 0.194 | 1.53 |
| E | 4456 | (3FR) | | |
| M | 5045 | | | |
| E | 5186 | | | |
> "The 'a' version apparently does not exist for solar system Model II; behavior
> was exhibited similar to that of periodic orbit M5-3."

### M6-9 (p.187)
| E | 506 | (S1SR)(L1SR)(3FR) | | 1.25, 10.8 |
| E | 2630 | −108 | 0.193 | 1.89, 1.89, 1.62 |
| M | 2772 | +34 | 0.314 | ∞ |
| E | 3360 | +622 | 0.194 | 1.53 |
| E | 4456 | (3FR) | | |
| M | 5045 | | | |
| E | 5186 | | | |
> "The 'a' version apparently does not exist, because of the similarity to M6-8
> above."

### M6-10 (p.187)
| E | 486 | (L1SR)(3FR)(L1SR) | | 0.92, 3.71 |
| E | 2651 | −87 | 0.175 | 0.92 |
| M | 2778 | +40 | 0.315 | 13.5 |
| E | 3360 | +622 | 0.201 | 1.55 |
| E | 4456 | (3FR) | | |
| M | 5038 | | | |
| E | 5165 | | | |
> "The 'a' version does not exist, because it intersects Earth twice in 64
> years."

### M6-11 (p.188)
| E | 564 | (3FR)(HR)(2FR) | | 2.18, 2.52 |
| E | 2573 | −166 | 0.286 | 1.90 |
| M | 2744 | +6 | 0.315 | 2.42 |
| E | 3360 | +622 | 0.176 | 1.09 |
| E | 4456 | (3FR) | | |
| M | 5072 | | | |
| E | 5243 | | | |
> (no trailing existence comment for M6-11)

**Appendix E row count: 15 periodic orbits (M4-1; M5-1..6; M6-1..11), each a
7-encounter collapsed listing → 105 encounter rows total.**

---

## 4. APPENDIX F — Eccentric Inclined member list (VERBATIM)

Full per-encounter listings (every E and M encounter). Columns:
**planet | date(JD last-5-digits, hyphen form) | speed(EMOS) | turn angle(deg) |
passing dist(local planet radii) | trajectory type**. These are long tables; I
transcribe each orbit's FIRST encounter row, LAST row, the Mars (M) encounter
rows in full (the load-bearing Mars-side V∞/passing-distance data), and the
total row count per orbit. The complete intermediate Earth-direct-return rows are
present in the source pp.191-203 and were read; reproducing all ~330 rows
verbatim here is deferred (flagged below) — the Mars rows + endpoints are the
catalogue-relevant signature; the Earth FR/IP rows are direct-return filler at
near-constant speed.

> Honesty note: I did NOT individually re-key every one of the ~330
> Earth-direct-return rows into this note. Every TABLE was readable (clean
> typescript, all numerals legible); the decision to capture Mars rows +
> endpoints + counts rather than all filler rows is an editorial scope choice,
> not a legibility failure. If a future pass needs the full per-row Earth FR/IP
> listing it is on report pp.191-203 and is fully legible.

### M4-1a (pp.191-193) — 8 encounter "blocks", ~57 rows; first 4-0538, last 8-7334
Mars (M) encounters, verbatim (date | speed | turn° | pass.dist | type):
- 4-1164 | 0.349 | 7.5 | 1.71 | IP
- 4-3495 | 0.353 | 6.1 | 2.10 | IP
- 4-6635 | 0.282 | 10.6 | 1.80 | IP
- 4-9718 | 0.326 | 5.9 | 2.55 | IP
- 5-2896 | 0.347 | 11.2 | 1.12 | IP
- 5-5945 | 0.301 | 4.6 | 3.82 | IP
- 5-9141 | 0.392 | 6.2 | 1.66 | IP
- 6-2178 | 0.278 | 3.7 | 5.64 | IP
- 6-5374 | 0.405 | 3.5 | 2.82 | IP
- 6-8420 | 0.256 | 3.4 | 7.36 | IP
- 7-1600 | 0.393 | 4.4 | 2.35 | IP
- 7-4673 | 0.245 | 3.6 | 7.63 | IP
- 7-7822 | 0.369 | 5.9 | 1.98 | IP
- 8-0936 | 0.258 | 7.4 | 3.18 | IP
- 8-4044 | 0.341 | 6.2 | 2.19 | IP
- 8-7194 | 0.301 | 7.8 | 2.20 | IP
First row: E | 4-0538 | 0.249 | 48.8 | 1.63 | S1SR. Last row: E | 8-7334 | 0.249
(speed only; terminal repeat row).

### M5-1a (p.194) — first 4-0447, last 5-2147
Mars encounters: 4-2018 | 0.252 | 3.9 | 6.57 | IP; 4-4227 | 0.328 | 9.4 | 1.53 |
IP; 4-5881 | 0.390 | 9.6 | 1.06 | IP; 4-8176 | 0.389 | 9.0 | 1.13 | IP;
4-9830 | 0.322 | 9.2 | 1.63 | IP; 5-2035 | 0.249 | 4.1 | 6.41 | IP.
First: E | 4-0447 | 0.188 | 79.2 | 1.14 | FR. Last: E | 5-2147 | 0.188. ~31 rows.

### M5-1b (p.195) — first 4-1233, last 5-2932
Mars: 4-2821 | 0.261 | 9.8 | 2.32 | IP; 4-4995 | 0.286 | 9.2 | 2.04 | IP;
4-6656 | 0.350 | 5.7 | 2.27 | IP; 4-8937 | 0.377 | 10.3 | 1.04 | IP;
5-0589 | 0.362 | 9.7 | 1.22 | IP; 5-2851 | 0.308 | 4.0 | 4.32 | IP.
First: E | 4-1233 | 0.208 | 55.2 | 1.90 | FR. Last: E | 5-2932 | 0.208. ~31 rows.

### M5-1c (p.196) — first 4-0235, last 5-1934
Mars: 4-1942 | 0.370 | 6.8 | 1.69 | IP; 4-3599 | 0.296 | 10.5 | 1.65 | IP;
4-4780 | 0.252 | 5.4 | 4.68 | IP; 4-7462 | 0.268 | 3.3 | 6.86 | IP;
4-9694 | 0.343 | 9.5 | 1.37 | IP; 5-1347 | 0.389 | 10.1 | 1.00 | IP.
First: E | 4-0235 | 0.174 | 75.8 | 1.46 | FR. Last: E | 5-1934 | 0.174. ~31 rows.

### M5-1d (p.197) — first 4-1025, last 5-2724
Mars: 4-2710 | 0.390 | 9.8 | 1.03 | IP; 4-4363 | 0.338 | 8.8 | 1.53 | IP;
4-6589 | 0.261 | 4.2 | 5.57 | IP; 4-8274 | 0.252 | 7.1 | 3.50 | IP;
5-0457 | 0.302 | 9.4 | 1.79 | IP; 5-2114 | 0.374 | 7.6 | 1.48 | IP.
First: E | 4-1025 | 0.176 | 76.9 | 1.39 | FR. Last: E | 5-2724 | 0.176. ~31 rows.

### M5-1e (p.198) — first 3-9662, last 5-1360
Mars: 4-1202 | 0.318 | 3.8 | 4.27 | IP; 4-3468 | 0.365 | 10.1 | 1.14 | IP;
4-5121 | 0.374 | 10.2 | 1.07 | IP; 4-7399 | 0.342 | 5.1 | 2.67 | IP;
4-9061 | 0.281 | 11.1 | 1.72 | IP; 5-1234 | 0.264 | 7.3 | 3.07 | IP.
First: E | 3-9662 | 0.214 | 59.9 | 1.55 | FR. Last: E | 5-1360 | 0.215. ~31 rows.

### M5-2a (p.199) — first 4-0539, last 5-2238
Mars: 4-1968 | 0.266 | 8.3 | 2.67 | IP; 4-4245 | 0.312 | 1.8 | 9.38 | IP;
4-5853 | 0.398 | 2.4 | 4.33 | IP; 4-8203 | 0.399 | 1.9 | 5.33 | IP;
4-9812 | 0.307 | 2.1 | 8.37 | IP; 5-2082 | 0.258 | 9.0 | 2.59 | IP.
First: E | 4-0539 | 0.281 | 45.7 | 1.42 | FR (note an HR appears as the 2nd
trajectory type, 4-0904). Last: E | 5-2238 | 0.281. ~31 rows.

### M5-2b (p.200) — first 4-1340, last 5-3039
Mars: 4-2793 | 0.250 | 4.2 | 6.10 | IP; 4-5017 | 0.272 | 1.8 | 12.6 | IP;
4-6618 | 0.380 | 6.5 | 1.67 | IP; 4-8957 | 0.372 | 2.4 | 4.93 | IP;
5-0569 | 0.350 | 2.4 | 5.66 | IP; 5-2894 | 0.343 | 10.0 | 1.31 | IP.
First: E | 4-1340 | 0.262 | 43.8 | 1.73 | FR (HR at 4-1705). Last: E | 5-3039 |
0.262. ~31 rows.

### M5-2c (p.201) — first 4-0240, last 5-1939
Mars: 4-1974 | 0.392 | 3.9 | 2.71 | IP; 4-3580 | 0.282 | 1.9 | 11.1 | IP;
4-5818 | 0.246 | 4.6 | 5.83 | IP; 4-7411 | 0.297 | 10.6 | 1.63 | IP;
4-9712 | 0.328 | 2.2 | 6.91 | IP; 5-1322 | 0.389 | 1.9 | 5.89 | IP.
First: E | 4-0240 | 0.176 | 78.8 | 1.31 | FR (HR at 4-2484, 4-6353). Last:
E | 5-1939 | 0.176. ~31 rows.

### M5-2d (p.202) — first 4-1019, last 5-2718
Mars: 4-2734 | 0.392 | 1.9 | 5.58 | IP; 4-4345 | 0.323 | 2.5 | 6.31 | IP;
4-6641 | 0.287 | 13.6 | 1.32 | IP; 4-8237 | 0.246 | 5.3 | 4.96 | IP;
5-0477 | 0.287 | 1.4 | 14.2 | IP; 5-2081 | 0.396 | 4.5 | 2.27 | IP.
First: E | 4-1019 | 0.179 | 78.1 | 1.29 | FR (HR at 4-3252, 4-7158, 5-1022).
Last: E | 5-2718 | 0.179. ~31 rows.

### M5-2e (p.203) — first 3-9740, last 5-1439
Mars: 4-1160 | 0.354 | 9.0 | 1.39 | IP; 4-3487 | 0.355 | 2.2 | 5.92 | IP;
4-5099 | 0.367 | 2.2 | 5.62 | IP; 4-7436 | 0.374 | 6.5 | 1.76 | IP;
4-9041 | 0.268 | 2.2 | 10.4 | IP; 5-1264 | 0.252 | 3.0 | 8.42 | IP.
First: E | 3-9740 | 0.261 | 57.8 | 1.11 | FR (HR at 4-0105, 4-8312). Last:
E | 5-1439 | 0.261. ~31 rows.

**Appendix F row count: 11 orbits (M4-1a; M5-1a..e; M5-2a..e). M4-1a ≈ 57 rows;
each M5-* version ≈ 31 rows → ~370 encounter rows total in source. Fully
transcribed above: all endpoints + all 70 Mars-encounter rows (the
catalogue-load-bearing signature); intermediate Earth FR/IP rows read but not
re-keyed (see honesty note).**

---

## 5. Per-orbit v4.2 provenance flags (NO writebacks)

Applies to all rows. Common to every orbit:
- **ephemeris / source_ephemeris**: Rall's OWN model — patched-conic
  (Ch.1.2.1) on a planetary model. Appendix E = circular-coplanar **Model I.B**
  (correct a, P; e=i=0). Appendix F = eccentric-inclined **Model III.B**
  (constant mean elements from a published series, "B" periodicity modification).
  **NEVER a DE / numerically-integrated kernel.** Tag literally:
  `source_ephemeris: patched-conic mean-element (Rall TE-34, 1969)`.
- **center**: heliocentric (Sun-centred) for all interplanetary arcs; flyby
  quantities (passing distance, turn angle) are planetocentric in **local
  planetary radii** (Earth radii at E, Mars radii at M). No catalogue-center
  ambiguity introduced.
- **model_assumption**: Appendix E orbits → `circular-coplanar`. Appendix F
  orbits → `analytic-ephemeris` (mean-element, NOT real DE).
- **tof_days_bounds**: **DERIVABLE (do NOT compute here)** from column-2 dates.
  Appendix E dates are JD−2440000 integers → leg ToF = difference of consecutive
  dates in days, e.g. M4-1 first E→M leg = 1194−1030 = 164 d (derivable). Total
  cycle span derivable from first→last date. Appendix F per-leg ToFs derivable
  from consecutive JD tokens once the JD-reconstruction ambiguity (§1) is
  resolved. **Flag: all ToFs are derived values; the date-encoding ambiguity in
  Appendix F must be resolved against Appendix E anchors before computing F-side
  ToFs.**

Per-orbit specifics:
| orbit | appendix | model | exists? (verbatim flag) |
|---|---|---|---|
| M4-1 / M4-1a | E & F | CC + ecc-incl III.B | "a" version exists |
| M5-1 / M5-1a..e | E & F | CC + ecc-incl III.B | all versions exist |
| M5-2 / M5-2a..e | E & F | CC + ecc-incl III.B | all versions exist (HR not optimally arranged, p.190) |
| M5-3 | E only | CC | "a" does NOT exist in CC region (the non-existence case discussed in Ch.2) |
| M5-4 | E only | CC | existence undetermined (no convergence near 180°) |
| M5-5 | E only | CC | "a" does NOT exist (intersects Mars twice in 32 yr) |
| M5-6 | E only | CC | "a" does NOT exist (intersects Earth) |
| M6-1 | E only | CC | all versions exist |
| M6-2 | E only | CC | "a" exists |
| M6-3 | E only | CC | "a" exists (occasional rearrangement needed) |
| M6-4, M6-5, M6-6, M6-7 | E only | CC | "a" exists for Model II |
| M6-8, M6-9 | E only | CC | "a" apparently does NOT exist (Model II; like M5-3) |
| M6-10 | E only | CC | "a" does NOT exist (intersects Earth twice in 64 yr) |
| M6-11 | E only | CC | (no existence comment given) |

---

## 6. Catalogue-eligibility recommendation per orbit (NO writebacks)

These ARE genuine repeating periodic Earth-Mars orbits (the defining problem
class), so the existing/converged ones are **catalogue-eligible** subject to the
standard provenance tags and a **superseded-by-modern-work** flag (1969
patched-conic mean-element; modern Byrnes/Russell/Casalino families supersede the
dynamics, but the dates/identifiers remain primary-source members).

- **M4-1 (+ M4-1a)** — ELIGIBLE. The headline 4-synodic-period family member
  (oldest in our problem space). Strongest candidate: exists in BOTH circular
  coplanar (App E) AND eccentric inclined (App F). multi-arc (2 E-M-E round trips
  + Earth direct returns), `circular-coplanar` (E) / `analytic-ephemeris` (F).
- **M5-1a..e, M5-2a..e** — ELIGIBLE. Five converged eccentric-inclined versions
  each; full per-encounter listings available (App F). Mark "version a/b/c/d/e"
  as distinct catalogue rows or as one family with 5 members — these are
  re-phasings (different oppositions), not different families.
- **M5-3, M5-5, M5-6, M6-8, M6-9, M6-10** — NOT eligible as flown orbits (author
  states the "a"/eccentric version does NOT exist / intersects a planet). Record
  as **negative results** only if the catalogue tracks non-existent schemes;
  otherwise exclude. Do NOT promote.
- **M5-4** — DEFER (existence undetermined; no convergence near 180°).
- **M6-1, M6-2, M6-3** — ELIGIBLE (CC exists; eccentric not listed in App F → no
  per-encounter member data, only the collapsed App-E signature). 6-synodic
  period, longer/less-efficient than M4/M5; lower priority.
- **M6-4..M6-7** — MARGINAL (exist only for Model II, the exactly-periodic
  eccentric model, not the realistic Model III). Tag `analytic-ephemeris (Model
  II, exactly-periodic)`; lower confidence.
- **M6-11** — DEFER (no existence verdict printed).

**Golden-test note:** Appendix E/F values are computed (patched-conic) outputs of
Rall's OWN program, NOT independently-sourced — per the repo rule
(`feedback_golden_tests_sourced_only`) they are **member-data provenance anchors,
not golden EXPECTED values**. They may seed/initialise our corrector
(continuation start points, exactly as the report intends — App E dates are "the
basic starting points in the search for eccentric inclined periodic orbits"),
but must not become a circular golden.

---

## 7. Unreadable-tables ledger (honest, per table)

- Appendix E (pp.182-188): **all rows fully legible.** No occlusion, no ambiguous
  glyphs. The "∞" symbols in passing-distance columns (M5-3, M6-4, M6-8, M6-9)
  are the report's own (a passing distance going to infinity = the relevant
  direct-return locus does not constrain that flyby); transcribed verbatim.
- Appendix F (pp.191-203): **all rows fully legible.** Editorial scope choice
  (Mars rows + endpoints re-keyed, Earth FR/IP filler not) is documented in §4 —
  NOT a legibility failure.
- **One genuine ambiguity, flagged not silently resolved**: the Appendix-F JD
  last-five-digits-with-hyphen reconstruction rule vs Appendix-E's clean
  "JD−2440000" (see §1). The printed TOKENS are unambiguous; the JD they
  reconstruct to has a documented ambiguity. Resolution path: anchor Appendix-F
  tokens to the Appendix-E "JD−2440000" integers (same orbit, same epoch family).
  Treat any JD conversion as a derived value to be confirmed.

---

## 8. Summary

- Appendix E: **15 circular-coplanar periodic orbits** (M4-1; M5-1..6; M6-1..11),
  105 collapsed encounter rows. Date convention: **JD − 2440000**, nearest day.
- Appendix F: **11 eccentric-inclined orbits** (M4-1a; M5-1a..e; M5-2a..e),
  ~370 full encounter rows (all 70 Mars rows + all endpoints transcribed). Date
  convention: **last-5-digits-of-JD, hyphen-separated** (ambiguous reconstruction
  — flagged).
- All quantities in the report's own units: speed = **EMOS** (Earth Mean Orbital
  Speed, ≈29.77 km/s); passing distance = **local planetary radii**; turn angle =
  **degrees**.
- Provenance for any future catalogue row: patched-conic mean-element (Rall
  TE-34, 1969), heliocentric, multi-arc; NEVER a DE kernel; superseded-by-modern
  -work flag applies.
- ET-vs-UT for the 1969 JD usage is **unstated** → conversion caveat recorded.

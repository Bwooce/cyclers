# Digest: Porter & Canup 2023, "Orbits and Masses of the Small Satellites of Pluto" (#960, wanted-list row 60)

S. B. Porter and R. M. Canup (Southwest Research Institute), The Planetary Science Journal 4:120 (14 pp), 2023 July,
doi 10.3847/psj/acde77. Received 2022 May 13, revised 2023 May 15, accepted 2023 May 22, published 2023 July 11.
Open access, CC BY 4.0.
- Proposed filename: `cyclers_pdf/papers/porter-canup-2023-orbits-masses-small-satellites-pluto-psj-4-120-doi-10.3847-psj-acde77.pdf`.
  - Publisher PDF, 14 pp., with a text layer, md5 99f33bd5a57ded820339be48efb36f79. PDF page n = journal page n.
- **How it was read:** the whole text layer (`pdftotext -layout`), start to end. Pages 10, 11 and 12 (Figure 1, Table 3,
  Table 4, the mass text) were also read as 200 dpi page images, cell by cell. Every number transcribed agreed in both.
  Tables 1 and 2 are excerpts of machine-readable tables (five rows from each end) and are not transcribed. The figures
  were looked at only for the numbers printed in them (Figure 1 prints the mass summary).
- **Companion `porter-canup-2023-tables.yaml`**: Table 3, Table 4, the masses (Figure 1 and text) and the text-only values.
  **Checks `porter-canup-2023-compare.py` / `.out`** (ours): mass and element comparison with Brozović & Jacobson 2024
  (B-J), the project registry, and Table 3 against `plu060.bsp`.

## 0. Verdict

**An independent, clean MCMC fit of the four small moons. It agrees with PLU060 (B-J 2024) within the quoted errors
for Nix and Hydra, and it gives tighter upper bounds for Styx and Kerberos. It does not fit Pluto or Charon, and it
prints no resonance argument.** So it is a cross-check for the #998 small-moon lane, not a source for Pluto or Charon
constants and not a source for Laplace-argument numbers.

- **Masses (GM, km^3/s^2, 1-sigma).** P&C: Nix (1.7455 +- 0.3472)e-3 and Hydra (2.0114 +- 0.2018)e-3 (Figure 1). Styx
  <0.026e-3 and Kerberos <0.046e-3 (Figure 1; text and abstract round to 0.03e-3 and 0.05e-3).
- **Against B-J 2024 (ours, `.out`):** Nix is +0.246e-3 higher (0.39 combined sigma), Hydra +0.001e-3 (0.00 sigma). The
  PLU060 kernel gives Nix 1.4962e-3 and Hydra 2.0080e-3. So Hydra agrees to 0.3 percent. Nix differs by 17 percent
  but well inside both errors. P&C is more precise than B-J for both moons (0.35 and 0.20 against 0.52 and 0.27).
- **Styx and Kerberos.** The P&C 1-sigma bounds (0.026e-3, 0.046e-3) are below the values that B-J HELD in PLU060
  (0.04e-3, 0.06e-3) by factors 1.5 and 1.3. B-J's own bounds (0.33e-3, 0.17e-3) are 13 and 3.7 times looser. All
  these are bounds from the posterior tail, not measurements.
- **Project registry** (`core/satellites.py`: Nix 0.0015, Hydra 0.0020, Styx 0.0003, Kerberos 0.0002):
  - Nix and Hydra agree with P&C to 0.71 and 0.06 P&C sigma. No change is needed.
  - Styx 0.0003 is 11.5 times the P&C bound and 7.5 times the B-J held value. Kerberos 0.0002 is 4.4 times the P&C
    bound and 3.3 times the held value. P&C is a second, independent source that says the registry limits (which are
    the B-J 1-sigma bounds, stored as if they were values) are too large. This repeats the proposal in the B-J digest:
    note it at the entries, or store a held value with the bound as metadata. It does not matter for massless targets.
- **Orbits.** Table 4 (mean osculating, 1000 yr) matches B-J Table 10 inside B-J's min/max range for all four moons
  (a differs by +9, -2, +8, +11 km). Period ratios agree to 1.2e-3 (Styx) and 4e-4 or better (others), which is the
  rounding of B-J's four-digit periods. See sec. 2.
- **Warning for consumers: the Table 3 small-moon states behave as PLUTO-centred, not barycentric.** The caption says
  "relative to the Pluto system barycenter". Against `plu060.bsp` at the printed epoch (ET 490276868 s), the kernel
  state relative to the barycentre differs from Table 3 by 2127-2138 km (about 5 percent of the orbit radius, 200 to
  600 times the printed sigma). The kernel state relative to Pluto (NAIF 999) differs by only 13.1, 5.2, 5.1 and 6.9 km
  (Styx, Nix, Kerberos, Hydra; 1.2, 1.4, 0.7, 1.7 times the printed position sigma). The velocity differences also
  equal Pluto's barycentric speed (24.3 m/s) in the barycentric comparison. The Table 3 Pluto and Charon rows are
  barycentric (they match PLU060 to 3.0 and 3.7 km; PLU055 against PLU060). This is our reading of the numbers; the
  paper's text says the astrometry is relative to the barycentre. Anyone who seeds a run from Table 3 must add Pluto's
  state (the Pluto row) or take the states from the kernel. This is a PROPOSAL-level warning for #998; nothing in the
  project is known to use Table 3.
  - The velocity match is weaker than the position match: after the Pluto shift the kernel differs by 11, 28, 48 and 43
    mm/s (Styx, Nix, Kerberos, Hydra), against printed 1-sigma of 25, 6, 9 and 4 mm/s. The two solutions agree in
    position at the 1-sigma level at this epoch and disagree in velocity by up to about 10 sigma (Hydra). Two different
    fits with different force models (see sec. 1) can do that; we have not chased it.
- **Resonances.** P&C print no resonance argument, no Laplace-like combination and no Charon-resonance argument. They
  give period ratios to Pluto-Charon only, and one qualitative remark (sec. 2). The Showalter-Hamilton test of B-J
  (Styx-Nix-Hydra argument circulates, about 140 yr) has no P&C counterpart.
- Wanted list: row 60 is closed by this paper (see sec. 4 for the rows it cites).
- No catalogue row cites this paper. No catalogue change is proposed.

## 1. Content

**Data** (sec. 2). Absolute astrometry (registered to Gaia DR2, not relative to Pluto, which avoids the Pluto
albedo-feature offset):
- HST 2005-2019: ACS/WFC, ACS/HRC, WFPC2 (2005-2007) and WFC3 (2010-2019; 107 orbits; 98 percent of images in F350LP).
  Pointing from Gaia stars, 1-11.8 mas for WFC3, 0.8-2.0 mas for the best images. Smeared Tiny Tim PSFs, a Pluto-Charon
  light model subtracted for F350LP and the WFPC2 F555W images. Uncertainty includes 100 km of barycentre-location
  error. Minimum exposure times: Nix and Hydra 30 s, Kerberos 90 s, Styx 120 s. Only 55 of 128 HST images gave a useful
  Styx point.
- New Horizons LORRI 4x4 images, 2015 February to July (OPNAV, U_HAZ, SCI and SATPHOT sequences), refit to Gaia DR2
  (hundreds of stars per image). Resolution 3000 to 300 km per pixel. Uncertainty includes 200 km of spacecraft-position
  error. The LORRI 1x1 images (the last week) were not used. S/N below 3 was dropped.
- Tables 1 and 2 print the format (dRA cos dec, d dec, sigma, correlation, S/N); the full sets are machine-readable.

**Model and fit** (sec. 3).
- Runge-Kutta-Nystrom adaptive integrator (PyNBody, Brankin 1989), Pluto, Charon and four small moons only. **No solar
  or planetary perturbations** (following Youdin 2012). Pluto and Charon masses and states are **PLU055, held fixed**.
  PLU055's own numbers are not printed, except its small-moon masses (Nix 2.15e-3, Hydra 3.66e-3, Kerberos 0.454e-3,
  Styx 2e-20).
- 28 parameters: ln(GM) and a barycentric state (ICRF, at the epoch) for each small moon. The epoch is just after the
  flyby: 2015-07-16 00:01:08 TDB (ours: 36 s of leap seconds plus 32.184 s gives 68.184 s from the printed UTC; consistent).
- emcee, 300 walkers; "over a week" on a 32-core machine; 11,400 retained cases (text, p.10). Three starting families
  (PLU055 masses; 1.0 g/cm^3; 0.5 g/cm^3 which relaxed to about 0.3), then merged.
- **Priors, which shape the result:** density 0.1-3.0 g/cm^3 using the Porter et al. 2021 volumes (a hard bound, not a
  fitted quantity), and a stability gate (sqrt(1-e^2) cos I within 1e-3 of the start; 10 yr back, then 1000 yr forward).
  The paper says the bound was needed to stop the chain exploring negative masses and says an unbounded run was
  not feasible with the data and compute then. Without the density gate the PLU055 start drifted to Hydra falling and
  Kerberos rising (to unrealistic densities). So the Nix and Hydra masses are partly data and partly this prior;
  they are still consistent with B-J, who used no density prior (B-J held Styx and Kerberos at 1 g/cm^3 instead).
- A minor posterior mode near Nix 0.6e-3 and Hydra 1.0e-3 (0.8 percent of cases) is mentioned and disregarded.
- Nix and Hydra masses are positively correlated through their mutual perturbation (Figure 1). Kerberos and Styx
  posteriors are monotone, falling away from zero, which is why the paper gives 1-sigma upper bounds only.

**Results** (sec. 3 and Table 4; transcribed in the yaml).
- Masses and densities: Nix 1.030 +- 0.203, Hydra 1.220 +- 0.151, Kerberos <2.03, Styx <2.20 g/cm^3 (Figure 1).
- Positional uncertainty at the epoch: 4 km (Nix, Hydra), 8 km (Kerberos), 11 km (Styx).
- Reduced chi^2 of the nominal solution: 0.376 (Nix), 0.331 (Hydra), 0.596 (Kerberos), 0.676 (Styx). The HST Styx points
  are not overfit (1.566); the other HST sets are (0.05-0.5). Against PLU055 the HST chi^2 is 4.0 (Nix) and 4.3 (Styx).
- Table 4 (mean osculating over 1000 yr, every 0.1 yr, whole cloud): T/T_PC, a (km), e, I (deg) for each moon.
  Kerberos e = 0.009901 +- 0.000006, below the Youdin 2012 limit of 0.02. Styx and Nix I are very small (0.041, 0.025
  deg). Hydra 0.282 and Kerberos 0.421 deg are the high-inclination pair.
- Discussion (sec. 4): Nix and Hydra are less dense than Charon (1.70 g/cm^3), about 1.1 g/cm^3; this implies ice-rich
  moons with low porosity (below about 30-35 percent) and supports formation from the Charon-forming impact debris
  (Canup 2011). Ice-fraction equations (1)-(4) and Figure 3 are for Nix only.

## 2. Comparison with Brozović & Jacobson 2024 and the project registry (ours; `.out`)

| Item | P&C 2023 | B-J 2024 | Registry / kernel |
|---|---|---|---|
| Nix GM (1e-3) | 1.7455 +- 0.3472 | 1.50 +- 0.52 | registry 1.5; kernel 1.4962 |
| Hydra GM (1e-3) | 2.0114 +- 0.2018 | 2.01 +- 0.27 | registry 2.0; kernel 2.0080 |
| Kerberos GM (1e-3) | <0.04568 | <0.17 (held 0.06) | registry 0.2; kernel 0.0604 |
| Styx GM (1e-3) | <0.02602 | <0.33 (held 0.04) | registry 0.3; kernel 0.0405 |
| Nix density (g/cm^3) | 1.030 +- 0.203 | 0.88 +- 0.31 | |
| Hydra density | 1.220 +- 0.151 | 1.21 +- 0.19 | |
| Pluto, Charon GM | not fitted (PLU055) | 869.3 +- 0.4, 106.1 +- 0.3 | project: 106.1 with 975.5 |

**Mass notes.**
- Sum of Nix and Hydra GM: P&C 3.757e-3, B-J 3.510e-3, kernel 3.504e-3; P&C is 7 percent higher, almost all in Nix.
- The P&C Nix value is within 0.4 sigma of B-J. Our reading is that the two fits agree and the combination
  (P&C more precise for Nix and Hydra) is the better estimate; we did not combine them and propose no change to the
  kernel-based values.
- The abstract says Nix and Hydra are "60% and 63%" of the Brozović 2015 means (3.0e-3 and 3.2e-3). With the printed
  1.74 and 2.01 the ratios are 58 and 63 percent; 60 percent uses the abstract's rounded 1.8.

**Elements (Table 4 against B-J Table 10, mean osculating).**

| Moon | a P&C (km) | a B-J (km) | diff | ratio P&C | ratio B-J (P/6.387221) | e P&C / B-J | I P&C / B-J (deg) |
|---|---|---|---|---|---|---|---|
| Styx | 43171.97 | 43163 | +8.97 | 3.269409 | 3.2706 | 0.02477 / 0.025 | 0.0411 / 0.049 |
| Nix | 49337.76 | 49340 | -2.24 | 3.994111 | 3.9939 | 0.01535 / 0.015 | 0.0249 / 0.022 |
| Kerberos | 58279.64 | 58272 | +7.64 | 5.127685 | 5.1274 | 0.00990 / 0.010 | 0.4210 / 0.422 |
| Hydra | 65186.17 | 65175 | +11.17 | 6.065662 | 6.0652 | 0.00857 / 0.009 | 0.2819 / 0.280 |

- Every P&C a lies inside B-J's min/max range (a range over time; P&C's +- is the spread over the MCMC cloud, so
  the two uncertainties are not comparable).
- These P&C ratios are the B-J Table 10 (mean osculating) form. The "registry ratios" 3.989 and 6.067 of the #998
  note are the same kind of number. In the B-J Table 11 (precessing ellipse, true mean motion) form the ratios are
  3.156, 3.891, 5.037, 5.981. The P&C text (p.11) says Nix and Hydra are "slightly closer to the 4:1 and 6:1 mean
  motion resonances" than earlier solutions and that Styx and Kerberos are farther: Nix 3.994111 and Hydra 6.065662
  against B-J's 3.9939 and 6.0652 (our values from B-J's printed periods). The statement is about the osculating-mean
  ratio only.
- **No Laplace-like or Charon-resonance argument is printed.** For the record only: the combination 3nS - 5nN + 2nH
  built from the P&C Table 4 ratios (with n_C = 360/6.387221) is -0.255 +- 0.0006 deg/day, and from the B-J Table 10
  periods -0.276 deg/day. B-J's value from the true mean motions (Table 11) is -7.3e-3 deg/day. The 40-times difference
  shows that the osculating-mean "period" is not a mean motion. Do NOT use P&C's Table 4 ratios to evaluate a
  resonance argument; use the kernel or B-J Table 11.
- Pluto-Charon period quoted in sec. 1 is 6.3872273 +- 0.0000003 d (Buie 2012), against B-J Table 10 6.387221 d (a
  mean osculating value, different definition) and the kernel header n = 56.36253 deg/day (6.387221 d). The 1e-6
  relative difference is below anything the #998 lane needs.

## 3. Internal inconsistencies in the paper (kept visible)

- Nix GM: text 1.74 +- 0.35 (p.10); Figure 1 1.74550 +- 0.34724; abstract 1.8 +- 0.4; B-J Table 8 prints 1.75 +- 0.35.
  The B-J value is the Figure 1 value rounded. Nix density 1.031 +- 0.204 (text) against 1.030 +- 0.203 (Figure 1);
  Hydra 1.220 +- 0.150 (text) against 0.151 (Figure 1).
- Styx density bound: Figure 1 <2.20; text (p.11) <2.1; abstract <2.1 (for both). Kerberos <2.03 (Figure 1), <2.0 (text).
- Table 3 caption says barycentric states; the numbers behave Pluto-centred (sec. 0).
- Typos: "0.1-3.0 g m-3" for g/cm^3 (p.9), "0.5 kg cm-3" (p.9), "Sytx" (abstract), "Ephermeris" (p.8).
- Table 4 note says the ratio is "period of the satellite over the orbital period of Pluto-Charon"; the sec. 1 values
  (3.269, 3.994, 5.128, 6.066) are those same numbers, rounded.
- P&C say Kenyon and Bromley 2022 give densities "1.4" (p.3) and "<1.4" (p.11).

## 4. Citation mining

About 72 references. `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md` hold none of the Pluto-system ones.

| Reference | Held? | Wanted list |
|---|---|---|
| Brozović, Showalter, Jacobson & Buie 2015, Icarus 246:317 (PLU043; the only prior full solution) | NOT held | row 28 (open) |
| Showalter & Hamilton 2015, Nature 522:45 | NOT held | row 27 (open) |
| Youdin, Kratter & Kenyon 2012, ApJ 755:17 (circumbinary stability, Nix and Hydra upper limits; P&C copy its no-solar-perturbation model) | NOT held | no row; proposed (low priority, #998 background) |
| Kenyon & Bromley 2019, AJ 158:69 (1 Gyr stability simulations; a different paper from Kenyon & Bromley AJ 157:79, the "Pluto-Charon Sonata" lead in the B-J digest) and 2022, AJ 163:238 (masses and densities of all four moons from stability) | NOT held | no row; proposed (low priority) |
| Bromley & Kenyon 2021, AJ 161:25; Woo & Lee 2020, AJ 159:277 (long-term dynamics of the Pluto system, named by P&C as the place for it) | NOT held | no row; low priority |
| Walsh & Levison 2015, AJ 150:11 (formation near resonances) | NOT held | no row; formation, not needed |
| Giuliatti Winter et al. 2013, 2014 (stable orbits inside Charon) | NOT held | no row; not needed |
| Porter et al. 2021, "The Pluto System After New Horizons" ch. (volumes: Nix D 36.5, Hydra D 36.2 km) | NOT held | no row; the B-J digest already notes these sizes |
| Buie, Tholen & Grundy 2012, AJ 144:15 (Pluto-Charon orbit, period 6.3872273 d) | NOT held | no row; observational |
| Tholen et al. 2008, AJ 135:777; Stern et al. 2015 (Science 350) and 2006; Weaver et al. 2016 (Science 351) and 2006; Canup 2005, 2011; Canup, Kratter & Neveu 2021 | NOT held | no row; observational or formation |
| Lauer et al. 2018 (New Horizons hazard search); Lindegren et al. 2018 (Gaia DR2); Foreman-Mackey et al. 2013 (emcee); Annex et al. 2020 (spiceypy) | NOT held | no row; method |

The two rows that matter for #998 are rows 27 and 28, already open. P&C do not cite Jacobson et al. 2015 (the PLU055
DPS abstract); they cite only the NAIF `plu055.cmt` file (footnote 1), which the project could read if needed.
The B-J 2024 paper (held) is not cited here, because P&C predate it; B-J cite P&C.

## 5. Proposals (PROPOSALS only; nothing in the repository was edited)

1. For the registry Styx and Kerberos GM: record next to the stored limit that P&C 2023 gives lower 1-sigma bounds
   (0.026e-3 and 0.046e-3) and that PLU060 held 0.04e-3 and 0.06e-3. A flyby-body use of these two moons should take the held values.
2. If any #998 script reads Table 3 states from this paper, add Pluto's state first (sec. 0), or use the kernel.
3. Nix GM: keep 1.5e-3 (kernel); P&C 1.75e-3 is a second value 0.39 sigma away. Not worth a change.
4. For resonance work cite B-J Table 11 and the kernel; do not take ratios from P&C Table 4.

*Proposed filing: `cyclers_pdf/papers/porter-canup-2023-orbits-masses-small-satellites-pluto-psj-4-120-doi-10.3847-psj-acde77.pdf`. Tables: `data/sources/porter-canup-2023-tables.yaml`; check script and output filed beside the PDF.*

*Filed as `cyclers_pdf/papers/porter-canup-2023-orbits-masses-small-satellites-pluto-psj-4-120-doi-10.3847-psj-acde77.pdf`. Tables: `data/sources/porter-canup-2023-tables.yaml`; comparison script filed beside the PDF.*

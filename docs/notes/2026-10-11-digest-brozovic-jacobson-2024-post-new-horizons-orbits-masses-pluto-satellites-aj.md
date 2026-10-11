# Digest: Brozović & Jacobson 2024, "Post-New-Horizons Orbits and Masses for the Satellites of Pluto" (#960, wanted-list row 61)

M. Brozović and R. A. Jacobson (JPL), The Astronomical Journal 167:256 (12 pp), 2024 June, doi 10.3847/1538-3881/ad39f0.
Received 2023 September 21, revised 2024 March 25, accepted 2024 April 1, published 2024 May 8. Gold open access,
CC BY 4.0, fetched from iopscience.
- Proposed filename: `cyclers_pdf/papers/brozovic-jacobson-2024-post-new-horizons-orbits-masses-satellites-pluto-aj-167-256-doi-10.3847-1538-3881-ad39f0.pdf`.
  - Publisher PDF, 12 pp., with a text layer, md5 f9d5b237871b376a2fd3afdda733dcdc (1,380,833 bytes). PDF page n = journal page n.
- **How it was read:** the whole text layer (`pdftotext -layout`). Tables 6, 7 and 10 were read on 220-250 dpi page
  crops, and Tables 8, 9 and 11 on 130 dpi page images, cell by cell. A scripted search then found all 151 numeric
  tokens of the transcribed tables in the text layer.
- **Companion `brozovic-jacobson-2024-tables.yaml`**: Tables 6-11, the text-only values (pole, resonance
  combinations, rates) and the plu060.bsp header GMs (marked as kernel values, not printed in the paper).
  Checks: `checks/check_pluto.py` / `.out` (masses, mu, Kepler, Table 7 geometry, Jbara mu), `checks/check_plu060_kernel.py` /
  `.out` (Table 7 against plu060.bsp, kernel comment area), `checks/check_yaml_textlayer.py` / `.out`. All values in
  these files are ours.

## 0. Verdict

**This is the paper that documents JPL's PLU060 solution, the satellite ephemeris in the project's pinned kernel
`plu060.bsp`.** It is a joint fit of Charon, Styx, Nix, Kerberos and Hydra to the astrometry from 1965 to 2019 (plates,
speckle, mutual events, HST, New Horizons LORRI, Gaia FPR). It gives the masses, the epoch state vectors, two kinds
of mean elements and the Charon pole. It also finds that the Styx-Nix-Hydra argument circulates (period about 140 yr).
So the Showalter-Hamilton (2015) Laplace-like resonance is not supported with the current masses.

Kernel identity, checked by us: plu060.bsp gives every Table 7 state at 2013-01-01 00:00 TDB (ET 410270400 s) to
1e-8 m and 2e-9 mm/s. The kernel comment names "PLU060", is dated BVS/NAIF 2024-04-03 (SATEPHGEN 21-Mar-2024), and
covers 1800-2199. Its pole constants (133.0073 deg, -6.2449 deg, 0.0036 and 2.5e-5 deg/century) and Charon mean motion
(56.36253 deg/day) are this paper's. The NAIF `spk/satellites/` listing (2026-10-11) has only `plu060.bsp` for
Pluto. Nothing newer exists, so the project kernel is current. On naming: the series runs PLU043 (Brozović et al.
2015, used for New Horizons navigation), PLU055 (Jacobson et al. 2015 DPS abstract; used by Porter & Canup 2023),
and PLU060 (this paper). There is no "PLU06x" successor. This paper does not supersede the project's kernel: it is
the published source of it.

**For the project (PROPOSALS only):**
- **Constants: the project's Pluto-Charon mu is a mix of two sources. The difference is small.**
  - The project uses GM_Charon 106.1 (this paper) with the system GM 975.5 (`satellites.py PRIMARIES["Pluto"]`,
    `pluto_charon_realeph.py`). The 975.5 comes from DE440 GM9 and is also Brozović et al. 2015's value. That
    gives mu = 106.1/975.5 = 0.1087647 (`pluto_charon_kk_sweep.PC_MU`).
  - This paper: mu = 106.1/(869.3 + 106.1) = 0.1087759. The kernel's full-precision GMs give 106.10114/975.42726 =
    0.1087740.
  - The project's mu is 0.0000093 below the kernel mu (relative 8.5e-5) and 0.000011 below the printed mu (relative
    1.0e-4). The 1-sigma of mu from GM_Charon alone is 3.1e-4. So the difference is 0.04 sigma.
  - The implied Pluto GM is 869.4 against the printed 869.3 +- 0.4.
  - Proposal: when a constant is next edited, take the Pluto-Charon pair from the kernel header (GM_P 869.326,
    GM_C 106.101) or from Table 8. Keep DE440 GM9 975.5 for the heliocentric body. This blocks nothing: the
    difference is far below the 12 percent barycentric offset of the #998 one-centre model.
- **Registry period question (#1022 sec. 3, #1022b, #1022c): this paper settles it.** The JPL SSD page pairs two
  different element sets.
  - The page's a values (43200, 49300, 58300, 65200 km) are Table 10's MEAN OSCULATING a rounded to 100 km
    (43163, 49340, 58272, 65175).
  - Its periods (20.16, 24.85, 32.17, 38.20 d) are Table 11's PRECESSING-ELLIPSE P.
  - Each table agrees with itself:
    - Table 10: Kepler P from its a, with the kernel GM_P + GM_C, matches its own P to 0.05 percent (ours).
    - Table 11: P = 360/n. Its a is NOT the Kepler radius of n. Kepler gives a period 0.87, 0.67, 0.45 and
      0.36 percent too long. This is about the size of the quadrupole estimate in #1022 sec. 3.
  - So the registry's Kepler-from-a mean motions reproduce the osculating-mean periods (20.91, 25.49, 32.78,
    38.76 d against Table 10's 20.89, 25.51, 32.75, 38.74 d). The "registry ratios" 3.989 and 6.067 in the #998
    note are the paper's osculating ratios 3.99 and 6.07.
  - The #998 plu060 fitted periods (20.16195, 24.85472, 32.16798, 38.20202 d) are the Table 11 periods to the
    printed digits. B-J sec. 4.1 says so itself: "the definition of elements matters".
  - Proposal: cite Table 10 / Table 11 in #1022 sec. 3 as the sourced explanation, which replaces "remainder
    unexplained". For any registry option (a), the sourced period is Table 11's P (or n), not a fit of ours.
- **Styx and Kerberos GM in the registry.** The registry stores 0.0003 and 0.0002. These are the rounded Table 8
  1-sigma upper bounds (<0.33e-3, <0.17e-3). PLU060 was integrated with HELD values 0.04e-3 and 0.06e-3 (Table 8
  note; kernel 4.045e-5 and 6.038e-5). The registry values are 7.5 and 3.3 times the held values (ours).
  - Proposal: add a note at the registry entries, or store the held values with the bound as metadata.
  - This does not matter for #998 (massless targets). It overstates the bending of any lane that uses Styx or
    Kerberos as a flyby body.
  - Radii: registry Nix 18.0 and Hydra 18.5 km (the SSD page values). Table 9 uses diameters 36.5 and 36.2 km
    (R 18.25 and 18.1), and Kerberos D 12 km against the registry's 6.0 km radius. These are small; noted only.
- **#998 model support.** Table 11's free e is below 1e-5 for Styx, Nix and Kerberos, and 0.006 for Hydra. The
  osculating e of 0.009-0.025 in Table 10 is forced by the binary, not free. Mean inclinations to Charon's plane:
  Styx 0.049, Nix 0.022, Kerberos 0.42, Hydra 0.28 deg. Both support the coplanar circular seed model of #998. The
  #998 normal angles (0.006, 0.010, 0.348, 0.260 deg, over 3 years) agree in order. The Table 7 states give 0.052,
  0.022, 0.424, 0.280 deg (ours).
- **Literature-gate anchor (#1052, row 4).** The anchor `Brozovic et al.` cites "Icarus 246:317-329 (2015);
  Brozovic et al., AJ 163:241 (2022)". Crossref (2026-10-11): AJ 163:241 (2022), doi 10.3847/1538-3881/ac617f, is
  Brozović & Jacobson, "Orbits of the Irregular Satellites of Uranus and Neptune". It is a real paper about the
  wrong system. This paper does not cite it.
  - Proposal: replace the second reference with AJ 167:256 (2024), doi 10.3847/1538-3881/ad39f0. The label
    {ephemeris} fits this paper exactly. Do not edit `literature_check.py` here.
- **The Showalter-Hamilton anchor (#1052, row 3, D3).** This paper finds the Styx-Nix-Hydra argument circulating
  with the present masses (Fig. 9, period about 140 yr). It also finds Charon-Kerberos and Charon-Hydra only near
  fourth- and fifth-order resonances (residuals -0.03 and +0.03 deg/day). It does not find a cycler or a
  trajectory. That supports D3 (resonant to ephemeris) for the moons' own orbits.
- Wanted list: row 61 is closed by this paper. Two new rows are proposed (sec. 4).
- No catalogue row cites this paper. No catalogue change is proposed.

## 1. Content

**Data** (sec. 2.1, Tables 1-5). Charon: 1965-2019, 28 data rows in Table 1, residual rms from about 100 mas (plates) to 1.4 mas
(Gaia, with the Buratti 2017 albedo correction). Small moons: mostly HST astrometry as re-reduced by Porter & Canup
(2023), 2002-2019, rms at or below 1-3 mas. New Horizons LORRI residuals of 19-42 km (Fig. 4). Mutual-event light curves
(Fig. 1).

**Model** (sec. 2.2). Barycentric, ICRF, TDB, after Peters (1981). Point-mass Sun, Jupiter, Saturn, Uranus and Neptune
(Table 6 GMs; positions from DE440). Gauss-Jackson integration. Over 200 parameters, including biases. The fit
showed no sensitivity to Styx or Kerberos masses. Making them massless changes Hydra by about 3 km along track over
2005-2020.

**Masses** (sec. 3.2, Table 8).
- GM_sys 975.4 +- 0.2 (formal 0.09), Charon 106.1 +- 0.3, Pluto 869.3 +- 0.4.
- mc/mp = 0.1220 +- 0.0003; ours 106.1/869.3 = 0.12205.
- Nix (1.50 +- 0.52)e-3 and Hydra (2.01 +- 0.27)e-3.
- Kerberos below 0.17e-3 and Styx below 0.33e-3, held in the fit at 0.06e-3 and 0.04e-3.
- Densities (Table 9) reproduce from the printed sizes (ours: 1.853, 1.705, 0.883, 1.212 g/cm^3).

**Mean elements** (sec. 3.3, Tables 10, 11). Over 400 yr (1800-2200), sampled daily. Two definitions are given and
they disagree on purpose.
- Table 10, mean osculating: Styx 43163 km / 20.89 d, Nix 49340 / 25.51, Kerberos 58272 / 32.75, Hydra 65175 / 38.74.
- Table 11, precessing ellipse: 42410 / 20.16, 48689 / 24.85, 57749 / 32.17, 64719 / 38.20.
- Charon (Plutocentric) a 19595.764 km, e 1.61e-4, P 6.387221 / 6.387222 d, n 56.362530 deg/day.
- Rates: Kerberos Omega_dot -0.11 deg/day (nodal period about 3300 d, about 9 yr). Hydra varpi_dot 0.070,
  Omega_dot -0.071 deg/day (about 5100 d, about 14 yr).

**Pole and plane** (sec. 3.4). Charon pole R.A. = 133.0073 + 0.0036 T, decl. = -6.2449 + 2.5e-5 T (deg, T in Julian
centuries from J2000), from 5000 yr of integration, good to about 2 millidegrees. The solar torque precesses the pole
at 3.79e-3 deg/century (period 9.5 Myr). Ours: the Table 7 Charon state gives an osculating pole of 133.0077,
-6.2443 deg at 2013.

**Uncertainties** (sec. 3.5, Fig. 8). Covariance mapping to 2050. Mean-motion 1-sigma: Charon 9e-7, Styx 2e-5, Nix
4e-6, Kerberos 7e-6, Hydra 2e-6 deg/day. PLU043 against PLU060 at the 2015 flyby: up to 244 km along track (Styx).

**Near resonances** (sec. 4.1). The period ratios to Charon "as printed":

| Element set | Styx | Nix | Kerberos | Hydra |
|---|---|---|---|---|
| Table 10 (osculating mean) | 3.27 | 3.99 | 5.13 | 6.07 |
| Table 11 (precessing ellipse) | 3.16 | 3.89 | 5.04 | 5.98 |

Ours (ratios of the printed P): 3.2706, 3.9939, 5.1274, 6.0652 and 3.1563, 3.8906, 5.0366, 5.9807. The introduction
calls the sequence "near-resonant orbits (1:3:4:5:6)". Showalter-Hamilton's 3:4:5:6 is the Table 10 form. In the
Table 11 (true mean-motion) form, Styx is 5 percent off 3:1 and Nix 3 percent off 4:1.
- nC - 5 nK + 4 Omega_dot_K = -0.03 deg/day: ours -0.0325.
- nC - 6 nH + 4 varpi_dot_H + Omega_dot_H = +0.03 deg/day: ours +0.0275.
- 3 nS - 5 nN + 2 nH = -7.3e-3 deg/day. Ours from the rounded rates is -7.0e-3; the rounding allows several e-3.
- Phi circulates with a period of about 140 yr (ours: 141 yr). It passes 180 deg in 2010-2012 "by coincidence".

## 2. Witnesses

- Witness C (decides): the page images, cell by cell. Every cell agreed with the text layer.
- Witness A: the text layer. All 151 numeric tokens of Tables 6, 7, 8, 10 and 11 were found in it (scripted).
  Glyph losses: "L" for "...", "v" for varpi-dot and "W" for Omega-dot (Table 11), and the "x" of "x 10^-3". In the
  sec. 3.2 text, sigma_GM prints as "sGM".

## 3. Cross-checks (ours; `checks/`)

- Table 7 against plu060.bsp: identical to 1e-8 m (above). Pluto (999) is 2131.3 km from the barycentre at that epoch.
  Charon's |r| is 17462.6 km and |r|/(1 - mu) = 19593.9 km.
- Kernel GMs: Charon 106.1011388, Pluto 869.3261226, system 975.4308664. The system includes the small moons, which sum
  to 3.605e-3. The kernel's own DE440 header prints GM9 = 975.5.
- Charon: Kepler P from a = 19595.764 with GM_P + GM_C = 975.42726 is 6.3872216 d, against 6.387221 printed. With
  975.4 it would be 6.38731, so the printed 975.4 is a rounding.
- Table 6 Sun GM 132713233263.35 equals the DE440 Sun plus Mercury, Venus, EMB and Mars (difference 793222.07, sum
  793222.08). The paper does not state this.
- Table 11: 360/n reproduces P (20.162, 24.855, 32.169, 38.200).
- **Mismatches inside the paper, kept visible:**
  - The abstract says "Kerberos and Styx have 0.4 deg and 0.3 deg inclinations". Tables 10 and 11 and sec. 3.3 say
    Kerberos and Hydra. The abstract is wrong.
  - Nix GM: the text on p.8 gives 0.0014 +- 0.0005; Table 8 and the abstract give (1.50 +- 0.52)e-3; the kernel has
    1.496e-3.
  - Charon formal sigma: the text gives 0.14; Table 8 gives {0.12}.
  - Brozović 2015 system GM: the text on p.6 says 975.2 (formal and inflated); Table 8 prints 975.5.
  - Porter & Canup Nix: the text gives 0.0017 +- 0.00035 (and 0.0003); Table 8 gives (1.75 +- 0.35)e-3.
  - Sec. 4.1 says "Nix, Hydra, and Kerberos reside in a Laplace-like resonance", but the combination it quotes is
    Styx-Nix-Hydra.
  - Table 3, Hydra 2019, prints 9 points for Da cos d and 13 for Delta delta (checked on the page image; probably a typo, since every other row pairs equal counts).

## 4. Citation mining

There are about 75 references. `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md` hold none of the
Pluto-system ones; Acton 1996 (SPICE) is held. Relevant rows:

| Reference | Held? | Wanted list |
|---|---|---|
| Showalter & Hamilton 2015, Nature 522:45 | NOT held | row 59 (open). Needed less now: this paper shows the claimed resonance circulates |
| Brozović, Showalter, Jacobson & Buie 2015, Icarus 246:317 | NOT held | row 60 (open). PLU043; superseded by this paper for constants; still the anchor's first citation |
| Porter & Canup 2023, PSJ 4:120, "Orbits and Masses of the Small Satellites of Pluto", doi 10.3847/psj/acde77 (Crossref 2026-10-11) | NOT held | **no row. Proposed new row:** the independent MCMC fit of the four small moons (Nix and Hydra masses, Styx and Kerberos bounds); PSJ is open access |
| Youdin, Kratter & Kenyon 2012, ApJ 755:17, "Circumbinary chaos: using Pluto's newest moon to constrain the masses of Nix and Hydra", doi 10.1088/0004-637x/755/1/17 | NOT held | **no row. Proposed, low priority:** circumbinary stability of the small-moon zone; background for #998 |
| Kenyon & Bromley 2014, AJ 147:8; Walsh & Levison 2015, AJ 150:11; Canup 2011, AJ 141:35; Canup et al. 2021 | NOT held | no row; origin and formation, not needed |
| Stern et al. 2018, ARA&A 56:357 | NOT held | no row; review (radii) |
| Park et al. 2021, AJ 161:105 (DE440/441) | NOT held (CORPUS_INDEX only mentions DE440) | no row; the kernel header carries what is needed |
| Peters 1981, A&A 104:37; Jacobson & Owen 2004, AJ 128:1412 | NOT held | no row; method |
| Ward 1975; Ward & Hamilton 2004, AJ 128:2501 | NOT held | no row; pole precession, not needed |
| Jacobson et al. 2015 (DPS abstract, PLU055) | NOT held | no row |
| Acton 1996, P&SS 44:65 | HELD (`acton` file) | n/a |
| Astrometry sources: Harrington & Christy, Buie 2006/2012, Tholen 2008 (two), Null, Olkin, Sicardy 2006/2011, Weaver 2006, Pelletier 2016, Thompson 2016, Tanga 2023, Herald 2020, Young, Beletic, Hege, Baier, Bonneau, Hetterich, Howell, Fuentes & Holman, Beauvalet 2013 | NOT held | none; observational, no row proposed |

Crossref lead (not cited here, not in the wanted list): Kenyon & Bromley 2019, AJ 157:79, "A Pluto-Charon Sonata:
The Dynamical Architecture of the Circumbinary Satellite System", doi 10.3847/1538-3881/aafa72. It covers the
stability of the region the #998 small moons live in. Proposal: one low-priority row beside Youdin 2012.

*Proposed filing: `cyclers_pdf/papers/brozovic-jacobson-2024-post-new-horizons-orbits-masses-satellites-pluto-aj-167-256-doi-10.3847-1538-3881-ad39f0.pdf`. Tables: `data/sources/brozovic-jacobson-2024-tables.yaml`.*

*Filed as `cyclers_pdf/papers/brozovic-jacobson-2024-post-new-horizons-orbits-masses-satellites-pluto-aj-167-256-doi-10.3847-1538-3881-ad39f0.pdf`. Tables: `data/sources/brozovic-jacobson-2024-tables.yaml`; check scripts filed beside the PDF.*

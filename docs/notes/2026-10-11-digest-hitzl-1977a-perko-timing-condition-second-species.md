# Digest: Hitzl 1977a, "On the Perko Timing Condition for Second Species Orbits in the Restricted Problem" (#960 batch 43)

D. L. Hitzl (Lockheed Palo Alto Research Laboratory), Astronomy & Astrophysics 54:47-55 (1977), ADS 1977A&A....54...47H.
Received May 21, 1976. No DOI: Crossref has no record (bibliographic query checked 2026-10-11; it returns only the
Celest. Mech., Acta Astronaut. and AIAA 76-158 papers of the same programme).
- Proposed filename: `cyclers_pdf/papers/hitzl-1977a-perko-timing-condition-second-species-orbits-restricted-problem-aa-54-47-ads-1977AA-54-47H.pdf`.
  - File the `ocrmypdf --force-ocr -l eng` copy, md5 bf3545358351cb44bcb3f23b452c02be (715 kB).
  - Source: image-only ADS JBIG2 scan, 9 pp. (pp.47-55), md5 ac40bb1ecfcec0742f64be57fb71ade0, fetched 2026-10-11 from
    articles.adsabs.harvard.edu with curl.
  - Page images survived the OCR. Dark-pixel fractions at 100 dpi, original against OCR copy: p.1 0.073/0.080,
    p.5 0.050/0.055, p.7 0.055/0.061 (`render_check.out`). The copy is slightly darker, with no blank pages.
- **How it was read:** all nine pages on page images (110 dpi), Table 1 on a 400 dpi crop. The OCR text layer was
  used only to find things.
- **Companion yaml:** `hitzl-1977a-tables.yaml`, Table 1 (12 rows), with the equations and a `conventions:` block.
- Checks (all values in them are ours): `recompute_table1.py/.out`, `twobody_check.py/.out`,
  `crosscheck_paperI.py/.out`, `paperI_miss.py/.out`, witness file `witness-comparison-hitzl1977a.tsv`,
  `witness_counts.out`, `tess_table1_psm6.txt`.
- Naming: Bruno 1994 calls this paper "Hitzl (1977a)". The corpus already calls the Celest. Mech. paper
  "hitzl-henon-1977" and the AIAA J paper "hitzl-1977". This digest uses "Hitzl 1977a" only for this A&A paper.

## 0. Verdict

**This paper shows that Perko's timing condition for second-species generating orbits (mu = 0) is a special case
of Hénon's 1968 consecutive-collision timing condition.** It adds one table of 12 example orbits and 13 orbit plots.
- The proof (pp.48-50): Hénon's condition F0(tau, eta) = 0 (eq. 9), rewritten in (a, e) as H0 = 0 (eq. 14), equals
  Perko's F_mk(a, e) = 0 (eq. 16) when tau = m pi, eta = k pi and the signs are sigma0 sigma1 sigma2 = - + -
  (eqs. 17-24), or - - + via Perko's second case (eqs. 25-27).
- New closed form (eq. 34): C on a consecutive-collision orbit as a function of (tau, eta) only. Critical orbits are
  the extrema of C on a characteristic, the condition G* = 0 (eqs. 35-37). This is the same condition as Hitzl &
  Hénon 1977 (Celest. Mech. 15:421) and 1977b (Acta Astronaut.), both held.
- Fig. 16: the four limiting shapes of the Moon passage (direct and retrograde "angular points") at mu = 0 and mu > 0,
  linked to Arenstorf's Types A and B.

**What it gives the project: little that is new.**
- The theory and the critical orbits are already in the held Hitzl-Hénon pair and in Hénon 1968.
- Table 1 is 12 example orbits at mu = 0, three of them critical. We reproduced every row (sec. 3). The three critical
  rows are less accurate than Hitzl & Hénon 1977, because they come from 5-digit (tau, eta).
- **Cross-link worth recording:** the captions tie two Huang 1962 orbits to Hénon families.
  - Fig. 7 (C23, retrograde): "A closely related orbit of this type appears in Huang (1962), as Figure 1".
  - Fig. 8 (C23(1), direct, critical): "An orbit of this type is presented by Huang as Figure 2".
  - Huang 1962 is digested in this same batch (`../huang/digest-huang-1962.md`). It found that Huang's retrograde
    2/3 family is the project's `#997` F4 2/3-r family, and the direct one continues to Causey 1964's 2/3 order 1.
    The C23 characteristic has a ≈ (2/3)^(2/3) = 0.763 near tau/pi = 2, eta/pi = 3. So Hitzl names the mu = 0
    generating families of both Huang 2/3 families: retrograde ↔ C23 (sigma1 = -1), direct ↔ C23(1) (sigma1 = +1).
  - This is Hitzl's caption claim ("closely related", "of this type"). We did not continue a Huang orbit to mu = 0.
    Huang's orbits pass the Moon at 0.08-0.20, so they are not near-collision orbits; the link is by shape and
    resonance, not by a computed continuation.
- **Proposals:**
  1. Record the Huang 1962 ↔ C23 attribution in the `#997` lineage notes. It gives a mu = 0 generating orbit
     (published, Hénon 1968 family) for both 2/3 families.
  2. Use Hitzl & Hénon 1977 (Celest. Mech.) for any critical-orbit golden, not this Table 1.
  3. No new code. Our eq. (7)/(9)/(33) implementation reproduces the table, and the held digests already carry it.
- For the Bruno 1994 citation (`#944` X2, `#948` R4 generating families, wanted-list row 28): this paper supplies
  attribution and the Perko-Hénon equivalence. It adds no family data beyond Hénon 1968 and Hitzl-Hénon 1977.

## 1. Content

**Setup** (pp.47-48, page images). Earth M1 (mass 1) at the origin of inertial axes; the Moon M2 on x = cos t,
y = sin t (eq. 2); collisions at t = ±tau at points P and Q (Fig. 1, from Hénon 1968). Particle on a Kepler ellipse,
eqs. (3)-(4), with the switch functions sigma0 (apse side), sigma1 (direct +1, retrograde -1), sigma2 (perigee +1,
apogee -1 at t = 0), eq. (1).
- Eq. (7): a = (1 - sigma c_tau c_eta)/s_eta², e = (sigma2 c_eta - sigma0 c_tau)/(1 - sigma c_tau c_eta), sigma = sigma0 sigma2.
- Eq. (9), Hénon's timing condition: F0 = sigma_eta (1 - sigma c_tau c_eta)^{1/2}[eta(1 - sigma c_tau c_eta) - s_eta(c_eta - sigma c_tau)] - tau s_eta³ = 0.
- Only sigma = sigma0 sigma2 enters (9) and (14); sigma1 does not (p.49).
- "solutions exist for all values of tau > 0.51500 ... and eta > 0" (p.48).

**Perko's condition** (p.49, eq. 16): F_mk(a, e) = sigma1 cos⁻¹[(1 - a(1 - e²))/e] - m pi + a^{3/2} k pi -
a^{3/2}[cos⁻¹((1 - a)/(ae)) + e(1 - ((1 - a)/(ae))²)^{1/2}] = 0, with a = (m/k)^{2/3} (eq. 18).
- Perko's figure has P and Q reversed relative to Hénon's, which forces sigma0 = -1 (eq. 17).
- Eqs. (19)-(24) map the terms one by one. The only possibility is sigma0 = -1, sigma1 = +1, sigma2 = -1.
- Perko p.203 adds a second case (eq. 25): F_mk + (m/k) pi = 0, giving signs - - + with 2l - 1 = m/k (eq. 27).
- Example: m = k = 2, l = 1 gives an orbit of family A1 (Figs. 3, 4).

**Table 1 columns** (p.51): n* is half the number of x-axis crossings in one period; x0 (eq. 28) at t = E = 0;
x1 (eq. 29) at E = pi; C by eq. (33) C = 2 sigma1 (a(1 - e²))^{1/2} + 1/a; T = 2 tau; "Critical Yes/No".

**Critical orbits** (pp.52-54). At mu = 0 a consecutive-collision orbit has k = ±infinity. An extremum of C along a
family must be a critical orbit with k = +1. So for small mu > 0 the jump in k becomes a finite stable interval
between two critical orbits, k = +1 and k = -1. Eq. (34) gives C(tau, eta) in closed form; eqs. (35)-(37) give
G* = 0. "a total of 179 critical second species orbits have been determined (Hitzl and Hénon, 1976a)" (p.54).
G* = 0 was also found independently from a first-order matching stability analysis (Hitzl and Hénon 1976b).

**Fig. 5 caption:** A0(-1) "had been found previously by Hénon and Guyot (1970) ... to be the termination of two
retrograde critical periodic families m1 and m2. These families exist over the full interval 0 ≤ mu ≤ 1 (m1) and
over the two symmetrical intervals 0 ≤ mu < mu0 = 0.327..., 1 - mu0 < mu ≤ 1 (m2)".

**Fig. 14 caption:** B1(0) is "the degenerate case when the particle describes exactly the same orbit as the Moon".

**Printed oddities** (kept as printed in the yaml).
- A0(-1) T = 1.34113 (and x0 -0.16294) follows the unrounded tau/pi = 0.2134482407 of Hitzl & Hénon 1977, not the printed 0.21345
  (which gives 1.34115). C (-0.39913 printed, -0.39914 from the rounded inputs) likewise.
- C23(1) and C24(2): a, e, x0, x1, C, T follow the 5-digit (tau, eta) and differ from Hitzl & Hénon 1977 by up to
  5.7e-4 (C23(1) x0: -0.52736 here, -0.52679 there). This is input rounding, not a misprint; the critical point is
  sensitive.
- C12: the printed (tau, eta) do not satisfy eq. (9) well (sec. 3). Unresolved.

## 2. Witnesses (Table 1)

- 96 numeric cells (tau/pi, eta/pi, a, e, x0, x1, C, T for 12 rows). Names, Fig., n*, signs and the Critical column
  were checked by eye against the tesseract output and agree.
- **A, page image** (400 dpi crop, decides).
- **B, tesseract 5** at 400 dpi, psm 6, on the same crop (not the ocrmypdf text layer): 94 agree after mapping the
  em-dash minus to "-". 2 differ:
  - C24 T: tesseract "1223946" drops the point; image 12.23946.
  - B1 T: tesseract "5.10777"; image 5.70777. T = 2 tau from the printed tau/pi 0.90842 gives 5.70777, so the image
    reading stands.

## 3. Checks (all ours)

**(a) Table 1 from its own (tau, eta)** (`recompute_table1.out`). Eqs. (7), (28), (29), (33) and T = 2 tau from the
printed tau/pi, eta/pi and signs reproduce 69 of 72 cells to the printed digits. The 3 others are A0(-1) x0, C and
T (sec. 1, rounding of tau: ours -0.16295, -0.39914, 1.34115 against printed -0.16294, -0.39913, 1.34113). Eq. (34) gives the same C as eq. (33) in all 12 rows.

**(b) Timing condition (9) at the printed (tau, eta).** |F0| is 4e-8 to 3e-5 in 10 rows, consistent with 5-digit
rounding. Two rows are larger:
- B1: F0 = 1.2e-4. The root at the printed tau is eta/pi = 1.332251 (printed 1.33223).
- C12: F0 = 6.0e-5. Here dF0/deta is small (-0.057), so the root at the printed tau is eta/pi = 2.085458 (printed
  2.08512), and the root at the printed eta is tau/pi = 0.966098 (printed 0.96601). The printed a, e, T follow
  0.96601. A digit transposition (0.96610) would fit the timing condition, but nothing else in the row supports it.
  Unresolved.

**(c) Integration** (`twobody_check.out`). At mu = 0 a CR3BP integration "at the paper's mu" is the two-body problem,
so we integrated the Kepler arc (Earth mass 1, DOP853, rtol 1e-12) from the apse state of eqs. (3)-(4) to t = tau and
compared it with the Moon at (cos tau, sin tau). We did this for all 12 rows, in the inertial frame and again in the
mu = 0 rotating frame.
- Moon miss distance: 4e-7 (B1~B1(0)), 1e-6 to 8e-5 in eight rows, 1.4e-4 (B1), 1.8e-3 (C12) and 3.0e-3 (C23(1)).
- With Hitzl & Hénon 1977's 10-digit (tau, eta) the three critical rows miss by 7e-11 (A0(-1)), 2e-8 (C23(1)) and
  1e-9 (C24(2)) (`paperI_miss.out`). So the C23(1) miss is the rounding of this table, and the orbit is right.
- Jacobi C from the rotating-frame state equals the printed C in 11 of 12 rows to the printed digits (A0(-1): rounding,
  as above).
- n*: interior x-axis crossings on (0, tau) plus one reproduce the printed n* in 12 of 12 rows.

**(d) Against the held Hitzl & Hénon 1977 (Celest. Mech. 15:421) and Hitzl 1977 (AIAA J).**
- A0(-1): every field equals Hitzl & Hénon Table II.
- C23(1), C24(2): equal to the AIAA J Table 1 values in tau/pi, eta/pi, signs and n* (and in C24(2) x0, and both x1
  to 1e-5); a, e, C, T and C23(1) x0 differ by the rounding effect of sec. 1. The AIAA J digest says "C23(1) and C24(2) were not in the digest excerpt and were not
  compared"; this closes that gap: the Celest. Mech. 10-digit values give a 0.76343, e 0.30997, x0 -0.52679,
  x1 1.00007, C 2.97130, T 12.53191 for C23(1), equal to the AIAA J table.

## 4. Citation mining

References (p.55), checked with `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md`:

| Reference | Held? |
|---|---|
| Arenstorf 1963, Amer. J. Math. 85:27 | HELD (`arenstorf-1963-periodic-solutions-...-amer-j-math-85-27-...`) |
| Arenstorf 1966, AAS/AAS Astronautics Symposium, Berkeley, pp.103-112 | NOT held. Low priority (review). |
| Barrar 1965, AJ 70:3 | NOT held, not on the wanted list. Low priority (existence of periodic orbits, theory). |
| Breakwell & Perko 1966, Prog. Astronaut. Aeronaut. 17:159 | HELD as the AIAA 65-689 preprint (`breakwell-perko-1965-...`) |
| Breakwell & Perko 1974, Celest. Mech. 9:437 | HELD |
| Bruno 1973, "Research on the Restricted Three Body Problem. III. Properties of the Solutions for mu = 0", Acad. Sci. Moscow | HELD in translation as `brjuno-1978b-...-celest-mech-18-51-...` (index line 742); same title, not checked that it is the same text |
| Darwin 1911, Scientific Papers IV | NOT held. Low priority (historic, mu = 1/11). |
| Guillaume 1973 (Tapley & Szebehely eds.), pp.161-174 | NOT held (the held Guillaume 1973 is Celest. Mech. 8:199, a different paper). Low priority. |
| Guillaume 1974, A&A 37:209 | NOT held. Low priority. |
| Guillaume 1975, Celest. Mech. 11:213 and 11:449 | HELD (both) |
| Hénon 1968, Bull. Astron. 3:377 | HELD |
| Hénon 1969, A&A 1:223 | HELD |
| Hénon & Guyot 1970, Giacaglia ed., pp.349-374 | HELD |
| Hitzl 1975a, A&A 40:147 | NOT held, not on the wanted list. Hitzl's own critical-orbit stability index paper; the Celest. Mech. digest cites it for "k = ±1". Medium-low priority; ADS free scan (1975A&A....40..147H, not verified). |
| Hitzl & Hénon 1976a (Celest. Mech., "to appear") | HELD as `hitzl-henon-1977-...-celest-mech-15-421-...` |
| Hitzl & Hénon 1976b ("in preparation") | HELD as `hitzl-henon-1977b-...-acta-astronaut-4-1019-...` |
| Huang 1962, AJ 67:304 | Fetched this batch (`../huang/`), wanted-list row 23 |
| Kevorkian & Lancaster 1968, AJ 73:791 | HELD |
| Markellos 1974a (Celest. Mech. 9:365), 1974b (Celest. Mech. 10:87), 1975 (Astrophys. Space Sci. 36:245) | NOT held, not on the wanted list. Numerical families at mu > 0; medium-low priority for second-species family data. |
| Perko 1974, SIAM J. Appl. Math. 27:200 | HELD |

- Also found on Crossref: Hitzl & Hénon, AIAA paper 76-158 (14th Aerospace Sciences Meeting), doi 10.2514/6.1976-158,
  "Periodic orbits in the restricted problem - Critical generating orbits for mu equals 0". Not held; it is the
  conference form of the Celest. Mech. paper. Low priority.
- No must-have lead. Proposal: one low-priority wanted row for Hitzl 1975a and Markellos 1974-75.

*Wanted-list row 28 (Bruno 1994 second-species sources): this paper is the "Hitzl, D. L. (1977a), A&A 54:47-55" item.*

*Filed as `cyclers_pdf/papers/hitzl-1977a-perko-timing-condition-second-species-orbits-restricted-problem-aa-54-47-ads-1977AA-54-47H.pdf`. Check scripts, outputs, witness tables and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/hitzl-1977a-tables.yaml`.*

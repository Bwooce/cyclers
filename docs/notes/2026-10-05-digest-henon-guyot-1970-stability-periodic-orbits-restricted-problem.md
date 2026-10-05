# Digest: Henon & Guyot 1970, "Stability of Periodic Orbits in the Restricted Problem" (#960)

M. Henon (Observatoire de Nice) and M. Guyot (Faculte des Sciences de Nice), "Stability of Periodic Orbits
in the Restricted Problem", in G. E. O. Giacaglia (ed.), Periodic Orbits, Stability and Resonances, Reidel,
Dordrecht, 1970, pp. 349-374, doi 10.1007/978-94-010-3323-7_33.
- Filed as `cyclers_pdf/papers/henon-guyot-1970-stability-periodic-orbits-restricted-problem-giacaglia-349-doi-10.1007-978-94-010-3323-7_33.pdf`.
  26 pages, text layer, md5 38635f4a2fdc28f7caa0cc083bdf73dd. Book page = PDF page + 348.
- I read all pages from the text layer. The table pages (pp.367-373) were rendered as page images and
  compared with the transcription.
- Row counts per family match the images exactly. The values were spot-checked on every table page.

## 0. Verdict

This is the source of the critical orbits (stability index a = +-1) of the six main symmetric families
f, g, h, i, l and m of the planar circular restricted problem, for all mass ratios mu in [0, 1]. It has a
tabulated (mu, x0, x1, C) for every critical-orbit family that the paper follows, including the Earth-Moon
value mu = 0.01215 (and 0.98785).
- Tables are printed for families m, l, h and i only.
- Families f and g follow by the symmetry mu -> 1 - mu, x -> -x (Eq. (12), p.355). Family f is the image of
  h, and g is the image of i.
- The abstract's threshold, recorded exactly: "One particular result is that if the relative mass of one
  of the bodies is less than 0.0477, all retrograde orbits around that body are stable" (p.349).
- The body text gives it as Eq. (24), "0 <= mu < 0.047 7 ...": "essentially all retrograde periodic orbits
  around the lighter body M2 are stable; the interval of stability extends from the orbits of infinitely
  small dimensions to the orbit of ejection with the heavy body M1 and farther, being limited only by the
  double-periodic critical orbit h26" (p.366).
- The value is 1 - 0.9523..., the maximum of branch h2 (p.359).
- The stability is planar only. The paper says vertical instability can occur (comment (a), p.366).
- In the printed discussion (p.374), Jefferys reports that for the Earth-Moon case he finds the retrograde
  orbit unstable "just before C = 2.9". Henon suggests the stable set may be very small and so missed.

## 1. Method (READ pp.350-354)

- Symmetric periodic orbits start perpendicular on the x-axis.
- A critical orbit satisfies two conditions: periodicity, and stability index a = +-1. With three unknowns
  (mu, x0, C), the critical orbits form one-parameter families as mu varies.
- The stability index and the orbit "type" are computed from a half orbit (Deprit & Price 1965). Types are
  conserved along a critical family.
- Starting points (p.350):
  - mu = 1/2 (Henon 1965b) and mu = 0 (Henon 1969).
  - mu = 1/11 (Shearing 1960) and mu = 0.01215 (Broucke 1968), by interpolation.
- The Earth-Moon value used is mu = 1/82.30 (IAU 1966) (p.367).
- Coordinates: origin at the barycentre. C = x^2 + y^2 + 2(1 - mu)/r1 + 2 mu/r2 - xdot^2 - ydot^2
  (Eq. (25), p.367).
- x0 and x1 are the two x-axis crossings, with x0 taken with ydot0 > 0 (p.367).
- Family definitions (Stromgren notation, p.350):
  - f: retrograde about the second body.
  - g: direct about the second body.
  - h: retrograde about the first body.
  - i: direct about the first body.
  - l: direct (fixed axes) about both bodies.
  - m: retrograde about both bodies.

## 2. Stability results (READ pp.355-366)

- **Family m (p.355):**
  - m1 extends over the whole range 0 <= mu <= 1.
  - m2-m3 is one family that runs from mu = 0 up to mu0 = 0.327... and back.
  - For mu0 < mu < 1 - mu0 there is only m1: "periodic orbits larger than m1 are stable, periodic orbits
    smaller than m1 are unstable".
  - For mu < mu0: stable until m3, unstable between m3 and m2, stable between m2 and m1, unstable after m1.
- **Family l (pp.356-357):** stable until l1, unstable between l1 and l'1, stable between l'1 and l2, mostly
  unstable after l2.
- **Family h (pp.358-360):** a detailed stability sequence over mu intervals, with these breakpoints: 0.783,
  0.844, 0.892, 0.909, 0.921 and 0.9523.
  - h4 and h5 pass through ejection orbits at mu = 0.094 and 0.102.
  - h2 has a maximum at mu = 0.9523 and becomes h'2.
- **Family i (pp.360-361):** i1, i2 and i14 are three branches of one family. i16 is new. i2, i3 and i16 were
  not followed to their ends.
- **mu -> 0 limit (pp.362-365):** the critical orbits terminate at consecutive-collision (second-species)
  orbits where C is extremal. Example, Eq. (17), p.364: family A0 at tau/pi = 0.21345..., x0 = -0.16294...,
  C = -0.39913... is the common end of m1 and m2. These are the same numbers as A0(-1) in Hitzl & Henon
  1977a and in the m1/m2-m3 rows below. This is a cross-source control.

## 3. Tables (READ pp.367-373; transcribed from the text layer, checked against page images)

Columns: mu, x0, x1, C. Type and a are as printed.

**Family m1** (type 1, a = +1; p.367), 12 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | -0.16294 | 1.00000 | -0.39913 |
| 0.01215 | -0.17509 | 0.99085 | -0.41292 |
| 0.05 | -0.21283 | 0.96198 | -0.45065 |
| 0.1 | -0.26219 | 0.92357 | -0.49167 |
| 0.15 | -0.31073 | 0.88519 | -0.52485 |
| 0.2 | -0.35824 | 0.84699 | -0.55158 |
| 0.25 | -0.40457 | 0.80896 | -0.57280 |
| 0.3 | -0.44961 | 0.77103 | -0.58923 |
| 0.35 | -0.49331 | 0.73307 | -0.60142 |
| 0.4 | -0.53571 | 0.69489 | -0.60983 |
| 0.45 | -0.57689 | 0.65628 | -0.61475 |
| 0.5 | -0.61702 | 0.61702 | -0.61636 |

**Family m2-m3** (type 5, a = -1; p.368), 22 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | -0.16294 | 1.00000 | -0.39913 |
| 0.01215 | -0.17814 | 0.99090 | -0.41304 |
| 0.05 | -0.22581 | 0.96278 | -0.45260 |
| 0.1 | -0.29019 | 0.92703 | -0.49976 |
| 0.15 | -0.35717 | 0.89397 | -0.54476 |
| 0.2 | -0.42844 | 0.86531 | -0.59216 |
| 0.25 | -0.50764 | 0.84441 | -0.64997 |
| 0.3 | -0.60723 | 0.84103 | -0.74080 |
| 0.32 | -0.66945 | 0.85573 | -0.81271 |
| 0.32 | -0.77688 | 0.91698 | -0.96696 |
| 0.3 | -0.82094 | 0.95526 | -1.03878 |
| 0.25 | -0.87466 | 1.01047 | -1.12755 |
| 0.2 | -0.90748 | 1.04636 | -1.17672 |
| 0.15 | -0.93207 | 1.07098 | -1.20412 |
| 0.1 | -0.95239 | 1.08506 | -1.21122 |
| 0.05 | -0.97095 | 1.08392 | -1.18846 |
| 0.02 | -0.98313 | 1.06644 | -1.14125 |
| 0.01215 | -0.98718 | 1.05588 | -1.11688 |
| 0.01 | -0.98838 | 1.05189 | -1.10801 |
| 0.005 | -0.99187 | 1.03923 | -1.08059 |
| 0.001 | -0.99637 | 1.01914 | -1.03871 |
| 0.0 | -1.00000 | 1.00000 | -1.00000 |

**Family l1** (type 5, a = -1; p.368), 11 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | -2.08008 | 2.08008 | 3.36525 |
| 0.05 | -2.09831 | 2.11319 | 3.38674 |
| 0.1 | -2.10987 | 2.13214 | 3.40126 |
| 0.15 | -2.11783 | 2.14323 | 3.41160 |
| 0.2 | -2.12357 | 2.14928 | 3.41903 |
| 0.25 | -2.12775 | 2.15178 | 3.42427 |
| 0.3 | -2.13072 | 2.15159 | 3.42773 |
| 0.35 | -2.13261 | 2.14924 | 3.42968 |
| 0.4 | -2.13344 | 2.14503 | 3.43027 |
| 0.45 | -2.13315 | 2.13912 | 3.42957 |
| 0.5 | -2.13186 | 2.13186 | 3.42772 |

**Family l'1** (type 6, a = -1; p.368), 11 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | -2.08008 | 2.08008 | 3.36525 |
| 0.05 | -2.04966 | 2.06755 | 3.36492 |
| 0.1 | -2.04230 | 2.07052 | 3.37109 |
| 0.15 | -2.04855 | 2.08047 | 3.38040 |
| 0.2 | -2.06056 | 2.09193 | 3.39028 |
| 0.25 | -2.07436 | 2.10264 | 3.39953 |
| 0.3 | -2.08811 | 2.11181 | 3.40769 |
| 0.35 | -2.10101 | 2.11925 | 3.41461 |
| 0.4 | -2.11272 | 2.12502 | 3.42025 |
| 0.45 | -2.12310 | 2.12925 | 3.42465 |
| 0.5 | -2.13186 | 2.13186 | 3.42772 |

**Family l2** (type 4, a = +1; p.369), 14 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | -1.52943 | 1.64537 | 3.14812 |
| 0.0001 | -1.52982 | 1.64573 | 3.14827 |
| 0.001 | -1.53269 | 1.64937 | 3.14958 |
| 0.01215 | -1.57029 | 1.68205 | 3.16435 |
| 0.05 | -1.65322 | 1.75430 | 3.19898 |
| 0.1 | -1.72191 | 1.81058 | 3.22886 |
| 0.15 | -1.77013 | 1.84706 | 3.25007 |
| 0.2 | -1.80673 | 1.87226 | 3.26601 |
| 0.25 | -1.83547 | 1.88982 | 3.27819 |
| 0.3 | -1.85832 | 1.90165 | 3.28744 |
| 0.35 | -1.87643 | 1.90884 | 3.29426 |
| 0.4 | -1.89049 | 1.91205 | 3.29894 |
| 0.45 | -1.90092 | 1.91169 | 3.30169 |
| 0.5 | -1.90796 | 1.90796 | 3.30259 |

**Family h1-h5** (type 6, a = -1; p.369), 31 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | -1.00000 | 1.00000 | -1.00000 |
| 0.0001 | -0.96674 | 0.96517 | -0.92885 |
| 0.001 | -0.93266 | 0.92481 | -0.84287 |
| 0.01215 | -0.87126 | 0.82223 | -0.61216 |
| 0.05 | -0.84380 | 0.69530 | -0.31908 |
| 0.1 | -0.85145 | 0.58877 | -0.07657 |
| 0.2 | -0.89737 | 0.42245 | 0.29428 |
| 0.3 | -0.95683 | 0.27803 | 0.61444 |
| 0.4 | -1.02163 | 0.14284 | 0.91770 |
| 0.5 | -1.08874 | 0.01232 | 1.21689 |
| 0.6 | -1.15644 | -0.11592 | 1.51875 |
| 0.7 | -1.22346 | -0.24343 | 1.82756 |
| 0.8 | -1.28967 | -0.37051 | 2.14592 |
| 0.85 | -1.32533 | -0.43153 | 2.30762 |
| 0.88 | -1.35360 | -0.46231 | 2.40178 |
| 0.90319 | -1.40000 | -0.46588 | 2.45970 |
| 0.90749 | -1.50000 | -0.39328 | 2.41200 |
| 0.89636 | -1.60000 | -0.29722 | 2.30836 |
| 0.88 | -1.68993 | -0.20533 | 2.19927 |
| 0.85 | -1.80314 | -0.08547 | 2.05293 |
| 0.8 | -1.92776 | 0.05124 | 1.90069 |
| 0.7 | -2.08333 | 0.23225 | 1.78721 |
| 0.6 | -2.18664 | 0.36700 | 1.81697 |
| 0.5 | -2.26086 | 0.48378 | 1.93355 |
| 0.4 | -2.31049 | 0.59226 | 2.10902 |
| 0.3 | -2.33348 | 0.69662 | 2.32413 |
| 0.2 | -2.32368 | 0.79889 | 2.55845 |
| 0.15 | -2.30516 | 0.84964 | 2.67192 |
| 0.05 | -2.25483 | 0.94857 | 2.84539 |
| 0.02 | -2.23731 | 0.97384 | 2.88470 |
| 0.0 | -2.17445 | 1.00000 | 2.97099 |

**Family h2** (type 5, a = -1; p.370), 27 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | -1.00000 | 1.00000 | -1.00000 |
| 0.001 | -1.00374 | 0.97840 | -0.95625 |
| 0.01215 | -1.01365 | 0.91713 | -0.82807 |
| 0.05 | -1.03082 | 0.81033 | -0.59694 |
| 0.1 | -1.04925 | 0.70651 | -0.36756 |
| 0.2 | -1.08574 | 0.53306 | 0.02221 |
| 0.3 | -1.12481 | 0.37774 | 0.37784 |
| 0.4 | -1.16661 | 0.23088 | 0.72057 |
| 0.5 | -1.21036 | 0.08848 | 1.05941 |
| 0.6 | -1.25487 | -0.05191 | 1.39965 |
| 0.7 | -1.29844 | -0.19237 | 1.74521 |
| 0.8 | -1.33848 | -0.33528 | 2.09970 |
| 0.9 | -1.37422 | -0.48081 | 2.46527 |
| 0.93 | -1.39435 | -0.51553 | 2.57073 |
| 0.95 | -1.47163 | -0.48045 | 2.59580 |
| 0.95212 | -1.60000 | -0.36550 | 2.48910 |
| 0.9472 | -1.80000 | -0.17184 | 2.19396 |
| 0.93272 | -2.00000 | 0.02108 | 1.61888 |
| 0.92 | -2.06283 | 0.07641 | 1.23426 |
| 0.9 | -2.08665 | 0.08378 | 0.75754 |
| 0.875 | -2.06050 | 0.02318 | 0.27853 |
| 0.84471 | -1.96472 | -0.20000 | -0.37650 |
| 0.85105 | -1.90643 | -0.40000 | -0.72124 |
| 0.9 | -1.86091 | -0.68863 | -1.10441 |
| 0.95 | -1.85793 | -0.86497 | -1.30251 |
| 0.98 | -1.86470 | -0.94989 | -1.39014 |
| 1.0 | -1.87212 | -1.00000 | -1.43948 |

**Family h3-h4** (type 1, a = +1; p.370), 24 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | -1.00000 | 1.00000 | -1.00000 |
| 0.005 | -1.11941 | 0.97765 | -0.93879 |
| 0.01215 | -1.15359 | 0.95711 | -0.88523 |
| 0.05 | -1.22162 | 0.87475 | -0.67931 |
| 0.1 | -1.26165 | 0.78461 | -0.46076 |
| 0.2 | -1.31176 | 0.62506 | -0.07851 |
| 0.3 | -1.35353 | 0.47843 | 0.27486 |
| 0.4 | -1.39560 | 0.33922 | 0.61697 |
| 0.5 | -1.44099 | 0.20565 | 0.95555 |
| 0.6 | -1.49265 | 0.07815 | 1.29475 |
| 0.7 | -1.56137 | -0.03655 | 1.63566 |
| 0.75 | -1.61930 | -0.07745 | 1.80427 |
| 0.78 | -1.70624 | -0.06795 | 1.90029 |
| 0.77947 | -1.80000 | -0.01082 | 1.90007 |
| 0.75713 | -1.90000 | 0.07343 | 1.85733 |
| 0.7 | -2.02811 | 0.20064 | 1.80366 |
| 0.6 | -2.16040 | 0.35425 | 1.82467 |
| 0.5 | -2.24633 | 0.47790 | 1.93760 |
| 0.4 | -2.30191 | 0.58944 | 2.11125 |
| 0.3 | -2.32823 | 0.69529 | 2.32540 |
| 0.2 | -2.32061 | 0.79836 | 2.55910 |
| 0.1 | -2.27807 | 0.89999 | 2.77159 |
| 0.05 | -2.25385 | 0.94882 | 2.84550 |
| 0.0 | -2.17445 | 1.00000 | 2.97099 |

**Family h26-h29** (type 6, a = -1; p.371), 19 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 1.0 | -1.00000 | -1.00000 | -1.00000 |
| 0.99975 | -1.20942 | -0.99438 | -1.17920 |
| 0.999 | -1.28384 | -0.98674 | -1.22556 |
| 0.98785 | -1.51254 | -0.93470 | -1.30099 |
| 0.96 | -1.70389 | -0.85561 | -1.27869 |
| 0.9 | -1.92156 | -0.72878 | -1.14731 |
| 0.8 | -2.15084 | -0.55670 | -0.86747 |
| 0.7 | -2.31884 | -0.40511 | -0.54758 |
| 0.6 | -2.45343 | -0.26484 | -0.19645 |
| 0.5 | -2.56396 | -0.13280 | 0.18277 |
| 0.4 | -2.65476 | -0.00888 | 0.58826 |
| 0.3 | -2.73080 | 0.10175 | 1.01723 |
| 0.25 | -2.76901 | 0.14487 | 1.23727 |
| 0.2 | -2.82741 | 0.15822 | 1.45091 |
| 0.18007 | -2.90000 | 0.11346 | 1.51810 |
| 0.225 | -3.11721 | -0.11037 | 1.37454 |
| 0.3 | -3.23311 | -0.25005 | 1.31770 |
| 0.4 | -3.32840 | -0.38232 | 1.38068 |
| 0.5 | -3.39134 | -0.49454 | 1.53549 |

**Family h27-h28** (type 1, a = +1; p.371), 17 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 1.0 | -1.87212 | -1.00000 | -1.43948 |
| 0.98785 | -1.89689 | -0.97133 | -1.41094 |
| 0.95 | -1.97660 | -0.88654 | -1.31600 |
| 0.9 | -2.07241 | -0.78484 | -1.18090 |
| 0.8 | -2.23829 | -0.60572 | -0.88476 |
| 0.7 | -2.37838 | -0.44701 | -0.55823 |
| 0.6 | -2.49813 | -0.30147 | -0.20370 |
| 0.5 | -2.60039 | -0.16609 | 0.17743 |
| 0.4 | -2.68751 | -0.04111 | 0.58392 |
| 0.3 | -2.76597 | 0.06596 | 1.01303 |
| 0.25 | -2.81226 | 0.10134 | 1.23237 |
| 0.21 | -2.88429 | 0.09072 | 1.40242 |
| 0.20544 | -3.00000 | -0.00190 | 1.42291 |
| 0.25 | -3.13882 | -0.14787 | 1.34421 |
| 0.3 | -3.21866 | -0.24166 | 1.31950 |
| 0.4 | -3.32150 | -0.37938 | 1.38149 |
| 0.5 | -3.38747 | -0.49343 | 1.53589 |

**Family i14-i1-i2** (type 1, a = +1; pp.371-372), 25 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 1.0 | -1.00000 | -1.00000 | 3.00000 |
| 0.999 | -0.97400 | -1.03180 | 3.04043 |
| 0.98785 | -0.93059 | -1.06654 | 3.18267 |
| 0.95 | -0.85591 | -1.08116 | 3.38219 |
| 0.9 | -0.77734 | -1.07028 | 3.50964 |
| 0.8 | -0.63440 | -1.02348 | 3.61469 |
| 0.7 | -0.49590 | -0.96208 | 3.63127 |
| 0.6 | -0.35602 | -0.88944 | 3.60595 |
| 0.5 | -0.20993 | -0.80096 | 3.56641 |
| 0.45 | -0.13142 | -0.74454 | 3.55695 |
| 0.425 | -0.08721 | -0.70401 | 3.57685 |
| 0.425 | -0.08234 | -0.69044 | 3.60868 |
| 0.45 | -0.11344 | -0.69419 | 3.67170 |
| 0.5 | -0.18136 | -0.72145 | 3.73896 |
| 0.6 | -0.31723 | -0.78754 | 3.80842 |
| 0.7 | -0.45198 | -0.85784 | 3.81637 |
| 0.8 | -0.58757 | -0.92775 | 3.75519 |
| 0.85 | -0.65631 | -0.96038 | 3.68884 |
| 0.88 | -0.69646 | -0.97673 | 3.63254 |
| 0.88 | -0.67780 | -0.95994 | 3.63281 |
| 0.85 | -0.61146 | -0.91996 | 3.69155 |
| 0.8 | -0.50372 | -0.85143 | 3.76668 |
| 0.7 | -0.34389 | -0.74693 | 3.85662 |
| 0.6 | -0.20143 | -0.65049 | 3.88505 |
| 0.5 | -0.06445 | -0.55686 | 3.85893 |

**Family i3** (type 6, a = -1; p.372), 9 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 1.0 | -1.00000 | -1.00000 | 3.00000 |
| 0.999 | -0.94110 | -1.00576 | 3.03933 |
| 0.98785 | -0.85721 | -1.00160 | 3.18451 |
| 0.96 | -0.76868 | -0.97871 | 3.36304 |
| 0.9 | -0.64393 | -0.92367 | 3.57568 |
| 0.8 | -0.47982 | -0.82911 | 3.76293 |
| 0.7 | -0.33334 | -0.73416 | 3.85381 |
| 0.6 | -0.19435 | -0.63989 | 3.88267 |
| 0.5 | -0.05924 | -0.54713 | 3.85685 |

**Family i15** (type 5, a = -1; p.372), 17 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | 0.48075 | -0.48075 | 3.46681 |
| 0.01215 | 0.46955 | -0.49175 | 3.44884 |
| 0.05 | 0.43421 | -0.53001 | 3.36348 |
| 0.1 | 0.33385 | -0.57652 | 3.24709 |
| 0.14 | 0.24566 | -0.61819 | 3.22186 |
| 0.2 | 0.13094 | -0.68907 | 3.23020 |
| 0.3 | -0.03448 | -0.80347 | 3.29238 |
| 0.4 | -0.18401 | -0.90259 | 3.37064 |
| 0.5 | -0.32456 | -0.98734 | 3.44585 |
| 0.6 | -0.45913 | -1.05983 | 3.50790 |
| 0.7 | -0.58956 | -1.12048 | 3.54626 |
| 0.8 | -0.71746 | -1.16676 | 3.54372 |
| 0.9 | -0.84529 | -1.18821 | 3.46028 |
| 0.95 | -0.91138 | -1.17640 | 3.34980 |
| 0.98785 | -0.96694 | -1.12685 | 3.17001 |
| 0.998 | -0.98753 | -1.07317 | 3.05864 |
| 1.0 | -1.00000 | -1.00000 | 3.00000 |

**Family i16** (type 6, a = -1; p.373), 10 rows

| mu | x0 | x1 | C |
|---|---|---|---|
| 0.0 | 0.48075 | -0.48075 | 3.46681 |
| 0.01215 | 0.46095 | -0.48325 | 3.47484 |
| 0.05 | 0.40409 | -0.49738 | 3.48717 |
| 0.1 | 0.33464 | -0.52457 | 3.48362 |
| 0.15 | 0.26724 | -0.55901 | 3.46028 |
| 0.2 | 0.19477 | -0.60715 | 3.40808 |
| 0.23 | 0.13726 | -0.65804 | 3.35350 |
| 0.245 | 0.09334 | -0.70342 | 3.31118 |
| 0.24674 | 0.08000 | -0.71832 | 3.29403 |
| 0.24249 | 0.07000 | -0.73046 | 3.26738 |

Notes on the tables:
- Families h26-h29, h27-h28, i14-i1-i2, i3 and i16 end with "..." in print: they were not followed further.
- Label spellings follow the page images. The text layer misreads them as "II", "1'1", "12", "ha6-ha9",
  "ha7-ha8" and "it5".
- In family m2-m3 the printed "0.01" row comes between 0.01215 and 0.005.

## 4. Use as controls

- **`#949` X4 (Pluto-Charon and Titan, DA enumerator):**
  - Pluto-Charon has mu of about 0.109. Interpolate in the tables between 0.1 and 0.15. For example, m1 at
    mu = 0.1 has x0 = -0.26219, x1 = 0.92357, C = -0.49167.
  - The critical orbits give stability boundaries that any computed family member must reproduce.
  - The 0.0477 retrograde-stability result covers Titan (2.4e-4) and every planet-moon pair.
- **`#956` R9 (exterior-realm, L2):**
  - The l and m families are the large "around both bodies" orbits. The l1, l'1 and l2 values at
    mu = 0.01215 or 0.1 bound the stable windows.
  - For Earth-Moon, the l2 row is mu = 0.01215, x0 = -1.57029, x1 = 1.68205, C = 3.16435.
- **Any stability claim:** a claimed stable or unstable symmetric orbit of these families at a given mu
  must sit on the correct side of the tabulated critical orbits. Planar only.
- **Cross-source:** the mu = 0 end points match Hitzl & Henon 1977a's A0(-1). The m2-m3 end at mu = 0 is
  the retrograde circular orbit of radius 1 (x0 = -1, C = -1).

## 5. Citation mining (references p.373)

Held:
- Broucke 1968, TR 32-1168.
- Henon 1968, Bull. Astron. 3:377.

Not held, in priority order. These are journal articles of 1965-1970; DOIs were not checked (ADS
bibcodes exist).
1. Henon, M. (1965b), "Exploration numerique du probleme restreint II. Masses egales, stabilite des orbites
   periodiques", Ann. Astrophys. 28:992-1007. The stability index and orbit types. Also wanted by the Hitzl
   digests.
2. Henon, M. (1969), "Numerical exploration of the restricted problem. V. Hill's case: periodic orbits and
   their stability", Astron. Astrophys. 1:223. The mu = 0 Hill-case critical orbits.
3. Henon, M. (1965a), Ann. Astrophys. 28:499 (part I); Henon 1966a, b, Bull. Astron. (3) 1 (parts III, IV);
   Henon 1970, A&A (vertical stability, then in press).
4. Shearing, G. (1960), thesis, Manchester (mu = 1/11 families). No DOI.
5. Deprit, A. & Price, J. F. (1965), AJ 70:836 (half-orbit stability); Deprit & Henrard (1969), AJ 74:308;
   Colombo, Franklin & Munford (1968), AJ 73:111; Guillaume (1969), A&A 3:57; Message (1970), same volume
   p.19; Bulirsch & Stoer (1966); Stromgren (1935), Copenhagen Publ. 100; Jackson (1913), MNRAS 74:62;
   Arnold 1963; Moser 1962.

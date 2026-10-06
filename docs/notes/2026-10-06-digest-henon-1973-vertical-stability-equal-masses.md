# Digest: Hénon 1973, "Vertical Stability of Periodic Orbits in the Restricted Problem. I. Equal Masses" (#960 batch 27)

M. Hénon (Observatoire de Nice), Astronomy & Astrophysics 28:415-426 (1973), ADS 1973A&A....28..415H.
Received 20 July 1973.
- Filed as `cyclers_pdf/papers/henon-1973-vertical-stability-periodic-orbits-restricted-problem-I-equal-masses-aa-28-415-ads-1973AA-28-415H.pdf`.
  - This is an `ocrmypdf --force-ocr` copy of the image-only ADS scan. The source scan (fetched by
    fetch-sonnet) is md5 c788dfa42a7c3d2870e8acd4bdea8fa2, 12 pp.
- **Companion `...-tables.txt`: Tables 1-2 transcribed from 400 dpi images.**
  - Three witnesses (tesseract with a whitelist, Vision, and the image, which decides).
  - Table 1 has 141 rows and 412 numeric cells. Tesseract differed in 67 (55 of them lost minus signs)
    and Vision in 2; all were settled on the image.
  - Checks: consecutive numbering in every family; c/l/m symmetry (x1 = -x0); ejection rows at +-0.5;
    L_1 and L_3 C values recomputed; the Table 2 symmetry required by eq. (45).
- Wanted list: removed in batch 27.

## 0. Verdict

**The best published numeric substitute for the absent Bruno & Varin Table 8 (family h, mu = 1/2): six
digits, against four in Hénon 1965b.** Hénon's x0 and x1 are the crossings, and his C_H = C - 0.25 in
the Bruno-Varin convention.

Family h rows of Table 1 (Name, x0, x1, C, y0dot, y1dot, Type):

| Orbit | x0 | x1 | C | Type |
|---|---|---|---|---|
| h1 | -1.088744 | 0.012321 | 1.216891 | D |
| h2 | -1.210359 | 0.088477 | 1.059408 | A |
| h3 | -1.440989 | 0.205646 | 0.955552 | C |
| h1v | -1.666728 | 0.292111 | 1.042853 | A_v |
| h2v | -1.919377 | 0.367855 | 1.344878 | A_v |
| h4 | -2.246340 | 0.477902 | 1.937595 | C |
| h5 | -2.260877 | 0.483782 | 1.933543 | D |
| h_2 6 | -2.563959 | -0.132801 | 0.182771 | D |
| h_2 7 | -2.600391 | -0.166092 | 0.177425 | C |
| h_2 3v | -2.808075 | -0.291588 | 0.297945 | A_v |
| h_2 4v | -3.213599 | -0.428523 | 1.091642 | A_v |
| h_2 8 | -3.387546 | -0.493456 | 1.535892 | C |
| h_2 9 | -3.391440 | -0.494564 | 1.535464 | D |

- The horizontal critical orbits h1-h_2 9 agree with Hénon 1965b's four-digit values. For example, h1:
  -1.0875 / 1.2175 there, against -1.088744 / 1.216891 here.
- The vertical critical orbits h1v, h2v, h_2 3v and h_2 4v are Bruno-Varin k = 4, 5, 11 and 12
  (KIAM 64/2005 sec. 7).
- **Vertical stability:** a_v = vertical index, stable when |a_v| < 1.
  - The near-circular orbits around one or both primaries (families f, g, h, i, l, m) are 3-D stable in
    extended regions. Instability generally appears first in the plane.
  - Exception: family m has a large region that is horizontally stable but vertically unstable.
  - Each vertical critical orbit is an intersection with a family of 3-D periodic orbits.
- **Use:** `#944`, `#956` and any spatial-stability claim at mu = 1/2. The 52 vertical critical orbits
  are the bifurcation points to 3-D families (Bray & Goudas-type).

## 1. Notes

- There is a short Celest. Mech. 8:269 note with the same title (doi 10.1007/BF01231427). It is not held
  and is not needed now.
- Oddities recorded by the transcriber:
  - "i18v" is printed without its subscript 2;
  - one stray mark on i2v;
  - i_3 26 is also an orbit of family h described three times (text p.10).

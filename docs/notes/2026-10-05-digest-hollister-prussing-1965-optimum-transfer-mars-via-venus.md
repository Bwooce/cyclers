# Digest: Hollister & Prussing 1965, "Optimum Transfer to Mars via Venus" (#960)

W. M. Hollister & J. E. Prussing (MIT Experimental Astronomy Laboratory), "Optimum Transfer to Mars via
Venus", AIAA Paper 65-700, AIAA/ION Astrodynamics Specialist Conference, Monterey, 16-17 September 1965,
doi 10.2514/6.1965-700 (wanted-list CONFIRMED).
- Filed as `cyclers_pdf/papers/hollister-prussing-1965-optimum-transfer-mars-via-venus-aiaa-65-700-doi-10.2514-6.1965-700.pdf`.
  15 pages, text layer (two columns interleaved in extraction), md5 98ac068248b5f891d320e750218703de.
- The journal version, Astronautica Acta 12(2):169-179 (1966), is not held. It stays on the wanted list
  for attribution only: nothing here suggests the content differs (unverified).
- I read all pages from the text layer, and the figure pages (Figs. 1-9) as thumbnails. There are no
  tables.

## 0. Verdict

This is a one-shot mission study: Earth-Mars round trips with ONE Venus flyby (outbound or homebound),
optionally with a thrust impulse at Venus, for dates 1970-1990. Nothing is periodic. **No literal
collision for `#942` R1(c).**
- **Planet constants (the #942 question):** none are stated. The paper gives no ephemeris source, planet
  elements or periods. It works in EMOS units (Earth mean orbital speed) and normalised flyby-planet units
  (lengths in planet radii, velocities in circular surface speed). For the Hollister group's elements, use
  the Rall thesis listing (1960 mean elements).

## 1. Content (READ)

- **Availability rule (pp.2-3, Fig. 1):**
  - A low-energy Venus-assisted transfer to Mars needs the Sun-Venus-Mars alignment about halfway
    through the trip, and on the correct side relative to Earth.
  - Venus "is available on one leg during every Mars opposition period"; during "every third Mars
    opposition period it is possible to utilize Venus both going and returning".
- **Thrusted flyby (pp.4-16):**
  - Given inbound and outbound v_inf vectors, find the smallest impulse that connects the two hyperbolas
    without passing below the surface. This is a planar (2-D) optimisation.
  - The common-peripoint impulse (the scalar difference of the two periapsis speeds) is within 3 % of
    the optimum in the typical case (p.13, Fig. 4).
  - When the common peripoint is below the surface, approximate the constrained optimum along the less
    energetic boundary hyperbola (p.16).
  - Computed with a MAD program on the IBM 7094.
- **Results (pp.17-20, Figs. 5-9):**
  - Contours of 0.2 EMOS launch speed for direct and via-Venus trips, 1970-1990.
  - Venus is available both ways in 1971, outbound only in 1973, homebound only in 1975, and both again
    in 1978 (p.17).
  - Whenever a pure Venus flyby would go below the surface, a neighbouring direct trip is cheaper.
  - Where a pure flyby is possible, thrust at Venus saves only "a few hundred feet per second".

## 2. Gate relevance

- **`#942` R1(c):** no collision. It is background on V-M geometry (the Sun-Venus-Mars alignment
  criterion) by the R1(b) author.
- **#942 H-M control:** no planet constants. Nothing to add beyond Rall's thesis listing.

## 3. Positive controls

- None numeric. The results are contour plots only.

## 4. Citation mining (references 1-15)

Held: Sohn 1964 [6] (filed today). Not held: Battin 1964 [4]; Ross 1963 [10] (in the wanted list).

Not held:
1. [5] Hollister, W. M. (1963), "The Mission for a Manned Expedition to Mars", Sc.D. thesis, MIT. Already
   in the wanted list.
2. [8] Hollister, W. M. (1964), "Mars Transfer via Venus", AIAA/ION Astrodynamics Guidance and Control
   Conference, AIAA 64-647. Not in the wanted list; DOI not checked. Low priority (predecessor of this
   paper).
3. [7] Deerwester (1965), AIAA 65-89; [10] Ross (1963); [9] Gobetz (1963), AIAA J 1(9); [12] Lee (1964);
   [13] Luidens & Kappraff (1965), NASA TN D-2605; [14] Titus (1965); [15] Fimple (1962). Background.
4. [1] Crocco 1956; [2] Lockheed 1962; [3] Ford Aeronutronic EMPIRE 1962. Background (Crocco is in the
   wanted list).

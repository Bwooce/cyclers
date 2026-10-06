# Digest: Lantukh & Russell 2012, "Automated Inclusion of n-pi Transfers in Gravity-Assist Flyby Tour Design" (#960, #943)

D. V. Lantukh & R. P. Russell (UT Austin), "Automated Inclusion of n-pi Transfers in Gravity-Assist Flyby
Tour Design", AIAA 2012-4749, AIAA/AAS Astrodynamics Specialist Conference, Minneapolis, August 2012,
doi 10.2514/6.2012-4749.
- Filed as `cyclers_pdf/papers/lantukh-russell-2012-automated-inclusion-n-pi-transfers-gravity-assist-flyby-tour-design-aiaa-2012-4749-doi-10.2514-6.2012-4749.pdf`.
  20 pages, text layer, md5 123cd1e4a14f3f0c725353ba9e5d5c38. Supplied by the owner.
- I read the abstract, sec. I and sec. IV (examples) from the text layer, and searched the whole text
  for moon and cycler content.

## 0. Verdict

**No `#943` gc collision.**
- A pathfinding method that strings same-body n-pi transfers (resonant even-n pi and half-rev odd-n pi)
  into sequences on the v_inf sphere. It works as a "gap" in a Lambert-based tour search.
- All examples (A-J, sec. IV, Tables 6-7) are in normalised units around one flyby body.
- The paper never names Ganymede, Callisto or Europa. "Jupiter" appears only in reference titles. There
  is no cycler and no two-moon sequence.
- My prior-art note had suspected "Jovian n-pi sequences" (`2026-10-06-943-gc-prior-art-search.md`,
  row 8). That suspicion is cleared.

## 1. Content (READ)

- n-pi transfers have extra degrees of freedom over Lambert. The v_inf-sphere elements are:
  - even-n pi: circles (resonant M:N);
  - odd-n pi: points.
- Sequences are enumerated with transfer-angle sets (Table 3; e.g. 2:3 then 1:1). Degrees of freedom are
  chosen by optimisation. Algorithm outlines are in Tables 1-2.
- Examples A-J vary the inbound angle (pi/6, 2pi/3, pi/3), the number of sequences and the turn
  constraints. Run times are on a single core.

## 2. Gate relevance

- X1: tooling background only (same-body n-pi sequencing, related to Russell & Ocampo 2005, held).

## 3. Positive controls

- None (normalised examples).

## 4. Citation mining

The references overlap the held corpus: Russell & Ocampo 2005, Russell & Strange 2009, McConaghy et al.
2005, Wolf & Smith 1995 (held); Uphoff, Roberts & Friedman 1974 (AIAA 74-781; the 1976 JSR version is in
the wanted list). Nothing new for X1.

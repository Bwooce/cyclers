# Digest — Rauwolf, Friedlander & Nock 2002, "A Mars Cycler Architecture Utilizing Low-Thrust Propulsion" (AIAA 2002-5046)

**Acquired 2026-06-22** (the lone #416 acquisition gap; user-supplied PDF). AIAA/AAS
Astrodynamics Specialist Conf., Monterey CA, 5–8 Aug 2002. SAIC (Rauwolf,
Friedlander) + Global Aerospace (Nock). NIAC-funded. **PDF needs filing into
`cyclers_pdf/papers/` + CORPUS_INDEX registration.**

Purpose here: the missing independent primary for the #415 ΔV-band Axis-B
(real-ephemeris cycler maintenance ΔV). All numbers verbatim with table/figure
citations.

## Cycler MAINTENANCE ΔV (the load-bearing data) — impulsive, conic-optimized, 200 km min flyby

- **Aldrin Up Escalator (outbound), Table 1 / Fig 3** — 3 mid-course ΔV over the
  15-yr (≈7-cycle) sequence (2011–2026 opportunities):
  ΔV1 = **0.211 km/s** (18-Apr-2017), ΔV2 = **0.679 km/s** (20-Jun-2019),
  ΔV3 = **0.671 km/s** (18-Aug-2021). **Total = 1.561 km/s / 15 yr.**
  Quote: "the Aldrin Cycler requires a modest orbit correction on **3 out of 7
  orbits** to maintain the proper orbit orientation." "a 200-km flyby altitude
  constraint has been imposed." "Mathematically correct, but impractical
  subsurface flybys can eliminate the need for corrections."
- **Aldrin Down Escalator (inbound), Table 2** — ΔV1 = **0.046 km/s**,
  ΔV2 = **0.626 km/s**, ΔV3 = **0.933 km/s**. **Total = 1.605 km/s / 15 yr.**
- BASIS: km/s, total over 15 yr (~7 cycles), 3 discrete impulses; conic
  (2-body + flyby) optimizer (SAIC CHEBY2); per Escalator (one Astrotel).

## SEP execution of the same corrections (band-4 / low-thrust boundary)

- The same 3 corrections flown by SEP (Xenon ion, **Isp 5000 s**, 69% power
  efficiency, 8 kg/kW, 150 kWe baseline): "propellant mass is about **3 metric
  tons over a 15-year cycle** … **propellant mass fraction of 4.0%** compared to
  38.2% if chemical propulsion was used" (Table 3/4 text). SEP burn ~430 days of
  the orbit (Fig 4), centered on apoapsis.
- This confirms the "~1500 m/s/15-yr SEP maintenance" figure (previously via Howe
  2025's secondhand citation) — it is the **same 1.561 km/s**, executed by SEP,
  not an extra cost.

## Encounter V∞ (corroborates high-V∞ Aldrin)

- Up Escalator (Table 1): Earth V∞ 5.7–5.9 km/s, Mars V∞ 7.1–11.5 km/s.
  Down (Table 2): Earth 5.5–6.2, Mars 7.2–11.9 km/s. Intro: Aldrin "encounter
  velocity magnitudes are as large as **11.9 km/s**."

## Taxi ΔV (NOT cycler maintenance — crew-transfer / establishment axis)

- Tables 7–8: Up taxi Earth-departure **~2.5 km/s**, Mars-departure **5–10 km/s**
  ("PATs for ΔV over 3.4 km/s; separate engines for ΔV > 6.6 km/s"). Down taxi
  similar (Mars 5.5–9.5 km/s). These are taxi/transfer costs, kept OUT of the
  maintenance band (cf. McConaghy 2006 taxi DSMs; the catalogue's
  establishment/`v_infinity_leveraging_dv_kms` axis).

## Band impact (folds into `2026-06-22-dv-band-definitions.md`)

**The Aldrin (powered, TR<1) maintenance ΔV is now THREE independent primaries:**
- Byrnes-Longuski-Aldrin 1993: 1.73 (out) / 2.04 (in) km/s / 15 yr.
- **Rauwolf 2002: 1.561 (out) / 1.605 (in) km/s / 15 yr.**
- Russell 2004: 0 m/s at the best 7-cycle launch window, ~3.3 km/s averaged over
  21 windows ("spiky").
All agree at **~1.5–2 km/s / 15 yr (≈7 cycles), 3 of 7 orbits**, with Russell
explaining the best-vs-average spread. This is NOT a conflict with Russell's
< 1/< 10/< 300 m/s tiers — those tiers are for the **near-ballistic/ballistic
class (TR ≥ 1)**, whereas Aldrin is the **powered class (TR = 0.86)**. So the two
bands are now each multi-sourced:
- **powered cycler:** ~1.5–2 km/s/15 yr (Byrnes 1993, Rauwolf 2002, Russell avg).
- **ballistic/near-ballistic:** < 1/< 10/< 300 m/s/7-cycle (Russell 2004),
  magnitude corroborated by McConaghy 2006 S1L1 ~10 m/s/30 yr.

Net: the Axis-B maintenance band is no longer single-sourced — the powered-class
magnitude is firmly triangulated; the near-ballistic class still leans on Russell
for the *cutoffs* but its *magnitude* is independently corroborated.

## References
- Rauwolf, Friedlander, Nock, AIAA 2002-5046 (this paper).
- `docs/notes/2026-06-22-dv-band-definitions.md` (#415 bands),
  `docs/notes/2026-06-17-digest-byrnes-longuski-aldrin-1993.md`,
  `docs/notes/2026-06-07-russell-2004-continuation-deepdive.md`.

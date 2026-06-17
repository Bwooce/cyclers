# Conte-Spencer 2018 — Mission analysis for Earth to Mars-Phobos distant Retrograde Orbits

**Deep-read verdict note, 2026-06-17 AET.** Supersedes the file-only entry in
`2026-06-17-mars-cycler-wave-digest.md`.

## Header

- **Title:** *Mission analysis for Earth to Mars-Phobos distant Retrograde Orbits*
- **Authors:** Davide Conte, David B. Spencer (Penn State, Dept. Aerospace Eng.)
- **Venue:** Acta Astronautica 151 (2018) 761-771
- **DOI:** 10.1016/j.actaastro.2018.06.049
- **Pages:** 11

## What the paper actually is

A trajectory-design paper for **LEO → Mars-Phobos DRO insertion** —
i.e. depart LEO, do an interplanetary Lambert/porkchop transfer to Mars
SOI, then perform a four-burn Concept of Operations (ConOps) sequence to
insert into a periodic Distant Retrograde Orbit around Phobos in the
**Mars-Phobos CR3BP** (mass ratio λ = 1.661 × 10⁻⁸, semimajor axis
9376 km, primary-period 7.657 h). DRO amplitudes A_x = 15 to 300 km
considered; below 15 km collides with Phobos, above 300 km is too far for
proximity ops.

The Mars-Phobos DRO is the same family as the Mars-Phobos CR3BP libration
orbits used in the Wallace-NASA work (their Wallace 2012 cite is reference
[2] of this paper).

Quantitative results:

- **Table 3** (worked example, 2035-06-27 departure / 2036-01-15 arrival,
  A_x = 300 km DRO):
  - Arrival v_∞ = 2.543 km/s, inclination i = 20.93°, phasing angle α = 15°
  - ΔV1 = 0.7481 km/s, ΔV2 = 0.1503, ΔV3 = 0.3479, ΔV4 = 0.2854 km/s
  - ΔV_tot = 1.5316 km/s, Mars-SOI-to-DRO TOF = 1.18 d
- **Table 4** (LEO-to-Mars-Phobos-DRO minimum ΔV per synodic period
  2020-2042, A_x = 30 km, h_LEO = 300 km):
  - Each row gives ΔV_tot (km/s), TOF (days), departure date, arrival date.
  - Minimum 2035-2036 ΔV = 5.2120 km/s, TOF 204 d (dep 27 Jun 2035,
    arr 15 Jan 2036).
  - Best over the 2020-2042 window is the 2030-2031 syn period:
    ΔV_tot = 5.7098 km/s — wait, the bolded minimum is the 2041-2042
    season at 5.1856 km/s, dep 20 Oct 2041 / arr 31 Aug 2042.
- **Table 5** (Lunar DRO A_x = 61500 km → Mars-Phobos DRO A_x = 30 km, same
  syn periods):
  - Best is **2035-2036 at ΔV_tot = 2.7317 km/s**, dep 23 Jun 2035 / arr
    14 Jan 2036, TOF 205 d. Lunar DRO is ~2/5 the ΔV of direct LEO-to-LMO,
    reducing required Earth-launch dry mass.

CR3BP analysis uses JPL DE430 + DE431 ephemerides for the interplanetary
portion; CR3BP for Phobos-vicinity only.

## Catalogue / KNOWN_CORPUS relevance

**This is NOT a cycler paper.** A Mars-Phobos DRO is a *parking orbit* —
the spacecraft arrives, inserts, and stays. There is no repeated Mars
encounter, no Earth-return free encounter, no resonance with the Mars-Earth
synodic period. Cyclic / quasi-cycler classifications do not apply.

**Possible relevance to existing Wallace Mars-Phobos anchor:** the current
KNOWN_CORPUS has a `Wallace Mars-Phobos CR3BP rendezvous trajectory (NASA TM)`
anchor that explicitly extends authorship to Genova/Mars Express and
cites Phobos BCR4BP context (#334). Conte-Spencer 2018 cites Wallace 2012
as their reference [2] (`M.S. Wallace, J.S. Parker, N.J. Strange, D. Grebow,
Orbital operations for Phobos and Deimos exploration, AIAA/AAS 2012-5067`).

**Recommendation: add Conte & Spencer as a secondary citation on the
existing Wallace anchor**, NOT a new anchor. Both papers cover the same
family of Mars-Phobos CR3BP DRO / libration orbits. The combined
authorship list becomes ("Wallace", "Parker", "Strange", "Grebow", "Conte",
"Spencer"). Keywords already cover ("Mars-Phobos CR3BP rendezvous",
"Mars-Phobos libration orbit") — perhaps add "Mars-Phobos distant
retrograde orbit" and "DRO insertion".

### Catalogue impact

No new row. The Mars-Phobos system is not currently in the catalogue
(catalogue scope expansion 2026-06-15 includes `cycler / quasi_cycler /
precursor_mga / mga_tour` for Mars-Phobos-like systems, but a DRO insertion
itinerary doesn't fit any of these four classes). The paper is **reference
for the Wallace anchor**.

If a `precursor_mga` row were to be added for Earth-to-Mars-Phobos staging
(LEO → Mars-Phobos DRO as a precursor to surface-Phobos operations),
Conte-Spencer Table 4/5 would provide the V0 reproducible data (synodic
period, ΔV_tot, TOF, departure/arrival dates). The `precursor_mga` class
admits "inserts_into a cycler" — Mars-Phobos DRO doesn't insert into a
cycler, so this is a stretch and likely not catalogue-worthy.

## Errata vs the pre-read survey

The pre-read survey filed this as "*Tangentially related to cyclers
(Mars-Phobos DRO is a parking orbit, not a cycler). File-only unless #313
Mars-Phobos work reopens.*" The deep-read **confirms** this verdict —
correctly identified as a parking-orbit DRO paper, no cycler relevance.
The only minor correction is that there IS an existing KNOWN_CORPUS Wallace
Mars-Phobos anchor that could be strengthened by adding Conte-Spencer as
a co-citation (the survey did not mention this anchor existed).

## Action items for parent

1. **No new KNOWN_CORPUS anchor.** Update the existing Wallace Mars-Phobos
   anchor to add Conte-Spencer as a co-citation:
   - Add to `authors`: `"Conte", "Spencer"`
   - Add to `keywords`: `"Mars-Phobos distant retrograde orbit"`,
     `"Mars-Phobos DRO"`
   - Update `citation`: append `"; Conte & Spencer (2018), Acta Astro
     151:761-771, DOI 10.1016/j.actaastro.2018.06.049"`
   - The Wallace anchor's `context` already mentions Phobos-Deimos and
     BCR4BP; no change needed.
2. **No new catalogue row.** Mars-Phobos DRO insertion is not a
   cycler / quasi_cycler / precursor_mga / mga_tour.
3. **#313 Mars-Phobos reactivation context (file-only):** if #313 is
   reactivated, Conte-Spencer is a reference for LEO-to-Mars-Phobos-DRO
   ΔV budget and ConOps. Their Table 4/5 windows (2035-2036, 2041-2042)
   are the cheapest of the 2020-2042 era.
4. **Note for the lunar-DRO branch (Kakoi 2014 / #316):** Conte-Spencer's
   Table 5 explicitly studies *Lunar DRO → Mars-Phobos DRO* as a 2.73 km/s
   transfer (2035-2036). This is the same Earth-Moon-to-Mars staging
   concept as #316 / Kakoi 2014, with Phobos DRO as the destination. If
   #316's framework needs an external check, Table 5 is a sourced
   independent data point.

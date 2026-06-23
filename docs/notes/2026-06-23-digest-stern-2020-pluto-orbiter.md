# Digest — Stern, Tapley, Finley & Scherrer 2020, "Pluto Orbiter–Kuiper Belt Explorer: Mission Design for the Gold Standard"

**Journal of Spacecraft and Rockets**, Article in Advance, 2020. DOI 10.2514/1.A34658.
Received 2019-11-19; accepted 2020-03-01. SwRI + Ponderosa Labs (T. Finley) + JHU/APL.
Corpus file (private `cyclers_pdf/papers/`):
`stern-tapley-finley-scherrer-2020-pluto-orbiter-kuiper-belt-explorer-jsr-a34658.pdf`.
Processed 2026-06-23 (#429 follow-on, user-supplied).

## What it is
A flagship-class **mission-design study** (proof of concept, notional 2028 launch) for a
combined **Pluto orbiter + Kuiper-Belt explorer** powered by RTG-sourced electric
propulsion (NEXT-C gridded ion, Isp 4195 s). Three legs: (1) EP+gravity-assist
heliocentric cruise to Pluto; (2) a **Pluto-system orbital tour that maneuvers almost
entirely on Charon gravity assists**; (3) Charon-powered escape back into the KB and a
low-thrust transfer to a second dwarf planet (Quaoar selected from 67 KBO candidates).
Not flown → a study, not a catalogue-admissible flown mission.

## Values relevant to us (flyby-altitude floors / #429)

**Pluto design floor = 100 km (SOURCED).** Section III tour spec 3: "Pluto periapse
altitudes are between 100 and 500 km to perform in situ atmospheric measurements." This
is a *design* minimum flyby/periapse altitude for an atmospheric-dip pass — exactly the
class of value we want for a flyby floor. It **matches our existing code floor**
(`constants.py` dwarf-planet `safe_alt_km` engineering default 100 km), so it upgrades
that constant's provenance **convention → sourced** with NO model change. (Note: the
code comment said "no atmosphere" for the 100 km default; Pluto *does* have a thin
atmosphere — the value is right, the rationale now sourced.)

**Charon = primary GA maneuvering body, but NO numeric periapsis altitude here.** Specs
4/8 use Charon as the primary orbit-maneuvering tool and for system escape; the tour
"targets close approaches of Charon" and a "final swingby" sets apoapse to 227,000 km —
but the paper gives no Charon flyby periapsis number. The detailed Charon-GA trajectory
geometry is deferred to **Finley et al. [5]** (the #279 Finley Pluto-tour paper). So
**Charon remains a flyby-floor gap** after this paper.

**Satellite observation distance ≠ flyby altitude.** Spec 2 ("observe each of the five
satellites ≥6 times at radial distance <10,000 km") is an *observation* range, NOT a
flyby altitude — not ingested as a floor.

**Table 2 — Pluto-system physical parameters (cross-source).** Pluto GM 869.339, Charon
GM 106.25 km³/s²; "Dimension" column is **diameter** despite the "Spherical radius"
label (Pluto 2376 = 2×1188.3 km; Charon 1212 = 2×606 km — same diameter/radius trap the
Stone-Miner 1986 table had). Consistent with our code radii (Pluto 1188.3, Charon 606.0);
do NOT read 2376/1212 as radii. Styx GM 0.001, Nix 0.003, Kerberos 0.0011, Hydra 0.0032
(triaxial-ellipsoid bodies) — negligible-mass, screen self-prunes.

## Other content (context, not ingested)
- Tour: total ΔV < 400 m/s (Charon-GA-dominated), ≤3 yr, hydrazine monoprop (legs too
  short for EP); 8 phases (capture → satellite-flyby setup → flyby phases → Pluto
  close-in atmospheric → Charon close-flyby → escape).
- Post-Pluto: Pluto→Quaoar low-thrust rendezvous (Table 3), ~7 yr cruise; 14 of 67 KBOs
  reachable for flybys, Quaoar/Ixion viable for orbit. Table 4 lists feasible secondary
  KBO flyby targets (Arrokoth/2014 MU69 among them).
- Mass: 2732 kg launch, NEXT-C ion ×2, 4 RTGs; Table 5/6/7 budgets.

## Catalogue / task implications
- **#429 gap partially closed:** Pluto planet now has a **sourced** design floor (100 km),
  upgrading the convention. Charon + small Pluto moons remain gaps (periapsis numbers
  live in Finley [5]).
- **Reference [5] resolves the #279 "does the Finley paper exist?" question — IT DOES.**
  Full citation from Stern's reference list: **Finley, T., Barth, E., Howett, C.,
  Zangari, A., Tapley, M., Scherrer, J., and Stern, A., "An Orbital Tour of Pluto and
  Its Moons," *Journal of Spacecraft and Rockets* (to be published)** — a real,
  separate, Finley-first-author JSR paper (in-press 2020 → likely published 2020-2021).
  This corrects the #279 disposition note ("no separate Finley publication appears to
  exist"); acquisition re-opened as a precise target — it holds the Charon GA periapsis
  altitudes. Reference **[4] Brozović, Showalter, Jacobson & Buie 2015, *Icarus*
  246:317-329** (DOI 10.1016/j.icarus.2014.03.015) is the Pluto-system masses source
  behind Table 2 and our `satellites.py` Charon/Nix/Hydra GMs.
- **Not catalogue-admissible** as a cycler/tour row: it is a notional study, not flown,
  and the tour is a one-off orbital tour (not a repeating cycler). Recorded as corpus +
  flyby-floor provenance only.
- The Charon-GA Pluto-system tour and Pluto↔Quaoar low-thrust leg are of interest to the
  moon-tour (#306) and low-thrust (#309) lanes as a *reference*, not an anchor.

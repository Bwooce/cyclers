# Forward-citation sweep 2 — the next wave (2026-06-13)

**Task:** #215 — successor to `docs/notes/2026-06-11-forward-citation-sweep.md` (sweep 1). All of
sweep 1's wave-2 acquisitions (Ross AAS 25-621, Liang CGE, Wittal IAC-22, Cuevas 2023, Saloglu
pair, Ellison) are mined; this sweep hunts the next wave and re-checks open acquisitions.
**Method:** WebSearch only, ~26 queries. The Semantic Scholar citations API was again unreachable
from this environment (fetch tooling denied), so — same caveat as sweep 1 — coverage is
search-engine recall, not exhaustive citation-graph traversal.
**Triage lens (unchanged):** HIT only for new sourced members/ICs, real-ephemeris reproduction or
maintenance ΔV, new construction methods, or E-M / moon-tour extension.

---

## Acquisition shortlist — FETCH-NOW (all free)

1. **Braik, A. & Ross, S.D., "Orbital Networks in the Three-Body Problem," arXiv:2605.31543
   (submitted 2026-05-29).** https://arxiv.org/abs/2605.31543 — **first forward citation of
   AAS 25-621 found anywhere** (same lab). Reachable-set framework: finite-ΔV / finite-TOF
   reachable-set overlaps infer accessibility between representative CR3BP periodic-orbit families
   on a common Jacobi level, assembled into a weighted orbital network whose nodes explicitly
   include **cycler** families alongside LPO, resonant, and DPO families. Method-grade HIT: this
   is family-to-family accessibility quantification — directly relevant to our family-selection /
   basin problems and to connecting the Ross E-M cycler families into the wider cislunar family
   set. May also restate or extend the AAS 25-621 family representatives (check for new tuples).
2. **Şaloğlu, K. & Taheri, E., JAS version-of-record is now OPEN ACCESS.**
   DOI 10.1007/s40295-025-00528-0, *J. Astronaut. Sci.* 72:54 (2025). Direct free PDF:
   https://link.springer.com/content/pdf/10.1007/s40295-025-00528-0.pdf (arXiv mirror
   arXiv:2501.01583). Sweep-1 target re-check that flipped from gated to free — fetch the typeset
   VoR to supersede/diff the held preprint (page-citable, final numbers).
3. **JPL Three-Body Periodic Orbits API — access route resolved.** Documented at
   https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html ; endpoint
   `GET https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=<system>&family=<family>[&branch|&libr]`
   (v1.0). Families include `resonant`, `dro`, `dpo`, `halo`, `lyapunov`, etc., for earth-moon and
   other systems. Machine-readable, JPL-sourced ICs — a sourced-rows feed for the catalogue's
   CR3BP lanes and an independent cross-check source for our own continuation output.
4. **"Identifying Fixed Points in the Three-Body Problem Using a High-Order Transfer Map,"
   arXiv:2509.12671 (Sept 2025).** https://arxiv.org/abs/2509.12671 — exhaustive periodic-orbit
   search in the Earth-Moon CRTBP via differential-algebra high-order transfer maps over Poincaré
   sections (planar to 9 revs, spatial to 4). Recovers DRO/Lyapunov families and reports **a
   previously undocumented family with hybrid DRO–Lyapunov character**. Method HIT + potential new
   sourced family for the CR3BP lane.
5. **Der, G.J., "An Elegant State Transition Matrix," *J. Astronaut. Sci.* 45(4), 1997** (also
   AIAA 96-3660). Free PDF: https://link.springer.com/content/pdf/10.1007/BF03546398.pdf —
   universal-variables Keplerian STM that eliminates the secular terms of the Goodyear/Battin
   forms and is the direct free substitute for the still-gated **Shepperd 1985** (Celest. Mech.
   35:129–144, DOI 10.1007/BF01227666). Unblocks the analytic-STM acquisition without buying
   Shepperd.
6. **Pellegrini, E., PhD dissertation, UT Austin 2017, "Multiple-shooting differential dynamic
   programming with applications to spacecraft trajectory optimization."**
   https://repositories.lib.utexas.edu/handle/2152/61764 — free route to the Pellegrini & Russell
   2016 STM-accuracy content (JGCD DOI 10.2514/1.G001920 remains gated); the dissertation covers
   STM quality/computation in the multiple-shooting context.
7. **Andreu 1998 quasi-bicircular thesis — candidate free copy located.** The UB Dynamical Systems
   Group preprint archive (maia.ub.edu/dsg) hosts 1998 preprints; search surfaced
   `http://www.maia.ub.es/dsg/1998/9801mangel_e.ps.gz` ("9801" + "M. Angel" matches M.À. Andreu,
   "The Quasi-Bicircular Problem"). **Verify the exact file on the DSG index before relying on
   it** (filename came via search-engine snippet; PostScript, will need conversion).
8. **Liang, G., Niu, J., Yang, H., Li, S., Li, H., "Utilization of Solar Gravity Perturbation in
   Moon-Aided Jovian Capture," *Space: Science & Technology*, DOI 10.34133/space.0285 (July
   2025).** Open access: https://spj.science.org/doi/pdf/10.34133/space.0285 — same NUAA group as
   the CGE triple-cycler paper; SGP + multi-moon-aided capture into the Jovian system in both
   CR3BP and high-fidelity models. Moon-tour-lane feed (capture into tour/cycler geometry); check
   whether it cites/extends the CGE families. Lower priority than 1–4.

---

## Per-target findings

### Target 1 — forward citations of pillar sources

| Citing work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| Braik & Ross, "Orbital Networks in the Three-Body Problem," arXiv:2605.31543 | 2026, arXiv | Cites Ross & Roberts-Tsoukkas 2025; reachable-set orbital network over CR3BP families incl. cycler/resonant/DPO/LPO nodes | **HIT (method)** | free (arXiv) |
| Liang, Niu, Yang, Li & Li, *Space: Sci. & Technol.*, 10.34133/space.0285 | 2025 | Same group as Liang CGE; SGP-assisted multi-moon Jovian capture; moon-tour lane context | marginal-HIT (lane feed) | free (SPJ OA) |

No other 2025–2026 citing works with trajectory numbers were found for Russell & Ocampo,
McConaghy/Longuski, Byrnes/Longuski/Aldrin 1993, or Jones/Hernandez/Jesick AAS 17-577 beyond what
sweep 1 already recorded (see no-change section).

### Target 2 — new 2025–2026 publications by lane

| Work | Year / venue | Relevance | Verdict | Access |
|---|---|---|---|---|
| "Identifying Fixed Points in the Three-Body Problem Using a High-Order Transfer Map," arXiv:2509.12671 | 2025, arXiv | Exhaustive E-M CRTBP PO search (DA/HOTM); previously undocumented hybrid DRO–Lyapunov family | **HIT** (method + possible new family) | free |
| Beolchi, Pontani, Howell, Pozzi, Swei, Fantino, "Low-Energy Round-Trip Trajectories to Near-Earth Objects Using Low Thrust," arXiv:2603.27683 | 2026 (Mar), arXiv | Sun-Earth CR3BP manifolds + heliocentric 2BP patched; mass-producible round trips, no inner-loop optimization. Round-trip lane method, not cyclers | BACKGROUND | free |
| "Round-Trip Mars Missions in the 2031 Window," preprints.org 202510.1072 (v2) | 2025 (Oct), preprint, CC BY | Claims 153-day round-trip via plane geometry of asteroid 2001 CA21; heavily media-hyped, **not peer-reviewed** — treat with strong skepticism; no sourced rows for us | BACKGROUND (skeptical) | free (preprints.org) |
| Quarta, A.A., "Round-Trip Heliocentric Trajectories for Continuous-Thrust CubeSats," *AM&S*, 10.1007/s42496-025-00297-x | 2025, Springer | Min-time round trips under propellant constraint; CubeSat SEP; powered round-trip lane only | BACKGROUND | gated |
| "Discovery of 10,059 new three-dimensional periodic orbits of general three-body problem," arXiv:2508.08568 | 2025, arXiv | General (non-restricted) TBP periodic orbits — wrong regime for mission catalogue | BACKGROUND | free |
| "Generative Design of Periodic Orbits in the Restricted Three-Body Problem," arXiv:2408.03691 | 2024, arXiv | ML-generative PO design; genome-roadmap adjacent (cf. Ozaki NN surrogate, sweep 1 bonus) | BACKGROUND | free |
| "Connecting Earth and Moon via the L1 Lagrangian point," arXiv:2502.11694 | 2025, arXiv | E-M transfer via L1; transfer not cycler | BACKGROUND | free |
| "Automated Tour Design in the Saturnian System," arXiv:2210.14996 (CMDA 2023, 10.1007/s10569-023-10179-8) | 2023 | Grid-based dynamic-programming moon-tour optimization Titan→Enceladus; Saturnian moon-tour lane backfill (not in sweep 1) | marginal-HIT (lane backfill) | free (arXiv) |

### Target 3 — open-acquisition re-checks

| Item | Status this sweep |
|---|---|
| Ross & Roberts-Tsoukkas 2026 journal VoR | **Not yet published.** No journal version found in JGCD/CMDA/JAS as of 2026-06; the VSGC-hosted manuscript (Apr 2026, https://vsgc.odu.edu/wp-content/uploads/2026/04/Roberts-Tsoukkas_Michael_Cycler-Journal-Paper.pdf) remains the latest text. Re-check next sweep. |
| Liang CGE JGCD VoR (10.2514/1.G008387) | **Still gated** at AIAA ARC (JGCD 48(1):146–155). The ResearchGate author full-text (pub. 383986230) remains the only free copy. No change. |
| Şaloğlu JAS VoR (10.1007/s40295-025-00528-0) | **FLIPPED TO FREE** — Springer open access (shortlist #2). |
| Shepperd 1985 (Celest. Mech. 35) | Still gated (Springer; NTRS 19850050934 is metadata-only). **Free substitute found:** Der 1997 (shortlist #5). |
| Pellegrini & Russell 2016 (10.2514/1.G001920) | JGCD VoR still gated. **Free substitute found:** Pellegrini's UT Austin dissertation (shortlist #6). |
| Friedlander 1986 (AIAA 86-2009, 10.2514/6.1986-2009) | Still gated at AIAA ARC; ResearchGate metadata-only. No free copy. |
| Andreu bicircular thesis | Candidate free copy on UB DSG preprint server (shortlist #7, verify). |
| JPL Three-Body Periodic Orbit Catalog | Access route resolved: documented public API (shortlist #3). |

### Target 4 — errata-ledger corrections

Searched for published errata/corrigenda for: Shakouri/Kiani/Pourtakdoust 2019 (shape-based
multiple-impulse), Ellison/Conway/Englander/Ozimek (analytic gradients), Liang CGE 2024,
Şaloğlu & Taheri iso-impulse, Ross AAS 25-621. **No published correction found for any of them.**
`/errata` page statuses unchanged. (Secondary note: the Şaloğlu VoR going OA means our
ledger entry can now be checked against the typeset final text, not just the arXiv preprint.)

---

## No-change section (clean negatives — do not redo next sweep without a new method)

- **Ross & Roberts-Tsoukkas journal VoR**: not yet published anywhere (checked JGCD, CMDA, JAS,
  general search). Re-check ~Q4 2026.
- **Forward citations of AAS 25-621**: exactly one found (Braik & Ross arXiv:2605.31543, same
  lab). No external citers yet — the paper is <1 yr old.
- **Forward citations of Liang CGE 2024**: none found beyond the group's own follow-on
  (10.34133/space.0285). No external citers with numbers.
- **Russell & Ocampo 2004/2005**: no new 2025–2026 citing work with cycler numbers beyond
  sweep 1's table.
- **McConaghy/Longuski (S1L1 / two-synodic)**: no new 2025–2026 citers; Wilde JSR line unchanged.
- **Byrnes/Longuski/Aldrin 1993**: no new quantitative citers (scite.ai report page exists at
  https://scite.ai/reports/cycler-orbit-between-earth-and-Wz4aVd but was not retrievable
  in-session; a future sweep with API access could mine it).
- **Jones/Hernandez/Jesick AAS 17-577 / VEM triple cyclers**: no new citers beyond the Liang
  lineage already held.
- **Earth-Mars cycler lane 2025–2026**: no new publication with sourced member rows. The only
  new E-M-adjacent items are round-trip papers (Target-2 table), none cycler-grade.
- **Saturnian / Enceladus-Titan moon-tour cyclers 2025–2026**: nothing new; latest remains the
  2023 automated-tour-design line (arXiv:2210.14996, listed above as backfill).
- **2026 AIAA/AAS Space Flight Mechanics Meeting (Orlando, Jan 2026)**: no cycler paper surfaced
  via search; the program (space-flight.org/docs/2026_winter/) is not search-indexed at paper
  level. Re-check when proceedings/NTRS postings appear.
- **IAC-25 (Sydney, Sep–Oct 2025)**: no cycler papers surfaced in indexed proceedings.
- **Wittal (KSC) follow-on work**: no post-IAC-22 cycler publication found.
- **Fornari & Pontani 2020 (10.1007/s42496-020-00050-6)**: still gated, no free copy; Springer
  explicitly offers no shared link. Pontani's 2025–2026 output is round-trip/NEO work (see
  arXiv:2603.27683), not new cycler families.
- **Hughes et al. JSR fast free returns (10.2514/1.A33293)**: still gated; no Purdue/NTRS free
  copy found.
- **Rogers et al. Acta Astronautica 112 (2015)**: still gated; the free Purdue AAC copy covers
  only the held conference version (AIAA 2012-4746).

## Counts

- Queries run: ~26 (WebSearch only; Semantic Scholar API + direct fetch denied in-session —
  the sweep-1 completeness caveat still stands).
- **FETCH-NOW (free): 8** (shortlist above).
- **ASK-USER (paywalled, worth it): 0 new** — the gated carryover list from sweep 1 (Rogers 2015,
  Fornari & Pontani 2020, Miguel 2023, Hughes JSR, p:q resonant Acta 170, Wilde JSR, Pascarella
  Acta 203) is unchanged, and two of its members (Shepperd, Pellegrini & Russell) now have free
  substitutes that remove the purchase case.
- **BACKGROUND: 7** (Target-2 table).
- Errata corrections found: 0.

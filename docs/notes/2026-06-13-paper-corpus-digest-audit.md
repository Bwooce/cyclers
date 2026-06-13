# Paper-corpus digest audit (2026-06-13)

**Purpose.** Classify every filed paper in the reference corpus (~67 PDFs) by digestion
status so genuinely-undigested cycler/astrodynamics-relevant papers can get per-paper
digest todos. Classification is from filename + whether the project already references the
paper (mining notes in `docs/notes/`, completed tasks recorded in `data/OUTSTANDING.md`,
catalogue citations in `data/catalogue.yaml`, or the in-flight digest todos #230-233).
PDFs were not re-opened.

**Method.** For each paper, multiple keys were grepped (author surname, year, venue id,
arXiv id, topic keywords) across `docs/notes/`, `data/OUTSTANDING.md`, and
`data/catalogue.yaml`. A single substantive author/topic match in a mining note or a
catalogue `orbit_source`/`source_quotes` citation counts as DIGESTED. Bibliography-only
mentions (a paper appearing solely as a `[n]` reference inside another paper's note) do
**not** count as digestion.

**Counts:** A (DIGESTED) = 49 · B (BACKGROUND-TRIAGED) = 9 · C (UNDIGESTED / needs
backfill) = 5 · plus 4 IN-PROGRESS-TODO (#230-233). Total = 67.

---

## BACKFILL TODOS NEEDED (bucket C only)

These are cycler/astrodynamics-relevant, have no mining note, no completed task, and no
catalogue citation (or only a bibliography-reference mention). Each needs a per-paper
digest todo:

1. **Fan 2025 — Electric-sail multi-target trajectory design (Aerospace 12:00196).**
   Proposed scope: triage-skim — is the multi-target E-sail sequencing / propellantless
   continuous-thrust method relevant to cycler establishment or maintenance, or
   out-of-scope (E-sail-only)? Record verdict + any reusable sequencing idea.

2. **Szebehely 1967 — *Theory of Orbits: The Restricted Problem of Three Bodies* (book).**
   Proposed scope: complete the already-flagged #185 CR3BP periodic-orbit / libration
   foundational mining (Jacobi integral, periodic-orbit families, stability) as an
   independent cross-check reference for the CR3BP lane.

3. **Silvestrini & Lavagna 2022 — Deep-learning ANN for spacecraft GNC (Drones 6:270).**
   Background, needs one-line triage. Proposed scope: one-line REFERENCE-ONLY /
   OUT-OF-SCOPE verdict (ANN closed-loop GNC, not cycler-data) — no triage note exists yet.

4. **Li, Topputo & Baoyin 2019 — Neural time-optimal orbit-raising EP-to-GEO
   (arXiv:1909.08768).** Background, needs one-line triage. Proposed scope: one-line
   verdict (low-thrust ML guidance; likely REFERENCE-ONLY) — no triage note exists yet.

5. **Singh & Junkins 2022 — Stochastic learning of the extremal field for low-thrust
   guidance (Sci. Rep.).** Background, needs one-line triage. Proposed scope: one-line
   verdict (low-thrust extremal-field ML guidance; likely REFERENCE-ONLY) — no triage
   note exists yet.

> Note: **Putnam & Braun 2005** (entry-system options) appears only as bibliography
> reference [3] inside the Donahue-Duggan note — never actually digested. It is
> EDL/aerocapture, tangential to cycler trajectory mining; if a backfill todo is wanted it
> is a one-line OUT-OF-SCOPE triage, but it is borderline-tangential rather than clearly
> relevant, so it is held in bucket B (acquired-as-background) rather than promoted to C.

---

## A. DIGESTED (49)

Each has a dedicated mining note, a completed task in `data/OUTSTANDING.md`, and/or a
catalogue citation.

| Paper | Evidence |
|---|---|
| Russell 2004 dissertation | `2026-06-07-russell-2004-dissertation-method-mining.md`, `-continuation-deepdive.md`, `-member-tables-transcription.md`; catalogue `orbit_source` (9.7k hits) |
| Hollister & Rall 1970 (NASA-CR periodic orbits) | `2026-06-07-hollister-rall-1970-periodic-orbits-mining.md`, `-rall-1970-appendices-transcription.md`; catalogue |
| Hollister & Menning 1970 (Earth-Venus periodic swingby, JSR 7-10) | catalogue family `hollister-menning-1970-ev-orbit-01..15`; OUTSTANDING (E-V family individuated) |
| McConaghy 2004 dissertation | `2026-06-10-mcconaghy-2004-dissertation-mining.md`, `-table71-reproduction.md`; catalogue |
| Rogers 2012 — V∞-leveraging cyclers (AIAA 2012-4746) | catalogue `orbit_source: rogers-2012-t1`, extensive `source_quotes` (371 hits) |
| Jones, Hernandez & Jesick 2017 — low-excess-speed VEM triple cyclers (AAS 17-577) | `2026-06-05-jones-aas17-577-vem-mining.md`, `2026-06-07-jones-aas17-577-method-deepdive.md`; catalogue |
| Liang 2024 — Callisto-Ganymede-Europa triple cyclers (JGCD) | `2026-06-11-liang-2024-cge-triple-cyclers-mining.md`, `2026-06-13-liang-abc-reproduction.md`, `-member-d-nbody.md`; catalogue |
| Genova & Aldrin 2015 — Earth-Moon cycler (AAS-15) | `2026-06-10-genova-aldrin-2015-mining.md` |
| Genova & Aldrin 2015 — free-return Earth-Moon cycler (AAS-15) | `2026-06-10-genova-aldrin-2015-mining.md` (free-return arc covered) |
| Genova 2016 — Phobos/Deimos PADME (AIAA 2016-5681) | `2026-06-05-v42-backfill-sweep.md` §7-context (PADME) |
| Ross, Roberts & Tsoukkas 2025 — stable ballistic Earth-Moon cyclers (AAS 25-621) | `2026-06-11-ross-roberts-tsoukkas-2025-mining.md`, `2026-06-12-ross-adoption-results.md`, `2026-06-13-ross-v2-assessment.md`, `-v2-longspan-evidence.md`; catalogue |
| Roberts-Tsoukkas — VSGC multi-orbiter cyclers (student summary) | covered by ross-roberts-tsoukkas mining + forward-citation sweeps |
| Campagnola & Russell 2009 — endgame Part A (VILM/leveraging graph, AAS 09-224) | `2026-06-05-endgame-tisserand-mining.md` |
| Campagnola & Russell 2009 — endgame Part B (multibody T-P graph, AAS 09-227) | `2026-06-05-endgame-tisserand-mining.md` |
| Vasile & Campagnola 2009 — low-thrust MGA Europa DFET (JBIS) | `2026-06-07-vasile-campagnola-dfet-method-mining.md`, `2026-06-05-vasile-tables-retranscription.md` |
| Vasile, Summerer & De Pascale 2005 — Earth-Mars evolutionary branching | `2026-06-05-v42-backfill-sweep.md`, `2026-06-07-external-algorithms-survey.md` |
| Ceriotti 2010 — global optimisation MGA (Glasgow PhD) | `2026-06-07-ceriotti-2010-mga-global-opt-mining.md` |
| Hiraiwa 2026 — lobe dynamics cislunar transfers (arXiv:2602.17444) | `2026-06-07-hiraiwa-lobe-dynamics-method-mining.md`, `2026-06-05-vasile-hiraiwa-scan.md` |
| Takao 2025 — Saturn Trojan 2019 UO14 mission analysis | `2026-06-07-takao-2025-mpga-1dsm-mining.md` (promoted from marginal triage) |
| Zhou, Armellin 2025 — single-impulse reachable set polynomials (arXiv:2502.11280) | `2026-06-07-zhou-2025-da-reachable-sets-mining.md` (promoted from marginal triage) |
| Cuevas del Valle 2023 — optimal Floquet stationkeeping CR3BP (Aerospace) | `2026-06-11-cuevas-del-valle-2023-floquet-mining.md` |
| Cuevas del Valle 2026 — fuel-optimal CR3BP rendezvous MPC (EuroGNC) | `2026-06-10-cuevas-del-valle-2026-cr3bp-mpc-mining.md`, `2026-06-13-cuevas-2026-l1-halo-seed-run.md` |
| Şaloğlu 2023 — infinitely many iso-impulse trajectories (JGCD) | `2026-06-10-saloglu-2023-iso-impulse-mining.md` |
| Şaloğlu 2025 — iso-impulse 3D classification/feasibility (arXiv) | `2026-06-10-saloglu-2025-iso-impulse-3d-mining.md` |
| Ellison 2018 — analytic-gradient bounded-impulse two-sided shooting (JGCD) | `2026-06-10-ellison-2018-analytic-gradients-mining.md` |
| Iorfida 2016 — geometric perpendicular-thrust optimization (JGCD) | `2026-06-10-iorfida-2016-perpendicular-thrust-mining.md` |
| Shakouri 2019 — shape-based multiple-impulse coplanar maneuvers (arXiv) | `2026-06-10-shakouri-2019-shape-based-mining.md` |
| Vallado 1991 — methods of astrodynamics (USAFA TR-91-6) | `2026-06-10-vallado-1991-tr916-mining.md` |
| Guzman et al. 2002 — primer-vector optimization survey (IAC-02-A.6.09) | `2026-06-07-guzman-2002-primer-survey-mining.md` |
| Beeson, Englander & Hughes 2015 — EMTG/GMAT low-thrust toolchain (AAS 15-278) | `2026-06-07-beeson-2015-emtg-gmat-toolchain-mining.md` |
| Englander & Englander 2014 — tuning MBH (ISSFD24-S7-3) | `2026-06-07-englander-2014-mbh-tuning-mining.md` |
| Ozimek 2019 — LinX low-thrust MGA optimization (AAS 19-348) | `2026-06-07-ozimek-linx-aas19-348-mining.md` |
| Hughes, Edelman & Longuski 2014 — fast Mars free returns via Venus GA (AIAA 2014-4109) | `2026-06-07-hughes-2014-fast-mars-free-returns-mining.md` |
| Okutsu & Longuski 2002 — Mars free returns via Venus GA (JSR 39-1) | `2026-06-07-okutsu-tito-free-returns-mining.md` |
| Tito, MacCallum, Carrico 2013 — manned Mars free-return 2018 (IEEE Aerospace) | `2026-06-07-okutsu-tito-free-returns-mining.md` |
| Landau & Longuski 2006 — human Mars trajectories Pt 1 impulsive (JSR) | `2026-06-04-agrawal-landau-howe-mining.md`; catalogue (Landau) |
| Landau & Longuski 2009 — comparative assessment human Mars technologies | `2026-06-04-agrawal-landau-howe-mining.md` |
| Agrawal 2022 — orbital logistics architecture sustainable Mars (Purdue PhD) | `2026-06-04-agrawal-landau-howe-mining.md` |
| Howe 2025 — tackling Mars cycler design head-on (ICES-2025-555) | `2026-06-04-agrawal-landau-howe-mining.md`; OUTSTANDING (Howe) |
| Donahue & Duggan 2022 — Boeing Mars-2033 human flyby (IAC-22) | `2026-06-07-donahue-duggan-2022-mars2033-flyby-mining.md` |
| Wittal 2022 — Earth-Moon cycler lunar logistics (IAC-22-C1.6.6) | `2026-06-11-wittal-2022-iac-mining.md` |
| Patel 2019 — Earth-Mars cycler vehicle conceptual design (FIT thesis) | `s1l1-target-topology-mining.md`, `2026-06-05-v42-backfill-sweep.md` §7; OUTSTANDING H.5 |
| Pascarella et al. 2022 — Pony Express low-thrust E-M cycler (AAS-22-015) | OUTSTANDING H.3 (task #38, two-paper split) |
| Sanchez-Net et al. 2022 — Cycler orbits / Pony Express near-ballistic (JSR) | OUTSTANDING H.4 (task #38); catalogue `sanchez-net-2022-*` entries |
| CCSDS 2023 — Orbit Data Messages 502.0-B-3 blue book | `2026-06-05-ccsds-odm-502-mining.md` |
| Rickman — intro orbital mechanics (NASA-NESC slides) | `2026-06-10-rickman-nesc-slides-triage.md` |
| Zhang 2024 — neural angle-only OD Earth-Moon libration (Remote Sensing) | `2026-06-07-zhang-2024-neural-od-mining.md` |
| Zhang, Acciarini 2026 — pretrained approximators low-thrust cost/reachability (arXiv:2605.26790) | `2026-06-07-ml-surrogate-investigation.md` Paper 1 |
| Zhang 2026 — neural porkchop low-thrust asteroid rendezvous (Astronautics) | `2026-06-07-ml-surrogate-investigation.md` Paper 2 |

---

## B. BACKGROUND-TRIAGED (9)

Acquired as part of an algorithms / ML / GNC research sweep and recorded as
background/method-reference, NOT cycler-data, with a triage note citing them.

| Paper | Triage note + verdict |
|---|---|
| Venigalla, Englander & Scheeres 2020 — low-thrust missed-thrust margin (AAS 20-438) | `2026-06-07-marginal-papers-triage.md` #1 — OUT-OF-SCOPE |
| Sinha & Beeson 2025 — initial-guess generation robust low-thrust (arXiv:2501.06694) | `2026-06-07-marginal-papers-triage.md` #2 — REFERENCE-ONLY |
| Hu, Yang, Li 2024 — robust low-thrust GA via RL (AIAA G009427) | `2026-06-07-marginal-papers-triage.md` #5 — OUT-OF-SCOPE |
| Blender & Singh 2025 — uncertainty-aware GBDT continuous-thrust guidance (AAS 25-524) | `2026-06-07-marginal-papers-triage.md` #6 — OUT-OF-SCOPE |
| Viavattene & Ceriotti 2021 — ANN multiple-NEA rendezvous continuous thrust (JSR) | `2026-06-07-ml-surrogate-investigation.md` Paper 3 — reference |
| Ozaki 2022 — NN surrogate global cycler search (arXiv:2111.11858) | `2026-06-11-ml-surrogate-trio-triage.md` — blueprint DEFERRED (training-data floor) |
| Leifsson 2022 — global surrogate modeling NN uncertainty (ICCS) | `2026-06-11-ml-surrogate-trio-triage.md` — background-only |
| Wu 2024 — physics-informed ML review condition-monitoring (ESWA) | `2026-06-11-ml-surrogate-trio-triage.md` — background-only |
| Putnam & Braun 2005 — entry-system options human return Moon/Mars (AIAA 2005-5915) | reference-only mention in `2026-06-07-donahue-duggan-...md` [3]; EDL/aerocapture, tangential — held as acquired-background |

---

## C. GENUINELY UNDIGESTED (5)

See "BACKFILL TODOS NEEDED" at the top for proposed per-paper digest scopes.

1. Fan 2025 — Electric-sail multi-target trajectory design (Aerospace 12:00196) —
   propellantless continuous-thrust multi-target sequencing; relevance to cycler
   establishment/maintenance unassessed.
2. Szebehely 1967 — *Theory of Orbits* (book) — foundational CR3BP reference; mining
   already flagged pending as #185, never executed.
3. Silvestrini & Lavagna 2022 — deep-learning ANN spacecraft GNC (Drones 6:270) — no
   triage note; needs one-line verdict.
4. Li, Topputo & Baoyin 2019 — neural time-optimal orbit-raising EP-to-GEO
   (arXiv:1909.08768) — no triage note; needs one-line verdict.
5. Singh & Junkins 2022 — stochastic learning of extremal field for low-thrust guidance
   (Sci. Rep.) — no triage note; needs one-line verdict.

---

## IN-PROGRESS-TODO (4) — already on digest todos #230-233

Identified and queued in `2026-06-13-forward-citation-sweep-2.md`; not undigested, not
yet fully mined. NOT to be re-filed.

| Paper | Todo |
|---|---|
| Braik & Ross 2026 — orbital networks in the three-body problem (arXiv:2605.31543) | #230 |
| "Identifying fixed points in the three-body problem using a high-order transfer map" (arXiv:2509.12671) | #231 |
| Şaloğlu & Taheri 2025 — iso-impulse 3D, JAS version-of-record (s40295-025-00528-0) | #232 (VoR diff vs held preprint) |
| Shepperd 1985 — universal Keplerian state-transition matrix (Celest. Mech. 35) | #233 |

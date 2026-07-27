# Digest: Kumar & Moreno 2025 (AAS 25-677, preprint)

**Paper:** "Networks of Periodic Orbits in the Earth–Moon System Through a Regularized and
Symplectic Lens." 2025 AAS/AIAA Astrodynamics Specialist Conference (Preprint), 20 pages.
**Authors:** Bhanu Kumar (James Van Loo Postdoctoral Assistant Professor / NSF Postdoctoral
Research Fellow, Dept. of Mathematics, University of Michigan), Agustin Moreno (Junior Professor,
Institute for Mathematics, Heidelberg University).
**PDF source:** direct S3-hosted conference-proceeding link,
`https://s3.amazonaws.com/amz.xcdsystem.com/A464D031-C624-C138-7D0E208E29BC4EDD_abstract_File25843/PreprintPaperUpload_677_0726010120.pdf`.
**Filed:** `kumar-moreno-2025-networks-periodic-orbits-earth-moon-regularized-symplectic-aas-25-677.pdf`
(private `cyclers_pdf` repo).
**Acquired/digested:** 2026-07-27 (`#728`).
**Text layer:** native LaTeX-produced text-layer PDF confirmed via `pdffonts` (embedded/subsetted
Type-1 Computer Modern + NimbusRoman fonts) and `pdftotext -layout` (clean extraction, 2175 lines).
**No OCR needed.**

**Why in corpus:** flagged by `#598`'s citation-mining pass (`data/OUTSTANDING.md` line ~15278) as
"directly relevant to `#570`'s cycler-network schema" — this task (`#728`) executes that
acquisition + the mandated cross-check against `#570`'s actual schema.

## 1. What the paper actually does

Uses numerical continuation + Kustaanheimo-Stiefel (KS) regularization (through the Earth and Moon
coordinate singularities) + a "symplectic toolkit" (Conley-Zehnder indices, Floer numerical
invariants, via a new public MATLAB tool) to conduct a bifurcation-network study of four planar,
symmetric Earth-Moon CR3BP periodic-orbit families originally catalogued by Broucke's 1968 JPL
report [1]: small circular orbits that become (1) lunar prograde, (2) Earth prograde, (3) lunar
retrograde, (4) Earth retrograde orbits. µ = 1.215058e-2 (Earth-Moon), standard CR3BP normalization
(Eq. 1-3, p.2).

**Headline result 1 — Broucke's 1968 conjecture confirmed (§"Lunar Prograde Orbits..."):** Broucke's
two lunar-prograde families H1 (low/distant prograde orbits stretching toward L1) and H2 (highly
eccentric distant-prograde orbits stretching toward L2) were conjectured by Broucke (p.71 of his
1968 report) to join into a single family if continued far enough. Kumar & Moreno report they found
this junction, "57 years later" — H1's final orbit (Broucke orbit no. 162) connects continuously to
H2's orbit no. 1 (a DPO) via a family passing once through the Earth singularity, then the Moon
singularity, spanning Jacobi constant C = -0.736 to 2.965 (Fig. 4, p.9). They explicitly distinguish
this from Lara & Russell [14]'s prior spatial (non-planar, z≠0) link between double covers of H1/H2 —
the new orbits are planar and symmetric like H1/H2 themselves, genuinely the same continuous family.

**Headline result 2 — bifurcation network of the unified H1-H2 family (Fig. 5-6, pp.9-11):** the
unified H1-H2 family was continued to a fold bifurcation at C = -1.282; along the way, 13
planar-to-spatial (out-of-plane) bifurcations were found (Floquet-multiplier-confirmed), of which
three (at C = 3.136, 3.026, 3.125 — critical H1 LPO/H2 DPO/H2 LPO orbits) generate NEW spatial orbit
families that connect H1/H2 to the well-known L1/L2 Halo-orbit double covers (period-doubling
bifurcations of Halos). Floer-invariant matching is verified across every bifurcation in the diagram
except one segment (CZ index 4↔5 red curve). The authors note prior work found *pieces* of this
picture — Howell & Campbell [16] computed some H2-to-L2-Halo orbits (but not the H2 connection);
Franz & Russell [4] found some H2-DPO-to-L1-Halo orbits (but not the full family or Halo connection,
since it requires KS regularization through the Moon singularity, which their study didn't use).

**Headline result 3 — infinite resonance chains (§ "1:2N and 1:2N+1 chains," p.17, Fig. 14-18):**
Broucke's A1 family (small circular Earth retrograde orbits) is shown to chain, via alternating Moon
and Earth singularity passages that each spawn a new "loop" and flip prograde/retrograde sense, into
a sequence of 1:2, 1:4, 1:6, ... exterior resonant periodic orbits — apparently continuing "ad
infinitum" (no natural termination found even after extensive continuation). Broucke's C family
(small circular lunar retrograde orbits) similarly chains into 1:1, 1:3, 1:5, ... (odd) resonances.
Orbit period vs. Jacobi constant plots (Fig. 16-17) show clear period jumps at each resonance
transition, always occurring as the family passes near/through the Moon singularity at x = 1-µ.
Bifurcations of these two families' single covers are also mapped (Fig. 18): only ~5 spatial
bifurcations found for C, 1 for A1 (far fewer than the prograde-family case), with one bifurcation
(C = -1.300) producing a spatial bridge directly connecting the geocentric C-family to a near-collision
retrograde A1-family orbit (Fig. 20).

**Method note (§"The Symplectic Toolkit," p.6 onward, not deep-read line-by-line but scanned):**
the CZ-index/Floer-invariant machinery is used purely as a *bifurcation detector and family-linkage
validator* (matching Floer numbers across a bifurcation is a necessary condition for the families on
either side to be part of a consistent global picture) — not as a search/discovery method in its own
right. The actual orbit-family discovery is standard numerical continuation + Floquet-multiplier
bifurcation detection; KS regularization is what lets continuation survive passage through the Earth
or Moon singularity, which is the mechanism that produces every "surprising" new connection in this
paper (all three headline results above hinge on a family passing through a coordinate singularity
that non-regularized continuation would have stopped at).

## 2. Cross-check — does this bear on `#570`'s cycler-network schema?

**Read `data/cycler_network.schema.json` in full before answering.** `#570`'s schema (v1.0) registers
`data/cycler_networks.yaml`, a **derived registry of SETS of already-catalogued cyclers sharing
phasing/downlink cadence** — modeled directly on Sanchez Net et al. 2022's Earth-Mars Cycler Orbit
(EMCO) *fleet* concept. Its fields are: `member_cycler_ids` (references into `data/catalogue.yaml`),
`downlink_cadence` (a per-encounter schedule of dates + M-E transit durations), and
`per_member_taxi_insertion_cost_kms` (an Earth-insertion ΔV surrogate per member, explicitly NOT an
inter-cycler transfer cost — the schema's own top-level description flags that as future M5/M6-tier
work). It ships genuinely empty (`[]`) — no real network is populated.

**Answer: no, this paper's "network" is a completely different mathematical object from `#570`'s
schema, and offers nothing this project's schema or discovery methods need or are missing.**

- Kumar & Moreno's "network of periodic orbits" is a **bifurcation graph**: nodes are individual
  periodic-orbit families (each a 1-parameter continuum of CR3BP periodic orbits, parameterized by
  Jacobi constant), edges are *bifurcation events* (a Floquet multiplier crossing 1, producing a new
  family branching off an old one) with symplectic-invariant (CZ index / Floer number) bookkeeping
  attached to each edge for validation. This is a structural/topological map of *how single-spacecraft
  periodic-orbit families relate to each other within one CR3BP system* — a mathematical object
  living entirely inside `data/catalogue.yaml`'s conceptual domain of "one orbit" (or rather, one
  orbit *family*), not a relation *between distinct catalogued cyclers*.
- `#570`'s schema's "network" is an operational/fleet concept: a SET of independently-flyable cycler
  vehicles (each already a complete, separately-catalogued trajectory) coordinated by a shared
  downlink/phasing schedule, directly modeling how a mission designer would field multiple cyclers as
  one constellation. The two concepts share only the English word "network"; nothing in Kumar &
  Moreno's schema, method, or scope (member IDs, downlink cadence, insertion cost) maps onto their
  bifurcation-graph object, and nothing in their bifurcation-graph maps onto a set of independently
  flyable, already-catalogued cyclers.
- **Does the paper describe/catalogue any SPECIFIC Earth-Moon periodic-orbit network member worth a
  future cross-check against this project's own `catalogue.yaml` rows?** This project's catalogue has
  no Earth-Moon `cycler`/`quasi_cycler` rows built from Broucke-lineage lunar-prograde/retrograde
  families or Halo-adjacent bifurcation branches (a quick domain check: the project's Earth-Moon work
  is concentrated in the QBCP/CCR4BP torus-corrector arc — `#538`/`#544`/`#611`-`#620` — and in the
  cislunar MMR/heteroclinic Kumar-Rawat-Rosengren-Ross lineage already in corpus, `#597`/`#598`/`#621`
  — neither of which is Broucke H1/H2/A1/C-family territory). No specific orbit in this paper (e.g.
  the H1-H2 junction orbits, or the H2-to-L2-Halo bridge family) corresponds to any existing catalogue
  row, so there is nothing to cross-check today. It IS a candidate future search target in its own
  right (their new H1-H2-Halo bridge family, Fig. 6, is flagged by the authors themselves as
  potentially useful for "cislunar PNT spacecraft or even transfers," p.19) but that would be a new
  discovery task, not a cross-check against existing rows — flagging for the coordinating session's
  awareness, not acting on it here (out of this task's scope).
- **Method takeaway worth carrying forward (not a schema gap, a search-method observation):** KS
  regularization through the Earth/Moon coordinate singularity is the load-bearing technique behind
  every new connection this paper finds. A grep of `src/cyclerfinder/core/` and `src/cyclerfinder/search/`
  during this digest pass found no KS-regularization module in this codebase's own CR3BP periodic-orbit
  continuation machinery — the project's continuation code (e.g. the various resonant-family/niching-GA
  searches, `#580`-`#591`) integrates the standard singular equations directly. This paper is a concrete
  demonstration that un-regularized continuation genuinely misses real, connected orbit families whose
  link passes through a primary. This is a plausible future capability-gap candidate (worth a
  standalone scoping task if the project ever wants to search for orbit families that graze a primary),
  not something to build reactively here.

**Conclusion for `#570`:** no correction or extension to the schema is warranted. The paper is
topically Earth-Moon-relevant and methodologically interesting (regularization + bifurcation-network
mapping), but it is not the kind of "network" `#570` was scoped to represent, and it describes no
specific orbit currently in this project's catalogue.

## 3. Mandatory citation-mining pass

All 19 references read; cross-checked against `docs/notes/CORPUS_INDEX.md` (2026-07-27) and the
`cyclers_pdf/papers/` filename listing directly.

**Already in corpus** (no action needed):
- [18] Kumar, Rawat, Rosengren, Ross 2024, IAC-24-C1.9.5 — `kumar-rawat-rosengren-ross-2024-interior-mmr-heteroclinic-earth-moon-IAC-24-C1.9.5.pdf`.
- [19] Rawat, Kumar, Rosengren, Ross 2025, AAS 25-569 — `rawat-kumar-rosengren-ross-2025-exterior-mmr-earth-moon-heteroclinic-AAS-25-569.pdf`.
- [3] Doedel et al. — the *Doedel-Keller-Kernevez* AUTO-foundations papers are in corpus, but this
  paper's specific [3] (Doedel, Romanov, Paffenroth, Keller, Dichmann, Galán-Vioque, Vanderbauwhede
  2007, "Elemental Periodic Orbits Associated With the Libration Points in the CR3BP," IJBC 17(8))
  is a DIFFERENT, more specific paper than what's filed — flagged below as a genuine gap, not
  double-counted as "already in corpus."

**Genuinely new candidates, flagged NOT acquired** (per task instructions):
1. **Broucke, R. 1968**, "Periodic Orbits in the Restricted Three-Body Problem with Earth-Moon
   Masses," JPL Tech. Rep. No. 32 [ref 1]. **HIGH priority** — the foundational source catalogue for
   ALL FOUR orbit families this paper studies (H1, H2, A1, C); already independently flagged as a
   corpus gap by the sibling `#725` Casoliva digest's own citation-mining pass (same reference).
   Two independent papers now flag this same 1968 JPL report — the strongest gap signal in this pass.
2. **Doedel, Romanov, Paffenroth, Keller, Dichmann, Galán-Vioque & Vanderbauwhede 2007**, "Elemental
   Periodic Orbits Associated With the Libration Points in the Circular Restricted 3-Body Problem,"
   IJBC 17(8):2625-2677 [ref 3]. HIGH priority — the L1-L5 libration-point-orbit bifurcation-network
   paper this project's own `#570` scoping note explicitly calls "often referred to by astrodynamicists"
   (per this paper's own Conclusion, p.19) and which this paper positions itself as a "stepping stone"
   sequel to; distinct from the Doedel-Keller-Kernevez AUTO-methodology papers already in corpus.
3. **Franz, C.J. & Russell, R.P. 2022**, "Database of Planar and Three-Dimensional Periodic Orbits and
   Families Near the Moon," J. Astronaut. Sci. 69(6):1573-1612, DOI `10.1007/s40295-022-00361-9`
   [ref 4]. HIGH priority — "recently developed a database of millions of POs near the Moon using a
   grid search" per this paper's own intro; a large, directly comparable lunar-PO census this project's
   own lunar-region search work could benchmark against.
4. **Frauenfelder, U., Koh, D. & Moreno, A. 2023**, "Symplectic methods in the numerical search of
   orbits in real-life planetary systems," SIAM J. Appl. Dyn. Syst. 22(4):3284-3319 [ref 5]. Medium
   priority — the foundational paper for the "symplectic toolkit" (CZ indices, Floer invariants) this
   paper's whole bifurcation-detection method builds on.
5. **Moreno, A., Aydin, C., von Koert, O., Frauenfelder, U. & Koh, D. 2024**, "Bifurcation Graphs for
   the CR3BP via Symplectic Methods," J. Astronaut. Sci. 71(6):51 [ref 9]. Medium priority — a
   directly related prior bifurcation-graph paper by an overlapping author set.
6. **Frauenfelder, U. & Moreno, A. 2023**, "On GIT quotients of the symplectic group, stability and
   bifurcations of periodic orbits (with a view towards practical applications)," J. Symplectic Geom.
   21(4):723-773 [ref 10]. Low-medium priority — more foundational symplectic-methods theory, likely
   heavier on pure math than directly actionable for this project.
7. **Aydin, C. & Batkhin, A. 2025**, "Studying network of symmetric periodic orbit families of the
   Hill problem via symplectic invariants," Celest. Mech. Dyn. Astron. 137(2):12 [ref 15]. Medium
   priority — explicitly cited as the impetus for this paper's own search (the Hill-problem analogue of
   the LPO-to-Halo bridge families found here); a genuine sibling result in a simpler (Hill) model.
8. **Howell, K. & Campbell, E. 1999**, "Three-dimensional periodic solutions that bifurcate from halo
   families in the circular restricted three-body problem," Adv. Astronaut. Sci. 102:891-910 [ref 16].
   Low-medium priority — computed some of the same H2-to-L2-Halo orbits from the Halo side; a partial
   precedent for this paper's own Fig. 6 result.
9. **Lara, M. & Russell, R.P. 2006**, "On the family 'g' of the restricted three-body problem,"
   Monografias de la Real Academia de Ciencias de Zaragoza 30:51-66 [ref 14]. Low priority — the prior
   (non-planar) H1-H2-double-cover link this paper explicitly distinguishes itself from.

None of the above were acquired this pass (per task scope: flag only). Item 1 (Broucke 1968) is now
independently flagged by two different digest passes (`#725` and this one) — the strongest
acquisition-priority candidate to come out of this task.

## Summary answer (for the dispatching session)

**Does this paper's periodic-orbit-network method offer anything `#570`'s schema or this project's
discovery methods don't already have?** No direct schema overlap — `#570`'s "network" (a fleet of
independently-catalogued cyclers sharing downlink cadence) and this paper's "network" (a bifurcation
graph linking periodic-orbit families within one CR3BP system) are different mathematical objects
sharing only a name; no correction to `#570` is warranted. The one genuinely reusable idea is a
*search-method* observation, not a schema one: KS regularization through a primary's coordinate
singularity is demonstrably load-bearing for finding real orbit-family connections this project's own
un-regularized continuation code would structurally miss — a candidate future capability gap, not
acted on here. No specific orbit in the paper corresponds to an existing catalogue row.

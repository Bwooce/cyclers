# Mars-cycler wave acquisition digest (2026-06-17)

**Trigger:** user uploaded 17 papers to cyclers_pdf 2026-06-17 16:39 AET
(commits `069e376`, `6b4ca1a`, `ae070aa` on `Bwooce/cyclers_pdf:main`).
Filed atomically at cyclers_pdf commit `b04d944` per the
`feedback_cyclers_pdf_filing_pattern` standard.

**Scope of this note:** one-paragraph digest per paper + recommended
follow-on action. Where a paper has substantial catalogue / corpus
leverage, a tracked task is created. Where it's reference-only, the
survey entry is the digest.

## Highest-leverage acquisitions (own digest tasks created)

### Russell-Ocampo 2003 (AAS-03-145, preprint of JGCD 2004 27(3):321-335)

**`russell-ocampo-2003-systematic-method-earth-mars-cyclers-direct-return-trajectories-AAS-03-145.pdf`**

The preprint of the paper that anchors the SnLm catalogue's `russell-ocampo-*`
rows. Abstract reports: **24 purely ballistic cyclers with periods of two to
four synodic periods, 92 ballistic cyclers with periods of five or six
synodic periods, and hundreds of near-ballistic cyclers** — a method for
systematic Earth-Mars cycler construction via Lambert + free-return + direct
returns + Earth-generated gravity-assist.

*Action:* this RESOLVES the Russell-Ocampo wishlist item on #116. Deep-read
needed to determine whether the published tables enable any V0→V1 promotion
for catalogue rows currently V0-by-publication-gap (per the
`project_validation_ceiling` memory rule). Tracked as new task.

### Byrnes-Longuski-Aldrin 1993 (JSR 30(3), DOI 10.2514/3.25519)

**`byrnes-longuski-aldrin-1993-cycler-orbit-earth-mars-jsr-doi-10.2514-3.25519.pdf`**

The canonical Aldrin cycler paper. The current KNOWN_CORPUS Aldrin anchor
cites Byrnes-McConaghy-Longuski AIAA 2002-4420 — this is the 1993 progenitor
JSR paper that should be cited alongside. Reading it may also tighten the
catalogue's Aldrin V0 row provenance.

*Action:* deep-read for KNOWN_CORPUS citation strengthening + Aldrin row
provenance audit. Tracked as new task.

### Nock-Friedlander 1987 (Acta Astronautica 15(6-7), DOI 10.1016/0094-5765(87)90189-5)

**`nock-friedlander-1987-elements-mars-transportation-system-acta-astro-doi-10.1016-0094-5765(87)90189-5.pdf`**

"Elements of a Mars Transportation System" — RESOLVES the Friedlander 1986
wishlist entry on #116 (the actual paper is 1987 June; common references say
"Friedlander 1986" or "1986/87" interchangeably). This is the predecessor
paper for the Aldrin cycler concept; likely cites the original
Friedlander-Niehoff-Byrnes-Longuski AIAA-86-2009 work.

*Action:* deep-read to confirm whether this *is* the AIAA-86-2009 archival
version OR a separate Acta Astro paper. May yield a new KNOWN_CORPUS anchor.
Tracked under the Aldrin-canonical task.

### McConaghy-Longuski-Byrnes 2004 (JSR 41(4), DOI 10.2514/1.11939)

**`mcconaghy-longuski-byrnes-2004-analysis-class-earth-mars-cycler-trajectories-jsr-doi-10.2514-1.11939.pdf`**

"Analysis of a Class of Earth-Mars Cycler Trajectories" — the SnLm class
analysis. Cited in the catalogue's `russell-ocampo-*` rows but mostly as
"see also". A deep-read may yield per-member orbital elements for V0→V1
promotion.

*Action:* deep-read paired with Russell-Ocampo 2003. Tracked under the SnLm
task.

### Rogers-Hughes-Longuski-Aldrin 2015 (Acta Astro 112, DOI 10.1016/j.actaastro.2015.03.002)

**`rogers-hughes-longuski-2015-establishing-cycler-trajectories-earth-mars-acta-astro-doi-10.1016-j.actaastro.2015.03.002.pdf`**

"Establishing cycler trajectories between Earth and Mars" — likely covers
the *insertion* / *establishment* problem (i.e. how a spacecraft enters a
cycler in the first place, the `precursor_mga` class from v4.7 scope
expansion). May surface anchors for the catalogue's `precursor_mga` rows.

*Action:* deep-read for `precursor_mga` class anchoring + KNOWN_CORPUS
extension. Tracked as new task.

## Catalogue-strengthening tier (single bulk task created)

### McConaghy-Russell-Longuski 2005 (JSR 42(4), DOI 10.2514/1.8123)

**`mcconaghy-russell-longuski-2005-standard-nomenclature-earth-mars-cycler-trajectories-jsr-doi-10.2514-1.8123.pdf`**

"Toward a Standard Nomenclature for Earth-Mars Cycler Trajectories" — the
SnLm naming convention paper. The catalogue uses `n.m.k±j` notation
consistent with this paper. Deep-read may surface dual-naming clarifications
or canonical IDs to align catalogue.

### McConaghy-Landau-Yam-Longuski 2006 (JSR 43(2), DOI 10.2514/1.15215)

**`mcconaghy-landau-yam-2006-notable-two-synodic-period-earth-mars-cycler-jsr-doi-10.2514-1.15215.pdf`**

"Notable Two-Synodic-Period Earth-Mars Cycler" — focused on S2 cyclers.
Catalogue has `mcconaghy-2006-em-k2` row (cross_validated tier). Worth
checking whether the published values match the catalogue exactly.

## Reference / methodology / off-scope (file-only)

* **Adamo 2025** — `adamo-2025-spanning-earth-mars-chasm-synodic-resonant-waypoints-AIAA-houston-LnL.pdf`. AIAA Houston Section Lunch & Learn talk. Synodic-resonant *waypoints* methodology rather than peer-reviewed family. *File-only.*

* **Conte-Spencer 2018** — `conte-spencer-2018-mission-analysis-earth-to-mars-phobos-dro-acta-astro-doi-10.1016-j.actaastro.2018.06.049.pdf`. Earth-to-Mars-Phobos distant retrograde orbit mission analysis. Tangentially related to cyclers (Mars-Phobos DRO is a parking orbit, not a cycler). *File-only unless #313 Mars-Phobos work reopens.*

* **de la Fuente Marcos & de la Fuente Marcos 2018** — `de-la-fuente-marcos-2018-geometric-characterization-arjuna-orbital-domain-arxiv-1410.4104v2.pdf`. Arjuna-class Earth co-orbital NEOs (1:1 resonance temporary trojans/horseshoes/quasi-satellites). Directly relevant for #308 asteroid-leveraging cycler search. *Tagged for #308 reactivation context if it ever re-fires.*

* **Gravier-Marchal-Culp 1972** — `gravier-marchal-culp-1972-optimal-trajectories-earth-mars-true-planetary-orbits-jota-doi-10.1007-bf00932349.pdf`. The 1972 JOTA paper on optimal Earth-Mars in their true planetary orbits. Foundational, predates cycler concept. *File-only as historical reference.*

* **Hintz 2023** — `hintz-2023-orbital-mechanics-astrodynamics-techniques-tools-space-missions-springer-textbook.pdf`. Springer textbook. *Reference resource.*

* **Jesick 2019** — `jesick-2019-mars-trojan-orbits-continuous-earth-mars-communication-jas-doi-10.1007-s40295-019-00195-y.pdf`. Mars Trojan orbits for continuous Earth-Mars communication. Sun-Mars co-orbital L4/L5 family. *File-only unless circumstellar-trojan work fires.*

* **Kakoi-Howell-Folta 2014** — `kakoi-howell-folta-2014-access-mars-from-earth-moon-libration-orbits-acta-astro-doi-10.1016-j.actaastro.2014.06.010.pdf`. Earth-Moon libration to Mars access. Sun-Earth↔Earth-Moon manifold-transit topic; same area as #316 (which is completed). *File-only as supporting reference.*

* **Montenbruck-Markgraf 2004** — `montenbruck-markgraf-2004-gps-sensor-impact-point-prediction-sounding-rockets-jsr-doi-10.2514-1.1962.pdf`. **Off-topic for cyclers** — GPS receiver paper for sounding rockets that happened to be in the same JSR issue as the McConaghy 2004 cycler paper. Archived for completeness but no cycler relevance.

* **Parker 2007** — `parker-2007-low-energy-ballistic-lunar-transfers-phd-thesis-cu-boulder.pdf`. PhD thesis on low-energy ballistic lunar transfers (WSB / BCT / weak-stability boundary). Relevant for #316 / cislunar BCT lit-pass. *File-only; reference resource for any cislunar BCT discussion.*

* **Russell 2012** — `russell-2012-survey-spacecraft-trajectory-design-strongly-perturbed-environments-jgcd-doi-10.2514-1.56813.pdf`. JGCD survey paper on trajectory design in strongly perturbed environments (Russell's career-survey article). Could yield methodology anchors for `KNOWN_CORPUS` Strange/Campagnola/Russell entry. *Light-read tagged for future KNOWN_CORPUS extension.*

## #116 wishlist status post-wave

**Resolved this wave:**
- ✅ Friedlander 1986 → Nock-Friedlander 1987 Acta Astro (likely same context)
- ✅ Russell & Ocampo 2006 → Russell-Ocampo 2003 AAS-03-145 (preprint of JGCD 2004; the "2006" wishlist year was a citation-style misreading; the paper is 2004 published / 2003 preprint)

**Still outstanding on #116:**
- ❌ Acton 1996 (SPICE design paper)
- ❌ Pellegrini & Russell 2016 (STM continuation)
- ❌ Taheri & Junkins 2020 (low-thrust optimization)
- ❌ JPL 3BP catalog (Stuart, Russell catalogue infrastructure)

## Tasks created

- **#351** (proposed) — SnLm catalogue strengthening via Russell-Ocampo 2003 + McConaghy 2004/2005/2006 + Rogers 2015 deep-read. Goal: identify any catalogue V0→V1 promotion opportunities and KNOWN_CORPUS citation strengthening.
- **#352** (proposed) — Aldrin-canonical deep-read (Byrnes-Longuski-Aldrin 1993 + Nock-Friedlander 1987). Goal: confirm Aldrin V0 rows' provenance, strengthen KNOWN_CORPUS Aldrin anchor citation.

# Citation-integrity audit — body/system mis-citation sweep (#482/#483/#485)

**Date:** 2026-06-26
**Trigger:** a citation hallucination found 2026-06-26 — a Jovian "Galilean
Io-Europa-Ganymede (IEG) triple cycler" claim was at risk of being cited to
Jones-Hernandez-Jesick **AAS 17-577**, which is actually *"Low Excess Speed
Triple Cyclers of **Venus, Earth, and Mars**"* (a heliocentric interplanetary
paper, confirmed from its title page). Root cause: **citation-by-concept-
collision** — "triple cycler" spans dynamical systems — that **propagates**
because a hallucinated citation, once written into our notes, is read back as
fact. See memory `feedback_ground_citations_against_content`.

This note records the structural defense built and the corpus-wide audit run.

## What was built

* **#483 (ground truth):** a sourced `system` field on
  `cyclerfinder.search.literature_check.CorpusAnchor`, extracted from each
  paper's **title/abstract** (the lowest-hallucination, near-verbatim surface),
  with a `system_grounded` property that derives the label from the
  title-sourced `primary` when no explicit override is set. Standard labels:
  `heliocentric`, `earth-moon`, `jovian`, `saturnian`, `uranian`, `neptunian`,
  `pluto-charon`, `mars-system`, `solar-system`. The two concept-collision
  anchors (Jovian IEG vs heliocentric VEM) and the cross-system Restrepo-Russell
  database carry an **explicit** `system` with a grounding comment.
* **#485 (the guard that can't be forgotten):** a frozen-census-style ratchet
  `tests/data/test_citation_integrity.py`. For every catalogue row whose
  `first_published` work is also pinned in `KNOWN_CORPUS`, it asserts the row's
  **claimed bodies ⊆ the cited work's bodies** and the systems agree.
  * **Positive control (MUST flag):** a `{Io, Europa, Ganymede}` jovian claim
    cited to the `{Venus, Earth, Mars}` heliocentric VEM work fails — verified
    the ratchet fires with the exact extra-bodies diagnostic.
  * **Known-good (MUST pass, non-vacuous):** a Galilean IEG claim cited to the
    correct Jovian IEG anchor passes and is asserted to strong-link (not pass by
    matching nothing).
  * **Non-vacuity guard:** asserts ≥20 catalogue rows actually strong-link.

### The body/system link logic (and its scoped coverage)

A catalogue row is linked to a `KNOWN_CORPUS` anchor only by a **strong link**:
the anchor's author surnames are a **subset** of the row's `first_published`
author surnames **and** the systems agree. A single shared surname is NOT enough
— "Longuski" is on Earth-Mars, VEM, and Uranian-tour works alike; an early naive
"any shared surname" pass produced 41 false positives (all the cross-system
surname collisions + the body-token spelling/Earth-Moon-primary artifacts below)
before tightening to author-set-subset. Bodies are normalised before the subset
test: single-letter planet codes (`J`→Jupiter, `S`→Saturn, `U`→Uranus,
`N`→Neptune, `Me`→Mercury, `E`→Earth) are spelled out, the **primary** body is
dropped (a moon-system anchor lists only the secondaries), and minor bodies
(Gaspra/Ida and other asteroid flyby targets) are out of system-integrity scope.

**Coverage caveat (honest):** the ratchet checks **author-consistent** citations
— a row citing a work whose author set it correctly names. The Galilean
hallucination is exactly this class (right authors, wrong system/bodies). A
mis-citation that ALSO names the wrong authors would not strong-link and is out
of this ratchet's reach; that is the #484 follow-on (structured citation keys).

## Audit result — corpus-wide sweep

| metric | value |
|---|---|
| catalogue rows | 319 |
| rows strong-linked to ≥1 KNOWN_CORPUS anchor | 48 |
| body/system mis-citations found | **0** |
| KNOWN_CORPUS anchors, all with a grounded system | 53/53 |

Anchors exercised by the sweep (rows → anchor): Rogers-2015 establishment (20),
McConaghy dissertation (9), McConaghy-Longuski-Byrnes 2004 (7), Roberts-Tsoukkas
& Ross (5), **Jones VEM (4)**, **Liang CGE (4)**, Braik-Ross (3),
**Hernandez/Jones/Jesick IEG (1)**, Tito 2018 (1), D'Amario-Bright-Wolf (1).

**No body/system mis-citation exists in the committed catalogue.** The #468/#470
Jovian/Saturnian moon-tour wave and the #345 mga_tour wave were scrutinised
specifically; all are body/system-consistent.

### The Galilean case — resolved, body/system-CORRECT

The triggering concern does **not** manifest as a literal mis-citation in any
committed artifact. Investigation found:

* The catalogue row `hernandez-2017-jovian-ieg-triple-family` cites the **real,
  distinct Jovian** paper *"One Class of Io-Europa-Ganymede Triple Cyclers"*
  (Hernandez, Jones & Jesick; Adv. Astronaut. Sci. 162, pp. 973-984; Semantic
  Scholar `7e1de63096852b5422107ffc23a9312ea3de54f3`). Its `first_published.note`
  already explicitly distinguishes it from AAS 17-577. **Correct.**
* `data/admission_proposals_468.jsonl` `matched_source` cites the same Jovian
  IEG work, not AAS 17-577. **Correct.**
* Every literal `AAS 17-577` reference in code/notes
  (`grep` across `src/`, `scripts/`, `docs/notes/`, `data/`) attaches to the
  **VEM** (Venus-Earth-Mars) paper. **Correct.** The Jovian paper, where given a
  number, is cited as **AAS 17-462** (`2026-06-11-forward-citation-sweep.md`).

So the propagated **literal** 577→Jovian error never landed in committed
artifacts — the risk the memory flags was that the Jovian citation was attached
by **concept match** rather than by opening the source, even though it happened
to land on the right real paper. It is now **grounded**: the `KNOWN_CORPUS` IEG
and VEM anchors carry explicit, grounded `system` labels + grounding comments,
and `docs/notes/2026-06-26-470-moontour-admission-wave-v2-verdict.md` carries a
grounding addendum (respectful-errata framing: our agent's concept-match, not
the authors' error).

### #480 implication

**A CORRECT Galilean moon-cycler reference DOES exist** and is already in the
corpus: *"One Class of Io-Europa-Ganymede Triple Cyclers"* (Hernandez, Jones &
Jesick, 2017), the Jovian-moon companion to the heliocentric VEM paper. So #480
(IEG real-ephemeris reproduction) has a genuine published target — it does **not**
need to be punted to #482 research for a reference. (Whether our reproduction
matches the paper's specific members is the separate #480 modelling question;
the #470 verdict already shows our single-window conic close is not V2-stable in
the circular-coplanar moon-ephemeris model.)

## Follow-ons (scoped, NOT done this run)

### #484 — structured citation keys (replace free-text citations)

**Problem.** Citations live as free text in code comments, note prose, and
catalogue `first_published`/`corroborating_sources` blocks. Free text is where
concept-collision hallucinations hide and propagate, and the #485 ratchet can
only check the author-consistent subset.

**Design sketch.** Introduce a single **citation registry** keyed by a stable id
(e.g. `jones-2017-vem-577`, `hernandez-2017-ieg-462`), each entry carrying the
grounded `{authors, year, title, venue, doi/url, system, bodies}` (the #483
surfaces, promoted to first-class). Then:
1. Replace free-text citations in `KNOWN_CORPUS`, catalogue blocks, and note
   prose with **references to keys** (`cite: jones-2017-vem-577`).
2. A loader resolves keys → grounded entries; a lint asserts every `cite:` key
   exists and every catalogue row's bodies/system are contained by its key's.
3. The #485 ratchet then runs over **every** keyed citation (not just
   author-consistent strong-links), closing the coverage caveat above.
Migration is incremental: seed the registry from the existing grounded
`KNOWN_CORPUS` anchors (already `{authors, system, body_set, citation, doi}`),
then backfill catalogue/notes by key.

### #486 — provenance tags: verified-vs-inherited

**Problem.** We cannot today tell, for a given citation, whether it was
**verified** (someone opened the source and confirmed bodies/system) or
**inherited** (copied from our own notes/code as fact — the propagation
vector). The memory rule "distrust inherited citations" has no machine surface.

**Design sketch.** Add a per-citation `provenance` tag with two values:
`verified` (grounded against the source's title/abstract or on-disk PDF, with a
`grounded_on` date + a `grounded_against` pointer to the OCR/text/PDF) and
`inherited` (present but not yet ground-truthed). New citations default to
`inherited`; a citation is promoted to `verified` only by an explicit grounding
step. A lint can then (a) refuse to let an `inherited` citation back a V-tier
promotion or a golden EXPECTED value, and (b) report the `inherited` backlog for
triage. Pairs with #484 (the registry is where the tag lives) and with
`feedback_ground_citations_against_content` (makes "ground before cite" a
checkable gate, not a discipline). The #483 grounding comments added this run
(IEG, VEM, Restrepo-Russell) are the first `verified` entries by example.

## Artifacts

* `src/cyclerfinder/search/literature_check.py` — #483 `system` field +
  `system_grounded` + grounded IEG/VEM/Restrepo-Russell anchors.
* `tests/data/test_citation_integrity.py` — #485 ratchet (57 tests, green).
* `docs/notes/2026-06-26-470-moontour-admission-wave-v2-verdict.md` — grounding
  addendum.

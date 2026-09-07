# #875 — Novelty policy decision (owner, 2026-09-07)

**Status:** DECIDED, written into `docs/spec.md` §16.4 (new "Novelty policy" bullet) and §16.5
(attribution fields). Registered by `#864` (`docs/notes/2026-09-05-864-project-feasibility-future-review.md`,
sec. 10) as the load-bearing owner decision of the Sep-Nov 2026 roadmap. Taken in chat on 2026-09-07;
this note is the decision record.

## The question

Three kinds of result the roadmap can produce sit between `known-class-member` and `candidate-novel`
in the §16.4 vocabulary. The project's own precedents (`#577` 0/36 Galilean ruling, 2026-07-12; `#817`'s
reading of it, 2026-08-10) leaned toward "known-class" for all three. Nothing in the spec distinguished
"the source authors excluded this" from "the source authors did not enumerate this", or "a theorem
guarantees objects of this type exist" from "this object is published".

| # | Case | Concrete instance | Precedent before this decision |
|---|---|---|---|
| (i) | Member of a class the source authors EXCLUDED on stated practicality grounds | Jones-Hernandez-Jesick 2017 searched VEM triple cyclers at 1-2 synodic periods, ≤6 flybys; 3-/5-synodic and >6-flyby classes excluded by design (`#867`'s target) | none explicit; `#577` by analogy |
| (ii) | Known architecture at a system the authors never treated | Russell-Strange 2009 one-working-node moon cycler at Uranus / Neptune (`#870`, `#874`) | `#817`: R-S-class member at a new pair is known-class |
| (iii) | Theorem-generic object | Two symmetric periodic orbits accumulating on `#781`'s Neptune-Triton homoclinic of the Miceli-Bosanac 4:5 saddle (`#868`) | torus rows: `known-class-member` because KAM guarantees existence |

## The decision

1. **(i) Author-excluded classes are `candidate-novel`.** The row attributes the source method
   (`corroborating_sources`) and quotes the exclusion (`notes`). Owner: "yes this is novel, it
   should attribute but it's genuinely new."
2. **(ii) Known architecture at a never-treated system is `candidate-novel`.** The row attributes the
   architecture's source and `notes` must explain HOW it was re-applied (system, flyby body, passive
   target, what changed physically). Owner: "yes this is novel, it should attribute and explain how
   it was reapplied." This supersedes `#817`'s reading for new-SYSTEM cases. The `#577` 0/36 ruling
   stands: those were same-system members of Russell-Strange's own Jovian class.
3. **(iii) Theorem-generic objects are `candidate-novel` if not a member of any published family.**
   Owner asked "if it's published then it's not novel?" — correct, and it applies at two levels:
   literal collision (the orbit or its family is in the source's tabulated members) ⇒
   `known-reproduction`/`known-class-member`; otherwise a theorem asserting that orbits of a type
   exist is a scope statement, not a computation — it gives neither the orbit, its period, nor its
   encounter geometry (inside/outside the moon's SOI, i.e. whether it is cycler-class at all). The
   row cites the theorem and attributes the parent object's source; wording is "first computed",
   never "first predicted". Torus rows around published orbits stay `known-class-member` (they
   characterise a neighbourhood of a known object, not a new orbit). Owner: "agree with your
   recommendation for iii."
4. **Literal collision always comes first.** Matcher, family-membership check (incl. bifurcated
   sub-families), and `literature_check.py` all run before any of (i)-(iii) applies.

## What is published vs ours in case (iii), for the record

- **Published:** the Neptune-Triton 4:5 saddle resonant periodic orbit itself (Miceli & Bosanac 2026,
  JAS 73:11, DOI 10.1007/s40295-025-00545-z, supplementary row; 63 other planar POs in the same
  files). That paper builds tours from motion primitives; it computes neither homoclinic connections
  nor cycler-class orbits.
- **Ours, adjudicated novel 2026-08-08 (`#781`):** two on-axis homoclinic connections of that saddle.
- **Ours, closed 2026-09-05 (`#864` refutation pass; `#868` to write back):** two symmetric periodic
  orbits near those homoclinics — PRIMARY x=1.16933872, ydot0_sign=-1, hc=5, T=68.7485 (2.26 base
  periods), closest Triton approach 0.58x SOI; SECONDARY x=-1.38561105, ydot0_sign=+1, hc=4, T=86.898
  (2.86 base periods), 0.46x SOI. Cycler-class under the `#855` rule.
- **Still required before the tag:** the family-membership (collision) check against all 64 Miceli
  rows and their bifurcation sub-families; reading Campagnola 2014 and AIAA 2024-1280; a live
  `literature_check.py` run.

## Consequences for the roadmap

- `#867` (Jones excluded VEM classes) and `#870` (Uranian one-working-node) are novelty-bearing, not
  census-only. The `#864` note's "~1 in 10" Uranian novelty figure and its "policy-gated" hedges on
  items 3, 6, 7, 10 are superseded; the V-tier and census value of those items is unchanged.
- `#868` writes the two Neptune-Triton rows at V1 as `candidate-novel` ONLY after the collision check
  and literature check pass; otherwise `known-class-member`/`known-reproduction` with the citation.
- `#865` applies the vocabulary to the 8 `source: discovered` rows that carry no `our_status` today.
- The website `/about/` "our_status" paragraph now states the three cases in plain language
  (cyclers.space repo, same day).
- `#577`, `#578`, `#817` bullets are NOT rewritten; this note and spec §16.4 are the superseding
  record, cited from the `#875` bullet.

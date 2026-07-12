# #579 — literature_check.py Antoniadou anchor mislabel fix + #287/#301 re-audit

Date: 2026-07-12
Issue: #579 (P1, corpus-accuracy, do first per Fable's priority order)

## 1. The bug

`src/cyclerfinder/search/literature_check.py` (~line 1522-1546) carried a
`CorpusAnchor` named `"Antoniadou-Voyatzis spatial resonant periodic orbits
in CR3BP (2018)"`, `authors=("Antoniadou", "Voyatzis")`, citing
arXiv:1811.09442 with `doi=None`. That arXiv id is actually **Antoniadou &
LIBERT** 2019 (MNRAS 483(3):2923-2940, DOI 10.1093/mnras/sty3195) — Voyatzis
is not an author. The correctly-cited twin anchor already existed in
`src/cyclerfinder/genome/known_corpus_3d.py` (fixed by #459 / commit
`0b1528f`, 2026-06-25); the `literature_check.py` copy was never updated.

The bad anchor's scope comment additionally claimed coverage of "1:1" MMR
resonance and was cited by name as "the anchor for #287's 3D Braik-Ross
(1,1) family extension (likely rediscovery target)." The paper actually
behind arXiv:1811.09442 (Antoniadou & Libert 2019) covers MMRs 3/2, 2/1,
5/2, 3/1, 4/1, 5/1 at mu=0.001 (Jupiter-mass planetary, no Earth-Moon ICs)
and explicitly excludes asymmetric/isolated spatial families — it does
**not** cover 1:1 resonance at all.

## 2. Fix applied

Relabel-only (no second anchor added — see §3 for why). In
`literature_check.py`:

- `name` → `"Antoniadou & Libert spatial resonant periodic orbits in the
  RTBP (2019)"`
- `authors` → `("Antoniadou", "Libert")`
- `citation` → full MNRAS citation with the correct DOI, explicitly stating
  the MMR list and the "does NOT cover 1:1" fact, and explicitly stating
  the anchor does not on its own establish #287's (1,1) family as a
  rediscovery.
- `doi` → `"10.1093/mnras/sty3195"` (was `None`).
- Scope comment (`period_band_tu` block) rewritten: MMR list corrected to
  3/2, 2/1, 5/2, 3/1, 4/1, 5/1 (NO 1:1); explicit note that the matcher's
  structural fields (`primary`, `body_set`, `topology_label`,
  `period_band_tu`) cannot express MMR-level distinctions, so a 1:1
  candidate is *not* mechanically excluded by this anchor even though it
  is out of the paper's real scope — only the comment records that. Also
  records that no anchor in the corpus currently covers restricted-problem
  (CR3BP) Earth-Moon 1:1 spatial resonance.
- The two doc-comment cross-references at `literature_check.py` ~line 107
  and ~313 (`CandidateSignature.period_band_tu` / `CorpusAnchor.period_band_tu`
  docstrings, both naming "#301 ... Antoniadou-Voyatzis 2018 anchor") were
  updated to the corrected name.

`docs/notes/CORPUS_INDEX.md` and the 2013-paper digest were already
self-corrected on 2026-07-12 per the task's own note (one-anchor-not-two
finding); no further digest edit was needed.

## 3. Why the 2014 IAU proceedings paper was NOT used to add 1:1 coverage

The user located and filed **Antoniadou, Voyatzis & Varvoglis (2014), "1/1
resonant periodic orbits in three dimensional planetary systems," Proc. IAU
9:82-83, DOI 10.1017/S1743921314007893** — topically the closest possible
match to the mislabeled anchor's "1:1" claim (same first author, "1/1
resonance," "three dimensional planetary systems" in the title).

It is filed and digested
(`docs/notes/2026-07-12-digest-antoniadou-voyatzis-varvoglis-2014-1-1-resonant.md`,
`docs/notes/CORPUS_INDEX.md` line ~345,
`antoniadou-voyatzis-varvoglis-2014-1-1-resonant-periodic-orbits-3d-planetary-systems-iau-proc-310-82-doi-10.1017-S1743921314007893.pdf`
in the private corpus). Read in full (it is a 2-page proceedings note).

**It does not validly fix the anchor.** Its own text states "we utilize the
spatial general TBP in a rotating frame of reference" — this is the
**general** three-body problem (both bodies massive), the same model class
as its 2013 sibling paper (`antoniadou-voyatzis-2013-2-1-resonant-...pdf`,
already digested and already determined "not a codebase capability fit").
The buggy anchor's title claims coverage "in the **Restricted** Three-Body
Problem" — the model class cyclerfinder's actual candidates (one massless
third body) live in. A general-TBP paper's 1:1-resonance families (two
massive planets, ρ = m2/m1 < 0.0205 for any spatial bifurcation to exist at
all) are not a substitute literature anchor for a restricted-problem
Earth-Moon CR3BP candidate — citing it that way would be exactly the kind
of model-class mismatch this task was opened to catch. So the anchor stays
a single relabel with no 1:1 coverage claim, and the 2014 paper is recorded
here as a *checked-and-rejected* candidate, not silently ignored.

The paper it itself cites for the genuinely restricted-problem-adjacent
work, Antoniadou & Voyatzis (2014, Ap&SS 349:657, DOI
10.1007/s10509-013-1679-8), is not yet acquired; if a real restricted-
problem 1:1-resonance reference is ever needed, that citation chain remains
open (not chased further here — out of #579's scope).

**Net result: no anchor in the corpus (literature_check.py KNOWN_CORPUS or
genome/known_corpus_3d.py KNOWN_CORPUS_3D) covers Earth-Moon (or any
restricted-three-body) 1:1 spatial resonant periodic orbits.** This is an
honest open literature gap, not papered over.

## 4. Re-audit: #287 Braik-Ross (1,1) "likely rediscovery" verdict

Traced the chain: `docs/notes/2026-06-16-287-3d-aldrin-scoping-spike.md` →
`docs/notes/2026-06-16-299-3d-lit-check-bifurcation.md` →
`docs/notes/2026-06-16-301-3d-subfamily-validation.md`.

### #287 itself (the spike)

The spike doc is already careful: **"Verdict: spike SUCCEEDED... likely
rediscovery... not a novelty claim... pending the full literature-check
chain."** It explicitly defers to a follow-up literature check and records
"NO catalogue writeback. NO novelty claim." No fix needed here — the spike
never asserted a completed rediscovery verdict.

### #299 Part A (the actual "clean rediscovery" verdict)

This is where the bug bites. #299 ran all 265 members of the base (1,1)
family (`data/family_296_3d_em_11.jsonl`, T_TU ≈ 9.3-10.3, inside the
anchor's period_band) through a literature check and reported:

> All 265 members hit the **Antoniadou-Voyatzis 2018** anchor with
> confidence 0.85 ... **Verdict: clean rediscovery ... No novelty claims
> warranted.**

That verdict leaned on the (now-corrected) anchor's claimed "1:1, 2:1, 3:2"
scope — a claim that was never true of the paper actually behind
arXiv:1811.09442 (Antoniadou & Libert 2019 does not cover 1:1 at all). The
matcher itself (`_candidate_anchors` in `literature_check.py`) only filters
on `primary` / `body_set` / `topology_label` / `period_band_tu` — it has no
field that discriminates MMR ratio, so relabeling the anchor does not by
itself change what the matcher would flag; it only removes the false
*documentation* claiming the underlying paper covers this topology. The
substance (does a real published corpus entry cover Earth-Moon 1:1 spatial
resonance) was never actually checked against a real anchor's real content.

**Independent corroboration that this verdict was wrong:** the later,
more careful `genome/known_corpus_3d.py` module (#434, 2026-06-16, later
fixed 2026-06-25 by #459) states in its own module docstring:

> "The planar (1,1) Braik-Ross root is `k_z = 0`; its #287 3D extension is
> the **novel-frontier target the campaign is testing for, NOT one of
> these anchors**."

That is the corpus's own later, correctly-cited authorship concluding the
opposite of #299's verdict — the #287 (1,1) 3D family is explicitly framed
as *not* covered by the (correctly cited) Antoniadou-Libert 2019 spatial
anchor.

**Conclusion: #299 Part A's "clean rediscovery, no novelty claims
warranted" verdict for the base 265-member (1,1) family is a false
NOT-novel, and is hereby REOPENED.** Its literature status is honestly
**inconclusive / open**, not "published." Nothing was ever admitted to
`data/catalogue.yaml` from this lineage (#299/#301 both explicitly recorded
"NO catalogue writeback" pending V0-V5 gauntlet adaptation for 3D orbits,
which never landed), so there is no catalogue row to correct or retract.
This reopening is scoped to a future task (a real literature search for a
restricted-problem Earth-Moon or mu-general 1:1 spatial resonant orbit
catalog, or acquisition of the still-missing Antoniadou-Voyatzis 2014
Ap&SS companion paper's restricted-problem-adjacent content) — not
attempted here, out of #579's scope.

### #301 Part (the k=3..6 Neimark-Sacker sub-families)

**Unaffected by this bug, and its "0/145 novelty-claimable" verdict
stands.** #301 explicitly and *correctly* EXCLUDED the AV/Antoniadou anchor
for all four sub-families via the `period_band_tu` filter (their T_TU ∈
[20,44] sits outside the anchor's [0,15] band regardless of which paper is
behind it). The actual "likely-rediscovery" call for those 145 members
rests on **other** anchors — the Braik-Ross / Roberts-Tsoukkas-Ross /
Kumar-Rosengren-Ross "Earth-Moon CR3BP family network" framework papers —
not on the mislabeled anchor. Those anchors are unaffected by this fix.
No re-audit action needed for #301's own verdict.

### Comments elsewhere still saying "Antoniadou-Voyatzis 2018"

`src/cyclerfinder/search/cr3bp_general_periodic_3d.py`,
`cr3bp_3d_family_tracer.py`, and `src/cyclerfinder/genome/qp_tori.py` (plus
`tests/genome/test_qp_tori.py`) still narrate "Antoniadou-Voyatzis 2018" in
prose docstrings/comments as shorthand for the anchor these #287/#299/#301
lineage docs point at. These are historical narration of what was believed
at the time (matching the doc trail above), not anchor definitions, and are
left as-is — the substantive fix (the anchor itself, the two
`literature_check.py` docstrings that actively describe matcher behaviour,
and this note) is what future readers need to find the corrected story.
`genome/spatial_novelty_prefilter.py` already cites Antoniadou & Libert
2019 correctly (post-#459) and needed no change.

## 5. Summary

| Item | Before | After |
|---|---|---|
| Anchor name | "Antoniadou-Voyatzis spatial resonant periodic orbits in CR3BP (2018)" | "Antoniadou & Libert spatial resonant periodic orbits in the RTBP (2019)" |
| Authors | (Antoniadou, Voyatzis) | (Antoniadou, Libert) |
| DOI | None | 10.1093/mnras/sty3195 |
| Scope comment "1:1" claim | present (wrong) | removed; MMR list corrected to 3/2, 2/1, 5/2, 3/1, 4/1, 5/1 |
| 2014 IAU 1:1 paper | not considered | acquired, read, and explicitly rejected as a fix (general-TBP model-class mismatch) |
| #287/#299 (1,1) base-family verdict | "clean rediscovery, no novelty claims warranted" | REOPENED as inconclusive/open — no real anchor covers this topology; no catalogue row affected (none was ever written back) |
| #301 sub-family verdict (k=3-6) | "0/145 novelty-claimable" | unchanged — rests on other, unaffected anchors |

# #1052 Rule U (untopologied signatures) and the Aldrin single-body fix — PRE-REGISTRATION (2026-10-11), committed before any code

**This is a LEAD RULING ON GATE SEMANTICS (2026-10-11), to be reported to the owner in the next
question batch.** The lead approved both changes "because both move only in the conservative
direction (published -> inconclusive, never toward not-found), and inconclusive is the owner's
decision point".

## Problem (the #1045 A1 review, finding 1; #1052 entry)

A signature with NO topology label reads "published" through anchors that are not cycler sources.
Baseline at HEAD (d4ff1a80), all untopologied and unlabelled:

| Signature | Today | Anchor |
|---|---|---|
| Saturn (Titan,) | published 0.95 | Davis-Phillips-McCarthy (Titan tulip / Enceladus NRHO orbiters) |
| Neptune (Triton,) | published 0.95 | Miceli-Bosanac (motion primitives) |
| Mars (Phobos,) | published 0.95 | Wallace (Mars-Phobos rendezvous) |
| Sun (Jupiter,) | published 0.95 | Koon-Lo-Marsden-Ross (no label) |
| Uranus (Umbriel, Titania) | published 0.95 | Pergola (manifold transfers) |
| Sun (E,), (M,), (E, E) | published 0.85 | Aldrin Earth-Mars cycler (one-body subsets) |
| Pluto (Charon, Nix) | published | Howett (no label); the Pluto lane |

The untopologied production paths (`signature_from_review_entry`, `saturn_uranus_campaign.py`)
produce repeated-moon cycler candidates.

Why relabelling cannot close it:
- G1' caps only anchors with a tour label and no cycler-class label.
- Davis-Phillips-McCarthy's sourced labels {tulip, halo, nrho} are cycler classes.
- Wallace and KLMR have no label. Adding one activates the #349 topology filter for TOPOLOGIED
  signatures, which can move them toward not-found.

## Rule U (untopologied signatures only)

For a signature whose `topology_label` is EMPTY, a hit synthesised from an anchor can make it
"published" only if that anchor's `topology_label` contains "repeated-moon", the cycler class.
- The hit of every other anchor is scored and then capped at `TOUR_ONLY_CAP` (below
  `MATCH_THRESHOLD`), so it gives at most "inconclusive".
- "Every other anchor" means: no label at all, or only resonant / halo / nrho / tulip / transfer /
  tour / ephemeris / quasi-satellite / retrograde-satellite / binary-coorbital / model-transition
  labels.
- When the best hit is capped, the notes name the F7 anchors as before (#1045 A1-2).
- Signatures WITH a topology label are untouched (the #349 filter governs them).
- Real web hits (no anchor identity) are untouched.
- G1' becomes its special case. Every anchor G1' caps (a tour label and no cycler-class label) lacks
  "repeated-moon", so rule U caps it too.

## The Aldrin single-body fix (the two parts land together or not at all)

(a) Tag `n_bodies_scope = 2` on the Earth-Mars cycler anchors at the Sun whose body set is {E, M} and
whose topology contains "repeated-moon":
- Aldrin;
- Russell-Ocampo / McConaghy SnLm;
- Genova-Aldrin;
- the three McConaghy anchors;
- Pisarevsky-Kogan-Guelman;
- Rall-Hollister;
- Rall;
- Rogers-Hughes-Longuski-Aldrin.

Source: each is an Earth-Mars cycler or periodic-orbit family with encounters at both planets (the
anchors' own citations and digests).

(b) EXTEND F7 (`_different_architecture_same_system`): it also fires on the "n-bodies" exclusion
reason when the signature visits a SINGLE distinct body, at body-set inclusion (the signature's body
is in the anchor's body set). The result is "inconclusive", naming the anchor.

**Without (b), the tag (a) alone would create a false-fresh risk.** It would exclude the anchors for
(E,) / (M,) / (E, E), and with no other match those signatures could fall to "not-found". That is
why (a) and (b) land together or not at all.

The extension is restricted to SINGLE-body signatures. Two-body signatures that the n-bodies tag
excludes from a three-body anchor (e.g. ev-A's (E, V) against Jones VEM) are NOT escalated, so the
ev candidates keep "not-found".

## Expected results (written before running; stop rule below)

Positives: rows that read "published" today via a repeated-moon anchor, untopologied, with the
working label the source allows. Each must stay "published".

| Primary | Positive signature | Today |
|---|---|---|
| Sun (E-M) | (E, M), unlabelled (aldrin-classic-em-k1-outbound) | published, Aldrin (Byrnes-Longuski-Aldrin) |
| Earth | ("Moon",), resonance 7:3, working "one" (Casoliva 7-3b, as #1047/#1062) | published, Casoliva |
| Jupiter | (Ganymede, Callisto, Ganymede), working "one" (R-S GanCal#1) | published, R-S |
| Saturn | (Titan, Enceladus, Titan), working "one" (R-S TitEnc#631) | published, R-S |

**Neptune, Mars, Uranus and Sun-Jupiter have NO repeated-moon anchor in the corpus, so they have
no positive.** Under rule U every untopologied signature there reads at most "inconclusive". That is
the intended outcome, not a failure.

Negatives:
- Jupiter (Io, Callisto), working "one", topology {"repeated-moon"}: known-architecture-new-system,
  unchanged.
- The pump-tour negative (`test_1045` G1 tour case, Uranus {"mga-tour"}): unchanged.

The problem cases, expected after the change:
- Saturn (Titan,), Neptune (Triton,), Mars (Phobos,), Sun (Jupiter,) and Uranus (Umbriel, Titania)
  -> inconclusive.
- Sun (E,), (M,) and (E, E) -> inconclusive, naming an Earth-Mars cycler anchor.
- Pluto (Charon, Nix) -> inconclusive (reported, not judged; the Pluto lane needs no relabel).

Byte-identical (topologied signatures only):
- `data/942_943_litcheck_scope.json` (26 controls, the negative, 5 candidates);
- `data/1025_litcheck_973_members.json` (42 members);
- `data/1047_literature/gate_results.json` (18 Earth-Moon entries and their controls).

Probe (`scripts/probe_1052_rule_u.py`, HEAD baseline module against the working tree):
- Signatures: every catalogue row with a sequence (primary inferred, unlabelled) and every 1-4-body
  subset (with repeat forms) of every anchor's body set.
- Topologies: none, repeated-moon, resonant, halo, nrho, tulip, mga-tour, pump-tour, transfer,
  quasi-satellite, binary-coorbital.
- **STOP RULE:** no signature in ANY topology may move toward "not-found" or
  "known-architecture-new-system".
  - Topologied signatures must not change at all.
  - Untopologied ones may move only published -> inconclusive.

Tests:
- Pinned tests (tests/search/, < 2 s): the positives, the problem cases and the single-body
  extension, including an (E, V) two-body case unchanged.
- The F13 ratchet and all the literature tests green; ruff and mypy clean.

## Also registered

The remaining #1052 items (ratchet evasions: a Call with a status argument; NamedExpr; nits) become
backlog #1066 (lead ruling).

# Digest: Howell & Pernicka 1988, "Numerical Determination of Lissajous Trajectories in the Restricted Three-Body Problem"

**Source:** K.C. Howell & H.J. Pernicka, *Celestial Mechanics* 41:107-124 (1987, published 1988).
DOI `10.1007/BF01238756`. Purdue University. Filed at
`cyclers_pdf/papers/howell-pernicka-1988-numerical-determination-lissajous-trajectories-celest-mech-41-107-doi-10.1007-BF01238756.pdf`
(md5 `b7ae7fcc7777d7eb28406aab46260518`), text-layer (image-scan-derived but OCR'd with a real text
layer, 35,883 chars — well above the 10 char/page floor, no re-OCR needed), 18 pages.

**Acquisition context:** `#730` backlog item 97 (registered 2026-08-07), flagged as the deepest
technical source cited by Parker, Davis & Born 2010's own corrector description (item 96, ref [21]
in that paper) — itself the technique now being applied by `#782`'s reopened Saturn-Titan chain
closure. User-supplied upload.

## What this paper actually contains — full corrector algorithm, worked in detail

This is the genuine, fully-derived multi-patchpoint differential-correction algorithm the later
papers (Parker/Davis/Born 2010, and by extension `#782`'s own current work) cite but don't re-derive.
**Two-level iterative structure:**

1. **Level 1 — position continuity.** The trajectory is divided into segments (target points at the
   start/end of each interval — roughly half-revolution spacing in the paper's own Lissajous
   application). Each segment is corrected independently via a **state-transition-matrix-based linear
   correction**: for segment `o->p`, propagate `X_o` (with STM) to get `p*`; the position mismatch
   `(δx_p, δy_p, δz_p)` relates to the velocity/time adjustment via
   `δX_p ≈ Φ(t_p*,t_o) δX_o + (∂X/∂t)|_p* δ(t_p*-t_o)` (their Eq. 7) — a 3-equation, 4-unknown
   underdetermined linear system (`δẋ_o, δẏ_o, δż_o, δ(t_p*-t_o)`), solved via the **minimum-Euclidean-
   norm pseudoinverse solution** `ξ = L^T(LL^T)^{-1}β` (Eq. 8-9) rather than a square Newton step.
   Segments are corrected one at a time this way until each is individually continuous in position.
2. **Level 2 — velocity-discontinuity (Δv) reduction at the patch points.** Once every segment is
   individually position-continuous, a finite velocity discontinuity Δv remains at each patch point
   (segments were corrected independently, so segment-two's arrival velocity at `p` generally
   disagrees with segment-one's). This level treats the CHANGE in each patch point's Δv as a function
   of small changes in ALL target-point positions and ALL segment times simultaneously (STM-derived
   sensitivities propagated backward/forward from each patch point), assembled into one linear system
   solved to drive every Δv toward zero at once — this is the actual "simultaneous" multi-point
   correction, not sequential per-segment correction.
3. The outer loop repeats: re-integrate all segments at the newly-adjusted target positions/times
   (Level 1 again), re-assess remaining Δv's (Level 2 again), until all Δv's fall below tolerance.
4. **Seeding**: first guesses for target-point states come from an analytic (3rd-order or lower)
   Lissajous approximation in-family, or — notably, directly relevant to `#782`'s own out-of-family
   chain-orbit problem — **a continuation method when beyond the analytic approximation's validity
   range** ("a continuation method was successfully employed to obtain first guesses for the target
   point state vectors" when analytic guesses failed).

**Relation to this project's own existing `cr3bp_multiple_shooting.correct_multiple_shooting`**: not
read/cross-checked in this digest pass (out of scope for a corpus filing task) — `#782`'s own agent
should make this comparison directly if it wants to sharpen its corrector, since this paper's
minimum-norm-pseudoinverse Level-1 step + separate Level-2 simultaneous Δv-reduction is a specific,
fully-specified algorithm that may or may not match what this project's existing multiple-shooting
utility already implements.

## Registration

Filed in `cyclers_pdf`, `CORPUS_INDEX.md` and `#730` backlog item 97 updated to ACQUIRED in the same
session. No catalogue/code changes. Directly relevant to `#782` (in progress at time of this digest) —
flagged to that task directly, not left for it to discover independently.

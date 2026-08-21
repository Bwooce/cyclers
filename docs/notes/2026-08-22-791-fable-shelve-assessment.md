# Fable's assessment of `#791`'s path forward (after both Stage 0 screens failed)

**Requested 2026-08-22** after `#791`'s Stage 0 checkpoint (pairwise Tisserand screen, then a
user-directed chained/energy-consistency follow-up) both found the Galilean moon-tour sequence
alphabet has a universal "hub" flyby energy (V∞=4-5 km/s links all 12 possible moon pairs at
once), so neither energy-based filter can discriminate at all
(`docs/notes/2026-08-22-791-stage0-tisserand-pruning-checkpoint.md`).

## Recommendation: (c) shelve `#791` as scoped, redirect to `#790`

Options on the table were (a) accept Stage 1's full un-narrowed cost (~10^6-10^7 cells), (b) build
a non-energy (phasing/resonance) filter, (c) shelve and redirect to `#790`.

**(a) rejected**: the hub-energy finding is not just a failed filter — it predicts Stage 1's own
downstream sequence-ranking mechanism would be noise-dominated even if funded. `joint_cell.py`'s
evaluator discriminates PER-CELL (`#501`'s 213/3072≈7% feasible), but Stage 1's actual decision
output is a PER-SEQUENCE ranking ("shoot the top decile"). Since every sequence in this alphabet
shares the same energetic character, whatever would separate sequences must come from phasing
statistics alone at 512-4k Sobol samples in 13-25 dimensions — too sparse to trust a ranking by.
Paying the full ~300-2800 CPU-hours would produce an untrustworthy ranking feeding a Stage 2 whose
novelty ceiling was already near zero (`#858`'s 0/36 Galilean precedent).

**(b) rejected**: sanity-checked concretely and it collapses into one of two failure modes. The
truly-cheap algebraic version (period commensurability) has almost no discriminating power for a
mirror-image reason to the energy screen: Io/Europa/Ganymede sit in the Laplace resonance
(n_Io − 3n_Eu + 2n_Ga = 0), so their relative geometry is already periodic and any commensurate tour
period is near-commensurate with all three — this reduces to roughly a 1-bit "does the sequence
contain Callisto" filter on 3269 sequences. The honest version (per-leg feasible-ToF windows chained
with accumulated moon phase) IS a coarse patched-conic evaluation — i.e. it's Stage 1 wearing a
Stage 0 costume, not a cheap gate, and would take days to build, duplicating `joint_cell.py`.

**(c) recommended**: shelve `#791` as scoped — the Stage 0 STOP rule fired twice on independently
verified physics, the same disposition logic `#861` used for Resonant Atlas. `#790` (itinerary
enumeration) was already ranked ahead by `#858` on expected value before today's result, and today
widens that gap further since `#791` lost the one mechanism (cheap pruning) that could have made it
cost-competitive. `#790`'s path is `#858` Sec. 4.3's own prerequisite: the Jacobi-constrained
multiple-shooting corrector (the `#687` extension) + a symmetric-only mini-pilot.

## Optional small follow-up (registered separately as `#863`, NOT this campaign's next funded slice)

A "density inversion" experiment `#501` never ran: pick ~5-10 sequences by judgment (short
inner-three Laplace-resonance tours, plus one or two literature-adjacent Callisto patterns) and run
*deep* per-sequence density — 32k+ cells each, `#501`'s Liang Member D positive control re-run
first, Stage 2 n-body shooting now legitimately unblocked since `#849` landed. This attacks the
actual identified flaw (`#858` Sec. 5.1's "sparse depth", not breadth) and would produce
empty-region stamps with real conditional meaning, unlike a stamp from an already-judged-too-sparse
method. Costs ~1-3 days. Still ranks below `#790` on expected value — idle-machine-time task, not a
priority dispatch.

## Disposition

`#791` SHELVED as scoped (mirrors `#861`/Resonant Atlas). `#863` registered as the optional
concentrated-density idle-time follow-up. Effort redirect recommendation: `#790` via `#858` Sec.
4.3's corrector-build path. No catalogue writeback; no code changes from this assessment itself.

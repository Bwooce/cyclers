# Roberts-Tsoukkas & Ross 2026 journal — "Stable Prograde Earth-Moon Multi-Orbiter Cyclers via Three-Body Dynamics" (discovery-capability digest)

> **SUPERSEDED 2026-06-13 (#251).** This note was written under a WebFetch
> DENIAL and ASSUMED the PDF held transcribable higher-μ (μ, C, T) tuples and
> printed ICs. WebFetch was later re-enabled and the full PDF was fetched: it
> contains **NO numeric data at all** (no tables, C, T, σ, or ICs) — the
> binary-star families exist only as figures. The "USER FETCH the PDF" follow-on
> is DISCHARGED (fetched; nothing to transcribe) and the "ADOPT gated on fetch"
> verdict collapses to **no adoption possible from this paper** (qualitative
> priors only). See `docs/notes/2026-06-13-roberts-tsoukkas-2026-multi-orbiter-mining.md`
> for the full-text verdict. The HIT-by-coverage hopes below were aspirational
> and are NOT borne out by the manuscript.

**Date:** 2026-06-13
**Task:** #250 discovery-capability sweep (head-start candidate 4).
**Source (free):** M. Roberts-Tsoukkas & S. D. Ross, "Stable Prograde
Earth-Moon Multi-Orbiter Cyclers via Three-Body Dynamics," journal manuscript
hosted at the Virginia Space Grant Consortium
(https://vsgc.odu.edu/wp-content/uploads/2026/04/Roberts-Tsoukkas_Michael_Cycler-Journal-Paper.pdf,
posted Apr 2026). This is the journal extension of **AAS 25-621**, already
mined and ADOPTED (5 `ross-rt-em-cycler-*-2025` V1 rows; see
`docs/notes/2026-06-11-ross-roberts-tsoukkas-2025-mining.md`).
**Writeback: NONE** (digest only; sourced-golden discipline still applies to
anything later wired).
**Access caveat:** read via WebSearch metadata + abstract only — `WebFetch`
was DENIED in this environment and `WebSearch` cannot extract the PDF's
numeric tables. The new-member tuples below are described from the abstract /
summary text, NOT transcribed digit-for-digit. **This note flags the paper for
a user PDF fetch (#116) so the new (μ, C, T) tuples can be transcribed under
the same rescan discipline applied to the AAS version.**

---

## Verdict: HIT (topology/coverage extension) — but ADOPT is gated on a user PDF fetch

The AAS 25-621 version gave us 5 Earth-Moon (k1,k2) families at the physical
lunar mass ratio (μ = 1.2150584e-2). The journal version's stated new
contributions (abstract + project summary, multiple independent search hits):

1. **Extension to HIGHER mass parameter, up to the binary-star range.** "The
   study aimed to expand previous results on multi-body cyclers to systems with
   higher values of the mass parameter." "Several cycler trajectories were
   discovered in [the] binary-star-system mass-parameter range." This is a
   genuine **coverage extension** — cycler families at μ values our catalogue
   has zero rows for (everything we hold is Earth-Moon-scale or smaller). It is
   *not* a new genome (still the symmetric fixed-Jacobi PCR3BP corrector +
   arc-length continuation of AAS 25-621), but it populates a region of
   (μ, C, T) space the project has never catalogued.
2. **"Universal existence of stable subfamilies."** "Substantial evidence was
   discovered in favor of the existence of a stable subfamily for every cycler
   family." The AAS version observed (p.14) that every (k1,k2) family it
   *examined* had a stable region near C^max; the journal elevates this to a
   near-universal structural claim across mass parameters. If true, it is a
   strong **prioritizer prior**: any (k1,k2) topology, at any μ, is predicted to
   contain a stable, ballistic, catalogue-grade member near its C^max — i.e. a
   guaranteed-non-empty region for the discovery daemon to seed.
3. **"Resonant and highly maneuverable regimes."** The summary mentions the
   systematic family-design method "revealing resonant and highly maneuverable
   regimes that offer new opportunities for cislunar mission design" — likely a
   fuller treatment of the synodic-resonance structure the AAS version only
   touched for (1,1)/(3,2) (2:3 and 2:5 synodic).
4. **"Multi-orbiter" framing.** Re-brands (k1,k2)-cyclers as multi-orbiter
   cyclers (k1 Earth orbits + k2 Moon orbits per period). Possibly more (k1,k2)
   combinations than the AAS version's five ("we explored several, but certainly
   not all," AAS p.14) — to be confirmed against the PDF tables.

## Discovery vs anchor

**Both, leaning discovery-via-coverage.** The Earth-Moon members duplicate /
supersede our 5 adopted V1 rows (anchor — diff for any digit revisions or the
two flagged AAS table self-inconsistencies C_(2,1)/C_(3,1) being resolved). The
**higher-μ / binary-star members are net-new catalogue territory** — a coverage
extension into a μ regime the project has never represented, found by a method
we have already validated (the Ross corrector, ADOPTED 2026-06-12, see
`docs/notes/2026-06-13-ross-adoption-results.md`). That is the discovery payoff:
our toolkit can *immediately* reproduce these members (same model, same
corrector) the moment the (μ, C, T) tuples are transcribed.

## Reproducible data we could adopt

- **New (μ, C^stable, T^stable) tuples** at higher mass parameters — the
  catalogue-grade columns, recoverable to a full state via the same 1-D
  symmetric solve used for the AAS rows (IC = (x0, 0, 0, ẏ0(x0,C))). **Cannot
  transcribe without the PDF** (WebFetch denied).
- Resolution of the AAS Tables-3/4 C_(2,1) (Δ 2.23e-6) and C_(3,1) (Δ 7.72e-3)
  bound-column self-inconsistencies — check whether the journal fixes them
  (our two `kind: conflict` data_gaps were explicitly parked "for the 2026
  journal").
- Possibly printed initial conditions (the AAS version printed none — the
  single biggest gap; the AAS mining note's §10 explicitly hoped the journal
  would print ICs).
- Any new (k1,k2) combinations beyond the AAS five.

## Method capability

No new genome — same fixed-Jacobi symmetric corrector + half-period Barden
monodromy + pseudo-arc-length continuation we already adopted. The new
capability is **demonstrated generality across μ**: the method is shown to find
stable cycler subfamilies at arbitrary mass parameter, which validates pointing
our (already-built) Ross corrector at *any* binary system, not just Earth-Moon.
That is directly actionable for Track A coverage and Track B prioritization.

## Proposed follow-on

1. **USER FETCH the PDF** (#116) — it is free but only the user can pull the
   numeric tables here (WebFetch denied). Highest value: the higher-μ
   (μ, C^stable, T^stable) tuples and any printed ICs.
2. On fetch, run the standard rescan discipline (two independent transcription
   passes, Decimal arithmetic cross-checks) as was done for the AAS version.
3. Then reproduce each higher-μ member with the adopted Ross corrector
   (same-model golden: EXPECTED = published (μ, C, T); recovered x0 is derived
   provenance, never a golden) → proposed new V1 rows in a μ regime the
   catalogue currently has none.
4. Feed the "universal stable subfamily near C^max" claim to the prioritizer
   (Track B): treat "(k1,k2) family, any μ, near C^max" as a predicted
   non-empty seed region for the discovery daemon, subsuming the relevant
   negative-registry entries under the new method-capability rule.

## Dedup / supersession

Supersedes (does not duplicate) `2026-06-11-ross-roberts-tsoukkas-2025-mining.md`
as the citable version-of-record once published in a journal. The 5 Earth-Moon
rows stay; the higher-μ members are additive. Forward-citation status: still
effectively one external citer (Braik & Ross arXiv:2605.31543, same lab,
already mined).

# Discovery Program — spec (2026-06-13)

**Premise (the reset):** the project's deliverable is **novel cyclers found and
validated**, not a literature census. Paper validation was instrumental — it earned
a **trusted, known-good toolkit**. That toolkit is now largely mature, so the
weight shifts from mining/validation to **discovery**. Validation does not go away:
it is the *quality bar* every novel candidate must clear (the unchanged V0–V5
gauntlet). It is simply no longer the scoreboard.

## What we already have (the earned toolkit — treat as done)
CR3BP corrector + Tier-2 propagator; multi-rev Lambert; FBS analytic-gradient legs
+ multi-arc chains + flyby coupling (#226); n-body harness (REBOUND/IAS15); GMAT V4
high-fidelity lane; continuation driver; MBH; self-seeding; reachable-set screening
(Braik-Ross energy-preserving #236 + Zhou-Armellin impulsive #239); the
method-versioned **empty-region registry**; the SILVER→gauntlet pipeline; the
solver-outcome logger (#210).

## What the discovery attempts so far tell us
Forge Phase 6, the multi-family sweeps (Jovian/Saturnian/heliocentric corners), and
the #218 continuation campaign returned **mostly empty** — captured honestly in the
negative-results registry. Read correctly, this is not "the space is empty"; it is
**"the current genome already searched the regions it can represent, and they're
sparse or already published."** The discovery frontier is therefore gated by
**capability (what topologies the genome can express)** and **coverage of the
unexplored space**, NOT by algorithm trust.

## The three discovery tracks

### Track A — Richer genome (HIGHEST leverage; tracked as #254)
A more capable genome opens regions the current one cannot even represent.
**Corollary (the coupling): a novel cycler usually arrives WITH a new algorithm** —
you find cyclers you couldn't before because you built/adopted the capability to
represent and close them (Liang's CGE members ⇐ repeated-moon multi-rev genome). So
the **capability-paper sweep (#250) is a primary discovery input, not a downranked
validation chore**: papers bring algorithms, algorithms make new cyclers findable.
Priorities, in rough order of evidenced payoff:
- **Repeated-moon multi-rev topologies.** The Liang CGE class *reopened* the Jovian
  empty region precisely because it lives in a topology our swept (zero-rev,
  single-encounter) genome does not cover (the §6b new-sourced-data arm, conditional
  on this genome). This is the most concrete "known to contain cyclers we can't yet
  represent" case.
- **Broken-plane / inclined** genome (the M-3D work lifted coplanarity for transfers;
  extend to the discovery genome).
- **Low-thrust** genome (Sims-Flanagan machinery exists; no sourced powered rows —
  but discovery need not be sourced).
Each new genome capability triggers the empty-region registry's §6b
capability-subsumption re-sweep of the regions it now reaches.

### Track B — The prioritizer (where to look)
Discovery search must be pointed, not blind (sweeping empty space faster finds
nothing). The reachable-set scorer is exactly this:
- **#249** — ungate the Braik-Ross scorer (recover the 4 unstable members so
  C32-dominance is testable).
- **#239 merge** — fold the Zhou impulsive reachable set in as the energy-changing
  complement (documented two-tier "cheap heading screen → budgeted impulsive
  confirm" path).
- **#162 family-selection** feeds the same question.
Output: a ranked list of accessible, non-empty regions/families to seed — so Track C
spends compute where cyclers plausibly are.

### Track C — Discovery-campaign daemon (how to run it; tracked as #253)
Novel-cycler search is **days of compute**, which is why one-shot agents HANG on it
(they background-and-poll; subagents get no completion signal — see
[[feedback_long_agents_commit_incrementally]]). Long compute belongs in a **daemon**,
like the surrogate-corpus daemon: a detached, resumable, quota-proof process that
sweeps a prioritized region, dedups candidates against the catalogue + the
negative-registry, and auto-feeds every SILVER survivor to the unchanged gauntlet.
This dissolves the agent-failure problem by moving long compute to where it belongs.
The corpus daemon is the working template to generalize.

## Cross-cutting
- **Gauntlet unchanged.** Every novel candidate clears the same V0–V5 bar. The
  validation work pays off *here*: trusted algorithms mean a discovery survives review.
- **ML false-positive flagger** (proposed): trained on our own past-bug signatures
  (ToF artifact, 63 s epoch offset, shared-cross-check non-independence) to flag
  suspicious "closures" for independent re-check. This matters MORE in discovery than
  validation — there is no published anchor to catch a false novel claim. Smallest,
  most-defensible ML item; build before the speculative surrogate.
- **ML surrogate / KKT amplifier (#240)** is DOWNSTREAM of Track B: it accelerates
  search but only pays off pointed at non-empty space, and only if the target lane is
  chosen (the corpus is CR3BP-corrector data, not the Ozaki flyby-cycler target — the
  logger needs no change, only a flyby-cycler data-gen driver). Defer until A+B exist.

## Task re-prioritization
- **Elevate (discovery-enabling):** #248 (multi-arc convergence — SOLVED for
  mcconaghy, 0.0987 km/s), #249 + #239-merge (the prioritizer), #250 (capability-paper
  sweep — DONE, no new genome paper exists), and the **Track-A genome effort #254**
  (repeated-moon multi-rev first).
- **New:** the discovery-campaign daemon **#253** (Track C) generalizing the corpus daemon; the
  ML false-positive flagger.
- **Downrank:** pure paper-validation/acquisition (#116) now matters only when a paper
  brings a NEW algorithm or topology, not just another anchor. The FBS-default ladder
  (#245/#246) is low-stakes housekeeping (FBS is already usable opt-in).
- **Keep deferred:** #240 (behind A+B), #115 (no published V∞).

## Success metric
**One novel cycler found AND passed the gauntlet** — distinct from, and the point of,
the 277-row literature census. Everything above is in service of that.

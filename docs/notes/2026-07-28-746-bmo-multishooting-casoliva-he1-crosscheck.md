# Finding: does BMO 2009's multiple-shooting extension apply to Casoliva's He1? (#746)

**Task:** `#746`, follow-up cross-check flagged by `#742` (item 2) on the `#725` Casoliva
digest. Question: does Barrabés/Mondelo/Ollé (BMO) 2009's own multiple-shooting
augmentation of their homoclinic-connection continuation system apply to — and would it
materially change or improve — Casoliva et al.'s own He1 family construction (113.6-day
connection, Casoliva's Eq. (20), presented in single-shooting form)?

**Sources read in full (not just the two digests):**
- `casoliva-...-2010-two-classes-cycler-trajectories-earth-moon-jgcd-33-5-...pdf` (+ `.txt`
  sidecar) — full text, §V.A/B/C (Homoclinic Connections and Their Relation to Cyclers,
  Continuation of Homoclinic Connections, Results and Discussion), Table 4/5/6, Fig. 8-12.
- `casoliva-...-2008-families-cycler-trajectories-earth-moon-AIAA-2008-6434.txt` — full text,
  grepped for "shoot".
- `barrabes-mondelo-olle-2009-numerical-continuation-homoclinic-connections-...pdf`, converted
  to text with `pdftotext -layout` for this task (no prior `.txt` sidecar existed) — full text,
  Abstract, Introduction, §2.2 (the system-of-equations/multiple-shooting section), §3.1
  (L1/L2, Earth-Moon), §3.2 (L3, Sun-Jupiter), Acknowledgments/References.

## Verdict: NO — clean negative. BMO's own paper scopes the multiple-shooting motivation to
## its L3 (Sun-Jupiter) case, not the L1 Earth-Moon case He1 belongs to, and Casoliva's own
## numbers/converegence evidence are consistent with that scoping. It also doesn't matter for
## this project's own code, because we do not have a literal port of Casoliva/BMO's continuation
## method at all — our homoclinic-connection tool is architecturally different.

### 1. BMO 2009 states no numeric threshold anywhere — only a qualitative, section-scoped one

Exhaustive grep of the full BMO text for "condition number", "ill-condition", "precision",
"digits", "accuracy", "single shoot" found **zero** occurrences of any explicit numeric
threshold (no stated integration-time cutoff in days or nondimensional units, no condition-number
bound). The paper's only statement tying multiple shooting to a concrete cause is in §2.2, right
after presenting the single-shooting system (their eq. (7)):

> "In our setting, the integration times `T^u`, `T^s` may become large (**this is the case in
> section 3.2**). In order to avoid loss of precision, we have also used a multiple shooting
> version of (7)."

Section 3.2 is explicitly the **L3, Sun–Jupiter horseshoe-orbit case** (`µ = µ_SJ = 9.53875e-4`,
stated at the top of §3.2) — a different libration point AND a different mass ratio from
Casoliva's He1 (L1, Earth-Moon, `µ_EM ≈ 0.01215`). Section 3.1 (L1/L2, Earth-Moon — the section
that computes the **same** `He_j`/`Hm_j`/`Hi_j` family labels Casoliva reuses verbatim, confirmed
at §3.1.2: "we have labelled the corresponding families `He_j`, j = 1, 2, 3, 4") is never named as
a case needing multiple shooting anywhere in the paper. (The Introduction's more general sentence,
"[t]he instability due to the hyperbolic character of all the p.o. considered is coped with using
a multiple shooting strategy," reads as a summary of the paper's toolset as a whole, not a claim
that every individual family in §3.1 required it — §2.2's more specific statement, tied
concretely to §3.2, is the one with an actual named locus.)

So there is no literal "threshold in days" to compare Casoliva's 113.6-day connection against —
BMO's own paper simply never publishes one. The comparison that BMO's own text supports is
qualitative: **is Casoliva's He1 in BMO's own L1/L2 Earth-Moon regime (§3.1, not flagged as
needing multiple shooting), or in BMO's own L3 Sun-Jupiter regime (§3.2, the one explicitly
flagged)?** It is squarely in the former — He1 is literally one of BMO's own §3.1.2 families.

### 2. Casoliva's own numbers are unremarkable relative to the p.o. period, not "large"

Casoliva's Table 4 (2010, p.1635) reports, for the closest-approach He1 connection at
`h = −1.45016232260699`:

- `T` (Lyapunov p.o. period) = 6.706878522271349 (nondim) = 29.1640 days (Table 5 caption).
- `T^u` = 27.6856812343605 (nondim), `T^s` = −27.33008060916397 (nondim).

Converting via the paper's own stated period-to-days ratio (29.1640 days / 6.706879 nondim units
≈ 4.3484 days/nondim unit, consistent with the Earth-Moon sidereal-month-based nondimensionalization),
`T^u` ≈ 120.4 days and `|T^s|` ≈ 118.8 days — **roughly 4.1-4.4 multiples of the p.o.'s own
period**, not an extreme multiple. (The separately reported "connection flight time" of
113.6319 days, periselene 1→19, is a different derived quantity — physical elapsed time along
the connection counting periselene passages — not literally `T^u + |T^s|`; both numbers describe
the same connection at a broadly consistent ~110-120-day scale, well short of anything BMO's own
L3 discussion implies by "large".) By contrast, BMO's own L3/Sun-Jupiter case involves
horseshoe-orbit manifold tubes that loop repeatedly around islands of forbidden motion near the
tangency energy `h_t`, generating (per BMO's own text) "an infinite number" of higher-order
homoclinics as energy approaches `h_t` — a qualitatively different, much longer/more-convoluted
regime than a single ~4-period Earth-Moon Lyapunov connection.

### 3. Casoliva's own text shows clean single-shooting convergence, no caveat

Both Casoliva papers were grepped in full for "shoot" — the only hit in either is the unrelated
reference-list title "Shoot the Moon 3D" (Parker & Lo 2006). Neither paper ever mentions single-
vs multiple-shooting as a design choice.

Casoliva's own §V.C text reports CPU times for the He1 family's successive connections in Fig. 8
(ordered by increasing energy): 121, 196, 397, and 1619 s, and states explicitly: **"The
continuation procedure slows down as energy increases due to the decreasing perigee distance"**
— i.e., the paper's own stated reason for growing cost is close-approach geometry (perigee/periselene
distance shrinking, a near-singularity effect), not growing integration time or conditioning
loss from `T^u`/`T^s`. The closest-approach example (Fig. 9, Table 4, `h = −1.450162`) took 5h20min
of CPU — a large number in absolute terms, but the paper offers this without any convergence
caveat, loosened tolerance, or failure/retry language anywhere nearby. This is consistent with —
not contradicting — BMO's own scoping: Casoliva's single-shooting corrector converged cleanly in
exactly the regime (L1/L2 Earth-Moon) BMO's own paper never flags as needing the augmentation.

### 4. Conclusion on the paper cross-check

Casoliva's Eq. (20) is presented in single-shooting form, and that is the methodologically
correct choice for their own regime by BMO's own paper's stated scoping (large-`T` multiple
shooting motivated concretely by §3.2/L3/Sun-Jupiter, not §3.1/L1/Earth-Moon where He1 lives).
The `#742` digest's characterization of this as an "omission" in Casoliva's paper is accurate as
a literal absence-of-mention, but this pass finds the omission does not indicate a hidden problem
— Casoliva's own numbers (moderate `T^u`/`T^s` at ~4 p.o. periods, monotone increasing but modest
CPU cost explicitly attributed to perigee geometry, no convergence caveats) are exactly what one
would expect from a single-shooting corrector that did not need the augmentation. **No correction
to `#725`'s golden-numeric-table digest is needed** — the 113.6319-day connection flight time and
surrounding Table 4/5/6 values are confirmed accurate and are not called into question by any
single/multiple-shooting concern.

### 5. Does this matter for this project's own code?

`grep -rniE "homoclinic|multiple_shoot|multi_shoot|He1|BMO|barrabes" src/cyclerfinder/search/`
turns up:
- `src/cyclerfinder/search/ccr4bp_heteroclinic_search.py` — this project's actual
  homoclinic/heteroclinic manifold-tube-intersection search tool (used for the `#701`/`#702`
  Uranus Umbriel-Titania CCR4BP torus-homoclinic discovery). **Its own docstring explicitly
  states it does NOT reuse `cr3bp_multiple_shooting.py`**: "the natural unknowns here are exactly
  the two continuous phases + two continuous flow times, a smaller, differently-shaped problem
  that does not fit that module's node/segment-time contract without a bigger adaptation than a
  from-scratch 4-unknown least-squares call." It refines a single candidate connection via
  `scipy.optimize.least_squares` over `(theta2_u, t_u, theta2_s, t_s)`, guarded by an independent-
  integrator (Radau vs DOP853) ghost-guard and a fixed `ref_vec` eigenvector-sign anchor
  (`#702`'s fix for a spurious lobe-mismatch artifact) — not a Newton/QR-pivoted continuation of
  an over-determined system like BMO's eq. (7)/Casoliva's eq. (20).
- `src/cyclerfinder/search/cr3bp_multiple_shooting.py`, `src/cyclerfinder/genome/multi_shooting.py`,
  `src/cyclerfinder/genome/family_switch.py` — this project's actual multiple-shooting
  infrastructure. Per `family_switch.py`'s own docstring, it exists for a **different** problem:
  closing PERIODIC orbits through an N-node periodicity loop (helping convergence specifically
  "near NRHO bifurcations"), not for continuing homoclinic-CONNECTION families in energy à la
  BMO. It is available in the codebase but architecturally unrelated to, and not wired into, the
  homoclinic-connection search path.

**Conclusion: this project does not have a literal implementation of Casoliva/BMO's
predictor-corrector homoclinic-connection-family continuator (single- or multiple-shooting) at
all.** `ccr4bp_heteroclinic_search.py` solves a differently-shaped problem (single-candidate
least-squares refine, not family continuation via an over-determined Newton system), so the
question "would BMO's multiple-shooting extension improve OUR corrector" doesn't have a direct
answer — there is no analogous corrector in this codebase to improve. The `#701`/`#702` known
failure mode in this project's homoclinic-search tooling was a `ref_vec` eigenvector-sign/lobe
anchoring artifact (a spurious ~1 km Radau/DOP853 disagreement), not a long-integration-time
precision-loss issue of the kind BMO's multiple shooting addresses.

## Recommendation

No code change and no correction to `#725`/`#742`'s digests. This is a clean negative on both
halves of the question (paper-level: BMO's own threshold doesn't reach Casoliva's He1 regime;
project-level: we have no corrector that would be affected either way). The one live, deferred
idea worth flagging for a future task (not registered here, per this task's scope): `#725`'s own
digest already noted Casoliva's He1/Hm1/Hm2 golden numeric table (Table 4-6) is "the first sourced
numeric target" for a future Barrabés-Mondelo-Ollé-style homoclinic-connection continuator if this
project ever builds one from scratch (as opposed to `ccr4bp_heteroclinic_search.py`'s existing
least-squares-refine approach). *If* such a continuator is ever built, BMO's multiple-shooting
augmentation (their eqs. following (7), §2.2) should be ported alongside the single-shooting
baseline as a documented option for other, longer-`T` regimes (e.g. an L3-analogue or higher-energy
`He2`/`He3`/`He4` members if their `T^u`/`T^s` turn out to grow — not verified in this pass, out
of scope) — but nothing in this cross-check shows that need exists today, and no new task number
is being self-assigned for it.

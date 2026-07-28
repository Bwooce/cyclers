# `#756`: relaxed-period-criterion search for Anderson-Lo 2011's 5:6-LO

**Task:** `#756`, dispatched 2026-07-28 as a direct continuation of `#755`
(see `docs/notes/2026-07-28-755-jupiter-europa-3-4-lo-5-6-lo-targeted-search.md`
for the full prior search history and its own "Reviewer verdict" section,
which is the reason this task exists — this note does not repeat that
background). `#755`'s own coordinating-session review confirmed 3:4-LO
(target `|λ|=1036.116088`) on the weight of evidence despite a real 2.1%
period offset, reasoning from Anderson & Lo 2011's own p.171 text that CR3BP
resonance is only *approximate* for strongly unstable families far from the
two-body integrable limit — and flagged that `#755`'s own 5:6-LO search
(target `|λ|=4445.387515`) may have implicitly filtered out its own answer
by favoring `period_over_2π` near the integer `q=6` while scanning. This
task's job: redo the 5:6-LO search with that period bias explicitly removed.

**Code:** `src/cyclerfinder/search/jovian_resonant_families.py` (extended) +
`tests/search/test_jovian_resonant_families.py` (32 tests, up from 25, all
passing). Raw search checkpoint (159 candidates, JSONL, one record per
converged orbit found):
`data/found/756_jupiter_europa_5_6_lo_relaxed_period/candidates.jsonl`.

---

## Method: rank by eigenvalue closeness first, period second

Checked the actual scan code first, as instructed: `survey_candidates`
itself does **not** filter or sort by period at all — it grids the
corrector's own crossing-miss function `g(x0) = xdot` at a fixed
`half_crossings` and Jacobi constant, brackets every sign change, and
Newton-converges each bracket via
`cr3bp_periodic.correct_symmetric_fixed_jacobi`, returning **every**
converged candidate regardless of its resulting period. The implicit
period bias in `#755`'s own search was therefore not a code-level filter,
but a *search-time judgment call*: which `(x0, half_crossings)` windows to
spend compute on, and which of the resulting candidates to report as "the"
closest miss. `#755`'s own narrative ("the single closest *simultaneous*
approach") reveals this — it was looking for candidates good on *both*
axes at once, rather than doing a two-stage search (rank by eigenvalue,
then separately assess period as corroboration).

This task instead ran a wide `(x0, ydot0_sign, half_crossings)` sweep —
`x0 ∈ [-1.9, 1.7]` (**both** sides of the secondary, not just the
`x0 ≈ -1.42` hotspot `#753`/`#755` had already focused on), both `ydot0`
signs, `half_crossings ∈ {2,4,5,6,7,8}`, at `ANDERSON_LO_C_FLYBY` directly
— logging **every** converged candidate (159 total) to a JSONL checkpoint,
then ranking purely by `|λ_recovered − 4445.387515| / 4445.387515`, with
period and Europa-approach distance checked only *afterward*, as
corroboration signals, exactly mirroring how `3:4-LO`'s confirmation
actually worked (near-machine-precision eigenvalue match, corroborated
*separately* by shape and close-approach altitude, not by period).

New module additions (kept, not throwaway): `europa_closest_approach()` —
a proper, tested, reusable "close flyby" corroboration check, extracting
the ad hoc trajectory-sampling `#755` did inline for 3:4-LO into a real
function — and `_756_RELAXED_SEARCH_NEAR_MISSES` /
`recover_756_near_miss()`, the three best-eigenvalue-match candidates
found by this task's sweep, preserved with full provenance (not just
mentioned in this note) per this project's negative-results-registry
discipline.

---

## Result: still NOT CONFIRMED — and the relaxed criterion changes nothing

| Candidate | `x0` | `ydot0_sign` | `half_crossings` | Recovered `\|λ\|` | Rel. err | `period/2π` | Europa closest approach |
|---|---|---|---|---|---|---|---|
| **Pre-existing (`#753`)** | 0.81360506 | +1 | 2 | 4533.602947 | **1.98%** | 16.11 (2.7x) | 0.0183 nondim (~12,300 km) |
| `756-nm-1` (new best) | 1.370014 | −1 | 5 | 3760.867 | 15.4% | 25.02 (4.2x) | 0.0423 nondim (~28,400 km) |
| `756-nm-2` | −1.384961 | +1 | 6 | 6864.864 | 54.4% | 28.04 (4.7x) | 0.0241 nondim (~16,200 km) |
| `756-nm-3` | −1.721452 | +1 | 7 | 1421.336 | 68.0% | 13.02 (2.2x) | 0.0199 nondim (~13,400 km) |
| *(reference)* confirmed 3:4-LO | −1.430408 | +1 | 6 | 1036.116117 | 2.8e-8 | 4.086 (2.1% off `q=4`) | **0.00245 nondim (~1,641 km)** |

The core finding, in one sentence: **relaxing the period-proximity
criterion, exactly as the `#755` ruling recommended, does not change the
verdict.** No candidate found this task — across a search domain
substantially wider than `#755`'s own (both sides of the secondary, six
`half_crossings` values, both `ydot0` signs) — beats the pre-existing
`#753` candidate's 1.98% eigenvalue match, and critically:

1. **No candidate has a plausible period.** Every candidate examined
   (the pre-existing one included) has `period_over_2π` 2.2x–4.7x the
   naive `q=6` value — not a "few percent to maybe 10-20% off" the way
   `3:4-LO`'s confirmed candidate was only 2.1% off `q=4`. This is a
   qualitatively different, much weaker situation, closer to `#753`'s
   *original* rejected near-misses (which had periods bearing no relation
   to any plausible `q`) than to `3:4-LO`'s striking case.
2. **No candidate makes a close Europa approach.** `europa_closest_approach()`
   (new, tested function) shows the confirmed `3:4-LO` orbit passes
   `~1,641 km` from Europa — a genuine close flyby, matching the paper's
   own attributed instability mechanism ("much of the unstable
   characteristics of this trajectory arise from the fact that it is
   performing a close flyby of Europa," p.177-178). Every 5:6-LO candidate
   examined here — the pre-existing best *and* all three new near-misses —
   stays `12,000–28,000 km` away, an order of magnitude farther, with no
   qualitative close-flyby signature at all.
3. **The best eigenvalue match available anywhere (1.98%, the pre-existing
   candidate) is five to six orders of magnitude less precise than
   3:4-LO's (2.8e-8).** It is not the kind of "essentially exact,
   corroborated on other axes" match that justified a reviewer judgment
   call for 3:4-LO — it is a modest near-miss, in the same rough tier as
   `#753`'s original 5:6-NO miss (27% eigenvalue error), not a striking one.

So this is not a case of "eigenvalue matches but period doesn't" (the
3:4-LO pattern) — it is a case where **neither** criterion is close, and
the two best-available signals (eigenvalue proximity, Europa-approach
distance) point at *different* candidates entirely (the pre-existing
0.8136 seed has the best eigenvalue but a merely-ordinary approach
distance; none of the new near-misses improve on either axis
simultaneously). There is no candidate this task can even present as "a
candidate found, reviewer judgment needed" the way 3:4-LO was — the
evidence itself does not clear that bar.

---

## What was searched (for future reference / anti-catalogue value)

* `x0 ∈ [-1.9, 1.7]` at `C = ANDERSON_LO_C_FLYBY = 2.99163956830415`
  (covering both the `x0 ≈ -1.42` "outer negative" hotspot `#753`/`#755`
  already found rich, AND the `x0 ≈ +1.3` to `+1.7` "outer positive" side,
  which `#755`'s own search covered less finely).
* `ydot0_sign ∈ {+1, -1}`.
* `half_crossings ∈ {2, 4, 5, 6, 7, 8}` (2 matches the pre-existing
  candidate's own complexity; 4-8 matches 3:4-LO's successful
  `half_crossings=6` and its near neighbors).
* Fine-grid zooms (`n_grid` up to 6000) around every promising region found
  by the coarser passes, mirroring `#755`'s own successful 3:4-LO strategy.
* Direct (non-bracketed) convergence sweeps from ~60 seed guesses around the
  pre-existing `x0=0.8136` candidate, confirming it sits in an isolated,
  well-converged basin (not a fluke of one lucky seed) — and surfacing a
  real, minor methodological caveat: `survey_candidates`' own bracket-based
  scan (using a fixed scan-only `t_hi`) did *not* rediscover this exact
  point via sign-flip detection at `half_crossings=2` (the scan's `g(x0)`
  stayed one-signed across a window containing the known root), while
  direct Newton convergence from a nearby guess reliably finds it. This is
  consistent with the low-`half_crossings` regime being more sensitive to
  the scan-tool's own integration horizon than the corrector's internal one
  (`t_hi = 1.25 * period_guess`) — a real tooling nuance, not a data error;
  ground truth for this candidate is the direct `recover_table1_candidate`
  reproduction, which is exact and stable.
* Confirmed the fractal/chaotic sensitivity `#753`/`#755` already
  documented is not confined to the one previously-explored hotspot: eigen-
  values ranging from `~1` (trivial) to `~1.9×10^8` were found within
  `Δx0 < 0.001` windows in *multiple* separate regions across the full
  `x0` range searched.

None of this searched territory is now "confirmed empty" in the strong
sense — the search domain is continuous and this task, like `#753`/`#755`
before it, sampled a large but necessarily finite set of `(x0,
half_crossings)` combinations. But three independent tasks' combined effort
(grid+bisection at `C_flyby` directly across a very wide `x0` domain, both
signs, `half_crossings` 1-40 across all three tasks combined, continuation-
descent from higher `C`, and the two-body flyby-vector-rotation seed
construction) have not found a 5:6-LO candidate that is both eigenvalue-
close and independently corroborated. This is reported as a genuine,
well-evidenced, **continued non-confirmation** — not fudged, not forced.

---

## Verification

* `uv run ruff check` / `ruff format --check` on both changed files: clean.
* `uv run mypy src tests` (project canonical invocation): clean, 821 files.
* `uv run pytest tests/search/test_jovian_resonant_families.py -q`: 32/32
  pass.
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: pass (run
  before committing the `OUTSTANDING.md` update).

---

## Opinion on `#754` (Task B) — not a decision

`#754` needs both `3:4-LO` and `5:6-LO`'s manifolds for its Table 2/3 gate
(or a re-scope). `3:4-LO` stands reviewer-confirmed per `#755`'s own
verdict. `5:6-LO` remains genuinely, thoroughly unconfirmed — and unlike
`3:4-LO`, there is no close call here for a reviewer to weigh: the
best-available candidate (the pre-existing 1.98%-eigenvalue-error one)
lacks the corroborating evidence (period plausibility, close Europa
approach) that made `3:4-LO`'s judgment call reasonable, and this task's
substantially widened relaxed-period search did not change that picture.
My opinion: `#754` still does not clear `#753`'s own bar for dispatch as
originally scoped (needs both rows' manifolds). The two live options
remain what `#753`/`#755` already identified: (a) further, differently-
targeted search for `5:6-LO` (a genuinely new seed strategy — e.g. actually
carrying out the paper's own multi-patchpoint flyby-vector-rotation
refinement `#755` built the first piece of but didn't have time to
complete — rather than more grid+bisection at `C_flyby` directly, which
three tasks have now applied at length), or (b) re-scoping `#754` around
only the two confirmed families (`5:6-LI`, `3:4-LO`). This is my
assessment for the user to weigh, not a decision — `#754` stays HELD.

# KKT pseudo-body database amplifier — design spec (#210, second half)

**Status:** DESIGN ONLY. Implementation deferred until the #210 outcome-logger
has accrued real solves and we decide to build the surrogate. No code here.

**Goal.** Multiply the #210 outcome-log corpus by ~one order of magnitude *per
real solve*, so the Ozaki surrogate's training floor (~7×10⁶ samples,
arXiv:2111.11858 / JGCD 2022) is reachable from ~7×10⁵ actual trajectory
optimizations instead of 7×10⁶. The amplifier manufactures additional valid
(inputs → cost) training tuples from a *single* converged optimal
Earth-body-Earth (E-B-E) / DSM-cell solution, with no extra optimization.

## 1. The idea (Ozaki contribution 3, in our terms)

A converged optimal E-B-E leg is a spacecraft trajectory `x*(t)` on `[t0, tf]`
that performs one flyby of the target body at epoch `t*`, where the body and the
spacecraft share position `r* = r_sc(t*) = r_body(t*)`. The optimization (min
total DSM ΔV subject to the patched-conic + flyby constraints) satisfies the
Karush-Kuhn-Tucker conditions at the optimum.

**Key observation:** the spacecraft trajectory `x*(t)` is *agnostic to which body
it flies by* — it only needs SOME body to be at `r*` at `t*`. So for ANY
"pseudo-body" whose Keplerian orbit passes through `r*` at epoch `t*`, the SAME
trajectory `x*(t)` is a *feasible* solution of that pseudo-body's E-B-E problem,
with the SAME cost (ΔV, ToF, emergent v∞). If, in addition, the pseudo-body's
velocity at `t*` is chosen so the KKT conditions still hold (the inserted body
does not perturb the optimum — see §2), then `x*(t)` is also the *optimal*
solution, and `(inputs, œ_pseudo) → (ΔV*, tf*, v∞*)` is a valid training tuple.

One real solve → many pseudo-bodies along (and consistent with) `x*` → many free
tuples.

## 2. KKT-consistency condition (what makes a pseudo-body admissible)

A pseudo-body placed at `r*` at `t*` defines a one-parameter-family of Keplerian
orbits through that point (varying the velocity `v_body(t*)`). Not all preserve
optimality:

- **Position match (feasibility):** `r_body,pseudo(t*) = r*`. Necessary for the
  flyby to occur on `x*`.
- **Flyby-geometry / KKT match (optimality):** the gravity-assist or DSM flyby
  in our MGA-1DSM model constrains the *relative* velocity v∞ = v_sc(t*) −
  v_body(t*). For the unpowered gravity-assist case the flyby turns v∞ at fixed
  magnitude; the KKT optimum fixes the incoming/outgoing v∞ split. A pseudo-body
  is admissible iff its `v_body(t*)` reproduces the SAME v∞ vector(s) the parent
  optimum used — i.e. `v_body,pseudo(t*) = v_sc(t*) − v∞*`. With position and the
  flyby v∞ both matched, the parent's costates/multipliers transfer unchanged →
  KKT holds → `x*` is optimal for the pseudo-body too.

So each *parent solve* yields pseudo-bodies by: (a) sampling alternative *flyby
epochs* `t*_k` along `x*` where the spacecraft is at a gravity-assist-feasible
state, and/or (b) at a fixed `t*`, the admissible `(r*, v_body)` is essentially
pinned by the KKT match — so the amplification comes mainly from (a) sampling
multiple feasible flyby points on the parent trajectory, each defining a distinct
pseudo-body orbit (distinct œ = [a,e,i,Ω,ω,M] back-computed from `(r*, v_body)`).

(Open design question to resolve in implementation: how many genuinely-distinct
admissible flyby points a single E-B-E arc affords in OUR MGA-1DSM model vs
Ozaki's — this sets the realized amplification factor. Verify empirically, §4;
do NOT assume 10× — measure it.)

## 3. Where it hooks in

- **Input:** the #210 outcome log must capture, for each REAL converged solve,
  enough to reconstruct `x*`: the genome/inputs already logged PLUS the converged
  flyby state(s) `(t*, r*, v_sc(t*), v∞*)` and the cost outputs. → ACTION:
  confirm/extend the #210 logger schema to include the flyby state(s); if absent,
  add them (logger is the dependency).
- **Amplifier:** a post-processor `surrogate/pseudo_body_amplifier.py` that reads
  real-solve records, generates admissible pseudo-bodies (§2), back-computes
  `œ_pseudo`, and emits amplified tuples tagged `source: "kkt-pseudo"` (distinct
  from `source: "real-solve"`).
- **Output:** the same JSONL corpus, with a provenance tag so real vs synthetic
  is always separable (never mix silently).

## 4. Validation gate (reproduce-before-trust, applied to synthetic data)

The KKT construction is a CLAIM ("pseudo-body cost == parent cost"). It must be
empirically verified, not trusted:

- For a random sample (e.g. 5%) of generated pseudo-bodies, actually run the real
  E-B-E / DSM solver on `œ_pseudo` and confirm the converged cost matches the
  parent's to tolerance. A mismatch means the admissibility condition (§2) is
  wrong or our flyby model differs from Ozaki's — STOP and fix, do not ship a
  corpus built on an unverified amplifier.
- A property test: for an admissible pseudo-body, `r_body,pseudo(t*)` equals `r*`
  and `v_sc(t*) − v_body,pseudo(t*)` equals `v∞*` to machine precision (the
  construction's own invariants).

## 5. Discipline boundary (hard)

- Amplified tuples are **TRAINING DATA ONLY** — never catalogue rows, never
  goldens, never validation evidence. A pseudo-body is a synthetic construct, not
  a real body; the surrogate trained on them only *proposes* search candidates,
  which still pass V0→V5 to earn anything. (See the orbit-closure / golden
  discipline.)
- Provenance tag mandatory: every corpus record carries `real-solve` vs
  `kkt-pseudo` so a future trainer can weight/audit them and we can measure the
  surrogate's accuracy on real-only holdout.

## 6. Scope / sequence

1. (#210, in progress) outcome logger captures real solves — DEPENDENCY.
2. Extend logger schema to include flyby state(s) if not already.
3. Build `pseudo_body_amplifier.py` + the §4 verification gate + property tests.
4. Measure the realized amplification factor (§2 open question) on a real batch.
5. Only then assess distance to the ~7×10⁶ floor and whether to train.

Do NOT build 3-5 until the logger is accruing real solves and the user calls it —
this spec exists so the path is captured, not to trigger immediate work.

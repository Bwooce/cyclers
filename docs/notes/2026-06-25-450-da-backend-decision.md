# #450 DA/HOTM backend decision record

**Date:** 2026-06-25
**Issue:** #450 (DA/HOTM global multi-rev CR3BP fixed-point enumeration lane)
**Status:** DECIDED — pure-Python truncated Taylor-map deliverable backend.

## Environment probe (recorded)

Probed once on 2026-06-25 (`uv run python -c "import importlib.util ..."`):

| module | available |
|---|---|
| `daceypy` | **NO** |
| `dace` | **NO** |
| `pyaudi` | **NO** |
| `mosek` | **NO** |
| `cvxpy` | **NO** |

No differential-algebra library and no convex/commercial optimizer is installable
in this environment today (matches the design-draft §0 / §8.1 verified state).

## The three options (design-draft §8.1 / §10.1)

1. **(a) Vendor DACEyPy + obtain a MOSEK academic licence** — paper-exact, fastest,
   licence cost. REJECTED: MOSEK is commercial; DACEyPy/DACE/pyaudi are not present
   and vendoring a DA library + a commercial solver is out of scope and against the
   project ground rules (NO MOSEK, NO cvxpy-commercial, NO DACEyPy/dace/pyaudi).
2. **(b) Pure-Python truncated Taylor-map fallback** — a single-revolution section
   map expanded to a chosen Taylor order about a section reference, composed to
   `Pⁿ` by truncated polynomial composition, with a non-commercial fixed-point
   root-finder (`scipy.optimize` / polynomial homotopy). NO MOSEK. **CHOSEN.**
3. **(c) Ship on the sampling backend, defer DA entirely** — proves the lane but
   leaves the paper's acceleration layer unbuilt. Superseded by (b): the
   Taylor-map IS buildable in pure Python and is the deliverable.

## Decision (USER, 2026-06-25)

**Chosen deliverable backend: option (b), a pure-Python truncated Taylor-map
`DASectionMap`.** No MOSEK, no cvxpy-commercial, no DACEyPy / dace / pyaudi at any
point.

Build order (de-risks the Taylor backend behind a validated oracle):

1. Build the cheap **`SamplingSectionMap` FIRST** (plan Tasks 3–6) as the
   *validation oracle* and the Png' recovery proof. This is the brute-force
   float-propagator realization of the method's geometry (the `reachable_impulsive.py`
   precedent), recovering Png' now at the cost of the paper's speed.
2. Build the pure-Python Taylor-map **`DASectionMap` to the SAME `SectionMap`
   interface** (the former "Task 7 deferred DA acceleration" is **PROMOTED to an
   in-scope deliverable**), and validate it against BOTH:
   - the sourced **Png' golden** (arXiv:2509.12671 printed tables), and
   - the **`SamplingSectionMap` oracle** (assert the two backends agree on the
     recovered fixed points / section geometry to a tolerance).

The validation path (Png' recovery + base-family triangulation) does NOT depend on
the DA backend choice: the sampling backend carries the capability proof, and the
Taylor backend is validated against it as a swappable seam.

## Honesty flag (carried from the design draft)

The METHOD is fully digested (`2026-06-13-high-order-transfer-map-2509.12671-mining.md`);
this is an *implementation* choice, not method invention. The paper's DA/ADS/SCOP
acceleration (DACEyPy + MOSEK) is environment-gated and we do NOT reproduce it; the
pure-Python Taylor map is the non-commercial substitute for its single-rev
polynomial-map + composition core. If the pure-Python Taylor backend cannot recover
Png' within the validation tolerance after a genuine escalation effort (higher
order, finer reference placement, better root-finder, domain splitting), that is a
REAL finding to report honestly — the sampling backend carries the validation in
that case, and a characterized honest negative on the Taylor backend is an
acceptable outcome. A faked pass is not.

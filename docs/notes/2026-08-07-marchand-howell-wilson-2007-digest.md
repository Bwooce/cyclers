# Digest: Marchand, Howell & Wilson 2007, "Improved Corrections Process for Constrained Trajectory Design in the n-Body Problem"

**Source:** B.G. Marchand, K.C. Howell (Purdue) & R.S. Wilson (JPL), *Journal of Spacecraft and
Rockets* 44(4):884-897 (2007). DOI `10.2514/1.27205`. Filed at
`cyclers_pdf/papers/marchand-howell-wilson-2007-improved-corrections-process-constrained-trajectory-design-n-body-jsr-44-884-doi-10.2514-1.27205.pdf`
(md5 `259351abcee68ba3f46b61c47ae79d4e`), text-layer, 14 pages.

**Acquisition context:** user-supplied, in response to a search for the unobtainable Wilson 2003 JPL
internal memo (`#730` backlog item 98, "Derivation of differential correctors used in GENESIS mission
design," JPL IOM 312.I-03-002, no DOI, not publicly distributed). **This paper's own third author is
the SAME Roby S. Wilson** — this is very plausibly the published, generalized, peer-reviewed
successor to that internal memo's own methodology, not merely a related paper. Directly on-point,
arguably a BETTER source than the memo itself (peer-reviewed, DOI-bearing, publicly obtainable).

## What this paper actually contains — a generalized constrained multi-arc corrector

Extends the classic multi-patchpoint corrector (Howell & Pernicka 1988, `#730` item 97, already
digested this session) to handle **arbitrary algebraic interior/exterior constraints**, not just
position/velocity continuity — including, explicitly, **periodicity as a constraint without requiring
prior knowledge of the solution's symmetry**. This is the single most directly relevant feature for
`#782`'s own problem (closing a periodic resonant-chain orbit with no known symmetry axis assumed in
advance).

Four-step design strategy: (1) model the trajectory as a series of arcs (analytical, numerical, or
conic; established from a three/four-body model) to establish general characteristics and timing;
(2) formulate the specific constraints and their partials; (3) run the corrections process to enforce
position/velocity continuity while satisfying the constraints (also yielding preliminary maneuver
requirements where applicable); (4) transition the converged solution to a full ephemeris model.

Core mechanics: the trajectory is decomposed into segments/patch points (their own explicit guidance:
"at least four patch points per revolution" near libration points, though "two patch points may be
sufficient" in simpler cases — a concrete, sourced patch-density heuristic `#782` could use directly).
Velocity discontinuities `ΔV_k` at each patch point are driven to zero via a Newton process built on
the state transition matrix `Φ(t,t0)` (their Eq. 5-6, the standard variational-equation STM
propagator this project's own code already implements) plus a **noncontemporaneous variation**
formulation (explicitly flagged by the authors as beneficial for deriving the corrections process —
worth checking whether this project's own multiple-shooting implementation already uses this
formulation or the simpler contemporaneous one). Constraint partials `∂ψ_k/∂X` are assembled
per-patch-point into the same linear system as the continuity conditions, so periodicity (or any other
algebraic constraint) is enforced SIMULTANEOUSLY with continuity, not as a separate outer loop.

## Registration

Filed in `cyclers_pdf`. `CORPUS_INDEX.md` and `#730` backlog item 98 to be updated to reflect this as
a strong substitute for the unobtainable JPL IOM. Directly relayed to `#782` (in progress at time of
this digest) rather than left for it to discover independently — its periodicity-as-constraint
formulation and patch-density heuristic are both immediately actionable for that task's own problem.

# Triage digest: Fu, Peng, Gong and Shi 2025, teardrop hovering formation along the NRHO

Date 2026-10-04. Triage only: the paper was not requested and is probably outside the project's scope. This
note establishes what it does, transcribes every printed table, and gives a scope verdict.

Citation: S. Fu, Y. Peng, S. Gong and P. Shi, "Design and Continuation of Nonlinear Teardrop Hovering Formation
along the Near Rectilinear Halo Orbit", arXiv:2504.11771v1 [math.OC], 16 April 2025, 18 PDF pages (17 numbered
pages plus a highlights page); preprint submitted to Aerospace Science and Technology. Authors at Beihang
University (Beijing). Filed in the private paper corpus as
fu-peng-gong-shi-2025-design-continuation-nonlinear-teardrop-hovering-formation-NRHO-arxiv-2504.11771v1.pdf.

Evidence labels: READ means read from the printed page in this session (all pages read as images); INFERRED
means my reading. Page numbers are the paper's own (1 to 17).

## 1. What it does (READ)

Abstract (p. 1): "This short communication is devoted to the design and continuation of a teardrop hovering
formation along the Near Rectilinear Halo orbit and provides further insights into future on-orbit services in
the cislunar space. First, we extend the concept of the teardrop hovering formation to scenarios along the Near
Rectilinear Halo orbit in the Earth-Moon circular restricted three-body problem. Then, we develop two methods
for designing these formations based on the nonlinear model for relative motion. The first method addresses the
design of the teardrop hovering formations with relatively short revisit distances, while the second method
continues hovering trajectories from short to longer revisit distances. In particular, new continuation method
is developed to meet the design requirements of this new scenario."

It is a formation-flying and proximity-operations paper. A deputy spacecraft hovers near a chief spacecraft on
the 9:2 southern L2 NRHO (the Gateway baseline orbit): the deputy must return to a fixed relative position rho
after one NRHO period (a 1:1 "teardrop" formation), with one impulse per revisit period. The design variable is
the initial relative velocity; the quantities reported are the revisit distance (1 km up to 50 km) and the
impulse (of order 1e-4 to 8e-2 m/s). It does not find or catalogue periodic orbits, cyclers or flyby
trajectories; the NRHO is taken as given.

## 2. Model (READ)

Earth-Moon circular restricted three-body problem, spatial, massless spacecraft, rotating frame, unit system
LU = 384405 km (the table gives 3.84405000e5 km), TU = T_EM / (2 pi). Equations of motion (1) to (3), Jacobi
constant (4). Integrator: MATLAB ode113 with absolute and relative tolerances 1e-13. Relative motion between the
deputy and the chief (an NRHO propagated in the same CR3BP) is propagated in the full nonlinear model (section
2.3, p. 6), with the linear model (state transition matrix, eqs (5) and (6)) used only for first guesses.

## 3. The "continuation method" (READ)

Section 3.3 (pp. 8-9). The parameter continued is the revisit distance rho (the length of the target relative
position), a parameter that appears in the boundary conditions of the problem, not in the dynamics. The
solution at step k is the initial relative velocity delta-v(t_0) that satisfies the revisit constraint psi = 0
(eq. 9: the relative position after one NRHO period equals the initial relative position). Step: rho^k = rho^{k-1}
+ Delta rho with Delta rho = 0.1 km, direction (alpha, beta) fixed (eq. 18). Predictor: the first-order Taylor
expansion of the constraints (19) gives a linear system A Delta(delta-v(t_0)) = b (20), where for this 1:1
formation A is a 3x3 block of the monodromy-type state transition matrix of the absolute hovering trajectory
(entries Phi_14 to Phi_36, eq. 21) and b is built from the same entries and the change in the initial relative
position (22). Prediction by least squares: Delta(delta-v(t_0)) = (A^T A)^{-1} A^T b (23). Corrector: the
predicted velocity is the initial guess for a differential correction, which MATLAB fmincon performs as a
nonlinear programme (minimise the norm of psi). Termination: constraint norm above tolerance 1e-9, or more than
499 steps (p. 9, Fig. 5 flow chart).

Is it new? The authors present it as new ("a new form of linear predictor", Remark 4, p. 10; highlights: "New
continuation method is developed to meet the design requirements of this scenario"). In substance it is a
predictor-corrector natural-parameter continuation with a tangent-style linear predictor from a first-order
sensitivity (state transition matrix) of the constraint, and a least-squares solve. The introduction (p. 2)
describes it as differing from earlier multi-body continuation work because the continued parameter "is part of
the constraints" rather than a parameter of the dynamics or states, and cites it as an alternative to
the natural parameter continuation of Singh et al. (ref. 33). There is no pseudo-arclength, no fold handling,
no stability or bifurcation analysis, no step-size control (fixed 0.1 km step), and no use of a Jacobian of the
full boundary-value problem including time or period. (INFERRED: this is an application of standard
predictor-corrector continuation with a first-order predictor, with the problem-specific part being that the
continued parameter enters the boundary condition.) The paper cites for the predictor Fu, Wu and Gong,
"Analytical nonlinear predictor for trajectory continuation in multibody models", J. Guid. Control Dyn. (2025)
(ref. 29, listed as 1-10, doi 10.2514/1.G008863), and Oshima (ref. 28) for continuation of periodic orbits;
neither is held.

## 4. Printed tables (READ, transcribed)

Table 1 (p. 4), "Parameter Setting for the Earth-Moon PCR3BP":
- mu = 1.21506683e-2 (Earth-Moon mass parameter, dimensionless)
- T_EM = 2.24735067e6 s (Earth-Moon period)
- R_E = 6378.145 km (mean Earth radius)
- R_M = 1737.100 km (mean Moon radius)
- LU = 3.84405000e5 km (length unit)
- TU = 3.75676968e5 s (time unit)

Table 2 (p. 5), "Initial States And Period Of The Considered NRHO" (the 9:2 NRHO):
- x_iNRHO = 0.987581435006489 LU
- z_iNRHO = 0.005276210630165 LU
- v_iNRHO = 2.120240531159090 LU/TU
- T_NRHO = 4 pi / 9 TU
(The text, p. 4, calls the NRHO "axisymmetric" (sic; the orbit is symmetric about the x-z plane). The initial
state is on the x-z plane with only a y-velocity; the paper labels the velocity component v_iNRHO. The
printed value of T_NRHO, 4 pi / 9 TU, is as printed; the paper does not give a numeric value in days beyond
"approximately 6-7 days" on p. 10, which is consistent with 4 pi / 9 TU times 4.348 d per TU.)

Table 3 (p. 11), "The relative states of the hovering trajectory with the minimum impulse at t = t_0":
- delta x(t_0) = 0 LU
- delta y(t_0) = -2.60142297836917e-6 LU
- delta z(t_0) = 0 LU
- delta u(t_0) = -3.2643727501816e-5 LU/TU
- delta v(t_0) = -1.98390221419e-7 LU/TU
- delta w(t_0) = 5.33425501523417e-4 LU/TU
Associated text (p. 10): "The minimum impulse required is 7.333e-4 m/s"; the minimum-impulse pair (alpha,
beta) is not tabulated, and rho = 1 km for this case.

Table 4 (p. 14), "Parameter settings for simulations" (fmincon): TolX = 1e-11, TolFun = 1e-11, TolCon = 1e-11,
MaxIter = 50000, MaxFunEvals = 50000.

Table 5 (p. 14), "Optimization variables ranges setting in the fmincon command": delta u(t_0), delta v(t_0),
delta w(t_0) each within the linear-model value +- 1.5 LU/TU.
(The Appendix text, p. 14, refers to "Table 4" for the parameter settings and "Table 5" for the ranges; the
captions printed are as transcribed here.)

No other tables exist in the paper.

## 5. Results (READ, short)

Short revisit distances (rho = 1 km): a grid over (alpha, beta) with step pi/100 (p. 10), impulse map in Fig. 6
(range roughly 0 to 80 m/s on the colour bar). The best case needs about 7.333e-4 m/s per revisit; its initial
relative state is almost along the eigenvector of the monodromy matrix with eigenvalue 1 (p. 11: "the impulse
calculated by the linear model is almost 0"), which the authors call a near-natural formation. Over ten revisit
periods the linear-model design drifts by tens of km (Fig. 10(a), up to about 85 km); the nonlinear design
drifts by order 1e-3 km (Fig. 10(b)). Long revisit distances: continuation from 1 km in steps of 0.1 km to 50 km
(Fig. 11, 12); the nonlinear impulses grow to about 0.07 m/s near rho = 10 km (Fig. 13(a), read from the
figure; approximate).

## 6. Verdict

Scope: OUT OF SCOPE for a catalogue of cycler, quasi-cycler and related periodic trajectories: it designs
near-periodic relative motion (a hovering formation) around a given NRHO and contains no new periodic orbit, no
cycler and no multi-body flyby sequence. The NRHO itself (Table 2) is a standard Earth-Moon CR3BP halo orbit
already catalogued in the periodic-orbit literature (Franz and Russell 2022 database, ref. 23, not held).

Does the continuation method offer anything for continuing a periodic orbit in a model parameter (#884, #890)?
No (INFERRED). It continues a constrained trajectory in a parameter of the boundary conditions with a first-order
sensitivity predictor and a MATLAB nonlinear-programme corrector, with a fixed step and no arclength
parametrisation, step control, fold handling or Floquet analysis; the project's continuation in a model parameter
(mass ratio) needs a predictor that stays well-posed at folds and a corrector on the periodic-orbit boundary
value problem including the period, which this paper does not treat. The one transferable idea, using the
state transition matrix block to build a least-squares tangent predictor, is standard textbook
sensitivity-based continuation and not new to this paper. Do not add to the catalogue and do not act on this paper.

# Digest: Guillaume 1975, "The restricted problem: an extension of Breakwell-Perko's matching theory"

P. Guillaume, "The restricted problem: an extension of Breakwell-Perko's matching theory", Celestial Mechanics 11:449-467 (1975), DOI 10.1007/BF01650284 (received 29 January 1974). Facultes Universitaires Notre-Dame de la Paix, Namur. Filed in the private paper corpus as
`guillaume-1975-restricted-problem-extension-breakwell-perko-matching-theory-celest-mech-11-449-doi-10.1007-BF01650284.pdf` (19 PDF pages; PDF page n is journal page 448 + n; a scan with an OCR layer; every page was read on a 130 dpi page image and nothing below was taken from the text layer).
Evidence tags: READ (journal page) is what the printed page says; COMPUTED is my own arithmetic of 2026-10-04; INFERRED is my reading across sources.
Companions read for this note: `docs/notes/2026-10-04-digest-perko-1976-second-species-O-mu-near-moon.md`, `docs/notes/2026-10-04-digest-guillaume-1973-periodic-symmetric-solutions-small-mu.md`, `docs/notes/2026-10-04-digest-henon-1968-consecutive-collision-orbits.md`, `docs/notes/2026-10-04-digest-bruno-1981-periodic-flybys-of-the-moon.md`, `docs/notes/2026-10-04-digest-henon-1997-generating-families.md`, `docs/notes/2026-10-04-digest-hitzl-henon-1977-critical-generating-orbits.md`.

## 0. What the paper is, and what it is not

This is the SECOND of two 1975 papers by Guillaume in the same volume. The first ("I" in the text, footnote p.449) is Guillaume, Celest. Mech. 11:213-254 (not held); Perko cites it as "1975a" and this one as "1975b". The paper is a derivation, not an application: it contains no table, no orbit, no value of mu, and no numerical example. Its content is the asymptotic matching of an outer (Earth-centred Kepler plus perturbations) and an inner (Moon-centred hyperbola plus perturbations) expansion of a second species solution of the planar-or-spatial circular restricted problem, with the orders of magnitude re-chosen so that the matching works when the incoming and outgoing Moon-relative velocities of the basic orbit are the same (V1' = V1, a zero zeroth-order turn). Breakwell and Perko's theory (1965 thesis; 1966 Progress in Astronautics and Aeronautics 17:160-182) assumed an O(mu) neighbourhood and could not describe that case (p.449).

It does NOT contain the coefficients of Guillaume 1973's local forms (the f1 f2 = mu hyperbola and the cubic). Section 5 explains why.

## 1. Setting and notation (READ pp.449, 451-452)

- Sidereal frame centred at the Earth E. r = ES, R = EM (Moon position), rho = MS = r - R, mu = m_M/(m_E + m_M). Equations (1.1), p.451:
  r'' = -r/r^3 + mu g(r, t),  g(r, t) = -[ rho/rho^3 + R/R^3 - r/r^3 ],  R'' = -R/R^3.
  At mu = 0 these are the two-body equations with Keplerian conics. (Units: Earth-Moon distance and the Earth gravitational parameter scaled to one, so that R'' = -R/R^3 is a unit circular orbit; the paper does not restate units; INFERRED from R'' = -R/R^3.)
- Basic orbit r^(0)(t): a Keplerian conic that meets the Moon at time t1, r^(0)(t1) = R(t1), with relative velocity V1 = r^(0)'(t1) - R'(t1) not equal to 0 (hypothesis H1, p.452). It does not meet the Moon in [t0, t1[. Its data at t0 are r0, r0'.
- The actual solution r(t, mu) has initial data r0 + dr0, r0' + dr0' (H2: dr0, dr0' at most O(mu^(1/3))) and a closest approach to the Moon at time t_p (H3: t_p - t1 less than O(mu^(1/3))). tau_p = t_p - t.
- Moon-relative state of the outer Kepler conic r(t) at t_p: rho_p, rho_p' (written rho_P in the paper). g1(rho, t) = -[ (rho + R)/|rho + R|^3 - R/R^3 ] is the Moon-removed perturbation (p.452).
- s1 = +1 or -1 as the basic orbit is an input conic (t0 < t1) or an output conic (t0 > t1) (p.450).
- Hyperbola data (p.455): F eccentric anomaly, n* mean motion, v_inf and Delta the velocity at infinity and the position vector of the asymptote with respect to the Moon (the arrival asymptote for s1 = +1, the departure asymptote for s1 = -1; Fig. 1, p.455). e* the eccentricity.

## 2. The outer expansion (READ pp.452-454)

r(t, mu) = r(t) + mu r^(1)(t) + mu^2 r^(2)(t) + r^(3)(t, mu)   (1.2), valid for t in [t0, t_p - K1 mu^(1/3)].

1. r(t) is a Keplerian conic with the same initial data as r(t, mu) at t0 (so it carries the dr0 variation); rho(t) = r(t) - R(t) expands in tau_p = O(mu^(1/3)) as
   rho(t) = rho_p - tau_p rho_p' + (tau_p^2/2!) g1(rho_p, t_p) - (tau_p^3/3!) g1'(rho_p, t_p) + (tau_p^4/4!) g1''(rho_p, t_p) + O(mu^(5/3))   (1.3).
2. mu r^(1)(t) is the first-order perturbation, the integral from t0 to t of Phi_rv(t, t') g(r(t'), t') dt' (Phi_rv the 3 by 3 matrix of partials of r(t) with respect to r'(t')). The part that is singular at the Moon is isolated:
   mu r_s^(1)(t) = -mu integral_{t0}^{t} (t - t') rho_s(t')/rho_s^3(t') dt'
   = -(mu s1/|rho_p'|^3) ln( gamma(t)/(2|rho_p'|^2) ) rho_p' - mu ( (rho_p' x rho_p) x rho_p' )/( |rho_p'|^3 gamma(t) )   (1.4),
   with gamma(t) = |rho_p'| |rho_s(t)| - s1 rho_p' . rho_s(t) and rho_s(t) = rho_p - tau_p rho_p' (the straight-line approximation of the relative path). (The printed vector product reads "(rho_p' X rho_p) X rho_p'" with the cross written as an X.)
   The remainder mu r_b^(1)(t) = mu r_bp - tau_p mu r_bp' + O(mu^(5/3))   (1.5), where r_bp and r_bp' are the regular parts, limits as t tends to t1:
   r_bp = lim { integral_{t0}^{t} Phi_rv^(0)(t, t') g(r^(0)(t'), t') dt' + (mu s1 V1/V1^3) ln|t1 - t| } + O(mu^(4/3)),
   r_bp' = lim { integral_{t0}^{t} Phi_vv^(0)(t, t') g(r^(0)(t'), t') dt' - (mu s1 V1/V1^3)(t1 - t)^(-1) } + O(mu^(4/3)).
   (Phi^(0) are the partials of the basic orbit; as printed the second term carries V1 over V1^3 with the vector V1 on top; I transcribed it as printed.) The expansion is valid only if rho_p . rho_p at most O(mu^(1/3)) (H4, p.454), equivalently (t1 - t_p) = V1 . dr1/V1^2 + O(mu^(1/3)) (H4'), with dr1 = r(t1) - r^(0)(t1).
3. mu^2 r^(2)(t) at most O(mu^(5/3)) (1.6) and r^(3)(t, mu) at most O(mu^(7/3)) (1.7) for |tau_p| at least O(mu^(1/3)), "with methods quite similar to Perko's (1965)". Perko assumed |t1 - t| at least O(mu^(1/2)) and |t_p - t1| at most O(mu); the outer region here (|t_p - t| at least O(mu^(1/3))) is strictly contained in his for small mu (p.454).

## 3. The inner expansion (READ pp.454-457)

Near the Moon, rho'' = -mu rho/rho^3 - (1 - mu)(r/r^3 - R/R^3) (p.454); the first term dominates for small rho. rho(t, mu) = rho^(0) + rho^(1) + rho^(2)   (1.8).

1. rho^(0) is the Moon-centred hyperbola osculating rho(t, mu) at t_p (rho'' = -mu rho/rho^3), written (1.9) to (1.12), pp.455:
   rho^(0)(t, mu) = -(s1 mu/v_inf^3) e* sh|F| v_inf + ( s1 mu v_inf/v_inf^3 + Delta )( 1 - e*^(-1) exp(-|F|) )   (1.9)
   |rho^(0)(t, mu)| = (mu/v_inf^2)( e* ch|F| - 1 )   (1.10)
   s1 tau_p v_inf^3/mu = -n* s(F) tau_p = e* sh|F| - |F|   (1.11)
   e*^2 = 1 + Delta^2 v_inf^4/mu^2   (1.12)
   (sh, ch are sinh, cosh; s(f) = +1 or -1 as f is positive or negative; the printed dagger footnote defines s(f); in (1.9) the bold v_inf in the second bracket is the vector, as read from the image.)
   Assumptions (p.455): v_inf = O(mu^0) (H5); O(mu^(2/3)) at most Delta at most O(mu^(1/3)) (H6), equivalent to O(mu^(2/3)) at most mu e* = O(Delta) at most O(mu^(1/3)) (H6') and to O(mu^(2/3)) at most rho(t_p, mu) at most O(mu^(1/3)) (p.457). For |tau_p| at most O(mu^(1/3)): mu e* exp|F| at most O(mu^(1/3)) (p.456).
2. rho^(1) = (1 - mu) integral from t_p to t of Phi_{rho rho'}(t, t') g1[rho^(0)(t', mu), t'] dt' (p.456). Perko's split Phi = A1 + A2 with A1 = (t - t') I3 and A2 at most O(mu^(4/3)/rho^(0)); A2 terms are at most O(mu^(5/3)) and dropped. With rho^(0)(t') replaced by Delta - tau_p' v_inf (negligible change):
   rho^(1)(t, mu) = integral (t - t') g1[Delta - tau_p' v_inf, t'] dt' = (tau_p^2/2!) g1(Delta, t_p) - (tau_p^3/3!) g1'(Delta, t_p) + (tau_p^4/4!) g1''(Delta, t_p) + O(mu^(5/3))   (1.13).
3. rho^(2) is discussed "in great detail by Perko (1965)" and is not reproduced. The paper notes that without it, the condition (H6) is necessary (p.457).

## 4. The matching (READ pp.457-462) and what changes against Breakwell-Perko

### 4.1 The first-stage matching (hypotheses 0.1, 0.2), pp.457-459

Constants K1 and Kp can be chosen so that the outer interval [t0, t_p - K1 mu^(1/3)] and the inner interval [t_p - Kp mu^(1/3), t_p] overlap. In the transition interval both expansions are written in a common independent variable (tau_p or F; F is preferred to avoid inverting Kepler's equation). Keeping terms larger than O(mu), outer and inner reduce to Delta - tau_p v_inf + O(mu) and rho_p - tau_p rho_p' + O(mu); comparing the tau_p^0 and tau_p^1 (= O(mu^(1/3))) coefficients gives (1.14):
   rho_p = Delta + O(mu),  rho_p' = v_inf + O(mu^(2/3)).
(The matching equations "imply" H4, H5, H6 as consequences of H1 and H2, with the word "imply" too strong, since all are used to build the expansions, p.457.) The terms (tau_p^k/k!) g1^(k-2)(Delta, t_p) and the same with rho_p are identified to O(mu^(5/3)) (p.458). The finer analysis expresses mu r_s^(1) in the inner variables (Delta, v_inf, F, e*) using (1.14):
   -rho_s(t) = Delta - tau_p v_inf + O(mu) = rho^(0)(t, mu) + O(mu) = (mu e* ch|F|/v_inf^2) + O(mu),
   -rho_p' . rho_s = -v_inf^2 tau_p + O(mu) = -(s1 mu e* sh|F|/v_inf) + O(mu),  |rho_p'| = v_inf + O(mu^(2/3)),
   -gamma(t) = (mu e* exp|F|/v_inf)(1 + O(mu^(2/3))),
   -ln( gamma(t)/(2|rho_p'|^2) ) = ln( mu e*/(2 v_inf^3) ) + |F| + O(mu^(2/3)),
   mu r_s^(1)(t) = -(mu s1 v_inf/v_inf^3){ ln( mu e*/(2 v_inf^3) ) + |F| } - Delta/(e* exp|F|) + O(mu^(5/3)),
(vector v_inf in the first term, as read.) Dropping everything at most O(mu^(5/3)), including the term -s1 mu v_inf/(v_inf^3 e* exp|F|) of rho^(0), the difference of outer and inner must vanish identically in tau_p; the coefficient of tau_p^0 and tau_p^1 give, p.459, the MATCHING EQUATIONS (1.15):
   rho_p + mu r_bp = Delta + (s1 mu v_inf/v_inf^3)[ 1 + ln( mu e*/(2 v_inf^3) ) ] + O(mu^(5/3)),
   rho_p' + mu r_bp' = v_inf + O(mu^(4/3)).
(Here v_inf in the numerator is the vector and v_inf^3 the cube of its magnitude; the text also prints V_inf with a capital in (2.8), a typesetting variation.) They are identical to Breakwell and Perko's equations (0.3) of this paper's introduction (their form (1.32) in paper I),
   dr1^(0) + V1 (t_p - t1) + mu r_b1 = Delta + (s1 mu/v_inf^3)(1 + ln(mu e/(2 v_inf^3))) v_inf + O(mu^2),  V1 + dr1'^(0) + mu r_b1' = v_inf + O(mu^(3/2))   (0.3, p.450),
in the form (0.4) where rho_P = dr1^(0) + V1 (t_p - t1) and rho_P' stand for the Moon-relative state of the outer conic extrapolated to t_p. Only the error terms differ: O(mu^2) and O(mu^(3/2)) in Breakwell-Perko against O(mu^(5/3)) and O(mu^(4/3)) here (p.450). Note (0.3) prints e inside the logarithm where (1.15) prints e*; the two are the same eccentricity (the paper's own variation in notation).

### 4.2 The general theory (Section 2), pp.458-462

Let delta be the order of tau_p in the transition interval, |K1| delta at most |tau_p| at most |Kp| delta (2.1), and epsilon the order of Delta, t1 - t_p, dr0, dr0'. The paper shows in that interval (p.458):
   (1) mu r_b^(1)(t) = mu r_bp - tau_p mu r_bp' + O(mu delta^2);
   (2) mu^2 r_b^(2)(t) at most O(mu^2);
   (3) r^(3)(t, mu) at most O(mu^3/delta^2);
   (4) inner perturbations beyond first order at most O(delta^5).
With mu delta^2 as the first order to be neglected: (1)+(2) need delta at least O(mu^(1/2)); (4) needs delta at most O(mu^(1/3)). Hence
   O(mu^(1/2)) at most delta at most O(mu^(1/3))   (2.2).
Breakwell and Perko took the lower limit (so mu delta^2 = mu^2) and had to retain more terms; this paper took the upper limit (mu delta^2 = mu^(5/3)), a "limiting case" in which the bounded second-order term is negligible (p.460). For delta = O(mu^alpha), alpha above 1/3 (delta smaller than mu^(1/3)), the singular second-order term mu^2 r_s^(2) = O(mu^2/delta) exceeds mu delta^2 and must be kept; that is the content of Section 2. The orders of dr0, dr0', Delta and t1 - t_p must satisfy
   O(delta^2) at most epsilon at most O(delta^4/mu)   (stated p.451, derived (2.3), p.461).
Special cases (p.451): delta = O(mu^(1/3)) gives epsilon = O(mu^(1/3)); delta = O(mu^(1/2)) gives epsilon = O(mu). The error terms of (0.4) become O(mu delta^2) and O(mu delta) respectively instead of O(mu^2) and O(mu^(3/2)).

The singular second-order term (p.460-461), with Gamma(rho) = -I3/rho^3 + 3 rho rho^T/rho^5:
   mu^2 r_s^(2)(t) = mu integral_{t0}^{t} (t - t') Gamma[rho_s(t')] [ mu r_s^(1)(t') + mu r_bp ] dt',
   which, using rho_s(t') = |rho_p'||tau_p'|(1 + O(epsilon/delta)) and mu r_s^(1)(t') = -(mu s1/|rho_p'|^3) ln|tau_p'| + O(mu epsilon/delta), integrates to (2.4)
   mu^2 r_s^(2)(t) = -(mu^2 s1/(V1^6 |tau_p|)) [ ln|tau_p| + 3/2 ] V1 + (mu^2/(2 V1^5 |tau_p|)) [ -V1^2 r_bp + 3 (V1 . r_bp) V1 ] + O(mu^2 epsilon/delta^2).
   (Printed V1^6 and V1^5 are powers of the speed; the first term has the vector V1 as the last factor, as read.)
   Requiring the omitted term at most O(mu delta^2) gives epsilon at most O(delta^4/mu), (2.3).
Second-stage matching (terms not at most O(mu delta)): the identification gives (2.5)
   rho_p + mu r_bp = Delta + (s1 mu v_inf/v_inf^3)[ 1 + ln(mu e*/(2 v_inf^3)) ] + O(mu delta),  rho_p' = v_inf + O(mu).
Inner variables of gamma, mu r_s^(1) and mu^2 r_s^(2) through (2.5) are printed as (2.6) and (2.7), p.462 (long expansions in |F|, e* exp|F| and the logarithm; every term is read but they are intermediate and I do not transcribe them; the full matching is "of an easy algebraic type but too long to be reproduced", p.462). The result, p.462, is that the matching reproduces (1.15) with the error terms replaced:
   rho_p + mu r_bp = Delta + (s1 mu v_inf/V_inf^3)[ 1 + ln( mu e*/(2 V_inf^3) ) ] + O(mu delta^2),   rho_p' + mu r_bp' = v_inf + O(mu delta)   (2.8).
"The result is remarkable: we find again Equations (1.15)" with only the error orders changed.

### 4.3 What it extends beyond Breakwell-Perko 1965/1966

(a) The transition interval is O(delta) with mu^(1/2) at most delta at most mu^(1/3), where Breakwell-Perko fixed delta = mu^(1/2) (the O(mu) neighbourhood of the Moon). (b) The size of the allowed variations: dr0, dr0', Delta, t1 - t_p up to O(epsilon) with epsilon up to delta^4/mu (mu^(1/3) at the upper limit), where Breakwell-Perko had O(mu). This is the point: with Delta = O(mu^nu) rather than O(mu) the zeroth-order turn can be zero (V1' = V1) while a nonzero turn of order mu^(1 - nu) arises from the hyperbola. (c) The matching equations keep the same form: the claim of Section 2 is that the first-order matching conditions (1.15) hold unchanged on the whole range, only the error terms change (O(mu^2), O(mu^(3/2)) become O(mu delta^2), O(mu delta)). (d) A consequence the paper states (p.457): the matching equations make H4 to H6 consequences of H1 and H2. Perko's Part C of the Perko digest uses exactly this theory (his "eq. 5, from Guillaume's eq. 2.8", with the "epsilon^2" footnote citing p.451).

Limitations stated or visible: the appendix assumes rho_p . rho_p' = o(delta) (A1.7, A1.8; "perhaps rather restrictive but it does not introduce incompatibilities", p.464), the outer region needs r^(0) bounded away from the Earth (A1.9 and following, p.464), and mu small enough that the Earth-centred conic stays away from the Earth (p.464). The remainder estimates for r^(3) and rho^(2) are quoted from Perko (1965 thesis, not held), not proved here.

## 5. Do the coefficients of Guillaume 1973's local forms appear here? No.

Checked against the 1973 digest (sections 4.1 and 4.2 there) and against this paper's content:
- The 1973 equation (2), f1 f2 = mu + O(mu^(3/2)), with f1 = 0 and f2 = 0 the curves A0 and I_r in the (x0, C) plane, and the cubic (4'') with its coefficients alpha, beta, alpha', beta', gamma, lambda, nu, a, b, delta, are statements about the SYMMETRIC conjunction conditions y = xdot = 0 written for perturbed initial conditions (x0 + dx0, ..., C + dC) in the SYNODIC frame. This paper never uses the synodic frame, the Jacobi constant, or the symmetry conditions; it works in the sidereal Earth-centred frame for a single arc ending at the Moon, with dr0, dr0' arbitrary. No equation here has the form f1 f2 = mu and the word cubic does not occur.
- Perko's Part A (p.419, in the Perko digest) puts the second species bifurcations under "Guillaume (1975a), pp. 253-254", i.e. the last pages of paper I (Celest. Mech. 11:213-254), not this paper. So the application of the matching to the bifurcations at (x0, C) = (-1, -1) and the n = -2 ellipse, with the printed coefficients, is expected at I, pp.253-254. That paper is NOT held. (INFERRED from Perko's citation; not verified.)
- What this paper does supply is the matching equations (2.8) in the form the symmetry equations of the 1973 local analysis need as input: if one writes rho_p, rho_p' through the outer conic's variations (dr0, dr0', t_p - t1, as in (0.4)) and eliminates Delta, v_inf by (2.8) and (1.12), the 1973 coefficients follow by algebra at each junction. I did not carry the algebra out; it is the content of paper I and of Perko 1981 (eq. 9 there). So the 1973 digest's open item 1 ("the coefficients missing from eqs. 3'' and 4''") is NOT settled; what changes is where to look (paper I, pp.253-254) and that the matching these equations rest on is now in hand.
- One structural fact that this paper does settle for the 1973 "O(mu^(4/3))" error order (1973 digest section 4.2, flagged there as unexplained): (1.15) has error O(mu^(5/3)) in the position equation and O(mu^(4/3)) in the velocity equation under hypotheses H1 to H6. The 4/3 in the cubic's (4'') is plausibly this same velocity-equation error (INFERRED, not stated by the author; it supports, not proves, the reading in the 1973 digest that the cubic scaling is mu^(1/3)).

## 6. The relation r_p = mu (e - 1)/v^2 with sin(turn/2) = 1/e: unchanged

(1.10) at F = 0 (periapsis, tau_p = 0): |rho^(0)| = (mu/v_inf^2)(e* - 1). With (1.12), e*^2 - 1 = (Delta v_inf^2/mu)^2, so cot(turn/2) = sqrt(e*^2 - 1) = Delta v_inf^2/mu, equivalent to tan(turn/2) = mu/(Delta v_inf^2) and sin(turn/2) = 1/e* (the turn angle itself is not named in the paper; this is the standard identity, and it is the same as Perko's eq. 3.8 in the Perko digest). Hence r_p = mu (e* - 1)/v_inf^2 with sin(turn/2) = 1/e* is exact for the osculating inner hyperbola (1.9) to (1.12).

COMPUTED check that I ran on the printed formulas: take Delta perpendicular to v_inf, build the vector (1.9) in components (along v_inf: (s1 mu/v_inf^2)(1 - e^(-1) exp(-|F|) - e sh|F|); perpendicular: Delta (1 - e^(-1) exp(-|F|))) and compare its norm to (1.10) at F = 0, 0.7, -1.3, 2.0, for s1 = +1 and -1, for (v_inf, Delta) = (0.5, 0.02), (1.0, 0.1), (0.3, 0.005) in the normalised units at mu = 0.0121529529: agreement to 1e-16 absolute in all 24 cases. The same cases give turn 135.273 deg (e* = 1.08133, r_p = 0.0039535), 13.858 deg (e* = 8.28899, r_p = 0.0885828) and 175.759 deg (e* = 1.000685, r_p = 9.2538e-5), with 2 arctan(mu/(Delta v^2)) equal to 2 arcsin(1/e*) to 1e-14 degrees. So (1.9), (1.10) and (1.12) are mutually consistent as transcribed, which also checks my transcription of (1.9).

ANSWER: the matching changes only the map from the Earth-side initial conditions (dr0, dr0', t_p - t1, with the logarithmic regular part r_bp) to the hyperbola parameters (Delta, v_inf, and via (1.11) the perilune time); it does not change the turn-gate relation. This confirms the Perko digest's conclusion. Two additions, from this paper only:
1. v_inf is not |rho_p'| exactly: (2.8) gives rho_p' + mu r_bp' = v_inf + O(mu delta), so the Moon-relative speed of the outer conic at the extrapolated periapsis differs from the hyperbola's v_inf by O(mu) (the r_bp' term), which has to be carried when the gate is evaluated from outer-conic states at mu = 0.0121529529 (an O(mu) = 1.2 percent shift of a speed of order one, COMPUTED size, not a coefficient).
2. Range of validity of the first-stage matching for Delta: H6 says Delta between O(mu^(2/3)) and O(mu^(1/3)), so the turn from tan(turn/2) = mu/(Delta v^2) is between O(mu^(1/3)) (Delta at the lower end) and O(mu^(2/3)) (Delta at the upper end). COMPUTED scales at mu = 0.0121529529 and 384,400 km: mu^(2/3) = 0.05286 (20,319 km), mu^(1/2) = 0.11024 (42,376 km), mu^(1/3) = 0.22991 (88,378 km), mu = 4,672 km, mu |ln mu| = 0.0536. These are scales with unknown O(1) coefficients, not boundaries. The general theory of Section 2 relaxes the bound on Delta through epsilon (O(delta^2) at most epsilon at most O(delta^4/mu)); the paper does not restate H6 for it (INFERRED that Delta follows epsilon).

## 7. Printed numbers usable as sourced tests

There are no printed numerical values (no tables, no example orbit). The printed content that can be tested is identities and exponents; every item is a published statement, not a value produced by project code.

| Item | Expected (as printed) | Page | Test form |
|---|---|---|---|
| Hyperbola magnitude (1.10) against vector form (1.9) | equal for all F, s1 | 455 | Identity. Verified by the arithmetic in section 6 (COMPUTED here, to 1e-16); make it a permanent unit test of a hyperbola helper with random (v_inf, Delta, F, s1). Source for the expected side: the printed (1.10), not the project's own function. |
| Eccentricity (1.12) | e*^2 = 1 + Delta^2 v_inf^4/mu^2 | 455 | Test against the gate: e* from (1.12) gives sin(turn/2) = 1/e* and r_p = mu (e* - 1)/v_inf^2 to rounding; ties the matching's hyperbola to `verify/turn_gate.py` `required_periapsis_alt_km`. |
| Kepler equation (1.11) | s1 tau_p v_inf^3/mu = e* sh|F| - |F| | 455 | Identity for the time on the hyperbola; test by integrating the two-body hyperbola numerically and comparing t(F). |
| Admissible transition width (2.2) | O(mu^(1/2)) at most delta at most O(mu^(1/3)) | 459 | Structural: the exponents. |
| Admissible variation size (2.3) | O(delta^2) at most epsilon at most O(delta^4/mu); delta = mu^(1/3) gives epsilon = mu^(1/3); delta = mu^(1/2) gives epsilon = mu | 451, 461 | Exponent check of (2.3): the two corner cases hold exactly (COMPUTED: delta^4/mu at delta = mu^(1/3) is mu^(1/3), and at mu^(1/2) it is mu). |
| Matching error orders | (1.15): O(mu^(5/3)), O(mu^(4/3)); (2.8): O(mu delta^2), O(mu delta); Breakwell-Perko (0.3): O(mu^2), O(mu^(3/2)) | 450, 459, 462 | Convergence slopes on log-log in mu of the residual of (1.15) against a full CR3BP integration of a flyby with v_inf, Delta fixed: expected slope 5/3 for position, 4/3 for velocity in the stated range; "at least" rather than equality, so a steeper slope is not a failure. |
| Validity hypotheses | H4: rho_p . rho_p at most O(mu^(1/3)); H5: v_inf = O(1); H6: Delta in [mu^(2/3), mu^(1/3)] | 454, 455 | Hypothesis checks before applying (1.15): the same kind of check the Hitzl-Henon digest proposed for the gate's validity numbers. |
| Outer-inner regions overlap | t_p - t order mu^(1/3) at the end points | 457 | Structural. |

A control that no longer rests on an invented number: the exact Kepler hyperbola with the project's own `required_periapsis_alt_km`; the paper's printed relation (1.10) at F = 0, i.e. r_p = mu (e* - 1)/v_inf^2, with e* from (1.12) is an independent source.

## 8. Reconciliation against project code

- Turn gate (`src/cyclerfinder/verify/turn_gate.py`): demanded_turn = arccos(v_in . v_out/(|v_in||v_out|)); available bend = 2 asin(1/(1 + r_p v^2/GM)); required altitude r_p = GM/v^2 (1/sin(demanded/2) - 1) minus R (module docstring lines 20-36; `required_periapsis_alt_km` at line 283 and `available_bend_rad` at line 278). With e = 1 + r_p v^2/GM, sin(turn/2) = 1/e, this is exactly (1.10) at F = 0 with GM = mu and v = v_inf, so the code agrees with the printed relation; nothing in this paper requires a change (section 6). The #888 wiring (`tests/verify/test_888_turn_gate_wiring.py`, `src/cyclerfinder/data/validation/moontour_turn.py`) calls the gate on V-infinity vectors at the encounters and is unaffected.
- Second species and matching: a grep of `src/cyclerfinder` and `scripts` for second-species or matching code finds only docstrings (`search/earth_moon_resonant_families.py`, `search/earth_moon_class1_resonant_connections.py`, the literature-check query string) saying the project does not implement Casoliva's second-species method; no outer or inner expansion, hyperbola matching map or Breakwell-Perko code exists. Nothing to reconcile for the matching; the paper is input to a new module.
- Not adopted anywhere else: (1.12) and (1.11) have no counterpart in the project; `core/flyby.py` has the bend relation only.

## 9. Techniques applicable to the project's problems

Everything below is INFERRED use; nothing was built or run.

**#899 (second species continuation to the Earth-Moon mass; reproduction of Casoliva 2010 with the Henon enumerator; read data/OUTSTANDING.md entry, including the 2026-10-04 controls and first build step).** The OUTSTANDING entry makes #899 a reproduction of Casoliva et al. 2010 section IV.C (second-species seeds at mu = 1e-6 from Barrabes and Gomez's matched maps, three-step continuation to the Earth-Moon mass), with the nine catalogued Casoliva rows as the control. This paper contributes the in/out map that a native seed generator would use, and the validity statement that tells the continuation where the seed is trustworthy.
1. Build `hyperbola_match(v_inf, Delta, mu, s1)` returning e*, perilune distance, perilune time and the turn from (1.9) to (1.12), with the outer-conic side supplied as (rho_p, rho_p') extrapolated to t_p. Test: the identities of section 7 on random inputs (positive control, expected side from the printed (1.10), not from the helper), plus a round trip (v_inf, Delta) to turn to `required_periapsis_alt_km` to r_p.
2. Seed generator for one flyby at small mu: choose the outer conic (r0, r0') with an O(mu^nu) intended Delta, evaluate (0.4) / (2.8) to produce (Delta, v_inf), then the outgoing conic from the hyperbola's outgoing asymptote (same v_inf, deflected by the turn about the Delta direction); r_bp is the regular part of the first-order integral (1.5), computed once per basic orbit by quadrature with the logarithm subtracted. Control: integrate the CR3BP (Sun-less planar, the project's existing integrator) from the seed at mu = 1e-6 and 1e-5 and compare the matched pair of conics with the numerical flyby; expected residual slopes in mu of 5/3 and 4/3 (section 7 row 6; an observed slope below those is a bug, one above is allowed). Positive control for the whole generator: reproduce the 2008 conference seeds at 1e-6 listed in the OUTSTANDING entry (the Casoliva seeds), by the criterion that the seed closes under the project's corrector, not by comparing digits.
3. Continuation guard: record nu_eff = ln(r_p/a_moon)/ln(mu) and (mu |ln mu|)/v_inf^3 along the continuation (this extends the regime classifier of the Perko digest, Part C follow-up 1). This paper's range statements mean the first-order matching applies for nu_eff between 1/3 and 2/3 with the delta bounds of (2.2); outside it, use Perko's Part A (large turn) or a regular perturbation. At the Earth-Moon mass the scales are in section 6; a continuation that reaches a passage under 20,000 km (mu^(2/3)) leaves the H6 range at the lower end, which is the "impact" end of the Casoliva continuation (the OUTSTANDING entry records lunar impact as the usual end of the fixed-period continuation).
4. The Casoliva three-step strategy (raise the periselene by Jacobi constant at fixed mass) can be read as moving Delta upward through the H6 window; a test that nu_eff decreases monotonically in the second step is a cheap sanity check on a successful continuation.
5. Not provided: the f1 f2 = mu coefficients (section 5), so the local-hyperbola fit remains a data-driven fit, as in the 1973 digest follow-up 3.

**#906 (turn-gate hardening; small turns near resonance are indeterminate, not rejected).** The paper shows that the gate relation is exact and that the matching, not the gate, carries the O(mu) uncertainty. Concrete implementation, consistent with the OUTSTANDING amendment (return "indeterminate: near resonance or outside first-order validity", never reject on that ground alone):
1. In the "indeterminate" branch, report the three validity numbers: mu |ln mu|/v_inf^3, r_p/sqrt(mu) (as in the amendment), and the new one from H6: r_p against mu^(2/3) and mu^(1/3), i.e. nu_eff of section 9 above. A turn below the O(mu) uncertainty of the v_inf speed (item 1 of section 6, an O(mu) shift of the speed) cannot be classified as zero from the vectors alone.
2. Test (positive control from this paper, not from the code): for the three (v_inf, Delta) cases in section 6, build a synthetic encounter with v_in, v_out equal in magnitude and the printed turn (135.273, 13.858, 175.759 deg, which are COMPUTED from the printed (1.12) and the identity), and confirm the gate returns required periapsis equal to (1.10)'s r_p to 1e-12. A second test: perturb |v_out| by 1.2 percent (mu) and check the gate flags the magnitude mismatch separately (the module already reports it) and the demanded turn changes by less than the amount implied by the O(mu) speed shift.
3. Zero-turn junction: the paper shows a zero zeroth-order turn is exactly the case it was written for (V1' = V1); a positive-definite turn of order mu^(1 - nu) arises. So "indeterminate" must not be "turn is zero": the gate should print the expected scale mu/(Delta v^2) for the observed Delta.

**#890 and #895 (Titania-Oberon, Uranian; methods only).** The matching is for one small mass; for two moons the outer expansion is the same Kepler (Uranus-centred) conic with perturbations from both moons, and the inner hyperbola about each moon follows (1.9) to (1.12) with that moon's mu. The method transfers if the two moons' encounter windows do not overlap (|t_p - t| of order mu_moon^(1/3) much less than the synodic gap). The mass ratios of Titania and Oberon are about 1e-5 to 1e-4 of Uranus (the project's `body_constants`, not checked here), so mu^(1/3) is about 0.02 to 0.05; the validity ranges of section 6 scale accordingly. Concrete use for #890: run the identity tests on the Uranian mu values as a unit test of the same helper, and use nu_eff at each encounter of the patched-conic closure to say which regime (large turn, small turn, regular perturbation) the closure sits in. For #895: the one-moon Titania-only and Oberon-only second-species maps are the controls for any two-moon construction, since each should reproduce the matching with the other moon removed. No printed two-moon data exist in this paper.

## 10. Recommended follow-ups (not registered)

1. Obtain Guillaume, Celest. Mech. 11:213-254 (1975, paper I), pp.253-254 especially: it should carry the second species bifurcation applications and the coefficients of the 1973 local forms (section 5). Check CORPUS_INDEX for it before ordering.
2. Implement `hyperbola_match` (section 9, item 1) with the section 7 identity tests; this is a prerequisite of both the #899 seed generator and the #906 validity numbers.
3. Measure the residual slopes of (1.15) against a numerical CR3BP flyby at mu = 1e-6, 1e-5, 1e-4, 1e-3 (expected 5/3 and 4/3 or steeper); a clean slope is a control for the seed generator and a check that the transcription here is correct.
4. Add nu_eff and the H6 comparison to the Perko digest's regime-classifier follow-up.
5. Derive the 1973 f1 f2 = c mu hyperbola by combining (2.8) with the Earth-side variational matrices for the Moon-to-Moon case; compare c with the digitised Fig. 4 of the 1973 paper and with Perko 1981 eq. 9.
6. Check the printed (2.6) and (2.7) term by term against a symbolic expansion before using them; they were read but not transcribed or cross-checked here.
7. Obtain Perko 1965 thesis and Perko 1976 Rocky Mountain J. Math. 6:130-145 for the rho^(2) and r^(3) remainder estimates the paper cites without proof.

## 11. Illegible or ambiguous items

- (2.6) and (2.7), p.462: legible but not transcribed (intermediate).
- (1.9) and (2.8): the bold/italic distinction between the vector v_inf and the speed V_inf is read from the image; the capital V_inf in (2.8) and the lower-case v_inf elsewhere I take to be the same speed.
- (0.3): the logarithm prints e where (1.15) prints e*; taken as the same quantity.
- p.457: "rho^(0)(t_1 mu)" in the matching paragraph is a missing comma on the page (t, mu).

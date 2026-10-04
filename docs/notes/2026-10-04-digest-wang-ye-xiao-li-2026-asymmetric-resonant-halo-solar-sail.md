# Digest: Wang, Ye, Xiao & Li (2026), "Design and analysis of asymmetric resonant Halo orbits for solar sail spacecraft" (太阳帆航天器非对称共振晕轨道设计与分析)

Wang Dawei (王大维), Ye Dong (叶东, corresponding author, yed@hit.edu.cn), Xiao Yan (肖岩), Li Haiyang (李海洋). State Key
Laboratory of Micro-Spacecraft Rapid Design and Intelligent Cluster, Harbin Institute of Technology; Deep Space Exploration Lab,
Beijing. Acta Aeronautica et Astronautica Sinica (航空学报), Vol. 47, No. 6, article 332576, 25 March 2026 issue; DOI
10.7527/S1000-6893.2025.32576. Received 2025-07-16, revised 2025-08-21, accepted 2025-10-09, online 2025-10-28. 27 pages (332576-1
to 332576-27; the last page is the English title page and abstract). The text is IN CHINESE; captions and the title block are
bilingual.
Filed in the private paper corpus as
wang-ye-xiao-li-2026-asymmetric-resonant-halo-orbits-solar-sail-acta-aeronautica-astronautica-sinica-47-6-332576-doi-10.7527-S1000-6893.2025.32576-chinese.pdf

Digested 2026-10-04 from the page images (all 27 pages) with the PDF text layer used to check digits and to take the Chinese
quotations. Marks: READ (seen on the page; page numbers are the article page suffixes 1 to 27), TRANSLATED (the English is MY
translation of the Chinese quoted just before it; not an official translation), COMPUTED (our arithmetic) or INFERRED (our reading).
Context: tasks #893 and #884; `core/er3bp.py`; the Sun-Mercury elliptic problem digest `2026-06-25-digest-peng-2017-sun-mercury-ERTBP.md`.

## 0. What the paper is

Abstract, READ (p1) in Chinese: 基于太阳-水星椭圆限制性三体问题这一动力学模型，提出3种递进式的太阳帆指向律，使太阳帆平面由最初正对太阳，过渡到在引力平面内偏离太阳方向，再到最终偏离引力平面。利用多段打靶法进行轨道延拓计算，获得了具有不同共振比的太阳帆共振晕轨道。 TRANSLATED: "Based on the Sun-Mercury elliptic restricted three-body problem as the dynamical model, three progressive solar-sail steering laws are proposed, taking the sail plane from initially facing the Sun, to deviating from the Sun direction inside the gravitational plane, to finally deviating out of the gravitational plane. Using multiple shooting for orbit continuation, solar-sail resonant halo orbits with different resonance ratios are obtained."

Chinese: 结果表明，前2种策略下得到的共振晕轨道在旋转坐标系中关于XZ平面呈对称结构，而第3种策略下获得的轨道，不关于任意平面或直线对称，表现出复杂的空间三维结构，同时保持良好的周期性。 TRANSLATED: "The results show that the resonant halo orbits obtained under the first two strategies are symmetric about the XZ plane in the rotating frame, while the orbits obtained under the third strategy are not symmetric about any plane or line, show a complex three-dimensional spatial structure, and at the same time keep good periodicity."

Chinese: 最后在高精度星历模型下对轨道稳定性进行分析，结果显示所构建的太阳帆共振晕轨道在无需轨道保持的情况下，可在约4~6个水星轨道周期内（约352~528 d）维持轨道形状且不发散。 TRANSLATED: "Finally the orbit stability is analysed in a high-precision ephemeris model; the constructed solar-sail resonant halo orbits can keep their shape without diverging for about 4 to 6 Mercury orbital periods (about 352 to 528 d) with no orbit maintenance."

The paper prints no initial condition, period, Jacobi constant or stability index for any orbit. Its tables are the model parameters
(table 1) and three continuation-limit tables (tables 2 to 4). All orbits appear only in figures.

Translation slip to note: the English abstract (p27) says the third law varies the sail's "pitch angle"; the Chinese abstract and
body call it the azimuth angle (方位角, beta) throughout. They refer to the same angle (INFERRED).

## 1. The model exactly as printed

### 1.1 Frame and units (READ, p2, section 1, eq. 1, fig. 1)

Chinese: 该坐标系的原点位于太阳（第一引力体）与水星（第二引力体）的质心，OX轴指向水星，OY轴位于水星的轨道面内垂直于OX轴，指向水星速度方向，OZ轴构成右手系且指向水星的角动量方向。 TRANSLATED: "The origin of this frame is at the barycentre of the Sun (first primary) and Mercury (second primary); OX points to Mercury, OY lies in Mercury's orbital plane perpendicular to OX and points in the direction of Mercury's velocity, and OZ completes a right-handed system and points along Mercury's angular momentum."

Barycentric pulsating-rotating frame OXYZ. Unit length "DU" is the instantaneous Sun-Mercury distance,
`DU = a(1 - e^2)/(1 + e cos f)` (eq. 1; a semi-major axis, e eccentricity, f Mercury's true anomaly); unit time `TU = DU^2/h` with h
Mercury's angular momentum (so that f is the independent variable and d f/d t = 1 in the normalised time); unit mass the sum of the
two masses, `mu = m2/(m1 + m2)`; Sun mass 1 - mu at `[-mu, 0, 0]`, Mercury at `[1 - mu, 0, 0]`, both fixed in the frame (p3).
The same convention as the project's (primary at -mu, secondary at 1 - mu, pulsating, true anomaly).

### 1.2 Equations of motion with sail (READ, p3, eqs. 2 to 4)

With dots denoting derivatives with respect to the independent variable f (the text says "独立变量f(真近点角)用来表征时间", TRANSLATED: "the independent variable f, the true anomaly, is used to characterise time"):

```
x'' =  2 y' + d(omega)/dx + a_SRP,x
y'' = -2 x' + d(omega)/dy + a_SRP,y        (eq. 2)
z'' =           d(omega)/dz + a_SRP,z

omega(x,y,z,f) = Omega~(x,y,z,f)/(1 + e cos f)
Omega~ = (1/2)(x^2 + y^2) + (1 - mu)/r1 + mu/r2 + (1/2) mu (1 - mu) - (1/2) e z^2 cos f      (eq. 3)
r1 = sqrt((x + mu)^2 + y^2 + z^2),   r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2)                  (eq. 4)
```

The note on the acceleration (READ, p3): 注意这里的加速度同样经过无量纲化处理 ("note that this acceleration is also non-dimensionalised"). The scaling
of `a_SRP` with the pulsating length and time units is NOT printed. COMPUTED (algebra): the sail-free part is the same function as in
`core/er3bp.py` apart from the additive constant `mu(1 - mu)/2`, which does not enter the gradient; the signs of the Coriolis terms
are the project's. Quote on the model: 式(3)说明椭圆限制性三体问题是一个非自治系统 ("eq. (3) shows that the ERTBP is a non-autonomous system", TRANSLATED).

Symmetry used (p3, eq. 5): `e cos(f +- pi) = -e cos f`. Chinese: 这一变换在某些情况下是非常方便的，例如，若航天器初始时刻对应的水星的真近点角为π，则可以将偏心率取相反数，同时将积分的初始时刻设为零，而积分得到的航天器轨迹与原先相同。 TRANSLATED: "This transformation is very convenient in some cases: for example, if the Mercury true anomaly at the spacecraft's initial time is pi, the eccentricity can be negated and the initial time of integration set to zero, and the integrated trajectory is the same as the original." (Valid for the sail-free equations by inspection; with a sail the same transformation is used by the authors; the dependence of `a_SRP` on f through d and the attitude is not discussed. INFERRED.)

### 1.3 Parameters (READ, p6, table 1)

| Parameter | Value (as printed) |
|---|---|
| Sun gravitational parameter (km^3 s^-2) | 1.327 1 x 10^11 |
| Mercury gravitational parameter (km^3 s^-2) | 22 032 |
| Mercury orbital eccentricity | 0.205 6 |
| Mercury orbital semi-major axis (km) | 5.971 x 10^7 |

Slip, stated factually: the semi-major axis in table 1 reads 5.971 x 10^7 km, while fig. 27's caption text on p21 gives
"(=5.791 x 10^7)", which is the usual value (57.91 million km). COMPUTED: with 5.791 x 10^7 km and the Sun plus Mercury
gravitational parameter, the Keplerian period is about 88.0 d, matching the 87.969 d the paper uses (p22); 5.971 x 10^7 km would
give about 92 d. So the table entry is a digit transposition; the computation evidently used 5.791 x 10^7 km (INFERRED).
COMPUTED: `mu = 22032/(1.3271e11 + 22032) = 1.660e-7`; the paper does not print mu.

### 1.4 Sail acceleration model (READ, p3 to p4, eqs. 6 to 12, figs. 2, 3)

Optical model of Rios-Reyes & Scheeres (ref. 24): flat sail, four forces (specular reflection, diffuse reflection, absorption,
thermal emission) with the Sun direction unit vector `s-hat` and the sail normal `n-hat`, cone angle `alpha = arccos(n-hat . s-hat)`,
alpha in [0, pi/2], and an azimuth angle beta in [-pi, pi] (eqs. 9, 10, fig. 3). Sail frame o_s x_s y_s z_s: origin at the sail, x_s along the
Sun direction, y_s perpendicular to x_s in the Sun-Mercury-sail plane (the "gravitational plane") pointing away from Mercury, z_s
completing the right-handed triad. In that frame `s = [1 0 0]^T`, `n = [cos(alpha), sin(alpha) cos(beta), sin(alpha) sin(beta)]^T`.

Acceleration (eq. 7): `a_SRP = (S0/(c d^2)) sigma cos(alpha) [ (C1 cos(alpha) + C2) n-hat + C3 s-hat ]`, `sigma = A/m` (the
area-to-mass ratio, "面质比", the sail's lightness measure; used in m^2/kg), with `S0 = 1366 kg/s^3` (solar constant), c the speed of
light printed "2.98 x 10^5 km/s", `d = d_s/AU` the sail-Sun distance in AU, AU = 1.496 x 10^8 km. Coefficients (eq. 8):

```
C1 = 2 rho s
C2 = B_f rho (1 - s) + (B_f eps_f - B_b eps_b)/(eps_f + eps_b) (1 - rho)
C3 = 1 - rho s
```

Optical values used (p4, from NEA Scout, ref. 5): `C1 = 1.711, C2 = 0.002, C3 = 0.145`. COMPUTED: C1 = 2 rho s = 1.711 gives
rho s = 0.8555 and so C3 = 1 - 0.8555 = 0.1445, consistent with the printed 0.145. No characteristic acceleration is printed; the
sail is described by sigma only.

Slips, stated respectfully (INFERRED): the speed of light is printed as 2.98 x 10^5 km/s (the usual value is 2.998 x 10^5 km/s, a
difference of 0.6 percent in the SRP acceleration if it was used as printed); eq. 6 writes the specular reflectance as delta in the
diffuse term and s in the specular term and in eq. 8, with s never defined in the text; and in the specular force of eq. 6 the page
image shows a squared cosine factor which would not reproduce eq. 7 (where C1 multiplies one cos(alpha) inside the bracket and one
outside), so one of the two is probably a typesetting slip. The text-layer extraction does not settle the exponent; this is not
confirmed.

Sail position is assumed to be anywhere in the frame (not on Mercury's side only); the sail force is the full 3-D vector in the
OXYZ frame via `s = C^s_r s^s`, `n = C^s_r n^s` (eqs. 11, 12, with `C` the rotation from the sail frame to OXYZ). No shadowing by
Mercury is in the equations; eclipse is checked afterwards (section 4.5 below).

### 1.5 The three steering laws (READ, p4, section 3, fig. 4)

Chinese: 目前已有文献中的太阳帆指向律主要有2种：一种是正对太阳，另一种是与其它引力体的运动相关联，如地月连线指向律。本文所采用的太阳帆指向律可以通过前文所定义的2个太阳帆姿态角统一表示：α≡α*，β≡β*。 TRANSLATED: "In the existing literature there are mainly two sail steering laws: one faces the Sun, the other is tied to the motion of another gravitating body, such as the Earth-Moon line law. The steering laws used here can be expressed uniformly with the two sail attitude angles defined above: alpha = alpha*, beta = beta*."

1. Sun-pointing law: alpha* = 0 (beta irrelevant).
2. Cone-varying law: alpha* in [0, pi/2], beta* in {0, pi}; the normal lies in the gravitational plane at angle alpha from the Sun line.
3. Azimuth-varying law: alpha* in [0, pi/2], beta* in [-pi, pi]; the normal is neither toward the Sun nor parallel to the
   gravitational plane; alpha fixed non-zero and beta not 0 or 180 deg.

The sail loading sigma and, for laws 2 and 3, alpha are held fixed in each run; the attitude angles are constants in the frame of
section 1.4 (so they rotate with the Sun-sail line) (INFERRED from the text; no time-varying schedule is printed).

## 2. What is computed

### 2.1 Orbit family and periodicity (READ, p4 to p7, section 4, eq. 13, figs. 5 to 7)

Resonant halo orbits of the Sun-Mercury problem about L1 and L2. The ERTBP is periodic with period 2 pi in f, so periodic orbits exist
only when the orbit period is an integer ratio of Mercury's period. Chinese (p2): 在椭圆限制性三体问题中，周期轨道不再连续分布，而是以离散的共振轨道形式存在，其轨道周期与第二主天体的轨道周期呈整数比。 TRANSLATED: "In the ERTBP periodic orbits are no longer continuously distributed but exist as discrete resonant orbits, whose orbital period is an integer ratio to the orbital period of the second primary."

Resonance ratio notation p:q with p/q the halo orbit frequency relative to Mercury's: the CR3BP halo period is `2 pi q/p`
(in time units where Mercury's period is 2 pi); in the ERTBP the orbit closes after `2 pi q` in f. Example (p7): "5:2 orbit, CRTBP period 4pi/5,
ERTBP period 4pi". COMPUTED check against fig. 5's dashed lines: Mercury's period 87.969 d divided by 4, 3.5, 3, 2.5 gives 21.99, 25.13,
29.32, 35.19 d, matching the dashed resonance lines at about 22, 25, 29 and 35 d for 4:1, 7:2, 3:1 and 5:2. Chinese on geometry
(p5): 对于共振比为p∶q的共振晕轨道，p/q的值越大，轨道离拉格朗日点越远，离第二引力体越近。 TRANSLATED: "For a resonant halo orbit of ratio p:q, the larger p/q, the farther the orbit from the Lagrange point and the nearer the second primary."

Resonances computed (tables 2 to 4): 3:1 (also written 6:2), 4:1, 5:2, 6:2, 7:2, each at L1 and L2, with two initial-epoch groups.

### 2.2 Initial-condition classification: perihelion and aphelion groups (READ, p5 to p6, section 4.2)

Chinese: 首先本文所采用的共振晕轨道的初值均位于脉动旋转坐标系XZ平面内，其对应的初始时刻仅有2个，一个是水星位于近日点（Perihelion）时，另一个是水星位于远日点（Aphelion）时。 TRANSLATED: "First, the initial values of the resonant halo orbits used here all lie in the XZ plane of the pulsating rotating frame, and there are only two corresponding initial times: one when Mercury is at perihelion, the other when it is at aphelion."

Labels: `IC[p:q, L1|L2, P|A]` (perihelion group P, aphelion group A). Only three states are non-zero for the symmetric
orbits (x0, z0, ydot0 at y0 = xdot0 = zdot0 = 0, INFERRED from fig. 9 and p17 text). By eq. 5 the aphelion group is integrated from f = 0 with
NEGATIVE eccentricity; Chinese (p6): 因此下文若出现负偏心率，仅代表其轨道初值属于远日点组（即初始时刻水星位于远日点），而实际天体轨道的偏心率始终为正。 TRANSLATED: "Therefore any negative eccentricity below only means the orbit's initial value belongs to the aphelion group (that is, Mercury is at aphelion at the initial time); the eccentricity of the real body's orbit is always positive." Mercury's semi-major axis and eccentricity set the dimensional scale.

### 2.3 What "asymmetric" means and how such orbits are found (READ, p4, p5, p10, p16)

Under the sun-pointing and cone-varying laws the sail halo orbits are symmetric about the rotating-frame XZ plane. Under the
azimuth-varying law they are not symmetric about any plane or line yet remain periodic. Chinese (p5): 在变方位角指向律下，共振晕轨道不关于任何平面对称，但仍然保持其周期性。 TRANSLATED: "Under the azimuth-varying law the resonant halo orbit is not symmetric about any plane, but it still keeps its periodicity."

The mechanism (p16, section 5.4): starting from the beta = pi cone-varying orbit and moving beta, "the sail-plane normal always lies on a circular cone about the Sun-light direction". Because the symmetry is lost the initial state is no longer in the XZ plane and needs all six components (p16): 此时轨道初值不能仅用3个非零状态量表示（x、z和ydot），而是需要6个状态量才能完整表示. TRANSLATED: "the initial value can then no longer be expressed with only three non-zero states (x, z and ydot), and all six states are needed". Symmetry property found (p16 to p18): orbits with the same sigma and alpha and opposite beta are mirror images about the XZ plane; the states x, z, ydot are even functions of beta and y, xdot, zdot are odd functions (fig. 23, p18); orbits with beta in the negative half drift to y < 0, positive beta to y > 0.

Orbits are found by multiple shooting (p5, section 4.1, eq. 13): N equal segments, N + 1 nodes of 6 states, segment time `dt = T/N`; bidirectional
integration to the segment midpoint (continuity `X_i^{+dt/2} - X_{i+1}^{-dt/2} = 0`, plus `X_{N+1} - X_1 = 0`); MATLAB ode113 and fmincon
interior-point with a constant (zero) objective, `StepTolerance 1e-8`, `ConstraintTolerance 1e-10`; fsolve with trust-region-dogleg
(`FunctionTolerance 1e-11`, `StepTolerance 1e-8`) is an alternative and takes within 3 to 20 percent of the time. No value of N is
printed for the CR3BP/ERTBP stage; N = 2p for the ephemeris stage (p21, 2p + 1 nodes).

### 2.4 Continuation procedure (READ, p6, section 4.3)

Four steps: (1) CR3BP resonant halo orbit with mu fixed, then continuation in Mercury's eccentricity to the real value (the sail-free
ERTBP orbit); (2) from that, sun-pointing law, continuation in sigma until the solver fails or the orbit hits Mercury's surface
(radius 2440 km) or sigma reaches 8 m^2/kg (set above LightSail 2's 6.4 m^2/kg, ref. 7); (3) from step 2, cone-varying law at fixed
sigma, beta = 0 or pi, continue alpha until failure or alpha = pi/2; (4) from step 3's orbit, azimuth-varying law at fixed sigma and
alpha, continue beta until failure or until beta covers [-pi, pi]. Pseudo-arclength continuation in e (from Spreen's thesis, ref. 25,
p7) is used for the eccentricity step. The authors note that the full three-parameter (sigma, alpha, beta) database was not built:
"当面质比较大时，变锥角与变方位角指向律下的共振晕轨道容易不收敛" (TRANSLATED: "when the area-to-mass ratio is large, resonant halo orbits under the cone-varying and azimuth-varying laws easily fail to converge").

### 2.5 Sail-free ERTBP results: bifurcations in eccentricity (READ, p6 to p8, figs. 6 to 10)

CRTBP resonant halos about L1 and L2 (fig. 6: L1 and L2 3:1, 4:1, 5:2, 7:2, the Sun-Mercury halos are nearly symmetric about Mercury because mu is tiny).
Continuation in Mercury's eccentricity (fig. 7) gives, from each CR3BP orbit, one perihelion-group and one aphelion-group orbit.
Period-doubling bifurcations along the L1-3:1 and L2-3:1 families (fig. 8: eigenvalues of the monodromy matrix in the complex plane,
coloured by eccentricity). Chinese (p7): 所谓倍周期分叉，是指若某一周期轨道的单值矩阵的3对特征根中有一对在复平面实轴上的-1点处碰撞，则存在另一个双倍周期的轨道族与原先的轨道族在分叉点处重合。 TRANSLATED: "A period-doubling bifurcation means that if, of the three pairs of eigenvalues of a periodic orbit's monodromy matrix, one pair collides at the point -1 on the real axis of the complex plane, then there is another family of doubled period that meets the original family at the bifurcation point."

Chinese (p7): L1-3∶1共振晕轨道在偏心率分别为-0.0179，0.0089，0.04处发生了倍周期分叉；L2-3∶1共振晕轨道在偏心率分别为-0.0189，0.0094，0.044处发生了倍周期分叉，可见L1与L2侧的3∶1共振晕轨道的分叉点是很接近的。 TRANSLATED: "The L1-3:1 resonant halo orbits undergo period-doubling bifurcations at eccentricities -0.0179, 0.0089 and 0.04; the L2-3:1 resonant halo orbits at -0.0189, 0.0094 and 0.044; so the bifurcation points of the L1-side and L2-side 3:1 orbits are very close." (Negative values are the aphelion group, section 2.2.)

The doubled-period branches are 6:2 orbits, named `IC[L1/L2, 6:2, A1/A2/A3]` for the bifurcations at the first, second and third eccentricity respectively (fig. 8 panels a,b / c,d / e,f). Continuing each in eccentricity to the real value (+-0.2056): chinese (p7 to p8) 结果表明图8中仅图8（c）~图8（f）4处倍周期分叉对应的6∶2共振晕轨道可以延拓到-0.2056，而图8（a）和图8（b）2处在延拓到0.037与0.038时开始不收敛。 TRANSLATED: "The results show that only the 6:2 resonant halo orbits corresponding to the four period-doubling bifurcations of figs. 8(c) to 8(f) can be continued to -0.2056, while the two of figs. 8(a) and 8(b) begin not to converge when continued to 0.037 and 0.038." So in the Sun-Mercury system there are four 6:2 halo orbits (A2 and A3 at L1 and L2), all of the aphelion group (p8, fig. 10 shows the L2 A3 orbit at e = -0.2056, with initial-condition curves in fig. 9).

### 2.6 Sail results, sun-pointing law (READ, p8 to p13, table 2, figs. 11 to 15)

Continuation in sigma stops by non-convergence or collision with Mercury; no orbit reached the 8 m^2/kg ceiling. Largest sigma:
IC[L2, 3:1, A], 6.65 m^2/kg. Chinese (p10): 从延拓结果来看，不论哪种类型的共振晕轨道都没有达到预设的面质比的上限8 m²/kg，其中最大面质比为IC[L2,3∶1,A]的6.65 m²/kg。 TRANSLATED: "From the continuation results, no type of resonant halo orbit reached the preset upper limit of 8 m^2/kg; the maximum is 6.65 m^2/kg for IC[L2,3:1,A]." L2-side orbits continue to larger sigma than L1 ones (they are easier to converge), because the radiation always arrives from the Sun direction, which breaks the L1/L2 near-symmetry. Orbits move toward the Sun (toward -X) as sigma grows; the minimum distance to Mercury's surface decreases with sigma (fig. 15, with sharp kinks when the time of minimum distance jumps between 88 d, 59 d and 118 d for L1-5:2 P). The count of distance extrema per orbit changes with sigma (L1-5:2 P: 5 minima and 5 maxima at small sigma, 7 and 7 above 2 m^2/kg, p10 to p11).

### 2.7 Sail results, cone-varying law (READ, p13 to p15, table 3, figs. 16 to 20)

Run at a fixed sigma from the sun-pointing orbit; alpha continued from 0 toward 90 deg. At alpha = 90 deg the radiation force vanishes, so the
orbit returns to the sail-free ERTBP orbit, and the end points of the beta = 0 and beta = pi branches coincide although the paths
differ (p13, p14). The two beta cases drive the orbit in opposite directions: Chinese (p13): 随着锥角的逐渐增大，轨道整体朝着远离太阳的方向移动，这一点与正对太阳指向律下的延拓结果刚好相反。 TRANSLATED: "As the cone angle grows, the whole orbit moves in the direction away from the Sun; this is exactly the opposite of the continuation results under the sun-pointing law."

### 2.8 Sail results, azimuth-varying law (READ, p16 to p19, table 4, figs. 21 to 23)

Run at sigma = 1 m^2/kg with alpha held. The set of azimuths over which the orbit exists depends on alpha (table 4). Chinese (p17):
当锥角较小时，方位角可以在[-π，π]范围内存在；而随着锥角逐渐变大，方位角仅在0与π附近存在（注意-π与π实际对应同一个方向）；随着锥角继续增大到接近π/2，方位角的存在范围又恢复到[-π，π]。 TRANSLATED: "When the cone angle is small the azimuth can exist over [-pi, pi]; as the cone angle grows the azimuth exists only near 0 and pi (note that -pi and pi are the same direction); as the cone angle grows further toward pi/2 the range of azimuth recovers to [-pi, pi]." Converged orbits are easier to find when alpha tends to 0 or pi/2 and when beta tends to 0 or pi (sail normal nearly parallel to the gravitational plane). Examples: fig. 21 (L2-5:2, alpha 30 deg, beta -20 deg) and fig. 22 (beta = +20 deg), mirror images about XZ.

### 2.9 Stability (READ, p20 to p24, section 6, figs. 27 to 35)

No linear stability analysis (no Floquet multipliers for the sail orbits) is printed; "stability" is tested by integration in a
high-fidelity ephemeris model without station-keeping. Model (p21 to p22): JPL DE435 (ref. 26 cites the SPICE toolkit paper; the
text says DE435); gravitation of the Sun, Mercury, Venus and Earth; SRP as in section 1.4; Mercury-centred ecliptic J2000 frame,
converted to a Mercury-Sun rotating non-pulsating frame (X to the Sun, Z along Mercury's angular momentum). The orbit is re-corrected
by the multiple-shooting method of section 2.3 with 2p + 1 nodes (N = 2p), equally spaced in time. Start epochs: perihelion group
2026-02-19 10:41:41 UTC; aphelion group 2026-04-04 10:19:38 UTC (p22 to p23). Integration times: L1-3:1 for 5 Mercury periods, L2-3:1 for 6, L2-4:1 and L2-6:2 for 4, L2-5:2 and L1-5:2 (law 2 and 3 examples) for 4 (Mercury period 87.969 d, p22). Fig. 27 shows Mercury's osculating elements 2025-01-01 to 2030-01-01 against the idealised ERTBP values: semi-major axis offset about -500 to -1200 km, eccentricity offset about 2.3 to 5.1 x 10^-5, inclination about 7.0032 to 7.0035 deg, ascending node about 48.2935 to 48.2995 deg (read off the plots).

Results (p22 to p24): sail orbits of the sun-pointing law (L1-3:1 at 4 m^2/kg, L2-3:1 at 4, L2-4:1 at 3, L2-6:2 A3 at 1 m^2/kg),
cone-varying (L2-5:2 at alpha 40 deg, beta 0 and 180 deg, sigma 1) and azimuth-varying (L1-5:2, alpha 63 deg, beta +-80 deg, sigma 1, 4 periods).
Chinese (p22): 可以明显看出L1-3∶1、L2-3∶1及L2-4∶1晕轨道在最后一个周期出现发散趋势，而L2-6∶2晕轨道由于仅积分了2个周期...并没有出现明显的发散趋势。 TRANSLATED: "It is clearly seen that the L1-3:1, L2-3:1 and L2-4:1 halo orbits show a diverging trend in the last period, while the L2-6:2 halo orbit, having been integrated for only 2 periods (here a period means the orbit's own period, not a Mercury period), shows no obvious divergence." The ephemeris orbits are only approximately symmetric about the Mercury-centred XZ plane. Chinese (p24): 在变方位角指向律下，所有的晕轨道均可以在星历模型下保持其轨道形状至少4个水星轨道周期，至多可以保持6个。 TRANSLATED: "Under the azimuth-varying law all the halo orbits keep their orbit shape in the ephemeris model for at least 4 Mercury orbital periods and at most 6." COMPUTED: 4 x 87.969 = 351.9 d and 6 x 87.969 = 527.8 d, matching the "352 to 528 d" of the abstract.

Eclipse (p18 to p20, figs. 24 to 26): a Mercury-shadow test by projecting the orbit on the YZ plane of the non-pulsating frame; the sun-pointing
L1-5:2 P orbit at 2.5 m^2/kg (alpha 0) has an eclipse; the cone-varying orbit (alpha 10 deg, beta 180 deg, same sigma) avoids it. The
paper presents this as a use of the cone-varying law, not as a computation of eclipse durations.

## 3. Printed tables, transcribed digit by digit

All three are READ from the page image and checked against the text layer. The paper's own typos in table titles ("restlts" in table 2,
"vanying" in table 3, "azimath" in table 4) are in the original.

**Table 1 (p6)**, parameters of the ERTBP: see section 1.3.

**Table 2 (p10)**, continuation results under the sun-pointing steering law. Columns: orbit type; maximum area-to-mass ratio
(m^2/kg); reason the continuation ended. 不收敛 = did not converge; 碰撞水星 = collided with Mercury (TRANSLATED).

| Orbit type | Max sigma (m^2/kg) | End reason |
|---|---|---|
| IC[L1,3:1,P] | 2.160 | did not converge |
| IC[L1,3:1,A] | 2.450 | did not converge |
| IC[L1,4:1,P] | 0.245 | did not converge |
| IC[L1,4:1,A] | 0.240 | did not converge |
| IC[L1,5:2,P] | 2.570 | collided with Mercury |
| IC[L1,5:2,A] | 2.000 | did not converge |
| IC[L1,6:2,A2] | 0.460 | did not converge |
| IC[L1,6:2,A3] | 0.485 | did not converge |
| IC[L1,7:2,P] | 0.910 | did not converge |
| IC[L1,7:2,A] | 0.900 | did not converge |
| IC[L2,3:1,P] | 4.020 | collided with Mercury |
| IC[L2,3:1,A] | 6.650 | collided with Mercury |
| IC[L2,4:1,P] | 3.640 | collided with Mercury |
| IC[L2,4:1,A] | 4.200 | did not converge |
| IC[L2,5:2,P] | 4.100 | collided with Mercury |
| IC[L2,5:2,A] | 5.080 | did not converge |
| IC[L2,6:2,A2] | 1.260 | did not converge |
| IC[L2,6:2,A3] | 1.980 | did not converge |
| IC[L2,7:2,P] | 0.435 | collided with Mercury |
| IC[L2,7:2,A] | 0.230 | collided with Mercury |

**Table 3 (p14)**, continuation results under the cone-varying steering law. Columns: orbit type; azimuth beta; sigma (m^2/kg);
maximum cone angle reached (deg).

| Orbit type | beta | sigma (m^2/kg) | Max alpha (deg) |
|---|---|---|---|
| IC[L1,3:1,A] | 0 | 2 | 14.1 |
| IC[L1,3:1,A] | pi | 2 | 90 |
| IC[L1,5:2,P] | 0 | 2 | 90 |
| IC[L1,5:2,P] | pi | 2 | 90 |
| IC[L2,4:1,A] | 0 | 3 | 26.5 |
| IC[L2,4:1,A] | pi | 3 | 90 |
| IC[L2,5:2,A] | 0 | 3.5 | 90 |
| IC[L2,5:2,A] | pi | 3.5 | 90 |
| IC[L2,6:2,A3] | 0 | 1 | 90 |
| IC[L2,6:2,A3] | pi | 1 | 90 |
| IC[L2,7:2,P] | 0 | 0.3 | 90 |
| IC[L2,7:2,P] | pi | 0.3 | 90 |

**Table 4 (p19)**, continuation results under the azimuth-varying steering law, sigma = 1 m^2/kg. Columns: orbit type; cone angle
(deg); range of azimuth for which the orbit exists (deg). Braces in the original join the three sub-intervals.

| Orbit type | Cone angle (deg) | Azimuth range (deg) |
|---|---|---|
| IC[L1,5:2,P] | 1 | [-180, 180] |
| | 8 | [-180, 180] |
| | 9 | {[-180, -103], [-69.8, 69.8], [103, 180]} |
| | 50 | {[-180, -154.4], [-27.4, 27.4], [154.4, 180]} |
| | 62 | {[-180, -104.8], [-74.4, 74.4], [104.8, 180]} |
| | 63 | [-180, 180] |
| | 80 | [-180, 180] |
| IC[L2,5:2,P] | 1 | [-180, 180] |
| | 8 | [-180, 180] |
| | 9 | {[-180, -103], [-75.9, 75.9], [103, 180]} |
| | 50 | {[-180, -146], [-32.4, 32.4], [146, 180]} |
| | 63 | {[-180, -103.6], [-75.6, 75.6], [103.6, 180]} |
| | 64 | [-180, 180] |
| | 80 | [-180, 180] |

Numbers in figures (READ off captions and colour bars, not tables): fig. 8 colour-bar eccentricity ranges about -0.0172 to -0.0188
(panel a), -0.0180 to -0.0200 (b), 0.00875 to 0.00910 (c; panel d's bar is printed "8 900 to 9 900" without the decimal prefix,
INFERRED to be the same quantity scaled by 10^-6, about 0.0089 to 0.0099), 0.038 to 0.042 (e) and 0.0425 to 0.045 (f). Fig. 9 annotates
the bifurcation values -0.0178 (L1; the text says -0.0179) and -0.0189 (L2), 0.0089 and 0.0094, 0.04 and 0.044.

## 4. Slips and inconsistencies noted (stated factually, INFERRED)

- Table 1 semi-major axis 5.971 x 10^7 km against 5.791 x 10^7 km on p21 (section 1.3).
- Orbit-group labels differ between table 3, the figure captions and the text for the same figures: table 3 lists IC[L2,4:1,A] and
  the p14 and fig. 16 and 17 captions say aphelion group, but fig. 20(c) is labelled IC[L2,4:1,P] and the p14 text says
  "IC[L2,4:1,P] ... not converging at sigma 3 for beta = 0"; table 3 and 4 list IC[L1,5:2,P] while the Chinese captions of
  figs. 18 and 19 and 34 say aphelion group; the Chinese caption of fig. 32 says perihelion group (近日点组) while the text on p23 says
  aphelion (远日点组) and table 3 lists IC[L2,5:2,A]; fig. 20(a) is titled "IC[L1, 3:3, A]". One label in each pair is probably a slip; the paper does not say which.
- The figure-9 annotation -0.0178 (L1-A1) against -0.0179 in the text.
- The speed of light and the specular-force exponent (section 1.4).
- The English abstract's "pitch angle" for the azimuth angle.

## 5. Positive controls for the project

The paper prints NO periodic orbit state, period or closure value for any orbit, with or without a sail. No orbit can be rebuilt from
it for closure testing. Candidate checks (all derived from printed numbers, none are initial-state goldens):

1. Sail-free Sun-Mercury elliptic problem (`mu` COMPUTED 1.660e-7 from table 1; `e = 0.2056`): the L1-3:1 and L2-3:1 resonant halo
   families, seeded by the project's own CR3BP halo corrector at period `2 pi/3`, should show a monodromy eigenvalue at -1 (period
   doubling) at e = -0.0179, 0.0089, 0.04 (L1) and e = -0.0189, 0.0094, 0.044 (L2), where negative e means start at f = pi from the
   +e orbit by eq. 5. The doubled-period 6:2 branches born at the second and third should continue to e = +-0.2056 and the first
   should fail near e = 0.037 (L1) and 0.038 (L2). This tests the elliptic model's stability and continuation behaviour with printed
   numbers but not closure; the precision is that of the printed digits (three to four decimals in e).
2. For a closure test with printed 15-digit states in the same system, use Peng, Bai & Xu (2017) table 2 (the ME-halo states,
   digested in `2026-06-25-digest-peng-2017-sun-mercury-ERTBP.md`); this paper is a sail extension of that line of work and does not
   reprint it.

## 6. What this means for the project

Scope. Sail-augmented, non-ballistic periodic orbits about the Sun-Mercury L1 and L2 points in the elliptic problem. They are not
cyclers, not quasi-cyclers and not resonant orbits with flybys (the orbits are near the libration point with no planetary
flyby; resonance means a rational ratio of halo period to Mercury's period). Outside the catalogue's four classes. The sail terms
also put them outside any ballistic criterion the project applies.

Overlap with open tasks. Only the sail-free sub-result (section 2.5) touches the project: it is a second published check on the
`core/er3bp.py` model and on elliptic-problem continuation (#893), at a mass ratio some five orders of magnitude from Earth-Moon,
and the pulsating-frame equation (eq. 3) agrees with the project's potential. It offers nothing for the #884 sun-forced or
bicircular work.

Reusable methods (INFERRED): the (P, A) initial-epoch classification with the eq. 5 trick of using negative eccentricity for the
aphelion group; the shooting formulation with bidirectional midpoint matching; the staged continuation (e, then sigma, then alpha,
then beta) with stop criteria; the observation that a lost symmetry requires a full six-state initial condition; and the proxy test of "natural
stability" by integrating a periodic-model orbit in an ephemeris model without control for 4 to 6 primary periods.

## 7. Peng & Xu (2015) and the corrupt held copy

The paper cites Peng & Xu (2015, CMDA 123:279) as ref. 13 together with Peng, Bai & Xu (2017) as ref. 14 (p2: Peng et al. studied the
influence of the primary's gravitational constant and the secondary's eccentricity on resonant halo orbits in the ERTBP). It reprints no
table, initial condition or result from either. It fills none of the Peng & Xu (2015) gap.

## 8. References cited (as printed; translations of Chinese titles are mine)

Elliptic problem and resonant halo orbits:
- [12] Hou X Y, Liu L. On motions around the collinear libration points in the elliptic restricted three-body problem [J]. Monthly
  Notices of the Royal Astronomical Society, 2011, 415(4): 3552-3560.
- [13] Peng H, Xu S J. Stability of two groups of multi-revolution elliptic halo orbits in the elliptic restricted three-body problem
  [J]. Celestial Mechanics and Dynamical Astronomy, 2015, 123(3): 279-303.
- [14] Peng H, Bai X L, Xu S J. Continuation of periodic orbits in the Sun-Mercury elliptic restricted three-body problem [J].
  Communications in Nonlinear Science and Numerical Simulation, 2017, 47: 1-15.
- [25] Spreen E M Z. Trajectory design and targeting for applications to the exploration program in cislunar space [D]. West
  Lafayette, Indiana: Purdue University, 2021 (cited for the pseudo-arclength continuation in eccentricity).
- [9] (circular problem, for context) Saleem V, Ram K. Families of periodic orbits about Lagrangian points L1, L2 and L3 with continuation method [J]. Planetary and
  Space Science, 2022, 217: 105491.

Solar-sail orbits in restricted problems (including the four-body and ephemeris checks):
- [15] McInnes C R. Solar sailing: Orbital mechanics and mission applications [J]. Advances in Space Research, 2003, 31(8): 1971-1980.
- [16] Bookless J, McInnes C. Control of Lagrange point orbits using solar sail propulsion [J]. Acta Astronautica, 2008, 62(2/3): 159-176.
- [17] Heiligers J, Hiddink S, Noomen R, et al. Solar sail Lyapunov and Halo orbits in the Earth-Moon three-body problem [J]. Acta
  Astronautica, 2015, 116(11/12): 25-35.
- [18] Heiligers J, Macdonald M, Parker J S. Extension of Earth-Moon libration point orbits with solar sail propulsion [J].
  Astrophysics and Space Science, 2016, 361(7): 241.
- [19] An Shiyu (安诗宇), Liu Ming (刘明), Li Huayi (李化义), et al. Design and applications of novel periodic orbits with solar sail
  in Earth-Moon system [J]. Acta Aeronautica et Astronautica Sinica, 2025, 46(4): 230828 (in Chinese).
- [20] Chujo T, Takao Y. Synodic resonant halo orbits of solar sails in restricted four-body problem [J]. Journal of Spacecraft and
  Rockets, 2022, 59(6): 2129-2147.
- [21] Gong S P, Li J F. Solar sail periodic orbits in the elliptic restricted three-body problem [J]. Celestial Mechanics and
  Dynamical Astronomy, 2015, 121(2): 121-137.
- [22] Huang J, Biggs J D, Cui N G. Families of halo orbits in the elliptic restricted three-body problem for a solar sail with
  reflectivity control devices [J]. Advances in Space Research, 2020, 65(3): 1070-1082.
- [23] Gamez Losada F, Heiligers J. New solar-sail orbits for observation of the Earth and Moon [J]. Journal of Guidance, Control,
  and Dynamics, 2021, 44(12): 2155-2171.
- [24] Rios-Reyes L, Scheeres D J. Generalized model for solar sails [J]. Journal of Spacecraft and Rockets, 2005, 42(1): 182-184.
- [26] Acton C H. Ancillary data services of NASA's navigation and ancillary information facility [J]. Planetary and Space Science,
  1996, 44(1): 65-70.

Context references: [1] Zhao P Y, Wu C C, Li Y M, Chinese Journal of Aeronautics 2023, 36(5): 125-144 (solar sailing review);
[2] Garwin 1958; [3] Tsuda et al. 2013 (IKAROS); [4] Johnson et al. 2011 (NanoSail-D); [5] Lockett et al. 2020 (NEA Scout);
[6] Ridenoure et al. 2016 (LightSail); [7] Spencer et al. 2021 (LightSail 2); [8] Takao et al. 2023 (OKEANOS); [10] Zhang Chen and Zhang Hao
2023; [11] Zhang Chen 2024 (distant retrograde orbit transfers, in Chinese).
The paper cites no separate reference for model-transition (CR3BP to ephemeris) studies beyond [20] and [23], which
verify their orbits in ephemeris models.

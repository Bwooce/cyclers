# Digest: Canalias 2007 PhD Thesis
**Source:** Canalias 2007 PhD thesis, "Homoclinic and heteroclinic connections between planar Lyapunov orbits..."
**Subject:** SE <-> EM manifold connections, heteroclinic bifurcations

## 1. Jacobi Convention
The Jacobi constant is defined in Equation 4.4:
`C(x, y, x_dot, y_dot) = -(x_dot^2 + y_dot^2) + 2 \Omega(x, y)`
where `\Omega(x, y) = 1/2(x^2 + y^2) + (1-\mu)/r_1 + \mu/r_2`
This expands to:
`C = (x^2 + y^2) + 2(1-\mu)/r_1 + 2\mu/r_2 - (x_dot^2 + y_dot^2)`
This exactly matches our `cr3bp.py` convention. There is no `\mu(1-\mu)` offset.

## 2. SE L1/L2 Lyapunov Families & Bifurcation Values
Canalias doesn't tabulate individual starting `[x, y, xdot, ydot]` ICs for the Lyapunov orbits in the text, but provides the limits of the usable Jacobi constants and corresponding amplitudes (Table 4.1):
- SE L1: `C \in [3.0007222915, 3.00090098]`, `Ax \in [1081.9, 359019.1]` km
- SE L2: `C \in [3.00072105, 3.0008969275]`, `Ax \in [548.3, 366122.8]` km
- EM L1: `C \in [3.149305, 3.20034403]`, `Ax \in [100.6, 12661.6]` km
- EM L2: `C \in [3.14445, 3.184163]`, `Ax \in [50.3, 17258.1]` km

Per-family bifurcation Jacobi values:
- **SE L1 Homoclinic (Ho11):** `C = 3.00088389`
- **SE L1->L2 Heteroclinic Bifurcation:** `C = 3.000863625` (Confirmed)

## 3. Section and Coupling Convention
The Sun-Earth and Earth-Moon problems are coupled by deploying the EM system inside the SE system. 
- `\beta` is the relative position angle between the Sun-Earth x-axis and the Earth-Moon x-axis. 
- Poincare section `S` is defined by the angle `\phi_{SE}` (e.g., 90 degrees) with respect to the SE x-axis.
- Intersection is determined by observing where the integrated SE and EM manifolds meet on `S`.

## 4. Cross-System Patch Coordinates
Table 5.3 specifies maneuvers (in m/s) at the intersecting points between the unstable manifold of the EM L2 Lyapunov orbit (`Ax = 17123.2` km) and the stable manifold of the SE L2 Lyapunov orbit (`Ax = 330102.8` km). Exact coordinate values are deferred to the attached DVD (`inters.dat`, `sistemn.dat`), but representative delta-v costs and phases are given in the tables.

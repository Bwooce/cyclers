# Full-body registry expansion (#260)

Date: 2026-06-14

Goal: make every meaningful solar-system body available for solution generation,
with SOURCED data only. Two registries grew:

- `src/cyclerfinder/core/satellites.py` — planet-centric moon systems
  (`PRIMARIES` system GMs + `SATELLITES`).
- `src/cyclerfinder/core/constants.py` — heliocentric bodies (`PLANETS`), plus
  the inclined-element map in `src/cyclerfinder/core/ephemeris.py`.

The architecture is registry-driven: every consumer (Tisserand/V∞ screen,
ephemeris, flyby floor, discovery enumeration) keys off these tables, so adding a
sourced entry makes a body resolvable everywhere. No catalogue rows were written.

## Caveat (completeness, not advocacy)

Most additions are LOW-MASS → weak gravity assists → poor flyby/cycler nodes.
Phobos/Deimos and the small irregulars (Amalthea, Hyperion, Nix, Hydra) have
GM < 1 km³/s². The dwarf planets/large asteroids are tiny relative to the
planets. Triton is large but RETROGRADE and inclined. They are added so the
Tisserand/V∞ screen has the full body set to self-prune from. The genuinely
interesting additions are the five regular Uranian moons
(Miranda/Ariel/Umbriel/Titania/Oberon) and the Pluto–Charon pair (Charon is
~11% of the system GM → a real binary).

---

## PART A — satellite systems (satellites.py)

### Primary system GMs added (`PRIMARIES`), JPL DE440 planetary constants
(ssd.jpl.nasa.gov astro_par / gm_de440.tpc; Pluto via Brozovic et al. 2015),
accessed 2026-06-14:

| Primary | system GM (km³/s²) |
|---|---|
| Mars    | 4.282837521e4 |
| Uranus  | 5.7945564e6 |
| Neptune | 6.836527100580e6 |
| Pluto   | 9.755e2 (Pluto+Charon system) |

(Earth/Jupiter/Saturn were already present.)

### Moons added (`SATELLITES`)

All from JPL SSD satellite physical parameters (`sats/phys_par/`, ref code in
parentheses) for GM + mean radius, and JPL SSD mean orbital elements
(`sats/elem/`) for semi-major axis. Accessed 2026-06-14. `mean_motion_deg_day`
is DERIVED (Kepler III) from `sma_km` + the primary GM, never hand-copied.
`safe_alt_km` is an engineering default (10 km for tiny irregulars, 100 km
otherwise; convention, not sourced physics).

| Moon | Primary | GM (km³/s²) | mean R (km) | a (km) | ref |
|---|---|---|---|---|---|
| Phobos   | Mars    | 0.0007087  | 11.08  | 9 375     | MAR097 |
| Deimos   | Mars    | 0.0000962  | 6.2    | 23 457    | MAR097 |
| Miranda  | Uranus  | 4.3        | 235.8  | 129 846   | URA111 |
| Ariel    | Uranus  | 83.5       | 578.9  | 190 929   | URA111 |
| Umbriel  | Uranus  | 85.1       | 584.7  | 265 986   | URA111 |
| Titania  | Uranus  | 226.9      | 788.9  | 436 298   | URA111 |
| Oberon   | Uranus  | 205.3      | 761.4  | 583 511   | URA111 |
| Triton   | Neptune | 1428.49546 | 1352.6 | 354 800   | NEP097 |
| Proteus  | Neptune | 2.58342    | 208.0  | 117 600   | NEP097 |
| Charon   | Pluto   | 106.1      | 606.0  | 19 600    | PLU060 |
| Nix      | Pluto   | 0.0015     | 18.0   | 49 300    | PLU060 |
| Hydra    | Pluto   | 0.0020     | 18.5   | 65 200    | PLU060 |
| Amalthea | Jupiter | 0.16456    | 83.5   | 181 400   | JUP365 |
| Iapetus  | Saturn  | 120.51511  | 734.3  | 3 561 700 | SAT441 |
| Hyperion | Saturn  | 0.37049    | 135.0  | 1 481 500 | SAT441 |

### Omitted for lack of a sourced value

- **Nereid (Neptune)** — JPL SSD lists its GM as 0.0 (mass not determined).
  Per the sourcing discipline we omit it rather than guess from a size estimate.

---

## PART B — heliocentric planetoids (constants.py PLANETS)

Two-letter codes (single letters are taken by the major planets).
`mean_motion_deg_day` is DERIVED from `sma_au` + MU_SUN (Kepler III).
`safe_alt_km` = 100 km engineering default. `ecc` lives in `PLANETS`;
inclination/Ω live in the inclined-element map (`ephemeris.py`), per the
existing convention (the live coplanar `PLANETS` table is never mutated). The
non-coplanar inclinations (Pluto 17°, Eris 44°, Makemake 29°, Haumea 28°)
matter for the 3-D Tisserand screen.

### GM (`mu_km3_s2`)

- **Pluto, Ceres, Vesta, Pallas** — JPL DE440 `gm_de440.tpc` (Park et al.,
  AJ 2021), direct GM:
  - Pluto (system) 975.5; Ceres 62.6288886444; Vesta 17.2882328792;
    Pallas 13.6658781460.
- **Eris, Makemake, Haumea** — JPL/SBDB publishes no GM, only a mass from
  satellite dynamics, so GM = G·M with G = 6.67430e-20 km³ kg⁻¹ s⁻²
  (CODATA 2018). MASS-derived → less precise (flagged inline):
  - Eris: M = (1.638±0.014)e22 kg (Brown & Schaller 2007, Science 316:1585)
    → GM ≈ 1093.25.
  - Makemake: M = (2.69±0.20)e21 kg (Bamberger et al. 2025, from moon
    S/2015 (136472) 1) → GM ≈ 179.54.
  - Haumea: M = (3.952±0.011)e21 kg (Proudfoot et al. 2024) → GM ≈ 263.77.

### Orbital elements (`sma_au` / `ecc` / inc / Ω)

- **Pluto** — Standish & Williams "Keplerian Elements for Approximate Positions
  of the Major Planets" (JPL SSD), **Table 2a (3000 BC–3000 AD)** — the only
  Standish table that retains Pluto (Table 1 dropped it post-IAU-2006).
  Verified against two independent copies (the JPL PDF and the soniakeys/aprx
  implementation): a0=39.48686035, e0=0.24885238, I0=17.14104260,
  L0=238.96535011, ϖ0=224.09702598, Ω0=110.30167986. Pluto's phasing angles
  (ϖ/L0) are therefore sourced; the SBDB-element bodies leave them at 0.0.
- **Ceres, Eris, Makemake, Haumea, Vesta, Pallas** — JPL Small-Body Database
  osculating heliocentric elements (`full-prec=1`), all at the common epoch
  **JD 2461200.5 (TDB)**, accessed 2026-06-14. Osculating (not a mean-element
  table), so ϖ/L0 are left at the 0.0 default — only sma/ecc/inc feed the screen.

| Body | code | a (AU) | e | i (deg) | Ω (deg) |
|---|---|---|---|---|---|
| Pluto    | Pl | 39.48686035        | 0.24885238          | 17.14104260       | 110.30167986 |
| Ceres    | Ce | 2.765552595034094  | 0.07969229514816586 | 10.58802780183462 | 80.24862682043221 |
| Eris     | Er | 67.93394687853566  | 0.4382385347971672  | 43.9258279471791  | 36.00477044417249 |
| Makemake | Mk | 45.57093317300052  | 0.1588889953992523  | 29.02785603743067 | 79.2948338209406 |
| Haumea   | Ha | 43.06029023650952  | 0.1944430148898797  | 28.20847393040364 | 121.7860561329425 |
| Vesta    | Ve | 2.361365965127599  | 0.09020374382834395 | 7.143925545058711 | 103.701293265032 |
| Pallas   | Pa | 2.769559010737709  | 0.2307000995648547  | 34.93279321851542 | 172.8866193357694 |

### Radii (`radius_eq_km`, flyby periapsis floor only)

Pluto 1188.3 (New Horizons / IAU 2015 WGCCRE), Ceres 469.7 (Dawn / WGCCRE),
Eris 1163 (Sicardy et al. 2011), Makemake 715 (Ortiz et al. 2012), Haumea ~780
(volume-equivalent, Ortiz et al. 2017; markedly triaxial), Vesta 262.7
(Russell et al. 2012), Pallas 256.5 (Marsset et al. 2020).

### DE440 ephemeris note

DE440 carries the 8 planets + the Pluto barycenter, but NOT the main-belt / TNO
small bodies. So Ceres/Eris/Makemake/Haumea/Vesta/Pallas are excluded from the
astropy real-ephemeris backend map (they remain fully usable via the
circular/inclined backends the Tisserand screen uses); Pluto stays resolvable.

---

## Tests

- `tests/core/test_satellites_registry.py` — extended: every new moon present +
  keyed by full name; mean motion derived; **each system enumerates ONLY its own
  moons** (Mars→{Phobos,Deimos}, Uranus→its five, etc.) and the slices are
  mutually disjoint.
- `tests/core/test_primaries.py` — extended: Mars/Uranus/Neptune/Pluto system
  GMs present + in expected bands.
- `tests/core/test_constants_planetoids.py` — new: sourced spot-checks
  (EXPECTED = published value), Kepler-III mean motion, V∞ ceiling finite,
  code-collision guard, major-planet byte-identical guard, enumerate/Tisserand
  behavioural smoke.
- `tests/core/test_constants_phasing_angles.py`,
  `tests/verify/test_ephemeris_crosscheck.py` — updated for the now-superset
  registry (was "exactly 8 planets").

Full suite green (`uv run pytest`, exit 0; only the pre-existing S1L1/Aldrin
model-boundary XFAILs remain). `ruff check`/`ruff format --check` clean;
`mypy src tests` clean.

## Decisions under ambiguity

- **Pluto appears in BOTH registries.** In `PRIMARIES`/`SATELLITES` it is the
  central body of the Pluto–Charon system (system GM 975.5). In `PLANETS` it is
  the heliocentric body (also GM 975.5 — DE440 publishes only the system value
  and Charon is ~11% of it; the system GM is what perturbs a heliocentric
  flyby). Documented inline in both files.
- **Mass→GM conversion** only for Eris/Makemake/Haumea (no published GM); used
  CODATA-2018 G and cited the mass source. All other GMs are direct DE440.
- **Osculating vs mean elements.** Standish (mean) is unavailable for the
  dwarf planets except Pluto, so SBDB osculating at a single common epoch was
  used and the time-phasing angles (which a single osculating epoch cannot
  supply consistently) were left at 0.0 rather than fabricated.

# Planet-centric moon-tour cycler support — design draft (task #76 / #117)

**Status:** design draft (no code, no task plan). DOCS-ONLY.
**Supersedes/extends:** `docs/superpowers/specs/2026-06-02-planet-centric-moon-tours-scoping.md`
(80-line scoping note — this draft is its brainstorming follow-on and resolves
its open questions).
**On-ramp note mined for inputs:** `docs/notes/2026-06-05-endgame-tisserand-mining.md`
(Campagnola & Russell "Endgame Problem" Parts 1 & 2, transcribed anchors).
**Scope fences (locked, see §9 non-goals):** M-3D (heliocentric inclination,
executing now); M-ED (heliocentric real-ephemeris multi-arc, landed); the Forge
(verdict gauntlet). This draft does **not** touch any of those; it adds a *new
central body* axis orthogonal to all three.

---

## 0. Headline finding (one tier is mostly built, the other is genuinely new)

The brief's premise — "make the 6 non-heliocentric rows computable and open
planet-centric search" — splits cleanly into **two tiers with very different
build cost**, because the existing machinery is *centre-parametric in some
layers and Sun-hardwired in others*. Read-only audit on 2026-06-06:

**Already centre-agnostic (works around any primary if fed the right μ/states):**

- **Lambert is centre-agnostic by signature.** `lambert(...)` takes
  `mu: float = MU_SUN_KM3_S2` as an explicit argument
  (`src/cyclerfinder/core/lambert.py:335`); the whole Stumpff/`_t_of_z` core
  threads `mu` through (`core/lambert.py:122-167`). Passing `mu = μ_Jupiter`
  makes it a Jovicentric solver with no code change — only the *default* is
  heliocentric.
- **The Tisserand layer is already written in canonical (dimensionless) form**
  and is the literal subject of the on-ramp note. `vinf_to_tisserand` /
  `tisserand_to_vinf` use `T = 3 − V∞²·a_p/μ` (`core/tisserand.py:88-110`);
  `linkable` intersects two bodies' constant-T contours (`core/tisserand.py:306-411`);
  `linkable_3d` adds the `cos(i)` term (`core/tisserand.py:471-680`). The note
  confirms these are exactly the Endgame Part-2 identity `T = 3 − v∞² ≈ J`
  (note lines 230-232) and that `linkable`'s contour-intersection *is* the
  multibody patch-point predicate (note lines 463-465). **The math generalises;
  only the body-SMA lookup `_a_p_km` is Sun-relative** (`core/tisserand.py:78-80`,
  reads `PLANETS[body].sma_au * AU_KM`).
- **The registry is the single source of truth and is designed to be extended.**
  `SUPPORTED_BODIES = tuple(PLANETS.keys())` with the explicit contract "adding a
  (sourced) `PlanetData` entry to `PLANETS` is all that is needed to make a new
  body resolvable everywhere" (`core/constants.py:229-236`). Task #75 made compute
  key off this table rather than hardcoding V/E/M.

**Sun-hardwired (blocks tier 1 until patched):**

- **The ballistic corrector assumes heliocentric throughout.** `ballistic_correct`
  / `_vinf_nodes` call `ephem.state(body, t)` and `lambert(r1, r2, tof)` **with
  no `mu` argument** (`search/correct.py:115`), so Lambert silently defaults to
  `MU_SUN`. `_max_bend_deg` reads `PLANETS[body].mu_km3_s2` (`search/correct.py:182-193`)
  — fine *if* the registry holds the moon. The frame assumption is implicit-Sun
  via the ephemeris, not via Lambert. **Verdict: the corrector is *patched-conic
  around any primary* in principle, but needs (a) the moon in the registry and
  (b) `mu=` plumbed into the `lambert()` call.** Small, surgical.
- **`PlanetData.sma_au` is intrinsically heliocentric** (field doc
  `core/constants.py:138-139`); a moon has no `sma_au`. So moons cannot be
  jammed into `PLANETS` as-is — they need a sibling record with `sma_km`
  about-the-primary (see §3).
- **The ephemeris is hard-Sun.** Both backends are heliocentric by construction:
  `_circular_inplane_state` uses `sma_au * AU_KM` (`core/ephemeris.py:117`); the
  astropy backend literally subtracts the Sun's barycentric posvel
  (`core/ephemeris.py:263-264`). A moon-tour propagation needs **planet-centred
  states** (moon position about its primary), which neither backend provides.
- **`vinf_kms_at_encounters` / the canonical V∞ signature is heliocentric V∞.**
  For Jovicentric tier-1 this is fine (V∞ about Jupiter is a real conic
  quantity). For the CR3BP rows it is *model-mismatched* — the rows already say
  so and store `vinf: null` with Jacobi as the real invariant
  (`data/catalogue.yaml:7462-7472`).

**Conclusion:** the **patched-conic moon-system tier (Jovian/Saturnian) is
reachable with a moon registry + a centred ephemeris + plumbing `mu=` into the
corrector** — mostly-existing machinery. The **Earth-Moon CR3BP tier needs
genuinely new dynamics** (CR3BP propagator, Jacobi constant, no conic V∞).
These are **two distinct tiers** and the recommendation (§8) is to ship tier 1
first.

---

## 1. The two-tier split (deliverables + which rows each unlocks)

### Tier 1 — Patched-conic moon systems (circular-coplanar / linked-conic)

Same dynamical model as the heliocentric catalogue (Kepler conics + impulsive
gravity-assist V∞ rematch), but with the **central body = a planet** and the
**flyby bodies = its moons**. The on-ramp's VILM / Tisserand machinery is the
linked-conic feasibility layer (note §"Paper 1", lines 44-63).

**Rows unlocked:**
- `hernandez-2017-jovian-ieg-triple-family` — `model_assumption: circular-coplanar`,
  `primary: Jupiter`, bodies Io/Europa/Ganymede, exploits the 1:2:4 Laplace
  resonance (`data/catalogue.yaml:7930-8089`). **Cleanest tier-1 target** — the
  row is *already* tagged circular-coplanar and explicitly patched-conic.
- `russell-strange-2009-jovian-multimoon-family` — `circular-coplanar`,
  `primary: Jupiter`, Galilean-multi (`data/catalogue.yaml:8094-8226`).
- `russell-strange-2009-saturnian-multimoon-family` — tagged `cr3bp` but the
  SCHEMA-MISMATCH note (`data/catalogue.yaml:8321-8336`) says it is **mixed**:
  Titan cyclers are patched-conic-modellable; Enceladus/Mimas/Tethys tend toward
  low-energy / manifold (CR3BP-natural) regimes. **Tier 1 unlocks the Titan
  sub-family; the midsize-moon members defer to tier 2.** This row should
  arguably be re-tagged `circular-coplanar` for its Titan members on ingest, per
  the existing note.

**Deliverables (tier 1):**
1. `SATELLITES` registry (§3): per-moon μ, radius, SMA-about-primary, mean
   motion-about-primary, sourced from JPL SSD.
2. A planet-centred circular ephemeris: a moon backend analogous to
   `_circular_inplane_state` but scaled in km about the primary (§3, §4).
3. `mu=` plumbed into `search/correct.py`'s `lambert()` call + a `primary`/`mu`
   parameter on `ballistic_correct` so the corrector closes a moon-tour chain
   around Jupiter/Saturn (§4).
4. `_a_p_km` in `core/tisserand.py` taught to resolve a moon's about-primary SMA
   (§4) so `linkable`/`linkable_3d` prune Jovicentric/Saturnian moon pairs.
5. **(Optional, the on-ramp's headline)** a VILM module implementing Endgame
   Part-1 Eq. (13) quadrature ΔV-floor + Eq. (9) V̄∞-efficiency root, validated
   against the mined anchors (§6). This is the *feasibility/search* layer; the
   corrector is the *verification* layer. Either can ship without the other.

### Tier 2 — CR3BP (Earth-Moon and low-energy moon regimes)

Genuinely new dynamics: a rotating-frame circular restricted three-body
propagator, the **Jacobi constant** as the conserved invariant (there is **no
conic V∞**), and validation against published periodic-orbit families. Mirrors
the JPL three-body periodic-orbit catalog (spec §16.7.4-§16.7.6,
`docs/spec.md:937-968`).

**Rows unlocked:**
- `arenstorf-em-figure8-1963` — CR3BP figure-8, `mass_ratio: 0.01215`
  (`data/catalogue.yaml:7490-7504`).
- `genova-aldrin-2015-em-3petal-cycler` — CR3BP 3:1-lunar-resonance, **already
  carries `trajectory.center: Earth` + 2 segments** (`data/catalogue.yaml:7659-7704`).
- `wittal-2022-em-cycler-family` — CR3BP family seed, ~90° inclination to lunar
  plane (`data/catalogue.yaml:7795-7853`).
- The midsize-moon (Enceladus/Mimas/Tethys) members of the Russell-Strange
  Saturnian row (low-energy regime, per its note).

**Deliverables (tier 2):**
1. A planar (then spatial) CR3BP propagator + `mass_ratio`-parametrised EOM
   (note lines 276-279: Endgame Part-2 Eq. (2)/(4) give the EOM + Jacobi
   constant form).
2. A Jacobi-constant computation and a **CR3BP differential corrector** (the
   §16.7.5 dispatch row: "CR3BP differential corrector (M8+)",
   `docs/spec.md:953`).
3. The `orbit_elements.cr3bp{}` identity populated from a backfill source (the
   rows' `data_gaps` all point at the JPL periodic-orbit catalog, e.g.
   `data/catalogue.yaml:7505-7520`).
4. A non-keplerian validation path that asserts `(jacobi_constant, period_nd,
   stability_index)` rather than a V∞ multiset (§6).

---

## 2. What "V∞ signature" and "period" mean per tier

The brief's hard question: extend the canonical signature/matching without
breaking heliocentric rows.

| Concept | Heliocentric (today) | Tier 1 (patched-conic moon) | Tier 2 (CR3BP) |
|---|---|---|---|
| Energy invariant | V∞ about the Sun | **V∞ about the primary** (Jovicentric V∞ is a real conic quantity; note line 230 `T=3−v∞²`) | **Jacobi constant** (no conic V∞ exists) |
| `vinf_kms_at_encounters` | populated | populated (planet-relative) | `null` + model-mismatch note (already done, `catalogue.yaml:7462-7472`) |
| Period | 2-body synodic `{pair,k,years}` | moon-pair synodic / Laplace-resonance multiple (e.g. Io 1.77 d, `catalogue.yaml:7949-7953`) | `period_nd` (dimensionless rotating-frame period) + `lunit/tunit` de-norm |
| Identity block | `orbit_elements{a,e}` | `orbit_elements{periapse_km, apoapse_km, reference_frame: planetcentric-inertial, center}` (spec §16.7.2, `docs/spec.md:896-907`) | `orbit_elements.cr3bp{jacobi_constant, period_nd, stability_index, state_nd, mass_ratio, lunit_km, tunit_s}` (spec §16.7.4) |

**How the canonical signature extends without breaking heliocentric rows:**

The schema already solved this. Per spec §16.2 and §16.7, **`model_assumption`
is a *pool pre-filter*, not a signature input** (`docs/spec.md:551`): "M7 *may*
additionally pre-filter the matcher pool by `model_assumption` … but this is a
*pool filter*, not a signature input." So:

- The **canonical signature stays exactly as-is for heliocentric rows** (V∞
  multiset + `(a,e)` + period). No field changes; no golden moves.
- A new finder hit is **first bucketed by `(model_assumption, primary)`**, then
  matched only against the same-bucket pool. A Jovicentric V∞ never compares to
  a heliocentric V∞ because they are different buckets — and a CR3BP Jacobi
  never compares to any V∞ at all (spec §16.7.5 dispatch table,
  `docs/spec.md:949-953`).
- Tier 1's V∞ is the *same kind of scalar* as heliocentric V∞ (conic excess
  speed), just about a different μ — the matcher code is unchanged; only the
  pool partition key gains `primary`.
- Tier 2's signature is `(jacobi_constant, period_nd, stability_index)` — a
  **separate dispatch branch** that never enters the V∞ comparator. This is the
  "incommensurable signature" the spec already mandates (`docs/spec.md:276`).

**Net: heliocentric rows are byte-identical; tier 1 reuses the V∞ comparator
under a new bucket key; tier 2 adds a parallel non-keplerian comparator.**

---

## 3. Moon/body registry design

A moon cannot live in `PLANETS` because `PlanetData.sma_au` is intrinsically
heliocentric (`core/constants.py:138-139`). Two viable shapes:

**Shape A (recommended): a sibling `SATELLITES` registry keyed by primary.**

```text
@dataclass(frozen=True)
class SatelliteData:
    name: str                # "Europa"
    code: str                # "Eu" (2-letter; 1-letter space is V/E/M-exhausted)
    primary: str             # "Jupiter" — links to a primary's μ
    mu_km3_s2: float         # moon GM
    radius_eq_km: float
    sma_km: float            # SMA ABOUT THE PRIMARY (not AU, not Sun-relative)
    mean_motion_deg_day: float   # about the primary, derived from sma_km + primary μ
    safe_alt_km: float

SATELLITES: dict[str, SatelliteData]   # keyed by moon code
PRIMARIES:  dict[str, float]           # primary code -> μ_primary (Jupiter, Saturn, Earth)
```

This mirrors the existing `PLANETS`/`SUPPORTED_BODIES` contract exactly: add a
sourced record, and the registry-keyed machinery (Tisserand `_a_p_km`, the
corrector's `_max_bend_deg`, a centred ephemeris) resolves the moon. The
about-primary `mean_motion_deg_day` is **derived at import** from `sma_km` +
`μ_primary` via Kepler III, exactly like `_mean_motion_deg_day`
(`core/constants.py:149-159`), so the table stays internally consistent rather
than hand-copied.

**Shape B: overload `PlanetData` with an optional `primary`/`sma_km`.** Rejected
— it pollutes the heliocentric record with always-null fields and risks a
heliocentric caller reading `sma_au` on a moon. Shape A keeps the moon model in
its own type.

**Sourcing discipline (non-negotiable, per global golden-tests rule + the v4.2
backfill checklist):** every numeric field traces to a *published* source, never
to our own compute.

- **Primary sourcing source:** JPL SSD (`http://ssd.jpl.nasa.gov/`) — the
  Endgame Part-1 Table 3 footnote (note line 351) cites exactly this for its
  moon physical data. **Use the original JPL tables, not the paper's
  transcription, for the *registry values*** — the paper's table is a *golden
  anchor* (independent check), so it must stay independent of the registry it
  validates (golden-discipline circularity rule).
- **Cross-check anchor:** the mined Endgame Part-1 Table 3 (note lines 337-352)
  gives μ̃, ã, Ṽ for Io/Europa/Ganymede/Callisto + Enceladus/Tethys/Dione/Rhea/
  Titan. Registry value sourced from JPL → must reproduce the paper's ã/Ṽ to a
  documented tolerance. This is a registry-construction golden, separate from §6.
- **Earth-Moon mass ratio** for tier 2 is already in the catalogue rows
  (`mass_ratio: 0.01215`, sourced Genova & Aldrin 2015, `catalogue.yaml:7494-7497`).
- **`source_ephemeris` / `center` / `tof_days_bounds` backfill:** per the v4.2
  backfill checklist, every new registry entry and every newly-computable row
  must record its center and source ephemeris.

---

## 4. Minimal code-touch map (tier 1)

Surgical changes that make tier-1 rows computable, with the centre-agnostic
layers reused as-is:

1. **`core/constants.py`** — add `SatelliteData` + `SATELLITES` + `PRIMARIES`
   (§3). Heliocentric `PLANETS` untouched.
2. **`core/ephemeris.py`** — add a `_CentredCircularBackend` (moon on its
   mean-motion circle *about the primary*, km-scaled) selected by a new
   `center=` / `model=` option. The astropy backend can later return real moon
   states (astropy resolves e.g. "europa" via the body-name map,
   `core/ephemeris.py:60-63`); the circular backend is enough for the
   circular-coplanar rows.
3. **`core/tisserand.py:78-80`** — `_a_p_km` resolves a moon code to
   `SATELLITES[m].sma_km` (about-primary). `vinf_to_tisserand`/`tisserand_to_vinf`
   already take μ as `MU_SUN_KM3_S2` *implicitly via the planet path* — they need
   a `mu=` (the primary's μ) so the canonical-units identity is about the right
   centre. **This is the only Tisserand change**; `linkable`/`linkable_3d`
   contour logic is centre-blind once `_a_p_km` and μ are right.
4. **`search/correct.py:115`** — pass `mu=μ_primary` into the `lambert()` call;
   add a `primary: str` / `mu_central: float` parameter to `ballistic_correct`
   and `_vinf_nodes`. `_max_bend_deg` already keys off the per-body registry
   (`correct.py:182-193`) — point it at `SATELLITES` for moon codes. **No
   residual-math change**: V∞-continuity + closure is centre-agnostic.

**Verification of the brief's claim "is correct.py centre-agnostic?":**
*Almost.* Its residual/bend logic is centre-agnostic; its two Sun-couplings are
(a) the implicit `MU_SUN` Lambert default at `correct.py:115`, and (b) the
heliocentric ephemeris it is handed. Fix both (pass `mu`, hand it a centred
ephemeris) and it closes a patched-conic chain around any primary. The brief's
"it should be patched-conic around any primary" is achievable with these two
edits.

---

## 5. T-P graph + VILM as the tier-1 feasibility/search layer

The on-ramp note is explicit that the two Endgame graphs *are* our Tisserand
machinery generalised (note lines 444-465):

- **Leveraging Graph (Part 1):** V∞ level sets + ΔV as a third axis. The
  ΔV-min **quadrature** (Eq. 13, note lines 153-162) is a closed-form ΔV floor
  for any VILM tour leg — usable as an **admissible lower bound for A*-style
  pruning** inside a Forge search and as a fast feasibility screen (note line
  458). The n:m_K± taxonomy (note lines 98-130) classifies each VILM leg.
- **T-P graph (Part 2):** Tisserand level sets that **extend past the
  linked-conic feasible region** (`T > 3` ⇒ ballistic intermoon transfers that
  linked-conics calls impossible, note lines 209-215). The **patch point** =
  intersection of two moons' T_M level sets (Eq. 13, note lines 242-243) — the
  multibody generalisation of our `linkable` contour intersection (note lines
  463-465).

**How it slots in:** tier 1's search layer = `linkable`/`linkable_3d` (already
the M4 pruning gate) prunes infeasible moon pairs at a given V∞; the **VILM
quadrature gives the ΔV cost** on the surviving pairs; the corrector
(`ballistic_correct`, §4) closes the chosen chain. The T-P extension (`T > 3`
ballistic region) is a **tier-1.5 enhancement** — it needs CR3BP reachability,
so it bridges to tier 2; ship the linked-conic VILM first.

---

## 6. Validation gates — which anchors gate what

Golden discipline: the EXPECTED side traces to a *published* source, never to
our own compute (global rule + project memory). The mined anchors (note
lines 329-442) are the gates. **Two suspect cells are explicitly excluded as
goldens** (note lines 412, 491-493): Endgame Part-2 Table 1 Titan J_L4@100km
"766.5/776.4" (MAX<MIN inversion). Do not use either as an EXPECTED value.

**Tier 1 gates:**

| Gate | Anchor | What it validates |
|---|---|---|
| Registry construction | Part-1 Table 3 ã/Ṽ (note 337-352) | `SATELLITES` SMA/velocity reproduce JPL-sourced values |
| VILM efficiency root | Part-1 Table 3 **V̄∞ E/I** (note 337-354) | a future `min_vinf_for_vilm(moon)` (Eq. 9 root) |
| ΔV-min quadrature | Part-1 Table 1 **ΔV_min** (no-GA, note 358-381) + Table 2 (with-GA, note 383-397) | Eq. (13) quadrature ΔV-floor per moon pair |
| Worked scalar | Part-1 Europa 3-VILM: 154 m/s / 46 d (note 436-438) | end-to-end VILM endgame ΔV |
| Ballistic paradox | Part-2 Table 1 (note 399-411, **minus** the flagged cell) | insertion-ΔV near-invariance under arrival angle |

**Tier 2 gates:**

| Gate | Anchor | What it validates |
|---|---|---|
| Jacobi at patch | Part-2 Ganymede/Europa **J=3.0052 / 3.0023** (note 424-428) | CR3BP Jacobi-constant computation |
| Patch-point r_a/r_p | Part-2 `(694641, 1021834) km` (note 426) | T_M1=T_M2 intersection (Eq. 13) |
| Periodic-orbit identity | JPL three-body catalog `(jacobi, period_nd, stability)` | `orbit_elements.cr3bp{}` backfill for Arenstorf/Genova/Wittal |

**Model-fidelity caveat to carry forward (note lines 286-293):** the *same*
endgame can read 5-10% different ΔV between linked-conic and CR3BP, and
ballistic (Δv=0) intermoon transfers exist in CR3BP that linked-conics declares
impossible. Any tier-1 ΔV is therefore an upper-bound-ish linked-conic
reference — the validator must use a tolerance band, not equality, and must not
"reject" a CR3BP row for disagreeing with a linked-conic anchor by ≤10%.

---

## 7. Fidelity-ladder + gauntlet integration

**What the rungs mean planet-centric.** The Axis-B ladder is
`circular-coplanar → analytic-ephemeris → real-de440`
(`src/cyclerfinder/data/provenance.py:95`; `verify/fidelity.py:15-25`). For a
**moon tour the ladder is the same three rungs but about the primary**:

- `circular-coplanar` → the `_CentredCircularBackend` (§4) + the resonance
  construction. This is the rung tier-1 rows live at (their `model_assumption`).
- `analytic-ephemeris` → still the documented **extension point** (no in-house
  backend today, `fidelity.py:22,167-169`); for moons this would be an analytic
  satellite theory.
- `real-de440` → astropy moon states (astropy resolves Galilean moons; the
  body-name map auto-derives from the registry, `core/ephemeris.py:60-63`). This
  rung's persistence check (does V∞ stay stable circular→real?) is the same
  Axis-B logic, now about Jupiter.

**Tier 2 does not fit the V∞ ladder** — its tracked quantity is the Jacobi
constant. The ladder generalises as: rung 1 = analytic CR3BP (the published
`mass_ratio`), rung 2/3 = ephemeris-perturbed CR3BP (real lunar ephemeris +
solar third-body, which is what Genova & Aldrin's "modest phasing maneuvers"
compensate for, `catalogue.yaml:7754-7761`). The **persistence axis tracks
Jacobi-constant drift** rather than V∞ drift — a parallel `PersistenceReport`
on a different scalar, same `_moves_toward_band` logic
(`verify/fidelity.py:393-411`).

**Gauntlet.** The Axes A-D combiner (`verify/gauntlet.py`) folds axis reports
into a verdict; it "never invents a value." Planet-centric changes are confined
to **what each axis computes**, not the combiner:
- Axis A (code-path agreement) — VILM quadrature vs corrector ΔV on the same
  moon pair is a *second code path* for tier 1 (exactly the kind of crosscheck
  Axis A wants).
- Axis B (fidelity persistence) — V∞ (tier 1) or Jacobi (tier 2) across rungs.
- Axis C (provenance) — the mined anchors as corroboration.
- Axis D (falsification) — a deliberately-bogus moon-pair / wrong-μ guard.

The §16.7.5 validation dispatch (`docs/spec.md:944-953`) already routes by
`cycler_class`: `single-ellipse`→constructor, `multi-arc`→multi-leg solver,
`non-keplerian`→CR3BP corrector. Tier 1 rows are `non-keplerian` today but are
*really* patched-conic multi-arc moon tours — **on ingest, re-tag the Jovian
patched-conic rows so the dispatch sends them to the multi-leg solver, not the
CR3BP corrector** (the Saturnian note already anticipates this split).

---

## 8. Architecture options + recommendation

**Option 1 — Tier 1 only, registry + plumbing (recommended first slice).**
Add `SATELLITES`, a centred circular ephemeris, plumb `mu=` into the corrector
and Tisserand. Unlocks the 2 Jovian rows + the Titan Saturnian sub-family. No
new dynamics. ~4 surgical edits (§4) on centre-agnostic substrate.
*Cost: low. Risk: low. Coverage: 3 of 6 rows (partial on the Saturnian).*

**Option 2 — Tier 1 + VILM feasibility module.** Option 1 plus the Endgame
Part-1 quadrature/efficiency module (Eq. 9, 13), validated against the mined
anchors, wired as a Forge search/feasibility layer with the ΔV-floor as an
admissible bound. *Cost: medium. Risk: low (closed-form, well-anchored).
Coverage: same rows as Opt 1, but adds the search capability the brief calls
"open planet-centric cycler search as a scope."*

**Option 3 — Full two-tier (Option 2 + CR3BP).** Adds the CR3BP propagator,
Jacobi constant, periodic-orbit corrector, non-keplerian validation. Unlocks all
6 rows. *Cost: high (new dynamics, multiple-shooting, basin-aware seeding — note
lines 304-307 warn gradient optimizers get trapped). Risk: high.*

**Recommendation: Option 2, sequenced as Option 1 → VILM add-on → (later)
Option 3.** Ship the patched-conic tier first; it is mostly registry + plumbing
on already-centre-agnostic Lambert/Tisserand machinery, unlocks the most rows
per unit effort, and the on-ramp note's anchors are linked-conic (so they
validate tier 1 directly). CR3BP is a separate, larger milestone (M8-class) that
should not gate the patched-conic win.

### Phased sketch

- **Phase 1 (tier 1 substrate):** `SATELLITES`/`PRIMARIES` registry + sourcing +
  registry-construction golden (Part-1 Table 3). Centred circular ephemeris.
- **Phase 2 (tier 1 compute):** plumb `mu=`/`primary` into `correct.py` +
  `_a_p_km`/μ into `tisserand.py`; close `hernandez-2017-jovian-ieg` end-to-end;
  re-tag Jovian rows `circular-coplanar`; gauntlet dispatch fix.
- **Phase 3 (VILM feasibility):** Eq. (9) V̄∞ root + Eq. (13) quadrature module,
  validated against Part-1 Tables 1-3 (excluding flagged cells); n:m_K± taxonomy.
- **Phase 4 (tier 2, separate milestone):** planar CR3BP propagator + Jacobi +
  periodic-orbit corrector; backfill `cr3bp{}` from JPL catalog; non-keplerian
  validation; unlock Arenstorf/Genova/Wittal + Saturnian midsize moons.

---

## 9. Non-goals (locked)

- **No heliocentric changes.** M-3D, M-ED, the Forge verdict logic, and every
  heliocentric golden stay byte-identical. The canonical signature is unchanged
  (§2).
- **No interplanetary→moon-capture chains** (a true "cruise to Jupiter *then*
  tour the moons" multi-patched-conic stack). Single central body per cycler, as
  the 2026-06-02 scoping note recommended. Deferred.
- **No asteroids / small-body tours** (the scoping note's separate track).
- **No CR3BP in tier 1.** The `T > 3` ballistic-transfer region (Part-2) and the
  Earth-Moon rows are tier 2 / Phase 4.
- **No invariant-manifold machinery.** The on-ramp explicitly avoids manifolds
  (note lines 319-325); if tier 2 ships, it uses patch-point + multiple-shooting,
  not manifold intersection.
- **No silent fixing of the two flagged suspect anchor cells** (note lines
  491-493) — they are excluded as goldens, preserved verbatim.

---

## 10. Open questions for the user

1. **Tier-1-first, yes?** Confirm shipping the patched-conic Jovian/Saturnian
   tier (Option 2) before any CR3BP work — i.e. accept that 3 of 6 rows
   (Arenstorf, Genova, Wittal) stay citation-only until a later milestone.
2. **VILM module in-scope or defer?** Is the Endgame Part-1 quadrature/efficiency
   module (the "open planet-centric *search*" half) wanted in this milestone, or
   is making the rows *computable* (corrector only) enough for now?
3. **Re-tag the mis-flagged rows?** OK to change `russell-strange-2009-saturnian`
   and the Jovian rows' `model_assumption`/dispatch to `circular-coplanar` for
   their patched-conic members (the Saturnian note already recommends this), or
   keep them as-is and handle the split purely in the matcher?
4. **Moon body-code convention?** 1-letter codes are V/E/M-exhausted. Adopt
   2-letter moon codes (Io, Eu, Ga, Ca, Ti, En…)? `data/README.md` reportedly
   reserved Saturnian codes already (`catalogue.yaml:8338-8340`) — confirm the
   scheme.
5. **CR3BP backfill source.** All tier-2 rows point `data_gaps` at the JPL
   three-body periodic-orbit catalog for `(jacobi, period_nd, stability)`. Is
   pulling those values (when tier 2 lands) acceptable as the golden source, or
   do you want the original Arenstorf/Genova/Wittal papers re-mined first?

---

## Approval (2026-06-06)

User-approved with all recommendations accepted: (Q1) Tier-1-first — ship the
patched-conic Jovian/Saturnian tier; Arenstorf/Genova/Wittal stay
citation-only until the CR3BP milestone; (Q2) the VILM quadrature/efficiency
module is IN scope (Option 2); (Q3) re-tag the mis-flagged Jovian/Saturnian
rows to their patched-conic reality (the Jovian rows' `non-keplerian` tag
would misroute the gauntlet dispatch); (Q4) adopt 2-letter moon body codes
(Io, Eu, Ga, Ca, Ti, En, …), reconciled against the codes already reserved in
the catalogue; (Q5) the JPL three-body periodic-orbit catalog is the accepted
golden source for Tier-2 `(jacobi, period_nd, stability)` backfill when that
tier lands. Next step: implementation plan (Tier 1 + VILM).

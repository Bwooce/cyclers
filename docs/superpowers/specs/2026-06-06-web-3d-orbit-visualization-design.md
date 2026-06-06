# Web orbit-visualization design draft (cyclers.space)

**Status:** DESIGN DRAFT — no code. Revisits the spec §15 out-of-scope line
("later add … a 3D orbit viewer").
**Date:** 2026-06-06
**Task:** #130
**Scope of this doc:** rendering-approach options + recommendation, data
pipeline (seeds-not-tracks), a phased scope ladder, site-principles compliance,
in-UI honesty markers, open questions. It is a planning artefact only; it
prescribes nothing that ships without a follow-up implementation plan.

---

## 0. Survey — what exists today (cite file:line)

The site is deliberately minimal and this design must not break that.

- **Stack:** Astro 6 static, zero runtime JS frameworks. `dependencies` are
  `astro`, `typescript`, `@astrojs/check`; `devDependencies` are `js-yaml` +
  types (`cyclers.space/package.json:18-27`). No React/Vue/Svelte, no Three.js,
  no WebGL anything today.
- **Design principles, verbatim:** *"legible, keyboard-navigable, no animation,
  no colour-only meaning"* (`cyclers.space/src/styles/global.css:1-2`), with a
  full dark-mode palette via `@media (prefers-color-scheme: dark)`
  (`global.css:21-35`). There is currently **no** `prefers-reduced-motion`
  block and **no** `@keyframes`/`transition`/`animation` rule anywhere in the
  stylesheet — the "no animation" rule is enforced by simply not having any.
- **How interactivity is done today:** one vanilla-TS `<script>` island inlined
  in `CatalogueTable.astro:120-223` — planet filter + column sort, no framework,
  progressive-enhancement style (it bails gracefully if the table is absent,
  `CatalogueTable.astro:138-144`). This is the established pattern any viz must
  follow: a small inlined/island script, no dependencies, JS-off safe.
- **Detail page:** `cyclers.space/src/pages/cycler/[id].astro` — statically
  generated per entry via `getStaticPaths()` (`[id].astro:35-40`). It already
  branches on `cycler_class`: single-ellipse renders heliocentric `(a,e,
  perihelion, aphelion, inclination)` (`[id].astro:202-214`); multi-arc renders
  invariants + per-segment legs, **no single (a,e)** (`[id].astro:216-238`);
  non-keplerian renders the CR3BP identity block only (`[id].astro:240-266`).
  This branch structure is exactly where a per-class viz hooks in.
- **Build pipeline:** `node scripts/sync-catalogue.mjs` runs at `predev` /
  `prebuild` (`package.json:9-11`), copying the upstream `data/catalogue.yaml`
  and `windows.json` into `src/data/`. Any build-time sampling step would attach
  here.
- **Data already on the page object** (`cyclers.space/src/lib/types.ts`):
  `orbit_elements` carries `a_au, e, inclination_deg` (`types.ts:68-72`) plus
  the v4 orientation children `raan_deg, arg_periapsis_deg, true_anomaly_deg,
  epoch_iso8601` (spec.md:485-488 — present in the schema, frequently `null`).
  `trajectory.segments[]` carry per-leg `a_au?, e?` (`types.ts:103`, may be
  null). `cr3bp` carries `jacobi_constant, period_nd, state_nd, lunit/tunit`
  (spec.md:503-512). `windows.json` carries **real** matched encounter dates,
  per-window `c3_km2_s2` and first-leg ToF (header at
  `cyclers.space/src/data/windows.json:1-12`).
- **Upstream viz already computes a trajectory plot** —
  `cyclers/src/cyclerfinder/viz/plots.py:162-214` (`trajectory()`): plots planet
  orbit *circles* (radius = encounter heliocentric distance,
  `plots.py:177-179`), then samples each leg by **two-body propagation** from the
  departure state — 80 samples via `propagate(r0, v0, dt)` (`plots.py:188-197`)
  — plus encounter `o` markers and a Sun `*`. Crucially it is a **2D XY
  (ecliptic-projected) plot** (`plots.py:208-209` label "x/y (km)"), it requires
  a *built* `OptimisationResult` (full state vectors), and matplotlib is a
  lazy-imported optional `[viz]` extra (`viz/__init__.py`). It is **not** wired
  to the catalogue YAML and the site never invokes it.

**Key takeaway from the survey:** the upstream plot proves the math (planet
circles + Kepler-sampled legs) but operates on *built state vectors*, which the
public catalogue rows generally do **not** carry. The site has only
literature-shape `(a,e,i,Ω,ω,ν?)` — often partial. So the viz must reconstruct
geometry from sparse Kepler elements, not replay stored tracks. This is the
"seeds-not-tracks" constraint and it shapes everything below.

---

## 1. Rendering-approach options

Three honest options, with the trade-off that matters most for *this* data:
**heliocentric cycler orbits are near-ecliptic.** Venus i = 3.39°, Mars i =
1.85°, Earth ≈ 0° (M-3D design,
`specs/2026-06-05-m-3d-inclination-lift-design.md:1,36-38`). All V/E/M cyclers
sit at i ≤ ~3.4°. **A top-down ecliptic view loses almost nothing visually** —
the out-of-plane excursion of a 1.5 AU leg at 2° is ~0.05 AU, a few pixels. 3D
becomes genuinely informative only for (a) the small z-structure the M-3D work
makes meaningful, and (b) future non-coplanar or planet-centric (Moon-tour)
families. For the bread-and-butter heliocentric catalogue, 3D is mostly
decoration.

### Option A — Three.js (or react-three-fiber)

- **What:** full WebGL scene, OrbitControls, real perspective camera.
- **Pros:** the obvious "3D orbit viewer" the spec gestures at; trivially
  handles future high-inclination / planet-centric cases; lots of reference
  material.
- **Cons:** bundle is **~150 KB+ min (≈600 KB unminified)** even tree-shaken,
  which is larger than the *entire current site's* JS. Pulls the site's first
  real heavyweight dependency. Canvas is opaque to screen readers and keyboard
  users unless you build a parallel a11y layer. Tension with "no animation" —
  the whole point of a WebGL canvas is the render loop. Dark-mode means managing
  a second set of materials/clear-colors. **Lazy-load discipline is mandatory:**
  it can only load on a detail page, only on explicit user action, and must be
  zero bytes on every other route.

### Option B — Lighter WebGL wrapper (e.g. regl / OGL / raw WebGL)

- **What:** thin WebGL layer, hand-rolled camera + line drawing.
- **Pros:** ~10–30 KB; real 3D when we want it; far smaller than Three.js.
- **Cons:** we write the camera, picking, axis gizmo, and label placement
  ourselves — non-trivial and exactly the kind of bespoke maintenance burden a
  one-maintainer static site should avoid. Same a11y/animation/canvas problems
  as A with less ecosystem support. Honestly a worse spot than either A or C: it
  has the WebGL downsides without Three's batteries.

### Option C — 2.5D inline SVG (top-down ecliptic, optional z-exaggeration) ★

- **What:** server-or-client-rendered `<svg>` of the ecliptic-plane projection:
  Sun at origin, planet orbits as ellipses/circles, each cycler leg as a sampled
  Kepler polyline, encounter markers as distinct shapes. Inclination shown
  *optionally* via a small "edge-on" companion panel (a second SVG, X–Z
  projection) with an explicit **exaggeration factor** label, rather than a true
  perspective tilt — this keeps z honest instead of letting a casual tilt imply
  precision we don't have.
- **Pros:** **fits the research-minimal aesthetic** (it is line art, like a
  paper figure). **~0 KB of library** — SVG is native; the only JS is a small
  island for pan/zoom/hover/keyboard, in the same style as the existing sort
  island. **Keyboard- and screen-reader-friendly:** every leg/marker is a real
  DOM node with `<title>`/`aria-label`, focusable, tab-navigable — impossible on
  a WebGL canvas without rebuilding it. **Static by default** (no render loop →
  trivially satisfies "no animation"). **Dark mode for free** via
  `currentColor` / CSS variables already defined in `global.css`. **JS-off safe:**
  a build-time-rendered SVG shows the orbit even with no JS; the island only
  *adds* pan/zoom. Crisp at any DPI; trivially exportable as a static figure.
- **Cons:** not "true 3D" — a tilt/rotate gesture would be faked
  (z-exaggerated companion view, not a free camera). For genuinely 3D future
  families (steep planet-centric tours) it would need the edge-on panel to carry
  more weight, or a later upgrade to Option A *scoped to those pages only*.

### Recommendation

**Ship Option C (2.5D inline SVG).** For near-ecliptic heliocentric orbits the
top-down view is the honest, faithful representation, and SVG uniquely satisfies
every site principle (no animation, keyboard, screen-reader, no colour-only
meaning, dark mode, JS-off, near-zero bundle) that a WebGL canvas fights. Treat
"3D" as an **optional edge-on companion panel with a labelled exaggeration
factor**, not a perspective camera. Keep Option A in reserve, *page-scoped and
lazy-loaded*, for a future high-inclination / planet-centric family where the
top-down projection genuinely fails — but do not adopt it for the catalogue's
overwhelmingly near-coplanar majority. The "3D viewer" in spec §15 is best read
as "an orbit viewer that is honest about z," and 2.5D SVG delivers that better
than 3D WebGL does.

---

## 2. Data pipeline — honouring seeds-not-tracks

**Principle:** the catalogue stores *identities* (sparse Kepler seeds), not
sampled trajectories. The viz must reconstruct geometry from
`(a,e,i,Ω,ω,ν?)` and must **not** invent a track where the data is null. Two
sub-decisions:

### 2a. Client-side Kepler sampling vs build-time polylines

| | Client-side sampling | Build-time polylines in JSON |
|---|---|---|
| Bundle | small Kepler sampler in the island (~a few KB of TS) | per-orbit point arrays embedded in page data |
| Page weight | tiny per page | grows with sample density × entries |
| Truth coupling | recomputes from elements on the page (single source) | a second derived artefact to keep in sync |
| Re-projection (tilt/zoom) | trivial — resample/reproject on the fly | needs enough stored points or re-sampling anyway |

**Recommendation:** **client-side Kepler sampling.** Sampling a single ellipse
from `(a,e)` (and orienting by `i,Ω,ω` when present) is a few dozen lines and a
handful of KB — far cheaper than the alternative and it keeps the *only* source
of truth the elements already on the page. This mirrors what `plots.py:182-197`
does, but from elements instead of state vectors (no Lambert/`propagate`
needed for a closed single ellipse — straightforward `(a,e) → r(ν)` polar
sampling; the upstream `propagate` approach is only needed for partial legs with
explicit departure state, which the public rows lack). Build-time sampling is a
fallback only if a future class needs server-side ephemeris we won't ship to the
client.

### 2b. What renders per `cycler_class`

- **single-ellipse** (`[id].astro:202-214`): **full render.** One closed
  heliocentric ellipse from `orbit_elements.(a_au,e)`; orient by
  `inclination_deg / raan_deg / arg_periapsis_deg` *if present*, else draw the
  in-plane ellipse and **label it coplanar-idealized** (§5). Overlay planet
  orbit circles (V/E/M at their semi-major axes, as `plots.py:177-179` does).
  Encounter markers placed from `windows.json` real dates where available, else
  from the published V∞-at-encounters geometry; mark which.
- **multi-arc** (`[id].astro:216-238`): **per-segment render with honest gaps.**
  Draw one arc per `trajectory.segments[]` entry **that has `(a_au,e)`**; for
  segments where those are `null` (common — `types.ts:103` makes them optional)
  **draw nothing and show an explicit gap legend entry** ("segment E→M:
  elements not published"). Never interpolate across a null segment — silent
  interpolation is precisely the provenance violation the brand forbids.
  Cycle-level identity (aphelion ratio / turn ratio) stays in the existing dl,
  not faked into geometry.
- **non-keplerian** (`[id].astro:240-266`): **not renderable from current
  public data.** CR3BP identity is `(jacobi, period_nd, stability)` plus a
  rotating-frame `state_nd` that the public rows do **not** carry as a sampled
  path (and propagating a CR3BP orbit client-side is out of scope). **Say so in
  the UI:** a placeholder panel — "This is a rotating-frame (CR3BP) periodic
  orbit; a faithful render requires numerical propagation in the synodic frame,
  not yet available on the site." No fake ellipse. (If a future build emits a
  sampled `state` path, this becomes renderable; design the placeholder to be
  swappable.)

---

## 3. Scope ladder

- **Phase 1 — per-cycler detail-page view (the obvious win).** A 2.5D SVG on
  `cycler/[id].astro`, class-aware per §2b, lazy/enhancement-only. Self-contained,
  high value, lowest risk. **Recommend shipping this alone first.**
- **Phase 2 — catalogue-wide explorer.** Many orbits overlaid in one view,
  driven by the existing filter/sort island (`CatalogueTable.astro:120-223`).
  Useful but heavier (layout, occlusion, performance with N orbits, a
  legend-design problem). Defer until Phase 1 has shipped and earned it.
- **Phase 3 — Tisserand explorer.** Defer hard. This is the (a,e)/V∞ Tisserand
  graph, conceptually tied to the upstream T-P (Tisserand–Poincaré) graph and
  the M-3D `tisserand_feasible` work
  (`specs/2026-06-05-m-3d-inclination-lift-design.md:264,333`) — a different
  visualization (parameter space, not configuration space) with its own design.
  It does not block, and should not be bundled with, the orbit viewer.

**Recommendation:** Phase 1 only, first. Re-evaluate 2 and 3 after it lands.

---

## 4. Site-principles compliance (the gating checklist)

Every item below is a hard requirement, traced to `global.css:1-2` and the
existing patterns.

- **"No animation":** static-by-default. No render loop, no autoplay rotation,
  no `@keyframes`. Any motion (rotate/zoom/time-scrub) is **user-initiated
  only** and instantaneous-on-input, not continuous.
- **`prefers-reduced-motion`:** add a `@media (prefers-reduced-motion: reduce)`
  block (the stylesheet has none today) that disables even the user-initiated
  transitions, snapping to end-state. Default to reduced-motion-safe behaviour.
- **Keyboard controls:** pan/zoom/focus-next-element fully keyboard-operable,
  same affordance level as the sortable headers (`CatalogueTable.astro:205-218`,
  which already wire `tabIndex` + Enter/Space). Every leg and marker focusable.
- **No colour-only meaning:** distinguish legs by **line style + marker shape +
  text label**, never colour alone (so the dotted/solid/dashed convention from
  `plots.py:179,199` carries over). Encounters get distinct glyphs (○ △ □) plus
  labels.
- **Dark-mode palettes:** draw with `currentColor` / the existing CSS variables
  (`--fg`, `--accent`, `--border`, `--warn-*`), so the viz inverts with the page
  automatically (`global.css:21-35`). No hard-coded hex.
- **Bundle discipline:** the viz script and any sampler load **only on
  `cycler/[id].astro`**, ideally only after user intent (e.g. an "explore"
  toggle / on first focus), and contribute **zero bytes** to home, catalogue,
  launch-windows, about. With Option C there is no library to lazy-load at all —
  just a small island, matching the existing inlined-script model.
- **JS-disabled fallback:** the page must be fully usable with JS off. The
  existing static content (the dl tables, legs list) is the fallback; the viz is
  **enhancement-only**. Build-time-rendered SVG (no-JS) is acceptable and
  preferred; the island only layers on interaction. Mirrors the table island's
  graceful bail (`CatalogueTable.astro:138-144`).

---

## 5. Honesty markers in the UI (the brand)

Provenance is the site's whole brand (the detail page already surfaces
`model_assumption` text, validation level, per-field source quotes —
`[id].astro:135-144,406-422`). The viz must not undercut that.

- **Label every rendered orbit by model fidelity**, reusing the existing
  `model_assumption` vocabulary (`types.ts:11`): *circular-coplanar
  idealization* vs *analytic-ephemeris* vs *cr3bp*. A coplanar-idealized orbit
  must say so on the figure (e.g. a caption badge "idealized: coplanar, planets
  on circles" echoing `[id].astro:142`), not silently imply real geometry.
- **Distinguish idealized vs real encounter placement.** If markers come from
  `windows.json` real DE440-matched dates (`windows.json` header :1-12), label
  them "real ephemeris (DE440)"; if from idealized signature geometry, label
  "idealized." Never blend the two without saying which is which.
- **Gaps are shown, never interpolated.** A null-element multi-arc segment
  renders as an explicit "not published" legend entry, not a guessed curve
  (§2b). A non-keplerian orbit renders as a "not renderable from current data"
  placeholder (§2b), not a stand-in ellipse.
- **No precision the data lacks.** The optional edge-on/z view carries an
  explicit exaggeration-factor label so a casual tilt can't imply we know the
  out-of-plane structure better than i ≤ 3.4° actually constrains.

---

## 6. Open questions for the user (verbatim)

1. Do we accept "2.5D SVG with an honest, labelled edge-on z-panel" as
   satisfying the spec §15 "3D orbit viewer," or do you specifically want a true
   WebGL/Three.js perspective camera (accepting the ~150 KB bundle and the
   keyboard/screen-reader cost)?
2. For Phase 1, is shipping **single-ellipse and multi-arc only** acceptable,
   with non-keplerian (CR3BP) explicitly shown as "not yet renderable," or must
   CR3BP render before we ship anything?
3. Should encounter markers prefer the **real `windows.json` DE440 dates** when
   available (and fall back to idealized geometry otherwise), or always render
   the idealized signature geometry for consistency across rows?
4. Is **client-side Kepler sampling** acceptable (keeps elements as the single
   source of truth), or do you want build-time sampled polylines committed as a
   reviewable, diffable artefact despite the truth-duplication?
5. When orientation angles (`raan_deg`, `arg_periapsis_deg`, `true_anomaly_deg`)
   are `null` (the common case), is drawing the in-plane ellipse with a
   "coplanar-idealized" label the right default, or should we suppress the
   render entirely until oriented elements exist?
6. Should the viz be visible by default on the detail page, or behind an
   explicit "show orbit" toggle (stronger bundle/scroll discipline, weaker
   discoverability)?

---

## Approval (2026-06-06)

User-approved with all recommendations accepted — AND one binding correction
to the design's premises:

- **The "no animation" principle is NOT a user requirement.** It originated
  in the site-bootstrap agent's stylesheet comment (cyclers.space
  `src/styles/global.css:1-2`), not from the user, and must not be treated
  as policy. Animation that genuinely serves the product — e.g. time-scrubbed
  or animated planetary/spacecraft motion along the rendered orbits — is ON
  the table for the final product. `prefers-reduced-motion` respect remains
  as accessibility good practice (the only hard constraint of that family).
  The 2.5D SVG recommendation stands on its other merits (near-ecliptic
  physics, bundle size, DOM accessibility, dark-mode, JS-off safety) — and
  SVG/CSS animation of markers along paths is fully compatible with it.
- Answers as recommended: (Q1) 2.5D SVG + labelled edge-on z-panel satisfies
  the §15 goal for the heliocentric catalogue; Three.js reserved for a
  future genuinely-3D family; (Q2) Phase 1 ships single-ellipse + multi-arc,
  CR3BP as an honest "not yet renderable" placeholder; (Q3) encounter
  markers prefer real windows.json DE440 dates, idealized fallback,
  provenance-labelled; (Q4) client-side Kepler sampling — elements stay the
  single source of truth; (Q5) null orientation angles render the in-plane
  ellipse with a "coplanar-idealized" label; (Q6) visible by default on
  detail pages, lazy where it costs.
- These decisions are EXPECTED TO EVOLVE against the final product — the
  design notes are kept as the record of why, not as immutable law. The
  global.css principles comment should be amended when the viz ships.

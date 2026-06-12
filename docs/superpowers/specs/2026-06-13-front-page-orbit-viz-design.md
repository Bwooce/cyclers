# Front-page hero orbit visualization — design (task #227)

Status: USER-APPROVED design, recorded before implementation. Implementation
lives in the cyclers.space repo (Astro static site); this spec is the upstream
record.

## Goal

Every independently reproduced orbit (validation V1+) in the catalogue, on the
front page, in 3D: a hybrid poster + live canvas, progressively updated at
build time. The hero count ("N independently reproduced orbits — and
counting") is computed from the live filter, never hard-coded.

## 1. Data

Build-time filter of the site's synced catalogue (`src/data/catalogue.yaml`)
on `validation_level in {V1, V2, V3}` (and any future V4/V5). At spec time:
26 rows — V3:2, V2:6 (Aldrin outbound powered + 5 Ross EM CR3BP), V1:18
(incl. 3 Liang Jovian).

Group by system (the row's `primary`, default "Sun") into three scenes:

### Heliocentric scene (primary = Sun; 18 rows)

- `aldrin-classic-em-k1-outbound` (V2) + `-inbound` (V1): single-ellipse rows
  with sourced (a, e) = (1.60, 0.393) — true Kepler curves drawn. NOTE: both
  rows publish the same (a, e) and no Ω/ω, so under the coplanar idealization
  the two curves coincide; the caption says so explicitly rather than
  inventing a relative orientation.
- 16 Russell multi-arc rows (2× V3, 14× V1): NO per-segment (a, e) is
  published (0/236 in the catalogue), so no full curve is drawable. Each row
  DOES carry a sourced max-aphelion (`orbit_elements.aphelion_au`, 1.54–2.22
  AU): rendered as faint reference rings, labelled "sourced max-aphelion ring
  only — per-arc conics unpublished", colour-keyed by tier (V3 vs V1) and
  summarised as count lines in the legend. This reuses the honest multi-arc
  treatment already shipped on the detail pages (commit 128da12).
- Planets: Earth + Mars true J2000 ellipses from the synced Standish elements
  (the existing `planet-elements.json` single source).
- Time-true animation: Earth/Mars and the two Aldrin craft move by the
  existing Kepler clock (`kepler-time.stateAt`) on an idealized phase clock
  (M0 = 0, no epoch asserted) — caption states "idealized phase clock (no
  epoch); encounter timing not asserted" (the per-row encounter re-phasing of
  the detail pages cannot satisfy two craft in one shared scene, so the hero
  does not claim it).

### Earth–Moon scene (primary = Earth; 5 Ross rows, V2)

Geometry lives in `orbit_elements.cr3bp`, not heliocentric ellipses. Every
row carries the full CR3BP identity: sourced (μ, C, T_nd, stability) plus
`state_nd` (DERIVED upstream by the fixed-Jacobi corrector from the sourced
(μ, C) — publication gap recorded in `data_gaps`, #216). That data DOES
support a faithful render: a planar-CR3BP propagation of `state_nd` for one
sourced period in the rotating frame is mechanical reproduction of the row's
own data, not invented geometry. So:

- A small pure-TS PCR3BP integrator (RK4, fixed step) propagates each row at
  render time (client for the gallery, build for the poster). Honesty
  guards: the propagation reports Jacobi drift and closure residual; unit
  tests assert C(state_nd) matches the row's SOURCED `jacobi_constant`
  (sourced-expected only — no golden values our own code computed) and that
  the orbit closes to tolerance at the sourced T.
- Scene frame: Earth–Moon rotating frame, unit = Earth–Moon distance; Earth
  at (−μ, 0), Moon at (1−μ, 0). Caption names the frame, the model
  (PCR3BP), and the state_nd provenance (derived upstream, publication gap).
- Time-true animation: the integrator's time-tagged samples, converted
  TU→days via the row's sourced `tunit_s` (periods 44.8–84.5 d), drive the
  markers — time-true in the rotating frame.

### Jovian scene (primary = Jupiter; 3 Liang rows, V1)

The rows publish V∞ multisets, transit times, and the sequence — but NOT the
moon orbital radii nor per-arc conic elements; reconstructing the idealized
circular-moon-orbits + conic-arcs picture would require external constants
and re-derivation (not straightforward). Per the approved design's honesty
rule ("else badge/count"), the Jovian scene is a badge panel: each row named
with its tier, sequence, V∞ multiset and transit times, plus the explicit
line "geometry not reconstructible from the catalogue's invariants alone —
no curve drawn". No fake curves, ever.

Any future V1+ row whose `primary` is none of Sun/Earth/Jupiter falls into a
generic badge scene (same honesty rule), so the hero never silently drops a
reproduced row and the count always equals the filter.

## 2. Live canvas (gallery)

New sibling module `src/lib/three-gallery.ts` reusing the existing pieces
(`three-axis.toThree`, `kepler-time`, palette/starfield/dispose patterns from
`three-view.ts`; the new `cr3bp-propagate.ts`):

- Auto-cycles the scenes ~8 s each, with manual prev/next buttons (manual
  nav stops the auto-cycle); slow camera orbit per scene.
- Orbits of a group drawn together with distinct colours; legend below the
  canvas carries per-curve validation-tier badges (V3/V2/V1) + fidelity.
- Per-scene honesty caption BELOW the canvas (the 3c79bd9 caption-below
  binding) naming each curve's fidelity.
- prefers-reduced-motion: no auto-cycle, no camera drift, no marker
  animation — static curves + manual prev/next only.
- The scene parameters (Kepler elements, ring radii, CR3BP tuples — a few
  KB, never sampled polylines) are serialised by the build into an inline
  JSON island (the existing clockConfig pattern); the gallery chunk
  regenerates the geometry from them client-side.

## 3. Poster — CHOSEN: build-time SVG montage (no browser in CI)

Of the two approved options (Playwright-CI screenshot vs build-time SVG),
the SVG montage is chosen: it needs no browser in CI, is deterministic,
regenerates from the same data every deploy by construction (an Astro static
endpoint), and cannot rot separately from the data.

- `src/lib/poster-svg.ts`: pure builder — three-panel montage (heliocentric
  top-down / Earth–Moon rotating frame / Jovian badge card), 1200×630,
  self-contained styling, title + live count + per-panel fidelity captions.
- `src/pages/poster.svg.ts`: static endpoint emitting it at `/poster.svg`
  each build. `src/pages/poster.astro`: a `/poster/` route inlining the same
  SVG, suitable for manual screenshots.
- No CI workflow change needed. (og-image PNG would need a rasterizer; out
  of scope for the SVG option — noted, not silently skipped.)

## 4. Front page

`src/components/HeroViz.astro`, mounted at the top of `src/pages/index.astro`:

- The poster `<img src="/poster.svg">` loads instantly (static asset).
- A play button swaps in the lazy-loaded gallery (dynamic
  `import("../lib/three-gallery")` — the existing View-in-3D pattern; zero
  WebGL bytes until click). A close control restores the poster.
- prefers-reduced-motion and mobile default to the poster (the gallery is
  click-only everywhere; under reduced-motion it mounts in static mode).
- Count line computed from the live filter, e.g. "26 independently
  reproduced orbits — and counting", linking to the validation-levels
  explainer.

## File plan (cyclers.space)

- `src/lib/hero-data.ts` — V1+ filter, system grouping, per-row render plan
  (kepler-ellipse | aphelion-ring | cr3bp | badge) with fidelity strings.
- `src/lib/cr3bp-propagate.ts` — pure planar-CR3BP RK4 + Jacobi constant.
- `src/lib/hero-scenes.ts` — JSON-serialisable SceneSpec[] (build-time).
- `src/lib/poster-svg.ts`, `src/pages/poster.svg.ts`, `src/pages/poster.astro`.
- `src/lib/three-gallery.ts` — lazy gallery renderer.
- `src/components/HeroViz.astro`, edits to `src/pages/index.astro`,
  `src/styles/global.css`.
- Tests: `hero-data`, `cr3bp-propagate` (sourced-expected Jacobi check),
  `hero-scenes`, `poster-svg` under `src/lib/__tests__/`.

## Honesty rules (binding, project law)

- Never render geometry the data doesn't support; representative renders are
  labelled as such; rows without honest geometry appear as badges/counts.
- Per-curve fidelity captions, below the canvas, always visible with it.
- The count comes from the live filter; nothing hard-codes 26.
- The public repo never references the private papers repo or local paths.

# Viz Phase 2b — the flying 3D camera (cyclers.space)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:subagent-driven-development`
> or `superpowers:executing-plans`. Checkbox steps; strict TDD where a test is
> possible (write failing test → run **red** → minimal impl → run **green** →
> commit), and a **built-HTML + browser** verification step where the artefact is
> WebGL canvas / camera behaviour that a unit test cannot reach. Work on `main` —
> **do NOT branch** (project rule). **Do NOT commit while this plan is authored
> under the docs-only mandate**; the commit messages below are the messages to
> *use* when the implementation phase runs.
>
> **This plan targets the SEPARATE repo `/home/bruce/dev/cyclers.space`** (its own
> git root, verified), not the `cyclers` modelling repo this plan file lives in.
> Every path below is relative to `cyclers.space/` unless prefixed otherwise.
> Toolchain there: Astro 6 + TypeScript, npm. Lint/type gate before **every**
> commit: `npm run astro check` (the `@astrojs/check` + `tsc` gate the project
> already ships). No test runner is installed today; **Task 0.0 adds a minimal
> `vitest` dev-dep** for the pure-TS unit tests this plan needs (the camera
> physics/keyframe maths), and browser checks use the Playwright MCP against
> `npm run preview`.
>
> This plan is the task-level expansion of the **approved** design
> `cyclers/docs/superpowers/specs/2026-06-07-viz-phase2-timetrue-flying-camera-design.md`
> (read it **including its Approval section** — §4 is the camera, and the seven
> Approval resolutions are binding). Where this plan and the design could diverge,
> **the design + its Approval win.** Phase 2a has **shipped** (see Survey); 2b is
> the §7 Option-1 "camera as its own plan" follow-up.

---

## Goal

Add an **opt-in, lazy-loaded Three.js 3D camera** to the per-cycler orbit view
that consumes the **same shared clock** (`kepler-time.ts`) the shipped 2a SVG
already runs on, so the 3D scene and the 2D SVG can never disagree about where
any body is at instant *t* (one clock, two renderers). The 3D view is **strictly
additive**: the SVG remains the accessible source of truth; Three.js ships **zero
bytes** until the user clicks "View in 3D" (lazy on intent, not merely on route).

Three camera modes, smallest-shippable-slice first:

1. **Orbit-cam** (the default) — slow user-or-auto orbit around the whole system,
   trajectory fully inked, **paused, time at the first encounter** on open.
2. **Chase-cam** (power-user) — rides the spacecraft, looks along velocity;
   **reduced-motion forces orbit-cam** (motion-sickness guard).
3. **Guided tour** (second button) — a scripted keyframe path over the didactic
   beats (departure → flyby → aphelion → return phasing), keyframes **derived
   per-cycler-class from the geometry**, never hardcoded per row.

Honest labelling carried into 3D (model badges, gap segments NOT drawn with
invented geometry, the CR3BP "not renderable" equivalent). The user's note stands:
tilt shipped but is physically subtle (V/E/M all ≤ 3.4° of ecliptic) — **the
camera is where the 3D experience actually lands**, and animation is fully
permitted (user-controlled, reduced-motion-aware).

---

## Survey — what 2a actually shipped (verified against live code, cite file:line)

2a is real and shipped in `cyclers.space`. 2b builds **on top of it** and must not
regress it.

- **The shared clock module** — `src/lib/kepler-time.ts` (277 lines). Pure TS,
  framework-free, the **one home for time → position**. 2b MUST consume it; it
  writes no new physics. Live exports 2b reuses (verified):
  - `interface Vec3 {x,y,z}` (`kepler-time.ts:21-25`) — **3D already present**;
    the SVG drops z, the 3D scene will *use* it.
  - `interface KeplerElements` (`kepler-time.ts:28-37`) — `a,e,i_deg,lan_deg,
    argp_deg,M0_deg,n_deg_per_day?,t_epoch_day?`.
  - `stateAt(el, t): Vec3` (`kepler-time.ts:111-121`) — the per-frame call for
    every moving body. **Returns position only** (no velocity).
  - `samplePath(el, nSamples=240): Vec3[]` (`kepler-time.ts:136-144`) — the inked
    trajectory / planet orbit line geometry, time-uniform.
  - `periodDays(el)` (`:124-127`), `distance(a,b)` (`:147-152`),
    `proximitySeries(...)` (`:260-277`), `planetToElements(record)` (`:181-192`),
    `rephaseToEncounter(...)` (`:203-229`), `isoToJ2000Days(iso)` (`:233-235`).
- **The mount point + the synchronisation contract** — `src/components/OrbitView.astro`
  (661 lines). 2b reuses, in particular:
  - **`clockConfig`** (`OrbitView.astro:256-266`), serialised into an inline
    `<script type="application/json" data-orbit-config="orbit-svg-${id}">`
    (`:377`). Fields: `regime, t0, t1, craft (KeplerElements), planets[{code,el}],
    scale, cx, cy, bodies[]`. **This is the synchronisation contract:** the 2a
    island reads it (`:517-525`) to drive Kepler-true motion; **2b's camera reads
    the SAME JSON** so both renderers share `craft`, `planets`, `t0/t1`. 2b adds
    NO new clock state — it parses the existing block.
  - The clock-regime + provenance machinery: `clockRegime`
    (`:93`), `clockT0/clockT1` (`:144-145`), `encMarks` time-true markers
    (`:162-168`), `proximity` series (`:149-154`), `clockLabel`/`encProvenance`
    (`:169-176`), `fidelityBadge` (`:217-223`). 2b's on-canvas captions reuse
    these strings verbatim (provenance is the brand).
  - The 2a island (`<script>` `:487-661`) already: parses `clockConfig`, has
    `setTime(t)` driving `stateAt` for craft+planets, a `proxLive` live region,
    a `reduced` check (`:538`), play/scrub, and `initTilt` (`:629-653`). 2b adds
    a sibling lazy module; it does **not** rewrite this island.
- **Planet elements (sourced)** — `src/data/planet-elements.json` (95 lines):
  all 8 planets, J2000 osculating `(a,e,i,Ω,ϖ,L0,n)`, citation "J2000 osculating
  ellipse (Standish & Williams Table 1)". Consumed via
  `src/lib/orbit.ts` `PLANETS` (`orbit.ts:41-43`), `samplePlanetEllipse`
  (`orbit.ts:46-50`), `PLANET_GEOMETRY_CITATION` (`orbit.ts:32`). 2b reuses these
  for the 3D planet orbit lines and markers.
- **The tilt implementation (stays SVG-side)** — `initTilt`
  (`OrbitView.astro:629-653`) applies an SVG `matrix(1 0 0 k 0 …)` Y-scale to the
  `[data-tilt-group]` (`:309`); control at `:393-401`. This is the **2.5D
  understudy** and remains pure-SVG, zero-bundle, for the 99% case. 2b does **not**
  touch it; the camera is the opt-in upgrade *next to* it.
- **The detail-page mount** — `src/pages/cycler/[id].astro:269` renders
  `<OrbitView entry={entry} />`. 2b's "View in 3D" button lives inside
  `OrbitView.astro` (so it travels with the component), gated to
  `cls === "single-ellipse"` initially (multi-arc/non-keplerian get honest
  "not available in 3D" notes — Slice 1 scope).
- **Bundle baseline (the zero-cost proof)** — `package.json` `dependencies` are
  exactly `@astrojs/check`, `astro`, `typescript` (no runtime framework, **no
  WebGL**). This is the invariant Task 1.1 measures against. Build = `astro build`;
  `predev`/`prebuild` run `scripts/sync-catalogue.mjs`.
- **A11y + theming primitives 2b must honour** — `src/styles/global.css`:
  `@media (prefers-reduced-motion: reduce)` block (`global.css:511-514`, already
  scoped to `.orbit-*`), and the dark-mode CSS variables under
  `@media (prefers-color-scheme: dark)` (`global.css:23+`, `--fg/--bg/--accent/…`).
  The 3D scene reads `matchMedia` for both at init (two material sets;
  no auto-motion under reduced-motion).

### Axis-convention decision (document against kepler-time.ts) — BINDING

`kepler-time.ts` works in **heliocentric J2000-ecliptic**, right-handed, Sun at
origin, units **AU**: `stateAt` returns `Vec3{x,y,z}` where the **ecliptic plane
is z = 0**, +z is ecliptic north (`perifocalToEcliptic`, `kepler-time.ts:89-103`,
`z = sinw·sinI·xp + cosw·sinI·yp`). The SVG today projects top-down by taking
`(x, y)` and flipping to screen with `cy − y·scale` (`orbit.ts:193-197`).

Three.js's default camera looks down **−Z** with **+Y up**. Mapping the physics
frame directly into Three's would put the ecliptic on Three's XY plane and the
camera would look at its edge — wrong. **Decision: map physics → Three as a
single fixed swap so the ecliptic is Three's ground plane (XZ) and +Y is ecliptic
north:**

```
three.x =  phys.x      // ecliptic x  → Three x
three.y =  phys.z      // ecliptic north (+z) → Three up (+y)
three.z = -phys.y      // ecliptic y  → Three -z  (preserves right-handedness)
```

One pure function `toThree(v: Vec3): THREE.Vector3` (Task 1.2) is the **only**
place this swap lives — the camera, the orbit lines, and the markers all route
through it, exactly as the SVG routes through `toSvgPath`. Rationale recorded in
the function's doc comment. This keeps "slow at aphelion" reading correctly from
an orbit-cam looking down on the ecliptic, and makes the ≤3.4° tilt finally
*visible* as real out-of-plane height (the y-axis), which the SVG could only fake
via the exaggerated edge-on panel. Scale: reuse `clockConfig.scale` semantics
but in 3D we work in **AU directly** (no px/AU) and frame the camera from the
trajectory's aphelion radius (Task 1.4) — document that 3D drops the SVG's px
scale and uses AU world units.

---

## Architecture

### One clock, two renderers — the synchronisation contract (design §4.4)

```
            src/data/planet-elements.json   (sourced, Standish Table 1)
                          │
                   src/lib/orbit.ts  (PLANETS, samplePlanetEllipse)
                          │
   OrbitView.astro build  │  computes clockConfig {regime,t0,t1,craft,planets,scale,cx,cy,bodies}
                          ▼
        <script type="application/json" data-orbit-config="orbit-svg-${id}">   ← THE CONTRACT
                ├────────────────────────────────┬───────────────────────────────
                ▼ (2a, shipped)                   ▼ (2b, this plan)
        SVG island setTime(t)             three-view.ts setTime(t)
        stateAt(craft/planets, t)         stateAt(craft/planets, t) → toThree()
        → SVG attrs                        → THREE.Vector3 positions
                          \________ both import kepler-time.ts ________/
```

2b adds **one new lazy module** `src/lib/three-view.ts` (the scene + camera
controller) and **one minimal controls helper** `src/lib/three-controls.ts`
(hand-rolled damped orbit/chase controls — NOT `three/examples` OrbitControls, to
keep the lazy chunk small, design §4.1). It reads the **existing** `clockConfig`
JSON; it introduces no second source of truth. A shared `t` cursor (a tiny
in-page event/exported setter) keeps the SVG marker and the 3D camera on the same
instant when both are visible, so switching to 3D "continues the same instant in
time" (design §4.2).

### Lazy-on-intent loading (design §4.1, non-negotiable)

`three` is added to `package.json` `dependencies` (Task 1.0) but **never
statically imported** by any `.astro` frontmatter or any eagerly-loaded script.
The only entry point is, inside a click handler on the "View in 3D" button:

```ts
const THREE = await import("three");
const { mountThreeView } = await import("../lib/three-view");
```

Astro/Vite code-splits `three` + `three-view` + `three-controls` into a chunk
fetched **only** on that click. Task 1.1 proves the before/after page weight:
the detail page's initial JS payload is unchanged; the three chunk appears in the
network log **only after** the click. `three-controls` is hand-rolled so we never
pull the examples bundle.

### Camera modes (design §4.2)

| Mode | Frame | Look-at | Motion | Default? |
|---|---|---|---|---|
| **orbit-cam** | world (Sun-centred) | system centroid | slow azimuth orbit (user/auto), damped | **YES** — paused, inked, t = first encounter |
| **chase-cam** | spacecraft position | along velocity (finite-diff of `stateAt`) | follows craft | power-user; **disabled under reduced-motion** |
| **guided-tour** | scripted | scripted | keyframe lerp over beats | second button |

Velocity for chase-cam look-along: `kepler-time.ts` `stateAt` returns position
only (`:109-110` notes "velocity is not needed … a finite difference suffices").
2b computes look direction as `normalize(stateAt(craft, t+δ) − stateAt(craft, t))`
with a small δ (e.g. `period/4000` days) — the same trick the proximity readout
implies. **Flag (refactor candidate, do NOT do here):** if chase-cam jitter from
finite-diff proves visible, the honest fix is a `velocityAt` export in
`kepler-time.ts` (closed-form from the vis-viva / perifocal velocity). Recorded
in the Self-review; not in 2b scope.

### Guided tour keyframes — derived from geometry, not hardcoded per row (design §4.2)

A `tourKeyframes(cfg)` pure function (Task 3.0) builds `(t, cameraPose, caption)`
tuples **per cycler-class** from the shared config:

- **departure** — `t = clockT0`; camera frames Earth + spacecraft start node
  (`stateAt(craftEl, t0)`, Earth via `planets['E']`); caption = "Earth departure".
- **flyby** — `t` = the proximity-series **minimum** for the primary visited
  non-Earth body (reuse `proximitySeries` minima already computed at 2a build,
  re-exposed in `clockConfig` — see Task 1.0b); camera dollies in to the
  encounter; caption reuses `encProvenance`.
- **aphelion** — `t` where `stateAt(craft,t)` radius is max over `[t0,t1]`
  (scan); camera pulls back; caption "long aphelion arc — time visibly slows
  (Kepler's second law)". This is the payoff beat.
- **return phasing** — `t` near `clockT1`; camera frames Earth + converging
  spacecraft; caption "Earth-return phasing".

Keyframes are **data computed from the geometry**, so any single-ellipse row gets
a correct tour with no per-row authoring. Multi-arc/non-keplerian: no tour
(honest — Slice 1 excludes them from 3D entirely).

### Honesty layer carried into 3D (design §5)

- **Model badges** — the on-canvas caption overlay shows `fidelityBadge`,
  `PLANET_GEOMETRY_CITATION`, and `clockLabel` (same strings as the SVG
  figcaption `OrbitView.astro:380-390`).
- **Gaps stay gaps** — Slice 1 ships single-ellipse only. When multi-arc lands
  (out of 2b scope, flagged), null-`tof_days` legs are drawn as **nothing** with
  an explicit "elements not published — not shown in 3D" note, mirroring the SVG
  gap discipline (`OrbitView.astro:450-459`). The camera never flies an unknown
  leg.
- **CR3BP placeholder equivalent** — non-keplerian rows show NO "View in 3D"
  button (a rotating-frame orbit is not renderable from heliocentric elements,
  mirroring `OrbitView.astro:281-293`); instead a one-line note "3D view not
  available for rotating-frame (CR3BP) orbits."
- **Idealized vs anchored** — the clock-regime label (`clockLabel`) rides into
  the 3D caption unchanged; an idealized phase clock never displays a fabricated
  date in 3D any more than in 2D.

### Accessibility (design §4.3) — designed, not waved at

- **Strictly additive** — the SVG (focusable DOM nodes, `<title>`/`aria-label`,
  the dl tables) remains the source of truth. WebGL unavailable/declined ⇒ the
  page is the SVG, unchanged.
- **Focus management** — the "View in 3D" button is a normal focusable control;
  on activate, focus moves into the canvas (`tabindex="0"`, visible focus ring via
  `:focus-visible`); **`Esc` returns focus to the SVG** and hides the canvas. A
  "Close 3D" button is the visible equivalent. Entering/leaving is announced.
- **Keyboard camera control** — full key map (table below), no mouse-only
  affordance, with an on-canvas key-help overlay toggled by `?`.
- **Live region** — an off-canvas `aria-live="polite"` element announces camera
  mode + current didactic beat + proximity ("Approaching Mars flyby;
  spacecraft–Mars distance 0.02 AU"), driven by the same `kepler-time.ts` state,
  reusing the existing `[data-prox-live]` pattern (`OrbitView.astro:423`).
- **prefers-reduced-motion** — **no auto-orbit, no tour autoplay; chase-cam
  disabled** (forces orbit-cam); all transitions snap (camera pose set
  instantly, no lerp). The button label notes "(reduced-motion: manual step)".
  Consistent with `global.css:511-514`.
- **Dark mode** — scene reads `prefers-color-scheme` at init; two material sets
  (clear-color, line, label) mirroring the CSS variables (`global.css:23+`).

### Keyboard binding table (BINDING — design §4.3, designed concretely)

| Key | Action | Modes | reduced-motion |
|---|---|---|---|
| `←` / `→` | orbit camera azimuth − / + | orbit-cam | step (snap) |
| `↑` / `↓` | orbit camera elevation + / − | orbit-cam | step (snap) |
| `+` / `=` | dolly in (zoom toward target) | all | step |
| `-` / `_` | dolly out | all | step |
| `[` | step time backward (one beat / Δt) | all | yes (primary) |
| `]` | step time forward | all | yes (primary) |
| `Space` | play / pause time | all | **inert** (no autoplay); announces "manual step only" |
| `1` | switch to orbit-cam | all | yes |
| `2` | switch to chase-cam | (ignored under reduced-motion → stays orbit-cam, announces) | n/a |
| `T` | start / stop guided tour | all | **inert** (no autoplay); `[`/`]` walk beats instead |
| `?` | toggle on-canvas key-help overlay | all | yes |
| `Esc` | exit 3D, return focus to the SVG | all | yes |

Pointer drag (orbit), wheel (dolly), and pinch mirror the arrow/dolly keys for
sighted mouse/touch users; every pointer action has a keyboard equivalent above.

---

## Smallest shippable slice (BINDING)

**Slice 1 = orbit-cam only, single-ellipse rows only.** Concretely: the "View in
3D" button (lazy import on click) → a Three scene with Sun, the visited planets'
sourced J2000 **orbit lines** (`samplePlanetEllipse` → `toThree`), the cycler's
**inked trajectory** (`samplePath(craft)` → `toThree`), planet + spacecraft
**markers** driven by `stateAt` on the shared clock, **paused at t = first
encounter**, an **orbit-cam** with keyboard + damped pointer controls, the
honesty caption overlay, dark-mode + reduced-motion handling, focus management
(Esc returns to SVG), and the zero-bytes-until-click bundle proof. **No chase-cam,
no tour, no multi-arc.** This is independently shippable and delivers "the 3D
idea" the user is missing. Slices 2 (chase) and 3 (tour) layer on without
touching Slice 1's contract.

---

## Phasing (independently shippable)

| Slice | Theme | Tasks | three dep | new lazy chunk | browser check |
|---|---|---|---|---|---|
| **0** | Test harness + the toThree axis fn (pure, no WebGL) | 0.0–0.1 (2) | no | no | no |
| **1** | Orbit-cam MVP: lazy load, scene, inked trajectory, shared-clock markers, a11y, bundle proof | 1.0–1.7 (8) | yes | yes | yes |
| **2** | Chase-cam (look-along-velocity; reduced-motion forces orbit-cam) | 2.0–2.2 (3) | yes | yes | yes |
| **3** | Guided tour (geometry-derived keyframes, live-region narration) | 3.0–3.2 (3) | yes | yes | yes |

Slice 0 is foundational (pure-TS, unit-testable, gates the axis convention before
any WebGL). Slice 1 is the smallest shippable 3D experience. **Total: 16 tasks.**

---

## Slice 0 — test harness + the axis-convention function (pure TS)

### Task 0.0 — add a unit-test runner (vitest) for the pure-TS 2b maths

**Files:** `package.json` (devDependencies + a `test` script);
`vitest.config.ts`; `src/lib/__tests__/three-axis.test.ts` (placeholder that
imports nothing yet, asserting the runner works).

The project ships no test runner. The 2b physics/keyframe maths (axis swap, tour
keyframes, camera framing) are pure functions and MUST be unit-tested. Add
`vitest` (+ `jsdom` only if a later DOM test needs it; Slice-0/3 maths are pure,
so plain node env first). `astro check` stays the type/lint gate.

#### Failing test — `src/lib/__tests__/three-axis.test.ts`

```ts
import { describe, it, expect } from "vitest";

describe("vitest harness", () => {
  it("runs", () => {
    expect(1 + 1).toBe(2);
  });
});
```

Run: `npm test` → **red** (no vitest). Add `vitest` to devDependencies, a
`"test": "vitest run"` script, `vitest.config.ts` (node env). `npm install`.
Run → **green**. `npm run astro check`. Commit:

```
viz-2b: add vitest harness for pure-TS camera maths (viz phase 2b slice 0)
```

### Task 0.1 — `toThree` axis convention (the single frame swap)

**Files:** create `src/lib/three-axis.ts`; test `src/lib/__tests__/three-axis.test.ts`.

`toThree(v: Vec3): {x,y,z}` implementing the BINDING swap
(`three.x=phys.x, three.y=phys.z, three.z=-phys.y`) with the rationale in the doc
comment (ecliptic = Three ground plane XZ, +Y = ecliptic north, right-handed
preserved). Returns a plain `{x,y,z}` (so the function is testable with **zero**
three import — three is loaded lazily; the maths is not). `three-view.ts` later
wraps it in `new THREE.Vector3(...)`.

#### Failing test — `src/lib/__tests__/three-axis.test.ts`

```ts
import { describe, it, expect } from "vitest";
import { toThree } from "../three-axis";
import type { Vec3 } from "../kepler-time";

describe("toThree axis convention (ecliptic z=0 -> Three ground plane XZ)", () => {
  it("maps ecliptic north (+z) to Three up (+y)", () => {
    const v: Vec3 = { x: 0, y: 0, z: 1 };
    expect(toThree(v)).toEqual({ x: 0, y: 1, z: 0 });
  });
  it("maps ecliptic +x to Three +x", () => {
    expect(toThree({ x: 2, y: 0, z: 0 })).toEqual({ x: 2, y: 0, z: 0 });
  });
  it("maps ecliptic +y to Three -z (right-handed)", () => {
    expect(toThree({ x: 0, y: 3, z: 0 })).toEqual({ x: 0, y: 0, z: -3 });
  });
});
```

Run → **red** → impl `toThree` (no three import) → **green**. `astro check`.
Commit:

```
viz-2b: toThree ecliptic->Three axis convention, the single frame swap (viz phase 2b slice 0)
```

---

## Slice 1 — orbit-cam MVP

### Task 1.0 — add `three` as a dependency (NOT imported anywhere yet)

**Files:** `package.json` (dependencies).

Add `"three": "^0.180.0"` and `"@types/three"` (devDependencies) to
`package.json`; `npm install`. **No import is added** — this task only makes the
package resolvable for the lazy `await import("three")` in Task 1.3. Verify (a
guard) that no `.astro`/`.ts` statically imports three yet.

- [ ] `npm install` succeeds; `three` resolvable.
- [ ] `grep -rn "from \"three\"\|from 'three'\|import(\"three\")" src/` shows
  **no** static `from "three"` (the dynamic `import("three")` arrives in 1.3).
- [ ] `npm run astro check` passes.

Commit:

```
viz-2b: add three + @types/three (lazy-only; no static import) (viz phase 2b slice 1)
```

### Task 1.0b — re-expose the encounter/aphelion clock data in `clockConfig`

**Files:** `src/components/OrbitView.astro` (extend the `clockConfig` object only).

The tour (Slice 3) and the orbit-cam default ("t = first encounter") need the
encounter times the SVG already computes. The 2a build computes `encounterTimes`
(`OrbitView.astro:126-136`) and `proximity` minima (`:149-154`) but does **not**
put them in `clockConfig`. Add two fields (additive — the 2a island ignores
unknown keys, verified at `:521`): `encounterTimes: number[]` and
`proximityMinima: {body,t,d_au}[]` (from `proximity.map(s => ({body:s.body,
t:s.minimum.t, d_au:s.minimum.d_au}))`). **No behaviour change to 2a.**

#### Verification (no new unit test needed — type-checked + browser-confirmed)

- [ ] `npm run astro check` passes (the `clockConfig` object stays typed).
- [ ] Built HTML: `npm run build` then inspect a single-ellipse detail page's
  `data-orbit-config` JSON contains `encounterTimes` and `proximityMinima`.
- [ ] The 2a SVG island still plays/scrubs unchanged (browser check, Task 1.7).

Commit:

```
viz-2b: surface encounterTimes + proximity minima in clockConfig for the camera (viz phase 2b slice 1)
```

### Task 1.1 — bundle baseline + the zero-cost-when-unused proof

**Files:** create `docs/viz-2b-bundle.md` in `cyclers.space` (a short measured
record — this is a code-repo artefact, allowed; it is NOT the plan doc). *(If
shared-doc concurrency forbids adding a file at run time, record the numbers in
the commit body instead and flag it.)*

Measure BEFORE any three import lands on the page: `npm run build`, record the
detail page's shipped JS (the `dist/_astro/*.js` the `cycler/[id]` page
references) and total. This is the baseline the post-Task-1.3 network check
compares against. The **proof obligation** (closed in Task 1.7): the detail page's
initial payload is byte-identical pre/post-button, and the `three` chunk loads
**only** after the click.

- [ ] Baseline JS bytes for `cycler/[id]` recorded.
- [ ] Note states the acceptance: three chunk (~150 KB tree-shaken core, design
  §4.1) is acceptable **only** behind the click; documented in the button's
  `title`.

Commit:

```
viz-2b: record detail-page JS bundle baseline (zero-cost-when-unused proof) (viz phase 2b slice 1)
```

### Task 1.2 — the "View in 3D" button + lazy mount shell (no scene yet)

**Files:** `src/components/OrbitView.astro` (button markup, gated to
`single-ellipse`; a sibling lazy-init script); a `<div class="orbit-3d" hidden>`
canvas host.

Add, only when `cls === "single-ellipse" && craftElems`: a `<button
class="orbit-3d-btn">View in 3D</button>` with `title` noting the ~150 KB cost,
a hidden `<div class="orbit-3d" data-orbit-3d="orbit-svg-${id}">` host, and a
small **eagerly-loaded** init script that ONLY wires the click handler — the
handler body does `await import("../lib/three-view")` then `mountThreeView(host,
cfg)`. For non-single-ellipse: render the honest note (no button). This task
ships the button + the dynamic-import wiring with a **stub** `three-view.ts` that
just sets `host.textContent = "3D loading…"` (proves the lazy path before the
scene exists).

#### Verification (browser, via Playwright MCP against `npm run preview`)

- [ ] On a single-ellipse page the button renders; on a CR3BP page it does not
  (the note does).
- [ ] Clicking the button reveals the host and triggers the dynamic import
  (network log shows the `three-view` chunk fetched **on click**, not before).
- [ ] `npm run astro check` passes.

Commit:

```
viz-2b: View-in-3D button + lazy import shell (stub scene) (viz phase 2b slice 1)
```

### Task 1.3 — the scene: Sun, sourced planet orbit lines, inked trajectory

**Files:** create `src/lib/three-view.ts` (`mountThreeView(host, cfg)`); it does
`const THREE = await import("three")` internally (or receives it — decide for
chunking; prefer importing inside so the button handler stays tiny).

Build the scene from the parsed `clockConfig`: Sun sphere at origin; for each
`cfg.planets[].code` a **line** from `samplePlanetEllipse(code)` → `toThree` (AU
world units); the cycler **inked trajectory** from `samplePath(cfg.craft)` →
`toThree`; a faint ecliptic grid + starfield for orientation (design §4.2 chase
guard, also helps orbit-cam). Dark/light material sets chosen from
`matchMedia("(prefers-color-scheme: dark)")`. Renderer sized to the host;
`renderer.dispose()` on close. **No camera controller yet** (fixed look-down
camera) — that is Task 1.4.

#### Verification (browser + a thin pure-TS unit where possible)

- [ ] Pure unit (no WebGL): a `buildOrbitLinePoints(cfg)` helper extracted from
  the scene builder returns the expected point count and routes every point
  through `toThree` (assert one known planet point maps correctly). Test in
  `src/lib/__tests__/three-view-geometry.test.ts`.
- [ ] Browser: clicking "View in 3D" shows Sun + planet ellipses + the cycler
  curve; the ellipses are visibly **eccentric** for Mars (not circles).
- [ ] No console errors; `astro check` passes.

Commit:

```
viz-2b: Three scene — Sun, sourced planet orbit lines, inked trajectory (viz phase 2b slice 1)
```

### Task 1.4 — orbit-cam controller (damped pointer + keyboard) framed to aphelion

**Files:** create `src/lib/three-controls.ts` (hand-rolled spherical orbit
controls with damping); `three-view.ts` wires it.

Hand-rolled (NOT `three/examples`): spherical `(radius, azimuth, elevation)`
around the system centroid, pointer-drag → azimuth/elevation, wheel → radius,
critically-damped easing toward targets. Initial `radius` framed from the
trajectory's **aphelion** (`max |stateAt(craft,t)|` over `[t0,t1]`, AU) so the
whole orbit fits. Keyboard bindings per the table: `←→↑↓` orbit, `+-` dolly. The
camera up-vector is Three +Y (ecliptic north) so "looking down on the ecliptic"
is the default pose.

#### Verification

- [ ] Pure unit: `frameRadiusAU(cfg)` returns ~aphelion (assert against a
  hand-computed `a(1+e)` for a fixture craft). `cameraPoseFromSpherical(r,az,el)`
  returns expected Vector3-like coords (test the maths without three by returning
  plain `{x,y,z}` and constructing the Vector3 at the call site).
- [ ] Browser: drag orbits smoothly (damped), wheel zooms, arrow keys orbit,
  `+/-` dolly; the full orbit is in frame on open.

Commit:

```
viz-2b: hand-rolled damped orbit-cam controller framed to aphelion (viz phase 2b slice 1)
```

### Task 1.5 — shared-clock markers, paused at first encounter; SVG<->3D time sync

**Files:** `three-view.ts` (per-frame marker update + the `setTime` hookup).

Spacecraft + planet **markers** as small spheres; each frame (or on time-step)
`stateAt(el, t)` → `toThree` → marker position, the SAME calls the SVG island
makes. **Default t = first encounter** (`cfg.encounterTimes[0]` if present, else
`cfg.t0`), **paused**. Wire a shared `t` cursor so stepping time in 3D (`[`/`]`,
`Space`) also moves the SVG craft marker and vice-versa (a tiny exported
`setOrbitTime(svgId, t)` the 2a island also calls, or a `CustomEvent` on the
host) — switching to 3D continues the same instant (design §4.2). `Space`
play/pause **only when not reduced-motion**; under reduced-motion `Space` is inert
and `[`/`]` step.

#### Verification

- [ ] Pure unit: a `markerWorldPos(el, t)` = `toThree(stateAt(el,t))` round-trips
  a fixture (craft at `t0` equals the SVG's `stateAt(craft,t0)` mapped).
- [ ] Browser: on open the scene is paused with the spacecraft at the first
  encounter; `]` advances time and BOTH the 3D craft and the SVG craft move
  together; aphelion visibly takes more wall-time per AU than periapsis when
  playing (Kepler's second law — the payoff).

Commit:

```
viz-2b: shared-clock markers paused at first encounter + SVG<->3D time sync (viz phase 2b slice 1)
```

### Task 1.6 — a11y: focus, key-help overlay, live region, reduced-motion, dark mode, honesty caption

**Files:** `three-view.ts` (canvas `tabindex`, key handlers, `Esc` exit, `?`
overlay, `aria-live` wiring, caption overlay); `OrbitView.astro` (the off-canvas
live region + caption strings passed in `cfg`); `src/styles/global.css` (focus
ring, `.orbit-3d` reduced-motion rules).

Implement the full a11y surface from Architecture: canvas `tabindex="0"` +
`:focus-visible` ring; `Esc` hides the canvas and returns focus to the SVG; `?`
toggles an on-canvas key-help overlay listing the binding table; an off-canvas
`aria-live="polite"` element announces mode + beat + proximity (reuse the
`[data-prox-live]` pattern); under reduced-motion no auto-orbit/no autoplay and
all camera moves snap; dark/light material sets from `matchMedia`. The honesty
caption overlay shows `fidelityBadge` + `PLANET_GEOMETRY_CITATION` + `clockLabel`
(pass these strings into `cfg` from the build — they already exist at
`OrbitView.astro:217-223, 173-176, 32`).

#### Verification

- [ ] Browser (Playwright MCP): Tab reaches the button; activating moves focus
  into the canvas; `Esc` returns focus to the SVG; `?` shows the overlay; the
  live region text updates on time-step; the caption shows the model/clock
  provenance strings.
- [ ] Emulate `prefers-reduced-motion: reduce`: no auto-orbit, `Space` inert,
  camera snaps; `[`/`]` still step.
- [ ] Emulate dark mode: scene uses the dark material set.

Commit:

```
viz-2b: 3D a11y — focus, key-help, live region, reduced-motion, dark mode, honesty caption (viz phase 2b slice 1)
```

### Task 1.7 — Slice-1 gate: bundle proof, regression, full browser pass

**Files:** none new; close `docs/viz-2b-bundle.md`; type/lint gate.

- [ ] **Zero-cost proof:** `npm run build`; the `cycler/[id]` initial JS is
  byte-identical to the Task-1.1 baseline; the `three`/`three-view`/`three-controls`
  chunk is fetched **only after** the "View in 3D" click (network log).
- [ ] **No 2a regression:** the SVG renders, the play/scrub island works, the
  tilt slider works, the proximity sparkline + live readout work — all unchanged.
- [ ] **Other routes untouched:** home / catalogue / launch-windows / about ship
  zero new bytes (no `three` chunk referenced).
- [ ] `npm run astro check` passes; `npm test` green.

Commit:

```
viz-2b: Slice 1 gate — orbit-cam MVP, zero-cost bundle proof, no 2a regression (viz phase 2b slice 1)
```

---

## Slice 2 — chase-cam

### Task 2.0 — look-along-velocity via finite-difference of stateAt

**Files:** `three-view.ts` (a `chaseLookDir(craft, t)` helper); test
`src/lib/__tests__/chase-look.test.ts`.

`chaseLookDir(craft, t) = normalize(stateAt(craft, t+δ) − stateAt(craft, t))`,
δ = `periodDays(craft)/4000`, mapped through `toThree`. Pure-testable (no three):
return a plain unit `{x,y,z}`.

#### Failing test

```ts
import { describe, it, expect } from "vitest";
import { chaseLookDir } from "../three-view-chase"; // pure helper module
import type { KeplerElements } from "../kepler-time";

const craft: KeplerElements = { a: 1.6, e: 0.2, i_deg: 0, lan_deg: 0, argp_deg: 0, M0_deg: 0 };

describe("chaseLookDir", () => {
  it("returns a unit vector", () => {
    const d = chaseLookDir(craft, 10);
    const m = Math.hypot(d.x, d.y, d.z);
    expect(m).toBeCloseTo(1, 6);
  });
  it("points roughly prograde (changes smoothly along the orbit)", () => {
    const a = chaseLookDir(craft, 0);
    const b = chaseLookDir(craft, 1);
    expect(a).not.toEqual(b); // direction evolves
  });
});
```

Run → **red** → impl `chaseLookDir` in a pure helper `three-view-chase.ts` →
**green**. `astro check`. Commit:

```
viz-2b/chase: look-along-velocity via finite-diff of stateAt (viz phase 2b slice 2)
```

### Task 2.1 — chase-cam mode + reduced-motion forces orbit-cam

**Files:** `three-view.ts`, `three-controls.ts` (mode switch `1`/`2`).

Chase-cam: camera position = `toThree(stateAt(craft,t))` + a small trailing
offset opposite `chaseLookDir`; look-at along `chaseLookDir`; the fixed
starfield/grid stays for orientation (motion-sickness mitigation). Key `2` →
chase, `1` → orbit. **Under reduced-motion, `2` is ignored** (stays orbit-cam)
and the live region announces "chase-cam disabled under reduced-motion".

#### Verification

- [ ] Pure unit: `chaseCameraPose(craft,t)` returns position+lookAt plain coords;
  position equals `toThree(stateAt(craft,t))` plus the offset.
- [ ] Browser: `2` enters chase (rides the craft, looks ahead); grid/stars give
  a stable reference; `1` returns to orbit-cam.
- [ ] Reduced-motion emulation: `2` does nothing; announcement fires.

Commit:

```
viz-2b/chase: chase-cam mode; reduced-motion forces orbit-cam (viz phase 2b slice 2)
```

### Task 2.2 — Slice-2 gate

**Files:** none new; type/lint + browser pass.

- [ ] Chase ↔ orbit switching is smooth; chase disabled under reduced-motion.
- [ ] No Slice-1 regression; bundle still zero-cost-until-click.
- [ ] `npm run astro check`; `npm test` green.

Commit:

```
viz-2b/chase: Slice 2 gate — chase-cam with motion-sickness guard (viz phase 2b slice 2)
```

---

## Slice 3 — guided tour

### Task 3.0 — geometry-derived tour keyframes (pure)

**Files:** create `src/lib/three-tour.ts` (`tourKeyframes(cfg)`); test
`src/lib/__tests__/three-tour.test.ts`.

`tourKeyframes(cfg) -> {t, pose, caption, beat}[]` built from the config:
departure (`t=t0`), flyby (`cfg.proximityMinima[0].t`), aphelion (scanned max
radius), return phasing (`~t1`). **Derived, not hardcoded per row.** Captions:
"Earth departure", the `encProvenance` flyby string (pass it in `cfg`), "long
aphelion arc — time visibly slows", "Earth-return phasing".

#### Failing test

```ts
import { describe, it, expect } from "vitest";
import { tourKeyframes } from "../three-tour";

const cfg = {
  t0: 0, t1: 1000,
  craft: { a: 1.6, e: 0.2, i_deg: 0, lan_deg: 0, argp_deg: 0, M0_deg: 0 },
  encounterTimes: [50],
  proximityMinima: [{ body: "M", t: 50, d_au: 0.01 }],
  planets: [{ code: "E", el: { a: 1, e: 0.0167, i_deg: 0, lan_deg: 0, argp_deg: 102.9, M0_deg: 0 } }],
  bodies: ["M"],
} as any;

describe("tourKeyframes (geometry-derived)", () => {
  it("emits the four didactic beats in time order", () => {
    const k = tourKeyframes(cfg);
    expect(k.map((x) => x.beat)).toEqual(["departure", "flyby", "aphelion", "return"]);
    for (let i = 1; i < k.length; i++) expect(k[i].t).toBeGreaterThanOrEqual(k[i - 1].t);
  });
  it("places the flyby beat at the proximity minimum", () => {
    expect(tourKeyframes(cfg).find((x) => x.beat === "flyby")!.t).toBe(50);
  });
  it("places aphelion at max heliocentric radius in [t0,t1]", () => {
    const ap = tourKeyframes(cfg).find((x) => x.beat === "aphelion")!.t;
    expect(ap).toBeGreaterThan(0);
    expect(ap).toBeLessThanOrEqual(1000);
  });
});
```

Run → **red** → impl `tourKeyframes` (reuse `stateAt`/`distance` for the aphelion
scan) → **green**. `astro check`. Commit:

```
viz-2b/tour: geometry-derived tour keyframes (departure/flyby/aphelion/return) (viz phase 2b slice 3)
```

### Task 3.1 — tour playback (keyframe lerp) + live-region narration; reduced-motion = step

**Files:** `three-view.ts` (tour runner), wire key `T`.

`T` starts/stops the tour: lerp camera pose + advance `t` between keyframes,
updating the live region with each beat's caption + live proximity. **Under
reduced-motion the tour does NOT autoplay**: `T` arms it and `[`/`]` walk beat to
beat with snapped poses (design §4.3). The tour is a clock script over the shared
timeline (reproducible, inspectable).

#### Verification

- [ ] Browser: `T` flies through departure → flyby (camera pulls in, proximity
  ~0) → aphelion (pulls back, time slows) → return; live region narrates each
  beat.
- [ ] Reduced-motion: `T` arms; `[`/`]` step beats with snapped camera; no
  autoplay.

Commit:

```
viz-2b/tour: keyframe playback + live-region narration; reduced-motion steps beats (viz phase 2b slice 3)
```

### Task 3.2 — Slice-3 gate + full 2b verification

**Files:** none new; type/lint + full browser regression.

- [ ] All three modes work; tour narrates; reduced-motion path correct across all
  three.
- [ ] Honesty preserved: badges/clock-regime/citation visible in 3D; CR3BP rows
  show no button; (multi-arc 3D remains out of scope — confirm the button is
  gated to single-ellipse and a note covers the rest).
- [ ] Bundle: zero-cost-until-click holds; no 2a regression; other routes
  untouched.
- [ ] `npm run astro check`; `npm test` green.

Commit:

```
viz-2b/tour: Slice 3 gate — guided tour + full 2b a11y/bundle/honesty verification (viz phase 2b slice 3)
```

---

## Self-review (author's pre-flight)

- **Contract honoured.** 2b adds **no** new physics: every position comes from
  `kepler-time.ts` `stateAt`/`samplePath`, every planet from the sourced
  `planet-elements.json` via `orbit.ts`. The 3D scene reads the **existing**
  `clockConfig` JSON (the synchronisation contract) — one clock, two renderers.
- **Zero-cost-when-unused is proven, not asserted** — baseline (1.1) → click-only
  chunk (1.2/1.3) → gate (1.7) with a network-log check. `three-controls` is
  hand-rolled to avoid the examples bundle.
- **Axis convention is one function** (`toThree`, Task 0.1), unit-tested before
  any WebGL, documented against `kepler-time.ts`'s ecliptic z=0 / +z-north frame.
- **A11y is designed, not waved at** — focus management (Esc → SVG), full keyboard
  map, key-help overlay, live-region narration, reduced-motion (no autoplay,
  chase disabled, snaps), dark mode (two material sets); the SVG remains the
  accessible source of truth throughout.
- **Honesty carried in** — model badge + clock-regime + citation on the canvas;
  CR3BP shows no button; gaps never invented; idealized clock never shows a date.
- **Refactor flags (do NOT do in 2b):**
  1. **`velocityAt` in `kepler-time.ts`** — chase-cam uses a finite-difference of
     `stateAt` because the module returns position only (`:109-110`). If chase
     jitter is visible, the honest fix is a closed-form `velocityAt` export
     upstream; flagged, not done here (would touch the shared module other agents
     and the SVG depend on).
  2. **`clockConfig` lacks `encounterTimes`/`proximityMinima`** — added additively
     in Task 1.0b (the 2a island ignores unknown keys); a deeper refactor to make
     `clockConfig` a typed shared interface in `kepler-time.ts` is desirable but
     out of 2b scope.
  3. **Multi-arc + non-keplerian 3D** — Slice 1 gates the button to
     single-ellipse; multi-arc 3D (per-segment arcs + honest gaps) is a future
     slice, flagged, not built.
- **Per-slice shippability** — Slice 1 (orbit-cam) is a complete, useful release
  on its own; 2 and 3 layer on without changing Slice 1's contract.
```

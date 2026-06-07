# Viz Phase 2 — time-true motion, real planet ellipses, a flying 3D camera (cyclers.space)

**Status:** DESIGN DRAFT — no code. Extends the approved Phase-1 viz
(`specs/2026-06-06-web-3d-orbit-visualization-design.md`, task #130; implemented
as task #132 in `cyclers.space/src/lib/orbit.ts` + `src/components/OrbitView.astro`).
**Date:** 2026-06-07
**Task:** #138
**Scope of this doc:** the physics correction the user asked for (time-true
Kepler motion on a shared clock), real planet ellipses (sourced osculating
elements, replacing reference circles), the 2.5D tilt morph, and the opt-in
Three.js flying camera — plus the honesty layer, sourcing requirements, phasing
options, and open questions. Planning artefact only: it prescribes nothing that
ships without a follow-up implementation plan and the user's answers to §8.

---

## 0. Survey — what Phase 1 actually built (cite file:line)

The implemented Phase-1 viz is real and shipped; this design corrects it, it
does not start from scratch.

- **Geometry is exact but TIME-FREE.** `cyclers.space/src/lib/orbit.ts:52-80`
  (`sampleEllipse`) samples the heliocentric ellipse by walking **true anomaly
  ν uniformly** (`nu = (2 * Math.PI * k) / n`, `orbit.ts:68`) and applying the
  exact polar form `r(ν) = p / (1 + e·cos ν)` (`orbit.ts:69`). This draws the
  correct *shape* but carries **no clock**: it has no mean anomaly, no Kepler
  solve, no time. There is no `M → E → ν` anywhere in the file (confirmed: the
  module's only angle input is ν itself).
- **The animation walks path-length, not Kepler-true speed.** The play/scrub
  island in `OrbitView.astro:335-342` (`setPhase`) places every marker at a
  **fractional arc length** via `getPointAtLength((t % 1) * len)`
  (`OrbitView.astro:338`). Arc-length-uniform motion is *not* Kepler motion: it
  is too slow near periapsis and too fast near aphelion (the exact inverse of
  the truth). The auto-play advances all markers by the same `phase`
  (`OrbitView.astro:362`, `~12.5 s per revolution`) — so the spacecraft and
  every planet share a *path-fraction* clock, not a *time* clock. **This is the
  user's observation (1b/1c) confirmed in code:** the planets and spacecraft
  have no physical time relationship, and the spacecraft itself moves at the
  wrong speed along its own ellipse.
- **Planets are circles.** `orbit.ts:33-37` (`PLANETS`) carries `sma_au` only;
  `sampleCircle` (`orbit.ts:114-121`) draws a pure circle at that radius;
  `OrbitView.astro:67` renders them as reference circles. Mars's e=0.0934 is
  not represented — confirming the user's observation (1a). The caption says
  "planet orbit (circle at sma)" (`OrbitView.astro:259`).
- **Encounter markers assume circular planets.** `idealEncounters`
  (`orbit.ts:137-158`) finds where the cycler ellipse crosses the planet's
  *circular* radius `planet.sma_au` (`orbit.ts:148`). On a real planet ellipse
  the crossing radius is no longer a single number, so this logic must change
  (§3).
- **Three.js was reserved, not killed.** The Phase-1 approval
  (`2026-06-06-web-3d-orbit-visualization-design.md:311-338`) explicitly keeps
  Option A "in reserve, page-scoped and lazy-loaded, for a future genuinely-3D
  family," and — binding correction — records that **"no animation" is NOT a
  user requirement**; time-scrubbed motion is on the table, only
  `prefers-reduced-motion` is a hard constraint
  (`...:316-325`). The current site already honours this: `global.css:3-4`
  documents user-controlled reduced-motion-aware motion and `global.css:463-465`
  has the `prefers-reduced-motion: reduce` block. So Phase 2 does **not** need
  to re-litigate animation; it needs to make the animation *physically true*.
- **The data we have for planets** (`cyclers/src/cyclerfinder/core/constants.py`):
  `PlanetData` (`constants.py:120-173`) carries `sma_au`, `ecc` (sourced,
  Standish & Williams Table 1 e_0 column, `constants.py:151-159`), and *defines*
  `inc_deg`/`lan_deg` fields whose **sourced J2000 values are recorded in the
  docstrings/comments** for all eight planets (e.g. Venus inc=3.39467605,
  lan=76.67984255 at `constants.py:217-220`; Mars inc=1.84969142,
  lan=49.55953891 at `constants.py:244-245`) but are left at the coplanar
  default `0.0` in the live circular backend (`constants.py:171-173, 215`).
- **What is MISSING for time-true positions** (verified by grep over
  `constants.py` and `cyclerfinder/ephemeris/`): there is **no
  longitude-of-perihelion ϖ** and **no mean-longitude-at-epoch L0** anywhere.
  Both are in the same already-cited Standish & Williams Table 1 (the `varpi_0`
  / `L_0` columns) and are **required** to know *where each planet is at a given
  date*. Without them we can draw a real ellipse but cannot place the planet on
  it at time t. **These two angles per planet must be sourced (§4).**
- **The spacecraft side has two regimes** (`cyclers.space/src/lib/types.ts`):
  epoch-free idealized rows carry `orbit_elements` with `epoch_iso8601`
  *usually null* (`types.ts:80`) and segments with optional `tof_days`
  (`types.ts:90,107`) and `tof_days_bounds` (`types.ts:118`); real-window rows
  expose real DE440-matched encounter dates through `encounterWindowsFor`
  (`cyclers.space/src/lib/catalogue.ts:44-53`) reading `next_encounters_iso` +
  `vinf_actual_kms`. So the spacecraft has an honest *real* clock only when a
  window is matched; otherwise it has only relative ToFs and its own period.

**Key takeaway:** Phase 1 nailed the *shapes* and the provenance discipline.
Phase 2's job is the *clock* — one shared time axis that drives Kepler-true
planet and spacecraft motion so the cycler's defining property (an encounter is
a coincidence in space **and** time) becomes visible — plus upgrading planets
from circles to sourced ellipses, and finally answering the user's "hard to get
the 3D idea" with a real (opt-in, lazy) camera.

---

## 1. The shared-clock physics layer — `kepler-time.ts`

The single new module Phase 2 needs. It introduces *time* to a stack that today
only knows *shape*. It is pure TS, framework-free, a few hundred lines, shared
by **both** renderers (SVG and Three.js) so the physics has exactly one home.

### 1.1 What it computes

```
M(t)  = M0 + n · (t − t_epoch)          // mean anomaly from the clock
E     : E − e·sin E = M                  // Kepler's equation (Newton solve)
ν(E)  = 2·atan2(√(1+e)·sin(E/2), √(1−e)·cos(E/2))
r(ν)  = a(1 − e²)/(1 + e·cos ν)
x⃗     = R3(−Ω)·R1(−i)·R3(−ω) · [r cosν, r sinν, 0]   // perifocal → ecliptic
```

- **`solveKepler(M, e)`** — Newton–Raphson on `E − e sin E − M`, seeded with
  `E0 = M + e sin M`, ~4 iterations to 1e-12 for e ≤ 0.21 (Mercury is the worst
  case in our table; cycler ellipses go higher but still elliptic). Falls back
  to a bisection guard for robustness. This is the exact step `orbit.ts` lacks.
- **`stateAt(elements, t)`** — returns the ecliptic position (and, cheaply, the
  velocity for the encounter-proximity readout in §6) for any body or trajectory
  leg described by `(a, e, i, Ω, ω, M0, t_epoch, n)`. `n` (mean motion) is
  derived from `a` via Kepler III when not supplied — exactly as
  `constants.py:176-186` (`_mean_motion_deg_day`) does upstream, keeping the two
  stacks consistent.
- **`samplePath(elements, t0, t1, n_samples)`** — replaces the ν-uniform
  `sampleEllipse`: samples *time* uniformly and converts each t to a position.
  The same point set, but now each rendered dot corresponds to an equal *time*
  step, which is what makes "slow at aphelion, fast at periapsis" visible.
  `sampleEllipse` is kept for the static no-JS shape (the closed curve is the
  same set of points regardless of clock), so Phase 1's JS-off fallback is
  untouched.

### 1.2 Per-leg time windows for multi-arc

A multi-arc cycler is not one ellipse on one clock. Each
`trajectory.segments[]` (`types.ts:90,107`) carries a `tof_days` (a **sourced
transit ToF**, when present). The module composes a piecewise clock: leg k
occupies `[Σtof_{<k}, Σtof_{<k} + tof_k]`; within a leg the spacecraft moves on
that leg's `(a,e)` (`types.ts:110`) by the Kepler solve above. Where a leg's
`tof_days` is **null**, the clock has a **gap** — the module returns an explicit
"unknown duration" marker and the renderer draws nothing across it (mirrors the
Phase-1 null-segment discipline at `OrbitView.astro:267-276`). We never invent a
transit time to make the animation continuous.

### 1.3 Two honest spacecraft clocks (the epoch-free problem)

This is the crux of doing time-true motion honestly given §0's two regimes:

- **Real-window rows (anchored clock).** When `encounterWindowsFor` returns
  matched DE440 dates (`catalogue.ts:44-53`), anchor `t_epoch` to a chosen real
  window's first-encounter date. Now the spacecraft *and* the planets (which
  have a real ephemeris clock from §4) run on the **same calendar time**, and
  the encounter markers land where the planet truly is on that date. This is the
  full payoff: real coincidence in space and time.
- **Epoch-free idealized rows (phase clock).** Most rows have
  `epoch_iso8601 = null` (`types.ts:80`) — there is no real t0. Honest
  convention: run an **idealized phase clock** whose unit is the cycler's own
  synodic/repeat period, **explicitly labelled "idealized phasing (no real
  epoch) — encounters shown at the geometry's own period."** Planets are placed
  by their *relative* mean motions (real n ratios from §4) but with an arbitrary
  common t=0, so the *relative* phasing (the thing a cycler is *about*) is true
  even though the absolute date is not. The label must make clear: the dance is
  real, the calendar is a placeholder. We never print a fake date.

**Provenance rule:** the clock regime is surfaced in the caption exactly like
Phase 1 surfaces `model_assumption` (`OrbitView.astro:120-126, 234-240`):
`real-window-anchored (DE440)` vs `idealized phase clock (no epoch)`.

---

## 2. Real planet ellipses (sourced elements)

Replace the reference circles with the planets' true J2000 osculating ellipses.

- **Data:** extend `PLANETS` in `orbit.ts:33-37` from `{sma_au}` to the full
  `(a, e, i, Ω, ϖ, L0)` per body. `a` and `e` already exist upstream
  (`constants.py:195-202` sma; `constants.py:214,231,242,260,…` ecc). `i` and Ω
  exist as sourced values *in the upstream comments* (`constants.py:217-220` etc)
  and must be promoted to live fields. **ϖ and L0 must be added (§4).** From
  these, `ω = ϖ − Ω` and `M0 = L0 − ϖ` — the standard Standish reduction.
- **Provenance label change:** the caption flips from "planet orbit (circle at
  sma)" (`OrbitView.astro:259`) to **"J2000 osculating ellipse (Standish &
  Williams Table 1)."** This is a *strengthening* of provenance — we move from
  an idealization to sourced geometry, with the citation already used elsewhere
  in `orbit.ts:20-26`.
- **Sourced JSON for the site.** Same discipline as the catalogue sync
  (`scripts/sync-catalogue.mjs`, design §0). Add a small, reviewable
  `planet-elements.json` emitted from `constants.py` (single source of truth
  upstream) carrying the 6 elements + the citation string per body, copied into
  `src/data/` at `predev`/`prebuild`. The site never hand-copies numbers; it
  consumes the synced artefact. (Alternative: inline the 8×6 table directly in
  `orbit.ts` with the citation comment, as the current `sma_au` values already
  are at `orbit.ts:33-37`. Recommend the JSON-from-upstream path so the numbers
  trace to `constants.py` and can't drift — but for V/E/M only, inlining is
  defensible and lighter.)

### Encounter markers when planets aren't circles

`idealEncounters` (`orbit.ts:137-158`) currently solves "where does the cycler
ellipse reach the planet's circular radius." With a real planet ellipse the
encounter is no longer a fixed radius. Two honest replacements, in order of
fidelity:

1. **Time-true (preferred when there's a clock).** The encounter marker is
   simply *the planet's position at the encounter time* (real date for
   anchored rows; phase-clock time for idealized rows). No radius-crossing math
   — the marker is wherever the planet actually is, and the
   encounter-proximity indicator (§6) confirms the spacecraft is there too. This
   is strictly more honest than radius-crossing and falls straight out of §1.
2. **Geometry-only fallback (no clock / static figure).** When we only have the
   shape, mark the cycler ellipse's *minimum-distance point* to the planet's
   ellipse (or the radius band [peri,apo] crossing), labelled "idealized
   crossing geometry." This generalizes the current single-radius logic to the
   eccentric-planet case and stays explicitly idealized.

The radial-crossing markers therefore **move** relative to Phase 1 (the
user's observation 1a): on Mars's real ±9% ellipse the crossing is at a
different place than on the circle, which is precisely the misplacement the user
flagged.

---

## 3. The tilt morph (2.5D) — keep, and make it the camera's understudy

Phase 1 ships a *separate* edge-on companion panel (`OrbitView.astro:244-255`,
`sampleEdgeOn` `orbit.ts:89-111`) with an explicit z-exaggeration label. Phase 2
adds a **continuous top-down → edge-on projection slider** on the *same* SVG:
one parameter `tilt ∈ [0°, 90°]` interpolating the projection matrix (the SVG
`y` becomes `y·cos(tilt) − z·sin(tilt)` with the existing exaggeration applied
to z). At 0° it is exactly today's top-down view; at 90° it is the edge-on view;
in between it reads as a tilting plane.

**Honest assessment (the user asked for this explicitly):** the tilt morph is a
**complement, not a replacement** for the camera, and it does **not** fully
deliver "the 3D idea."

- It *does* cheaply communicate "this is a tilted plane in space," it stays in
  pure SVG (zero bundle), it keeps every leg/marker a focusable DOM node, and it
  inverts with dark mode and respects reduced-motion for free — all the
  Phase-1 virtues. For the near-coplanar V/E/M majority (i ≤ 3.4°,
  Phase-1 design §1) it is genuinely *enough* to show that the orbits are nearly
  flat.
- It does **not** give parallax, depth ordering, or the felt sense of *being in
  the system* that the user means by "the 3D idea" and "a flying camera." A
  single-axis tilt of line art is still line art; it cannot convey the long
  aphelion arc swinging out past Mars the way a moving perspective camera can.

**Recommendation:** ship the tilt morph (it's small and improves the default,
no-WebGL experience), but treat it as the **understudy** that every user gets,
with the Three.js camera (§4) as the opt-in upgrade for users who click "View in
3D." Tilt does *not* let us defer the camera — the user explicitly wants the
camera — but it does mean the 99%-case page never needs WebGL.

---

## 4. The flying camera (Three.js, opt-in, lazy)

This re-opens the Option-A reserve from the Phase-1 approval
(`2026-06-06-…:331-334`), scoped exactly as that approval permits: **page-scoped,
lazy-loaded, zero bytes on every other route.**

### 4.1 Bundle discipline (non-negotiable)

- The site today has **zero runtime framework / WebGL** (deps:
  `astro`, `typescript`, `@astrojs/check` only —
  `cyclers.space/package.json` dependencies block). Three.js must not change
  that for any page that isn't actively showing 3D.
- **Route-level lazy `import()` on explicit user action.** A "View in 3D"
  button on `cycler/[id].astro` triggers `const THREE = await import("three")`
  (+ a hand-rolled minimal controls module, not the full examples bundle). Zero
  bytes ship on home, catalogue, launch-windows, about, or even the detail page
  until the button is pressed. This is stricter than "lazy on the route"; it's
  "lazy on intent," matching the Phase-1 enhancement-only model
  (`OrbitView.astro:291` ships the controls `hidden`, JS reveals them).
- Tree-shaken core three is ~150 KB min (Phase-1 design Option A). Acceptable
  *only* behind the click. Document the byte cost in the button's title.

### 4.2 Camera modes

Three modes, sharing `kepler-time.ts` (§1) so the 3D scene and the SVG run the
**same clock** — switching to 3D continues the same instant in time:

- **Chase-cam** — rides the spacecraft, looking ahead along velocity. Best for
  "what does the journey feel like." Risk: disorienting; mitigate with a faint
  fixed star-field / ecliptic grid for reference.
- **Orbit-cam** — slow user-or-auto orbit around the whole system with the full
  trajectory inked in (the user's "slow orbit around the system with the
  trajectory inked in"). Best for comprehension of the *shape* in 3D. This is
  the safest default.
- **Guided tour** — a **scripted camera path hitting the didactic moments**:
  (1) launch / Earth departure, (2) the Mars flyby (camera pulls in to the
  encounter, proximity indicator §6 at ~0), (3) the long aphelion arc (camera
  pulls back, time visibly slows — the Kepler payoff), (4) the Earth-return
  phasing (camera frames Earth and spacecraft converging). Each beat is a
  keyframed (time, camera-pose, caption) tuple; the tour is just a clock script
  over the shared timeline, so it is reproducible and inspectable.

**Recommended default path:** **Orbit-cam, paused, trajectory fully inked,
time at the first encounter** — a legible static 3D figure the moment 3D opens,
with Play and a "Guided tour" button as the two next actions. Chase-cam is the
power-user/"feel it" mode, not the default (too easy to get lost on first open).

### 4.3 Accessibility in a WebGL context (designed, not waved at)

The canvas is opaque to AT — this is the real cost of §4 and must be paid:

- **The 3D view is strictly additive.** The accessible source of truth remains
  the Phase-1 SVG (focusable DOM nodes, `<title>`/`aria-label`,
  `OrbitView.astro:197-205`) and the dl tables. 3D is never the *only* way to
  get any information. If WebGL is unavailable/declined, the SVG is the page.
- **Keyboard:** the canvas gets `tabindex="0"` and a documented key map —
  arrows orbit, `+`/`-` dolly, `[`/`]` step time, `Space` play/pause, `T`
  guided tour, `Esc` returns focus to the SVG. A visible focus ring on the
  canvas and an on-canvas key-help overlay (toggled with `?`). No mouse-only
  affordance.
- **A live region** (`aria-live="polite"`, off-canvas DOM) announces the camera
  mode and the current didactic beat ("Approaching Mars flyby; spacecraft–Mars
  distance 0.02 AU") so a screen-reader user gets the *narrative* the sighted
  user gets from the tour, driven by the same `kepler-time.ts` state.
- **`prefers-reduced-motion`:** **the camera does not auto-orbit and the tour
  does not auto-advance** under reduced-motion (consistent with
  `global.css:463-465` and the Phase-1 hard rule). The user steps time/camera
  explicitly; all transitions snap. The button label notes "(reduced-motion:
  manual step)."
- **Dark mode:** the scene reads `prefers-color-scheme` (and the CSS variables
  Phase 1 uses) at init to pick clear-color, line, and label materials — a
  light scene on light pages, dark on dark — mirroring `global.css:23-35`. Two
  small material sets, switched once at init.

### 4.4 Sharing the physics

Both renderers import `stateAt`/`samplePath` from `kepler-time.ts`. The SVG
calls them at build-time (static curve) and in its scrub island; the Three.js
scene calls them per frame for the moving bodies and once for the inked
trajectory. There is **one** Kepler solver and **one** clock convention; the 3D
view can never disagree with the 2D view about where anything is.

---

## 5. Honesty layer

Provenance is the brand (`[id].astro:135-144`, Phase-1 design §5). Every mode
is labelled by **two independent axes**, never blended:

- **Clock axis:** `real-window-anchored (DE440)` vs `idealized phase clock (no
  epoch)` (§1.3). An idealized clock never displays a fabricated calendar date.
- **Geometry axis:** `planets: J2000 osculating ellipse (Standish Table 1)`
  (the Phase-2 default) vs the row's own `model_assumption`
  (`circular-coplanar` / `analytic-ephemeris` / `cr3bp`, `types.ts:11`) for the
  *spacecraft* curve. The two can differ (real planet ellipses around an
  idealized-coplanar cycler) and the caption says so.
- **Encounter markers** keep the Phase-1 real-vs-idealized split
  (`OrbitView.astro:79-81`), upgraded per §2 (time-true marker when a clock
  exists; idealized minimum-distance otherwise).
- **Gaps stay gaps** — null-`tof_days` legs and CR3BP rows render as explicit
  "not shown" notes (Phase-1 `OrbitView.astro:267-276,144-156`), in 2D and 3D
  alike. The camera never flies across an unknown leg.

### The pedagogical payoff — encounter-proximity indicator

The reason time-true motion matters: surface a **spacecraft↔planet distance
readout** that **shrinks toward ~0 at true encounters**. This is the single most
important new affordance — it makes the cycler's defining property (coincidence
in space *and* time) *visible and measurable*:

- A small panel / on-figure label showing live `|r_sc − r_planet|` in AU (and
  km at close approach) for the currently-relevant body, computed from
  `kepler-time.ts` velocities/positions.
- A sparkline of that distance over the timeline, with minima marked — the dips
  ARE the encounters. On a real-window-anchored row the minima sit on the DE440
  dates; on an idealized phase clock they sit at the geometry's own period
  (labelled as such).
- This is the answer to the user's "a cycler IS a phasing phenomenon": the
  indicator is the proof that walking markers at uniform path-speed (Phase 1)
  was hiding the whole point — with the true clock, the distance actually goes
  to zero *at the same time* the spacecraft reaches the planet's orbit.

---

## 6. Site-principles compliance (gating checklist)

Unchanged from Phase-1 §4 except where Phase 2 adds surface:

- **Motion is user-initiated; reduced-motion respected** — already in
  `global.css:3-4,463-465`; extend the reduced-motion block to the 3D canvas
  (no auto-orbit, no auto-tour) and the tilt morph (snap, no transition).
- **No colour-only meaning** — carries over (dash + glyph + text,
  `OrbitView.astro:259-264`); in 3D, line style + label + marker shape, never
  hue alone.
- **Keyboard** — SVG already focusable; the camera gets the full key map (§4.3).
- **Dark mode** — SVG via `currentColor`/vars (Phase 1); 3D via dual material
  sets (§4.3).
- **Bundle** — SVG + tilt + `kepler-time.ts` are a few KB, ship on the detail
  page as today; **Three.js ships zero bytes until the "View in 3D" click**
  (§4.1). Home/catalogue/launch/about: unchanged, zero new bytes.
- **JS-off safe** — the build-time static SVG (Phase 1) is the floor; tilt and
  3D are enhancements that simply don't appear without JS.

---

## 7. Phasing / scoping options + recommendation

**Option 1 — Physics-first, camera later (RECOMMENDED).**
Phase 2a: `kepler-time.ts` + real planet ellipses + time-true SVG animation +
encounter-proximity indicator + tilt morph. Phase 2b (separate plan, after 2a
lands): the Three.js flying camera. *Rationale:* the user's correctness
observations (1a/1b/1c) are the substance and are deliverable with **zero new
dependencies**; they also de-risk the camera (which just reuses
`kepler-time.ts`). Ship the truth first, the spectacle second.

**Option 2 — Camera-first.**
Build the Three.js scene now because "hard to get the 3D idea" is the loudest
pain. *Rejected:* a flying camera over *wrong* motion (path-speed, circular
planets) would animate a falsehood beautifully — exactly the provenance sin the
brand forbids. The clock must be true before we make it cinematic.

**Option 3 — Everything in one Phase 2.**
One big drop: physics + ellipses + tilt + camera + proximity. *Rejected for
scope/risk:* couples a zero-dependency correctness win to a 150 KB WebGL
feature with a hard a11y surface; if the camera slips, the correctness fix
slips with it. Split per Option 1.

**Recommendation: Option 1.** Ship Phase 2a (time-true physics, real ellipses,
proximity indicator, tilt) as a dependency-free correctness release; follow with
Phase 2b (the opt-in lazy camera) as its own plan once 2a is validated.

---

## 8. Open questions for the user (verbatim)

1. Do you want me to source **ϖ (longitude of perihelion) and L0 (mean
   longitude at epoch)** for all eight planets from Standish & Williams Table 1
   and add them to `constants.py` (upstream single source) + a synced
   `planet-elements.json` for the site — or, for now, only V/E/M, inlined in
   `orbit.ts` with the citation, deferring the other five until a non-V/E/M
   family needs them?
2. For epoch-free rows, do you accept the **"idealized phase clock (no real
   epoch)"** convention — real relative mean-motion ratios, arbitrary common
   t=0, encounters at the geometry's own period, never a fabricated calendar
   date — as the honest default, with real-window rows anchored to a chosen
   DE440 date?
3. Phasing: **Option 1 (physics-first 2a, camera-later 2b)** as recommended, or
   do you want the camera in the same release despite the bundle/a11y cost?
4. Camera default: **Orbit-cam, paused, trajectory inked, time at first
   encounter** as the opening view (chase-cam as a power-user mode), or do you
   want chase-cam or the guided tour to be what loads first?
5. For real-window rows, which **encounter date anchors t=0** — the earliest
   matched window, the lowest-V∞ window, or a user-selectable one?
6. Is the **encounter-proximity indicator** (live spacecraft↔planet distance +
   minima sparkline) worth building in 2a as the headline pedagogical feature,
   or is plain time-true motion enough for the first cut?
7. Planet ellipse sourcing path: **synced JSON emitted from `constants.py`**
   (numbers trace upstream, can't drift) or **inlined in `orbit.ts`** (lighter,
   like the current `sma_au` table) — your call on the trade between traceability
   and simplicity?

---

## Approval (2026-06-07)

User-approved with all recommendations accepted: (Q1) source ϖ and L0 for
ALL EIGHT planets from Standish & Williams Table 1 into constants.py
(upstream single source, same citation discipline) + the synced
planet-elements.json for the site; (Q2) the "idealized phase clock (no real
epoch)" convention adopted — never a fabricated calendar date; real-window
rows anchor to a DE440 date; (Q3) Option 1 physics-first — 2a now, the
camera as its own 2b plan; (Q4) camera default = orbit-cam, paused,
trajectory inked, time at first encounter; (Q5) the lowest-mismatch window
anchors t=0 (most physically meaningful), user-selectable later; (Q6) the
encounter-proximity indicator IS built in 2a — the pedagogical headline;
(Q7) synced JSON emitted from constants.py — traceability over lightness.

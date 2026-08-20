# #857 — representing `data/manifold_connections.yaml` on cyclers.space: design recommendation

**Date**: 2026-08-21
**Scope**: design proposal only. No file in either repo (`Bwooce/cyclers` or `Bwooce/cyclers.space`)
is modified by this note. This follows the `#838` -> `#856` process precedent (design note first,
separate user-approved implementation task after). Sources read in full for this design:
`data/manifold_connection.schema.json` (v1.0) + `data/manifold_connections.yaml` (2 entries) +
`docs/notes/2026-08-21-838-connection-schema-design-recommendation.md` in the cyclerfinder repo;
and in the cyclers.space repo: `README.md`, `package.json`, `.gitignore`,
`scripts/sync-catalogue.mjs`, `.github/workflows/refresh-windows.yml`,
`src/pages/cycler/[id].astro`, `src/pages/catalogue/index.astro`,
`src/components/CatalogueTable.astro`, `src/pages/index.astro`, `src/lib/errata.ts`, and the
loader/sanitizer surface of `src/lib/catalogue.ts`.

---

## 1. The question, restated

The cyclerfinder repo just gained a dedicated registry, `data/manifold_connections.yaml`
(schema v1.0, `#838` design / `#856` implementation): verified extrinsic manifold connections
(heteroclinic/homoclinic transfers) between two orbits, where at least one endpoint is a
catalogued row and the other may be an inline `uncatalogued` descriptor (the half-catalogued base
case). The public website renders `data/catalogue.yaml` but has no representation of this
registry — nor of its sibling `data/cycler_networks.yaml` (`#570`), which ships genuinely empty
and so has never needed one. How should the connections registry appear on cyclers.space:
what sync mechanism, what per-row display, what catalogue-table indicator, whether a dedicated
page, and how to surface the `evidence_class` distinction without visually implying a validation
promotion the schema explicitly forbids?

## 2. Recommendation (read this first)

**Extend the existing `prebuild` fetch script (`scripts/sync-catalogue.mjs`) with a fourth
soft-fail block for `manifold_connections.yaml` (the committed `errata.yaml` pattern, not the
gitignored `catalogue.yaml` pattern), add a small `src/lib/connections.ts` loader mirroring
`errata.ts`, and render a new "Manifold connections" section on the two-to-three touched
`/cycler/<id>/` detail pages — connections render inline there, with no per-connection page.**
Defer both the `/connections/` index page and any catalogue-table indicator until the registry
grows past its current 2 near-identical entries (`#840`/`#854`/`#839` are queued, so likely
soon); the loader should be written list-first so the index page is a cheap follow-on.
`evidence_class` IS surfaced, as a neutral-styled label with plain-language wording, under an
explicit "does not change either orbit's validation level" caption.

## 3. Reasoning, point by point

### 3.1 Sync mechanism — extend `sync-catalogue.mjs`, one block, nothing new invented

The premise that the site's data arrives by a one-off manual copy is not what the code shows.
The actual mechanism, confirmed in `cyclers.space/package.json` and `scripts/sync-catalogue.mjs`:

- `predev`/`prebuild` npm hooks run `node scripts/sync-catalogue.mjs` before **every** dev
  session and every build. The script fetches from
  `https://raw.githubusercontent.com/Bwooce/cyclers/main/...` with per-file env-var overrides.
- It already handles **four** files under **two** policies:
  - `src/data/catalogue.yaml` — **gitignored** (`.gitignore`: "canonical catalogue is fetched at
    build time ... never committed"), **hard-fail** on non-200 ("Refusing to build with
    stale/missing catalogue"), with a content sniff (`body.includes("- id:")`).
  - `src/data/planet-elements.json`, `src/data/errata.yaml`, and
    `public/data/sampled/<id>.json` — **committed** copies refreshed by the sync, **soft-fail**
    (keep the committed copy on fetch failure; hard-fail only if the committed copy is also
    missing). Rationale stated in-file: small files, offline-reproducible build.
- Separately, `.github/workflows/refresh-windows.yml` (weekly cron) curls `catalogue.yaml` and
  recomputes `windows.json`. The README's "syncs on a schedule" claim refers to this; whether
  the cron is currently firing is a CI question (GitHub Actions minutes are exhausted for
  Aug 2026), but it is orthogonal — the prebuild fetch runs on every build regardless, so the
  2026-07-28 mtime on `src/data/catalogue.yaml` just dates the last local build, not a manual
  copy step.

So the answer is the boring one: **add a fourth fetch block** to `sync-catalogue.mjs` for
`data/manifold_connections.yaml` -> `src/data/manifold_connections.yaml`, using the
**committed + soft-fail** policy (`errata.yaml` is the exact template, ~15 lines, same
`- id:` sniff). Reasons for committed-soft-fail over gitignored-hard-fail: the file is ~5 KB
(errata-sized, not catalogue-sized), and a fetch outage should not brick the whole site build
over a supplementary registry.

Two deliberate non-syncs:

- **The schema file does not need syncing.** The site performs no JSON-Schema validation of any
  synced file today (there is no schema copy of `catalogue.schema.json` or `errata.schema.json`
  in the site repo either); structural validation is upstream's job
  (`src/cyclerfinder/data/validate_connections.py` exists and gates the registry in
  cyclerfinder's own test suite). The site's `types.ts`-style interface is its own contract, as
  it is for the catalogue.
- **`cycler_networks.yaml` stays unsynced.** It is empty; when it first populates, it can reuse
  this exact block pattern (and possibly share a page — see 3.4).

One ordering constraint worth stating: the raw-URL fetch only works once `#856`'s registry is
pushed to `main` of the public `Bwooce/cyclers` repo. It currently sits in the working tree /
recent commits; the site work should land after that push.

### 3.2 Per-row display — yes: a "Manifold connections" section on `/cycler/<id>/`, inline, no per-connection page

**Show it.** The detail page is exactly where a reader of `vaquero-31-c254-em-cycler-2013`
should learn "a verified ballistic transfer arrives here from a 2:1 family member at C=2.54."
The page template already has the precedents needed:

- **Lookup pattern**: `errataForRow(entry.id)` in `lib/errata.ts` filters a small synced YAML
  registry by row id at render time; `connectionsForRow(entry.id)` is its twin (filter entries
  where any endpoint's `row_ref === id`).
- **Placement**: as its own `<h2>` section, after the "CR3BP orbit identity" block (all touched
  rows are `cycler_class: non-keplerian`, so that section is always present for them) and near
  `<OrbitView>` — i.e. with the row's dynamical identity, not with citations. Section renders
  only when `connectionsForRow` is non-empty, the same conditional-section style every other
  block on the page uses.
- **Cross-row links**: the `inserts_into` rendering (page lines ~219–231) is the exact template
  for a `row_ref` endpoint: resolve via `getEntryById`, link `/cycler/<id>/` when found, fall
  back to `<code>{id}</code>` + a muted "(unresolved)" note when not — never a broken link,
  never a build failure on skew between the (hard-fetched) catalogue and the (soft-fail,
  possibly older committed) connections file.

**Both endpoint shapes, concretely** (the ordered pair is directed: index 0 = unstable-manifold
origin, index 1 = stable-manifold destination):

- `row_ref` endpoint, and it is the row being displayed: render as **"this orbit"** (bold, no
  self-link), with its role stated — "departs along its unstable manifold" / "arrives on its
  stable manifold". The two live entries have the catalogued row on *opposite* sides
  (destination at C=2.54, origin at C=2.66), so the role wording must come from the endpoint
  index, never be hardcoded.
- `row_ref` endpoint, other row: name-link via `getEntryById`, unresolved fallback as above.
  (No both-catalogued entry exists today, but `#839`/`#840` could produce one; the rendering
  handles it for free.)
- `uncatalogued` endpoint: render the `family` label as plain text with a muted qualifier
  "(not a catalogued orbit)", plus its key identifying numbers (C, nondimensional period) in a
  muted sub-line. Nothing to link to — and no pretend-link. The full state vector /
  `derivation` text can sit in a `title` attribute or a `<details>` fold rather than the main
  flow.

**No per-connection page.** A connection entry is one screen-paragraph of content
(direction, kind, model, C, ΔV, evidence label, caveat); a dedicated
`/connection/<id>/` route would be a page for two paragraphs, and would create a second URL
namespace to maintain for a registry of 2. Inline rendering on each touched row's page (yes,
the same connection renders on up to two pages — it is tiny and generated from one source)
plus, later, a `/connections/` index (3.4) covers every navigation need. Each rendered
connection should carry an `id=` anchor (`#conn-em-vaquero-hetero-...`) so a future index page
can deep-link to it.

**What to show per connection** (and what not to): kind + "ballistic" when `dv_kms === 0`
(that is the headline claim — a free transfer), system/model, Jacobi constant, the two
endpoints with roles, the evidence label (3.5), and `round_trip_note` **whenever present** — it
is an honesty field ("reverse direction NOT demonstrated at this C") and hiding it would
overclaim. The raw `evidence` battery (`full_state_gap`, `ghost_distance_*`, `radau_gap`, ...)
and the `connection` geometry block stay out of the default view — a `<details>`
"verification numbers" fold is the right home if shown at all; the site's convention for deep
provenance is "click through to the source", and `provenance.data`/`provenance.notes` paths can
render as links into the GitHub tree the way the catalogue table already links sources.

**Sanitization is mandatory, and the registry is a stress case for it.** `lib/catalogue.ts`'s
`sanitizeCatalogueText` exists precisely because upstream prose is written against the internal
`#NNN` tracker, and the connections registry's free text is *dense* with it
(`identity_evidence`: "#828 independent re-run: ...", `evidence_class`: "... #822 verification
battery ... #828 ...", `round_trip_note`: "... (#840, registered, ...)"). The loader must apply
the sanitizer to every free-text field (`identity_evidence`, `derivation`, `model_note`,
`evidence_class`, `round_trip_note`, `family`), exactly as `errata.ts` does per-field at load
time. The implementer should eyeball the sanitized output of both live entries — leading-"#828:"
sentence-openers are a shape the existing regexes may not leave grammatical.

### 3.3 Catalogue-table indicator — no new column; defer even the badge

The table (`CatalogueTable.astro`) already renders **13–14 columns** (Name, Class, Bodies,
Struct, Period, Identity, V∞, ΔV band, Returns, conditional Validity, Data, Defined,
Validation, Source) inside an `overflow-x:auto` wrapper, plus a five-control filter bar. A
"Connections" column would be an em-dash for ~397 of ~400 rows. Clear no.

The cheaper option is an inline badge in the Name cell — the `discovery-badge`
("Discovered here") precedent shows the mechanism costs almost nothing. But I recommend
**deferring it too**, for two reasons: (a) at 2 connections touching 3 rows, a table-level
affordance optimizes discovery of something a reader can't yet do anything with (there is no
index page to pivot to); and (b) every existing badge in that cell/column region encodes a
*status* (`discovery-badge`, `our-status-badge`, `vlevel`) — adding a connection marker into
that visual family is exactly the "looks like a tier" risk 3.5 warns about. Revisit alongside
the `/connections/` page: when that ships, a small neutral glyph in the Name cell linking to
`/connections/#conn-<id>` becomes genuinely useful.

### 3.4 A dedicated `/connections/` page — defer, but build the loader list-first

Cost side: a new page + probably no new component (a static table, no filter island — the
interactivity in `CatalogueTable.astro` is ~200 lines of script the connections list will not
need for years), plus a nav/README entry. Small, but not free. Benefit side today: an index of
**two rows that differ only in Jacobi constant and direction**. That page would be padding.

But the deferral should be short-fused and cheap to end, because the growth is already queued
in the upstream ledger: `#840` (reverse-direction band-edge demos -> 2 more directed entries),
`#854` (Kumar digit-grade entry -> first second-`evidence_class` entry), `#839` (possible
C=3.13 connection touching a third row). Concretely:

- Write `src/lib/connections.ts` with `loadConnections(): ManifoldConnection[]` as the primary
  accessor and `connectionsForRow()` derived from it — so `/connections/` is a one-file
  follow-on, not a refactor.
- Trigger for building the page: any of — registry reaches ~5 entries; first `DIGIT-GRADE`
  entry lands (the page then has a story to tell, not just rows); or a second system/family
  pair appears.
- Related future need, noted not designed: `data/cycler_networks.yaml` (`#570`) is the same
  species (cross-row relation registry, currently empty). If both populate, a single
  "Relations" or "Connections & networks" page may serve better than two stub pages — one more
  reason not to lock in a `/connections/`-shaped URL this week.

### 3.5 Surfacing `evidence_class`, without faking a validation tier

**Surface it.** This site already treats its public audience as caring about epistemic
provenance — that is its differentiator: V0–V5 badges with an explainer page, `our_status`
badges, an `/errata/` page with confidence taxonomies, and the standing "`—` means the value is
not in the accessible source and is deliberately left blank rather than guessed" footer. The
self-consistency-vs-digit-grade distinction is the same genre of honesty, and the pending
`#854` entry is the whole point: two connections touching the same node with genuinely
different evidentiary weight. Hiding it would make the site *less* honest than its own data.

Mechanics: `evidence_class` is schema-level free text, so the loader classifies by prefix —
`/^SELF-CONSISTENCY/i` -> kind `self-consistency`, `/^DIGIT-GRADE/i` -> kind `digit-grade`,
else `other` — and renders a short plain-language label with the full sanitized text as the
`title`/expandable detail:

- *self-consistency*: "verified by this project's own numerical battery (no published transfer
  state exists to compare against)"
- *digit-grade*: "reproduces a state published in the literature"
- *other*: show the sanitized text itself, so a future third class degrades to visible rather
  than to mislabeled.

**Anti-implication rules** (the schema's "provenance/audit only — not a promotion gate" clause,
translated to CSS):

1. The connection section and evidence label must **not** reuse `.vlevel`, `.our-status-badge`,
   `.dv-band-badge`, or `.orbit-class-badge` classes or their palettes. Neutral outline chip,
   its own class (`.conn-evidence`), no green/tier coloring.
2. The section opens with a one-line muted caption, stating the rule in reader language:
   *"Recorded transport evidence between orbits. A connection does not change either orbit's
   validation level — see the endpoints' own V-levels."* (Optionally anchor an explainer
   paragraph on `/about/` next to `#validation-levels`, the page's established pattern for
   every other badge.)
3. Nothing about a connection appears adjacent to the row's `validation: V1` line in the page
   header — it stays in its own section, below the identity blocks.
4. `round_trip_note` renders whenever present (see 3.2): the registry records *directed*
   demonstrations, and the site must not let "connection" read as "round trip".

## 4. Concrete sketch

**`scripts/sync-catalogue.mjs`** — append a block (verbatim errata pattern):

```js
// --- manifold_connections.yaml (committed; soft-fail on a stale remote) ------
const CONN_URL = process.env.CONNECTIONS_URL ??
  "https://raw.githubusercontent.com/Bwooce/cyclers/main/data/manifold_connections.yaml";
const CONN_OUT = "src/data/manifold_connections.yaml";
// fetch; sniff body.includes("- id:"); write on success;
// on failure keep committed copy, hard-fail only if committed copy missing.
```

**`src/lib/connections.ts`** (mirrors `errata.ts`):

```ts
import yaml from "js-yaml";
import rawYaml from "../data/manifold_connections.yaml?raw";
import { sanitizeCatalogueText } from "./catalogue";

export type ConnectionEndpoint =
  | { row_ref: string; identity_evidence: string; model_note?: string | null }
  | { uncatalogued: { family: string; jacobi_constant?: number | null;
      period_nd?: number | null; derivation: string; /* ... */ } };

export interface ManifoldConnection {
  id: string;
  kind: "heteroclinic" | "homoclinic";
  model: { type: string; system: string; mass_ratio: number };
  jacobi_constant?: number | null;
  endpoints: [ConnectionEndpoint, ConnectionEndpoint]; // [0]=Wu origin, [1]=Ws destination
  evidence_class: string;
  round_trip_note?: string | null;
  dv_kms?: number | null;
  provenance: { task_refs: string[]; data: string; module: string; notes?: string[] | null };
}

export function loadConnections(): ManifoldConnection[] { /* parse, sanitize free text, cache */ }
export function connectionsForRow(rowId: string): ManifoldConnection[] {
  return loadConnections().filter((c) =>
    c.endpoints.some((e) => "row_ref" in e && e.row_ref === rowId));
}
export type EvidenceKind = "self-consistency" | "digit-grade" | "other";
export function evidenceKind(c: ManifoldConnection): EvidenceKind { /* prefix match */ }
```

**`src/pages/cycler/[id].astro`** — new section after the CR3BP-identity block:

```astro
{rowConnections.length > 0 && (
  <>
    <h2>Manifold connections</h2>
    <p class="muted">
      Recorded transport evidence between orbits (verified heteroclinic/homoclinic
      transfers). A connection does not change either orbit's validation level —
      see the endpoints' own V-levels. <a href="/about/#manifold-connections">What
      this means →</a>
    </p>
    <ul class="bare">
      {rowConnections.map((c) => {
        const [origin, dest] = c.endpoints;
        return (
          <li id={`conn-${c.id}`}>
            {/* headline: kind + ballistic + energy */}
            <strong>{c.kind}</strong>
            {c.dv_kms === 0 && <span> · ballistic (0 km/s)</span>}
            {c.jacobi_constant != null && <span> · C = {c.jacobi_constant}</span>}
            <span class="muted"> · {c.model.system}</span>
            {/* endpoints, role-labelled from position, three shapes each:
                this-row (bold "this orbit"), other row_ref (link or unresolved
                <code> fallback), uncatalogued (family text + "(not a catalogued
                orbit)" + muted C/period sub-line) */}
            <div>{renderEndpoint(origin, "departs along its unstable manifold", entry.id)}
              {" → "}
              {renderEndpoint(dest, "arrives on its stable manifold", entry.id)}</div>
            <div><span class="conn-evidence" title={c.evidence_class}>
              {EVIDENCE_LABEL[evidenceKind(c)]}
            </span></div>
            {c.round_trip_note && <div class="muted">{c.round_trip_note}</div>}
            <details><summary class="muted">verification numbers</summary>
              {/* evidence battery table + provenance links, optional */}
            </details>
          </li>
        );
      })}
    </ul>
  </>
)}
```

`.conn-evidence` styling: neutral outline chip, no tier palette. Deferred follow-ons, in order:
`/about/#manifold-connections` explainer anchor; `/connections/` index page (trigger conditions
in 3.4); Name-cell glyph in `CatalogueTable.astro` linking to it.

## 5. What I am NOT certain about — owner sanity checks before approving

1. **Sync timing/skew.** Committed-soft-fail means the connections copy can lag the hard-fetched
   catalogue within one build if the raw fetch fails; the unresolved-`row_ref` fallback (3.2)
   makes that safe, but the owner may prefer hard-fail-both for strict consistency. Also confirm
   the intended landing order: `#856`'s registry must be pushed to public `main` before the
   site's prebuild fetch can see it, and whether the weekly cron is actually firing (CI minutes
   exhausted this month) is worth a one-time check — the prebuild path works regardless.
2. **Sanitizer coverage.** The registry's free text opens sentences with task refs ("#828
   independent re-run: ..."), a shape `sanitizeCatalogueText`'s current rules may leave
   ungrammatical after stripping. Eyeball both live entries' sanitized output; a small rule
   addition may be needed. Same question for `provenance.data`/`module` repo paths — link to
   GitHub, render as plain code, or omit from the public view entirely.
3. **How much raw evidence to expose.** I put the numeric battery behind `<details>`; the owner
   may prefer omitting it entirely (click through to GitHub for numbers) or promoting one
   headline number (e.g. `full_state_gap`). Any shown number needs a caption a non-specialist
   can parse.
4. **Defer-the-page judgment.** If the owner expects `#840`/`#854` to land within weeks anyway,
   building `/connections/` now (one static page over the already-needed loader) is defensible;
   my deferral leans on the 2-near-identical-entries content problem, not on cost.
5. **Naming/anchors.** "Manifold connections" as the public-facing section title assumes the
   audience tolerates the term (the site already says "heteroclinic" nowhere today; `/about/`
   would need to carry the definition). An alternative reader-first title: "Verified transfers
   to/from other orbits".
6. **`cycler_networks.yaml`.** Explicitly out of scope here, but the shared-page question in 3.4
   (one relations page vs two) is a fork the owner may want to pre-decide before any
   `/connections/` URL ships publicly.

## 6. Summary

| Question | Recommendation | Confidence |
|---|---|---|
| Sync | One new soft-fail committed block in `scripts/sync-catalogue.mjs` (the `errata.yaml` template); no schema sync; no new mechanism — the "manual copy" premise is outdated, prebuild fetch already exists | High |
| Per-row display | New conditional "Manifold connections" section on `/cycler/<id>/`, inline, no per-connection page; endpoint renderer handles this-row / other-`row_ref` (link or unresolved fallback, per the `inserts_into` precedent) / `uncatalogued` (labelled text, no fake link); `round_trip_note` always shown | High |
| Table indicator | No column (13–14 columns already, ~397/400 rows would be em-dash); defer even the Name-cell badge until `/connections/` exists to link to | High on no-column; medium on deferring the badge |
| Dedicated page | Defer; build `lib/connections.ts` list-first so the page is a one-file follow-on; triggers: ~5 entries, first digit-grade entry, or second system; consider a combined relations page with the (empty) networks registry | Medium-high |
| `evidence_class` / honesty | Surface it as a neutral outline chip with plain-language labels (prefix-classified, raw text as fallback/detail); never reuse tier/status badge styling; explicit "does not change validation level" caption; sanitize all `#NNN` task refs | High |

No file in either repo has been modified by this task. Implementation, if approved, is a
separate task in the cyclers.space repo (sync block + loader + detail-page section), landing
after `#856`'s registry is pushed to the public repo's `main`.

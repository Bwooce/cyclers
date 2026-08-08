# Digest: Vaquero (2013) Ch. 4.4 — "Resonance Transition in the Earth-Moon System"

**Task:** `#787`, registered 2026-08-08 during `#786`'s own mandatory literature-novelty gate (not
dispatched at registration time). Vaquero's 2013 Purdue dissertation is already fully acquired in
this project's corpus (`#765`'s own primary source — full citation, filing path, and md5 are in
`docs/notes/CORPUS_INDEX.md`, not repeated here); `#765`/`#780`/`#783` only ever deep-read her
Saturn-Titan chapter (Sec. 4.3.1). **This task deep-reads her Earth-Moon chapter, Sec. 4.4
("Resonance Transition in the Earth-Moon System"), pp.133-172, in full** — not just the three
sub-sections (4.4.1, 4.4.2, 4.4.7) `#786`'s own dispatch note flagged from a first skim.

**Scope of this pass:** full read of Sec. 4.4.1 through 4.4.7 (lines 4717-5865 of the OCR text
sidecar, PDF pages 132-181, printed pp.133-172), text-layer first (born-digital PDF, confirmed
by `#765`; no OCR artifacts found this pass either — clean extraction throughout). One direct
page-image (vision) read: Fig. 4.44 (PDF page 188, printed p.173), the ONLY place in the whole
chapter where a per-member IC axis (`x0`) exists for the Sec. 4.4.7 cycler families — load-bearing
for the "is this a catalogue gap" question, so read directly rather than trusted from the OCR
caption text alone. Ch. 3's Table 3.5 (printed p.89, referenced from 4.4.1) was cross-read from
the text layer only (not vision-verified) — it is not load-bearing for any code/catalogue claim
in this task, only a cross-reference for the resonance-labeling convention check below, and its
extraction is internally self-consistent (16 rows, sane monotonic-ish x/C values, no garbling).

---

## Headline findings (read this first)

1. **Sec. 4.4.1's genuine 1:2/2:3 homoclinic-and-heteroclinic connections at `C=2.8284`** (the
   finding `#786` already needed) are real, correctly summarized by that task, and are FIGURE/
   PROSE-ONLY — no digit-grade IC table exists anywhere in Ch. 4.4 for this or any other orbit in
   the chapter (confirmed by direct grep for `2.8284` and `2.5945`, and by reading every `Table
   4.x` in the section's line range — see §7 below).
2. **Sec. 4.4.7's "Earth-Moon Periodic Cyclers" are NOT a fully-stable family** — `#786`'s own
   note text ("separately catalogues STABLE `2:1`/`3:1` resonant cyclers") overstates Vaquero's
   own printed criterion and has been corrected in that note (see §5 below). Both families span a
   stability spectrum, confirmed directly in the Fig. 4.44 vision read.
3. **No digit-grade IC table exists for the 4.4.7 cyclers either** — Fig. 4.44's own `x0` axis is
   the ONLY per-member numeric handle Vaquero ever prints, and it carries no printed value labels
   (a genuinely non-digitizable graphical source, same gap class as the 3:4 resonant-chain family
   in her own Sec. 4.3.1, already flagged by `#765`).
4. **A specific, already-catalogued orbit (Braik-Ross 2026's R31-S) sits at almost exactly the
   Jacobi constant Vaquero prints as her 3:1 family's own upper bound (`C=3.1294` vs. `3.13`)** —
   but a direct geometric cross-check (§4 below) finds it does NOT satisfy Vaquero's own printed
   cycler-selection criteria (perigee, close lunar approach). Best current read: coincidental
   proximity in `C`, not the same orbit — but this is the single most concrete, actionable
   cross-check candidate for any future writeback/reproduction task on this family.
5. **The Casoliva/Vaquero-lineage Earth-Moon resonant-cycler CLASS has zero catalogue rows** —
   despite extensive prior digest work (`#725`, `#780`, `#786`), nothing from either paper's own
   p:q-resonance-labeled family has ever been written into `data/catalogue.yaml`. The
   Braik-Ross/Ross-Roberts-Tsoukkas manifold-tube lineage IS well represented (18 rows). This is
   the accurately-scoped gap, not "the 4.4.7 family specifically" (see §6).

---

## 1. Sec. 4.4 intro + 4.4.1 (pp.133-138) — planar homoclinic/heteroclinic connections

Motivating context (p.133): Earth-Moon has only one major satellite (unlike Jupiter/Saturn) and a
much larger mass parameter (`µ≈0.0122` vs. Saturn-Titan's `µ≈2.3658e-4`), so lunar gravity affects
periodic-orbit stability far more strongly — flagged by Vaquero herself as a design *opportunity*,
not just a complication.

**4.4.1 (pp.134-138):** takes one member each from the planar 1:2 (exterior, p<q) and 2:3
(exterior, p<q) resonant families (Figs. 4.22(a)-(b); families themselves plotted earlier in Fig.
3.17(b)-(g)), both at `C=2.8284`. Manifolds computed with a 30 km offset, 10,000 fixed points/orbit,
propagated 40 nondim time units (174 days) each direction. One-sided Poincaré maps (`x`-`ẋ` at
`y=0`, ~7-year background integration) show multiple homoclinic AND heteroclinic intersections
(Fig. 4.23) at this energy. Three connections are worked as illustrative examples (p.136-138, Fig.
4.24): a heteroclinic 1:2→2:3 transfer, a 1:2 homoclinic (which shadows a 5:4 INTERIOR resonant
orbit through the interior region — the only place a 5:4 orbit is mentioned in the whole
dissertation), and a 2:3 homoclinic (which shadows a 3:4 EXTERIOR resonant orbit, staying
exterior throughout). All three are corrected via multiple shooting from manifold-arc initial
guesses (continuity in position/velocity at patch points, no maneuvers) — the same general method
lineage as her own Sec. 4.3.1 3:4-orbit homoclinic connection and periodic resonant-chain work
(already `#765`'s anchor). Vaquero explicitly flags (p.137) that these are candidates for
periodic "resonant chains" analogous to her Saturn-Titan work, but does not construct one here.

**No numeric table anywhere in 4.4.1** — everything is figure + prose. The only digit-grade
number in the whole subsection is the shared `C=2.8284`.

## 2. Sec. 4.4.2 (pp.139-146) — 3D (spatial) extensions

Extends the same connection machinery to 3D using a 4D Poincaré-map representation (`x, ẋ, z` in
position/velocity subspace + directional arrows encoding `ż` at each crossing — a technique
credited to prior higher-dimensional-map literature, refs [83]-[88]). Worked example: an unstable
'northern'/'southern' 3D 1:2 resonant orbit and a 3D 2:3 resonant orbit, BOTH at `C=2.5945` (a
DIFFERENT, lower energy than 4.4.1's planar `C=2.8284` — not the same orbit family member). A
homoclinic-type 3D connection (1:2 self) and a heteroclinic-type 3D connection (1:2→2:3, manifolds
propagated 217 days) are both constructed (Figs. 4.28-4.29). Vaquero notes (p.146) both northern
and southern variants exist with analogous connections — not separately worked. Again, no numeric
IC table; `C=2.5945` is the only sourced number.

## 3. Sec. 4.4.3-4.4.6 (pp.147-169) — ΔV transfers from LEO to libration-point orbits

This block (planar 4.4.3, 3D 4.4.4, higher-fidelity-ephemeris 4.4.5, LPO-tour 4.4.6) is a
DIFFERENT design problem — using resonant-orbit manifolds purely as intermediate TRANSFER
MECHANISMS between a 180-km LEO and Earth-Moon libration-point orbits (L1-L5), not studying the
resonant orbits' own dynamics. Not central to this project's cycler-discovery scope, but captured
for completeness since the task asked for a full-section read:

- The `4:3` interior resonant family (Fig. 4.31, used as the general-purpose "tour the system"
  intermediate orbit for these transfers) is **entirely unstable** in the Earth-Moon system
  (`|λu|_max = 2513.2`) vs. mostly linearly STABLE in Saturn-Titan (`|λu| = 1.4704` for the few
  unstable members) — Vaquero's own direct illustration (p.151-152) of how much larger `µ`
  destabilizes resonant families. This 4:3 family is explicitly NOT one of her "cyclers" (no
  stability/perigee/period filter applied — it is chosen purely because its manifolds tour the
  whole Earth-Moon system).
- **Table 4.2** (p.158): ΔV/TOF for LEO→L1/L2/L3/L5 transfers via conic+resonant+LPO manifold
  arcs. `ΔV1` (LEO departure) 3.10-3.14 km/s always; totals `ΔVT` 3.30-3.67 km/s; TOF 2.89-31.74
  days. **Table 4.3** (p.161): a locally-optimized LEO-L5 variant, `ΔVT` 3.67→3.40 km/s, TOF
  27.68→23.78 days (fmincon, `TolX=TolCon=TolFun=1e-12`, exitflag=1). **Table 4.4** (p.164): a 3D
  LEO→L4-axial-orbit transfer via a 3:2 axial resonant family + L2 axial-orbit manifolds, `ΔVT =
  3.27` km/s, TOF 22.54 days. **Table 4.5** (p.166): the same LEO-L5 optimal transfer re-run in an
  ephemeris model with solar gravity (year 2020) — `ΔVT = 3.36` km/s vs. `3.40` km/s in the CR3BP,
  epoch-dependent (summer months costlier than winter for the year sampled).
- **4.4.6 LPO tour** (p.168-169, Fig. 4.42): LEO→L1→L2→L5→L4→L3→LEO tour using a 4:3 resonant
  manifold arc for the L2→L5 leg. `ΔV1=3.13` km/s (LEO departure), `ΔV2:7 = [374.12, 58.11, 49.21,
  80.44, 79.41, 29.59]` m/s (6 intermediate maneuvers), `ΔV8=634.17` m/s (L3→LEO return), total
  TOF between LPOs 165 days (excludes station-keeping time on each LPO itself).

None of these numbers are load-bearing for cycler discovery (this is a ΔV-transfer design study,
not a periodic-orbit family census) — recorded here only per the task's "read the whole section"
instruction, not adopted as any project constant.

## 4. Sec. 4.4.7 (pp.169-172, Figs. 4.43-4.44) — "Earth-Moon Periodic Cyclers"

**Motivation and criteria (p.169-170):** explicitly framed around lunar-infrastructure logistics
(telecom, nav, human outpost support), citing prior cycler literature (Byrnes/Longuski/Aldrin
[98], McConaghy/Russell/Longuski [100], Russell/Ocampo x3 [101]-[103], Russell/Strange [104] for
Earth-Mars/planetary-moon cyclers generally) and **explicitly adopting Mondelo & Villac's own
criteria** — i.e. refs [105]/[106], **which ARE Casoliva, Mondelo, Villac, Mease, Barrabés & Ollé's
own 2008/2010 papers**, already in this project's corpus (`#725`/`#780`). Vaquero's own four
selection criteria (quoted structure, p.169-170):

- **Close Earth approach:** insertion from a circular LEO with `180 km ≤ r-rE ≤ 35,786 km` (GEO);
  only INTERIOR resonant orbits (`p>q`) qualify by construction.
- **Close Moon approach:** either surface contact or a small-amplitude L1/L2 LPO connection
  (cislunar or circumlunar coverage both considered).
- **Period:** Earth→Moon time-of-flight `≤ 7` days.
- **Stability:** **"only resonant cyclers that are stable OR POSSESS SMALL UNSTABLE MODES are
  considered"** (direct quote, p.170) — NOT a strict stability requirement.

**Two families selected** (Fig. 4.43): planar **2:1** resonant cyclers (interior, circumlunar,
connects naturally to L2) and planar **3:1** resonant cyclers (interior, cislunar, connects
naturally to L1). Sourced numeric ranges (prose only, p.171):

| Family | Earth→Moon TOF range | Jacobi constant range |
|---|---|---|
| 2:1 | 4.91 d (`C=2.66`) to 6.39 d (`C=1.98`) | `C ∈ [1.98, 2.66]` |
| 3:1 | 4.90 d (`C=2.54`) to 5.04 d (`C=3.13`) | `C ∈ [2.54, 3.13]` |

A cost-free (unstable-to-unstable, same-`C`) transfer between the two families exists at
**`C ∈ [2.54, 2.66]` approximately** (the overlap of the two families' own ranges) — the source
text prints this range reversed, "`2.66 < C < 2.54`" (p.172); flagged here as a likely
typesetting slip (respectful-errata-framing discipline, not corrected in any project constant —
nothing in this project depends on it), since `[2.54, 2.66]` is exactly the two families' own
range overlap and the ONLY internally-consistent reading. A separate free connection to an L2
Lyapunov orbit at the 3:1 family's own boundary energy is also noted (no `C` value given).

**Fig. 4.44 vision read (PDF p.188/printed p.173) — the load-bearing check for this task.**
Direct image read confirms:
- Axes: `x0` (nd) on the horizontal, `~0.6` to `~1.2+`; Jacobi constant on the vertical, `~2.0` to
  `~3.4`. Two horizontal reference lines mark `C=C_L1` (`≈3.19`) and `C=C_L2` (`≈3.0`). A
  "Free-Transfers" band is drawn at roughly `C≈2.5-2.7`, consistent with the prose range above.
- 3:1 (blue circle, per legend) and 2:1 (red square, per legend) members are plotted as discrete
  dots along two continuous-looking curves in `(x0, C)` space, each with an attached in-plane
  stability arrow (`ν2D`, length-coded, gray) and an out-of-plane stability color (`ν3D`, via the
  right-hand colorbar spanning `0.5` to `>2.5`).
- **The `ν3D` colorbar spans across the Barden `|ν|=2` instability threshold within the plotted
  range of both families** — visibly orange/red arrows appear on the higher-`C` portion of the
  3:1 curve (near `C_L1`/`C_L2`), confirming these families are NOT uniformly stable, directly
  corroborating the "stable or small unstable modes" criterion read in §5 below.
- **No numeric value is printed on the plot anywhere** — `x0` per member must be read off the axis
  by eye; there is no accompanying table. This is a genuinely non-digitizable graphical source for
  IC purposes (same class as the Sec. 4.3.1 resonant-chain family Vaquero also never tabulates).

## 5. Correction to `#786`'s own framing (this task's own required precision check)

`#786`'s note (`docs/notes/2026-08-08-786-earth-moon-class1-resonant-connections.md`) originally
described Sec. 4.4.7 as cataloguing "STABLE `2:1`/`3:1` resonant cyclers." This is corrected (see
the note's own inline update, made this task): the source's OWN printed criterion is "stable or
possess small unstable modes," and the source's own free-transfer sentence — "a free-transfer
between two periodic orbits may exist if BOTH orbits are UNSTABLE and possess the same value of
Jacobi constant" — is stated in the SAME paragraph discussing the 2:1/3:1 cost-free transfer,
proving unstable members exist within these very families. The Fig. 4.44 vision read (§4 above)
confirms this directly and independently of the prose. The dispatch note's own "STABLE... family"
framing (`data/OUTSTANDING.md`'s `#787` bullet) is corrected by this same finding.

## 6. Cross-check against Casoliva's Table 3 (`#780`) and this project's own catalogue

**p:q convention check (precision item — the task's own explicit instruction).** Vaquero defines
`p:q` resonance (Eq. 3.1, p.83-84) as: spacecraft (body B) completes `p` revolutions in the time
the moon (body A) completes `q` revolutions; `p<q` = exterior, `p>q` = interior. Casoliva's own Eq.
(9) (`q·T_M = p·T_s`) states explicitly: "the spacecraft traverses its (inertial) elliptic orbit
`p` times, while the moon completes `q` revolutions" — **the SAME role assignment** (`p`=
spacecraft, `q`=moon), cross-validated against Casoliva's own Table 1 formula `a_s=(q/p)^(2/3)`
(e.g. `p=7,q=3` → `0.5686`, matching the table's printed `0.5684`). **This means Vaquero's `1:2`,
`2:3`, `2:1`, `3:1` labels and Casoliva's `1-2`, `2-1`, `3-2`, `7-3` labels use the identical digit
order and can be compared directly, no relabeling needed.** (A prose-gloss error inverting this
role assignment was found and corrected in `docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-
families-digest.md` this task — confirmed NOT to have propagated into
`search/earth_moon_resonant_families.py`: `Table3Row.p`/`.q` are opaque metadata from Casoliva's
own designation digits, and `satisfies_resonance`/`exists_in_em_system` are verbatim
footnote-flag transcriptions from her own printed table, never derived from a p/q formula in this
project's own code — direct read of the module confirms this.)

**Direct numeric overlap check, resonance label by resonance label:**

- **1:2 exterior** (Vaquero, `C=2.8284`, 4.4.1) vs. Casoliva's `1-2c/d/e` (`C_J = 1.5692 / 2.5803 /
  2.7630`). None match `2.8284` exactly; `1-2e` (`2.7630`) is the closest, ~2% away — same
  resonance family, evidently a different specific member, not confirmed identical.
- **2:3 exterior** (Vaquero, `C=2.8284`) — Casoliva's Table 3 has NO `2-3`/`3-2`-labeled exterior
  row at all (her only related row, `3-2c`, is `p=3,q=2`, `p>q`, INTERIOR — a genuinely different
  orbit type, not the same resonance reversed).
- **4:3 interior** (Vaquero, unstable `|λu|=2513.2`, 4.4.3 — NOT a "cycler" per her own criteria)
  — no `4-3` row in Casoliva's Table 3 at all.
- **2:1 interior cyclers** (Vaquero, `C∈[1.98,2.66]`) vs. Casoliva's `2-1a/b` (`C_J=0.4887/1.1964`)
  — **no overlap whatsoever**; Casoliva's own tabulated `2-1` members sit at much lower `C` than
  Vaquero's entire stated cycler range. Different points on presumably the same underlying
  continuous family (both methods can in principle continue a `2:1` family across a wide `C`
  range), but NOT confirmed identical — no direct reproduction attempted (out of scope, this is a
  digest task).
- **3:1 interior cyclers** (Vaquero, `C∈[2.54,3.13]`) — Casoliva's Table 3 has NO `3-1`-labeled row
  at all.
- **3:2 axial resonant family** (Vaquero, 3D, 4.4.4, used for the L4 transfer) vs. Casoliva's
  `3-2c` (planar, `C_J=0.7089`) — different dynamical class entirely (axial/out-of-plane vs.
  planar); not comparable even setting `C` aside.

**Verdict: no exact numeric identity found between any Vaquero Sec. 4.4 orbit and any Casoliva
Table 3 row.** Same resonance LABELS recur (1:2/2:1/3:2), but the specific tabulated members sit
at different `C` values or are simply absent from one paper's table — consistent with both papers
sampling different points of what MAY be shared underlying families (both use two-body/elliptical-
orbit-seeded continuation as their starting method), but this is not established without a direct
family-continuation reproduction, which this task does not attempt.

**Cross-check against `data/catalogue.yaml` — the actually load-bearing finding of this section.**
`grep -i "casoliva|mondelo|villac"` against the full catalogue returns **zero rows**. Despite three
prior tasks digesting/gating Casoliva's own paper (`#725` digest, `#780` gate module, `#786`
connection-search), **the Casoliva/Vaquero-lineage p:q-resonant Earth-Moon cycler CLASS has never
been written into the catalogue at all** — it exists only as `#780`'s own `earth_moon_resonant_
families.py` gate module (a validation harness, not a catalogue writeback). By contrast, the
Braik-Ross/Ross-Roberts-Tsoukkas manifold-tube-intersection lineage IS well represented (18
catalogue rows: `ross-rt-em-cycler-*` (1,1)/(2,1)/(3,1)/(3,2)/(3,3), `braik-ross-c21/c32-*-
corridor-*`, `braik-ross-planar-r{21,31,52}-s-corridor`, `braik-ross-common-energy-*`). **This is
the accurately-scoped catalogue gap** — not "the 4.4.7 family specifically" as `#786`'s own
dispatch framing suggested, but the entire two-body/elliptical-seed p:q-resonance-generation
METHOD lineage (Casoliva AND Vaquero both use it) vs. the well-represented manifold-tube-
intersection lineage.

**The R31-S near-coincidence (specific, checked directly).** `data/catalogue.yaml`'s
`braik-ross-planar-r31-s-corridor` row (Braik & Ross 2026 Table 2, already sourced+catalogued) has
`jacobi_constant=3.1294` — within `0.001` of Vaquero's own printed 3:1-family upper bound `C=3.13`.
A direct geometric check on whether this is the SAME orbit (Keplerian estimate, indicative not
decisive on a genuine CR3BP orbit):

- `state_nd = [x0=-0.8081272738, 0, 0, 0, ẏ0=0.1389495551, 0]`, `µ=0.0121505843`, primary (Earth)
  at `x=-µ`, secondary (Moon) at `x=1-µ` (this project's own convention, `core/cr3bp.py`).
- `r1` (Earth distance) `= |x0+µ| = 0.79598` LU `= 305,973` km. Small `ẏ0` (`≈142` m/s in
  physical units) is consistent with this being an apoapsis crossing.
- `x0` sits on the OPPOSITE side of Earth from the Moon (`x0=-0.808` vs. Earth at `-0.012` vs.
  Moon at `+0.988`) — this specific crossing point is `≈690,000` km from the Moon, nowhere near a
  close lunar approach.
- Kepler's-3rd-law estimate for a 3:1 resonance (`p=3,q=1`): `a_s=(q/p)^(2/3)=0.4807` LU (matches
  Casoliva's own Table 1 printed value exactly). If `x0` is apoapsis, implied periapsis `≈2a-r1 =
  0.1654` LU `≈ 63,585` km.
- **This periapsis estimate (`≈63,600` km) already exceeds Vaquero's own printed GEO ceiling
  (`42,164` km) for a cycler's Earth approach** — and the direct crossing point itself is nowhere
  near the Moon. On this estimate, **R31-S does not appear to satisfy Vaquero's own printed
  cycler-selection criteria**, despite the near-identical `C`.

**Read this as an indicative, not decisive, geometric estimate — weaker than it may look.** The
`0.4807`/`63,600` km figures assume a single shared two-body ellipse, but R31-S is a genuine 3:1
CR3BP orbit: it makes THREE Earth loops per period (`tof_days_bounds` on the catalogue row is the
full `27.2517`-day period, one loop is `~9.08` days), and a CR3BP resonant orbit's three perigees
need not be equal — the true minimum perigee across all three loops could sit materially below
this single-ellipse estimate. It also does not rule out a same-side close lunar approach at some
OTHER point on the orbit (only the one printed apse crossing, `~690,000` km from the Moon, was
checked). A genuine two-body-ellipse periapsis this far above the GEO ceiling (`1.5x`) makes a
same-orbit identity look unlikely, but this paragraph should NOT be read as a settled verdict —
`#798` (registered below) proposes the actual check (direct propagation of the full orbit, or a
family-continuation comparison). Separately, Vaquero's own `C=3.13` is itself a printed-to-3-sig-
fig range endpoint, not necessarily a literal match target — per this project's own
`[[feedback_published_rounded_values_are_display]]` discipline. Direction, not certainty: same
`C`, apparently different single-apse geometry — most likely a coincidental proximity in energy,
not confirmed same or different orbit without the fuller check. This is the single most concrete,
falsifiable target for
a future task that wants to resolve the question properly (full CR3BP family continuation +
direct IC comparison, not a Keplerian estimate). `braik-ross-planar-r21-s-corridor` (`C=3.1294`)
is more simply excluded — its `C` sits entirely outside Vaquero's printed 2:1 range (`[1.98,
2.66]`), no geometric check needed.

## 7. Table check (confirms no digit-grade orbit-state table exists in Sec. 4.4)

Direct grep for every `Table 4.x` reference within lines 4717-5865 of the text sidecar: **Tables
4.2, 4.3, 4.4, 4.5 all appear** (§3 above) — every one of them is a maneuver-cost/TOF table for
the LEO→LPO transfer studies (4.4.3-4.4.6), NOT an orbit-state (IC/period/Jacobi) table. No
`Table 4.6` or later exists before Sec. 4.5 begins. The nearest orbit-state table is Ch.3's **Table
3.5** (printed p.89, referenced from 4.4.1 for the general p:q family definitions) — but its
"highlighted" 1:2 member (`C=1.753856`) and 2:3 member (`C=2.500500`) are DIFFERENT specific
family members than the `C=2.8284` orbits actually used in Sec. 4.4.1's connections (Table 3.5 is
a general family-catalog reference, not the specific IC source for the connection work).

## 8. Citation-mining pass (Sec. 4.4's own background/related-work references)

Per `[[feedback_corpus_document_policy]]` step 3 — checked every reference cited across Sec. 4.4
(refs [41], [83]-[106]) against `docs/notes/CORPUS_INDEX.md` and the acquisition backlog:

- **[105]/[106] (Casoliva et al. 2008/2010)** — ARE this project's own already-corpus Mondelo &
  Villac cycler-criteria papers (`#725`), directly cited by Vaquero as the source of her own 4.4.7
  selection criteria. No new acquisition; strengthens the case that a direct Casoliva/Vaquero
  cross-reproduction is a coherent, well-motivated future task (§6).
- **[41] (Parker & Lo, "Unstable Resonant Orbits near Earth...", AIAA 2004-22819)** — same title
  as the already-corpus `lo-parker-2004-...-AIAA-2004-5304.pdf` (`#730` item 30); author-order and
  paper-number differences are consistent with a citation-formatting variant of the SAME paper,
  already reconciled by a prior task per that corpus entry's own note. No action needed.
- **[89] (Vaquero & Howell 2013, Acta Astronaut. — the JSR/Acta companion to this dissertation)**
  and **[92]/[93] (Vaquero & Howell conference papers, AAS-13-334 / IAA-AAS-DyCoSS1-05-09)** —
  [89] is already backlog item 85 (`#730` list); [92]/[93] are conference-paper variants of the
  same author/topic lineage, not independently flagged (same content class as the dissertation
  itself per this project's own `#764` scoping precedent).
- **[90] (Villac & Scheeres 2003, "Escaping Trajectories in the Hill Three-Body Problem")** —
  already backlog item 61 (HIGH priority, 3-way corroborated before this task); this is a FOURTH
  independent citation of it, reinforcing its priority further (not re-registered as a new item).
- **[91] (Davis & Howell 2011, Acta Astronaut. 69:1038-1049)** — already in corpus, digested for
  `#683`. No action.
- **[104] (Russell & Strange, "Planetary Moon Cycler Trajectories," AAS-07-118)** — the 2007
  conference precursor of the already-corpus, already-digested Russell & Strange 2009 JGCD journal
  version (`docs/notes/2026-06-30-digest-russell-strange-2009-planetary-moon-cyclers.md`). Same
  content-class precedent as `#764`; not independently acquired.
- **Four genuinely new, low-priority candidates** (LEO→L3/L4/L5 preliminary transfer-design
  papers, background context for Sec. 4.4.3-4.4.6, not resonant-orbit/cycler methodology papers)
  — registered as backlog items 100-103 (see below; the backlog file's item numbers are NOT
  strictly ordered by position, so 94-99 were already taken by earlier entries — checked directly
  before numbering), no DOI search performed this pass (honest "not searched" per the backlog's
  own methodology-notes convention):
  - [43] Perozzi, E. & Di Salvo, A., "Novel Spaceways for Reaching the Moon: An Assessment for
    Exploration," *Celest. Mech. Dyn. Astron.* 102(1-3):207-218 (2008).
  - [94] Davis, K.E., Born, G.H., Deilami, M., Larsen, A. & Butcher, E.A., "Transfers to
    Earth-Moon L3 Halo Orbits," AIAA 2012-4593.
  - [95] Larsen, A. et al., "Optimal Transfers with Guidance to the Earth-Moon L1 and L3
    Libration Points using Invariant Manifolds: A Preliminary Study," AIAA 2012-4667.
  - [97] Salazar, F.J.T., de Melo, C.F., Macau, E.E.N. & Winter, O.C., "Three-Body Problem, Its
    Lagrangian Points and How to Exploit Them Using an Alternative Transfer to L4 and L5,"
    *Celest. Mech. Dyn. Astron.* 114(1-2):201-213 (2012).

## 9. Verdict

**DIGESTED (this task's own object).** Sec. 4.4 of Vaquero 2013 is now fully read and recorded.
Nothing in this chapter is adopted as a sourced golden constant in any project code this task
(pure digest, no code/catalogue touched, per the dispatch's own explicit scope limit). The
chapter's own content is reproduction of Vaquero's already-published work, not novel by
construction — no literature-novelty framing applies to anything in this note.

**Two concrete follow-on tasks registered** in `data/OUTSTANDING.md` (`#797`, `#798`) — see below
— neither dispatched.

---

## Verification

This is a pure digest/documentation task (no code, no catalogue edits, per the dispatch's own
explicit "do not touch `data/catalogue.yaml`" instruction). Ran the mandatory `OUTSTANDING.md`
structural ratchets after editing that file:

```
uv run pytest tests/data/test_outstanding_structure.py tests/data/test_outstanding_header_body_consistency.py -q
```

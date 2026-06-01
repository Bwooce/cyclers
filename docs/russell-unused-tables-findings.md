# Russell 2004 — unused-table extraction & coverage audit

Research-only findings. Source: `docs/refs/russell-2004-dissertation.pdf`
(Russell 2004, "Global Search and Optimization for Free-Return Earth-Mars
Cyclers"). Catalogue snapshot read for reference only: `data/catalogue.yaml`
(a concurrent agent is editing it; segment ids/line numbers below are an
approximate snapshot).

This document does two things:

- **Part A** — transcribes Tables 3.5–3.8 (PDF page 99 = printed page 84),
  which give per-encounter cumulative epochs and heliocentric flyby-velocity
  components for four worked-example cyclers, and produces a mechanical
  FILL PLAN that converts the cumulative times into per-leg `tof_days`.
- **Part B** — audits Tables 4.7 & 4.8 (PDF page 139 = printed page 124)
  against the catalogue.

All numbers are transcribed verbatim from the visual PDF layout (Read tool,
`pages: "99"` / `"139"` / `"90"`), cross-checked against `pdftotext`. Where a
value is ambiguous in the scan it is flagged `[unreadable in scan]`. Nothing
is rounded, inferred, or invented.

---

## Interpretation of the Δvx / Δvy / Δvz columns (READ FIRST)

The four tables label three velocity-component rows `Δvx (km/s)`,
`Δvy (km/s)`, `Δvz (km/s)` at each encounter. **These are NOT propulsive
maneuver burns.** They are the **heliocentric velocity components of the
spacecraft at / across each ballistic gravity-assist flyby** (equivalently,
the velocity change imparted *by the flyby's gravity assist*, i.e. the
bend in the heliocentric velocity vector — not a rocket burn).

Justification:

1. **These are ballistic cyclers.** The intro paragraph above the tables
   (PDF p99) says the tables "have sufficient data to simulate one complete
   cycle plus the first leg of the second cycle for each described cycler.
   The ecliptic is the x-y axis plane and the initial position of the Earth
   is always on the x-axis." This describes a *kinematic state description*
   in an ecliptic Cartesian frame, not a burn schedule. Russell's whole
   Chapter 3 search produces *ballistic* (zero-ΔV) free-return cyclers; a
   propulsive-ΔV interpretation would contradict the ballistic premise.

2. **Table 3.3 ("Coordinates for heliocentric flyby velocities", PDF p90 =
   printed p75)** defines, for the flyby velocity diagrams (Figures 3.8a/3.8b),
   the latitude/longitude of the **incoming heliocentric excess-velocity
   vectors** `v∞1-`, `v∞2-`, `v∞3-`:

   | vector | Fig 3.8a latitude | Fig 3.8a longitude (west) | Fig 3.8b latitude | Fig 3.8b longitude (west) |
   |--------|------------------|---------------------------|-------------------|----------------------------|
   | v∞1-   | φGR              | 0                         | φGR               | λ                          |
   | v∞2-   | φFR              | 0                         | φFR               | λb+λ                       |
   | v∞3-   | φFR              | λa                        | φFR               | λb+λ                       |

   The entire flyby section is framed in terms of **heliocentric velocity
   vectors and turning angles** (Eqs. 3.5–3.8 give *geocentric turning angles*
   ωc, ωb that bend the heliocentric velocity). So the per-encounter
   (x,y,z) triples in Tables 3.5–3.8 are the resolved Cartesian components of
   exactly these heliocentric flyby velocities at each encounter.

3. **The footnote on the first column** of every one of the four tables reads
   "initial v∞ with respect to Earth" (Tables 3.5/3.6/3.8) — i.e. the first
   encounter's triple is explicitly the hyperbolic-excess velocity at the
   first Earth departure, a v∞ vector, not a ΔV burn. (Table 3.7's first
   Earth column footnote is the same "initial v∞ with respect to Earth"; its
   Mars columns carry footnote b: "0.008 AU from Mars (cycler aphelion)".)

**Bottom line for the catalogue:** if these triples are ever ingested, they
must go into velocity / `vinf`-style fields (or a heliocentric-state block),
**never** into `maneuvers[].dv_kms`. The cyclers remain ballistic
(`delta_v_kms: 0`). Do not let these be mislabeled as propulsive ΔV.

A second sourced datum these tables provide is the **Mars position vector at
t0**, footnoted `r_mars at t0 = [...] AU` under each table (values transcribed
per table below).

---

## PART A — Tables 3.5–3.8 (PDF page 99 / printed page 84)

Verbatim lead-in (PDF p99): *"Details about the flyby maneuvers for a few of
the discussed solutions are provided in Table 3.5 - Table 3.8. They have
sufficient data to simulate one complete cycle plus the first leg of the second
cycle for each described cycler. The ecliptic is the x-y axis plane and the
initial position of the Earth is always on the x-axis."*

Footnote markers used on this page: **a** = "initial v∞ with respect to Earth";
**b** (Table 3.7 only) = "0.008 AU from Mars (cycler aphelion)".

The **time (days)** row is a **cumulative epoch** at each encounter (starts at
0). Per-leg time-of-flight is therefore `t[i+1] − t[i]`.

### Table 3.5 — Cycler 2.5.1.+0

Caption (verbatim): **Table 3.5: Cycler 2.5.1.+0**

| Encounter:    | 1 (Earth) | 2 (Mars) | 3 (Earth) | 4 (Earth) | 5 (Earth) | 6 (Earth) | 7 (Mars) |
|---------------|-----------|----------|-----------|-----------|-----------|-----------|----------|
| Location      | Earth     | Mars     | Earth     | Earth     | Earth     | Earth     | Mars     |
| time (days)   | 0         | 94       | 652       | 1018      | 1200      | 1565      | 1659     |
| Δvx (km/s)    | 6.50 ᵃ    | 0        | −5.19     | 1.40      | −1.40     | −5.29     | 0        |
| Δvy (km/s)    | 4.35 ᵃ    | 0        | −1.41     | −6.12     | 6.12      | −0.98     | 0        |
| Δvz (km/s)    | 0 ᵃ       | 0        | 4.55      | 3.20      | 3.20      | 4.55      | 0        |

Footnote: `r_mars at t0 = [1.41  0.57  0] AU`.
(Encounter sequence: Earth–Mars–Earth–Earth–Earth–Earth–Mars = one cycle
E→M→…→M plus the first leg of cycle 2.)

**Catalogue match:** id `russell-ocampo-2.5.1+0` (snapshot ~line 1274).
- `sequence_canonical: "E-E-M-M"`, period 4.27 yr (k=2). 4 Earth flybys.
- Segments: `out-em` (tof_days = **94**, set), `ret-me` (tof_days **null**),
  `loop-ee` (tof_days **null**, single representative slot).
- data_gaps: `ret-me.tof_days` (unknown), `loop-ee` (derive — "to be derived
  by the multi-rev Lambert solver"), `segments` topology (uncertain),
  `out-em.a_au` (unknown).

**FILL PLAN (2.5.1.+0).** Cumulative epochs → per-leg ToF:

| leg (encounter→encounter) | arithmetic        | ToF (days) | maps to catalogue segment |
|---------------------------|-------------------|-----------:|----------------------------|
| Earth→Mars (outbound)     | 94 − 0            | **94**     | `out-em` (already 94 — confirms) |
| Mars→Earth (return)       | 652 − 94          | **558**    | `ret-me` (closes `ret-me.tof_days` gap) |
| Earth→Earth (loop 1)      | 1018 − 652        | **366**    | `loop-ee` (loop arc #1) |
| Earth→Earth (loop 2)      | 1200 − 1018       | **182**    | (loop arc #2 — needs a 2nd loop segment) |
| Earth→Earth (loop 3)      | 1565 − 1200       | **365**    | (loop arc #3 — needs a 3rd loop segment) |
| Earth→Mars (cycle-2 leg)  | 1659 − 1565       | **94**     | start of next cycle (matches 94-day outbound) |

Closes: `ret-me.tof_days = 558 d`. The single `loop-ee` slot becomes
**three** Earth-Earth loop arcs (366 / 182 / 365 d) — consistent with the
data_gaps note "4 Earth flyby(s) … implying 3 intermediate Earth-Earth
loop(s)." So the catalogue's provisional single `loop-ee` slot should expand
to 3 loop segments with these ToFs. Total cycle (E→M→E→E→E→E) = 1565 d, plus
the 94-d cycle-2 outbound leg = 1659 d total tabulated.

### Table 3.6 — Cycler 3.1.2.+1

Caption (verbatim): **Table 3.6: Cycler 3.1.2.+1**

| Encounter:    | 1 (Earth) | 2 (Mars) | 3 (Earth) | 4 (Earth) | 5 (Earth) | 6 (Mars) |
|---------------|-----------|----------|-----------|-----------|-----------|----------|
| Location      | Earth     | Mars     | Earth     | Earth     | Earth     | Mars     |
| time (days)   | 0         | 181      | 1083      | 1265      | 2348      | 2529     |
| Δvx (km/s)    | 0.71 ᵃ    | 0        | −0.09     | −1.48     | −1.28     | 0        |
| Δvy (km/s)    | 3.32 ᵃ    | 0        | −3.59     | −3.27     | 0.62      | 0        |
| Δvz (km/s)    | 0 ᵃ       | 0        | 3.39      | 3.39      | 0         | 0        |

Footnote: `r_mars at t0 = [1.15  0.99  0] AU`.
(Sequence: Earth–Mars–Earth–Earth–Earth–Mars.)

**Catalogue match:** id `russell-ocampo-3.1.2+1` (snapshot ~line 1129).
- `sequence_canonical: "E-E-M-M"`, period 6.41 yr (k=3). 3 Earth flybys.
  Note already records "Table 3.6 gives complete-cycle time 2529 d ≈ 6.93 yr".
- Segments: `out-em` (tof_days = **181**, set), `ret-me` (**null**),
  `loop-ee` (**null**).
- data_gaps: `ret-me.tof_days`, `loop-ee` (derive), `segments` (uncertain),
  `out-em.a_au`.

**FILL PLAN (3.1.2.+1).**

| leg                       | arithmetic     | ToF (days) | maps to segment |
|---------------------------|----------------|-----------:|------------------|
| Earth→Mars (outbound)     | 181 − 0        | **181**    | `out-em` (already 181 — confirms) |
| Mars→Earth (return)       | 1083 − 181     | **902**    | `ret-me` |
| Earth→Earth (loop 1)      | 1265 − 1083    | **182**    | `loop-ee` (loop arc #1) |
| Earth→Earth (loop 2)      | 2348 − 1265    | **1083**   | (loop arc #2 — needs 2nd loop segment) |
| Earth→Mars (cycle-2 leg)  | 2529 − 2348    | **181**    | next cycle outbound (matches 181-day outbound) |

Closes: `ret-me.tof_days = 902 d`. The `loop-ee` slot becomes **two**
Earth-Earth loop arcs (182 / 1083 d) — consistent with "3 Earth flyby(s) …
implying 2 intermediate Earth-Earth loop(s)." Full cycle (E→M→E→E→E) = 2348 d.

### Table 3.7 — Cycler 4.3.1.-5

Caption (verbatim): **Table 3.7: Cycler 4.3.1.-5**

| Encounter:    | 1 (Earth) | 2 (Mars ᵇ) | 3 (Earth) | 4 (Earth) | 5 (Mars ᵇ) |
|---------------|-----------|------------|-----------|-----------|------------|
| Location      | Earth     | Mars ᵇ     | Earth     | Earth     | Mars ᵇ     |
| time (days)   | 0         | 268        | 2583      | 3131      | 3399       |
| Δvx (km/s)    | −1.24 ᵃ   | 0          | 0.18      | 2.42      | 0          |
| Δvy (km/s)    | 2.84 ᵃ    | 0          | −3.24     | −2.16     | 0          |
| Δvz (km/s)    | 0 ᵃ       | 0          | 3.09      | 3.09      | 0          |

Footnotes: ᵃ "initial v∞ with respect to Earth"; ᵇ "0.008 AU from Mars
(cycler aphelion)". Footnote: `r_mars at t0 = [0.93  1.20  0] AU`.
(Sequence: Earth–Mars–Earth–Earth–Mars.)

**Catalogue match:** id `russell-ocampo-4.3.1-5` (snapshot ~line 947).
- `sequence_canonical: "E-E-M-M"`, period 8.54 yr (k=4). 2 Earth flybys.
  Note already records "Table 3.7 lists complete cycle times: 3399 d ≈ 9.30 yr"
  and "Table 3.7 has the full flyby velocity data for one cycle."
- Segments: `out-em` (tof_days = **268**, set), `ret-me` (**null**),
  `loop-ee` (**null**).
- data_gaps: `ret-me.tof_days`, `loop-ee` (derive), `segments` (uncertain),
  `out-em.a_au`.

**FILL PLAN (4.3.1.-5).**

| leg                       | arithmetic     | ToF (days) | maps to segment |
|---------------------------|----------------|-----------:|------------------|
| Earth→Mars (outbound)     | 268 − 0        | **268**    | `out-em` (already 268 — confirms) |
| Mars→Earth (return)       | 2583 − 268     | **2315**   | `ret-me` |
| Earth→Earth (loop 1)      | 3131 − 2583    | **548**    | `loop-ee` (the single loop arc) |
| Earth→Mars (cycle-2 leg)  | 3399 − 3131    | **268**    | next cycle outbound (matches 268-day outbound) |

Closes: `ret-me.tof_days = 2315 d` and `loop-ee.tof_days = 548 d`. Here the
single `loop-ee` slot is **exactly right** (2 Earth flybys → 1 loop), so this
cycler's loop gap fully closes with one value. Full cycle (E→M→E→E) = 3131 d.
Note Mars encounters here are at aphelion (0.008 AU from Mars), consistent with
the AR=0.99 "doesn't quite reach Mars" note.

### Table 3.8 — Cycler 4.5.2.-2

Caption (verbatim): **Table 3.8: Cycler 4.5.2.-2**

| Encounter:    | 1 (Earth) | 2 (Mars) | 3 (Earth) | 4 (Earth) | 5 (Earth) | 6 (Earth) | 7 (Mars) |
|---------------|-----------|----------|-----------|-----------|-----------|-----------|----------|
| Location      | Earth     | Mars     | Earth     | Earth     | Earth     | Earth     | Mars     |
| time (days)   | 0         | 191      | 1109      | 1474      | 1657      | 2022      | 3131     |
| Δvx (km/s)    | −0.71 ᵃ   | 0        | 3.38      | −3.29     | 3.29      | −1.80     | 1.29 → 0 |
| Δvy (km/s)    | 3.34 ᵃ    | 0        | −2.86     | −0.75     | 0.75      | −4.04     | 0.62 → 0 |
| Δvz (km/s)    | 0 ᵃ       | 0        | −0.50     | −2.91     | −2.91     | −0.50     | 0 → 0    |

Footnote: ᵃ "initial v∞ with respect to Earth"; `r_mars at t0 = [1.03 1.12 0] AU`.

**Column-alignment note for Table 3.8 (the one to double-check on re-ingest):**
the visual grid has **8 column slots** in the velocity rows but the header
Location/time rows show 7 encounters (Earth Mars Earth Earth Earth Earth Mars
with cumulative times 0/191/1109/1474/1657/2022/3131). The Δv rows read across
as: `−0.71 / 0 / 3.38 / −3.29 / 3.29 / −1.80 / 1.29 / 0` (x),
`3.34 / 0 / −2.86 / −0.75 / 0.75 / −4.04 / 0.62 / 0` (y),
`0 / 0 / −0.50 / −2.91 / −2.91 / −0.50 / 0 / 0` (z) — i.e. **8 entries** per
row vs **7** time columns. The final `… 0 / 0 / 0` triple is the closing
zero-vector at the terminal Mars encounter (the Mars encounters carry the
all-zero v-triple in every table, cf. Tables 3.5–3.7). The trailing
`1.29 / 0.62 / 0` is the penultimate (6th Earth, t=2022) value and the
final all-zero triple aligns to the terminal Mars (t=3131). Treat the
**7-encounter / 7-time mapping** as authoritative; the apparent 8th column is
the terminal Mars zero-vector. Re-verify visually before any edit.
(Encounter sequence as printed: Earth–Mars–Earth–Earth–Earth–Earth–Mars.)

**Catalogue match:** id `russell-ocampo-4.5.2-2` (snapshot ~line 3958).
- `sequence_canonical: "E-E-M-M"`, period 8.54 yr (k=4). 5 Earth flybys.
- Segments: `out-em` (tof_days = **191**, set), `ret-me` (**null**),
  `loop-ee` (**null**).
- data_gaps: `ret-me.tof_days`, `loop-ee` (derive), `segments` (uncertain),
  `out-em.a_au`.

**FILL PLAN (4.5.2.-2).**

| leg                       | arithmetic      | ToF (days) | maps to segment |
|---------------------------|-----------------|-----------:|------------------|
| Earth→Mars (outbound)     | 191 − 0         | **191**    | `out-em` (already 191 — confirms) |
| Mars→Earth (return)       | 1109 − 191      | **918**    | `ret-me` |
| Earth→Earth (loop 1)      | 1474 − 1109     | **365**    | `loop-ee` (loop arc #1) |
| Earth→Earth (loop 2)      | 1657 − 1474     | **183**    | (loop arc #2) |
| Earth→Earth (loop 3)      | 2022 − 1657     | **365**    | (loop arc #3) |
| Earth→Mars (cycle-2 leg)  | 3131 − 2022     | **1109**   | next-cycle leg (see caveat) |

Closes: `ret-me.tof_days = 918 d`. The single `loop-ee` slot becomes **three**
Earth-Earth loop arcs (365 / 183 / 365 d) — but note the data_gaps text says
"5 Earth flyby(s) … implying 4 intermediate Earth-Earth loop(s)." This table
only resolves **3** Earth-Earth loops between the two Mars encounters shown
(t=191 and the next Mars), so the 4-loop topology in the catalogue note is NOT
fully confirmed by Table 3.8 — the table covers "one complete cycle plus the
first leg of the second cycle," and the terminal segment 2022→3131 (1109 d)
spans into Mars again. **Flag:** reconcile the 3-resolved-loops in Table 3.8
against the catalogue's "4 loops" assumption before expanding segments; the
last 1109-d leg may itself contain an unflagged encounter or be the cycle-2
M→E return rather than a clean E→M outbound. Recommend the multi-rev Lambert
solver / topology re-check for this one rather than a purely mechanical fill.

### Part A summary of recoverable leg-times

For each cycler, the **outbound `out-em`** value is already in the catalogue
and is *confirmed* by Table 3.x; the new recoverable values are the **return
leg** and the **loop legs**:

| cycler        | `ret-me` (new) | loop legs (new)        | new leg-times recovered |
|---------------|---------------:|------------------------|------------------------:|
| 2.5.1.+0      | 558 d          | 366 / 182 / 365 d      | 4 |
| 3.1.2.+1      | 902 d          | 182 / 1083 d           | 3 |
| 4.3.1.-5      | 2315 d         | 548 d                  | 2 |
| 4.5.2.-2      | 918 d          | 365 / 183 / 365 d (caveat) | 4 |

**Total: 13 previously-null leg-times become source-recoverable** across the
four cyclers (4 + 3 + 2 + 4). All trace directly to Tables 3.5–3.8 cumulative
epochs via `t[i+1] − t[i]`. The four `out-em` values are independently
re-confirmed. Caveat: 4.5.2.-2's loop topology (3 resolved vs 4 assumed) needs
a topology re-check before mechanical expansion.

---

## PART B — Tables 4.7 & 4.8 coverage audit (PDF page 139 / printed page 124)

Captions (verbatim): **Table 4.7: Original list of ballistic cyclers** and
**Table 4.8: Rearranged list of ballistic cyclers**. Columns (both):
`v∞ (km/s)`, `Δt (yr)`, `Cycler Family`, `Characteristic Information`.
The tables do **not** continue onto p140 (p140 is text: §4.7.3 + §4.8 lead-in).

### What Tables 4.7/4.8 actually are

These are **methodological illustrations of the post-processing step (§4.7.2)**,
NOT a master enumeration of named cyclers. Verbatim context (PDF p139, §4.7.2):

> "After the algorithm has run for the entire range of v∞, a long list is
> compiled of the ballistic cyclers… **Table 4.7 illustrates an example list.**
> When v∞ is incremented by a small amount, the same cyclers will reappear
> with slightly different values for Δt. The similar cyclers are then placed
> into common bins **as illustrated in Table 4.8.**"

The rows contain only placeholder ellipses ("…") and synthetic example data:
v∞ ∈ {5.32, 5.33, 5.34, 5.35} km/s, Δt ∈ {±0.001…0.026} yr, generic
`Cycler Family` labels **A / B / C**, and `Characteristic Information` = "…"
(all ellipsis). They demonstrate the *binning + zero-Δt interpolation* method.

The single concrete result derived on the page is, verbatim:

> "From Table 4.8, cycler family C experiences a switching sign for Δt over the
> range shown from v∞ of 5.32 to 5.35 km/s. Thus, **cycler family C has a
> ballistic cycler for the true Earth-Mars synodic period with v∞ = 5.335
> km/s.**"

### Are 4.7/4.8 a superset needing new catalogue entries?

**No.** Tables 4.7/4.8 are an illustrative excerpt of the *method*, fully
superseded by the actual result tables. §4.8 ("RESULTS", PDF p139–140):
*"Table 4.9 through Table 4.12 present all of the ballistic solutions"*; Table
4.13 adds the near-ballistic three-synodic cyclers. §5.5 (PDF ~p172) states the
complete enumeration explicitly: **"1 one-synodic period parent cycler and 5
two-synodic period parent cyclers… The remaining seventy-one parent cyclers
have a three-synodic period simple model repeat time"** = **77 parent cyclers
total**, all named with the shorthand convention (e.g. `4.991gG2(#83)`,
`5.399G1(#1)` = Aldrin) and presented in Tables 4.9–4.13.

The catalogue's `russell-ch4-*` ids use exactly this Table-4.9–4.13 shorthand.

### Coverage table

| 4.7/4.8 row / item | what it is | catalogued? |
|--------------------|------------|-------------|
| Header / columns (v∞, Δt, Family, Char. Info) | schema of the illustrative method list | N/A — method illustration, not a cycler |
| Placeholder rows v∞=5.32, families A/B/C, Δt=0.020/−0.008/0.003, etc. | synthetic example data (all "Characteristic Information" = "…") | N/A — not real named cyclers; no data to ingest |
| Placeholder rows v∞=5.33 (A/B/C), 5.34 (A/B/C), 5.35 (A/B/C) | synthetic example data | N/A — same |
| **Derived result: family C → ballistic cycler at v∞ = 5.335 km/s, true E-M synodic period** | the one concrete cycler the §4.7.2 example produces | **Represented** — corresponds to the ~5.3 km/s three-synodic ballistic entries in Tables 4.9–4.13 already catalogued (`russell-ch4-5.30gGf3`, `russell-ch4-5.30ggF3`); the 5.335 value is the *interpolated* member of that v∞≈5.3 family. Not a distinct unrepresented cycler. |

### Verdict (Part B)

Tables 4.7/4.8 are **a reorganized illustrative excerpt of the post-processing
method, fully covered conceptually by Tables 4.9–4.13** (the real enumeration
of 77 parent cyclers). They contain **no named cyclers and no ingestible data**
— only placeholder/synthetic example rows. The one derived datum (family C,
v∞ = 5.335 km/s) is the interpolated representative of the v∞≈5.3 km/s
three-synodic ballistic family that is **already represented** in the catalogue
(`russell-ch4-5.30gGf3` / `russell-ch4-5.30ggF3`).

**Unrepresented entries from 4.7/4.8: 0.** No candidate new entries.
(There are currently 14 `russell-ch4-*` entries in the catalogue snapshot;
the full Russell enumeration is 77 parent cyclers in Tables 4.9–4.13 — any
gap is against 4.9–4.13, which is out of scope here, NOT against 4.7/4.8.)

---

## Files referenced (absolute paths)

- Source PDF: `/home/bruce/dev/cyclers/docs/refs/russell-2004-dissertation.pdf`
  (Tables 3.5–3.8 = PDF p99; Table 3.3 = PDF p90; Tables 4.7/4.8 = PDF p139)
- Catalogue (reference snapshot, do not edit): `/home/bruce/dev/cyclers/data/catalogue.yaml`
  - `russell-ocampo-2.5.1+0` (~line 1274)
  - `russell-ocampo-3.1.2+1` (~line 1129)
  - `russell-ocampo-4.3.1-5` (~line 947)
  - `russell-ocampo-4.5.2-2` (~line 3958)
  - `russell-ch4-5.30gGf3` (~line 35661), `russell-ch4-5.30ggF3` (~line 36306)

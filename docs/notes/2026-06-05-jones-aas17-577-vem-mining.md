# Jones / Hernandez / Jesick 2017 — VEM Triple Cycler Mining

Mined 2026-06-05. Resolves OUTSTANDING.md item **D** (Jones 2017 VEM
triple cyclers — full member list).

**Source (cite exactly this, no file path):**
Jones, Hernandez & Jesick, AAS 17-577, 2017 (NTRS 20190028464).
"Low Excess Speed Triple Cyclers of Venus, Earth, and Mars," AAS/AIAA
Astrodynamics Specialist Conference, Stevenson, WA, Aug 20-24 2017.
12 pages (incl. references). All authors JPL/Caltech.

> NOTE ON TABLE QUALITY: all four tables (Table 1 itineraries, Tables 2-3
> flyby summaries, Table 4 transit characteristics) rendered cleanly from
> the PDF. No digits were ambiguous; every cell below is transcribed
> verbatim. Where a cell is genuinely blank in the source it is shown as
> "-" (matching the paper) or noted explicitly.

---

## 1. The member tables (verbatim)

The paper gives **three** fully-tabulated cycler solutions plus the
itinerary enumeration table:

- **Table 1** (p.3): enumeration of all permitted itinerary families
  (not numeric solutions — the family list).
- **Table 2** (p.10): one EMEVVE *outbound* family cycler, full
  encounter-by-encounter listing over two repeat periods.
- **Table 3** (p.10): one MEEVEM *inbound* family cycler, full
  encounter-by-encounter listing over two repeat periods.
- **Table 4** (p.11): an *example optimal inbound+outbound pair* starting
  2020 (transit characteristics: per-encounter v∞, transit ToF, taxi Δv).

### 1.1 Table 1 (p.3) — "Triple cycler itinerary combinations"

The complete enumerated itineraries (max 6 flybys per cycle, 1 or 2
synodic periods per cycle, low-v∞ near-Hohmann transit constraint):

| Outbound | Inbound |
|---|---|
| EMEVE  | MEVEM  |
| EMEEVE | MEEVEM |
| EMEVVE | MEVVEM |
| EMEVEE | MEVEEM |
| EMMEVE | MEVEMM |

Verbatim caption: **"Table 1: Triple cycler itinerary combinations"**.

Note: column labels are Outbound / Inbound. Each outbound itinerary
begins with a near-Hohmann **E→M** arc; each inbound begins with a
near-Hohmann **M→E** arc. (p.3: "each outbound cycler begins with a
near-Hohmann Earth-to-Mars arc (Mars-to-Earth for inbound)".)

This list is the direct answer to OUTSTANDING.md item D's "member list"
at the *family* level. The catalogue's `vem-emeeve-3syn` row uses the
sequence "E-M-E-E-V-E" = EMEEVE, which **is** a listed outbound family
here (row 2 of Table 1). Confirmed.

### 1.2 Table 2 (p.10) — EMEVVE family cycler, two repeat periods

Verbatim caption: **"Table 2: Flyby summary for EMEVVE family cycler over
two repeat periods"**. Columns: Flyby body | Date | Excess speed, km/sec |
Periapsis altitude, km | Flight Time, days.

| Flyby body | Date | Excess speed, km/sec | Periapsis altitude, km | Flight Time, days |
|---|---|---|---|---|
| Earth | 07-Aug-2022 | 4.72 | 100   | -   |
| Mars  | 12-Jun-2023 | 2.50 | 4,164 | 309 |
| Earth | 01-Oct-2025 | 5.81 | 3,814 | 842 |
| Venus | 26-Apr-2030 | 7.00 | 684   | 1668 |
| Venus | 08-Feb-2034 | 7.00 | 1,985 | 1384 |
| Earth | 22-May-2035 | 4.21 | 1,998 | 468 |
| Mars  | 05-Feb-2036 | 2.79 | 1,754 | 259 |
| Earth | 15-Aug-2038 | 5.04 | 3,213 | 922 |
| Venus | 26-Dec-2039 | 4.49 | 3,319 | 498 |
| Venus | 17-Apr-2045 | 4.49 | 836   | 1939 |
| Earth | 05-Mar-2048 | 5.67 | 100   | 1053 |

Body sequence read top-to-bottom: E-M-E-V-V-E-M-E-V-V-E. That is two
EMEVVE cycles sharing the middle Earth (EMEVVE + EMEVVE overlapping at
the repeat boundary). First row's Flight Time is "-" (cycle start;
verbatim). "Flight Time" column = days **since the previous row's
encounter** (cumulative-to-previous, not cumulative-to-start; e.g.
07-Aug-2022 → 12-Jun-2023 = 309 d). Transit legs (E→M / M→E) per p.9 text
for the *outbound* example are **309 and 259 days** — matches rows
Mar/Mars 309 and Feb-2036 259.

### 1.3 Table 3 (p.10) — MEEVEM family cycler, two repeat periods

Verbatim caption: **"Table 3: Flyby summary for MEEVEM family cycler over
two repeat periods"**. Same columns as Table 2.

| Flyby body | Date | Excess speed, km/sec | Periapsis altitude, km | Flight Time, days |
|---|---|---|---|---|
| Mars  | 24-Jun-2022 | 3.85 | 100    | -    |
| Earth | 19-Mar-2023 | 3.48 | 7,831  | 268  |
| Earth | 28-Aug-2024 | 3.42 | 967    | 528  |
| Venus | 13-Mar-2028 | 5.16 | 29,777 | 1293 |
| Earth | 02-Oct-2032 | 3.84 | 2,545  | 1664 |
| Mars  | 08-Apr-2035 | 3.12 | 249    | 918  |
| Earth | 17-Nov-2035 | 2.98 | 3,719  | 223  |
| Earth | 23-Apr-2039 | 2.92 | 2,224  | 1253 |
| Venus | 23-Oct-2039 | 4.31 | 12,349 | 183  |
| Earth | 11-Aug-2045 | 4.97 | 7,484  | 2119 |
| Mars  | 21-Jan-2048 | 2.42 | 100    | 893  |

Body sequence top-to-bottom: M-E-E-V-E-M-E-E-V-E-M = two MEEVEM cycles.
p.9 text: this *inbound* example transit legs are **268 and 223 days**
(matches rows 19-Mar-2023 268 and 17-Nov-2035 223). First row Flight Time
"-" verbatim.

### 1.4 Table 4 (p.11) — example optimal inbound+outbound pair (2020)

Verbatim caption: **"Table 4: Example transit characteristics for triple
cycler transportation architecture"**. This is a two-panel table:
Outbound (left) and Inbound (right). Columns each side:
Body | Date | v∞, km/sec | Transit, days | Taxi Δv, km/sec.
The Transit column is populated only on transit legs (E→M outbound,
M→E inbound); blanks shown as empty verbatim.

**OUTBOUND panel:**

| Body | Date | v∞, km/sec | Transit, days | Taxi Δv, km/sec |
|---|---|---|---|---|
| Earth | 30-Jun-2020 | 4.76 |     | 4.22 |
| Mars  | 05-May-2021 | 2.90 | 309 | 2.24 |
| Earth | 26-Sep-2022 | 4.88 |     | 4.28 |
| Mars  | 08-May-2023 | 2.59 | 224 | 2.09 |
| Earth | 28-Sep-2024 | 3.59 |     | 3.82 |
| Mars  | 03-Aug-2025 | 2.55 | 309 | 2.07 |
| Earth | 05-Oct-2026 | 3.70 |     | 3.85 |
| Mars  | 05-Aug-2027 | 2.99 | 304 | 2.29 |
| Earth | 02-Jan-2029 | 4.20 |     | 4.02 |
| Mars  | 03-Sep-2029 | 3.76 | 244 | 2.72 |
| Earth | 21-Dec-2030 | 3.58 |     | 3.81 |
| Mars  | 11-Sep-2031 | 3.66 | 264 | 2.66 |
| Earth | 14-Apr-2033 | 4.26 |     | 4.04 |
| Mars  | 04-Dec-2033 | 3.93 | 234 | 2.82 |

**INBOUND panel:**

| Body | Date | v∞, km/sec | Transit, days | Taxi Δv, km/sec |
|---|---|---|---|---|
| Mars  | 24-Jun-2022 | 3.85 |     | 2.77 |
| Earth | 19-Mar-2023 | 3.48 | 268 | 3.78 |
| Earth | 25-Jul-2024 | 2.94 |     | 2.26 |
| Earth | 15-May-2025 | 2.95 | 294 | 3.64 |
| Earth | 26-Aug-2026 | 3.16 |     | 2.37 |
| Earth | 11-Jun-2027 | 3.47 | 289 | 3.78 |
| Earth | 05-Oct-2028 | 2.79 |     | 2.18 |
| Earth | 10-Aug-2029 | 4.71 | 309 | 4.21 |
| Earth | 03-Dec-2030 | 2.73 |     | 2.15 |
| Earth | 25-Jul-2031 | 4.61 | 234 | 4.17 |
| Mars  | 14-Feb-2033 | 2.48 |     | 2.04 |
| Earth | 21-Sep-2033 | 3.12 | 219 | 3.76 |
| Earth | 08-Apr-2035 | 2.35 |     | (none) |
| Earth | 17-Nov-2035 | 2.97 | 223 | 3.64 |

> CAUTION on Table 4 inbound: the inbound panel mixes Mars and Earth
> rows and the body labels do not form a single clean MEVEM string —
> this is the matching algorithm "mixing cycler families" (p.7: "The
> matching also permits the mixing of cycler families") across 7 crewed
> round-trip missions, **not** one cycler's itinerary. Table 4 is an
> *architecture* table (which transports a crew rendezvouses with, per
> opportunity), NOT a single-cycler member. Per row labels read
> verbatim; the row "08-Apr-2035 Earth 2.35" has no Taxi Δv entry in the
> source (blank), shown "(none)". Treat Table 4 rows as transit-leg
> v∞/Δv samples for the architecture, not as one cycler signature.

### 1.5 Worked-example numbers in the prose (p.10)

Verbatim mission walk-through (the 2020 outbound example from Table 4):
> "a crew would launch from Earth in late-June 2020 and expend 4.22
> km/sec of Δv to rendezvous with a transport vehicle already on the
> cycler. After a 309-day transit, the crew use 2.24 km/sec of taxi
> vehicle Δv to capture at Mars, where they will remain until June 2022
> (415-day stay). The crew will then expend 2.77 km/sec to rendezvous
> with an inbound transport, returning them to Earth in 268-days. The
> 992-day mission will complete with a 3.78 km/sec capture burn at Earth."

Architecture totals (p.10): "Transporting at every opportunity would
require a total of twelve cycling transport vehicles (6 inbound and 6
outbound). A total of seven round-trip crewed missions may be extracted
and analyzed from Table 4." Transit leg flight times across the
architecture: **"vary from 219 to 309 days"**; taxi vehicle Δv **"between
2.04 and 4.28 km/sec"**.

---

## 2. Method / model assumption

The model is **patched-conic on real planetary ephemeris**, NOT
circular-coplanar. Verbatim (p.5, Methodology):

> "A broad search algorithm is developed to identify near-ballistic
> cycler solutions using approximate dynamics. A zero-sphere-of-influence
> patched conic gravity model is used with the real planetary ephemeris,
> and Lambert's problem is solved to determine legs connecting consecutive
> encounters."

So **`model_assumption: analytic-ephemeris`** for any catalogue row
derived from Tables 2/3/4 (the broad-search and optimized results use the
real ephemeris). Note the Tisserand graph in Fig.2 *is* circular-coplanar
("a simplified circular coplanar Solar System model", p.4) but that is
only the qualitative explanatory tool — the actual solutions are
ephemeris-based.

Key method facts (all verbatim / closely paraphrased with page):
- **Lambert broad search** over a discrete grid of flight times; flight
  time is the primary search variable; revolutions 1→max enumerated, plus
  fast/slow (type 1 / type 2) arcs; only prograde; integer-π multiples
  excluded; cites Russell [12] for Lambert solutions. (p.5)
- **B-plane / powered-flyby targeting machinery**: turning angle
  δ = ∠(v∞⁻, v∞⁺) (Eq.1); periapsis radius r_p solved iteratively from
  Eq.2; periapsis speeds Eq.3; B-plane unit vectors Ŝ, T̂, R̂ Eq.4;
  B-plane angle θ_B Eq.5; Luidens [14] used for SOI propagation time.
  (pp.5-6)
- **Flyby feasibility**: altitude window **100 km < r_p − R_planet <
  100,000 km** (Eq.7); v∞ matching tolerance **Δv∞^max between 100 and
  200 m/sec** ("velocity increments below 200 m/sec are permitted since
  experience has shown these can be differentially corrected ... to be
  entirely ballistic"). (pp.5-7)
- **Seed-leg constraint**: seed legs are near-Hohmann, deviating **no
  more than 50 days** from the Hohmann transfer flight time; seed v∞ at
  Earth/Mars below **v∞^max = 5 km/sec**. (p.7)
- **Optimization in true ephemeris** (the "OPTIMIZED RESULTS"): two-step
  homotopy continuation, Step 1 Sun + planets, Step 2 adds all planets +
  Earth's moon; control-point/break-point model; **SNOPT** SQP optimizer
  [15]; continuity tolerance 1.0E-3 km position, 1.0E-6 km/sec velocity;
  flyby altitudes 100 km–100,000 km enforced. (pp.7-8)

---

## 3. Definitions (verbatim, p.3)

> "**Cycle:** Portion of trajectory with flight time equal to an integer
> number of T_syn, and that starts and ends at the same body (Earth or
> Mars in this work)."
> "**Repeat period:** The flight time of a single cycle."
> "**Cycler:** Trajectory that completes one or more cycles."

Synodic period (p.3):
> "The time it takes to repeat a given angular alignment of the three
> planets (the synodic period T_syn) is approximately 6.4 years. This is
> about three Earth-Mars synodic periods, and the three planets inertially
> align approximately every 32 years, or 5 T_syn."

**This settles the catalogue's open `period.basis` ambiguity (data_gaps
on `jones-2017-vem-triple-family`).** "Synodic period" in this paper =
**T_syn ≈ 6.4 yr** (the three-body VEM beat), NOT the 2.135-yr E-M
synodic. Therefore:
- "two synodic period" cycler = **2 × 6.4 = 12.8 yr** repeat period.
  Confirmed by p.9: "Recall that the repeat period T is 12.8 years."
- The E-M opportunity spacing is the *separate* 2.13 yr number: p.3 "The
  opportunities open every 2.13 years."

So the catalogue's `jones-2017-vem-triple-family` row, which currently
interprets "two synodic period" as **k=2 E-M synodic = 4.27 yr**, is
**INCORRECT per this source**. The paper's two-synodic cyclers have a
**12.8-yr** repeat period. (This is exactly the alternate reading flagged
in that row's `data_gaps[period.basis]` — the source now resolves it
toward 12.8 yr. NOTE: I am NOT editing catalogue.yaml; flagging for the
owner.)

"Low excess speed" threshold: transit-leg v∞ at Earth/Mars **below 5
km/sec** (search constraint v∞^max = 5 km/sec, p.7; abstract: "average
transit leg excess speed below 5 km/sec").

**Search scope** (p.3): "attention is restricted to families with 1 or 2
synodic periods in a cycle, and a maximum of six flybys per cycle." And
the key negative result (p.8 Broad Search Results):
> "No feasible solutions were found (of any family) with a repeat period
> of one synodic period (6.4 years). ... In contrast, thousands of
> feasible two-synodic period cyclers were obtained."

So **all feasible discovered cyclers are 2-synodic = 12.8-yr repeat
period**; the 1-synodic (6.4-yr) family is *empty* (only subsurface-flyby
near-solutions). This directly impacts both catalogue placeholder rows
(see §6).

---

## 4. Family structure

Solutions are classified by (a) integer number of synodic periods per
cycle, (b) flyby itinerary (order of bodies), and (c) initialization
year / opportunity. Verbatim (p.3): "Cycler solutions are categorized
into families based on the integer number of synodic periods in a cycle,
and the itinerary of flybys (the order of bodies encountered)."

The 10 enumerated families are Table 1 (§1.1). Outbound vs inbound is the
top-level split (outbound = favorable E→M legs; inbound = favorable M→E).
p.8: itineraries with six flybys and consecutive Earth or Venus
encounters (EMEVVE and EMEEVE) "seemed to exhibit the best overall
characteristics."

Implication for catalogue: this supports **family-level rows** (one per
Table 1 itinerary that we choose to capture) plus **specific-solution
rows** for the three tabulated example cyclers (Tables 2, 3). Both
existing placeholder rows map onto Table 1 families (EMEEVE and the
generic VEM family).

---

## 5. Relationship to prior art

- **EMV/triple-cycler concept**: p.2 "To the authors' knowledge, the
  literature contains only two references to triple cyclers,[6,7] and
  none which involve the planets." Refs [6] Moir & Barr (interstellar
  cycling) and [7] Lynam & Longuski "Laplace-resonant triple-cyclers for
  missions to Jupiter" (Acta Astronautica 2011). So this paper claims to
  be the **first planetary (VEM) triple cycler** — confirms catalogue
  novelty note. Abstract: "computed here for the first time."
- **Companion paper**: [8] Hernandez, Jones & Jesick, "Families of
  Io-Europa-Ganymede Triple Cycler," same 2017 conference — the Jovian
  sibling. (This is OUR `vem-emeeve` neighbour material — the IEG triple
  cycler is a *different* paper; do not conflate.)
- **Earth-Mars cyclers**: cites Byrnes/Longuski/Aldrin 1993 [1] (Aldrin),
  McConaghy/Landau/Yam/Longuski 2006 [2] (S1L1 "Notable two-synodic
  Earth-Mars cycler"), Russell & Ocampo 2004 [3], Pisarevsky 2008 [4].
- **S1L1 comparison (directly relevant to our S1L1 work)**: p.8 "one of
  the best Earth-Mars cyclers (the S1L1) has maximum excess speed
  exceeding 7 km/sec." p.10: "one of the best S1L1 Earth-Mars cyclers in
  Ref. [2] has flight times between 115 to 223 days, but the taxi vehicle
  Δv can be up to 5.33 km/sec at Earth and 5.64 km/sec at Mars." — these
  are *quoted from McConaghy 2006*, useful cross-checks for our S1L1 rows
  (note 5.33/5.64 km/sec Δv values, and the 7+ km/sec max v∞).
- **Venus-flyby benefit**: cites Vanderveen 1969 [9] (triple-planet
  ballistic flybys of Mars and Venus) and Okutsu & Longuski 2002 [10]
  (Mars free returns via Venus gravity assist).
- **Nothing here directly substantiates the `vem-emeeve-3syn` row's
  6.4-yr / 3×E-M = 4×E-V "beat" claim as a *discovered cycler*** — the
  6.4-yr (1-synodic) family produced **no feasible solutions** (§3). The
  6.4-yr number is real (it IS T_syn) but as a *repeat period* it is
  empty in this paper. See §6.

---

## 6. Candidate catalogue rows

Two fully-specified single cyclers + the architecture pair. Model
assumption for all = **analytic-ephemeris** (real ephemeris patched
conic; the optimized cases are ballistic in full n-body dynamics).

### Row candidate A — EMEVVE outbound cycler (Table 2)
- sequence: E-M-E-V-V-E (repeating); itinerary family **EMEVVE** (outbound)
- cycler_class: **multi-arc** (6 flybys/cycle, multiple gravity-assist arcs)
- repeat period T = **12.8 yr** (2 synodic periods of 6.4 yr)
- per-encounter v∞ (km/s), in order, cycle 1: E 4.72, M 2.50, E 5.81,
  V 7.00, V 7.00→[boundary]; cycle 2 continues E 4.21, M 2.79, E 5.04,
  V 4.49, V 4.49, E 5.67
- transit-leg ToF (E→M / M→E): **309 d and 259 d** (p.9 text)
- flyby periapsis altitudes (km): see Table 2 (100; 4,164; 3,814; 684;
  1,985; 1,998; 1,754; 3,213; 3,319; 836; 100)
- epoch: starts 07-Aug-2022 (Table 2). Epoch-specific (real ephemeris).
- start year of the figure version: 2022 (Fig.6, "two repeat periods")

### Row candidate B — MEEVEM inbound cycler (Table 3)
- sequence: M-E-E-V-E-M (repeating); itinerary family **MEEVEM** (inbound)
- cycler_class: **multi-arc**
- repeat period T = **12.8 yr**
- per-encounter v∞ (km/s) in order: M 3.85, E 3.48, E 3.42, V 5.16,
  E 3.84, M 3.12, E 2.98, E 2.92, V 4.31, E 4.97, M 2.42
- transit-leg ToF (M→E / E→M): **268 d and 223 d** (p.9 text)
- flyby periapsis altitudes (km): Table 3 (100; 7,831; 967; 29,777;
  2,545; 249; 3,719; 2,224; 12,349; 7,484; 100)
- epoch: starts 24-Jun-2022 (Table 3)

### Row candidate C — Table 4 architecture pair (2020) — NOT a single cycler
- Use as an *architecture* reference / transit-leg sample only, not a
  single-cycler signature (see §1.4 caution). Could seed a "Jones 2017
  transport architecture" note row but should not be a single-ellipse or
  single multi-arc cycler entry.

### Corrections to existing placeholder rows (flag for owner — NOT edited)
1. `jones-2017-vem-triple-family`: `period.years = 4.27` and `period.k=2`
   (E-M synodic) are **wrong per source**. The paper's two-synodic
   cyclers have **T = 12.8 yr** (2 × 6.4-yr T_syn). The `data_gaps`
   `period.basis` uncertainty is now resolved toward 12.8 yr.
   `model_assumption` should change `circular-coplanar →
   analytic-ephemeris`. `cycler_class` should be `multi-arc`.
2. `vem-emeeve-3syn`: `period.years = 6.41` (1 synodic period / 3×E-M
   beat) corresponds to a family that produced **NO feasible solutions**
   in this paper (p.8). The 6.4-yr value is T_syn (correct as a beat) but
   is **not a realized repeat period** here. Sequence "E-M-E-E-V-E" =
   EMEEVE is a valid Table 1 outbound family, but any *realized* EMEEVE
   cycler in this paper repeats at 12.8 yr, not 6.4 yr. `model_assumption`
   likewise should be `analytic-ephemeris` not `circular-coplanar`.
   `cycler_class` should be `multi-arc`.

---

## 7. Candidate golden anchors

A golden-rediscovery test needs a (sequence, period, per-leg ToF, v∞ set)
tuple traceable to the published source. Two qualify cleanly:

**Golden anchor 1 — EMEVVE outbound (Table 2 + p.9):**
- sequence EMEVVE, repeat period T = 12.8 yr
- transit legs: E→M 309 d, M→E 259 d
- transit-leg v∞ multiset (the near-Hohmann legs): Earth 4.72 / Mars 2.50
  (cycle-1 E→M); Earth 4.21 / Mars 2.79 (cycle-2). All < 5 km/s except
  the Venus-adjacent Earth flybys (5.81, 5.04, 5.67) which are
  gravity-assist legs, not transit legs.
- full per-encounter v∞ + altitude set = Table 2 (11 rows).
- Source: AAS 17-577 Table 2 (p.10) + p.9 text. EXPECTED side is
  fully source-traced (satisfies golden-discipline: no self-computed
  values).

**Golden anchor 2 — MEEVEM inbound (Table 3 + p.9):**
- sequence MEEVEM, repeat period T = 12.8 yr
- transit legs: M→E 268 d, E→M 223 d
- transit-leg v∞: Mars 3.85 / Earth 3.48 (M→E cycle 1); Earth 2.98 /
  Mars 3.12 region cycle 2 (see Table 3 ordering).
- full per-encounter v∞ + altitude set = Table 3 (11 rows).
- Source: AAS 17-577 Table 3 (p.10) + p.9 text.

**Weaker anchor — Table 4 mission walk-through (p.10):**
- 2020 outbound: launch Δv 4.22 km/s; 309-d transit; Mars capture taxi Δv
  2.24 km/s; 415-day Mars stay; inbound rendezvous Δv 2.77 km/s; 268-d
  return; Earth capture 3.78 km/s; total mission 992 days. Fully
  source-traced but it's an architecture mission profile (Δv-based), not a
  cycler orbital signature — use only as a Δv cross-check.

**S1L1 cross-check anchors (for our existing S1L1 rows):**
- p.8: S1L1 max excess speed > 7 km/sec.
- p.10: S1L1 flight times 115–223 d; taxi Δv up to 5.33 km/s at Earth,
  5.64 km/s at Mars (quoted from McConaghy 2006 [2]).

---

## 8. Honest "not extracted / ambiguous" list

- **No orbital elements (a, e, peri, apo, i) are tabulated** anywhere in
  this paper. The solutions are given purely as encounter dates +
  v∞ + periapsis altitudes + transit ToF. So catalogue
  `orbit_elements.*` for these rows must stay **null** (cannot be
  source-derived without re-running Lambert — which would be circular and
  violates golden-discipline). This is the single biggest gap.
- **Per-leg (a,e) of individual arcs** likewise not given.
- **Turning angles δ and B-plane angles θ_B** for the example cyclers are
  not tabulated (the machinery is defined Eqs.1-5 but no per-flyby δ /
  θ_B values are listed). So `flyby_mechanics.turning_angle_deg` not
  extractable; `min_altitude_km` IS extractable (periapsis altitude
  columns of Tables 2/3).
- **Total count of cyclers**: "thousands" (abstract, p.8) — no exact
  integer. `fleet_size` (12 transport vehicles for the architecture: 6
  inbound + 6 outbound, p.10) is the only hard count.
- **Table 4 inbound panel body ordering** does not form a clean single
  itinerary string (family-mixed architecture, §1.4) — do NOT infer a
  single MEVEM sequence from it.
- **Δv for the broad-search Tables 2/3 cyclers** is not given per-row;
  Δv numbers (taxi/launch/capture) appear only in Table 4 and the p.10
  walk-through, which are the optimized architecture, a different set of
  trajectories from Tables 2/3.
- **Fig.2 Tisserand graph** is circular-coplanar (qualitative only); do
  NOT treat any number read off it as a solution value.
- Epoch-independence: the abstract says low-v∞ solutions exist
  "independent of encounter epoch," but each tabulated cycler has a
  specific real-ephemeris start date; "epoch-independent" means the
  *family* property, not that any single cycler is epoch-free.

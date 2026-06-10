# McConaghy 2004 PhD dissertation — depth-mining for per-member cycler data (task #183)

**Date:** 2026-06-10
**Source:** Thomas Troy McConaghy, *Design and Optimization of Interplanetary
Spacecraft Trajectories*, PhD thesis, Purdue University, December 2004 (advisor
J. Longuski). UMI/ProQuest no. 3166673.
File: `/home/bruce/dev/cyclers_pdf/papers/mcconaghy-2004-design-optimization-interplanetary-spacecraft-trajectories-purdue-phd.pdf`
(PRIVATE repo — never linked from the public tree; see memory `cyclers_pdf private repo`).
Read with the Read tool (rendered pages). **PDF-page↔printed-page offset is +18 for
Chapters 5–7** (PDF p.147 = printed p.129); cited values below use the **printed page**
shown in the page header.
**Writeback: NONE.** This note is a structured extraction + honest ceiling assessment
+ PROPOSED updates only, held to main-session review. No edit to `data/catalogue.yaml`
or `validate.py`. Golden discipline: every value below carries a chapter/table/page
citation; nothing is a value our own code produced.

---

## 0. The single most important finding (read this first)

**The dissertation DOES contain per-MEMBER, reproducible, real-ephemeris cycler data —
but ONLY for the ballistic S1L1 cycler (Chapter 7, Tables 7.1–7.5).** For S1L1 it gives
four full DE405 itineraries (calendar date, V∞, closest-approach distance, and leg ToF
at *every* Earth and Mars encounter over 33 years). This is exactly the per-member
state needed to lift S1L1 toward V1+ and it directly bears on **#94** ("blocked:
McConaghy full text") — #94 can be UNBLOCKED on the source-availability axis.

For **every other family** (the `nUr`/`nLr`/`nSr` period classes, the `nPr` 1–6-synodic
cyclers, the n=7 15-year cyclers, the (2,1) two-synodic families incl. S1S2/U0L1/L2U0/
L3U0) the dissertation gives **only one representative summary member per descriptor**,
in the **simplified circular-coplanar** model (Tables 5.4, 5.5, 5.8) — the same kind of
summary data Russell published. **These do NOT move the validation ceiling.** There is
no per-member ephemeris table for any family except S1L1.

**Bearing on the ~200 `russell-ocampo` rows: NONE.** McConaghy does not reproduce
Russell's n.m.k catalogue and gives no per-member data for those rows. The ceiling for
the ocampo rows is **confirmed**, not lifted.

---

## 1. Document map (where the cycler content lives)

Table of contents (printed p. v–ix). Cycler material is Chapters 5–7:

- **Ch 5 "Earth-Mars Cyclers in a Simple Dynamical Model"** (printed 97–131).
  Defining characteristic `T = nS` (Eq. 5.2, p.99). Simplified model = circular
  coplanar Earth/Mars, S = 2 1/7 yr. Period-class taxonomy Tables 5.1–5.3
  (pictograms only, no numbers). **Table 5.4** (p.112) the most-promising 1≤n≤6
  `nPr` cyclers. **Table 5.5** (p.113) the n=7 (15-yr) cyclers. §5.6 the (2,1)
  two-synodic cyclers with an intermediate Earth encounter; **Tables 5.6–5.8**
  (p.121–124).
- **Ch 6 "A Nomenclature for Earth-Mars Cyclers"** (printed 132–144). The
  Russell-McConaghy joint nomenclature (`n d1 d2 … dK`, leg descriptors
  g/f/h = generic/full-rev/half-rev; Table 6.1 p.140). **Table 6.2** (p.143) labels
  for well-known cyclers, incl. the exact S1L1 label.
- **Ch 7 "Earth-Mars Cyclers in a More Realistic Dynamical Model"** (printed 145–159).
  The "ephemeris model" = JPL **DE405** for Earth/Mars; SNOPT + continuation. §7.2 the
  ballistic S1L1 cycler in the ephemeris model. **Tables 7.1–7.5** — the per-encounter
  itineraries.

---

## 2. The decisive material — S1L1 (Chapter 7 + Tables 5.8 / 6.2)

### 2.1 Exact S1L1 identification (Table 6.2, printed p.143; Table 5.8, p.124)

- **Cycler label (nomenclature):** Ballistic S1L1 cycler =
  **`2g(2.8277, 657.97 deg, U) … g(1.4508, 522.29 deg, L)`** (Table 6.2, p.143,
  references [88,89]). Read: n=2 (two-synodic repeat), K=2 generic Earth-Earth legs;
  leg 1 generic g(tf=2.8277 yr, θ=657.97°, ε=U); leg 2 generic g(tf=1.4508 yr,
  θ=522.29°, ε=L). **This pins S1L1 exactly and confirms it is a MULTI-ARC cycler
  (two generic Earth-Earth return legs), corroborating memory `S1L1 nomenclature`
  and `S1L1 real-eph closure blocker` (S1L1 is multi-arc, two generic-return arcs).**
  θ>360° on both legs ⇒ each leg is a multi-revolution Earth-Earth transfer.
- **Family / member (simple model, Table 5.8, p.124):** the ΔV-minimising member of the
  S1L1 family is **S1L1(τ=2.8277)** with, per leg (leg1, leg2):
  - ΔV per flyby: **0.00 km/s** (ballistic — required turn ≤ max possible)
  - aphelion: **1.64, 1.22 AU**
  - period: **1.49, 1.07 yr**
  - Earth V∞: **4.7, 4.7 km/s**
  - Mars V∞: **5.0, NA km/s** (only the first leg crosses Mars' orbit)
  - Mars-orbit crossing times: 0.42, 0.92, 1.91, 2.41 yr
  - §5.6 text (p.123): "only S1L1(2.8277) has an acceptably low V∞ at Earth and Mars
    (4.7 and 5.0 km/s)." S1L1's ΔV-minimum is exactly ballistic (Fig. 5.11, p.127:
    the V-shaped ΔV(τ) touches 0 at τ=2.8277).

This is a **same-model corroboration** of the catalogue's existing `russell-ch4` S1L1
descriptor (Russell 2004 Table 4.9: aphelion 1.64, g(1.4612 yr), G(2.8096 yr),
v_inf_E=4.99, v_inf_M=5.10 — per `docs/notes/2026-06-08-self-seeding-results.md` §1).
The two sources agree on aphelion (1.64) and Earth/Mars V∞ to ~0.3 km/s; the τ values
differ slightly (2.8277 vs 2.8096) because the two groups used slightly different
synodic-period / model conventions. McConaghy's leg-2 aphelion 1.22 AU is sub-Mars,
which is *why* "only the first leg crosses Mars" — consistent with the catalogue's
multi-arc framing.

### 2.2 The per-member ephemeris itineraries — Tables 7.1–7.5 (printed p.149–159)

These are the data the ocampo/SnLm rows lack and that S1L1 now HAS. The ephemeris model
uses **DE405** (printed p.145). The cycler needs **4 vehicles** (two outbound, two
inbound). Each table lists, for every encounter in a 33-year (24- or 27-body) itinerary:
encounter label (Earth-k / Mars-k), **calendar date (mm/dd/yyyy)**, **V∞ (km/s)**,
**closest-approach distance (km)**, and **leg TOF (days)**.

**Table 7.1 — Outbound Cycler Vehicle 1 (24 encounters, p.149).** Verbatim extract
(date | V∞ km/s | CA km | leg TOF d):

- Earth-1  08/13/2005 | 4.01 | — | —
- Mars-2   02/27/2006 | 3.02 | 4,816 | 198
- Earth-3  06/09/2008 | 6.89 | 20,130 | 833
- Earth-4  12/03/2009 | 6.90 | 31,110 | 541
- Mars-5   06/06/2010 | 4.31 | 17,710 | 186
- Earth-6  08/24/2012 | 6.42 | 26,490 | 809
- Earth-7  02/14/2014 | 6.43 | 41,520 | 539
- Mars-8   07/03/2014 | 7.14 | 12,190 | 138
- Earth-9  12/09/2016 | 4.01 | 27,730 | 890
- Earth-10 05/22/2018 | 4.03 | 19,920 | 530
- Mars-11  09/15/2018 | 6.47 | 11,580 | 115
- Earth-12 04/06/2021 | 4.61 | 22,990 | 934
- Earth-13 09/20/2022 | 4.59 | 14,780 | 532
- Mars-14  05/01/2023 | 2.77 | 7,601 | 223
- Earth-15 07/02/2025 | 7.08 | 23,860 | 793
- Earth-16 12/26/2026 | 7.09 | 35,120 | 542
- Mars-17  06/14/2027 | 5.27 | 13,840 | 170
- Earth-18 09/21/2029 | 5.80 | 26,850 | 830
- Earth-19 03/12/2031 | 5.80 | 37,520 | 537
- Mars-20  07/15/2031 | 7.85 | 8,802 | 125
- Earth-21 01/15/2034 | 4.21 | 24,870 | 915
- Earth-22 06/28/2035 | 4.20 | 2,756 | 529
- Mars-23  11/12/2035 | 5.87 | 1,770 | 137
- Earth-24 05/06/2038 | 7.23 | — | 906

**Table 7.2 — Outbound Cycler Vehicle 2 (p.150):** 24 encounters, launch Earth-1
09/30/2007 (V∞ 4.52) … Earth-24 05/11/2040 (V∞ 11.05). Mars-2 V∞ 3.00, CA 6,601, TOF
231 d. (Full table available in source; first/last rows quoted to characterise.)

**Table 7.3 — Inbound Cycler Vehicle 1 (p.151):** 24 encounters, Earth-1 04/01/2005
(V∞ 3.33) … Earth-24 12/13/2037 (V∞ 9.12). Mars-2 10/05/2007 V∞ 7.25, CA 12,140,
TOF 918 d.

**Table 7.4 — Inbound Cycler Vehicle 2 (p.152):** 24 encounters, Earth-1 05/29/2007
(V∞ 4.29) … Earth-24 04/07/2040 (V∞ 3.40). Mars-2 10/17/2009 V∞ 5.83, CA 4,553,
TOF 872 d.

**Table 7.5 — Outbound Vehicle 1, an itinerary near the joint optimum (27 encounters,
p.159).** This is the "best balanced" example (boxed maxima: 7.09 km/s at Earth
[Earth-16], 222.6 d max Earth-Mars TOF, 7.70 km/s at Mars [Mars-20]). Launch Earth-1
09/09/2005 (V∞ 5.32) … Earth-27 11/13/2042 (V∞ 7.87). Selected short Earth→Mars legs:
Mars-2 03/03/2006 V∞ 3.00 CA 9,961 TOF 175; Mars-14 05/01/2023 V∞ 2.77 CA 7,593 TOF
223 (boxed). Average over the itinerary (text p.154): taxi-Earth V∞ 5.50, Mars V∞ 5.06,
Earth-Mars TOF 162 d.

### 2.3 The feasibility envelope (Figs 7.2–7.5, text p.153–154)

S1L1 is not a point solution — it is a 2-parameter family of feasible itineraries
(no altitude constraint active; N eqs in N+2 unknowns ⇒ 2 DOF). For outbound vehicle 1,
**2,549** feasible (launch, arrival) pairs were sampled; **launch-date range 128 d,
arrival-date range 74 d** (Fig. 7.2). Across all itineraries (text p.154):
- smallest max taxi-Earth V∞ = **7.0864 km/s**
- smallest max Earth-Mars TOF = **222.6 d**
- smallest max Mars V∞ = **7.6897 km/s**
These are SOURCED scalar bounds that any reproduction of S1L1 must respect.

### 2.4 Why this is genuinely reproducible (V1+ enabling)

The Ch-7 tables give absolute calendar dates + a named ephemeris (DE405) + per-encounter
V∞ and closest-approach. A reproduction run can: (a) take Earth-1's date and V∞ as the
initial condition, (b) propagate against DE405 (or DE440 with a noted Δ), (c) check that
each subsequent encounter date/V∞/closest-approach matches the table within a stated
band. That is a closed, sourced, per-member target — the defining thing the V0 rows
lack. (Contrast: the ocampo rows give n.m.k + a summary V∞ and NO per-member state.)

---

## 3. The summary-only material — Chapters 5 (and what it cannot lift)

### 3.1 Table 5.4 — most-promising nPr cyclers, 1≤n≤6 (printed p.112, simple model)

Columns: aphelion radius AU | V∞ Earth km/s | V∞ Mars km/s | shortest transfer time d |
required turn deg | max possible turn deg. **One representative member per descriptor**,
circular-coplanar. Verbatim:

| nPr | aph AU | V∞ E | V∞ M | transfer d | req turn° | max turn° |
|---|---|---|---|---|---|---|
| 1L1 (Aldrin) | 2.23 | 6.54 | 9.75 | 146 | 84 | 72 |
| 2L2 | 2.33 | 10.06 | 11.27 | 158 | 134 | 44 |
| 2L3 (Byrnes Case 1) | 1.51 | 5.65 | 3.05 | 280 | 135 | 82 |
| 3L4 | 1.89 | 11.78 | 9.68 | 189 | 167 | 35 |
| 3L5 | 1.45 | 7.61 | 2.97 | 274 | 167 | 62 |
| 3S5 | 1.52 | 12.27 | 5.45 | 134 | 167 | 33 |
| 4S5 | 1.82 | 11.23 | 8.89 | 88 | 167 | 38 |
| 4S6 | 1.53 | 8.51 | 4.07 | 157 | 167 | 54 |
| 5S4 | 2.49 | 10.62 | 12.05 | 75 | 134 | 41 |
| 5S5 | 2.09 | 9.08 | 9.87 | 89 | 134 | 50 |
| 5S6 | 1.79 | 7.51 | 7.32 | 111 | 135 | 62 |
| 5S7 | 1.54 | 5.86 | 3.67 | 170 | 135 | 79 |
| 5S8 | 1.34 | 4.11 | 0.71 | 167 | 136 | 103 |
| 6S4 | 2.81 | 7.93 | 12.05 | 87 | 83 | 59 |
| 6S5 | 2.37 | 6.94 | 10.44 | 97 | 84 | 68 |
| 6S6 | 2.04 | 5.96 | 8.69 | 111 | 84 | 78 |
| 6S7 | 1.78 | 4.99 | 6.66 | 133 | 85 | 90 (ballistic) |
| 6S8 | 1.57 | 4.02 | 3.90 | 179 | 85 | 104 (ballistic) |
| 6S9 | 1.40 | 3.04 | 1.21 | 203 | 86 | 120 (ballistic) |

Notes (from table footnotes / text p.111): aphelion superscript `c` = below Mars
(1.52 AU); the 6S7/6S8/6S9 are the ballistic n=6 cyclers (req turn < max turn); required
Earth flyby periapse altitudes 1,402 / 5,408 / 13,836 km respectively. **1L1 = Aldrin**
[81]; **2L3 = Byrnes Case 1** [86].

This is genuinely the **same data class as the catalogue's existing simple-model V∞
anchors** (it corroborates 1L1's V∞ 6.54/9.75 already in the `aldrin-classic` row, which
cites Russell 2004 Table 3.4 with 6.5/9.7). It is a single exemplar per descriptor, not
a per-member ephemeris table — **cannot lift V0→V1**.

### 3.2 Table 5.5 — n=7 (15-year) cyclers (printed p.113)

Columns: r (revs/15yr) | period (15/r) yr | aphelion-radius range AU (for Rp∈(0,1)) |
yr between Earth encounters | yr between Mars encounters. r=10 = **VISIT 2** [84];
r=12 = **VISIT 1** [84]. These are PARAMETRIC ranges (aphelion given as an interval
because perihelion Rp is a free DOF), not pinned members — text p.114, Eqs 5.10–5.15
give a,e,ω analytically from (r, Rp). Useful for the VISIT-cycler rows as a
cross-check on period/aphelion bounds, but it does not pin a unique member (Rp free).

### 3.3 Table 5.8 — notable two-synodic (2,1) cyclers (printed p.124, simple model)

The per-family ΔV-minimising member (one each) for S1S2, S1S1, S1L1, U0L1, L1L1(Aldrin
& min), L2U0, L3U0. Columns per leg (leg1,leg2): ΔV/flyby | aphelion AU | period yr |
Earth V∞ | Mars V∞ | Mars-crossing times yr. Verbatim:

| cycler | ΔV | aph (l1,l2) | period (l1,l2) | V∞E (l1,l2) | V∞M (l1,l2) |
|---|---|---|---|---|---|
| S1S2(2.4885) | 0.41 | 1.83, 1.21 | 1.40, 0.71 | 13.9, 13.7 | 10.2, NA |
| S1S1(2.9124) | 0.90 | 1.62, 1.07 | 1.50, 0.95 | 3.7, 3.1 | 4.7, NA |
| **S1L1(2.8277)** | **0.00** | **1.64, 1.22** | **1.49, 1.07** | **4.7, 4.7** | **5.0, NA** |
| U0L1(2.7540) | 0.00 | 3.20, 1.54 | 2.93, 1.18 | 11.3, 11.3 | 14.0, 5.4 |
| L1L1(2 1/7) Aldrin | 1.41 | 2.23, 2.23 | 2.02, 2.02 | 6.5, 6.5 | 9.8, 9.8 |
| L1L1(2.1604) | 1.19 | 2.24, 2.22 | 2.03, 2.02 | 6.9, 6.2 | 9.9, 9.6 |
| L2U0(2.5408) | 0.00 | 1.36, 2.20 | 1.08, 1.94 | 8.8, 8.8 | NA, 10.3 |
| L3U0(2.7531) | 1.00 | 1.31, 2.29 | 0.80, 1.82 | 15.0, 15.6 | NA, 13.5 |

Three ballistic (2,1) cyclers (ΔV=0): **S1L1(2.8277), U0L1(2.7540), L2U0(2.5408)** —
text p.123. Of these only S1L1 has low V∞ at both planets, hence it (alone) is carried
to the Chapter-7 ephemeris analysis. **Table 5.7** (p.122) gives the τ-intervals where
each family's ΔV<2.5 km/s (e.g. S1L1: τ 2.794–2.860; first leg only crosses Mars).
Again: one exemplar per family, simple model — **summary, not per-member ephemeris.**

---

## 4. HONEST validation-ceiling assessment

| dissertation content | data class | does it move the ceiling? |
|---|---|---|
| **Ch 7 Tables 7.1–7.5 (S1L1)** | per-ENCOUNTER DE405 itineraries (date+V∞+CA+TOF), 4 vehicles, 33 yr | **YES — for S1L1 only.** Reproducible per-member state. Enables V1+ work on S1L1. |
| Ch 7 Figs 7.2–7.5 / text scalars | sourced feasibility bounds (launch/arrival ranges, max-V∞/TOF) | Supporting constraints for an S1L1 reproduction; not a member by themselves. |
| Ch 6 Table 6.2 | exact nomenclature labels (incl. S1L1) | Confirms topology/descriptor; pins multi-arc structure. Metadata, not a member. |
| Ch 5 Tables 5.4 / 5.5 / 5.8 | ONE simple-model exemplar per descriptor (circular-coplanar) | **NO.** Same class as Russell's summary V∞ tables. |
| (anything for `russell-ocampo` n.m.k rows) | **absent** | **NO — ceiling confirmed for the ~200 ocampo rows.** |

**Verdict:** McConaghy is **NOT** a blanket lift for the V0 backlog. It is a **targeted,
decisive lift for the SINGLE S1L1 cycler** (the catalogue's already-special multi-arc
row) plus same-model corroboration for the simple-model families. Like Russell, it
publishes summary tables + exemplars for the bulk of families; UNLIKE Russell, it adds
genuine per-member ephemeris itineraries for the one cycler (S1L1) it carries into the
realistic model. The ocampo n.m.k summary rows get **no** new per-member data — their
validation ceiling stands (memory `validation ceiling`: V0 is a publication gap; for
the ocampo rows McConaghy does not close that gap).

---

## 5. PROPOSED catalogue updates (NO writeback — for main-session review)

All sourced; nothing computed by us. To be applied (or rejected) by the main session.

1. **#94 — unblock the source-availability axis.** #94 is recorded "blocked: McConaghy
   full text." The full text is in hand and Tables 7.1–7.5 + 6.2 + 5.8 provide what #94
   needed. PROPOSE: change #94 from "blocked: McConaghy full text" to "source acquired;
   S1L1 per-member DE405 itineraries available (McConaghy 2004 Tables 7.1–7.5)." Whether
   S1L1 can be *reproduced* to V1+ is a separate compute task (see #6 below) — the
   *blocker* (missing source) is removed.

2. **S1L1 `russell-ch4` row — add corroborating source + exact label.** PROPOSE adding to
   the existing S1L1 row's `corroborating_sources` / `notes`:
   - exact nomenclature label `2g(2.8277,657.97°,U)…g(1.4508,522.29°,L)` (McConaghy 2004
     Table 6.2, p.143) — pins the multi-arc (two generic Earth-Earth leg) topology;
   - simple-model member S1L1(2.8277): aphelion 1.64/1.22 AU, period 1.49/1.07 yr,
     Earth V∞ 4.7/4.7, Mars V∞ 5.0, ΔV/flyby 0.00 (McConaghy 2004 Table 5.8, p.124) —
     an INDEPENDENT same-model corroboration of Russell 2004 Table 4.9 (aph 1.64,
     V∞E 4.99, V∞M 5.10). Note the small τ convention difference (2.8277 vs 2.8096).
   - These are corroboration of EXISTING fields, not new golden targets derived from our
     code.

3. **S1L1 — record the ephemeris (DE405) itinerary as a per-member sub-entry / data
   pointer.** PROPOSE capturing Table 7.1 (outbound vehicle 1) — and optionally 7.2–7.5
   — as a sourced per-encounter trajectory block (the segment dates + V∞ + closest
   approach + TOF are all in §2.2 above). This is the concrete artefact that could carry
   S1L1 to V1+ once a reproduction run confirms it. Mark `source_ephemeris: DE405`,
   `model_assumption: analytic-ephemeris`, `cycler_class: multi-arc`.

4. **Aldrin row corroboration (minor).** Table 5.4 (1L1: V∞E 6.54, V∞M 9.75, transfer
   146 d, req turn 84° vs max 72°) independently corroborates the existing
   `aldrin-classic` row's 6.5/9.7 V∞ and the 84°/72° turn-deficit (already cited there to
   McConaghy 2002 / Russell 2004). Optional: add Table 5.4 p.112 as a corroborating
   citation. Low value — the row is already V2.

5. **VISIT 1 / VISIT 2 rows (if present).** Table 5.5 (p.113) identifies VISIT 2 = r=10
   (period 1.5 yr, aphelion 1.62–2.62) and VISIT 1 = r=12 (period 1.25 yr, aphelion
   1.32–2.32), and Table 6.2 gives labels `7f(3:2,φ,0°)^5` and `7f(5:4,φ,0°)^3`. These
   are parametric (Rp free) — usable as period/aphelion-bound cross-checks only, NOT as
   pinned members.

6. **NO ocampo writeback.** Confirm the ~200 `russell-ocampo` rows stay V0; McConaghy
   adds nothing per-member for them. Feed this as a sourced negative into the
   negative-results registry (memory `negative results registry`): "McConaghy 2004
   surveyed — no per-member data for Russell n.m.k ocampo rows; ceiling confirmed."

7. **Follow-on compute task (not this task).** Reproduce S1L1 Table 7.1 against DE440
   (note the DE405→DE440 Δ) on the independent REBOUND/IAS15 arbiter: propagate from
   Earth-1 (08/13/2005, V∞ 4.01) and check each subsequent encounter date / V∞ /
   closest-approach within a stated band. A PASS would lift S1L1 from V0/V_descriptor
   toward V1+. This is the actionable next step the source now makes possible.

---

## 6. Honesty ledger

- Every numeric value above is transcribed from a rendered dissertation page and tagged
  with table + printed page. Nothing here is a value our code computed (golden
  discipline: EXPECTED side traces only to the published source).
- The "decisive" claim is scoped: per-member reproducible data exists **for S1L1 only**
  (Ch 7). For all other families and for the ocampo rows the dissertation is
  summary/exemplar-only — stated plainly, as the brief requires ("a clear 'no depth
  here' is a valuable, honest result"). The honest net is: **one row's worth of genuine
  lift (S1L1, the multi-arc special case), ceiling otherwise confirmed.**
- S1L1 reproduction is NOT claimed here — only that the SOURCE is now sufficient. The
  arbiter run (item 7) is required before any V-level change is asserted.
- No writeback performed. No tolerance set. No fabrication: where a full table was not
  transcribed (7.2–7.4) the first/last/representative rows are quoted and the source
  table cited for the remainder.

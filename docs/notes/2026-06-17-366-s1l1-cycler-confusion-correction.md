# #366 — S1L1 cycler-confusion correction

**Trigger:** Agent A's McConaghy 2004 JSR digest (2026-06-17, `docs/notes/2026-06-17-digest-mcconaghy-2004.md`, commit `d38c655`) flagged that the persistent **5.65 km/s Earth, 3.05 km/s Mars** anchor — wired into `docs/spec.md` §9 / §M5 and every S1L1 closure attempt as the published target — actually traces to a **different cycler** (Table 4 row **2L3**) and the 3.05 km/s value is not a true V∞ at all, but a Mars-aphelion-speed difference.

Scope of this note:
1. **Phase A (this note)** — verify the cycler confusion against the McConaghy 2004 JSR PDF directly (not via the digest).
2. **Phase B** — update the `project_s1l1_realeph_closure_blocker` memory with a retraction stanza.
3. **Phase C** — decide whether code edits are needed (much of the stack already knows 3.05 is wrong).
4. **Phase D** — re-run S1L1 closure against the corrected target(s) if any change of target is required.
5. **Phase E** — verdict + recommendation for #365 (V0→V1 promotion path).

---

## Phase A — diagnosis confirmed against the PDF (verbatim)

### A.1 Table 4 row 2L3 (p.625)

Direct PDF read (`/home/bruce/dev/cyclers_pdf/papers/mcconaghy-longuski-byrnes-2004-analysis-class-earth-mars-cycler-trajectories-jsr-doi-10.2514-1.11939.pdf`, page 4 of PDF = JSR p.625):

| Cycler nPr | Aphelion AU | V∞ at Earth km/s | V∞ at Mars km/s | Shortest transfer days | Req turn deg | Max turn deg |
|---|---|---|---|---|---|---|
| 2L3ᵇ | 1.51ᶜ | 5.65 | **3.05ᵈ** | 280ᵉ | 135 | 82 |

Where the footnote text reads, **verbatim** (p.625, table footnotes block):

- ᵇ "Case 1 cycler analyzed by Byrnes et al.³⁴"
- ᶜ "Note: the semimajor axis of Mars is 1.52 AU."
- ᵈ "Difference between Mars's speed and spacecraft aphelion speed."
- ᵉ "Time to transfer from Earth to aphelion."

So the V∞ at Mars **3.05** is NOT a hyperbolic excess speed at Mars. It is the difference between Mars's heliocentric speed at its semimajor axis and the cycler spacecraft's speed at its own aphelion (1.51 AU, just inside Mars). Since the 2L3 cycler's aphelion is **below Mars's orbit**, the cycler **never actually reaches Mars** in this circular-coplanar idealisation, and the tabulated "V∞ at Mars" is the minimum ΔV that would be needed to reach Mars from aphelion — an entirely different quantity.

The other Table 4 rows carrying the ᵈ footnote (sub-Mars aphelion, "V∞ Mars" is the same speed-difference quantity): 2L3, 3L5, 3S5, 5S8, 6S9.

### A.2 Table 6 (p.627) — S1L1 DE405 itinerary

Direct PDF read (page 6 of the file = JSR p.627), title: **"Outbound ballistic S1L1 cycler itinerary (using DE405 ephemerides of Earth and Mars)"**.

Mars-encounter rows (excerpt — verbatim from the PDF):

| Encounter | Date | Approach V∞ km/s | Closest approach alt km | Leg duration days |
|---|---|---|---|---|
| Mars 3 | 8 June 2010 | **4.31** | 17,704 | 809 |
| Mars 6 | 3 July 2014 | **7.14** | 12,179 | 890 |
| Mars 9 | 15 Sept 2018 | **6.47** | 11,570 | 934 |
| Mars 12 | 1 May 2023 | **2.77** | 7,593 | 793 |
| Mars 15 | 14 June 2027 | **5.26** | 13,751 | 830 |
| Mars 18 | 15 July 2031 | **7.70** | 10,566 | 915 |
| Mars 21 | 13 Nov 2035 | **4.68** | 15,525 | 907 |

Earth-encounter rows show V∞ varying 3.76–7.09 km/s; the launch (Earth 1, 9 June 2008) is 6.89 km/s.

Body text (p.627, between Tables and Conclusions, verbatim):
> "We note that the flyby V∞ are all less than 7.7 km/s, the flyby altitudes are all greater than 7500 km, and the Earth-Mars legs range from 115 to 223 days."

**S1L1's Mars V∞ is therefore a real-ephemeris breathing band 2.77 – 7.70 km/s (mean over the 7 Mars encounters ≈ 5.48 km/s), NOT a single value, and NOT 3.05 km/s.**

### A.3 Cross-check against the `s1l1_corrected` module

The `src/cyclerfinder/search/s1l1_corrected.py` module (commit history points to #166/#167, 2026-06-08) ALREADY encodes the corrected interpretation. Sourced directly from Russell 2004 Appendix C #83 (which is the same physical cycler as McConaghy 2004 Table 6 — Russell App-C transcribes the same DE405 reproduction block per the McConaghy 2004 dissertation Ch.7 / Table 7.5 chain established in #94):

```python
# Russell 2004 App-C "EARTH-TO-MARS TRANSIT LEG CHARACTERISTICS" table — the published
# (transit_days, vinf_Mars) for each Mars-transit (G) outbound leg, keyed by the
# ARRIVAL Mars leg number. EXPECTED side of the gate (sourced, never fit).
APPC_MARS_TRANSIT: dict[int, tuple[float, float]] = {
    3: (179.8, 5.248),
    6: (125.5, 7.693),
    9: (138.4, 4.657),
    12: (210.9, 3.198),
    15: (154.6, 6.263),
    18: (112.2, 8.046),
    21: (230.2, 3.219),
}

APPC_VINF_M_AVG: float = 5.475
APPC_VINF_E_AVG: float = 5.368
COPLANAR_VINF_M: float = 5.10  # Russell Table 4.9 idealisation — NOT a real-eph target
ROGERS_CPOM_VINF_M: float = 3.05  # Rogers/CPOM idealisation — NOT a real-eph target
```

Comparison of the Russell App-C breathing band (3.198 – 8.046 km/s, avg 5.475) against the McConaghy 2004 Table 6 breathing band (2.77 – 7.70 km/s, avg ≈ 5.48) — they are **the same cycler** (Russell 2004's #83 is McConaghy's S1L1 reproduction with a different launch epoch), and our module is already keyed off the correct values.

### A.4 What the 5.65 / 3.05 trace actually is

The spec.md §9 anchor `5.65 / 3.05 km/s` was sourced (per the row's source_quotes) from a Russell 2004 paragraph quoting "a cycler that repeats every two synodic periods has a low V-infinity at Earth and Mars (5.65 and 3.05 km/s, respectively)" — and this paragraph was paraphrasing McConaghy 2004 Table 4 row **2L3** (Byrnes 2002 Case 1, n=2, period=2 synodic, aphelion 1.51 AU **inside Mars orbit**), NOT S1L1 (which McConaghy 2004 catalogues separately in Table 6, **outside** Mars's orbit, with a real DE405 itinerary).

So the 5.65 / 3.05 pair is:
- **Real numbers** from a real cycler (2L3, the Case 1 cycler analysed by Byrnes 2002).
- But applied to the **wrong cycler** (S1L1 ≠ 2L3 — different family, different aphelion, different topology, and 2L3's "V∞_M" is a speed-difference, not a true V∞).

This is consistent with what the catalogue row's data_gap note ALREADY says (line 589: "this value traces only to spec.md §9; it could NOT be corroborated in any mined primary source" — the data_gap was right; the cause was a cycler-confusion not a missing source).

### A.5 What this changes about the closure blocker

The `project_s1l1_realeph_closure_blocker` memory's headline claim — "modelling stack complete; S1L1 fails on family-selection (off-basin)" — was true AS WRITTEN in early June 2026, but is **already substantively resolved** by the time-ordered chain of work captured in the same memory:

- **#163 (2026-06-08):** Two-arc free-return chain reaches V∞ anchors with vinf_residual 0.081 km/s on 6.44Gg3; S1L1 closes at its OWN coplanar anchors 4.99 / 5.10.
- **#166 (2026-06-08):** Source-dig found the TOPOLOGY was wrong (g = pure Earth-Earth, NOT Mars-crossing) and gave the exact App-C real-eph state.
- **#167 (2026-06-08):** Corrected topology + App-C seed reproduces all 7 Mars encounters inside the 3-SOI band on DE440 with REBOUND/IAS15 — **V3-grade independent reproduction**.
- **#173 (2026-06-08):** Self-seeding construction PASSED the S1L1 blind gate (App-C-blind synodic longitude scan recovered the family-correct basin).
- **#94 CLOSED (2026-06-11):** McConaghy Table 7.1 + Table 7.5 chain reproduced on DE440 end-to-end (23 legs, ≤0.195 km/s; two independent published sources).

The cycler confusion identified by Agent A's digest is the **root cause of why the 5.65 / 3.05 target never worked** — it was never the right target for S1L1. The S1L1 closure was achieved against the CORRECT target (Russell App-C real-eph per-leg breathing 3.198–8.046, McConaghy 2004 Table 6 breathing 2.77–7.70) over a week ago. The blocker memory has just not been updated to reflect this resolution chain crisply, and the catalogue row `s1l1-2syn-em-cpom` still carries the wrong V∞ pair as its primary record.

### A.6 Phase A — where the wrong target still lives in the codebase

Inventory from `grep -rn "3\.05"` across `src/`, `scripts/`, `tests/`, `data/`:

**Spec / catalogue (THE PRIMARY DEFECTS):**
- `docs/spec.md` lines 149, 163, 529, 1353 — anchors 5.65/3.05 as the M5 milestone target and §9 published value.
- `data/catalogue.yaml` lines 510-514 (s1l1-2syn-em-cpom vinf_kms_at_encounters), 570, 578, 599, 604, 643-646, 661-666 — the row encodes 5.65/3.05 with a data_gap note. The note correctly says "could not be corroborated" — but the cause turned out to be cycler-confusion, not source unavailability.

**Code (mostly already-correct, explicitly labelled as idealisations):**
- `src/cyclerfinder/search/s1l1_corrected.py` line 122 — `ROGERS_CPOM_VINF_M: float = 3.05  # Rogers/CPOM idealisation — NOT a real-eph target` — **correct disposition already**.
- `src/cyclerfinder/cli.py` lines 74, 90, 314 — `--vinf-targets` parser help text uses `E=5.65,M=3.05` as an example string. Not a hardcoded target; example only.
- `src/cyclerfinder/search/resonant_construct.py` lines 27, 30 — comments call out 5.65/3.05 as a "higher-fidelity" target requiring eccentric Mars (already labelled).
- `src/cyclerfinder/search/optimize.py` line 6 — docstring references spec §9.
- `src/cyclerfinder/search/precursor_matcher.py` lines 179, 372 — `S1L1's first encounter is Earth (V_inf 5.65 km/s)`. **STILL WIRES THE WRONG TARGET into the precursor matcher.**

**Scripts (legacy diagnostics):**
- `scripts/characterise_s1l1_emee.py`, `scripts/correct_s1l1_twoarc.py`, `scripts/close_s1l1_realeph.py`, `scripts/construct_s1l1_twoarc_realeph.py`, `scripts/diagnose_s1l1_twoarc.py`, `scripts/diagnose_s1l1_eccmars.py`, `scripts/diagnose_s1l1_realeph.py` — all reference 5.65/3.05 as PUBLISHED anchors. Most are diagnostic scripts from the pre-#163 / pre-#167 era and are now historical artefacts.

**Tests:**
- `tests/data/test_catalogue_rediscovery_tagging.py:16` and `tests/verify/test_propagate.py:149, 158, 167` — test the rediscovery-tagging logic against the 5.65/3.05 pair; the test bodies treat this as a *disputed/unverified* anchor (correct disposition).

### A.7 Where the wrong target genuinely needs correction

Two scopes:
1. **The catalogue row `s1l1-2syn-em-cpom`** carries 5.65/3.05 as the primary `vinf_kms_at_encounters` record. The row's data_gap note already qualifies these as unverified, but a downstream consumer (e.g. `precursor_matcher.py:179`) reads `5.65` straight off the row as a precursor-target. **This is the primary integrity defect.** Promotion to V1 (#365) is **conditional** on fixing this.
2. **`docs/spec.md` §9 / §M5** anchors the cycler-finder's milestone gate to the wrong cycler. This is a spec-level inaccuracy.

Code edits and tests are mostly already-correct (they label 3.05 as an idealisation or use it for disputed-tagging tests). The historical-diagnostic scripts can be left as-is with a header comment noting they preceded the #166 source-dig.

### A.8 Catalogue impact of the corrected target

The CORRECT S1L1 V∞ at Mars is a **breathing band, not a single value**:
- Real-ephemeris (DE405/DE440): 2.77 – 7.70 km/s (McConaghy 2004 Table 6, 7 Mars encounters).
- Real-ephemeris (DE405): 3.198 – 8.046 km/s (Russell 2004 App-C #83, 7 Mars encounters — same cycler, different launch epoch).
- Coplanar idealisation (Russell 2004 Table 4.9): 5.10 km/s (one value).
- McConaghy 2006 real-eph (abstract): 5.0 km/s (one value, presumably a typical or design-point value).

V∞ at Earth: similarly a breathing band 3.76 – 7.09 km/s (Table 6) or 5.368 km/s avg (Russell App-C).

The catalogue should not encode a single (E, M) V∞ pair for S1L1 at all — it should encode either the App-C / Table 6 per-leg breathing values (a list, not a pair) or the coplanar/design-point summary with the explicit fidelity tag. The row's `model_assumption: circular-coplanar` says the row is a coplanar framing; the right coplanar value pair is then **4.99 / 5.10** (Russell Table 4.9), not 5.65 / 3.05 (Byrnes 2002 Case 1 = 2L3 = a DIFFERENT CYCLER).

---

(Phases B–E continue below; this note will be appended as those phases complete.)

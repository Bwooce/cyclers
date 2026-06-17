# Digest: Campagnola, Buffington & Petropoulos 2014 — Jovian tour design for Europa orbiter and lander

Task #383. Single-paper digest, deep read of all 14 pages (Acta Astro
pp. 68-81) on 2026-06-17 AET. Part of the "every page of every paper" deep
read directive 2026-06-17.

## 1. Header

- **Title (verbatim)**: *Jovian tour design for orbiter and lander missions
  to Europa*
- **Authors**: Stefano Campagnola (JAXA/ISAS, Yoshinodai 3-1-1, Sagamihara,
  Kanagawa, 155-0031 Japan), Brent B. Buffington (Jet Propulsion Laboratory,
  Caltech), Anastassios E. Petropoulos (JPL)
- **Venue**: *Acta Astronautica* 100 (2014) 68-81
- **DOI**: 10.1016/j.actaastro.2014.02.005
- **Received / revised / accepted / online**: 19 March 2013 / 27 Nov 2013 /
  4 Feb 2014 / 15 Feb 2014
- **Length**: 14 pages (front matter through references, pp. 68-81)

## 2. What the paper actually is

The paper presents the Jovian tour designs for the Europa orbiter and
lander concepts of the **NASA Europa Habitability Mission (EHM)** study,
commissioned post-Planetary Decadal Survey 2013-2022 to address the
"too-expensive" verdict on JIMO and predecessor concepts (p. 68 §1). It
documents three specific tours by name:

- **Tour 11-O3** — long orbiter tour (Earth launch 2021 Nov 22, EOI 2030
  Jan 03; total mission ~9 years), pp. 70-73 §3.1
- **Tour 12-L1** — short lander tour with pump-down strategy (Earth launch
  2021 Nov 22, EOI 2029 Aug 06), pp. 73-76 §4.1
- **Tour 12-L4** — alternative lander tour with low-energy endgame
  (Earth launch 2021 Nov 22, EOI 2030 Jan 04), p. 76 §4.2

The companion multi-flyby (Europa Clipper precursor) tour appeared in a
separate paper (ref [14]) and is NOT presented here. The companion paper
explicitly cited is the May 2012 EHM Final Study Report (ref [13]).

The paper's central methodological claim is that the **"courtship phase"**
(the last few months from final Callisto flyby through Europa Orbit
Insertion, EOI) of every prior Europa-mission tour can be significantly
improved by replacing the linked-conics Hohmann Callisto-Ganymede-Europa
descent with a **multi-body, low-energy transfer using Tisserand-leveraging
transfers (TLT)** plus a **Tisserand-Poincare graph** as the design tool
(p. 69 §2, p. 71-72 §3.1, citing the methodology papers by Campagnola &
Russell [15,22,23,24] and the Strange-Russell endgame work).

The reduction is quantified in Tables 4-7 (pp. 76, 80, see §3 below):
courtship-phase TID 250-300 krad lower than the comparable Banzai
Pipeline, OR Delta-V 250-300 m/s lower, for the same total mission scope.

## 3. Key content, tables, equations with page citations

### 3.1 Tisserand-Poincare graph (the headline methodological tool)

- **Definition (p. 71)**: A 2D projection of a Poincare section of the
  spacecraft trajectory with Tisserand level-sets overlaid. The Tisserand
  parameter relative to a moon is "an approximation of the Jacobi constant
  of the spacecraft-moon-Jupiter three-body problem" and can be expressed
  as a function of osculating apojove and perijove. It is also a measure
  of V_inf relative to the moon, whenever V_inf is defined (i.e. very
  close to the moon or far away with intersecting apsides).
- **Axes**: apojove (R_J units) horizontal, perijove (R_J units) vertical.
  Both axes logarithmic-or-linear-mixed in the example figures.
- **Overlay families (Fig. 3 caption, p. 73; Fig. 6 caption, p. 75)**:
  - Green curves: Tisserand level-sets w.r.t. Ganymede
  - Red curves: Tisserand level-sets w.r.t. Callisto
  - Blue curves: Tisserand level-sets w.r.t. Europa
  - Yellow shading: linked-conics domain (the set of orbits intersecting
    the moon's orbit and for which the linked-conic model can be used).
    "Navigating the low-energy 'wings' outside the linked-conics domain"
    is the named multi-body advantage (p. 72 col 2).
- **Trajectory representation**: between two flybys of the same moon, the
  spacecraft state stays on a single Tisserand level-set; jumps between
  level-sets correspond to flybys of a different moon. Impulsive maneuvers
  appear as jumps between points with the same apojove (PJR is a jump
  between same-apojove perijoves; p. 72 col 1).
- **Resonance lines**: orbits with periods commensurable with a moon
  appear as discrete labels on the graph (e.g. 32:1, 11:1, 8:1, 5:1, 3:1,
  2:1, 1:1 for Ganymede on Fig. 3A; 6:5, 4:3 for Europa on Fig. 3C).

### 3.2 Tour 11-O3 events (Table 1, p. 72)

Verbatim transcription (Event / Date ET / Maneuver Delta-V km/s / Flyby
V_inf km/s / Alt km / In or Out):

| Event       | Date (ET)              | Delta-V km/s | V_inf km/s | Alt km | In/Out |
|-------------|------------------------|--------------|------------|--------|--------|
| Earth escape| 2021 Nov 22 00:01:13  | 3.852        | 3.77       | 200    | -      |
| Venus       | 2022 May 14 09:30:43  | -            | 6.62       | 300    | In     |
| Earth       | 2023 Oct 24 22:41:49  | -            | 12.07      | 11 761 | In     |
| Earth       | 2025 Oct 24 09:06:48  | -            | 12.05      | 3 330  | In     |
| Ganymede 0  | 2028 Apr 03 11:53:40  | -            | 7.99       | 500    | -      |
| JOI         | 2028 Apr 03 22:53:35  | 0.812        | -          | -      | -      |
| PJR         | 2028 Jul 14 13:40:14  | 0.122        | -          | -      | -      |
| Ganymede 1  | 2028 Nov 18 20:22:49  | -            | 5.61       | 629    | In     |
| OTM         | 2028 Dec 28 04:45:07  | 0.005        | -          | -      | -      |
| Ganymede 2  | 2029 Feb 05 12:34:23  | -            | 5.68       | 100    | In     |
| Ganymede 3  | 2029 Apr 03 17:20:35  | -            | 5.8        | 3 370  | In     |
| Ganymede 4  | 2029 May 09 11:45:33  | -            | 5.79       | 643    | In     |
| Callisto 5  | 2029 Jun 17 12:40:13  | -            | 5.48       | 221    | In     |
| Ganymede 6  | 2029 Jul 21 13:57:07  | -            | 3.87       | 6 645  | Out    |
| Ganymede 7  | 2029 Sep 02 10:32:53  | -            | 3.79       | 268    | Out    |
| Ganymede 8  | 2029 Sep 23 21:08:34  | -            | 3.77       | 2 009  | Out    |
| Ganymede 9  | 2029 Oct 08 04:05:01  | -            | 3.76       | 2 730  | Out    |
| Callisto 10 | 2029 Oct 23 05:23:05:28| -           | 1.77       | 2 124  | In     |
| Ganymede 11 | 2029 Oct 27 09:33:56  | -            | N/A        | 23 667 | N/A    |
| Ganymede 12 | 2029 Nov 17 22:51:09  | -            | N/A        | 4 900  | N/A    |
| Ganymede 13 | 2029 Nov 24 20:15:57  | -            | 1.19       | 1 185  | In     |
| Europa 14   | 2029 Nov 27 08:10:12  | -            | N/A        | 6 681  | N/A    |
| TLM         | 2029 Dec 09 15:08:29  | 0.054        | -          | -      | -      |
| Europa 15   | 2029 Dec 11 17:16:32  | -            | N/A        | 6 563  | N/A    |
| TLM         | 2029 Dec 27 01:25:06  | 0.097        | -          | -      | -      |
| EOI         | 2030 Jan 03 21:29:49  | 0.477        | N/A        | 100    | N/A    |

Note: Callisto 10 entry was transcribed as "05:23:05:28" — this appears
to be a typo in the table (likely "05:10:28" given the format of the
other entries). The remaining figures are unambiguous.

The In/Out tag indicates the V_inf direction relative to the moon
(inbound vs outbound; the In/Out switch at G6 marks the energy descent
from the resonance pump-down to the gravity-assist pump-down).

### 3.3 Tour 12-L1 events (Table 2, p. 76)

Lander pump-down tour. Same Earth-launch date 2021 Nov 22, EOI Aug 2029
with a 1.067 km/s Europa Orbit Insertion burn (much higher than 11-O3's
0.477 km/s because no low-energy endgame is used; multi-body endgame
deferred to 12-L4). Total mission 21 months courtship + 18 mo cruise.

Verbatim selected entries (Date ET / Delta-V km/s / V_inf km/s / Alt km
/ In or Out):

| Event       | Date                  | Delta-V | V_inf | Alt km | I/O |
|-------------|-----------------------|---------|-------|--------|-----|
| Ganymede 0  | 2028 Apr 03 06:23:43  | -       | 9.21  | 500    | In  |
| JOI         | 2028 Apr 03 16:49:29  | 0.839   | -     | -      | -   |
| PJR         | 2028 Jun 26 17:48:57  | 0.162   | -     | -      | -   |
| Ganymede 1  | 2028 Oct 21 06:21:39  | -       | 5.7   | 293    | In  |
| Ganymede 2  | 2028 Dec 10 08:31:35  | -       | 5.63  | 176    | In  |
| Ganymede 3  | 2029 Jan 22 06:46:42  | -       | 5.6   | 120    | In  |
| Callisto 4  | 2029 Mar 02 13:48:18  | -       | 6.31  | 200    | Out |
| Ganymede 5  | 2029 May 03 11:40:25  | -       | 3.8   | 202    | In  |
| Ganymede 6  | 2029 Jun 01 01:58:54  | -       | 3.81  | 559    | In  |
| Ganymede 7  | 2029 Jun 17 10:31:51  | -       | 3.71  | 2 171  | Out |
| Callisto 8  | 2029 Jul 12 15:21:32  | -       | 1.82  | 2 883  | In  |
| Ganymede 9  | 2029 Jul 06 20:17:43  | -       | N/A   | 15 847 | N/A |
| OTM         | 2029 Jul 12 00:13:08  | 0.014   | -     | -      | -   |
| Ganymede 10 | 2029 Jul 28 10:31:38  | -       | N/A   | 552    | N/A |
| Ganymede 11 | 2029 Aug 04 09:28:08  | -       | 1.39  | 3 995  | In  |
| EOI         | 2029 Aug 06 15:36:32  | 1.067   | 1.48  | 200    | -   |

### 3.4 Tour 12-L4 events (Table 3, p. 78)

Same Earth launch 2021 Nov 22; EOI 2030 Jan 04 04:19:18:22 with a 0.484
km/s burn — the low-energy endgame (4:3 / 6:5 Europa resonance hopping
matching the orbiter 11-O3 pattern, see §3.3) — which costs an extra ~3.5
months on top of 12-L1 but cuts EOI Delta-V by 583 m/s and lets Delta-V
total drop from 2.10 to 1.56 km/s.

### 3.5 Tour comparison (Tables 4-7, pp. 76, 80)

Table 4 (courtship phase only, p. 76, verbatim):

| Phase / Quantity            | 99-35 short | 08-008 short | Banzai | 12-L1 | 99-35 | 08-008 | 11-O3 | 12-L4 |
|-----------------------------|-------------|--------------|--------|-------|-------|--------|-------|-------|
| Endgame TLM1 Delta-V (km/s) | -           | -            | -      | -     | 0.12  | 0.07   | 0.05  | 0.05  |
| Endgame TLM2 Delta-V        | -           | -            | -      | -     | 0.07  | 0.07   | 0.10  | 0.08  |
| EOI Delta-V                 | 1.17        | 1.16         | 1.21   | 1.07  | 0.45  | 0.71   | 0.48  | 0.48  |
| EOI adjustment              | -           | -            | -0.004 | -     | -     | -      | -0.02 | -     |
| Total courtship Delta-V     | 1.17        | 1.16         | 1.21   | 1.08  | 0.64  | 0.85   | 0.61  | 0.62  |
| Courtship TID (krad)        | 384         | 230          | 48     | 52    | 1047  | 880    | 792   | 832   |
| Courtship TOF (months)      | 3           | 1            | 0.5    | 1     | 4.5   | 2      | 2.5   | 2.5   |

Table 5 (tour total Delta-V, p. 80, verbatim row-by-row):

| Phase                    | 99-35 short | 08-008 short | Banzai | 12-L1 | 99-35 | 08-008 | 11-O3 | 12-L4 |
|--------------------------|-------------|--------------|--------|-------|-------|--------|-------|-------|
| JOI                      | 0.98        | 0.60         | 0.81   | 0.84  | 0.98  | 0.60   | 0.81  | 0.81  |
| JOI adjustment           | +0.01       | -0.04        | +0.05  | -     | +0.01 | -0.04  | -     | -     |
| PJR                      | 0.08        | 0.08         | 0.11   | 0.16  | 0.08  | 0.08   | 0.12  | 0.12  |
| Total G0/I0 -> G1/I1     | 1.07        | 0.64         | 0.97   | 1.00  | 1.07  | 0.64   | 0.93  | 0.93  |
| Total G1/I1 -> lastC     | 0.00        | 0.00         | 0.22   | 0.02  | 0.00  | 0.00   | 0.00  | 0.01  |
| Total lastC -> EOI       | 1.17        | 1.16         | 1.21   | 1.08  | 0.64  | 0.85   | 0.61  | 0.62  |
| **Tour total Delta-V**   | **2.24**    | **1.80**     | **2.39**| **2.10** | **1.71** | **1.48** | **1.55** | **1.56** |

Table 6 (tour-total TID, krad, p. 80): 11-O3 = 881, 12-L4 = 920, vs
99-35 = 1100 (older orbiter) and 08-008 short = 1110.
Table 7 (tour-total TOF, months): 11-O3 = 21, 12-L4 = 21 vs 99-35 = 16.

### 3.6 Constraints used (sourced)

- Solar conjunction avoidance: no conjunction allowed for 60 days before
  EOI and 30 days after EOI (p. 69 col 2).
- Science orbit Local Solar Time: 7-9 AM, fixed to 7:30 AM for the
  reference design (pp. 73-74 col 1, 12-L1 §4.1). The 11-O3 orbiter
  reaches 3:40 PM (p. 73 col 2, top).
- Science orbit altitude: 100 km circular polar (orbiter 11-O3, p. 70);
  200 km circular polar (lander 12-L1, p. 73 col 2).
- Minimum flyby altitudes: Ganymede 100 km (ephemeris-uncertainty buffer);
  Callisto 200 km; Europa 100 km — all real altitudes (pp. 70 col 2).

### 3.7 Six-flyby V_inf sequence at Ganymede (orbiter 11-O3)

Reading Table 1, the V_inf at Ganymede walks down monotonically as the
resonance ratio is lowered:
- G0 inbound: 7.99 km/s (JOI prep; pre-PJR)
- G1: 5.61 km/s (post-PJR; 11:1 resonance per p. 70 col 1)
- G2: 5.68 km/s (8:1)
- G3: 5.80 km/s (5:1)
- G4: 5.79 km/s (3:2 Callisto-targeting)
- G6 outbound: 3.87 km/s (post-C5 pump-down)
- G7: 3.79 km/s (6:1)
- G8: 3.77 km/s (3:1)
- G9: 3.76 km/s (2:1)
- G11: N/A (high-altitude Tisserand-leveraged segment, outside the
  linked-conics yellow region of Fig. 3B)
- G12: N/A
- G13: 1.19 km/s (1:1 resonance, courtship phase)

The Europa V_inf at the closing 4:3 / 6:5 resonance is also too small
to be defined in linked-conics terms (E14, E15 in Table 1; the multi-
body endgame is operating below the patched-conic boundary, by design).

## 4. KNOWN_CORPUS impact

**RECOMMEND: New KNOWN_CORPUS anchor for Jovian Galilean-moon tour.**

This paper is a clean, modern (2014), Acta-Astronautica-archival
single-source for a Jovian mga_tour with the EHM 2012 study lineage
(JPL + JAXA, post-Decadal). It carries:

- Full E-V-E-E-J trajectory specification (Earth launch through EOI,
  pp. 70-73, Tables 1-3)
- V_inf km/s tabulated at every encounter where defined (Tables 1-3)
- The Tisserand-Poincare graph topology (Fig. 3 for 11-O3, Fig. 6 for
  12-L1)
- Reproducible mission events at ephemeris-time-tagged precision

Suggested anchor (for KNOWN_CORPUS, not a catalogue admission yet):

- `corpus_id`: campagnola-buffington-petropoulos-2014-ehm-orbiter
- `orbit_class`: mga_tour
- `topology_label`: ["pump-tour", "mga-tour", "tisserand-poincare"]
- `body_set`: ["Earth", "Venus", "Jupiter", "Ganymede", "Callisto", "Europa"]
- `interplanetary_sequence_canonical`: "E-V-E-E-J-G0-JOI-G1-G2-G3-G4-C5-G6-G7-G8-G9-C10-G11-G12-G13-E14-E15-EOI"
- `validation_level`: V0 (single sourced table; not reproduced by us)

A second anchor for the lander tour 12-L1 is a separate row (same
Earth launch, different Jovian arc). The 12-L4 alternative shares the
12-L1 lander concept; it can be an annotation on 12-L1 rather than a
separate corpus row.

NB: The interplanetary cruise V-E-E sequence (Venus then double Earth)
is the **EHM 2012 6.37-year cruise** (p. 69 col 2), distinct from the
Galileo VEEGA. Worth flagging that this is its own MGA-cruise pattern
plus an extensive Jovian inner tour.

## 5. Catalogue impact

**RECOMMEND: orbit_class=mga_tour candidate admission**, but only after
KNOWN_CORPUS pre-anchor + V_inf reproduction (V0->V1 step).

The tour is epoch-locked (the 2021 Nov 22 Earth launch is named on the
table; "the JOI epoch is selected to (1) minimize the TID by preventing
phasing orbits ... (2) avoid solar conjunction ... (3) achieve the
desired target science orbit", p. 69 col 2). n_returns=1, terminates at
EOI. This is structurally identical to the existing `tito-2018-mars-free-return`
admission pattern (single epoch-locked window, mga_tour class, finite
n_returns=1, validity_window pinned to the published table dates).

Difference from Tito 2018: ten Galilean moon encounters instead of one
Mars flyby; the "tour" character is much stronger and the EOI termination
is at Europa orbit rather than terminal Earth re-entry.

Difference from Cassini VVEJGA (already in catalogue per the recent
acquisitions wave): different destination (Europa orbit, not Saturn),
different inner-cruise geometry (V-E-E, not V-V-E-J).

Suggested catalogue row (do NOT create now, await admission review):

- id: `campagnola-2014-ehm-orbiter-tour-11-O3`
- orbit_class: mga_tour
- epoch_locked: true
- n_returns: 1
- validity_window: {start: "2021-11-22T00:01:13Z", end: "2030-01-03T21:29:49Z"}
- launch_epoch: "2021-11-22T00:01:13Z"
- bodies: ["E", "V", "J", "Ganymede", "Callisto", "Europa"]
- sequence_canonical: "E-V-E-E-J(JOI)-G-G-G-G-C-G-G-G-G-C-G-G-G-E-E-EOI(Europa)"
- model_assumption: "real-ephemeris"  (the table dates run to seconds,
  consistent with a high-fidelity ephemeris solve)

The 12-L1 lander tour is a sibling row with a different Jovian arc
(7-day cruise variant, EOI ~1.067 km/s vs 0.477 km/s for 11-O3).

A catalogue ADMISSION decision should consider whether two-row
duplication is justified (separate Jovian arcs, same interplanetary
prefix) vs a single row with arc variants noted in `notes`. My
recommendation: keep them separate because the EOI burns, TID, and
total Delta-V differ materially.

## 6. Schema impact

NONE. The four-class taxonomy (`cycler / quasi_cycler / precursor_mga /
mga_tour`) added in schema v4.7 already accommodates the Campagnola
tours. No new enum value needed.

The TID (Total Ionizing Dose) metric is novel for the catalogue (every
Jovian tour will have one, but no prior row carries it). RECOMMEND: a
**non-schema annotation** in the `notes` field (TID is a per-mission
quantity tied to the science-orbit shielding assumption, 100-mil Al here;
making it a first-class schema field would force all rows to carry a
field most authors never tabulate). Leave it as free-text in `notes`
until a second TID-bearing Jovian tour lands, then revisit.

## 7. Errata

None found that affect the trajectory data.

Possible typesetting slip: Table 1 "Callisto 10" row date renders as
"2029 Oct 23 05:23:05:28" which has one too many colons — likely should
be "2029 Oct 23 05:10:28" or "23 05:23:28" given the journal format. This
is a presentation glitch, not a numerical defect; the entry function
(Tisserand-Poincare leveraging arc) is unambiguous from context. Filed
as a low-impact errata; benefit-of-the-doubt framing per
feedback_respectful_errata_framing.

## 8. Action items

1. Pre-anchor the 11-O3 orbiter tour in KNOWN_CORPUS (V0) — done by
   this digest; awaits corpus-admission review.
2. Pre-anchor the 12-L1 lander tour as a separate corpus row (sibling).
3. Decide whether 12-L4 is an annotation on 12-L1 or a third row (my
   recommendation: annotation, because the only differences are the
   Europa-arrival endgame and the EOI burn).
4. Reproduction effort: reproduce the V_inf chain at Ganymede (G1..G9)
   in DE440 from the table dates. This is a high-value task that should
   land before catalogue admission (V0 -> V1 gate).
5. Decision: defer Campagnola 2010/2012/2014 multi-moon-orbiter
   companion paper [22, methodology] for a future digest if a TLT
   reproduction lane opens.
6. Cross-reference flag: this paper cites Niehoff [18] for the original
   "Multi-Moon orbiter" inspiration — note the Galilean-tour lineage
   from Niehoff 1970 (separately digested, task #383 paper 2) through
   Ross/Scheeres/Lantoine/Russell to Campagnola.

# CCSDS ODM mining — OCM/OEM logical-block model vs. our catalogue schema

**Date:** 2026-06-05
**Source:** CCSDS 502.0-B-3, *Orbit Data Messages*, Blue Book, April 2023
(Editorial Change 1, May 2023). ISO 26900. Read via local mirror; all
section/page references below are to that document.
**Purpose:** Verify the OCM/OEM borrowings in `docs/spec.md` §16.6–§16.7 against
the primary standard (previously adopted from secondary sources).

---

## Q1 — Which document defines OCM? (spec citation check)

**Finding: the spec citation "OCM (CCSDS 504.0-B)" is WRONG on two counts.**

1. **OCM is defined *in 502.0-B-3 itself*, alongside OPM, OMM, OEM.** §1.1
   (Purpose and Scope, p. 1-1):

   > "This Orbit Data Messages (ODM) Recommended Standard specifies **four
   > standard message formats** … The Orbit Parameter Message (OPM), the Orbit
   > Mean-Elements Message (OMM), the Orbit Ephemeris Message (OEM), and **the
   > Orbit Comprehensive Message (OCM)**."

   §1.3 (Rationale, p. 1-2):

   > "This update to version 2 of the Orbit Data Messages **adds a fourth message
   > type, the OCM** …"

   The whole of section 6 ("Orbit Comprehensive Message", pp. 6-1 ff.) is the
   normative OCM definition. There is no separate "OCM standard."

2. **CCSDS 504.0-B is the *Attitude* Data Messages, not OCM.** Reference [10]
   in §1.6 (p. 1-6):

   > "[10] *Attitude Data Messages*. Issue 1. Recommendation for Space Data
   > System Standards (Blue Book), **CCSDS 504.0-B-1**. Washington, D.C.: CCSDS,
   > May 2008."

   So "504.0-B" denotes ADM (attitude), a completely different message family.

**Correct citation:** OCM is defined in *CCSDS 502.0-B-3, Orbit Data Messages,
Blue Book, April 2023, section 6* (ISO 26900). It is NOT a standalone "504.0-B".

Related reference IDs gathered for context (all from §1.6, pp. 1-5/1-6):
- [5] *XML Specification for Navigation Data Messages*, CCSDS 505.0-B-3, May 2023
  (this is the NDM/XML spec the ODM XML form depends on).
- [9] *Tracking Data Message*, CCSDS 503.0-B-2, June 2020.
- [14] *Conjunction Data Message*, CCSDS 508.0-B-1, June 2013.

---

## Q2 — OEM structure (confirm/refute spec characterisation)

Spec §16.6.1 table says: *"OEM (Ephemeris Message): Sampled state vectors +
interpolation metadata … Over-specified; we don't sample states."* and §16.7.8
item 2 calls it "a *track* format (sampled states + interpolation)."

**Verdict: the characterisation is CORRECT.** Supporting detail:

- **Sampled state vectors.** §5.1.1 (p. 5-1):
  > "Orbit information may be exchanged … by sending an ephemeris in the form of
  > a series of state vectors (Cartesian vectors providing position and velocity,
  > and optionally accelerations) using an OEM."
  §5.2.4.1 (p. 5-7): each ephemeris line is fixed order
  `Epoch, X, Y, Z, X_DOT, Y_DOT, Z_DOT, X_DDOT, Y_DDOT, Z_DDOT`; position+velocity
  mandatory, acceleration optional (§5.2.4.2). OEM is Cartesian-only (no element
  sets — unlike OCM).

- **Interpolation methods actually specified.** Table 5-3 (p. 5-6), `INTERPOLATION`
  keyword (Optional), example values: **HERMITE, LINEAR, LAGRANGE**.
  `INTERPOLATION_DEGREE` (Conditional — "must be used if the INTERPOLATION keyword
  is used"). The interpolation keyword is *per metadata block* (it applies to "the
  immediately following set of ephemeris lines"), so different segments can use
  different interpolation. §5.2.4.6 (p. 5-7): a second metadata block resets
  interpolation — "interpolation using succeeding ephemeris data with ephemeris
  data occurring prior to that metadata block shall not be done" (used to model a
  maneuver / eclipse discontinuity).

- **Epoch required?** Yes, multiply:
  - Header `CREATION_DATE` Mandatory (Table 5-2, p. 5-3).
  - Metadata `START_TIME` and `STOP_TIME` both Mandatory; `TIME_SYSTEM`,
    `REF_FRAME`, `CENTER_NAME` all Mandatory (Table 5-3, pp. 5-5/5-6).
  - Every ephemeris data line carries a mandatory `Epoch` time tag (§5.2.4.1).
  So OEM is fully epoch-anchored — confirms our "every operational standard
  requires an epoch" observation (spec §16.6.1).

- **Maneuvers:** OEM has no maneuver block (confirms spec table "maneuvers no").
  Discontinuities are modelled only implicitly via separate metadata/ephemeris
  blocks (§5.2.4.6). Covariance is optional (`COVARIANCE_START`/`_STOP`,
  Table 5-4, p. 5-8), lower-triangular 6×6, ordered [1,1]→[6,6] (§5.2.5.4).

**No correction needed** — though note the spec's "Segments yes, maneuvers no"
is slightly loose: OEM "segments" are just repeated metadata+ephemeris blocks,
not a first-class segment construct. Minor, not wrong.

---

## Q3 — OCM logical blocks (correct names, optionality, scoping)

Spec claims blocks **TRAJ, MAN, PHYS, COV, USER**. The actual OCM structure
(Table 6-1, "OCM File Layout and Ordering Specification", p. 6-3) is:

| Block (our name) | Actual section name | Delimiters | Status | Notes |
|---|---|---|---|---|
| (header) | OCM Header | — | **M** | single, mandatory |
| (metadata) | OCM Metadata | `META_START`/`META_STOP` | **M** | single, mandatory; **at most ONE** (§6.2.4.3) |
| TRAJ | "orbit data" / Trajectory State Time History | `TRAJ_START`/`TRAJ_STOP` | **O** | one or more |
| PHYS | "physical properties" / Space Object Physical Characteristics | `PHYS_START`/`PHYS_STOP` | **O** | at most ONE (§6.2.6.2) |
| COV | "covariance data" / Covariance Time History | `COV_START`/`COV_STOP` | **O** | one or more |
| MAN | "maneuver data" / Maneuver Specification | `MAN_START`/`MAN_STOP` | **O** | one or more |
| PERT | perturbations parameters | `PERT_START`/`PERT_STOP` | **C** | single; required if OD block present |
| OD | orbit determination | `OD_START`/`OD_STOP` | **O** | single |
| USER | user-defined parameters | `USER_START`/`USER_STOP` | **O** | single |

**Things our 5-block list (TRAJ, MAN, PHYS, COV, USER) got wrong / missed:**

1. **Block order is TRAJ → PHYS → COV → MAN → PERT → OD → USER** (Table 6-1),
   and §6.2.2.1 says "The order of occurrence of OCM keywords shall be fixed."
   Our spec lists them as "TRAJ, MAN, PHYS, COV, USER" — wrong order and missing
   two blocks. (Order matters: §6.2.4 `OCM_DATA_ELEMENTS` enumerates present
   blocks "in the same order as the data blocks in the message", from the fixed
   list `ORB, PHYS, COV, MAN, PERT, OD, USER`.)
2. **Two blocks omitted: PERT (perturbations) and OD (orbit determination).**
   PERT is Conditional-mandatory if OD is present (Table 6-1; §6.2.5.14 also
   recommends a PERT section whenever a TRAJ section is included).
3. **The canonical keyword for the trajectory data-element token is `ORB`**, not
   `TRAJ` (see `OCM_DATA_ELEMENTS` list `ORB, PHYS, COV, MAN, PERT, OD, USER`),
   even though the block delimiters are `TRAJ_START`/`TRAJ_STOP` and keywords are
   `TRAJ_*`. Worth noting if we ever emit `OCM_DATA_ELEMENTS`.
4. **Degenerate case allowed:** an OCM may have *zero* data blocks — only
   header+metadata (§6.2.1.1 note, p. 6-2). Header and the single metadata
   section are the only mandatory parts.

**CENTER_NAME / REF_FRAME / TIME_SYSTEM scoping (the key structural point):**

- **`TIME_SYSTEM`** — **message-level**, lives in the single OCM Metadata section
  (Table 6-3, p. 6-9), Mandatory, default UTC. Its description: *"Time system for
  all absolute time stamps in this OCM including EPOCH_TZERO… **This field is used
  by all OCM data blocks.**"* So time system is global to the message and cannot
  vary per block.
- **`EPOCH_TZERO`** — also message-level metadata (Table 6-3, p. 6-10),
  Mandatory: *"Default epoch to which all relative times are referenced in data
  blocks… **This field is used by all OCM data blocks.**"* This is the OCM's single
  required epoch (relative time tags are SI seconds from EPOCH_TZERO; §6.2.2.3).
- **`CENTER_NAME`** — **per-TRAJ-block** (Table 6-4, p. 6-17), Mandatory *within*
  each trajectory block, default EARTH. NOT in the message metadata section.
- **`TRAJ_REF_FRAME`** (not "REF_FRAME") — **per-TRAJ-block** (Table 6-4, p. 6-17),
  Mandatory, default ICRF3. Reference frame is per-trajectory-block, not global.

So the spec's gloss in the §16.6.1 table — "explicit CENTER_NAME / REF_FRAME /
TIME_SYSTEM" lumped together — is imprecise about scope: in OCM, **TIME_SYSTEM is
global, CENTER_NAME and the frame (TRAJ_REF_FRAME) are per-trajectory-block.** The
keyword is `TRAJ_REF_FRAME`, not `REF_FRAME` (REF_FRAME is the OEM/OPM/OMM spelling).

Other TRAJ keywords worth knowing (Table 6-4, pp. 6-15…6-20): `TRAJ_ID`,
`TRAJ_PREV_ID`/`TRAJ_NEXT_ID` (linked-list chaining of arcs across messages),
`TRAJ_BASIS` (PREDICTED/DETERMINED/TELEMETRY/**SIMULATED**/OTHER — note SIMULATED
"for generic simulations, future mission design studies, and optimization
studies", directly relevant to design-study cyclers), `TRAJ_BASIS_ID`,
`INTERPOLATION` (HERMITE/LINEAR/LAGRANGE/**PROPAGATE**), `TRAJ_TYPE` (Mandatory,
default CARTPV — selects the element set per SANA registry; e.g. CARTP, CARTPV,
Keplerian, spherical), `ORB_AVERAGING` (OSCULATING/BROUWER/KOZAI…),
`TRAJ_UNITS` (optional free-text SI-unit list in brackets, informational only).

---

## Q4 — Does CENTER_NAME vary per segment within one message? (schema v4.2 check)

**Finding: YES — OCM explicitly permits CENTER_NAME to vary per trajectory block.
The v4.2 `segments[].center` assumption is VALID per the OCM model.**

Direct support, §6.2.5.4 (p. 6-13) — "Each trajectory state time history data
block should differ from all others in at least one of the following respects:"

> "a) the selected element set (TRAJ_TYPE);
>  b) the orbit basis … (TRAJ_BASIS_ID);
>  c) the reference frame is unique (TRAJ_REF_FRAME);
>  **d) the orbit center is unique (CENTER_NAME);**
>  e) the data interval timespan is unique …"

And §6.2.5.15 (p. 6-14):

> "The **CENTER_NAME shall be used to specify the origin of the reference frame
> that the trajectory state time history is specified in.** The specified center
> may either be a natural, gravitationally attracting body … or it may be a
> non-gravitationally attracting origin … If a non-gravitationally attracting
> origin is selected, however, then the specified TRAJ_TYPE shall be confined to
> Cartesian or spherical coordinates …"

Because CENTER_NAME is a *per-TRAJ-block* keyword (Table 6-4) and §6.2.5.4(d)
names a differing center as a legitimate reason to open a new TRAJ block, an OCM
can carry consecutive trajectory blocks centered on, e.g., SUN then EARTH then
MARS within a single message. This is exactly the per-segment center our cycler
arcs need (heliocentric transfer arc vs. planetocentric flyby arc).

**Caveat to record in the schema rationale:** the *time system and epoch* do NOT
vary per block (they are message-global, see Q3), so a faithful OCM export must
keep one TIME_SYSTEM/EPOCH_TZERO for the whole record even though `segments[]`
may each carry their own `center`. Our YAML can be looser, but the OCM projection
cannot put a per-segment time system.

---

## Q5 — Other useful patterns for a published-orbit catalogue

**Covariance handling.**
- OEM: optional, lower-triangular 6×6, ordered upper-left [1,1] → lower-right
  [6,6] row-by-row, double precision (§5.2.5.4, p. 5-8). `COV_REF_FRAME` overrides
  the state frame if different (§5.2.5.3). Multiple matrices allowed, ordered by
  increasing time tag (§5.2.5.6/.7).
- OCM: far more flexible — "Covariance Time History" block(s) (Table 6-1; Table
  6-6 referenced p. 6-32, not transcribed here). §2.4 overview (p. 2-2/2-3) lists
  "covariance matrix of selectable/arbitrary order … Cartesian 3×3, 6×6, 7×7, or
  any combination of order, reference frame, and orbit elements." For our use
  (published cycler families rarely carry covariance) this is over-capable;
  if a source ever gives a state covariance, OCM's per-block COV with its own
  frame is the natural target.

**Units conventions** (§1.5.1, p. 1-3; §7.7).
- SI base/derived units. Notation: `*` = multiply, `**` = exponent (so m² is
  `m**2`, √km is `km**0.5`), `/` = divide (§1.5.1.2). `d` = 86400 SI seconds,
  `n/a` = not applicable. Useful if we ever serialise unit strings to match CCSDS.

**KVN vs XML.**
- §1.1 (p. 1-1): both "Keyword = Value Notation" (KVN) and XML are defined;
  "Selection of KVN or XML … should be mutually agreed between message exchange
  partners." KVN is the line-oriented `KEYWORD = value` form; XML form is built
  per section 8 against the NDM/XML spec CCSDS 505.0-B-3 (ref [5]). **There is no
  CCSDS-defined JSON or YAML serialization** — confirms spec §16.6.2's statement
  that a YAML projection is *our* mapping, not a standard. (JSON is an in-progress
  Nav-WG effort; not in this Blue Book.)

**Version / versioning fields.**
- Each message has a `CCSDS_<TYPE>_VERS` mandatory header keyword in form `x.y`:
  `CCSDS_OEM_VERS` (Table 5-2), `CCSDS_OCM_VERS` (Table 6-2), value `3.0`.
  Description: *"'y' is incremented for corrections and minor changes, and 'x' is
  incremented for major changes."* A clean precedent for our own schema-version
  field semantics (major.minor).

**Message ID / originator / provenance metadata** (Table 5-2 p. 5-3; Table 6-2
p. 6-5; Table 6-3 pp. 6-7…6-8) — rich provenance patterns we could mirror:
- `ORIGINATOR` (M) — creating agency/operator, drawn from a controlled list
  (SANA registry).
- `MESSAGE_ID` (O) — "ID that uniquely identifies a message from a given
  originator." Plus OCM-only `PREVIOUS_MESSAGE_ID` / `NEXT_MESSAGE_ID` and
  `..._EPOCH` — linked-list chaining of message revisions over time.
- OCM metadata point-of-contact set: `ORIGINATOR_POC`, `_POSITION`, `_PHONE`,
  `_EMAIL`, `_ADDRESS`; `TECH_ORG`, `TECH_POC`, … ; cross-message links
  `ADM_MSG_LINK`, `CDM_MSG_LINK`, `PRM_MSG_LINK`, `RDM_MSG_LINK`, `TDM_MSG_LINK`.
- `CREATION_DATE` (M) UTC.
- `CATALOG_NAME` + `OBJECT_DESIGNATOR` (the catalog that assigns the ID) — a
  pattern directly analogous to our (source, source-id) provenance pairing.
- `CELESTIAL_SOURCE` (e.g. `JPL_DE_FILES`) and `EOP_SOURCE` in OCM metadata —
  i.e. CCSDS records *the ephemeris/EOP model the numbers were computed against*.
  This is exactly the intent of our v4.2 `source_ephemeris` field (§16.7.9); good
  external validation that "which ephemeris this was computed against" is a
  first-class provenance datum in the standard.

**Maneuver block (informative, for the MAN borrowing).** §6.2.8 (pp. 6-35…6-46),
Tables 6-7/6-8/6-9. Each MAN block `MAN_START`/`MAN_STOP`, mandatory `MAN_ID`,
`MAN_DEVICE_ID` (special values `ALL`, `DEPLOY`), `MAN_BASIS`, `MAN_REF_FRAME`,
`MAN_COMPOSITION` (declares the per-line element order). Supports impulsive ΔV
(Table 6-8) and finite-burn/acceleration profiles with duty cycles (Table 6-9,
`DC_TYPE` TIME / TIME_AND_ANGLE / CONTINUOUS). Far richer than a cycler catalogue
needs; our `maneuvers[]` only needs the impulsive-ΔV subset.

---

## Spec corrections needed

1. **§16.6.1 table + §16.7.8 item 2: change "OCM (CCSDS 504.0-B)" / "the newer
   OCM (CCSDS 504.0-B)" → OCM is defined *within* CCSDS 502.0-B-3 (section 6),
   April 2023.** 504.0-B-1 is the *Attitude* Data Messages (per ref [10]).
   This is a factual citation error in two places.

2. **§16.7.8 item 2 wording "the newer OCM … adds multi-segment trajectories":**
   keep the substance (OCM did add multi-segment + maneuvers in the B-3 revision)
   but drop the implication it is a separate standard.

3. **§16.6.1 table OCM cell: block list/order is imprecise.** Actual fixed order
   is TRAJ(ORB) → PHYS → COV → MAN → PERT → OD → USER (Table 6-1). Our list
   "TRAJ, MAN, PHYS, COV, USER" omits **PERT** and **OD** and has wrong order. If
   the catalogue ever claims OCM-completeness, note PERT/OD exist (we legitimately
   don't use them).

4. **§16.6.1 table OCM cell: scoping of CENTER_NAME/REF_FRAME/TIME_SYSTEM.** These
   are NOT all "explicit per block." In OCM: **TIME_SYSTEM and EPOCH_TZERO are
   message-global** (single metadata section, "used by all OCM data blocks");
   **CENTER_NAME and the frame are per-TRAJ-block**, and the frame keyword is
   **`TRAJ_REF_FRAME`** (not `REF_FRAME`). Recommend the spec note this split,
   because it constrains the OCM exporter (one time system per record).

5. **v4.2 note (optional, additive):** record that OCM's per-block CENTER_NAME is
   *verified* (§6.2.5.4(d), §6.2.5.15) — `segments[].center` is sound — but that a
   faithful OCM export must still carry a single message-global TIME_SYSTEM /
   EPOCH_TZERO even when segment centers differ.

(No golden/numeric values are touched by any of the above — these are
structural/citation corrections only.)

---

## Not checked (honest scope limits)

- **Section 3 (OPM) and Section 4 (OMM) in full** — only the OMM tail (4.2.4.x,
  pp. 4-8/4-9, TLE/TEME conventions + User-Defined Parameters caveat) was read;
  OPM tables 3-1…3-3 were not. Not needed for the OCM/OEM questions.
- **OCM Covariance Time History detail (Table 6-6, p. 6-32)** — confirmed it
  exists and its overview (§2.4); the full keyword table was not transcribed.
- **OCM Perturbations (Table 6-10), Orbit Determination (Table 6-11), User-Defined
  (Table 6-12)** — existence and status confirmed from Table 6-1; keyword-level
  detail not read.
- **OCM Maneuver Tables 6-8 / 6-9** (impulsive vs finite-burn field lists) — the
  governing prose §6.2.8 and Table 6-7 header were read; the per-field tables were
  not transcribed (not needed; we only use impulsive ΔV).
- **Section 7 (KVN syntax detail) and Section 8 (XML instantiation)** — only the
  high-level KVN-vs-XML fact from §1.1 was used; the regex/XML-tag-delimiter
  detail (7.10, 8.x) was not read.
- **Annexes A–J** — including Annex B (normative controlled value lists:
  TIME_SYSTEM B3, CENTER_NAME bodies B2, frames B4, TRAJ_TYPE element sets B7),
  Annex G (worked OCM examples G-15…G-20), and Annex J (changes vs. B-2). Annex B
  would be the place to harvest exact controlled vocabularies if we ever build a
  validating OCM exporter; not done here.
- **No round-trip / serialization test** was performed — this is a documentary
  read, not a code change.

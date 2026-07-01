# Task #510: NAIF public SPK kernel archive survey vs catalogue gaps

**Date:** 2026-07-01
**Task:** survey NAIF's public satellite-SPK archive (`https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/`)
and map it against what the catalogue currently lacks. RESEARCH ONLY — no code,
no catalogue edits, no downloads, no commits.

---

## 0. What's already in-hand

| Kernel | Location | Size | Used by |
|---|---|---|---|
| `jup365.bsp` | `~/dev/references/kernels/jup365.bsp` | 1.1 GB | Galilean cycler real-eph (#480), Liang et al. Member D SPICE row |
| `ura111.bsp` | `~/GMAT/R2022a/data/planetary_ephem/spk/uranian/` | 162 MB | #331/#332/#335 Uranian V3/V4-strict SILVER gauntlet |
| DE440 (planets) | astropy cache | ~114 MB | `spice_kernels.py` cross-check, all V3+ rows |
| `naif0012.tls` (LSK) | astropy cache | tiny | ET/UTC conversion |
| DE405/DE421/DE424 | GMAT install | — | GMAT-native, unused for satellite work |

`~/dev/references/kernels/` contains only `jup365.bsp`; the Uranian kernel lives
under the GMAT install tree by a different convention (`scripts/install_uranian_spice.sh`).
No Saturn, Neptune, or Pluto satellite kernel exists anywhere in the environment.

## 1. Catalogue gap scan

Grepped `data/catalogue.yaml` for `validation_level`, `bodies:`, and
real-ephemeris/SPICE mentions:

- **Pluto-Charon**: exactly one row, `ross-rt-pc-cycler-32-2026` (V2-ballistic,
  CR3BP-periodic). #506 (2026-07-01) already characterized this row as
  V3-gated on Gate (a) "no Pluto-Charon SPK kernel in-repo" — this survey's
  §2.5 resolves that gate.
- **Saturn/Titan-Enceladus**: 22 `russell-strange-2009-titenc-*` rows + 1
  family-seed row, ALL at V0. These are capped at V0 by the **#76 Tier-1/Tier-2
  split** (patched-conic ideal-model; CR3BP Tier-2 reproduction was deferred,
  not attempted) — **not** by kernel availability. A Saturn satellite kernel
  would only become relevant once/if a Tier-2 CR3BP or real-eph pass is
  actually run on one of these rows; right now there is no gated row to unblock.
- **Uranus**: `heaton-longuski-2003-uranian-tour-u00-01` (V0, patched-conic) +
  the #327/#331/#332/#335 SILVER lineage (V4, already using URA111). No
  additional Uranus row is currently kernel-gated.
- **Neptune/Triton**: zero catalogue rows reference Neptune or Triton at all
  (only physical parameters in `core/satellites.py` for screening). No benefit
  from any Neptune kernel today.
- **Binary asteroids (Didymos-Dimorphos)**: `#308` (asteroid-leveraging search)
  references NEAs generically (Eros, Bennu, Itokawa, Didymos as individual
  bodies for flyby leveraging), not binary-asteroid mutual dynamics. No
  catalogue row models a binary-asteroid system as a two-primary cycler
  substrate. No benefit today.

## 2. NAIF archive survey

Source: `https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/`
(+ `a_old_versions/` subdirectory) and `aa_summaries.txt` (per-file body list
+ time-coverage index) fetched 2026-07-01. All filenames/sizes/coverage below
are read directly off the live NAIF listing, not inferred.

### 2.1 Jupiter/Galilean — already have, current is fine

| Kernel | Size | Bodies | Coverage | Note |
|---|---|---|---|---|
| **jup365.bsp** (in-repo) | 1.1 GB | Io, Europa, Ganymede, Callisto + Amalthea/Thebe/Adrastea/Metis | 1600–2200 | current best for Galilean cyclers; long span, in use |
| jup347.bsp | 879 MB | Jupiter + 75 moons (incl. irregulars) | 1799–2200 | superset of bodies, shorter span — no benefit unless irregular Jovian moons enter scope |
| jup348.bsp / jup349.bsp | 57–93 MB | small subsets of newly-numbered moons | 1799–2199 | irrelevant to cyclers |

**No re-acquisition warranted.** jup365 already covers the only bodies the
catalogue's Galilean/Liang-et-al rows use, over a longer span than the
alternatives.

### 2.2 Saturn/Titan + icy moons

| Kernel | Size | Bodies | Coverage | Priority |
|---|---|---|---|---|
| **sat441.bsp** | 631 MB | Mimas, Enceladus, Tethys, Dione, Rhea, **Titan**, Hyperion, **Iapetus**, Phoebe, Helene, Telesto, Calypso, Methone, Polydeuces + Saturn | 1749–2250 | **medium** — this is the file `core/satellites.py` already cites (`SAT441`) for Iapetus/Hyperion GM; it covers all 7 classical moons the Russell-Strange Titan-Enceladus rows need |
| sat415.bsp | 595 MB | ring-region minor moons (Janus, Epimetheus, Atlas, Prometheus, Pandora, Pan, Methone, Pallene, Anthe, Aegaeon) | 1949–2050 | not relevant to Titan-Enceladus rows |
| sat455–459.bsp | 68–278 MB each | small/irregular outer moons (Ymir, Paaliaq, Tarvos, S/YYYY_S_## series) | 1749/1750–2150 | not relevant |
| sat441xl_part-1/2.bsp | 1.9 GB each | extreme-precision extended-length version of sat441 | huge span | not warranted (same overkill pattern as ura111xl) |

sat441 is mass-consistent with the values `satellites.py` already uses
(SAT441-sourced GM for Iapetus/Hyperion) — no model-mismatch risk if adopted.

**What it would unblock:** nothing *today*. The 22 `russell-strange-2009-titenc-*`
rows are V0 by the #76 Tier-1/Tier-2 gap (no CR3BP reproduction attempted yet),
not by kernel absence. `sat441.bsp` would become directly useful the moment a
Tier-2 CR3BP or real-eph pass is scoped for any Titan-Enceladus row — worth
pre-fetching only if that work is imminent.

### 2.3 Uranus — already have URA111; a newer comprehensive kernel exists but isn't needed yet

| Kernel | Size | Bodies | Coverage | Priority |
|---|---|---|---|---|
| **ura111.bsp** (in-repo, GMAT) | 162 MB | Ariel, Umbriel, Titania, Oberon, Miranda, Uranus | 1900–2099 | in use; #335's V4-strict SILVER gauntlet already runs on it |
| ura184_part-3.bsp | 369 MB | same 5 classical moons + Uranus + 8 irregular moons + S/2023_U_1 (newly discovered) | 1600–2399 | **low** — extends span 199 yr → 799 yr and adds the 2023-discovered moon; current #335 gauntlet only tested launch epochs 2000/2030/2050, all inside URA111's 1900–2099 window |
| ura115.bsp | 1.4 GB | (not itemized; large single-file classical+irregular kernel) | — | not warranted, same overkill pattern |
| ura111xl-\*.bsp (6 files) | ~2 GB each | per-body extreme-length versions | 12001 BC – 17000 AD | explicitly rejected already (`install_uranian_spice.sh` comment) — navigation-precision files, not needed |

**Delta if upgraded:** ura184_part-3 would only matter if a future Uranian row
needs an epoch outside 1900–2099 (e.g. the #335 note flags "URA111 SPICE
kernel coverage ends 2099" as a known boundary — if a `validity_window` needs
to extend past 2099, this is the file to fetch). Not recommended now — no
row currently needs it.

### 2.4 Neptune/Triton — no catalogue benefit

| Kernel | Size | Bodies | Coverage | Priority |
|---|---|---|---|---|
| nep097.bsp | 100 MB | Triton, Neptune | 1600–2399 | **none** — zero catalogue rows reference Neptune/Triton |
| nep104.bsp | 318 MB | Halimede, Psamathe, Sao, Laomedeia, Neso + 2 newly-numbered moons | 1600–2399 | none |
| nep095.bsp / nep101xl\*.bsp | 1.5–2.0 GB | extended/XL variants | huge | none, not warranted regardless |

No action. Recorded for completeness only — do not manufacture urgency here.

### 2.5 Pluto system (Pluto-Charon + Nix/Hydra/Kerberos/Styx) — the one that matters

| Kernel | Size | Bodies (NAIF IDs) | Coverage | Priority |
|---|---|---|---|---|
| **plu060.bsp** | **129 MB** | Charon (901), Nix (902), Hydra (903), Kerberos (904), Styx (905), Pluto (999), Earth, Sun | **1800–2199** | **HIGH — directly unblocks #506** |

**This is a correction to #506's estimate.** #506 (2026-07-01, same day) named
`plu058.bsp (~3.8 MB)` as the file to fetch. That filename/size is stale —
the *current* NAIF Pluto-system kernel is **plu060.bsp at 129 MB** (the
larger size reflects real multi-body dynamics for all 5 moons, not just
Charon, vs. the older 2-body-only plu058 solution). `plu058.bsp` does not
appear in the current top-level listing; it would need to be pulled from
`a_old_versions/` if ever wanted for a historical-reproducibility reason, but
there is no reason to prefer it over plu060.

**Mass/GM consistency check against `core/satellites.py`:** the satellites.py
comment for the Pluto system already reads `# JPL SSD satellite physical
(phys_par, ref PLU060) + mean elements` — i.e. the GM values already baked
into the catalogue's Charon/Nix/Hydra screening (GM_Charon=106.1,
GM_Nix=0.0015, GM_Hydra=0.0020 km³/s²) are **already sourced from PLU060**.
Fetching `plu060.bsp` is therefore not just compatible but literally
same-source as the mass data already in the repo — zero model-mismatch risk,
unlike the eccentricity-driven CR3BP/real-eph mismatch #506 flagged for Gate
(b) (which is a physics issue independent of which kernel file is used).

**What it unblocks:** `plu060.bsp` resolves #506's Gate (a) ("no PLU SPK
in-repo or GMAT install") completely. It does **not** resolve Gate (b) — #506
already established that V3-as-naive-real-eph-propagation is the wrong test
for this CR3BP-periodic-orbit row regardless of kernel availability (the
Charon eccentricity model-mismatch argument is orthogonal to whether the
kernel exists). Installing plu060.bsp would only become useful if a future
task pursues #506's stated harder alternative: differentially correcting the
CR3BP IC into the real-eph model to find the real-eph analog orbit, which
needs the kernel as an input regardless of methodology.

### 2.6 Binary asteroids / comet systems — lower priority, no current gap

- **Didymos-Dimorphos**: NAIF hosts a DART-mission-specific SPICE archive
  (`pub/naif/pds/pds4/dart/dart_spice/`) with Didymos-system ephemeris SPKs,
  but these are narrow-window mission-reconstruction kernels centered on the
  2022-09-26 DART impact epoch, not a multi-decade generic binary-asteroid
  ephemeris analogous to the planetary-satellite series. No catalogue row
  currently models Didymos-Dimorphos as a two-primary binary-cycler substrate
  (only individual-NEA flyby leveraging under #308). No acquisition
  recommended; flagged only in case #308/#494's binary-cycler genome scope
  ever extends to a real binary-asteroid pair.
- No other binary-NEA SPK products were found in the generic_kernels archive
  (binary asteroids are not systematically covered the way planetary
  satellites are — most published binary-asteroid ephemerides exist as
  Keplerian mean elements in the literature, not SPK kernels).

## 3. Download URL pattern

```
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/<filename>
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/<filename>  (superseded versions)
```
Per-file summaries (body lists + coverage intervals) for every file in the
directory are in `aa_summaries.txt` in the same directory — worth fetching
once and keeping as a reference rather than re-querying per file.

## 4. Priority acquisition list (top of the list, in order)

| # | Kernel | Size | Unblocks | Priority |
|---|---|---|---|---|
| 1 | **plu060.bsp** | 129 MB | #506 Gate (a): Pluto-Charon V2→V3 infra gate — clears the "no kernel" blocker outright; mass-consistent with existing `satellites.py` PLU060 citations | **HIGH** |
| 2 | sat441.bsp | 631 MB | No row gated *today*; would unblock a future Tier-2 CR3BP/real-eph pass on any of the 22 Russell-Strange Titan-Enceladus V0 rows once #76 Tier-2 work is scoped | MEDIUM (pre-fetch only if that work is imminent) |
| 3 | ura184_part-3.bsp | 369 MB | No row gated today; only useful if a future Uranian `validity_window` needs to extend past URA111's 2099 boundary | LOW |
| 4 | nep097.bsp | 100 MB | No catalogue row references Neptune/Triton at all | LOW / parking only |

**Recommendation: fetch only #1 (plu060.bsp, 129 MB) now.** It is the sole
kernel in this survey that resolves a concretely identified, already-documented
catalogue gate (#506). Items 2–4 are legitimate future needs but have no
gated row waiting on them today — fetching them now would be acquiring
infrastructure ahead of demand, which the project record (#506, #335) shows
tends to sit unused for weeks until the consuming task actually fires.

## 5. Cumulative size (feeds #509 CI-cache planning)

**Currently in environment:** jup365.bsp (1.1 GB, `~/dev/references/kernels/`)
+ ura111.bsp (162 MB, GMAT tree) + DE440 (~114 MB, astropy cache) ≈ **1.4 GB**.

**If the full "one current/best kernel per system" set surveyed here were
fetched** (plu060 + sat441 + ura184_part-3 + nep097, on top of what's already
in-hand):

| Item | Size |
|---|---|
| jup365.bsp (have) | 1.1 GB |
| ura111.bsp (have) | 162 MB |
| plu060.bsp (new) | 129 MB |
| sat441.bsp (new) | 631 MB |
| ura184_part-3.bsp (new, optional upgrade) | 369 MB |
| nep097.bsp (new, no current benefit) | 100 MB |
| **Total** | **≈ 2.5 GB** |

**If only the recommended top-1 (plu060.bsp) is added:** ≈ 1.5 GB cumulative
(vs. today's 1.4 GB) — a 129 MB delta.

**Do not extrapolate to the "everything in the archive" figure** — the XL/
extended-length per-body kernels (`ura111xl-*`, `sat441xl_part-*`,
`nep097xl-*`, `jup310xl` etc.) run 1–2 GB *each* and sum to tens of GB
archive-wide; none of them are warranted for cycler-catalogue work (the
project has already explicitly rejected this tier once, in
`install_uranian_spice.sh`'s comment on why URA111 was chosen over
ura116xl/ura111xl). For #509's CI-cache-limit purposes, the realistic
working set is the ≈2.5 GB "one kernel per system, current version" table
above, not the full archive.

## 6. Discipline anchors

- `feedback_check_dont_guess` — every filename/size/coverage claim in this
  doc was read directly off the live NAIF listing and `aa_summaries.txt`,
  not inferred from memory or the #506 doc's (stale) plu058 reference.
- `feedback_ground_citations_against_content` — the plu058→plu060 correction
  in §2.5 is exactly this discipline: don't trust an inherited filename/size
  without checking the live source.
- `project_negative_results_registry` — §2.4 (Neptune) and §2.6 (binary
  asteroids) are stated as plain "no current benefit" rather than manufactured
  urgency, per the task's explicit instruction.

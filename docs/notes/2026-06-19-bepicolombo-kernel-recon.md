# BepiColombo SPK kernel recon for `mga_tour` admission (#345 / #399)

**Date:** 2026-06-19
**Purpose:** Reconnaissance for admitting a BepiColombo `mga_tour` row. Identify the
reconstructed NAIF SPK kernel, confirm the spacecraft NAIF ID empirically, and
extract the planetocentric hyperbolic-excess velocity (V∞) at each of the 9
gravity-assist flybys via `cyclerfinder.verify.mission_spk.vinf_at_flyby`, with a
published closest-approach-altitude cross-check (sourced-only).

This is a recon deliverable only — the catalogue row itself is written by the
parent/concurrent agent, not here.

---

## Recommended kernel

- **Filename:** `bc_mpo_fcp_00226_20181020_20270407_v01.bsp`
- **Base URL:** `http://spiftp.esac.esa.int/data/SPICE/BEPICOLOMBO/kernels/spk/`
- **Size:** 61,204,480 bytes (~61 MB) — full-fidelity, NOT a coarse summary kernel.
- **Confirmed NAIF spacecraft ID:** **-121** (BepiColombo MPO — the Mercury
  Planetary Orbiter; during cruise it carries the MTM + MMO composite). Confirmed
  empirically: `spkobj` returns exactly `[-121]`.
- **Coverage (from `spkcov` + LSK):** single continuous interval
  **2018-10-20T02:12:20 .. 2027-04-07T00:00:00** (TDB→UTC). Spans launch and all 9
  flybys with margin — one kernel suffices, no per-flyby kernel swapping needed.

`fcp` = flight-dynamics control product; the `00226` delivery is the latest
reconstructed full-cruise solution at the time of recon. The coarse
`bc_mmo_scp_cruise_20181016_20251205_v02.bsp` (7 KB, body -68) is NOT usable
(spkezr errors); avoided as instructed. The per-swingby `bc_mpo_fcp_*SwingbyMTP*`
kernels are predicted planning products and were not used.

### How it was selected
- Directory listing (`curl … | grep bsp`): 306+ `bc_mpo_cog_*`, 267 `bc_mpo_fcp_*`,
  62 `bc_mpo_fmt_*`, plus MMO/MTM kernels.
- The `cog` (centre-of-gravity) per-delivery kernels are small (58 KB span 8 yr →
  too coarse to resolve a periapsis); the `fmt` kernels are ~520 KB (also coarse).
- The numbered `fcp` cruise kernel is 61 MB over the same span → real fidelity.
  Empirically it resolves every periapsis to <1.1% of the published CA altitude.

---

## The 9 flybys — extracted vs published

Extractor call (identical for all): `vinf_at_flyby(spk, -121, <body>, <epoch>,
mission="BepiColombo", window_minutes=2880.0, n_samples=41)`. Periapsis refinement
on (default). Body GM/radius sourced from `cyclerfinder.core.constants.PLANETS`.

| Flyby | Body | Sourced epoch (UTC, CA) | Extracted V∞ (km/s) | vis-viva std (km/s) | Extracted CA alt (km) | Published CA alt (km) | %diff | CA offset (min) |
|-------|------|-------------------------|---------------------|---------------------|-----------------------|-----------------------|-------|-----------------|
| Earth     | Earth   | 2020-04-10 04:25      | 3.9932 | 0.00337 | 12685.9 | 12677  | +0.07 | −0.04 |
| Venus-1   | Venus   | 2020-10-15 03:58      | 7.9791 | 0.00707 | 10721.9 | 10720  | +0.02 | +0.53 |
| Venus-2   | Venus   | 2021-08-10 13:51      | 8.0221 | 0.00710 |   552.6 |   552  | +0.10 | +0.90 |
| Mercury-1 | Mercury | 2021-10-01 23:34:41   | 6.6586 | 0.01046 |   198.8 |   199  | −0.10 | +0.01 |
| Mercury-2 | Mercury | 2022-06-23 09:44      | 6.4231 | 0.02647 |   198.1 |   200  | −0.97 | +0.37 |
| Mercury-3 | Mercury | 2023-06-19 19:34      | 3.6401 | 0.01664 |   234.8 |   236  | −0.51 | +0.42 |
| Mercury-4 | Mercury | 2024-09-04 21:48      | 2.7483 | 0.01431 |   163.3 |   165  | −1.01 | +0.82 |
| Mercury-5 | Mercury | 2024-12-01 14:23      | 2.6121 | 0.01288 | 37625.6 | 37626  | −0.00 | +0.69 |
| Mercury-6 | Mercury | 2025-01-08 05:59      | 1.8921 | 0.01019 |   295.0 |   295  | +0.00 | −0.14 |

All 9 periapses refined (`closest_approach_refined=True` for every row). The
vis-viva window std is ≤0.027 km/s everywhere (sub-1% of V∞ → converged
asymptotic V∞). The refined CA epoch sits within ~1 minute of the published time
in all cases, confirming the sourced epochs are accurate to the minute.

**Max |%diff| vs published CA altitude: 1.01% (Mercury-4).** All others ≤1%.

### Notes on individual flybys
- **Mercury-5 (2024-12-01)** is the post-anomaly high flyby: a solar-electric
  thruster power shortfall (discovered Sept 2024) raised the targeted altitude from
  the planned ~37 km to **37,626 km**. It is a genuine gravity-assist flyby (not the
  terminal capture), correctly extracted; the high altitude is physical, not an
  extraction defect. The extractor reproduces it to −0.00%.
- **No flyby failed.** Every row converged cleanly on the first attempt with the
  default `window_minutes=2880` (48 h) window; no epoch correction or window
  widening was needed.
- The cosmos.esa.int flyby index lists some flybys one calendar day later (e.g.
  "Mercury 4: September 5th 2024", "Mercury 6: January 9 2025") — this is a local-
  timezone (CEST) vs UTC display difference, not a discrepancy. ESA's own
  press release for Mercury-4 states 21:48 UTC / 23:48 CEST, matching the UTC
  epoch used here.

---

## Launch and Mercury Orbit Insertion

- **Launch:** 2018-10-20, **01:45 UTC**, on an **Ariane 5 ECA** (flight VA245),
  from Kourou. (Wikipedia, citing ESA.)
- **Mercury Orbit Insertion (MOI):** originally December 2025, **delayed to
  November 2026** following the September 2024 thruster-power anomaly. MOI is the
  TERMINAL CAPTURE, not a clean flyby — no V∞ extracted for it, as instructed.

---

## Sources (published values — sourced-only discipline)

- ESA / Wikipedia consolidated flyby table (dates, CA UTC times, CA altitudes;
  launch date+rocket; MOI date): <https://en.wikipedia.org/wiki/BepiColombo>
- ESA flyby portal (flyby roster + dates):
  <https://www.cosmos.esa.int/web/bepicolombo-flyby>
- ESA press release, fourth Mercury flyby (21:48 UTC / 23:48 CEST, 4 Sep 2024,
  ~165 km): <https://www.esa.int/Science_Exploration/Space_Science/BepiColombo/Fourth_Mercury_flyby_begins_BepiColombo_s_new_trajectory>
- ESA, first Venus flyby (15 Oct 2020, ~10,720 km):
  <https://www.esa.int/Science_Exploration/Space_Science/BepiColombo/BepiColombo_flies_by_Venus_en_route_to_Mercury>
- ESA, ready for double Venus flyby (Venus-2 context, 10 Aug 2021, ~550 km):
  <https://www.esa.int/Science_Exploration/Space_Science/ESA_gets_ready_for_double_Venus_flyby>
- ESA, third Mercury flyby (19 Jun 2023, ~236 km):
  <https://www.esa.int/Science_Exploration/Space_Science/BepiColombo/BepiColombo_braces_for_third_Mercury_flyby>
- Kernel archive (SPK directory):
  <http://spiftp.esac.esa.int/data/SPICE/BEPICOLOMBO/kernels/spk/>

Extracted V∞ and CA values are DERIVED by our own code from the reconstructed SPK
at the sourced epochs — they are not published and are not auto-admittable.

---

## Ready to admit? — verdict

**Yes.** A single reconstructed kernel
(`bc_mpo_fcp_00226_20181020_20270407_v01.bsp`, NAIF -121) cleanly resolves all 9
gravity-assist flybys: every periapsis refines, every vis-viva V∞ converges
(std ≤0.027 km/s), and every extracted closest-approach altitude matches the
ESA/published value to within 1.01% (max), confirming both the kernel fidelity and
the sourced epochs. The flyby sequence (1×Earth, 2×Venus, 6×Mercury) with V∞
ranging 8.02 km/s (Venus-2) down to 1.89 km/s (Mercury-6, the final pre-capture
flyby) is a textbook `mga_tour`. Recommend admission with the kernel + NAIF ID +
extracted V∞ tuple above, citing the sources listed; the anomalous Mercury-5 high
flyby should carry a one-line note that its 37,626 km altitude is the physical
post-thruster-anomaly value, not an extraction artefact.

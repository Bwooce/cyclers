# Task #514: NAIF kernel-freshness checker

**Date:** 2026-07-01
**Task:** stop discovering superseded NAIF SPK kernels by accident. #506
(same day) assumed the Pluto-Charon satellite kernel was `plu058.bsp`;
#510's survey found that filename no longer exists in NAIF's public archive
at all -- the current kernel is `plu060.bsp` (129 MB vs. #506's stale
`~3.8 MB` estimate). That was a one-off manual catch during a survey. This
task makes the check mechanical and repeatable.

## What it checks

`scripts/check_kernel_freshness.py` holds an extensible manifest
(`KERNEL_MANIFEST`) of every NAIF satellite SPK kernel the project currently
relies on, per `docs/notes/2026-07-01-510-naif-spk-kernel-survey.md`:

| System | Pinned filename | Local path convention |
|---|---|---|
| jupiter | `jup365.bsp` | `~/dev/references/kernels/jup365.bsp` |
| uranus | `ura111.bsp` | `~/GMAT/R2022a/data/planetary_ephem/spk/uranian/ura111.bsp` |
| pluto | `plu060.bsp` | `~/dev/references/kernels/plu060.bsp` |

Saturn (`sat441.bsp`) and Neptune (`nep097.bsp`) are not yet tracked --
#510 established no catalogue row is currently gated on either, so there is
nothing pinned to go stale. Adding a system when one becomes needed is a
one-entry addition to the dict; the parsing/comparison logic is generic over
`KernelSpec`, not per-body.

For each tracked system the script fetches NAIF's live directory listing
(`https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/`, an
Apache `mod_autoindex` HTML page) via `urllib.request` -- the same HTTP
mechanism `src/cyclerfinder/verify/spice_kernels.py` already uses for the
leapseconds kernel; no new HTTP dependency was added (the checker has zero
third-party imports, stdlib only). It parses every `.bsp` row into
`{filename: size_bytes}`, matches candidates in the SAME numbered series as
the pinned filename (e.g. `jup3\d{3}\.bsp`, `plu0\d{3}\.bsp`,
`ura\d{3}(_part-N)?\.bsp`), and picks the highest version number as
"current upstream". The series patterns are anchored so they never match
the XL/extended-length navigation kernels (`ura111xl-701.bsp`,
`ura116xl.bsp`, etc.) -- those are a different kernel class the project has
already explicitly rejected (see `scripts/install_uranian_spice.sh`'s
comment on why URA111 was chosen over the XL releases).

Two independent axes are reported per kernel:

* **freshness** -- pinned-filename version vs. highest version found in the
  live listing. This is a pure manifest-vs-remote comparison; it needs NO
  local file on disk, so it works unmodified in CI (which never holds the
  multi-GB binaries -- they are never committed). This is the ONLY thing
  that drives a nonzero exit code / a filed GitHub issue.
* **local presence** -- is the pinned kernel actually sitting at its
  documented local path on the machine running the check? Purely
  informational (mirrors `ensure_jup365_kernel()`'s local-path-only,
  no-auto-download-by-default convention). Absence alone never fails the
  exit code.

Combining the two gives three reported statuses per kernel:

* `OK` -- pinned filename is current upstream and present locally.
* `STALE` -- a newer kernel exists upstream in the same series (regardless
  of local presence). **Drives exit code 1 and a filed/updated GitHub
  issue.**
* `MISSING` -- pinned filename is current upstream but not found locally.
  Exit code stays 0 (nothing to fix upstream-side; the developer just
  hasn't fetched it).

A network failure while fetching the NAIF listing is reported separately
and gives exit code 2 (distinguished from staleness so a transient NAIF
outage never gets misreported as "kernel superseded").

## Running it manually

```sh
uv run python scripts/check_kernel_freshness.py            # all tracked systems
uv run python scripts/check_kernel_freshness.py --system pluto
uv run python scripts/check_kernel_freshness.py --fetch    # download any STALE kernel's
                                                             # current upstream version
```

`--fetch` is opt-in only and never runs unattended (mirrors the
auto-download-refusal convention already established for `jup365.bsp` and
the leapseconds kernel): it downloads the new file into the SAME local
directory the stale pin lives in, under the new upstream filename -- it
never overwrites the old file and never edits `KERNEL_MANIFEST` for you.
Verified live against the real NAIF archive while building this (2026-07-01):

```
jupiter: OK -- jup365.bsp is current and present locally
pluto: OK -- plu060.bsp is current and present locally
uranus: STALE -- pinned ura111.bsp (v111), current upstream is
  ura184_part-1.bsp (v184, 2040109466 bytes, delta +1870180762 bytes)
```

This exactly reproduces #510's own finding (`ura184_part-3.bsp` "newer
alternative but not urgent") -- the checker correctly flags Uranus as
STALE by design; whether to actually re-acquire ura184 in place of ura111
remains a human judgement call (per #510, URA111's 1900-2099 coverage
already spans every epoch any current Uranian row needs), which is exactly
why `--fetch` requires an explicit human decision rather than
auto-upgrading.

## Test coverage

`tests/scripts/test_check_kernel_freshness.py` -- no network access; all
HTTP is mocked. The fixture HTML is a verbatim excerpt of NAIF's real
directory listing (captured 2026-07-01 while drafting the checker), so the
parsing regex is exercised against the real row format. Covers: `.bsp` row
parsing, non-`.bsp` companion-file exclusion (`.cmt`/`.mrg`), XL-variant
exclusion, the full OK/STALE/MISSING matrix, the #506/#510 "STALE
regardless of local presence" scenario specifically, and `run()`'s
exit-code contract including a zero-local-files CI simulation.

```
uv run pytest tests/ -k kernel_freshness -q
# 18 passed
uv run ruff check scripts/check_kernel_freshness.py tests/scripts/test_check_kernel_freshness.py
# All checks passed!
```

## Scheduled workflow

`.github/workflows/kernel-freshness.yml` -- monthly cron (`17 4 1 * *`,
04:17 UTC on the 1st of each month) plus `workflow_dispatch` for manual
runs. No precedent scheduled workflow existed in this repo (`.github/workflows/`
only had `ci.yml`, push/PR-triggered); the schedule + `workflow_dispatch` +
concurrency-group shape follows the pattern already established in the
`cyclers.space` companion repo's `refresh-windows.yml` (weekly cron +
manual dispatch), adapted for this repo's purpose: this workflow never
commits or pushes anything (no kernel binary is ever committed), so on
staleness it files or updates a GitHub issue via `gh issue create`/`gh issue
comment` instead of `refresh-windows.yml`'s commit-and-push pattern. The
job installs nothing beyond `actions/setup-python` (the checker has zero
third-party dependencies) and needs no kernel files present -- it runs the
same network-only manifest-vs-remote comparison described above.

## Discipline anchors

- `feedback_check_dont_guess` -- the manifest's series-comparison logic was
  verified against the LIVE NAIF listing (curl + regex dry-run) before being
  encoded, not inferred from the #510 survey doc's summary tables alone.
- `feedback_ground_citations_against_content` -- this task exists precisely
  because #506 inherited a stale filename without checking the live source;
  the checker is the mechanical version of that discipline.
- `feedback_golden_tests_sourced_only` -- the test fixture is real NAIF
  listing HTML (not invented), so the parsing test is checking OUR code
  against a REAL row shape, not a shape we made up to make the test pass.

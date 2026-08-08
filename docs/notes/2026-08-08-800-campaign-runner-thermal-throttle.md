# #800 — thermal/duty-cycle throttle for `campaign_runner.py`

**Trigger**: this Mac now also runs the self-hosted GitHub Actions CI runner
(set up the same day, CI minutes exhausted for August 2026). Any of the
registered multi-week campaigns (`#789`-`#792`) dispatched with
`CampaignRunnerConfig`'s default `n_workers=-1` would peg all 8 logical cores
(4P+4E) continuously for days-to-weeks on a laptop, stacked on top of
whatever CI jobs land — a real sustained-heat concern the user flagged
directly ("we need to somehow detune the multiweek campaigns to not run this
laptop too hot").

## What was added

`src/cyclerfinder/search/campaign_runner.py`:

- `CampaignRunnerConfig.pause_seconds_per_batch: float = 0.0` — sleeps this
  long between batches (skipped after the final batch). A sensor-independent
  duty-cycle knob: even with `n_workers` already capped below
  `os.cpu_count()`, a deliberate breather between checkpoint batches reduces
  average sustained duty cycle over a long unattended run.
- `CampaignRunnerConfig.thermal_backoff_seconds: float = 0.0` — if > 0, polls
  a new `_os_thermal_throttled()` helper between batches and sleeps this long
  extra whenever macOS itself reports it is already limiting CPU speed for
  heat. Opt-in (default 0.0), so it's a no-op on Linux/other platforms
  without `pmset`.
- `_os_thermal_throttled() -> bool`: runs `pmset -g therm`, parses the
  `CPU_Speed_Limit`/`CPU_Scheduler_Limit` percentage fields (< 100 = actively
  throttling). Fails open (`False`) on any error — missing binary,
  non-macOS, unparseable output, timeout — so it can never hang a campaign.
  Checked live on this machine: `pmset -g therm` exists and returns cleanly
  (currently "no thermal warning level recorded" — the normal idle state);
  this is a reactive signal (fires once macOS is *already* throttling), not
  a proactive temperature reading — Apple Silicon exposes no unprivileged
  live temperature sensor. The proactive half of the mitigation is
  `pause_seconds_per_batch` and, per-dispatch, capping `n_workers` well below
  8 (recommend 4-5 for a laptop running unattended for days).
- `_sleep(seconds)`: a one-line indirection so tests can stub the campaign's
  own pauses without patching process-wide `time.sleep` — discovered during
  testing that doing so also intercepts loky's internal polling sleep calls
  (245k+ calls recorded in one test run), making the executor spin instead
  of pausing.

`tests/search/test_campaign_runner.py`: 8 new tests (18 total, all green) —
pause-between-not-after-last-batch, thermal-backoff-fires-when-throttled,
no-backoff-when-not-throttled, and 5 parametrized/edge cases for
`_os_thermal_throttled`'s `pmset` output parsing (nominal, throttled on
either field, garbage output, missing binary).

## Not done (deliberately out of scope here)

- No change to `CampaignRunnerConfig.n_workers`'s own default (`-1` = all
  cores) — that default is fine for the existing fast toy-worker tests and
  short campaigns; the mitigation for real long campaigns is an explicit
  lower `n_workers` passed at dispatch time (a per-campaign call-site
  decision, not a library-wide default change), combined with these two new
  pause knobs.
- No proactive temperature-based throttle (no unprivileged live sensor
  exists on Apple Silicon without `sudo powermetrics`) — the reactive
  `pmset` check plus the unconditional duty-cycle pause is the practical
  ceiling without asking for elevated privileges.

## Recommended defaults for real campaign dispatches (`#789`-`#792`)

`CampaignRunnerConfig(n_workers=4, pause_seconds_per_batch=30.0,
thermal_backoff_seconds=120.0, ...)` — leaves 4 of 8 logical cores free for
the CI runner + interactive use, takes a 30s breather every
`checkpoint_batch_size` cells, and backs off an extra 2 minutes if macOS
itself ever reports active thermal limiting.

## Verification

`uv run pytest tests/search/test_campaign_runner.py -q` — 18/18 pass.
`uv run ruff check .` / `ruff format --check .` clean on both changed files.
`uv run mypy src tests` — clean (839 files; one `attr-defined` error from a
test monkeypatching `campaign_runner.subprocess` instead of the already-
imported `subprocess` module directly, fixed).

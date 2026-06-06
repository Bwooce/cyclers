# cyclerfinder CLI quickstart

`cyclerfinder` is a thin, typed driver over the existing enumerator / optimiser
/ discover / gauntlet stack. It does no new physics: every subcommand is a shell
over an already-shipped function. The default body set is the **M8 Venus-Earth-
Mars (VEM) anchor** at the 6.4-yr beat (spec §8 line 152), so the headline VEM
campaign runs with zero required flags.

## Install

```sh
uv sync                 # base CLI (no plotting)
uv sync --extra viz     # adds matplotlib for the `viz` subcommand
```

The console script `cyclerfinder` is installed by `[project.scripts]`. You can
also run it as `python -m cyclerfinder.cli`.

## Subcommands

```text
cyclerfinder [--version] [-v|--verbose] <subcommand> ...
```

### enumerate — list cells + Tisserand feasibility

```sh
cyclerfinder enumerate --feasible-only            # VEM defaults
cyclerfinder enumerate --bodies V,E,M --l-max 6 --k-max 3 --format json
```

`--period` accepts `beat` (natural beat, `period_basis=None`) or an anchor pair
like `E-M` (sets `Cell.period_basis`). Output: `--format table|json|csv`.

### solve — optimise a single cell

```sh
# Build the cell from flags:
cyclerfinder solve --bodies E,M --sequence E-M-E --k 2 \
    --revs 0,0 --branch single,single --fidelity idealized --format json

# Or pass a cell id (round-trips through Cell.id):
cyclerfinder solve --cell-id 'EM|E-M-E|k2|r00|bss' --fidelity idealized
```

Optimiser outputs land under a `computed` key (closure residual, max V∞,
maintenance ΔV, taxi cost) — never presented as sourced anchors. Convergence is
carried in `constraints_satisfied`; a non-converged cell is honest data, not a
CLI failure (exit stays 0).

`--fidelity ephemeris --mode maintenance` runs the real-DE440 maintenance path
(slow). `--mode ballistic` is **exposed but M-ED-gated** (see below).

### discover — ledger-backed campaign with V-level gates

```sh
cyclerfinder discover --ledger out/vem.jsonl                  # VEM defaults
cyclerfinder discover --bodies E,M --k 2 --ledger out/em.jsonl --max-cells 50
```

Tallies cells searched / solved and per-V-level counts; the ledger is always
written. `--enable-v3` wires the ballistic V3 closure gate (**M-ED-gated**).

### report — campaign summary (markdown + JSON)

```sh
cyclerfinder report --ledger out/vem.jsonl --out out/vem --with-verdicts
```

Writes `<out>.md` and/or `<out>.json` (`--format md|json|both`). Each candidate
is partitioned into **disjoint** `sourced` (catalogue-traceable) and `computed`
(optimiser / ledger) blocks — the golden-discipline split is structural, not a
convention. `--with-verdicts` attaches the gauntlet `VerdictTier` (GOLD / SILVER
/ BRONZE / REJECTED) **verbatim**; the report never upgrades a tier.

### viz — visualisations (requires the `[viz]` extra)

```sh
cyclerfinder viz beat --bodies V,E,M --out out/beat.png
cyclerfinder viz porkchop --bodies E,M --tof-min 100 --tof-max 400 \
    --epoch-range 2032-01-01:2039-01-01 --out out/porkchop.png
cyclerfinder viz trajectory --cell-id 'EM|E-M-E|k2|r00|bss' --out out/traj.png
```

matplotlib is a lazy-imported optional extra (Agg backend, file output). Without
it the `viz` subcommand exits **3** with an install hint; the base CLI never
imports matplotlib.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | OK |
| 2 | usage error (bad flags / unknown body / no subcommand) |
| 3 | missing `[viz]` extra |
| 4 | not implemented (interim, no longer reachable) |
| 5 | no candidates / empty ledger |

## Honesty contract (sourced vs computed)

The CLI is a *driver*, not a discovery claim (spec §11.3, §17; M8-Core golden
discipline). `solve` and `report` keep **sourced** facts (catalogue period /
sequence / sourced V∞ multisets) strictly separate from **computed** results
(our optimiser's V∞ / ΔV / closure residual). A value our own code produced is
never presented as a sourced or golden anchor. The `report` artifact carries
each candidate's gauntlet verdict tier verbatim so the trust ladder (spec §14,
V0-V5) stays visible.

## Worked VEM campaign walkthrough

```sh
cyclerfinder enumerate --feasible-only
cyclerfinder discover --ledger out/vem.jsonl
cyclerfinder report --ledger out/vem.jsonl --out out/vem --with-verdicts
cyclerfinder viz beat --bodies V,E,M --out out/beat.png
```

## M-ED-gated paths

The ballistic corrector and the V3 ballistic-closure gate are deliverables of
the concurrent **M-ED** milestone (`docs/superpowers/plans/2026-06-05-m-ed-
ballistic-corrector.md`). `solve --mode ballistic` and `discover --enable-v3`
are wired and accept their flags, but their convergence is **exposed-but-inert**
until M-ED lands; the corresponding convergence smoke tests are skipped with a
recorded reason.

# Vendored neural-low-thrust-approximator pretrained weights

These CSV files are pretrained model weights for the Zhang–Acciarini–Izzo–Baoyin–Topputo
2026 neural low-thrust trajectory cost / reachability surrogate, used by
`cyclerfinder.search.neural_reach_prefilter` (#276) as a tier-0 prefilter.

## Source
- **Upstream:** https://github.com/zhong-zh15/neural-low-thrust-approximator
- **Paper:** Zhang, Acciarini, Izzo, Baoyin, Topputo (2026), "Pretrained Approximators
  for Low-Thrust Trajectory Cost and Reachability," arXiv:2605.26790.
- **Local copy:** `cyclers_pdf/papers/zhang-acciarini-2026-pretrained-approximators-lowthrust-cost-reachability-arxiv-2605.26790.pdf`

## Vendored content
- `eigen_model_large/` — ΔV surrogate (~2.1 MB; 10-layer MLP, CSV per layer).
- `eigen_model_tmin_large/` — minimum-ToF surrogate (~2.1 MB).
- `LICENSE.upstream` — the upstream Mozilla Public License 2.0 (MPL-2.0) text.

## License
The vendored weights are distributed under **MPL-2.0** (verbatim, per
`LICENSE.upstream`). MPL-2.0 is a file-level copyleft: modifications to these
files retain MPL-2.0; the surrounding cyclerfinder code is unaffected.

If you redistribute cyclerfinder with these weights, include `LICENSE.upstream`
and this NOTICE. Modifications to the weights themselves must be MPL-2.0 and
must indicate the changes.

## Why vendored
The weights are 4.2 MB total, the upstream is MPL-2.0 (redistribution-permissible),
and vendoring removes the env-var / clone-upstream dance from every user of the
prefilter. Tests no longer skip when the upstream is absent.

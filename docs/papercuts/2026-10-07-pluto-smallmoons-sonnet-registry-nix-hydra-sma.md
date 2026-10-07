# Registry semi-major axes for Nix and Hydra disagree with plu060 and the published period ratios

- Date: 2026-10-07
- Agent: pluto-smallmoons-sonnet
- Seen before: no

What happened: `src/cyclerfinder/core/satellites.py` gives Nix a = 49300 km and Hydra a = 65200 km (comment: JPL SSD mean elements, ref PLU060, accessed 2026-06-14). Through Kepler III with the registry Pluto-system GM (975.5) these give period ratios to Charon of 3.989 (Nix) and 6.067 (Hydra). The published ratios (and the lead's brief) are 3.89 and 5.98. Fitting plu060.bsp directly (2030-01-01 plus 3 yr, 3000 samples, J2000, relative to the Pluto-system barycentre, in-plane angle fit) gives mean periods Charon 6.38722 d, Styx 20.16195, Nix 24.85472, Kerberos 32.16798, Hydra 38.20202 d, so ratios 3.1566, 3.8913, 5.0363, 5.9810. The Kepler a from those periods is Nix 48481 km (registry +1.7 percent) and Hydra 64569 km (registry +1.0 percent); the kernel mean barycentric radius is 48689 and 64718 km. Charon (19600 vs 19596) agrees. Styx and Kerberos are not in the registry at all.
Workaround: #998 builds its ideal model from the kernel-fitted periods and records the registry values beside them.
Suggested fix: verify the Pluto small-moon registry entries against plu060 and the published elements (the lead is registering a task); add Styx and Kerberos. Check whether the #320 Pluto rows (which used the registry sma) change if corrected.

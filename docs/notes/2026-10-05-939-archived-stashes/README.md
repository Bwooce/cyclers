# Archived stash branches (#939, 2026-10-05)

Two stash commits had been pushed as branches from the other machine ("amdnuc return"). With the owner's
approval, the branches were deleted on 2026-10-05 after their diffs were saved here.

- `cyclers-amdnuc-return-wip-main-b44c34ed.patch`: `origin/amdnuc-return/wip-main` in Bwooce/cyclers,
  stash of uncommitted work on 2026-10-04 on top of 1d3d7165 (#884). Touches
  `scripts/screen_884_sun_forced_em_cyclers.py` and `scripts/screen_886_titan_rhea_torus_check.py`.
  Main has since rewritten both scripts (b1777003, 0c11834f); `stage_posthoc`/`run_pool` from this
  stash do not exist on main.
- `cyclers.space-amdnuc-return-stash-0-1fc5628.patch`: `origin/amdnuc-return/stash-0` in
  Bwooce/cyclers.space, the 2026-06-20 "tmp-410" stash on top of 0aa44a7 (errata.yaml,
  hero-gallery.ts, global.css). Its untracked-files part (`src/lib/hero-legend.ts` and its test) is in
  `cyclers.space-amdnuc-return-stash-0-untracked-7d4a355.patch`.

Apply with `git apply --3way <patch>` if anything here is needed.

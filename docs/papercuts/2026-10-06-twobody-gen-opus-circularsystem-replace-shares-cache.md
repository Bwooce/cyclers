# dataclasses.replace on CircularSystem shares its flyby-body cache

- Date: 2026-10-06
- Agent: twobody-gen-opus
- Seen before: no

What happened: `CircularSystem` (src/cyclerfinder/search/two_working_body.py) caches `FlybyBody` records in a dataclass field `_fb`. `dataclasses.replace(circ, flyby_overrides=...)` copies the field by reference, so the new system returned the original body constants from the shared cache and ignored the override. A Callisto-GM sweep (#943, note 6.21) showed a gate ratio that did not change with GM; the instrument check caught it before any verdict.
Workaround: pass `_fb={}` to `replace`.
Suggested fix: declare `_fb` with `field(default_factory=dict, init=False, compare=False)` and create it in `__post_init__`, or key the cache on the overrides, so `replace` always starts with an empty cache.

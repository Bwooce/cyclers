"""#791 Stage 0 -- Tisserand/energy-graph pruning of the Galilean moon-tour sequence alphabet.

Per `#858`'s review (`docs/notes/2026-08-21-858-campaigns-789-790-791-review.md` Sec. 5.4):
"Stage 0 (hours): Tisserand/energy-graph pruning of the sequence alphabet itself -- pure algebra,
before any cell is evaluated. Report the surviving sequence count; if 10^3 sequences survive at
length <= 8, the pruning is too weak to proceed." This is that pass, and ONLY that pass -- no cell
(epoch/ToF/n_rev) is evaluated here, no `campaign_runner` dispatch, no n-body shooting. Stage 1
(patched-conic enumeration of survivors) is a separate, later task, gated on this report.

Alphabet: the four Galilean moons (Io, Europa, Ganymede, Callisto), matching every prior campaign
in this lane (`#318`/`#501`/`releg_moontour.py`). Saturnian moons need `joint_cell.py`'s
`primary="Jupiter"` hardcode addressed first (`#858` Sec. 5.3) -- out of scope here.

Sequence definition: closed tours (``seq[0] == seq[-1]``), adjacency-distinct (no immediate
same-moon repeat -- a flyby of X immediately followed by another flyby of X with no intervening
leg is degenerate, mirroring `sequence.enumerate_cells`'s own invariant), lengths 3-8 encounters
(`#791`'s own registration language; matches `scan_501_broadened_joint_search.py`'s tag convention
where e.g. "EGE" = 3 encounters / 2 legs).

Pruning gates, two of them (the 2026-08-22 run adds the second after the first proved too weak):

1. **Pairwise** (original Stage 0): for every leg, does `tisserand.linkable` find a shared V_inf
   in the project's own established Jovian probe band (`releg_moontour._VINF_PROBE_KMS`) at which
   the two moons' constant-V_inf contours intersect -- each leg free to use its OWN V_inf?
2. **Chained** (energy-consistency, this run's addition): does a SINGLE common V_inf in the same
   probe band link EVERY leg of the sequence simultaneously? This is the physically correct
   ballistic constraint, not an arbitrary tightening: a flyby rotates V_inf direction but preserves
   its MAGNITUDE (`releg_moontour.py`'s own docstring, "Continuity model" section, already encodes
   this for the powered lane), so an unpowered patched-conic tour holds V_inf quasi-constant across
   the whole cycle -- a sequence whose legs only individually link at DIFFERENT V_inf values could
   never actually be flown ballistically end-to-end.

Both are deliberately narrower than `moon_prune.moon_leg_admissible` (which also folds in a VILM
Delta-v floor and a bend check) -- this stays the pure ENERGETIC (Tisserand) screen Fable's Stage 0
spec asked for, no ΔV-budget assumption baked in this early (a budget belongs to Stage 1/2's actual
cell evaluation).

The 6 sequences `#501` already swept and stamped to `empty_regions.jsonl`
(`scan_501_broadened_joint_search.py`) are excluded from the "new" survivor count by
capability-subsumption (skip, not re-sweep) but are still individually reported.

Output: ``data/found/791_stage0_sequence_pruning/results.jsonl`` (one row per enumerated closed
sequence) + a summary printed to stdout / captured in the companion results note.
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from cyclerfinder.core.satellites import PRIMARIES  # noqa: E402
from cyclerfinder.search.moon_prune import _a_range_au  # noqa: E402
from cyclerfinder.search.releg_moontour import _VINF_PROBE_KMS  # noqa: E402
from cyclerfinder.search.tisserand import linkable  # noqa: E402

MOONS: tuple[str, ...] = ("Io", "Europa", "Ganymede", "Callisto")
PRIMARY = "Jupiter"
L_MIN, L_MAX = 3, 8

# #501's own 6 already-stamped sequences (scan_501_broadened_joint_search.py SEQUENCE_CONFIGS)
# plus #318's own earlier CGCEC smoke-test sequence (jovian.CGCEC) -- the registration's own
# "7 hand-picked encounter sequences" count (data/OUTSTANDING.md's #791 bullet), corrected here:
# the #501 driver module itself only lists 6; CGCEC is the 7th, found by checking #318's own
# smoke-test script rather than trusting the registration's summary.
ALREADY_STAMPED: frozenset[tuple[str, ...]] = frozenset(
    {
        ("Europa", "Ganymede", "Europa"),
        ("Ganymede", "Callisto", "Ganymede"),
        ("Europa", "Ganymede", "Callisto", "Europa"),
        ("Io", "Europa", "Io"),
        ("Io", "Europa", "Ganymede", "Io"),
        ("Europa", "Ganymede", "Callisto", "Ganymede", "Europa"),
        ("Callisto", "Ganymede", "Callisto", "Europa", "Callisto"),  # jovian.CGCEC, #318 smoke
    }
)


def enumerate_closed_sequences(
    moons: tuple[str, ...], l_min: int, l_max: int
) -> list[tuple[str, ...]]:
    """Closed tours (seq[0]==seq[-1]), adjacency-distinct, encounter-count in [l_min, l_max]."""
    out: list[tuple[str, ...]] = []
    for length in range(l_min, l_max + 1):
        free_len = length - 1
        for body in itertools.product(moons, repeat=free_len):
            if any(body[i] == body[i + 1] for i in range(free_len - 1)):
                continue
            if body[-1] == body[0]:
                continue  # closing leg would repeat the last free moon
            out.append((*body, body[0]))
    return out


def build_pair_linkable_table(moons: tuple[str, ...]) -> dict[tuple[str, str], frozenset[float]]:
    """For every ordered pair of distinct moons, the set of probe-band V_inf that link them."""
    mu = PRIMARIES[PRIMARY]
    table: dict[tuple[str, str], frozenset[float]] = {}
    for a in moons:
        for b in moons:
            if a == b:
                continue
            a_range = _a_range_au(a, b)
            table[(a, b)] = frozenset(
                float(v)
                for v in _VINF_PROBE_KMS
                if linkable(a, b, float(v), a_range_au=a_range, mu=mu)
            )
    return table


def sequence_tisserand_survives(
    seq: tuple[str, ...], table: dict[tuple[str, str], frozenset[float]]
) -> bool:
    """Pairwise screen: every leg individually linkable, legs free to use different V_inf."""
    return all(table[(seq[i], seq[i + 1])] for i in range(len(seq) - 1))


def sequence_chained_survives(
    seq: tuple[str, ...], table: dict[tuple[str, str], frozenset[float]]
) -> frozenset[float]:
    """Chained screen: the set of V_inf that link EVERY leg simultaneously (empty = fails)."""
    common: frozenset[float] | None = None
    for i in range(len(seq) - 1):
        vs = table[(seq[i], seq[i + 1])]
        common = vs if common is None else (common & vs)
        if not common:
            return frozenset()
    return common if common is not None else frozenset()


def main() -> int:
    out_path = REPO_ROOT / "data" / "found" / "791_stage0_sequence_pruning" / "results.jsonl"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    table = build_pair_linkable_table(MOONS)

    # Positive control: the known-closing IEGI family (Io-Europa-Ganymede-Io,
    # tests/search/test_moon_prune.py::test_prune_keeps_the_known_closing_ieg_family) must survive
    # BOTH screens -- closing is impossible if the legs are not even energetically linkable at a
    # single common V_inf, so a real closer failing here would mean a gate is broken.
    pc_seq = ("Io", "Europa", "Ganymede", "Io")
    pc_pairwise = sequence_tisserand_survives(pc_seq, table)
    pc_chained = sequence_chained_survives(pc_seq, table)
    if not pc_pairwise or not pc_chained:
        print(
            f"POSITIVE CONTROL FAILED: {pc_seq} pairwise={pc_pairwise} chained={bool(pc_chained)}"
        )
        return 1

    sequences = enumerate_closed_sequences(MOONS, L_MIN, L_MAX)
    n_total = len(sequences)
    n_stamped = sum(1 for s in sequences if s in ALREADY_STAMPED)

    rows = []
    n_pairwise_total = 0
    n_pairwise_new = 0
    n_chained_total = 0
    n_chained_new = 0
    for seq in sequences:
        pairwise = sequence_tisserand_survives(seq, table)
        chained_vs = sequence_chained_survives(seq, table)
        chained = bool(chained_vs)
        stamped = seq in ALREADY_STAMPED
        if pairwise:
            n_pairwise_total += 1
            if not stamped:
                n_pairwise_new += 1
        if chained:
            n_chained_total += 1
            if not stamped:
                n_chained_new += 1
        rows.append(
            {
                "sequence": list(seq),
                "length": len(seq),
                "n_legs": len(seq) - 1,
                "already_stamped": stamped,
                "tisserand_survives": pairwise,
                "chained_tisserand_survives": chained,
                "chained_common_vinf_kms": sorted(chained_vs),
            }
        )

    with out_path.open("w") as f:
        for row in rows:
            f.write(json.dumps(row) + "\n")

    print(f"positive_control ({'-'.join(pc_seq)}): pairwise=SURVIVES chained=SURVIVES")
    print(f"total enumerated closed sequences (length {L_MIN}-{L_MAX}): {n_total}")
    print(f"already-stamped (skipped from 'new' counts): {n_stamped}")
    print(f"[pairwise] surviving (total/new): {n_pairwise_total}/{n_pairwise_new}")
    print(f"[chained]  surviving (total/new): {n_chained_total}/{n_chained_new}")
    print(f"[chained]  pruned (energy-inconsistent): {n_total - n_chained_total}")
    print(f"wrote {len(rows)} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# ruff: noqa: E501  (catalogue YAML and evidence text carry long lines)
# STATUS (2026-10-12): READY. Window 4 lands this script alone (owner: gc4-19 held). Dry-run clean alone on 5a84d622 (ruff, mypy, tests/data, writeback/census/rediscovery tests exit 0) and with gc4-19 on ac4a90d7 and 3ccbf5e7.
"""Catalogue window 4, Earth-Venus part (#1025 rung-d note secs. 11.6-11.9; lead orders 2026-10-11/12).

Rebuilt 2026-10-12 from the committed record after the session scratchpad was wiped. It reproduces
twobody-ext-opus's last window-4 script (md5 4d125826, rebuilt by replaying its recorded edits) with
these changes: the repo root is a required argument, a re-run refuses with a message before any
edit is prepared, the timeout limit in the no-chain text is read from the tool log, the two
retried rows cite sec. 11.9 in their validation_level comment, and ruff runs through uv.

Usage: uv run python scripts/window4/apply_1025_ev_window4.py <repo-root>
Every edit asserts exactly one match of its anchor (catalogue edits: inside the row's own block);
files are written only after every edit has been prepared. No row is added or removed.

Edits:
- data/catalogue.yaml:
  - ev5-4, ev5-5, ev5-6, ev5-7, ev5-8, ev5-9: validation_level V1 -> V2 (comment), a V2-ballistic
    notes paragraph;
  - ev4-1, ev4-4, ev4-6, ev4-7: data_gaps path "validation_level.v2" (not judgeable, L = 38;
    the 10-cycle diagnostic);
  - ev5-4, ev5-5: data_gaps path "validation_level.v2" (the e0 no-chain and the e1 band failure).
- src/cyclerfinder/data/validate.py: _V1_1025_EV's last sentence; new _V2_1025_EV; the 6
  _LEVEL_EVIDENCE entries re-keyed V1 -> V2.
- tests/data/test_schema_v45_fields.py: the 6 map entries V1 -> V2; the V2-count assert N -> N+6
  (N read from the file, not typed).
All numbers in the row text come from data/1025_ev_realeph/v2_10cyc/judge_*.json.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path()  # set from argv[1] in main()
V2 = ["ev5-4", "ev5-5", "ev5-6", "ev5-7", "ev5-8", "ev5-9"]
K4 = ["ev4-1", "ev4-4", "ev4-6", "ev4-7"]
NOTE = "docs/notes/2026-10-11-1025-ev-realeph-rung-d.md"
DATA = "data/1025_ev_realeph/v2_10cyc"

edits: dict[Path, str] = {}


def rid(m: str) -> str:
    return f"earth-venus-two-working-body-cycler-{m}-2026"


def text(rel: str) -> str:
    p = ROOT / rel
    return edits.get(p, p.read_text())


def rep(rel: str, old: str, new: str) -> None:
    s = text(rel)
    n = s.count(old)
    if n != 1:
        sys.exit(f"ABORT: {rel}: anchor found {n} times (need 1):\n{old[:300]}")
    edits[ROOT / rel] = s.replace(old, new)


def rep_row(m: str, old: str, new: str) -> None:
    """Replace inside one catalogue row's block only."""
    s = text("data/catalogue.yaml")
    head = f"- id: {rid(m)}\n"
    if s.count(head) != 1:
        sys.exit(f"ABORT: row {rid(m)} found {s.count(head)} times")
    i = s.index(head)
    j = s.find("\n- id: ", i + 1)
    j = len(s) if j < 0 else j + 1
    blk = s[i:j]
    if blk.count(old) != 1:
        sys.exit(f"ABORT: {m}: anchor found {blk.count(old)} times in its row:\n{old[:300]}")
    edits[ROOT / "data/catalogue.yaml"] = s[:i] + blk.replace(old, new) + s[j:]


def judge(m: str) -> dict[str, Any]:
    j: dict[str, Any] = json.loads((ROOT / DATA / f"judge_{m}.json").read_text())
    return j


def chains(j: dict[str, Any]) -> list[tuple[dict[str, Any], str, dict[str, Any]]]:
    """(epoch record, route, judged record) for every reached route."""
    return [
        (e, t["route"], t["judged"]) for e in j["epochs"] for t in e["routes"] if t.get("reached")
    ]


def rng(v: list[float], fmt: str) -> str:
    lo, hi = min(v), max(v)
    return f"{lo:{fmt}}" if f"{lo:{fmt}}" == f"{hi:{fmt}}" else f"{lo:{fmt}}-{hi:{fmt}}"


def lit(s: str, indent: int, first: str = "", width: int = 99) -> str:
    """``s`` as implicitly concatenated string literals, each line <= width (no escapes needed)."""
    assert '"' not in s and "\\" not in s
    words = s.split(" ")
    lines: list[str] = []
    cur = ""
    pad = " " * indent
    for w in words:
        cand = (cur + " " + w) if cur else w
        room = width - indent - 2 - (len(first) if not lines else 0) - 1
        if len(cand) > room and cur:
            lines.append(cur + " ")
            cur = w
        else:
            cur = cand
    lines.append(cur)
    out = [f'{pad}{first if i == 0 else ""}"{ln}"' for i, ln in enumerate(lines)]
    return "\n".join(out) + "\n"


def srng(v: list[float], fmt: str) -> str:
    return f"{min(v):+{fmt}} to {max(v):+{fmt}}"


def already_applied() -> str | None:
    """A reason string if any part of this patch is already in the tree."""
    v = (ROOT / "src/cyclerfinder/data/validate.py").read_text()
    if "_V2_1025_EV" in v:
        return "_V2_1025_EV is already in src/cyclerfinder/data/validate.py"
    cat = (ROOT / "data/catalogue.yaml").read_text()
    for m in V2 + K4:
        head = f"- id: {rid(m)}\n"
        if cat.count(head) != 1:
            return None  # the anchors below report this
        i = cat.index(head)
        j = cat.find("\n- id: ", i + 1)
        blk = cat[i : (len(cat) if j < 0 else j)]
        if 'path: "validation_level.v2"' in blk:
            return f"{m} already has a validation_level.v2 data_gaps entry"
        if m in V2 and "\n  validation_level: V2 " in blk:
            return f"{m} is already V2 in data/catalogue.yaml"
    return None


def main() -> None:
    global ROOT
    if len(sys.argv) != 2:
        sys.exit("usage: uv run python scripts/window4/apply_1025_ev_window4.py <repo-root>")
    ROOT = Path(sys.argv[1]).resolve()
    if not (ROOT / "data" / "catalogue.yaml").is_file():
        sys.exit(f"ABORT: {ROOT} is not a cyclers checkout (no data/catalogue.yaml)")
    why = already_applied()
    if why is not None:
        sys.exit(
            f"REFUSED: the Earth-Venus window-4 patch is already applied ({why}); nothing written"
        )

    # --- the promoted rows ---------------------------------------------------------------------------
    v2_lines: dict[str, str] = {}
    for m in V2:
        j = judge(m)
        assert j["judgeable"] and j["verdict"] == "PASS" and j["L"] == 1, (m, j["verdict"])
        good = [(e, r, x) for e, r, x in chains(j) if x["pass"]]
        worst = [x["gate"]["worst_ratio"] for _, _, x in good]
        refly = [x["refly_max_km"] for _, _, x in good]
        junc = [x["junction_kms"] for _, _, x in good]
        dmax = [max(x["c4"]["Dprime_d"]) for _, _, x in good]
        dband = [x["c4"]["band_D_d"] for _, _, x in good]
        wmax = [max(x["c4"]["Wprime_kms"]) for _, _, x in good]
        wband = [x["c4"]["band_W_kms"] for _, _, x in good]
        sd = [x["c4"]["slope_D_d_per_cycle"] for _, _, x in good]
        sw = [x["c4"]["slope_W_kms_per_cycle"] for _, _, x in good]
        ep_pass = ", ".join(
            f"JD {e['epoch_jd']:.1f}" + ("" if r == "S1" else f" ({r})") for e, r, _ in good
        )
        whys = []
        for e in j["epochs"]:
            if e["pass"]:
                continue
            judged = [t["judged"] for t in e["routes"] if t.get("reached")]
            if judged:
                c4 = judged[0]["c4"]
                whys.append(
                    f"epoch {e['epoch']} fails the date band (D' {max(c4['Dprime_d']):.2f} / "
                    f"{c4['band_D_d']:.2f} d)"
                )
                continue
            timed = []
            for rt in ("S1", "S2", "S3"):
                log = ROOT / DATA / m / rt / f"e{e['epoch']}" / "tool.log"
                if log.exists() and log.read_text().startswith("TIMEOUT"):
                    mt = re.match(r"TIMEOUT after (\d+) s", log.read_text())
                    if mt is None:
                        sys.exit(f"ABORT: unrecognised TIMEOUT line in {log}")
                    timed.append(f"{rt} timed out at {mt.group(1)} s under load")
            whys.append(
                f"epoch {e['epoch']} has no 10-cycle chain"
                + (f" ({'; '.join(timed)}; not retried)" if timed else "")
            )
        ep_none = "; ".join(whys)
        n = j["passes"]
        orphan = (
            "\n      The JD 2462505.4 chain came from a run that overran the 420-s limit (sec. 11.6,"
            "\n      disclosure 2); without it the row is 3/5, still PASS."
            if m == "ev5-9"
            else "\n      The JD 2467752.5 chain is from the approved retry (sec. 11.9; S2 converged in"
            "\n      under a minute at low load; at 420 s under load it had timed out)."
            if m in ("ev5-4", "ev5-5")
            else ""
        )
        para = (
            "\n    V2-ballistic (#1025; spec §14 owner amendment 2026-10-08, same-phase drift wording of\n"
            f"    2026-10-11; {NOTE} secs. 11.1-11.6, 11.9):\n"
            f"    PASS at {n} of 5 Standish epochs, 10 continuous cycles each (80 yr; L = 1 at 0.02 rev;\n"
            "    first and last cycle excluded).\n"
            f"    - Passing epochs: {ep_pass}. Not passing: {ep_none}.{orphan}\n"
            f"    - Junction continuity <= {max(junc):.1e} km/s; gate worst {rng(worst, '.3f')}; "
            f"re-fly <= {max(refly):.1e} km.\n"
            f"    - Same-phase drift: dates {rng(dmax, '.2f')} d against {rng(dband, '.2f')}-d bands; "
            f"|V_inf|\n      {rng(wmax, '.3f')} km/s against {rng(wband, '.2f')}.\n"
            f"    - Slope diagnostics (handed to V3): dates {srng(sd, '.3f')} d per cycle; |V_inf|\n"
            f"      {srng(sw, '.4f')} km/s per cycle.\n"
            "    - Controls: Hollister-Menning 1970 orbit 1 passes and a broken chain fails (k = 2,\n"
            "      L = 5); no published k = 5 member exists to control L = 1.\n"
            "    - V2 is on Standish mean elements; the chains run to about 2110-2142, past the\n"
            "      elements' 1800-2050 fit span (a test of the multi-cycle structure). The DE440 standing\n"
            "      is unchanged (see data_gaps).\n"
        )
        rep_row(
            m,
            "\n    Literature: literature_check v3-A1 (#1025)",
            para + "\n    Literature: literature_check v3-A1 (#1025)",
        )
        rep_row(
            m,
            "  validation_level: V1   # spec §14 V1 (#1025 writeback): lamberthub Izzo+Gooding agreement "
            "and the #973 independent DOP853 re-fly; see validate._LEVEL_EVIDENCE. NOT V2.\n",
            "  validation_level: V2   # spec §14 V2-ballistic (owner amendment 2026-10-08; same-phase drift "
            f"wording 2026-10-11; #1025): {n} of 5 Standish epochs, 10 cycles (80 yr), L = 1; see "
            f"validate._LEVEL_EVIDENCE and {NOTE} "
            + ("secs. 11.6, 11.9.\n" if m in ("ev5-4", "ev5-5") else "sec. 11.6.\n"),
        )
        v2_lines[m] = (
            f"Row: {m} V2 {n}/5 (final wording; Standish, 10 cycles = 80 yr, {DATA}/); gate worst "
            f"{rng(worst, '.3f')}; re-fly <= {max(refly):.1e} km; dates {rng(dmax, '.2f')} d against "
            f"{rng(dband, '.2f')} d; slopes dates {srng(sd, '.3f')} d, |V_inf| {srng(sw, '.4f')} km/s "
            "per cycle. "
        )

    # --- data_gaps: k = 4 (not judgeable) and ev5-4/5 (open) -----------------------------------------
    gap_anchor = "  legs: []\n"
    for m in K4:
        j = judge(m)
        assert not j["judgeable"] and j["L"] == 38 and j["c123_epochs"] == 5, m
        worst = [x["gate"]["worst_ratio"] for _, _, x in chains(j)]
        note = (
            "V2-ballistic is NOT JUDGEABLE under the spec §14 wording: the Earth-Venus k = 4 geometric "
            "repeat is L = 38 at 0.02 rev (cyclerfinder.verify.geometric_repeat), so the same-phase "
            "test needs >= 76 cycles (about 486 yr). Diagnostic (not a V2 test): ballistic 10-cycle "
            "(64-yr) Standish chains close at 5 of 5 epochs with junction continuity < 1e-6 km/s, the "
            f"gate passing (worst {rng(worst, '.3f')}) and the DOP853 re-fly <= 9.5e-4 km; same-phase "
            "values left blank."
        )
        rep_row(
            m,
            gap_anchor,
            '  - path: "validation_level.v2"\n'
            '    kind: "not-applicable"\n'
            f'    note: "{note}"\n'
            f'    source_hint: "{NOTE} secs. 11.1, 11.7-11.8; {DATA}/judge_{m}.json"\n'
            '    todo_ref: "#1025"\n' + gap_anchor,
        )
    for m in ("ev5-4", "ev5-5"):
        j = judge(m)
        assert j["verdict"] == "PASS" and j["passes"] == 3, m
        assert [e["outcome"] for e in j["epochs"]][:2] == ["NO 10-CYCLE CHAIN", "FAIL"], m
        note = (
            "V2-ballistic PASS at 3 of 5 epochs (spec §14 final wording, L = 1, 10-cycle 80-yr Standish "
            "chains). Two epochs do not pass: at JD 2463078.5 no 10-cycle chain exists from the ladder "
            "(S1 does not converge; S2 and S3 are not reached, in 39 and 15 s, no timeout; sec. 11.9), "
            "and at JD 2465416.0 the chain fails only the same-phase date band (D' 4.78 against 2.73 d; "
            "gate 0.791). Whether a ballistic 10-cycle chain exists at JD 2463078.5 is unknown (a fold "
            "of the ladder is not a negative)."
        )
        rep_row(
            m,
            gap_anchor,
            '  - path: "validation_level.v2"\n'
            '    kind: "unknown"\n'
            f'    note: "{note}"\n'
            f'    source_hint: "{NOTE} secs. 11.6, 11.9; {DATA}/judge_{m}.json"\n'
            '    todo_ref: "#1025"\n' + gap_anchor,
        )

    # --- validate.py -------------------------------------------------------------------------------
    rep(
        "src/cyclerfinder/data/validate.py",
        '    "gate-passing) and the DE440 blend (data/1025_ev_realeph/). NOT V2: no multi-lap "\n'
        '    "bounded-drift campaign has been run. "\n'
        ")\n",
        lit(
            "gate-passing) and the DE440 blend (data/1025_ev_realeph/). V2 was judged separately "
            f"({NOTE} sec. 11); a row that stays V1 states the result in its data_gaps. ",
            4,
        )
        + ")\n",
    )
    rep(
        "src/cyclerfinder/data/validate.py",
        "# #1047 rows A, B, D (owner decision 2026-10-11",
        "# #1025 Earth-Venus k = 5 V2-ballistic promotions (lead order 2026-10-11, window 4; the final\n"
        "# spec §14 wording, L computed at 0.02 rev).\n"
        "_V2_1025_EV: str = (\n"
        + lit(
            "spec §14 V2-ballistic under the owner amendment of 2026-10-08 for exact-commensurate "
            "patched-conic multi-arc rows, with the same-phase drift wording of 2026-10-11 (#1007 "
            "secs. 8.4b, 8.17; L from cyclerfinder.verify.geometric_repeat at 0.02 rev: Earth-Venus "
            f"k = 5 -> L = 1), judged in {NOTE} secs. 11.1-11.6 "
            "(scripts/v2_ev10_1025.py, which imports final_record and final_controls from "
            "scripts/v2_samephase_1007.py unchanged; chains from scripts/ev_realeph_1025.py run "
            f"--n-cycles 10; {DATA}/): >= 3 of 5 epochs with a continuous ballistic Standish J2000 "
            "mean-element patched-conic chain of 10 cycles (80 yr, >= 2L), junction |V_inf| "
            "continuity < 1e-6 km/s, the demanded-turn gate passing at every flyby, DOP853 re-fly "
            "< 1 km, and same-phase date and |V_inf| drift (cycle i+1 against i, i = 1..7) inside the "
            "band (dates <= 5 % of the shortest Lambert leg, |V_inf| <= 10 % of the smallest "
            "departure |V_inf|). Slopes are diagnostics, not gates. Controls, run first with the same "
            "code: Hollister-Menning 1970 orbit 1 (k = 2, L = 5) passes; a broken chain fails; no "
            "published k = 5 member exists, so L = 1 is controlled through the code path only. The "
            "chains run to about 2110-2142, past the 1800-2050 fit span of the mean elements: a test "
            "of the multi-cycle structure. INDEPENDENCE gate: real-ephemeris body states, not the "
            "ideal model the row was found in. shared with primary path: the chain tool's Lambert "
            "solver, the full-rev seed rule and the cycle structure. Also carries the row's V1 "
            "evidence (_V1_1025_EV). ",
            4,
        )
        + ")\n\n"
        "# #1047 rows A, B, D (owner decision 2026-10-11",
    )
    for m in V2:
        s = text("src/cyclerfinder/data/validate.py")
        old_head = f'    ("{rid(m)}", "V1"): (\n        _V1_1025_EV + "Row: {m} ('
        if s.count(old_head) != 1:
            sys.exit(f"ABORT: validate.py entry for {m} found {s.count(old_head)} times")
        rep(
            "src/cyclerfinder/data/validate.py",
            old_head,
            f'    ("{rid(m)}", "V2"): (\n'
            f"        _V2_1025_EV\n"
            + lit(v2_lines[m], 8, first="+ ")
            + f'        + _V1_1025_EV\n        + "Row: {m} (',
        )

    # --- tests/data/test_schema_v45_fields.py ------------------------------------------------------
    for m in V2:
        rep(
            "tests/data/test_schema_v45_fields.py",
            f'        "{rid(m)}": "V1",\n',
            f'        # window 4 (#1025 sec. 11.6): V1 -> V2\n        "{rid(m)}": "V2",\n',
        )
    s = text("tests/data/test_schema_v45_fields.py")
    pat = re.compile(r'    assert sum\(1 for lvl in byid\.values\(\) if lvl == "V2"\) == (\d+)\n')
    found = pat.findall(s)
    if len(found) != 1:
        sys.exit(f"ABORT: V2-count assert found {len(found)} times")
    n_v2 = int(found[0])
    edits[ROOT / "tests/data/test_schema_v45_fields.py"] = pat.sub(
        f"    # {n_v2} -> {n_v2 + 6} (2026-10-12, window 4, #1025 secs. 11.6, 11.9): + ev5-4, ev5-5,\n"
        "    # ev5-6, ev5-7, ev5-8, ev5-9 V2-ballistic.\n"
        f'    assert sum(1 for lvl in byid.values() if lvl == "V2") == {n_v2 + 6}\n',
        s,
    )

    for p, s in edits.items():
        p.write_text(s)
        print("edited", p.relative_to(ROOT))
    # ruff format the two Python files (the commit hook would refuse them otherwise)
    subprocess.run(
        [
            "uv",
            "run",
            "ruff",
            "format",
            "src/cyclerfinder/data/validate.py",
            "tests/data/test_schema_v45_fields.py",
        ],
        check=True,
        cwd=ROOT,
    )
    print(f"OK: all edits applied (V2 count {n_v2} -> {n_v2 + 6})")


if __name__ == "__main__":
    main()

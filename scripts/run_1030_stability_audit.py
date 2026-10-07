"""#1030: audit orbit_elements.cr3bp.stability_index on every catalogue row that stores one.

For each row with a stored stability_index, state_nd, period_nd and mass_ratio, integrate the
full-period 6 x 6 monodromy at the row's own mu (core.cr3bp.propagate, stm_mode="fixed_path")
and compute:
- planar rows (z = zdot = 0): k_par = tr(M4) - 2 (in-plane, x y xdot ydot block) and
  k_perp = tr(Mz) (vertical z zdot block). Each is stable iff it lies in [-2, 2]. The Barden half
  index is nu = k_par / 2 (stable iff |nu| <= 1);
- spatial rows: the two nontrivial reciprocal pairs of the 6 x 6 monodromy, b_i = lambda + 1/lambda;
- the closure |x(T) - x(0)| of the stored state. A row that closes worse than 1e-8 (a rounded
  published state) is re-closed at its stored C with the symmetric corrector when it starts at a
  perpendicular crossing, and judged on the re-closed orbit; otherwise it is not judged;
- the spectral radius max |lambda| of the 6 x 6 monodromy (another convention in use).

The stored value is then identified as k_par, k_perp, nu = k_par/2, k_perp/2, or one of the
spatial b_i, by the smallest relative difference. The row's notes, name and field comments are
searched for "stable" wording (a word "stable" or "STABLE" not preceded by "un"/"UN").
Read only: data/catalogue.yaml is not modified. Output: data/1030_stability_audit/audit.json.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from cyclerfinder.core.cr3bp import CR3BPSystem, jacobi_constant, propagate
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "1030_stability_audit" / "audit.json"
CLOSURE_OK = 1e-8
STABLE_RE = re.compile(r"(?<![a-zA-Z])(?<!un)(?<!UN)(?<!Un)(stable|STABLE|Stable)(?![a-zA-Z])")


def stable_phrases(text: str) -> list[str]:
    out = []
    for m in STABLE_RE.finditer(text):
        a, b = max(0, m.start() - 60), min(len(text), m.end() + 40)
        out.append(text[a:b].replace("\n", " "))
    return out


def raw_row_text(lines: list[str], start: int) -> str:
    """The row's raw YAML text (comments included), from its '- id:' line to the next row."""
    end = start + 1
    while end < len(lines) and not lines[end].startswith("- id: "):
        end += 1
    return "".join(lines[start:end])


def main() -> None:
    preflight_search(
        task_no=1030,
        region_id="stability-index-audit-catalogue-cr3bp-rows",
        method=MethodCapability(
            genome="catalogue CR3BP rows with a stored stability_index",
            corrector="none (monodromy of the stored state)",
            capability_tags=frozenset({"cr3bp"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=50,
    )
    text = (REPO / "data" / "catalogue.yaml").read_text()
    lines = text.splitlines(keepends=True)
    starts = {
        ln[len("- id: ") :].strip(): i for i, ln in enumerate(lines) if ln.startswith("- id: ")
    }
    rows = yaml.safe_load(text)
    out: list[dict[str, Any]] = []
    for r in rows:
        c = ((r.get("orbit_elements") or {}).get("cr3bp")) or {}
        k_st = c.get("stability_index")
        if k_st is None:
            continue
        st, per, mu = c.get("state_nd"), c.get("period_nd"), c.get("mass_ratio")
        rec: dict[str, Any] = {"id": r["id"], "stored": float(k_st)}
        raw = raw_row_text(lines, starts[r["id"]])
        rec["stable_wording"] = stable_phrases(raw)
        if not (st and per and mu):
            rec["note"] = "no state/period/mass_ratio; not computed"
            out.append(rec)
            continue
        s0 = np.asarray(st, dtype=float)
        sysm = CR3BPSystem(mu=float(mu), primary="", secondary="", l_km=1.0, t_s=1.0)
        try:
            arc = propagate(sysm, s0, float(per), with_stm=True, stm_mode="fixed_path")
        except RuntimeError as exc:
            rec["note"] = f"propagation failed: {exc}"
            out.append(rec)
            continue
        m = np.asarray(arc.stm)
        rec["closure"] = float(np.linalg.norm(arc.state_f - s0))
        planar = bool(abs(s0[2]) < 1e-12 and abs(s0[5]) < 1e-12)
        rec["planar"] = planar
        symmetric = bool(planar and abs(s0[1]) < 1e-12 and abs(s0[3]) < 1e-12)
        rec["evaluated_on"] = "stored state"
        if rec["closure"] > CLOSURE_OK:
            if not symmetric:
                rec["note"] = "stored state does not close to 1e-8; not at a symmetric crossing"
                out.append(rec)
                print(f"{r['id']:48s} NOT JUDGED closure={rec['closure']:.1e}")
                continue
            # A rounded published state: re-close it at the stored C with the project's
            # symmetric corrector, then judge the corrected orbit (also give Barden's nu).
            # the row's stored C, not the C of its rounded state: for ross-rt-em-cycler-21 the
            # 10-digit state is 2e-11 off in C, more than its stable window is wide (#1031)
            jc = float(c.get("jacobi_constant") or jacobi_constant(s0, float(mu)))
            orb = correct_symmetric_fixed_jacobi(
                sysm, float(s0[0]), jc, float(per), ydot0_sign=float(np.sign(s0[4])), tol=1e-12
            )
            if not orb.converged:
                rec["note"] = "stored state does not close and the symmetric corrector failed"
                out.append(rec)
                continue
            s0 = np.array([orb.x0, 0.0, 0.0, 0.0, orb.ydot0, 0.0])
            arc = propagate(sysm, s0, orb.period, with_stm=True, stm_mode="fixed_path")
            m = np.asarray(arc.stm)
            rec["evaluated_on"] = (
                f"re-closed at stored C (dx0 {orb.x0 - st[0]:.2e}, dT {orb.period - per:.2e})"
            )
            rec["closure_reclosed"] = float(np.linalg.norm(arc.state_f - s0))
        rec["spectral_radius"] = float(np.max(np.abs(np.linalg.eigvals(m))))
        cands: dict[str, float] = {}
        if planar:
            m4 = m[np.ix_([0, 1, 3, 4], [0, 1, 3, 4])]
            mz = m[np.ix_([2, 5], [2, 5])]
            k_par = float(np.trace(m4) - 2.0)
            k_perp = float(np.trace(mz))
            rec.update({"k_par": k_par, "k_perp": k_perp})
            cands = {
                "k_par": k_par,
                "k_perp": k_perp,
                "nu=k_par/2": k_par / 2,
                "k_perp/2": k_perp / 2,
            }
            rec["inplane_stable"] = bool(abs(k_par) <= 2.0)
            rec["vertical_stable"] = bool(abs(k_perp) <= 2.0)
        else:
            eig = np.linalg.eigvals(m)
            order = np.argsort(np.abs(eig - 1.0))
            nontriv = eig[order[2:]]
            b = sorted({round(float((e + 1 / e).real), 9) for e in nontriv}, key=abs, reverse=True)
            rec["b_spatial"] = b
            for i, v in enumerate(b[:2]):
                cands[f"b{i + 1}"] = v
                cands[f"b{i + 1}/2"] = v / 2
            rec["spatial_stable"] = bool(all(abs(v) <= 2.0 for v in b[:2]))
        cands["spectral_radius"] = rec["spectral_radius"]
        rel = {k: abs(v - k_st) / max(abs(k_st), 1e-12) for k, v in cands.items()}
        best = min(rel, key=lambda k: rel[k])
        rec["stored_is"] = (
            best if rel[best] < 1e-3 else f"none (closest {best}, rel {rel[best]:.2e})"
        )
        rec["stored_rel_diff"] = rel[best]
        out.append(rec)
        print(
            f"{r['id']:48s} stored={k_st:+.5g} is={rec['stored_is']:28s} "
            + (
                f"k_par={rec['k_par']:+.4g} k_perp={rec['k_perp']:+.4g}"
                if planar
                else f"b={rec['b_spatial'][:2]}"
            )
            + f" closure={rec['closure']:.1e} stable_words={len(rec['stable_wording'])}"
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()

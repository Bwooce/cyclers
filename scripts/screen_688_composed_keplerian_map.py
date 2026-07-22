"""#688 -- composed two-map Keplerian screen (Stage A of #686's CCR4BP plan).

Composes the positive-controlled RS07 exterior Keplerian map (`#500`) into an
alternating Jupiter-Europa / Jupiter-Ganymede system and searches for periodic
alternating itineraries phase-locked to the moons' ~2:1 Laplace commensurability.

SCREEN-GRADE ONLY: each map ignores the other moon while in control.  A negative
here is NOT registry-grade (no `empty_regions.jsonl` stamp); a positive is a SEED
for Stage B, not a result.  No catalogue writeback.

Run:  uv run python scripts/screen_688_composed_keplerian_map.py
Outputs -> data/found/688_composed_keplerian_map/result.json
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.genome.composed_moon_map import (  # noqa: E402
    _GM_JUPITER,
    ComposedMoonMap,
    ComposedState,
    apoapsis_norm,
    closure_residual,
    eccentricity_from_tisserand,
    moon_config,
    periapsis_norm,
    resonance_semimajor,
)

OUT_DIR = ROOT / "data" / "found" / "688_composed_keplerian_map"
C_J_REF = 3.0  # RS07 running-example Tisserand constant


def _f(x: float) -> float:
    return float(x)


# ---------------------------------------------------------------------------
# 1. Single-moon sanity (PC0 applied to Europa & Ganymede individually)
# ---------------------------------------------------------------------------


def single_moon_sanity(cm: ComposedMoonMap) -> dict[str, Any]:
    """Confirm each moon's map reproduces RS07 structure: f(0)~0, odd, 1:2 fixed point."""
    out: dict[str, Any] = {}
    a_res_12 = resonance_semimajor(1, 2)  # 2^(2/3) = 1.5874, the 1:2 resonance
    k_res = -0.5 / a_res_12
    for moon, kmap in cm.maps.items():
        f0 = kmap.kick(0.0)
        f_pos = kmap.kick(0.05)
        f_neg = kmap.kick(-0.05)
        odd_resid = abs(f_pos + f_neg) / (abs(f_pos) + abs(f_neg) + 1e-15)
        # 1:2 fixed point: one step returns (omega->0 mod 2pi, K unchanged)
        w1, k1 = kmap.step(0.0, k_res, u=0.0)
        out[moon] = {
            "mu": _f(cm.cfg[moon].mu),
            "sma_km": _f(cm.cfg[moon].sma_km),
            "period_days": _f(cm.cfg[moon].period_s / 86400.0),
            "f_at_0": _f(f0),
            "odd_residual": _f(odd_resid),
            "fixed_point_dK": _f(abs(k1 - k_res)),
            "fixed_point_omega": _f(w1),
            "sign_ok": bool(f_pos < 0.0 < f_neg),  # RS07 sign convention
        }
    return out


# ---------------------------------------------------------------------------
# 2. Regime-overlap analysis (the structural obstruction, quantified)
# ---------------------------------------------------------------------------


def regime_analysis(cm: ComposedMoonMap) -> dict[str, Any]:
    """Quantify whether the two exterior maps' valid encounter regimes overlap.

    (a) Periapsis-shell disjointness: a Europa encounter pins r_peri ~ Europa's
        orbit; a Ganymede encounter pins r_peri ~ Ganymede's orbit; the exterior
        map conserves periapsis, so no single orbit alternates.
    (b) The CCR4BP paper's key Jupiter-Ganymede resonances (3:2, 4:3) are
        INTERIOR (a<1) -> outside the RS07 exterior map entirely.
    """
    eur = cm.cfg["Europa"]
    gan = cm.cfg["Ganymede"]
    # Europa exterior-encounter shell in PHYSICAL km, and its image in Ganymede units.
    eur_peri_lo_km = 1.0 * eur.sma_km
    eur_peri_hi_km = cm.shell_hi * eur.sma_km
    gan_peri_lo_km = 1.0 * gan.sma_km
    gan_peri_hi_km = cm.shell_hi * gan.sma_km
    # Does the Europa encounter shell (in km) ever reach the Ganymede shell?
    shells_overlap = eur_peri_hi_km >= gan_peri_lo_km

    # Paper resonances: is each interior or exterior?
    paper_res = {
        "Europa_3:4": resonance_semimajor(3, 4),
        "Europa_5:4": resonance_semimajor(5, 4),
        "Ganymede_3:2": resonance_semimajor(3, 2),
        "Ganymede_4:3": resonance_semimajor(4, 3),
    }
    res_class = {
        k: {"a_norm": _f(v), "regime": "exterior" if v > 1.0 else "interior"}
        for k, v in paper_res.items()
    }
    return {
        "europa_encounter_shell_km": [_f(eur_peri_lo_km), _f(eur_peri_hi_km)],
        "ganymede_encounter_shell_km": [_f(gan_peri_lo_km), _f(gan_peri_hi_km)],
        "shells_overlap": bool(shells_overlap),
        "europa_peri_in_ganymede_units": _f(eur.sma_km / gan.sma_km),
        "sma_ratio_gan_over_eur": _f(gan.sma_km / eur.sma_km),
        "period_ratio_gan_over_eur": _f(gan.period_s / eur.period_s),
        "paper_resonance_regime": res_class,
    }


# ---------------------------------------------------------------------------
# 3. Composed itinerary search (physical-preserving patch)
# ---------------------------------------------------------------------------


def _state_at_europa_resonance(
    cm: ComposedMoonMap, p: int, q: int, varpi: float, c_j: float
) -> ComposedState:
    a_norm = resonance_semimajor(p, q)
    k = -0.5 / a_norm
    e = eccentricity_from_tisserand(k, c_j)
    a_phys = a_norm * cm.cfg["Europa"].sma_km
    return ComposedState(a_phys_km=a_phys, e=e, varpi_rad=varpi, t_s=0.0)


def itinerary_search(cm: ComposedMoonMap) -> dict[str, Any]:
    """Scan alternating [(Europa,N_E),(Ganymede,N_G)] itineraries for closure.

    Reports the best physical-element closure residual AND -- decisively -- the
    fraction of Ganymede-segment passages that were self-consistent encounters
    (r_peri in the exterior shell).  If that fraction is ~0, the Ganymede map is
    dynamically inert after a Europa patch and no genuine alternation occurs.
    """
    best: dict[str, Any] | None = None
    gan_valid_total = 0
    gan_steps_total = 0
    eur_valid_total = 0
    eur_steps_total = 0
    n_itins = 0

    # Seed from exterior Europa resonances the CCR4BP literature names.
    europa_res = [(3, 4), (2, 3), (3, 5), (1, 2)]
    for p, q in europa_res:
        for varpi in [i * math.pi / 6 for i in range(12)]:
            start = _state_at_europa_resonance(cm, p, q, varpi, C_J_REF)
            if math.isnan(start.e):
                continue
            for n_e in (2, 3, 4, 6):
                for n_g in (1, 2, 3):
                    itin = [("Europa", n_e), ("Ganymede", n_g), ("Europa", n_e), ("Ganymede", n_g)]
                    end, recs = cm.run_itinerary(start, itin)
                    n_itins += 1
                    if math.isnan(end.e):
                        continue
                    for r in recs:
                        if r.moon == "Ganymede":
                            gan_steps_total += 1
                            gan_valid_total += int(r.encounter_valid)
                        else:
                            eur_steps_total += 1
                            eur_valid_total += int(r.encounter_valid)
                    resid = closure_residual(start, end)
                    score = resid["da_rel"] + resid["de"] + 0.1 * resid["dvarpi_rad"]
                    if best is None or score < best["score"]:
                        best = {
                            "score": _f(score),
                            "europa_resonance": f"{p}:{q}",
                            "varpi_rad": _f(varpi),
                            "n_e": n_e,
                            "n_g": n_g,
                            "residual": {k: _f(v) for k, v in resid.items()},
                            "gan_encounters_valid_in_itin": sum(
                                1 for r in recs if r.moon == "Ganymede" and r.encounter_valid
                            ),
                            "gan_steps_in_itin": sum(1 for r in recs if r.moon == "Ganymede"),
                        }
    return {
        "n_itineraries_tried": n_itins,
        "ganymede_selfconsistent_fraction": _f(
            gan_valid_total / gan_steps_total if gan_steps_total else 0.0
        ),
        "europa_selfconsistent_fraction": _f(
            eur_valid_total / eur_steps_total if eur_steps_total else 0.0
        ),
        "ganymede_steps_total": gan_steps_total,
        "europa_steps_total": eur_steps_total,
        "best_closure": best,
    }


# ---------------------------------------------------------------------------
# 4. Seed geometry: Laplace-compatible resonance pairs + transfer geometry
# ---------------------------------------------------------------------------


def seed_geometry(cm: ComposedMoonMap) -> dict[str, Any]:
    """Enumerate Europa/Ganymede resonance pairs and their transfer geometry.

    For each (Europa exterior resonance, Ganymede resonance) pair, report the
    physical apsides and the minimum impulsive dv to connect the two orbits at
    a shared radius (if their [peri,apo] ranges overlap), plus the moon-phase
    (Laplace) super-period commensurability.  This is the SEED geometry Stage B
    would target -- explicitly not a result.
    """
    eur = cm.cfg["Europa"]
    gan = cm.cfg["Ganymede"]
    europa_res = [(3, 4), (2, 3), (3, 5), (5, 7)]
    ganymede_res = [(3, 2), (4, 3), (5, 4)]  # interior side (a<1)
    pairs: list[dict[str, Any]] = []
    for pe, qe in europa_res:
        ae = resonance_semimajor(pe, qe)
        ke = -0.5 / ae
        ee = eccentricity_from_tisserand(ke, C_J_REF)
        if math.isnan(ee):
            continue
        eur_peri_km = periapsis_norm(ae, ee) * eur.sma_km
        eur_apo_km = apoapsis_norm(ae, ee) * eur.sma_km
        for pg, qg in ganymede_res:
            ag = resonance_semimajor(pg, qg)
            kg = -0.5 / ag
            eg = eccentricity_from_tisserand(kg, C_J_REF)
            if math.isnan(eg):
                continue
            gan_peri_km = periapsis_norm(ag, eg) * gan.sma_km
            gan_apo_km = apoapsis_norm(ag, eg) * gan.sma_km
            # Overlap of the two orbits' radial ranges (necessary for a ballistic patch).
            lo = max(eur_peri_km, gan_peri_km)
            hi = min(eur_apo_km, gan_apo_km)
            overlap = hi >= lo
            dv_km_s = float("nan")
            if overlap:
                r_patch = 0.5 * (lo + hi)  # a shared radius both orbits reach
                # circular-restricted vis-viva speeds at r_patch on each orbit
                a_e_km = ae * eur.sma_km
                a_g_km = ag * gan.sma_km
                v_e = math.sqrt(_GM_JUPITER * (2.0 / r_patch - 1.0 / a_e_km))
                v_g = math.sqrt(_GM_JUPITER * (2.0 / r_patch - 1.0 / a_g_km))
                dv_km_s = abs(v_e - v_g)  # coarse coplanar speed-magnitude bound
            # Laplace phasing: spacecraft:Europa = pe:qe, spacecraft:Ganymede = pg:qg.
            # A super-period closing both requires a common multiple of the segment times.
            pairs.append(
                {
                    "europa_res": f"{pe}:{qe}",
                    "ganymede_res": f"{pg}:{qg}",
                    "europa_a_km": _f(ae * eur.sma_km),
                    "ganymede_a_km": _f(ag * gan.sma_km),
                    "europa_apsides_km": [_f(eur_peri_km), _f(eur_apo_km)],
                    "ganymede_apsides_km": [_f(gan_peri_km), _f(gan_apo_km)],
                    "radial_overlap": bool(overlap),
                    "patch_dv_km_s": _f(dv_km_s),
                    "ganymede_regime": "exterior" if ag > 1.0 else "interior",
                }
            )
    # rank by dv among overlapping pairs
    overlapping = [p for p in pairs if p["radial_overlap"] and not math.isnan(p["patch_dv_km_s"])]
    overlapping.sort(key=lambda d: d["patch_dv_km_s"])
    return {
        "n_pairs": len(pairs),
        "n_overlapping": len(overlapping),
        "pairs": pairs,
        "best_overlapping_by_dv": overlapping[:5],
    }


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main() -> None:
    t0 = time.time()
    eur = moon_config("Europa")
    gan = moon_config("Ganymede")
    cm = ComposedMoonMap(eur, gan, c_j_ref=C_J_REF, a_ref_norm=1.35, shell_hi=1.30)

    result: dict[str, Any] = {
        "task": "#688 composed two-map Keplerian screen (Stage A of #686 CCR4BP)",
        "screen_grade_only": True,
        "c_j_ref": C_J_REF,
        "single_moon_sanity": single_moon_sanity(cm),
        "regime_analysis": regime_analysis(cm),
        "itinerary_search": itinerary_search(cm),
        "seed_geometry": seed_geometry(cm),
    }
    result["elapsed_s"] = _f(time.time() - t0)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "result.json"
    out_path.write_text(json.dumps(result, indent=2))

    # -- concise console summary --
    ra = result["regime_analysis"]
    it = result["itinerary_search"]
    sg = result["seed_geometry"]
    print(f"[#688] composed Keplerian-map screen  ({result['elapsed_s']:.1f}s)")
    print(f"  shells_overlap (Europa vs Ganymede encounter shells): {ra['shells_overlap']}")
    print(
        f"  paper Ganymede resonances regime: "
        f"3:2={ra['paper_resonance_regime']['Ganymede_3:2']['regime']}, "
        f"4:3={ra['paper_resonance_regime']['Ganymede_4:3']['regime']}"
    )
    print(f"  itineraries tried: {it['n_itineraries_tried']}")
    print(
        f"  Ganymede self-consistent-encounter fraction: "
        f"{it['ganymede_selfconsistent_fraction']:.4f} "
        f"(Europa: {it['europa_selfconsistent_fraction']:.4f})"
    )
    if it["best_closure"]:
        bc = it["best_closure"]
        print(
            f"  best closure: Eur {bc['europa_resonance']} N_E={bc['n_e']} N_G={bc['n_g']} "
            f"da_rel={bc['residual']['da_rel']:.3e} de={bc['residual']['de']:.3e} "
            f"gan_valid_encounters={bc['gan_encounters_valid_in_itin']}"
        )
    print(f"  seed pairs: {sg['n_pairs']} ({sg['n_overlapping']} radially overlapping)")
    if sg["best_overlapping_by_dv"]:
        b = sg["best_overlapping_by_dv"][0]
        print(
            f"  cheapest ballistic-ish patch: Eur {b['europa_res']} <-> Gan {b['ganymede_res']} "
            f"dv~{b['patch_dv_km_s']:.3f} km/s ({b['ganymede_regime']} Ganymede)"
        )
    print(f"  wrote {out_path}")


if __name__ == "__main__":
    main()

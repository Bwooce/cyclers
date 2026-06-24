"""Discovery probe (registry-clean): leveraging-VILM endgame at Uranus/Neptune moon
systems. Only the NON-leveraging repeated-moon lane swept these; the phase-full VILM
endgame (solve_endgame) did not -> distinct method capability, not a re-sweep.
REPORT-ONLY: prints findings/survivors per topology; no data-file writes (the
empty-region registry emit is done under review after this completes)."""

from datetime import UTC, datetime

from cyclerfinder.data.discover_novel import _tour_topologies, discover_endgame_moon

URANUS_TOURS = [
    ("Ariel", "Umbriel", "Titania", "Ariel"),
    ("Miranda", "Ariel", "Umbriel", "Miranda"),
    ("Umbriel", "Titania", "Oberon", "Umbriel"),
    ("Ariel", "Umbriel", "Oberon", "Ariel"),
]
NEPTUNE_TOURS = [("Triton", "Proteus", "Triton")]


def run(center, tours):
    specs = []
    for seq in tours:
        specs.extend(_tour_topologies(seq, period_ks=(1, 2)))
    print(f"\n=== {center}: {len(tours)} tours -> {len(specs)} topology specs ===", flush=True)
    survivors = []
    for f in discover_endgame_moon(
        topologies=specs,
        center=center,
        target_vinf_floor_kms=6.0,
        dv_budget_kms=4.0,
        n_epochs=16,
        span_days=8.0,
        vinf_cap=14.0,
    ):
        survivors.append(f)
        seq = getattr(f, "sequence", None) or getattr(
            getattr(f, "candidate", None), "sequence", None
        )
        print(
            f"  SURVIVOR: powered={getattr(f, 'powered', None)} seq={seq} "
            f"vinf={getattr(f, 'vinf_kms', getattr(f, 'best_vinf_kms', '?'))} "
            f"dv={getattr(f, 'total_dv_kms', getattr(f, 'dv_kms', '?'))}",
            flush=True,
        )
    print(f"  {center}: {len(survivors)} survivor(s)", flush=True)
    return survivors


print(f"START {datetime.now(UTC).isoformat()}", flush=True)
allsurv = run("Uranus", URANUS_TOURS) + run("Neptune", NEPTUNE_TOURS)
print(f"\nTOTAL survivors across Uranus+Neptune endgame: {len(allsurv)}", flush=True)
print(f"END {datetime.now(UTC).isoformat()}", flush=True)

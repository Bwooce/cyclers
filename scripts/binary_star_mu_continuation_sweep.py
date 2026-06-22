import json
import logging

from cyclerfinder.core.cr3bp import CR3BPSystem
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi
from cyclerfinder.search.mu_continuation import continue_in_mu, scan_c_family_at_mu

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

seeds = [
    {
        "id": "ross-rt-em-cycler-31-2025",
        "k1": 3,
        "k2": 1,
        "mu_start": 1.2150584270572e-2,
        "c_start": 3.161784147013429,
        "period_start": 14.78849241668140,
        "x0_guess": -0.3209891696,
        "ydot0_sign": -1.0,
        "mu_target": 0.3,
    },
    {
        "id": "ross-rt-em-cycler-11-2025",
        "k1": 1,
        "k2": 1,
        "mu_start": 1.2150584270572e-2,
        "c_start": 3.151175879508174,
        "period_start": 10.29206921007976,
        "x0_guess": -0.7682140805,
        "ydot0_sign": -1.0,
        "mu_target": 0.5,
    },
]


def main():
    negatives = []

    for seed_info in seeds:
        logging.info(f"=== Sweeping {seed_info['id']} to mu = {seed_info['mu_target']} ===")
        system_start = CR3BPSystem(
            mu=seed_info["mu_start"], primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8
        )

        logging.info("Recovering the seed...")
        seed_orbit = correct_symmetric_fixed_jacobi(
            system=system_start,
            x0_guess=seed_info["x0_guess"],
            jacobi=seed_info["c_start"],
            period_guess=seed_info["period_start"],
            ydot0_sign=seed_info["ydot0_sign"],
            half_crossings=None,
        )

        if not seed_orbit.converged:
            logging.error("Failed to recover seed.")
            continue

        logging.info(
            f"Seed recovered: x0={seed_orbit.x0}, period={seed_orbit.period}, half_crossings guess established."
        )

        # We need to determine the half_crossings index that `correct_symmetric_fixed_jacobi` auto-determined
        # We can just look at `seed_orbit.t_half` and how it corresponds to the crossings? Wait, `continue_in_mu` requires `half_crossings` as an int.
        # How do we get the half_crossings index?
        # Let's count the crossings up to t_half + eps
        import numpy as np

        from cyclerfinder.search.cr3bp_periodic import _xaxis_crossings

        state0 = np.array([seed_orbit.x0, 0.0, 0.0, 0.0, seed_orbit.ydot0, 0.0])
        times, _ = _xaxis_crossings(
            system_start, state0, seed_orbit.t_half * 1.1, with_stm=False, rtol=1e-12, atol=1e-12
        )
        idx = int(np.argmin(np.abs(times - seed_orbit.t_half))) + 1
        logging.info(f"Seed half_crossings index determined as: {idx}")

        logging.info(f"Continuing {seed_info['id']} to mu = {seed_info['mu_target']} ...")
        branch = continue_in_mu(
            seed=seed_orbit,
            mu_start=seed_info["mu_start"],
            half_crossings=idx,
            ydot0_sign=seed_info["ydot0_sign"],
            mu_target=seed_info["mu_target"],
            ds_max=1e-2,
            period_jump_frac=0.3,  # safeguard against branch switching
        )

        logging.info(f"Branch stopped with reason: {branch.stop_reason}")
        if not branch.members:
            logging.error("No members on branch.")
            continue

        final_member = branch.members[-1]
        logging.info(f"Final reached mu: {final_member.mu}")

        target_reached = abs(final_member.mu - seed_info["mu_target"]) < 1e-6
        if target_reached:
            logging.info("Target mu reached. Checking topology...")
            topo = winding_topology(final_member.mu, final_member.state0, final_member.period)
            logging.info(
                f"Topology at target: k1={topo.k1}, k2={topo.k2}, prograde={topo.prograde}"
            )

            if topo.k1 == seed_info["k1"] and topo.k2 == seed_info["k2"] and topo.prograde:
                logging.info("Topology matches target cycler family! Checking stability...")
                if final_member.stable:
                    logging.info(f"STABLE WINDOW FOUND at target mu! nu = {final_member.nu}")
                    # Novelty check required!
                else:
                    logging.info(
                        f"Member at target mu is unstable (nu = {final_member.nu}). Scanning C-family..."
                    )
                    members = scan_c_family_at_mu(
                        mu=seed_info["mu_target"],
                        x0_guess=final_member.x0,
                        c_center=final_member.jacobi,
                        period_guess=final_member.period,
                        half_crossings=idx,
                        ydot0_sign=seed_info["ydot0_sign"],
                        dc=5e-4,
                        n_each=20,
                    )
                    stable_members = [m for m in members if m.stable]
                    if stable_members:
                        logging.info(
                            f"Found {len(stable_members)} stable members in C-scan! Novelty check required."
                        )
                    else:
                        logging.info("No stable members found in C-scan.")
                        negatives.append(
                            {
                                "task": 315,
                                "method": "mu-continuation",
                                "seed_id": seed_info["id"],
                                "target_mu": seed_info["mu_target"],
                                "result": "Topology reached, but C-scan found no stable window.",
                            }
                        )
            else:
                logging.info(
                    f"Topology mismatch. Expected ({seed_info['k1']},{seed_info['k2']}), got ({topo.k1},{topo.k2})."
                )
                negatives.append(
                    {
                        "task": 315,
                        "method": "mu-continuation",
                        "seed_id": seed_info["id"],
                        "target_mu": seed_info["mu_target"],
                        "result": f"Target mu reached but topology changed to ({topo.k1},{topo.k2}) prograde={topo.prograde}.",
                    }
                )
        else:
            logging.info(
                f"Target mu NOT reached. Branch stopped at mu = {final_member.mu} due to {branch.stop_reason}"
            )
            negatives.append(
                {
                    "task": 315,
                    "method": "mu-continuation",
                    "seed_id": seed_info["id"],
                    "target_mu": seed_info["mu_target"],
                    "result": f"Branch terminated early at mu={final_member.mu:.6f} due to {branch.stop_reason}.",
                }
            )

    if negatives:
        with open("data/empty_regions.jsonl", "a") as f:
            for neg in negatives:
                f.write(json.dumps(neg) + "\n")
        logging.info(f"Appended {len(negatives)} negative results to data/empty_regions.jsonl")


if __name__ == "__main__":
    main()

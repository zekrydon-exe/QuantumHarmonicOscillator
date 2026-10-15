import numpy as np
from scipy.linalg import eigh

from src.hamiltonian import create_grid, build_hamiltonian


def main():
    domain_limits = [4.0, 6.0, 8.0, 10.0]

    target_dq = 0.04
    n_states = 5

    results = []

    for q_max in domain_limits:

        q_min = -q_max

        total_intervals = round(
            (q_max - q_min) / target_dq
        )

        n_points = total_intervals - 1

        q, dq = create_grid(
            q_min,
            q_max,
            n_points
        )

        H = build_hamiltonian(q, dq)

        energies, wavefunctions = eigh(
            H,
            subset_by_index=(0, n_states - 1)
        )

        exact_energies = np.arange(n_states) + 0.5

        absolute_errors = np.abs(
            energies - exact_energies
        )

        results.append({
            "q_max": q_max,
            "n_points": n_points,
            "dq": dq,
            "energies": energies,
            "absolute_errors": absolute_errors,
        })

    print("DOMAIN-CONVERGENCE STUDY")
    print("=" * 80)

    print(f"Target grid spacing: {target_dq}")
    print(f"Number of states: {n_states}")

    for result in results:

        print("\n" + "-" * 80)

        print(
            f"Domain: [{-result['q_max']}, "
            f"{result['q_max']}]"
        )

        print(
            f"Interior points: {result['n_points']} | "
            f"Grid spacing: {result['dq']:.6f}"
        )

        print("-" * 80)

        print(
            f"{'n':>3} "
            f"{'Numerical Energy':>18} "
            f"{'Exact Energy':>15} "
            f"{'Absolute Error':>18}"
        )

        print("-" * 80)

        exact_energies = np.arange(n_states) + 0.5

        for n in range(n_states):

            numerical = result["energies"][n]
            exact = exact_energies[n]

            error = result["absolute_errors"][n]

            print(
                f"{n:3d} "
                f"{numerical:18.10f} "
                f"{exact:15.8f} "
                f"{error:18.10e}"
            )

    print("\n\nCHANGE IN ENERGY BETWEEN DOMAINS")
    print("=" * 80)

    print(
        f"{'n':>3} "
        f"{'4 -> 6':>18} "
        f"{'6 -> 8':>18} "
        f"{'8 -> 10':>18}"
    )

    print("-" * 65)

    for n in range(n_states):

        changes = []

        for i in range(len(results) - 1):

            energy_1 = results[i]["energies"][n]
            energy_2 = results[i + 1]["energies"][n]

            changes.append(abs(energy_2 - energy_1))

        print(
            f"{n:3d} "
            f"{changes[0]:18.10e} "
            f"{changes[1]:18.10e} "
            f"{changes[2]:18.10e}"
        )


if __name__ == "__main__":
    main()
import numpy as np
from scipy.linalg import eigh

from src.hamiltonian import create_grid, build_hamiltonian


def main():
    q_min = -4.0
    q_max = 4.0

    grid_sizes = [49, 99, 199, 399]
    n_states = 5

    results = []

    for n_points in grid_sizes:

        q, dq = create_grid(q_min, q_max, n_points)

        H = build_hamiltonian(q, dq)

        # Calculate only the five lowest eigenvalues.
        energies, wavefunctions = eigh(
            H,
            subset_by_index=(0, n_states - 1)
        )

        # Exact dimensionless harmonic-oscillator energies.
        quantum_numbers = np.arange(n_states)
        exact_energies = quantum_numbers + 0.5

        # Calculate absolute and relative errors.
        absolute_errors = np.abs(energies - exact_energies)

        relative_errors = (
            absolute_errors / np.abs(exact_energies)
        )

        results.append({
            "n_points": n_points,
            "dq": dq,
            "energies": energies,
            "absolute_errors": absolute_errors,
            "relative_errors": relative_errors,
        })

    # --------------------------------------------------
    # 3. Print the energy convergence table
    # --------------------------------------------------

    print("GRID CONVERGENCE STUDY")
    print("=" * 85)
    print(f"Domain: [{q_min}, {q_max}]")
    print(f"Number of states analysed: {n_states}")

    for result in results:

        print("\n" + "-" * 85)
        print(
            f"Grid points: {result['n_points']:4d} | "
            f"Grid spacing: {result['dq']:.6f}"
        )

        print("-" * 85)

        print(
            f"{'n':>3} "
            f"{'Numerical Energy':>18} "
            f"{'Exact Energy':>15} "
            f"{'Absolute Error':>18} "
            f"{'Relative Error (%)':>20}"
        )

        print("-" * 85)

        exact_energies = np.arange(n_states) + 0.5

        for n in range(n_states):

            numerical = result["energies"][n]
            exact = exact_energies[n]

            absolute_error = result["absolute_errors"][n]
            relative_error = result["relative_errors"][n] * 100

            print(
                f"{n:3d} "
                f"{numerical:18.10f} "
                f"{exact:15.8f} "
                f"{absolute_error:18.10e} "
                f"{relative_error:20.8f}"
            )

    # --------------------------------------------------
    # 4. Estimate the convergence order
    # --------------------------------------------------

    print("\n\nOBSERVED CONVERGENCE ORDER")
    print("=" * 85)

    print(
        f"{'n':>3} "
        f"{'49 -> 99':>15} "
        f"{'99 -> 199':>15} "
        f"{'199 -> 399':>15}"
    )

    print("-" * 55)

    for n in range(n_states):

        observed_orders = []

        for i in range(len(results) - 1):

            h1 = results[i]["dq"]
            h2 = results[i + 1]["dq"]

            e1 = results[i]["absolute_errors"][n]
            e2 = results[i + 1]["absolute_errors"][n]

            if e1 > 0 and e2 > 0:

                order = np.log(e1 / e2) / np.log(h1 / h2)

            else:
                order = np.nan

            observed_orders.append(order)

        print(
            f"{n:3d} "
            f"{observed_orders[0]:15.6f} "
            f"{observed_orders[1]:15.6f} "
            f"{observed_orders[2]:15.6f}"
        )


if __name__ == "__main__":
    main()
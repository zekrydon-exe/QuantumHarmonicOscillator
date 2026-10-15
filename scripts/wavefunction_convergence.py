import numpy as np

from src.harmonic_oscillator import solve_harmonic_oscillator
from src.analytical_solution import analytical_wavefunction


def main():
    q_min = -8.0
    q_max = 8.0

    grid_sizes = [99, 199, 399, 799]
    n_states = 5

    results = []

    for n_points in grid_sizes:

        q, dq, energies, wavefunctions = solve_harmonic_oscillator(
            q_min=q_min,
            q_max=q_max,
            n_points=n_points,
            n_states=n_states,
        )

        exact_wavefunctions = np.column_stack([
            analytical_wavefunction(q, n)
            for n in range(n_states)
        ])

        relative_errors = []

        for n in range(n_states):

            numerical = wavefunctions[:, n]
            analytical = exact_wavefunctions[:, n]

            # Align the arbitrary eigenvector sign.
            overlap = np.trapezoid(
                numerical * analytical,
                q,
            )

            if overlap < 0:
                numerical = -numerical

            # Calculate the relative L2 error.
            difference = numerical - analytical

            error_norm = np.sqrt(
                np.trapezoid(np.abs(difference)**2, q)
            )

            analytical_norm = np.sqrt(
                np.trapezoid(np.abs(analytical)**2, q)
            )

            relative_error = error_norm / analytical_norm

            relative_errors.append(relative_error)

        results.append({
            "n_points": n_points,
            "dq": dq,
            "energies": energies,
            "relative_errors": np.array(relative_errors),
        })

    print("WAVEFUNCTION CONVERGENCE STUDY")
    print("=" * 85)

    print(f"Domain: [{q_min}, {q_max}]")
    print(f"Number of states: {n_states}")

    print("\nRELATIVE L2 ERRORS")
    print("-" * 85)

    print(
        f"{'Grid points':>12} "
        f"{'dq':>10} "
        + " ".join(
            f"{'n=' + str(n):>15}"
            for n in range(n_states)
        )
    )

    print("-" * 85)

    for result in results:

        print(
            f"{result['n_points']:12d} "
            f"{result['dq']:10.5f} "
            + " ".join(
                f"{error:15.8e}"
                for error in result["relative_errors"]
            )
        )

    print("\nOBSERVED CONVERGENCE ORDERS")
    print("-" * 85)

    print(
        f"{'State':>8} "
        f"{'99 -> 199':>15} "
        f"{'199 -> 399':>15} "
        f"{'399 -> 799':>15}"
    )

    print("-" * 60)

    for n in range(n_states):

        orders = []

        for i in range(len(results) - 1):

            e1 = results[i]["relative_errors"][n]
            e2 = results[i + 1]["relative_errors"][n]

            h1 = results[i]["dq"]
            h2 = results[i + 1]["dq"]

            order = np.log(e1 / e2) / np.log(h1 / h2)

            orders.append(order)

        print(
            f"{n:8d} "
            f"{orders[0]:15.6f} "
            f"{orders[1]:15.6f} "
            f"{orders[2]:15.6f}"
        )


if __name__ == "__main__":
    main()
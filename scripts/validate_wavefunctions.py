import numpy as np

from src.harmonic_oscillator import solve_harmonic_oscillator
from src.analytical_solution import (
    analytical_energy,
    analytical_wavefunction,
)


def main():
    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=-8.0,
        q_max=8.0,
        n_points=399,
        n_states=5,
    )

    exact_wavefunctions = np.column_stack([
        analytical_wavefunction(q, n)
        for n in range(5)
    ])

    exact_energies = np.array([
        analytical_energy(n)
        for n in range(5)
    ])

    print("WAVEFUNCTION NORMALISATION")
    print("=" * 70)

    print(
        f"{'n':>3} "
        f"{'Numerical norm':>20} "
        f"{'Analytical norm':>20}"
    )

    print("-" * 50)

    for n in range(5):

        numerical_norm = np.trapezoid(
            np.abs(wavefunctions[:, n])**2,
            q,
        )

        analytical_norm = np.trapezoid(
            np.abs(exact_wavefunctions[:, n])**2,
            q,
        )

        print(
            f"{n:3d} "
            f"{numerical_norm:20.12f} "
            f"{analytical_norm:20.12f}"
        )

    for n in range(5):

        overlap = np.trapezoid(
            wavefunctions[:, n] * exact_wavefunctions[:, n],
            q,
        )

        if overlap < 0:
            wavefunctions[:, n] *= -1

    print("\nWAVEFUNCTION ERROR ANALYSIS")
    print("=" * 70)

    print(
        f"{'n':>3} "
        f"{'Overlap':>18} "
        f"{'Relative L2 error':>22}"
    )

    print("-" * 50)

    for n in range(5):

        difference = (
            wavefunctions[:, n] - exact_wavefunctions[:, n]
        )

        error_norm = np.sqrt(
            np.trapezoid(np.abs(difference)**2, q)
        )

        exact_norm = np.sqrt(
            np.trapezoid(
                np.abs(exact_wavefunctions[:, n])**2,
                q,
            )
        )

        relative_error = error_norm / exact_norm

        overlap = np.trapezoid(
            wavefunctions[:, n] * exact_wavefunctions[:, n],
            q,
        )

        print(
            f"{n:3d} "
            f"{overlap:18.10f} "
            f"{relative_error:22.10e}"
        )

    print("\nENERGY VALIDATION")
    print("=" * 70)

    for n in range(5):

        print(
            f"n = {n}: "
            f"numerical = {energies[n]:.10f}, "
            f"analytical = {exact_energies[n]:.10f}"
        )


if __name__ == "__main__":
    main()
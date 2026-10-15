import numpy as np
from scipy.linalg import eigh

from src.hamiltonian import create_grid, build_hamiltonian


def main():
    q_min = -4.0
    q_max = 4.0
    n_points = 99

    q, dq = create_grid(q_min, q_max, n_points)

    H = build_hamiltonian(q, dq)

    energies, wavefunctions = eigh(H)

    n_states_to_show = 5

    energies = energies[:n_states_to_show]
    wavefunctions = wavefunctions[:, :n_states_to_show]

    residuals = H @ wavefunctions - wavefunctions * energies

    residual_norms = np.linalg.norm(residuals, axis=0)

    print("\nEIGENPAIR VALIDATION")
    print("-" * 40)
    print("Maximum residual norm:", np.max(residual_norms))

    quantum_numbers = np.arange(len(energies))
    analytical_energies = quantum_numbers + 0.5

    print("QUANTUM HARMONIC OSCILLATOR")
    print("=" * 65)

    print(f"Domain: [{q_min}, {q_max}]")
    print(f"Number of interior grid points: {n_points}")
    print(f"Grid spacing: {dq:.6f}")

    print("\nENERGY COMPARISON")
    print("-" * 65)

    print(
        f"{'n':>3} "
        f"{'Numerical':>15} "
        f"{'Analytical':>15} "
        f"{'Absolute Error':>15} "
        f"{'Relative Error (%)':>18}"
    )

    print("-" * 65)

    for n, (numerical, analytical) in enumerate(
        zip(energies, analytical_energies)
    ):
        absolute_error = abs(numerical - analytical)

        relative_error = (
            absolute_error / abs(analytical)
        ) * 100

        print(
            f"{n:3d} "
            f"{numerical:15.8f} "
            f"{analytical:15.8f} "
            f"{absolute_error:15.8f} "
            f"{relative_error:18.6f}"
        )

    print("\nEIGENVECTOR INFORMATION")
    print("-" * 40)

    print("Hamiltonian shape:", H.shape)
    print("Eigenvalue array shape:", energies.shape)
    print("Eigenvector matrix shape:", wavefunctions.shape)


if __name__ == "__main__":
    main()
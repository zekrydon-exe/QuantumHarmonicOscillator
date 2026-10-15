import csv
from pathlib import Path

import numpy as np

from src.analytical_solution import analytical_energy
from src.hamiltonian import build_hamiltonian
from src.harmonic_oscillator import solve_harmonic_oscillator


def main():
    # Computational parameters
    q_min = -8.0
    q_max = 8.0

    n_points = 399
    n_states = 5

    # Create the results directory
    results_dir = Path("results")
    results_dir.mkdir(parents=True, exist_ok=True)

    # Solve the numerical problem
    q, dq, energies, wavefunctions = solve_harmonic_oscillator(
        q_min=q_min,
        q_max=q_max,
        n_points=n_points,
        n_states=n_states,
    )

    # Calculate analytical reference energies
    exact_energies = np.array([
        analytical_energy(n)
        for n in range(n_states)
    ])

    # Calculate energy errors
    absolute_errors = np.abs(energies - exact_energies)

    relative_errors = absolute_errors / np.abs(exact_energies)

    # Display the calculation parameters
    print("QUANTUM HARMONIC OSCILLATOR")
    print(f"Domain: [{q_min}, {q_max}]")
    print(f"Interior grid points: {n_points}")
    print(f"Grid spacing: {dq:.6f}")
    print(f"Number of states: {n_states}")

    # Display energy comparison
    print("\nENERGY COMPARISON")

    print(
        f"{'n':>3} "
        f"{'Numerical':>15} "
        f"{'Analytical':>15} "
        f"{'Absolute Error':>18} "
        f"{'Relative Error (%)':>20}"
    )

    for n in range(n_states):
        print(
            f"{n:3d} "
            f"{energies[n]:15.10f} "
            f"{exact_energies[n]:15.10f} "
            f"{absolute_errors[n]:18.10e} "
            f"{100 * relative_errors[n]:20.8f}"
        )

    # Validate wavefunction normalisation
    normalisation_integrals = np.array([
        np.trapezoid(
            np.abs(wavefunctions[:, n]) ** 2,
            q,
        )
        for n in range(n_states)
    ])

    print("\nWAVEFUNCTION NORMALISATION")
    print(f"{'n':>3} {'Normalisation Integral':>28}")

    for n, norm in enumerate(normalisation_integrals):
        print(f"{n:3d} {norm:28.12f}")

    # Construct the Hamiltonian for eigenpair validation
    H = build_hamiltonian(q, dq)

    # Convert continuum-scaled wavefunctions into discrete eigenvectors
    eigenvectors = wavefunctions * np.sqrt(dq)

    # Calculate eigenpair residuals
    residuals = (
        H @ eigenvectors
        - eigenvectors * energies
    )

    residual_norms = np.linalg.norm(
        residuals,
        axis=0,
    )

    print("\nEIGENPAIR VALIDATION")
    print(f"{'n':>3} {'Residual Norm':>22}")

    for n, residual in enumerate(residual_norms):
        print(f"{n:3d} {residual:22.12e}")

    print(
        "Maximum residual norm: "
        f"{np.max(residual_norms):.12e}"
    )

    # Export energy comparison
    energy_file = results_dir / "energy_comparison.csv"

    with energy_file.open("w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "quantum_number",
            "numerical_energy",
            "analytical_energy",
            "absolute_error",
            "relative_error_percent",
        ])

        for n in range(n_states):
            writer.writerow([
                n,
                energies[n],
                exact_energies[n],
                absolute_errors[n],
                relative_errors[n] * 100,
            ])

    print(f"\nEnergy results saved to: {energy_file}")

    # Export validation metrics
    validation_file = results_dir / "validation_metrics.csv"

    with validation_file.open("w", newline="") as file:
        writer = csv.writer(file)

        writer.writerow([
            "quantum_number",
            "normalisation_integral",
            "eigenpair_residual_norm",
        ])

        for n in range(n_states):
            writer.writerow([
                n,
                normalisation_integrals[n],
                residual_norms[n],
            ])

    print(f"Validation metrics saved to: {validation_file}")
    print("Calculation completed.")


if __name__ == "__main__":
    main()
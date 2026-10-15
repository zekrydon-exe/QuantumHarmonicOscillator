import numpy as np
import matplotlib.pyplot as plt

from src.harmonic_oscillator import solve_harmonic_oscillator
from src.analytical_solution import analytical_wavefunction


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

    for n in range(5):

        overlap = np.trapezoid(
            wavefunctions[:, n] * exact_wavefunctions[:, n],
            q,
        )

        if overlap < 0:
            wavefunctions[:, n] *= -1

    fig, axes = plt.subplots(
        3, 2,
        figsize=(12, 12),
        sharex=True,
    )

    axes = axes.flatten()

    for n in range(5):

        axes[n].plot(
            q,
            exact_wavefunctions[:, n],
            label="Analytical",
            linestyle="--",
        )

        axes[n].plot(
            q,
            wavefunctions[:, n],
            label="Numerical",
            linestyle="-",
        )

        axes[n].set_title(f"State n = {n}")
        axes[n].set_ylabel(r"$\psi_n(q)$")
        axes[n].grid(True, alpha=0.3)
        axes[n].legend()

    axes[4].set_xlabel(r"Dimensionless coordinate $q$")
    axes[5].axis("off")

    fig.suptitle(
        "Quantum Harmonic Oscillator: Wavefunctions",
        fontsize=15,
    )

    fig.tight_layout()

    fig.savefig(
        "results/wavefunctions.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    fig, axes = plt.subplots(
        3, 2,
        figsize=(12, 12),
        sharex=True,
    )

    axes = axes.flatten()

    for n in range(5):

        numerical_density = np.abs(wavefunctions[:, n])**2
        analytical_density = np.abs(exact_wavefunctions[:, n])**2

        axes[n].plot(
            q,
            analytical_density,
            label="Analytical",
            linestyle="--",
        )

        axes[n].plot(
            q,
            numerical_density,
            label="Numerical",
            linestyle="-",
        )

        axes[n].set_title(f"State n = {n}")
        axes[n].set_ylabel(r"$|\psi_n(q)|^2$")
        axes[n].grid(True, alpha=0.3)
        axes[n].legend()

    axes[4].set_xlabel(r"Dimensionless coordinate $q$")
    axes[5].axis("off")

    fig.suptitle(
        "Quantum Harmonic Oscillator: Probability Densities",
        fontsize=15,
    )

    fig.tight_layout()

    fig.savefig(
        "results/probability_densities.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()


if __name__ == "__main__":
    main()
import numpy as np
from scipy.linalg import eigh

from src.hamiltonian import create_grid, build_hamiltonian


def solve_harmonic_oscillator(
    q_min: float,
    q_max: float,
    n_points: int,
    n_states: int = 5,
):
    """
    Parameters
    q_min : float
        Left boundary of the computational domain.
    q_max : float
        Right boundary of the computational domain.
    n_points : int
        Number of interior grid points.
    n_states : int
        Number of lowest-energy states to calculate.

    Returns
    q : numpy.ndarray
        Interior spatial grid.
    dq : float
        Grid spacing.
    energies : numpy.ndarray
        Lowest numerical energy eigenvalues.
    wavefunctions : numpy.ndarray
        Continuum-scaled numerical wavefunctions.
        Each column represents one eigenstate.
    """

    if isinstance(n_states, (bool, np.bool_)) or not isinstance(
        n_states, (int, np.integer)
    ):
        raise TypeError("n_states must be an integer.")

    if isinstance(n_points, (bool, np.bool_)) or not isinstance(
        n_points, (int, np.integer)
    ):
        raise TypeError("n_points must be an integer.")

    if n_states < 1 or n_states > n_points:
        raise ValueError(
            "n_states must be between 1 and n_points."
        )

    q, dq = create_grid(q_min, q_max, n_points)

    H = build_hamiltonian(q, dq)

    energies, eigenvectors = eigh(
        H,
        subset_by_index=(0, n_states - 1),
    )

    # Convert discrete eigenvector normalisation into
    # the corresponding uniform-grid wavefunction scaling.
    wavefunctions = eigenvectors / np.sqrt(dq)

    return q, dq, energies, wavefunctions
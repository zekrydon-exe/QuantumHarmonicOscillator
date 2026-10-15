import numpy as np


def create_grid(q_min: float, q_max: float, n_points: int):
    """
    Parameters
    q_min : float
        Left boundary of the computational domain.
    q_max : float
        Right boundary of the computational domain.
    n_points : int
        Number of interior grid points.

    Returns
    q : numpy.ndarray
        Array containing the interior grid points.
    dq : float
        Uniform grid spacing.
    """

    if q_min >= q_max:
        raise ValueError("q_min must be smaller than q_max.")

    if n_points < 1:
        raise ValueError("n_points must be at least 1.")

    dq = (q_max - q_min) / (n_points + 1)

    q = np.linspace(
        q_min + dq,
        q_max - dq,
        n_points
    )

    return q, dq

def potential_energy(q: np.ndarray) -> np.ndarray:
    """
    Parameters
    q : numpy.ndarray
        Dimensionless spatial grid.

    Returns
    numpy.ndarray
        Potential-energy matrix.
    """

    potential = 0.5 * q**2

    return np.diag(potential)

def kinetic_energy(n_points: int, dq: float) -> np.ndarray:
    """
    Parameters
    n_points : int
        Number of interior grid points.
    dq : float
        Uniform grid spacing.

    Returns
    numpy.ndarray
        Kinetic-energy matrix.
    """

    if n_points < 1:
        raise ValueError("n_points must be at least 1.")

    if dq <= 0:
        raise ValueError("dq must be positive.")

    diagonal = np.full(n_points, 1.0 / dq**2)

    off_diagonal = np.full(
        n_points - 1,
        -1.0 / (2.0 * dq**2)
    )

    T = np.diag(diagonal)

    T += np.diag(off_diagonal, k=1)
    T += np.diag(off_diagonal, k=-1)

    return T

def build_hamiltonian(q: np.ndarray, dq: float) -> np.ndarray:
    """
    Parameters
    q : numpy.ndarray
        Interior dimensionless spatial grid.
    dq : float
        Uniform grid spacing.

    Returns
    numpy.ndarray
        Hamiltonian matrix.
    """

    n_points = len(q)

    T = kinetic_energy(n_points, dq)
    V = potential_energy(q)

    H = T + V

    return H
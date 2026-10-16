import numpy as np

def create_grid(q_min: float, q_max: float, n_points: int):
    """
    Create a uniformly spaced grid excluding the boundaries.

    Parameters
    ----------
    q_min : float
        Left boundary of the computational domain.
    q_max : float
        Right boundary of the computational domain.
    n_points : int
        Number of interior grid points.

    Returns
    -------
    q : numpy.ndarray
        Array containing the interior grid points.
    dq : float
        Uniform grid spacing.

    Raises
    ------
    TypeError
        If the boundaries are not real numbers or n_points is not an integer.
    ValueError
        If the boundaries are non-finite, improperly ordered, or
        n_points is less than 1.
    """

    if isinstance(q_min, (bool, np.bool_)) or not isinstance(
        q_min, (int, float, np.integer, np.floating)
    ):
        raise TypeError("q_min must be a real number.")

    if isinstance(q_max, (bool, np.bool_)) or not isinstance(
        q_max, (int, float, np.integer, np.floating)
    ):
        raise TypeError("q_max must be a real number.")

    if not np.isfinite(q_min) or not np.isfinite(q_max):
        raise ValueError("Grid boundaries must be finite.")

    if q_min >= q_max:
        raise ValueError("q_min must be smaller than q_max.")

    if isinstance(n_points, (bool, np.bool_)) or not isinstance(
        n_points, (int, np.integer)
    ):
        raise TypeError("n_points must be an integer.")

    if n_points < 1:
        raise ValueError("n_points must be at least 1.")

    dq = (q_max - q_min) / (n_points + 1)

    q = np.linspace(
        q_min + dq,
        q_max - dq,
        n_points,
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
    Construct the kinetic-energy matrix using the central
    finite-difference approximation.

    Parameters
    ----------
    n_points : int
        Number of interior grid points.
    dq : float
        Uniform grid spacing.

    Returns
    -------
    numpy.ndarray
        Kinetic-energy matrix.
    """

    if (
        isinstance(n_points, bool)
        or not isinstance(n_points, (int, np.integer))
        or n_points < 1
    ):
        raise ValueError("n_points must be a positive integer.")

    if (
        isinstance(dq, bool)
        or not isinstance(dq, (int, float, np.integer, np.floating))
        or not np.isfinite(dq)
        or dq <= 0
    ):
        raise ValueError("dq must be a finite positive number.")

    diagonal = np.full(n_points, 1.0 / dq**2)

    off_diagonal = np.full(
        n_points - 1,
        -1.0 / (2.0 * dq**2)
    )

    return (
        np.diag(diagonal)
        + np.diag(off_diagonal, k=1)
        + np.diag(off_diagonal, k=-1)
    )

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
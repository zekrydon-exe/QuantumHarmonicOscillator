import numpy as np
from scipy.special import eval_hermite, factorial


def analytical_energy(n: int) -> float:
    """
    Parameters
    n : int
        Quantum number, starting from zero.

    Returns
    float
        Exact dimensionless energy.
    """

    if not isinstance(n, (int, np.integer)) or isinstance(n, bool):
        raise TypeError("The quantum number n must be an integer.")

    if n < 0:
        raise ValueError("The quantum number n must be non-negative.")

    return n + 0.5


def analytical_wavefunction(
    q: np.ndarray,
    n: int
) -> np.ndarray:
    """
    Parameters
    q : numpy.ndarray
        Dimensionless spatial coordinates.
    n : int
        Quantum number, starting from zero.

    Returns
    numpy.ndarray
        Analytical wavefunction evaluated at q.
    """

    if not isinstance(n, (int, np.integer)) or isinstance(n, bool):
        raise TypeError("The quantum number n must be an integer.")

    if n < 0:
        raise ValueError("The quantum number n must be non-negative.")

    q = np.asarray(q, dtype=float)

    hermite_polynomial = eval_hermite(n, q)

    normalisation = (
        np.pi**(-0.25)
        / np.sqrt(2.0**n * factorial(n))
    )

    wavefunction = (
        normalisation
        * hermite_polynomial
        * np.exp(-q**2 / 2.0)
    )

    return wavefunction
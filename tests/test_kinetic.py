import numpy as np
import pytest

from src.hamiltonian import create_grid, kinetic_energy


def test_kinetic_matrix_has_correct_shape():
    q, dq = create_grid(-4.0, 4.0, 9)

    T = kinetic_energy(len(q), dq)

    assert T.shape == (9, 9)


def test_kinetic_matrix_is_symmetric():
    q, dq = create_grid(-4.0, 4.0, 9)

    T = kinetic_energy(len(q), dq)

    assert np.allclose(T, T.T)


def test_kinetic_matrix_diagonal_elements():
    q, dq = create_grid(-4.0, 4.0, 9)

    T = kinetic_energy(len(q), dq)

    expected_diagonal = 1.0 / dq**2

    assert np.allclose(np.diag(T), expected_diagonal)


def test_kinetic_matrix_off_diagonal_elements():
    q, dq = create_grid(-4.0, 4.0, 9)

    T = kinetic_energy(len(q), dq)

    expected_off_diagonal = -1.0 / (2.0 * dq**2)

    assert np.allclose(
        np.diag(T, k=1),
        expected_off_diagonal,
    )

    assert np.allclose(
        np.diag(T, k=-1),
        expected_off_diagonal,
    )


def test_kinetic_matrix_is_tridiagonal():
    q, dq = create_grid(-4.0, 4.0, 9)

    T = kinetic_energy(len(q), dq)

    tridiagonal_T = (
        np.diag(np.diag(T))
        + np.diag(np.diag(T, k=1), k=1)
        + np.diag(np.diag(T, k=-1), k=-1)
    )

    assert np.allclose(T, tridiagonal_T)


@pytest.mark.parametrize("n_points", [0, -1, 1.5, True, "5", None])
def test_kinetic_energy_rejects_invalid_n_points(n_points):
    with pytest.raises(ValueError):
        kinetic_energy(n_points, 0.1)


@pytest.mark.parametrize(
    "dq",
    [0, -0.1, np.nan, np.inf, -np.inf, True, "0.1", None],
)
def test_kinetic_energy_rejects_invalid_dq(dq):
    with pytest.raises(ValueError):
        kinetic_energy(5, dq)
import numpy as np

from src.hamiltonian import create_grid, potential_energy


def test_potential_matrix_has_correct_shape():
    q, dq = create_grid(-4.0, 4.0, 9)

    V = potential_energy(q)

    assert V.shape == (9, 9)


def test_potential_matrix_is_symmetric():
    q, dq = create_grid(-4.0, 4.0, 9)

    V = potential_energy(q)

    assert np.allclose(V, V.T)


def test_potential_energy_values_are_correct():
    q, dq = create_grid(-4.0, 4.0, 9)

    V = potential_energy(q)

    expected_diagonal = 0.5 * q**2

    assert np.allclose(np.diag(V), expected_diagonal)


def test_potential_matrix_is_diagonal():
    q, dq = create_grid(-4.0, 4.0, 9)

    V = potential_energy(q)

    assert np.allclose(V, np.diag(np.diag(V)))


def test_potential_energy_is_nonnegative_and_minimum_at_origin():
    q, dq = create_grid(-4.0, 4.0, 9)

    V = potential_energy(q)

    potential_values = np.diag(V)

    assert np.all(potential_values >= 0.0)
    assert np.isclose(potential_values[4], 0.0)
    assert np.isclose(np.min(potential_values), 0.0)
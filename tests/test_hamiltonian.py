import numpy as np

from src.hamiltonian import (
    create_grid,
    kinetic_energy,
    potential_energy,
    build_hamiltonian,
)


def test_hamiltonian_has_correct_shape():
    q, dq = create_grid(-4.0, 4.0, 9)

    H = build_hamiltonian(q, dq)

    assert H.shape == (9, 9)


def test_hamiltonian_is_symmetric():
    q, dq = create_grid(-4.0, 4.0, 9)

    H = build_hamiltonian(q, dq)

    assert np.allclose(H, H.T)


def test_hamiltonian_equals_kinetic_plus_potential():
    q, dq = create_grid(-4.0, 4.0, 9)

    T = kinetic_energy(len(q), dq)
    V = potential_energy(q)
    H = build_hamiltonian(q, dq)

    assert np.allclose(H, T + V)


def test_hamiltonian_diagonal_elements():
    q, dq = create_grid(-4.0, 4.0, 9)

    H = build_hamiltonian(q, dq)

    expected_diagonal = 1.0 / dq**2 + 0.5 * q**2

    assert np.allclose(np.diag(H), expected_diagonal)


def test_hamiltonian_off_diagonal_elements():
    q, dq = create_grid(-4.0, 4.0, 9)

    H = build_hamiltonian(q, dq)

    expected_off_diagonal = -1.0 / (2.0 * dq**2)

    assert np.allclose(
        np.diag(H, k=1),
        expected_off_diagonal,
    )

    assert np.allclose(
        np.diag(H, k=-1),
        expected_off_diagonal,
    )


def test_hamiltonian_is_tridiagonal():
    q, dq = create_grid(-4.0, 4.0, 9)

    H = build_hamiltonian(q, dq)

    tridiagonal_H = (
        np.diag(np.diag(H))
        + np.diag(np.diag(H, k=1), k=1)
        + np.diag(np.diag(H, k=-1), k=-1)
    )

    assert np.allclose(H, tridiagonal_H)
import numpy as np

from src.hamiltonian import (
    create_grid,
    kinetic_energy,
    potential_energy,
    build_hamiltonian,
)


q, dq = create_grid(-4.0, 4.0, 9)

T = kinetic_energy(len(q), dq)
V = potential_energy(q)
H = build_hamiltonian(q, dq)

print("Hamiltonian matrix:")
print(H)

print("\nMatrix shape:", H.shape)

print("Hamiltonian is symmetric:", np.allclose(H, H.T))

print(
    "Hamiltonian equals T + V:",
    np.allclose(H, T + V)
)

print(
    "Off-diagonal elements correct:",
    np.allclose(
        np.diag(H, k=1),
        -1.0 / (2.0 * dq**2)
    )
)

print(
    "Matrix is tridiagonal:",
    np.allclose(H, np.diag(np.diag(H)))
    or np.allclose(
        H,
        np.diag(np.diag(H))
        + np.diag(np.diag(H, k=1), k=1)
        + np.diag(np.diag(H, k=-1), k=-1)
    )
)
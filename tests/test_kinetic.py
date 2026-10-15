import numpy as np

from src.hamiltonian import create_grid, kinetic_energy


q, dq = create_grid(-4.0, 4.0, 9)

T = kinetic_energy(len(q), dq)

print("Kinetic-energy matrix:")
print(T)

print("\nMatrix shape:", T.shape)
print("Main diagonal:", np.diag(T))

print("Upper diagonal:", np.diag(T, k=1))
print("Lower diagonal:", np.diag(T, k=-1))

print("\nIs the matrix symmetric?", np.allclose(T, T.T))

expected_diagonal = 1.0 / dq**2
expected_off_diagonal = -1.0 / (2.0 * dq**2)

print(
    "Diagonal elements correct?",
    np.allclose(np.diag(T), expected_diagonal)
)

print(
    "Off-diagonal elements correct?",
    np.allclose(np.diag(T, k=1), expected_off_diagonal)
    and np.allclose(np.diag(T, k=-1), expected_off_diagonal)
)
import numpy as np

from src.hamiltonian import create_grid, potential_energy


q, dq = create_grid(-4.0, 4.0, 9)

V = potential_energy(q)

print("Potential-energy matrix:")
print(V)

print("\nMatrix shape:", V.shape)
print("Diagonal elements:", np.diag(V))

print("\nIs the matrix symmetric?", np.allclose(V, V.T))
print("Are all potential energies non-negative?",
      np.all(np.diag(V) >= 0))

print("Is the potential minimum at the origin?",
      np.isclose(V[4, 4], 0.0))
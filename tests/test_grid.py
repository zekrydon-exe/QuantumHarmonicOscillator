from src.hamiltonian import create_grid


q_min = -4.0
q_max = 4.0
n_points = 9

q, dq = create_grid(q_min, q_max, n_points)

print("Number of grid points:", len(q))
print("Grid spacing:", dq)
print("Grid coordinates:", q)
print("Grid includes origin:", 0.0 in q)
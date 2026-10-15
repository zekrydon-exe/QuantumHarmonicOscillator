# Quantum Harmonic Oscillator

A computational quantum mechanics project that solves the one-dimensional
quantum harmonic oscillator using the finite-difference method and validates
the numerical results against analytical solutions.

## 1. Overview

The quantum harmonic oscillator is a fundamental model in quantum mechanics.
It provides an exactly solvable system for studying quantised energy levels
and wavefunctions.

This project numerically solves the dimensionless time-independent
Schrödinger equation, calculates the lowest energy eigenvalues and
eigenfunctions, and compares them with the exact analytical solutions.

The project also investigates numerical accuracy through grid convergence,
domain convergence, wavefunction validation, and eigenpair residual analysis.

## 2. Scientific Model

The dimensionless Schrödinger equation is

\[
-\frac{1}{2}\frac{d^2\psi(q)}{dq^2}
+\frac{1}{2}q^2\psi(q)
=\varepsilon\psi(q)
\]

The exact dimensionless energy eigenvalues are

\[
\varepsilon_n=n+\frac{1}{2},
\qquad n=0,1,2,\ldots
\]

The numerical solution is obtained by discretising the second derivative
using a central finite-difference approximation and constructing a
tridiagonal Hamiltonian matrix.

The resulting matrix eigenvalue problem is solved using SciPy.

For the production calculation, the parameters are:

| Parameter | Value |
|---|---:|
| Dimensionless domain | [-8, 8] |
| Number of interior grid points | 399 |
| Grid spacing | 0.04 |
| Number of eigenstates | 5 |
| Finite-difference method | Central difference |
| Eigenvalue solver | `scipy.linalg.eigh` |

The detailed mathematical derivation is available in
[`theory/theory.md`](theory/theory.md).

## 3. Project Structure

```text
02_Quantum_Harmonic_Oscillator/
├── src/
│   ├── __init__.py
│   ├── analytical_solution.py
│   ├── hamiltonian.py
│   └── harmonic_oscillator.py
├── scripts/
│   ├── __init__.py
│   ├── solve_eigenvalues.py
│   ├── run_calculation.py
│   ├── plot_wavefunctions.py
│   ├── validate_wavefunctions.py
│   ├── convergence_study.py
│   ├── domain_convergence.py
│   └── wavefunction_convergence.py
├── tests/
├── theory/
│   └── theory.md
├── report/
├── results/
│   ├── energy_comparison.csv
│   ├── validation_metrics.csv
│   ├── wavefunctions.png
│   └── probability_densities.png
├── requirements.txt
└── README.md
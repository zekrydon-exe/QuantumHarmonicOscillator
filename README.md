# Quantum Harmonic Oscillator

A computational quantum mechanics project that solves the one-dimensional quantum harmonic oscillator using the finite-difference method and validates the numerical results against analytical solutions.

The project calculates low-lying energy eigenvalues and eigenfunctions, compares them with exact reference solutions, and investigates numerical accuracy through grid convergence, domain convergence, wavefunction validation, and eigenpair residual analysis.

## 1. Scientific Model

The dimensionless time-independent Schrödinger equation is

$$
-\frac{1}{2}\frac{d^2\psi(q)}{dq^2}
+\frac{1}{2}q^2\psi(q)
=\varepsilon\psi(q).
$$

The exact dimensionless energy eigenvalues are

$$
\varepsilon_n=n+\frac{1}{2},
\qquad n=0,1,2,\ldots
$$

The numerical solution is obtained by approximating the second derivative using a central finite-difference formula and constructing a tridiagonal Hamiltonian matrix. The resulting real, symmetric matrix eigenvalue problem is solved using `scipy.linalg.eigh`.

The detailed mathematical model and numerical derivation are documented in [`theory/theory.md`](theory/theory.md).

## 2. Production Configuration

The main calculation uses the following parameters:

| Parameter | Value |
|---|---:|
| Dimensionless domain | $[-8,8]$ |
| Number of interior grid points | 399 |
| Grid spacing | 0.04 |
| Number of eigenstates | 5 |
| Finite-difference approximation | Central difference |
| Eigenvalue solver | `scipy.linalg.eigh` |

The dimensionless energy values are reported in units of $\hbar\omega$, where $\omega$ is the oscillator angular frequency.

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
│   └── technical_report.md
├── results/
│   ├── energy_comparison.csv
│   ├── validation_metrics.csv
│   ├── wavefunctions.png
│   └── probability_densities.png
├── requirements.txt
└── README.md
```

The `src/` directory contains the core scientific implementation. The `scripts/` directory contains executable calculations, plotting, validation, and convergence studies. The `tests/` directory contains automated tests, `theory/` contains the mathematical documentation, `report/` contains the technical report, and `results/` stores generated numerical data and figures.

Python cache directories and virtual environments are intentionally excluded from version control.

## 4. Requirements

The project uses Python and these direct dependencies:

- **NumPy** for arrays and numerical operations.
- **SciPy** for linear algebra and Hermite polynomials.
- **Matplotlib** for scientific visualisation.
- **pytest** for automated testing.

The dependency versions are recorded in [`requirements.txt`](requirements.txt).

## 5. Installation

From the project directory, create a virtual environment if you do not already have one:

```bash
python3 -m venv .venv
```

Activate it on Linux:

```bash
source .venv/bin/activate
```

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

## 6. Run the Main Calculation

Run the calculation from the project root:

```bash
python -m scripts.run_calculation
```

The script calculates the lowest five energy eigenvalues and wavefunctions, compares numerical energies with the analytical reference values, evaluates wavefunction normalisation, calculates eigenpair residual norms, and exports numerical data.

The generated CSV files are:

- `results/energy_comparison.csv`
- `results/validation_metrics.csv`

The energy CSV contains the quantum number, numerical energy, analytical energy, absolute error, and relative error percentage. The validation CSV contains normalisation integrals and eigenpair residual norms.

## 7. Generate Plots

Run:

```bash
python -m scripts.plot_wavefunctions
```

The script generates plots of the numerical wavefunctions and their probability densities in the `results/` directory.

The probability density for state $n$ is

$$
\rho_n(q)=|\psi_n(q)|^2.
$$

The wavefunction can take positive or negative values, while its probability density is non-negative.

## 8. Validation and Convergence

The project includes several independent checks.

### Energy validation

Numerical energies are compared with the exact values

$$
\varepsilon_n=n+\frac{1}{2}.
$$

Absolute and relative errors quantify the differences between numerical and analytical energies.

### Wavefunction validation

The numerical wavefunctions are checked for normalisation and compared with analytical wavefunctions using their overlap and relative $L^2$ error. An overall sign difference between two eigenfunctions is physically immaterial, so the comparison aligns their signs before evaluating the difference.

### Eigenpair residuals

For a numerical eigenpair, the residual is

$$
\mathbf{r}_n
=
\mathbf{H}\boldsymbol{\psi}_n
-
\varepsilon_n\boldsymbol{\psi}_n.
$$

A small residual indicates that the computed eigenpair satisfies the discretised matrix eigenvalue equation accurately. It does not, by itself, establish that the discretised result is equally accurate relative to the continuous Schrödinger equation.

### Grid convergence

Run:

```bash
python -m scripts.convergence_study
```

This study varies the grid resolution over the fixed domain $[-4,4]$. The observed ground-state energy convergence is consistent with the second-order accuracy expected from the central finite-difference approximation.

### Domain convergence

Run:

```bash
python -m scripts.domain_convergence
```

This study varies the domain while maintaining a target grid spacing of approximately 0.04. The energies for the first five states become effectively unchanged between the domains $[-8,8]$ and $[-10,10]$ at this resolution.

### Wavefunction convergence

Run:

```bash
python -m scripts.wavefunction_convergence
```

This study compares numerical wavefunctions with analytical solutions as the grid is refined.

### Additional eigenvalue calculation

The project also includes:

```bash
python -m scripts.solve_eigenvalues
```

This script provides a separate entry point for calculating the low-lying energy eigenvalues.

## 9. Automated Tests

Run the test suite from the project root:

```bash
python -m pytest tests/
```

The recorded test run completed with **61 tests passing**. The tests check selected behaviour of the analytical solution and numerical implementation. Passing tests support implementation reliability, but they do not replace physical validation and convergence analysis.

## 10. Recorded Numerical Results

For the production configuration $q\in[-8,8]$, with 399 interior grid points and grid spacing 0.04, the calculated energies were:

| Quantum number $n$ | Numerical energy | Exact energy | Absolute error |
|---:|---:|---:|---:|
| 0 | 0.4999499950 | 0.5000000000 | $5.0005\times10^{-5}$ |
| 1 | 1.4997499550 | 1.5000000000 | $2.5005\times10^{-4}$ |
| 2 | 2.4993498249 | 2.5000000000 | $6.5018\times10^{-4}$ |
| 3 | 3.4987495446 | 3.5000000000 | $1.2505\times10^{-3}$ |
| 4 | 4.4979490539 | 4.5000000000 | $2.0509\times10^{-3}$ |

The numerical wavefunctions have normalisation integrals equal to 1.000000000000 to the displayed precision. Their relative $L^2$ errors range from approximately $7.36\times10^{-5}$ for the ground state to $1.10\times10^{-3}$ for the fourth excited state.

The maximum eigenpair residual norm recorded for this calculation is approximately $1.81\times10^{-13}$. This indicates high accuracy in solving the discretised matrix problem, while the larger energy discrepancies relative to the exact spectrum primarily reflect finite-difference discretisation.

These results are generated by the project scripts and are stored in the `results/` directory.

## 11. Reproducibility

To reproduce the calculation:

1. Use a Python environment compatible with the dependencies in `requirements.txt`.
2. Install the dependencies.
3. Run `python -m scripts.run_calculation` from the project root.
4. Inspect the generated CSV files.
5. Run the plotting, validation, and convergence scripts as needed.
6. Run `python -m pytest tests/` to execute the automated tests.

The production parameters, mathematical formulation, numerical method, and validation procedures are documented so that the calculation can be repeated and independently assessed.

## 12. Possible Extensions

Possible future extensions include:

- Increasing grid resolution to reduce discretisation error.
- Comparing finite differences with alternative numerical methods.
- Calculating additional quantum states.
- Extending the framework to other one-dimensional potentials.
- Applying the same numerical workflow to other quantum-mechanical systems.

## 13. Project Context

This project is part of a broader computational quantum mechanics sequence. It builds on the numerical treatment of a particle in a one-dimensional infinite potential well and develops the methods needed to solve and validate more general Schrödinger eigenvalue problems.
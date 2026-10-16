# Numerical Study of the Quantum Harmonic Oscillator

## 1. Abstract

This project investigates the quantum harmonic oscillator by numerically solving the time-independent Schrödinger equation using the finite-difference method. The dimensionless Hamiltonian is discretised on a uniformly spaced one-dimensional grid, producing a real symmetric matrix whose eigenvalues and eigenvectors approximate the energy levels and wavefunctions.

The numerical energies are compared with the analytical energy spectrum. The computed wavefunctions are examined through normalisation, orthogonality, parity and probability-density properties. Numerical reliability is further assessed using eigenpair residuals, automated tests and grid-convergence studies.

For the standard calculation, the computational domain is $q\in[-8,8]$, with 399 interior grid points and grid spacing $\Delta q=0.04$. The five lowest dimensionless energy levels are calculated. The numerical results agree closely with the analytical spectrum, with remaining discrepancies primarily associated with finite-difference discretisation.

## 2. Introduction

The quantum harmonic oscillator is a fundamental model in quantum mechanics. It describes a particle subject to a potential energy that increases quadratically with displacement from an equilibrium position. Despite its simple potential, the system provides a useful framework for studying quantisation, stationary states, wavefunctions, probability densities and numerical solutions of the Schrödinger equation.

The harmonic oscillator is particularly suitable for computational investigation because its energy eigenvalues and eigenfunctions are known analytically. These exact solutions provide reference values against which numerical methods can be evaluated.

### 2.1 Objectives

The principal objectives are:

1. Formulate the one-dimensional quantum harmonic oscillator in dimensionless form.
2. Discretise the Schrödinger equation using a finite-difference approximation.
3. Construct the corresponding Hamiltonian matrix.
4. Calculate the lowest energy eigenvalues and eigenvectors numerically.
5. Compare the numerical energy spectrum with the analytical solution.
6. Examine normalisation, orthogonality and parity of the numerical wavefunctions.
7. Investigate numerical accuracy through residual calculations and convergence studies.
8. Develop a reproducible Python implementation supported by automated tests.

## 3. Theoretical Background

### 3.1 Harmonic oscillator potential

For a particle of mass $m$ moving in a harmonic potential, the potential energy is

$$
V(x)=\frac{1}{2}m\omega^2x^2,
$$

where $\omega$ is the angular frequency of the oscillator.

The time-independent Schrödinger equation is

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}
+\frac{1}{2}m\omega^2x^2\psi(x)
=E\psi(x).
$$

Here, $\psi(x)$ is the stationary-state wavefunction and $E$ is the corresponding energy eigenvalue.

### 3.2 Dimensionless formulation

Introduce the dimensionless coordinate

$$
q=\sqrt{\frac{m\omega}{\hbar}}\,x
$$

and dimensionless energy

$$
\varepsilon=\frac{E}{\hbar\omega}.
$$

The Schrödinger equation becomes

$$
-\frac{1}{2}\frac{d^2\psi(q)}{dq^2}
+\frac{1}{2}q^2\psi(q)
=\varepsilon\psi(q).
$$

The dimensionless Hamiltonian is

$$
\hat{H}=-\frac{1}{2}\frac{d^2}{dq^2}+\frac{1}{2}q^2,
$$

and its eigenvalue equation is

$$
\hat{H}\psi_n(q)=\varepsilon_n\psi_n(q).
$$

This representation removes the explicit dependence on particle mass and oscillator frequency from the numerical problem.

### 3.3 Analytical energy spectrum

The exact dimensionless energy eigenvalues are

$$
\varepsilon_n=n+\frac{1}{2},\qquad n=0,1,2,\ldots
$$

The dimensional energies are

$$
E_n=\hbar\omega\left(n+\frac{1}{2}\right).
$$

The ground-state energy is $\varepsilon_0=\frac{1}{2}$, and consecutive dimensionless levels are separated by

$$
\varepsilon_{n+1}-\varepsilon_n=1.
$$

This equally spaced spectrum provides a reference for validating the numerical eigenvalues.

### 3.4 Analytical wavefunctions

The normalised analytical eigenfunctions are

$$
\psi_n(q)=\frac{H_n(q)e^{-q^2/2}}{\pi^{1/4}\sqrt{2^n n!}},
$$

where $H_n(q)$ is the physicists' Hermite polynomial of degree $n$.

For the ground state, $H_0(q)=1$, giving

$$
\psi_0(q)=\frac{1}{\pi^{1/4}}e^{-q^2/2}.
$$

These analytical wavefunctions provide reference solutions for examining the numerical eigenstates.

### 3.5 Normalisation and orthogonality

A physically admissible bound-state wavefunction satisfies

$$
\int_{-\infty}^{\infty}|\psi_n(q)|^2\,dq=1.
$$

Different energy eigenstates are orthogonal:

$$
\int_{-\infty}^{\infty}\psi_m^*(q)\psi_n(q)\,dq=\delta_{mn},
$$

where $\delta_{mn}$ is the Kronecker delta. On a finite numerical grid, these integrals are approximated over the computational domain, so accuracy depends on grid spacing and domain size.

## 4. Numerical Methodology

### 4.1 Computational grid

The coordinate $q$ is replaced by a uniformly spaced grid containing $N$ interior points between $q_{\min}$ and $q_{\max}$. The spacing is

$$
\Delta q=\frac{q_{\max}-q_{\min}}{N+1}.
$$

The interior coordinates are

$$
q_i=q_{\min}+i\Delta q,\qquad i=1,2,\ldots,N.
$$

The boundary points are excluded from the array of interior points, corresponding to zero wavefunction values at the computational boundaries. For the standard calculation,

$$
q_{\min}=-8,\qquad q_{\max}=8,\qquad N=399,
$$

so

$$
\Delta q=\frac{16}{400}=0.04.
$$

### 4.2 Finite-difference approximation

The second derivative is approximated by the central finite-difference formula

$$
\left.\frac{d^2\psi}{dq^2}\right|_{q_i}
\approx\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{(\Delta q)^2}.
$$

Substitution into the dimensionless Schrödinger equation gives

$$
-\frac{1}{2}\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{(\Delta q)^2}
+\frac{q_i^2}{2}\psi_i
=\varepsilon\psi_i.
$$

Rearranging,

$$
-\frac{\psi_{i-1}}{2(\Delta q)^2}
+\left[\frac{1}{(\Delta q)^2}+\frac{q_i^2}{2}\right]\psi_i
-\frac{\psi_{i+1}}{2(\Delta q)^2}
=\varepsilon\psi_i.
$$

This expression defines the matrix elements of the discretised Hamiltonian.

### 4.3 Kinetic-energy matrix

The kinetic-energy operator is

$$
\hat{T}=-\frac{1}{2}\frac{d^2}{dq^2}.
$$

Its finite-difference matrix has diagonal elements

$$
T_{ii}=\frac{1}{(\Delta q)^2}
$$

and nearest-neighbour off-diagonal elements

$$
T_{i,i+1}=T_{i+1,i}=-\frac{1}{2(\Delta q)^2}.
$$

All other elements are zero. The kinetic-energy matrix is therefore tridiagonal and symmetric.

### 4.4 Potential-energy matrix

The dimensionless potential energy is

$$
V(q)=\frac{q^2}{2}.
$$

In the position-grid representation, the potential matrix is diagonal:

$$
V_{ij}=
\begin{cases}
\dfrac{q_i^2}{2}, & i=j,\\
0, & i\ne j.
\end{cases}
$$

### 4.5 Hamiltonian matrix

The complete discretised Hamiltonian is

$$
H=T+V.
$$

Its diagonal elements are

$$
H_{ii}=\frac{1}{(\Delta q)^2}+\frac{q_i^2}{2},
$$

while the nearest-neighbour elements are

$$
H_{i,i+1}=H_{i+1,i}=-\frac{1}{2(\Delta q)^2}.
$$

The resulting matrix is real, symmetric and tridiagonal.

### 4.6 Numerical eigenvalue problem

The discretised Schrödinger equation becomes

$$
H\mathbf{v}_n=\varepsilon_n\mathbf{v}_n,
$$

where eigenvalues approximate dimensionless energies and eigenvectors represent the corresponding discrete states. The implementation uses the symmetric eigensolver `scipy.linalg.eigh` and requests only the required lowest-energy states.

The discrete eigenvectors are Euclidean-normalised. The continuum-scaled wavefunctions are obtained from

$$
\psi_n(q_i)=\frac{v_{i,n}}{\sqrt{\Delta q}}.
$$

Thus, the discrete normalisation condition

$$
\sum_i|v_{i,n}|^2=1
$$

is consistent with the uniform-grid approximation

$$
\sum_i|\psi_n(q_i)|^2\Delta q\approx1.
$$

## 5. Software Implementation

The project uses Python, NumPy, SciPy, Matplotlib and pytest. Its modules construct the Hamiltonian, calculate analytical solutions, solve the eigenvalue problem, perform convergence studies and validate results.

### 5.1 Main source modules

**`src/hamiltonian.py`**

- `create_grid()` constructs the uniformly spaced interior grid and returns its spacing.
- `kinetic_energy()` constructs the finite-difference kinetic-energy matrix.
- `potential_energy()` constructs the diagonal harmonic potential matrix.
- `build_hamiltonian()` combines the kinetic and potential matrices.

Input validation rejects invalid grid parameters and invalid kinetic-energy matrix inputs.

**`src/analytical_solution.py`**

- `analytical_energy(n)` calculates the exact dimensionless energy.
- `analytical_wavefunction(q, n)` evaluates the analytical wavefunction using Hermite polynomials.

**`src/harmonic_oscillator.py`**

`solve_harmonic_oscillator()` constructs the grid and Hamiltonian, calculates the requested lowest-energy eigenstates, and scales wavefunctions for integration on the coordinate grid.

### 5.2 Supporting scripts

The project includes scripts for calculations, plotting, validation and convergence:

- `scripts/convergence_study.py`: investigates energy dependence on grid resolution.
- `scripts/domain_convergence.py`: investigates the effect of computational-domain size.
- `scripts/wavefunction_convergence.py`: investigates wavefunction convergence.
- `scripts/validate_wavefunctions.py`: supports wavefunction validation.
- `scripts/plot_wavefunctions.py`: generates wavefunction plots.
- `scripts/solve_eigenvalues.py`: calculates low-lying eigenvalues.
- `scripts/run_calculation.py`: runs the main calculation workflow.

Generated outputs include energy-comparison and validation CSV files, together with wavefunction and probability-density plots.

### 5.3 Automated testing

The pytest suite covers grid construction and input validation, matrix properties, analytical solutions, numerical energies, wavefunction properties, and physical consistency. Automated tests help identify implementation errors and make subsequent modifications easier to verify.

## 6. Results and Discussion

### 6.1 Energy eigenvalues

The five lowest dimensionless energy levels were calculated for domain $[-8,8]$, 399 interior points and spacing $\Delta q=0.04$. The exact spectrum is

$$
\varepsilon_n=n+\frac{1}{2}.
$$

**Table 1. Numerical and analytical energy eigenvalues.**

| Quantum number $n$ | Numerical energy | Analytical energy | Absolute error |
|---:|---:|---:|---:|
| 0 | 0.4999499950 | 0.5000000000 | $5.0005\times10^{-5}$ |
| 1 | 1.4997499550 | 1.5000000000 | $2.5005\times10^{-4}$ |
| 2 | 2.4993498249 | 2.5000000000 | $6.5018\times10^{-4}$ |
| 3 | 3.4987495446 | 3.5000000000 | $1.2505\times10^{-3}$ |
| 4 | 4.4979490539 | 4.5000000000 | $2.0509\times10^{-3}$ |

The numerical energies are close to the analytical values. The error generally becomes larger for higher states in this calculation, reflecting their greater sensitivity to the finite-difference approximation. Error magnitudes depend on grid spacing, domain size and state.

### 6.2 Energy-level spacing

The analytical levels are equally spaced:

$$
\varepsilon_{n+1}-\varepsilon_n=1.
$$

The numerical spectrum reproduces this behaviour approximately. Deviations arise from the discretised Hamiltonian and are checked by the project's physics-validation tests.

### 6.3 Wavefunctions and probability densities

The probability density for state $n$ is

$$
\rho_n(q)=|\psi_n(q)|^2.
$$

Expected qualitative properties include:

- The ground-state wavefunction has even parity.
- Successive eigenstates alternate between even and odd parity.
- Probability density is non-negative and symmetric about the origin.
- Higher-energy states have additional nodes.

These properties provide checks on the numerical eigenstates and their physical interpretation.

### 6.4 Interpretation of energy errors

The central finite-difference approximation to the second derivative has a leading truncation error of order $(\Delta q)^2$ for sufficiently smooth solutions. At fixed domain and sufficiently small spacing, numerical eigenvalues are therefore expected to exhibit second-order convergence under usual conditions.

The finite computational domain introduces a separate error because the wavefunction is constrained to vanish at the boundaries. For domain $[-8,8]$, the low-lying oscillator states are strongly localised away from the boundaries, so this effect is small for the states examined. The convergence studies investigate discretisation and boundary effects separately.

## 7. Numerical Validation

### 7.1 Energy comparison

The numerical energies are compared with the exact spectrum

$$
\varepsilon_n=n+\frac{1}{2}.
$$

The absolute error is

$$
\Delta\varepsilon_n=
\left|\varepsilon_n^{\mathrm{numerical}}-\varepsilon_n^{\mathrm{exact}}\right|.
$$

### 7.2 Wavefunction normalisation

The continuum normalisation condition is

$$
\int_{-\infty}^{\infty}|\psi_n(q)|^2\,dq=1.
$$

The numerical integral is evaluated over the finite computational domain:

$$
I_n\approx\int_{q_{\min}}^{q_{\max}}|\psi_n(q)|^2\,dq.
$$

The computed normalisation integrals are approximately one for the states tested.

### 7.3 Orthogonality

The numerical overlap is evaluated as

$$
S_{mn}\approx\int_{q_{\min}}^{q_{\max}}\psi_m(q)\psi_n(q)\,dq.
$$

For distinct states, the expected overlap is zero. The automated tests verify that overlaps between the first five numerical wavefunctions remain below the specified tolerance of $10^{-3}$.

### 7.4 Eigenpair residuals

A numerical eigenpair should satisfy

$$
H\mathbf{v}_n\approx\varepsilon_n\mathbf{v}_n.
$$

The residual vector and its Euclidean norm are

$$
\mathbf{r}_n=H\mathbf{v}_n-\varepsilon_n\mathbf{v}_n,
\qquad
r_n=\|\mathbf{r}_n\|_2.
$$

For the standard calculation, the maximum residual norm is approximately

$$
\max_n r_n\approx1.81\times10^{-13}.
$$

This indicates that the computed eigenpairs satisfy the discretised matrix eigenvalue problem to high numerical precision. A small residual alone does not establish accuracy relative to the exact continuous problem; energy comparisons and convergence studies are also necessary.

### 7.5 Automated test results

The automated suite checks grid dimensions, spacing and boundary exclusion; invalid inputs; Hamiltonian symmetry and tridiagonal structure; kinetic- and potential-energy matrix elements; analytical solutions; energy accuracy; normalisation and orthogonality; eigenpair residuals; energy spacing; parity; probability-density properties; and solver input validation.

The recorded test run completed with **61 passing tests**. This confirms that all tests collected in that run passed, but does not imply that every possible input or numerical configuration has been tested.

## 8. Convergence Analysis

Convergence analysis examines how numerical results change as the grid is refined or the computational domain is enlarged. These studies address distinct errors from discretisation and finite-domain boundaries.

### 8.1 Grid-spacing convergence

The grid-convergence study uses fixed domain $[-4,4]$ and interior-point counts

$$
N\in\{49,99,199,399\}.
$$

For state $n$, the energy error is

$$
e_n(\Delta q)=
\left|\varepsilon_n^{\mathrm{numerical}}-\varepsilon_n^{\mathrm{exact}}\right|.
$$

For spacings $\Delta q_1$ and $\Delta q_2$, the observed order is estimated by

$$
p=\frac{\ln\left(e_n(\Delta q_1)/e_n(\Delta q_2)\right)}
{\ln(\Delta q_1/\Delta q_2)}.
$$

If the leading error is proportional to $(\Delta q)^2$, the observed order should approach two as the grid is refined. The estimate can be unreliable when errors are very small or other error sources dominate. The domain $[-4,4]$ is retained for this study, although higher states can be more sensitive to its boundaries. The order should therefore be taken from the script output rather than assumed to equal two.

### 8.2 Computational-domain convergence

A separate study investigates changing the computational boundaries while maintaining appropriate grid resolution. Increasing the domain moves the imposed zero-wavefunction boundaries farther into regions where bound-state wavefunctions are small.

The standard calculation uses $[-8,8]$. Checks at spacing approximately $0.04$ found the low-lying energies for domains $[-8,8]$ and $[-10,10]$ effectively identical at the reported precision. This supports the standard domain for the states examined, but the conclusion should not automatically be extended to arbitrarily high-energy states.

### 8.3 Wavefunction convergence

Energy convergence alone does not fully establish eigenstate accuracy. Wavefunction comparisons across grids must account for different coordinate grids and the arbitrary overall sign of each real eigenvector. Analytical wavefunctions provide an additional reference. Energy, domain and wavefunction convergence studies provide complementary assessments of numerical accuracy.

## 9. Limitations and Possible Improvements

### 9.1 Finite computational domain

The calculation uses a finite interval instead of the entire real line and constrains the wavefunction to vanish at its boundaries. Increasing the domain can reduce this error, but increases computational cost if spacing is held constant.

### 9.2 Finite-difference discretisation

The central finite-difference formula introduces error that depends on grid spacing. Higher-order finite-difference schemes could be investigated to determine whether comparable accuracy can be achieved with fewer grid points.

### 9.3 Wavefunction comparison

Numerical eigenvectors are defined only up to an overall sign for a real-valued, non-degenerate state. Comparisons should account for this sign freedom. Wavefunctions on different grids should be interpolated onto a common grid or compared using another appropriate numerical measure.

### 9.4 Numerical integration

Normalisation and orthogonality use numerical integration over the finite domain. Their accuracy depends on grid resolution, integration method and domain size.

### 9.5 Computational efficiency

The implementation constructs a dense Hamiltonian even though the finite-difference Hamiltonian is tridiagonal. For the grid sizes used, the dense implementation is convenient and straightforward to validate. For substantially larger grids, a sparse representation and suitable sparse eigensolver could reduce memory use and computational cost.

### 9.6 Further investigations

Possible extensions include:

- Comparing dense and sparse eigensolvers.
- Investigating higher-order finite-difference schemes.
- Quantifying wavefunction errors against analytical solutions.
- Extending the calculation to other one-dimensional potentials.
- Investigating how numerical error changes with quantum number.
- Separating discretisation and boundary effects in a more systematic error analysis.

## 10. Conclusion

This project demonstrates a numerical solution of the one-dimensional quantum harmonic oscillator using the finite-difference method. The dimensionless Schrödinger equation was discretised on a uniform grid, and the resulting real symmetric Hamiltonian was solved using a numerical eigensolver.

The five lowest numerical energies were compared with the analytical spectrum for domain $[-8,8]$ and 399 interior grid points. The calculated wavefunctions were examined through normalisation, orthogonality, parity and probability-density properties. Eigenpair residuals were evaluated, with the maximum residual norm approximately $1.81\times10^{-13}$ for the five states considered.

Automated tests check mathematical properties and invalid input parameters. Grid-spacing, computational-domain and wavefunction-convergence studies provide additional methods for assessing reliability. The project illustrates why numerical validation requires more than a successful eigensolver: analytical comparisons, physical consistency checks and convergence studies are necessary to assess accuracy.

## 11. Reproducibility

The project uses Python, NumPy, SciPy, Matplotlib and pytest.

### 11.1 Environment setup

Create and activate a virtual environment if needed, then install dependencies:

```bash
python -m pip install -r requirements.txt
```

### 11.2 Running the numerical calculation

From the project root, run:

```bash
python -m scripts.run_calculation
```

### 11.3 Running convergence studies

Run the grid-convergence study:

```bash
python -m scripts.convergence_study
```

The domain- and wavefunction-convergence studies can be run similarly if their scripts support module execution.

### 11.4 Running automated tests

Execute the test suite:

```bash
python -m pytest tests/ -v
```

The recorded test run completed with **61 passing tests**.

### 11.5 Generated results

The `results/` directory stores generated numerical results, including:

- `energy_comparison.csv`
- `validation_metrics.csv`
- `wavefunctions.png`
- `probability_densities.png`

Regenerate these outputs whenever the numerical implementation or calculation parameters change.

## 12. References

1. Griffiths, D. J., and Schroeter, D. F. *Introduction to Quantum Mechanics*. 3rd edition. Cambridge University Press, 2018.
2. Levine, I. N. *Quantum Chemistry*. 7th edition. Pearson, 2013.
3. Atkins, P., and Friedman, R. *Molecular Quantum Mechanics*. 4th edition. Oxford University Press, 2005.
4. NumPy Developers. *NumPy Documentation*. https://numpy.org/doc/
5. SciPy Developers. *SciPy Documentation*. https://docs.scipy.org/doc/scipy/
6. pytest Development Team. *pytest Documentation*. https://docs.pytest.org/
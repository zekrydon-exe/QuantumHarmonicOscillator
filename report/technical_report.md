# Numerical Study of the Quantum Harmonic Oscillator

## 1. Abstract

This project investigates the quantum harmonic oscillator by numerically
solving the time-independent Schrödinger equation using the finite-difference
method. The dimensionless Hamiltonian is discretised on a uniformly spaced
one-dimensional grid, producing a real symmetric matrix whose eigenvalues
and eigenvectors approximate the energy levels and wavefunctions of the
quantum system.

The numerical energies are compared with the analytical energy spectrum.
The computed wavefunctions are examined through normalisation,
orthogonality, parity and probability-density properties. Numerical
reliability is further assessed using eigenpair residuals, automated
tests and grid-convergence studies.

For the standard calculation, the computational domain is [-8, 8], with
399 interior grid points and a grid spacing of 0.04. The five lowest
dimensionless energy levels are calculated. The numerical results agree
closely with the analytical spectrum, while the remaining discrepancies
are primarily associated with the finite-difference discretisation.

## 2. Introduction

The quantum harmonic oscillator is a fundamental model in quantum
mechanics. It describes a particle subject to a potential energy that
increases quadratically with displacement from an equilibrium position.

Despite the simplicity of its potential, the system provides a useful
framework for studying quantisation, stationary states, wavefunctions,
probability densities and numerical solutions of the Schrödinger equation.

The harmonic oscillator is particularly suitable for computational
investigation because its energy eigenvalues and eigenfunctions are known
analytically. These exact solutions provide reference values against which
numerical methods can be evaluated.

### 2.1 Objectives

The principal objectives of this project are:

1. Formulate the one-dimensional quantum harmonic oscillator in
   dimensionless form.
2. Discretise the Schrödinger equation using a finite-difference
   approximation.
3. Construct the corresponding Hamiltonian matrix.
4. Calculate the lowest energy eigenvalues and eigenvectors numerically.
5. Compare the numerical energy spectrum with the analytical solution.
6. Examine the normalisation, orthogonality and parity of the
   numerical wavefunctions.
7. Investigate numerical accuracy through residual calculations and
   grid-convergence studies.
8. Develop a reproducible Python implementation supported by automated
   tests.

## 3. Theoretical Background

### 3.1 The harmonic oscillator potential

For a particle of mass m moving in a harmonic potential, the potential
energy is

$$
V(x) = \frac{1}{2}m\omega^2x^2,
$$

where omega is the angular frequency of the oscillator.

The time-independent Schrödinger equation is

$$
-\frac{\hbar^2}{2m}
\frac{d^2\psi(x)}{dx^2}
+
\frac{1}{2}m\omega^2x^2\psi(x)
=
E\psi(x).
$$

Here, psi(x) is the stationary-state wavefunction and E is the
corresponding energy eigenvalue.

### 3.2 Dimensionless formulation

Introduce the dimensionless coordinate

$$
q = \sqrt{\frac{m\omega}{\hbar}}\,x
$$

and the dimensionless energy

$$
\varepsilon = \frac{E}{\hbar\omega}.
$$

The Schrödinger equation becomes

$$
-\frac{1}{2}\frac{d^2\psi(q)}{dq^2}
+
\frac{1}{2}q^2\psi(q)
=
\varepsilon\psi(q).
$$

Equivalently, the dimensionless Hamiltonian is

$$
\hat{H}
=
-\frac{1}{2}\frac{d^2}{dq^2}
+
\frac{1}{2}q^2.
$$

The corresponding eigenvalue equation is

$$
\hat{H}\psi_n(q)
=
\varepsilon_n\psi_n(q).
$$

This dimensionless representation removes the explicit dependence on
the particle mass and oscillator frequency from the numerical problem.

### 3.3 Analytical energy spectrum

The exact dimensionless energy eigenvalues are

$$
\varepsilon_n = n + \frac{1}{2},
\qquad n=0,1,2,\ldots
$$

The corresponding dimensional energies are

$$
E_n = \hbar\omega\left(n+\frac{1}{2}\right).
$$

The ground-state energy is therefore

$$
\varepsilon_0 = \frac{1}{2}.
$$

The energy difference between consecutive states is constant:

$$
\varepsilon_{n+1}-\varepsilon_n=1.
$$

This equally spaced spectrum provides an important reference for
validating the numerical eigenvalues.

### 3.4 Analytical wavefunctions

The normalised analytical eigenfunctions are

$$
\psi_n(q)
=
\frac{H_n(q)e^{-q^2/2}}
{\pi^{1/4}\sqrt{2^n n!}},
$$

where H_n(q) denotes the physicists' Hermite polynomial of degree n.

For the ground state, H_0(q) = 1, giving

$$
\psi_0(q)
=
\frac{1}{\pi^{1/4}}e^{-q^2/2}.
$$

These analytical wavefunctions provide reference solutions for examining
the numerical eigenstates.

### 3.5 Normalisation and orthogonality

A physically admissible bound-state wavefunction must satisfy the
normalisation condition

$$
\int_{-\infty}^{\infty}
|\psi_n(q)|^2\,dq=1.
$$

Different energy eigenstates of this Hamiltonian are orthogonal:

$$
\int_{-\infty}^{\infty}
\psi_m^*(q)\psi_n(q)\,dq
=
\delta_{mn},
$$

where delta_mn is the Kronecker delta.

For a numerical calculation on a finite grid, these integrals are
approximated over the computational domain. The accuracy of the
approximation depends on the grid spacing and the chosen domain.

## 4. Numerical Methodology

### 4.1 Computational grid

The continuous coordinate q is replaced by a uniformly spaced grid
containing N interior points between the boundaries q_min and q_max.

The grid spacing is defined as

$$
\Delta q
=
\frac{q_{\max}-q_{\min}}{N+1}.
$$

The interior grid coordinates are

$$
q_i=q_{\min}+i\Delta q,
\qquad i=1,2,\ldots,N.
$$

The boundary values are excluded from the array of interior points.
This corresponds to imposing zero wavefunction values at the
computational boundaries.

For the standard calculation,

$$
q_{\min}=-8,\qquad q_{\max}=8,
\qquad N=399.
$$

Consequently,

$$
\Delta q=\frac{16}{400}=0.04.
$$

### 4.2 Finite-difference approximation

The second derivative is approximated using the central
finite-difference formula:

$$
\frac{d^2\psi}{dq^2}\bigg|_{q_i}
\approx
\frac{
\psi_{i+1}-2\psi_i+\psi_{i-1}
}{(\Delta q)^2}.
$$

Substituting this expression into the dimensionless Schrödinger
equation gives

$$
-\frac{1}{2}
\frac{
\psi_{i+1}-2\psi_i+\psi_{i-1}
}{(\Delta q)^2}
+
\frac{q_i^2}{2}\psi_i
=
\varepsilon\psi_i.
$$

Rearranging,

$$
-\frac{1}{2(\Delta q)^2}\psi_{i-1}
+
\left[
\frac{1}{(\Delta q)^2}
+
\frac{q_i^2}{2}
\right]\psi_i
-
\frac{1}{2(\Delta q)^2}\psi_{i+1}
=
\varepsilon\psi_i.
$$

This expression defines the matrix elements of the discretised
Hamiltonian.

### 4.3 Kinetic-energy matrix

The kinetic-energy operator is

$$
\hat{T}
=
-\frac{1}{2}\frac{d^2}{dq^2}.
$$

Its finite-difference matrix has the diagonal elements

$$
T_{ii}=\frac{1}{(\Delta q)^2}
$$

and the nearest-neighbour off-diagonal elements

$$
T_{i,i+1}
=
T_{i+1,i}
=
-\frac{1}{2(\Delta q)^2}.
$$

All other elements are zero.

Therefore, the kinetic-energy matrix is tridiagonal and symmetric.

### 4.4 Potential-energy matrix

The dimensionless potential energy is

$$
V(q)=\frac{q^2}{2}.
$$

In the position-grid representation, the potential-energy matrix is
diagonal:

$$
V_{ij}
=
\begin{cases}
\dfrac{q_i^2}{2}, & i=j,\
0, & i\ne j.
\end{cases}
$$

### 4.5 Hamiltonian matrix

The complete discretised Hamiltonian is constructed as

$$
H=T+V.
$$

Its diagonal elements are

$$
H_{ii}
=
\frac{1}{(\Delta q)^2}
+
\frac{q_i^2}{2},
$$

while the nearest-neighbour elements are

$$
H_{i,i+1}
=
H_{i+1,i}
=
-\frac{1}{2(\Delta q)^2}.
$$

The resulting matrix is real, symmetric and tridiagonal.

### 4.6 Numerical eigenvalue problem

The discretised Schrödinger equation becomes the matrix eigenvalue
problem

$$
H\mathbf{v}_n
=
\varepsilon_n\mathbf{v}_n,
$$

where the eigenvalues approximate the dimensionless energies and
the eigenvectors represent the corresponding discrete states.

The implementation uses the symmetric eigensolver provided by
`scipy.linalg.eigh`. Only the requested lowest-energy eigenstates
are calculated.

The discrete eigenvectors are normalised according to the
Euclidean norm. To obtain wavefunctions suitable for integration
on a uniform coordinate grid, the implementation scales each
eigenvector by

$$
\psi_n(q_i)
=
\frac{v_{i,n}}{\sqrt{\Delta q}}.
$$

This scaling makes the discrete normalisation condition

$$
\sum_i |v_{i,n}|^2=1
$$

consistent with the uniform-grid approximation

$$
\sum_i |\psi_n(q_i)|^2\Delta q\approx 1.
$$

## 5. Software Implementation

The project is implemented in Python using NumPy, SciPy and
Matplotlib. The code is organised into modules for constructing
the Hamiltonian, calculating analytical solutions, solving the
numerical eigenvalue problem, performing convergence studies and
validating the results.

### 5.1 Main source modules

#### `src/hamiltonian.py`

This module implements the numerical grid and the discretised
Hamiltonian components.

Its principal functions are:

- `create_grid()`: constructs a uniformly spaced grid of interior
  points and returns the grid spacing.
- `kinetic_energy()`: constructs the finite-difference
  kinetic-energy matrix.
- `potential_energy()`: constructs the diagonal potential-energy
  matrix for the harmonic oscillator.
- `build_hamiltonian()`: combines the kinetic and potential
  matrices to obtain the Hamiltonian.

Input validation is included in the grid-generation function to
reject invalid boundaries and grid-point counts.

#### `src/analytical_solution.py`

This module provides the analytical reference solutions.

- `analytical_energy(n)` calculates the exact dimensionless
  energy for quantum number n.
- `analytical_wavefunction(q, n)` evaluates the analytical
  harmonic-oscillator wavefunction using Hermite polynomials.

These functions allow the numerical results to be compared
against known solutions.

#### `src/harmonic_oscillator.py`

The `solve_harmonic_oscillator()` function performs the main
numerical calculation.

The function constructs the grid, builds the Hamiltonian and
uses `scipy.linalg.eigh` to calculate the requested lowest-energy
eigenstates.

The returned wavefunctions are scaled by the square root of the
grid spacing to make them suitable for continuum-style
integration on the numerical grid.

### 5.2 Supporting scripts

The project contains scripts for running calculations, plotting
wavefunctions, validating numerical results and investigating
convergence.

The convergence-related scripts include:

- `scripts/convergence_study.py`: investigates the dependence
  of the calculated energy levels on the number of grid points.
- `scripts/domain_convergence.py`: investigates the effect of
  changing the computational domain.
- `scripts/wavefunction_convergence.py`: investigates numerical
  wavefunction convergence.
- `scripts/validate_wavefunctions.py`: supports wavefunction
  validation.
- `scripts/plot_wavefunctions.py`: generates wavefunction plots.
- `scripts/solve_eigenvalues.py`: calculates numerical energy
  eigenvalues.
- `scripts/run_calculation.py`: provides a script for running
  the calculation workflow.

The generated results include energy-comparison and validation
CSV files, together with plots of wavefunctions and probability
densities.

### 5.3 Automated testing

The project uses pytest to automate checks of the numerical
implementation.

The tests cover grid construction, input validation, matrix
properties, analytical solutions, numerical energy levels,
wavefunction properties and physical consistency.

This testing framework helps identify errors in the numerical
implementation and makes subsequent modifications easier to
verify.

## 6. Results and Discussion

### 6.1 Energy eigenvalues

The five lowest dimensionless energy levels were calculated
using the standard configuration:

- Computational domain: [-8, 8]
- Number of interior grid points: 399
- Grid spacing: 0.04
- Number of states: 5

The exact energy spectrum is

$$
\varepsilon_n=n+\frac{1}{2}.
$$

The numerical and analytical energies are compared in
Table 1.

**Table 1. Numerical and analytical energy eigenvalues.**

| Quantum number n | Numerical energy | Analytical energy | Absolute error |
|---:|---:|---:|---:|
| 0 | 0.499949994999 | 0.5 | 5.0005 × 10⁻⁵ |
| 1 | 1.499749954978 | 1.5 | 2.5005 × 10⁻⁴ |
| 2 | 2.499349824875 | 2.5 | 6.5018 × 10⁻⁴ |
| 3 | 3.498749544572 | 3.5 | 1.2505 × 10⁻³ |
| 4 | 4.497949053887 | 4.5 | 2.0509 × 10⁻³ |

The numerical energies are slightly lower than the corresponding
analytical values.

The absolute error increases with the quantum number over the
five states considered. This is consistent with the increasing
sensitivity of higher-energy states to the finite-difference
approximation.

The errors remain small relative to the corresponding energy
values. However, their magnitude depends on the grid spacing,
computational domain and state being considered.

### 6.2 Energy-level spacing

The analytical harmonic-oscillator spectrum has equally spaced
energy levels:

$$
\varepsilon_{n+1}-\varepsilon_n=1.
$$

The numerical energy levels reproduce this behaviour
approximately. The deviations from exact equal spacing arise
from the numerical approximation of the Hamiltonian.

The energy-spacing property is tested automatically by the
project's physics-validation tests.

### 6.3 Wavefunctions and probability densities

The numerical eigenvectors provide approximations to the
stationary-state wavefunctions.

The probability density associated with the nth state is

$$
\rho_n(q)=|\psi_n(q)|^2.
$$

The wavefunction and probability-density plots generated by the
project allow the spatial distributions of the states to be
examined.

The expected qualitative properties are:

- The ground-state wavefunction has even parity.
- Successive eigenstates alternate between even and odd parity.
- The probability density is non-negative.
- The probability density of each stationary state is symmetric
  about the origin.
- Higher-energy states exhibit additional nodes.

These properties provide useful checks on the numerical
eigenstates and their physical interpretation.

### 6.4 Interpretation of the energy errors

The finite-difference approximation replaces a continuous
differential operator with a discrete matrix. Consequently,
the numerical eigenvalues are approximations to the exact
eigenvalues.

For sufficiently smooth solutions, the central finite-difference
approximation to the second derivative has a leading truncation
error of order (Δq)².

For a fixed computational domain and sufficiently small grid
spacing, this leads to second-order convergence of the numerical
eigenvalues under the usual conditions.

The finite computational domain introduces a separate source
of error because the wavefunction is constrained to vanish at
the boundaries. For the standard domain [-8, 8], the low-lying
harmonic-oscillator wavefunctions are strongly localised away
from the boundaries, making this effect small for the states
examined.

The relative contributions of discretisation and boundary errors
are investigated separately in the convergence studies.

## 7. Numerical Validation

Several independent checks are used to assess the numerical
solution.

### 7.1 Energy comparison

The numerical energy eigenvalues are compared against the exact
spectrum:

$$
\varepsilon_n=n+\frac{1}{2}.
$$

The absolute error is defined as

$$
\Delta\varepsilon_n
=
\left|
\varepsilon_n^{\mathrm{numerical}}
-
\varepsilon_n^{\mathrm{exact}}
\right|.
$$

The values in Table 1 quantify the accuracy of the first five
numerical energy levels.

### 7.2 Wavefunction normalisation

The continuum normalisation condition is

$$
\int_{-\infty}^{\infty}|\psi_n(q)|^2\,dq=1.
$$

The numerical approximation is evaluated using trapezoidal
integration over the computational grid:

$$
I_n
\approx
\int_{q_{\min}}^{q_{\max}}|\psi_n(q)|^2\,dq.
$$

The computed normalisation integrals are approximately equal
to one for the states tested.

### 7.3 Orthogonality

Different eigenstates of the real symmetric Hamiltonian are
orthogonal. Their numerical overlap is evaluated as

$$
S_{mn}
\approx
\int_{q_{\min}}^{q_{\max}}
\psi_m(q)\psi_n(q)\,dq.
$$

For distinct states, the expected overlap is zero.

The automated tests verify that the overlaps between the first
five numerical wavefunctions remain below the specified
tolerance of 10⁻³.

### 7.4 Eigenpair residuals

A numerical eigenpair should satisfy

$$
H\mathbf{v}_n
\approx
\varepsilon_n\mathbf{v}_n.
$$

The residual vector is defined as

$$
\mathbf{r}_n
=
H\mathbf{v}_n-\varepsilon_n\mathbf{v}_n.
$$

Its Euclidean norm is

$$
r_n=\|\mathbf{r}_n\|_2.
$$

For the standard calculation, the maximum residual norm
reported for the five calculated states is approximately

$$
\max_n r_n
\approx 1.81\times10^{-13}.
$$

This small residual indicates that the numerical eigenpairs
satisfy the discretised matrix eigenvalue problem to high
numerical precision.

A small residual does not, by itself, establish that the
discretised solution is close to the exact continuous solution.
The energy comparisons and convergence studies are therefore
also necessary.

### 7.5 Automated test results

The project's automated test suite checks:

- Grid dimensions, spacing, symmetry and boundary exclusion.
- Rejection of invalid grid parameters.
- Hamiltonian symmetry and tridiagonal structure.
- Correct kinetic- and potential-energy matrix elements.
- Analytical energy and wavefunction properties.
- Numerical energy accuracy.
- Wavefunction normalisation and orthogonality.
- Eigenpair residuals.
- Energy spacing and wavefunction parity.
- Probability-density properties.
- Rejection of invalid solver parameters.

The verified test run completed with **47 passing tests**.

This result confirms that all tests collected in that run passed.
It does not imply that every possible input or numerical
configuration has been tested.

## 8. Convergence Analysis

Convergence analysis examines how numerical results change as
the grid is refined or the computational domain is enlarged.

Both investigations are important because grid discretisation
and finite-domain boundaries introduce distinct numerical errors.

### 8.1 Grid-spacing convergence

The grid-convergence study uses the fixed computational domain
[-4, 4] and the following numbers of interior grid points:

$$
N\in\{49,99,199,399\}.
$$

For each grid, the numerical eigenvalues are compared with the
analytical values.

The absolute error for state n is

$$
e_n(\Delta q)
=
\left|
\varepsilon_n^{\mathrm{numerical}}
-
\varepsilon_n^{\mathrm{exact}}
\right|.
$$

For two calculations with spacings Δq₁ and Δq₂, the observed
convergence order is estimated by

$$
p
=
\frac{
\ln\left(e_n(\Delta q_1)/e_n(\Delta q_2)\right)
}{
\ln(\Delta q_1/\Delta q_2)
}.
$$

If the errors are dominated by a leading term proportional to
(Δq)², the observed order should approach two as the grid is
refined.

The estimate can be unreliable when the errors are very small
or when different error sources dominate at different
resolutions.

The domain [-4, 4] is deliberately retained in this study.
Although it is useful for investigating grid refinement, the
higher-energy states may be more sensitive to the boundary
conditions than they are in the larger standard domain.

The numerical convergence orders should therefore be taken
from the output of the convergence script rather than assumed
to be exactly two.

### 8.2 Computational-domain convergence

A separate domain-convergence study investigates the effect
of changing the computational boundaries while maintaining
an appropriate grid resolution.

Increasing the domain moves the imposed zero-wavefunction
boundaries farther into regions where bound-state wavefunctions
are small.

For the harmonic oscillator, this reduces the influence of
the artificial boundaries on the low-lying states.

The domain-convergence script is used to compare the resulting
energy eigenvalues and assess whether the chosen computational
domain is sufficiently large.

The standard calculation uses [-8, 8]. Previous checks at a
grid spacing of approximately 0.04 found the low-lying energies
for domains [-8, 8] and [-10, 10] to be effectively identical
at the reported precision.

This supports the use of [-8, 8] for the standard calculation.
The conclusion applies to the states and numerical precision
examined; it should not automatically be extended to arbitrarily
high-energy states.

### 8.3 Wavefunction convergence

Energy convergence alone does not completely establish the
accuracy of the calculated eigenstates.

The wavefunction-convergence study examines how numerical
wavefunctions change as the grid is refined. Meaningful
comparisons must account for the fact that wavefunctions are
defined on different grids and that an eigenvector can be
multiplied by an overall sign without changing its physical
state.

The analytical wavefunctions provide an additional reference
for assessing the spatial structure of the numerical solutions.

Together, the energy, domain and wavefunction convergence
studies provide complementary assessments of numerical
accuracy.

## 9. Limitations and Possible Improvements

Although the numerical results agree closely with the analytical
solution, several limitations should be considered.

### 9.1 Finite computational domain

The numerical calculation is performed over a finite interval
rather than the entire real line.

The wavefunction is constrained to vanish at the computational
boundaries. This approximates the behaviour of the bound-state
wavefunctions at large distances but introduces a potential
source of error.

Increasing the domain can reduce this error, although it also
increases the computational cost if the grid spacing is kept
constant.

### 9.2 Finite-difference discretisation

The second derivative is approximated using a central
finite-difference formula. This introduces discretisation error
that depends on the grid spacing.

A finer grid generally improves the approximation, provided
that other numerical errors do not dominate.

Higher-order finite-difference schemes could be investigated
to determine whether comparable accuracy can be achieved with
fewer grid points.

### 9.3 Wavefunction comparison

Numerical eigenvectors are defined only up to an overall sign
for a real-valued, non-degenerate eigenstate. Consequently,
direct comparisons between wavefunctions should account for
possible sign differences.

Furthermore, wavefunctions calculated on different grids must
be interpolated onto a common grid or compared using another
appropriate numerical measure.

### 9.4 Numerical integration

Normalisation and orthogonality are evaluated using numerical
integration over the finite computational domain.

The resulting integrals are approximations to their continuous
counterparts. Their accuracy depends on the grid resolution,
integration method and domain size.

### 9.5 Computational efficiency

The current implementation constructs a dense Hamiltonian
matrix even though the finite-difference Hamiltonian is
tridiagonal.

For the grid sizes used in this project, the dense
implementation is convenient and straightforward to validate.
For substantially larger grids, a sparse-matrix representation
and a suitable sparse eigensolver could reduce memory
requirements and computational cost.

### 9.6 Further investigations

Possible extensions include:

- Comparing dense and sparse eigensolvers.
- Investigating higher-order finite-difference schemes.
- Quantifying wavefunction errors against analytical solutions.
- Extending the calculation to other one-dimensional
  potentials.
- Investigating how the numerical error changes with the
  quantum number.
- Developing a more systematic error analysis that separates
  discretisation and boundary effects.

## 10. Conclusion

This project demonstrates a numerical solution of the
one-dimensional quantum harmonic oscillator using the
finite-difference method.

The dimensionless Schrödinger equation was discretised on a
uniform grid, and the resulting real symmetric Hamiltonian
was solved using a numerical eigensolver.

The five lowest numerical energy eigenvalues were compared
with the analytical spectrum. For the standard calculation
on the domain [-8, 8], using 399 interior grid points, the
numerical energies were close to the exact values.

The calculated wavefunctions were examined through
normalisation, orthogonality, parity and probability-density
properties. The eigenpair residuals were also evaluated,
with the maximum residual norm approximately equal to
1.81 × 10⁻¹³ for the five states considered.

Automated tests were used to check the mathematical
properties of the numerical implementation and to detect
invalid input parameters.

Grid-spacing, computational-domain and wavefunction
convergence studies provide additional methods for assessing
the reliability of the numerical results.

Overall, the project illustrates how a continuous quantum
mechanical eigenvalue problem can be transformed into a
finite-dimensional matrix problem. It also demonstrates why
numerical validation requires more than a successful
eigensolver: analytical comparisons, physical consistency
checks and convergence studies are necessary to assess
the accuracy of the calculated results.

## 11. Reproducibility

The project is implemented in Python and uses NumPy, SciPy,
Matplotlib and pytest.

### 11.1 Environment setup

Create and activate a virtual environment if one is not
already configured.

Install the dependencies listed in `requirements.txt`:

```bash
python -m pip install -r requirements.txt

```

### 11.2 Running the numerical calculation

From the project root, run the calculation script:

```bash
python -m scripts.run_calculation
```

This assumes that `scripts/run_calculation.py` supports execution as a module.

### 11.3 Running the convergence study

Run the grid-convergence study:

```bash
python -m scripts.convergence_study
```

The domain-convergence and wavefunction-convergence studies can be
run similarly, provided their scripts support module execution.

### 11.4 Running the automated tests

Execute the complete test suite:

```bash
python -m pytest -v
```

The verified test run completed with 47 passing tests.

### 11.5 Generated results

The project stores generated numerical results in the `results/`
directory. The principal outputs include:

- `energy_comparison.csv`
- `validation_metrics.csv`
- `wavefunctions.png`
- `probability_densities.png`

Regenerate these outputs whenever the numerical implementation
or calculation parameters change.

## 12. References

1. Griffiths, D. J., and Schroeter, D. F. *Introduction to Quantum
   Mechanics*. 3rd edition. Cambridge University Press, 2018.

2. Levine, I. N. *Quantum Chemistry*. 7th edition.
   Pearson, 2013.

3. Atkins, P., and Friedman, R. *Molecular Quantum Mechanics*.
   4th edition. Oxford University Press, 2005.

4. NumPy Developers. *NumPy Documentation*.
   https://numpy.org/doc/

5. SciPy Developers. *SciPy Documentation*.
   https://docs.scipy.org/doc/scipy/

6. pytest Development Team. *pytest Documentation*.
   https://docs.pytest.org/

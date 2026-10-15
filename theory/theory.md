# Quantum Harmonic Oscillator

## 1. Theoretical Model

### 1.1 Introduction

The quantum harmonic oscillator is a fundamental model in quantum mechanics. It describes a particle moving in a quadratic potential and provides an exactly solvable system against which numerical methods can be tested.

The model is relevant to molecular vibrations, quantised vibrational energy levels, and the study of more complicated quantum systems.

This project solves the time-independent Schrödinger equation numerically and compares the calculated energies and wavefunctions with their analytical counterparts.

### 1.2 Potential Energy

The potential energy of a one-dimensional harmonic oscillator is

$$
V(x)=\frac{1}{2}kx^2,
$$

where $x$ is the displacement from equilibrium and $k$ is the force constant. The potential is parabolic, with a minimum at $x=0$. The corresponding restoring force is linear in displacement.

### 1.3 Time-Independent Schrödinger Equation

For a particle of mass $m$ moving in one dimension, the time-independent Schrödinger equation is

$$
\hat{H}\psi(x)=E\psi(x),
$$

where $\hat{H}$ is the Hamiltonian operator, $\psi(x)$ is the stationary-state wavefunction, and $E$ is the energy eigenvalue.

The Hamiltonian consists of kinetic and potential energy operators:

$$
\hat{H}
=
-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}
+\frac{1}{2}kx^2.
$$

The Schrödinger equation is therefore

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}
+\frac{1}{2}kx^2\psi(x)
=E\psi(x).
$$

The first term represents kinetic energy and the second represents potential energy.

### 1.4 Dimensionless Formulation

Define the angular frequency

$$
\omega=\sqrt{\frac{k}{m}},
$$

and the dimensionless coordinate

$$
q=\sqrt{\frac{m\omega}{\hbar}}\,x.
$$

The dimensionless energy is

$$
\varepsilon=\frac{E}{\hbar\omega}.
$$

Using these definitions, the Schrödinger equation becomes

$$
-\frac{1}{2}\frac{d^2\psi(q)}{dq^2}
+\frac{1}{2}q^2\psi(q)
=\varepsilon\psi(q).
$$

The dimensionless Hamiltonian is

$$
\hat{H}
=
-\frac{1}{2}\frac{d^2}{dq^2}
+\frac{1}{2}q^2,
$$

and the eigenvalue equation is

$$
\hat{H}\psi_n(q)=\varepsilon_n\psi_n(q).
$$

This formulation removes explicit dependence on the particle mass, force constant, and reduced Planck constant from the dimensionless eigenvalue problem.

### 1.5 Analytical Energy Eigenvalues

The exact dimensionless energy eigenvalues are

$$
\varepsilon_n=n+\frac{1}{2},
\qquad n=0,1,2,\ldots
$$

The corresponding physical energies are

$$
E_n=\hbar\omega\left(n+\frac{1}{2}\right).
$$

The ground-state energy is non-zero:

$$
\varepsilon_0=\frac{1}{2}.
$$

This is the zero-point energy. Adjacent levels are equally spaced:

$$
\varepsilon_{n+1}-\varepsilon_n=1,
$$

or, in physical units, by $\hbar\omega$.

### 1.6 Analytical Wavefunctions

The normalised analytical wavefunctions are

$$
\psi_n(q)
=
\frac{H_n(q)e^{-q^2/2}}
{\pi^{1/4}\sqrt{2^n n!}},
$$

where $H_n(q)$ is the physicists' Hermite polynomial of degree $n$.

The first five Hermite polynomials are

$$
H_0(q)=1,
$$

$$
H_1(q)=2q,
$$

$$
H_2(q)=4q^2-2,
$$

$$
H_3(q)=8q^3-12q,
$$

$$
H_4(q)=16q^4-48q^2+12.
$$

States with even $n$ have even parity, and states with odd $n$ have odd parity. The $n$-th oscillator eigenfunction has $n$ nodes.

### 1.7 Normalisation and Orthogonality

A normalised bound-state wavefunction satisfies

$$
\int_{-\infty}^{\infty}|\psi_n(q)|^2\,dq=1.
$$

Different energy eigenstates are orthogonal:

$$
\int_{-\infty}^{\infty}
\psi_m^*(q)\psi_n(q)\,dq
=\delta_{mn},
$$

where $\delta_{mn}$ is the Kronecker delta.

These properties provide criteria for validating the numerical wavefunctions.

## 2. Numerical Methodology

### 2.1 Computational Domain and Grid Discretisation

The analytical oscillator is defined over the infinite domain

$$
-\infty<q<\infty.
$$

Numerical calculations use a finite interval

$$
q_{\min}\leq q\leq q_{\max}.
$$

For $N$ interior points on a uniform grid, with the endpoints treated as boundaries, the spacing is

$$
h=\frac{q_{\max}-q_{\min}}{N+1}.
$$

The interior grid points are

$$
q_i=q_{\min}+ih,
\qquad i=1,2,\ldots,N.
$$

The wavefunction is represented by its values at the grid points:

$$
\psi(q_i)\approx\psi_i.
$$

The finite-domain calculation effectively imposes zero wavefunction values at the endpoints. This is a good approximation when the bound-state wavefunctions are negligible there. The adequacy of the chosen domain is examined separately through domain convergence.

### 2.2 Central Finite-Difference Approximation

The dimensionless Schrödinger equation is

$$
-\frac{1}{2}\frac{d^2\psi(q)}{dq^2}
+\frac{1}{2}q^2\psi(q)
=\varepsilon\psi(q).
$$

At grid point $q_i$, the central finite-difference approximation is

$$
\left.\frac{d^2\psi}{dq^2}\right|_{q_i}
\approx
\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{h^2}.
$$

Substituting this approximation into the Schrödinger equation gives

$$
-\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{2h^2}
+\frac{1}{2}q_i^2\psi_i
=\varepsilon\psi_i.
$$

Rearranging,

$$
-\frac{1}{2h^2}\psi_{i-1}
+
\left(\frac{1}{h^2}+\frac{1}{2}q_i^2\right)\psi_i
-\frac{1}{2h^2}\psi_{i+1}
=\varepsilon\psi_i.
$$

The central finite-difference approximation has a leading truncation error of order $h^2$ for sufficiently smooth functions.

### 2.3 Construction of the Hamiltonian Matrix

The discretised equation is written as a matrix eigenvalue problem:

$$
\mathbf{H}\boldsymbol{\psi}
=\varepsilon\boldsymbol{\psi}.
$$

The diagonal elements are

$$
H_{ii}=\frac{1}{h^2}+\frac{1}{2}q_i^2,
$$

and the nearest-neighbour off-diagonal elements are

$$
H_{i,i+1}=H_{i+1,i}=-\frac{1}{2h^2}.
$$

All other elements are zero. The Hamiltonian is therefore tridiagonal. Its diagonal contains the finite-difference kinetic-energy contribution and the potential energy, while its off-diagonal elements arise from the kinetic-energy operator.

The matrix is real and symmetric, so its eigenvalues are real and its eigenvectors can be chosen orthogonal under the discrete Euclidean inner product.

### 2.4 Numerical Eigenvalue Solution

The numerical eigenpairs satisfy

$$
\mathbf{H}\boldsymbol{\psi}_n
=\varepsilon_n\boldsymbol{\psi}_n.
$$

The eigenvalues are ordered from lowest to highest energy. The first eigenvalue corresponds to the ground state, and subsequent eigenvalues correspond to excited states.

The implementation uses `scipy.linalg.eigh` to solve the real symmetric matrix eigenvalue problem and requests only the lowest five states for the production calculation.

The numerical energies are compared with

$$
\varepsilon_n^{\mathrm{exact}}=n+\frac{1}{2}.
$$

### 2.5 Numerical Wavefunction Scaling

The eigensolver returns eigenvectors normalised according to the discrete Euclidean norm:

$$
\sum_{i=1}^{N}|v_i|^2=1.
$$

For a uniform grid, the corresponding continuum integral is approximated by a weighted sum. The eigenvectors are therefore scaled as

$$
\psi_i=\frac{v_i}{\sqrt{h}}.
$$

This gives wavefunctions with the appropriate continuum scaling. Their normalisation is then checked using numerical integration over the computational grid.

## 3. Numerical Validation and Convergence

### 3.1 Energy Eigenvalue Validation

The absolute energy error is

$$
\Delta\varepsilon_n
=
\left|
\varepsilon_n^{\mathrm{numerical}}
-\varepsilon_n^{\mathrm{exact}}
\right|.
$$

The relative energy error is

$$
\delta_n
=
\frac{
\left|
\varepsilon_n^{\mathrm{numerical}}
-\varepsilon_n^{\mathrm{exact}}
\right|
}{
\left|\varepsilon_n^{\mathrm{exact}}\right|
}.
$$

The percentage relative error is

$$
\delta_n(\%)=100\delta_n.
$$

Absolute error measures the magnitude of the discrepancy, while relative error expresses it relative to the exact energy.

### 3.2 Wavefunction Normalisation

The normalisation integral is

$$
I_n=\int_{q_{\min}}^{q_{\max}}|\psi_n(q)|^2\,dq.
$$

It is evaluated numerically using the trapezoidal rule. A value close to unity indicates that the wavefunction is normalised over the selected domain. Normalisation alone does not establish that the wavefunction has the correct shape, so comparison with analytical wavefunctions is also performed.

### 3.3 Wavefunction Comparison

An eigenfunction can be multiplied by an overall sign without changing its physical state. The numerical wavefunction is therefore sign-aligned with the analytical wavefunction before calculating the difference.

The absolute overlap is

$$
S_n=
\left|
\int
\psi_n^{\mathrm{analytical}}(q)
\psi_n^{\mathrm{numerical}}(q)\,dq
\right|.
$$

For matching normalised wavefunctions, the overlap is unity.

The relative $L^2$ error is

$$
e_n=
\frac{
\left\|
\psi_n^{\mathrm{numerical}}
-\psi_n^{\mathrm{analytical}}
\right\|_2
}{
\left\|\psi_n^{\mathrm{analytical}}\right\|_2
}.
$$

A smaller value indicates closer agreement.

### 3.4 Eigenpair Residuals

The residual vector is

$$
\mathbf{r}_n
=
\mathbf{H}\boldsymbol{\psi}_n
-\varepsilon_n\boldsymbol{\psi}_n.
$$

Its Euclidean norm is

$$
R_n=\|\mathbf{r}_n\|_2.
$$

A small residual indicates that the computed eigenpair satisfies the discretised matrix equation to high numerical accuracy. It does not necessarily mean that the eigenvalue is equally close to the exact eigenvalue of the continuous differential equation. The residual assesses the matrix problem; the analytical energy error assesses the approximation to the continuous problem.

### 3.5 Grid Convergence

The grid convergence study keeps the domain fixed at $[-4,4]$ and refines the grid. For a sufficiently smooth wavefunction and an adequately large fixed domain, the leading finite-difference error is expected to scale as

$$
\Delta\varepsilon_n\propto h^2.
$$

Halving the grid spacing should therefore reduce the leading discretisation error by approximately a factor of four when other error sources are small.

An observed convergence order can be estimated from two calculations:

$$
p\approx
\frac{
\ln\left(|\Delta\varepsilon_{n,1}|/
|\Delta\varepsilon_{n,2}|\right)
}{
\ln(h_1/h_2)
}.
$$

This estimate requires non-zero errors and is most meaningful when the same error regime dominates both calculations.

The ground-state results show observed orders close to two as the grid is refined, consistent with the central finite-difference method. Higher-state energy errors can be more affected by the finite domain, so their convergence should be interpreted with that limitation in mind.

### 3.6 Domain Convergence

A finite computational domain introduces an approximation because the exact problem extends over all space. If the boundaries are too close to the centre, the wavefunction may not have decayed sufficiently before reaching them.

The domain study compares

$$
[-4,4],\quad[-6,6],\quad[-8,8],\quad[-10,10],
$$

while maintaining a target grid spacing of approximately $0.04$.

The first five energies become effectively unchanged between $[-8,8]$ and $[-10,10]$. This supports using $[-8,8]$ for the production calculation at the selected resolution. It does not eliminate finite-difference discretisation error.

### 3.7 Sources of Numerical Error

The main sources of numerical error are:

1. **Finite-difference discretisation:** the second derivative is approximated on a discrete grid.
2. **Finite computational domain:** the infinite domain is replaced by a bounded interval.
3. **Numerical integration:** continuous integrals are approximated numerically.
4. **Floating-point arithmetic:** computations use finite-precision numbers.

Grid refinement and domain enlargement address different sources of error and should be investigated independently.

## 4. Results and Discussion

### 4.1 Energy Eigenvalues

The production calculation uses the domain $[-8,8]$, 399 interior grid points, a grid spacing of $0.04$, and five states.

| Quantum number $n$ | Exact energy | Numerical energy | Absolute error |
|---:|---:|---:|---:|
| 0 | 0.5000000000 | 0.4999499950 | $5.0005\times10^{-5}$ |
| 1 | 1.5000000000 | 1.4997499550 | $2.5005\times10^{-4}$ |
| 2 | 2.5000000000 | 2.4993498249 | $6.5018\times10^{-4}$ |
| 3 | 3.5000000000 | 3.4987495446 | $1.2505\times10^{-3}$ |
| 4 | 4.5000000000 | 4.4979490539 | $2.0509\times10^{-3}$ |

The numerical values are close to the analytical energy levels. The remaining discrepancies are primarily due to the finite-difference approximation. The absolute error increases across these five states, illustrating that higher-energy states can be more sensitive to the discretisation.

### 4.2 Wavefunction Validation

The numerical normalisation integrals for all five states are 1.000000000000 to the displayed precision.

The measured relative $L^2$ errors are:

| Quantum number $n$ | Relative $L^2$ error |
|---:|---:|
| 0 | $7.36\times10^{-5}$ |
| 1 | $2.09\times10^{-4}$ |
| 2 | $4.18\times10^{-4}$ |
| 3 | $7.12\times10^{-4}$ |
| 4 | $1.10\times10^{-3}$ |

The small errors indicate strong agreement between the numerical and analytical wavefunctions. The errors increase with quantum number for the states examined, consistent with the greater spatial structure of higher states.

### 4.3 Eigenpair Residuals

For the production calculation, the maximum residual norm among the five eigenpairs was approximately

$$
R_{\max}\approx1.81\times10^{-13}.
$$

This indicates that the computed eigenpairs satisfy the discretised matrix eigenvalue equation to high numerical accuracy. The residual is much smaller than the energy discrepancies relative to the continuous analytical solution because the residual measures the matrix eigensolver's solution of the discretised problem.

### 4.4 Grid Convergence

The grid convergence study used the fixed domain $[-4,4]$ and the following resolutions:

| Interior grid points | Grid spacing | Ground-state absolute error |
|---:|---:|---:|
| 49 | 0.16 | $8.0074\times10^{-4}$ |
| 99 | 0.08 | $1.9958\times10^{-4}$ |
| 199 | 0.04 | $4.9511\times10^{-5}$ |
| 399 | 0.02 | $1.2009\times10^{-5}$ |

The ground-state error decreases by approximately a factor of four when the grid spacing is halved. The observed orders are close to two, as expected for the central finite-difference approximation.

The higher-state errors are affected by both grid spacing and finite-domain effects. Consequently, they do not always show the same clean monotonic behaviour as the ground state.

### 4.5 Domain Convergence

At a target grid spacing of approximately $0.04$, the domain study compared $[-4,4]$, $[-6,6]$, $[-8,8]$, and $[-10,10]$.

The differences between the energy values for $[-8,8]$ and $[-10,10]$ are negligible at the displayed numerical precision. This indicates that the selected domain is sufficiently large for the first five states at this resolution.

The remaining energy discrepancies are therefore dominated primarily by finite-difference discretisation rather than the finite domain.

### 4.6 Overall Interpretation

The numerical method reproduces the low-energy spectrum and wavefunctions of the quantum harmonic oscillator with good accuracy. Energy comparisons, wavefunction normalisation, analytical wavefunction comparison, eigenpair residuals, and convergence studies provide complementary checks.

The results illustrate that a reliable computational quantum mechanics calculation requires more than obtaining eigenvalues from a numerical solver. The discretisation, spatial domain, wavefunctions, and convergence behaviour must also be examined.

## 5. Project Structure and Reproducibility

### 5.1 Project Organisation

The project is organised into separate components:

| Component | Purpose |
|---|---|
| `src/hamiltonian.py` | Defines the grid and constructs the Hamiltonian matrix |
| `src/analytical_solution.py` | Implements analytical energies and wavefunctions |
| `src/harmonic_oscillator.py` | Solves the numerical eigenvalue problem |
| `scripts/` | Contains calculation, plotting, validation, and convergence scripts |
| `tests/` | Contains automated tests |
| `results/` | Stores generated numerical data and figures |
| `theory/` | Contains theoretical documentation |
| `report/` | Contains the technical report |
| `README.md` | Gives the project overview and execution instructions |
| `requirements.txt` | Lists direct Python dependencies |

### 5.2 Dependencies

The project uses:

- **NumPy** for arrays and numerical operations.
- **SciPy** for linear algebra and Hermite polynomials.
- **Matplotlib** for visualisation.
- **pytest** for automated testing.

Install the dependencies from the project root using

```bash
python -m pip install -r requirements.txt
```

The versions in `requirements.txt` record the direct dependency versions used for the project.

### 5.3 Main Calculation

Run the production calculation from the project root:

```bash
python -m scripts.run_calculation
```

The script writes `results/energy_comparison.csv` and `results/validation_metrics.csv`.

### 5.4 Visualisation and Analysis Scripts

Generate plots:

```bash
python -m scripts.plot_wavefunctions
```

Validate the numerical wavefunctions:

```bash
python -m scripts.validate_wavefunctions
```

Run the grid convergence study:

```bash
python -m scripts.convergence_study
```

Run the domain convergence study:

```bash
python -m scripts.domain_convergence
```

Run the wavefunction convergence study:

```bash
python -m scripts.wavefunction_convergence
```

Calculate the low-lying eigenvalues using the separate eigenvalue script:

```bash
python -m scripts.solve_eigenvalues
```

These commands should be executed from the project root so that module imports and relative paths resolve correctly.

### 5.5 Automated Tests

Run the test suite:

```bash
python -m pytest tests/
```

The recorded test run completed with 10 tests passing. The tests check selected behaviour of the analytical and numerical implementation; they complement rather than replace the physical validation and convergence studies.

### 5.6 Reproducibility Configuration

The production configuration is:

| Parameter | Value |
|---|---:|
| Dimensionless domain | $[-8,8]$ |
| Interior grid points | 399 |
| Grid spacing | 0.04 |
| Number of states | 5 |
| Differential operator approximation | Central finite difference |
| Eigenvalue solver | `scipy.linalg.eigh` |

To reproduce the calculation, install the listed dependencies, run the main calculation, inspect the generated CSV files, regenerate the plots if needed, and execute the validation scripts and automated tests.
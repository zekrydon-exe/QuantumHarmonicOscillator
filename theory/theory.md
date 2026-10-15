# Quantum Harmonic Oscillator

## 1. Theoretical Model

### 1.1 Introduction

The quantum harmonic oscillator is a fundamental model in quantum mechanics. It describes a particle moving in a quadratic potential and provides an exactly solvable system against which numerical methods can be tested.

The model is relevant to molecular vibrations, quantised vibrational energy levels, and the study of more complicated quantum systems.

In this project, we solve the time-independent Schrödinger equation numerically and compare the calculated energies and wavefunctions with their analytical counterparts.

### 1.2 Potential Energy

The potential energy of a one-dimensional harmonic oscillator is

\[
V(x) = \frac{1}{2}kx^2
\]

where:

- \(x\) is the displacement from equilibrium.
- \(k\) is the force constant.
- \(V(x)\) is the potential energy.

The potential is quadratic, producing a parabolic potential-energy curve. The potential has its minimum at \(x=0\), and the restoring force increases linearly with displacement.

### 1.3 Time-Independent Schrödinger Equation

For a particle of mass \(m\) moving in one dimension, the time-independent Schrödinger equation is

\[
\hat{H}\psi(x) = E\psi(x)
\]

where \(\hat{H}\) is the Hamiltonian operator, \(\psi(x)\) is the stationary-state wavefunction, and \(E\) is the corresponding energy eigenvalue.

The Hamiltonian consists of kinetic and potential energy operators:

\[
\hat{H}
=
-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}
+
\frac{1}{2}kx^2
\]

Therefore, the Schrödinger equation becomes

\[
-\frac{\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}
+
\frac{1}{2}kx^2\psi(x)
=
E\psi(x)
\]

The first term represents the kinetic energy, while the second represents the potential energy.

### 1.4 Dimensionless Formulation

To simplify the equation, we introduce the dimensionless coordinate

\[
q = \sqrt{\frac{m\omega}{\hbar}}x
\]

where the angular frequency of the oscillator is

\[
\omega = \sqrt{\frac{k}{m}}
\]

We also define the dimensionless energy

\[
\varepsilon = \frac{E}{\hbar\omega}
\]

Using these definitions, the Schrödinger equation becomes

\[
-\frac{1}{2}\frac{d^2\psi(q)}{dq^2}
+
\frac{1}{2}q^2\psi(q)
=
\varepsilon\psi(q)
\]

Equivalently, the dimensionless Hamiltonian is

\[
\hat{H}
=
-\frac{1}{2}\frac{d^2}{dq^2}
+
\frac{1}{2}q^2
\]

and the eigenvalue equation is

\[
\hat{H}\psi_n(q)=\varepsilon_n\psi_n(q)
\]

The dimensionless formulation removes explicit dependence on the particle mass, force constant, and reduced Planck constant. It also allows the numerical results to be compared directly with the exact dimensionless energy levels.

### 1.5 Analytical Energy Eigenvalues

The one-dimensional quantum harmonic oscillator has the exact dimensionless energy eigenvalues

\[
\boxed{\varepsilon_n=n+\frac{1}{2}}
\]

where

\[
n=0,1,2,\ldots
\]

The corresponding physical energies are

\[
E_n=\hbar\omega\left(n+\frac{1}{2}\right)
\]

The ground-state energy is non-zero:

\[
\varepsilon_0=\frac{1}{2}
\]

This is the zero-point energy. It demonstrates that the oscillator cannot have zero energy in its ground state.

Adjacent energy levels are equally spaced:

\[
\varepsilon_{n+1}-\varepsilon_n=1
\]

In physical units, the spacing is \(\hbar\omega\).

### 1.6 Analytical Wavefunctions

The normalised analytical wavefunctions are

\[
\psi_n(q)
=
\frac{1}{\sqrt{2^n n!\sqrt{\pi}}}
H_n(q)e^{-q^2/2}
\]

where:

- \(n\) is the quantum number.
- \(H_n(q)\) is the physicists' Hermite polynomial of degree \(n\).
- \(e^{-q^2/2}\) is the Gaussian factor.
- The prefactor ensures normalisation.

The first few Hermite polynomials are

\[
H_0(q)=1
\]

\[
H_1(q)=2q
\]

\[
H_2(q)=4q^2-2
\]

\[
H_3(q)=8q^3-12q
\]

\[
H_4(q)=16q^4-48q^2+12
\]

The wavefunctions alternate in parity: states with even \(n\) are even functions, while states with odd \(n\) are odd functions.

The \(n\)-th state has \(n\) nodes, excluding the asymptotic approach to zero at large distances.

### 1.7 Normalisation and Orthogonality

A physically acceptable bound-state wavefunction must be normalised:

\[
\int_{-\infty}^{\infty}
|\psi_n(q)|^2\,dq=1
\]

Different energy eigenstates are orthogonal:

\[
\int_{-\infty}^{\infty}
\psi_m^*(q)\psi_n(q)\,dq
=
\delta_{mn}
\]

where \(\delta_{mn}\) is the Kronecker delta.

These properties provide useful criteria for validating the numerically calculated wavefunctions.

## 2. Numerical Methodology

### 2.1 Computational Domain and Grid Discretisation

The analytical harmonic oscillator is defined over the infinite domain

\[
-\infty < q < \infty
\]

However, numerical calculations require a finite computational domain. We therefore restrict the dimensionless coordinate to

\[
q_{\min} \leq q \leq q_{\max}
\]

The computational domain is discretised into a uniformly spaced grid. For a grid containing \(N\) interior points, the grid spacing is

\[
h=\frac{q_{\max}-q_{\min}}{N+1}
\]

The interior grid points are

\[
q_i=q_{\min}+ih,\qquad i=1,2,\ldots,N
\]

The wavefunction is represented by its values at these grid points:

\[
\psi(q_i)\approx\psi_i
\]

The wavefunction is assumed to approach zero at the boundaries of a sufficiently large computational domain. This approximation is appropriate for bound states whose probability densities become negligible near the domain boundaries.

### 2.2 Finite-Difference Approximation

The dimensionless time-independent Schrödinger equation is

\[
-\frac12\frac{d^2\psi(q)}{dq^2}
+\frac12q^2\psi(q)
=\varepsilon\psi(q)
\]

The second derivative is approximated using the central finite-difference formula:

\[
\frac{d^2\psi}{dq^2}\bigg|_{q_i}
\approx
\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{h^2}
\]

This approximation replaces the continuous second derivative with a combination of wavefunction values at neighbouring grid points.

Substituting it into the Schrödinger equation gives

\[
-\frac{\psi_{i+1}-2\psi_i+\psi_{i-1}}{2h^2}
+\frac12q_i^2\psi_i
=\varepsilon\psi_i
\]

Rearranging,

\[
-\frac{1}{2h^2}\psi_{i-1}
+
\left(\frac{1}{h^2}+\frac12q_i^2\right)\psi_i
-
\frac{1}{2h^2}\psi_{i+1}
=
\varepsilon\psi_i
\]

This equation relates the wavefunction at a grid point to its values at the two neighbouring points.

### 2.3 Construction of the Hamiltonian Matrix

The discretised Schrödinger equation can be expressed as a matrix eigenvalue problem:

\[
\mathbf{H}\boldsymbol{\psi}
=
\varepsilon\boldsymbol{\psi}
\]

Here, \(\mathbf{H}\) is the discretised Hamiltonian matrix, \(\boldsymbol{\psi}\) is an eigenvector containing the wavefunction values at the grid points, and \(\varepsilon\) is the corresponding energy eigenvalue.

The diagonal elements of the Hamiltonian are

\[
H_{ii}=\frac{1}{h^2}+\frac12q_i^2
\]

The elements immediately above and below the diagonal are

\[
H_{i,i+1}=H_{i+1,i}=-\frac{1}{2h^2}
\]

All other elements are zero.

Consequently, the Hamiltonian has a tridiagonal structure:

\[
\mathbf{H}=
\begin{pmatrix}
\frac{1}{h^2}+\frac12q_1^2 & -\frac{1}{2h^2} & 0 & \cdots & 0\\
-\frac{1}{2h^2} & \frac{1}{h^2}+\frac12q_2^2 & -\frac{1}{2h^2} & \cdots & 0\\
0 & -\frac{1}{2h^2} & \frac{1}{h^2}+\frac12q_3^2 & \ddots & \vdots\\
\vdots & \vdots & \ddots & \ddots & -\frac{1}{2h^2}\\
0 & 0 & \cdots & -\frac{1}{2h^2} & \frac{1}{h^2}+\frac12q_N^2
\end{pmatrix}
\]

The diagonal contains contributions from both kinetic and potential energy. The off-diagonal elements arise from the finite-difference approximation to the kinetic energy operator.

The matrix is real and symmetric. Its eigenvalues are therefore real, and its eigenvectors corresponding to distinct eigenvalues are orthogonal under the standard discrete inner product.

### 2.4 Numerical Eigenvalue Solution

The Hamiltonian eigenvalue problem is solved using a numerical eigensolver:

\[
\mathbf{H}\boldsymbol{\psi}_n
=
\varepsilon_n\boldsymbol{\psi}_n
\]

The eigenvalues are ordered from lowest to highest energy. The first eigenvalue corresponds to the ground state, while subsequent eigenvalues correspond to excited states.

For this project, only the lowest five eigenstates are required. Computing this subset avoids calculating the entire spectrum when only a small number of low-energy states is of interest.

The numerical eigenvalues are compared against the exact dimensionless energies:

\[
\varepsilon_n^{\mathrm{exact}}=n+\frac12
\]

The resulting differences quantify the accuracy of the numerical approximation.

### 2.5 Numerical Wavefunction Scaling

The eigensolver returns eigenvectors normalised according to the discrete Euclidean norm:

\[
\sum_{i=1}^{N}|v_i|^2=1
\]

A continuum-normalised wavefunction should instead satisfy

\[
\int_{-\infty}^{\infty}|\psi(q)|^2\,dq=1
\]

Using a uniform grid, the integral is approximated by a weighted sum. Therefore, the discrete eigenvectors are scaled according to

\[
\psi_i=\frac{v_i}{\sqrt{h}}
\]

This scaling makes the eigenvectors suitable approximations to continuum-normalised wavefunctions.

The resulting wavefunctions are checked numerically by evaluating

\[
\int |\psi_n(q)|^2\,dq
\]

using numerical integration over the computational grid.

## 3. Numerical Validation and Convergence

A numerical solution should not be considered reliable merely because an eigensolver returns eigenvalues and eigenvectors. The results must be tested against known solutions and checked for numerical consistency.

The quantum harmonic oscillator is particularly suitable for this purpose because its exact energies and analytical wavefunctions are known.

### 3.1 Energy Eigenvalue Validation

The exact dimensionless energy eigenvalues are

\[
\varepsilon_n^{\mathrm{exact}}=n+\frac12
\]

The numerical energies are obtained by diagonalising the discretised Hamiltonian. Their accuracy can be assessed using the absolute error:

\[
\Delta\varepsilon_n
=
\left|
\varepsilon_n^{\mathrm{numerical}}
-
\varepsilon_n^{\mathrm{exact}}
\right|
\]

The relative error is defined as

\[
\delta_n
=
\frac{
\left|
\varepsilon_n^{\mathrm{numerical}}
-
\varepsilon_n^{\mathrm{exact}}
\right|
}{
\left|\varepsilon_n^{\mathrm{exact}}\right|
}
\]

The relative error expressed as a percentage is

\[
\delta_n(\%)=100\delta_n
\]

These quantities measure the difference between the calculated energies and the analytical reference values.

The absolute error indicates the magnitude of the discrepancy, while the relative error expresses that discrepancy relative to the magnitude of the exact energy.

### 3.2 Wavefunction Normalisation

A wavefunction must satisfy the normalisation condition

\[
\int_{-\infty}^{\infty}|\psi_n(q)|^2\,dq=1
\]

Because the numerical wavefunctions are defined on a finite grid, the integral is evaluated numerically:

\[
I_n\approx
\int_{q_{\min}}^{q_{\max}}|\psi_n(q)|^2\,dq
\]

The numerical integration is performed using the trapezoidal rule.

A value of \(I_n\) close to unity indicates that the calculated wavefunction satisfies the normalisation condition over the chosen domain.

However, normalisation alone does not establish that a wavefunction is correct. An incorrectly shaped wavefunction can also be normalised. Therefore, the numerical wavefunctions are additionally compared with the analytical solutions.

### 3.3 Wavefunction Comparison

The analytical and numerical wavefunctions can differ by an overall sign without representing physically different states. If \(\psi_n(q)\) is an eigenfunction, then \(-\psi_n(q)\) is also an eigenfunction with the same energy.

Consequently, the numerical wavefunction is aligned in sign with the analytical wavefunction before comparison.

The overlap between the two normalised wavefunctions is

\[
S_n=
\left|
\int
\psi_n^{\mathrm{analytical}}(q)
\psi_n^{\mathrm{numerical}}(q)\,dq
\right|
\]

For perfectly matching normalised wavefunctions, the overlap is unity.

The relative \(L^2\) error provides another measure of agreement:

\[
e_n=
\frac{
\left\|
\psi_n^{\mathrm{numerical}}
-
\psi_n^{\mathrm{analytical}}
\right\|_2
}{
\left\|
\psi_n^{\mathrm{analytical}}
\right\|_2
}
\]

A smaller value indicates closer agreement between the numerical and analytical wavefunctions.

### 3.4 Eigenpair Residuals

The numerical eigenvalues and eigenvectors should satisfy the discretised eigenvalue equation

\[
\mathbf{H}\boldsymbol{\psi}_n
=
\varepsilon_n\boldsymbol{\psi}_n
\]

The residual vector is defined as

\[
\mathbf{r}_n
=
\mathbf{H}\boldsymbol{\psi}_n
-
\varepsilon_n\boldsymbol{\psi}_n
\]

Its Euclidean norm is

\[
R_n=\|\mathbf{r}_n\|_2
\]

A small residual indicates that the computed eigenpair satisfies the discretised matrix equation to high numerical accuracy.

This test assesses the numerical eigensolver and the consistency of the computed eigenpairs with the Hamiltonian matrix.

An important distinction is that a small residual does not necessarily imply that the numerical energy is close to the exact physical energy. The residual measures agreement with the discretised problem, whereas the analytical energy error measures agreement with the continuous quantum harmonic oscillator.

### 3.5 Grid Convergence

The second derivative is approximated using a central finite-difference formula. This approximation has a leading truncation error proportional to \(h^2\), where \(h\) is the grid spacing.

For a sufficiently smooth wavefunction and a fixed, adequately large computational domain, the leading discretisation error is therefore expected to scale approximately as

\[
\Delta\varepsilon_n\propto h^2
\]

When the grid spacing is reduced, the numerical approximation to the differential operator generally improves.

A grid convergence study evaluates the calculated energies using progressively finer grids while keeping the computational domain fixed.

For a second-order method, halving the grid spacing should reduce the leading discretisation error by approximately a factor of four, provided that other sources of error are sufficiently small.

The observed convergence rate can be estimated from two calculations:

\[
p\approx
\frac{
\ln\left(|\Delta\varepsilon_{n,1}|/
|\Delta\varepsilon_{n,2}|\right)
}{
\ln(h_1/h_2)
}
\]

where \(h_1\) and \(h_2\) are two grid spacings and the corresponding errors are non-zero.

A value of \(p\) close to two is consistent with second-order convergence.

### 3.6 Domain Convergence

The analytical harmonic oscillator extends over an infinite spatial domain, whereas the numerical calculation uses finite boundaries.

Restricting the domain introduces a separate source of error. If the boundaries are too close to the oscillator's centre, the wavefunction may not have decayed sufficiently before reaching them.

A domain convergence study investigates this effect by increasing the domain size while keeping the grid spacing approximately constant.

For this project, the domains examined are

\[
[-4,4],\quad[-6,6],\quad[-8,8],\quad[-10,10]
\]

The grid spacing is maintained at approximately \(0.04\).

The results are compared to determine whether increasing the domain changes the calculated low-energy eigenvalues significantly.

Once the energies become insensitive to further increases in the domain size, the remaining error is expected to be dominated by other numerical approximations, particularly the finite-difference discretisation.

### 3.7 Sources of Numerical Error

The principal sources of error in this calculation are:

1. **Finite-difference discretisation:** The continuous second derivative is approximated on a discrete grid.
2. **Finite computational domain:** The infinite spatial domain is replaced by a bounded interval.
3. **Numerical integration:** Normalisation and wavefunction comparison involve approximating continuous integrals.
4. **Floating-point arithmetic:** Numerical operations are performed with finite precision.

These sources of error should be distinguished when interpreting the results. In particular, improving the grid spacing and increasing the domain size address different limitations and should be investigated separately.

## 4. Results and Discussion

### 4.1 Energy Eigenvalues

The lowest five eigenvalues were calculated by diagonalising the discretised Hamiltonian over the dimensionless domain \([-8,8]\), using 399 interior grid points and a grid spacing of \(0.04\).

The numerical eigenvalues were compared with the exact harmonic oscillator energies:

\[
\varepsilon_n^{\mathrm{exact}}=n+\frac12
\]

The calculated values are summarised below.

| Quantum number \(n\) | Exact energy | Numerical energy |
|---:|---:|---:|
| 0 | 0.500000 | 0.499950 |
| 1 | 1.500000 | 1.499750 |
| 2 | 2.500000 | 2.499350 |
| 3 | 3.500000 | 3.498750 |
| 4 | 4.500000 | 4.497949 |

The numerical energies are close to the analytical values. The discrepancies arise primarily from the finite-difference approximation to the kinetic energy operator.

The absolute energy error increases across the five states. This is consistent with the greater sensitivity of higher-energy states to the discretisation, although the precise error also depends on the grid spacing and computational domain.

### 4.2 Wavefunction Validation

The calculated wavefunctions were checked for normalisation and compared with the corresponding analytical solutions.

For all five states, the numerical normalisation integrals were approximately unity to the displayed precision.

The relative \(L^2\) errors obtained from the wavefunction comparison were:

| Quantum number \(n\) | Relative \(L^2\) error |
|---:|---:|
| 0 | \(7.36\times10^{-5}\) |
| 1 | \(2.09\times10^{-4}\) |
| 2 | \(4.18\times10^{-4}\) |
| 3 | \(7.12\times10^{-4}\) |
| 4 | \(1.10\times10^{-3}\) |

The small relative errors indicate strong agreement between the numerical and analytical wavefunctions.

The error increases with the quantum number for the five states examined. Higher-energy states exhibit more spatial structure, including additional nodes, which makes accurate representation on a finite grid more demanding.

The numerical and analytical wavefunctions therefore provide complementary validation: the energy comparison assesses the spectrum, while the wavefunction comparison assesses the calculated eigenstates.

### 4.3 Eigenpair Residuals

The numerical eigenpairs were checked using the residual

\[
\mathbf{r}_n
=
\mathbf{H}\boldsymbol{\psi}_n
-
\varepsilon_n\boldsymbol{\psi}_n
\]

The residual norms were calculated using the same discretised Hamiltonian that was used to obtain the eigenvalues.

Small residual norms indicate that the computed eigenvalues and eigenvectors satisfy the matrix eigenvalue equation to high numerical accuracy.

For the production calculation, the maximum residual norm among the five
computed eigenpairs was approximately

\[
R_{\max}\approx 1.81\times10^{-13}
\]

This small value indicates that the computed eigenpairs satisfy the discretised Hamiltonian eigenvalue equation to high numerical accuracy. It does not imply that the discretisation error relative to the exact continuous problem is equally small.

### 4.4 Grid Convergence

The grid convergence study examined the effect of changing the grid resolution while maintaining the domain \([-4,4]\).

For the ground state, the absolute energy errors decreased as follows:

| Grid points | Grid spacing | Ground-state absolute error |
|---:|---:|---:|
| 49 | 0.16 | \(8.01\times10^{-4}\) |
| 99 | 0.08 | \(2.00\times10^{-4}\) |
| 199 | 0.04 | \(4.95\times10^{-5}\) |
| 399 | 0.02 | \(1.20\times10^{-5}\) |

The ground-state error decreases by approximately a factor of four when the grid spacing is halved. This is consistent with the second-order accuracy expected from the central finite-difference approximation.

The behaviour of higher-energy states must be interpreted with additional care because their errors can be influenced by the finite computational domain as well as the grid spacing.

### 4.5 Domain Convergence

The domain convergence study compared the first five energy eigenvalues over the intervals \([-4,4]\), \([-6,6]\), \([-8,8]\), and \([-10,10]\), maintaining a grid spacing of approximately \(0.04\).

The results changed appreciably between the smaller domains, particularly for the higher-energy states. However, the calculated energies became effectively unchanged between \([-8,8]\) and \([-10,10]\).

This indicates that the finite-domain effect is negligible for the first five states at the selected resolution when the domain is \([-8,8]\).

The remaining discrepancy from the exact energies is therefore primarily associated with the finite-difference discretisation rather than the choice of domain.

### 4.6 Interpretation

The results demonstrate that the finite-difference method can reproduce the low-energy spectrum and wavefunctions of the quantum harmonic oscillator with good accuracy.

The combination of energy comparisons, wavefunction validation, eigenpair residuals, grid convergence, and domain convergence provides a more complete assessment than any single test.

The study also illustrates an important principle in computational quantum mechanics: numerical accuracy depends on both the discretisation of the differential operator and the representation of the spatial domain.

Although the calculated results agree closely with the analytical solution, further improvements in grid resolution would be required to reduce the remaining discretisation errors.

## 5. Project Structure and Reproducibility

### 5.1 Project Organisation

The project separates the mathematical model, numerical implementation, validation procedures, and generated results into distinct components.

The main directories and files serve the following purposes:

| Component | Purpose |
|---|---|
| `src/` | Contains the core numerical implementation |
| `src/hamiltonian.py` | Defines the computational grid and constructs the Hamiltonian matrix |
| `src/analytical_solution.py` | Implements the exact energy eigenvalues and analytical wavefunctions |
| `src/harmonic_oscillator.py` | Solves the discretised eigenvalue problem and returns the numerical eigenstates |
| `scripts/` | Contains scripts for calculations, validation, plotting, and convergence studies |
| `tests/` | Contains automated tests for the implementation |
| `results/` | Stores generated plots and numerical data |
| `theory/` | Contains theoretical documentation |
| `report/` | Reserved for the final project report |
| `README.md` | Provides an overview of the project and instructions for use |
| `requirements.txt` | Lists the Python dependencies required by the project |

### 5.2 Software and Dependencies

The implementation uses Python and the following scientific computing libraries:

- **NumPy:** Numerical arrays and mathematical operations.
- **SciPy:** Linear algebra routines, including the symmetric eigensolver, and evaluation of Hermite polynomials.
- **Matplotlib:** Visualisation of wavefunctions and probability densities.
- **pytest:** Automated testing of the implementation.

The Python environment should contain compatible versions of these dependencies.

Install the dependencies using

```bash
python -m pip install -r requirements.txt
# SQUID Inverse Holography

## Mathematical Formulation of Inverse Matterwave Holography using SQUID Arrays

In the forward matterwave holography problem (developed in the companion document), we compute the probability density at a target plane given a known magnetic flux distribution across a SQUID array. The **inverse holographic mapping** problem reverses this: given a desired 2D particle deposition pattern — a target probability density $P_{\text{target}}(x,y)$ — compute the magnetic flux required at each SQUID loop to sculpt a coherent matterwave beam into that pattern.

The inverse problem decomposes into two mathematically distinct stages:

1. **Phase retrieval**: Finding the continuous phase map $\phi_{\text{grid}}(x,y)$ at the SQUID array plane that, under forward Fresnel propagation, produces the desired target intensity.
2. **Hardware inversion**: Translating that phase map into discrete drive currents for a physical SQUID array, accounting for mutual inductance, discretization, and device constraints.

This separation is both conceptually clean and practically necessary — the physics of wave propagation and the physics of superconducting circuits are independent problems that couple only through the Aharonov–Bohm phase relation.

---

### 1. The Target Constraint and the Free Phase Degree of Freedom

Let the target (deposition) surface be located at distance $z$ from the SQUID array. The desired deposition pattern is specified as a normalized probability density $P_{\text{target}}(x,y) \geq 0$.

The amplitude of the matterwave at the target plane is strictly determined by this density:

$$A_{\text{target}}(x,y) = \sqrt{P_{\text{target}}(x,y)}$$

However, the **phase** at the target plane, $\phi_{\text{target}}(x,y)$, is physically unconstrained: particle deposition (or any intensity-based detection) depends only on $|\psi|^2$, not on the phase of $\psi$. Two wavefields with identical amplitude but different phase profiles produce the same deposition pattern.

This unconstrained target phase is not a nuisance — it is the **essential free degree of freedom** that makes the inverse problem solvable. Without it, the problem would be overconstrained: a phase-only modulator has $N^2$ real parameters (one flux per pixel), while a complex target field with fixed phase would impose $2N^2$ real constraints (amplitude and phase at each target point). The free target phase reduces the effective constraint count to $N^2$ (amplitude only), matching the number of degrees of freedom.

The target wavefunction is therefore:

$$\psi_{\text{target}}(x,y) = A_{\text{target}}(x,y)\, e^{i\phi_{\text{target}}(x,y)}$$

where $\phi_{\text{target}}$ is unknown and must be determined as part of the solution.

---

### 2. The Propagation Operators

Before describing the retrieval algorithm, we define the forward and backward propagation operators that it relies on. Both are derived from the paraxial free-particle Schrödinger equation and are mathematically identical to the Fresnel diffraction integrals of scalar optics (see the companion forward-propagation document).

**Forward propagator** $\mathcal{P}_{+z}$: maps the field at $z = 0$ to the field at $z$:

$$\mathcal{P}_{+z}\{\psi\}(x,y) = \frac{e^{ikz}}{i\lambda z} \iint \psi(x',y')\, \exp\left[\frac{i\pi}{\lambda z}\left((x-x')^2 + (y-y')^2\right)\right] dx'\, dy'$$

**Backward propagator** $\mathcal{P}_{-z}$: the adjoint (conjugate) of the forward propagator, mapping the field at $z$ back to $z = 0$:

$$\mathcal{P}_{-z}\{\psi\}(x,y) = \frac{e^{-ikz}}{-i\lambda z} \iint \psi(x',y')\, \exp\left[-\frac{i\pi}{\lambda z}\left((x-x')^2 + (y-y')^2\right)\right] dx'\, dy'$$

These operators are exact inverses: $\mathcal{P}_{-z}\{\mathcal{P}_{+z}\{\psi\}\} = \psi$ for any square-integrable $\psi$. However, this exact inversion requires the **complete complex field** — both amplitude and phase. Since we know only $A_{\text{target}}$ and not $\phi_{\text{target}}$, a single back-propagation does not yield the solution. Instead, we use the operators iteratively within a constrained optimization loop.

---

### 3. Iterative Phase Retrieval

#### 3.1 The Gerchberg–Saxton Algorithm

The classical approach to this class of problem is the Gerchberg–Saxton (GS) algorithm, which alternates between the two planes, enforcing the known physical constraint at each:

- **At $z = 0$ (SQUID grid):** The modulator is phase-only, so the amplitude must equal the uniform incident beam amplitude $A_{\text{in}} = \sqrt{\rho_0}$.
- **At $z$ (target plane):** The amplitude must match $A_{\text{target}}(x,y) = \sqrt{P_{\text{target}}(x,y)}$.

Let $n$ denote the iteration index. Initialize with a random phase guess $\phi_{\text{target}}^{(0)}(x,y)$ (drawn uniformly from $[0, 2\pi)$ at each point).

**Step A — Back-propagate to the grid plane:**

$$\psi_{\text{grid}}^{(n)}(x,y) = \mathcal{P}_{-z}\left\{A_{\text{target}}(x,y)\, e^{i\phi_{\text{target}}^{(n-1)}(x,y)}\right\}$$

**Step B — Enforce the grid constraint (phase-only modulation):**

Extract the phase: $\phi_{\text{grid}}^{(n)}(x,y) = \arg\!\left(\psi_{\text{grid}}^{(n)}(x,y)\right)$.

Replace the amplitude with the known incident amplitude:

$$\hat{\psi}_{\text{grid}}^{(n)}(x,y) = A_{\text{in}}\, e^{i\phi_{\text{grid}}^{(n)}(x,y)}$$

**Step C — Forward-propagate to the target plane:**

$$\psi_{\text{target}}^{(n)}(x,y) = \mathcal{P}_{+z}\left\{\hat{\psi}_{\text{grid}}^{(n)}(x,y)\right\}$$

**Step D — Enforce the target constraint (desired intensity):**

Extract the phase: $\phi_{\text{target}}^{(n)}(x,y) = \arg\!\left(\psi_{\text{target}}^{(n)}(x,y)\right)$.

Replace the amplitude with the target amplitude:

$$\psi_{\text{target}}^{(n)}(x,y) = A_{\text{target}}(x,y)\, e^{i\phi_{\text{target}}^{(n)}(x,y)}$$

Repeat from Step A until convergence.

#### 3.2 Convergence Criteria

The algorithm requires a quantitative measure of solution quality. Common error metrics include:

**Normalized root-mean-square error (NRMSE):**

$$\varepsilon = \frac{\left[\iint \left(|\psi_{\text{target}}^{(n)}|^2 - P_{\text{target}}\right)^2 dx\, dy\right]^{1/2}}{\left[\iint P_{\text{target}}^2\, dx\, dy\right]^{1/2}}$$

**Diffraction efficiency** (fraction of total flux landing in the desired pattern region $\Omega$):

$$\eta = \frac{\iint_\Omega |\psi_{\text{target}}^{(n)}|^2\, dx\, dy}{\iint |\psi_{\text{target}}^{(n)}|^2\, dx\, dy}$$

**Intensity uniformity** within the target region (relevant for flat-top patterns):

$$U = 1 - \frac{\max(|\psi|^2) - \min(|\psi|^2)}{\max(|\psi|^2) + \min(|\psi|^2)} \quad \text{over } \Omega$$

The iteration terminates when the chosen metric crosses a predefined threshold (e.g., $\varepsilon < 0.01$) or when successive iterations produce negligible improvement (stagnation detection).

#### 3.3 Limitations and Improved Variants

The basic GS algorithm has well-documented limitations:

**Local optima.** GS performs alternating projections between two constraint sets. While each projection reduces (or maintains) the error with respect to its own constraint, the algorithm is not guaranteed to find the globally optimal phase. It converges to a local optimum that depends on the initial random phase guess. In practice, the algorithm is run from multiple independent random initializations, and the best result is selected.

**Stagnation.** For complex target patterns — particularly those with sharp edges, high dynamic range, or isolated bright features surrounded by dark regions — GS can stagnate far from an acceptable solution. The amplitude replacement step (Step D) is aggressive: it discards all amplitude information from the propagated field, which can slow convergence or trap the algorithm.

**Speckle and non-uniformity.** Because the target phase is unconstrained, the converged solution typically exhibits rapid random phase variations across the target plane. These produce granular intensity fluctuations (speckle) within the nominally bright regions of the pattern.

Several improved algorithms address these issues:

- **Weighted Gerchberg–Saxton (iterative Fourier transform algorithm, IFTA):** Instead of hard amplitude replacement at the target, the amplitude is nudged toward $A_{\text{target}}$ using a weighted feedback rule: $A^{(n)} \to A^{(n)} + \alpha(A_{\text{target}} - A^{(n)})$, where $\alpha < 1$. This gentler update improves uniformity at the cost of slower convergence.

- **Adaptive-additive method:** Adjusts the feedback weight $\alpha$ spatially and per-iteration based on the local error, applying stronger corrections where the intensity deviates most from the target.

- **Gradient descent / adjoint optimization:** Treats the phase retrieval as a continuous optimization problem with a differentiable cost function (e.g., mean squared error between $|\psi|^2$ and $P_{\text{target}}$). The gradient of the cost with respect to each phase pixel can be computed analytically through the Fresnel integral, enabling gradient-based optimizers (L-BFGS, Adam, etc.) to search the parameter space more effectively than alternating projections. This approach naturally accommodates additional regularization terms (e.g., smoothness penalties on the phase to reduce speckle, or constraints on the maximum phase gradient to respect hardware limits).

---

### 4. Phase-to-Flux Translation

#### 4.1 The Aharonov–Bohm Mapping

Once the phase retrieval algorithm converges, it delivers a grid-plane phase map $\phi_{\text{grid}}(x,y)$. This phase represents the *total* phase of the modulated beam at $z = 0$. The AB phase shift that the SQUID array must impose is the difference between this retrieved phase and the phase of the unmodulated incident beam:

$$\Delta\phi(x,y) = \phi_{\text{grid}}(x,y) - \phi_{\text{in}}(x,y)$$

For a normally incident plane wave ($\phi_{\text{in}} = \text{const}$ across the array), this simplifies to $\Delta\phi = \phi_{\text{grid}}$ up to a global constant (which has no physical effect on the intensity pattern).

The Aharonov–Bohm relation then gives the required local magnetic flux:

$$\Phi(x,y) = \frac{\hbar}{q}\,\Delta\phi(x,y)$$

where $q$ is the charge of the beam particles (e.g., $q = e$ for electrons). This expression uses the general charge $q$ rather than the electron charge specifically, since the framework applies to any coherent charged-particle beam.

#### 4.2 Phase Wrapping and Flux Quantization

The AB phase shift is periodic: a flux change of one flux quantum $\Phi_0 = h/q$ produces a phase shift of $2\pi$, which is physically indistinguishable from zero. The required flux can therefore always be reduced modulo $\Phi_0$:

$$\Phi_{\text{reduced}}(x,y) = \Phi(x,y) \mod \Phi_0 = \frac{\hbar}{q}\left[\Delta\phi(x,y) \mod 2\pi\right]$$

This is practically important: the phase retrieval algorithm may produce phase values spanning many multiples of $2\pi$, but the physical flux at each SQUID need never exceed a single flux quantum. For electrons, $\Phi_0 = h/e \approx 4.14 \times 10^{-15}\,\text{Wb}$, which corresponds to modest drive currents in a superconducting loop with typical inductances.

Phase wrapping also has implications for the smoothness of the drive current distribution. A slowly varying phase $\Delta\phi(x,y)$ that crosses a $2\pi$ boundary produces a sharp discontinuity in $\Phi_{\text{reduced}}$, which translates to an abrupt jump in drive current between neighboring pixels. The hardware must be able to accommodate these transitions without cross-talk artifacts (see Section 5).

---

### 5. Hardware Inversion: From Flux to Drive Currents

#### 5.1 The Mutual Inductance Model

The final step translates the continuous flux distribution into discrete drive currents for a physical SQUID array of $N$ loops.

Discretize the target flux onto the array grid: $\Phi_i \equiv \Phi(x_i, y_i)$ for $i = 1, \ldots, N$. In a coupled superconducting circuit, the flux threading the $i$-th loop is not simply proportional to its own drive current — it receives contributions from every other loop via mutual inductance:

$$\Phi_i = \sum_{j=1}^{N} M_{ij}\, I_j$$

where $M_{ij}$ is the mutual inductance between loops $i$ and $j$ (for $i \neq j$) and $M_{ii} \equiv L_i$ is the self-inductance of loop $i$.

In matrix form:

$$\boldsymbol{\Phi} = \mathbf{M}\, \mathbf{I}$$

The inductance matrix $\mathbf{M}$ is symmetric ($M_{ij} = M_{ji}$ by reciprocity) and positive-definite (since the magnetic energy stored in any nonzero current distribution, $\frac{1}{2}\mathbf{I}^T \mathbf{M}\, \mathbf{I} > 0$, must be positive). Positive-definiteness guarantees that $\mathbf{M}$ is invertible, so the required drive currents are:

$$\mathbf{I} = \mathbf{M}^{-1}\, \boldsymbol{\Phi}$$

#### 5.2 Structure and Scaling of the Inductance Matrix

For a regular planar array, the mutual inductance between two loops depends primarily on their separation $r_{ij}$ and falls off as $M_{ij} \sim r_{ij}^{-3}$ for loops small compared to their spacing. This gives $\mathbf{M}$ a characteristic structure: it is **diagonally dominant** (self-inductance $L_i$ is much larger than any single off-diagonal term) and **approximately banded** (distant loops contribute negligibly).

This structure has computational implications. For a large array with $N = 10^4$ or more pixels:

- **Direct inversion** ($O(N^3)$) is expensive but feasible for moderate $N$ and need only be performed once for a given array geometry.
- **Iterative solvers** (conjugate gradient, GMRES) exploit the diagonal dominance and sparsity pattern to solve $\mathbf{M}\,\mathbf{I} = \boldsymbol{\Phi}$ in $O(N^2)$ or better per iteration, with rapid convergence.
- **Truncated approximation:** Setting $M_{ij} = 0$ for $|r_{ij}| > r_{\text{cut}}$ yields a sparse matrix that can be factored efficiently. The truncation error is small if $r_{\text{cut}}$ spans several pixel pitches.

In all cases, the matrix $\mathbf{M}$ depends only on the array geometry and can be precomputed and stored.

#### 5.3 Nonlinear Device Effects

The linear model $\Phi_i = \sum M_{ij} I_j$ assumes that each SQUID operates as a simple linear inductance. In practice, real SQUIDs exhibit:

- **Periodic flux-to-voltage transfer functions:** A SQUID's response is periodic in $\Phi_0$, so the operating point must be chosen to lie within a single monotonic branch of the transfer function, or the nonlinearity must be calibrated and compensated.
- **Flux trapping:** Ambient magnetic flux or thermal cycling can trap vortices in the superconducting material, producing fixed offsets in the flux of individual loops. These offsets shift the effective phase origin of each pixel and must be measured and subtracted during calibration.
- **Hysteresis and noise:** Practical SQUIDs have finite flux noise (typically $\sim 10^{-6}\,\Phi_0 / \sqrt{\text{Hz}}$ for state-of-the-art devices) and can exhibit hysteretic behavior near critical current. Both effects set a floor on the phase accuracy achievable at each pixel.

A full hardware model would replace the linear matrix equation with a nonlinear system $\boldsymbol{\Phi} = \mathbf{f}(\mathbf{I})$ and solve it via Newton iteration or similar methods, using the linear solution as an initial guess.

---

### 6. Discretization and Sampling Constraints

The continuous flux distribution $\Phi(x,y)$ emerging from the phase retrieval algorithm must be sampled onto the finite pixel grid of the SQUID array. This discretization imposes fundamental limits on the holographic output.

#### 6.1 Spatial Bandwidth

If the array has pixel pitch $d$ (center-to-center spacing), the maximum spatial frequency the phase screen can represent is:

$$f_{\max} = \frac{1}{2d}$$

by the Nyquist–Shannon sampling theorem. This limits the maximum diffraction angle:

$$\theta_{\max} = \arcsin(\lambda f_{\max}) = \arcsin\!\left(\frac{\lambda}{2d}\right)$$

and therefore the **field of view** at the target plane:

$$\text{FOV} \approx 2z\tan\theta_{\max} \approx \frac{\lambda z}{d}$$

for small angles. Any features in the target pattern that require spatial frequencies beyond $f_{\max}$ are lost — they either alias into the representable band or are simply absent from the reconstructed pattern.

#### 6.2 Space-Bandwidth Product

The total information capacity of the holographic system is characterized by its **space-bandwidth product** (SBP):

$$\text{SBP} = \left(\frac{D}{d}\right)^2 = N$$

where $D$ is the total array aperture and $N$ is the number of SQUID pixels. The SBP equals the number of independent complex degrees of freedom the system can control, and it sets an upper bound on the complexity (number of resolvable spots) in the holographic output. A target pattern requiring more than $N$ resolvable features cannot be faithfully reproduced by the array.

---

### 7. Diffraction Efficiency and Energy Budget

Not all incident particles contribute to the desired pattern. The **diffraction efficiency** $\eta$ quantifies the fraction of the total beam flux that lands within the target region:

$$\eta = \frac{\iint_\Omega P(x,y,z)\, dx\, dy}{\iint P(x,y,z)\, dx\, dy}$$

where $\Omega$ is the spatial region containing the desired pattern and the denominator integrates over the entire target plane.

For a phase-only modulator, the theoretical upper bound on $\eta$ depends on the target pattern:

- **Simple periodic structures** (gratings, arrays of spots): $\eta$ can approach 100%, since the phase profile can be designed to direct nearly all flux into the desired orders.
- **Complex images** with fine detail and high dynamic range: $\eta$ typically falls to 50–80%, with the remainder scattered into the zero-order spot, the conjugate (twin) image, and higher-order diffraction artifacts.
- **Isolated small features** (e.g., a single focused spot): $\eta$ can be very high, since this is equivalent to a lens — a well-known optimal use case for phase-only elements.

The lost flux is not absorbed (the modulator is lossless) but is redistributed to undesired locations on the target plane. Understanding the efficiency budget is essential for predicting exposure times in a deposition application: if $\eta = 0.6$, the required beam dose is $1/0.6 \approx 1.67$ times the dose that would be needed with a perfect (amplitude-and-phase) modulator.

---

### 8. Summary of the Inverse Pipeline

The complete inverse holographic computation proceeds as follows:

1. **Specify the target** $P_{\text{target}}(x,y)$: the desired particle deposition pattern at distance $z$.

2. **Retrieve the phase** $\phi_{\text{grid}}(x,y)$: using the Gerchberg–Saxton algorithm or an improved variant, iteratively propagate between the grid and target planes, enforcing the phase-only constraint at the grid and the amplitude constraint at the target, until the error metric converges.

3. **Convert phase to flux** $\Phi(x,y) = (\hbar/q)\,\Delta\phi(x,y)$: apply the Aharonov–Bohm relation, reducing modulo $\Phi_0$ to keep the flux within a single quantum.

4. **Discretize and sample** onto the $N$-pixel SQUID grid, respecting the Nyquist limit set by the pixel pitch.

5. **Invert the hardware** $\mathbf{I} = \mathbf{M}^{-1}\boldsymbol{\Phi}$: solve the coupled inductance problem to obtain the drive current for each SQUID loop, accounting for mutual inductance, device nonlinearities, and calibration offsets.

6. **Apply the currents** to the SQUID array and expose the target to the modulated matterwave beam.

Each stage introduces its own error sources — algorithmic (phase retrieval stagnation, speckle), discretization (aliasing, finite SBP), and hardware (flux noise, crosstalk calibration, nonlinearity). A complete system design must budget these errors jointly to ensure the final deposition pattern meets the required fidelity.
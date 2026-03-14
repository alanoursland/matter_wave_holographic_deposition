## Quantum Substrate Deposition & Topological Caging

**Mathematical Framework for Spatial Floquet Engineering and Lieb Lattice Localization**

---

### 1. Coherent Beam Preparation

The atomic beam is modeled as a massive bosonic field $\psi(\mathbf{r})$ governed by the time-independent Schrödinger equation. Global phase synchronization of the atomic ensemble is described by the **Kuramoto order parameter**, $r$, defining the complex coherence factor $\mathcal{C}$:


$$r e^{i\theta} = \frac{1}{N} \sum_{j=1}^{N} e^{i\theta_j}$$


The initial wave packet is modulated by a spatial coherence envelope determined by the beam coherence length $L_c$:


$$\psi_{0}(\mathbf{r}) \propto \exp\left(-\frac{r^2}{2\sigma^2}\right) \exp(i k_0 x) \exp\left(-\frac{r^2}{2L_c^2}\right)$$

---

### 2. Spatial Phase Imprinting (Aharonov-Bohm Analogue)

The substrate imprints a geometric phase $\phi(\mathbf{r})$ via a synthetic gauge field. For a hexagonal vortex lattice, the phase accumulation around a singularity is given by:


$$\phi(\mathbf{r}) = \sum_{i} s_i \arctan\left(\frac{y-y_i}{x-x_i}\right) \cdot \text{env}(\mathbf{r})$$


where $s_i = \pm 1$ is the vortex sign and $\text{env}(\mathbf{r})$ is the **cosine-squared apodization envelope** applied to eliminate boundary discontinuities at the substrate perimeter $L$:


$$\text{env}(x) = \begin{cases} \cos^2\left(\frac{\pi}{2} \frac{|x| - (L/2 - m)}{m}\right) & |x| > L/2 - m \\ 1.0 & \text{otherwise} \end{cases}$$

---

### 3. Spatial Floquet Dressing

The transverse potential $V(\mathbf{r}) \propto |\phi(\mathbf{r})|$ acts as a periodic drive in the spatial domain. The state is decomposed into **Floquet sidebands** $|n\rangle$ with energies $E_n = n \hbar \omega$. The effective Hamiltonian in the Floquet basis is:


$$H_{F} = \sum_{n} n\hbar\omega |n\rangle\langle n| + \frac{V(\mathbf{r})}{2} \sum_{n} (|n\rangle\langle n+1| + |n+1\rangle\langle n|)$$


The local occupancy of the $n$-th sideband, $|c_n(\mathbf{r})|^2$, determines the deposition probability for a species resonant with that specific energy level.

---

### 4. Selective Adsorption Filter

Atomic deposition is modeled as a Lorentzian binding filter centered on a resonant sideband $n_{res}$:


$$\psi_{ads}(\mathbf{r}) = \sum_{n} \left[ \frac{\Gamma^2}{(n - n_{res})^2 + \Gamma^2} \right] c_n(\mathbf{r}) \psi_{prop}(\mathbf{r})$$


This allows for **multi-species deterministic placement**, where different atomic species A and B are directed to inter-vortex ($n=0$) or vortex-core ($n=\pm 2$) regions respectively.

---

### 5. Topological Caging on the Lieb Lattice

The deposited pattern is loaded onto a 2D Lieb lattice. To achieve Aharonov-Bohm (AB) caging, a **staggered Landau gauge** is implemented to ensure a net flux of $\Phi = \pi$ per plaquette:


$$H_{Lieb} = -J \sum_{\langle i,j \rangle} e^{i\theta_{ij}} \hat{c}_i^\dagger \hat{c}_j + \sum_{i} \delta\epsilon_i \hat{n}_i$$

* **Gauge Condition:** $\sum_{plaquette} \theta_{ij} = \pi$.
* **Disorder Robustness:** On-site disorder $\delta\epsilon \in [-W/2, W/2]$ is suppressed by the topological flat bands, maintaining pattern fidelity $\mathcal{F}$:

$$\mathcal{F}(t) = |\langle \psi(0) | \psi(t) \rangle|^2$$



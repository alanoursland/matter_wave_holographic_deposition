# Hologram Basics: A Mathematical Foundation

## 1. The Core Problem: Recording Phase

All optical detectors — film, CCD sensors, the human eye — respond to **intensity**, the time-averaged energy flux of a light wave. Intensity is proportional to the squared magnitude of the electric field and is always a real, positive number. This means that when you photograph a scene, you capture the brightness at each point but irretrievably lose the **phase** — the information about where each wave was in its oscillation cycle at the moment of detection.

Phase is what encodes three-dimensional structure. When you look at a real object, light waves scattered from its near surfaces arrive at your eye with different phases than light from its far surfaces, because they have traveled different distances. Your visual system interprets these phase relationships (via interference at the retina and binocular disparity) as depth. A photograph collapses all of that into a flat intensity map.

Holography is the technique for encoding phase information into an intensity pattern, so that a conventional recording medium can capture the full complex wavefield — amplitude *and* phase — and later reconstruct it.

## 2. Mathematical Representation of the Wavefield

We describe a monochromatic (single-frequency) wave at a recording plane $(x, y)$ using **complex amplitude**:

$$U(x,y) = A(x,y)\,e^{i\phi(x,y)}$$

where $A(x,y) \geq 0$ is the real amplitude and $\phi(x,y)$ is the phase. The physical wave is understood to be $\mathrm{Re}\{U\, e^{-i\omega t}\}$, but because all waves in the system share the same temporal frequency $\omega$, the time dependence cancels in every calculation that follows and we can work entirely with the spatial part $U(x,y)$.

A holographic system requires two mutually coherent waves at the recording plane:

- **The Object Wave** $O(x,y) = o(x,y)\,e^{i\phi_o(x,y)}$: the complex field scattered or diffracted by the target. This is the wavefield we want to record and later reconstruct.

- **The Reference Wave** $R(x,y)$: a known, controlled wave, typically a plane wave or spherical wave. For a plane wave arriving at angle $\theta$ to the recording plane normal, the reference takes the form:

$$R(x,y) = r\, e^{i 2\pi f_x x}$$

where $r$ is a constant real amplitude and $f_x = \frac{\sin\theta}{\lambda}$ is the **spatial carrier frequency** — the number of wavefront cycles per unit length across the recording plane.

## 3. Coherence Requirements

The interference central to holography demands that the object and reference waves maintain a stable phase relationship. This imposes two practical constraints:

**Temporal coherence.** The light source must have a sufficiently narrow spectral bandwidth $\Delta\lambda$ that the coherence length $l_c \approx \lambda^2 / \Delta\lambda$ exceeds the maximum optical path difference between the reference and object arms. For a He-Ne laser ($\lambda = 632.8$ nm, $\Delta\lambda \sim 0.001$ nm), $l_c$ can exceed 100 meters. An ordinary white-light source, by contrast, has $l_c$ on the order of micrometers — far too short for most holographic geometries.

**Spatial coherence.** The source must approximate a point emitter (or a single spatial mode) so that wavefronts from different parts of the source interfere constructively at the recording plane. Extended sources wash out fringe visibility. In practice, this is satisfied by using a single-mode laser.

Without adequate coherence, the interference fringes that encode phase information have reduced contrast or vanish entirely, and the hologram fails.

## 4. The Recording Process: Interference

When the object and reference waves superimpose, the total field at the recording plane is $U_{\text{total}} = O + R$. The recording medium captures the intensity:

$$I(x,y) = |O + R|^2 = (O + R)(O + R)^*$$

Expanding:

$$\boxed{I(x,y) = |O|^2 + |R|^2 + OR^* + O^*R}$$

This is the fundamental equation of holographic recording. Each term has a distinct physical role:

| Term | Name | Physical meaning |
|------|------|-----------------|
| $\|O\|^2$ | Object self-interference | Intensity of the object wave alone — a blurred halo with no useful phase content. |
| $\|R\|^2$ | Reference bias | A uniform (DC) background if $R$ is a plane wave. |
| $OR^*$ | Object term | The key term. It encodes the full complex amplitude of $O$, modulated onto the carrier frequency of $R$. |
| $O^*R$ | Conjugate term | The complex conjugate of the object wave, modulated onto the opposite sideband. |

**Why does an intensity pattern encode phase?** Consider the $OR^*$ term explicitly. If $O = o\,e^{i\phi_o}$ and $R = r\,e^{i2\pi f_x x}$, then:

$$OR^* = o\, r\, e^{i(\phi_o - 2\pi f_x x)}$$

The real part of this (which contributes to the total intensity) is $o\, r\, \cos(\phi_o - 2\pi f_x x)$. This is a set of fringes at spatial frequency $f_x$, but **locally shifted** by the object phase $\phi_o(x,y)$. Where the object phase advances, the fringes shift one way; where it retards, they shift the other. The phase has been converted into a spatial pattern of fringe positions — a form that intensity-sensitive media can record.

## 5. Resolution Requirements of the Recording Medium

The interference pattern from the previous section contains spatial frequencies up to $f_x + B$, where $B$ is the spatial bandwidth of the object wave. The recording medium must be able to resolve these frequencies, which means its resolution (in line pairs per mm) must satisfy:

$$f_{\text{medium}} \geq f_x + B = \frac{\sin\theta}{\lambda} + B$$

For visible light ($\lambda \approx 500$ nm) at a moderate off-axis angle ($\theta = 30°$), the carrier frequency alone is $f_x = \sin 30° / 0.5\,\mu\text{m} = 1000\,\text{lp/mm}$. Conventional photographic film resolves perhaps 100–200 lp/mm. Holographic recording therefore requires specialized high-resolution media — silver halide emulsions with grain sizes below 50 nm, photopolymers, or photorefractive crystals — capable of resolving several thousand line pairs per millimeter.

This resolution requirement is a direct physical consequence of the mathematical structure of the interference equation and represents one of the primary engineering constraints in holography.

## 6. The Amplitude Transmittance

After exposure and processing, the recording medium acts as a diffractive element. Assuming a **linear** recording regime — meaning the exposure falls within the medium's linear response range — the amplitude transmittance of the developed hologram is:

$$t(x,y) = t_b + \beta\, I(x,y)$$

where $t_b$ is the background transmittance and $\beta$ is a real constant determined by the medium's sensitivity and the exposure time.

The linear assumption is important and limited. Real media saturate at high intensities and have a threshold at low intensities. If the object-to-reference intensity ratio is too high, or if the total exposure leaves the linear regime, intermodulation terms ($|O|^2$ cross-products) generate spurious diffraction orders that degrade image quality. In practice, the reference beam is typically made 3–10 times brighter than the object beam to keep the modulation depth within the linear range.

For the derivations that follow, we adopt the idealized linear case and set $t(x,y) \propto I(x,y)$ for clarity.

## 7. The Reconstruction Process: Diffraction

To reconstruct the stored wavefield, the hologram is illuminated by a **reconstruction wave** $C(x,y)$. The field transmitted through the hologram is the product of the illuminating wave and the transmittance:

$$U_{\text{trans}}(x,y) = C(x,y) \cdot t(x,y) \propto C(x,y) \cdot I(x,y)$$

Substituting the expanded intensity:

$$U_{\text{trans}} = C\,|O|^2 + C\,|R|^2 + C\,OR^* + C\,O^*R$$

### The Exact Reconstruction Case: $C = R$

If the reconstruction wave is identical to the original reference wave ($C = R$), we obtain:

$$U_{\text{trans}} = R\,|O|^2 + R\,|R|^2 + R\,OR^* + R\,O^*R$$

The critical simplification occurs in the third term. Since $RR^* = |R|^2 = r^2$ (a real positive constant):

$$R \cdot OR^* = |R|^2\, O = r^2\, O$$

The full reconstructed field is therefore:

$$\boxed{U_{\text{trans}} = R(|O|^2 + |R|^2) \;+\; r^2\, O \;+\; R^2\, O^*}$$

This field contains three physically distinct waves:

**Term 1: Zero-Order Beam** — $R(|O|^2 + |R|^2)$

This is the reconstruction beam transmitted through the hologram, modulated by the sum of the object and reference intensities. It propagates along the direction of $R$ and carries no useful phase information about the object. It is essentially a slightly attenuated version of the illuminating beam with a weak halo from $|O|^2$.

**Term 2: Virtual Image** — $r^2\, O$

This is the prize. Because $r^2$ is a positive real scalar, this term is an exact copy of the original object wave $O(x,y)$, differing only in brightness. An observer looking through the hologram into this diffracted wave sees a three-dimensional virtual image of the original object, indistinguishable from the original wavefield in both amplitude and phase. The image exhibits full parallax: moving your head reveals different perspectives, exactly as with the real object.

**Term 3: Real (Conjugate) Image** — $R^2\, O^*$

This is the complex conjugate of the object wave, multiplied by $R^2$. The conjugation reverses the sign of the phase, which means the wave converges rather than diverges. It forms a **real image** — one that can be projected onto a screen — on the opposite side of the hologram from the virtual image. However, because the phase is conjugated, the depth is inverted: points that were originally far appear near and vice versa. This **pseudoscopic** image is a characteristic artifact of single-exposure holography.

## 8. The Off-Axis Solution: Separating the Beams

### The Problem with In-Line Holography

In Dennis Gabor's original 1948 configuration, the reference and object waves propagated along the same axis ($\theta = 0$, so $f_x = 0$). While this successfully demonstrated the principle, all three reconstructed beams — zero-order, virtual image, and conjugate image — propagated in the same direction, superimposed on one another. The observer saw the desired image contaminated by the unfocused conjugate and the bright zero-order background.

### The Off-Axis Solution (Leith and Upatnieks, 1962)

The introduction of a nonzero reference angle $\theta$ assigns each term a distinct propagation direction via its spatial frequency content. With $R = r\,e^{i2\pi f_x x}$:

- **Virtual Image Term:** $R\,R^*\,O = r^2\, O$ — propagates along the original object direction (centered at the object's spatial frequencies).

- **Conjugate Image Term:** $R^2\, O^* = r^2\, e^{i4\pi f_x x}\, O^*$ — propagates at an angle corresponding to spatial frequency $2f_x$, i.e., roughly twice the reference angle.

- **Zero-Order Term:** $R(|O|^2 + |R|^2)$ — propagates along the reference direction (centered at $f_x$).

If the reference angle $\theta$ is chosen large enough that $f_x$ exceeds the spatial bandwidth $B$ of the object wave, these three terms occupy non-overlapping regions in the spatial frequency domain. In the far field (or at a lens's Fourier plane), they separate angularly, and the pure object wave $O$ can be observed or isolated without contamination.

The minimum reference angle required for clean separation is:

$$\sin\theta_{\min} \geq 3\lambda B$$

where $B$ is the highest spatial frequency in the object wave. This condition ensures that the sidebands of the three terms do not overlap.

## 9. Summary

Holography converts phase into intensity via interference with a known reference wave, records the resulting fringe pattern, and reconstructs the original wavefield by diffraction. The mathematical structure guarantees exact reconstruction (up to a scalar) when the recording is linear and the reconstruction wave matches the reference. Off-axis geometry separates the desired image from parasitic terms, making practical holography possible.

The key constraints are: mutual coherence of the source, sufficient resolution of the recording medium to capture the carrier frequency, and a linear recording regime to avoid intermodulation artifacts. Within these constraints, the reconstructed wave is a faithful reproduction of the original complex field — amplitude, phase, and all.
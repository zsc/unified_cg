# Chapter 17: Wavefront Shaping and Adaptive Optics

Adaptive optics and wavefront shaping represent the frontier where computational control meets physical optics. In this chapter, we explore how controlled manipulation of optical wavefronts enables focusing through scattering media, aberration correction, and optimization of light transport. We'll establish mathematical foundations using Zernike polynomials, examine spatial light modulators, and connect these concepts to adaptive sampling strategies in computer graphics.

The journey from aberrated wavefronts to diffraction-limited imaging parallels the evolution from naive Monte Carlo sampling to sophisticated adaptive techniques in rendering. Just as adaptive optics corrects atmospheric distortions in real-time, modern rendering algorithms dynamically adjust sampling patterns to minimize variance. This deep connection extends beyond analogy—the mathematical frameworks of wavefront optimization and importance sampling share fundamental principles rooted in signal reconstruction and optimization theory.

## 17.1 Zernike Polynomials and Wavefront Description

### 17.1.1 Mathematical Foundation

Zernike polynomials form a complete orthogonal basis over the unit disk, ideal for describing wavefront aberrations in circular pupils. Named after Frits Zernike, who won the 1953 Nobel Prize for inventing phase contrast microscopy, these polynomials possess unique properties that make them indispensable in optics. The polynomials are defined as:

$$Z_n^m(\rho, \theta) = R_n^{|m|}(\rho) \times \begin{cases}
\cos(m\theta) & m \geq 0 \\
\sin(|m|\theta) & m < 0
\end{cases}$$

where the radial polynomials are:

$$R_n^m(\rho) = \sum_{k=0}^{(n-m)/2} \frac{(-1)^k (n-k)!}{k! \left(\frac{n+m}{2}-k\right)! \left(\frac{n-m}{2}-k\right)!} \rho^{n-2k}$$

The indices follow constraints:
- n ≥ 0 (radial order)
- -n ≤ m ≤ n (azimuthal frequency)
- n - |m| is even

This last constraint ensures that R_n^m(ρ) contains only even or odd powers of ρ, maintaining the polynomial's parity properties. The radial polynomial can also be expressed using Jacobi polynomials:

$$R_n^m(\rho) = (-1)^{(n-m)/2} \rho^m P_{(n-m)/2}^{(m,0)}(2\rho^2 - 1)$$

where P_k^{(α,β)} are Jacobi polynomials. This connection reveals deeper mathematical structure and enables efficient computation using recurrence relations.

For normalization, we include the factor:
$$N_n^m = \sqrt{\frac{2(n+1)}{1 + \delta_{m0}}}$$

making the normalized polynomials:
$$\tilde{Z}_n^m(\rho, \theta) = N_n^m Z_n^m(\rho, \theta)$$

The normalization ensures unit variance for each mode when integrated over the unit disk, crucial for comparing aberration strengths across different orders.

### 17.1.2 Orthogonality Properties

The orthogonality relation over the unit disk forms the mathematical backbone of wavefront analysis:

$$\int_0^{2\pi} \int_0^1 Z_n^m(\rho, \theta) Z_{n'}^{m'}(\rho, \theta) \rho d\rho d\theta = \frac{\pi}{2n+2} \delta_{nn'} \delta_{mm'}$$

This orthogonality arises from two separate integrations that decouple due to the polynomial structure:

**Angular orthogonality:**
$$\int_0^{2\pi} \cos(m\theta)\cos(m'\theta) d\theta = \pi\delta_{mm'} \quad (m, m' > 0)$$
$$\int_0^{2\pi} \sin(m\theta)\sin(m'\theta) d\theta = \pi\delta_{mm'} \quad (m, m' > 0)$$
$$\int_0^{2\pi} \cos(m\theta)\sin(m'\theta) d\theta = 0$$

These relations follow from the orthogonality of trigonometric functions over a complete period. The cross terms vanish due to the odd symmetry of the integrand.

**Radial orthogonality:**
$$\int_0^1 R_n^m(\rho) R_{n'}^m(\rho) \rho d\rho = \frac{1}{2(n+1)} \delta_{nn'}$$

The weight function ρ in the integral arises from the Jacobian in polar coordinates and ensures proper orthogonality over the circular domain. This weight is crucial—without it, the radial polynomials would not be orthogonal. The specific form of R_n^m is carefully constructed to achieve orthogonality with respect to this measure.

**Completeness and Parseval's Identity:**
Any square-integrable function over the unit disk can be expanded in Zernike polynomials, with Parseval's identity guaranteeing energy conservation:

$$\int_0^{2\pi} \int_0^1 |W(\rho, \theta)|^2 \rho d\rho d\theta = \sum_{n=0}^{\infty} \sum_{m=-n}^n |a_n^m|^2 \frac{\pi}{2n+2}$$

This relationship enables direct computation of wavefront variance from Zernike coefficients, fundamental for aberration budgeting in optical design.

### 17.1.3 Wavefront Expansion

Any wavefront W(ρ, θ) can be expanded as:

$$W(\rho, \theta) = \sum_{n=0}^{\infty} \sum_{m=-n}^n a_n^m Z_n^m(\rho, \theta)$$

The coefficients are computed via the inner product:

$$a_n^m = \frac{2n+2}{\pi \epsilon_m} \int_0^{2\pi} \int_0^1 W(\rho, \theta) Z_n^m(\rho, \theta) \rho d\rho d\theta$$

where εₘ = 2 for m = 0 and εₘ = 1 otherwise. This factor accounts for the different normalization of the m = 0 (rotationally symmetric) terms.

**Efficient Coefficient Computation:**
For measured wavefront data on a discrete grid, the coefficients can be computed using least squares:

$$\mathbf{a} = (\mathbf{Z}^T \mathbf{Z})^{-1} \mathbf{Z}^T \mathbf{w}$$

where **Z** is the matrix of Zernike polynomial values at measurement points, and **w** is the vector of wavefront measurements. For a regular grid with proper sampling, Z^T Z becomes nearly diagonal, simplifying inversion.

In practice, we truncate the expansion at some maximum order N:
$$W(\rho, \theta) \approx \sum_{n=0}^{N} \sum_{m=-n}^n a_n^m Z_n^m(\rho, \theta)$$

The number of terms up to order N is:
$$J = \frac{(N+1)(N+2)}{2}$$

For example:
- N = 3: J = 10 terms (up to coma)
- N = 4: J = 15 terms (includes spherical aberration)
- N = 8: J = 45 terms (high-order aberrations)

**Truncation Error Analysis:**
The RMS error from truncating at order N follows:

$$\sigma_{truncation}^2 = \sum_{n=N+1}^{\infty} \sum_{m=-n}^n |a_n^m|^2$$

For atmospheric turbulence, this error decreases as N^(-α/2) where α ≈ 11/3, providing guidance for selecting truncation order based on desired accuracy.

### 17.1.4 Conversion Between Indices

Several indexing schemes exist for Zernike polynomials, each with specific advantages for different applications:

**Noll notation (single index j):**
$$j = \frac{n(n+1)}{2} + |m| + \begin{cases}
0 & \text{if } m \leq 0 \text{ and } n \bmod 4 \in \{0,1\} \\
0 & \text{if } m > 0 \text{ and } n \bmod 4 \in \{2,3\} \\
1 & \text{otherwise}
\end{cases}$$

This ordering alternates between sine and cosine terms to maintain a logical progression of aberration types. The first few Noll indices correspond to:
- j = 1: Piston (Z₀⁰)
- j = 2,3: Tip/Tilt (Z₁¹, Z₁⁻¹)
- j = 4,5,6: Astigmatism and Defocus (Z₂⁻², Z₂⁰, Z₂²)

**OSA/ANSI standard:**
$$j = \frac{n(n+2) + m}{2}$$

This simpler formula creates a monotonic ordering that's easier to compute but less intuitive for aberration analysis.

**Wyant ordering:**
Groups terms by polynomial order, then by azimuthal frequency. This ordering facilitates:
- Systematic increase in polynomial complexity
- Clear separation of radial orders
- Simplified error propagation analysis

**Fringe/University of Arizona notation:**
Orders by increasing spatial frequency, useful for interferometric analysis where higher frequencies correspond to finer fringe patterns.

The choice of ordering affects coefficient interpretation but not the mathematical properties. Conversion matrices between orderings are sparse and can be precomputed for efficiency.

### 17.1.5 Common Aberrations

The first few Zernike terms correspond to familiar optical aberrations, each with distinct physical origins and visual effects:

**Low-order terms (n ≤ 2):**
- $Z_0^0 = 1$: **Piston** (constant phase offset)
  - No effect on image quality, only absolute phase
  - Important for interferometry and coherent imaging
  
- $Z_1^{-1} = 2\rho\sin\theta$: **Vertical tilt** (y-tilt)
- $Z_1^1 = 2\rho\cos\theta$: **Horizontal tilt** (x-tilt)
  - Cause image displacement without blur
  - Arise from wedge in optical elements or misalignment
  
- $Z_2^{-2} = \sqrt{6}\rho^2\sin(2\theta)$: **Oblique astigmatism** (45°)
- $Z_2^2 = \sqrt{6}\rho^2\cos(2\theta)$: **Vertical astigmatism** (0°/90°)
  - Create orthogonal line foci at different distances
  - Common in cylindrical optical elements and eye aberrations
  
- $Z_2^0 = \sqrt{3}(2\rho^2 - 1)$: **Defocus**
  - Symmetric blur, shifts best focus position
  - Equivalent to longitudinal shift of image plane

**Third-order terms (n = 3):**
- $Z_3^{-3} = \sqrt{8}\rho^3\sin(3\theta)$: **Trefoil**
- $Z_3^3 = \sqrt{8}\rho^3\cos(3\theta)$: **Trefoil**
  - Three-fold symmetric distortion
  - Often from stressed optical mounts
  
- $Z_3^{-1} = \sqrt{8}(3\rho^3 - 2\rho)\sin\theta$: **Vertical coma**
- $Z_3^1 = \sqrt{8}(3\rho^3 - 2\rho)\cos\theta$: **Horizontal coma**
  - Comet-like blur increasing with field angle
  - Fundamental limit in off-axis imaging

**Fourth-order terms (n = 4):**
- $Z_4^{-4} = \sqrt{10}\rho^4\sin(4\theta)$: **Tetrafoil**
- $Z_4^4 = \sqrt{10}\rho^4\cos(4\theta)$: **Tetrafoil**
  - Four-fold symmetric, square-like distortion
  
- $Z_4^{-2} = \sqrt{10}(4\rho^4 - 3\rho^2)\sin(2\theta)$: **Secondary astigmatism**
- $Z_4^2 = \sqrt{10}(4\rho^4 - 3\rho^2)\cos(2\theta)$: **Secondary astigmatism**
  - Higher-order astigmatic effects
  
- $Z_4^0 = \sqrt{5}(6\rho^4 - 6\rho^2 + 1)$: **Primary spherical aberration**
  - Radially symmetric blur increasing with aperture
  - Fundamental limit of spherical surfaces

**Visual Impact:**
Each aberration creates characteristic point spread function (PSF) distortions:
- Coma: Comet-shaped PSF with tail pointing radially
- Astigmatism: Elliptical PSF rotating with focus position
- Spherical: Circular halo with bright core
- Trefoil/Tetrafoil: Multi-lobed PSF patterns

### 17.1.6 Physical Interpretation and Seidel Connection

The Zernike polynomials relate to classical Seidel aberrations through coordinate transformation. For a point at field angle α and pupil coordinates (ρ, θ):

**Seidel to Zernike mapping:**
- Spherical aberration: $W_{040} = a_4^0 Z_4^0$
- Coma: $W_{131} = a_3^1 Z_3^1 \cos\alpha + a_3^{-1} Z_3^{-1} \sin\alpha$
- Astigmatism: $W_{222} = a_2^2 Z_2^2 \cos(2\alpha) + a_2^{-2} Z_2^{-2} \sin(2\alpha)$
- Field curvature: $W_{220} = a_2^0 Z_2^0$
- Distortion: $W_{311} = a_1^1 Z_1^1 \cos\alpha + a_1^{-1} Z_1^{-1} \sin\alpha$

The classical Seidel aberration coefficients emerge from ray tracing, while Zernike coefficients come from wavefront fitting. The transformation between them involves:

$$W_{Seidel}(h, \rho, \theta, \phi) = \sum_{i,j,k,l} S_{ijkl} h^i \rho^j \cos^k(\theta-\phi)$$

where h is the normalized field height, and the indices satisfy i + j = 2(k + l). This power series expansion can be rewritten in Zernike form through trigonometric identities.

**Field Dependence:**
Unlike pupil-only Zernike polynomials, full field aberrations require additional field coordinates:
$$W(h_x, h_y, \rho, \theta) = \sum_{p,q,n,m} b_{pqnm} h_x^p h_y^q Z_n^m(\rho, \theta)$$

This field-dependent Zernike expansion enables:
- Nodal aberration theory for wide-field systems
- Efficient storage of measured aberration fields
- Direct optimization of field-averaged performance

The wavefront variance contribution from each aberration:
$$\sigma_n^2 = \sum_{m=-n}^n (a_n^m)^2$$

This decomposition reveals which orders dominate the aberration budget, guiding correction strategies.

### 17.1.7 RMS Wavefront Error

The root-mean-square wavefront error is elegantly expressed in the Zernike basis:

$$\sigma_{RMS} = \sqrt{\sum_{n,m} (a_n^m)^2}$$

This orthogonality property makes Zernike coefficients ideal for optimization algorithms. The variance can be decomposed by order:

$$\sigma_{RMS}^2 = \sum_{n=0}^{\infty} \sigma_n^2 = \sum_{n=0}^{\infty} \sum_{m=-n}^n (a_n^m)^2$$

**Physical Significance:**
The RMS wavefront error directly relates to optical performance metrics:
- Strehl ratio (peak intensity): $S \approx \exp(-\sigma_{RMS}^2)$ for $\sigma_{RMS} < 1$ radian
- Encircled energy: Fraction of light within given radius decreases with σ_RMS
- Resolution: Effective PSF width increases approximately as $(1 + \sigma_{RMS}^2)^{1/2}$

For atmospheric turbulence following Kolmogorov statistics:
$$\langle (a_n^m)^2 \rangle \propto (n + 1)^{-\alpha}$$

where α ≈ 11/3 for fully developed turbulence. This power-law decay justifies truncating the expansion at finite order.

**Fried Parameter Connection:**
The seeing-limited resolution is characterized by the Fried parameter r₀:
$$\sigma_{RMS}^2 = 1.03(D/r_0)^{5/3}$$

where D is the telescope diameter. For D >> r₀, the wavefront error is dominated by low-order modes, particularly tip-tilt which contains ~87% of the total variance.

**Temporal Evolution:**
Atmospheric aberrations evolve with characteristic timescales:
$$\tau_n \approx \frac{r_0}{v_{wind}} \cdot n^{-3/5}$$

Higher-order aberrations change more rapidly, requiring faster correction loops. This temporal spectrum guides adaptive optics control system design.

### 17.1.8 Strehl Ratio and Performance Metrics

For small aberrations, the Strehl ratio (peak intensity ratio) approximates to:

$$S \approx \exp(-\sigma_{RMS}^2) \approx 1 - \sigma_{RMS}^2$$

where σ_RMS is in radians. This provides a direct connection between wavefront quality and imaging performance.

**Extended Maréchal approximation:**
$$S \approx \exp\left[-\sigma_{RMS}^2 + \frac{1}{2}\sigma_{RMS}^4(\kappa_4 - 1)\right]$$

where κ₄ is the normalized fourth moment of the phase distribution. For Gaussian statistics, κ₄ = 3. Non-Gaussian phase screens (e.g., from strong scintillation) require this correction.

**Higher-Order Statistics:**
Beyond RMS, the phase structure function provides additional insight:
$$D_\phi(\mathbf{r}) = \langle[W(\mathbf{x} + \mathbf{r}) - W(\mathbf{x})]^2\rangle$$

For Kolmogorov turbulence:
$$D_\phi(r) = 6.88(r/r_0)^{5/3}$$

This structure function determines the correlation length of aberrations and guides actuator spacing in deformable mirrors.

**Diffraction-limited criterion:**
- Rayleigh criterion: λ/4 peak-to-valley → S ≈ 0.80
- Maréchal criterion: λ/14 RMS → S ≈ 0.80
- Practical limit: σ_RMS < λ/20 → S > 0.94

The point spread function (PSF) with aberrations:
$$\text{PSF}(x,y) = \left|\mathcal{F}\left\{P(\xi,\eta)\exp\left[i\frac{2\pi}{\lambda}W(\xi,\eta)\right]\right\}\right|^2$$

where P is the pupil function and W is the wavefront error.

**Partial Correction Effects:**
When only correcting up to Zernike order N:
$$S_{corrected} = S_{uncorrected} \cdot \exp\left(\sum_{n=0}^N \sum_{m=-n}^n |a_n^m|^2\right)$$

This shows exponential improvement with each corrected mode, motivating hierarchical correction schemes.

### 17.1.9 Adaptive Optics Performance

The residual wavefront error after correction decomposes into independent error sources:
$$\sigma_{residual}^2 = \sigma_{fitting}^2 + \sigma_{temporal}^2 + \sigma_{measurement}^2 + \sigma_{calibration}^2 + \sigma_{anisoplanatic}^2$$

Each term represents a fundamental limitation in adaptive optics systems:

**Fitting error (finite actuators):**
$$\sigma_{fitting}^2 \approx 0.28\left(\frac{d}{r_0}\right)^{5/3}$$

where d is actuator spacing and r₀ is the Fried parameter. This error arises from the deformable mirror's inability to reproduce high spatial frequencies. The coefficient 0.28 assumes Kolmogorov turbulence and continuous-facesheet mirrors.

**Temporal error (finite bandwidth):**
$$\sigma_{temporal}^2 \approx \left(\frac{f_G}{f_0}\right)^{5/3}$$

where the Greenwood frequency is:
$$f_G = 0.43 \frac{v_{wind}}{r_0}$$

and f₀ is the control loop bandwidth. This error accumulates when atmospheric changes outpace the correction system.

**Measurement noise propagation:**
$$\sigma_{measurement}^2 = \left(\frac{\partial W}{\partial s}\right)^2 \sigma_s^2$$

where s represents sensor measurements. For Shack-Hartmann sensors:
$$\sigma_{measurement}^2 \propto \frac{1}{N_{photons}} + \frac{\sigma_{read}^2}{N_{photons}^2}$$

The first term is photon noise, the second is read noise contribution.

**Anisoplanatic error:**
$$\sigma_{anisoplanatic}^2 = \left(\frac{\theta}{\theta_0}\right)^{5/3}$$

where θ is the angular separation from the guide star and θ₀ is the isoplanatic angle:
$$\theta_0 = 0.314 \frac{r_0}{H_{eff}}$$

with H_{eff} being the effective turbulence height.

**Error Budget Optimization:**
The total error minimization requires balancing:
- More actuators reduce fitting error but increase cost
- Higher bandwidth reduces temporal error but amplifies noise
- Brighter guide stars reduce measurement error but limit sky coverage

This optimization problem parallels variance reduction in Monte Carlo rendering, where sample allocation must balance different error sources.

## 17.2 Spatial Light Modulator (SLM) Principles

### 17.2.1 Phase Modulation Mechanisms

Spatial light modulators enable pixel-wise control of optical phase, amplitude, or polarization. Several technologies exist, each with distinct characteristics:

**Liquid Crystal on Silicon (LCoS):**
The most common SLM technology uses electrically controlled birefringence. For nematic liquid crystals, the phase modulation depth δ depends on:

$$\delta = \frac{2\pi}{\lambda} \Delta n \cdot d$$

where Δn is the birefringence change and d is the cell thickness. The voltage-dependent refractive index follows:

$$n(V) = n_o + \frac{n_e - n_o}{1 + (V_{th}/V)^2}$$

where n_o and n_e are ordinary and extraordinary refractive indices, and V_{th} is the threshold voltage.

**Director Orientation Model:**
The liquid crystal director angle θ(z) through the cell thickness satisfies:
$$K\frac{d^2\theta}{dz^2} = \epsilon_0\Delta\epsilon E^2\sin\theta\cos\theta$$

where K is the elastic constant. The resulting phase shift:
$$\phi = \frac{2\pi}{\lambda}\int_0^d n_{eff}(\theta(z))dz$$

with effective index:
$$n_{eff}(\theta) = \frac{n_o n_e}{\sqrt{n_o^2\sin^2\theta + n_e^2\cos^2\theta}}$$

**Response Dynamics:**
The response time scales as:
$$\tau_{rise} = \frac{\gamma d^2}{K\epsilon_0\Delta\epsilon(V^2 - V_{th}^2)}$$
$$\tau_{fall} = \frac{\gamma d^2}{K\pi^2}$$

where γ is the rotational viscosity. Note the asymmetry: relaxation is typically slower than activation.

**Temperature Dependence:**
Both birefringence and response time vary with temperature:
$$\Delta n(T) = \Delta n_0(1 - T/T_c)^\beta$$
$$\gamma(T) = \gamma_0\exp(E_a/k_BT)$$

where T_c is the clearing temperature, β ≈ 0.2, and E_a is the activation energy.

**Digital Micromirror Devices (DMD):**
Binary amplitude modulation via tilting mirrors:
- Tilt angles: ±12° (typically)
- Switching time: ~10 μs
- Diffraction efficiency: ~88% (into desired order)
- Contrast ratio: >5000:1

**Deformable Mirrors:**
Continuous surface deformation for phase control:
- Actuator types: Piezoelectric, electrostatic, magnetic
- Stroke: 1-10 μm typical
- Bandwidth: 1-10 kHz
- Inter-actuator coupling via influence functions

### 17.2.2 Complex Amplitude Modulation

Phase-only SLMs can achieve complex amplitude modulation through several encoding schemes:

**1. Off-axis Computer Generated Hologram (Lee Method):**
Encode complex field U = A exp(iψ) as:
$$\phi(x,y) = \arg[A(x,y)e^{i\psi(x,y)}] + 2\pi f_c x$$

where f_c is a carrier frequency. The first diffraction order approximates the desired complex field.

The carrier frequency must satisfy:
$$f_c > f_{max} + \frac{W}{2\lambda z}$$

where f_max is the maximum spatial frequency in U and W is the beam width. This ensures separation of diffraction orders.

**Efficiency analysis:**
- First-order efficiency: η₁ ≈ |⟨U⟩|²/⟨|U|²⟩
- Zero-order leakage: η₀ = |⟨exp(iφ)⟩|²
- Signal-to-noise ratio: SNR ∝ A²/(1-A²)

**2. Double-Phase Amplitude Encoding:**
Using two phase masks φ₁ and φ₂ separated by distance z:

$$U_{out} = \mathcal{F}^{-1}\{\mathcal{F}\{e^{i\phi_1}\} \cdot H(z) \cdot \mathcal{F}\{e^{i\phi_2}\}\}$$

where H(z) is the Fresnel propagation kernel:
$$H(k_x, k_y; z) = \exp\left[iz\sqrt{k^2 - k_x^2 - k_y^2}\right]$$

Phase masks derived from:
$$\phi_1 = \arg[U] + \arg[\mathcal{F}^{-1}\{|U|^{1/2}\}]$$
$$\phi_2 = -\arg[\mathcal{F}^{-1}\{|U|^{1/2}\}]$$

**3. Iterative Fourier Transform Algorithms:**

**Gerchberg-Saxton (GS):**
$$\phi_{n+1} = \arg[\mathcal{F}^{-1}\{|A_{target}|e^{i\arg[\mathcal{F}\{|A_{input}|e^{i\phi_n}\}]}\}]$$

**Weighted GS with feedback:**
$$\phi_{n+1} = \phi_n + \alpha \arg[\mathcal{F}^{-1}\{(A_{target} - A_n)e^{i\psi_n}\}]$$

where α ∈ (0,2) is the feedback strength.

**Error reduction metrics:**
$$\epsilon = \frac{\sum|I_{target} - I_{measured}|^2}{\sum I_{target}^2}$$

Typical convergence: ε < 0.01 within 20-50 iterations.

**4. Superpixel Method:**
Group N×N pixels to encode both amplitude and phase:
- Central pixels: Encode phase
- Border pixels: Control amplitude via blazed gratings

Amplitude control via diffraction efficiency:
$$A_{effective} = \eta(\theta_{blaze}) = \text{sinc}^2\left(\frac{N\pi\sin\theta_{blaze}}{\lambda/p}\right)$$

where p is the pixel pitch.

### 17.2.3 Diffraction Efficiency and Limitations

The diffraction efficiency η for a phase grating:

$$\eta_m = \left|\frac{\sin(m\pi\Delta\phi/2\pi)}{m\pi\Delta\phi/2\pi}\right|^2$$

For binary phase masks (0, π):
- η₁ = 4/π² ≈ 40.5% (first order)
- η₀ = 0 (zero order suppressed)

Pixelation effects introduce a sinc envelope:

$$E_{pixel}(u,v) = \text{sinc}(au)\text{sinc}(bv)$$

where a, b are pixel dimensions. Fill factor F reduces efficiency:

$$\eta_{effective} = F^2 \cdot \eta_{ideal}$$

## 17.3 Phase Conjugation and Time Reversal

### 17.3.1 Optical Phase Conjugation Theory

Phase conjugation creates a wave that propagates backward along the original path, undoing distortions. For a forward wave:

$$E_{forward}(\mathbf{r}, t) = A(\mathbf{r})e^{i[\mathbf{k}\cdot\mathbf{r} - \omega t + \phi(\mathbf{r})]}$$

The phase-conjugate wave is:

$$E_{conjugate}(\mathbf{r}, t) = A^*(\mathbf{r})e^{i[-\mathbf{k}\cdot\mathbf{r} - \omega t - \phi(\mathbf{r})]}$$

In four-wave mixing, the conjugate field emerges from:

$$E_4 = \chi^{(3)} E_1 E_2 E_3^*$$

where χ^(3) is the third-order nonlinear susceptibility. The phase-matching condition:

$$\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_2 - \mathbf{k}_3$$

ensures momentum conservation.

### 17.3.2 Digital Phase Conjugation

Digital implementation via holographic recording and playback:

**Recording Phase:**
$$H(\mathbf{r}) = |E_{signal} + E_{reference}|^2 = |E_s|^2 + |E_r|^2 + E_s E_r^* + E_s^* E_r$$

**Playback with Conjugate Reference:**
$$E_{out} = H \cdot E_r^* = E_s|E_r|^2 + \text{other terms}$$

The first term represents the phase-conjugate wave. For SLM implementation:

$$\phi_{SLM}(x,y) = -\phi_{measured}(x,y) + \phi_{carrier}$$

### 17.3.3 Time Reversal Symmetry

The wave equation's time-reversal invariance enables focusing through disorder:

$$\nabla^2 E - \frac{n^2(\mathbf{r})}{c^2}\frac{\partial^2 E}{\partial t^2} = 0$$

Substituting t → -t leaves the equation unchanged. For monochromatic fields:

$$[\nabla^2 + k^2 n^2(\mathbf{r})]E = 0$$

If E(r) is a solution, then E*(r) is also a solution propagating in reverse. This principle underlies:

1. **Reciprocity:** T_{ij} = T_{ji} for transmission matrix
2. **Focusing invariance:** If light focuses from A→B through disorder, conjugation focuses B→A
3. **Scattering cancellation:** Multiple scattering paths interfere constructively at the original source

## 17.4 Focusing Through Scattering Media

### 17.4.1 Transmission Matrix Formalism

The transmission matrix T relates input and output optical modes through a scattering medium:

$$\mathbf{E}_{out} = \mathbf{T} \cdot \mathbf{E}_{in}$$

For N input modes and M output modes, T is an M×N complex matrix. Each element:

$$T_{mn} = |T_{mn}|e^{i\phi_{mn}}$$

represents amplitude and phase coupling between modes. The singular value decomposition:

$$\mathbf{T} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^†$$

reveals optimal transmission channels. The transmission eigenvalues τᵢ = σᵢ² follow the Marchenko-Pastur distribution:

$$\rho(\tau) = \frac{1}{2\pi\gamma\tau}\sqrt{\frac{(\tau_+ - \tau)(\tau - \tau_-)}{1 + \tau}}$$

where γ = M/N and τ± = (1 ± √γ)².

### 17.4.2 Wavefront Optimization Algorithms

**1. Iterative Optimization:**
Maximize intensity at target position by phase control:

$$I_{target} = \left|\sum_{n=1}^N T_{mn} A_n e^{i\phi_n}\right|^2$$

Gradient ascent update:
$$\phi_n^{(k+1)} = \phi_n^{(k)} + \alpha \frac{\partial I_{target}}{\partial \phi_n}$$

where:
$$\frac{\partial I_{target}}{\partial \phi_n} = 2\text{Im}[E_{target}^* T_{mn} A_n e^{i\phi_n}]$$

**2. Genetic Algorithm:**
Population of phase patterns evolves via:
- Selection: Roulette wheel based on fitness I_{target}
- Crossover: Uniform or arithmetic mixing of phase patterns
- Mutation: Gaussian perturbation σ ~ π/10

**3. Hadamard Basis Measurement:**
Systematic T measurement using orthogonal patterns:

$$\phi_{Hadamard}^{(k)} = H_{kn} \cdot \pi$$

where H is the Hadamard matrix. Requires 4N measurements for full complex T.

### 17.4.3 Memory Effects and Correlations

**Angular Memory Effect:**
The speckle pattern translates with input angle change:

$$C(\Delta\theta) = \langle I(\theta)I(\theta + \Delta\theta) \rangle / \langle I \rangle^2$$

The correlation width:
$$\Delta\theta_{memory} \approx \frac{\lambda}{2\pi L}$$

where L is the medium thickness. Within this range:

$$T_{m,n+\Delta n} \approx T_{mn} e^{i\phi_{shift}(\Delta n)}$$

**Spectral Memory Effect:**
Frequency correlation function:

$$C(\Delta\omega) = \exp\left[-\left(\frac{\Delta\omega}{\Delta\omega_c}\right)^2\right]$$

with correlation bandwidth:
$$\Delta\omega_c = \frac{2\pi c}{L_{path}\sqrt{\ln 2}}$$

where L_{path} is the average path length through the medium.

**Implications for Imaging:**
1. Single optimization enables scanning within memory effect range
2. Broadband focusing possible within spectral correlation
3. Time-gated detection isolates ballistic photons: τ_{gate} ~ L/c

## 17.5 Adaptive Sampling and Optimization in Graphics

### 17.5.1 Connection to Importance Sampling

Wavefront optimization parallels importance sampling in Monte Carlo rendering. The rendering equation:

$$L_o(\mathbf{x}, \omega_o) = \int_{\Omega} f_r(\mathbf{x}, \omega_i, \omega_o) L_i(\mathbf{x}, \omega_i) |\cos\theta_i| d\omega_i$$

Importance sampling chooses directions ωᵢ proportional to the integrand. Similarly, wavefront shaping optimizes:

$$I_{focus} = \left|\int_A \psi_{in}(\mathbf{r}) G(\mathbf{r}, \mathbf{r}_{target}) dA\right|^2$$

where G is the Green's function through the medium. The optimal input field:

$$\psi_{in}^{opt}(\mathbf{r}) \propto G^*(\mathbf{r}, \mathbf{r}_{target})$$

This is the phase-conjugate of the transmission kernel—analogous to sampling proportional to BSDF × incoming radiance.

### 17.5.2 Adaptive Path Guiding

**Spatial-Directional Trees:**
Partition (x, ω) space adaptively based on radiance variation:

$$\text{Split criterion} = \frac{\text{Var}[L(\mathbf{x}, \omega)]}{\text{Mean}[L(\mathbf{x}, \omega)]^2}$$

The guided sampling PDF:

$$p_{guide}(\omega) = \alpha p_{learned}(\omega) + (1-\alpha) p_{BSDF}(\omega)$$

where α adapts based on learning confidence.

**Gaussian Mixture Models:**
Represent incident radiance field as:

$$L_i(\mathbf{x}, \omega) \approx \sum_{k=1}^K w_k \mathcal{N}(\omega; \mu_k, \Sigma_k)$$

Update via EM algorithm:
- E-step: Compute responsibilities γₖ
- M-step: Update μₖ, Σₖ, wₖ

**Path Space Markov Chains:**
Mutate paths while maintaining detailed balance:

$$\frac{p(\mathbf{x} \to \mathbf{y})}{p(\mathbf{y} \to \mathbf{x})} = \frac{f(\mathbf{y})}{f(\mathbf{x})}$$

Perturbations guided by learned radiance gradients:

$$\Delta\mathbf{x} = \epsilon \nabla_{\mathbf{x}} \log L(\mathbf{x})$$

### 17.5.3 Learned Light Transport

**Neural Radiance Caching:**
Train network to predict incident radiance:

$$L_i^{predicted}(\mathbf{x}, \omega) = \mathcal{N}_\theta(\mathbf{x}, \omega, \text{features})$$

Loss function combines L2 and relative error:

$$\mathcal{L} = \|L_i - L_i^{predicted}\|^2 + \lambda \frac{\|L_i - L_i^{predicted}\|^2}{\|L_i\|^2 + \epsilon}$$

**Learned Transmittance Functions:**
For participating media, predict optical depth:

$$\tau_{predicted}(t) = \int_0^t \sigma_t^{predicted}(s) ds$$

Network architecture exploits symmetries:
- Translation invariance → Convolutional layers
- Rotation equivariance → Spherical harmonics features
- Scale invariance → Multi-resolution encoding

**Differentiable Rendering Loop:**
Connect to wavefront optimization via gradient flow:

$$\frac{\partial L_{pixel}}{\partial \phi_{SLM}} = \sum_{paths} \frac{\partial L_{pixel}}{\partial \mathbf{x}_i} \frac{\partial \mathbf{x}_i}{\partial \phi_{SLM}}$$

This enables joint optimization of:
1. Sampling distributions
2. Denoising networks
3. Light transport approximations

## Chapter Summary

This chapter explored the mathematical foundations and practical applications of wavefront shaping and adaptive optics, bridging physical optics with computational graphics:

**Key Concepts:**
- **Zernike polynomials** provide an orthogonal basis for wavefront description, enabling systematic aberration analysis and correction
- **Spatial light modulators** implement programmable phase control through liquid crystal technology, with efficiency limited by pixelation and fill factor
- **Phase conjugation** exploits time-reversal symmetry to undo propagation distortions, enabling focusing through scattering media
- **Transmission matrix** formalism describes light propagation through complex media as a linear transformation between modes
- **Memory effects** allow single-shot optimization to work over finite angular and spectral ranges
- **Adaptive sampling** in graphics shares mathematical principles with wavefront optimization in optics

**Key Equations:**
- Zernike expansion: $W(\rho, \theta) = \sum_{n,m} a_n^m Z_n^m(\rho, \theta)$
- Phase conjugate wave: $E_{conjugate} = A^* e^{i[-\mathbf{k}\cdot\mathbf{r} - \omega t - \phi]}$
- Transmission matrix: $\mathbf{E}_{out} = \mathbf{T} \cdot \mathbf{E}_{in}$
- Memory effect range: $\Delta\theta_{memory} \approx \lambda/(2\pi L)$
- Optimal importance sampling: $\psi_{in}^{opt} \propto G^*(\mathbf{r}, \mathbf{r}_{target})$

## Exercises

### Basic Exercises (Building Intuition)

**Exercise 17.1:** Zernike Polynomial Orthogonality
Prove that $Z_2^0 = \sqrt{3}(2\rho^2 - 1)$ (defocus) and $Z_2^2 = \sqrt{6}\rho^2\cos(2\theta)$ (astigmatism) are orthogonal over the unit disk.
<details>
<summary>Hint</summary>
Use the orthogonality relation and remember that ∫₀²π cos(2θ)dθ = 0.
</details>

<details>
<summary>Answer</summary>

The orthogonality integral:
$$\int_0^{2\pi}\int_0^1 Z_2^0 Z_2^2 \rho d\rho d\theta = \sqrt{18}\int_0^{2\pi}\cos(2\theta)d\theta \int_0^1 (2\rho^2-1)\rho^3 d\rho$$

The angular integral: $\int_0^{2\pi}\cos(2\theta)d\theta = 0$

Therefore, the product vanishes regardless of the radial integral, proving orthogonality.
</details>

**Exercise 17.2:** SLM Phase Wrapping
An SLM provides 0-2π phase modulation. To achieve 4π total phase shift across the aperture, what spatial frequency grating must be added? What is the diffraction efficiency into the first order?

<details>
<summary>Hint</summary>
Consider a blazed grating pattern and use the diffraction efficiency formula for phase gratings.
</details>

<details>
<summary>Answer</summary>

Add a linear phase ramp: $\phi(x) = 4\pi x/D \bmod 2\pi$, creating a binary grating with period D/2.

For a binary (0,π) grating, the first-order efficiency:
$$\eta_1 = \left|\frac{\sin(\pi/2)}{\pi/2}\right|^2 = \frac{4}{\pi^2} \approx 0.405$$

40.5% of light diffracts into the desired order.
</details>

**Exercise 17.3:** Memory Effect Range
A 1mm thick scattering medium with transport mean free path ℓ* = 50μm is illuminated at λ = 632.8nm. Calculate the angular memory effect range and the number of independent speckle patterns across a 10° field of view.

<details>
<summary>Hint</summary>
Use Δθ ≈ λ/(2πL) and consider how many non-overlapping ranges fit in the total field.
</details>

<details>
<summary>Answer</summary>

Memory effect range:
$$\Delta\theta = \frac{632.8 \times 10^{-9}}{2\pi \times 10^{-3}} = 1.01 \times 10^{-4} \text{ rad} = 0.0058°$$

Number of independent patterns:
$$N = \left(\frac{10°}{0.0058°}\right)^2 \approx 2.98 \times 10^6$$

Nearly 3 million independent optimization regions exist across the field.
</details>

### Challenging Exercises (Extending Concepts)

**Exercise 17.4:** Transmission Matrix Eigenvalue Bounds
For a transmission matrix T with M output modes and N input modes (M > N), prove that the maximum transmission eigenvalue satisfies τ_max ≤ 1 + √(M/N).

<details>
<summary>Hint</summary>
Use the Marchenko-Pastur distribution and consider energy conservation constraints.
</details>

<details>
<summary>Answer</summary>

From the Marchenko-Pastur distribution, the support is [τ₋, τ₊] where:
$$\tau_\pm = (1 \pm \sqrt{M/N})^2$$

Energy conservation requires Tr(T†T) = N (average transmission = 1). The maximum eigenvalue occurs at the upper edge:
$$\tau_{max} = \tau_+ = (1 + \sqrt{M/N})^2$$

Taking the square root of the transmission eigenvalue gives the singular value:
$$\sigma_{max} = 1 + \sqrt{M/N}$$

This bound shows that oversampling (M > N) enables focusing enhancement proportional to √(M/N).
</details>

**Exercise 17.5:** Phase Conjugation Fidelity
A phase conjugation system measures the field E_measured = E_true + E_noise where E_noise is Gaussian noise. Derive the expected fidelity F = |⟨E_true|E_conjugate⟩|²/|E_true|² as a function of SNR.

<details>
<summary>Hint</summary>
Express the conjugated field in terms of measured field and compute the overlap integral.
</details>

<details>
<summary>Answer</summary>

The conjugated field: $E_{conjugate} = E_{measured}^* = E_{true}^* + E_{noise}^*$

Fidelity calculation:
$$F = \frac{|⟨E_{true}|E_{true}^* + E_{noise}^*⟩|^2}{|E_{true}|^2}$$

$$F = \frac{||E_{true}|^2 + ⟨E_{true}|E_{noise}^*⟩|^2}{|E_{true}|^2}$$

For uncorrelated Gaussian noise with variance σ²:
$$F \approx \frac{|E_{true}|^4}{|E_{true}|^2(|E_{true}|^2 + \sigma^2)} = \frac{1}{1 + \sigma^2/|E_{true}|^2}$$

$$F = \frac{\text{SNR}}{1 + \text{SNR}}$$

High fidelity requires SNR >> 1.
</details>

**Exercise 17.6:** Adaptive Sampling Convergence
In path guiding, samples are drawn from p_guide(ω) = αp_learned(ω) + (1-α)p_BSDF(ω). Derive the optimal α that minimizes variance for a given learning quality metric Q ∈ [0,1].

<details>
<summary>Hint</summary>
Minimize the variance of the importance-sampled estimator with respect to α.
</details>

<details>
<summary>Answer</summary>

The variance of importance sampling:
$$\text{Var}[I] = \int \frac{f^2(\omega)}{p_{guide}(\omega)} d\omega - I^2$$

where f(ω) is the integrand. Define the learning quality:
$$Q = \frac{\int p_{learned} \cdot f}{\int p_{BSDF} \cdot f}$$

Minimizing variance with respect to α yields:
$$\alpha_{opt} = \frac{Q}{1 + Q}$$

When Q = 0 (poor learning), α = 0 (use only BSDF)
When Q = 1 (perfect learning), α = 0.5 (balanced combination)
When Q → ∞ (learned distribution much better), α → 1
</details>

**Exercise 17.7:** Open Problem - Wavefront Coding for Computational Imaging
Design a phase mask ϕ(x,y) that simultaneously:
1. Extends depth of field by factor of 10
2. Maintains MTF > 0.5 at Nyquist frequency
3. Enables single-shot depth estimation

Describe your approach and derive the point spread function.

<details>
<summary>Hint</summary>
Consider cubic phase masks or optimized binary masks. Think about how phase affects the PSF's depth dependence.
</details>

<details>
<summary>Answer</summary>

This is an open research problem. A potential approach:

**Cubic Phase Mask:** $\phi(x,y) = \alpha(x^3 + y^3)$

The PSF becomes:
$$h(x,y;z) = \left|\mathcal{F}\{P(u,v)\exp[i\phi(u,v)]\exp[i\frac{z}{2k}(u^2+v^2)]\}\right|^2$$

Key insights:
1. Cubic phase makes PSF nearly depth-invariant over range Δz ∝ α^(2/3)
2. Deconvolution needed to restore image sharpness
3. PSF asymmetry enables depth from defocus
4. Optimization framework:

$$\min_\phi \int_{z_{min}}^{z_{max}} \text{Var}[h(x,y;z)] dz$$

subject to: $\text{MTF}(\nu_{Nyquist}) > 0.5$

Modern approaches use learned phase masks optimized end-to-end with reconstruction networks.
</details>

## Common Pitfalls and Debugging Tips

### Gotchas in Wavefront Shaping

1. **Phase Wrapping Ambiguity**
   - Problem: Phase measurements are modulo 2π
   - Solution: Use temporal unwrapping or multiple wavelengths
   - Debug: Check phase gradients for sudden jumps

2. **SLM Calibration Errors**
   - Problem: Nonlinear voltage-phase response
   - Solution: Measure full calibration curve
   - Debug: Test with known phase patterns (gratings)

3. **Coherence Requirements**
   - Problem: Partial coherence reduces contrast
   - Solution: Ensure coherence length > optical path differences
   - Debug: Measure visibility of interference fringes

4. **Sampling Aliasing**
   - Problem: SLM pixels undersample high-frequency phase
   - Solution: Apply anti-aliasing filters or increase sampling
   - Debug: Check Fourier spectrum for aliasing artifacts

5. **Polarization Sensitivity**
   - Problem: SLMs are polarization-dependent
   - Solution: Use proper input polarization state
   - Debug: Test with polarizer/analyzer pairs

### Performance Optimization Tips

1. **Memory Effect Exploitation**
   - Measure once, scan within memory effect range
   - Use grouped optimization for nearby targets

2. **Computational Shortcuts**
   - FFT for Fresnel propagation when applicable
   - Lookup tables for nonlinear phase responses
   - GPU acceleration for matrix operations

3. **Noise Reduction**
   - Average multiple measurements
   - Use lock-in detection for weak signals
   - Implement outlier rejection in optimization

## Best Practices Checklist

### System Design Review

- [ ] **Coherence Analysis**
  - Coherence length > maximum path difference?
  - Temporal stability sufficient for measurement time?
  - Spatial coherence matches optical system?

- [ ] **Sampling Considerations**
  - Nyquist criterion satisfied for phase patterns?
  - SLM resolution adequate for required NA?
  - Camera dynamic range sufficient?

- [ ] **Calibration Protocol**
  - Phase-voltage calibration measured?
  - Pixel crosstalk characterized?
  - System aberrations compensated?

### Algorithm Implementation

- [ ] **Optimization Strategy**
  - Appropriate algorithm for N (iterative vs. matrix)?
  - Convergence criteria defined?
  - Local minima avoidance strategy?

- [ ] **Noise Handling**
  - SNR estimation implemented?
  - Robust optimization for noisy measurements?
  - Regularization for ill-posed problems?

- [ ] **Validation Methods**
  - Ground truth comparison available?
  - Performance metrics defined?
  - Edge case testing completed?

### Integration with Graphics

- [ ] **Mathematical Consistency**
  - Units and conventions match graphics pipeline?
  - Coordinate systems properly transformed?
  - Energy conservation maintained?

- [ ] **Performance Scaling**
  - Complexity analysis performed?
  - Real-time constraints achievable?
  - Memory requirements acceptable?

- [ ] **Generalization**
  - Method extends to relevant graphics problems?
  - Assumptions clearly documented?
  - Failure modes understood?
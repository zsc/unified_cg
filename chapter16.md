# Chapter 16: Coherence Theory

This chapter introduces the fundamental concepts of optical coherence, bridging the gap between deterministic wave optics and statistical descriptions of light fields. We develop the mathematical framework for describing partial coherence, establishing key theorems that govern how coherence properties transform through optical systems and propagate through space. These concepts are essential for understanding advanced rendering effects involving interference and diffraction with realistic light sources.

## 16.1 Temporal Coherence and Spectra

Temporal coherence describes the correlation of a light wave with itself at different time delays. For a scalar optical field E(t), we define the temporal coherence function:

**Γ(τ) = ⟨E*(t)E(t+τ)⟩**

where ⟨·⟩ denotes time averaging and * indicates complex conjugation.

### 16.1.1 Coherence Time and Length

The normalized degree of temporal coherence is:

**γ(τ) = Γ(τ)/Γ(0)**

The coherence time τ_c is typically defined as the time delay at which |γ(τ)| falls to 1/e or 1/2:

**τ_c = ∫₀^∞ |γ(τ)|² dτ**

The coherence length is then:

**l_c = c·τ_c**

where c is the speed of light.

### 16.1.2 Spectral Width Relationship

By the Wiener-Khinchin theorem (detailed in Section 16.3), the temporal coherence function is the Fourier transform of the power spectral density S(ω):

**Γ(τ) = ∫_{-∞}^∞ S(ω)e^{iωτ} dω**

For a Gaussian spectrum with FWHM Δν:

**S(ω) = S₀ exp[-(ω-ω₀)²/(2σ²)]**

where σ = 2πΔν/(2√(2ln2)), the coherence time is:

**τ_c ≈ 0.44/Δν**

This fundamental relationship shows that narrow spectral lines produce long coherence times.

### 16.1.3 Examples of Coherence Times

- Laser (single mode): τ_c ~ 10⁻³ s, l_c ~ 300 km
- LED: τ_c ~ 10⁻¹⁴ s, l_c ~ 3 μm  
- Sunlight: τ_c ~ 10⁻¹⁵ s, l_c ~ 0.3 μm
- Sodium lamp: τ_c ~ 10⁻¹² s, l_c ~ 0.3 mm

## 16.2 Spatial Coherence and Young's Experiment

Spatial coherence describes the correlation between light fields at different spatial points at the same time. Young's double-slit experiment provides the canonical framework for understanding spatial coherence.

### 16.2.1 Mutual Intensity and Visibility

Consider two points P₁ and P₂ illuminated by a source. The mutual intensity is:

**J₁₂ = ⟨E*(P₁,t)E(P₂,t)⟩**

In Young's experiment, the intensity at observation point P is:

**I(P) = I₁ + I₂ + 2Re[J₁₂ exp(ikΔ)]**

where Δ is the path difference. The visibility of fringes is:

**V = (I_max - I_min)/(I_max + I_min) = 2√(I₁I₂)|γ₁₂|/(I₁ + I₂)**

For equal intensities (I₁ = I₂), V = |γ₁₂|.

### 16.2.2 Complex Degree of Coherence

The complex degree of coherence is defined as:

**γ₁₂ = J₁₂/√(J₁₁J₂₂)**

where J₁₁ = I₁ and J₂₂ = I₂ are the intensities at points 1 and 2.

Properties:
- |γ₁₂| ≤ 1 (equality for fully coherent light)
- γ₁₁ = γ₂₂ = 1
- γ₁₂ = γ₂₁*

### 16.2.3 Extended Sources

For an extended incoherent source, the degree of coherence depends on the source size and observation geometry. For a uniformly illuminated circular source of diameter D at distance z, observing points separated by d:

**|γ₁₂| ≈ |2J₁(πDd/λz)/(πDd/λz)|**

where J₁ is the first-order Bessel function. The first zero occurs at:

**d_coh ≈ 1.22λz/D**

This defines the transverse coherence length.

## 16.3 Mutual Coherence Function and Wiener-Khinchin Theorem

The mutual coherence function provides a complete statistical description of partially coherent fields.

### 16.3.1 General Definition

For a scalar field U(r,t), the mutual coherence function is:

**Γ(r₁,r₂,τ) = ⟨U*(r₁,t)U(r₂,t+τ)⟩**

The equal-time mutual coherence function (mutual intensity) is:

**J(r₁,r₂) = Γ(r₁,r₂,0)**

### 16.3.2 Cross-Spectral Density

Taking the Fourier transform with respect to τ:

**W(r₁,r₂,ω) = ∫_{-∞}^∞ Γ(r₁,r₂,τ)e^{-iωτ} dτ**

This cross-spectral density function satisfies:

**Γ(r₁,r₂,τ) = ∫_{-∞}^∞ W(r₁,r₂,ω)e^{iωτ} dω/(2π)**

### 16.3.3 Wiener-Khinchin Theorem

For stationary fields, when r₁ = r₂ = r:

**W(r,r,ω) = S(r,ω)**

is the power spectral density. The theorem states:

**⟨|U(r,t)|²⟩ = ∫_{-∞}^∞ S(r,ω) dω/(2π)**

This connects the time-averaged intensity to the integral of the spectrum.

### 16.3.4 Quasi-Monochromatic Approximation

For narrow-band light centered at ω₀:

**Γ(r₁,r₂,τ) ≈ J(r₁,r₂)e^{iω₀τ}g(τ)**

where g(τ) is a slowly varying envelope function. This simplifies many calculations while retaining essential physics.

## 16.4 van Cittert-Zernike Theorem

The van Cittert-Zernike theorem is one of the most important results in coherence theory, establishing how spatial coherence emerges from incoherent sources through propagation.

### 16.4.1 Statement of the Theorem

For a planar incoherent source with intensity distribution I_s(ξ,η) in the source plane, the complex degree of coherence in an observation plane at distance z is:

**γ₁₂ = ∫∫ I_s(ξ,η) exp[ik(r₁-r₂)·(ξ,η)/z] dξdη / ∫∫ I_s(ξ,η) dξdη**

In the paraxial approximation, this becomes:

**γ₁₂ = ℱ{I_s(ξ,η)}|_{(u,v)=(x₁-x₂)/λz, (y₁-y₂)/λz} / ∫∫ I_s(ξ,η) dξdη**

where ℱ denotes the 2D Fourier transform.

### 16.4.2 Physical Interpretation

The theorem reveals that:
1. The degree of coherence is the normalized Fourier transform of the source intensity distribution
2. Larger sources produce smaller coherence areas
3. The coherence function inherits the symmetry of the source

### 16.4.3 Examples and Applications

**Circular Source:**
For a uniform circular source of radius a:

**γ₁₂ = 2J₁(2πa|r₁-r₂|/λz)/(2πa|r₁-r₂|/λz)**

The coherence radius (first zero) is:
**ρ_c = 0.61λz/a**

**Rectangular Source:**
For a rectangular source of dimensions L_x × L_y:

**γ₁₂ = sinc(πL_x(x₁-x₂)/λz) × sinc(πL_y(y₁-y₂)/λz)**

**Double Star:**
For two point sources separated by angle θ:

**γ₁₂ = cos(πθ|r₁-r₂|/λ)**

This is the basis of stellar interferometry.

### 16.4.4 Generalized van Cittert-Zernike Theorem

For non-planar geometries and arbitrary propagation distances, the generalized form uses the Green's function:

**J(r₁,r₂) = ∫∫ I_s(r_s)G*(r₁,r_s)G(r₂,r_s) d²r_s**

where G is the appropriate Green's function for the geometry.

## 16.5 Propagation of Partially Coherent Light

Understanding how coherence properties change during propagation is crucial for accurate modeling of optical systems.

### 16.5.1 Propagation Law for Mutual Coherence

The mutual coherence function satisfies a pair of wave equations:

**∇₁²J(r₁,r₂) + k²J(r₁,r₂) = 0**
**∇₂²J(r₁,r₂) + k²J(r₁,r₂) = 0**

where ∇ᵢ² operates on coordinate rᵢ. These are known as the Wolf equations.

### 16.5.2 Propagation Through Free Space

For propagation from plane z=0 to plane z, using the Fresnel approximation:

**J(x₁,y₁,x₂,y₂;z) = (k/2πz)² ∫∫∫∫ J₀(ξ₁,η₁,ξ₂,η₂) × exp[ik/2z((x₁-ξ₁)² + (y₁-η₁)² - (x₂-ξ₂)² - (y₂-η₂)²)] dξ₁dη₁dξ₂dη₂**

### 16.5.3 Coherence Mode Representation

Any partially coherent field can be decomposed into coherent modes:

**J(r₁,r₂) = Σₙ λₙ φₙ*(r₁)φₙ(r₂)**

where λₙ are eigenvalues and φₙ are orthonormal eigenfunctions satisfying:

**∫ J(r₁,r₂)φₙ(r₂) d²r₂ = λₙφₙ(r₁)**

This is the coherent mode representation, analogous to principal component analysis.

### 16.5.4 Schell-Model Sources

A important class of sources satisfies:

**J(r₁,r₂) = √[I(r₁)I(r₂)] μ(r₁-r₂)**

where μ is the coherence function depending only on r₁-r₂. These Schell-model sources maintain this form during propagation:

**J_z(r₁,r₂) = √[I_z(r₁)I_z(r₂)] μ_z(r₁-r₂)**

### 16.5.5 Connection to Volume Rendering

In the context of volume rendering with partially coherent illumination, the radiance at each point becomes:

**L(x,ω) = ∫∫ W(x,x',ω)σ_s(x')p(x',ω'→ω) dω' d³x'**

where W(x,x',ω) is the cross-spectral density of the illumination. This generalizes the standard volume rendering equation to include coherence effects, essential for accurate modeling of:
- Laser scanning microscopy
- Optical coherence tomography
- Holographic displays
- Interferometric imaging

## Chapter Summary

This chapter established the mathematical framework for optical coherence, connecting statistical properties of light fields to observable interference phenomena. Key concepts include:

1. **Temporal Coherence**: Related to spectral bandwidth through Fourier transform relationships. Coherence time τ_c ≈ 1/Δν determines the path length difference over which interference can occur.

2. **Spatial Coherence**: Quantified by the complex degree of coherence γ₁₂, directly observable through fringe visibility in interference experiments.

3. **Mutual Coherence Function**: Γ(r₁,r₂,τ) provides complete second-order statistical description of partially coherent fields, with cross-spectral density W as its frequency domain counterpart.

4. **van Cittert-Zernike Theorem**: Establishes that spatial coherence in the far field is the Fourier transform of the source intensity distribution, fundamental for understanding coherence from extended sources.

5. **Propagation Laws**: Wolf equations govern how coherence properties evolve during propagation, maintaining the wave nature while incorporating statistical effects.

The volume rendering equation generalizes to:
**L(x,ω) = ∫∫ W(x,x',ω)σ_s(x')p(x',ω'→ω) dω' d³x'**

incorporating coherence effects through the cross-spectral density.

## Exercises

### Exercise 16.1: Coherence Time Calculation
A helium-neon laser has a spectral linewidth of 1 GHz. Calculate:
a) The coherence time τ_c
b) The coherence length l_c
c) The maximum path difference for which interference fringes have visibility > 0.5

*Hint: Use the relationship τ_c ≈ 0.44/Δν for Gaussian spectra.*

<details>
<summary>Solution</summary>

a) τ_c = 0.44/Δν = 0.44/(10⁹ Hz) = 4.4 × 10⁻¹⁰ s = 0.44 ns

b) l_c = c·τ_c = (3 × 10⁸ m/s)(4.4 × 10⁻¹⁰ s) = 0.132 m = 13.2 cm

c) For Gaussian spectrum, visibility V = exp(-π²τ²Δν²/2ln2)
   Setting V = 0.5: τ = √(2ln2·ln2)/πΔν ≈ 0.37/Δν
   Maximum path difference = c·τ = 0.37c/Δν = 11.1 cm
</details>

### Exercise 16.2: Young's Double Slit with Extended Source
A sodium lamp (λ = 589 nm) with circular aperture of diameter 2 mm illuminates Young's double slits separated by 0.5 mm at a distance of 1 m. Calculate:
a) The degree of coherence at the slits
b) The fringe visibility
c) The slit separation for which visibility drops to zero

*Hint: Apply the van Cittert-Zernike theorem for a circular source.*

<details>
<summary>Solution</summary>

a) Using van Cittert-Zernike for circular source:
   γ₁₂ = 2J₁(πDd/λz)/(πDd/λz)
   where D = 2 mm, d = 0.5 mm, z = 1 m, λ = 589 nm
   
   πDd/λz = π(2×10⁻³)(0.5×10⁻³)/(589×10⁻⁹)(1) = 5.33
   γ₁₂ = 2J₁(5.33)/5.33 ≈ 2(-0.327)/5.33 = -0.123

b) For equal intensity slits: V = |γ₁₂| = 0.123

c) First zero when πDd/λz = 3.83 (first zero of J₁)
   d = 3.83λz/πD = 3.83(589×10⁻⁹)(1)/π(2×10⁻³) = 0.36 mm
</details>

### Exercise 16.3: Wiener-Khinchin Application
A light source has a Lorentzian spectrum: S(ω) = S₀Γ²/[(ω-ω₀)² + Γ²]
Derive:
a) The temporal coherence function Γ(τ)
b) The coherence time
c) Compare with a Gaussian spectrum of same FWHM

*Hint: Use contour integration for the Fourier transform.*

<details>
<summary>Solution</summary>

a) Γ(τ) = ∫ S(ω)e^{iωτ} dω = S₀Γ² ∫ e^{iωτ}/[(ω-ω₀)² + Γ²] dω
   
   Using residue theorem with pole at ω = ω₀ + iΓ (for τ > 0):
   Γ(τ) = 2πiS₀Γ²·e^{i(ω₀+iΓ)τ}/(2iΓ) = πS₀Γe^{iω₀τ}e^{-Γτ}
   
   Normalizing: γ(τ) = e^{iω₀τ}e^{-Γ|τ|}

b) τ_c = ∫₀^∞ |γ(τ)|² dτ = ∫₀^∞ e^{-2Γτ} dτ = 1/(2Γ)
   
   For FWHM = 2Γ: τ_c = 1/FWHM

c) Gaussian: τ_c = 0.44/Δν
   Lorentzian: τ_c = 1/(2πΔν) ≈ 0.16/Δν
   The Lorentzian has shorter coherence time due to extended wings.
</details>

### Exercise 16.4: Coherence Mode Decomposition
For a Schell-model source with Gaussian intensity I(x) = I₀exp(-x²/w₀²) and Gaussian coherence μ(Δx) = exp(-Δx²/2σ_c²), find the first three coherent modes.

*Hint: Use Hermite-Gaussian functions as basis.*

<details>
<summary>Solution</summary>

The eigenvalue equation: ∫ J(x,x')φₙ(x') dx' = λₙφₙ(x)

For this Gaussian Schell-model, eigenfunctions are Hermite-Gaussians:
φₙ(x) = (2^n n!√π σ)^{-1/2} Hₙ(x/σ) exp(-x²/2σ²)

where σ⁴ = w₀²σ_c²/2

Eigenvalues: λₙ = λ₀(σ_c²/(σ_c² + w₀²))^n
where λ₀ = I₀√(2πσ_c²w₀²/(σ_c² + w₀²))

First three modes:
- n=0: φ₀(x) = (πσ²)^{-1/4} exp(-x²/2σ²), λ₀
- n=1: φ₁(x) = (πσ²)^{-1/4} √(2)x/σ exp(-x²/2σ²), λ₁
- n=2: φ₂(x) = (πσ²)^{-1/4} (2x²/σ² - 1)/√2 exp(-x²/2σ²), λ₂
</details>

### Exercise 16.5: Propagation of Coherence (Challenge)
A partially coherent beam has initial mutual intensity J₀(x₁,x₂) = exp(-|x₁-x₂|²/2σ₀²)exp(-(x₁²+x₂²)/4w₀²). Find J(x₁,x₂,z) after propagating distance z.

*Hint: Use the Fresnel propagation integral for mutual intensity.*

<details>
<summary>Solution</summary>

Using Fresnel propagation:
J(x₁,x₂,z) = (k/2πz)² ∫∫ J₀(ξ₁,ξ₂) exp[ik(x₁-ξ₁)²/2z] exp[-ik(x₂-ξ₂)²/2z] dξ₁dξ₂

Substituting the Gaussian form and completing the square:
J(x₁,x₂,z) = A(z) exp(-|x₁-x₂|²/2σ²(z)) exp(-(x₁²+x₂²)/4w²(z))

where:
- w²(z) = w₀²(1 + z²/z_R²), z_R = πw₀²/λ
- σ²(z) = σ₀² + λ²z²/(4π²σ₀²)
- A(z) includes normalization factors

The beam maintains Gaussian-Schell form with evolving parameters.
</details>

### Exercise 16.6: Van Cittert-Zernike for Stellar Interferometry (Challenge)
Two telescopes separated by baseline B observe a binary star with angular separation θ and intensity ratio R. Derive the visibility curve V(B) and show how to extract θ and R.

*Hint: Model as two incoherent point sources.*

<details>
<summary>Solution</summary>

For two stars at angles ±θ/2 with intensities I₁, I₂:
γ₁₂ = [I₁exp(ikθB/2) + I₂exp(-ikθB/2)]/(I₁ + I₂)

Let R = I₂/I₁, then:
γ₁₂ = [exp(ikθB/2) + R·exp(-ikθB/2)]/(1 + R)
    = [(1+R)cos(kθB/2) + i(1-R)sin(kθB/2)]/(1 + R)

Visibility: V(B) = |γ₁₂| = √[cos²(πθB/λ) + ((1-R)/(1+R))²sin²(πθB/λ)]

Analysis:
- At B = 0: V = 1 (full coherence)
- First minimum at πθB/λ = π/2 gives θ = λ/2B_min
- Visibility at minimum: V_min = |1-R|/(1+R) gives R
- For equal stars (R=1): V = |cos(πθB/λ)|
</details>

### Exercise 16.7: Coherence in Volume Rendering
Derive the modified volume rendering equation for a laser-scanned microscope where the illumination has Gaussian spatial coherence with radius ρ_c.

*Hint: Start with the coherent volume rendering equation and average over the coherence function.*

<details>
<summary>Solution</summary>

Start with coherent illumination at point x':
L_coh(x,ω) = ∫ E(x')σ_s(x')p(x',ω'→ω)G(x',x) d³x'

For partially coherent illumination with mutual intensity J(x₁,x₂):
L(x,ω) = ∫∫ J(x₁,x₂)σ_s(x₁)σ_s*(x₂)p(x₁,ω'→ω)p*(x₂,ω'→ω)G(x₁,x)G*(x₂,x) d³x₁d³x₂

For Gaussian coherence: J(x₁,x₂) = I(x₁)exp(-|x₁-x₂|²/2ρ_c²)

In the limit ρ_c → 0 (incoherent):
L(x,ω) = ∫ I(x')|σ_s(x')|²|p(x',ω'→ω)|²|G(x',x)|² d³x'

For finite ρ_c, coherent volume effects occur within radius ρ_c, leading to speckle patterns with characteristic size ~λz/ρ_c at distance z.
</details>

### Exercise 16.8: Cross-Spectral Purity (Open-ended)
Investigate conditions under which a partially coherent field can be spectrally decomposed such that each frequency component is fully coherent. What physical sources satisfy this condition?

*Hint: Consider when W(r₁,r₂,ω) factorizes as U*(r₁,ω)U(r₂,ω).*

<details>
<summary>Solution</summary>

A field is cross-spectrally pure if:
W(r₁,r₂,ω) = U*(r₁,ω)U(r₂,ω)S(ω)

This requires the spectral degree of coherence:
μ(r₁,r₂,ω) = W(r₁,r₂,ω)/√[W(r₁,r₁,ω)W(r₂,r₂,ω)] = 1

for all ω where S(ω) ≠ 0.

Physical examples:
1. **Filtered thermal light**: Passing white light through narrow-band filter
2. **Mode-locked lasers**: Each longitudinal mode is coherent
3. **Stationary Schell-model sources**: Under specific propagation conditions

Counter-examples:
- Moving sources (Doppler broadening breaks spectral purity)
- Nonlinear processes (frequency mixing)
- Time-varying media

Implications for rendering: Cross-spectrally pure sources allow frequency-by-frequency calculation of coherence effects, greatly simplifying computations.
</details>

## Common Pitfalls and Errors

1. **Confusing Coherence Time and Coherence Length**
   - Error: Using τ_c where l_c is needed
   - Fix: Remember l_c = c·τ_c, coherence length has units of distance

2. **Misapplying van Cittert-Zernike Theorem**
   - Error: Using it for coherent or partially coherent sources
   - Fix: Theorem only applies to incoherent sources; check source properties first

3. **Incorrect Normalization of Degree of Coherence**
   - Error: Forgetting to normalize by √(I₁I₂)
   - Fix: Always compute γ₁₂ = J₁₂/√(J₁₁J₂₂)

4. **Sign Errors in Complex Coherence**
   - Error: Ignoring phase of γ₁₂ in interference calculations
   - Fix: Keep track of complex nature; phase affects fringe position

5. **Assuming Coherence is Preserved**
   - Error: Treating coherence as invariant during propagation
   - Fix: Use Wolf equations or appropriate propagation laws

6. **Mixing Temporal and Spatial Coherence**
   - Error: Using temporal coherence formulas for spatial problems
   - Fix: Identify whether dealing with time delays (temporal) or spatial separations

7. **Quasi-Monochromatic Approximation Misuse**
   - Error: Applying to broadband sources
   - Fix: Check Δν/ν₀ << 1 before using approximation

8. **Forgetting Statistical Averaging**
   - Error: Using instantaneous fields instead of ensemble averages
   - Fix: All coherence functions involve time or ensemble averaging ⟨·⟩

## Best Practices Checklist

### Design Review
- [ ] Identified whether dealing with temporal, spatial, or both types of coherence
- [ ] Verified source coherence properties (coherent/partially coherent/incoherent)
- [ ] Selected appropriate coherence measure (Γ, J, W, γ)
- [ ] Checked validity of approximations (paraxial, quasi-monochromatic)

### Mathematical Verification
- [ ] Coherence functions properly normalized
- [ ] Complex conjugates placed correctly
- [ ] Fourier transform pairs consistent (time ↔ frequency, space ↔ spatial frequency)
- [ ] Units dimensionally correct throughout

### Physical Constraints
- [ ] |γ₁₂| ≤ 1 everywhere
- [ ] Coherence functions Hermitian: Γ(r₁,r₂,τ) = Γ*(r₂,r₁,-τ)
- [ ] Power spectral density non-negative: S(ω) ≥ 0
- [ ] Energy conservation in propagation

### Computational Considerations
- [ ] Sampling sufficient for coherence scale (Nyquist for ρ_c)
- [ ] Numerical integration stable for oscillatory integrands
- [ ] Mode truncation error estimated for coherent mode expansion
- [ ] Statistical convergence verified for Monte Carlo approaches

### Experimental Validation
- [ ] Coherence length matches spectral width measurements
- [ ] Visibility measurements consistent with theoretical γ₁₂
- [ ] Propagation effects match Wolf equation predictions
- [ ] Van Cittert-Zernike validated for extended sources
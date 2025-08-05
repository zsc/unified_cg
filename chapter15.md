# Chapter 15: Scalar Wave Optics Foundations

In this chapter, we transition from geometric optics and volume rendering to wave optics, establishing the mathematical foundation for understanding light as an electromagnetic wave. We'll derive the scalar wave approximation from Maxwell's equations and explore fundamental diffraction phenomena that become crucial when wavelength-scale effects matter. This bridge from ray-based to wave-based descriptions enriches our understanding of light transport and sets the stage for advanced optical phenomena in computer graphics.

The transition from rays to waves fundamentally changes how we model light transport. Where geometric optics treats light as infinitesimal rays following straight paths, wave optics reveals that light spreads, diffracts around edges, and interferes with itself. These effects become essential when:
- Feature sizes approach the wavelength of light (~400-700nm)
- Coherent illumination is present (lasers, some LEDs)
- High-fidelity rendering of optical phenomena is required
- Microscale surface structures create visual effects

We'll see how the volume rendering equation naturally extends to include wave phenomena through the Green's function formalism, providing a unified framework that encompasses both ray and wave regimes.

## Mathematical Foundations and Context

Before diving into wave equations, let's establish the mathematical context. The transition from geometric to wave optics represents a fundamental shift in how we describe light propagation:

**Geometric Optics**: Light intensity I(x,ω) follows rays according to:
dI/ds = -σₜI along ray parameterized by s

**Wave Optics**: Complex field amplitude u(x,t) satisfies wave equations:
∇²u - (1/c²)∂²u/∂t² = 0

The connection emerges through the **eikonal approximation**. For highly oscillatory fields u = A exp(ikS), where k >> 1:
- Amplitude A varies slowly
- Phase S satisfies the eikonal equation: |∇S|² = n²
- Rays are orthogonal to surfaces of constant phase

This chapter explores what happens when k is finite, revealing diffraction, interference, and the wave nature of light.

## 15.1 From Maxwell's Equations to the Helmholtz Equation

### Vector Wave Equation

We begin with Maxwell's equations in a source-free, homogeneous medium:

∇ × **E** = -∂**B**/∂t  (Faraday's law)
∇ × **H** = ∂**D**/∂t   (Ampère-Maxwell law)
∇ · **D** = 0           (No free charges)
∇ · **B** = 0           (No magnetic monopoles)

For linear, isotropic media: **D** = ε**E** and **B** = μ**H**, where ε = ε₀εᵣ and μ = μ₀μᵣ.

These equations embody fundamental electromagnetic principles:
- **Faraday's law**: Time-varying magnetic fields induce electric fields
- **Ampère-Maxwell law**: Time-varying electric fields (displacement current) and conduction currents create magnetic fields
- **Gauss's law**: Electric field divergence relates to charge density
- **No monopoles**: Magnetic field lines are always closed loops

Taking the curl of Faraday's law:
∇ × (∇ × **E**) = -∇ × (∂**B**/∂t) = -∂(∇ × **B**)/∂t = -μ∂(∇ × **H**)/∂t

Using the vector identity ∇ × (∇ × **E**) = ∇(∇ · **E**) - ∇²**E** and noting that ∇ · **E** = 0 in source-free regions:

-∇²**E** = -μ∂(∇ × **H**)/∂t = -μ∂(∂**D**/∂t)/∂t = -με∂²**E**/∂t²

This yields the vector wave equation:

∇²**E** - με(∂²**E**/∂t²) = 0

Or in more compact form:
∇²**E** - (1/v²)(∂²**E**/∂t²) = 0

The wave velocity is v = 1/√(με) = c/n, where:
- c = 1/√(μ₀ε₀) ≈ 2.998×10⁸ m/s is the speed of light in vacuum
- n = √(εᵣμᵣ) ≈ √εᵣ is the refractive index (since μᵣ ≈ 1 for most optical materials)

An identical equation holds for the magnetic field **H**. These vector equations couple the three spatial components of the fields through:
1. **Boundary conditions** at material interfaces
2. **Transversality constraint**: ∇ · **E** = 0 implies **k** · **E** = 0 for plane waves
3. **Impedance relations**: Z = √(μ/ε) relates **E** and **H** magnitudes

### Mathematical Structure

The vector wave equation exhibits several key mathematical properties:

**Linearity**: Solutions can be superposed
If **E₁** and **E₂** are solutions, then α**E₁** + β**E₂** is also a solution

**Time-reversal symmetry**: Replace t → -t yields valid solutions
Forward and backward propagating waves are equally valid

**Gauge invariance**: In vacuum, we can choose ∇ · **A** = 0 (Coulomb gauge) or 
∂Φ/∂t + ∇ · **A** = 0 (Lorenz gauge), where **E** = -∇Φ - ∂**A**/∂t

**Energy conservation**: The Poynting vector **S** = **E** × **H** satisfies:
∂u/∂t + ∇ · **S** = 0
where u = (ε|**E**|² + μ|**H**|²)/2 is the electromagnetic energy density

### Scalar Wave Approximation

For many optical phenomena, we can approximate the vector field with a scalar field U(r,t). This approximation is valid when:
- The medium is homogeneous over wavelength scales (∇n·λ << n)
- Polarization effects are negligible (unpolarized or fixed polarization)
- The field varies slowly compared to wavelength (paraxial approximation)
- We're far from material interfaces where boundary conditions couple components

#### Rigorous Derivation

To derive the scalar approximation systematically, we start with the vector Helmholtz equation for each component. Consider a predominantly z-propagating wave with electric field:

**E** = E_x **x̂** + E_y **ŷ** + E_z **ẑ**

From Maxwell's equations, the transversality condition ∇ · **E** = 0 gives:
∂E_x/∂x + ∂E_y/∂y + ∂E_z/∂z = 0

For paraxial waves where ∂/∂z ~ ik (rapid phase variation) but transverse derivatives are small:
E_z ≈ -(1/ik)(∂E_x/∂x + ∂E_y/∂y)

This shows E_z << E_x, E_y for paraxial propagation, justifying focus on transverse components.

Each transverse component satisfies:
∇²E_⊥ - (1/v²)(∂²E_⊥/∂t²) = 0

For monochromatic fields with angular frequency ω:
E_⊥(r,t) = Re[e_⊥(r)e^(-iωt)]

The complex amplitude representation separates time and space:
- Temporal: e^(-iωt) with ∂²/∂t² → -ω²
- Spatial: e_⊥(r) contains all spatial variation

Substituting:
∇²e_⊥ + (ω²/v²)e_⊥ = 0

Defining the wavenumber k = ω/v = 2πn/λ, we get:

∇²u + k²u = 0

This is the **Helmholtz equation**, where u represents any scalar component of the field.

#### Validity Limits

The scalar approximation breaks down when:

1. **High numerical aperture** (NA > 0.6):
   - Vector effects: longitudinal fields become significant
   - Polarization coupling in tight focusing
   - Use vector diffraction theory (Richards-Wolf)

2. **Near material interfaces**:
   - Boundary conditions couple field components
   - Fresnel coefficients depend on polarization
   - Surface plasmons and guided modes

3. **Subwavelength structures**:
   - Near-field enhancement
   - Evanescent waves dominate
   - Full vector treatment required

4. **Birefringent media**:
   - Different propagation for orthogonal polarizations
   - Coupled wave equations
   - Jones or Mueller calculus needed

### Physical Interpretation

The Helmholtz equation describes monochromatic wave propagation where:
- k = 2πn/λ represents the spatial frequency (rad/m)
- The equation balances spatial curvature (∇²u) against phase accumulation (k²u)
- Solutions form a complete basis for arbitrary fields

#### Fundamental Solutions

1. **Plane waves**: u = A exp(i**k**·**r**)
   - Wavevector **k** = k(sin θ cos φ **x̂** + sin θ sin φ **ŷ** + cos θ **ẑ**)
   - |**k**| = k = 2πn/λ
   - Constant amplitude surfaces perpendicular to **k**
   - Energy flux along **k** direction
   - Basis for angular spectrum representation

   Verification: ∇²[exp(i**k**·**r**)] = -k²exp(i**k**·**r**) ✓

2. **Spherical waves**: u = (A/r)exp(±ikr)
   - Point source/sink at origin
   - ± for outgoing/incoming waves
   - Amplitude ∝ 1/r (energy conservation)
   - Intensity ∝ 1/r² (inverse square law)
   - Phase surfaces are concentric spheres

   In spherical coordinates with radial symmetry:
   ∇²u = (1/r²)d/dr(r²du/dr) = (A/r²)d/dr[r²d/dr((1/r)e^(ikr))]
   After calculation: ∇²u = -k²u ✓

3. **Gaussian beams**: u = (A₀/q(z))exp[ikz + ikr²/2q(z)]
   
   Complex beam parameter: q(z) = z - iz₀ where z₀ = πw₀²/λ
   
   Beam properties:
   - Beam width: w(z) = w₀√(1 + (z/z₀)²)
   - Wavefront radius: R(z) = z(1 + (z₀/z)²)
   - Gouy phase: ζ(z) = arctan(z/z₀)
   - Rayleigh range: z₀ (beam doubles in area)
   
   Near axis (r << w): Satisfies paraxial wave equation
   ∂²u/∂x² + ∂²u/∂y² + 2ik∂u/∂z = 0

4. **Bessel beams**: u = J₀(k_⊥r)exp(ik_z z)
   - Non-diffracting solution
   - k_⊥² + k_z² = k²
   - Infinite energy (not physically realizable)
   - Approximated by finite apertures

5. **Hermite-Gaussian modes**: u_{mn} = H_m(√2x/w)H_n(√2y/w)exp(-r²/w²)×[Gaussian beam factor]
   - H_m: Hermite polynomials
   - Orthogonal mode basis
   - Rectangular symmetry
   - Important for laser cavities

6. **Laguerre-Gaussian modes**: u_{pl} = (r/w)^|l| L_p^|l|(2r²/w²)exp(-r²/w²)exp(ilφ)×[Gaussian beam factor]
   - L_p^|l|: Associated Laguerre polynomials
   - Orbital angular momentum: l𝗁 per photon
   - Cylindrical symmetry
   - Optical vortices for l ≠ 0

### Connection to Volume Rendering

The Helmholtz equation naturally connects to our volume rendering framework through the Green's function formalism. This connection reveals how wave optics emerges from and extends the radiative transfer equation.

#### Green's Function Formulation

Consider the frequency-domain rendering equation with coherent illumination:

L(x,ω) = L₀(x,ω) + ∫ σₛ(x')p(x',ω'→ω)G(x,x')L(x',ω')dV'

The Green's function G(x,x') represents coherent propagation from x' to x and satisfies:

(∇² + k²)G(x,x') = -δ(x-x')

This inhomogeneous Helmholtz equation has the fundamental solution:
G(x,x') = exp(ik|x-x'|)/(4π|x-x'|)

Physical interpretation:
- Outgoing spherical wave from point source at x'
- Phase accumulation: k|x-x'|
- Amplitude decay: 1/|x-x'|
- Satisfies Sommerfeld radiation condition

#### Scale-Dependent Regimes

The parameter kr = k|x-x'| determines the propagation regime:

1. **Geometric optics limit** (kr >> 1):
   - Rapid phase oscillation: exp(ikr)
   - Stationary phase approximation applies
   - Only paths with ∇φ = 0 contribute
   - Green's function → δ-function along rays
   
   G(x,x') ≈ δ(s - |x-x'|)/|∂s/∂x'|
   
   where s parameterizes the ray from x' to x

2. **Wave regime** (kr ~ 1):
   - Phase and amplitude comparable
   - Diffraction effects significant
   - Interference between multiple paths
   - Full Green's function needed

3. **Near-field** (kr << 1):
   - Static field approximation
   - G(x,x') ≈ 1/(4π|x-x'|) (Coulomb-like)
   - Evanescent waves dominate
   - Non-propagating near-field coupling

#### Coherent vs Incoherent Rendering

The transition between coherent and incoherent rendering depends on source coherence:

**Coherent sources** (lasers, single-mode fibers):
- Fields add: E_total = E₁ + E₂ + ...
- Cross terms: |E_total|² = |E₁|² + |E₂|² + 2Re(E₁*E₂) + ...
- Interference fringes with visibility V = |γ₁₂|

**Partially coherent sources** (LEDs, thermal):
- Mutual coherence function: Γ(x₁,x₂) = ⟨E*(x₁)E(x₂)⟩
- Van Cittert-Zernike theorem relates source size to coherence
- Coherence length: l_c = λ²/(Δλ) 
- Coherence area: A_c = λ²R²/A_source

**Incoherent limit** (most rendering):
- Intensities add: I_total = I₁ + I₂ + ...
- No interference terms
- Ray optics sufficient
- Standard rendering equation applies

#### Extended Volume Rendering Equation

For partial coherence, the rendering equation generalizes to:

L(x₁,x₂,ω) = L₀(x₁,x₂,ω) + ∬ σₛ(x'₁)σₛ(x'₂)p₁p₂G*(x₁,x'₁)G(x₂,x'₂)L(x'₁,x'₂,ω)dx'₁dx'₂

This 6D equation reduces to:
- Standard rendering (x₁ = x₂, diagonal terms only)
- Coherent rendering (factorizable L)
- Speckle/interference (off-diagonal terms)

#### Practical Implications

1. **Multi-scale rendering**:
   - Geometric optics: λ << feature size
   - Wave corrections: λ ~ feature size  
   - Full wave solution: λ >> feature size

2. **Unified algorithms**:
   - Path tracing with phase tracking
   - Beam propagation methods
   - Hybrid ray-wave techniques

3. **New phenomena**:
   - Diffraction from edges/apertures
   - Interference in thin films
   - Speckle from rough surfaces
   - Focusing and caustics beyond geometric limit

## 15.2 Huygens-Fresnel Principle

### Historical Development

Christiaan Huygens (1678) proposed that each point on a wavefront acts as a source of secondary spherical wavelets. Augustin-Jean Fresnel (1815) added the principle of interference, explaining diffraction patterns through the coherent superposition of these wavelets.

### Mathematical Formulation

Consider a wavefront Σ at time t. The field at point P at time t + Δt is:

u(P) = (1/iλ) ∫∫_Σ u(Q) (e^(ikr))/r K(χ) dS

where:
- Q is a point on the wavefront Σ
- r = |P - Q| is the distance
- K(χ) is the obliquity factor
- χ is the angle between normal and P-Q direction

### Obliquity Factor

Fresnel originally proposed K(χ) = (1 + cos χ)/2, which:
- Equals 1 for forward propagation (χ = 0)
- Equals 0 for backward propagation (χ = π)
- Provides smooth angular dependence

The physical meaning:
- cos χ term: projection of wavelet onto observation direction
- Constant term: isotropic contribution
- Together: cardioid radiation pattern

This obliquity factor ensures:
1. No backward propagating waves (causality)
2. Maximum contribution in forward direction
3. Smooth variation preventing discontinuities
4. Energy conservation in the far field

### Kirchhoff's Rigorous Formulation

Gustav Kirchhoff (1882) derived the Huygens-Fresnel principle from the Helmholtz equation using Green's theorem:

u(P) = (1/4π) ∫∫_Σ [e^(ikr)/r ∂u/∂n - u ∂/∂n(e^(ikr)/r)] dS

For an aperture in an opaque screen with incident field u_inc:
- On aperture: u = u_inc, ∂u/∂n = ∂u_inc/∂n
- On screen: u = 0, ∂u/∂n = 0

This yields the **Kirchhoff diffraction formula**:

u(P) = (1/iλ) ∫∫_aperture u_inc(Q) (e^(ikr))/r (1 + cos χ)/2 dS

### Connection to Rendering

The Huygens-Fresnel principle parallels importance sampling in rendering:

| Wave Optics | Rendering |
|-------------|----------|
| Secondary sources | Sample points |
| Wavelet superposition | Monte Carlo integration |
| Obliquity factor K(χ) | Cosine weighting (N·L) |
| Coherent addition | Complex phasor sum |
| Intensity = |∑ fields|² | Radiance accumulation |

Key differences:
- Wave optics: Complex amplitudes with phase
- Rendering: Real-valued intensities
- Coherence introduces interference not present in incoherent rendering

This suggests extensions to rendering:
1. Complex-valued path tracing for coherent sources
2. Phase-aware importance sampling
3. Interference effects in material models

## 15.3 Fresnel Diffraction Integral

### Near-Field Geometry

Consider a planar aperture in the z=0 plane illuminated by a field u₀(x₀,y₀). The field at observation point P(x,y,z) is given by the Kirchhoff integral. For near-field diffraction, we must carefully expand the distance r that appears in both the amplitude and phase terms.

Let **r** = (x,y,z) be the observation point and **r₀** = (x₀,y₀,0) be a point in the aperture, then:

r = |**r** - **r₀**| = √[(x-x₀)² + (y-y₀)² + z²]

The key insight is that phase varies much more rapidly than amplitude:
- Phase variation: kr ~ 10⁶ rad/m for visible light
- Amplitude variation: 1/r changes slowly over wavelength scales

This allows different approximation orders for phase and amplitude terms.

### Fresnel Approximation

For z >> (x-x₀), (y-y₀), we use the binomial expansion:

r = z√[1 + ((x-x₀)² + (y-y₀)²)/z²]

Let ρ² = (x-x₀)² + (y-y₀)². Using (1+ε)¹/² ≈ 1 + ε/2 - ε²/8 + ... for ε << 1:

r ≈ z[1 + ρ²/2z² - ρ⁴/8z⁴ + ...]

For the phase term kr, we keep terms that contribute phase errors < π/2:
- First order: kz
- Second order: kρ²/2z (Fresnel term)
- Third order: -kρ⁴/8z³ (usually neglected)

For the amplitude term 1/r, we keep only the leading term:
1/r ≈ 1/z

This yields:
e^(ikr)/r ≈ (e^(ikz)/z) exp[ikρ²/2z] = (e^(ikz)/z) exp[ik/2z((x-x₀)² + (y-y₀)²)]

### Fresnel Diffraction Formula

Substituting into the Kirchhoff integral:

u(x,y,z) = (e^(ikz)/iλz) ∫∫_aperture u₀(x₀,y₀) exp[ik/2z((x-x₀)² + (y-y₀)²)] dx₀dy₀

Expanding the quadratic term:

u(x,y,z) = (e^(ikz)/iλz) exp[ik/2z(x² + y²)] × 
           ∫∫ u₀(x₀,y₀) exp[ik/2z(x₀² + y₀²)] exp[-ik/z(xx₀ + yy₀)] dx₀dy₀

### Validity Conditions

The Fresnel approximation validity depends on the phase error from neglected terms. The quartic term contributes a phase:

Φ₄ = -kρ₄/8z³

Requiring |Φ₄|_max < π/2:

kρ₄_max/8z³ < π/2

Substituting k = 2π/λ and solving:

z³ > ρ₄_max/4λ = [(x-x₀)² + (y-y₀)²]²_max/(4λ)

Define the **Fresnel number**:

F = a²/(λz)

where a is the characteristic aperture dimension. The approximation regimes:

1. **F >> 1**: Geometric shadow (ray optics)
2. **F ~ 1**: Fresnel diffraction (near field)
3. **F << 1**: Fraunhofer diffraction (far field)

Physical interpretation:
- F compares aperture area (a²) to diffraction area (λz)
- Large F: Many Fresnel zones visible, geometric limit
- Small F: Single Fresnel zone, pure diffraction

### Computational Methods

#### 1. Direct Integration
Numerical quadrature of the Fresnel integral:

u(x,y,z) = (e^(ikz)/iλz) ∬ u₀(x₀,y₀) exp[ik/2z((x-x₀)² + (y-y₀)²)] dx₀dy₀

- Complexity: O(N⁴) for N×N grids
- Accurate but computationally prohibitive
- Useful for irregular apertures or sparse sampling

#### 2. FFT Convolution Method
Rewrite the Fresnel integral as a convolution:

u(x,y,z) = C × [u₀(x,y)exp(ik(x²+y²)/2z)] ⊗ exp(ik(x²+y²)/2z)

where C = exp(ikz)/(iλz) and ⊗ denotes convolution.

Implementation:
1. Multiply input by quadratic phase (chirp)
2. FFT to frequency domain
3. Multiply by transfer function
4. Inverse FFT
5. Multiply by output chirp

- Complexity: O(N²log N)
- Requires careful sampling to avoid aliasing
- Zero-padding needed for accuracy

#### 3. Angular Spectrum Method
Propagate the field in the spatial frequency domain:

u(x,y,z) = ℱ⁻¹{ℱ{u₀(x₀,y₀)} × H(fₓ,fᵧ,z)}

where the transfer function:
H(fₓ,fᵧ,z) = exp[ikz√(1-(λfₓ)²-(λfᵧ)²)]

For (λfₓ)² + (λfᵧ)² < 1: Propagating waves
For (λfₓ)² + (λfᵧ)² > 1: Evanescent waves (exponential decay)

Advantages:
- Most efficient: O(N²log N)
- Exact within sampling limits
- Handles arbitrary propagation distances
- Natural treatment of evanescent waves

Sampling requirement:
Δx < λz/(2X) where X is the field extent

## 15.4 Fraunhofer Diffraction and Fourier Optics

### Far-Field Approximation

In the Fraunhofer (far-field) regime, we further approximate the Fresnel integral by assuming the observation distance z is so large that:

z >> k(x₀² + y₀²)_max/2

This allows us to move the quadratic phase term outside the integral:

u(x,y,z) = (e^(ikz)/iλz) exp[ik/2z(x² + y²)] × 
           ∫∫ u₀(x₀,y₀) exp[-ik/z(xx₀ + yy₀)] dx₀dy₀

### Fourier Transform Relationship

The integral is now a 2D Fourier transform of the aperture field:

u(x,y,z) = (e^(ikz)/iλz) exp[ik/2z(x² + y²)] × ℱ{u₀(x₀,y₀)}|_{fₓ=x/λz, fᵧ=y/λz}

For a plane wave incident on the aperture (u₀ = A(x₀,y₀) where A is the aperture function):

u(x,y,z) ∝ ℱ{A(x₀,y₀)}

**Key insight**: The far-field diffraction pattern is the Fourier transform of the aperture.

### Examples

1. **Rectangular aperture** A(x₀,y₀) = rect(x₀/a)rect(y₀/b):
   u(x,y) ∝ sinc(ax/λz)sinc(by/λz)

2. **Circular aperture** of radius a:
   u(r,θ) ∝ 2J₁(kar/z)/(kar/z)
   where J₁ is the Bessel function of the first kind.

3. **Double slit** with separation d:
   u(x) ∝ sinc(ax/λz)cos(πdx/λz)

### Angular Spectrum Representation

Any field can be decomposed into plane waves propagating in different directions:

u(x,y,z) = ∫∫ A(kₓ,kᵧ) exp[i(kₓx + kᵧy + kᵣz)] dkₓdkᵧ

where the z-component of the wavevector:
kᵣ = √(k² - kₓ² - kᵧ²)

Two cases arise:

1. **Propagating waves** (kₓ² + kᵧ² < k²):
   - kᵣ is real
   - Plane waves propagate without decay
   - Direction cosines: (α,β,γ) = (kₓ/k, kᵧ/k, kᵣ/k)
   - Physical angles: θₓ = arcsin(kₓ/k), θᵧ = arcsin(kᵧ/k)

2. **Evanescent waves** (kₓ² + kᵧ² > k²):
   - kᵣ = iκ where κ = √(kₓ² + kᵧ² - k²)
   - Exponential decay: exp(-κz)
   - Confined to near-field (z ~ 1/κ ~ λ)
   - Carry sub-wavelength information

The angular spectrum at z = 0:
A(kₓ,kᵧ) = (1/2π)² ∫∫ u(x,y,0) exp[-i(kₓx + kᵧy)] dxdy = ℱ{u(x,y,0)}

This representation provides:
- Complete description of the field
- Natural propagation: multiply by exp(ikᵣz)
- Direct connection to Fourier optics
- Basis for understanding resolution limits

### Connection to Rendering

The Fourier optics framework provides deep insights for rendering:

#### 1. Frequency Analysis of Materials
The BRDF acts as a transfer function in angular frequency space:

- **Spatial BRDF**: ρ(x,ω₀,ωᵢ)
- **Angular spectrum**: ρ̃(k,ω₀,ωᵢ) = ℱ_x{ρ(x,ω₀,ωᵢ)}
- **Bandwidth**: Determines required sampling rate

Mirror: ρ̃ ~ δ(k) (all frequencies)
Diffuse: ρ̃ ~ sinc(k) (low-pass)
Glossy: Intermediate bandwidth

#### 2. Sampling Theory Applications

**Nyquist-Shannon theorem** in rendering context:
- Spatial: Δx < 1/(2f_max) where f_max is highest spatial frequency
- Angular: Δω < π/k_max for BRDF sampling
- Temporal: Δt < 1/(2f_motion) for motion blur

**Practical implications**:
- Texture filtering: mipmap levels based on frequency content
- Shadow map resolution: determined by light frequency
- Importance sampling: concentrate samples where |ρ̃| is large

#### 3. Anti-aliasing as Filtering

Rendering pipeline in frequency domain:

1. **Scene spectrum**: S̃(k) = ℱ{scene geometry/materials}
2. **Sampling**: Multiplication by comb function
3. **Reconstruction**: Convolution with filter kernel
4. **Display**: Band-limited by pixel grid

Optimal anti-aliasing:
- Pre-filter to remove frequencies > Nyquist
- Common filters: Box (sinc), Gaussian (Gaussian), Lanczos (windowed sinc)
- Trade-off: Sharpness vs. aliasing

#### 4. Light Field Analysis

4D light field L(x,y,u,v) has 4D Fourier transform:
L̃(k_x,k_y,k_u,k_v)

Key insights:
- Lambertian surfaces: Energy concentrated at k_u = k_v = 0
- Specular surfaces: Energy along k_x = λk_u, k_y = λk_v
- Depth creates shearing in frequency domain
- Enables optimal sampling strategies

#### 5. Coherent Rendering Effects

Extending rendering equation for coherence:

L(x,ω) = L₀(x,ω) + ∫ ρ(x,ω'→ω)L(x,ω')V(x,x')G(x,x')dx'

where V(x,x') is the mutual coherence function:
V(x,x') = 〈E*(x)E(x')〉 / √(I(x)I(x'))

This enables:
- Laser speckle simulation
- Holographic displays
- Interference in thin films
- Coherent subsurface scattering

## 15.5 Diffraction-Limited Imaging Systems

### Point Spread Function (PSF)

An ideal imaging system maps each object point to a unique image point. However, diffraction limits this ideal behavior. The image of a point source is the **Point Spread Function (PSF)**.

For a circular aperture of diameter D and focal length f, the PSF in the image plane is:

h(r) = |ℱ{P(x,y)}|² = [2J₁(πDr/λf)/(πDr/λf)]²

where P(x,y) is the pupil function (1 inside aperture, 0 outside).

### Airy Disk and Resolution

The PSF for a circular aperture forms the **Airy pattern**:

h(r) = [2J₁(πDr/λf)/(πDr/λf)]²

Characteristics:
- Central bright disk (Airy disk) contains 83.8% of total energy
- First dark ring at r₁ = 1.22λf/D
- First bright ring: 7.2% of energy
- Second bright ring: 2.8% of energy
- Ring radii: r_n ≈ (n + 0.22)λf/D for n ≥ 1

The Airy disk radius (first zero):

r₀ = 1.22λf/D = 1.22λF#

where F# = f/D is the f-number.

**Energy distribution**:
- Within r₀: 83.8%
- Within 2r₀: 91.0%
- Within 3r₀: 93.8%

This concentration of energy in the central disk is why the Airy disk radius serves as a practical measure of resolution.

### Rayleigh Criterion

Two point sources are "just resolved" when the maximum of one Airy disk falls on the first minimum of the other:

θ_min = 1.22λ/D

This angular resolution limit is fundamental to all imaging systems.

### Coherent vs Incoherent Imaging

**Incoherent imaging** (typical for natural light):
- Intensities add: I_total = I₁ + I₂
- Image intensity = |Object|² ⊗ |PSF|²
- Linear in intensity

**Coherent imaging** (laser illumination):
- Fields add: U_total = U₁ + U₂
- Image field = Object ⊗ PSF
- Linear in complex amplitude
- Can exhibit interference effects

### Transfer Functions

**Optical Transfer Function (OTF)** for incoherent imaging:
OTF(f) = ℱ{|PSF|²}

**Modulation Transfer Function (MTF)**:
MTF(f) = |OTF(f)|

For a circular aperture:
MTF(ν) = (2/π)[arccos(ν) - ν√(1-ν²)] for ν ≤ 1
MTF(ν) = 0 for ν > 1

where ν = λf·f_spatial/D is the normalized spatial frequency.

### Implications for Computer Graphics

#### 1. Depth of Field and Diffraction Limits

The circle of confusion (CoC) has two contributions:
- **Geometric**: C_geom = D|z - z_f|/z_f (defocus)
- **Diffraction**: C_diff = 2.44λF# (Airy disk diameter)

Total CoC: C_total = √(C_geom² + C_diff²)

Consequences:
- Minimum CoC at optimal aperture: F# = √(|z - z_f|λ/(2.44z_f))
- Diffraction-limited for F# > 8-11 in visible light
- Hyperfocal distance: H = f²/(F#c) + f, where c includes diffraction

#### 2. Physically-Based Bokeh

Bokeh shape depends on:

**Geometric limit** (F# < 5.6):
- Shape matches aperture geometry
- Sharp edges from aperture blades
- Uniform intensity distribution

**Transition regime** (F# ~ 5.6-11):
- Diffraction softens edges
- Brightness varies: brighter center
- Convolution: Bokeh = Aperture ⊗ Airy

**Diffraction limit** (F# > 11):
- Circular regardless of aperture shape
- Airy pattern dominates
- Rings may be visible in high contrast

Implementation approach:
1. Compute geometric bokeh kernel
2. Convolve with wavelength-dependent Airy function
3. Sum over visible spectrum for color effects

#### 3. Wave-Optical Material Effects

**Glints and Sparkles**:
- Caused by coherent reflection from rough surfaces
- Each microfacet creates diffraction pattern
- Interference between nearby facets
- Statistics: I = |E₁ + E₂ + ...|u00b2 follows speckle statistics

Modeling approach:
- Heightfield h(x,y) with correlation length ξ
- Phase variation: φ = 2kh cosθ
- Speckle size: Δx ~ λz/ξ
- Implement as normal-mapped diffraction

**Iridescence**:
- Thin-film interference
- Structural color from periodic nanostructures
- Wavelength-dependent reflection
- Requires wave-based BRDF models

#### 4. Advanced Camera Models

**Beyond thin lens**:
1. **Wavefront aberrations**: Φ(x,y) = ∑ Z_n(x,y)
   - Zernike polynomials Z_n describe aberrations
   - PSF = |ℱ{P(x,y)exp(ikΦ(x,y))}|²
   - Spatially-varying blur kernels

2. **Chromatic effects**:
   - Longitudinal: focal length f(λ)
   - Lateral: magnification m(λ)
   - PSF varies with wavelength
   - Natural chromatic aberration

3. **Polarization**:
   - Fresnel coefficients depend on polarization
   - Polarizing filters in lens systems
   - Sky models with polarization

4. **Coherence effects**:
   - Partial coherence from extended sources
   - Coherence area: A_c ~ λ²R²/A_s
   - Affects contrast and resolution

## Summary

This chapter established the mathematical foundation for wave optics, transitioning from Maxwell's equations to practical diffraction formulas. Key concepts include:

1. **Helmholtz Equation**: ∇²u + k²u = 0 - the fundamental equation for monochromatic wave propagation
2. **Huygens-Fresnel Principle**: Each wavefront point acts as a secondary source
3. **Fresnel Diffraction**: Near-field with quadratic phase approximation
4. **Fraunhofer Diffraction**: Far-field reduces to Fourier transform
5. **Resolution Limits**: Diffraction fundamentally limits imaging resolution

The wave nature of light introduces phenomena beyond geometric optics:
- Interference and diffraction patterns
- Fundamental resolution limits (Rayleigh criterion)
- Coherence effects in imaging
- Frequency-domain analysis of optical systems

These concepts bridge to computer graphics through:
- Physical camera models with diffraction
- Wave-based material appearance (glints, iridescence)
- Fourier analysis of rendering algorithms
- Connection to volume rendering via Green's functions

## Exercises

### Basic Understanding (3 problems)

**Exercise 15.1**: Helmholtz Equation Solutions
Show that u(r) = (A/r)exp(ikr) is a solution to the Helmholtz equation in spherical coordinates. What physical wave does this represent?

*Hint*: Use the spherical Laplacian: ∇²u = (1/r²)d/dr(r²du/dr) for spherically symmetric functions.

<details>
<summary>Solution</summary>

For u(r) = (A/r)exp(ikr):

du/dr = A[(-1/r²)exp(ikr) + (ik/r)exp(ikr)] = (A/r)exp(ikr)[ik - 1/r]

r²du/dr = A·r·exp(ikr)[ik - 1/r] = A·exp(ikr)[ikr - 1]

d/dr(r²du/dr) = A[ik·exp(ikr)·[ikr - 1] + exp(ikr)·ik]
                = A·exp(ikr)[−k²r + 2ik]

∇²u = (A/r²)exp(ikr)[−k²r + 2ik] = (A/r)exp(ikr)[−k²]

Therefore: ∇²u + k²u = (A/r)exp(ikr)[−k² + k²] = 0 ✓

This represents an outgoing spherical wave from a point source.
</details>

**Exercise 15.2**: Fresnel Number
A plane wave (λ = 500nm) illuminates a circular aperture of radius a = 1mm. At what distance z does the Fresnel number F = a²/λz equal 1? What approximation regime is this?

*Hint*: The Fresnel approximation is valid for F ≳ 1, while Fraunhofer requires F << 1.

<details>
<summary>Solution</summary>

Given: λ = 500 × 10⁻⁹ m, a = 1 × 10⁻³ m

F = a²/λz = 1

Solving for z:
z = a²/λ = (10⁻³)² / (500 × 10⁻⁹) = 10⁻⁶ / (5 × 10⁻⁷) = 2 m

At z = 2m, we're at the transition between Fresnel (near-field) and Fraunhofer (far-field) regimes. For z < 2m, use Fresnel diffraction; for z >> 2m, Fraunhofer approximation is valid.
</details>

**Exercise 15.3**: Airy Disk Size
A camera lens has focal length f = 50mm and aperture diameter D = 25mm (f/2). Calculate the Airy disk radius for green light (λ = 550nm). How does this compare to typical pixel sizes?

*Hint*: The Airy disk radius is r₀ = 1.22λf/D.

<details>
<summary>Solution</summary>

Given: f = 50mm, D = 25mm, λ = 550nm = 550 × 10⁻⁹ m

r₀ = 1.22λf/D = 1.22 × (550 × 10⁻⁹) × (50 × 10⁻³) / (25 × 10⁻³)
   = 1.22 × 550 × 10⁻⁹ × 2
   = 1.342 × 10⁻⁶ m = 1.34 μm

Diameter = 2r₀ = 2.68 μm

Modern camera sensors have pixel sizes of 1-5 μm, so the Airy disk spans approximately 1-3 pixels. This shows that many cameras are near the diffraction limit, especially at small apertures.
</details>

### Advanced Problems (3 problems)

**Exercise 15.4**: Fourier Optics and Rendering
Show that the rendering equation in the Fourier domain becomes a convolution. Start with:
L₀(x,ω₀) = ∫ ρ(x,ω₀,ωᵢ)L(x,ωᵢ)(ω₀·n)dωᵢ

*Hint*: Take the 2D spatial Fourier transform and use the convolution theorem.

<details>
<summary>Solution</summary>

Taking the 2D Fourier transform over x:

ℱ{L₀(x,ω₀)} = ℱ{∫ ρ(x,ω₀,ωᵢ)L(x,ωᵢ)(ω₀·n)dωᵢ}

For spatially-invariant BRDF ρ(x,ω₀,ωᵢ) = ρ(ω₀,ωᵢ):

L̃₀(k,ω₀) = ∫ ρ(ω₀,ωᵢ)ℱ{L(x,ωᵢ)}(ω₀·n)dωᵢ
          = ∫ ρ(ω₀,ωᵢ)L̃(k,ωᵢ)(ω₀·n)dωᵢ

For textured surfaces where ρ varies with x:

L̃₀(k,ω₀) = ∫ [ρ̃(k,ω₀,ωᵢ) ⊗ L̃(k,ωᵢ)](ω₀·n)dωᵢ

This shows that spatial texture variations cause frequency-domain convolution, leading to blur and aliasing if not properly sampled.
</details>

**Exercise 15.5**: Kirchhoff Boundary Conditions
Derive the Kirchhoff diffraction formula from the Helmholtz equation using Green's theorem. Show why the boundary conditions on an opaque screen are problematic.

*Hint*: Use Green's function G = exp(ikr)/r and Green's theorem: ∫∫∫_V (ψ∇²φ - φ∇²ψ)dV = ∫∫_S (ψ∂φ/∂n - φ∂ψ/∂n)dS

<details>
<summary>Solution</summary>

Let u satisfy (∇² + k²)u = 0 and G = exp(ikr)/r satisfy (∇² + k²)G = -4πδ(r).

Applying Green's theorem with ψ = G and φ = u:

∫∫∫_V [G∇²u - u∇²G]dV = ∫∫_S [G∂u/∂n - u∂G/∂n]dS

Since ∇²u = -k²u and ∇²G = -k²G - 4πδ(r-r₀):

-4πu(r₀) = ∫∫_S [G∂u/∂n - u∂G/∂n]dS

u(P) = (1/4π) ∫∫_S [exp(ikr)/r ∂u/∂n - u ∂/∂n(exp(ikr)/r)]dS

Kirchhoff boundary conditions assume:
- On aperture: u = u_incident, ∂u/∂n = ∂u_incident/∂n
- On screen: u = 0, ∂u/∂n = 0

The problem: These conditions are inconsistent at the aperture edge where u must jump from u_incident to 0 discontinuously, violating the wave equation. This is the "Kirchhoff paradox" - the approximation works well in practice despite theoretical inconsistency.
</details>

**Exercise 15.6**: Volume Rendering Connection
Show how the volume rendering equation with scattering reduces to the Huygens-Fresnel principle in the appropriate limit. Consider:
L(x,ω) = ∫ σₛ(x')p(x',ω'→ω)G(x,x')L(x',ω')dx'

*Hint*: Consider a thin scattering layer and the Green's function for the Helmholtz equation.

<details>
<summary>Solution</summary>

For monochromatic light, the Green's function satisfies:
(∇² + k²)G(x,x') = -δ(x-x')

In free space: G(x,x') = exp(ik|x-x'|)/(4π|x-x'|)

For a thin scattering layer at z = 0 with σₛ(x') = σ₀δ(z')A(x',y'):

L(x,y,z) = ∫∫ σ₀A(x',y')p(θ)G(x,x')L₀(x',y')dx'dy'

For forward scattering p(θ) ≈ (1 + cos θ)/2 and incident field L₀:

L(x,y,z) = σ₀/(4π) ∫∫ A(x',y')L₀(x',y') × 
           [exp(ikr)/r][(1 + cos χ)/2]dx'dy'

Setting σ₀/(4π) = 1/(iλ) recovers the Huygens-Fresnel formula:

u(P) = (1/iλ) ∫∫ u₀(Q)[exp(ikr)/r]K(χ)dS

This shows that the Huygens-Fresnel principle emerges from volume scattering in the limit of a thin layer with appropriate scattering properties.
</details>

### Challenge Problems (2 problems)

**Exercise 15.7**: Computational Complexity
Compare the computational complexity of three methods for computing Fresnel diffraction patterns:
1. Direct numerical integration
2. FFT-based convolution  
3. Angular spectrum method

For an N×N sampling grid, derive the complexity and discuss trade-offs.

*Hint*: Consider both computational cost and memory requirements.

<details>
<summary>Solution</summary>

1. **Direct Integration**: 
   - For each output point (N² total), integrate over N² input points
   - Complexity: O(N⁴)
   - Memory: O(N²)
   - Accurate but prohibitively slow for large N

2. **FFT Convolution**:
   - Fresnel integral as convolution with chirp function
   - Steps: FFT input (O(N²log N)), multiply (O(N²)), inverse FFT (O(N²log N))
   - Complexity: O(N²log N)
   - Memory: O(N²)
   - Requires careful sampling to avoid aliasing

3. **Angular Spectrum**:
   - Propagate in Fourier domain: H(fx,fy) = exp[ikz√(1-(λfx)²-(λfy)²)]
   - Steps: FFT (O(N²log N)), multiply by H (O(N²)), inverse FFT (O(N²log N))
   - Complexity: O(N²log N)
   - Memory: O(N²)
   - Most efficient, handles evanescent waves correctly

Trade-offs:
- Direct: Most flexible (arbitrary geometries) but slowest
- FFT convolution: Fast but can have aliasing issues with quadratic phase
- Angular spectrum: Fastest and most accurate for planar geometries

For typical N = 1024: Direct takes ~10¹² operations vs ~10⁷ for FFT methods.
</details>

**Exercise 15.8**: Unified Framework
Develop a unified mathematical framework that encompasses both geometric ray tracing and wave optics. Show how to transition smoothly between regimes based on the Fresnel number.

*Hint*: Consider the stationary phase approximation and the relationship between rays and wavefronts.

<details>
<summary>Solution</summary>

**Unified Framework**: Wigner Distribution Function (WDF)

The WDF W(x,k) combines position and momentum (direction) information:

W(x,k,z) = ∫ u*(x - ξ/2,z)u(x + ξ/2,z)exp(-ik·ξ)dξ

Properties:
- Marginals give intensity and angular spectrum: ∫W dk = |u(x)|², ∫W dx = |ũ(k)|²
- Evolution: ∂W/∂z + (k/k₀)·∇W = 0 (free space)
- Reduces to ray density in geometric limit

**Regime Transition**:

Define normalized scale parameter: ε = λz/a² = 1/F

1. **Geometric Optics** (ε → 0, F → ∞):
   - WDF → ray phase space density
   - W(x,k) = ∑ᵢ δ(x - xᵢ(z))δ(k - kᵢ)
   - Ray tracing valid

2. **Fresnel Regime** (ε ~ 1, F ~ 1):
   - Quadratic phase approximation
   - W spreads in phase space
   - Use Fresnel integrals

3. **Fraunhofer Regime** (ε >> 1, F << 1):
   - Position-momentum uncertainty maximized
   - W(x,k) ≈ W₀(x)W̃₀(k)
   - Fourier optics applies

**Smooth Transition**:

Propagation operator: P(z) = exp[iz(k²/2k₀ + Φ(x,k,ε))]

where Φ interpolates:
- Φ → 0 as ε → 0 (geometric)
- Φ → higher-order terms as ε increases

This framework unifies:
- Ray tracing (ε → 0)
- Gaussian beam propagation (intermediate ε)
- Full wave optics (arbitrary ε)

The WDF provides a phase-space representation that smoothly transitions between particle-like rays and wave-like diffraction, controlled by the Fresnel number.
</details>

## Common Pitfalls and Errors

### Approximation Validity
1. **Scalar Approximation**: Invalid for:
   - Strong focusing (NA > 0.6)
   - Near-field of subwavelength features
   - Polarization-dependent effects

2. **Fresnel vs Fraunhofer**: 
   - Fresnel: F ≳ 1 (near-field)
   - Fraunhofer: F << 1 (far-field)
   - Transition region needs careful handling

3. **Paraxial Approximation**: Breaks down for:
   - Large angles (> 15-20°)
   - Wide-aperture systems
   - Off-axis points

### Numerical Considerations

1. **Sampling Requirements**:
   - Quadratic phase in Fresnel integral requires dense sampling
   - Nyquist criterion: Δx < λz/(2X) where X is field extent
   - Aliasing causes artificial fringes

2. **FFT Artifacts**:
   - Periodic boundary conditions create wraparound
   - Zero-padding needed for accurate convolution
   - Windowing functions reduce edge effects

3. **Phase Unwrapping**:
   - Computed phase limited to [-π, π]
   - Unwrapping algorithms needed for continuous phase
   - Sensitive to noise and undersampling

4. **Numerical Precision**:
   - Large k values cause precision loss
   - exp(ikr) oscillates rapidly for large r
   - Use differential propagation for long distances

## Best Practice Checklist

### Design Review Points

✓ **Physical Validity**
- [ ] Wavelength range specified
- [ ] Coherence properties defined
- [ ] Polarization effects considered
- [ ] Material dispersion included if needed

✓ **Approximation Choice**
- [ ] Fresnel number calculated
- [ ] Appropriate regime selected
- [ ] Error bounds estimated
- [ ] Edge cases identified

✓ **Numerical Implementation**
- [ ] Sampling rate meets Nyquist criterion
- [ ] FFT size includes padding
- [ ] Boundary conditions properly handled
- [ ] Precision adequate for phase calculations

✓ **Performance Optimization**
- [ ] Algorithm complexity analyzed
- [ ] Memory requirements estimated
- [ ] GPU acceleration considered
- [ ] Multi-scale methods evaluated

✓ **Validation Strategy**
- [ ] Analytical test cases verified
- [ ] Energy conservation checked
- [ ] Reciprocity maintained
- [ ] Comparison with geometric limit

✓ **Graphics Integration**
- [ ] Rendering pipeline compatibility
- [ ] Real-time constraints evaluated
- [ ] Level-of-detail strategy defined
- [ ] Perceptual importance assessed
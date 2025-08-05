# Chapter 1: Geometric Optics and Rendering Fundamentals

This chapter establishes the mathematical foundation for computer graphics through the lens of geometric optics. We develop the rendering equation as our central framework, introduce key radiometric concepts, and establish the path integral formulation that will unify all subsequent rendering techniques. By treating light transport as a high-dimensional integration problem, we set the stage for understanding point-based, image-based, and neural rendering methods as different approaches to solving the same fundamental equation.

## Learning Objectives

After completing this chapter, you will be able to:
1. Derive the rendering equation from first principles using energy conservation
2. Transform between different coordinate systems while preserving radiometric quantities
3. Analyze BRDF properties and verify physical plausibility
4. Apply Monte Carlo methods to estimate high-dimensional integrals with known error bounds
5. Express light transport as a path integral and connect it to volume rendering

## 1.1 Ray Tracing Basics and the Rendering Equation

### Light as Rays in Geometric Optics

In geometric optics, we model light propagation using rays—infinitesimal beams that travel in straight lines through homogeneous media. This approximation holds when wavelength λ << feature size, allowing us to ignore diffraction and interference. A ray is parameterized as:

**r**(t) = **o** + t**d**

where **o** ∈ ℝ³ is the origin, **d** ∈ ℝ³ is the direction (||**d**|| = 1), and t ≥ 0 is the parameter along the ray.

The ray equation emerges from the eikonal equation ∇S = n**k̂** in the limit λ → 0, where S is the phase and n is the refractive index. In inhomogeneous media, rays follow curved paths satisfying:

d/ds(n d**r**/ds) = ∇n

This reduces to straight lines when n is constant.

**Connection to Wave Optics**: The geometric optics approximation emerges from the WKB (Wentzel-Kramers-Brillouin) approximation of the wave equation. When we substitute ψ = A exp(ikS) into the Helmholtz equation and take k → ∞:

∇²ψ + k²n²ψ = 0 → (∇S)² = n² (eikonal equation)

The surfaces of constant phase S = const are wavefronts, and rays are orthogonal trajectories to these wavefronts. This connection becomes crucial when we extend to wave optics in later chapters.

**Ray Optics Validity**: The geometric optics approximation breaks down when:
1. Feature size ~ wavelength (diffraction becomes significant)
2. Near caustics where ray density → ∞
3. In the presence of sharp edges or discontinuities
4. For coherent phenomena requiring phase information

**Fermat's Principle**: Rays follow paths of stationary optical path length:

δ∫ n(**r**) ds = 0

This variational principle unifies ray behavior: straight lines in homogeneous media, Snell's law at interfaces, and curved paths in gradient-index media. It also connects to the principle of least action in physics and the geodesic equation in differential geometry.

### Radiometric Quantities

Before deriving the rendering equation, we must establish our radiometric framework. These quantities form a hierarchy, each building upon the previous:

**Radiant energy** Q measures total electromagnetic energy:
Q [J]

**Radiant flux (power)** Φ measures energy per unit time:
Φ = dQ/dt [W]

**Radiant intensity** I measures flux per unit solid angle from a point source:
I = dΦ/dω [W/sr]

**Irradiance** E measures flux incident per unit area:
E = dΦ/dA [W/m²]

**Radiant exitance** M measures flux leaving per unit area:
M = dΦ/dA [W/m²]

**Radiance** L measures flux per unit area per unit solid angle:
L = d²Φ/(dA cos θ dω) = d²Φ/(dA⊥ dω) [W/(m²·sr)]

where dA⊥ = dA cos θ is the projected area perpendicular to the ray direction.

Radiance is the fundamental quantity in rendering because:
1. It remains constant along rays in vacuum (radiance invariance theorem)
2. It's what cameras and eyes measure
3. All other radiometric quantities can be derived from it

**Photometric vs Radiometric Quantities**: While we focus on radiometric quantities (physical energy), rendering for human perception often uses photometric quantities:
- Luminous flux [lm] = Radiant flux [W] × luminous efficacy
- Luminance [cd/m²] = Radiance × photopic response V(λ)
- The CIE luminosity function V(λ) peaks at 555 nm (green)

**Spectral Radiance**: In reality, radiance varies with wavelength:
L(x, ω, λ) [W/(m²·sr·nm)]

For rendering, we typically use:
- RGB approximation: 3 samples of the spectrum
- Spectral rendering: N wavelength samples (typ. 10-100)
- Hero wavelengths: Stochastic sampling of spectrum

**Coherent vs Incoherent Addition**: Radiometry assumes incoherent light—intensities add directly. For coherent sources (lasers), we must track phase and add complex amplitudes:
I_total = |E₁ + E₂|² ≠ |E₁|² + |E₂|² (in general)

### The Rendering Equation

The rendering equation, introduced by Kajiya (1986), describes the equilibrium distribution of light in a scene. It emerges from power balance: at any surface point, outgoing power equals emitted plus reflected power.

At surface point **x** with normal **n**, the outgoing radiance L_o in direction **ω**_o satisfies:

L_o(**x**, **ω**_o) = L_e(**x**, **ω**_o) + ∫_Ω f_r(**x**, **ω**_i, **ω**_o) L_i(**x**, **ω**_i) (**ω**_i · **n**) dω_i

where:
- L_e(**x**, **ω**_o) is emitted radiance from **x** in direction **ω**_o
- f_r(**x**, **ω**_i, **ω**_o) is the BRDF [sr⁻¹]
- L_i(**x**, **ω**_i) is incident radiance at **x** from direction **ω**_i
- Ω is the hemisphere above **x** (where **ω** · **n** > 0)
- (**ω**_i · **n**) = cos θ_i accounts for projected area

The integral represents the scattering integral—summing contributions from all incident directions, weighted by the BRDF and cosine foreshortening.

### Energy Conservation and the Measurement Equation

Energy conservation constrains the BRDF. The directional-hemispherical reflectance (albedo) must satisfy:

ρ(**ω**_o) = ∫_Ω f_r(**x**, **ω**_i, **ω**_o) cos θ_i dω_i ≤ 1 for all **ω**_o

Equality holds for lossless surfaces. The white furnace test verifies energy conservation: in a uniformly lit environment (L_i = L_0), a closed surface should neither gain nor lose energy.

**Detailed Energy Balance**: For a surface element dA, conservation requires:

∫_Ω L_o(**x**, **ω**) cos θ dω dA = L_e dA + ∫_Ω L_i(**x**, **ω**) cos θ dω dA

In a closed system at thermal equilibrium, Kirchhoff's law relates emissivity to absorptivity:
ε(λ, θ) = α(λ, θ) = 1 - ρ(λ, θ)

**The Measurement Equation**: The measurement equation connects scene radiance to sensor response:

I_j = ∫_A ∫_Ω W_j(**x**, **ω**) L(**x**, **ω**) cos θ dω dA

where W_j(**x**, **ω**) is the importance (sensitivity) function for pixel j. This duality between radiance and importance enables bidirectional algorithms.

**Importance Transport**: Importance satisfies an adjoint equation:

W(**x**, **ω**) = W_e(**x**, **ω**) + ∫_Ω f_r(**x**, **ω**, **ω**') W(**x**, **ω**') cos θ' dω'

This symmetry leads to:
- Bidirectional path tracing
- Photon mapping (forward light, backward importance)
- Adjoint methods for gradient computation

For a pinhole camera with pixel j subtending solid angle Ω_j from the pinhole:

I_j = ∫_{Ω_j} L(**x**_lens, **ω**) cos⁴ θ dω

The cos⁴ θ term accounts for:
- cos θ: projected lens area
- cos³ θ: inverse square falloff and pixel foreshortening

**Finite Aperture Cameras**: For realistic cameras with aperture A_lens:

I_j = (1/A_lens) ∫_{A_lens} ∫_{A_pixel} L(**x**_lens → **x**_pixel) G(**x**_lens ↔ **x**_pixel) dA_pixel dA_lens

This leads to depth of field effects and requires careful sampling strategies.

### Operator Form and Neumann Series

The rendering equation admits an elegant operator formulation. Define the transport operator 𝒯:

(𝒯L)(**x**, **ω**) = ∫_Ω f_r(**x**, **ω**', **ω**) L(**x**, **ω**') (**ω**' · **n**) dω'

Then the rendering equation becomes:

L = L_e + 𝒯L

This is a Fredholm equation of the second kind. The solution via Neumann series:

L = (I - 𝒯)⁻¹L_e = ∑_{k=0}^∞ 𝒯^k L_e = L_e + 𝒯L_e + 𝒯²L_e + ...

Each term has physical meaning:
- L_e: Direct illumination (emission only)
- 𝒯L_e: Single-bounce illumination
- 𝒯²L_e: Two-bounce illumination
- 𝒯^k L_e: k-bounce illumination

The series converges when ||𝒯|| < 1, which occurs when max albedo < 1. This decomposition naturally leads to path tracing algorithms that sample paths of increasing length.

### Three-Point Form and Geometric Coupling

The rendering equation can be rewritten in three-point form, making the geometric coupling explicit:

L(**x** → **x**') = L_e(**x** → **x**') + ∫_M f_r(**x**'' → **x** → **x**') L(**x**'' → **x**) G(**x**'' ↔ **x**) dA(**x**'')

where the geometry factor is:

G(**x** ↔ **x**') = V(**x** ↔ **x**') cos θ cos θ' / ||**x** - **x**'||²

with:
- V(**x** ↔ **x**'): binary visibility function (1 if mutually visible, 0 otherwise)
- cos θ, cos θ': angles between surface normals and connecting line
- ||**x** - **x**'||²: squared distance for inverse square falloff

This form emphasizes that light transport couples all surface points, leading to the path integral formulation.

**Visibility Complexity**: The visibility function V(**x** ↔ **x**') makes the rendering equation non-linear and non-local:
- Discontinuous: Creates hard shadows and occlusion boundaries
- Expensive to evaluate: Requires ray-scene intersection
- Couples all geometry: Changes anywhere affect visibility everywhere

**Kernel Properties**: The transport kernel K(**x**'' → **x**) = f_r G V has important properties:
- Singular along **x** = **x**'' (requires careful regularization)
- Discontinuous at occlusion boundaries
- Satisfies reciprocity: K(**x** → **x**') = K(**x**' → **x**)

**Connection to Heat Equation**: Without visibility, the rendering equation resembles the heat equation with a non-local kernel. This analogy helps understand:
- Smoothing properties of multiple scattering
- Diffusion approximation for optically thick media
- Finite element and multigrid solution methods

## 1.2 Coordinate Systems and Transformations

### World, Camera, and Object Spaces

Rendering pipelines involve a hierarchy of coordinate systems, each optimized for specific calculations:

1. **Object space (Model space)**: Geometry defined in canonical form
   - Origin typically at object center or base
   - Axes aligned with natural symmetries
   - Simplifies modeling and animation

2. **World space**: Unified scene coordinates
   - All objects transformed to common frame
   - Lighting and physics calculations
   - Ray-object intersections

3. **Camera space (View space)**: Observer-centric coordinates
   - Origin at eye point
   - -z axis along view direction (OpenGL convention)
   - +z into screen (DirectX convention)
   - Simplifies projection and culling

4. **Clip space**: Post-projection homogeneous coordinates
   - 4D coordinates before perspective divide
   - View frustum becomes [-1,1]³ cube (NDC)

5. **Screen space (Raster space)**: Final 2D image coordinates
   - Integer pixel coordinates
   - Origin at top-left or bottom-left

### Homogeneous Coordinates and Transformations

Homogeneous coordinates unify translation and linear transformations. A 3D point **p** = (x, y, z) becomes **p̃** = (x, y, z, 1), while vectors use **ṽ** = (x, y, z, 0).

The general affine transformation matrix:

**M** = [**A** **t**]
      [**0** 1  ]

where **A** is 3×3 linear part and **t** is translation. Common transformations:

**Translation by (tx, ty, tz):**
[1  0  0  tx]
[0  1  0  ty]
[0  0  1  tz]
[0  0  0  1 ]

**Rotation around axis **a** by angle θ:**
**R** = cos θ **I** + (1 - cos θ) **a****a**^T + sin θ [**a**]_×

where [**a**]_× is the skew-symmetric cross-product matrix.

**Scale by (sx, sy, sz):**
[sx 0  0  0]
[0  sy 0  0]
[0  0  sz 0]
[0  0  0  1]

### Normal and Tangent Transformations

Normals must transform to remain perpendicular to surfaces. Given transformation **M** for points:

**n**' = (**M**^{-T})^{3×3} **n** / ||(**M**^{-T})^{3×3} **n**||

Proof: For tangent **t** on surface, **n** · **t** = 0. After transformation:
**n**' · **t**' = (**M**^{-T}**n**) · (**M****t**) = **n**^T **M**^{-1} **M** **t** = **n** · **t** = 0

For orthonormal tangent frames {**t**, **b**, **n**}:
- Forward: transform **t** and **b**, then **n** = **t** × **b**
- Or transform **n** as above, then reconstruct frame

**Area and Volume Elements**: Under transformation **M**, differential elements scale as:
- Length: dl' = ||**M****v**|| dl (for direction **v**)
- Area: dA' = |det(**M**)| ||(**M**^{-T})**n**|| dA
- Volume: dV' = |det(**M**)| dV

**Non-uniform Scaling Issues**: Non-uniform scaling breaks isotropy:
- Spheres → ellipsoids
- Isotropic BRDFs → anisotropic BRDFs
- Care needed for physically-based materials

**Handedness Preservation**: When det(**M**) < 0, the transformation flips orientation:
- Right-handed → left-handed coordinate system
- Normal directions must be flipped
- Critical for consistent front/back face determination

### Spherical and Solid Angle Parameterizations

Spherical coordinates provide natural parameterization for directions:

**ω** = (sin θ cos φ, sin θ sin φ, cos θ)

where:
- θ ∈ [0, π]: polar angle from +z axis
- φ ∈ [0, 2π]: azimuthal angle from +x axis

The Jacobian gives differential solid angle:

dω = |∂(ω_x, ω_y)/∂(θ, φ)| dθ dφ = sin θ dθ dφ

Total solid angle of hemisphere: ∫_Ω dω = 2π

Alternative parameterizations useful for sampling:

**Concentric disk mapping** (Shirley-Chiu):
(u, v) ∈ [-1, 1]² → (r, φ) → (x, y) on unit disk

**Octahedral mapping**:
Unit sphere → octahedron → unit square
Preserves area better than spherical coordinates

### Change of Variables in Integrals

The general change of variables formula for integrals:

∫_Ω f(**x**) d**x** = ∫_Ω' f(**x**(**u**)) |det(∂**x**/∂**u**)| d**u**

Critical for rendering:

**Solid angle to area**:
∫_Ω L(**x**, **ω**) cos θ dω = ∫_A L(**x**, **ω**(**x**')) G(**x** ↔ **x**') dA'

where G(**x** ↔ **x**') = V(**x** ↔ **x**') cos θ cos θ' / ||**x** - **x**'||²

**Hemisphere to disk** (for cosine-weighted sampling):
Map (θ, φ) → (r, φ) where r = sin θ
Then p(r, φ) = p(θ, φ) |∂(θ, φ)/∂(r, φ)| = p(θ, φ) / cos θ

**Measure Theory Foundation**: The change of variables formula has measure-theoretic underpinnings:
- Pushforward measure: μ'(A) = μ(f^{-1}(A))
- Radon-Nikodym derivative gives the Jacobian
- Critical for understanding Monte Carlo convergence

### Projective Transformations and Perspective

The perspective projection matrix maps view frustum to clip space:

**P** = [n/r   0     0          0     ]
       [0     n/t   0          0     ]
       [0     0     -(f+n)/(f-n)  -2fn/(f-n)]
       [0     0     -1         0     ]

where n, f are near/far planes, r, t are right/top at near plane.

After perspective divide by w:
- x_ndc = x_clip / w_clip ∈ [-1, 1]
- y_ndc = y_clip / w_clip ∈ [-1, 1]
- z_ndc = z_clip / w_clip ∈ [-1, 1]

Important properties:
- Lines remain lines (except through eye)
- Planes remain planes
- Depth precision is non-linear (more near than far)

### Barycentric Coordinates and Interpolation

For triangle with vertices **v**₀, **v**₁, **v**₂, barycentric coordinates (u, v, w) satisfy:

**p** = u**v**₀ + v**v**₁ + w**v**₂

with constraint u + v + w = 1. Computation via areas:

u = Area(**p**, **v**₁, **v**₂) / Area(**v**₀, **v**₁, **v**₂)

Properties:
- u, v, w ∈ [0, 1] iff **p** inside triangle
- Linear interpolation: f(**p**) = uf₀ + vf₁ + wf₂
- Perspective-correct interpolation requires 1/z correction

For perspective-correct attribute interpolation:
1. Interpolate a/z, b/z, c/z and 1/z in screen space
2. Recover attributes: a = (a/z)/(1/z)

### Differential Geometry and Local Frames

At each surface point, we construct a local frame for shading calculations:

**Tangent space basis**:
- **n**: surface normal (∂**p**/∂u × ∂**p**/∂v normalized)
- **t**: tangent (often ∂**p**/∂u normalized)
- **b**: bitangent (**n** × **t**)

**First Fundamental Form**: The metric tensor describes local surface geometry:
**I** = [E F]
      [F G]

where:
- E = ∂**p**/∂u · ∂**p**/∂u
- F = ∂**p**/∂u · ∂**p**/∂v
- G = ∂**p**/∂v · ∂**p**/∂v

Arc length: ds² = E du² + 2F du dv + G dv²
Area element: dA = √(EG - F²) du dv

**Second Fundamental Form**: Describes surface curvature:
**II** = [e f]
       [f g]

where e = **n** · ∂²**p**/∂u², etc.

Principal curvatures κ₁, κ₂ are eigenvalues of **II****I**^{-1}
- Mean curvature: H = (κ₁ + κ₂)/2
- Gaussian curvature: K = κ₁κ₂

**Transformation to/from world space**:
[**t**_world]   [t_x t_y t_z] [**t**_local]
[**b**_world] = [b_x b_y b_z] [**b**_local]
[**n**_world]   [n_x n_y n_z] [**n**_local]

This orthonormal matrix can be inverted by transpose.

**Anisotropic BRDF parameterization**:
Many BRDFs depend on angle relative to tangent:
- φ_h: azimuthal angle of half-vector in tangent space
- Enables modeling of brushed metals, fabrics, hair

**Parallel Transport**: When tracing rays on surfaces, tangent frames must be parallel transported:
- Maintains orientation consistency
- Preserves anisotropic appearance
- Related to geometric phase in optics

## 1.3 BRDF, BSDF, and BSSRDF

### Bidirectional Reflectance Distribution Function (BRDF)

The BRDF f_r quantifies the differential relationship between incident irradiance and reflected radiance:

f_r(**x**, **ω**_i, **ω**_o) = dL_o(**x**, **ω**_o) / dE_i(**x**, **ω**_i) = dL_o(**x**, **ω**_o) / (L_i(**x**, **ω**_i) cos θ_i dω_i) [sr⁻¹]

Physically, it represents the probability density (after normalization) that a photon from direction **ω**_i scatters into direction **ω**_o.

The BRDF can be decomposed into components:
f_r = f_d + f_s + f_g + ...

where f_d is diffuse, f_s is specular, f_g is glossy, etc. This decomposition aids importance sampling.

### Fundamental BRDF Properties

**Helmholtz Reciprocity:**
f_r(**x**, **ω**_i, **ω**_o) = f_r(**x**, **ω**_o, **ω**_i)

This follows from time-reversal symmetry of Maxwell's equations and the principle of detailed balance. It enables bidirectional path tracing and photon mapping.

**Energy Conservation:**
The directional-hemispherical reflectance must satisfy:

ρ(**ω**_i) = ∫_Ω f_r(**x**, **ω**_i, **ω**_o) cos θ_o dω_o ≤ 1 for all **ω**_i

For energy-conserving BRDFs, equality holds when absorption is zero. The hemispherical-hemispherical reflectance:

ρ_hh = (1/π) ∫_Ω ∫_Ω f_r(**x**, **ω**_i, **ω**_o) cos θ_i cos θ_o dω_i dω_o ≤ 1

**Non-negativity:**
f_r(**x**, **ω**_i, **ω**_o) ≥ 0

Negative values would imply energy absorption dependent on outgoing direction, violating causality.

**Measurability and Integrability:**
For Monte Carlo integration convergence:
f_r ∈ L²(Ω × Ω) (square-integrable)

### Classical BRDF Models

**Lambertian (Perfectly Diffuse):**
f_r = ρ_d/π

where ρ_d ∈ [0, 1] is the diffuse albedo. Energy-conserving by construction.

**Phong Model:**
f_r = (ρ_d/π) + ρ_s (n+2)/(2π) (**r** · **ω**_o)^n

where **r** = 2(**n** · **ω**_i)**n** - **ω**_i is the reflection direction. Not reciprocal!

**Blinn-Phong (Reciprocal):**
f_r = (ρ_d/π) + ρ_s (n+2)/(8π) (**n** · **h**)^n / max(cos θ_i, cos θ_o)

where **h** = (**ω**_i + **ω**_o)/||**ω**_i + **ω**_o|| is the half-vector.

**Cook-Torrance Microfacet Model:**
f_r = (ρ_d/π) + D(**h**)G(**ω**_i, **ω**_o)F(**ω**_i, **h**) / (4 cos θ_i cos θ_o)

where:
- D(**h**): Normal distribution function (e.g., GGX)
- G(**ω**_i, **ω**_o): Geometric attenuation (masking/shadowing)
- F(**ω**_i, **h**): Fresnel reflectance

### Extension to BSDF

The Bidirectional Scattering Distribution Function (BSDF) unifies reflection and transmission:

f_s(**x**, **ω**_i, **ω**_o) = {
  f_r(**x**, **ω**_i, **ω**_o) if **ω**_i · **n** and **ω**_o · **n** have same sign
  f_t(**x**, **ω**_i, **ω**_o) if **ω**_i · **n** and **ω**_o · **n** have opposite sign
}

For dielectric interfaces (e.g., glass), Snell's law governs refraction:
n_i sin θ_i = n_o sin θ_o

The Fresnel equations determine reflection/transmission probabilities:
F_r = ((n_i cos θ_i - n_o cos θ_o)/(n_i cos θ_i + n_o cos θ_o))² (s-polarized)

**Generalized Reciprocity for BTDF:**
Due to radiance compression/expansion across interfaces:

n_i² f_t(**x**, **ω**_i, **ω**_o) = n_o² f_t(**x**, **ω**_o, **ω**_i)

This accounts for the n² factor in radiance L/n² being invariant.

### BSSRDF for Subsurface Scattering

The Bidirectional Scattering Surface Reflectance Distribution Function generalizes the BRDF to non-local transport:

S(**x**_i, **ω**_i, **x**_o, **ω**_o) = dL_o(**x**_o, **ω**_o) / dΦ_i(**x**_i, **ω**_i) [m⁻²sr⁻¹]

Key differences from BRDF:
- Couples different surface points
- Units include inverse area
- No longer a pure material property (depends on geometry)

The rendering equation with BSSRDF:

L_o(**x**_o, **ω**_o) = L_e(**x**_o, **ω**_o) + ∫_A ∫_Ω S(**x**_i, **ω**_i, **x**_o, **ω**_o) L_i(**x**_i, **ω**_i) cos θ_i dω_i dA_i

**Diffusion Approximation:**
For highly scattering media, the BSSRDF can be approximated:

S(**x**_i, **ω**_i, **x**_o, **ω**_o) ≈ (1/π)F_t(**ω**_i)R(||**x**_i - **x**_o||)F_t(**ω**_o)

where R(r) is the diffusion profile and F_t is the Fresnel transmittance.

### Mathematical Constraints and Physical Plausibility

A physically valid BRDF must satisfy:

1. **Reciprocity**: f_r(**x**, **ω**_i, **ω**_o) = f_r(**x**, **ω**_o, **ω**_i)
   - Test: Render scene with swapped lights and cameras

2. **Energy Conservation**: ∀**ω**_i: ∫_Ω f_r(**x**, **ω**_i, **ω**_o) cos θ_o dω_o ≤ 1
   - Test: White furnace test (uniform illumination)

3. **Non-negativity**: f_r(**x**, **ω**_i, **ω**_o) ≥ 0
   - Violations cause energy absorption anomalies

4. **Smoothness**: f_r should be C⁰ continuous (C¹ preferred)
   - Discontinuities cause sampling difficulties

5. **Fresnel Behavior**: f_r → 1 as θ → π/2 for smooth surfaces
   - All surfaces become mirrors at grazing angles

### Anisotropic BRDFs

For materials with directional structure (brushed metal, fabric, hair), the BRDF depends on the azimuthal angle:

f_r(**x**, **ω**_i, **ω**_o, φ) 

where φ is the angle between the half-vector projection and tangent direction.

**Ward Anisotropic Model:**
f_r = (ρ_d/π) + ρ_s exp(-tan²θ_h(cos²φ/α_x² + sin²φ/α_y²)) / (4π α_x α_y √(cos θ_i cos θ_o))

where α_x, α_y control anisotropic roughness.

### Spatially Varying BRDFs (SVBRDFs)

Real materials exhibit spatial variation:
f_r(**x**, **ω**_i, **ω**_o) = f_r(u, v, **ω**_i, **ω**_o)

where (u, v) are texture coordinates. This enables:
- Texture mapping of material properties
- Measured BRDF data (BTF - Bidirectional Texture Function)
- Procedural material variation

## 1.4 Monte Carlo Integration in Rendering

### Expected Value and Variance

Monte Carlo integration estimates integrals using random sampling:

I = ∫_Ω f(**x**) d**x** ≈ (1/N) ∑_{i=1}^N f(**X**_i)/p(**X**_i)

where **X**_i ~ p(**x**) are samples from probability density p.

The estimator is unbiased: E[Î] = I

The variance is:
Var[Î] = (1/N) ∫_Ω (f(**x**)/p(**x**) - I)² p(**x**) d**x**

### Importance Sampling

Optimal sampling minimizes variance by matching p to |f|:

p*(**x**) = |f(**x**)| / ∫_Ω |f(**x**)| d**x**

For the rendering equation, good sampling strategies include:
- BRDF sampling: p(**ω**) ∝ f_r(**ω**_i, **ω**_o)
- Light sampling: p(**ω**) ∝ L_e
- Cosine sampling: p(**ω**) ∝ cos θ

### Multiple Importance Sampling (MIS)

When multiple sampling strategies are available, MIS combines them optimally:

Î = ∑_{i=1}^{n_f} w_f(**X**_{f,i}) f(**X**_{f,i})/p_f(**X**_{f,i}) + ∑_{j=1}^{n_g} w_g(**X**_{g,j}) f(**X**_{g,j})/p_g(**X**_{g,j})

The balance heuristic provides good weights:
w_f(**x**) = n_f p_f(**x**) / (n_f p_f(**x**) + n_g p_g(**x**))

### Russian Roulette

To create unbiased estimators with finite computation, Russian roulette randomly terminates paths:

L'_i = {
  L_i/q  with probability q
  0      with probability 1-q
}

This maintains E[L'_i] = L_i while bounding computation.

### Convergence Rates and Error Bounds

Monte Carlo convergence follows the Central Limit Theorem:

P(|Î - I| ≤ ε) ≈ 2Φ(ε√N/σ) - 1

where Φ is the normal CDF and σ² is the variance. The error decreases as O(1/√N), independent of dimension—crucial for high-dimensional light transport.

## 1.5 Path Integral Formulation

### Light Transport as Path Integration

We can reformulate the rendering equation as an integral over all possible light paths. A path of length k is:

**x̄** = **x**₀**x**₁...**x**_k

where **x**₀ is on a light source and **x**_k is on the camera sensor.

### Path Space and Measure

The path space Ω̄_k consists of all valid paths of length k. The measure for a path is:

dμ(**x̄**) = dA(**x**₀) ∏_{i=1}^{k} dA(**x**_i)

The contribution of a path is:

f(**x̄**) = L_e(**x**₀ → **x**₁) (∏_{i=1}^{k-1} f_s(**x**_{i-1} → **x**_i → **x**_{i+1}) G(**x**_i ↔ **x**_{i+1})) W(**x**_{k-1} → **x**_k)

where G is the geometry factor:

G(**x** ↔ **x**') = V(**x** ↔ **x**') cos θ cos θ' / ||**x** - **x**'||²

### Connection to Feynman Path Integrals

The path integral formulation resembles Feynman's approach to quantum mechanics:

I = ∑_{k=2}^∞ ∫_{Ω̄_k} f(**x̄**) dμ(**x̄**)

This infinite sum over all path lengths captures global illumination. Each term represents paths with k-1 bounces.

### Recursive Formulation and Neumann Series

The value at a point satisfies the recursive relation:

L(**x**, **ω**) = L_e(**x**, **ω**) + ∫_M f_s(**y** → **x** → **ω**) L(**y**, **x** - **y**) G(**y** ↔ **x**) dA(**y**)

This leads to the Neumann series solution:

L = L^{(0)} + L^{(1)} + L^{(2)} + ...

where L^{(k)} represents k-bounce illumination.

### Volume Rendering Equation Preview

The path integral naturally extends to participating media. For a volume with absorption σ_a, scattering σ_s, and phase function p:

L(**x**, **ω**) = ∫₀^∞ T(0,t) [σ_a L_e + σ_s ∫_{S²} p(**ω**', **ω**) L(**x**+t**ω**, **ω**') dω'] dt + T(0,∞) L_∞

where T(s,t) = exp(-∫_s^t σ_t(**x**+u**ω**) du) is transmittance.

This unified formulation will connect all rendering methods in subsequent chapters.

## Chapter Summary

This chapter established the mathematical foundation for computer graphics through geometric optics:

1. **The rendering equation** L_o = L_e + ∫ f_r L_i cos θ dω governs light transport
2. **Coordinate transformations** preserve radiometric quantities when properly applied
3. **BRDFs** must satisfy reciprocity, energy conservation, and non-negativity
4. **Monte Carlo methods** solve high-dimensional integrals with O(1/√N) convergence
5. **Path integrals** unify light transport as integration over all possible paths

These concepts form the basis for all rendering algorithms. The path integral formulation particularly enables our unified treatment of point-based, image-based, and neural rendering methods as different approaches to the same fundamental equation.

## Exercises

### Exercise 1.1: Radiance Along a Ray
Prove that radiance remains constant along a ray in vacuum. Start from the definition of radiance and use the inverse square law.

**Hint:** Consider two differential areas dA₁ and dA₂ along the ray and show that L₁ = L₂.

<details>
<summary>Solution</summary>

Consider differential areas dA₁ and dA₂ at distances r₁ and r₂ along a ray. The solid angles subtended are:

dω₁ = dA₂ cos θ₂ / r₁₂²
dω₂ = dA₁ cos θ₁ / r₁₂²

The flux leaving dA₁ toward dA₂ is:
dΦ = L₁ dA₁ cos θ₁ dω₁ = L₁ dA₁ cos θ₁ dA₂ cos θ₂ / r₁₂²

By energy conservation, this equals the flux arriving at dA₂:
dΦ = L₂ dA₂ cos θ₂ dω₂ = L₂ dA₂ cos θ₂ dA₁ cos θ₁ / r₁₂²

Therefore L₁ = L₂, proving radiance invariance.
</details>

### Exercise 1.2: BRDF Energy Conservation
Prove that a Lambertian BRDF f_r = ρ/π satisfies energy conservation for any albedo ρ ≤ 1.

**Hint:** Integrate over the hemisphere using spherical coordinates.

<details>
<summary>Solution</summary>

For Lambertian BRDF f_r = ρ/π, the directional-hemispherical reflectance is:

∫_Ω f_r cos θ_o dω_o = (ρ/π) ∫₀^{2π} ∫₀^{π/2} cos θ sin θ dθ dφ

= (ρ/π) · 2π · ∫₀^{π/2} cos θ sin θ dθ

= (ρ/π) · 2π · [sin² θ/2]₀^{π/2}

= (ρ/π) · 2π · (1/2) = ρ

Since ρ ≤ 1, energy conservation is satisfied.
</details>

### Exercise 1.3: Monte Carlo Variance
Derive the optimal importance sampling distribution for estimating ∫₀¹ √x dx and calculate the resulting variance.

**Hint:** The optimal PDF is proportional to |f(x)|.

<details>
<summary>Solution</summary>

For f(x) = √x on [0,1], the optimal PDF is:

p*(x) = √x / ∫₀¹ √x dx = √x / (2/3) = (3/2)√x

To sample: X = F⁻¹(U) where F(x) = ∫₀ˣ (3/2)√t dt = x^{3/2}

Therefore X = U^{2/3}

With this sampling, the estimator becomes:
f(X)/p*(X) = √X / ((3/2)√X) = 2/3

Since the estimator is constant, the variance is 0—perfect importance sampling eliminates variance.
</details>

### Exercise 1.4: Path Integral Convergence (Challenge)
Show that the Neumann series for the rendering equation converges when the maximum albedo in the scene is less than 1.

**Hint:** Use the operator norm ||𝒯|| ≤ ρ_max < 1.

<details>
<summary>Solution</summary>

The transport operator 𝒯 satisfies:

||(𝒯L)||_∞ ≤ ρ_max ||L||_∞

where ρ_max = max_{x,ω} ∫_Ω f_r(x,ω',ω) cos θ' dω' < 1

By induction: ||𝒯^k L_e||_∞ ≤ ρ_max^k ||L_e||_∞

The Neumann series converges:
||∑_{k=0}^∞ 𝒯^k L_e||_∞ ≤ ∑_{k=0}^∞ ρ_max^k ||L_e||_∞ = ||L_e||_∞/(1-ρ_max)

This proves convergence for ρ_max < 1.
</details>

### Exercise 1.5: Coordinate Transform Jacobian
Derive the Jacobian for converting the rendering equation from integration over solid angle to integration over surface area.

**Hint:** Start with dω = dA cos θ'/r².

<details>
<summary>Solution</summary>

Given points x and x' with connecting vector r = x' - x:

dω = dA' cos θ' / ||r||²

The rendering equation becomes:

L_o(x,ω_o) = L_e(x,ω_o) + ∫_M f_r(x,ω(x'),ω_o) L_i(x',−ω(x')) V(x↔x') (cos θ cos θ')/||r||² dA'

where:
- ω(x') = (x' - x)/||x' - x|| 
- cos θ = n(x) · ω(x')
- cos θ' = -n(x') · ω(x')
- V(x↔x') is the visibility function

The Jacobian is |∂ω/∂x'| = cos θ'/||r||².
</details>

### Exercise 1.6: BSDF Reciprocity with Refraction (Challenge)
Derive the generalized reciprocity relation for BSDF with refraction between media with refractive indices n₁ and n₂.

**Hint:** Use radiance in phase space L̃ = L/n².

<details>
<summary>Solution</summary>

Define phase space radiance L̃ = L/n². In a medium with index n:

L̃ is invariant along rays (generalized radiance theorem)

At an interface, power conservation requires:

L̃₁(x,ω_i) dω_i = L̃₂(x,ω_t) dω_t

Using Snell's law: n₁ sin θ_i = n₂ sin θ_t

The solid angle ratio is: dω_t/dω_i = (n₁/n₂)² cos θ_t/cos θ_i

For the BTDF: f_t(ω_i→ω_t) = dL₂/dE₁ = (n₂²/n₁²) dL̃₂/dẼ₁

By time-reversal symmetry of Maxwell's equations:

n₁² f_t(x,ω_i→ω_t) = n₂² f_t(x,ω_t→ω_i)
</details>

### Exercise 1.7: Russian Roulette Bias
Prove that Russian roulette with continuation probability q creates an unbiased estimator.

**Hint:** Use the law of total expectation.

<details>
<summary>Solution</summary>

Let L be the true value and L' be the Russian roulette estimator:

L' = {
  L/q  with probability q
  0    with probability 1-q
}

The expected value is:

E[L'] = q · (L/q) + (1-q) · 0 = L

Therefore E[L'] = L, proving the estimator is unbiased.

The variance is:
Var[L'] = E[L'²] - (E[L'])² = q(L/q)² - L² = L²/q - L² = L²(1-q)/q

Note that variance increases as q decreases, showing the bias-variance tradeoff.
</details>

### Exercise 1.8: Volume Rendering Convergence (Open Problem)
Consider the volume rendering equation with spatially varying extinction. Under what conditions does the path integral formulation converge? Discuss the relationship between optical thickness and convergence rate.

**Hint:** Consider the maximum optical thickness τ_max along any ray.

<details>
<summary>Solution</summary>

For volumes with bounded extinction σ_t ≤ σ_max and albedo α = σ_s/σ_t ≤ α_max < 1:

The k-th scattering term is bounded by:
||L^{(k)}||_∞ ≤ ||L_e||_∞ α_max^k

This ensures geometric convergence.

For optical thickness τ = ∫ σ_t ds along a ray:
- Low τ: Single scattering dominates, fast convergence
- High τ: Multiple scattering important, slower convergence
- τ → ∞: Diffusion regime, may need specialized methods

Open questions:
1. Optimal importance sampling for heterogeneous media
2. Convergence rate vs. spatial frequency of σ_t
3. Connection to transport mean free path in multiple scattering
</details>

## Common Pitfalls and Errors (Gotchas)

1. **Radiance vs. Irradiance Confusion**
   - Radiance has units W/(m²·sr), irradiance has W/m²
   - Always check dimensional consistency in equations
   - Remember: cameras measure irradiance (integrated radiance)

2. **Incorrect Normal Transformation**
   - Normals transform by inverse transpose, not the direct matrix
   - Failing to renormalize after transformation
   - Sign errors when transforming between coordinate systems

3. **BRDF Energy Conservation Violations**
   - Forgetting the cosine term in the integral
   - Not accounting for Fresnel effects at grazing angles
   - Incorrect normalization of analytical BRDF models

4. **Monte Carlo Bias from PDF Errors**
   - PDF must be > 0 wherever integrand is non-zero
   - Normalization errors in importance sampling
   - Forgetting Jacobian when changing variables

5. **Numerical Precision Issues**
   - Ray-sphere intersection can fail due to catastrophic cancellation
   - Small denominators in geometry factor G(x↔x')
   - Accumulated error in long ray paths

## Best Practices Checklist

### Design Review

- [ ] **Verify Physical Units**: Check all equations for dimensional consistency
- [ ] **Energy Conservation**: Ensure BRDFs satisfy ∫ f_r cos θ dω ≤ 1
- [ ] **Reciprocity**: Verify f_r(ω_i, ω_o) = f_r(ω_o, ω_i)
- [ ] **Coordinate System Consistency**: Document and verify all coordinate conventions
- [ ] **Sampling Strategy**: Match importance sampling to dominant contributions

### Implementation Review

- [ ] **Robust Ray-Primitive Intersection**: Use numerically stable algorithms
- [ ] **PDF Validation**: Assert p(x) > 0 and ∫ p(x) dx = 1
- [ ] **Variance Reduction**: Implement MIS for multiple sampling strategies
- [ ] **Floating Point Hygiene**: Check for NaN/Inf propagation
- [ ] **Termination Criteria**: Use Russian roulette for infinite recursion

### Debugging Checklist

- [ ] **White Furnace Test**: Uniform emission should produce uniform radiance
- [ ] **Reciprocity Test**: Swap light and camera, verify same result
- [ ] **Energy Audit**: Track total flux in = total flux out
- [ ] **Convergence Analysis**: Plot variance vs. sample count
- [ ] **Reference Comparison**: Validate against analytical solutions when available
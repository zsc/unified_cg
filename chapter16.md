# Chapter 16: Coherence Theory

This chapter introduces the fundamental concepts of optical coherence, bridging the gap between deterministic wave optics and statistical descriptions of light fields. We develop the mathematical framework for describing partial coherence, establishing key theorems that govern how coherence properties transform through optical systems and propagate through space. These concepts are essential for understanding advanced rendering effects involving interference and diffraction with realistic light sources.

## 16.1 Temporal Coherence and Spectra

Temporal coherence describes the correlation of a light wave with itself at different time delays. For a scalar optical field E(t), we define the temporal coherence function:

**Γ(τ) = ⟨E*(t)E(t+τ)⟩**

where ⟨·⟩ denotes time averaging and * indicates complex conjugation.

For a stationary random process, this can be written more explicitly as:

**Γ(τ) = lim_{T→∞} (1/T) ∫_{-T/2}^{T/2} E*(t)E(t+τ) dt**

The physical interpretation is straightforward: Γ(τ) measures how similar the field is to a time-delayed version of itself. For τ = 0, we have Γ(0) = ⟨|E(t)|²⟩, which is the time-averaged intensity.

### Statistical Foundation

The field E(t) is treated as a complex-valued random process. For ergodic processes, the time average equals the ensemble average:

**Γ(τ) = E[E*(t)E(t+τ)]**

where E[·] denotes expectation value. The autocorrelation function satisfies:

**Γ(-τ) = Γ*(τ)** (Hermitian symmetry)

This follows from stationarity: ⟨E*(t)E(t-τ)⟩ = ⟨E*(t+τ)E(t)⟩ = ⟨E(t)E*(t+τ)⟩*.

### Connection to Interferometry

In a Michelson interferometer with path difference Δ = cτ, the intensity at the output is:

**I_out = I₁ + I₂ + 2Re[√(I₁I₂)γ(τ)exp(iω₀τ)]**

where I₁, I₂ are intensities from the two arms and γ(τ) = Γ(τ)/Γ(0) is the normalized coherence function. The fringe visibility directly measures |γ(τ)|:

**V(τ) = 2√(I₁I₂)|γ(τ)|/(I₁ + I₂)**

For balanced arms (I₁ = I₂), we have V(τ) = |γ(τ)|, providing a direct experimental method to measure temporal coherence.

### 16.1.1 Coherence Time and Length

The normalized degree of temporal coherence is:

**γ(τ) = Γ(τ)/Γ(0)**

This complex-valued function satisfies several important properties:
- **γ(0) = 1** (perfect self-correlation at zero delay)
- **|γ(τ)| ≤ 1** (Schwarz inequality)
- **γ(-τ) = γ*(τ)** (Hermitian symmetry)
- **γ(τ) → 0** as **τ → ∞** for finite bandwidth sources

The coherence time τ_c can be defined in several equivalent ways:

1. **1/e definition**: Time at which |γ(τ)| = 1/e
2. **FWHM definition**: Full width at half maximum of |γ(τ)|²
3. **Integral definition**: **τ_c = ∫₀^∞ |γ(τ)|² dτ**

The integral definition is most fundamental as it represents the effective duration over which the field maintains correlation. The coherence length is then:

**l_c = c·τ_c**

where c is the speed of light in the medium (c = c₀/n for refractive index n).

### 16.1.2 Spectral Width Relationship

By the Wiener-Khinchin theorem (detailed in Section 16.3), the temporal coherence function is the Fourier transform of the power spectral density S(ω):

**Γ(τ) = ∫_{-∞}^∞ S(ω)e^{iωτ} dω**

**S(ω) = (1/2π) ∫_{-∞}^∞ Γ(τ)e^{-iωτ} dτ**

This Fourier transform relationship implies a fundamental uncertainty principle:

**Δω·Δτ ≥ K**

where K is a constant of order unity depending on how the widths are defined.

#### Derivation of the Uncertainty Principle

Starting from the Schwarz inequality for Fourier transforms:

**∫|f(t)|² dt · ∫|F(ω)|² dω ≥ (1/4π)|∫f(t) dt|²**

Applying this to f(t) = tΓ(t) and using integration by parts:

**⟨τ²⟩⟨ω²⟩ ≥ 1/4**

where ⟨τ²⟩ = ∫τ²|γ(τ)|² dτ/∫|γ(τ)|² dτ and ⟨ω²⟩ = ∫(ω-ω₀)²S(ω) dω/∫S(ω) dω.

This gives the minimum uncertainty product:

**Δτ_rms · Δω_rms ≥ 1/2**

For other definitions:
- FWHM: Δτ_FWHM · Δω_FWHM ≥ 4ln(2)/π ≈ 0.88
- 1/e width: Δτ_{1/e} · Δω_{1/e} ≥ 1

For specific spectral shapes:

**1. Gaussian Spectrum:**
**S(ω) = S₀ exp[-(ω-ω₀)²/(2σ²)]**

where σ = 2πΔν/(2√(2ln2)) for FWHM Δν, yielding:

**γ(τ) = exp(iω₀τ) exp(-σ²τ²/2)**

**τ_c = ∫₀^∞ exp(-σ²τ²) dτ = √π/(2σ) ≈ 0.44/Δν**

**2. Lorentzian Spectrum:**
**S(ω) = S₀Γ₀/π/[(ω-ω₀)² + Γ₀²]**

This gives:
**γ(τ) = exp(iω₀τ) exp(-Γ₀|τ|)**

**τ_c = 1/(2Γ₀) = 1/(2πΔν)**

**3. Rectangular Spectrum:**
**S(ω) = S₀ for |ω-ω₀| < πΔν, 0 otherwise**

**γ(τ) = exp(iω₀τ) sinc(πΔντ)**

**τ_c ≈ 1/Δν**

The Gaussian gives the minimum time-bandwidth product, while the rectangular spectrum produces the longest coherence time for a given bandwidth.

### 16.1.3 Examples of Coherence Times

Different light sources exhibit vastly different coherence properties:

**1. Stabilized He-Ne Laser (single mode):**
- Δν ~ 1 kHz - 1 MHz
- τ_c ~ 10⁻³ - 1 s
- l_c ~ 300 km - 300,000 km
- Applications: Interferometry, holography
- Spectral profile: Nearly Lorentzian due to cavity dynamics

**2. Semiconductor Laser Diode:**
- Δν ~ 1-10 MHz (single mode)
- τ_c ~ 10⁻⁷ - 10⁻⁸ s
- l_c ~ 30 - 300 m
- Applications: Fiber optic communications, CD/DVD readers
- Temperature dependence: Δν ∝ T³/² (carrier scattering)

**3. Light Emitting Diode (LED):**
- Δν ~ 10 THz (Δλ ~ 30 nm for visible)
- τ_c ~ 10⁻¹⁴ s
- l_c ~ 3 μm
- Applications: Display technology, illumination
- Spectrum shape: Approximately Gaussian from band-edge emission

**4. Filtered Sunlight (1 nm filter):**
- Δν ~ 1 THz
- τ_c ~ 10⁻¹² s  
- l_c ~ 0.3 mm
- Natural illumination reference
- Blackbody spectrum: S(ω) ∝ ω³/(exp(ℏω/k_BT) - 1)

**5. White Light (unfiltered):**
- Δν ~ 300 THz (400-700 nm range)
- τ_c ~ 10⁻¹⁵ s
- l_c ~ 0.3 μm
- Essentially incoherent for most applications
- Coherence function: γ(τ) ≈ rect(τ/τ_c)exp(iω_cτ)

**6. Sodium D-line (low pressure lamp):**
- Δν ~ 500 MHz (Doppler broadened)
- τ_c ~ 10⁻⁹ s
- l_c ~ 0.3 m
- Classic atomic spectral line
- Fine structure: D₁ (589.6 nm) and D₂ (589.0 nm) doublet

**7. Synchrotron Radiation:**
- Δν/ν ~ 10⁻³ (typical undulator)
- τ_c ~ 10⁻¹² s at λ = 1 nm
- l_c ~ 0.3 mm
- Highly directional, partially coherent beam

**8. Free Electron Laser (FEL):**
- Δν/ν ~ 10⁻⁴ - 10⁻³
- τ_c ~ 10⁻¹¹ - 10⁻¹² s (X-ray regime)
- l_c ~ 3-30 mm
- SASE process produces partial temporal coherence

### 16.1.4 Measurement of Coherence Time

Coherence time can be measured through several methods:

**1. Michelson Interferometry:**
Measure fringe visibility V as function of path difference Δ:
**V(Δ) = |γ(Δ/c)|**

The coherence length l_c is where V drops to 1/e or 1/2.

For quasi-monochromatic light:
**I(Δ) = I₀[1 + V(Δ)cos(k₀Δ + φ)]**

The visibility envelope V(Δ) directly traces |γ(τ)|. Key considerations:
- Mechanical stability: Δ stable to λ/20
- Equal arm balance: Compensate dispersion
- Detector response: Must resolve fringes

**2. Spectral Analysis:**
Measure power spectrum S(ω) directly using spectrometer, then calculate:
**τ_c = 2π/Δω_eff**

where Δω_eff is effective spectral width.

Effective width definitions:
- RMS width: Δω_rms = √(⟨ω²⟩ - ⟨ω⟩²)
- FWHM: Full width at half maximum
- Equivalent width: Δω_eq = ∫S(ω)dω/S_max

Resolution requirements:
- Spectrometer resolution δω << Δω
- Free spectral range > full spectrum width

**3. Intensity Correlation:**
For thermal light, measure intensity correlation function:
**g⁽²⁾(τ) = ⟨I(t)I(t+τ)⟩/⟨I(t)⟩²**

Related to field correlation by:
**g⁽²⁾(τ) = 1 + |γ(τ)|²**

This is the Hanbury Brown-Twiss effect. Implementation:
- Fast detectors: Response time << τ_c
- Correlator: Digital or analog
- Photon counting: For weak signals

**4. Fourier Transform Spectroscopy:**
Scan Michelson interferometer over large Δ range:
**I(Δ) = ∫S(ω)[1 + cos(ωΔ/c)]dω**

Fourier transform gives S(ω):
**S(ω) ∝ ℱ{I(Δ) - I_∞}**

Advantages:
- Multiplex advantage (Fellgett)
- High light gathering (Jacquinot)
- Accurate wavelength calibration

## 16.2 Spatial Coherence and Young's Experiment

Spatial coherence describes the correlation between light fields at different spatial points at the same time. Young's double-slit experiment provides the canonical framework for understanding spatial coherence.

The fundamental question of spatial coherence is: given two points P₁ and P₂ in space, how correlated are the optical fields at these points? This correlation determines whether interference fringes can be observed when light from these two points is combined.

### Historical Context and Modern Applications

Young's 1801 experiment not only demonstrated the wave nature of light but established the principle that coherence determines interference visibility. Modern applications include:

- **Stellar interferometry**: Measuring star diameters using baseline correlation
- **Optical coherence tomography**: Depth imaging using coherence gating
- **Lithography**: Partial coherence effects in pattern transfer
- **Quantum optics**: EPR correlations and Bell inequality tests

### 16.2.1 Mutual Intensity and Visibility

Consider two points P₁ and P₂ illuminated by a source. The mutual intensity is:

**J₁₂ = ⟨E*(P₁,t)E(P₂,t)⟩**

This complex quantity contains both amplitude and phase information about the correlation.

#### Mathematical Structure

The mutual intensity forms a Hermitian matrix:

**J = [J₁₁  J₁₂]**
**    [J₂₁  J₂₂]**

with properties:
- J₁₁, J₂₂ ≥ 0 (intensities)
- J₂₁ = J₁₂* (Hermitian)
- det(J) ≥ 0 (positive semi-definite)
- |J₁₂|² ≤ J₁₁J₂₂ (Schwarz inequality)

#### Young's Double Slit Analysis

In Young's experiment with slits at positions r₁ and r₂, the field at observation point P is:

**E(P,t) = K₁E(r₁,t-t₁) + K₂E(r₂,t-t₂)**

where Kᵢ are propagation factors and tᵢ = |P-rᵢ|/c are propagation times.

The intensity becomes:

**I(P) = |K₁|²I₁ + |K₂|²I₂ + 2Re[K₁*K₂J₁₂ exp(ikΔ)]**

where Δ = |P-r₂| - |P-r₁| is the path difference.

For the paraxial case with equal propagation factors:

**I(P) = I₁ + I₂ + 2√(I₁I₂)|J₁₂|cos(kΔ + φ₁₂)**

where φ₁₂ = arg(J₁₂) is the phase of mutual intensity.

#### Fringe Visibility

The intensity varies between:
- **I_max = I₁ + I₂ + 2√(I₁I₂)|γ₁₂|**
- **I_min = I₁ + I₂ - 2√(I₁I₂)|γ₁₂|**

The visibility is defined as:

**V = (I_max - I_min)/(I_max + I_min) = 2√(I₁I₂)|γ₁₂|/(I₁ + I₂)**

Special cases:
1. **Equal intensities** (I₁ = I₂ = I₀): **V = |γ₁₂|**
2. **Fully coherent light**: |γ₁₂| = 1, **V = 2√(I₁I₂)/(I₁ + I₂) ≤ 1**
3. **Incoherent light**: γ₁₂ = 0, **V = 0** (no fringes)

### 16.2.2 Complex Degree of Coherence

The complex degree of coherence is defined as:

**γ₁₂ = J₁₂/√(J₁₁J₂₂) = J₁₂/√(I₁I₂)**

#### Mathematical Properties

1. **Schwarz Inequality**: |γ₁₂| ≤ 1
   Proof: By Cauchy-Schwarz, |⟨E₁*E₂⟩|² ≤ ⟨|E₁|²⟩⟨|E₂|²⟩

2. **Hermitian Symmetry**: γ₁₂ = γ₂₁*
   Since J₁₂ = ⟨E₁*E₂⟩ = ⟨E₂*E₁⟩* = J₂₁*

3. **Self-coherence**: γ₁₁ = γ₂₂ = 1
   Every point is perfectly coherent with itself

4. **Phase Information**: γ₁₂ = |γ₁₂|exp(iφ₁₂)
   The phase φ₁₂ affects fringe position, not visibility

#### Physical Interpretation

- **|γ₁₂| = 1**: Fully coherent - perfect correlation
- **0 < |γ₁₂| < 1**: Partially coherent - reduced fringe visibility  
- **|γ₁₂| = 0**: Incoherent - no correlation, no fringes

#### Connection to Classical Correlation

The degree of coherence is analogous to the correlation coefficient in statistics:

**γ₁₂ = Cov(E₁,E₂)/(σ₁σ₂)**

where Cov is covariance and σᵢ are standard deviations.

### 16.2.3 Extended Sources

For extended sources, each point on the source contributes independently (for incoherent sources), and the total mutual intensity is the sum of contributions.

#### General Formulation

For an extended incoherent source with intensity distribution I_s(r_s), the mutual intensity at observation points r₁, r₂ is:

**J₁₂ = ∫∫∫ I_s(r_s) K*(r₁,r_s)K(r₂,r_s) d³r_s**

where K(r,r_s) is the propagation kernel from source point r_s to observation point r.

In the Fresnel approximation:
**K(r,r_s) = (exp(ikR)/R) × exp[ik(r-r_s)²/2R]**

where R = |r_z - r_{s,z}| is the axial distance.

#### Circular Source Analysis

Consider a uniformly illuminated circular incoherent source of diameter D at distance z from the observation plane. For two observation points P₁ and P₂ separated by distance d:

**J₁₂ = ∬_source I_s(ξ,η) exp[ik(r₂-r₁)·r_s/z] dξ dη**

In the paraxial approximation with r₂-r₁ = (d,0,0):

**γ₁₂ = (2/πa²) ∬_{|ρ|<a} exp(ikdξ/z) dξ dη / 1**

Converting to polar coordinates and integrating:

**|γ₁₂| = |2J₁(πDd/λz)/(πDd/λz)|**

where J₁ is the first-order Bessel function of the first kind.

#### Detailed Derivation

Starting with polar coordinates (ρ,θ) in the source plane:

**γ₁₂ = (1/πa²) ∫₀^{2π} ∫₀^a exp(ikdρcos(θ)/z) ρ dρ dθ**

The angular integral gives:
**∫₀^{2π} exp(ikdρcos(θ)/z) dθ = 2πJ₀(kdρ/z)**

where J₀ is the zeroth-order Bessel function. The radial integral becomes:

**γ₁₂ = (2/a²) ∫₀^a J₀(kdρ/z) ρ dρ**

Using the identity ∫₀^x tJ₀(t)dt = xJ₁(x):

**γ₁₂ = 2J₁(kda/z)/(kda/z) = 2J₁(πDd/λz)/(πDd/λz)**

This is the celebrated Airy pattern for coherence.

#### Coherence Radius

The transverse coherence radius ρ_c is defined where |γ₁₂| first reaches zero:

**πDd_coh/λz = 3.832** (first zero of J₁)

**d_coh = ρ_c ≈ 1.22λz/D = 0.61λz/a**

This is identical to the Airy disk radius in diffraction!

#### Other Source Geometries

**1. Rectangular Source (L_x × L_y):**
**γ₁₂ = sinc(πL_xΔx/λz) × sinc(πL_yΔy/λz)**

Coherence lengths: ρ_x = λz/L_x, ρ_y = λz/L_y

**2. Double Star (two point sources, separation θ):**
**γ₁₂ = [exp(iπθd/λ) + exp(-iπθd/λ)]/2 = cos(πθd/λ)**

First zero at d = λ/(2θ) - basis of Michelson stellar interferometry

**3. Gaussian Source (1/e radius w_s):**
**γ₁₂ = exp(-π²w_s²d²/λ²z²)**

Coherence radius (1/e): ρ_c = λz/(πw_s)

#### Angular Diameter and Coherence

For a source subtending angle θ_s = D/z:

**ρ_c ≈ λ/θ_s**

This fundamental relationship shows:
- Smaller angular sources → larger coherence areas
- Stars (microarcseconds) produce meter-scale coherence
- The Sun (θ_s ≈ 0.5°) gives ρ_c ≈ 0.07 mm at 500 nm

## 16.3 Mutual Coherence Function and Wiener-Khinchin Theorem

The mutual coherence function provides a complete second-order statistical description of partially coherent fields. This section develops the mathematical framework connecting coherence in time domain to spectral properties in frequency domain.

### 16.3.1 General Definition

For a scalar field U(r,t), the mutual coherence function is:

**Γ(r₁,r₂,τ) = ⟨U*(r₁,t)U(r₂,t+τ)⟩**

This four-point function (two spatial, one temporal, one ensemble) captures all second-order coherence properties.

#### Properties of Mutual Coherence Function

1. **Hermitian Symmetry**: Γ(r₁,r₂,τ) = Γ*(r₂,r₁,-τ)

2. **Positive Semi-definiteness**: For any functions f₁(r), f₂(r):
   ∬∬ f₁*(r₁)Γ(r₁,r₂,0)f₂(r₂) dr₁dr₂ ≥ 0

3. **Boundary Values**:
   - Γ(r,r,0) = I(r) (intensity)
   - |Γ(r₁,r₂,τ)| ≤ √[Γ(r₁,r₁,0)Γ(r₂,r₂,0)]

4. **Stationarity**: For stationary fields:
   Γ(r₁,r₂,τ) depends only on τ, not absolute time

The equal-time mutual coherence function (mutual intensity) is:

**J(r₁,r₂) = Γ(r₁,r₂,0) = ⟨U*(r₁,t)U(r₂,t)⟩**

The normalized mutual coherence function:

**γ(r₁,r₂,τ) = Γ(r₁,r₂,τ)/√[Γ(r₁,r₁,0)Γ(r₂,r₂,0)]**

### 16.3.2 Cross-Spectral Density

The cross-spectral density is the Fourier transform of mutual coherence with respect to time delay:

**W(r₁,r₂,ω) = ∫_{-∞}^∞ Γ(r₁,r₂,τ)e^{-iωτ} dτ**

**Γ(r₁,r₂,τ) = (1/2π) ∫_{-∞}^∞ W(r₁,r₂,ω)e^{iωτ} dω**

#### Physical Meaning

W(r₁,r₂,ω) represents the correlation between spectral components at frequency ω at points r₁ and r₂. It can be written as:

**W(r₁,r₂,ω) = ⟨Û*(r₁,ω)Û(r₂,ω)⟩**

where Û(r,ω) is the Fourier transform of U(r,t).

#### Properties of Cross-Spectral Density

1. **Hermitian**: W(r₁,r₂,ω) = W*(r₂,r₁,ω)

2. **Non-negative on diagonal**: W(r,r,ω) ≥ 0 (power spectrum)

3. **Spectral degree of coherence**:
   **μ(r₁,r₂,ω) = W(r₁,r₂,ω)/√[W(r₁,r₁,ω)W(r₂,r₂,ω)]**
   
   with |μ(r₁,r₂,ω)| ≤ 1

4. **Total intensity**: 
   **I(r) = (1/2π) ∫_{-∞}^∞ W(r,r,ω) dω**

### 16.3.3 Wiener-Khinchin Theorem

The Wiener-Khinchin theorem establishes the fundamental connection between autocorrelation and power spectrum.

#### Statement of the Theorem

For a stationary field at position r:

**Γ(r,r,τ) = ⟨U*(r,t)U(r,t+τ)⟩ = ∫_{-∞}^∞ S(r,ω)e^{iωτ} dω**

**S(r,ω) = (1/2π) ∫_{-∞}^∞ Γ(r,r,τ)e^{-iωτ} dτ**

where S(r,ω) = W(r,r,ω)/(2π) is the power spectral density.

#### Implications

1. **Energy Conservation**:
   **I(r) = ⟨|U(r,t)|²⟩ = Γ(r,r,0) = ∫_{-∞}^∞ S(r,ω) dω**

2. **Coherence-Bandwidth Product**:
   The theorem implies Δω·τ_c ~ 2π, a form of uncertainty principle

3. **White Light Limit**:
   For S(ω) = S₀ (constant), Γ(τ) = 2πS₀δ(τ) - perfect incoherence

#### Generalized Wiener-Khinchin

For non-stationary fields, we use the Wigner distribution:

**W_U(t,ω) = ∫ U*(t-τ/2)U(t+τ/2)e^{-iωτ} dτ**

The ensemble average gives:
**⟨W_U(t,ω)⟩ = ∫ Γ(t-τ/2,t+τ/2,τ)e^{-iωτ} dτ**

### 16.3.4 Quasi-Monochromatic Approximation

Many practical sources have narrow spectral width compared to center frequency: Δω/ω₀ << 1.

#### Mathematical Formulation

For narrow-band light centered at ω₀:

**U(r,t) = A(r,t)e^{-iω₀t}**

where A(r,t) is a slowly varying complex amplitude.

The mutual coherence function becomes:

**Γ(r₁,r₂,τ) = ⟨A*(r₁,t)A(r₂,t+τ)⟩e^{iω₀τ}**

**= J(r₁,r₂)e^{iω₀τ}g(τ)**

where:
- J(r₁,r₂) = ⟨A*(r₁,t)A(r₂,t)⟩ is the mutual intensity
- g(τ) is a slowly varying envelope with g(0) = 1

#### Validity Conditions

1. **Spectral condition**: Δω/ω₀ << 1
2. **Temporal condition**: Changes in A(t) occur on timescale >> 1/ω₀
3. **Propagation condition**: Δk·L << 2π (L is propagation distance)

#### Cross-Spectral Density

Under quasi-monochromatic approximation:

**W(r₁,r₂,ω) ≈ J(r₁,r₂)G(ω-ω₀)**

where G(ω) is the Fourier transform of g(τ), centered at ω = 0.

#### Applications

1. **Interference**: Phase differences kΔ ≈ k₀Δ use center wavenumber
2. **Diffraction**: Single frequency calculation suffices
3. **Coherence propagation**: Simplifies to monochromatic case
4. **Spectroscopy**: Natural linewidths of atomic transitions

## 16.4 van Cittert-Zernike Theorem

The van Cittert-Zernike theorem is one of the most important results in coherence theory, establishing how spatial coherence emerges from incoherent sources through propagation.

### 16.4.1 Statement of the Theorem

For a planar incoherent source with intensity distribution I_s(ξ,η) in the source plane, the complex degree of coherence in an observation plane at distance z is:

**γ₁₂ = ∫∫ I_s(ξ,η) exp[ik(r₁-r₂)·(ξ,η)/z] dξdη / ∫∫ I_s(ξ,η) dξdη**

In the paraxial approximation, this becomes:

**γ₁₂ = ℱ{I_s(ξ,η)}|_{(u,v)=(x₁-x₂)/λz, (y₁-y₂)/λz} / ∫∫ I_s(ξ,η) dξdη**

where ℱ denotes the 2D Fourier transform.

#### Rigorous Derivation

Starting from the mutual intensity for an incoherent source:

**J(r₁,r₂) = ∫∫ I_s(r_s) G*(r₁,r_s)G(r₂,r_s) d²r_s**

where G is the Green's function. In free space:

**G(r,r_s) = exp(ik|r-r_s|)/4π|r-r_s|**

For paraxial propagation with z >> |r_⊥|, |r_s|:

**|r-r_s| ≈ z + (r_⊥ - r_s)²/2z**

Substituting and simplifying:

**J(r₁,r₂) = (k/2πz)² exp[ik(|r₁|² - |r₂|²)/2z] × ∫∫ I_s(r_s) exp[ik(r₂-r₁)·r_s/z] d²r_s**

The normalized degree of coherence becomes:

**γ₁₂ = exp[ik(|r₁|² - |r₂|²)/2z] × ℱ{I_s}(k(r₂-r₁)/2πz) / I_total**

For points at equal distance from the axis, the quadratic phase cancels, giving the classic result.

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
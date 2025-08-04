# Chapter 7: Dynamic Neural Radiance Fields

The real world is inherently dynamic. While static Neural Radiance Fields (NeRF) have revolutionized 3D scene representation, extending them to handle temporal changes opens up possibilities for modeling everything from human performances to fluid dynamics. This chapter explores how to incorporate time as an additional dimension in neural radiance fields, transforming the volume rendering equation to handle dynamic scenes while maintaining the mathematical elegance of our unified framework.

Dynamic scenes present unique challenges: capturing fast motions without aliasing, maintaining temporal coherence, handling topological changes, and managing the exponential growth in data complexity when adding time as a fourth dimension. We'll see how careful mathematical formulation, combined with insights from physics and signal processing, leads to practical and efficient solutions.

## Learning Objectives

After completing this chapter, you will be able to:

1. **Formulate** the time-extended volume rendering equation for dynamic scenes
2. **Design** deformation fields that map between canonical and observation spaces
3. **Derive** optical and scene flow constraints from the volume rendering framework
4. **Analyze** the computational and memory tradeoffs of different temporal representations
5. **Implement** regularization strategies for temporally coherent reconstructions
6. **Evaluate** the convergence properties of dynamic neural radiance fields

## Prerequisites

This chapter builds upon:
- Chapter 6: Neural Radiance Fields (static formulation)
- Understanding of vector calculus (Jacobians, divergence theorem)
- Familiarity with optimization and regularization techniques
- Basic knowledge of continuum mechanics (helpful but not required)

## 7.1 Time-varying Radiance Field Representations

Dynamic radiance fields require us to model how both geometry (density $\sigma$) and appearance (color $\mathbf{c}$) evolve over time. The key insight is that while the mathematical framework remains elegant, the choice of temporal representation profoundly impacts both expressiveness and computational efficiency.

### 7.1.1 Extending the Volume Rendering Equation

The static volume rendering equation from Chapter 6:

$$C(\mathbf{r}) = \int_{t_n}^{t_f} T(t)\sigma(\mathbf{r}(t))\mathbf{c}(\mathbf{r}(t), \mathbf{d}) dt$$

where $T(t) = \exp\left(-\int_{t_n}^{t} \sigma(\mathbf{r}(s))ds\right)$

For dynamic scenes, we introduce time $\tau$ as an additional parameter (note: we use $\tau$ for scene time to distinguish from ray parameter $t$):

$$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T(t, \tau)\sigma(\mathbf{r}(t), \tau)\mathbf{c}(\mathbf{r}(t), \mathbf{d}, \tau) dt$$

The transmittance also becomes time-dependent:

$$T(t, \tau) = \exp\left(-\int_{t_n}^{t} \sigma(\mathbf{r}(s), \tau)ds\right)$$

This seemingly simple extension has profound implications:

1. **Increased Dimensionality**: The radiance field is now $f: \mathbb{R}^3 \times \mathbb{S}^2 \times \mathbb{R} \rightarrow \mathbb{R}^4$
2. **Temporal Coherence**: Neighboring time frames should produce similar radiance values
3. **Motion Blur**: Fast motions during exposure time require integration over $\tau$

For motion blur modeling, the observed color becomes:

$$C_{\text{blur}}(\mathbf{r}) = \frac{1}{\Delta\tau} \int_{\tau_0}^{\tau_0+\Delta\tau} C(\mathbf{r}, \tau) d\tau$$

where $\Delta\tau$ is the exposure time.

### 7.1.2 Discrete vs Continuous Time Modeling

The choice between discrete and continuous time representations involves fundamental tradeoffs between expressiveness, memory efficiency, and optimization complexity.

**Discrete Time Representation:**
For a sequence of $N$ time steps $\{\tau_i\}_{i=1}^N$, we can represent the radiance field as:

$$\mathcal{F}_{\text{discrete}} = \{f_{\theta_i}: \mathbb{R}^3 \times \mathbb{S}^2 \rightarrow \mathbb{R}^4\}_{i=1}^N$$

Advantages:
- Perfect reconstruction at sampled times
- No temporal aliasing at discrete points
- Independent optimization per frame

Disadvantages:
- Memory complexity: $O(N \cdot |\theta|)$
- No intermediate frame generation
- Temporal discontinuities possible

**Continuous Time Representation:**
A single network that takes time as input:

$$f_\theta: \mathbb{R}^3 \times \mathbb{S}^2 \times \mathbb{R} \rightarrow \mathbb{R}^4$$

The time dimension typically uses positional encoding:
$$\gamma(\tau) = [\sin(2^0\pi\tau), \cos(2^0\pi\tau), ..., \sin(2^{L-1}\pi\tau), \cos(2^{L-1}\pi\tau)]$$

Advantages:
- Memory complexity: $O(|\theta|)$ (independent of sequence length)
- Natural temporal interpolation
- Smooth motion trajectories

Disadvantages:
- Limited temporal resolution (bounded by network capacity)
- Potential temporal aliasing
- Harder to fit rapid changes

**Hybrid Approaches:**
Combine benefits by using time-conditioned hypernetworks:

$$f_\theta(\mathbf{x}, \mathbf{d}, \tau) = g_{\phi(\tau)}(\mathbf{x}, \mathbf{d})$$

where $\phi(\tau)$ generates time-specific network parameters.

### 7.1.3 Temporal Basis Functions

To balance expressiveness and efficiency, we can decompose the time-varying field using basis functions:

$$\sigma(\mathbf{x}, \tau) = \sum_{k=1}^K \sigma_k(\mathbf{x}) \cdot \phi_k(\tau)$$

$$\mathbf{c}(\mathbf{x}, \mathbf{d}, \tau) = \sum_{k=1}^K \mathbf{c}_k(\mathbf{x}, \mathbf{d}) \cdot \phi_k(\tau)$$

This separable representation reduces the problem to learning $K$ spatial fields and temporal basis functions.

**Common Basis Function Choices:**

1. **Fourier Basis** (for periodic motions):
   $$\phi_k(\tau) = \begin{cases}
   \sin(2\pi k\tau/T) & k \text{ odd} \\
   \cos(2\pi k\tau/T) & k \text{ even}
   \end{cases}$$
   
   Properties:
   - Global support (affects entire timeline)
   - Orthogonal basis
   - Natural for cyclic motions
   - Spectral interpretation via FFT

2. **B-spline Basis** (for local control):
   $$\phi_k(\tau) = B_n\left(\frac{\tau - \tau_k}{\Delta\tau}\right)$$
   
   where $B_n$ is the $n$-th order B-spline kernel.
   
   Properties:
   - Compact support: $\phi_k(\tau) = 0$ for $|\tau - \tau_k| > n\Delta\tau/2$
   - $C^{n-1}$ continuity
   - Local editing capability
   - Efficient evaluation

3. **Learned Basis** (data-driven):
   $$\phi_k(\tau) = \text{softmax}(\text{MLP}_k(\gamma(\tau)))_k$$
   
   Properties:
   - Adapts to data distribution
   - Can capture complex temporal patterns
   - Requires careful regularization

**Optimal Basis Selection:**
The choice depends on the expected motion characteristics:
- Periodic scenes → Fourier basis
- Localized changes → B-splines
- Complex, non-periodic → Learned basis

The truncation parameter $K$ controls the temporal bandwidth:
$$K \approx 2 \cdot f_{\text{max}} \cdot T$$

where $f_{\text{max}}$ is the maximum motion frequency and $T$ is the sequence duration.

### 7.1.4 Frequency Analysis of Temporal Changes

Understanding the frequency content of dynamic scenes is crucial for avoiding aliasing and ensuring adequate temporal resolution.

**Nyquist-Shannon Sampling Theorem:**
The fundamental requirement for alias-free reconstruction:

$$f_{\text{sample}} \geq 2 f_{\text{max}}$$

where $f_{\text{max}}$ is the maximum frequency of scene changes.

**Motion-Specific Analysis:**

1. **Periodic Motion** (e.g., rotating objects):
   For motion with angular frequency $\omega$:
   $$\mathbf{x}(\tau) = \mathbf{R}(\omega\tau)\mathbf{x}_0$$
   
   Minimum samples per period:
   $$N_{\text{min}} = \frac{4\pi}{\omega \Delta\tau} \geq 4$$

2. **Linear Motion** (constant velocity $\mathbf{v}$):
   $$\mathbf{x}(\tau) = \mathbf{x}_0 + \mathbf{v}\tau$$
   
   Frequency content depends on observation window:
   $$f_{\text{max}} \approx \frac{\|\mathbf{v}\|}{2\lambda_{\text{min}}}$$
   
   where $\lambda_{\text{min}}$ is the smallest spatial feature size.

3. **Accelerated Motion**:
   $$\mathbf{x}(\tau) = \mathbf{x}_0 + \mathbf{v}_0\tau + \frac{1}{2}\mathbf{a}\tau^2$$
   
   Instantaneous frequency:
   $$f(\tau) = \frac{\|\mathbf{v}_0 + \mathbf{a}\tau\|}{2\pi\lambda_{\text{min}}}$$

**Practical Sampling Strategy:**

Given a capture rate $f_{\text{capture}}$, the representable motion bandwidth is:

$$f_{\text{motion}} < \frac{f_{\text{capture}}}{2}$$

For neural representations with positional encoding level $L$:
$$f_{\text{neural}} \leq 2^{L-1}$$

The effective temporal resolution is:
$$f_{\text{effective}} = \min(f_{\text{motion}}, f_{\text{neural}})$$

**Anti-aliasing Strategies:**

1. **Temporal Supersampling**: Sample at $kf_{\text{sample}}$ and average
2. **Motion Blur Integration**: Explicitly model exposure time
3. **Adaptive Sampling**: Increase samples in high-motion regions

## 7.2 Deformation Field Modeling

Rather than learning a time-varying radiance field directly, we can decompose the problem into learning a static canonical representation and a time-varying deformation field. This approach leverages the insight that many dynamic scenes exhibit strong spatial-temporal correlations that can be captured through warping.

### 7.2.1 Forward Deformation Fields

The core idea is to model a deformation field that maps points from a canonical (reference) space to the observation space at each time:

$$\mathbf{x}_{\text{obs}} = \mathbf{W}(\mathbf{x}_{\text{can}}, \tau)$$

where:
- $\mathbf{x}_{\text{can}} \in \mathbb{R}^3$: Point in canonical space
- $\mathbf{x}_{\text{obs}} \in \mathbb{R}^3$: Corresponding point in observation space at time $\tau$
- $\mathbf{W}: \mathbb{R}^3 \times \mathbb{R} \rightarrow \mathbb{R}^3$: Deformation field

**Volume Rendering with Deformation:**

The key challenge is correctly transforming the volume rendering integral. Starting from:

$$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T(t, \tau)\sigma(\mathbf{r}(t), \tau)\mathbf{c}(\mathbf{r}(t), \mathbf{d}, \tau) dt$$

We substitute the warped coordinates and must account for the change of variables:

$$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T(t)\sigma_{\text{can}}(\mathbf{W}^{-1}(\mathbf{r}(t), \tau))\mathbf{c}_{\text{can}}(\mathbf{W}^{-1}(\mathbf{r}(t), \tau), \mathbf{d}_{\text{can}}) \left|\det J_{\mathbf{W}^{-1}}\right| dt$$

where:
- $J_{\mathbf{W}^{-1}} = \frac{\partial \mathbf{W}^{-1}}{\partial \mathbf{x}}$ is the Jacobian of the inverse deformation
- $\mathbf{d}_{\text{can}}$ is the transformed view direction

**Density Transformation:**
The density transforms according to:
$$\sigma_{\text{obs}}(\mathbf{x}, \tau) = \sigma_{\text{can}}(\mathbf{W}^{-1}(\mathbf{x}, \tau)) \left|\det J_{\mathbf{W}^{-1}}(\mathbf{x}, \tau)\right|$$

This ensures mass conservation: integrating density over any region gives the same total mass in both canonical and observation spaces.

### 7.2.2 Volume Preservation Constraints

Many physical deformations preserve volume locally (incompressibility), which provides powerful constraints for regularization.

**Incompressibility Condition:**
For volume-preserving deformations:

$$\det(J_{\mathbf{W}}) = 1$$

where $J_{\mathbf{W}} = \frac{\partial \mathbf{W}}{\partial \mathbf{x}_{\text{can}}}$ is the deformation gradient.

**Velocity Field Formulation:**
Define the velocity field as the time derivative of deformation:

$$\mathbf{v}(\mathbf{x}, \tau) = \frac{\partial \mathbf{W}(\mathbf{x}, \tau)}{\partial \tau}$$

The incompressibility constraint in Eulerian form:

$$\nabla \cdot \mathbf{v} = 0$$

**Practical Enforcement:**

1. **Soft Constraint** (penalty method):
   $$\mathcal{L}_{\text{volume}} = \lambda \int_{\Omega} (\det(J_{\mathbf{W}}) - 1)^2 d\mathbf{x}$$

2. **Projection Method**:
   Decompose velocity into incompressible and potential components:
   $$\mathbf{v} = \mathbf{v}_{\text{incomp}} + \nabla \phi$$
   
   Then solve the Poisson equation:
   $$\nabla^2 \phi = \nabla \cdot \mathbf{v}_{\text{initial}}$$
   
   And project: $\mathbf{v}_{\text{incomp}} = \mathbf{v}_{\text{initial}} - \nabla \phi$

3. **Stream Function** (2D) or **Vector Potential** (3D):
   For 2D: $\mathbf{v} = \nabla^\perp \psi = (\frac{\partial \psi}{\partial y}, -\frac{\partial \psi}{\partial x})$
   
   Automatically satisfies $\nabla \cdot \mathbf{v} = 0$

**Approximate Volume Preservation:**
For small deformations, linearize the constraint:
$$\det(J_{\mathbf{W}}) \approx 1 + \text{tr}(J_{\mathbf{W}} - I) = 1 + \nabla \cdot (\mathbf{W} - \mathbf{x})$$

Leading to the simpler constraint:
$$\nabla \cdot \mathbf{u} = 0$$

where $\mathbf{u} = \mathbf{W} - \mathbf{x}$ is the displacement field.

### 7.2.3 Regularization Strategies

Deformation fields have many degrees of freedom and require careful regularization to ensure physically plausible and stable solutions.

**1. Spatial Smoothness Regularization:**
Penalizes high-frequency spatial variations:

$$\mathcal{L}_{\text{smooth}} = \int_{\Omega} \int_{\mathcal{T}} \|\nabla^2 \mathbf{W}(\mathbf{x}, \tau)\|_F^2 d\mathbf{x} d\tau$$

Alternative first-order smoothness:
$$\mathcal{L}_{\text{grad}} = \int_{\Omega} \int_{\mathcal{T}} \|\nabla \mathbf{W} - I\|_F^2 d\mathbf{x} d\tau$$

**2. Rigidity Regularization:**
Encourages locally rigid transformations (rotation + translation):

$$\mathcal{L}_{\text{rigid}} = \int_{\Omega} \|J_{\mathbf{W}}^T J_{\mathbf{W}} - I\|_F^2 d\mathbf{x}$$

This is zero when $J_{\mathbf{W}}$ is a rotation matrix. For computational efficiency, use the polar decomposition:
$$J_{\mathbf{W}} = \mathbf{R}\mathbf{S}$$

where $\mathbf{R}$ is rotation and $\mathbf{S}$ is symmetric. Then penalize:
$$\mathcal{L}_{\text{rigid}} = \|\mathbf{S} - I\|_F^2$$

**3. Temporal Coherence:**
Ensures smooth motion over time:

$$\mathcal{L}_{\text{temporal}} = \int_{\Omega} \int_{\mathcal{T}} \left\|\frac{\partial^2 \mathbf{W}}{\partial \tau^2}\right\|^2 d\mathbf{x} d\tau$$

For discrete time steps:
$$\mathcal{L}_{\text{temporal}} = \sum_{i} \|\mathbf{W}(\mathbf{x}, \tau_{i+1}) - 2\mathbf{W}(\mathbf{x}, \tau_i) + \mathbf{W}(\mathbf{x}, \tau_{i-1})\|^2$$

**4. Elastic Energy:**
Based on continuum mechanics, penalizes deformation energy:

$$\mathcal{L}_{\text{elastic}} = \int_{\Omega} \frac{\lambda}{2}(\text{tr}(\mathbf{E}))^2 + \mu \text{tr}(\mathbf{E}^2) d\mathbf{x}$$

where $\mathbf{E} = \frac{1}{2}(J_{\mathbf{W}}^T J_{\mathbf{W}} - I)$ is the Green strain tensor, and $\lambda, \mu$ are Lamé parameters.

**5. As-Rigid-As-Possible (ARAP):**
$$\mathcal{L}_{\text{ARAP}} = \sum_{i} \sum_{j \in \mathcal{N}(i)} w_{ij} \|(\mathbf{W}(\mathbf{x}_i) - \mathbf{W}(\mathbf{x}_j)) - \mathbf{R}_i(\mathbf{x}_i - \mathbf{x}_j)\|^2$$

where $\mathbf{R}_i$ is the optimal rotation at point $i$.

**Combined Regularization:**
$$\mathcal{L}_{\text{reg}} = \lambda_1 \mathcal{L}_{\text{smooth}} + \lambda_2 \mathcal{L}_{\text{rigid}} + \lambda_3 \mathcal{L}_{\text{temporal}} + \lambda_4 \mathcal{L}_{\text{volume}}$$

The weights $\lambda_i$ control the relative importance and depend on the scene type:
- Fluid scenes: High $\lambda_1$, low $\lambda_2$
- Articulated objects: High $\lambda_2$, moderate $\lambda_3$
- Soft bodies: Balanced weights

### 7.2.4 Bijective Mapping Guarantees

Ensuring that deformation fields are invertible (bijective) is crucial for stable optimization and physically meaningful results.

**1. Residual Formulation:**
Parameterize deformation as identity plus small displacement:

$$\mathbf{W}(\mathbf{x}, \tau) = \mathbf{x} + \epsilon \cdot \mathbf{u}(\mathbf{x}, \tau)$$

**Invertibility Condition:**
By the Banach fixed-point theorem, if:
$$\|\nabla \mathbf{u}\|_{\text{op}} < \frac{1}{\epsilon}$$

then $\mathbf{W}$ is a diffeomorphism (smooth bijection).

**Proof Sketch:**
The inverse satisfies: $\mathbf{y} = \mathbf{W}^{-1}(\mathbf{y}) + \epsilon \mathbf{u}(\mathbf{W}^{-1}(\mathbf{y}))$

Define operator $T(\mathbf{x}) = \mathbf{y} - \epsilon \mathbf{u}(\mathbf{x})$. If $\epsilon\|\nabla \mathbf{u}\|_{\text{op}} < 1$, then $T$ is a contraction, guaranteeing unique fixed point.

**2. Neural ODE Formulation:**
Model deformation through continuous flow:

$$\frac{d\mathbf{x}}{d\tau} = \mathbf{v}(\mathbf{x}, \tau), \quad \mathbf{x}(0) = \mathbf{x}_0$$

**Advantages:**
- Invertibility guaranteed by ODE theory (if $\mathbf{v}$ is Lipschitz)
- Inverse computed by integrating backward: $\frac{d\mathbf{x}}{d\tau} = -\mathbf{v}(\mathbf{x}, -\tau)$
- Natural temporal continuity

**Implementation:**
$$\mathbf{W}(\mathbf{x}_0, \tau) = \mathbf{x}_0 + \int_0^\tau \mathbf{v}(\mathbf{x}(s), s) ds$$

Use numerical ODE solvers (e.g., Runge-Kutta) with adaptive step size.

**3. Diffeomorphic Registration:**
Ensure smooth bijection through regularized velocity fields:

$$\mathcal{L}_{\text{diffeo}} = \int_0^T \int_{\Omega} \|\mathbf{L}\mathbf{v}(\mathbf{x}, t)\|^2 d\mathbf{x} dt$$

where $\mathbf{L}$ is a differential operator (e.g., $\mathbf{L} = \alpha\nabla^2 + \beta I$).

**4. Barrier Methods:**
Add penalty that approaches infinity as Jacobian determinant approaches zero:

$$\mathcal{L}_{\text{barrier}} = \int_{\Omega} \phi(\det(J_{\mathbf{W}})) d\mathbf{x}$$

where $\phi(s) = -\log(s)$ for $s > 0$ or $\phi(s) = \frac{1}{s} - 1$ for $s > 0$.

**Practical Guidelines:**
- Initialize near identity: $\mathbf{W}(\mathbf{x}, 0) = \mathbf{x}$
- Monitor $\min_{\mathbf{x}} \det(J_{\mathbf{W}})$ during training
- Use graduated optimization: start with strong regularization, gradually relax
- For neural networks, use spectral normalization to control Lipschitz constant

## 7.3 Optical Flow and Scene Flow Constraints

### 7.3.1 3D Scene Flow Definition

Scene flow $\mathbf{s}: \mathbb{R}^3 \times \mathbb{R} \rightarrow \mathbb{R}^3$ describes the 3D motion field:

$$\mathbf{s}(\mathbf{x}, \tau) = \frac{\partial \mathbf{W}(\mathbf{x}, \tau)}{\partial \tau}$$

### 7.3.2 2D Optical Flow from Volume Rendering

The 2D optical flow in image space can be derived from scene flow through the rendering equation:

$$\mathbf{u}_{\text{2D}}(\mathbf{p}) = \frac{\int_{t_n}^{t_f} T(t)\sigma(\mathbf{r}(t)) w(\mathbf{r}(t)) \Pi(\mathbf{s}(\mathbf{r}(t))) dt}{\int_{t_n}^{t_f} T(t)\sigma(\mathbf{r}(t)) w(\mathbf{r}(t)) dt}$$

where:
- $\Pi$ is the projection operator from 3D to 2D
- $w(\mathbf{x}) = \frac{\partial C}{\partial \mathbf{x}}$ is the contribution weight

### 7.3.3 Flow Consistency Equations

**Brightness Constancy:**
$$\frac{\partial I}{\partial \tau} + \nabla I \cdot \mathbf{u}_{\text{2D}} = 0$$

**Volume Rendering Consistency:**
$$\frac{\partial C(\mathbf{r}, \tau)}{\partial \tau} + \int_{t_n}^{t_f} T(t)\sigma(\mathbf{r}(t), \tau) \nabla_{\mathbf{x}} \mathbf{c} \cdot \mathbf{s} dt = 0$$

### 7.3.4 Multi-view Flow Constraints

For multiple cameras with centers $\{\mathbf{o}_i\}$ and directions $\{\mathbf{d}_i\}$:

**Epipolar Constraint on Flow:**
$$(\mathbf{u}_{\text{2D}}^{(i)} - \mathbf{u}_{\text{2D}}^{(j)})^T \mathbf{F}_{ij} \mathbf{p} = 0$$

where $\mathbf{F}_{ij}$ is the fundamental matrix between views $i$ and $j$.

**Depth-Flow Consistency:**
$$z_i \mathbf{u}_{\text{2D}}^{(i)} = \mathbf{K}_i \mathbf{R}_i \mathbf{s} + \frac{\partial z_i}{\partial \tau} \mathbf{K}_i \mathbf{R}_i \mathbf{d}_i$$

## 7.4 Warp-based Motion Representation

### 7.4.1 Space-Time Warping Functions

General warping function:
$$\Psi: \mathbb{R}^3 \times \mathbb{R} \rightarrow \mathbb{R}^3 \times \mathbb{R}$$
$$(\mathbf{x}', \tau') = \Psi(\mathbf{x}, \tau)$$

**Separable Warps:**
$$\Psi(\mathbf{x}, \tau) = (\mathbf{W}(\mathbf{x}, \tau), \tau)$$

**Non-separable Warps:**
$$\Psi(\mathbf{x}, \tau) = (\mathbf{W}(\mathbf{x}, \tau), T(\tau))$$

### 7.4.2 Hierarchical Motion Decomposition

$$\mathbf{W}(\mathbf{x}, \tau) = \mathbf{W}_{\text{global}}(\tau) \circ \mathbf{W}_{\text{local}}(\mathbf{x}, \tau)$$

where:
- $\mathbf{W}_{\text{global}}$: Rigid body motion (6 DOF)
- $\mathbf{W}_{\text{local}}$: Non-rigid deformation

**Lie Group Representation for Global Motion:**
$$\mathbf{W}_{\text{global}}(\tau) = \exp(\sum_{i=1}^{6} \theta_i(\tau) \mathbf{G}_i)$$

where $\{\mathbf{G}_i\}$ are generators of SE(3).

### 7.4.3 Temporal Interpolation Schemes

**Linear Interpolation:**
$$\mathbf{W}(\mathbf{x}, \tau) = (1-\alpha)\mathbf{W}(\mathbf{x}, \tau_i) + \alpha\mathbf{W}(\mathbf{x}, \tau_{i+1})$$

**Spherical Linear Interpolation (for rotations):**
$$\mathbf{R}(\tau) = \mathbf{R}_i (\mathbf{R}_i^T \mathbf{R}_{i+1})^{\alpha}$$

**Cubic Hermite Splines:**
$$\mathbf{W}(\mathbf{x}, \tau) = \sum_{j=0}^{3} h_j(\alpha) \mathbf{c}_j(\mathbf{x})$$

where $h_j$ are Hermite basis functions.

### 7.4.4 Inverse Warping Computation

**Fixed-Point Iteration:**
$$\mathbf{x}_{n+1} = \mathbf{x}_0 - \mathbf{W}(\mathbf{x}_n, \tau) + \mathbf{x}_{\text{target}}$$

Convergence guaranteed when:
$$\|\nabla_{\mathbf{x}} \mathbf{W}\|_{\text{op}} < 1$$

**Neural Inverse Networks:**
$$\mathbf{W}^{-1} = f_{\phi}(\mathbf{x}, \tau)$$

with cycle consistency loss:
$$\mathcal{L}_{\text{cycle}} = \|\mathbf{W}(\mathbf{W}^{-1}(\mathbf{x}, \tau), \tau) - \mathbf{x}\|^2$$

## 7.5 Canonical Space Methods

### 7.5.1 Canonical Space Definition

The canonical space represents a reference configuration of the scene:

$$\mathcal{C} = \{(\sigma_{\text{can}}(\mathbf{x}), \mathbf{c}_{\text{can}}(\mathbf{x}, \mathbf{d})) : \mathbf{x} \in \Omega_{\text{can}}\}$$

Common choices:
- **Rest pose**: For articulated objects
- **Mean shape**: $\mathbf{x}_{\text{can}} = \frac{1}{|\mathcal{T}|}\int_{\mathcal{T}} \mathbf{W}^{-1}(\mathbf{x}, \tau) d\tau$
- **Learned canonical**: Jointly optimized with deformations

### 7.5.2 Deformation to Observation Space

The rendering equation in canonical space:

$$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T_{\text{can}}(t')\sigma_{\text{can}}(\mathbf{r}_{\text{can}}(t'))\mathbf{c}_{\text{can}}(\mathbf{r}_{\text{can}}(t'), \mathbf{d}_{\text{can}}) \left|\frac{dt'}{dt}\right| dt$$

where:
- $t' = t'(t, \tau)$ is the warped ray parameter
- $\mathbf{r}_{\text{can}}(t') = \mathbf{W}^{-1}(\mathbf{r}(t), \tau)$
- $\mathbf{d}_{\text{can}} = \frac{d\mathbf{r}_{\text{can}}}{dt'} / \|\frac{d\mathbf{r}_{\text{can}}}{dt'}\|$

### 7.5.3 Template Learning Strategies

**Joint Optimization:**
$$\min_{\theta, \phi} \sum_{i,j} \|C_{ij} - \hat{C}_{ij}(\theta, \phi)\|^2 + \lambda \mathcal{R}(\phi)$$

where:
- $\theta$: Canonical radiance field parameters
- $\phi$: Deformation field parameters
- $\mathcal{R}$: Regularization term

**Progressive Training:**
1. Initialize with static reconstruction at $\tau_0$
2. Gradually add temporal samples
3. Refine canonical space and deformations jointly

### 7.5.4 Handling Topological Changes

For scenes with changing topology (e.g., fluids), we extend the canonical space:

$$\sigma_{\text{can}}(\mathbf{x}) = \sigma_{\text{base}}(\mathbf{x}) + \sum_{k} \alpha_k(\tau) \sigma_{\text{residual}}^{(k)}(\mathbf{x})$$

where $\alpha_k(\tau) \in [0,1]$ controls the appearance/disappearance of components.

### 7.5.5 Regularization in Canonical Space

**Spatial Smoothness:**
$$\mathcal{L}_{\text{can-smooth}} = \int_{\Omega_{\text{can}}} \|\nabla \sigma_{\text{can}}\|^2 + \|\nabla \mathbf{c}_{\text{can}}\|^2 d\mathbf{x}$$

**Sparsity Prior:**
$$\mathcal{L}_{\text{sparse}} = \int_{\Omega_{\text{can}}} \log(1 + \sigma_{\text{can}}^2/\epsilon^2) d\mathbf{x}$$

**Deformation Regularization:**
$$\mathcal{L}_{\text{def-reg}} = \int_{\mathcal{T}} \int_{\Omega} \|\mathbf{W}(\mathbf{x}, \tau) - \mathbf{x}\|^2 \rho(\mathbf{x}) d\mathbf{x} d\tau$$

where $\rho(\mathbf{x})$ weights regions (e.g., higher for rigid parts).

## 本章小结

Dynamic Neural Radiance Fields extend the static formulation by incorporating time as an additional dimension. The key mathematical frameworks include:

1. **Time-Extended Volume Rendering Equation:**
   $$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T(t, \tau)\sigma(\mathbf{r}(t), \tau)\mathbf{c}(\mathbf{r}(t), \mathbf{d}, \tau) dt$$

2. **Deformation Field Formulation:**
   $$\mathbf{x}_{\text{obs}} = \mathbf{W}(\mathbf{x}_{\text{can}}, \tau)$$
   
   with volume preservation: $\det(J_{\mathbf{W}}) = 1$

3. **Scene Flow Consistency:**
   $$\mathbf{s}(\mathbf{x}, \tau) = \frac{\partial \mathbf{W}(\mathbf{x}, \tau)}{\partial \tau}$$

4. **Canonical Space Rendering:**
   $$C(\mathbf{r}, \tau) = \int T_{\text{can}}\sigma_{\text{can}}(\mathbf{W}^{-1}(\mathbf{r}(t), \tau))\mathbf{c}_{\text{can}} \left|\det J_{\mathbf{W}^{-1}}\right| dt$$

5. **Regularization Framework:**
   - Spatial smoothness: $\mathcal{L}_{\text{smooth}}$
   - Temporal coherence: $\mathcal{L}_{\text{temporal}}$
   - Rigidity constraints: $\mathcal{L}_{\text{rigid}}$

These formulations provide a principled way to model dynamic scenes while maintaining the elegance of the volume rendering framework.

## 练习题

### 基础题

**习题 7.1**: 时间采样分析
考虑一个以频率 $f = 2$ Hz 旋转的物体。如果我们以 30 FPS 采样，是否满足 Nyquist 条件？如果物体同时包含 5 Hz 的振动，需要什么采样率？

*提示*: Nyquist 频率 $f_s \geq 2f_{\max}$

<details>
<summary>答案</summary>

对于 2 Hz 旋转：
- Nyquist 要求: $f_s \geq 2 \times 2 = 4$ Hz
- 30 FPS 采样率满足条件 (30 > 4)

对于 5 Hz 振动：
- Nyquist 要求: $f_s \geq 2 \times 5 = 10$ Hz  
- 需要至少 10 FPS

组合运动：$f_{\max} = 5$ Hz，需要至少 10 FPS 采样率。

</details>

**习题 7.2**: 体积保持约束
给定 2D 变形场 $\mathbf{W}(x,y) = (x + 0.1\sin(y), y + 0.1\cos(x))$，验证这是否满足体积保持约束。

*提示*: 计算 Jacobian 行列式

<details>
<summary>答案</summary>

Jacobian 矩阵：
$$J = \begin{bmatrix}
\frac{\partial W_x}{\partial x} & \frac{\partial W_x}{\partial y} \\
\frac{\partial W_y}{\partial x} & \frac{\partial W_y}{\partial y}
\end{bmatrix} = \begin{bmatrix}
1 & 0.1\cos(y) \\
-0.1\sin(x) & 1
\end{bmatrix}$$

行列式：
$$\det(J) = 1 \times 1 - (0.1\cos(y))(-0.1\sin(x)) = 1 + 0.01\sin(x)\cos(y)$$

由于 $\det(J) \neq 1$，不满足体积保持约束。最大偏差为 0.01。

</details>

**习题 7.3**: 光流计算
给定场景流 $\mathbf{s}(x,y,z) = (0, 0, -v)$ (沿 z 轴的匀速运动) 和相机内参矩阵 $\mathbf{K} = \text{diag}(f, f, 1)$，推导 2D 光流。

*提示*: 使用投影关系 $u = fx/z$, $v = fy/z$

<details>
<summary>答案</summary>

3D 点投影：$(u, v) = (fx/z, fy/z)$

时间导数：
$$\frac{du}{dt} = f \frac{d}{dt}\left(\frac{x}{z}\right) = f \frac{\dot{x}z - x\dot{z}}{z^2}$$

代入 $\dot{x} = 0$, $\dot{z} = -v$：
$$\frac{du}{dt} = f \frac{xv}{z^2} = \frac{uv}{z}$$

类似地：
$$\frac{dv}{dt} = \frac{vv}{z}$$

2D 光流：$\mathbf{u}_{2D} = \frac{v}{z}(u, v)$ - 径向扩张模式。

</details>

### 挑战题

**习题 7.4**: 变形场设计
设计一个变形场 $\mathbf{W}(\mathbf{x}, \tau)$ 来表示心跳运动，满足：
1. 周期为 $T = 1$ 秒
2. 体积在 $\pm 10\%$ 范围内变化
3. 运动主要沿径向

*提示*: 考虑球坐标系和时变缩放

<details>
<summary>答案</summary>

球坐标表示：$\mathbf{x} = (r, \theta, \phi)$

径向缩放变形：
$$\mathbf{W}(r, \theta, \phi, \tau) = s(\tau) \cdot r \cdot \hat{\mathbf{r}}$$

其中缩放因子：
$$s(\tau) = 1 + 0.1\sin(2\pi \tau/T)$$

笛卡尔坐标：
$$\mathbf{W}(\mathbf{x}, \tau) = s(\tau) \cdot \mathbf{x}$$

体积变化：
$$\det(J_{\mathbf{W}}) = s(\tau)^3 \approx 1 + 0.3\sin(2\pi\tau/T)$$

满足 $\pm 10\%$ 体积变化（一阶近似）。

</details>

**习题 7.5**: 正则化分析
证明对于小变形 $\mathbf{W}(\mathbf{x}) = \mathbf{x} + \epsilon \mathbf{u}(\mathbf{x})$，如果 $\|\nabla \mathbf{u}\|_{op} < 1/\epsilon$，则变形是可逆的。

*提示*: 使用 Banach 不动点定理

<details>
<summary>答案</summary>

求逆映射需要解：$\mathbf{y} = \mathbf{x} + \epsilon \mathbf{u}(\mathbf{x})$

定义算子：$T(\mathbf{x}) = \mathbf{y} - \epsilon \mathbf{u}(\mathbf{x})$

对于两点 $\mathbf{x}_1, \mathbf{x}_2$：
$$\|T(\mathbf{x}_1) - T(\mathbf{x}_2)\| = \epsilon\|\mathbf{u}(\mathbf{x}_1) - \mathbf{u}(\mathbf{x}_2)\|$$

由中值定理：
$$\|\mathbf{u}(\mathbf{x}_1) - \mathbf{u}(\mathbf{x}_2)\| \leq \|\nabla \mathbf{u}\|_{op} \|\mathbf{x}_1 - \mathbf{x}_2\|$$

因此：
$$\|T(\mathbf{x}_1) - T(\mathbf{x}_2)\| \leq \epsilon \|\nabla \mathbf{u}\|_{op} \|\mathbf{x}_1 - \mathbf{x}_2\|$$

当 $\epsilon \|\nabla \mathbf{u}\|_{op} < 1$ 时，$T$ 是压缩映射，由 Banach 定理存在唯一不动点，即逆映射存在。

</details>

**习题 7.6**: 多视图一致性
推导两个视图之间的场景流一致性约束。给定两个相机的投影矩阵 $\mathbf{P}_1, \mathbf{P}_2$，以及 3D 场景流 $\mathbf{s}$。

*提示*: 利用极线几何和流的投影关系

<details>
<summary>答案</summary>

3D 点 $\mathbf{X}$ 在两个视图的投影：
- 视图 1: $\mathbf{x}_1 = \mathbf{P}_1 \mathbf{X}$
- 视图 2: $\mathbf{x}_2 = \mathbf{P}_2 \mathbf{X}$

场景流后：$\mathbf{X}' = \mathbf{X} + \mathbf{s}$

投影变化：
- $\mathbf{x}'_1 = \mathbf{P}_1(\mathbf{X} + \mathbf{s})$
- $\mathbf{x}'_2 = \mathbf{P}_2(\mathbf{X} + \mathbf{s})$

2D 流：
- $\mathbf{u}_1 = \mathbf{x}'_1 - \mathbf{x}_1 = \mathbf{P}_1 \mathbf{s}$
- $\mathbf{u}_2 = \mathbf{x}'_2 - \mathbf{x}_2 = \mathbf{P}_2 \mathbf{s}$

消除 $\mathbf{s}$：如果 $\mathbf{P}_1^+$ 是 $\mathbf{P}_1$ 的伪逆，则：
$$\mathbf{u}_2 = \mathbf{P}_2 \mathbf{P}_1^+ \mathbf{u}_1$$

这建立了两视图光流的线性关系。

</details>

**习题 7.7**: 时间插值优化
给定离散时间点 $\{t_i\}$ 的变形场样本 $\{\mathbf{W}_i\}$，设计一个保持 $C^2$ 连续性的插值方案，并分析计算复杂度。

*提示*: 考虑三次样条插值

<details>
<summary>答案</summary>

三次 Hermite 样条插值：

对于 $t \in [t_i, t_{i+1}]$，令 $s = (t - t_i)/(t_{i+1} - t_i)$：

$$\mathbf{W}(t) = h_{00}(s)\mathbf{W}_i + h_{10}(s)\mathbf{m}_i\Delta t + h_{01}(s)\mathbf{W}_{i+1} + h_{11}(s)\mathbf{m}_{i+1}\Delta t$$

其中 Hermite 基函数：
- $h_{00}(s) = 2s^3 - 3s^2 + 1$
- $h_{10}(s) = s^3 - 2s^2 + s$
- $h_{01}(s) = -2s^3 + 3s^2$
- $h_{11}(s) = s^3 - s^2$

切线 $\mathbf{m}_i$ 通过解三对角系统获得（保证 $C^2$ 连续）：

$$\mathbf{m}_{i-1} + 4\mathbf{m}_i + \mathbf{m}_{i+1} = \frac{3}{\Delta t}(\mathbf{W}_{i+1} - \mathbf{W}_{i-1})$$

复杂度：
- 预计算：$O(n)$ 解三对角系统
- 查询：$O(\log n)$ 二分查找 + $O(1)$ 插值
- 空间：$O(n)$ 存储切线

</details>

**习题 7.8**: 拓扑变化建模
设计一个数学框架来处理场景中物体的出现和消失（例如：倒水时水的体积变化）。要求保持时间连续性。

*提示*: 考虑使用软掩码和连续变化的密度场

<details>
<summary>答案</summary>

时变密度场分解：
$$\sigma(\mathbf{x}, t) = \sigma_{\text{static}}(\mathbf{x}) + \sum_k \alpha_k(t) \sigma_k(\mathbf{x})$$

软出现/消失函数：
$$\alpha_k(t) = \begin{cases}
0 & t < t_k^{\text{start}} \\
\frac{1}{2}(1 - \cos(\pi \frac{t - t_k^{\text{start}}}{\Delta t})) & t \in [t_k^{\text{start}}, t_k^{\text{start}} + \Delta t] \\
1 & t \in [t_k^{\text{start}} + \Delta t, t_k^{\text{end}}] \\
\frac{1}{2}(1 + \cos(\pi \frac{t - t_k^{\text{end}}}{\Delta t})) & t \in [t_k^{\text{end}}, t_k^{\text{end}} + \Delta t] \\
0 & t > t_k^{\text{end}} + \Delta t
\end{cases}$$

体积守恒（对于流体）：
$$\int_{\Omega} \sigma(\mathbf{x}, t) d\mathbf{x} = V_0 + \int_0^t Q(\tau) d\tau$$

其中 $Q(t)$ 是流入/流出率。

渲染方程保持不变，自然处理拓扑变化：
$$C(\mathbf{r}, t) = \int T(s,t) \sigma(\mathbf{r}(s), t) \mathbf{c}(\mathbf{r}(s), \mathbf{d}, t) ds$$

</details>

## 常见陷阱与错误 (Gotchas)

### 1. 时间混叠问题

**问题**: 快速运动导致时间欠采样，产生伪影
```
症状：闪烁、运动模糊、鬼影
原因：违反 Nyquist 采样定理
```

**解决方案**:
- 增加时间采样密度
- 使用运动模糊建模
- 时间超分辨率技术

### 2. 非双射变形

**问题**: 变形场产生自相交或折叠
```
症状：渲染伪影、优化不收敛
检测：det(J) < 0 或 det(J) >> 1
```

**解决方案**:
- 限制变形幅度：$\|\mathbf{W} - \mathbf{x}\| < \epsilon$
- 使用残差参数化
- 添加 Jacobian 正则化

### 3. 规范空间歧义

**问题**: 多个规范空间-变形对产生相同观察
```
例子：膨胀的小球 vs 静止的大球
影响：优化陷入局部最优
```

**解决方案**:
- 强先验约束（如最小变形）
- 多阶段优化策略
- 利用多视图约束

### 4. 计算复杂度爆炸

**问题**: 4D 表示导致内存和计算需求剧增
```
静态 NeRF: O(XYZ)
动态 NeRF: O(XYZT)
```

**解决方案**:
- 分解表示（如 K-Planes）
- 自适应采样
- 增量式训练

### 5. 时间插值伪影

**问题**: 离散时间采样之间的插值不自然
```
线性插值：运动不流畅
高阶插值：可能过拟合
```

**解决方案**:
- 物理约束的插值
- 学习式时间编码
- 运动先验正则化

### 6. 光流-几何不一致

**问题**: 2D 光流与 3D 场景流投影不匹配
```
原因：遮挡、透明度变化
表现：优化震荡
```

**解决方案**:
- 遮挡感知的流计算
- 鲁棒损失函数
- 分层优化策略

## 最佳实践检查清单

### 设计阶段
- [ ] **时间建模选择**
  - [ ] 分析场景运动特性（周期性？局部性？）
  - [ ] 估算最大运动频率
  - [ ] 选择合适的时间表示（连续/离散/混合）

- [ ] **变形场设计**
  - [ ] 确定是否需要拓扑变化支持
  - [ ] 选择参数化方式（前向/后向/双向）
  - [ ] 设计正则化策略

- [ ] **内存预算规划**
  - [ ] 估算 4D 表示的内存需求
  - [ ] 考虑分解或压缩方案
  - [ ] 规划增量训练策略

### 实现阶段
- [ ] **数值稳定性**
  - [ ] 检查 Jacobian 计算的数值稳定性
  - [ ] 实现梯度裁剪
  - [ ] 添加 epsilon 防止除零

- [ ] **采样策略**
  - [ ] 实现自适应时间采样
  - [ ] 平衡空间和时间分辨率
  - [ ] 考虑重要性采样

- [ ] **优化流程**
  - [ ] 从粗到细的训练策略
  - [ ] 交替优化规范空间和变形
  - [ ] 监控各项损失的平衡

### 验证阶段
- [ ] **质量指标**
  - [ ] 时间一致性度量
  - [ ] 运动平滑度评估
  - [ ] 多视图一致性检查

- [ ] **性能分析**
  - [ ] 渲染时间 vs 质量权衡
  - [ ] 内存使用监控
  - [ ] 训练收敛速度

- [ ] **鲁棒性测试**
  - [ ] 快速运动场景
  - [ ] 大变形情况
  - [ ] 稀疏视图条件

### 调试技巧
- [ ] **可视化工具**
  - [ ] 变形场可视化（箭头图/网格）
  - [ ] Jacobian 热力图
  - [ ] 时间切片对比

- [ ] **诊断输出**
  - [ ] 逐帧光流检查
  - [ ] 规范空间稳定性
  - [ ] 正则化项贡献分析

- [ ] **消融实验**
  - [ ] 移除时间维度验证静态情况
  - [ ] 固定变形场测试规范空间
  - [ ] 简化场景逐步增加复杂度

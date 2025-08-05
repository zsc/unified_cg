# 第5章：基于物理的渲染

本章深入探讨基于物理的渲染（PBR）的数学基础，将光与物质的相互作用表述为辐射传输方程。我们将建立统一的体积渲染框架来描述散射介质、次表面散射和参与介质，并严格证明能量守恒与互易性原理在渲染中的应用。通过将传统的表面渲染扩展到体积渲染，本章为后续的神经渲染方法奠定物理基础。

## 5.1 辐射传输方程

辐射传输方程（Radiative Transfer Equation, RTE）是描述光在参与介质中传播的基本方程。它统一了吸收、发射和散射过程，为基于物理的体积渲染提供了数学框架。

### 5.1.1 基本形式推导

考虑光在位置 **x** 沿方向 **ω** 传播时的辐射亮度 L(**x**, **ω**, t)。在微小距离 ds 内，辐射亮度的变化由四个物理过程决定：

1. **吸收（Absorption）**：介质吸收光能，减少辐射亮度
2. **外散射（Out-scattering）**：光被散射到其他方向
3. **发射（Emission）**：介质自身发光
4. **内散射（In-scattering）**：其他方向的光被散射到当前方向

**从微观到宏观的严格推导**

考虑体积元 dV = dA·ds，其中 dA 是垂直于传播方向的截面积。时间 dt 内通过该体积元的能量变化：

1. **入射能量**：L(**x**, **ω**) dA dω dt
2. **出射能量**：L(**x** + ds**ω**, **ω**) dA dω dt
3. **吸收损失**：σₐ(**x**) L(**x**, **ω**) dV dω dt
4. **散射损失**：σₛ(**x**) L(**x**, **ω**) dV dω dt
5. **发射增益**：σₐ(**x**) Lₑ(**x**, **ω**) dV dω dt
6. **散射增益**：σₛ(**x**) ∫_{S²} p(**ω'** → **ω**) L(**x**, **ω'**) dω' dV dω dt

能量守恒要求：
$$[L(\mathbf{x}, \boldsymbol{\omega}) - L(\mathbf{x} + ds\boldsymbol{\omega}, \boldsymbol{\omega})] dA = [-\sigma_t L + \sigma_a L_e + \sigma_s \int p L' d\omega'] dA \cdot ds$$

取极限 ds → 0，得到微分形式：

$$\frac{dL(\mathbf{x}, \boldsymbol{\omega})}{ds} = -\sigma_t(\mathbf{x}) L(\mathbf{x}, \boldsymbol{\omega}) + \sigma_a(\mathbf{x}) L_e(\mathbf{x}, \boldsymbol{\omega}) + \sigma_s(\mathbf{x}) \int_{S^2} p(\mathbf{x}, \boldsymbol{\omega}' \to \boldsymbol{\omega}) L(\mathbf{x}, \boldsymbol{\omega}') d\boldsymbol{\omega}'$$

其中：
- σₜ = σₐ + σₛ：消光系数（extinction coefficient）
- σₐ：吸收系数（absorption coefficient）
- σₛ：散射系数（scattering coefficient）
- Lₑ：发射项（emission term）
- p(**ω'** → **ω**)：相位函数（phase function）

**方向导数的几何解释**

射线参数化：**x**(s) = **x₀** + s**ω**，则：
$$\frac{dL}{ds} = \frac{\partial L}{\partial x_i} \frac{dx_i}{ds} = \frac{\partial L}{\partial x_i} \omega_i = \boldsymbol{\omega} \cdot \nabla L$$

这将沿射线的变化率转换为空间导数。

**辐射传输方程的标准形式**

$$(\boldsymbol{\omega} \cdot \nabla) L(\mathbf{x}, \boldsymbol{\omega}) + \sigma_t(\mathbf{x}) L(\mathbf{x}, \boldsymbol{\omega}) = \sigma_a(\mathbf{x}) L_e(\mathbf{x}, \boldsymbol{\omega}) + \sigma_s(\mathbf{x}) \int_{S^2} p(\mathbf{x}, \boldsymbol{\omega}' \to \boldsymbol{\omega}) L(\mathbf{x}, \boldsymbol{\omega}') d\boldsymbol{\omega}'$$

**时间相关形式**

对于非稳态问题，需要考虑时间导数：
$$\frac{1}{c} \frac{\partial L}{\partial t} + (\boldsymbol{\omega} \cdot \nabla) L + \sigma_t L = \sigma_a L_e + \sigma_s \int_{S^2} p L' d\boldsymbol{\omega}'$$

其中 c 是介质中的光速。

**物理量的尺度分析**

为了更好理解RTE中各项的相对重要性，考虑典型尺度：
- 特征长度 L₀：场景或介质的典型尺寸
- 平均自由程 ℓ = 1/σₜ：光子碰撞间的平均距离
- 光学厚度 τ = L₀/ℓ = σₜL₀：无量纲参数

根据τ的大小，介质分类为：
- **光学薄介质**（τ ≪ 1）：单次散射主导，可用Born近似
- **光学厚介质**（τ ≫ 1）：多次散射主导，适用扩散近似
- **中等光学厚度**（τ ≈ 1）：需要完整RTE求解

**RTE的算子形式**

定义传输算子 𝒯 和散射算子 𝒮：
$$\mathcal{T} = -(\boldsymbol{\omega} \cdot \nabla) - \sigma_t$$
$$\mathcal{S}[L](\mathbf{x}, \boldsymbol{\omega}) = \sigma_s(\mathbf{x}) \int_{S^2} p(\mathbf{x}, \boldsymbol{\omega}' \to \boldsymbol{\omega}) L(\mathbf{x}, \boldsymbol{\omega}') d\boldsymbol{\omega}'$$

RTE可写作：
$$\mathcal{T}[L] = -\sigma_a L_e - \mathcal{S}[L]$$

或重排为：
$$L = \mathcal{T}^{-1}[-\sigma_a L_e - \mathcal{S}[L]]$$

这形式便于迭代求解和理论分析。

### 5.1.2 边界条件与初始条件

对于有界域 Ω ⊂ ℝ³，需要指定边界条件。设 ∂Ω 为域的边界，**n** 为外法向量。

**入射边界条件**：对于 **x** ∈ ∂Ω 且 **ω** · **n** < 0：
$$L(\mathbf{x}, \boldsymbol{\omega}) = L_{inc}(\mathbf{x}, \boldsymbol{\omega})$$

**出射边界条件**：对于 **x** ∈ ∂Ω 且 **ω** · **n** > 0：
- 真空边界：无约束
- 反射边界：$$L(\mathbf{x}, \boldsymbol{\omega}) = \int_{\boldsymbol{\omega}' \cdot \mathbf{n} < 0} f_r(\mathbf{x}, \boldsymbol{\omega}', \boldsymbol{\omega}) L(\mathbf{x}, \boldsymbol{\omega}') |\boldsymbol{\omega}' \cdot \mathbf{n}| d\boldsymbol{\omega}'$$

**界面条件**：在两种介质的界面上，需要考虑菲涅尔方程。设折射率为 n₁ 和 n₂，则：
$$L_t = \frac{n_2^2}{n_1^2} T(\theta_i) L_i$$
$$L_r = R(\theta_i) L_i$$

其中 T 和 R 分别为透射率和反射率，满足 T + R = 1（能量守恒）。

**Robin边界条件**

对于半透明边界，使用混合边界条件：
$$L(\mathbf{x}, \boldsymbol{\omega}) + \ell_b (\boldsymbol{\omega} \cdot \mathbf{n}) \frac{\partial L}{\partial n} = g(\mathbf{x}, \boldsymbol{\omega})$$

其中 ℓ_b 是外推长度，与表面反射特性相关。

**周期边界条件**

对于周期性结构（如光子晶体）：
$$L(\mathbf{x} + \mathbf{L}, \boldsymbol{\omega}) = L(\mathbf{x}, \boldsymbol{\omega})$$

其中 **L** 是晶格周期向量。

### 5.1.3 解析解与数值方法

辐射传输方程的解析解仅在特殊情况下存在。我们首先考虑几种可解情况，然后讨论一般的数值方法。

**情况1：均匀介质无散射**（σₛ = 0）

沿射线参数化：**x**(s) = **x₀** + s**ω**，方程简化为常微分方程：
$$\frac{dL}{ds} + \sigma_a L = \sigma_a L_e$$

这是一阶线性ODE，使用积分因子μ(s) = e^{σₐs}：
$$\frac{d}{ds}[L(s)e^{\sigma_a s}] = \sigma_a L_e(s) e^{\sigma_a s}$$

积分得到通解：
$$L(s) = L(0)e^{-\sigma_a s} + \sigma_a \int_0^s L_e(s') e^{-\sigma_a(s-s')} ds'$$

这就是比尔-朗伯定律（Beer-Lambert law）的推广，包含了发射项。

**情况2：各向同性单次散射**

假设散射仅发生一次，且相位函数各向同性 p = 1/(4π)。使用Born近似：
$$L = L^{(0)} + L^{(1)} + \mathcal{O}(\sigma_s^2)$$

其中 L^{(0)} 是无散射解，L^{(1)} 是单次散射贡献：
$$L^{(1)}(\mathbf{x}, \boldsymbol{\omega}) = \int_0^s T(0,s') \frac{\sigma_s(\mathbf{x}')}{4\pi} \int_{S^2} L^{(0)}(\mathbf{x}(s'), \boldsymbol{\omega}') d\boldsymbol{\omega}' ds'$$

其中 T(s₁,s₂) = exp(-∫_{s₁}^{s₂} σₜ(s')ds') 为光学透射率。

**情况3：小角度散射近似**

当相位函数强烈前向散射（g → 1）时，使用Fokker-Planck近似：
$$\frac{\partial L}{\partial s} + \sigma_t L = \sigma_a L_e + \sigma_s L + \frac{\sigma_s(1-g)}{2} \nabla_\perp^2 L$$

其中 ∇_⊥² 是垂直于 **ω** 的拉普拉斯算子。

**积分方程形式**

将RTE改写为Volterra型积分方程，便于迭代求解：
$$L(\mathbf{x}, \boldsymbol{\omega}) = L_0(\mathbf{x}, \boldsymbol{\omega}) + \int_0^{\tau_b} e^{-\tau} Q(\mathbf{x}(\tau), \boldsymbol{\omega}) d\tau$$

其中：
- τ = ∫₀ˢ σₜ(s')ds' 是光学深度
- Q = σₐLₑ + σₛ∫p L dω' 是源项
- L₀ 是边界贡献

**特征线方法**

RTE沿特征线（光线）是一维问题。定义特征坐标：
$$\frac{d\mathbf{x}}{ds} = \boldsymbol{\omega}, \quad \frac{d\boldsymbol{\omega}}{ds} = 0$$

沿特征线，RTE简化为：
$$\frac{dL}{ds} + \sigma_t L = S(\mathbf{x}(s), \boldsymbol{\omega})$$

这可用标准ODE方法求解。

**数值方法详述**

1. **离散坐标法（Discrete Ordinates, S_N方法）**：
   
   将角度空间离散化为 N 个方向 {**ωᵢ**}，权重 {wᵢ}满足：
   $$\sum_{i=1}^N w_i = 4\pi, \quad \sum_{i=1}^N w_i \omega_{i,k} = 0, \quad \sum_{i=1}^N w_i \omega_{i,k}\omega_{i,l} = \frac{4\pi}{3}\delta_{kl}$$
   
   离散化的RTE系统：
   $$(\boldsymbol{\omega}_i \cdot \nabla) L_i + \sigma_t L_i = \sigma_a L_{e,i} + \sigma_s \sum_{j=1}^N w_j p_{ij} L_j$$
   
   常用正交集：Gauss-Legendre、Chebyshev等。

2. **球谐函数展开（P_N方法）**：
   
   将辐射亮度展开为球谐函数：
   $$L(\mathbf{x}, \boldsymbol{\omega}) = \sum_{l=0}^{\infty} \sum_{m=-l}^{l} L_{lm}(\mathbf{x}) Y_{lm}(\boldsymbol{\omega})$$
   
   其中球谐系数：
   $$L_{lm}(\mathbf{x}) = \int_{S^2} L(\mathbf{x}, \boldsymbol{\omega}) Y_{lm}^*(\boldsymbol{\omega}) d\boldsymbol{\omega}$$
   
   代入RTE，利用球谐函数的递推关系和正交性，得到耦合PDE系统。P₁近似保留l=0,1项，导出扩散方程。

3. **蒙特卡洛方法（Monte Carlo）**：
   
   路径积分表示：
   $$L(\mathbf{x}, \boldsymbol{\omega}) = \sum_{k=0}^{\infty} L^{(k)}(\mathbf{x}, \boldsymbol{\omega})$$
   
   其中 L^{(k)} 是k次散射贡献。使用俄罗斯轮盘终止无限级数：
   $$L \approx \sum_{k=0}^{K} \frac{L^{(k)}}{p_k}$$
   
   其中 pₖ 是继续追踪的概率。

   方差减少技术：
   - **重要性采样**：根据贡献大小调整采样分布
   - **分层采样**：将样本空间分层以减少聚集
   - **控制变量**：使用已知期望的相关变量

4. **有限元方法（FEM）**：
   
   弱形式：
   $$\int_{\Omega \times S^2} \psi \left[ (\boldsymbol{\omega} \cdot \nabla) L + \sigma_t L \right] d\mathbf{x} d\boldsymbol{\omega} = \int_{\Omega \times S^2} \psi S d\mathbf{x} d\boldsymbol{\omega}$$
   
   选择适当的测试函数ψ和基函数，构建稀疏线性系统。

   **流线上风Petrov-Galerkin（SUPG）**：
   对于对流占优问题，添加人工扩散项：
   $$\psi_{SUPG} = \psi + \tau (\boldsymbol{\omega} \cdot \nabla)\psi$$
   
   其中 τ 是稳定化参数。

5. **格子Boltzmann方法**：
   
   基于速度离散的动力学方程：
   $$f_i(t+\Delta t, \mathbf{x}+\boldsymbol{\omega}_i\Delta t) - f_i(t,\mathbf{x}) = \Omega_i$$
   
   其中 Ω_i 是碰撞算子，f_i 是分布函数。

**收敛性分析**

对于迭代方法，收敛条件是散射算子的谱半径小于1：
$$\rho(\mathcal{K}) = \sup_{\mathbf{x}} \frac{\sigma_s(\mathbf{x})}{\sigma_t(\mathbf{x})} < 1$$

收敛速度：
$$||L^{(n)} - L|| \leq C \rho^n ||L^{(0)} - L||$$

**误差估计**

对于离散化方法，总误差包括：
- **截断误差**：O(h^p)，其中h是网格尺寸，p是方法阶数
- **统计误差**（蒙特卡洛）：O(1/√N)，N是样本数
- **射线效应**（S_N方法）：在光学薄区域的非物理条纹

## 5.2 体积散射与相位函数

体积散射是光与介质中微粒相互作用的结果。相位函数描述了散射光的角度分布，是参与介质渲染的核心组件。

### 5.2.1 散射系数与吸收系数

**微观到宏观的推导**

考虑介质中密度为 ρ(x) 的散射粒子，每个粒子的散射截面为 σₛ，吸收截面为 σₐ。宏观系数为：

$$\sigma_s(\mathbf{x}) = \rho(\mathbf{x}) \cdot \sigma_{s,particle}$$
$$\sigma_a(\mathbf{x}) = \rho(\mathbf{x}) \cdot \sigma_{a,particle}$$

对于混合介质，总系数是各成分的线性组合：
$$\sigma_t = \sum_i \rho_i \sigma_{t,i}$$

**散射截面的物理含义**

散射截面定义为入射能量流与散射功率的比值：
$$\sigma_s = \frac{P_{scattered}}{I_{incident}}$$

对于球形粒子，几何截面 σ_geo = πr²，但散射截面可能大于或小于几何截面：
- **散射效率因子** Q_s = σ_s/σ_geo
- 小粒子（Rayleigh）：Q_s ∝ r⁴/λ⁴
- 大粒子（几何光学）：Q_s → 2（衍射贡献）

**单散射反照率（Single Scattering Albedo）**

定义单散射反照率 α 为散射与消光的比值：
$$\alpha = \frac{\sigma_s}{\sigma_t} = \frac{\sigma_s}{\sigma_a + \sigma_s}$$

物理意义：光子与介质相互作用时被散射（而非吸收）的概率。
- α = 0：纯吸收介质
- α = 1：纯散射介质（无吸收）
- 0 < α < 1：一般情况

**平均自由程（Mean Free Path）**

光子在介质中两次相互作用之间的平均距离：
$$\ell = \frac{1}{\sigma_t}$$

概率分布遵循指数分布：
$$P(s) = \sigma_t e^{-\sigma_t s}$$

相关长度尺度：
- **散射平均自由程**：ℓ_s = 1/σ_s
- **吸收平均自由程**：ℓ_a = 1/σ_a
- **传输平均自由程**：ℓ* = 1/σ_t' = 1/[σ_a + σ_s(1-g)]

**波长依赖性**

光学系数通常与波长相关：
- **Rayleigh散射**：σ_s ∝ λ⁻⁴（蓝天现象）
- **Mie散射**：复杂振荡行为
- **吸收**：取决于材质的电子跃迁和振动模式

### 5.2.2 相位函数的数学性质

相位函数 p(**ω'** → **ω**) 描述了光从方向 **ω'** 散射到方向 **ω** 的概率密度。

**归一化条件**
$$\int_{S^2} p(\boldsymbol{\omega}' \to \boldsymbol{\omega}) d\boldsymbol{\omega} = 1$$

**对称性**

1. **旋转不变性**：大多数介质中，相位函数仅依赖于散射角 θ = arccos(**ω'** · **ω**)：
   $$p(\boldsymbol{\omega}' \to \boldsymbol{\omega}) = p(\cos\theta)$$

2. **互易性**（某些情况下）：
   $$p(\boldsymbol{\omega}' \to \boldsymbol{\omega}) = p(\boldsymbol{\omega} \to \boldsymbol{\omega}')$$

3. **前后对称性**（某些粒子）：
   $$p(\cos\theta) = p(\cos(\pi - \theta))$$

**矩与各向异性参数**

第 n 阶矩：
$$\mu_n = \int_{-1}^{1} \cos^n\theta \cdot p(\cos\theta) d(\cos\theta)$$

各向异性参数 g（平均余弦）：
$$g = \mu_1 = \int_{S^2} (\boldsymbol{\omega}' \cdot \boldsymbol{\omega}) p(\boldsymbol{\omega}' \to \boldsymbol{\omega}) d\boldsymbol{\omega}$$

- g > 0：前向散射为主
- g = 0：各向同性散射
- g < 0：后向散射为主

**Legendre多项式展开**

相位函数可展开为Legendre多项式：
$$p(\cos\theta) = \sum_{l=0}^{\infty} \frac{2l+1}{4\pi} \chi_l P_l(\cos\theta)$$

其中展开系数：
$$\chi_l = 2\pi \int_{-1}^{1} p(\cos\theta) P_l(\cos\theta) d(\cos\theta)$$

注意：χ₀ = 1（归一化），χ₁ = g（各向异性参数）

**相位函数的信息熵**

定义Shannon熵量化散射的随机性：
$$H[p] = -\int_{S^2} p(\boldsymbol{\omega}) \ln p(\boldsymbol{\omega}) d\boldsymbol{\omega}$$

- 各向同性：H_max = ln(4π)
- 完全前向：H_min = 0
- 实际介质：0 < H < ln(4π)

### 5.2.3 常见相位函数模型

**1. 各向同性散射**
$$p(\cos\theta) = \frac{1}{4\pi}$$

特点：g = 0，数学最简单，适用于高度多次散射的情况。

采样：cosθ = 1 - 2ξ₁, φ = 2πξ₂

应用场景：
- 密集烟雾的多次散射
- 扩散近似的验证
- 理论分析的基准

**2. Rayleigh 散射**

适用于粒子尺寸远小于波长的情况（x = 2πr/λ ≪ 1）：
$$p_{Rayleigh}(\cos\theta) = \frac{3}{16\pi}(1 + \cos^2\theta)$$

各向异性参数：g = 0（对称散射）

偏振形式（考虑偏振态）：
$$p_{Rayleigh}(\cos\theta, \phi) = \frac{3}{16\pi} \begin{pmatrix} 1 + \cos^2\theta & \sin^2\theta \cos(2\phi) \\ \sin^2\theta \cos(2\phi) & 1 + \cos^2\theta \end{pmatrix}$$

采样方法：
1. 生成 u = 2ξ₁ - 1
2. 如果 u² > (1 + u²)ξ₂，拒绝并重试
3. 否则 cosθ = u

**3. Henyey-Greenstein 相位函数**

最常用的参数化相位函数：
$$p_{HG}(\cos\theta) = \frac{1}{4\pi} \frac{1 - g^2}{(1 + g^2 - 2g\cos\theta)^{3/2}}$$

优点：
- 单参数 g ∈ [-1,1] 控制各向异性
- 解析形式简单
- 满足归一化条件
- 存在解析采样方法

缺点：
- 无法准确模拟强后向散射峰
- 对某些材质（如生物组织）过于简化

勒让德展开：
$$p_{HG}(\cos\theta) = \frac{1}{4\pi} \sum_{l=0}^{\infty} (2l+1) g^l P_l(\cos\theta)$$

解析采样（见后续重要性采样部分）

**4. Schlick 相位函数**

Henyey-Greenstein 的有理近似，计算更高效：
$$p_{Schlick}(\cos\theta) = \frac{1}{4\pi} \frac{1 - k^2}{(1 + k\cos\theta)^2}$$

参数转换：k = 1.55g - 0.55g³（保持相同的平均余弦）

误差分析：|p_{Schlick} - p_{HG}|/p_{HG} < 0.02 对大部分角度

计算效率提升：避免了HG中的幂运算和平方根

**5. 双 Henyey-Greenstein（Double HG）**

混合前向和后向散射，更好地拟合实测数据：
$$p_{DHG}(\cos\theta) = \alpha p_{HG}(\cos\theta; g_1) + (1-\alpha) p_{HG}(\cos\theta; g_2)$$

参数约束确保物理合理性：
- α ∈ [0,1]：混合权重
- g₁ > 0：前向散射成分
- g₂ < 0：后向散射成分
- 平均各向异性：g = αg₁ + (1-α)g₂

典型参数（生物组织）：
- 皮肤：α ≈ 0.8, g₁ ≈ 0.9, g₂ ≈ -0.3
- 肌肉：α ≈ 0.9, g₁ ≈ 0.95, g₂ ≈ -0.5

**6. Fournier-Forand 相位函数**

海洋光学中常用，考虑颗粒大小分布：
$$p_{FF}(\cos\theta) = \frac{1}{4\pi(1-\delta)^2\delta^v} \left[ v(1-\delta) - (1-\delta^v) \right] + \frac{1-\delta^v}{4\pi(1-\delta)^2\delta^v} \frac{3\cos^2\theta - 1}{(\sin\theta/2)^2}$$

其中：
- δ = 4sin²(θ/2)/(3(n-1)²)
- v = (3-μ)/2
- μ：颗粒大小分布的Junge指数
- n：相对折射率

特点：
- 准确描述海洋悬浮粒子
- 包含强前向散射峰
- 小角度近似为衍射峰

**7. Mie 散射理论**

对于球形粒子的精确解，基于麦克斯韦方程：

散射振幅函数：
$$S_1(\theta) = \sum_{n=1}^{\infty} \frac{2n+1}{n(n+1)} [a_n \pi_n(\cos\theta) + b_n \tau_n(\cos\theta)]$$
$$S_2(\theta) = \sum_{n=1}^{\infty} \frac{2n+1}{n(n+1)} [a_n \tau_n(\cos\theta) + b_n \pi_n(\cos\theta)]$$

其中 aₙ, bₙ 是Mie系数，涉及球贝塞尔函数：
$$a_n = \frac{m\psi_n(mx)\psi'_n(x) - \psi_n(x)\psi'_n(mx)}{m\psi_n(mx)\xi'_n(x) - \xi_n(x)\psi'_n(mx)}$$

相位函数：
$$p_{Mie}(\cos\theta) = \frac{4\pi}{k^2\sigma_s} \left[ |S_1(\theta)|^2 + |S_2(\theta)|^2 \right]$$

计算考虑：
- 级数截断：n_max ≈ x + 4x^(1/3) + 2
- 数值稳定性：使用向下递推计算贝塞尔函数
- 效率优化：预计算并缓存系数

**8. SGGX 相位函数**

基于微表面理论，适用于纤维和毛发：
$$p_{SGGX}(\cos\theta) = \frac{D(\mathbf{h})}{4|\boldsymbol{\omega}' \cdot \mathbf{h}|}$$

其中 D(h) 是SGGX分布，h 是半向量。

SGGX分布：
$$D(\mathbf{h}) = \frac{1}{\pi \alpha_x \alpha_y} \frac{1}{(\mathbf{h} \cdot \mathbf{S}^{-1} \mathbf{h})^2}$$

S 是3×3协方差矩阵，描述纤维方向分布。

**相位函数的重要性采样**

**Henyey-Greenstein 采样**：

给定均匀随机数 ξ₁, ξ₂ ∈ [0,1]：
1. 采样散射角：
$$\cos\theta = \begin{cases}
\frac{1}{2g}\left[1 + g^2 - \left(\frac{1-g^2}{1-g+2g\xi_1}\right)^2\right] & g \neq 0 \\
1 - 2\xi_1 & g = 0
\end{cases}$$

2. 采样方位角：φ = 2πξ₂

3. 构造散射方向（局部坐标系）：
   - 构建正交基 {t, b, ω'}
   - sinθ = √(1 - cos²θ)
   - ω = sinθ·cosφ·t + sinθ·sinφ·b + cosθ·ω'

PDF值：p_{HG}(cosθ)

**拒绝采样法（通用方法）**：

对于复杂相位函数：
1. 找到包络函数 M·q(θ) ≥ p(θ)
2. 从 q(θ) 采样候选值
3. 以概率 p(θ)/(M·q(θ)) 接受

效率：η = 1/M，选择好的q(θ)最大化效率

**表格法采样**：

1. 预计算CDF：F(θᵢ) = ∫₀^{θᵢ} p(θ)sinθ dθ
2. 二分查找：F⁻¹(ξ)
3. 线性插值获得精确值

存储优化：非均匀采样，在变化剧烈区域加密

**混合重要性采样**：

结合相位函数和其他因素（如光源方向）：
$$p_{mix}(\boldsymbol{\omega}) = \alpha p_{phase}(\boldsymbol{\omega}) + (1-\alpha) p_{light}(\boldsymbol{\omega})$$

使用多重重要性采样（MIS）组合不同策略：
$$w_i = \frac{p_i^2}{\sum_j p_j^2} \quad \text{(balance heuristic)}$$

**相位函数的预计算与优化**

1. **球谐系数缓存**：预计算前N项Legendre系数
2. **查找表**：离散化角度空间，存储函数值
3. **解析近似**：使用有理函数逼近复杂相位函数
4. **GPU优化**：向量化计算，利用纹理缓存

## 5.3 通过扩散的次表面散射

次表面散射（Subsurface Scattering, SSS）描述光进入半透明材质后在内部散射并从不同位置出射的现象。当散射远强于吸收时，可用扩散近似简化辐射传输方程。

### 5.3.1 扩散方程推导

**P₁近似（球谐函数展开）**

将辐射亮度展开为球谐函数的前两项：
$$L(\mathbf{x}, \boldsymbol{\omega}) \approx \frac{1}{4\pi}\phi(\mathbf{x}) + \frac{3}{4\pi}\mathbf{E}(\mathbf{x}) \cdot \boldsymbol{\omega}$$

其中：
- φ(**x**) = ∫_{S²} L(**x**, **ω**) d**ω**：辐射通量（fluence）
- **E**(**x**) = ∫_{S²} L(**x**, **ω**)**ω** d**ω**：辐射通量矢量

**从RTE到扩散方程的严格推导**

将P₁展开代入RTE：
$$(\boldsymbol{\omega} \cdot \nabla)\left[\frac{\phi}{4\pi} + \frac{3}{4\pi}\mathbf{E} \cdot \boldsymbol{\omega}\right] + \sigma_t\left[\frac{\phi}{4\pi} + \frac{3}{4\pi}\mathbf{E} \cdot \boldsymbol{\omega}\right] = \text{源项}$$

对所有方向积分（零阶矩）：
$$\nabla \cdot \mathbf{E} + \sigma_a \phi = Q$$

对**ω**加权积分（一阶矩）：
$$\frac{1}{3}\nabla\phi + \sigma_t' \mathbf{E} = \mathbf{0}$$

其中利用了 ∫**ω**⊗**ω** d**ω** = (4π/3)**I**。

**扩散近似的条件**

1. 散射占主导：σₛ ≫ σₐ（即 α ≈ 1）
2. 近似各向同性：|**E**| ≪ φ
3. 远离边界和光源（距离 > 几个平均自由程）

在这些条件下，辐射通量矢量遵循Fick定律：
$$\mathbf{E}(\mathbf{x}) = -D \nabla \phi(\mathbf{x})$$

扩散系数：
$$D = \frac{1}{3\sigma_t'} = \frac{1}{3[\sigma_a + \sigma_s(1-g)]}$$

其中 σₜ' = σₐ + σₛ' 是约化消光系数，σₛ' = σₛ(1-g) 是约化散射系数。

**扩散方程**

将Fick定律代入连续性方程：
$$\nabla \cdot (-D\nabla\phi) + \sigma_a\phi = Q$$

对于均匀介质：
$$\nabla^2 \phi(\mathbf{x}) - \frac{\sigma_a}{D} \phi(\mathbf{x}) = -\frac{Q(\mathbf{x})}{D}$$

或写作：
$$\nabla^2 \phi - \kappa^2 \phi = -\frac{Q}{D}$$

其中：
- κ = √(σₐ/D) = √(3σₐσₜ')：有效衰减系数
- 扩散长度：L_D = 1/κ

**边界条件**

在介质-真空界面，需要考虑菲涅尔反射。

**精确边界条件**（Marshak）：
$$\phi(\mathbf{x}_s) = \int_{\boldsymbol{\omega} \cdot \mathbf{n} < 0} L(\mathbf{x}_s, \boldsymbol{\omega}) d\boldsymbol{\omega} = 2\pi \int_0^{1} L(\mathbf{x}_s, -\mu) \mu d\mu$$

**外推边界条件**（近似）：
$$\phi(\mathbf{x}_s) + 2AD(\mathbf{n} \cdot \nabla)\phi(\mathbf{x}_s) = 0$$

其中 A 是与菲涅尔反射相关的参数：
$$A = \frac{1 + F_{dr}}{1 - F_{dr}}$$

**菲涅尔反射率的计算**

漫反射菲涅尔反射率：
$$F_{dr} = \int_0^{\pi/2} R(\theta) \sin(2\theta) d\theta$$

对于相对折射率 η = n₂/n₁：
- η < 1：F_dr ≈ -1.440η⁻² + 0.710η⁻¹ + 0.668 + 0.0636η
- η > 1：F_dr ≈ -0.4399 + 0.7099/η - 0.3319/η² + 0.0636/η³

**扩散近似的有效性范围**

定义无量纲参数：
- 光学厚度：τ = σₜ'L（L是特征长度）
- 反照率：α' = σₛ'/σₜ'

扩散近似误差：
$$\epsilon \approx \frac{1}{\sqrt{3\tau(1-\alpha')}}$$

要求 ε < 0.1，需要：
- τ > 10 且 α' > 0.9
- 或 τ(1-α') > 3.3

### 5.3.2 偶极子近似

**点光源的格林函数**

对于无限介质中的点光源 Q(**x**) = δ(**x** - **x₀**)，扩散方程的解为：
$$\phi(\mathbf{x}) = \frac{1}{4\pi D} \frac{e^{-\kappa r}}{r}$$

其中 r = |**x** - **x₀**|。

**偶极子源配置**

为满足边界条件，使用镜像源方法：
1. 真实源：位于 **xᵣ** = **x₀** + zₘ**n**（深度 zₘ = 1/σₜ'）
2. 镜像源：位于 **xᵥ** = **x₀** - (zₘ + 2zᵦ)**n**

其中 zᵦ = 2AD 是外推边界距离。

**BSSRDF的扩散近似**

双向次表面散射分布函数（BSSRDF）定义出射辐射亮度与入射辐照度的关系：
$$S(\mathbf{x}_i, \boldsymbol{\omega}_i; \mathbf{x}_o, \boldsymbol{\omega}_o) = \frac{dL_o(\mathbf{x}_o, \boldsymbol{\omega}_o)}{dE_i(\mathbf{x}_i, \boldsymbol{\omega}_i)}$$

在扩散近似下：
$$S_d(\mathbf{x}_i, \mathbf{x}_o) = \frac{1}{4\pi} F_t(\eta, \boldsymbol{\omega}_i) R_d(||\mathbf{x}_i - \mathbf{x}_o||) F_t(\eta, \boldsymbol{\omega}_o)$$

扩散剖面 Rd(r) 为：
$$R_d(r) = \frac{\alpha'}{4\pi} \left[ \frac{e^{-\sigma_{tr}d_r}}{d_r^2} + \frac{e^{-\sigma_{tr}d_v}}{d_v^2} \right]$$

其中：
- dᵣ = √(r² + zᵣ²)：到真实源的距离
- dᵥ = √(r² + zᵥ²)：到虚拟源的距离
- σₜᵣ = √(3σₐσₜ')：有效传输系数
- α' = σₛ'/σₜ'：约化反照率

### 5.3.3 多层材质模型

**层状介质的解析解**

对于N层平行介质，每层有不同的光学性质{σₐ,ᵢ, σₛ,ᵢ, gᵢ}。在每层内：
$$\frac{d^2\phi_i}{dz^2} - \kappa_i^2 \phi_i = 0$$

通解形式：
$$\phi_i(z) = A_i e^{-\kappa_i z} + B_i e^{\kappa_i z}$$

**界面连续性条件**

在层间界面 z = zᵢ：
1. 通量连续：φᵢ(zᵢ) = φᵢ₊₁(zᵢ)
2. 通量导数连续：Dᵢ(dφᵢ/dz)|ᵢ = Dᵢ₊₁(dφᵢ₊₁/dz)|ᵢ

构成线性方程组，可解得各层系数。

**多极展开方法**

对于更复杂的几何，使用多极展开：
$$\phi(\mathbf{x}) = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} a_{nm} \psi_{nm}(\mathbf{x})$$

其中 ψₙₘ 是满足边界条件的基函数。

**快速多极方法（FMM）**

对于大量散射源，直接求和的复杂度为 O(N²)。FMM通过层次分解将复杂度降至 O(N log N)：

1. 将空间分解为八叉树
2. 对远场使用多极展开
3. 对近场直接计算

**量化多层SSS**

定义有效穿透深度：
$$d_{eff} = \frac{1}{\sigma_{tr}} = \frac{1}{\sqrt{3\sigma_a\sigma_t'}}$$

多层材质的有效BSSRDF可通过加权平均近似：
$$S_{multi} \approx \sum_i w_i S_i$$

权重 wᵢ 基于各层的光学厚度和位置。

## 5.4 参与介质

参与介质（Participating Media）是指光在其中传播时会发生吸收、发射和散射的介质，如烟雾、云、水等。本节探讨如何在体积渲染框架下处理这些介质。

### 5.4.1 均匀与非均匀介质

**均匀介质**

在均匀介质中，光学性质不随空间变化：σₐ、σₛ、g 为常数。

透射率的解析形式：
$$T(s) = e^{-\sigma_t s}$$

内散射积分简化为：
$$L_s = \sigma_s \int_0^s e^{-\sigma_t t} \int_{S^2} p(\boldsymbol{\omega}' \to \boldsymbol{\omega}) L_{in}(\mathbf{x}_0 + t\boldsymbol{\omega}, \boldsymbol{\omega}') d\boldsymbol{\omega}' dt$$

**非均匀介质**

光学性质随空间变化：σₐ(**x**)、σₛ(**x**)、g(**x**)。

透射率的一般形式：
$$T(s_1, s_2) = \exp\left(-\int_{s_1}^{s_2} \sigma_t(\mathbf{x}(s)) ds\right)$$

光学厚度（optical depth）：
$$\tau(s_1, s_2) = \int_{s_1}^{s_2} \sigma_t(\mathbf{x}(s)) ds$$

则 T = e^(-τ)。

**介质的数学分类**

1. **分段常数**：空间划分为有限个均匀区域
2. **连续变化**：σₜ(**x**) 是连续函数
3. **随机介质**：σₜ(**x**) 是随机场

对于随机介质，使用统计描述：
$$\langle T \rangle = \exp\left(-\int \langle\sigma_t\rangle ds\right) \cdot C(\tau)$$

其中 C(τ) 是相关修正因子。

### 5.4.2 光线行进算法

**Ray Marching 基础**

将射线离散为 N 步，步长 Δs：
$$L_{out} = L_{in} \prod_{i=0}^{N-1} T_i + \sum_{i=0}^{N-1} L_{s,i} \prod_{j=i+1}^{N-1} T_j$$

其中：
- Tᵢ = exp(-σₜ,ᵢ Δs)：第 i 步的透射率
- Lₛ,ᵢ：第 i 步的散射贡献

**自适应步长策略**

1. **基于光学厚度**：保持每步 Δτ = σₜ Δs 恒定
   $$\Delta s_i = \frac{\Delta \tau_{target}}{\sigma_t(\mathbf{x}_i)}$$

2. **基于误差估计**：使用Richardson外推
   $$\epsilon = |L_{N} - L_{2N}|$$
   
   当 ε > εₘₐₓ 时细分步长。

3. **重要性引导**：在重要区域（如光源附近）使用更小步长

**Woodcock Tracking（Delta Tracking）**

虚拟粒子方法，将非均匀介质转化为均匀介质处理：

1. 定义主消光系数：σ̄ₜ ≥ max σₜ(**x**)
2. 虚拟消光系数：σₙ(**x**) = σ̄ₜ - σₜ(**x**)
3. 采样自由程：s ~ σ̄ₜ exp(-σ̄ₜs)
4. 在采样点以概率 σₜ(**x**)/σ̄ₜ 发生真实相互作用

优点：无偏，处理高频变化介质
缺点：σ̄ₜ ≫ σₜ 时效率低

**Ratio Tracking**

Woodcock tracking 的推广，使用控制变量：
$$\bar{\sigma}_t(\mathbf{x}) = \mu(\mathbf{x}) \sigma_t^{max}$$

其中 0 < μ(**x**) ≤ 1 是空间变化的控制函数。

**Spectral Tracking**

对于波长相关的介质，独立追踪每个波长会导致相关性丢失。联合采样策略：
$$P(s, \lambda) = \sigma_t(\mathbf{x}(s), \lambda) \exp\left(-\int_0^s \sigma_t(\mathbf{x}(t), \lambda) dt\right) P(\lambda)$$

### 5.4.3 重要性采样策略

**距离采样**

给定透射率 T(s) = exp(-∫σₜds)，采样碰撞距离：

1. **均匀介质**：逆变换采样
   $$s = -\frac{\ln(1-\xi)}{\sigma_t}$$

2. **分段常数**：逐段采样
   ```
   τ_target = -ln(1-ξ)
   τ_sum = 0
   for each segment i:
       τ_i = σ_t,i * length_i
       if τ_sum + τ_i > τ_target:
           s = (τ_target - τ_sum) / σ_t,i
           break
       τ_sum += τ_i
   ```

3. **一般非均匀**：数值求解
   $$\int_0^s \sigma_t(\mathbf{x}(t)) dt = -\ln(1-\xi)$$

**相位函数采样**

1. **Henyey-Greenstein**：见5.2.3节
2. **Rayleigh**：使用拒绝采样或分解方法
3. **Tabulated**：构建CDF并二分查找

**联合重要性采样**

考虑光源方向的重要性：
$$p(\boldsymbol{\omega}) = \alpha p_{phase}(\boldsymbol{\omega}) + (1-\alpha) p_{light}(\boldsymbol{\omega})$$

MIS（Multiple Importance Sampling）权重：
$$w_{phase} = \frac{p_{phase}^2}{p_{phase}^2 + p_{light}^2}$$

**Equiangular Sampling**

对于点光源或小光源，沿射线的重要性呈现峰值。等角采样：
$$t(\xi) = D \tan\left((1-\xi)\theta_a + \xi\theta_b\right)$$

其中 D 是到光源的最近距离，[θₐ, θᵦ] 是角度范围。

**Joint Importance Sampling**

结合距离采样和等角采样：
$$p(t) = \frac{1}{2} p_{distance}(t) + \frac{1}{2} p_{equiangular}(t)$$

使用MIS组合两种策略的贡献。

**体积光子映射**

存储介质中的光子：
```
struct VolumePhoton {
    vec3 position
    vec3 direction  
    vec3 power
    float radius    // 用于密度估计
}
```

辐射亮度估计：
$$L(\mathbf{x}, \boldsymbol{\omega}) = \frac{1}{4\pi r^3} \sum_{i \in N(\mathbf{x},r)} p(\boldsymbol{\omega}_i \to \boldsymbol{\omega}) \Phi_i$$

## 5.5 能量守恒与互易性

能量守恒和互易性是基于物理的渲染的两个基本原理。它们不仅确保渲染结果的物理正确性，还提供了验证算法实现的有力工具。

### 5.5.1 能量守恒的数学表述

**BRDF的能量守恒**

对于任意入射方向 **ωᵢ**，反射的总能量不能超过入射能量：
$$\int_{\Omega} f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o \leq 1$$

定义方向反射率（directional albedo）：
$$\rho_d(\boldsymbol{\omega}_i) = \int_{\Omega} f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o$$

能量守恒要求：ρd(**ωᵢ**) ≤ 1 对所有方向成立。

**半球反射率（Hemispherical Albedo）**

平均所有入射方向：
$$\rho_{hh} = \frac{1}{\pi} \int_{\Omega} \int_{\Omega} f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_i \cdot \mathbf{n}) (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_i d\boldsymbol{\omega}_o$$

对于朗伯表面：ρhh = ρ（反射率常数）。

**相位函数的能量守恒**

已在5.2.2节讨论：
$$\int_{S^2} p(\boldsymbol{\omega}' \to \boldsymbol{\omega}) d\boldsymbol{\omega} = 1$$

**辐射传输方程的能量守恒**

对RTE在所有方向积分：
$$\int_{S^2} (\boldsymbol{\omega} \cdot \nabla) L d\boldsymbol{\omega} + \sigma_a \phi - \sigma_a \phi_e = 0$$

其中利用了散射项的守恒性。这导出辐射通量的守恒方程：
$$\nabla \cdot \mathbf{F} + \sigma_a (\phi - \phi_e) = 0$$

**渲染方程的能量分析**

定义传输算子 𝒯：
$$(\mathcal{T}L)(\mathbf{x}, \boldsymbol{\omega}) = \int_{\mathcal{M}} \int_{S^2} K(\mathbf{x}, \boldsymbol{\omega} \leftarrow \mathbf{x}', \boldsymbol{\omega}') L(\mathbf{x}', \boldsymbol{\omega}') d\boldsymbol{\omega}' d\mathbf{x}'$$

能量守恒要求算子范数 ||𝒯|| < 1，保证Neumann级数收敛：
$$L = L_e + \mathcal{T}L_e + \mathcal{T}^2L_e + ...$$

### 5.5.2 Helmholtz互易原理

**基本陈述**

光路可逆：交换光源和探测器位置，测得的辐射亮度相同。

数学表述：
$$f_r(\boldsymbol{\omega}_i \to \boldsymbol{\omega}_o) = f_r(\boldsymbol{\omega}_o \to \boldsymbol{\omega}_i)$$

**互易性的推广**

1. **BSSRDF互易性**：
   $$S(\mathbf{x}_i, \boldsymbol{\omega}_i; \mathbf{x}_o, \boldsymbol{\omega}_o) = S(\mathbf{x}_o, -\boldsymbol{\omega}_o; \mathbf{x}_i, -\boldsymbol{\omega}_i)$$

2. **相位函数互易性**：
   $$p(\boldsymbol{\omega}_i \to \boldsymbol{\omega}_o) = p(-\boldsymbol{\omega}_o \to -\boldsymbol{\omega}_i)$$

3. **传输互易性**：
   设 G(**x**, **y**) 为从 **y** 到 **x** 的格林函数，则：
   $$G(\mathbf{x}, \mathbf{y}) = G(\mathbf{y}, \mathbf{x})$$

**互易性的物理起源**

基于时间反演对称性和洛伦兹互易定理。在经典电磁理论中：
$$\int_V \mathbf{J}_1 \cdot \mathbf{E}_2 dV = \int_V \mathbf{J}_2 \cdot \mathbf{E}_1 dV$$

**非互易情况**

1. **磁光效应**：存在外磁场时
2. **非线性光学**：强光场下
3. **荧光和磷光**：涉及能级跃迁

**互易性在渲染中的应用**

1. **双向路径追踪**：连接光路时利用互易性
2. **光子映射**：光子传播和聚集查询的对称性
3. **算法验证**：检查BRDF实现是否满足互易性

### 5.5.3 白炉测试与验证

**白炉测试（White Furnace Test）**

在完全封闭的均匀发光环境中，所有点的辐射亮度应相等。

设置：
- 封闭场景，所有表面 ρ = 1（完美漫反射）
- 均匀发射 Lₑ = 1
- 无吸收介质

理论结果：L(**x**, **ω**) = 1 处处成立。

**数值验证方法**

1. **蒙特卡洛收敛性**：
   $$\text{Variance}[L_N] \propto \frac{1}{N}$$
   
   其中 N 是样本数。

2. **能量平衡检查**：
   $$\left| \frac{\text{Power}_{out} - \text{Power}_{in}}{\text{Power}_{in}} \right| < \epsilon$$

3. **互易性测试**：
   交换光源和相机，比较结果：
   $$\left| \frac{L_{forward} - L_{backward}}{L_{forward}} \right| < \epsilon$$

**体积介质的白炉测试**

对于参与介质，修改条件：
- 单散射反照率 α = 1（纯散射）
- 各向同性相位函数

稳态解：φ(**x**) = const

**灰炉测试（Gray Furnace Test）**

更一般的情况，表面反射率 ρ < 1：
$$L = \frac{L_e}{1 - \rho}$$

这提供了多次反射正确实现的验证。

**基于物理约束的调试**

1. **正定性**：L ≥ 0, T ∈ [0,1]
2. **归一化**：∫p dω = 1, ∫BRDF cos θ dω ≤ 1  
3. **对称性**：检查旋转不变性等
4. **极限情况**：
   - σₜ → 0：退化为真空
   - σₛ → 0：纯吸收（Beer定律）
   - g → 0：各向同性散射

**误差度量**

相对误差：
$$E_{rel} = \frac{||L_{computed} - L_{reference}||}{||L_{reference}||}$$

RMSE（均方根误差）：
$$E_{RMS} = \sqrt{\frac{1}{N} \sum_i (L_i^{comp} - L_i^{ref})^2}$$

感知误差（如SSIM）考虑人眼特性。

## 本章小结

本章建立了基于物理的渲染的数学框架，核心是辐射传输方程（RTE）：

$$(\boldsymbol{\omega} \cdot \nabla) L + \sigma_t L = \sigma_a L_e + \sigma_s \int_{S^2} p(\boldsymbol{\omega}' \to \boldsymbol{\omega}) L(\boldsymbol{\omega}') d\boldsymbol{\omega}'$$

关键概念：
- **光学系数**：吸收 σₐ、散射 σₛ、消光 σₜ = σₐ + σₛ
- **相位函数**：描述散射方向分布，如Henyey-Greenstein
- **扩散近似**：高散射介质的简化，导出BSSRDF
- **参与介质算法**：Ray Marching、Woodcock Tracking等
- **基本原理**：能量守恒和Helmholtz互易性

这些概念构成了体积渲染的统一框架，为后续章节中的神经渲染方法提供物理基础。

## 练习题

### 基础题

**练习 5.1** 证明对于纯吸收介质（σₛ = 0），辐射传输方程的解满足比尔-朗伯定律。

<details>
<summary>提示</summary>
将 σₛ = 0 代入RTE，得到一阶常微分方程。
</details>

<details>
<summary>答案</summary>
设 σₛ = 0，RTE简化为：
$$(\boldsymbol{\omega} \cdot \nabla) L + \sigma_a L = \sigma_a L_e$$

沿射线参数化，∂L/∂s + σₐL = σₐLₑ。

齐次解：L_h = Ce^(-σₐs)

特解：设L_p = A，代入得A = Lₑ

通解：L(s) = Lₑ + (L(0) - Lₑ)e^(-σₐs)

当Lₑ = 0时，L(s) = L(0)e^(-σₐs)，即比尔-朗伯定律。
</details>

**练习 5.2** 推导Henyey-Greenstein相位函数的归一化常数。

<details>
<summary>提示</summary>
使用球坐标系，积分∫p(cosθ)sinθdθdφ = 1。
</details>

<details>
<summary>答案</summary>
HG相位函数：
$$p_{HG}(\cos\theta) = C \frac{1 - g^2}{(1 + g^2 - 2g\cos\theta)^{3/2}}$$

归一化条件：
$$\int_0^{2\pi} d\phi \int_0^{\pi} p_{HG}(\cos\theta) \sin\theta d\theta = 1$$

令μ = cosθ：
$$2\pi C (1-g^2) \int_{-1}^{1} \frac{d\mu}{(1 + g^2 - 2g\mu)^{3/2}} = 1$$

计算积分：
$$I = \frac{2}{|g|} \left[ \frac{1}{\sqrt{1+g^2-2g\mu}} \right]_{-1}^{1} = \frac{2}{|g|} \left( \frac{1}{|1-g|} - \frac{1}{1+|g|} \right) = \frac{4}{1-g^2}$$

因此 C = 1/(4π)。
</details>

**练习 5.3** 证明扩散系数 D = 1/(3σₜ') 的物理意义。

<details>
<summary>提示</summary>
考虑各向同性点源的稳态扩散，利用Fick定律。
</details>

<details>
<summary>答案</summary>
从RTE的P₁近似出发，辐射通量矢量：
$$\mathbf{E} = \int_{S^2} L(\boldsymbol{\omega})\boldsymbol{\omega} d\boldsymbol{\omega}$$

在扩散近似下，L近似各向同性加小的各向异性项：
$$L \approx \frac{\phi}{4\pi} + \frac{3}{4\pi} \mathbf{E} \cdot \boldsymbol{\omega}$$

将此代入RTE并对**ω**加权积分，得到：
$$\nabla \phi + 3\sigma_t' \mathbf{E} = 0$$

即 **E** = -∇φ/(3σₜ')，对比Fick定律 **E** = -D∇φ，得 D = 1/(3σₜ')。

物理意义：扩散系数与约化平均自由程 1/σₜ' 成正比，因子1/3来自三维空间的各向同性平均。
</details>

### 挑战题

**练习 5.4** 推导双层介质（表层+基底）的BSSRDF解析表达式。

<details>
<summary>提示</summary>
使用层间的连续性条件和边界条件，求解耦合的扩散方程。
</details>

<details>
<summary>答案</summary>
设两层光学性质为 (D₁, κ₁) 和 (D₂, κ₂)，厚度为 d。

每层的通解：
- 层1：φ₁(z) = A₁e^(-κ₁z) + B₁e^(κ₁z)
- 层2：φ₂(z) = A₂e^(-κ₂(z-d))

边界条件：
1. z = 0（表面）：φ₁ + 2AD₁∂φ₁/∂z = 0
2. z = d（界面）：φ₁ = φ₂, D₁∂φ₁/∂z = D₂∂φ₂/∂z
3. z → ∞：φ₂ → 0

求解得到传递矩阵，最终BSSRDF为两个指数项的组合，反映两层的贡献。
</details>

**练习 5.5** 设计一个自适应的ray marching算法，基于局部光学厚度调整步长。

<details>
<summary>提示</summary>
目标是保持每步的光学厚度Δτ近似恒定，同时限制最大和最小步长。
</details>

<details>
<summary>答案</summary>
算法框架：
1. 设定目标光学厚度 Δτ_target（如0.1）
2. 在当前位置 x 估计 σₜ(x)
3. 计算步长：Δs = min(max(Δτ_target/σₜ(x), Δs_min), Δs_max)
4. 预测下一步的σₜ，如果变化剧烈则细分
5. 使用误差估计（如Richardson外推）验证精度

关键优化：
- 使用σₜ的梯度信息预测变化
- 在边界附近强制小步长
- 缓存已计算的σₜ值
- 考虑视觉重要性（近处用小步长）
</details>

**练习 5.6** 分析Woodcock tracking在极端情况下的效率，并提出改进方案。

<details>
<summary>提示</summary>
考虑σₜ_max ≫ σₜ_average的情况，分析期望的拒绝次数。
</details>

<details>
<summary>答案</summary>
Woodcock tracking的效率分析：

平均接受概率：
$$p_{accept} = \frac{\langle \sigma_t \rangle}{\sigma_t^{max}}$$

期望步数：N = 1/p_accept

极端情况：
1. 稀疏峰值：σₜ在小区域很大，其他地方接近0
   - 效率极低：N ≈ σₜ^max/σₜ^avg ≫ 1
   
2. 改进方案：
   - 空间划分：将空间分成多个区域，每个使用局部σₜ^max
   - Ratio tracking：使用空间变化的控制函数μ(x)
   - 混合方法：在高变化区域用ray marching，均匀区域用Woodcock

实现要点：
- 预计算σₜ的空间分布统计
- 自适应选择采样策略
- 使用重要性信息（如到光源距离）调整控制函数
</details>

**练习 5.7** 证明在多次散射主导的情况下，各向异性散射介质可以用等效的各向同性散射介质近似。

<details>
<summary>提示</summary>
使用相似性理论（similarity theory），考虑多次散射后的角度分布。
</details>

<details>
<summary>答案</summary>
相似性理论的核心思想：经过多次散射后，光的方向分布趋于各向同性。

定义约化系数：
- σₛ' = σₛ(1-g)：约化散射系数
- σₜ' = σₐ + σₛ'：约化消光系数

证明步骤：
1. 考虑n次散射后的角度分布，使用中心极限定理
2. 显示均方散射角 <θ²> ∝ n(1-g)
3. 定义扩散长度 ℓ* = 1/σₜ'，在距离 ≫ ℓ* 时分布接近各向同性

等效性：
- 原介质：(σₐ, σₛ, g)
- 等效介质：(σₐ, σₛ', 0)

适用条件：
- 光学厚度 τ ≫ 1
- 观察距离 ≫ ℓ*
- 单散射反照率 α ≈ 1

这解释了为什么扩散近似在次表面散射中如此有效。
</details>

## 常见陷阱与错误

1. **相位函数归一化错误**
   - 错误：忘记4π因子或使用错误的积分域
   - 正确：∫_{S²} p dω = 1，注意是在整个球面积分

2. **能量不守恒的BRDF**
   - 错误：BRDF参数设置导致 ∫f_r cosθ dω > 1
   - 检查：实现白炉测试验证

3. **扩散近似的误用**
   - 错误：在薄层或低散射介质使用扩散近似
   - 准则：光学厚度 τ > 10 且 α > 0.9

4. **数值不稳定**
   - 错误：在 σₜ ≈ 0 时直接计算 exp(-σₜs)/σₜ
   - 解决：使用泰勒展开或特殊处理

5. **采样偏差**
   - 错误：重要性采样的PDF与实际采样不匹配
   - 验证：统计采样分布，与理论PDF比较

6. **边界条件处理**
   - 错误：忽略菲涅尔效应或使用错误的外推距离
   - 影响：次表面散射结果偏暗或偏亮

## 最佳实践检查清单

### 算法实现
- [ ] RTE求解器通过白炉测试
- [ ] BRDF满足能量守恒和互易性
- [ ] 相位函数正确归一化
- [ ] 数值方法在极限情况下稳定

### 性能优化
- [ ] 使用空间数据结构加速介质查询
- [ ] 实现多种采样策略并用MIS组合
- [ ] 缓存重复计算（如光学厚度积分）
- [ ] 考虑LOD：远处使用简化模型

### 验证测试
- [ ] 与解析解对比（如单次散射）
- [ ] 检查能量守恒（输入功率=输出功率）
- [ ] 验证互易性（交换光源和相机）
- [ ] 收敛性分析（误差 ∝ 1/√N）

### 参数设置
- [ ] 物理合理的光学系数（σₐ, σₛ基于测量数据）
- [ ] 相位函数的g参数符合材质特性
- [ ] 步长设置平衡精度和性能
- [ ] 采样数量基于场景复杂度调整
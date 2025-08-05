# 第6章：神经辐射场基础（NeRF）

神经辐射场（Neural Radiance Fields, NeRF）代表了计算机图形学中场景表示的范式转变。通过使用神经网络对连续体积密度和颜色场进行编码，NeRF统一了几何和外观建模，同时提供了前所未有的视图合成质量。本章探讨NeRF的数学基础，将其置于我们统一体积渲染框架内，并建立与后续章节中更高级神经表示方法的联系。

## 学习目标

完成本章后，您将能够：
1. 推导神经网络如何表示连续辐射场
2. 分析位置编码在克服频谱偏差中的作用
3. 将体积渲染方程表述为函数逼近问题
4. 理解NeRF优化的梯度流动态
5. 设计适当的正则化策略以提高泛化能力

## 章节大纲

### 6.1 引言与动机

传统的场景表示方法，如体素网格、点云或网格，本质上是离散的。它们在固定分辨率下采样3D空间，导致内存需求与分辨率立方成正比的问题。具体而言，对于分辨率为 $N^3$ 的体素网格，存储复杂度为 $\mathcal{O}(N^3)$，这在高分辨率下变得不可行。神经辐射场通过使用神经网络作为连续函数逼近器来解决这一限制。

#### 6.1.1 从离散到连续的范式转变

考虑传统离散表示的数学形式。体素网格可表示为：

$$\rho_{\text{voxel}}(\mathbf{x}) = \sum_{i,j,k} \rho_{ijk} \cdot \mathbb{1}_{V_{ijk}}(\mathbf{x})$$

其中 $\mathbb{1}_{V_{ijk}}$ 是体素 $(i,j,k)$ 的指示函数，$\rho_{ijk}$ 是存储的属性值。这种表示存在几个根本限制：

1. **分辨率-内存权衡**：提高分辨率导致内存需求立方增长
2. **固定采样**：无法在采样点之间进行连续查询
3. **混叠效应**：离散采样导致高频细节丢失
4. **拓扑约束**：离散网格强加了固定的连接性结构
5. **各向异性限制**：轴对齐的体素难以表示任意方向的细节

从信息论角度，离散表示的基本限制可以通过采样定理理解。对于带宽为 $B$ 的信号，Nyquist-Shannon定理要求采样率至少为 $2B$。在3D空间中，这转化为：

$$N_{\min} \geq (2B)^3$$

其中 $N_{\min}$ 是最小体素数。这表明高频细节需要指数级增长的存储。

神经辐射场通过将场景表示为神经网络参数化的连续函数来克服这些限制：

$$F_\Theta: (\mathbf{x}, \mathbf{d}) \mapsto (\sigma, \mathbf{c})$$

其中 $\mathbf{x} \in \mathbb{R}^3$ 是空间位置，$\mathbf{d} \in \mathbb{S}^2$ 是视角方向，$\sigma \in \mathbb{R}_+$ 是体积密度，$\mathbf{c} \in [0,1]^3$ 是RGB颜色，$\Theta \in \mathbb{R}^p$ 是神经网络的 $p$ 个参数。

**连续表示的信息论优势**：

考虑Kolmogorov复杂度 $K(S)$，即描述场景 $S$ 所需的最短程序长度。对于许多自然场景：

$$K(S_{\text{neural}}) \ll K(S_{\text{voxel}})$$

这是因为神经网络可以利用场景中的规律性和对称性。例如，一个包含重复纹理的表面在体素表示中需要 $\mathcal{O}(N^2)$ 存储，但神经网络可以通过学习周期函数以 $\mathcal{O}(1)$ 复杂度表示。

#### 6.1.2 与体积渲染方程的自然对应

这种表示与我们在第3章中建立的统一体积渲染方程自然契合：

$$L(\mathbf{r}) = \int_{t_n}^{t_f} T(t) \sigma(\mathbf{r}(t)) \mathbf{c}(\mathbf{r}(t), \mathbf{d}) dt$$

其中透射率定义为：

$$T(t) = \exp\left(-\int_{t_n}^t \sigma(\mathbf{r}(s)) ds\right)$$

射线参数化为 $\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$，其中 $\mathbf{o}$ 是射线原点，$\mathbf{d}$ 是单位方向向量。

这个方程的深层数学结构揭示了几个关键洞察：

1. **指数衰减的物理意义**：透射率 $T(t)$ 满足微分方程：
   $$\frac{dT}{dt} = -\sigma(\mathbf{r}(t)) T(t)$$
   
   这是Beer-Lambert定律在参与介质中的表现，描述了光在传播过程中的指数衰减。

2. **测度论解释**：体积渲染积分可以理解为关于测度 $d\mu = \sigma(\mathbf{r}(t)) dt$ 的Lebesgue-Stieltjes积分：
   $$L(\mathbf{r}) = \int_{[t_n, t_f]} T(t) \mathbf{c}(\mathbf{r}(t), \mathbf{d}) d\mu(t)$$
   
   这提供了处理不连续密度场（如表面）的严格框架。

3. **算子理论视角**：定义渲染算子 $\mathcal{R}$：
   $$\mathcal{R}[\sigma, \mathbf{c}](\mathbf{r}) = \int_{t_n}^{t_f} \exp\left(-\int_{t_n}^t \sigma(\mathbf{r}(s)) ds\right) \sigma(\mathbf{r}(t)) \mathbf{c}(\mathbf{r}(t), \mathbf{d}) dt$$
   
   这是一个非线性算子，其Fréchet导数为：
   $$D\mathcal{R}[\sigma, \mathbf{c}][\delta\sigma, \delta\mathbf{c}] = \int_{t_n}^{t_f} T(t) \left[\delta\sigma(t) \mathbf{c}(t) + \sigma(t) \delta\mathbf{c}(t) - \sigma(t) \mathbf{c}(t) \int_{t_n}^t \delta\sigma(s) ds\right] dt$$
   
   这对理解梯度传播至关重要。

#### 6.1.3 神经表示的理论优势

神经网络表示提供了几个理论优势：

1. **压缩表示**：根据Kolmogorov-Arnold表示定理，任何连续函数可以表示为有限个单变量连续函数的组合。神经网络通过其层次结构自然实现这种分解。

   更精确地，Kolmogorov-Arnold定理指出：对于任意连续函数 $f: [0,1]^n \rightarrow \mathbb{R}$，存在连续函数 $\phi_q$ 和 $\psi_{pq}$ 使得：
   $$f(x_1, ..., x_n) = \sum_{q=0}^{2n} \phi_q\left(\sum_{p=1}^n \psi_{pq}(x_p)\right)$$
   
   神经网络的层次结构提供了这种分解的可学习实现。

2. **隐式正则化**：神经网络的参数化引入了隐式的平滑性先验，这可以通过神经切线核（NTK）理论来理解：
   $$K_{\text{NTK}}(\mathbf{x}, \mathbf{x}') = \lim_{m \to \infty} \frac{1}{m} \sum_{i=1}^m \frac{\partial f(\mathbf{x}; \Theta)}{\partial \Theta_i} \frac{\partial f(\mathbf{x}'; \Theta)}{\partial \Theta_i}$$
   
   其中 $m$ 是网络宽度。这个核决定了函数空间的归纳偏置。
   
   **定理（NTK的RKHS特性）**：在适当的初始化下，NTK诱导的再生核希尔伯特空间（RKHS）包含所有满足某种平滑性条件的函数。具体地，对于深度 $L$ 的网络：
   $$\|f\|_{\mathcal{H}_{\text{NTK}}}^2 \approx \int |\hat{f}(\mathbf{k})|^2 (1 + |\mathbf{k}|^2)^L d\mathbf{k}$$
   
   这表明深度网络隐式地施加了高阶平滑性约束。

3. **自适应分辨率**：不同于固定分辨率的离散表示，神经网络可以自适应地分配容量到场景的不同部分。信息论分析表明，对于熵为 $H$ 的场景，神经网络可以达到接近最优的压缩率 $R \approx H$。
   
   **率失真理论分析**：设场景的概率分布为 $p(S)$，失真度量为 $d(S, \hat{S})$。率失真函数：
   $$R(D) = \inf_{p(\hat{S}|S): \mathbb{E}[d(S,\hat{S})] \leq D} I(S; \hat{S})$$
   
   神经网络通过变分推断近似这个最优编码：
   $$\min_\Theta \mathbb{E}_{S \sim p(S)} [d(S, F_\Theta)] + \beta \cdot \text{Complexity}(\Theta)$$

4. **函数空间的几何结构**：神经网络参数空间诱导了场景空间的黎曼度量：
   $$g_{ij}(\Theta) = \mathbb{E}_{\mathbf{x}, \mathbf{d}} \left[\frac{\partial F_\Theta(\mathbf{x}, \mathbf{d})}{\partial \Theta_i} \frac{\partial F_\Theta(\mathbf{x}, \mathbf{d})}{\partial \Theta_j}\right]$$
   
   这个Fisher信息度量提供了理解优化景观和泛化的几何框架。

#### 6.1.4 NeRF在统一框架中的位置

在我们的统一计算机图形学框架中，NeRF占据了一个特殊位置，连接了几个重要概念：

1. **与基于点的渲染的联系**（第3章）：NeRF可以视为无限密集点云的极限情况，其中每个空间点都有定义的属性。
   
   数学上，点云表示：
   $$\rho_{\text{points}}(\mathbf{x}) = \sum_{i=1}^N w_i \delta(\mathbf{x} - \mathbf{x}_i)$$
   
   NeRF提供了平滑逼近：
   $$\rho_{\text{NeRF}}(\mathbf{x}) = \lim_{\epsilon \to 0} \sum_{i=1}^N w_i K_\epsilon(\mathbf{x} - \mathbf{x}_i)$$
   
   其中 $K_\epsilon$ 是宽度为 $\epsilon$ 的平滑核。

2. **与基于图像的渲染的联系**（第4章）：NeRF隐式编码了光场，可以生成任意视角的图像。
   
   光场 $L(\mathbf{x}, \mathbf{d})$ 与NeRF的关系：
   $$L(\mathbf{x}, \mathbf{d}) = \int_0^\infty T(t) \sigma(\mathbf{x} + t\mathbf{d}) \mathbf{c}(\mathbf{x} + t\mathbf{d}, -\mathbf{d}) dt$$
   
   这建立了4D光场与5D辐射场之间的积分变换。

3. **与物理渲染的联系**（第5章）：通过适当的修改，NeRF可以扩展到建模参与介质和次表面散射。
   
   扩展到完整的辐射传输方程：
   $$\frac{d L(\mathbf{x}, \mathbf{d})}{dt} = -\sigma_t L + \sigma_s \int_{\mathbb{S}^2} p(\mathbf{d}' \rightarrow \mathbf{d}) L(\mathbf{x}, \mathbf{d}') d\mathbf{d}' + \sigma_a L_e$$
   
   NeRF的简化假设相当于忽略散射项（$\sigma_s = 0$）。

4. **向波动光学的桥梁**（第15-20章）：连续表示为后续引入相位信息和波动效应提供了自然基础。
   
   从标量场到复值场的推广：
   $$F_\Theta: (\mathbf{x}, \mathbf{d}, \lambda) \mapsto (A e^{i\phi}, \mathbf{c}) \in \mathbb{C} \times \mathbb{R}^3$$
   
   这允许建模干涉、衍射等波动现象。

5. **与逆向渲染的联系**（第11-14章）：NeRF的可微性使其成为逆向问题的理想表示。
   
   逆向渲染可表述为变分问题：
   $$\min_\Theta \sum_{i} \|\mathcal{R}[F_\Theta](\mathbf{r}_i) - I_i\|^2 + \mathcal{R}_{\text{reg}}(\Theta)$$
   
   其中 $\mathcal{R}$ 是渲染算子，$I_i$ 是观测图像。

### 6.2 连续体积表示

#### 6.2.1 神经隐式表示基础

隐式表示将几何编码为水平集函数或占用场。对于NeRF，我们扩展这一概念以包含外观信息。核心思想是学习一个函数 $F_\Theta$，它隐式地编码场景的几何和外观属性。

考虑几种经典的隐式表示及其与NeRF的关系：

1. **有符号距离函数（SDF）**：
   $$\phi: \mathbb{R}^3 \rightarrow \mathbb{R}, \quad \text{其中} \quad |\nabla \phi| = 1 \text{ a.e.}$$
   
   表面定义为零水平集：$\mathcal{S} = \{\mathbf{x} : \phi(\mathbf{x}) = 0\}$
   
   **Eikonal方程约束**：SDF满足偏微分方程：
   $$|\nabla \phi(\mathbf{x})| = 1, \quad \forall \mathbf{x} \in \mathbb{R}^3 \setminus \mathcal{S}_{\text{skeleton}}$$
   
   其中 $\mathcal{S}_{\text{skeleton}}$ 是中轴骷髅。这保证了到表面的最短距离属性。

2. **占用场（Occupancy Field）**：
   $$o: \mathbb{R}^3 \rightarrow [0,1], \quad o(\mathbf{x}) = \mathbb{P}[\mathbf{x} \in \mathcal{V}]$$
   
   其中 $\mathcal{V}$ 是占用体积。
   
   **与SDF的关系**：通过sigmoid函数转换：
   $$o(\mathbf{x}) = \text{sigmoid}(-\alpha \phi(\mathbf{x})) = \frac{1}{1 + \exp(\alpha \phi(\mathbf{x}))}$$
   
   其中 $\alpha$ 控制过渡的锐度。

3. **NeRF的体积密度场**：
   $$\sigma: \mathbb{R}^3 \rightarrow \mathbb{R}_+$$
   
   这可以理解为微分不透明度：$\sigma(\mathbf{x}) = -\frac{d\log T}{dt}|_{\mathbf{x}}$
   
   **物理解释**：$\sigma(\mathbf{x})$ 表示单位长度内光子被吸收的概率密度。在微小距离 $dt$ 内：
   $$\mathbb{P}[\text{absorption}] = \sigma(\mathbf{x}) dt + o(dt)$$

4. **广义隐式神经表示**：
   $$F_\Theta: \mathcal{X} \times \mathcal{D} \rightarrow \mathcal{Y}$$
   
   其中：
   - $\mathcal{X} \subseteq \mathbb{R}^3$：空间域
   - $\mathcal{D} \subseteq \mathbb{S}^2$：方向域
   - $\mathcal{Y} \subseteq \mathbb{R}^{d_y}$：属性空间

NeRF通过联合表示几何和外观扩展了这些概念：
$$F_\Theta: \mathbb{R}^3 \times \mathbb{S}^2 \rightarrow \mathbb{R}_+ \times [0,1]^3$$

**定理（表示能力）**：对于紧支撑集 $K \subset \mathbb{R}^3$ 上的任意连续密度场 $\sigma^*$ 和颜色场 $\mathbf{c}^*$，以及任意 $\epsilon > 0$，存在一个有限参数的神经网络 $F_\Theta$ 使得：

$$\sup_{(\mathbf{x},\mathbf{d}) \in K \times \mathbb{S}^2} \|F_\Theta(\mathbf{x},\mathbf{d}) - (\sigma^*(\mathbf{x}), \mathbf{c}^*(\mathbf{x},\mathbf{d}))\| < \epsilon$$

*证明要点*：利用Stone-Weierstrass定理的神经网络版本。由于 $K \times \mathbb{S}^2$ 是紧致的Hausdorff空间，且神经网络形成的函数类分离点、包含常函数且在点乘下封闭，故在连续函数空间中稠密。

**复杂度分析**：

这种表示的优势包括：
- **内存效率**：存储需求 $\mathcal{O}(p)$ 与场景复杂度无关，仅取决于网络参数数量 $p$
- **连续性**：可在任意分辨率下查询，支持超分辨率渲染
- **可微性**：几乎处处可微，支持基于梯度的优化和逆向渲染
- **适应性**：网络容量自动分配到复杂区域

#### 6.2.2 从体素到连续函数

考虑离散体素表示与连续神经表示之间的数学关系。离散体素网格可以视为：

$$\sigma_{\text{voxel}}(\mathbf{x}) = \sum_{i,j,k} \sigma_{ijk} \cdot \mathbb{1}_{V_{ijk}}(\mathbf{x})$$

其中 $\mathbb{1}_{V_{ijk}}$ 是体素 $(i,j,k)$ 的指示函数。这种表示本质上是在函数空间中使用了阶跃基函数。

**函数空间理论视角**：

设 $\mathcal{F}$ 为紧支撑上的连续函数空间 $C_c(\mathbb{R}^3)$。不同的表示方法对应于 $\mathcal{F}$ 中不同的子空间：

1. **体素空间**：$\mathcal{F}_{\text{voxel}} = \text{span}\{\mathbb{1}_{V_{ijk}}\}$
2. **样条空间**：$\mathcal{F}_{\text{spline}} = \{f : f \in C^k, f|_{V_{ijk}} \in \mathcal{P}^k\}$
3. **神经空间**：$\mathcal{F}_{\text{neural}} = \{F_\Theta : \Theta \in \mathbb{R}^p\}$

我们可以将这种表示推广到更一般的基函数展开：

$$\sigma(\mathbf{x}) = \sum_{n=1}^N c_n \psi_n(\mathbf{x})$$

不同的基函数选择导致不同的表示：

1. **体素基**：$\psi_n = \mathbb{1}_{V_n}$ - 零阶不连续
   - 逼近误差：$\mathcal{O}(h)$，$h$ 为体素尺寸
   - 存储复杂度：$\mathcal{O}(N^3)$

2. **三线性插值**：$\psi_n = \prod_{i=1}^3 (1-|x_i - x_i^n|/h)_+$ - 一阶连续
   - 逼近误差：$\mathcal{O}(h^2)$
   - 存储复杂度：$\mathcal{O}(N^3)$

3. **B样条基**：$\psi_n = B^k((\mathbf{x} - \mathbf{x}_n)/h)$ - $k$阶连续
   - 逼近误差：$\mathcal{O}(h^{k+1})$
   - 紧支撑：$\text{supp}(B^k) = [-k/2, k/2]^3$

4. **径向基函数**：$\psi_n = \exp(-\|\mathbf{x} - \mathbf{x}_n\|^2/\sigma^2)$ - 无限阶平滑
   - 谱衰减：$\hat{\psi}(\mathbf{k}) \propto \exp(-\sigma^2\|\mathbf{k}\|^2)$
   - 正定性：$\sum_{i,j} c_i c_j \psi(\mathbf{x}_i - \mathbf{x}_j) \geq 0$

神经网络表示可以理解为使用自适应、可学习基函数的推广：

$$\sigma_{\text{neural}}(\mathbf{x}) = g\left(\sum_{l=1}^L W_l^{(L)} \phi_{l-1}\left(\cdots \phi_1(W_1^{(1)}\mathbf{x} + \mathbf{b}_1^{(1)})\cdots\right)\right)$$

其中 $g$ 是输出激活函数（如softplus），$\phi_l$ 是隐藏层激活函数。

**定理（基函数等价性）**：具有ReLU激活的深度神经网络定义了分片线性基函数的自适应划分，其中：
- 划分区域数量：$\mathcal{O}((\text{neurons})^{\text{depth}})$
- 每个区域内：函数是线性的
- 边界复杂度：由网络架构决定

*证明要点*：ReLU神经网络 $f(\mathbf{x}) = \max(0, W\mathbf{x} + b)$ 在超平面 $\{\mathbf{x}: W_i\mathbf{x} + b_i = 0\}$ 上不可微。这些超平面将 $\mathbb{R}^3$ 划分成多个凸多面体区域，每个区域内网络是仿射函数。

**压缩效率分析**：

这种表示通过神经网络的组合性质实现了更高效的场景编码，其压缩率可以通过率失真理论分析：

$$R(D) = \inf_{F_\Theta} \{H(\Theta) : \mathbb{E}[\|\sigma - F_\Theta(\cdot)\|^2] \leq D\}$$

其中 $H(\Theta)$ 是参数的熵，$D$ 是允许的失真。

**定理（神经压缩效率）**：对于具有Lipschitz常数 $L$ 的场景函数，神经网络表示的压缩率满足：
$$R(D) \leq C \cdot \log\left(\frac{L^3 \text{Vol}(K)}{D}\right)$$

其中 $\text{Vol}(K)$ 是场景体积，$C$ 是常数。这表明神经网络可以达到近乎最优的压缩率。

#### 6.2.3 MLP作为辐射场编码器

标准的NeRF架构使用多层感知器（MLP）作为函数逼近器。让我们详细分析其架构设计和数学性质。

**基本架构**：

$$\mathbf{h}_0 = \gamma(\mathbf{x})$$
$$\mathbf{h}_{l+1} = \phi(W_l \mathbf{h}_l + \mathbf{b}_l), \quad l = 0, ..., L-1$$

对于密度预测：
$$\mathbf{f}_\sigma = \mathbf{h}_L$$
$$\sigma = \text{softplus}(w_\sigma^T \mathbf{f}_\sigma + b_\sigma) = \log(1 + \exp(w_\sigma^T \mathbf{f}_\sigma + b_\sigma))$$

对于颜色预测：
$$\mathbf{f}_c = [\mathbf{h}_L, \gamma(\mathbf{d})]$$
$$\mathbf{c} = \text{sigmoid}(W_c \mathbf{f}_c + \mathbf{b}_c) = \frac{1}{1 + \exp(-W_c \mathbf{f}_c - \mathbf{b}_c)}$$

其中 $\gamma$ 是位置编码函数，$\phi$ 是ReLU激活函数：$\phi(x) = \max(0, x)$。

**架构的关键设计原则**：

1. **视角无关的几何**：密度 $\sigma$ 仅依赖于位置 $\mathbf{x}$，确保几何的视角一致性
   
   数学上，这意味着：
   $$\frac{\partial \sigma}{\partial \mathbf{d}} = 0$$
   
   这保证了从不同视角看到的几何是一致的。
   
2. **视角相关的外观**：颜色 $\mathbf{c}$ 依赖于位置和方向 $(\mathbf{x}, \mathbf{d})$，允许建模镜面反射等效果
   
   这可以建模复杂的BRDF：
   $$\mathbf{c}(\mathbf{x}, \mathbf{d}) \approx \int_{\Omega} f_r(\mathbf{x}, \mathbf{d}_i \rightarrow \mathbf{d}) L_i(\mathbf{x}, \mathbf{d}_i) (\mathbf{n} \cdot \mathbf{d}_i) d\mathbf{d}_i$$

3. **跳跃连接**：在第 $l^*$ 层重新注入位置编码
   $$\mathbf{h}_{l^*+1} = \phi(W_{l^*} [\mathbf{h}_{l^*}, \gamma(\mathbf{x})] + \mathbf{b}_{l^*})$$
   
   这改善了梯度流并保留了高频信息。
   
   **理论分析**：跳跃连接可以视为残差学习：
   $$F(\mathbf{x}) = F_{\text{low}}(\mathbf{x}) + F_{\text{high}}(\mathbf{x})$$
   
   其中 $F_{\text{low}}$ 由前 $l^*$ 层学习，$F_{\text{high}}$ 由后续层学习。

**网络容量分析**：

对于深度 $L$、宽度 $w$ 的全连接网络，其表达能力可以通过以下度量：

1. **参数数量**：$p = \mathcal{O}(Lw^2)$
   
   具体地：
   $$p = \sum_{l=0}^{L-1} (w_l \cdot w_{l+1} + w_{l+1}) \approx Lw^2 + Lw$$

2. **线性区域数量**：$N_{\text{regions}} = \mathcal{O}((w/L)^L)$
   
   **定理（Montufar et al., 2014）**：具有 $L$ 层、每层 $w$ 个神经元的ReLU网络最多可以表示：
   $$N_{\text{regions}} \leq \left(\prod_{l=1}^{L-1} \lfloor w/d \rfloor^d\right) \cdot \sum_{j=0}^d \binom{w}{j}$$
   
   其中 $d$ 是输入维度。

3. **Lipschitz常数**：$\text{Lip}(F_\Theta) \leq \prod_{l=1}^L \|W_l\|_2$
   
   这决定了函数的平滑性：
   $$|F_\Theta(\mathbf{x}_1) - F_\Theta(\mathbf{x}_2)| \leq \text{Lip}(F_\Theta) \cdot \|\mathbf{x}_1 - \mathbf{x}_2\|$$

4. **VC维**：$\text{VCdim} = \mathcal{O}(pL)$
   
   这决定了泛化能力。

**激活函数选择的理论依据**：

- **Softplus用于密度**：确保 $\sigma \geq 0$，且在零点处平滑
  $$\text{softplus}'(x) = \text{sigmoid}(x) \in (0, 1)$$
  $$\text{softplus}''(x) = \text{sigmoid}(x)(1 - \text{sigmoid}(x)) \geq 0$$
  
  这保证了密度函数的凸性。
  
- **Sigmoid用于颜色**：确保 $\mathbf{c} \in [0,1]^3$，物理上对应归一化的RGB值
  
  Sigmoid的对称性：$\text{sigmoid}(-x) = 1 - \text{sigmoid}(x)$
  
- **ReLU用于隐藏层**：计算效率高，且保持了分片线性性质
  
  ReLU的正同性：$\text{ReLU}(\alpha x) = \alpha \text{ReLU}(x), \forall \alpha \geq 0$

**初始化策略**：

适当的初始化对于避免梯度消失/爆炸至关重要：

1. **Xavier/He初始化**：
   $$W_{ij} \sim \mathcal{N}(0, \sqrt{2/n_{\text{in}}})$$
   
   这保持了前向传播中激活值的方差稳定。

2. **SIREN初始化**（用于周期激活函数）：
   $$W_{ij} \sim \mathcal{U}(-\sqrt{6/n_{\text{in}}}, \sqrt{6/n_{\text{in}}})$$

#### 6.2.4 密度与颜色的联合建模

密度和颜色的耦合通过共享的特征表示实现：

$$\mathbf{f} = F_{\text{backbone}}(\gamma(\mathbf{x}))$$
$$\sigma = F_{\text{density}}(\mathbf{f})$$
$$\mathbf{c} = F_{\text{color}}(\mathbf{f}, \gamma(\mathbf{d}))$$

这种分解确保了物理合理性：
- 几何（密度）与视角无关
- 外观（颜色）可以建模视角相关效果（如镜面反射）

**信息论解释**：

从信息论角度，这种分解可以视为学习一个充分统计量：
$$\mathbf{f} = T(\mathbf{x})$$

其中 $T$ 满足：
$$I(\mathbf{x}; \sigma, \mathbf{c} | \mathbf{d}) = I(T(\mathbf{x}); \sigma, \mathbf{c} | \mathbf{d})$$

这意味着 $\mathbf{f}$ 包含了从位置预测密度和颜色所需的所有信息。

**分离表示的优势**：

1. **参数共享**：减少参数数量，提高效率
   $$p_{\text{shared}} < p_{\text{separate}} = p_{\sigma} + p_{\mathbf{c}}$$

2. **一致性约束**：密度和颜色共享底层特征，保证了几何和外观的一致性

3. **梯度传播**：来自颜色预测的梯度也会影响底层特征，间接改善密度预测

**理论分析：解耦学习**：

设联合分布为 $p(\sigma, \mathbf{c} | \mathbf{x}, \mathbf{d})$，我们可以分解为：
$$p(\sigma, \mathbf{c} | \mathbf{x}, \mathbf{d}) = p(\sigma | \mathbf{x}) \cdot p(\mathbf{c} | \mathbf{x}, \mathbf{d}, \sigma)$$

由于密度与视角无关：
$$p(\sigma | \mathbf{x}, \mathbf{d}) = p(\sigma | \mathbf{x})$$

这种因式分解指导了网络架构设计。

**物理约束的编码**：

1. **能量守恒**：通过限制颜色输出范围
   $$\|\mathbf{c}\|_\infty \leq 1$$

2. **空间连续性**：通过网络的Lipschitz性质
   $$\|\mathbf{c}(\mathbf{x}_1) - \mathbf{c}(\mathbf{x}_2)\| \leq L_c \|\mathbf{x}_1 - \mathbf{x}_2\|$$

3. **视角一致性**：对于漫反射表面
   $$\mathbf{c}(\mathbf{x}, \mathbf{d}_1) \approx \mathbf{c}(\mathbf{x}, \mathbf{d}_2)$$

### 6.3 位置编码与频谱偏差

#### 6.3.1 神经网络的频谱偏差问题

标准神经网络表现出强烈的频谱偏差，倾向于学习低频函数。这种现象可以从多个理论角度理解。

**频谱偏差的傅里叶分析**：

考虑目标函数的傅里叶展开：
$$f(\mathbf{x}) = \int_{\mathbb{R}^3} \hat{f}(\mathbf{k}) e^{i\mathbf{k} \cdot \mathbf{x}} d\mathbf{k}$$

标准神经网络在训练初期主要捕获低频成分：
$$f_{\text{neural}}(\mathbf{x}; t) \approx \int_{|\mathbf{k}| < k_c(t)} \hat{f}(\mathbf{k}) e^{i\mathbf{k} \cdot \mathbf{x}} d\mathbf{k}$$

其中截止频率 $k_c(t)$ 随训练时间缓慢增长。

**神经切线核（NTK）视角**：

在无限宽度极限下，神经网络的训练动态由NTK控制：
$$\frac{\partial f(\mathbf{x}; \Theta_t)}{\partial t} = -\eta \int K(\mathbf{x}, \mathbf{x}') \nabla_f \mathcal{L}(\mathbf{x}') d\mathbf{x}'$$

其中核函数为：
$$K(\mathbf{x}, \mathbf{x}') = \lim_{w \to \infty} \frac{1}{w} \sum_{i=1}^w \frac{\partial f(\mathbf{x}; \Theta)}{\partial \Theta_i} \frac{\partial f(\mathbf{x}'; \Theta)}{\partial \Theta_i}$$

对于ReLU网络，NTK具有特定的频谱特性：
$$\hat{K}(\mathbf{k}, \mathbf{k}') = \mathcal{F}[K](\mathbf{k}, \mathbf{k}') \propto \exp(-c|\mathbf{k}|^\alpha)$$

其中 $\alpha > 0$，表明高频成分被指数抑制。

**学习动态的频率依赖性**：

定义频率 $\mathbf{k}$ 的学习率为：
$$\tau_\mathbf{k}^{-1} = \lambda_\mathbf{k} \cdot \eta$$

其中 $\lambda_\mathbf{k}$ 是NTK在频率 $\mathbf{k}$ 处的特征值。对于标准初始化：
$$\lambda_\mathbf{k} \propto |\mathbf{k}|^{-\beta}, \quad \beta > 0$$

这导致高频成分的学习时间尺度：
$$\tau_\mathbf{k} \propto |\mathbf{k}|^\beta$$

因此，学习频率 $\mathbf{k}$ 所需的迭代次数随 $|\mathbf{k}|$ 多项式增长。

#### 6.3.2 傅里叶特征映射

为了克服频谱偏差，NeRF使用位置编码将输入映射到高维空间。这种技术基于随机傅里叶特征的理论，但使用确定性的频率选择。

**标准位置编码**：

$$\gamma(\mathbf{x}) = \left[\mathbf{x}, \sin(2^0\pi\mathbf{x}), \cos(2^0\pi\mathbf{x}), ..., \sin(2^{L-1}\pi\mathbf{x}), \cos(2^{L-1}\pi\mathbf{x})\right]^T$$

展开后，对于 $\mathbf{x} = (x, y, z)$：
$$\gamma(\mathbf{x}) = [x, y, z, \sin(\pi x), \cos(\pi x), \sin(\pi y), \cos(\pi y), ..., \sin(2^{L-1}\pi z), \cos(2^{L-1}\pi z)]^T$$

输出维度：$d_{\gamma} = 3 + 6L$（原始坐标 + 每个频率级别的正弦余弦对）。

**一般化的傅里叶特征**：

更一般的形式使用任意频率矩阵：
$$\gamma(\mathbf{x}) = \left[\sin(\mathbf{B}\mathbf{x}), \cos(\mathbf{B}\mathbf{x})\right]^T$$

其中 $\mathbf{B} \in \mathbb{R}^{m \times 3}$ 是频率矩阵。不同的 $\mathbf{B}$ 选择导致不同的编码：

1. **NeRF标准编码**：$\mathbf{B} = \text{diag}(2^0\pi, 2^1\pi, ..., 2^{L-1}\pi) \otimes \mathbf{I}_3$
2. **随机傅里叶特征**：$\mathbf{B}_{ij} \sim \mathcal{N}(0, \sigma^2)$
3. **学习的频率**：$\mathbf{B}$ 作为可训练参数

**与核方法的联系**：

位置编码等价于在特征空间中使用线性模型：
$$f(\mathbf{x}) = \mathbf{w}^T \gamma(\mathbf{x})$$

这对应于原始空间中的核函数：
$$K(\mathbf{x}, \mathbf{x}') = \gamma(\mathbf{x})^T \gamma(\mathbf{x}')$$

对于傅里叶特征：
$$K(\mathbf{x}, \mathbf{x}') = \sum_{j=1}^m \cos(\mathbf{b}_j^T(\mathbf{x} - \mathbf{x}'))$$

这是平稳核，其频谱由 $\mathbf{B}$ 的行向量 $\mathbf{b}_j$ 决定。

#### 6.3.3 位置编码的数学分析

位置编码从根本上改变了神经网络的逼近性质。让我们严格分析其数学效果。

**频谱扩展定理**：

**定理 6.1**：设 $f: \mathbb{R}^3 \rightarrow \mathbb{R}$ 是带限函数，其傅里叶变换满足 $\text{supp}(\hat{f}) \subseteq \{\mathbf{k}: |\mathbf{k}|_\infty \leq K\}$。使用位置编码 $\gamma$ 且最大频率 $2^{L-1}\pi \geq K$ 的单层线性网络可以精确表示 $f$。

*证明概要*：由Nyquist-Shannon定理，$f$ 可以表示为：
$$f(\mathbf{x}) = \sum_{\mathbf{n} \in \mathbb{Z}^3} f(\mathbf{n}/2K) \text{sinc}(2K\mathbf{x} - \mathbf{n})$$

使用傅里叶级数展开和位置编码的完备性，存在权重 $\mathbf{w}$ 使得：
$$f(\mathbf{x}) = \mathbf{w}^T \gamma(\mathbf{x})$$

**组合核分析**：

位置编码诱导的核函数：
$$K_\gamma(\mathbf{x}, \mathbf{x}') = \gamma(\mathbf{x})^T \gamma(\mathbf{x}')$$

对于NeRF的标准编码：
$$K_\gamma(\mathbf{x}, \mathbf{x}') = \mathbf{x}^T\mathbf{x}' + \sum_{l=0}^{L-1} \sum_{i=1}^3 \cos(2^l\pi(x_i - x_i'))$$

这个核的频谱具有离散支撑：
$$\hat{K}_\gamma(\mathbf{k}) = \delta(\mathbf{k}) + \sum_{l=0}^{L-1} \sum_{\epsilon \in \{-1,1\}^3} \delta(\mathbf{k} - \epsilon 2^l\pi)$$

**收敛速率分析**：

**定理 6.2**：对于Lipschitz连续函数 $f$ 和使用位置编码的神经网络 $f_\Theta$，逼近误差满足：
$$\|f - f_\Theta\|_{L^2} \leq C \cdot 2^{-L} \cdot \text{Lip}(f) + \mathcal{O}(n^{-1/2})$$

其中 $n$ 是网络参数数量，$L$ 是编码的频率级别数。

这表明：
1. 第一项：由带限近似引起的误差，随 $L$ 指数衰减
2. 第二项：由有限网络容量引起的误差，随参数数量多项式衰减

#### 6.3.4 其他编码方案比较

除了傅里叶编码，还有其他方案：

1. **学习的编码**：
   $$\gamma(\mathbf{x}) = \text{MLP}_\text{small}(\mathbf{x})$$

2. **球谐函数编码**（用于方向）：
   $$\gamma(\mathbf{d}) = [Y_l^m(\theta, \phi)]_{l=0}^{L_\text{max}}$$

3. **哈希编码**（用于加速）：
   $$\gamma(\mathbf{x}) = \bigoplus_{l=1}^L \text{hash}_l(\lfloor 2^l \mathbf{x} \rfloor)$$

每种编码在表达能力、计算效率和内存需求之间有不同权衡。

### 6.4 作为函数逼近的体积渲染

#### 6.4.1 离散积分的连续化

体积渲染方程需要评估积分：

$$C(\mathbf{r}) = \int_{t_n}^{t_f} T(t) \sigma(\mathbf{r}(t)) \mathbf{c}(\mathbf{r}(t), \mathbf{d}) dt$$

实际中，我们使用数值积分：

$$\hat{C}(\mathbf{r}) = \sum_{i=1}^N T_i (1 - \exp(-\sigma_i \delta_i)) \mathbf{c}_i$$

其中：
- $T_i = \exp\left(-\sum_{j=1}^{i-1} \sigma_j \delta_j\right)$ 是离散透射率
- $\delta_i = t_{i+1} - t_i$ 是采样间隔

这种离散化引入了两种误差：
1. **积分误差**：$\mathcal{O}(\max_i \delta_i)$
2. **函数逼近误差**：神经网络逼近真实场景的误差

#### 6.4.2 神经网络的万能逼近性质

**定理（万能逼近）**：对于任意连续函数 $f: K \rightarrow \mathbb{R}$（$K$ 是紧集）和 $\epsilon > 0$，存在一个神经网络 $F_\Theta$ 使得：

$$\sup_{\mathbf{x} \in K} |f(\mathbf{x}) - F_\Theta(\mathbf{x})| < \epsilon$$

对于NeRF，这意味着存在网络参数 $\Theta^*$ 使得：

$$\sup_{\mathbf{x}, \mathbf{d}} |(\sigma^*(\mathbf{x}), \mathbf{c}^*(\mathbf{x}, \mathbf{d})) - F_{\Theta^*}(\mathbf{x}, \mathbf{d})| < \epsilon$$

然而，找到这样的 $\Theta^*$ 是非凸优化问题。

#### 6.4.3 采样策略与积分精度

采样点的选择显著影响渲染质量。考虑误差界：

$$|C(\mathbf{r}) - \hat{C}(\mathbf{r})| \leq \frac{(t_f - t_n)^2}{8N} \max_{t \in [t_n, t_f]} \left|\frac{d^2}{dt^2}[T(t)\sigma(t)\mathbf{c}(t)]\right|$$

这表明：
- 均匀采样的误差为 $\mathcal{O}(1/N)$
- 在高频区域需要更密集的采样

#### 6.4.4 层次采样与重要性采样

NeRF使用两阶段采样策略：

**粗采样**：
$$t_i^c \sim \text{Uniform}[t_n, t_f], \quad i = 1, ..., N_c$$

**细采样**：基于粗采样的权重分布
$$w_i = T_i(1 - \exp(-\sigma_i \delta_i))$$
$$\hat{w}_i = w_i / \sum_j w_j$$
$$t^f \sim \sum_{i=1}^{N_c} \hat{w}_i \cdot \text{Uniform}[t_i^c, t_{i+1}^c]$$

这种重要性采样减少了方差：

$$\text{Var}[\hat{C}_{\text{hierarchical}}] \leq \text{Var}[\hat{C}_{\text{uniform}}]$$

### 6.5 基于梯度的优化

#### 6.5.1 损失函数设计

NeRF的基本损失函数是渲染图像与真实图像之间的光度误差：

$$\mathcal{L}_{\text{photo}} = \sum_{\mathbf{r} \in \mathcal{R}} \|\hat{C}(\mathbf{r}) - C_{\text{gt}}(\mathbf{r})\|_2^2$$

其中 $\mathcal{R}$ 是训练射线集合。这可以扩展为：

$$\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{photo}} + \lambda_{\text{reg}} \mathcal{L}_{\text{reg}}$$

常见的正则化项包括：
- **密度稀疏性**：$\mathcal{L}_{\text{sparse}} = \sum_{\mathbf{x}} \sigma(\mathbf{x})$
- **TV正则化**：$\mathcal{L}_{\text{TV}} = \sum_{\mathbf{x}} \|\nabla \sigma(\mathbf{x})\|_1$

#### 6.5.2 梯度流分析

考虑单个采样点的梯度：

$$\frac{\partial \mathcal{L}}{\partial \sigma_i} = 2(\hat{C} - C_{\text{gt}})^T \frac{\partial \hat{C}}{\partial \sigma_i}$$

其中：
$$\frac{\partial \hat{C}}{\partial \sigma_i} = T_i \delta_i \mathbf{c}_i - \sum_{j>i} \frac{\partial \hat{C}}{\partial T_j} \delta_i$$

这显示了梯度通过透射率的传播，创建了长程依赖关系。

**梯度消失问题**：当 $T_i \approx 0$（被遮挡的点），梯度 $\frac{\partial \mathcal{L}}{\partial \sigma_i} \approx 0$。这可能导致遮挡区域的学习困难。

#### 6.5.3 优化器选择与超参数

标准选择是Adam优化器，具有以下考虑：

$$\mathbf{m}_t = \beta_1 \mathbf{m}_{t-1} + (1-\beta_1) \mathbf{g}_t$$
$$\mathbf{v}_t = \beta_2 \mathbf{v}_{t-1} + (1-\beta_2) \mathbf{g}_t^2$$
$$\Theta_{t+1} = \Theta_t - \alpha \frac{\mathbf{m}_t}{\sqrt{\mathbf{v}_t} + \epsilon}$$

典型超参数：
- 学习率：$\alpha = 5 \times 10^{-4}$，带指数衰减
- $\beta_1 = 0.9$，$\beta_2 = 0.999$
- 批大小：1024-4096条射线

#### 6.5.4 收敛性分析

在理想条件下，NeRF优化可以视为随机梯度下降。收敛速率取决于：

1. **Lipschitz常数**：$L = \sup_{\Theta} \|\nabla^2 \mathcal{L}(\Theta)\|$
2. **方差界**：$\mathbb{E}[\|\nabla \mathcal{L}_i - \nabla \mathcal{L}\|^2] \leq \sigma^2$

在凸情况下，收敛速率为：
$$\mathbb{E}[\mathcal{L}(\Theta_T)] - \mathcal{L}^* \leq \mathcal{O}\left(\frac{L}{\sqrt{T}} + \frac{\sigma^2}{T}\right)$$

然而，NeRF优化是高度非凸的，实际收敛依赖于：
- 良好的初始化（如SIREN初始化）
- 适当的学习率调度
- 充分的网络容量

### 6.6 正则化与先验

#### 6.6.1 过拟合与泛化

NeRF容易过拟合训练视图，特别是在稀疏视图设置下。考虑经验风险最小化：

$$\hat{\Theta} = \arg\min_\Theta \frac{1}{|\mathcal{D}|} \sum_{(\mathbf{r}, C) \in \mathcal{D}} \|\hat{C}(\mathbf{r}; \Theta) - C\|^2$$

泛化误差界可以通过Rademacher复杂度分析：

$$\mathbb{E}[\mathcal{L}_{\text{test}}] \leq \mathcal{L}_{\text{train}} + 2\mathcal{R}_n(\mathcal{F}) + \mathcal{O}\left(\sqrt{\frac{\log(1/\delta)}{n}}\right)$$

其中 $\mathcal{R}_n(\mathcal{F})$ 是函数类 $\mathcal{F}$ 的Rademacher复杂度。

#### 6.6.2 几何正则化

几何正则化鼓励物理合理的密度场：

1. **深度先验**：
   $$\mathcal{L}_{\text{depth}} = \sum_{\mathbf{r}} |\hat{d}(\mathbf{r}) - d_{\text{sensor}}(\mathbf{r})|$$
   其中 $\hat{d}(\mathbf{r}) = \sum_i t_i w_i$ 是期望深度

2. **表面法线一致性**：
   $$\mathcal{L}_{\text{normal}} = \sum_{\mathbf{x}} \left\|\nabla_{\mathbf{x}} \sigma(\mathbf{x}) - \mathbf{n}_{\text{gt}}(\mathbf{x})\right\|^2$$

3. **Eikonal正则化**（用于SDF变体）：
   $$\mathcal{L}_{\text{eikonal}} = \sum_{\mathbf{x}} (|\nabla_{\mathbf{x}} f(\mathbf{x})| - 1)^2$$

#### 6.6.3 外观正则化

外观正则化改善颜色预测：

1. **视角一致性**：
   $$\mathcal{L}_{\text{view}} = \sum_{\mathbf{x}, \mathbf{d}_1, \mathbf{d}_2} w(\mathbf{d}_1, \mathbf{d}_2) \|\mathbf{c}(\mathbf{x}, \mathbf{d}_1) - \mathbf{c}(\mathbf{x}, \mathbf{d}_2)\|^2$$
   其中 $w(\mathbf{d}_1, \mathbf{d}_2)$ 根据材质属性加权

2. **语义一致性**：
   $$\mathcal{L}_{\text{semantic}} = -\sum_{\mathbf{r}, k} s_k(\mathbf{r}) \log \hat{s}_k(\mathbf{r})$$
   其中 $\hat{s}_k$ 是渲染的语义标签

#### 6.6.4 稀疏性与平滑性先验

1. **密度稀疏性**：
   $$\mathcal{L}_{\text{sparse}} = \sum_{\mathbf{x}} \rho(\sigma(\mathbf{x}))$$
   
   常用的 $\rho$ 函数：
   - L1：$\rho(x) = |x|$
   - 熵：$\rho(x) = -x\log(x + \epsilon)$
   - Cauchy：$\rho(x) = \log(1 + x^2/\epsilon^2)$

2. **空间平滑性**：
   $$\mathcal{L}_{\text{smooth}} = \sum_{\mathbf{x}} \|\nabla^2 \sigma(\mathbf{x})\|_F^2$$
   
   这促进了二阶平滑性，减少了高频噪声。

3. **信息瓶颈正则化**：
   $$\mathcal{L}_{\text{IB}} = \beta I(Z; X) - I(Z; Y)$$
   
   其中 $Z$ 是潜在表示，平衡压缩与保真度。

### 6.7 本章小结

本章介绍了神经辐射场（NeRF）的数学基础，展示了如何使用神经网络表示连续的体积密度和颜色场。关键概念包括：

1. **连续表示**：神经网络作为从3D位置和视角方向到体积属性的连续映射
2. **位置编码**：通过傅里叶特征克服神经网络的频谱偏差
3. **体积渲染**：将积分方程离散化并通过层次采样提高效率
4. **优化**：理解梯度流、收敛性和非凸优化挑战
5. **正则化**：几何和外观先验以改善泛化

核心方程汇总：
- 神经辐射场：$F_\Theta: (\mathbf{x}, \mathbf{d}) \mapsto (\sigma, \mathbf{c})$
- 体积渲染：$C(\mathbf{r}) = \int_0^\infty T(t) \sigma(\mathbf{r}(t)) \mathbf{c}(\mathbf{r}(t), \mathbf{d}) dt$
- 位置编码：$\gamma(\mathbf{x}) = [\sin(2^0\pi\mathbf{x}), \cos(2^0\pi\mathbf{x}), ..., \sin(2^{L-1}\pi\mathbf{x}), \cos(2^{L-1}\pi\mathbf{x})]^T$
- 损失函数：$\mathcal{L} = \sum_{\mathbf{r}} \|\hat{C}(\mathbf{r}) - C_{\text{gt}}(\mathbf{r})\|^2 + \lambda \mathcal{L}_{\text{reg}}$

### 6.8 练习题

#### 基础题

**练习 6.1**：推导透射率的离散近似误差界。给定连续透射率 $T(t) = \exp(-\int_0^t \sigma(s)ds)$ 和离散近似 $\hat{T}_i = \exp(-\sum_{j=1}^{i-1} \sigma_j \delta_j)$，证明：
$$|T(t_i) - \hat{T}_i| \leq T(t_i) \cdot \mathcal{O}(\max_j \delta_j^2)$$

<details>
<summary>提示</summary>
使用泰勒展开和误差累积分析。考虑 $\exp(x) \approx 1 + x + x^2/2$ 的近似。
</details>

<details>
<summary>答案</summary>

从连续到离散的误差来源于黎曼和近似：
$$\int_{t_j}^{t_{j+1}} \sigma(s)ds \approx \sigma_j \delta_j$$

使用中值定理，存在 $\xi_j \in [t_j, t_{j+1}]$ 使得：
$$\int_{t_j}^{t_{j+1}} \sigma(s)ds = \sigma(\xi_j) \delta_j$$

因此误差为：
$$e_j = |\sigma(\xi_j) - \sigma_j| \delta_j \leq \sup_s |\sigma'(s)| \cdot \delta_j^2 / 2$$

累积误差通过指数函数传播：
$$|T(t_i) - \hat{T}_i| = T(t_i) |1 - \exp(\sum_{j<i} e_j)| \leq T(t_i) \cdot \sum_{j<i} e_j \leq T(t_i) \cdot \mathcal{O}(i \cdot \max_j \delta_j^2)$$
</details>

**练习 6.2**：分析位置编码对神经网络表达能力的影响。证明使用 $L$ 个频率级别的位置编码可以精确表示频率最高为 $2^{L-1}\pi$ 的带限函数。

<details>
<summary>提示</summary>
考虑傅里叶级数展开和奈奎斯特采样定理。
</details>

<details>
<summary>答案</summary>

任何带限函数 $f$ 可以表示为：
$$f(\mathbf{x}) = \sum_{|\mathbf{k}|_\infty \leq k_{\max}} \hat{f}_\mathbf{k} e^{i\mathbf{k} \cdot \mathbf{x}}$$

使用位置编码 $\gamma(\mathbf{x})$ 包含频率 $\{2^0\pi, 2^1\pi, ..., 2^{L-1}\pi\}$。

通过欧拉公式：
$$e^{i\mathbf{k} \cdot \mathbf{x}} = \cos(\mathbf{k} \cdot \mathbf{x}) + i\sin(\mathbf{k} \cdot \mathbf{x})$$

对于 $|\mathbf{k}|_\infty \leq 2^{L-1}\pi$，存在系数使得：
$$\cos(\mathbf{k} \cdot \mathbf{x}) = \sum_{l} a_l \cos(2^l\pi \mathbf{x}) + b_l \sin(2^l\pi \mathbf{x})$$

因此，配备位置编码的线性层可以精确表示所有频率不超过 $2^{L-1}\pi$ 的函数。
</details>

**练习 6.3**：推导层次采样的方差减少。设均匀采样的估计量方差为 $\text{Var}[\hat{C}_{\text{uniform}}]$，证明层次采样满足：
$$\text{Var}[\hat{C}_{\text{hierarchical}}] \leq \text{Var}[\hat{C}_{\text{uniform}}]$$

<details>
<summary>提示</summary>
使用条件期望的性质和方差分解定理。
</details>

<details>
<summary>答案</summary>

设 $\hat{C}_c$ 为粗采样估计，$\hat{C}_f$ 为基于粗采样的细采样估计。

通过方差分解：
$$\text{Var}[\hat{C}_f] = \mathbb{E}[\text{Var}[\hat{C}_f | \hat{C}_c]] + \text{Var}[\mathbb{E}[\hat{C}_f | \hat{C}_c]]$$

由于细采样基于粗采样的重要性分布：
$$\mathbb{E}[\hat{C}_f | \hat{C}_c] = \hat{C}_c$$

因此：
$$\text{Var}[\hat{C}_f] = \mathbb{E}[\text{Var}[\hat{C}_f | \hat{C}_c]] \leq \mathbb{E}[\text{Var}[\hat{C}_{\text{uniform}}]] = \text{Var}[\hat{C}_{\text{uniform}}]$$

等号成立当且仅当重要性采样完美匹配目标分布。
</details>

#### 挑战题

**练习 6.4**：设计一个自适应位置编码方案，根据局部场景复杂度调整频率。推导选择最优频率的准则。

<details>
<summary>提示</summary>
考虑局部傅里叶分析和信号的局部带宽。
</details>

<details>
<summary>答案</summary>

定义局部复杂度度量：
$$\mathcal{C}(\mathbf{x}) = \|\nabla \sigma(\mathbf{x})\|^2 + \|\nabla \mathbf{c}(\mathbf{x})\|^2$$

自适应频率选择：
$$\mathbf{B}(\mathbf{x}) = \mathbf{B}_0 \cdot (1 + \alpha \mathcal{C}(\mathbf{x}))$$

其中 $\mathbf{B}_0$ 是基础频率矩阵，$\alpha$ 控制自适应强度。

最优准则通过最小化重构误差导出：
$$\min_{\mathbf{B}} \mathbb{E}_{\mathbf{x}} \left[\|f(\mathbf{x}) - \hat{f}_{\mathbf{B}}(\mathbf{x})\|^2\right]$$

使用变分法，最优 $\mathbf{B}(\mathbf{x})$ 满足：
$$\mathbf{B}(\mathbf{x}) \propto \text{diag}(\lambda_1(\mathbf{x}), \lambda_2(\mathbf{x}), \lambda_3(\mathbf{x}))$$

其中 $\lambda_i$ 是局部Hessian矩阵 $\nabla^2 f(\mathbf{x})$ 的特征值。
</details>

**练习 6.5**：分析NeRF中的模糊-锐利权衡。给定有限数量的训练视图，推导重建质量与细节保留之间的理论界限。

<details>
<summary>提示</summary>
使用信息论和采样理论的工具。
</details>

<details>
<summary>答案</summary>

设场景的信息内容为 $I(S)$，$N$ 个视图提供的信息为 $I(V_1, ..., V_N)$。

通过数据处理不等式：
$$I(\hat{S}; S) \leq I(V_1, ..., V_N; S)$$

对于bandlimited场景，每个视图的信息量受限于：
$$I(V_i; S) \leq W \times H \times \log(1 + \text{SNR})$$

总信息界限：
$$I(\hat{S}; S) \leq N \cdot W \times H \times \log(1 + \text{SNR})$$

细节级别 $L$ 需要的信息：
$$I_L = \mathcal{O}(L^3 \log L)$$（对于3D场景）

因此可恢复的最大细节级别：
$$L_{\max} = \mathcal{O}\left((N \cdot W \times H / \log L)^{1/3}\right)$$

这建立了视图数量、分辨率和可达到细节之间的基本权衡。
</details>

**练习 6.6**：提出一种结合物理约束的NeRF变体。推导如何将能量守恒和互易性纳入损失函数。

<details>
<summary>提示</summary>
从渲染方程的物理约束开始。
</details>

<details>
<summary>答案</summary>

能量守恒约束：
$$\int_{\mathbb{S}^2} \mathbf{c}(\mathbf{x}, \mathbf{d}_{\text{out}}) \cos\theta_{\text{out}} d\mathbf{d}_{\text{out}} \leq \int_{\mathbb{S}^2} L_i(\mathbf{x}, \mathbf{d}_{\text{in}}) \cos\theta_{\text{in}} d\mathbf{d}_{\text{in}}$$

互易性约束（对于BRDF）：
$$f_r(\mathbf{x}, \mathbf{d}_i \rightarrow \mathbf{d}_o) = f_r(\mathbf{x}, \mathbf{d}_o \rightarrow \mathbf{d}_i)$$

将这些约束编码为损失项：

1. 能量守恒损失：
$$\mathcal{L}_{\text{energy}} = \sum_{\mathbf{x}} \max\left(0, \int_{\mathbb{S}^2} \mathbf{c}(\mathbf{x}, \mathbf{d}) d\mathbf{d} - 1\right)^2$$

2. 互易性损失：
$$\mathcal{L}_{\text{reciprocity}} = \sum_{\mathbf{x}, \mathbf{d}_i, \mathbf{d}_o} \|\mathbf{c}(\mathbf{x}, \mathbf{d}_i \rightarrow \mathbf{d}_o) - \mathbf{c}(\mathbf{x}, \mathbf{d}_o \rightarrow \mathbf{d}_i)\|^2$$

总损失：
$$\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{photo}} + \lambda_1 \mathcal{L}_{\text{energy}} + \lambda_2 \mathcal{L}_{\text{reciprocity}}$$

这确保学习的表示遵守物理定律，提高泛化能力。
</details>

**练习 6.7**：开放性问题：如何将NeRF扩展到处理波动光学效应（如衍射和干涉）？概述所需的数学框架修改。

<details>
<summary>提示</summary>
考虑从标量场到复值场的转变。
</details>

<details>
<summary>答案</summary>

扩展到波动光学需要几个关键修改：

1. **复值场表示**：
   $$F_\Theta: (\mathbf{x}, \mathbf{d}, \lambda) \mapsto (A(\mathbf{x}), \phi(\mathbf{x})) \in \mathbb{C}$$
   其中 $A$ 是振幅，$\phi$ 是相位，$\lambda$ 是波长。

2. **波传播积分**：
   替换射线积分为Huygens-Fresnel积分：
   $$U(\mathbf{x}) = \frac{1}{i\lambda} \int_\Sigma U(\mathbf{x}') \frac{\exp(ikr)}{r} \cos\theta d\Sigma$$

3. **相干叠加**：
   $$I(\mathbf{x}) = |U(\mathbf{x})|^2 = |\sum_i A_i \exp(i\phi_i)|^2$$

4. **修改的神经架构**：
   - 输出复数值（实部和虚部）
   - 保持相位连续性
   - 编码波长依赖性

5. **新的损失函数**：
   $$\mathcal{L}_{\text{wave}} = \sum_{\mathbf{x}} |I_{\text{predicted}}(\mathbf{x}) - I_{\text{observed}}(\mathbf{x})|^2$$
   
   加上相位正则化：
   $$\mathcal{L}_{\text{phase}} = \sum_{\mathbf{x}} \|\nabla \phi(\mathbf{x})\|^2$$

这种方法可以捕捉衍射图案、薄膜干涉和其他波动效应，但计算成本显著增加。
</details>

### 6.9 常见陷阱与错误

1. **位置编码频率选择不当**
   - **错误**：使用过高的频率导致训练不稳定
   - **正确**：根据场景尺度选择合适的最大频率，通常 $L = 10$ 对于米级场景

2. **采样点数量不足**
   - **错误**：使用过少的采样点（如 $N < 64$）导致积分误差
   - **正确**：使用层次采样，粗采样64点+细采样128点

3. **忽视梯度消失**
   - **错误**：不使用跳跃连接导致深层梯度消失
   - **正确**：在中间层注入位置编码，使用残差连接

4. **过度正则化**
   - **错误**：过强的稀疏性约束导致几何细节丢失
   - **正确**：逐步增加正则化权重，监控验证集性能

5. **视角采样偏差**
   - **错误**：训练时未均匀采样视角方向
   - **正确**：确保训练射线覆盖所有视角范围

### 6.10 最佳实践检查清单

#### 架构设计
- [ ] 使用8层或更深的MLP以获得足够容量
- [ ] 在第4层后注入视角方向
- [ ] 使用跳跃连接改善梯度流
- [ ] 密度输出使用softplus激活
- [ ] 颜色输出使用sigmoid激活

#### 位置编码
- [ ] 根据场景尺度选择合适的频率级别
- [ ] 对位置和方向使用不同的编码策略
- [ ] 考虑使用学习的或自适应的编码

#### 采样策略
- [ ] 实现层次采样以提高效率
- [ ] 使用重要性采样聚焦于表面附近
- [ ] 添加随机偏移以避免混叠
- [ ] 验证采样密度足以捕捉细节

#### 训练过程
- [ ] 使用指数衰减的学习率调度
- [ ] 批量采样射线而非图像块
- [ ] 监控训练和验证PSNR
- [ ] 定期可视化深度图和法线图

#### 正则化
- [ ] 从低正则化开始，逐步增加
- [ ] 使用多种正则化的组合
- [ ] 根据具体应用调整权重
- [ ] 验证正则化未过度平滑细节

#### 评估
- [ ] 在保留的测试视图上评估
- [ ] 报告PSNR、SSIM和LPIPS指标
- [ ] 检查新视角的泛化能力
- [ ] 分析失败案例和边缘情况
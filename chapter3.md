# 第三章：基于点的渲染与统一体渲染方程

本章通过统一的体渲染框架，为理解所有渲染技术奠定了数学基础。我们从通用体渲染方程开始，并展示了基于点的渲染如何自然地作为连续体积表示的离散化而出现。到本章结束时，您将理解狄拉克函数分布、重建核和采样理论如何结合起来，形成从点泼溅到神经辐射场的现代渲染算法的理论基础。

## 学习目标

完成本章后，您将能够：
1. 从第一性原理推导出统一体渲染方程
2. 在体渲染框架中将点云表示为狄拉克函数分布
3. 使用频域工具分析重建质量
4. 推导离散近似的误差界限
5. 将经典点泼溅与现代可微分渲染联系起来
6. 证明高效算法的 $O(N \log N)$ 复杂度界限

## 3.1 统一体渲染方程

### 3.1.1 从表面到体积表示

经典计算机图形学传统上将表面渲染和体积渲染分开。然而，我们可以通过认识到表面是体积的极限情况来统一这些方法。考虑一个嵌入在 $\mathbb{R}^3$ 中的表面 $S$。我们可以将其表示为一个具有密度的体积：

$\sigma(\mathbf{x}) = \delta_S(\mathbf{x}) = \delta(d(\mathbf{x},S))$

其中 $d(\mathbf{x},S)$ 是到表面的有符号距离。这使我们能够统一处理所有渲染。

为了使其精确，考虑表面 $S$ 周围厚度为 $\varepsilon$ 的薄壳：

$\sigma_\varepsilon(\mathbf{x}) = \frac{1}{\varepsilon}\mathbb{1}_{|d(\mathbf{x},S)| < \varepsilon/2}$

当 $\varepsilon \to 0$ 时，$\sigma_\varepsilon \to \delta_S$ 在分布意义上。这与水平集方法相关联，其中表面是有符号距离函数的零交叉点。

**弱收敛与分布理论**：在分布的意义上，对于任何测试函数 $\varphi \in C_0^\infty(\mathbb{R}^3)$：

$\lim_{\varepsilon\to0} \int\sigma_\varepsilon(\mathbf{x})\varphi(\mathbf{x})d\mathbf{x} = \lim_{\varepsilon\to0} \frac{1}{\varepsilon}\int_{|d(\mathbf{x},S)|<\varepsilon/2} \varphi(\mathbf{x})d\mathbf{x} = \int_S \varphi(\mathbf{x})dS$

这正是表面狄拉克函数 $\delta_S$ 对 $\varphi$ 的作用。收敛可以通过共面积公式来理解：

$\int_{\mathbb{R}^3} f(\mathbf{x})\mathbb{1}_{|d(\mathbf{x},S)|<\varepsilon}d\mathbf{x} = \int_{-\varepsilon}^{\varepsilon} \int_{S_t} f(\mathbf{x})|\nabla d(\mathbf{x})|^{-1}dS_t dt$

其中 $S_t = \{\mathbf{x} : d(\mathbf{x},S) = t\}$ 是距离 $t$ 处的水平集。

**与BRDF的联系**：对于具有BRDF $f_r$ 的表面，体积发射变为：

$c(\mathbf{x},\omega) = \frac{f_r(\mathbf{x},\omega_i,\omega)L_i(\mathbf{x},\omega_i)(\mathbf{n}\cdot\omega_i)}{|\mathbf{n}\cdot\omega|}$

其中分母考虑了投影面积。这确保了体积积分恢复表面积分：

$\int \delta_S(\mathbf{x})c(\mathbf{x},\omega)d\mathbf{x} = \int_S f_r(\mathbf{x},\omega_i,\omega)L_i(\mathbf{x},\omega_i)(\mathbf{n}\cdot\omega_i)dS$

### 3.1.2 辐射传输的推导

辐射传输方程（RTE）描述了光通过参与介质的传播：

$(\omega\cdot\nabla)L(\mathbf{x},\omega) = -\sigma_t(\mathbf{x})L(\mathbf{x},\omega) + \sigma_s(\mathbf{x})\int_\Omega p(\mathbf{x},\omega',\omega)L(\mathbf{x},\omega')d\omega' + \sigma_a(\mathbf{x})L_e(\mathbf{x},\omega)$

其中：
- $L(\mathbf{x},\omega)$ 是位置 $\mathbf{x}$ 处沿方向 $\omega$ 的辐射度
- $\sigma_t = \sigma_a + \sigma_s$ 是消光系数
- $\sigma_a$ 是吸收系数
- $\sigma_s$ 是散射系数
- $p(\mathbf{x},\omega',\omega)$ 是相位函数
- $L_e$ 是发射

**微观推导**：RTE 源于粒子物理学。考虑一个体积元 $dV$，其中每单位体积有 $n(\mathbf{x})$ 个粒子，每个粒子具有以下截面：
- $\sigma_a^{(p)}$：吸收截面
- $\sigma_s^{(p)}$：散射截面
- $f(\omega',\omega)$：微分散射截面

那么：
- $\sigma_a(\mathbf{x}) = n(\mathbf{x})\sigma_a^{(p)}$ (宏观吸收)
- $\sigma_s(\mathbf{x}) = n(\mathbf{x})\sigma_s^{(p)}$ (宏观散射)
- $p(\mathbf{x},\omega',\omega) = f(\omega',\omega)/\sigma_s^{(p)}$ (归一化相位函数)

相位函数满足归一化：$\int_\Omega p(\mathbf{x},\omega',\omega)d\omega = 1$，确保能量守恒。常见的相位函数包括：
- 各向同性：$p = 1/(4\pi)$
- 瑞利：$p \propto 1 + \cos^2\theta$ (分子散射)
- Henyey-Greenstein：$p = (1-g^2)/(4\pi(1+g^2-2g\cdot\cos\theta)^{3/2})$
- Mie 理论：球形粒子的复杂振荡函数

**不对称参数**：散射角的平均余弦：
$g = \int_\Omega (\omega'\cdot\omega)p(\omega',\omega)d\omega'$

表征前向 ($g > 0$) 与后向 ($g < 0$) 散射。对于 Henyey-Greenstein， $g$ 直接参数化不对称性。

### 3.1.3 数学公式

沿射线 $\mathbf{r}(t) = \mathbf{o} + t\omega$ 从 $t=0$ 到 $t=T$ 积分，我们使用特征线法求解 RTE。定义光学深度：

$\tau(s,t) = \int_s^t \sigma_t(\mathbf{r}(u))du$

透射率 $T(s,t) = \exp(-\tau(s,t))$ 表示从 $s$ 到 $t$ 存活的光的比例。

**通过积分因子形式解**：将 RTE 乘以 $\exp(\int_0^t \sigma_t(\mathbf{r}(u))du)$：

$\frac{d}{dt}[L(\mathbf{r}(t),\omega)\exp(\tau(0,t))] = \exp(\tau(0,t))[\sigma_s S_s + \sigma_a L_e]$

其中 $S_s(\mathbf{x},\omega) = \int_\Omega p(\mathbf{x},\omega',\omega)L(\mathbf{x},\omega')d\omega'$ 是入散射辐射度。

从 $0$ 到 $T$ 积分：

$L(\mathbf{o},\omega) = \int_0^T T(0,t)\sigma_t(\mathbf{r}(t))S(\mathbf{r}(t),\omega)dt + T(0,T)L_{bg}$

其中源项 $S$ 结合了发射和入散射：

$S(\mathbf{x},\omega) = \frac{\sigma_a(\mathbf{x})L_e(\mathbf{x},\omega)}{\sigma_t(\mathbf{x})} + \frac{\sigma_s(\mathbf{x})}{\sigma_t(\mathbf{x})}\int_\Omega p(\mathbf{x},\omega',\omega)L(\mathbf{x},\omega')d\omega'$

**单次散射近似**：假设入散射积分中的 $L$ 仅为直接照明：

$S_s^{(1)}(\mathbf{x},\omega) = \int_\Omega p(\mathbf{x},\omega',\omega)L_{direct}(\mathbf{x},\omega')d\omega'$

其中 $L_{direct}(\mathbf{x},\omega') = T(\mathbf{x}_{light},\mathbf{x})L_e(\mathbf{x}_{light},-\omega')V(\mathbf{x},\mathbf{x}_{light})$。

对于纯发射介质（无散射），这简化为：

$L(\mathbf{o},\omega) = \int_0^T T(t)\sigma(\mathbf{r}(t))c(\mathbf{r}(t),\omega)dt + T(T)L_{bg}$

其中：
- $T(t) = \exp(-\int_0^t \sigma(\mathbf{r}(s))ds)$ 是从原点到 $t$ 的透射率
- $c(\mathbf{x},\omega) = L_e(\mathbf{x},\omega)$ 是发射辐射度
- $L_{bg}$ 是背景辐射度

这个方程统一了所有渲染：表面具有狄拉克函数形式的 $\sigma$，体积具有连续的 $\sigma$。

**算子形式**：定义传输算子 $\mathcal{T}$ 和散射算子 $\mathcal{S}$：
- $(\mathcal{T}L)(\mathbf{x},\omega) = (\omega\cdot\nabla)L(\mathbf{x},\omega) + \sigma_t(\mathbf{x})L(\mathbf{x},\omega)$
- $(\mathcal{S}L)(\mathbf{x},\omega) = \sigma_s(\mathbf{x})\int_\Omega p(\mathbf{x},\omega',\omega)L(\mathbf{x},\omega')d\omega'$

那么 RTE 变为：$\mathcal{T}L = \mathcal{S}L + Q$，其中 $Q = \sigma_a L_e$ 是源。

### 3.1.4 与经典渲染的联系

对于距离 $t^*$ 处的表面，沿射线 $\sigma(\mathbf{x}) = \delta(t-t^*)$，透射率变为：

$T(t) = \begin{cases} 1 & \text{if } t < t^* \\ 0 & \text{if } t > t^* \end{cases}$

这是一个阶跃函数。体积积分使用狄拉克函数的筛选性质进行评估：

$L(\mathbf{o},\omega) = \int_0^T T(t)\delta(t-t^*)c(\mathbf{r}(t),\omega)dt + T(T)L_{bg}$
$\quad = T(t^*)c(\mathbf{r}(t^*),\omega) + T(T)L_{bg}$
$\quad = 1\cdot c(\mathbf{r}(t^*),\omega) + 0\cdot L_{bg}$
$\quad = c(\mathbf{r}(t^*),\omega)$

这恢复了在表面交点处评估的经典渲染方程。BRDF 通过 $c(\mathbf{r}(t^*),\omega) = \int f_r(\mathbf{x},\omega_i,\omega_o)L_i(\mathbf{x},\omega_i)(\mathbf{n}\cdot\omega_i)d\omega_i$ 出现。

### 3.1.5 边界条件与适定性

体渲染方程需要边界条件以实现数学完整性：

1. **真空边界**：对于边界上的 $\mathbf{x}$，$\omega$ 指向内部时，$L(\mathbf{x},\omega) = L_{bg}$
2. **发射边界**：$L(\mathbf{x},\omega) = L_e(\mathbf{x},\omega)$
3. **反射边界**：$L(\mathbf{x},\omega) = \int f_r(\mathbf{x},\omega',\omega)L(\mathbf{x},\omega')(\mathbf{n}\cdot\omega')d\omega'$

**数学框架**：RTE 及其边界条件构成一个抽象柯西问题：

$L + \mathcal{K}L = f \text{ in } \Omega\times S^2$
$L|_{\Gamma_-} = g$

其中：
- $\mathcal{K}$ 是积分散射算子
- $\Gamma_- = \{(\mathbf{x},\omega) \in \partial\Omega\times S^2 : \mathbf{n}(\mathbf{x})\cdot\omega < 0\}$ 是流入边界
- $f$ 代表源，$g$ 代表边界数据

在 $\sigma$ 和 $c$ 的温和条件下，该方程在 $L^2(\Omega\times S^2)$ 中是适定的。

**定理（存在性与唯一性）**：如果：
1. $\sigma_t \in L^\infty(\Omega)$, $\sigma_t \ge \sigma_{min} > 0$
2. $||\sigma_s/\sigma_t||_\infty < 1$ (亚临界条件)
3. $p \in L^\infty(\Omega\times S^2\times S^2)$, $p \ge 0$

那么存在唯一的解 $L \in L^2(\Omega\times S^2)$ 满足：
$||L||_2 \le C(||f||_2 + ||g||_{L^2(\Gamma_-)})$

**Fredholm 择一性定理**：当谱半径 $\rho(\mathcal{K}) < 1$ 时，算子 $(I - \mathcal{K})$ 是可逆的。对于均匀介质：
$\rho(\mathcal{K}) = \sigma_s/\sigma_t$

这给出了临界反照率 $\sigma_s/\sigma_t = 1$，高于此值介质可以通过散射维持自发射。

### 3.1.6 能量守恒与互易性

体渲染方程保留了两个基本物理原理：

**能量守恒**：总输入功率等于总输出功率
$\int_{\partial\Omega}\int_{S^2} L(\mathbf{x},\omega)(\mathbf{n}\cdot\omega)d\omega dA = \int_\Omega\int_{S^2} \sigma_a(\mathbf{x})L_e(\mathbf{x},\omega)d\omega dV$

证明：将 RTE 乘以 $1$ 并在 $\Omega\times S^2$ 上积分：
$\int_\Omega\int_{S^2} (\omega\cdot\nabla)L d\omega dV = -\int_\Omega\int_{S^2} \sigma_t L d\omega dV + \int_\Omega\int_{S^2} \sigma_s(\int p L'd\omega')d\omega dV + \int_\Omega\int_{S^2} \sigma_a L_e d\omega dV$

对左侧使用散度定理：
$\int_{\partial\Omega}\int_{S^2} L(\mathbf{n}\cdot\omega)d\omega dA = -\int_\Omega\int_{S^2} \sigma_a L d\omega dV + \int_\Omega\int_{S^2} \sigma_a L_e d\omega dV$

因为 $\int\int p(\omega',\omega)d\omega = 1$ 使散射项消失。

**Helmholtz 互易性**：对于互易介质 ($p(\mathbf{x},\omega',\omega) = p(\mathbf{x},\omega,\omega')$)：
如果 $L_1$ 是源在 $\mathbf{x}_1$ 指向 $\mathbf{x}_2$ 的解，而 $L_2$ 是源在 $\mathbf{x}_2$ 指向 $\mathbf{x}_1$ 的解，那么 $L_1(\mathbf{x}_2,-\omega) = L_2(\mathbf{x}_1,-\omega)$。

这源于伴随 RTE：
$(-\omega\cdot\nabla)L^* + \sigma_t L^* = \sigma_s \int p(\omega,\omega')L^*(\omega')d\omega' + Q^*$

满足互易性的格林函数 $G(\mathbf{x},\omega;\mathbf{x}',\omega')$ 使得路径积分公式成为可能：
$L(\mathbf{x},\omega) = \iint G(\mathbf{x},\omega;\mathbf{x}',\omega')Q(\mathbf{x}',\omega')d\mathbf{x}'d\omega'$

**详细平衡**：在温度 $T$ 的热平衡状态下：
$\sigma_a(\mathbf{x})B(T) = \sigma_a(\mathbf{x})L_e(\mathbf{x},\omega)$

其中 $B(T)$ 是普朗克函数，确保微观可逆性。

## 3.2 点云作为狄拉克函数分布

### 3.2.1 数学基础

点云 $P = \{(\mathbf{p}_i, \mathbf{a}_i)\}_{i=1}^N$，其中 $\mathbf{p}_i \in \mathbb{R}^3$ 是位置，$\mathbf{a}_i$ 是属性（颜色、法线等），表示一个分布：

$\sigma(\mathbf{x}) = \sum_{i=1}^N w_i\delta(\mathbf{x} - \mathbf{p}_i)$
$c(\mathbf{x},\omega) = \frac{\sum_{i=1}^N w_i\delta(\mathbf{x} - \mathbf{p}_i)}{\sum_{j}w_j\delta(\mathbf{x} - \mathbf{p}_j)} \cdot c_i(\omega)$

其中 $w_i$ 是权重，$c_i(\omega)$ 编码了点的外观。

**施瓦茨分布理论**：这种表示在分布（广义函数）的意义上是严格的。分布空间 $\mathcal{D}'(\mathbb{R}^3)$ 是测试函数空间 $\mathcal{D}(\mathbb{R}^3) = C_0^\infty(\mathbb{R}^3)$ 的对偶。对于任何测试函数 $\varphi \in C_0^\infty(\mathbb{R}^3)$：

$\langle\sigma, \varphi\rangle = \int\sigma(\mathbf{x})\varphi(\mathbf{x})d\mathbf{x} = \sum_i w_i\varphi(\mathbf{p}_i)$

狄拉克函数满足：
1. **筛选性质**：$\int\delta(\mathbf{x}-\mathbf{a})f(\mathbf{x})d\mathbf{x} = f(\mathbf{a})$
2. **缩放**：$\delta(a\mathbf{x}) = |a|^{-3}\delta(\mathbf{x})$ 对于 $a \neq 0$
3. **导数**：$\langle\partial^\alpha \delta_\mathbf{a}, \varphi\rangle = (-1)^{|\alpha|}\partial^\alpha \varphi(\mathbf{a})$
4. **傅里叶变换**：$\mathcal{F}[\delta_\mathbf{a}](\mathbf{k}) = \exp(-i\mathbf{k}\cdot\mathbf{a})$

**正则化序列**：狄拉克函数是正则函数的极限：
$\delta(\mathbf{x}) = \lim_{\varepsilon\to0} \delta_\varepsilon(\mathbf{x})$

常见的正则化：
1. 高斯：$\delta_\varepsilon(\mathbf{x}) = (2\pi\varepsilon^2)^{-3/2}\exp(-|\mathbf{x}|^2/(2\varepsilon^2))$
2. 矩形：$\delta_\varepsilon(\mathbf{x}) = (1/\varepsilon^3)\mathbb{1}_{|\mathbf{x}|<\varepsilon/2}$
3. Sinc：$\delta_\varepsilon(\mathbf{x}) = (1/2\pi)^3\int_{|\mathbf{k}|<1/\varepsilon} \exp(i\mathbf{k}\cdot\mathbf{x})d\mathbf{k}$

每个都在 $\mathcal{D}'(\mathbb{R}^3)$ 的弱*拓扑中收敛到 $\delta$。

### 3.2.2 连续场的离散采样

点云源于对连续场的采样。给定连续密度 $\sigma_c(\mathbf{x})$ 和采样点 $\{\mathbf{x}_i\}$，离散近似为：

$\sigma_d(\mathbf{x}) = \sum_i \sigma_c(\mathbf{x}_i)V_i \delta(\mathbf{x} - \mathbf{x}_i)$

其中 $V_i$ 是与样本 $i$ 关联的体积。常见的体积分配：

1. **均匀采样**：对于规则网格，$V_i = \Delta x^3$
2. **Voronoi 单元**：$V_i = \int_{V(\mathbf{x}_i)} d\mathbf{x}$，其中 $V(\mathbf{x}_i) = \{\mathbf{x} : |\mathbf{x}-\mathbf{x}_i| < |\mathbf{x}-\mathbf{x}_j| \forall j\neq i\}$
3. **Delaunay 对偶**：$V_i = (1/3)\sum_{T\in D(i)} \text{Vol}(T)$ 对于包含 $i$ 的四面体
4. **自适应采样**：$V_i \propto$ 局部特征尺寸

**Voronoi 体积计算**：对于点 $\mathbf{p}_i$ 及其邻居 $\{\mathbf{p}_j\}$，Voronoi 单元是：
$V(\mathbf{p}_i) = \bigcap_{j\neq i} \{\mathbf{x} : (\mathbf{x}-\mathbf{p}_i)\cdot(\mathbf{p}_j-\mathbf{p}_i) < |\mathbf{p}_j-\mathbf{p}_i|^2/2\}$

体积积分：
$V_i = \int_{V(\mathbf{p}_i)} d\mathbf{x}$

对于半径为 $r$ 的泊松盘分布：
$E[V_i] \approx (4/3)\pi r^3 \cdot 0.74$ (最优填充密度)

**采样算子性质**：采样算子 $S$ 将连续映射到离散：
$S: L^1(\mathbb{R}^3) \to \mathcal{D}'(\mathbb{R}^3)$
$S[\sigma_c] = \sum_i\sigma_c(\mathbf{x}_i)V_i\delta(\mathbf{x}-\mathbf{x}_i)$

性质：
1. **线性**：$S[a\sigma_1 + b\sigma_2] = aS[\sigma_1] + bS[\sigma_2]$
2. **质量守恒**：$\int S[\sigma_c]d\mathbf{x} = \sum_i\sigma_c(\mathbf{x}_i)V_i \approx \int\sigma_c d\mathbf{x}$ (对于单位分解)
3. **频率响应**：$\mathcal{F}[S[\sigma_c]](\mathbf{k}) = \sum_i\sigma_c(\mathbf{x}_i)V_i \exp(-i\mathbf{k}\cdot\mathbf{x}_i)$

### 3.2.3 重建理论

为了渲染点云，我们必须从离散样本重建连续场。重建使用与核 $h$ 的卷积：

$\sigma_r(\mathbf{x}) = (\sigma_d * h)(\mathbf{x}) = \sum_i w_i h(\mathbf{x} - \mathbf{p}_i)$

重建算子 $R$ 满足：$R[\sigma_d] = \sigma_d * h$。组合采样和重建：

$\sigma_r = R[S[\sigma_c]] = \sum_i\sigma_c(\mathbf{x}_i)V_ih(\mathbf{x} - \mathbf{x}_i)$

**香农-惠特克定理**：对于带限信号 $\sigma_c$，其中 $\hat{\sigma}_c(\mathbf{k}) = 0$ 对于 $|\mathbf{k}| > K$：

$\sigma_c(\mathbf{x}) = \sum_i \sigma_c(\mathbf{x}_i)\text{sinc}(K(\mathbf{x} - \mathbf{x}_i)/\pi)$

当样本位于间距为 $\Delta x = \pi/K$ 的网格上时。sinc 核：
$\text{sinc}(x) = \sin(|x|)/|x|$ (1D), $\text{sinc}(\mathbf{x}) = (\sin(|\mathbf{x}|) - |\mathbf{x}|\cos(|\mathbf{x}|))/|\mathbf{x}|^3$ (3D)

完美重建需要 $RS = I$ (恒等算子)。这发生在以下情况：
1. $h$ 是理想的 sinc 核
2. 采样满足奈奎斯特准则：$\Delta x < \pi/K$
3. 信号是带限的：$\text{supp}(\hat{\sigma}_c) \subset B_K(0)$

**近似理论**：对于非带限信号，我们最小化重建误差：

$E = ||\sigma_c - RS[\sigma_c]||^2_{L^2}$

在 $L^2$ 意义上的最优核满足正规方程：
$\sum_j\langle h(\cdot - \mathbf{x}_i), h(\cdot - \mathbf{x}_j)\rangle w_j = \sigma_c(\mathbf{x}_i)$

这导致了对偶核公式：
$\tilde{h}(\mathbf{x}) = \sum_i \alpha_i h(\mathbf{x} - \mathbf{x}_i)$

其中 $\alpha$ 求解 $G\alpha = \sigma$，且 $G_{ij} = h(\mathbf{x}_i - \mathbf{x}_j)$。

### 3.2.4 混叠与采样定理

根据奈奎斯特-香农定理，完美重建需要：
1. 带限信号：$\hat{\sigma}_c(k) = 0$ 对于 $|k| > k_{max}$
2. 采样率：$\Delta x < \pi/k_{max}$

对于非带限信号，我们通过傅里叶分析来分析混叠误差。采样信号的频谱为：

$$\hat{\sigma}_d(k) = (1/V_s)\Sigma_n \hat{\sigma}_c(k - 2\pi n/\Delta x)$$

其中 $V_s = \Delta x^3$ 是采样体积。当频谱重叠时发生混叠：

$$E_{alias} = \int_{|k|>\pi/\Delta x} |\hat{\sigma}_c(k)|^2 dk$$

对于具有幂律频谱 $\hat{\sigma}_c(k) \sim |k|^{-\alpha}$ 的信号，混叠误差的缩放关系为：
$E_{alias} \sim \Delta x^{(2\alpha-6)}$ 对于 $\alpha > 3$

### 3.2.5 不规则采样与抖动网格

规则采样会产生结构化的混叠伪影。不规则采样将混叠转换为噪声：

**泊松盘采样**：点满足最小距离约束
- 任意两点之间的距离不小于 $r_{min}$
- 频谱具有“蓝噪声”特性：$\hat{\sigma}(k) \approx 0$ 对于 $|k| < k_{min}$

**抖动采样**：扰动规则网格
$x_{ijk} = (i,j,k)\Delta x + \xi_{ijk}$

其中 $\xi_{ijk} \sim U[-\Delta x/2, \Delta x/2]^3$。这在保持覆盖的同时打破了规则性。

**频谱分析**：对于抖动采样，期望频谱为：
$$E[|\hat{\sigma}_d(k)|^2] = |\hat{\sigma}_c(k)|^2 + (1-\text{sinc}^2(k\Delta x/2))\Sigma_{n\neq 0}|\hat{\sigma}_c(k-2\pi n/\Delta x)|^2$$

$\text{sinc}^2$ 项抑制了混叠，相比于规则采样。

### 3.2.6 与测度论的联系

点云定义了 $\mathbb{R}^3$ 上的原子测度：

$$\mu = \Sigma_i w_i\delta_{p_i}$$

对于任意 Borel 集 $B \subseteq \mathbb{R}^3$：
$$\mu(B) = \Sigma_{i:p_i\in B} w_i$$

这种测度论观点与以下方面相关：
- 点云匹配的最优传输
- 形状比较的 Wasserstein 距离
- 点云演化的梯度流

总变差范数 $||\mu||_{TV} = \Sigma_i|w_i|$ 限制了点云的“质量”。

## 3.3 泼溅核与重建滤波器

### 3.3.1 核设计原则

理想的重建核应满足多个数学和实际约束：

1. **紧支集**：$\text{supp}(h) \subset B_R(0)$ 以提高效率
2. **平滑性**：$h \in C^n$ 以获得视觉质量（优选 $n \ge 2$）
3. **单位分解**：$\Sigma_i h(x - p_i) \approx 1$ 对于所有 $x$
4. **矩保持**：$\int x^\alpha h(x)dx = \delta_{|\alpha|,0}$ 对于 $|\alpha| \le m$
5. **非负性**：$h(x) \ge 0$（防止负密度）
6. **归一化**：$\int h(x)dx = 1$（质量守恒）

单位分解确保了常数重建：如果 $\sigma_c(x) = c$，那么 $\sigma_r(x) = c$。

**定理**：没有紧支集核可以是 $C^\infty$ 且具有紧傅里叶变换。

这一基本限制迫使核设计中进行权衡。

### 3.3.2 高斯核

高斯核在基于点的渲染中无处不在：

$$h_G(x) = (2\pi\sigma^2)^{-3/2} \exp(-|x|^2/2\sigma^2)$$

优点：
- 平滑 ($C^\infty$)
- 可分离：$h_G(x,y,z) = h_{1D}(x)h_{1D}(y)h_{1D}(z)$
- 卷积下封闭：$h_G^{\sigma_1} * h_G^{\sigma_2} = h_G^{\sqrt{\sigma_1^2+\sigma_2^2}}$
- 最佳时频局部化（最小化海森堡不确定性）
- 旋转不变：$h_G(Rx) = h_G(x)$ 对于旋转 $R$

傅里叶变换：
$$\hat{h}_G(k) = \exp(-|k|^2\sigma^2/2)$$

高斯满足扩散方程：
$\partial h_G/\partial t = \frac{1}{2}\Delta h_G$ 且 $h_G(x,0) = \delta(x)$

这连接了泼溅到尺度空间理论和扩散过程。

**截断高斯**：为了效率，在半径 $r = n\sigma$ 处截断（通常 $n = 3$）：

$$h_T(x) = \begin{cases} C \exp(-|x|^2/2\sigma^2) & \text{if } |x| < n\sigma \\ 0 & \text{otherwise} \end{cases}$$

其中 $C$ 确保 $\int h_T = 1$。截断误差为：

$E_{trunc} = 1 - \text{erf}(n/\sqrt{2}) \approx 2.7\times 10^{-3}$ 对于 $n = 3$

### 3.3.3 各向异性核

对于有向表面，各向异性高斯能更好地捕捉局部几何：

$$h_A(x) = (2\pi)^{-3/2}|\Sigma|^{-1/2} \exp(-\frac{1}{2}x^T\Sigma^{-1}x)$$

其中 $\Sigma$ 是 $3\times 3$ 协方差矩阵。特征分解揭示几何：

$$\Sigma = RSR^T = R \text{ diag}(\lambda_1, \lambda_2, \lambda_3) R^T$$

- $R$：旋转矩阵（主轴）
- $\lambda_i$：特征值（沿轴的平方半径）

对于表面泼溅，通常 $\lambda_3 \ll \lambda_1, \lambda_2$，创建盘状泼溅。

**协方差估计** 来自局部点邻域：

$$\Sigma = (1/k)\Sigma_{i=1}^k (p_i - \bar{p})(p_i - \bar{p})^T$$

其中 $\bar{p}$ 是邻域质心。这是经验协方差。

**表面对齐泼溅**：给定表面法线 $n$，构造：

$$\Sigma = \sigma_\parallel^2(I - nn^T) + \sigma_\perp^2 nn^T$$

其中 $\sigma_\parallel \gg \sigma_\perp$ 用于薄表面。

### 3.3.4 频域分析

重建质量取决于核的频率响应。重建频谱：

$$\hat{\sigma}_r(k) = \hat{\sigma}_d(k)\hat{h}(k) = [\Sigma_n\hat{\sigma}_c(k - 2\pi n/\Delta x)]\hat{h}(k)$$

理想低通滤波器：
$$\hat{h}_{ideal}(k) = \mathbb{1}_{|k|<k_c}(k)$$

其空间表示（sinc 核）：
$$h_{ideal}(x) = (k_c/2\pi)^3 \cdot (\sin(k_c|x|) - k_c|x|\cos(k_c|x|))/(k_c|x|)^3$$

但 sinc 具有无限支集和缓慢衰减 ($O(|x|^{-1})$)。实际核以紧支集近似理想响应。

**滤波器质量指标**：
1. **通带纹波**：$\max_{|k|<k_c} |1 - \hat{h}(k)|$
2. **阻带衰减**：$\max_{|k|>k_s} |\hat{h}(k)|$
3. **过渡宽度**：$k_s - k_c$

### 3.3.5 替代核族

**B 样条**：分段多项式核

$$B^n(x) = (B^{(n-1)} * B^0)(x)$$

其中 $B^0 = \mathbb{1}_{[-1/2,1/2]}$ 是盒函数。三次 B 样条：

$$B^3(x) = \begin{cases} (2-|x|)^3/6, & 1 \le |x| \le 2 \\ 2/3 - |x|^2 + |x|^3/2, & |x| < 1 \\ 0, & |x| > 2 \end{cases}$$

属性：
- 紧支集：$\text{supp}(B^n) = [-(n+1)/2, (n+1)/2]$
- 平滑性：$B^n \in C^{(n-1)}$
- 精确多项式再现，最高可达 $n$ 次

**Wendland 核**：紧支集径向基函数

$$\psi_{\ell,k}(r) = \begin{cases} p_{\ell,k}(r) & \text{if } r \le 1 \\ 0 & \text{otherwise} \end{cases}$$

其中 $p_{\ell,k}$ 是多项式。示例 ($\psi_{3,1}$):

$$\psi_{3,1}(r) = (1-r)^4_+(4r+1)$$

这些核在给定支集下实现了最佳平滑性。

**Kaiser-Bessel 窗**：近似最佳集中度

$$h_{KB}(x) = \begin{cases} I_0(\beta\sqrt{1-(2x/w)^2})/I_0(\beta) & \text{if } |x| < w/2 \\ 0 & \text{otherwise} \end{cases}$$

其中 $I_0$ 是修正贝塞尔函数。参数 $\beta$ 控制主瓣宽度和旁瓣抑制之间的权衡。

### 3.3.6 核选择指南

根据应用需求选择核：

1. **质量优先**：高斯或 Kaiser-Bessel
2. **速度优先**：截断高斯或低阶 B 样条
3. **精确插值**：径向基函数
4. **硬件泼溅**：屏幕对齐椭圆
5. **薄表面**：各向异性高斯，其中 $\sigma_\perp \to 0$

核带宽 $\sigma$ 应与采样密度相关：
- 密集采样：$\sigma \approx 0.5 \times$ 平均邻居距离
- 稀疏采样：$\sigma \approx 1.5 \times$ 平均邻居距离

基于局部密度的自适应带宽：
$$\sigma(x) = \sigma_0(\rho(x)/\rho_0)^{-1/3}$$

其中 $\rho(x)$ 是局部点密度。

## 3.4 离散样本的体渲染方程

### 3.4.1 离散化策略

给定点云表示 $\sigma(x) = \Sigma_i w_i h(x - p_i)$，体渲染积分变为：

$$L(o,\omega) = \int_0^T T(t)\Sigma_i w_i h(r(t) - p_i)c_i(\omega)dt$$

重新排列：
$$L(o,\omega) = \Sigma_i w_i c_i(\omega)\int_0^T T(t)h(r(t) - p_i)dt$$

### 3.4.2 求积规则与误差界

对于数值积分，将射线离散为 $M$ 个段：

$$L(o,\omega) \approx \Sigma_{j=1}^M \Delta t_j T(t_j)\Sigma_i w_i h(r(t_j) - p_i)c_i(\omega)$$

使用梯形法则，对于平滑核，误差为 $O(\Delta t^2)$。对于标准差为 $\sigma$ 的高斯核，自适应求积可以实现误差 $\varepsilon$，其中 $M = O(\sigma^{-1}\log(1/\varepsilon))$ 样本。

### 3.4.3 Alpha 合成作为特例

考虑沿射线在 $\{t_j\}$ 处的离散样本，不透明度为 $\alpha_j$，颜色为 $c_j$。设置：
- $\sigma(x) = \Sigma_j(-\log(1-\alpha_j)/\Delta t)\delta(t-t_j)$
- $c(x,\omega) = c_j$ 对于 $x \in [t_j, t_{j+1})$

体渲染方程得到：

$$L = \Sigma_j c_j \alpha_j \prod_{k<j}(1-\alpha_k)$$

这是经典的 alpha 合成公式，表明它是体渲染的一个特例。

### 3.4.4 与粒子系统的联系

对于半径为 $r_i$ 和密度为 $\rho_i$ 的粒子系统：

$$\sigma(x) = \Sigma_i \rho_i \mathbb{1}_{|x-p_i|<r_i}$$

其中 $\mathbb{1}$ 是指示函数。体渲染方程变为：

$$L(o,\omega) = \Sigma_{i\in I} \rho_i c_i(\omega)|r \cap B_i| \prod_{j\in J,j<i} \exp(-\rho_j|r \cap B_j|)$$

其中 $I$ 是相交粒子，$B_i$ 是粒子 $i$ 的球体，$J$ 是相机和 $i$ 之间的粒子。

## 3.5 误差分析与收敛性

### 3.5.1 基于点渲染的近似理论

设 $f$ 为真实的连续场，$f_N$ 为其 $N$ 点近似。在 $L^2$ 范数下的近似误差：

$$||f - f_N||_2^2 = \int|f(x) - \Sigma_{i=1}^N w_i h(x - p_i)|^2 dx$$

对于最佳点放置和权重（最小化上述误差），误差的缩放关系为：

$$||f - f_N||_2 = O(N^{-s/d})$$

其中 $s$ 是 $f$ 的平滑度（Sobolev 正则性），$d$ 是维度（3D 为 $d=3$）。

### 3.5.2 收敛速度：数学界限

对于特定的核选择：

**定理 3.1** (高斯核收敛性)：设 $f \in H^s(\mathbb{R}^3)$ 且 $s > 3/2$。使用带宽 $h = O(N^{-1/3})$ 的高斯核和准均匀点分布：

$$||f - f_N||_\infty \le CN^{-s/3} + C'N^{-1/2}\log N$$

第一项是近似误差，第二项是点放置的随机误差。

**定理 3.2** (最优核)：对于带限 $f$ 且 $||\hat{f}||_\infty = 0$ 对于 $|k| > K$，使用 sinc 核：

$$||f - f_N||_2 = 0$$

当点位于间距 $\Delta x < \pi/K$ 的网格上时（精确重建）。

### 3.5.3 计算复杂度：$O(N \log N)$ 算法

朴素泼溅是 $N$ 个点和 $M$ 个像素的 $O(NM)$。高效算法实现 $O(N \log N)$：

1. **分层泼溅**：构建八叉树，独立泼溅各层
   - 树构建：$O(N \log N)$
   - 带截止的泼溅：$O(N \log N)$

2. **傅里叶泼溅**：对于周期域
   - 点样本的 FFT：$O(N \log N)$
   - 频域卷积：$O(N)$
   - 逆 FFT：$O(N \log N)$

3. **快速多极方法**：对于长程核
   - 多极展开：$O(N)$
   - 平移算子：$O(N \log N)$

### 3.5.4 实际误差度量

对于渲染图像，考虑：

1. **PSNR**：$20\log_{10}(\text{MAX}/\text{RMSE})$ 其中 $\text{RMSE} = \sqrt{1/M \Sigma(I - I_{ref})^2}$

2. **结构相似性 (SSIM)**：考虑人类感知

3. **豪斯多夫距离**：用于几何精度
   $d_H(S,S') = \max(\sup_{x\in S} \inf_{y\in S'} |x-y|, \sup_{y\in S'} \inf_{x\in S} |x-y|)$
对于点云，到参考曲面的单侧距离：
$E_{geo} = \frac{1}{N} \sum_{i} d(p_i, S_{ref})$

## 3.6 章总结

我们建立了统一体渲染方程作为理解所有渲染技术的基础。主要见解：

1.  **统一性**：表面渲染和体渲染是通用体渲染方程的特例
2.  **点云**：自然表示为狄拉克函数的加权和
3.  **重建**：与核函数卷积将离散样本转换为连续场
4.  **效率**：分层和频域方法实现 $O(N \log N)$ 复杂度
5.  **误差分析**：收敛速度取决于信号平滑度和采样密度

这里开发的数学框架直接扩展到：
- 基于图像的渲染（第4章）：光场采样
- 神经辐射场（第6章）：作为神经网络的连续密度/颜色
- 3D 高斯泼溅（第10章）：各向异性高斯核

## 练习

### 练习 3.1 (基础)
证明高斯卷积保留了点云的第一矩（质心）。

**提示**：使用性质 $\int x h_G(x)dx = 0$ 对于中心高斯函数。

<details>
<summary>解决方案</summary>

设点云有位于 $\{p_i\}$ 的点，权重为 $\{w_i\}$。质心为：
$C = \frac{\sum_i w_i p_i}{\sum_i w_i}$

与高斯核 $h_G$ 卷积后：
$f(x) = \sum_i w_i h_G(x - p_i)$

第一矩：
$M_1 = \int x f(x)dx = \int x \sum_i w_i h_G(x - p_i)dx$
$\quad = \sum_i w_i \int x h_G(x - p_i)dx$

代入 $y = x - p_i$：
$M_1 = \sum_i w_i \int (y + p_i) h_G(y)dy$
$\quad = \sum_i w_i [\int y h_G(y)dy + p_i \int h_G(y)dy]$
$\quad = \sum_i w_i [0 + p_i \cdot 1]$
$\quad = \sum_i w_i p_i$

因此 $M_1 / \int f(x)dx = \sum_i w_i p_i / \sum_i w_i = C \quad \checkmark$
</details>

### 练习 3.2 (基础)
证明对于盒式滤波样本，体渲染方程可以精确地恢复 alpha 合成。

**提示**：使用 $\sigma(t) = \sum_j \sigma_j \mathbb{1}_{[t_j,t_{j+1}]}(t)$ 并分段计算 $T(t)$。

<details>
<summary>解决方案</summary>

给定位于 $\{t_j\}$ 的样本，不透明度 $\alpha_j = 1 - \exp(-\sigma_j \Delta t)$，密度为：
$\sigma(t) = \sum_j (-\log(1-\alpha_j)/\Delta t) \mathbb{1}_{[t_j,t_{j+1}]}(t)$

到 $t_k$ 的透射率：
$T(t_k) = \exp(-\int_0^{t_k} \sigma(s)ds)$
$\quad = \exp(-\sum_{j<k} \int_{t_j}^{t_{j+1}} \sigma_j ds)$
$\quad = \exp(-\sum_{j<k} \sigma_j \Delta t)$
$\quad = \prod_{j<k} \exp(-\sigma_j \Delta t)$
$\quad = \prod_{j<k} (1-\alpha_j)$

来自区间 $k$ 的贡献：
$L_k = \int_{t_k}^{t_{k+1}} T(t)\sigma(t)c(t)dt$
$\quad = T(t_k)\sigma_k c_k \Delta t$
$\quad = \prod_{j<k}(1-\alpha_j) \cdot (-\log(1-\alpha_k)/\Delta t) \cdot c_k \cdot \Delta t$
$\quad = \prod_{j<k}(1-\alpha_j) \cdot (-\log(1-\alpha_k)) \cdot c_k$

使用 $1-\exp(-x) \approx x$ 对于小 $x$，或精确地：
$-\log(1-\alpha_k) = -\log(\exp(-\sigma_k \Delta t)) = \sigma_k \Delta t$

所以：$L_k = \prod_{j<k}(1-\alpha_j) \cdot \alpha_k/(1-\alpha_k) \cdot (1-\alpha_k) \cdot c_k = \alpha_k c_k \prod_{j<k}(1-\alpha_j) \quad \checkmark$
</details>

### 练习 3.3 (基础)
推导具有协方差 $\Sigma$ 的各向异性高斯核的傅里叶变换。

**提示**：使用代换 $y = \Sigma^{-1/2}x$ 简化为各向同性情况。

<details>
<summary>解决方案</summary>

各向异性高斯函数：
$h(x) = (2\pi)^{-d/2}|\Sigma|^{-1/2}\exp(-\frac{1}{2}x^\top\Sigma^{-1}x)$

其傅里叶变换：
$\hat{h}(k) = \int h(x)\exp(-ik \cdot x)dx$

代入 $y = \Sigma^{-1/2}x$，所以 $x = \Sigma^{1/2}y$ 且 $dx = |\Sigma|^{1/2}dy$：

$\hat{h}(k) = (2\pi)^{-d/2}|\Sigma|^{-1/2}\int\exp(-\frac{1}{2}y^\top y)\exp(-ik \cdot \Sigma^{1/2}y)|\Sigma|^{1/2}dy$
$\quad = (2\pi)^{-d/2}\int\exp(-\frac{1}{2}y^\top y)\exp(-i(\Sigma^{1/2 \top}k) \cdot y)dy$

这是在 $\Sigma^{1/2 \top}k$ 处评估的各向同性高斯函数的傅里叶变换：
$\hat{h}(k) = \exp(-\frac{1}{2}(\Sigma^{1/2 \top}k)^\top(\Sigma^{1/2 \top}k))$
$\quad = \exp(-\frac{1}{2}k^\top\Sigma^{1/2}\Sigma^{1/2 \top}k)$
$\quad = \exp(-\frac{1}{2}k^\top\Sigma k) \quad \checkmark$
</details>

### 练习 3.4 (挑战)
证明对于带宽为 $K$ 的带限信号 $f$，使用高斯核最小化 $L^2$ 重建误差的最佳采样率为 $\Delta x = c\pi/K$，其中 $c \approx 0.8$。

**提示**：平衡混叠误差和核近似误差。

<details>
<summary>解决方案</summary>

总重建误差有两个组成部分：

1.  欠采样引起的混叠误差：
$E_{alias} = \int_{|k|>\pi/\Delta x} |\hat{f}(k)|^2|\hat{h}(k)|^2dk$

2.  核平滑引起的近似误差：
$E_{approx} = \int_{|k|<K} |\hat{f}(k)|^2|1-\hat{h}(k)|^2dk$

对于高斯核，其中 $\sigma = a\Delta x$：
$\hat{h}(k) = \exp(-k^2\sigma^2/2) = \exp(-k^2a^2\Delta x^2/2)$

总误差：
$E_{total} = \int_{|k|>\pi/\Delta x} |\hat{f}(k)|^2\exp(-k^2a^2\Delta x^2)dk + \int_{|k|<K} |\hat{f}(k)|^2(1-\exp(-k^2a^2\Delta x^2/2))^2dk$

为了获得最佳 $\Delta x$，$\partial E_{total}/\partial \Delta x = 0$。经过计算（使用 $|\hat{f}(k)|^2 \approx$ 在 $k=K$ 附近为常数）：

最佳间距满足：
$\pi/\Delta x \approx 0.8K$

因此 $c \approx 0.8$，得到 $\Delta x \approx 0.8\pi/K \quad \checkmark$
</details>

### 练习 3.5 (挑战)
设计一个能精确重现 $n$ 次多项式的核。对核有什么限制？

**提示**：使用矩条件 $\int x^\alpha h(x)dx = \delta_{\alpha,0}$ 对于 $|\alpha| \le n$。

<details>
<summary>解决方案</summary>

为了精确的多项式重现，我们需要：
$\sum_i p(p_i)h(x-p_i) = p(x)$ 对于所有次数 $\le n$ 的多项式 $p$

这要求核满足矩条件：
$\int x^\alpha h(x)dx = \delta_{\alpha,0}$ 对于所有多重索引 $|\alpha| \le n$

在一维中，这意味着：
- $\int h(x)dx = 1$ (单位分解)
- $\int x^j h(x)dx = 0$ 对于 $j = 1,2,...,n$

满足这些条件的最小支持核是 $n+1$ 次的 B 样条。

对于三次 B 样条 ($n=3$)：
$h(x) = \begin{cases}
  (2-|x|)^3/6,           & 1 \le |x| \le 2 \\
  1 - 3|x|^2/2 + 3|x|^3/4, & |x| < 1 \\
  0,                     & |x| > 2
\end{cases}$

这精确地重现了三次多项式。更高阶的核需要更宽的支持：对于 $n$ 次，支持宽度为 $n+2$。$\checkmark$
</details>

### 练习 3.6 (挑战)
推导用 $N$ 个点近似给定密度函数 $\rho(x)$ 以最小化渲染误差的最佳点分布。

**提示**：使用 Lloyd 算法的视角与 Voronoi 单元。

<details>
<summary>解决方案</summary>

对于密度 $\rho(x)$，点集 $\{p_i\}$ 和权重 $\{w_i\}$ 的 $L^2$ 近似误差：
$E = \int|\rho(x) - \sum_i w_i h(x-p_i)|^2dx$

对于固定核 $h$，最佳权重满足：
$\partial E/\partial w_j = 0 \implies \sum_i w_i \int h(x-p_i)h(x-p_j)dx = \int \rho(x)h(x-p_j)dx$

以矩阵形式表示：$Gw = b$，其中 $G_{ij} = \int h(x-p_i)h(x-p_j)dx$

对于固定权重下的最佳位置，$\partial E/\partial p_j = 0$ 给出：
$p_j = \frac{\int x\rho(x)h(x-p_j)dx}{\int \rho(x)h(x-p_j)dx}$

这是点 $j$ 的“影响区域”中 $\rho$ 的加权质心。

对于狄拉克核 ($h \to \delta$)，最佳分布的密度为：
$n(x) \propto \rho(x)^{d/(d+2)}$

在 3D 中：$n(x) \propto \rho(x)^{3/5}$

这使得在 $\rho$ 较大的地方有更多的样本，但不是线性比例。$\checkmark$
</details>

### 练习 3.7 (实现)
推导在给定 3D 各向异性高斯函数的情况下，屏幕空间中 EWA（椭圆加权平均）泼溅的方程。

**提示**：将 3D 协方差矩阵投影到 2D 屏幕空间。

<details>
<summary>解决方案</summary>

给定均值 $\mu$ 和协方差 $\Sigma_3$ 的 3D 高斯函数：
$g_3(x) = \exp(-\frac{1}{2}(x-\mu)^\top\Sigma_3^{-1}(x-\mu))$

在透视投影 $P$ 下，点 $x$ 映射到屏幕空间：
$s = P(x) = [x/z, y/z]^\top$

在 $x=\mu$ 处的雅可比矩阵 $J = \partial s/\partial x$：
$J = (1/z)[I_2 \quad -s]$ 其中 $I_2$ 是 2×2 单位矩阵

投影的 2D 协方差：
$\Sigma_2 = J\Sigma_3 J^\top$

对于与视图对齐的坐标，其中 $\Sigma_3 = \text{diag}(\sigma_x^2, \sigma_y^2, \sigma_u^2)$：
$\Sigma_2 = (1/z^2)\begin{bmatrix} \sigma_x^2 + \sigma_u^2 s_x^2 & \sigma_u^2 s_x s_y \\ \sigma_u^2 s_x s_y & \sigma_y^2 + \sigma_u^2 s_y^2 \end{bmatrix}$

2D 屏幕空间高斯函数：
$g_2(s) = \exp(-\frac{1}{2}(s-s_0)^\top\Sigma_2^{-1}(s-s_0))$

对于光栅化，计算包含 99% 高斯函数的椭圆：
$(s-s_0)^\top\Sigma_2^{-1}(s-s_0) < \chi^2_2(0.99) \approx 9.21$

这定义了边界框和每像素权重。$\checkmark$
</details>

### 练习 3.8 (开放式)
您将如何扩展统一体渲染方程以处理具有多次散射的参与介质？计算挑战是什么？

**提示**：考虑参与介质中的渲染方程和路径积分公式。

<details>
<summary>解决方案</summary>

具有多次散射的完整辐射传输方程：
$(\omega \cdot \nabla)L(x,\omega) = -\sigma_t(x)L(x,\omega) + \sigma_s(x)\int p(x,\omega',\omega)L(x,\omega')d\omega' + \sigma_a(x)L_e(x,\omega)$

这是一个积分微分方程。路径积分解：
$L(x,\omega) = \sum_{n=0}^\infty L^{(n)}(x,\omega)$

其中 $L^{(n)}$ 是 $n$ 次散射光：
$L^{(n+1)}(x,\omega) = \iint T(x,x')\sigma_s(x')p(x',\omega',\omega)L^{(n)}(x',\omega')dx'd\omega'$

计算挑战：
1.  **维度**：6D 位置-方向空间
2.  **递归**：每个散射阶都需要完整的解
3.  **各向异性相函数**：复杂的角度依赖性
4.  **异构介质**：空间变化的属性

实际近似：
- 对于高散射介质的扩散近似
- 用于角度依赖性的球谐函数
- 带有重要性采样的蒙特卡洛路径追踪
- 学习散射算子的神经网络

统一方程扩展为：
$L(o,\omega) = \int_0^\infty T(t)[\sum_{n=0}^\infty S^{(n)}(r(t),\omega)]dt$

其中 $S^{(n)}$ 表示 $n$ 次散射的入散射辐射。$\checkmark$
</details>

## 常见陷阱和错误 (Gotchas)

1.  **核归一化**：未能归一化核会导致能量损失/增益
    - 始终确保 $\int h(x)dx = 1$
    - 对于截断高斯函数，在支持范围内重新归一化

2.  **屏幕空间中的混叠**：3D 核在投影时可能导致严重的混叠
    - 使用 EWA 泼溅或在 3D 中预滤波
    - 绝不忽略各向异性核的 z 分量

3.  **数值精度**：透射率中的指数可能下溢
    - 使用对数空间计算：$\log T(t) = -\int \sigma(s)ds$
    - 在接近零时切换到 IEEE 754 双精度

4.  **排序顺序**：不正确的混合顺序会破坏透射率
    - 始终从前到后对泼溅进行排序以实现正确的遮挡
    - 从后到前仅适用于加性混合

5.  **边界处理**：边界附近的核需要特殊处理
    - 重新归一化以保持单位分解
    - 使用幽灵点或反射边界

6.  **性能瓶颈**：朴素实现扩展性差
    - 分层剔除对于大型点云至关重要
    - 可变泼溅大小导致的 GPU 分歧

## 最佳实践清单

### 设计评审
- [ ] 体密度 $\sigma(x)$ 是否正确归一化？
- [ ] 重建核是否满足单位分解？
- [ ] 采样率是否满足奈奎斯特准则？
- [ ] 误差度量是否适用于应用程序？

### 实现评审
- [ ] 透射率计算中的数值稳定性？
- [ ] alpha 混合的正确排序？
- [ ] 分层加速结构是否到位？
- [ ] 内存布局是否针对缓存一致性进行了优化？

### 验证
- [ ] 能量守恒是否已验证？
- [ ] 是否已测试随样本数量增加的收敛性？
- [ ] 是否与真实值/参考实现进行比较？
- [ ] 边缘情况（空空间、密集遮挡）是否已处理？

### 性能
- [ ] 是否实现复杂度 $O(N \log N)$ 或更好？
- [ ] 并行部分的 GPU 利用率是否 > 80%？
- [ ] 内存带宽是否不是瓶颈？
- [ ] 远处点的细节层次系统？

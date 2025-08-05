# 第12章：可微渲染 (Differentiable Rendering)

可微渲染是逆向渲染的计算基础，它使得我们能够通过梯度下降优化场景参数。本章深入探讨如何使渲染过程可微，特别是处理可见性不连续性带来的挑战。我们将从可微光线追踪的基础开始，逐步深入到边缘采样、软阴影等高级主题，并讨论实际实现中的偏差-方差权衡。

## 12.1 引言

渲染本质上是求解积分方程，但当我们需要对场景参数（几何、材质、光照）求导时，传统渲染算法面临根本性挑战。考虑像素强度的一般形式：

$$I(\boldsymbol{\theta}) = \int_{\Omega} f(\boldsymbol{x}, \boldsymbol{\theta}) \, d\boldsymbol{x}$$

其中 $\boldsymbol{\theta} \in \mathbb{R}^n$ 是场景参数，$\Omega$ 是积分域（如光线空间、路径空间），$f$ 是被积函数。可微渲染的核心问题是计算梯度 $\nabla_{\theta} I = \frac{\partial I}{\partial \boldsymbol{\theta}}$。

从测度论的角度，渲染积分可以更严格地表述为：

$$I(\boldsymbol{\theta}) = \int_{\mathcal{M}} L(\boldsymbol{x}, \boldsymbol{\omega}; \boldsymbol{\theta}) \, d\mu(\boldsymbol{x}, \boldsymbol{\omega})$$

其中 $\mathcal{M}$ 是相空间流形，$\mu$ 是其上的测度。当参数 $\boldsymbol{\theta}$ 变化时，不仅被积函数 $L$ 变化，测度 $\mu$ 和积分域 $\mathcal{M}$ 本身也可能变化，这导致了可微性的根本挑战。

从更一般的角度，渲染可以视为一个映射：

$$\mathcal{R}: \Theta \times \mathcal{S} \rightarrow \mathcal{I}$$

其中 $\Theta$ 是参数空间，$\mathcal{S}$ 是场景描述空间，$\mathcal{I}$ 是图像空间。可微渲染要求这个映射关于 $\theta$ 是可微的，即存在 Fréchet 导数：

$$\lim_{\|\delta\theta\| \to 0} \frac{\|\mathcal{R}(\theta + \delta\theta) - \mathcal{R}(\theta) - D\mathcal{R}(\theta)[\delta\theta]\|}{\|\delta\theta\|} = 0$$

**分层积分表述**

在实际渲染中，积分通常具有分层结构。例如，考虑直接光照和间接光照的分解：

$$I = I_{\text{direct}} + I_{\text{indirect}} = \int_{\mathcal{V}} \int_{\Omega} L_e V(\boldsymbol{x}, \boldsymbol{y}) G(\boldsymbol{x}, \boldsymbol{y}) \, d\omega_y d\boldsymbol{x} + \int_{\mathcal{P}_k} f_k(\bar{\boldsymbol{x}}) \, d\mu_k$$

其中 $\mathcal{V}$ 是可见表面，$\mathcal{P}_k$ 是长度为 $k$ 的路径空间。这种分层结构允许我们对不同类型的光传输应用不同的微分策略。

### 12.1.1 可微性的挑战

主要挑战来自三个层面：

1. **几何层面 - 可见性不连续**：当物体边缘移动时，积分域 $\Omega$ 随参数变化，导致分段连续性
2. **计算层面 - 高维积分**：路径追踪涉及高维积分（典型地 $d > 100$），直接微分计算代价为 $\mathcal{O}(d^2)$
3. **数值层面 - 稀有事件**：某些重要贡献（如焦散）概率密度极低，梯度估计方差高

对于参数依赖的积分域，广义 Leibniz 积分规则给出：

$$\frac{\partial}{\partial \theta} \int_{\Omega(\theta)} f(x, \theta) \, dx = \int_{\Omega(\theta)} \frac{\partial f}{\partial \theta} \, dx + \int_{\partial \Omega(\theta)} f \cdot v_n \, ds$$

其中 $v_n = \frac{\partial \boldsymbol{x}_{\text{boundary}}}{\partial \theta} \cdot \boldsymbol{n}$ 是边界的法向速度。第二项（边界项）在渲染中对应物体轮廓、阴影边界等几何不连续处的贡献。

**分布理论的严格处理**

更精确地，使用分布理论，可见性函数 $V$ 的导数包含 Dirac delta：

$$\frac{\partial V}{\partial \theta} = \sum_{i \in \text{silhouettes}} \delta(\boldsymbol{x} - \boldsymbol{x}_i(\theta)) \cdot v_{\perp,i}$$

这解释了为什么朴素的点采样方法会遗漏边缘梯度——采样到零测集的概率为零。

从 Sobolev 空间的角度，可见性函数 $V \in L^{\infty}(\Omega)$ 但 $V \notin W^{1,p}(\Omega)$ 对任意 $p \geq 1$，即 $V$ 不属于任何 Sobolev 空间。这导致标准的微分算子不适用，必须在弱意义或分布意义下理解导数。

**高维积分的维度诅咒**

对于长度为 $k$ 的光路，相空间维度为 $d = 3k$（忽略方向的约束）。梯度的 Hessian 矩阵规模为 $\mathcal{O}(d^2)$，导致：

1. **存储复杂度**：$\mathcal{O}(n^2 d^2)$，其中 $n$ 是参数数量
2. **计算复杂度**：每次 Hessian-向量乘积需要 $\mathcal{O}(nd)$ 操作
3. **数值条件数**：随维度指数增长，$\kappa(H) \sim e^{\alpha d}$

这促使我们寻找低秩近似或稀疏表示方法。

### 12.1.2 数学框架

可微渲染的数学基础建立在以下框架上：

**1. 路径积分表述**

渲染方程的路径空间形式：

$$I_j = \int_{\mathcal{P}} f_j(\bar{x}) d\mu(\bar{x})$$

其中 $\bar{x} = (\boldsymbol{x}_0, \boldsymbol{x}_1, ..., \boldsymbol{x}_k)$ 是长度为 $k$ 的路径，$\mu$ 是路径测度。梯度计算需要：

$$\nabla_\theta I_j = \int_{\mathcal{P}} \nabla_\theta f_j(\bar{x}) d\mu(\bar{x}) + \int_{\mathcal{P}} f_j(\bar{x}) \nabla_\theta \log p(\bar{x}|\theta) d\mu(\bar{x})$$

第二项来自测度的参数依赖性（重参数化梯度）。

路径测度的具体形式为：
$$d\mu(\bar{x}) = dA(\boldsymbol{x}_0) \prod_{i=1}^{k-1} G(\boldsymbol{x}_i, \boldsymbol{x}_{i+1}) dA(\boldsymbol{x}_i)$$

其中 $G(\boldsymbol{x}, \boldsymbol{y}) = \frac{V(\boldsymbol{x}, \boldsymbol{y})|\cos\theta_x||\cos\theta_y|}{|\boldsymbol{x} - \boldsymbol{y}|^2}$ 是几何项。当几何参数变化时，$G$ 和 $V$ 都会变化，导致复杂的梯度表达式。

**2. 伴随方法**

对于复杂系统，伴随方法提供高效的梯度计算。定义拉格朗日量：

$$\mathcal{L} = I(u, \theta) + \langle \lambda, G(u, \theta) \rangle$$

其中 $G(u, \theta) = 0$ 是约束（如渲染方程），$\lambda$ 是伴随变量。最优性条件给出：

$$\frac{dI}{d\theta} = -\left\langle \lambda, \frac{\partial G}{\partial \theta} \right\rangle$$

其中 $\lambda$ 满足伴随方程：$(\partial G/\partial u)^T \lambda = -\partial I/\partial u$。

在渲染中，$u = L(\boldsymbol{x}, \boldsymbol{\omega})$ 是辐射度场，约束 $G$ 是渲染方程：
$$G(L, \theta) = L - L_e - \mathcal{T}[L] = 0$$

其中 $\mathcal{T}$ 是传输算子。伴随方程变为：
$$(\mathcal{I} - \mathcal{T}^*)[\lambda] = -\frac{\partial I}{\partial L}$$

这里 $\mathcal{T}^*$ 是 $\mathcal{T}$ 的伴随算子，在物理上对应重要性传输。

**3. 变分原理**

渲染可以表述为变分问题：

$$I = \min_L \mathcal{F}[L]$$

其中 $\mathcal{F}$ 是某个泛函。参数变化引起的一阶变分：

$$\delta I = \int \frac{\delta \mathcal{F}}{\delta L} \delta L \, d\Omega$$

这提供了另一种计算梯度的视角。

具体地，对于路径追踪，泛函可以取为：
$$\mathcal{F}[L] = \frac{1}{2}\|L - L_e - \mathcal{T}[L]\|^2_{L^2(\mathcal{M})}$$

这导致梯度流方程：
$$\frac{\partial L}{\partial t} = -\frac{\delta \mathcal{F}}{\delta L} = (\mathcal{I} - \mathcal{T})^*(\mathcal{I} - \mathcal{T})[L] - (\mathcal{I} - \mathcal{T})^*[L_e]$$

### 12.1.3 可微渲染的应用

可微渲染在以下领域有重要应用：

- **3D重建**：从2D图像恢复3D几何和材质，解决逆问题 $\theta^* = \arg\min_\theta \|I_{\text{obs}} - I(\theta)\|^2$
- **场景理解**：推断光照、材质分解，通常表述为贝叶斯推断 $p(\theta|I) \propto p(I|\theta)p(\theta)$
- **内容创作**：通过优化生成满足约束的3D内容，如 $\min_\theta \mathcal{L}_{\text{style}}(I(\theta)) + \lambda \mathcal{R}(\theta)$
- **机器人视觉**：主动感知和场景操作，需要实时梯度计算
- **计算成像**：联合优化光学系统和计算管线

**逆问题的正则化**

逆渲染本质上是病态的（ill-posed），因为：
1. **非唯一性**：多个场景配置可能产生相同图像
2. **不稳定性**：小的观测噪声可能导致解的巨大变化
3. **不存在性**：某些图像可能无法由物理有效的场景产生

为此，我们需要正则化。常用方法包括：

- **Tikhonov 正则化**：$\mathcal{L}_{\text{reg}} = \|I_{\text{obs}} - I(\theta)\|^2 + \lambda\|\theta - \theta_0\|^2$
- **稀疏性约束**：$\mathcal{L}_{\text{sparse}} = \|I_{\text{obs}} - I(\theta)\|^2 + \lambda\|\theta\|_1$
- **流形约束**：限制 $\theta \in \mathcal{M}_{\text{valid}}$，其中 $\mathcal{M}_{\text{valid}}$ 是物理有效参数的流形

**优化景观的特性**

可微渲染的损失函数通常具有复杂的优化景观：

1. **多模态性**：存在多个局部最小值，对应不同的场景解释
2. **非凸性**：由于可见性变化和材质-光照歧义
3. **梯度稀疏性**：大部分参数对特定像素的梯度为零
4. **尺度敏感性**：不同参数类型（几何、材质、光照）的梯度尺度相差巨大

这些特性要求专门的优化策略，如自适应学习率、预条件器和多尺度方法。

## 12.2 可微光线追踪

可微光线追踪是可微渲染的基础。我们需要计算光线-物体交点及相关量对场景参数的导数。这涉及到隐函数定理、微分几何和数值稳定性等核心数学工具。

### 12.2.1 光线-表面交点的导数

考虑参数化光线 $\boldsymbol{r}(t) = \boldsymbol{o} + t\boldsymbol{d}$，其中 $\boldsymbol{o} \in \mathbb{R}^3$ 是原点，$\boldsymbol{d} \in \mathbb{S}^2$ 是单位方向。对于隐式表面 $F: \mathbb{R}^3 \times \mathbb{R}^n \rightarrow \mathbb{R}$，其中 $F(\boldsymbol{x}, \boldsymbol{\theta}) = 0$ 定义了表面，交点满足：

$$F(\boldsymbol{o} + t^*\boldsymbol{d}, \boldsymbol{\theta}) = 0$$

**隐函数定理的应用**

假设 $F$ 是 $C^1$ 级且 $\boldsymbol{d} \cdot \nabla_x F \neq 0$（光线不与表面相切），则存在隐函数 $t^*(\boldsymbol{\theta})$ 使得：

$$G(\boldsymbol{\theta}, t) \equiv F(\boldsymbol{o} + t\boldsymbol{d}, \boldsymbol{\theta}) = 0$$

对 $\boldsymbol{\theta}$ 求全微分：

$$\frac{\partial G}{\partial \boldsymbol{\theta}} + \frac{\partial G}{\partial t}\frac{\partial t^*}{\partial \boldsymbol{\theta}} = 0$$

解得：

$$\frac{\partial t^*}{\partial \boldsymbol{\theta}} = -\frac{\partial F/\partial \boldsymbol{\theta}}{\boldsymbol{d} \cdot \nabla_x F} = -\frac{\partial F/\partial \boldsymbol{\theta}}{\partial F/\partial t}$$

**几何解释与退化情况**

从几何角度，分母 $\boldsymbol{d} \cdot \nabla_x F$ 表示光线方向与表面法线的点积（缩放后）。当此值接近零时：

1. **掠射情况**：光线几乎平行于表面，$|\boldsymbol{d} \cdot \nabla_x F| \ll 1$
2. **曲率奇异**：在高曲率区域，法线变化剧烈
3. **参数奇异**：某些参数化在特定点退化（如球坐标的极点）

为处理退化，我们引入正则化：

$$\frac{\partial t^*}{\partial \boldsymbol{\theta}} = -\frac{\partial F/\partial \boldsymbol{\theta}}{\boldsymbol{d} \cdot \nabla_x F + \epsilon \|\nabla_x F\|}$$

其中 $\epsilon \sim 10^{-6}$ 提供数值稳定性，同时保持物理意义。

**交点位置的全导数**

交点 $\boldsymbol{x}^* = \boldsymbol{o} + t^*\boldsymbol{d}$ 的导数使用链式法则：

$$\frac{\partial \boldsymbol{x}^*}{\partial \boldsymbol{\theta}} = \frac{\partial t^*}{\partial \boldsymbol{\theta}} \boldsymbol{d} + \frac{\partial \boldsymbol{o}}{\partial \boldsymbol{\theta}} + t^* \frac{\partial \boldsymbol{d}}{\partial \boldsymbol{\theta}}$$

对于不同的参数化方式：
- 若 $\boldsymbol{\theta}$ 仅影响表面：$\frac{\partial \boldsymbol{o}}{\partial \boldsymbol{\theta}} = \frac{\partial \boldsymbol{d}}{\partial \boldsymbol{\theta}} = 0$
- 若 $\boldsymbol{\theta}$ 包含相机参数：需要考虑光线的变化

**高阶导数与曲率效应**

二阶导数揭示了曲率对交点的影响：

$$\frac{\partial^2 t^*}{\partial \boldsymbol{\theta}^2} = -\frac{1}{\boldsymbol{d} \cdot \nabla F}\left[\frac{\partial^2 F}{\partial \boldsymbol{\theta}^2} - \frac{\partial t^*}{\partial \boldsymbol{\theta}} \otimes \frac{\partial}{\partial \boldsymbol{\theta}}(\boldsymbol{d} \cdot \nabla F)\right]$$

这在优化中用于：
1. **Newton 方法**：需要 Hessian 信息
2. **曲率感知步长**：在高曲率区域减小步长
3. **置信区域方法**：二阶近似的有效范围

**数值稳定性考虑**

当 $\boldsymbol{d} \cdot \nabla F \approx 0$（掠射角度）时，导数计算不稳定。实际实现中使用：

$$\frac{\partial t^*}{\partial \boldsymbol{\theta}} = -\frac{\partial F/\partial \boldsymbol{\theta}}{\max(\boldsymbol{d} \cdot \nabla F, \epsilon)}$$

其中 $\epsilon \sim 10^{-6}$ 避免除零。

更精细的处理使用自适应正则化：
$$\epsilon_{\text{adaptive}} = \epsilon_0 \cdot \max(1, \|\nabla^2 F\| \cdot \Delta t)$$

其中 $\Delta t$ 是光线步长，$\|\nabla^2 F\|$ 估计局部曲率。

### 12.2.2 表面法线的微分

表面法线 $\boldsymbol{n} = \nabla F / |\nabla F|$ 的导数涉及归一化操作的微分。

**直接计算方法**

设 $\boldsymbol{g} = \nabla_x F$，则 $\boldsymbol{n} = \boldsymbol{g}/|\boldsymbol{g}|$。使用链式法则：

$$\frac{\partial \boldsymbol{n}}{\partial \boldsymbol{\theta}} = \frac{\partial}{\partial \boldsymbol{\theta}}\left(\frac{\boldsymbol{g}}{|\boldsymbol{g}|}\right) = \frac{1}{|\boldsymbol{g}|}\frac{\partial \boldsymbol{g}}{\partial \boldsymbol{\theta}} - \frac{\boldsymbol{g}}{|\boldsymbol{g}|^3}\boldsymbol{g}^T\frac{\partial \boldsymbol{g}}{\partial \boldsymbol{\theta}}$$

整理后：

$$\frac{\partial \boldsymbol{n}}{\partial \boldsymbol{\theta}} = \frac{1}{|\nabla F|}\left(\boldsymbol{I} - \boldsymbol{n} \otimes \boldsymbol{n}\right)\frac{\partial \nabla F}{\partial \boldsymbol{\theta}}$$

其中 $(\boldsymbol{I} - \boldsymbol{n} \otimes \boldsymbol{n})$ 是到切平面的投影算子。

**几何意义与 Weingarten 映射**

法线导数与微分几何中的 Weingarten 映射（形状算子）密切相关。对于参数曲面 $\boldsymbol{x}(u,v)$，Weingarten 映射定义为：

$$\mathcal{W} = -d\boldsymbol{n} = \begin{bmatrix} \frac{\partial \boldsymbol{n}}{\partial u} & \frac{\partial \boldsymbol{n}}{\partial v} \end{bmatrix}$$

主曲率 $\kappa_1, \kappa_2$ 是 $\mathcal{W}$ 的特征值。当参数 $\boldsymbol{\theta}$ 改变时，曲率的变化为：

$$\frac{\partial \kappa_i}{\partial \boldsymbol{\theta}} = \boldsymbol{e}_i^T \frac{\partial \mathcal{W}}{\partial \boldsymbol{\theta}} \boldsymbol{e}_i$$

其中 $\boldsymbol{e}_i$ 是主方向。这解释了为什么高曲率区域的法线对参数变化更敏感。

**Hessian 矩阵的计算**

混合偏导数 $\frac{\partial \nabla_x F}{\partial \boldsymbol{\theta}}$ 是一个 $3 \times n$ 矩阵：

$$\left[\frac{\partial \nabla_x F}{\partial \boldsymbol{\theta}}\right]_{ij} = \frac{\partial^2 F}{\partial x_i \partial \theta_j}$$

对于特定几何形状：
- **平面**: Hessian 为零，法线导数简化
- **球面**: $\nabla F = 2(\boldsymbol{x} - \boldsymbol{c})$，$\frac{\partial \nabla F}{\partial \boldsymbol{c}} = -2\boldsymbol{I}$
- **二次曲面**: Hessian 为常数矩阵

**几何意义**

法线导数描述了表面弯曲如何随参数变化。其在切平面内的分量与曲率变化相关：

$$\kappa_{\text{change}} = \boldsymbol{t}^T \frac{\partial \boldsymbol{n}}{\partial \boldsymbol{\theta}}$$

其中 $\boldsymbol{t}$ 是切向量，$\kappa_{\text{change}}$ 表示沿 $\boldsymbol{t}$ 方向的曲率变化。

### 12.2.3 参数化表面的处理

对于参数化表面 $\boldsymbol{x}: \mathbb{R}^2 \times \mathbb{R}^n \rightarrow \mathbb{R}^3$，即 $\boldsymbol{x}(u, v, \boldsymbol{\theta})$，交点满足：

$$\boldsymbol{x}(u^*, v^*, \boldsymbol{\theta}) = \boldsymbol{o} + t^*\boldsymbol{d}$$

这给出三个约束方程。

**约束系统的微分**

定义约束函数 $\boldsymbol{F}: \mathbb{R}^3 \times \mathbb{R}^n \rightarrow \mathbb{R}^3$：

$$\boldsymbol{F}(u, v, t, \boldsymbol{\theta}) = \boldsymbol{x}(u, v, \boldsymbol{\theta}) - \boldsymbol{o} - t\boldsymbol{d} = \boldsymbol{0}$$

对 $\boldsymbol{\theta}$ 求全微分：

$$\frac{\partial \boldsymbol{F}}{\partial (u, v, t)} \begin{bmatrix} \frac{\partial u^*}{\partial \boldsymbol{\theta}} \\ \frac{\partial v^*}{\partial \boldsymbol{\theta}} \\ \frac{\partial t^*}{\partial \boldsymbol{\theta}} \end{bmatrix} + \frac{\partial \boldsymbol{F}}{\partial \boldsymbol{\theta}} = \boldsymbol{0}$$

**Jacobian 矩阵**

Jacobian 矩阵为：

$$\boldsymbol{J} = \frac{\partial \boldsymbol{F}}{\partial (u, v, t)} = \begin{bmatrix} \frac{\partial \boldsymbol{x}}{\partial u} & \frac{\partial \boldsymbol{x}}{\partial v} & -\boldsymbol{d} \end{bmatrix}$$

这是一个 $3 \times 3$ 矩阵。假设 $\det(\boldsymbol{J}) \neq 0$（非退化情况），则：

$$\begin{bmatrix} \frac{\partial u^*}{\partial \boldsymbol{\theta}} \\ \frac{\partial v^*}{\partial \boldsymbol{\theta}} \\ \frac{\partial t^*}{\partial \boldsymbol{\theta}} \end{bmatrix} = -\boldsymbol{J}^{-1} \frac{\partial \boldsymbol{x}}{\partial \boldsymbol{\theta}}$$

**特殊情况处理**

1. **退化情况**: 当 $\det(\boldsymbol{J}) = 0$ 时，通常发生在：
   - 光线与表面相切
   - 参数化奇异点（如球面两极）
   
2. **过参数化**: 对于过参数化表面（如 NURBS），使用最小二乘：
   $$\Delta \boldsymbol{p} = (\boldsymbol{J}^T\boldsymbol{J})^{-1}\boldsymbol{J}^T \Delta \boldsymbol{r}$$

3. **约束优化**: 对于复杂参数化，可能需要迭代求解：
   $$\boldsymbol{p}_{k+1} = \boldsymbol{p}_k - \alpha \boldsymbol{J}_k^{-1} \boldsymbol{F}_k$$

### 12.2.4 反射和折射方向的导数

**反射方向的微分**

反射定律：$\boldsymbol{r} = \boldsymbol{i} - 2(\boldsymbol{i} \cdot \boldsymbol{n})\boldsymbol{n}$，其中 $\boldsymbol{i}$ 是入射方向，$\boldsymbol{n}$ 是法线。

使用乘积法则和链式法则：

$$\frac{\partial \boldsymbol{r}}{\partial \boldsymbol{\theta}} = \frac{\partial \boldsymbol{i}}{\partial \boldsymbol{\theta}} - 2\frac{\partial}{\partial \boldsymbol{\theta}}[(\boldsymbol{i} \cdot \boldsymbol{n})\boldsymbol{n}]$$

展开第二项：

$$\frac{\partial}{\partial \boldsymbol{\theta}}[(\boldsymbol{i} \cdot \boldsymbol{n})\boldsymbol{n}] = (\boldsymbol{i} \cdot \boldsymbol{n})\frac{\partial \boldsymbol{n}}{\partial \boldsymbol{\theta}} + \boldsymbol{n} \otimes \frac{\partial}{\partial \boldsymbol{\theta}}(\boldsymbol{i} \cdot \boldsymbol{n})$$

其中：
$$\frac{\partial}{\partial \boldsymbol{\theta}}(\boldsymbol{i} \cdot \boldsymbol{n}) = \frac{\partial \boldsymbol{i}}{\partial \boldsymbol{\theta}} \cdot \boldsymbol{n} + \boldsymbol{i} \cdot \frac{\partial \boldsymbol{n}}{\partial \boldsymbol{\theta}}$$

**折射方向的微分**

Snell 定律的向量形式：

$$\boldsymbol{t} = \frac{n_1}{n_2}\boldsymbol{i} + \left(\frac{n_1}{n_2}\cos\theta_i - \cos\theta_t\right)\boldsymbol{n}$$

其中 $\cos\theta_t = \sqrt{1 - \left(\frac{n_1}{n_2}\right)^2(1 - \cos^2\theta_i)}$。

对 $\boldsymbol{\theta}$ 求导：

$$\frac{\partial \boldsymbol{t}}{\partial \boldsymbol{\theta}} = \frac{n_1}{n_2}\frac{\partial \boldsymbol{i}}{\partial \boldsymbol{\theta}} + \frac{\partial}{\partial \boldsymbol{\theta}}\left[\left(\frac{n_1}{n_2}\cos\theta_i - \cos\theta_t\right)\boldsymbol{n}\right]$$

关键是计算 $\frac{\partial \cos\theta_t}{\partial \boldsymbol{\theta}}$：

$$\frac{\partial \cos\theta_t}{\partial \boldsymbol{\theta}} = \frac{(n_1/n_2)^2\sin\theta_i}{\cos\theta_t}\frac{\partial \cos\theta_i}{\partial \boldsymbol{\theta}}$$

**全内反射的处理**

当 $(n_1/n_2)\sin\theta_i > 1$ 时发生全内反射。在临界角附近，导数变得不稳定。实际实现中：

1. **平滑过渡**: 使用 sigmoid 函数平滑切换
2. **分支处理**: 分别计算折射和反射分支的梯度

## 12.3 边缘采样与重参数化

可见性不连续是可微渲染的核心挑战。当物体移动时，其轮廓边缘导致积分域的变化，产生 Dirac delta 函数形式的梯度。本节深入探讨如何数学上严格处理这些不连续性。

### 12.3.1 边缘积分理论

**Reynolds 传输定理的应用**

考虑依赖于参数 $\theta \in \mathbb{R}$ 的时变积分域 $\Omega(\theta) \subset \mathbb{R}^d$：

$$I(\theta) = \int_{\Omega(\theta)} f(\boldsymbol{x}, \theta) \, d\boldsymbol{x}$$

Reynolds 传输定理（也称为 Leibniz-Reynolds 定理）给出：

$$\frac{dI}{d\theta} = \int_{\Omega(\theta)} \frac{\partial f}{\partial \theta} \, d\boldsymbol{x} + \oint_{\partial\Omega(\theta)} f \, v_n \, ds - \int_{\Omega(\theta)} f \, \nabla \cdot \boldsymbol{v} \, d\boldsymbol{x}$$

其中：
- $v_n = \boldsymbol{v} \cdot \boldsymbol{n}$ 是边界的法向速度
- $\boldsymbol{v} = \frac{\partial \boldsymbol{x}}{\partial \theta}$ 是速度场
- $\boldsymbol{n}$ 是边界的外法向

在渲染中，通常假设被积函数不随位置移动（$\nabla \cdot \boldsymbol{v} = 0$），简化为：

$$\frac{dI}{d\theta} = \int_{\Omega(\theta)} \frac{\partial f}{\partial \theta} \, d\boldsymbol{x} + \oint_{\partial\Omega(\theta)} f \, v_n \, ds$$

**渲染中的具体形式**

对于像素积分 $I_p = \int_{\mathcal{A}_p} L(\boldsymbol{x}, \boldsymbol{\omega}) \, dA$，其中 $\mathcal{A}_p$ 是像素覆盖区域：

$$\frac{\partial I_p}{\partial \theta} = \int_{\mathcal{A}_p} \frac{\partial L}{\partial \theta} \, dA + \sum_{e \in \text{edges}} \int_{e \cap \mathcal{A}_p} L \cdot v_{\perp} \, dl$$

其中 $v_{\perp}$ 是边缘在图像平面上的垂直速度。

### 12.3.2 可见性函数的导数

**分布理论视角**

可见性函数 $V: \mathbb{R}^3 \times \mathbb{R}^3 \rightarrow \{0, 1\}$ 定义为：

$$V(\boldsymbol{x}, \boldsymbol{y}) = \begin{cases}
1 & \text{如果射线 } \overline{\boldsymbol{xy}} \text{ 无遮挡} \\
0 & \text{否则}
\end{cases}$$

在分布意义下，其导数为：

$$\frac{\partial V}{\partial \theta} = \sum_{i \in \mathcal{S}} \delta_{\Gamma_i} \cdot n_i(\theta)$$

其中：
- $\mathcal{S}$ 是所有轮廓边缘的集合
- $\Gamma_i$ 是第 $i$ 条边缘曲线
- $\delta_{\Gamma_i}$ 是沿 $\Gamma_i$ 的 Dirac 测度
- $n_i(\theta)$ 是法向速度函数

**几何解释**

对于光线 $\boldsymbol{r}(t) = \boldsymbol{x} + t(\boldsymbol{y} - \boldsymbol{x})$ 和遮挡物边缘 $\Gamma(\theta)$，定义：

- **签名距离**: $d(\boldsymbol{r}, \Gamma) = \min_{\boldsymbol{p} \in \Gamma} \|\boldsymbol{r} - \boldsymbol{p}\|$
- **投影速度**: $v_{\perp} = \boldsymbol{v}_{\Gamma} \cdot \boldsymbol{n}_{\perp}$

其中 $\boldsymbol{v}_{\Gamma} = \frac{\partial \Gamma}{\partial \theta}$ 是边缘速度，$\boldsymbol{n}_{\perp}$ 是垂直于光线的单位向量。

则：
$$\frac{\partial V}{\partial \theta} = \delta(d) \cdot v_{\perp} \cdot \|\boldsymbol{y} - \boldsymbol{x}\|^{-1}$$

最后一项来自射线参数化的 Jacobian。

### 12.3.3 边缘采样策略

**轮廓边缘的检测**

对于三角网格，轮廓边缘（silhouette edge）的数学定义为：

$$\mathcal{E}_{\text{sil}} = \{e \in \mathcal{E} : (\boldsymbol{n}_1 \cdot \boldsymbol{v})(\boldsymbol{n}_2 \cdot \boldsymbol{v}) < 0\}$$

其中：
- $\mathcal{E}$ 是所有边的集合
- $\boldsymbol{n}_1, \boldsymbol{n}_2$ 是共享边 $e$ 的两个三角面的法线
- $\boldsymbol{v}$ 是视线方向

对于光滑曲面，轮廓满足：
$$\boldsymbol{n}(\boldsymbol{x}) \cdot \boldsymbol{v}(\boldsymbol{x}) = 0$$

**边缘的参数化**

边缘可以参数化为曲线 $\boldsymbol{\gamma}: [0, 1] \rightarrow \mathbb{R}^3$。在图像平面上的投影为：

$$\boldsymbol{p}(t) = \Pi(\boldsymbol{\gamma}(t))$$

其中 $\Pi$ 是投影算子。边缘积分变为：

$$\int_{\text{edge}} f \, ds = \int_0^1 f(\boldsymbol{p}(t)) \left\|\frac{d\boldsymbol{p}}{dt}\right\| dt$$

**重要性采样**

为了高效采样，我们根据边缘对像素的贡献分配样本。定义重要性度量：

$$w(t) = \left\|\frac{d\boldsymbol{p}}{dt}\right\| \cdot |f(\boldsymbol{p}(t))| \cdot K(\boldsymbol{p}(t) - \boldsymbol{p}_{\text{pixel}})$$

其中 $K$ 是像素滤波核。采样密度正比于 $w(t)$。

**分层采样策略**

实际实现中使用分层方法：

1. **粗级选择**: 使用层次包围盒快速筛选潜在边缘
2. **精确检测**: 对候选边缘计算精确的轮廓条件
3. **自适应细分**: 根据曲率和投影长度细分边缘

### 12.3.4 重参数化技巧

**重参数化的数学基础**

重参数化将不连续的采样决策转换为连续参数。考虑离散选择：

$$I = \sum_{i} \mathbb{1}[\text{condition}_i] \cdot f_i$$

重参数化为：

$$I(\theta) = \sum_{i} w_i(\theta) \cdot f_i(\theta)$$

其中 $w_i$ 是连续权重函数。

**边缘积分的重参数化**

对于边缘积分，我们将变化的积分域转换为固定域：

$$\int_{\Gamma(\theta)} f \, ds = \int_0^1 f(\boldsymbol{\gamma}(t, \theta)) \left\|\frac{\partial \boldsymbol{\gamma}}{\partial t}\right\| dt$$

关键是找到合适的参数化 $\boldsymbol{\gamma}(t, \theta)$ 使得：
1. $t \in [0,1]$ 与 $\theta$ 无关
2. $\boldsymbol{\gamma}$ 关于 $\theta$ 可微

**具体重参数化方法**

1. **投影重参数化**
   
   将 3D 边缘投影到 2D 图像平面：
   $$\boldsymbol{p}(t, \theta) = \Pi(\boldsymbol{\gamma}(t, \theta))$$
   
   边缘积分变为：
   $$\int_{\text{edge}} L \, ds = \int_0^1 L(\boldsymbol{p}(t)) J(t, \theta) dt$$
   
   其中 $J$ 是 Jacobian 行列式。

2. **解析积分**
   
   对于特定情况，可以解析计算。例如，三角形遮挡的线性边缘：
   
   $$\int_{\text{linear edge}} f_0 + \boldsymbol{g} \cdot \boldsymbol{x} \, ds = f_0 L + \frac{L^2}{2} \boldsymbol{g} \cdot \boldsymbol{t}$$
   
   其中 $L$ 是边缘长度，$\boldsymbol{t}$ 是单位切向量。

3. **基于采样的方法**
   
   使用重要性采样和控制变量：
   
   $$\nabla_\theta I \approx \frac{1}{N} \sum_{i=1}^N \frac{f(\boldsymbol{x}_i)}{p(\boldsymbol{x}_i)} \nabla_\theta \log p(\boldsymbol{x}_i|\theta)$$
   
   这是 REINFORCE 算法在渲染中的应用。

## 12.4 可微阴影与可见性

阴影是可见性的特殊情况，涉及光源、遮挡物和接收表面三者的关系。可微阴影计算需要正确处理软阴影的连续变化。本节探讨如何使阴影计算在数学上可微且数值稳定。

### 12.4.1 硬阴影的可微化

硬阴影的可见性函数是二值的：

$$V_{\text{hard}}(\boldsymbol{x}, \boldsymbol{l}) = \begin{cases}
1 & \text{如果光线 } \boldsymbol{x} \to \boldsymbol{l} \text{ 无遮挡} \\
0 & \text{否则}
\end{cases}$$

其梯度在阴影边界处包含delta函数。为了可微化，我们使用软化近似：

$$V_{\text{soft}}(\boldsymbol{x}, \boldsymbol{l}) = \sigma\left(\frac{d_{\text{signed}}}{\epsilon}\right)$$

其中 $\sigma$ 是sigmoid函数，$d_{\text{signed}}$ 是到阴影边界的带符号距离，$\epsilon$ 控制软化程度。

### 12.4.2 面光源的软阴影

对于面光源 $A_L$，阴影计算涉及光源上的积分：

$$L_{\text{shadow}}(\boldsymbol{x}) = \int_{A_L} L_e(\boldsymbol{y}) V(\boldsymbol{x}, \boldsymbol{y}) G(\boldsymbol{x}, \boldsymbol{y}) \, dA_{\boldsymbol{y}}$$

其中 $G$ 是几何项。梯度计算需要考虑：

1. **光源边界**：当遮挡物移动时光源可见部分的变化
2. **遮挡物边界**：影响可见性函数的轮廓边缘

总梯度为：

$$\frac{\partial L_{\text{shadow}}}{\partial \theta} = \int_{A_L} L_e \frac{\partial V}{\partial \theta} G \, dA + \oint_{\partial A_V} L_e G v_n \, ds$$

其中 $A_V$ 是光源的可见部分，第二项是边缘贡献。

### 12.4.3 球面光源的高效采样

对于球面光源，我们可以解析计算某些边缘积分。考虑半径为 $r$ 的球面光源，从点 $\boldsymbol{x}$ 观察时的立体角为：

$$\Omega = 2\pi\left(1 - \cos\alpha\right)$$

其中 $\alpha$ 是半锥角。当球被部分遮挡时，可见立体角的梯度为：

$$\frac{\partial \Omega_{\text{vis}}}{\partial \theta} = \oint_{\mathcal{C}} \sin\phi \, v_{\phi} \, d\phi$$

其中 $\mathcal{C}$ 是遮挡轮廓在球面上的投影。

### 12.4.4 多重遮挡处理

当存在多个遮挡物时，可见性函数变为：

$$V(\boldsymbol{x}, \boldsymbol{y}) = \prod_{i=1}^N V_i(\boldsymbol{x}, \boldsymbol{y})$$

使用对数空间避免数值问题：

$$\log V = \sum_{i=1}^N \log V_i$$

梯度计算使用链式法则：

$$\frac{\partial V}{\partial \theta} = V \sum_{i=1}^N \frac{1}{V_i} \frac{\partial V_i}{\partial \theta}$$

这种形式在 $V_i \approx 0$ 时数值稳定。

## 12.5 梯度偏差与方差权衡

可微渲染中的梯度估计器在偏差和方差之间存在基本权衡。理解这种权衡对于设计高效的优化算法至关重要。

### 12.5.1 无偏梯度估计器

理想的梯度估计器满足：

$$\mathbb{E}[\nabla_{\theta} \hat{I}] = \nabla_{\theta} I$$

对于连续被积函数，标准蒙特卡洛估计器是无偏的：

$$\nabla_{\theta} I \approx \frac{1}{N} \sum_{i=1}^N \nabla_{\theta} f(\boldsymbol{x}_i, \theta)$$

但对于包含可见性的积分，边缘项使得无偏估计器的方差极高：

$$\text{Var}[\nabla_{\theta} \hat{I}] = \mathcal{O}(1/\sqrt{N})$$

收敛速度慢，因为需要采样到罕见的边缘事件。

### 12.5.2 有偏低方差估计器

通过引入系统性偏差，我们可以显著降低方差。常见方法包括：

1. **软化可见性**：
$$V_{\epsilon}(\boldsymbol{x}, \boldsymbol{y}) = \sigma(d/\epsilon)$$
偏差为 $\mathcal{O}(\epsilon)$，方差降低到 $\mathcal{O}(\epsilon^2)$

2. **有限差分**：
$$\nabla_{\theta} I \approx \frac{I(\theta + h) - I(\theta - h)}{2h}$$
偏差为 $\mathcal{O}(h^2)$，但避免了边缘采样

3. **平滑聚合**：
$$I_{\text{smooth}} = \int I(\theta') K(\theta - \theta') d\theta'$$
其中 $K$ 是平滑核

### 12.5.3 控制变量法

控制变量通过引入相关但易计算的辅助函数减少方差：

$$\nabla_{\theta} \hat{I}_{\text{CV}} = \nabla_{\theta} \hat{I} - \beta(\nabla_{\theta} \hat{C} - \nabla_{\theta} C)$$

其中 $C$ 是控制变量，$\beta$ 是最优系数：

$$\beta^* = \frac{\text{Cov}[\nabla_{\theta} \hat{I}, \nabla_{\theta} \hat{C}]}{\text{Var}[\nabla_{\theta} \hat{C}]}$$

常用控制变量包括：
- 简化场景的解析解
- 低分辨率渲染
- 预计算的辐照度缓存

### 12.5.4 自适应采样策略

根据梯度估计的不确定性动态分配样本：

1. **重要性度量**：
$$w(\boldsymbol{x}) = |\nabla_{\theta} f(\boldsymbol{x})| \cdot \sigma_{\text{local}}(\boldsymbol{x})$$

2. **样本分配**：
$$N_i \propto \sqrt{w_i \sum_j w_j}$$

3. **早停准则**：
$$\frac{\sigma[\nabla_{\theta} \hat{I}]}{|\nabla_{\theta} \hat{I}|} < \tau$$

### 12.5.5 混合估计器

结合多种方法的优点：

$$\nabla_{\theta} I_{\text{hybrid}} = \alpha \nabla_{\theta} I_{\text{edge}} + (1-\alpha) \nabla_{\theta} I_{\text{interior}}$$

其中：
- $I_{\text{edge}}$：使用边缘采样的无偏估计
- $I_{\text{interior}}$：忽略边缘的有偏估计
- $\alpha$：基于置信度的混合权重

## 12.6 实现框架

现代可微渲染框架提供了不同层次的抽象和优化。我们介绍主要框架的设计理念和使用场景。

### 12.6.1 Mitsuba 3

Mitsuba 3 是基于 Dr.Jit 的可微渲染器，支持多种后端（CPU、CUDA、LLVM）：

**核心特性**：
- **自动微分**：前向和反向模式
- **边缘采样**：内置轮廓积分
- **路径重放**：高效梯度计算

**梯度计算模式**：
1. **前向模式**：计算 $\frac{\partial I}{\partial \theta_i}$ 对单个参数
2. **反向模式**：计算 $\nabla_\theta L$ 对所有参数
3. **路径空间微分**：处理间接光照

**典型工作流**：参数化场景 → 渲染 → 计算损失 → 反向传播

### 12.6.2 PyTorch3D 与神经渲染

PyTorch3D 将可微渲染集成到深度学习框架：

**光栅化管线**：
- **可微光栅化**：软边缘近似
- **纹理采样**：双线性插值的梯度
- **混合**：可微 alpha 合成

**关键算法**：
$$I_p = \sum_{k=1}^K w_k(\boldsymbol{p}) C_k \prod_{j=1}^{k-1} (1 - w_j(\boldsymbol{p}))$$

其中 $w_k$ 是软化的覆盖权重。

### 12.6.3 JAX 生态系统

JAX 提供函数式编程范式的可微渲染：

**优势**：
- **JIT 编译**：XLA 后端优化
- **vmap 批处理**：自动向量化
- **纯函数式**：易于调试和组合

**典型实现模式**：
```
render_fn = jit(vmap(ray_trace, in_axes=(0, None)))
grad_fn = jit(grad(loss_fn, argnums=1))
```

### 12.6.4 性能优化策略

1. **内存管理**：
   - 检查点技术：trade计算换内存
   - 梯度累积：处理大场景
   - 稀疏表示：只存储非零梯度

2. **计算优化**：
   - 层次采样：先粗后细
   - 重要性缓存：重用采样信息
   - 并行化：光线和像素级并行

3. **数值稳定性**：
   - 对数空间计算：避免下溢
   - 梯度裁剪：防止爆炸
   - 正则化：平滑不连续性

### 12.6.5 框架选择指南

- **Mitsuba 3**：物理准确的逆向渲染
- **PyTorch3D**：神经网络集成，实时应用
- **JAX**：研究原型，自定义算法
- **Taichi**：GPU 编程，自定义内核
- **DIRT/Redner**：特定应用的轻量级方案

## 12.7 本章小结

可微渲染通过使渲染过程可微，为逆向渲染和场景优化提供了计算基础。关键概念包括：

1. **可微光线追踪**：使用隐函数定理计算交点和法线的导数
2. **边缘处理**：通过边缘采样或软化近似处理可见性不连续
3. **阴影计算**：软阴影的可微表示和多重遮挡的稳定处理
4. **偏差-方差权衡**：在梯度估计的准确性和效率之间平衡
5. **实现框架**：根据应用需求选择合适的工具

核心挑战是处理可见性导致的不连续性，解决方案包括边缘积分、重参数化和混合估计器。

## 12.8 练习题

### 基础题

**练习 12.1**：推导球体的光线-表面交点对球心位置的导数。
*提示*：使用隐函数定理，球体方程为 $|\boldsymbol{x} - \boldsymbol{c}|^2 - r^2 = 0$。

<details>
<summary>答案</summary>

设光线 $\boldsymbol{r}(t) = \boldsymbol{o} + t\boldsymbol{d}$，球体 $F(\boldsymbol{x}, \boldsymbol{c}) = |\boldsymbol{x} - \boldsymbol{c}|^2 - r^2 = 0$。

交点条件：$F(\boldsymbol{o} + t^*\boldsymbol{d}, \boldsymbol{c}) = 0$

使用隐函数定理：
$$\frac{\partial t^*}{\partial \boldsymbol{c}} = -\frac{\partial F/\partial \boldsymbol{c}}{\partial F/\partial t} = \frac{2(\boldsymbol{x}^* - \boldsymbol{c})}{2\boldsymbol{d} \cdot (\boldsymbol{x}^* - \boldsymbol{c})}$$

因此：
$$\frac{\partial \boldsymbol{x}^*}{\partial \boldsymbol{c}} = \frac{\partial t^*}{\partial \boldsymbol{c}} \boldsymbol{d} = \frac{(\boldsymbol{x}^* - \boldsymbol{c}) \otimes \boldsymbol{d}}{\boldsymbol{d} \cdot (\boldsymbol{x}^* - \boldsymbol{c})}$$
</details>

**练习 12.2**：证明软化可见性函数 $V_\epsilon(d) = \sigma(d/\epsilon)$ 的梯度在 $\epsilon \to 0$ 时收敛到 Dirac delta 函数。
*提示*：考虑 sigmoid 函数的导数性质。

<details>
<summary>答案</summary>

Sigmoid 函数 $\sigma(x) = 1/(1 + e^{-x})$ 的导数为：
$$\sigma'(x) = \sigma(x)(1 - \sigma(x))$$

对于 $V_\epsilon(d) = \sigma(d/\epsilon)$：
$$\frac{\partial V_\epsilon}{\partial d} = \frac{1}{\epsilon}\sigma'(d/\epsilon) = \frac{1}{\epsilon}\sigma(d/\epsilon)(1 - \sigma(d/\epsilon))$$

当 $\epsilon \to 0$：
- 若 $d > 0$：$\sigma(d/\epsilon) \to 1$，导数 $\to 0$
- 若 $d < 0$：$\sigma(d/\epsilon) \to 0$，导数 $\to 0$
- 若 $d = 0$：导数 $\sim 1/(4\epsilon)$

满足：$\int_{-\infty}^{\infty} \frac{\partial V_\epsilon}{\partial d} dd = 1$，因此收敛到 $\delta(d)$。
</details>

**练习 12.3**：计算三角形边缘在图像平面上的投影长度对顶点位置的导数。
*提示*：使用透视投影公式和链式法则。

<details>
<summary>答案</summary>

设三角形边缘端点为 $\boldsymbol{v}_1, \boldsymbol{v}_2$，透视投影：
$$\boldsymbol{p}_i = \frac{f}{z_i}[x_i, y_i]^T$$

投影长度：
$$L = |\boldsymbol{p}_2 - \boldsymbol{p}_1| = \sqrt{(p_{2x} - p_{1x})^2 + (p_{2y} - p_{1y})^2}$$

对 $\boldsymbol{v}_1$ 的导数：
$$\frac{\partial L}{\partial \boldsymbol{v}_1} = \frac{1}{L}(\boldsymbol{p}_2 - \boldsymbol{p}_1) \cdot \frac{\partial \boldsymbol{p}_1}{\partial \boldsymbol{v}_1}$$

其中：
$$\frac{\partial \boldsymbol{p}_1}{\partial \boldsymbol{v}_1} = \begin{bmatrix}
f/z_1 & 0 & -fx_1/z_1^2 \\
0 & f/z_1 & -fy_1/z_1^2
\end{bmatrix}$$
</details>

### 挑战题

**练习 12.4**：设计一个混合梯度估计器，结合边缘采样和有限差分，证明其收敛性。
*提示*：考虑不同区域使用不同方法，分析整体误差。

<details>
<summary>答案</summary>

混合估计器：
$$\nabla I_{\text{hybrid}} = \begin{cases}
\nabla I_{\text{edge}} & \text{如果 } d < \delta \\
\nabla I_{\text{FD}} & \text{否则}
\end{cases}$$

其中 $d$ 是到最近边缘的距离，$\delta$ 是阈值。

误差分析：
1. 边缘区域：$|\nabla I_{\text{edge}} - \nabla I| = \mathcal{O}(1/\sqrt{N})$（无偏）
2. 内部区域：$|\nabla I_{\text{FD}} - \nabla I| = \mathcal{O}(h^2)$（有偏）

总误差：
$$\mathbb{E}[|\nabla I_{\text{hybrid}} - \nabla I|^2] \leq P_{\text{edge}} \cdot \mathcal{O}(1/N) + P_{\text{interior}} \cdot \mathcal{O}(h^4)$$

选择 $\delta = \mathcal{O}(h^{2/3})$ 可以平衡两种误差。
</details>

**练习 12.5**：推导球面光源软阴影的解析梯度公式。
*提示*：使用立体角和球面几何。

<details>
<summary>答案</summary>

可见立体角：
$$\Omega_{\text{vis}}(\boldsymbol{x}, \boldsymbol{c}) = 2\pi(1 - \cos\alpha) \cdot V(\boldsymbol{x}, \boldsymbol{c})$$

其中 $\cos\alpha = \sqrt{1 - (r/d)^2}$，$d = |\boldsymbol{x} - \boldsymbol{c}|$。

对光源中心 $\boldsymbol{c}$ 的梯度：
$$\frac{\partial \Omega_{\text{vis}}}{\partial \boldsymbol{c}} = 2\pi V \frac{\partial(1-\cos\alpha)}{\partial \boldsymbol{c}} + 2\pi(1-\cos\alpha)\frac{\partial V}{\partial \boldsymbol{c}}$$

第一项（连续部分）：
$$\frac{\partial(1-\cos\alpha)}{\partial \boldsymbol{c}} = \frac{r^2(\boldsymbol{x} - \boldsymbol{c})}{d^3\sqrt{d^2 - r^2}}$$

第二项（边缘项）需要边缘积分处理。
</details>

**练习 12.6**：分析路径追踪中高阶反射的梯度计算复杂度，提出优化策略。
*提示*：考虑路径长度和梯度传播。

<details>
<summary>答案</summary>

$k$ 次反射的梯度复杂度：
- 直接计算：$\mathcal{O}(k \cdot |\theta|)$ 每条路径
- 反向模式AD：$\mathcal{O}(k)$ 计算，$\mathcal{O}(k \cdot |\theta|)$ 内存

优化策略：
1. **路径重用**：存储前向路径，反向只计算梯度
2. **截断近似**：忽略低贡献的高阶项
3. **检查点**：只存储关键状态，重计算中间值
4. **重要性采样**：优先计算高梯度路径

实现伪代码：
```
grad = 0
for path in important_paths:
    if contribution(path) > threshold:
        grad += backprop(path, checkpoint_interval=5)
```
</details>

## 12.9 常见陷阱与错误

1. **忽略边缘梯度**
   - 错误：只计算内部像素的梯度
   - 后果：优化陷入局部最小值
   - 解决：显式边缘采样或软化近似

2. **数值不稳定**
   - 错误：直接计算 $\log(V)$ 当 $V \approx 0$
   - 后果：梯度爆炸或 NaN
   - 解决：使用 `log1p` 或裁剪

3. **梯度偏差累积**
   - 错误：链式使用多个有偏估计器
   - 后果：系统性误差增大
   - 解决：混合无偏估计或误差补偿

4. **内存爆炸**
   - 错误：存储所有中间变量用于反向传播
   - 后果：大场景无法处理
   - 解决：检查点技术或梯度累积

5. **采样效率低**
   - 错误：均匀采样边缘
   - 后果：高方差，收敛慢
   - 解决：基于重要性的自适应采样

## 12.10 最佳实践检查清单

### 算法设计
- [ ] 明确处理可见性不连续（边缘采样/软化）
- [ ] 选择合适的偏差-方差权衡点
- [ ] 实现数值稳定的梯度计算
- [ ] 考虑多尺度/层次化方法

### 实现优化
- [ ] 使用自动微分而非手动推导
- [ ] 实现高效的边缘检测和采样
- [ ] 优化内存使用（检查点、稀疏表示）
- [ ] 并行化光线和梯度计算

### 验证测试
- [ ] 梯度检查（有限差分验证）
- [ ] 收敛性测试（不同采样率）
- [ ] 边界情况（退化几何、极端光照）
- [ ] 与解析解对比（简单场景）

### 应用集成
- [ ] 选择适合任务的框架
- [ ] 设计合理的损失函数
- [ ] 实现早停和正则化
- [ ] 监控优化过程（梯度范数、损失曲线）
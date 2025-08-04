# 第12章：可微渲染 (Differentiable Rendering)

可微渲染是逆向渲染的计算基础，它使得我们能够通过梯度下降优化场景参数。本章深入探讨如何使渲染过程可微，特别是处理可见性不连续性带来的挑战。我们将从可微光线追踪的基础开始，逐步深入到边缘采样、软阴影等高级主题，并讨论实际实现中的偏差-方差权衡。

## 12.1 引言

渲染本质上是求解积分方程，但当我们需要对场景参数（几何、材质、光照）求导时，传统渲染算法面临根本性挑战。考虑像素强度的一般形式：

$$I(\boldsymbol{\theta}) = \int_{\Omega} f(\boldsymbol{x}, \boldsymbol{\theta}) \, d\boldsymbol{x}$$

其中 $\boldsymbol{\theta}$ 是场景参数。可微渲染的核心问题是计算 $\frac{\partial I}{\partial \boldsymbol{\theta}}$。

### 12.1.1 可微性的挑战

主要挑战来自两个方面：

1. **可见性不连续**：当物体边缘移动时，积分域 $\Omega$ 随参数变化
2. **高维积分**：路径追踪涉及高维积分，直接微分计算代价高昂

对于连续被积函数，Leibniz积分规则给出：

$$\frac{\partial}{\partial \theta} \int_{\Omega(\theta)} f(x, \theta) \, dx = \int_{\Omega(\theta)} \frac{\partial f}{\partial \theta} \, dx + \int_{\partial \Omega(\theta)} f \cdot v_n \, ds$$

第二项是边界项，处理积分域的变化。

### 12.1.2 可微渲染的应用

可微渲染在以下领域有重要应用：

- **3D重建**：从2D图像恢复3D几何和材质
- **场景理解**：推断光照、材质分解
- **内容创作**：通过优化生成满足约束的3D内容
- **机器人视觉**：主动感知和场景操作

## 12.2 可微光线追踪

可微光线追踪是可微渲染的基础。我们需要计算光线-物体交点及相关量对场景参数的导数。

### 12.2.1 光线-表面交点的导数

考虑参数化光线 $\boldsymbol{r}(t) = \boldsymbol{o} + t\boldsymbol{d}$，其中 $\boldsymbol{o}$ 是原点，$\boldsymbol{d}$ 是方向。对于隐式表面 $F(\boldsymbol{x}, \boldsymbol{\theta}) = 0$，交点满足：

$$F(\boldsymbol{o} + t^*\boldsymbol{d}, \boldsymbol{\theta}) = 0$$

使用隐函数定理，交点参数 $t^*$ 对 $\boldsymbol{\theta}$ 的导数为：

$$\frac{\partial t^*}{\partial \boldsymbol{\theta}} = -\frac{\partial F/\partial \boldsymbol{\theta}}{\boldsymbol{d} \cdot \nabla F}$$

交点位置的导数：

$$\frac{\partial \boldsymbol{x}^*}{\partial \boldsymbol{\theta}} = \frac{\partial t^*}{\partial \boldsymbol{\theta}} \boldsymbol{d} + \frac{\partial \boldsymbol{o}}{\partial \boldsymbol{\theta}} + t^* \frac{\partial \boldsymbol{d}}{\partial \boldsymbol{\theta}}$$

### 12.2.2 表面法线的微分

表面法线 $\boldsymbol{n} = \nabla F / |\nabla F|$ 的导数需要考虑归一化：

$$\frac{\partial \boldsymbol{n}}{\partial \boldsymbol{\theta}} = \frac{1}{|\nabla F|} \left( \frac{\partial \nabla F}{\partial \boldsymbol{\theta}} - \boldsymbol{n} \otimes \boldsymbol{n} \cdot \frac{\partial \nabla F}{\partial \boldsymbol{\theta}} \right)$$

其中 $\otimes$ 表示外积。Hessian矩阵 $\frac{\partial \nabla F}{\partial \boldsymbol{\theta}}$ 包含二阶导数信息。

### 12.2.3 参数化表面的处理

对于参数化表面 $\boldsymbol{x}(u, v, \boldsymbol{\theta})$，我们需要同时追踪参数坐标的变化。设交点参数为 $(u^*, v^*)$，则：

$$\begin{bmatrix}
\frac{\partial u^*}{\partial \boldsymbol{\theta}} \\
\frac{\partial v^*}{\partial \boldsymbol{\theta}}
\end{bmatrix} = -\boldsymbol{J}^{-1} \frac{\partial \boldsymbol{r}}{\partial \boldsymbol{\theta}}$$

其中 $\boldsymbol{J}$ 是关于 $(u, v)$ 的Jacobian矩阵。

### 12.2.4 反射和折射方向的导数

反射方向 $\boldsymbol{r} = \boldsymbol{i} - 2(\boldsymbol{i} \cdot \boldsymbol{n})\boldsymbol{n}$ 的导数：

$$\frac{\partial \boldsymbol{r}}{\partial \boldsymbol{\theta}} = \frac{\partial \boldsymbol{i}}{\partial \boldsymbol{\theta}} - 2\left[ (\boldsymbol{i} \cdot \boldsymbol{n})\frac{\partial \boldsymbol{n}}{\partial \boldsymbol{\theta}} + \boldsymbol{n} \otimes \left(\frac{\partial \boldsymbol{i}}{\partial \boldsymbol{\theta}} \cdot \boldsymbol{n} + \boldsymbol{i} \cdot \frac{\partial \boldsymbol{n}}{\partial \boldsymbol{\theta}}\right) \right]$$

折射方向使用Snell定律的向量形式，其导数计算更复杂但遵循类似原理。

## 12.3 边缘采样与重参数化

可见性不连续是可微渲染的核心挑战。当物体移动时，其轮廓边缘导致积分域的变化，产生Dirac delta函数形式的梯度。

### 12.3.1 边缘积分理论

考虑依赖于参数 $\theta$ 的2D积分域 $\Omega(\theta)$：

$$I(\theta) = \int_{\Omega(\theta)} f(\boldsymbol{x}) \, d\boldsymbol{x}$$

使用Reynolds传输定理，导数为：

$$\frac{dI}{d\theta} = \int_{\Omega(\theta)} \frac{\partial f}{\partial \theta} \, d\boldsymbol{x} + \oint_{\partial\Omega(\theta)} f \, v_n \, ds$$

其中 $v_n$ 是边界的法向速度。第二项是边缘项，在渲染中对应物体轮廓的贡献。

### 12.3.2 可见性函数的导数

定义可见性函数 $V(\boldsymbol{x}, \boldsymbol{y})$：

$$V(\boldsymbol{x}, \boldsymbol{y}) = \begin{cases}
1 & \text{如果 } \boldsymbol{x} \text{ 和 } \boldsymbol{y} \text{ 互相可见} \\
0 & \text{否则}
\end{cases}$$

其导数包含Dirac delta函数：

$$\frac{\partial V}{\partial \theta} = \delta(\text{dist}_{\text{edge}}) \cdot v_{\perp}$$

其中 $\text{dist}_{\text{edge}}$ 是到轮廓边缘的距离，$v_{\perp}$ 是边缘的垂直速度。

### 12.3.3 边缘采样策略

为了处理delta函数，我们需要显式采样轮廓边缘。对于三角网格，轮廓边缘满足：

$$(\boldsymbol{n}_1 \cdot \boldsymbol{v}) \cdot (\boldsymbol{n}_2 \cdot \boldsymbol{v}) < 0$$

其中 $\boldsymbol{n}_1, \boldsymbol{n}_2$ 是相邻面的法线，$\boldsymbol{v}$ 是视线方向。

边缘上的积分可以参数化为：

$$\int_{\text{edge}} g(s) \, ds = \int_0^L g(s(t)) \left|\frac{ds}{dt}\right| dt$$

### 12.3.4 重参数化技巧

重参数化将不可微的采样过程转换为可微形式。对于边缘采样，我们使用：

1. **边缘参数化**：将3D边缘投影到2D图像平面
2. **重要性采样**：根据边缘对像素的贡献分配样本
3. **解析边缘积分**：对简单情况（如多边形遮挡）使用闭式解

例如，对于线性边缘和常数被积函数，边缘积分有解析解：

$$\int_{\text{edge}} f_0 \, ds = f_0 \cdot \text{length}(\text{edge}) \cdot \cos(\phi)$$

其中 $\phi$ 是边缘与投影方向的夹角。

## 12.4 可微阴影与可见性

阴影是可见性的特殊情况，涉及光源、遮挡物和接收表面三者的关系。可微阴影计算需要正确处理软阴影的连续变化。

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
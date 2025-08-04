# 第13章：材质与几何重建

本章探讨逆向渲染的核心问题：从观察图像中恢复场景的材质属性和几何形状。我们将建立在前面章节的可微渲染基础上，展示如何通过优化框架同时重建材质和几何。这些技术在计算机视觉、虚拟现实和数字资产创建中有着广泛应用。

## 学习目标

完成本章后，您将能够：
1. 从多视图图像估计BRDF和BSSRDF参数
2. 理解并推导shape-from-shading的变分公式
3. 分析多视图立体重建中的优化景观
4. 设计联合材质-几何优化的能量函数
5. 应用物理约束提高重建质量
6. 评估不同先验对重建结果的影响

## 13.1 BRDF/BSSRDF估计

### 13.1.1 逆向渲染方程

给定观察图像 $I(\mathbf{x})$，我们寻求恢复表面反射属性。渲染方程的逆问题可表述为：

$$\min_{\rho} \sum_{\mathbf{x}} \left\| I(\mathbf{x}) - \int_{\Omega} L(\mathbf{x}, \boldsymbol{\omega}_o) d\boldsymbol{\omega}_o \right\|^2$$

其中出射辐射度由BRDF $\rho$ 决定：

$$L(\mathbf{x}, \boldsymbol{\omega}_o) = \int_{\Omega} \rho(\mathbf{x}, \boldsymbol{\omega}_i, \boldsymbol{\omega}_o) L_i(\mathbf{x}, \boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n}) d\boldsymbol{\omega}_i$$

### 13.1.2 参数化BRDF模型

实际应用中，我们通常采用参数化BRDF模型：

**微表面模型**：
$$\rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \frac{F(\boldsymbol{\omega}_i, \mathbf{h}) G(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) D(\mathbf{h})}{4(\boldsymbol{\omega}_i \cdot \mathbf{n})(\boldsymbol{\omega}_o \cdot \mathbf{n})}$$

其中：
- $F$: Fresnel项，参数为折射率 $\eta$
- $G$: 几何衰减，参数为粗糙度 $\alpha$
- $D$: 法线分布，如GGX分布
- $\mathbf{h} = (\boldsymbol{\omega}_i + \boldsymbol{\omega}_o)/\|\boldsymbol{\omega}_i + \boldsymbol{\omega}_o\|$

### 13.1.3 梯度计算

对BRDF参数 $\boldsymbol{\theta} = [\alpha, \eta, \mathbf{k}_d, \mathbf{k}_s]$ 的梯度：

$$\frac{\partial \mathcal{L}}{\partial \boldsymbol{\theta}} = -2\sum_{\mathbf{x}} (I(\mathbf{x}) - \hat{I}(\mathbf{x})) \frac{\partial \hat{I}(\mathbf{x})}{\partial \boldsymbol{\theta}}$$

其中渲染图像对参数的导数通过链式法则计算：

$$\frac{\partial \hat{I}}{\partial \boldsymbol{\theta}} = \int_{\Omega} \frac{\partial \rho}{\partial \boldsymbol{\theta}} L_i (\boldsymbol{\omega}_i \cdot \mathbf{n}) d\boldsymbol{\omega}_i$$

### 13.1.4 BSSRDF估计

对于半透明材质，需要考虑次表面散射：

$$L_o(\mathbf{x}_o, \boldsymbol{\omega}_o) = \int_A \int_{\Omega} S(\mathbf{x}_o, \boldsymbol{\omega}_o, \mathbf{x}_i, \boldsymbol{\omega}_i) L_i(\mathbf{x}_i, \boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n}_i) d\boldsymbol{\omega}_i dA$$

偶极子近似下的BSSRDF：
$$S(\mathbf{x}_o, \mathbf{x}_i) = \frac{1}{\pi} F_t(\eta) R_d(\|\mathbf{x}_o - \mathbf{x}_i\|) F_t(\eta)$$

其中 $R_d(r)$ 为扩散剖面，依赖于吸收系数 $\sigma_a$ 和散射系数 $\sigma_s$。

## 13.2 形状从明暗恢复

### 13.2.1 经典形状从明暗

对于朗伯表面在正交投影下，亮度方程为：

$$I(x,y) = \rho \mathbf{n}(x,y) \cdot \mathbf{l}$$

其中表面法线与深度函数 $z(x,y)$ 的关系：

$$\mathbf{n} = \frac{(-\partial z/\partial x, -\partial z/\partial y, 1)}{\sqrt{1 + (\partial z/\partial x)^2 + (\partial z/\partial y)^2}}$$

### 13.2.2 变分公式

能量泛函包含数据项和正则化项：

$$E[z] = \int_{\Omega} (I(x,y) - \rho \mathbf{n}[z] \cdot \mathbf{l})^2 dxdy + \lambda \int_{\Omega} |\nabla z|^2 dxdy$$

欧拉-拉格朗日方程：

$$\frac{\partial}{\partial z}\left(\frac{(I - \rho \mathbf{n} \cdot \mathbf{l})^2}{\sqrt{1 + |\nabla z|^2}}\right) - \lambda \nabla^2 z = 0$$

### 13.2.3 多光源形状从明暗

给定 $m$ 个光照条件下的图像 $\{I_j\}_{j=1}^m$：

$$E[z] = \sum_{j=1}^m \int_{\Omega} (I_j - \rho \mathbf{n}[z] \cdot \mathbf{l}_j)^2 dxdy$$

最优性条件产生非线性PDE系统。

### 13.2.4 光度立体

当 $m \geq 3$ 时，可以线性求解法线：

$$\begin{bmatrix} I_1 \\ I_2 \\ \vdots \\ I_m \end{bmatrix} = \rho \begin{bmatrix} \mathbf{l}_1^T \\ \mathbf{l}_2^T \\ \vdots \\ \mathbf{l}_m^T \end{bmatrix} \mathbf{n}$$

通过最小二乘：$\mathbf{n} = \frac{1}{\rho}(\mathbf{L}^T\mathbf{L})^{-1}\mathbf{L}^T\mathbf{I}$

## 13.3 多视图立体重建中的优化

### 13.3.1 体积雕刻能量

定义占用场 $\phi: \mathbb{R}^3 \rightarrow \{0,1\}$，光度一致性能量：

$$E[\phi] = \sum_{i,j} \int_{\partial\Omega[\phi]} \|I_i(\pi_i(\mathbf{x})) - I_j(\pi_j(\mathbf{x}))\|^2 dS$$

其中 $\pi_i$ 为相机 $i$ 的投影函数。

### 13.3.2 变分水平集方法

使用隐式表面 $\{\mathbf{x}: \phi(\mathbf{x}) = 0\}$：

$$E[\phi] = \sum_{i,j} \int_{\mathbb{R}^3} \|I_i - I_j\|^2 \delta(\phi) |\nabla\phi| d\mathbf{x}$$

演化方程：
$$\frac{\partial \phi}{\partial t} = \delta(\phi) \left[ \sum_{i,j} (I_i - I_j)(\nabla I_i - \nabla I_j) \cdot \frac{\nabla\phi}{|\nabla\phi|} + \lambda \kappa \right]$$

### 13.3.3 PatchMatch立体

离散优化公式，为每个像素分配深度 $d_p$：

$$E(\{d_p\}) = \sum_p C(p, d_p) + \sum_{(p,q) \in \mathcal{N}} V(d_p, d_q)$$

其中：
- $C(p,d)$: 匹配代价
- $V(d_p,d_q)$: 平滑项

### 13.3.4 平面扫描体积

构建代价体积 $\mathcal{C}(x,y,d)$：

$$\mathcal{C}(x,y,d) = \frac{1}{|\mathcal{V}|} \sum_{i \in \mathcal{V}} \|I_{ref}(x,y) - I_i(\mathbf{H}_i(x,y,d))\|$$

其中 $\mathbf{H}_i$ 为深度 $d$ 处的单应变换。

## 13.4 联合材质-几何优化

### 13.4.1 耦合优化框架

同时优化几何 $\mathcal{G}$ 和材质 $\mathcal{M}$：

$$E[\mathcal{G}, \mathcal{M}] = E_{photo}[\mathcal{G}, \mathcal{M}] + \lambda_g E_{geom}[\mathcal{G}] + \lambda_m E_{mat}[\mathcal{M}]$$

### 13.4.2 交替优化

**几何步骤**（固定材质）：
$$\mathcal{G}^{(k+1)} = \arg\min_{\mathcal{G}} E[\mathcal{G}, \mathcal{M}^{(k)}]$$

**材质步骤**（固定几何）：
$$\mathcal{M}^{(k+1)} = \arg\min_{\mathcal{M}} E[\mathcal{G}^{(k+1)}, \mathcal{M}]$$

### 13.4.3 联合梯度下降

计算完整梯度：
$$\begin{bmatrix} \mathcal{G}^{(k+1)} \\ \mathcal{M}^{(k+1)} \end{bmatrix} = \begin{bmatrix} \mathcal{G}^{(k)} \\ \mathcal{M}^{(k)} \end{bmatrix} - \alpha \begin{bmatrix} \nabla_{\mathcal{G}} E \\ \nabla_{\mathcal{M}} E \end{bmatrix}$$

需要考虑Hessian的条件数。

### 13.4.4 分解歧义性

材质-几何-光照分解存在固有歧义：
- 尺度歧义：$(\alpha\mathcal{M}, \mathcal{L}/\alpha)$ 产生相同图像
- 低频歧义：光照与反照率的低频分量难以区分

## 13.5 物理约束与先验

### 13.5.1 材质物理约束

**能量守恒**：
$$\int_{\Omega} \rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o \leq 1$$

**Helmholtz互易性**：
$$\rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \rho(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)$$

**单调性约束**（对于粗糙度）：
$$\alpha_1 < \alpha_2 \Rightarrow D_{\alpha_1}(\mathbf{h}) > D_{\alpha_2}(\mathbf{h}), \forall \mathbf{h} \neq \mathbf{n}$$

### 13.5.2 几何先验

**最小表面积**：
$$E_{area}[\mathcal{S}] = \int_{\mathcal{S}} dS$$

**平均曲率流**：
$$E_{smooth}[\mathcal{S}] = \int_{\mathcal{S}} H^2 dS$$

其中 $H = (\kappa_1 + \kappa_2)/2$ 为平均曲率。

**体积保持**：
$$\int_{\Omega} \phi d\mathbf{x} = V_0$$

### 13.5.3 统计先验

**材质分布先验**（对数正态）：
$$p(\alpha) = \frac{1}{\alpha\sigma\sqrt{2\pi}} \exp\left(-\frac{(\ln\alpha - \mu)^2}{2\sigma^2}\right)$$

**空间相关性**（MRF）：
$$p(\mathcal{M}) \propto \exp\left(-\sum_{(i,j) \in \mathcal{N}} \psi(\mathcal{M}_i, \mathcal{M}_j)\right)$$

### 13.5.4 学习先验

使用神经网络编码先验：
$$E_{prior}[\mathcal{M}] = \|\mathcal{M} - G_\theta(z)\|^2$$

其中 $G_\theta$ 为预训练的材质生成器。

## 本章小结

本章介绍了材质与几何重建的核心技术：

**关键概念**：
- BRDF/BSSRDF参数估计的优化框架
- 形状从明暗的变分公式和数值解法
- 多视图立体的连续和离散优化方法
- 联合材质-几何优化的耦合问题
- 物理约束和统计先验的作用

**核心方程**：
1. **逆向渲染**：$\min_{\rho} \|I - \mathcal{R}[\rho]\|^2$
2. **形状从明暗**：$(I - \rho\mathbf{n}\cdot\mathbf{l})^2 + \lambda|\nabla z|^2$
3. **多视图立体**：$\sum_{i,j}\|I_i(\pi_i(\mathbf{x})) - I_j(\pi_j(\mathbf{x}))\|^2$
4. **联合优化**：$E_{photo} + \lambda_g E_{geom} + \lambda_m E_{mat}$

## 常见陷阱与错误

1. **局部最小值**：非凸优化容易陷入局部最优，需要良好初始化
2. **尺度歧义**：忽略归一化导致材质-光照分解不稳定
3. **过拟合**：参数过多时容易过拟合观察数据
4. **数值稳定性**：法线计算中分母接近零需要正则化
5. **边界处理**：形状边界的可见性变化需要特殊处理

## 最佳实践检查清单

- [ ] 选择合适的BRDF参数化以平衡表达能力和稳定性
- [ ] 使用多尺度优化策略避免局部最小值
- [ ] 加入物理约束确保结果合理性
- [ ] 考虑光照-材质-几何的耦合关系
- [ ] 验证优化算法的收敛性
- [ ] 评估重建结果的不确定性

## 练习题

### 练习 13.1：微表面BRDF的能量守恒
证明GGX微表面模型满足能量守恒约束。具体地，证明：
$$\int_{\Omega} \rho_{GGX}(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o \leq 1$$

**提示**：利用微表面理论中的遮蔽-阴影函数 $G$ 的性质。

<details>
<summary>解答</summary>

对于微表面BRDF：
$$\rho = \frac{F(\boldsymbol{\omega}_i, \mathbf{h}) G(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) D(\mathbf{h})}{4(\boldsymbol{\omega}_i \cdot \mathbf{n})(\boldsymbol{\omega}_o \cdot \mathbf{n})}$$

积分可以变换到半角向量空间：
$$\int_{\Omega} \rho (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o = \int_{\Omega_h} F G D \frac{(\mathbf{h} \cdot \boldsymbol{\omega}_i)}{(\boldsymbol{\omega}_i \cdot \mathbf{n})} d\mathbf{h}$$

利用Smith遮蔽函数的性质：
$$\int_{\Omega} G(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) \frac{(\mathbf{h} \cdot \boldsymbol{\omega}_i)}{(\boldsymbol{\omega}_i \cdot \mathbf{n})} D(\mathbf{h}) d\mathbf{h} = 1$$

由于 $F \leq 1$，因此：
$$\int_{\Omega} \rho (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o \leq 1$$
</details>

### 练习 13.2：形状从明暗的唯一性
考虑朗伯表面的形状从明暗问题。证明在以下条件下解不唯一：
1. 单一光源
2. 无边界条件

构造两个不同的表面产生相同的图像。

**提示**：考虑凸/凹歧义。

<details>
<summary>解答</summary>

给定亮度方程：$I(x,y) = \mathbf{n}(x,y) \cdot \mathbf{l}$

设光源方向 $\mathbf{l} = (0, 0, 1)$（顶光照明）。则：
$$I(x,y) = \frac{1}{\sqrt{1 + p^2 + q^2}}$$

其中 $p = \partial z/\partial x$, $q = \partial z/\partial y$。

两个解：
1. 凸解：$z_1(x,y) = \sqrt{R^2 - x^2 - y^2}$（上半球）
2. 凹解：$z_2(x,y) = -\sqrt{R^2 - x^2 - y^2}$（下半球）

两者产生相同的图像，因为法线的 $z$ 分量相同：
$$n_z = \frac{1}{\sqrt{1 + (x/z)^2 + (y/z)^2}} = \frac{|z|}{\sqrt{x^2 + y^2 + z^2}} = \frac{|z|}{R}$$
</details>

### 练习 13.3：光度立体的最优光源配置
给定 $n = 3$ 个光源，求使光度立体重建最稳定的光源方向配置。稳定性用条件数 $\kappa(\mathbf{L}^T\mathbf{L})$ 衡量。

**提示**：最小化条件数等价于最大化最小特征值。

<details>
<summary>解答</summary>

光源矩阵：
$$\mathbf{L} = \begin{bmatrix} \mathbf{l}_1^T \\ \mathbf{l}_2^T \\ \mathbf{l}_3^T \end{bmatrix}$$

条件数：$\kappa = \lambda_{max}/\lambda_{min}$

对于 $n = 3$，最优配置使 $\mathbf{L}^T\mathbf{L} = I$，即光源正交：
$$\mathbf{l}_i \cdot \mathbf{l}_j = \delta_{ij}$$

一个最优配置：
$$\mathbf{l}_1 = (1, 0, 0), \quad \mathbf{l}_2 = (0, 1, 0), \quad \mathbf{l}_3 = (0, 0, 1)$$

此时 $\kappa = 1$（最小可能值）。

实际中考虑阴影，通常选择：
$$\mathbf{l}_i = (\sin\theta\cos(2\pi i/3), \sin\theta\sin(2\pi i/3), \cos\theta)$$
其中 $\theta \approx 45°$。
</details>

### 练习 13.4：BSSRDF的扩散近似误差
评估偶极子近似对薄材质的误差。设材质厚度为 $d$，推导当 $d \ll l_t$（输运平均自由程）时的相对误差。

**提示**：比较精确解（平板几何）与偶极子近似。

<details>
<summary>解答</summary>

平板几何的精确反射率：
$$R_{exact} = \frac{1 - \exp(-2\tau)}{1 + \exp(-2\tau)}$$
其中 $\tau = d/l_t$ 为光学厚度。

偶极子近似：
$$R_{dipole} \approx \frac{A}{1 + A}$$
其中 $A = (1 + F_{dr})/(1 - F_{dr})$，$F_{dr}$ 为内部反射系数。

泰勒展开对小 $\tau$：
$$R_{exact} \approx \tau - \frac{2\tau^3}{3} + O(\tau^5)$$
$$R_{dipole} \approx \tau - \tau^2 + O(\tau^3)$$

相对误差：
$$\epsilon = \frac{|R_{dipole} - R_{exact}|}{R_{exact}} \approx \tau = \frac{d}{l_t}$$

因此误差与厚度成正比，薄材质时偶极子近似失效。
</details>

### 练习 13.5：联合优化的收敛性分析
考虑简化的联合材质-几何优化：
$$E(\alpha, z) = (I - \alpha h(z))^2 + \lambda z^2$$
其中 $h(z)$ 是 $z$ 的非线性函数。分析交替优化的收敛条件。

**提示**：构造Lyapunov函数。

<details>
<summary>解答</summary>

交替优化步骤：
1. 固定 $z$：$\alpha^{(k+1)} = I h(z^{(k)})/h^2(z^{(k)})$
2. 固定 $\alpha$：$z^{(k+1)} = \arg\min_z (I - \alpha^{(k+1)} h(z))^2 + \lambda z^2$

定义增广能量：
$$\tilde{E}(\alpha, z, \beta) = (I - \alpha h(z))^2 + \lambda z^2 + \mu(\alpha - \beta)^2$$

可以证明：
$$E(\alpha^{(k+1)}, z^{(k+1)}) \leq E(\alpha^{(k)}, z^{(k)})$$

收敛条件：
1. $h(z)$ Lipschitz连续：$|h(z_1) - h(z_2)| \leq L|z_1 - z_2|$
2. $\lambda > L^2 I^2/(4h_{min}^2)$

其中 $h_{min} = \min_z |h(z)|$。
</details>

### 练习 13.6：多视图立体的基线选择
推导多视图立体中深度估计误差与相机基线的关系。设深度 $z$，基线 $b$，视差误差 $\Delta d$。

**提示**：使用三角测量的误差传播。

<details>
<summary>解答</summary>

三角测量关系：
$$z = \frac{fb}{d}$$

其中 $f$ 为焦距，$d$ 为视差。

深度误差：
$$\Delta z = \frac{\partial z}{\partial d} \Delta d = -\frac{fb}{d^2} \Delta d = -\frac{z^2}{fb} \Delta d$$

相对误差：
$$\frac{\Delta z}{z} = \frac{z}{fb} \Delta d$$

结论：
1. 深度误差与深度平方成正比
2. 增大基线 $b$ 减小误差
3. 但基线过大导致遮挡和匹配困难

最优基线选择需平衡精度和匹配可靠性：
$$b_{opt} \approx 0.1 \cdot z_{avg}$$
</details>

### 练习 13.7：材质先验的信息论分析
使用最大熵原理推导粗糙度参数 $\alpha \in [0,1]$ 的先验分布。已知约束：$E[\alpha] = \mu$，$E[\log\alpha] = \nu$。

**提示**：拉格朗日乘数法求解约束最大熵问题。

<details>
<summary>解答</summary>

最大熵问题：
$$\max_{p(\alpha)} H[p] = -\int_0^1 p(\alpha) \log p(\alpha) d\alpha$$

约束条件：
1. $\int_0^1 p(\alpha) d\alpha = 1$
2. $\int_0^1 \alpha p(\alpha) d\alpha = \mu$
3. $\int_0^1 \log\alpha \cdot p(\alpha) d\alpha = \nu$

拉格朗日函数：
$$\mathcal{L} = H[p] + \lambda_0(1 - \int p) + \lambda_1(\mu - \int \alpha p) + \lambda_2(\nu - \int \log\alpha \cdot p)$$

变分得：
$$p(\alpha) = \frac{1}{Z} \exp(-\lambda_1 \alpha - \lambda_2 \log\alpha) = \frac{1}{Z} \alpha^{-\lambda_2} \exp(-\lambda_1 \alpha)$$

这是广义逆高斯分布的特例。当 $\lambda_2 = 1$ 时，退化为指数分布。
</details>

### 练习 13.8：物理约束的投影算法
设计投影算子 $\Pi$ 将任意函数投影到满足互易性的BRDF空间。定义合适的度量并证明投影的最优性。

**提示**：考虑 $L^2$ 度量下的投影。

<details>
<summary>解答</summary>

互易性约束：$\rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \rho(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)$

对任意函数 $f$，定义投影：
$$\Pi[f] = \arg\min_{\rho \in \mathcal{R}} \int_{\Omega^2} |f(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) - \rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o)|^2 d\boldsymbol{\omega}_i d\boldsymbol{\omega}_o$$

其中 $\mathcal{R}$ 为互易BRDF集合。

解为对称化：
$$\Pi[f](\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \frac{1}{2}[f(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) + f(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)]$$

证明最优性：设 $\rho \in \mathcal{R}$，则：
$$\|f - \rho\|^2 = \|\frac{f + f^T}{2} - \rho\|^2 + \|\frac{f - f^T}{2}\|^2 \geq \|f - \Pi[f]\|^2$$

因为 $\rho$ 与 $(f - f^T)/2$ 正交。
</details>
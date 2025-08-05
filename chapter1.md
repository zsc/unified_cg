# 第1章：几何光学与渲染基础

本章通过几何光学的视角建立计算机图形学的数学基础。我们将渲染方程作为核心框架，介绍关键的辐射度量概念，并建立路径积分表述，这将统一所有后续的渲染技术。通过将光传输视为高维积分问题，我们为理解基于点的渲染、基于图像的渲染和神经渲染方法作为解决同一基本方程的不同方法奠定了基础。

## 学习目标

完成本章后，您将能够：
1. 使用能量守恒从第一性原理推导渲染方程
2. 在保持辐射度量不变的情况下在不同坐标系之间转换
3. 分析BRDF属性并验证物理合理性
4. 应用蒙特卡洛方法估计具有已知误差界限的高维积分
5. 将光传输表达为路径积分并将其与体积渲染联系起来

## 1.1 光线追踪基础与渲染方程

### 几何光学中的光线

在几何光学中，我们使用光线来建模光的传播——在均匀介质中沿直线传播的无限细光束。当波长 $\lambda \ll$ 特征尺寸时，这种近似成立，使我们可以忽略衍射和干涉。光线参数化为：

$\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$

其中 $\mathbf{o} \in \mathbb{R}^3$ 是原点，$\mathbf{d} \in \mathbb{R}^3$ 是方向（$||\mathbf{d}|| = 1$），$t \ge 0$ 是沿光线的参数。

光线方程源自程函方程 $\nabla S = n\mathbf{\hat{k}}$ 在 $\lambda \to 0$ 极限下的结果，其中 $S$ 是相位，$n$ 是折射率。在非均匀介质中，光线遵循满足以下方程的曲线路径：

$\frac{d}{ds}\left(n \frac{d\mathbf{r}}{ds}\right) = \nabla n$

当 $n$ 为常数时，这简化为直线。

**与波动光学的联系**：几何光学近似源自波动方程的WKB（Wentzel-Kramers-Brillouin）近似。当我们将 $\psi = A \exp(ikS)$ 代入亥姆霍兹方程并取 $k \to \infty$ 时：

$\nabla^2\psi + k^2n^2\psi = 0 \to (\nabla S)^2 = n^2$ （程函方程）

等相位面 $S = \text{const}$ 是波前，光线是这些波前的正交轨迹。当我们在后续章节扩展到波动光学时，这种联系变得至关重要。

**光线光学的有效性**：几何光学近似在以下情况下失效：
1. 特征尺寸 $\sim$ 波长（衍射变得显著）
2. 靠近焦散面，光线密度 $\to \infty$
3. 存在尖锐边缘或不连续性
4. 需要相位信息的相干现象

**费马原理**：光线遵循光程稳定的路径：

$\delta\int n(\mathbf{r}) ds = 0$

这个变分原理统一了光线行为：在均匀介质中的直线、界面处的斯涅尔定律以及梯度折射率介质中的曲线路径。它还与物理学中的最小作用量原理和微分几何中的测地线方程相关联。

### 辐射度量

在推导渲染方程之前，我们必须建立辐射度量框架。这些量形成一个层次结构，每个都建立在前一个的基础上：

**辐射能** $Q$ 测量总电磁能量：
$Q \text{ [J]}$

**辐射通量（功率）** $\Phi$ 测量单位时间的能量：
$\Phi = \frac{dQ}{dt} \text{ [W]}$

**辐射强度** $I$ 测量点源单位立体角的通量：
$I = \frac{d\Phi}{d\omega} \text{ [W/sr]}$

**辐照度** $E$ 测量入射单位面积的通量：
$E = \frac{d\Phi}{dA} \text{ [W/m}^2\text{]}$

**辐射出射度** $M$ 测量离开单位面积的通量：
$M = \frac{d\Phi}{dA} \text{ [W/m}^2\text{]}$

**辐射率** $L$ 测量单位面积单位立体角的通量：
$L = \frac{d^2\Phi}{dA \cos \theta d\omega} = \frac{d^2\Phi}{dA_\perp d\omega} \text{ [W/(m}^2\text{·sr)]}$

其中 $dA_\perp = dA \cos \theta$ 是垂直于光线方向的投影面积。

辐射率是渲染中的基本量，因为：
1. 它在真空中沿光线保持恒定（辐射率不变性定理）
2. 它是相机和眼睛测量的量
3. 所有其他辐射度量都可以从它导出

**光度量与辐射度量**：虽然我们关注辐射度量（物理能量），但为人类感知进行渲染时通常使用光度量：
- 光通量 [lm] = 辐射通量 [W] $\times$ 光效能
- 亮度 [cd/m$^2$] = 辐射率 $\times$ 明视觉响应 $V(\lambda)$
- CIE光度函数 $V(\lambda)$ 在555 nm（绿色）处达到峰值

**光谱辐射率**：实际上，辐射率随波长变化：
$L(\mathbf{x}, \omega, \lambda) \text{ [W/(m}^2\text{·sr·nm)]}$

对于渲染，我们通常使用：
- RGB近似：光谱的3个样本
- 光谱渲染：N个波长样本（通常10-100）
- 主波长：光谱的随机采样

**相干与非相干叠加**：辐射度量假设非相干光——强度直接相加。对于相干光源（激光），我们必须跟踪相位并叠加复振幅：
$I_{\text{total}} = |E_1 + E_2|^2 \neq |E_1|^2 + |E_2|^2$ （一般情况下）

### 渲染方程

渲染方程由Kajiya（1986）提出，描述场景中光的平衡分布。它源于功率平衡：在任何表面点，输出功率等于发射功率加反射功率。

在具有法线 $\mathbf{n}$ 的表面点 $\mathbf{x}$ 处，沿方向 $\omega_o$ 的出射辐射率 $L_o$ 满足：

$L_o(\mathbf{x}, \omega_o) = L_e(\mathbf{x}, \omega_o) + \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) L_i(\mathbf{x}, \omega_i) (\omega_i \cdot \mathbf{n}) d\omega_i$

其中：
- $L_e(\mathbf{x}, \omega_o)$ 是从 $\mathbf{x}$ 沿方向 $\omega_o$ 的发射辐射率
- $f_r(\mathbf{x}, \omega_i, \omega_o)$ 是BRDF [sr$^{-1}$]
- $L_i(\mathbf{x}, \omega_i)$ 是在 $\mathbf{x}$ 处来自方向 $\omega_i$ 的入射辐射率
- $\Omega$ 是 $\mathbf{x}$ 上方的半球（其中 $\omega \cdot \mathbf{n} > 0$）
- $(\omega_i \cdot \mathbf{n}) = \cos \theta_i$ 考虑了投影面积

积分表示散射积分——对所有入射方向的贡献求和，由BRDF和余弦投影加权。

### 能量守恒与测量方程

能量守恒约束BRDF。方向-半球反射率（反照率）必须满足：

$\rho(\omega_o) = \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) \cos \theta_i d\omega_i \le 1$ 对所有 $\omega_o$

对于无损表面，等号成立。白炉测试验证能量守恒：在均匀照明环境中（$L_i = L_0$），封闭表面既不应获得也不应失去能量。

**详细能量平衡**：对于表面元素 $dA$，守恒要求：

$\int_\Omega L_o(\mathbf{x}, \omega) \cos \theta d\omega dA = L_e dA + \int_\Omega L_i(\mathbf{x}, \omega) \cos \theta d\omega dA$

在热平衡的封闭系统中，基尔霍夫定律将发射率与吸收率关联：
$\varepsilon(\lambda, \theta) = \alpha(\lambda, \theta) = 1 - \rho(\lambda, \theta)$

**测量方程**：测量方程将场景辐射率与传感器响应联系起来：

$I_j = \int_A \int_\Omega W_j(\mathbf{x}, \omega) L(\mathbf{x}, \omega) \cos \theta d\omega dA$

其中 $W_j(\mathbf{x}, \omega)$ 是像素 $j$ 的重要性（灵敏度）函数。辐射率和重要性之间的这种对偶性使双向算法成为可能。

**重要性传输**：重要性满足伴随方程：

$W(\mathbf{x}, \omega) = W_e(\mathbf{x}, \omega) + \int_\Omega f_r(\mathbf{x}, \omega, \omega') W(\mathbf{x}, \omega') \cos \theta' d\omega'$

这种对称性导致：
- 双向路径追踪
- 光子映射（正向光，反向重要性）
- 梯度计算的伴随方法

对于针孔相机，像素 $j$ 从针孔张成立体角 $\Omega_j$：

$I_j = \int_{\Omega_j} L(\mathbf{x}_{\text{lens}}, \omega) \cos^4 \theta d\omega$

$\cos^4 \theta$ 项考虑了：
- $\cos \theta$：投影透镜面积
- $\cos^3 \theta$：平方反比衰减和像素投影

**有限孔径相机**：对于具有孔径 $A_{\text{lens}}$ 的真实相机：

$I_j = \frac{1}{A_{\text{lens}}} \int_{A_{\text{lens}}} \int_{A_{\text{pixel}}} L(\mathbf{x}_{\text{lens}} \to \mathbf{x}_{\text{pixel}}) G(\mathbf{x}_{\text{lens}} \leftrightarrow \mathbf{x}_{\text{pixel}}) dA_{\text{pixel}} dA_{\text{lens}}$

这导致景深效果并需要仔细的采样策略。

### 算子形式与诺伊曼级数

渲染方程允许优雅的算子表述。定义传输算子 $\mathcal{T}$：

$(\mathcal{T}L)(\mathbf{x}, \omega) = \int_\Omega f_r(\mathbf{x}, \omega', \omega) L(\mathbf{x}, \omega') (\omega' \cdot \mathbf{n}) d\omega'$

则渲染方程变为：

$L = L_e + \mathcal{T}L$

这是第二类弗雷德霍姆方程。通过诺伊曼级数的解：

$L = (I - \mathcal{T})^{-1}L_e = \sum_{k=0}^\infty \mathcal{T}^k L_e = L_e + \mathcal{T}L_e + \mathcal{T}^2L_e + \dots$

每项都有物理意义：
- $L_e$：直接照明（仅发射）
- $\mathcal{T}L_e$：单次反弹照明
- $\mathcal{T}^2L_e$：两次反弹照明
- $\mathcal{T}^k L_e$：k次反弹照明

当 $||\mathcal{T}|| < 1$ 时级数收敛，这在最大反照率 $< 1$ 时发生。这种分解自然导致采样路径长度递增的路径追踪算法。

### 三点形式与几何耦合

渲染方程可以改写为三点形式，使几何耦合明确：

$L(\mathbf{x} \to \mathbf{x}') = L_e(\mathbf{x} \to \mathbf{x}') + \int_M f_r(\mathbf{x}'' \to \mathbf{x} \to \mathbf{x}') L(\mathbf{x}'' \to \mathbf{x}) G(\mathbf{x}'' \leftrightarrow \mathbf{x}) dA(\mathbf{x}'')$

其中几何因子是：

$G(\mathbf{x} \leftrightarrow \mathbf{x}') = V(\mathbf{x} \leftrightarrow \mathbf{x}') \frac{\cos \theta \cos \theta'}{||\mathbf{x} - \mathbf{x}'||^2}$

其中：
- $V(\mathbf{x} \leftrightarrow \mathbf{x}')$：二元可见性函数（相互可见为1，否则为0）
- $\cos \theta, \cos \theta'$：表面法线与连接线之间的角度
- $||\mathbf{x} - \mathbf{x}'||^2$：平方反比衰减的平方距离

这种形式强调光传输耦合所有表面点，导致路径积分表述。

**可见性复杂性**：可见性函数 $V(\mathbf{x} \leftrightarrow \mathbf{x}')$ 使渲染方程成为非线性和非局部的：
- 不连续：创建硬阴影和遮挡边界
- 计算昂贵：需要光线-场景相交
- 耦合所有几何：任何地方的变化都会影响所有地方的可见性

**核函数性质**：传输核 $K(\mathbf{x}'' \to \mathbf{x}) = f_r G V$ 具有重要性质：
- 沿 $\mathbf{x} = \mathbf{x}''$ 奇异（需要仔细正则化）
- 在遮挡边界处不连续
- 满足互易性：$K(\mathbf{x} \to \mathbf{x}') = K(\mathbf{x}' \to \mathbf{x})$

**与热方程的联系**：没有可见性时，渲染方程类似于具有非局部核的热方程。这个类比有助于理解：
- 多重散射的平滑性质
- 光学厚介质的扩散近似
- 有限元和多重网格求解方法

## 1.2 坐标系统与变换

### 世界、相机和物体空间

渲染管线涉及坐标系统层次结构，每个都针对特定计算进行优化：

1. **物体空间（模型空间）**：以规范形式定义的几何
   - 原点通常在物体中心或底部
   - 轴与自然对称性对齐
   - 简化建模和动画

2. **世界空间**：统一的场景坐标
   - 所有物体变换到公共框架
   - 光照和物理计算
   - 光线-物体相交

3. **相机空间（视图空间）**：以观察者为中心的坐标
   - 原点在眼点
   - -z轴沿视线方向（OpenGL约定）
   - +z进入屏幕（DirectX约定）
   - 简化投影和剔除

4. **裁剪空间**：投影后的齐次坐标
   - 透视除法前的4D坐标
   - 视锥体变为[-1,1]$^3$立方体（NDC）

5. **屏幕空间（光栅空间）**：最终的2D图像坐标
   - 整数像素坐标
   - 原点在左上角或左下角

### 齐次坐标与变换

齐次坐标统一了平移和线性变换。3D点 $\mathbf{p} = (x, y, z)$ 变为 $\mathbf{\tilde{p}} = (x, y, z, 1)$，而向量使用 $\mathbf{\tilde{v}} = (x, y, z, 0)$。

一般仿射变换矩阵：

$\mathbf{M} = \begin{bmatrix} \mathbf{A} & \mathbf{t} \\ \mathbf{0} & 1 \end{bmatrix}$

其中 $\mathbf{A}$ 是3$\times$3线性部分，$\mathbf{t}$ 是平移。常见变换：

**平移 ($t_x, t_y, t_z$)：**
$\begin{bmatrix} 1 & 0 & 0 & t_x \\ 0 & 1 & 0 & t_y \\ 0 & 0 & 1 & t_z \\ 0 & 0 & 0 & 1 \end{bmatrix}$

**绕轴 $\mathbf{a}$ 旋转角度 $\theta$：**
$\mathbf{R} = \cos \theta \mathbf{I} + (1 - \cos \theta) \mathbf{a}\mathbf{a}^T + \sin \theta [\mathbf{a}]_\times$

其中 $[\mathbf{a}]_\times$ 是反对称叉积矩阵。

**缩放 ($s_x, s_y, s_z$)：**
$\begin{bmatrix} s_x & 0 & 0 & 0 \\ 0 & s_y & 0 & 0 \\ 0 & 0 & s_z & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}$

### 法线和切线变换

法线必须变换以保持与表面垂直。给定点的变换 $\mathbf{M}$：

$\mathbf{n}' = \frac{(\mathbf{M}^{-T})^{3\times3} \mathbf{n}}{||(\mathbf{M}^{-T})^{3\times3} \mathbf{n}||}$

证明：对于表面上的切线 $\mathbf{t}$，$\mathbf{n} \cdot \mathbf{t} = 0$。变换后：
$\mathbf{n}' \cdot \mathbf{t}' = (\mathbf{M}^{-T}\mathbf{n}) \cdot (\mathbf{M}\mathbf{t}) = \mathbf{n}^T \mathbf{M}^{-1} \mathbf{M} \mathbf{t} = \mathbf{n} \cdot \mathbf{t} = 0$

对于正交标架 $\{\mathbf{t}, \mathbf{b}, \mathbf{n}\}$：
- 正向：变换 $\mathbf{t}$ 和 $\mathbf{b}$，然后 $\mathbf{n} = \mathbf{t} \times \mathbf{b}$
- 或如上变换 $\mathbf{n}$，然后重建标架

**面积和体积元素**：在变换 $\mathbf{M}$ 下，微分元素缩放为：
- 长度：$dl' = ||\mathbf{M}\mathbf{v}|| dl$（对于方向 $\mathbf{v}$）
- 面积：$dA' = |\det(\mathbf{M})| ||(\mathbf{M}^{-T})\mathbf{n}|| dA$
- 体积：$dV' = |\det(\mathbf{M})| dV$

**非均匀缩放问题**：非均匀缩放破坏各向同性：
- 球体 $\to$ 椭球体
- 各向同性BRDF $\to$ 各向异性BRDF
- 需要注意基于物理的材质

**手性保持**：当 $\det(\mathbf{M}) < 0$ 时，变换翻转方向：
- 右手系 $\to$ 左手系坐标系统
- 法线方向必须翻转
- 对一致的正面/背面判定至关重要

### 球面和立体角参数化

球面坐标为方向提供自然参数化：

$\omega = (\sin \theta \cos \phi, \sin \theta \sin \phi, \cos \theta)$

其中：
- $\theta \in [0, \pi]$：从+z轴的极角
- $\phi \in [0, 2\pi]$：从+x轴的方位角

雅可比给出微分立体角：

$d\omega = \left|\frac{\partial(\omega_x, \omega_y)}{\partial(\theta, \phi)}\right| d\theta d\phi = \sin \theta d\theta d\phi$

半球总立体角：$\int_\Omega d\omega = 2\pi$

用于采样的替代参数化：

**同心圆盘映射**（Shirley-Chiu）：
$(u, v) \in [-1, 1]^2 \to (r, \phi) \to (x, y)$ 在单位圆盘上

**八面体映射**：
单位球面 $\to$ 八面体 $\to$ 单位正方形
比球面坐标更好地保持面积

### 积分中的变量替换

积分的一般变量替换公式：

$\int_\Omega f(\mathbf{x}) d\mathbf{x} = \int_{\Omega'} f(\mathbf{x}(\mathbf{u})) \left|\det\left(\frac{\partial\mathbf{x}}{\partial\mathbf{u}}\right)\right| d\mathbf{u}$

对渲染至关重要：

**立体角到面积**：
$\int_\Omega L(\mathbf{x}, \omega) \cos \theta d\omega = \int_A L(\mathbf{x}, \omega(\mathbf{x}')) G(\mathbf{x} \leftrightarrow \mathbf{x}') dA'$

其中 $G(\mathbf{x} \leftrightarrow \mathbf{x}') = V(\mathbf{x} \leftrightarrow \mathbf{x}') \frac{\cos \theta \cos \theta'}{||\mathbf{x} - \mathbf{x}'||^2}$

**半球到圆盘**（用于余弦加权采样）：
映射 $(\theta, \phi) \to (r, \phi)$ 其中 $r = \sin \theta$
则 $p(r, \phi) = p(\theta, \phi) \left|\frac{\partial(\theta, \phi)}{\partial(r, \phi)}\right| = \frac{p(\theta, \phi)}{\cos \theta}$

**测度论基础**：变量替换公式具有测度论基础：
- 推前测度：$\mu'(A) = \mu(f^{-1}(A))$
- Radon-Nikodym导数给出雅可比
- 对理解蒙特卡洛收敛性至关重要

### 投影变换与透视

透视投影矩阵将视锥体映射到裁剪空间：

$\mathbf{P} = \begin{bmatrix}
n/r & 0 & 0 & 0 \\
0 & n/t & 0 & 0 \\
0 & 0 & -(f+n)/(f-n) & -2fn/(f-n) \\
0 & 0 & -1 & 0
\end{bmatrix}$

其中 $n, f$ 是近/远平面，$r, t$ 是近平面的右/顶。

透视除法后：
- $x_{ndc} = x_{clip} / w_{clip} \in [-1, 1]$
- $y_{ndc} = y_{clip} / w_{clip} \in [-1, 1]$
- $z_{ndc} = z_{clip} / w_{clip} \in [-1, 1]$

重要性质：
- 直线保持为直线（除了通过眼点的）
- 平面保持为平面
- 深度精度是非线性的（近处比远处精度高）

### 重心坐标与插值

对于具有顶点 $\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2$ 的三角形，重心坐标 $(u, v, w)$ 满足：

$\mathbf{p} = u\mathbf{v}_0 + v\mathbf{v}_1 + w\mathbf{v}_2$

约束条件 $u + v + w = 1$。通过面积计算：

$u = \frac{\text{Area}(\mathbf{p}, \mathbf{v}_1, \mathbf{v}_2)}{\text{Area}(\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2)}$

性质：
- $u, v, w \in [0, 1]$ 当且仅当 $\mathbf{p}$ 在三角形内
- 线性插值：$f(\mathbf{p}) = uf_0 + vf_1 + wf_2$
- 透视正确的插值需要 $1/z$ 校正

透视正确的属性插值：
1. 在屏幕空间插值 $a/z, b/z, c/z$ 和 $1/z$
2. 恢复属性：$a = (a/z)/(1/z)$

### 微分几何与局部标架

在每个表面点，我们构建用于着色计算的局部标架：

**切空间基**：
- $\mathbf{n}$：表面法线（$\frac{\partial\mathbf{p}}{\partial u} \times \frac{\partial\mathbf{p}}{\partial v}$ 归一化）
- $\mathbf{t}$：切线（通常 $\frac{\partial\mathbf{p}}{\partial u}$ 归一化）
- $\mathbf{b}$：副切线（$\mathbf{n} \times \mathbf{t}$）

**第一基本形式**：度量张量描述局部表面几何：
$\mathbf{I} = \begin{bmatrix} E & F \\ F & G \end{bmatrix}$

其中：
- $E = \frac{\partial\mathbf{p}}{\partial u} \cdot \frac{\partial\mathbf{p}}{\partial u}$
- $F = \frac{\partial\mathbf{p}}{\partial u} \cdot \frac{\partial\mathbf{p}}{\partial v}$
- $G = \frac{\partial\mathbf{p}}{\partial v} \cdot \frac{\partial\mathbf{p}}{\partial v}$

弧长：$ds^2 = E du^2 + 2F du dv + G dv^2$
面积元素：$dA = \sqrt{EG - F^2} du dv$

**第二基本形式**：描述表面曲率：
$\mathbf{II} = \begin{bmatrix} e & f \\ f & g \end{bmatrix}$

其中 $e = \mathbf{n} \cdot \frac{\partial^2\mathbf{p}}{\partial u^2}$ 等。

主曲率 $\kappa_1, \kappa_2$ 是 $\mathbf{II}\mathbf{I}^{-1}$ 的特征值
- 平均曲率：$H = (\kappa_1 + \kappa_2)/2$
- 高斯曲率：$K = \kappa_1\kappa_2$

**世界空间的变换**：
$\begin{bmatrix} \mathbf{t}_{\text{world}} \\ \mathbf{b}_{\text{world}} \\ \mathbf{n}_{\text{world}} \end{bmatrix} = \begin{bmatrix} t_x & t_y & t_z \\ b_x & b_y & b_z \\ n_x & n_y & n_z \end{bmatrix} \begin{bmatrix} \mathbf{t}_{\text{local}} \\ \mathbf{b}_{\text{local}} \\ \mathbf{n}_{\text{local}} \end{bmatrix}$

这个正交矩阵可以通过转置求逆。

**各向异性BRDF参数化**：
许多BRDF依赖于相对于切线的角度：
- $\phi_h$：半向量在切空间中的方位角
- 能够建模拉丝金属、织物、毛发

**平行传输**：在表面上追踪光线时，切标架必须平行传输：
- 保持方向一致性
- 保留各向异性外观
- 与光学中的几何相位相关

## 1.3 BRDF、BSDF和BSSRDF

### 双向反射分布函数（BRDF）

BRDF $f_r$ 量化入射辐照度和反射辐射率之间的微分关系：

$f_r(\mathbf{x}, \omega_i, \omega_o) = \frac{dL_o(\mathbf{x}, \omega_o)}{dE_i(\mathbf{x}, \omega_i)} = \frac{dL_o(\mathbf{x}, \omega_o)}{L_i(\mathbf{x}, \omega_i) \cos \theta_i d\omega_i} \text{ [sr}^{-1}\text{]}$

物理上，它表示来自方向 $\omega_i$ 的光子散射到方向 $\omega_o$ 的概率密度（归一化后）。

BRDF可以分解为分量：
$f_r = f_d + f_s + f_g + \dots$

其中 $f_d$ 是漫反射，$f_s$ 是镜面反射，$f_g$ 是光泽反射等。这种分解有助于重要性采样。

### BRDF基本性质

**亥姆霍兹互易性：**
$f_r(\mathbf{x}, \omega_i, \omega_o) = f_r(\mathbf{x}, \omega_o, \omega_i)$

这源自麦克斯韦方程的时间反演对称性和细致平衡原理。它使双向路径追踪和光子映射成为可能。

**能量守恒：**
方向-半球反射率必须满足：

$\rho(\omega_i) = \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) \cos \theta_o d\omega_o \le 1$ 对所有 $\omega_i$

对于能量守恒的BRDF，当吸收为零时等号成立。半球-半球反射率：

$\rho_{hh} = \frac{1}{\pi} \int_\Omega \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) \cos \theta_i \cos \theta_o d\omega_i d\omega_o \le 1$

**非负性：**
$f_r(\mathbf{x}, \omega_i, \omega_o) \ge 0$

负值将意味着依赖于出射方向的能量吸收，违反因果性。

**可测性和可积性：**
对于蒙特卡洛积分收敛：
$f_r \in L^2(\Omega \times \Omega)$ （平方可积）

### 经典BRDF模型

**朗伯（完全漫反射）：**
$f_r = \frac{\rho_d}{\pi}$

其中 $\rho_d \in [0, 1]$ 是漫反射率。构造上能量守恒。

**Phong模型：**
$f_r = \frac{\rho_d}{\pi} + \rho_s \frac{n+2}{2\pi} (\mathbf{r} \cdot \omega_o)^n$

其中 $\mathbf{r} = 2(\mathbf{n} \cdot \omega_i)\mathbf{n} - \omega_i$ 是反射方向。不满足互易性！

**Blinn-Phong（互易）：**
$f_r = \frac{\rho_d}{\pi} + \rho_s \frac{n+2}{8\pi} \frac{(\mathbf{n} \cdot \mathbf{h})^n}{\max(\cos \theta_i, \cos \theta_o)}$

其中 $\mathbf{h} = \frac{\omega_i + \omega_o}{||\omega_i + \omega_o||}$ 是半向量。

**Cook-Torrance微面元模型：**
$f_r = \frac{\rho_d}{\pi} + \frac{D(\mathbf{h})G(\omega_i, \omega_o)F(\omega_i, \mathbf{h})}{4 \cos \theta_i \cos \theta_o}$

其中：
- $D(\mathbf{h})$：法线分布函数（如GGX）
- $G(\omega_i, \omega_o)$：几何衰减（遮蔽/阴影）
- $F(\omega_i, \mathbf{h})$：菲涅尔反射率

### 扩展到BSDF

双向散射分布函数（BSDF）统一了反射和透射：

$f_s(\mathbf{x}, \omega_i, \omega_o) = \begin{cases}
f_r(\mathbf{x}, \omega_i, \omega_o) & \text{如果 } \omega_i \cdot \mathbf{n} \text{ 和 } \omega_o \cdot \mathbf{n} \text{ 同号} \\
f_t(\mathbf{x}, \omega_i, \omega_o) & \text{如果 } \omega_i \cdot \mathbf{n} \text{ 和 } \omega_o \cdot \mathbf{n} \text{ 异号}
\end{cases}$

对于电介质界面（如玻璃），斯涅尔定律控制折射：
$n_i \sin \theta_i = n_o \sin \theta_o$

菲涅尔方程确定反射/透射概率：
$F_r = \left(\frac{n_i \cos \theta_i - n_o \cos \theta_o}{n_i \cos \theta_i + n_o \cos \theta_o}\right)^2$ （s偏振）

**BTDF的广义互易性：**
由于界面上的辐射率压缩/扩展：

$n_i^2 f_t(\mathbf{x}, \omega_i, \omega_o) = n_o^2 f_t(\mathbf{x}, \omega_o, \omega_i)$

这解释了辐射率 $L/n^2$ 不变中的 $n^2$ 因子。

### 次表面散射的BSSRDF

双向散射表面反射分布函数将BRDF推广到非局部传输：

$S(\mathbf{x}_i, \omega_i, \mathbf{x}_o, \omega_o) = \frac{dL_o(\mathbf{x}_o, \omega_o)}{d\Phi_i(\mathbf{x}_i, \omega_i)} \text{ [m}^{-2}\text{sr}^{-1}\text{]}$

与BRDF的关键区别：
- 耦合不同的表面点
- 单位包含逆面积
- 不再是纯材质属性（依赖于几何）

带BSSRDF的渲染方程：

$L_o(\mathbf{x}_o, \omega_o) = L_e(\mathbf{x}_o, \omega_o) + \int_A \int_\Omega S(\mathbf{x}_i, \omega_i, \mathbf{x}_o, \omega_o) L_i(\mathbf{x}_i, \omega_i) \cos \theta_i d\omega_i dA_i$

**扩散近似：**
对于高散射介质，BSSRDF可以近似为：

$S(\mathbf{x}_i, \omega_i, \mathbf{x}_o, \omega_o) \approx \frac{1}{\pi}F_t(\omega_i)R(||\mathbf{x}_i - \mathbf{x}_o||)F_t(\omega_o)$

其中 $R(r)$ 是扩散轮廓，$F_t$ 是菲涅尔透射率。

### 数学约束与物理合理性

物理有效的BRDF必须满足：

1. **互易性**：$f_r(\mathbf{x}, \omega_i, \omega_o) = f_r(\mathbf{x}, \omega_o, \omega_i)$
   - 测试：交换光源和相机渲染场景

2. **能量守恒**：$\forall \omega_i: \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) \cos \theta_o d\omega_o \le 1$
   - 测试：白炉测试（均匀照明）

3. **非负性**：$f_r(\mathbf{x}, \omega_i, \omega_o) \ge 0$
   - 违反会导致能量吸收异常

4. **平滑性**：$f_r$ 应该是 $C^0$ 连续的（$C^1$ 更好）
   - 不连续性导致采样困难

5. **菲涅尔行为**：对于光滑表面，当 $\theta \to \pi/2$ 时 $f_r \to 1$
   - 所有表面在掠射角度都变成镜面

### 各向异性BRDF

对于具有方向结构的材料（拉丝金属、织物、毛发），BRDF依赖于方位角：

$f_r(\mathbf{x}, \omega_i, \omega_o, \phi)$

其中 $\phi$ 是半向量投影与切线方向之间的角度。

**Ward各向异性模型：**
$f_r = \frac{\rho_d}{\pi} + \frac{\rho_s \exp\left(-\tan^2\theta_h\left(\frac{\cos^2\phi}{\alpha_x^2} + \frac{\sin^2\phi}{\alpha_y^2}\right)\right)}{4\pi \alpha_x \alpha_y \sqrt{\cos \theta_i \cos \theta_o}}$

其中 $\alpha_x, \alpha_y$ 控制各向异性粗糙度。

### 空间变化BRDF（SVBRDF）

真实材料表现出空间变化：
$f_r(\mathbf{x}, \omega_i, \omega_o) = f_r(u, v, \omega_i, \omega_o)$

其中 $(u, v)$ 是纹理坐标。这使得以下成为可能：
- 材质属性的纹理映射
- 测量的BRDF数据（BTF - 双向纹理函数）
- 过程化材质变化

## 1.4 渲染中的蒙特卡洛积分

### 期望值与方差

蒙特卡洛积分使用随机采样估计积分：

$I = \int_\Omega f(\mathbf{x}) d\mathbf{x} \approx \frac{1}{N} \sum_{i=1}^N \frac{f(\mathbf{X}_i)}{p(\mathbf{X}_i)}$

其中 $\mathbf{X}_i \sim p(\mathbf{x})$ 是来自概率密度 $p$ 的样本。

估计器是无偏的：$E[\hat{I}] = I$

方差是：
$\text{Var}[\hat{I}] = \frac{1}{N} \int_\Omega \left(\frac{f(\mathbf{x})}{p(\mathbf{x})} - I\right)^2 p(\mathbf{x}) d\mathbf{x}$

### 重要性采样

最优采样通过匹配 $p$ 和 $|f|$ 来最小化方差：

$p^*(\mathbf{x}) = \frac{|f(\mathbf{x})|}{\int_\Omega |f(\mathbf{x})| d\mathbf{x}}$

对于渲染方程，好的采样策略包括：
- BRDF采样：$p(\omega) \propto f_r(\omega_i, \omega_o)$
- 光源采样：$p(\omega) \propto L_e$
- 余弦采样：$p(\omega) \propto \cos \theta$

### 多重重要性采样（MIS）

当有多个采样策略可用时，MIS将它们最优地组合：

$\hat{I} = \sum_{i=1}^{n_f} w_f(\mathbf{X}_{f,i}) \frac{f(\mathbf{X}_{f,i})}{p_f(\mathbf{X}_{f,i})} + \sum_{j=1}^{n_g} w_g(\mathbf{X}_{g,j}) \frac{f(\mathbf{X}_{g,j})}{p_g(\mathbf{X}_{g,j})}$

平衡启发式提供好的权重：
$w_f(\mathbf{x}) = \frac{n_f p_f(\mathbf{x})}{n_f p_f(\mathbf{x}) + n_g p_g(\mathbf{x})}$

### 俄罗斯轮盘赌

为了用有限计算创建无偏估计器，俄罗斯轮盘赌随机终止路径：

$L'_i = \begin{cases}
L_i/q & \text{以概率 } q \\
0 & \text{以概率 } 1-q
\end{cases}$

这保持 $E[L'_i] = L_i$ 同时限制计算量。

### 收敛速率与误差界限

蒙特卡洛收敛遵循中心极限定理：

$P(|\hat{I} - I| \le \varepsilon) \approx 2\Phi(\varepsilon\sqrt{N}/\sigma) - 1$

其中 $\Phi$ 是正态分布CDF，$\sigma^2$ 是方差。误差以 $O(1/\sqrt{N})$ 递减，与维度无关——这对高维光传输至关重要。

## 1.5 路径积分表述

### 光传输作为路径积分

我们可以将渲染方程重新表述为对所有可能光路径的积分。长度为k的路径是：

$\bar{\mathbf{x}} = \mathbf{x}_0\mathbf{x}_1\dots\mathbf{x}_k$

其中 $\mathbf{x}_0$ 在光源上，$\mathbf{x}_k$ 在相机传感器上。

### 路径空间与测度

路径空间 $\bar{\Omega}_k$ 包含所有长度为k的有效路径。路径的测度是：

$d\mu(\bar{\mathbf{x}}) = dA(\mathbf{x}_0) \prod_{i=1}^{k} dA(\mathbf{x}_i)$

路径的贡献是：

$f(\bar{\mathbf{x}}) = L_e(\mathbf{x}_0 \to \mathbf{x}_1) \left(\prod_{i=1}^{k-1} f_s(\mathbf{x}_{i-1} \to \mathbf{x}_i \to \mathbf{x}_{i+1}) G(\mathbf{x}_i \leftrightarrow \mathbf{x}_{i+1})\right) W(\mathbf{x}_{k-1} \to \mathbf{x}_k)$

其中G是几何因子：

$G(\mathbf{x} \leftrightarrow \mathbf{x}') = V(\mathbf{x} \leftrightarrow \mathbf{x}') \frac{\cos \theta \cos \theta'}{||\mathbf{x} - \mathbf{x}'||^2}$

### 与费曼路径积分的联系

路径积分表述类似于费曼对量子力学的方法：

$I = \sum_{k=2}^\infty \int_{\bar{\Omega}_k} f(\bar{\mathbf{x}}) d\mu(\bar{\mathbf{x}})$

这个对所有路径长度的无限和捕获了全局光照。每一项代表k-1次反弹的路径。

### 递归表述与诺伊曼级数

点上的值满足递归关系：

$L(\mathbf{x}, \omega) = L_e(\mathbf{x}, \omega) + \int_M f_s(\mathbf{y} \to \mathbf{x} \to \omega) L(\mathbf{y}, \mathbf{x} - \mathbf{y}) G(\mathbf{y} \leftrightarrow \mathbf{x}) dA(\mathbf{y})$

这导致诺伊曼级数解：

$L = L^{(0)} + L^{(1)} + L^{(2)} + \dots$

其中 $L^{(k)}$ 代表k次反弹照明。

### 体积渲染方程预览

路径积分自然扩展到参与介质。对于具有吸收 $\sigma_a$、散射 $\sigma_s$ 和相位函数 $p$ 的体积：

$L(\mathbf{x}, \omega) = \int_0^\infty T(0,t) \left[\sigma_a L_e + \sigma_s \int_{S^2} p(\omega', \omega) L(\mathbf{x}+t\omega, \omega') d\omega'\right] dt + T(0,\infty) L_\infty$

其中 $T(s,t) = \exp\left(-\int_s^t \sigma_t(\mathbf{x}+u\omega) du\right)$ 是透射率。

这个统一的表述将在后续章节中连接所有渲染方法。

## 本章小结

本章通过几何光学建立了计算机图形学的数学基础：

1. **渲染方程** $L_o = L_e + \int f_r L_i \cos \theta d\omega$ 控制光传输
2. **坐标变换** 在正确应用时保持辐射度量
3. **BRDF** 必须满足互易性、能量守恒和非负性
4. **蒙特卡洛方法** 以 $O(1/\sqrt{N})$ 收敛率求解高维积分
5. **路径积分** 将光传输统一为对所有可能路径的积分

这些概念构成了所有渲染算法的基础。路径积分表述特别使我们能够将基于点的、基于图像的和神经渲染方法统一处理为解决同一基本方程的不同方法。

## 练习题

### 练习 1.1：沿光线的辐射率
证明在真空中沿光线辐射率保持恒定。从辐射率的定义开始，并使用平方反比定律。

**提示：** 考虑沿光线的两个微分面积 $dA_1$ 和 $dA_2$，并证明 $L_1 = L_2$。

<details>
<summary>解答</summary>

考虑沿光线距离 $r_1$ 和 $r_2$ 处的微分面积 $dA_1$ 和 $dA_2$。所张的立体角为：

$d\omega_1 = dA_2 \cos \theta_2 / r_{12}^2$
$d\omega_2 = dA_1 \cos \theta_1 / r_{12}^2$

离开 $dA_1$ 朝向 $dA_2$ 的通量为：
$d\Phi = L_1 dA_1 \cos \theta_1 d\omega_1 = L_1 dA_1 \cos \theta_1 dA_2 \cos \theta_2 / r_{12}^2$

根据能量守恒，这等于到达 $dA_2$ 的通量：
$d\Phi = L_2 dA_2 \cos \theta_2 d\omega_2 = L_2 dA_2 \cos \theta_2 dA_1 \cos \theta_1 / r_{12}^2$

因此 $L_1 = L_2$，证明了辐射率不变性。
</details>

### 练习 1.2: BRDF能量守恒
证明朗伯BRDF $f_r = \rho/\pi$ 对于任何反照率 $\rho \le 1$ 都满足能量守恒。

**提示：** 使用球面坐标在半球上积分。

<details>
<summary>解答</summary>

对于朗伯BRDF $f_r = \rho/\pi$，方向-半球反射率为：

$\int_\Omega f_r \cos \theta_o d\omega_o = (\rho/\pi) \int_0^{2\pi} \int_0^{\pi/2} \cos \theta \sin \theta d\theta d\phi$

$= (\rho/\pi) \cdot 2\pi \cdot \int_0^{\pi/2} \cos \theta \sin \theta d\theta$

$= (\rho/\pi) \cdot 2\pi \cdot [\sin^2 \theta/2]_0^{\pi/2}$

$= (\rho/\pi) \cdot 2\pi \cdot (1/2) = \rho$

由于 $\rho \le 1$，能量守恒得到满足。
</details>

### 练习 1.3: 蒙特卡洛方差
推导估计 $\int_0^1 \sqrt{x} dx$ 的最优重要性采样分布，并计算由此产生的方差。

**提示：** 最优PDF与 $|f(x)|$ 成比例。

<details>
<summary>解答</summary>

对于 $[0,1]$ 上的 $f(x) = \sqrt{x}$，最优PDF为：

$p^*(x) = \sqrt{x} / \int_0^1 \sqrt{x} dx = \sqrt{x} / (2/3) = (3/2)\sqrt{x}$

要采样：$X = F^{-1}(U)$ 其中 $F(x) = \int_0^x (3/2)\sqrt{t} dt = x^{3/2}$

因此 $X = U^{2/3}$

使用这种采样，估计器变为：
$f(X)/p^*(X) = \sqrt{X} / ((3/2)\sqrt{X}) = 2/3$

由于估计器是常数，方差为0——完美的重要性采样消除了方差。
</details>

### 练习 1.4: 路径积分收敛性（挑战）
证明当场景中的最大反照率小于1时，渲染方程的诺伊曼级数收敛。

**提示：** 使用算子范数 $||\mathcal{T}|| \le \rho_{\text{max}} < 1$。

<details>
<summary>解答</summary>

传输算子 $\mathcal{T}$ 满足：

$||(\mathcal{T}L)||_\infty \le \rho_{\text{max}} ||L||_\infty$

其中 $\rho_{\text{max}} = \max_{x,\omega} \int_\Omega f_r(x,\omega',\omega) \cos \theta' d\omega' < 1$

通过归纳法：$||\mathcal{T}^k L_e||_\infty \le \rho_{\text{max}}^k ||L_e||_\infty$

诺伊曼级数收敛：
$||\sum_{k=0}^\infty \mathcal{T}^k L_e||_\infty \le \sum_{k=0}^\infty \rho_{\text{max}}^k ||L_e||_\infty = ||L_e||_\infty/(1-\rho_{\text{max}})$

这证明了当 $\rho_{\text{max}} < 1$ 时的收敛性。
</details>

### 练习 1.5: 坐标变换雅可比
推导将渲染方程从立体角积分转换为表面面积积分的雅可比。

**提示：** 从 $d\omega = dA \cos \theta'/r^2$ 开始。

<details>
<summary>解答</summary>

给定点 $\mathbf{x}$ 和 $\mathbf{x}'$ 以及连接向量 $\mathbf{r} = \mathbf{x}' - \mathbf{x}$：

$d\omega = dA' \cos \theta' / ||\mathbf{r}||^2$

渲染方程变为：

$L_o(\mathbf{x},\omega_o) = L_e(\mathbf{x},\omega_o) + \int_M f_r(\mathbf{x},\omega(\mathbf{x}'),\omega_o) L_i(\mathbf{x}',-\omega(\mathbf{x}')) V(\mathbf{x}\leftrightarrow\mathbf{x}') (\cos \theta \cos \theta')/||\mathbf{r}||^2 dA'$

其中：
- $\omega(\mathbf{x}') = (\mathbf{x}' - \mathbf{x})/||\mathbf{x}' - \mathbf{x}||$
- $\cos \theta = \mathbf{n}(\mathbf{x}) \cdot \omega(\mathbf{x}')$
- $\cos \theta' = -\mathbf{n}(\mathbf{x}') \cdot \omega(\mathbf{x}')$
- $V(\mathbf{x}\leftrightarrow\mathbf{x}')$ 是可见性函数

雅可比是 $|\partial\omega/\partial\mathbf{x}'| = \cos \theta'/||\mathbf{r}||^2$。
</details>

### 练习 1.6: BSDF折射互易性（挑战）
推导折射率分别为 $n_1$ 和 $n_2$ 的介质之间BSDF的广义互易关系。

**提示：** 使用相空间辐射率 $\tilde{L} = L/n^2$。

<details>
<summary>解答</summary>

定义相空间辐射率 $\tilde{L} = L/n^2$。在折射率为 $n$ 的介质中：

$\tilde{L}$ 沿光线不变（广义辐射率定理）

在界面处，功率守恒要求：

$\tilde{L}_1(\mathbf{x},\omega_i) d\omega_i = \tilde{L}_2(\mathbf{x},\omega_t) d\omega_t$

使用斯涅尔定律：$n_1 \sin \theta_i = n_2 \sin \theta_t$

立体角比率为：$d\omega_t/d\omega_i = (n_1/n_2)^2 \cos \theta_t/\cos \theta_i$

对于BTDF：$f_t(\omega_i\to\omega_t) = dL_2/dE_1 = (n_2^2/n_1^2) d\tilde{L}_2/d\tilde{E}_1$

根据麦克斯韦方程的时间反演对称性：

$n_1^2 f_t(\mathbf{x},\omega_i\to\omega_t) = n_2^2 f_t(\mathbf{x},\omega_t\to\omega_i)$
</details>

### 练习 1.7: 俄罗斯轮盘赌偏差
证明具有延续概率 $q$ 的俄罗斯轮盘赌创建了一个无偏估计器。

**提示：** 使用全期望定律。

<details>
<summary>解答</summary>

设 $L$ 为真实值，$L'$ 为俄罗斯轮盘赌估计器：

$L' = \begin{cases}
L/q & \text{以概率 } q \\
0 & \text{以概率 } 1-q
\end{cases}$

期望值为：

$E[L'] = q \cdot (L/q) + (1-q) \cdot 0 = L$

因此 $E[L'] = L$，证明了估计器是无偏的。

方差为：
$\text{Var}[L'] = E[L'^2] - (E[L'])^2 = q(L/q)^2 - L^2 = L^2/q - L^2 = L^2(1-q)/q$

请注意，方差随 $q$ 的减小而增加，这表明了偏差-方差权衡。
</details>

### 练习 1.8: 体积渲染收敛性（开放问题）
考虑具有空间变化消光的体积渲染方程。在什么条件下路径积分公式收敛？讨论光学厚度与收敛速率之间的关系。

**提示：** 考虑沿任何光线的最大光学厚度 $\tau_{\text{max}}$。

<details>
<summary>解答</summary>

对于具有有界消光 $\sigma_t \le \sigma_{\text{max}}$ 和反照率 $\alpha = \sigma_s/\sigma_t \le \alpha_{\text{max}} < 1$ 的体积：

第k次散射项有界：
$||L^{(k)}||_\infty \le ||L_e||_\infty \alpha_{\text{max}}^k$

这确保了几何收敛。

对于沿光线的光学厚度 $\tau = \int \sigma_t ds$：
- 低 $\tau$：单次散射占主导，收敛快
- 高 $\tau$：多次散射重要，收敛慢
- $\tau \to \infty$：扩散机制，可能需要专门方法

开放问题：
1. 异构介质的最优重要性采样
2. 收敛速率与 $\sigma_t$ 的空间频率关系
3. 多次散射中与传输平均自由程的联系
</details>

## 常见陷阱和错误（Gotchas）

1. **辐射率与辐照度混淆**
   - 辐射率单位为 W/(m$^2$·sr)，辐照度单位为 W/m$^2$
   - 始终检查方程中的量纲一致性
   - 记住：相机测量辐照度（积分辐射率）

2. **法线变换不正确**
   - 法线通过逆转置矩阵变换，而不是直接矩阵
   - 变换后未能重新归一化
   - 在坐标系之间变换时的符号错误

3. **BRDF能量守恒违规**
   - 忘记积分中的余弦项
   - 未考虑掠射角处的菲涅尔效应
   - 分析BRDF模型归一化不正确

4. **PDF错误导致的蒙特卡洛偏差**
   - PDF必须在被积函数非零的任何地方都大于0
   - 重要性采样中的归一化错误
   - 改变变量时忘记雅可比

5. **数值精度问题**
   - 光线-球体交点可能因灾难性抵消而失败
   - 几何因子 $G(\mathbf{x}\leftrightarrow\mathbf{x}')$ 中的小分母
   - 长光线路径中累积的误差

## 最佳实践清单

### 设计评审

- [ ] **验证物理单位**：检查所有方程的量纲一致性
- [ ] **能量守恒**：确保BRDF满足 $\int f_r \cos \theta d\omega \le 1$
- [ ] **互易性**：验证 $f_r(\omega_i, \omega_o) = f_r(\omega_o, \omega_i)$
- [ ] **坐标系一致性**：记录并验证所有坐标约定
- [ ] **采样策略**：将重要性采样与主要贡献匹配

### 实现评审

- [ ] **鲁棒的光线-图元交点**：使用数值稳定的算法
- [ ] **PDF验证**：断言 $p(x) > 0$ 且 $\int p(x) dx = 1$
- [ ] **方差减少**：为多种采样策略实现MIS
- [ ] **浮点卫生**：检查NaN/Inf传播
- [ ] **终止标准**：使用俄罗斯轮盘赌进行无限递归

### 调试清单

- [ ] **白炉测试**：均匀发射应产生均匀辐射率
- [ ] **互易性测试**：交换光源和相机，验证结果相同
- [ ] **能量审计**：跟踪总输入通量 = 总输出通量
- [ ] **收敛分析**：绘制方差与样本数的关系
- [ ] **参考比较**：在可用时与解析解进行验证


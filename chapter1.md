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

在几何光学中，我们使用光线来建模光的传播——在均匀介质中沿直线传播的无限细光束。当波长 λ << 特征尺寸时，这种近似成立，使我们可以忽略衍射和干涉。光线参数化为：

**r**(t) = **o** + t**d**

其中 **o** ∈ ℝ³ 是原点，**d** ∈ ℝ³ 是方向（||**d**|| = 1），t ≥ 0 是沿光线的参数。

光线方程源自程函方程 ∇S = n**k̂** 在 λ → 0 极限下的结果，其中 S 是相位，n 是折射率。在非均匀介质中，光线遵循满足以下方程的曲线路径：

$\frac{d}{ds}\left(n \frac{d\mathbf{r}}{ds}\right) = \nabla n$

当 n 为常数时，这简化为直线。

**与波动光学的联系**：几何光学近似源自波动方程的WKB（Wentzel-Kramers-Brillouin）近似。当我们将 ψ = A exp(ikS) 代入亥姆霍兹方程并取 k → ∞ 时：

$\nabla^2\psi + k^2n^2\psi = 0 \rightarrow (\nabla S)^2 = n^2$ （程函方程）

等相位面 S = const 是波前，光线是这些波前的正交轨迹。当我们在后续章节扩展到波动光学时，这种联系变得至关重要。

**光线光学的有效性**：几何光学近似在以下情况下失效：
1. 特征尺寸 ~ 波长（衍射变得显著）
2. 靠近焦散面，光线密度 → ∞
3. 存在尖锐边缘或不连续性
4. 需要相位信息的相干现象

**费马原理**：光线遵循光程稳定的路径：

$\delta\int n(\mathbf{r}) ds = 0$

这个变分原理统一了光线行为：在均匀介质中的直线、界面处的斯涅尔定律以及梯度折射率介质中的曲线路径。它还与物理学中的最小作用量原理和微分几何中的测地线方程相关联。

### 辐射度量

在推导渲染方程之前，我们必须建立辐射度量框架。这些量形成一个层次结构，每个都建立在前一个的基础上：

**辐射能** Q 测量总电磁能量：
Q [J]

**辐射通量（功率）** Φ 测量单位时间的能量：
$\Phi = \frac{dQ}{dt}$ [W]

**辐射强度** I 测量点源单位立体角的通量：
$I = \frac{d\Phi}{d\omega}$ [W/sr]

**辐照度** E 测量入射单位面积的通量：
$E = \frac{d\Phi}{dA}$ [W/m²]

**辐射出射度** M 测量离开单位面积的通量：
$M = \frac{d\Phi}{dA}$ [W/m²]

**辐射率** L 测量单位面积单位立体角的通量：
$L = \frac{d^2\Phi}{dA \cos \theta d\omega} = \frac{d^2\Phi}{dA_\perp d\omega}$ [W/(m²·sr)]

其中 $dA_\perp = dA \cos \theta$ 是垂直于光线方向的投影面积。

辐射率是渲染中的基本量，因为：
1. 它在真空中沿光线保持恒定（辐射率不变性定理）
2. 它是相机和眼睛测量的量
3. 所有其他辐射度量都可以从它导出

**光度量与辐射度量**：虽然我们关注辐射度量（物理能量），但为人类感知进行渲染时通常使用光度量：
- 光通量 [lm] = 辐射通量 [W] × 光效能
- 亮度 [cd/m²] = 辐射率 × 明视觉响应 V(λ)
- CIE光度函数 V(λ) 在555 nm（绿色）处达到峰值

**光谱辐射率**：实际上，辐射率随波长变化：
$L(\mathbf{x}, \omega, \lambda)$ [W/(m²·sr·nm)]

对于渲染，我们通常使用：
- RGB近似：光谱的3个样本
- 光谱渲染：N个波长样本（通常10-100）
- 主波长：光谱的随机采样

**相干与非相干叠加**：辐射度量假设非相干光——强度直接相加。对于相干光源（激光），我们必须跟踪相位并叠加复振幅：
$I_{\text{total}} = |E_1 + E_2|^2 \neq |E_1|^2 + |E_2|^2$ （一般情况下）

### 渲染方程

渲染方程由Kajiya（1986）提出，描述场景中光的平衡分布。它源于功率平衡：在任何表面点，输出功率等于发射功率加反射功率。

在具有法线 **n** 的表面点 **x** 处，沿方向 $\omega_o$ 的出射辐射率 $L_o$ 满足：

$L_o(\mathbf{x}, \omega_o) = L_e(\mathbf{x}, \omega_o) + \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) L_i(\mathbf{x}, \omega_i) (\omega_i \cdot \mathbf{n}) d\omega_i$

其中：
- $L_e(\mathbf{x}, \omega_o)$ 是从 **x** 沿方向 $\omega_o$ 的发射辐射率
- $f_r(\mathbf{x}, \omega_i, \omega_o)$ 是BRDF [sr⁻¹]
- $L_i(\mathbf{x}, \omega_i)$ 是在 **x** 处来自方向 $\omega_i$ 的入射辐射率
- Ω 是 **x** 上方的半球（其中 $\omega \cdot \mathbf{n} > 0$）
- $(\omega_i \cdot \mathbf{n}) = \cos \theta_i$ 考虑了投影面积

积分表示散射积分——对所有入射方向的贡献求和，由BRDF和余弦投影加权。

### 能量守恒与测量方程

能量守恒约束BRDF。方向-半球反射率（反照率）必须满足：

$\rho(\omega_o) = \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) \cos \theta_i d\omega_i \leq 1$ 对所有 $\omega_o$

对于无损表面，等号成立。白炉测试验证能量守恒：在均匀照明环境中（$L_i = L_0$），封闭表面既不应获得也不应失去能量。

**详细能量平衡**：对于表面元素 dA，守恒要求：

$\int_\Omega L_o(\mathbf{x}, \omega) \cos \theta d\omega dA = L_e dA + \int_\Omega L_i(\mathbf{x}, \omega) \cos \theta d\omega dA$

在热平衡的封闭系统中，基尔霍夫定律将发射率与吸收率关联：
$\varepsilon(\lambda, \theta) = \alpha(\lambda, \theta) = 1 - \rho(\lambda, \theta)$

**测量方程**：测量方程将场景辐射率与传感器响应联系起来：

$I_j = \int_A \int_\Omega W_j(\mathbf{x}, \omega) L(\mathbf{x}, \omega) \cos \theta d\omega dA$

其中 $W_j(\mathbf{x}, \omega)$ 是像素 j 的重要性（灵敏度）函数。辐射率和重要性之间的这种对偶性使双向算法成为可能。

**重要性传输**：重要性满足伴随方程：

$W(\mathbf{x}, \omega) = W_e(\mathbf{x}, \omega) + \int_\Omega f_r(\mathbf{x}, \omega, \omega') W(\mathbf{x}, \omega') \cos \theta' d\omega'$

这种对称性导致：
- 双向路径追踪
- 光子映射（正向光，反向重要性）
- 梯度计算的伴随方法

对于针孔相机，像素 j 从针孔张成立体角 $\Omega_j$：

$I_j = \int_{\Omega_j} L(\mathbf{x}_{\text{lens}}, \omega) \cos^4 \theta d\omega$

$\cos^4 \theta$ 项考虑了：
- $\cos \theta$：投影透镜面积
- $\cos^3 \theta$：平方反比衰减和像素投影

**有限孔径相机**：对于具有孔径 $A_{\text{lens}}$ 的真实相机：

$I_j = \frac{1}{A_{\text{lens}}} \int_{A_{\text{lens}}} \int_{A_{\text{pixel}}} L(\mathbf{x}_{\text{lens}} \rightarrow \mathbf{x}_{\text{pixel}}) G(\mathbf{x}_{\text{lens}} \leftrightarrow \mathbf{x}_{\text{pixel}}) dA_{\text{pixel}} dA_{\text{lens}}$

这导致景深效果并需要仔细的采样策略。

### 算子形式与诺伊曼级数

渲染方程允许优雅的算子表述。定义传输算子 𝒯：

$(\mathcal{T}L)(\mathbf{x}, \omega) = \int_\Omega f_r(\mathbf{x}, \omega', \omega) L(\mathbf{x}, \omega') (\omega' \cdot \mathbf{n}) d\omega'$

则渲染方程变为：

$L = L_e + \mathcal{T}L$

这是第二类弗雷德霍姆方程。通过诺伊曼级数的解：

$L = (I - \mathcal{T})^{-1}L_e = \sum_{k=0}^\infty \mathcal{T}^k L_e = L_e + \mathcal{T}L_e + \mathcal{T}^2L_e + ...$

每项都有物理意义：
- $L_e$：直接照明（仅发射）
- $\mathcal{T}L_e$：单次反弹照明
- $\mathcal{T}^2L_e$：两次反弹照明
- $\mathcal{T}^k L_e$：k次反弹照明

当 $||\mathcal{T}|| < 1$ 时级数收敛，这在最大反照率 < 1 时发生。这种分解自然导致采样路径长度递增的路径追踪算法。

### 三点形式与几何耦合

渲染方程可以改写为三点形式，使几何耦合明确：

$L(\mathbf{x} \rightarrow \mathbf{x}') = L_e(\mathbf{x} \rightarrow \mathbf{x}') + \int_M f_r(\mathbf{x}'' \rightarrow \mathbf{x} \rightarrow \mathbf{x}') L(\mathbf{x}'' \rightarrow \mathbf{x}) G(\mathbf{x}'' \leftrightarrow \mathbf{x}) dA(\mathbf{x}'')$

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

**核函数性质**：传输核 $K(\mathbf{x}'' \rightarrow \mathbf{x}) = f_r G V$ 具有重要性质：
- 沿 $\mathbf{x} = \mathbf{x}''$ 奇异（需要仔细正则化）
- 在遮挡边界处不连续
- 满足互易性：$K(\mathbf{x} \rightarrow \mathbf{x}') = K(\mathbf{x}' \rightarrow \mathbf{x})$

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
   - 视锥体变为[-1,1]³立方体（NDC）

5. **屏幕空间（光栅空间）**：最终的2D图像坐标
   - 整数像素坐标
   - 原点在左上角或左下角

### 齐次坐标与变换

齐次坐标统一了平移和线性变换。3D点 **p** = (x, y, z) 变为 **p̃** = (x, y, z, 1)，而向量使用 **ṽ** = (x, y, z, 0)。

一般仿射变换矩阵：

$\mathbf{M} = \begin{bmatrix} \mathbf{A} & \mathbf{t} \\ \mathbf{0} & 1 \end{bmatrix}$

其中 **A** 是3×3线性部分，**t** 是平移。常见变换：

**平移 (tx, ty, tz)：**
$\begin{bmatrix} 1 & 0 & 0 & t_x \\ 0 & 1 & 0 & t_y \\ 0 & 0 & 1 & t_z \\ 0 & 0 & 0 & 1 \end{bmatrix}$

**绕轴 **a** 旋转角度 θ：**
$\mathbf{R} = \cos \theta \mathbf{I} + (1 - \cos \theta) \mathbf{a}\mathbf{a}^T + \sin \theta [\mathbf{a}]_\times$

其中 $[\mathbf{a}]_\times$ 是反对称叉积矩阵。

**缩放 (sx, sy, sz)：**
$\begin{bmatrix} s_x & 0 & 0 & 0 \\ 0 & s_y & 0 & 0 \\ 0 & 0 & s_z & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}$

### 法线和切线变换

法线必须变换以保持与表面垂直。给定点的变换 **M**：

$\mathbf{n}' = \frac{(\mathbf{M}^{-T})^{3×3} \mathbf{n}}{||(\mathbf{M}^{-T})^{3×3} \mathbf{n}||}$

证明：对于表面上的切线 **t**，$\mathbf{n} \cdot \mathbf{t} = 0$。变换后：
$\mathbf{n}' \cdot \mathbf{t}' = (\mathbf{M}^{-T}\mathbf{n}) \cdot (\mathbf{M}\mathbf{t}) = \mathbf{n}^T \mathbf{M}^{-1} \mathbf{M} \mathbf{t} = \mathbf{n} \cdot \mathbf{t} = 0$

对于正交标架 {**t**, **b**, **n**}：
- 正向：变换 **t** 和 **b**，然后 $\mathbf{n} = \mathbf{t} \times \mathbf{b}$
- 或如上变换 **n**，然后重建标架

**面积和体积元素**：在变换 **M** 下，微分元素缩放为：
- 长度：$dl' = ||\mathbf{M}\mathbf{v}|| dl$（对于方向 **v**）
- 面积：$dA' = |\det(\mathbf{M})| ||(\mathbf{M}^{-T})\mathbf{n}|| dA$
- 体积：$dV' = |\det(\mathbf{M})| dV$

**非均匀缩放问题**：非均匀缩放破坏各向同性：
- 球体 → 椭球体
- 各向同性BRDF → 各向异性BRDF
- 需要注意基于物理的材质

**手性保持**：当 $\det(\mathbf{M}) < 0$ 时，变换翻转方向：
- 右手系 → 左手系坐标系统
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
$(u, v) \in [-1, 1]^2 \rightarrow (r, \phi) \rightarrow (x, y)$ 在单位圆盘上

**八面体映射**：
单位球面 → 八面体 → 单位正方形
比球面坐标更好地保持面积

### 积分中的变量替换

积分的一般变量替换公式：

$\int_\Omega f(\mathbf{x}) d\mathbf{x} = \int_{\Omega'} f(\mathbf{x}(\mathbf{u})) \left|\det\left(\frac{\partial\mathbf{x}}{\partial\mathbf{u}}\right)\right| d\mathbf{u}$

对渲染至关重要：

**立体角到面积**：
$\int_\Omega L(\mathbf{x}, \omega) \cos \theta d\omega = \int_A L(\mathbf{x}, \omega(\mathbf{x}')) G(\mathbf{x} \leftrightarrow \mathbf{x}') dA'$

其中 $G(\mathbf{x} \leftrightarrow \mathbf{x}') = V(\mathbf{x} \leftrightarrow \mathbf{x}') \frac{\cos \theta \cos \theta'}{||\mathbf{x} - \mathbf{x}'||^2}$

**半球到圆盘**（用于余弦加权采样）：
映射 $(\theta, \phi) \rightarrow (r, \phi)$ 其中 $r = \sin \theta$
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

其中 n, f 是近/远平面，r, t 是近平面的右/顶。

透视除法后：
- $x_{ndc} = x_{clip} / w_{clip} \in [-1, 1]$
- $y_{ndc} = y_{clip} / w_{clip} \in [-1, 1]$
- $z_{ndc} = z_{clip} / w_{clip} \in [-1, 1]$

重要性质：
- 直线保持为直线（除了通过眼点的）
- 平面保持为平面
- 深度精度是非线性的（近处比远处精度高）

### 重心坐标与插值

对于具有顶点 $\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2$ 的三角形，重心坐标 (u, v, w) 满足：

$\mathbf{p} = u\mathbf{v}_0 + v\mathbf{v}_1 + w\mathbf{v}_2$

约束条件 $u + v + w = 1$。通过面积计算：

$u = \frac{\text{Area}(\mathbf{p}, \mathbf{v}_1, \mathbf{v}_2)}{\text{Area}(\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2)}$

性质：
- $u, v, w \in [0, 1]$ 当且仅当 **p** 在三角形内
- 线性插值：$f(\mathbf{p}) = uf_0 + vf_1 + wf_2$
- 透视正确的插值需要 1/z 校正

透视正确的属性插值：
1. 在屏幕空间插值 a/z, b/z, c/z 和 1/z
2. 恢复属性：$a = (a/z)/(1/z)$

### 微分几何与局部标架

在每个表面点，我们构建用于着色计算的局部标架：

**切空间基**：
- **n**：表面法线（$\frac{\partial\mathbf{p}}{\partial u} \times \frac{\partial\mathbf{p}}{\partial v}$ 归一化）
- **t**：切线（通常 $\frac{\partial\mathbf{p}}{\partial u}$ 归一化）
- **b**：副切线（$\mathbf{n} \times \mathbf{t}$）

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
$\begin{bmatrix} \mathbf{t}_{world} \\ \mathbf{b}_{world} \\ \mathbf{n}_{world} \end{bmatrix} = \begin{bmatrix} t_x & t_y & t_z \\ b_x & b_y & b_z \\ n_x & n_y & n_z \end{bmatrix} \begin{bmatrix} \mathbf{t}_{local} \\ \mathbf{b}_{local} \\ \mathbf{n}_{local} \end{bmatrix}$

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

$f_r(\mathbf{x}, \omega_i, \omega_o) = \frac{dL_o(\mathbf{x}, \omega_o)}{dE_i(\mathbf{x}, \omega_i)} = \frac{dL_o(\mathbf{x}, \omega_o)}{L_i(\mathbf{x}, \omega_i) \cos \theta_i d\omega_i}$ [sr⁻¹]

物理上，它表示来自方向 $\omega_i$ 的光子散射到方向 $\omega_o$ 的概率密度（归一化后）。

BRDF可以分解为分量：
$f_r = f_d + f_s + f_g + ...$

其中 $f_d$ 是漫反射，$f_s$ 是镜面反射，$f_g$ 是光泽反射等。这种分解有助于重要性采样。

### BRDF基本性质

**亥姆霍兹互易性：**
$f_r(\mathbf{x}, \omega_i, \omega_o) = f_r(\mathbf{x}, \omega_o, \omega_i)$

这源自麦克斯韦方程的时间反演对称性和细致平衡原理。它使双向路径追踪和光子映射成为可能。

**能量守恒：**
方向-半球反射率必须满足：

$\rho(\omega_i) = \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) \cos \theta_o d\omega_o \leq 1$ 对所有 $\omega_i$

对于能量守恒的BRDF，当吸收为零时等号成立。半球-半球反射率：

$\rho_{hh} = \frac{1}{\pi} \int_\Omega \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) \cos \theta_i \cos \theta_o d\omega_i d\omega_o \leq 1$

**非负性：**
$f_r(\mathbf{x}, \omega_i, \omega_o) \geq 0$

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

$S(\mathbf{x}_i, \omega_i, \mathbf{x}_o, \omega_o) = \frac{dL_o(\mathbf{x}_o, \omega_o)}{d\Phi_i(\mathbf{x}_i, \omega_i)}$ [m⁻²sr⁻¹]

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

2. **能量守恒**：$\forall \omega_i: \int_\Omega f_r(\mathbf{x}, \omega_i, \omega_o) \cos \theta_o d\omega_o \leq 1$
   - 测试：白炉测试（均匀照明）

3. **非负性**：$f_r(\mathbf{x}, \omega_i, \omega_o) \geq 0$
   - 违反会导致能量吸收异常

4. **平滑性**：$f_r$ 应该是 $C^0$ 连续的（$C^1$ 更好）
   - 不连续性导致采样困难

5. **菲涅尔行为**：对于光滑表面，当 $\theta \rightarrow \pi/2$ 时 $f_r \rightarrow 1$
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

$P(|\hat{I} - I| \leq \varepsilon) \approx 2\Phi(\varepsilon\sqrt{N}/\sigma) - 1$

其中 $\Phi$ 是正态分布CDF，$\sigma^2$ 是方差。误差以 $O(1/\sqrt{N})$ 递减，与维度无关——这对高维光传输至关重要。

## 1.5 路径积分表述

### 光传输作为路径积分

我们可以将渲染方程重新表述为对所有可能光路径的积分。长度为k的路径是：

$\bar{\mathbf{x}} = \mathbf{x}_0\mathbf{x}_1...\mathbf{x}_k$

其中 $\mathbf{x}_0$ 在光源上，$\mathbf{x}_k$ 在相机传感器上。

### 路径空间与测度

路径空间 $\bar{\Omega}_k$ 包含所有长度为k的有效路径。路径的测度是：

$d\mu(\bar{\mathbf{x}}) = dA(\mathbf{x}_0) \prod_{i=1}^{k} dA(\mathbf{x}_i)$

路径的贡献是：

$f(\bar{\mathbf{x}}) = L_e(\mathbf{x}_0 \rightarrow \mathbf{x}_1) \left(\prod_{i=1}^{k-1} f_s(\mathbf{x}_{i-1} \rightarrow \mathbf{x}_i \rightarrow \mathbf{x}_{i+1}) G(\mathbf{x}_i \leftrightarrow \mathbf{x}_{i+1})\right) W(\mathbf{x}_{k-1} \rightarrow \mathbf{x}_k)$

其中G是几何因子：

$G(\mathbf{x} \leftrightarrow \mathbf{x}') = V(\mathbf{x} \leftrightarrow \mathbf{x}') \frac{\cos \theta \cos \theta'}{||\mathbf{x} - \mathbf{x}'||^2}$

### 与费曼路径积分的联系

路径积分表述类似于费曼对量子力学的方法：

$I = \sum_{k=2}^\infty \int_{\bar{\Omega}_k} f(\bar{\mathbf{x}}) d\mu(\bar{\mathbf{x}})$

这个对所有路径长度的无限和捕获了全局光照。每一项代表k-1次反弹的路径。

### 递归表述与诺伊曼级数

点上的值满足递归关系：

$L(\mathbf{x}, \omega) = L_e(\mathbf{x}, \omega) + \int_M f_s(\mathbf{y} \rightarrow \mathbf{x} \rightarrow \omega) L(\mathbf{y}, \mathbf{x} - \mathbf{y}) G(\mathbf{y} \leftrightarrow \mathbf{x}) dA(\mathbf{y})$

这导致诺伊曼级数解：

$L = L^{(0)} + L^{(1)} + L^{(2)} + ...$

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
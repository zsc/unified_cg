# 第18章：统计光学与散斑

当相干光从粗糙表面反射或透过随机介质传播时，会产生一种称为散斑的随机强度图案。这种现象在激光物理、光学成像和计算机图形学中都有重要应用。本章将建立光场的统计描述框架，分析散斑的形成机理和统计特性，并探讨其在成像系统中的影响。我们将看到，散斑不仅是相干成像的固有限制，也可以作为测量和成像的有力工具。

## 学习目标

完成本章学习后，您将能够：
1. 使用统计方法描述随机光场的一阶和二阶特性
2. 推导散斑强度的概率分布（Rayleigh分布和Rice分布）
3. 分析粗糙表面散射的统计模型，包括Kirchhoff近似
4. 计算成像系统中散斑的统计特性和相关函数
5. 设计散斑抑制方案和利用散斑进行测量
6. 将散斑现象与体积渲染方程建立联系

## 章节大纲

## 18.1 光场的统计描述

在许多实际情况中，光场的精确值是未知的或快速变化的，需要使用统计方法来描述。这种随机性可能源于光源的热涨落、传播介质的不均匀性或散射表面的粗糙度。统计光学为理解这些现象提供了强大的数学框架，将确定性的麦克斯韦方程与随机过程理论相结合。

### 18.1.1 随机光场的表征

考虑一个标量光场 $U(\mathbf{r}, t)$，它是空间位置 $\mathbf{r}$ 和时间 $t$ 的随机函数。我们使用复解析信号表示：

$U(\mathbf{r}, t) = A(\mathbf{r}, t) \exp[i\varphi(\mathbf{r}, t)]$

其中 $A(\mathbf{r}, t)$ 是振幅，$\varphi(\mathbf{r}, t)$ 是相位。

对于随机过程，我们定义系综平均：

$\langle U(\mathbf{r}, t) \rangle = \int U(\mathbf{r}, t) p(U) dU$

其中 $p(U)$ 是概率密度函数。

**平稳性和各态历经性**：

如果统计特性不随时间平移改变，则过程是平稳的：
$\langle U(\mathbf{r}, t) \rangle = \langle U(\mathbf{r}, t + \tau) \rangle$

更强的条件是广义平稳（弱平稳），要求二阶矩也平移不变：
$\langle U^*(\mathbf{r}, t)U(\mathbf{r}', t + \tau) \rangle = \Gamma(\mathbf{r}, \mathbf{r}', \tau)$

如果系综平均等于时间平均，则过程是各态历经的：
$\langle U(\mathbf{r}, t) \rangle = \lim_{T\to\infty} \frac{1}{T} \int_{0}^{T} U(\mathbf{r}, t) dt$

各态历经性允许我们用单一实现的时间平均代替系综平均，这在实验测量中极为重要。

**空间平稳性**：

类似地，如果统计特性不随空间平移改变：
$\langle U(\mathbf{r}, t) \rangle = \langle U(\mathbf{r} + \mathbf{\Delta r}, t) \rangle$

这在均匀介质和平移不变系统中常见。对于空间平稳过程：
$\langle U^*(\mathbf{r}, t)U(\mathbf{r} + \mathbf{\rho}, t) \rangle = \Gamma(\mathbf{\rho}, t)$

仅依赖于相对位置 $\mathbf{\rho}$。

**准单色近似**：

实际光源通常是准单色的，即：
$U(\mathbf{r}, t) = U_0(\mathbf{r}, t) \exp(-i\omega_0 t)$

其中 $U_0(\mathbf{r}, t)$ 是缓变包络（相对于载波频率 $\omega_0$），满足：
$|\partial U_0/\partial t| \ll \omega_0|U_0|$

这允许我们分离快速振荡的载波和缓变的包络，大大简化分析。

**解析信号表示**：

实际测量的光场是实函数，但使用复解析信号更方便：
$U^{+}(t) = \frac{1}{2\pi} \int_{0}^{\infty} \tilde{U}(\omega) \exp(-i\omega t) d\omega$

其中 $\tilde{U}(\omega)$ 是实场的傅里叶变换。解析信号的实部给出物理场：
$E(t) = 2\text{Re}[U^{+}(t)]$

**随机性的分类**：

1. **时间随机性**：源于光源的量子涨落或热涨落
   - 特征时间：相干时间 $\tau_c$
   - 频域表现：有限带宽 $\Delta\omega \approx 1/\tau_c$

2. **空间随机性**：源于波前畸变或散射
   - 特征长度：相干长度 $lc$
   - 角域表现：有限角度谱 $\Delta\theta \approx \lambda/lc$

3. **偏振随机性**：源于双折射或去偏振
   - 需要矢量理论描述
   - 使用Stokes参数或相干矩阵

### 18.1.2 一阶统计量

**平均强度**：
$I(\mathbf{r}, t) = \langle|U(\mathbf{r}, t)|^2\rangle = \langle U(\mathbf{r}, t)U^*(\mathbf{r}, t) \rangle$

这是最基本的统计量，直接对应于探测器测量的平均功率密度。对于各态历经过程：
$I(\mathbf{r}) = \lim_{T\to\infty} \frac{1}{T} \int_{0}^{T} |U(\mathbf{r}, t)|^2 dt$

**归一化一阶相关函数**：
$g^{(1)}(\mathbf{r}_1, t_1; \mathbf{r}_2, t_2) = \frac{\langle U^*(\mathbf{r}_1, t_1)U(\mathbf{r}_2, t_2) \rangle}{\sqrt{I(\mathbf{r}_1, t_1)I(\mathbf{r}_2, t_2)}}$

对于平稳随机过程，相关函数仅依赖于时间差 $\tau = t_2 - t_1$：
$g^{(1)}(\mathbf{r}_1, \mathbf{r}_2, \tau) = \frac{\langle U^*(\mathbf{r}_1, t)U(\mathbf{r}_2, t + \tau) \rangle}{\sqrt{I(\mathbf{r}_1)I(\mathbf{r}_2)}}$

**互强度函数**：

定义互强度为：
$J(\mathbf{r}_1, \mathbf{r}_2, \tau) = \langle U^*(\mathbf{r}_1, t)U(\mathbf{r}_2, t + \tau) \rangle$

它描述了两点间场的相关性，是相干理论的核心量。互强度满足波动方程：
$\left(\nabla_1^2 - \frac{1}{c^2} \frac{\partial^2}{\partial t^2}\right)J(\mathbf{r}_1, \mathbf{r}_2, \tau) = 0$

对于准单色光（$\Delta\omega/\omega_0 \ll 1$）：
$J(\mathbf{r}_1, \mathbf{r}_2, \tau) \approx J(\mathbf{r}_1, \mathbf{r}_2, 0) \exp(-i\omega_0\tau)$

**复相干度**：

归一化的互强度定义为复相干度：
$\gamma(\mathbf{r}_1, \mathbf{r}_2, \tau) = \frac{J(\mathbf{r}_1, \mathbf{r}_2, \tau)}{\sqrt{I(\mathbf{r}_1)I(\mathbf{r}_2)}}$

其模值 $|\gamma|$ 描述相干程度，相位 $\arg(\gamma)$ 描述两点间的平均相位差。

**相干度的物理意义**：

复相干度 $|\gamma^{(1)}|$ 的值域为 $[0, 1]$：
- $|\gamma^{(1)}| = 1$：完全相干
- $0 < |\gamma^{(1)}| < 1$：部分相干  
- $|\gamma^{(1)}| = 0$：非相干

相干度直接影响干涉条纹的可见度：
$V = \frac{I_{\text{max}} - I_{\text{min}}}{I_{\text{max}} + I_{\text{min}}} = \frac{2\sqrt{I_1I_2}}{I_1 + I_2} |\gamma^{(1)}|$

对于等强度光束（$I_1 = I_2$），$V = |\gamma^{(1)}|$。

**Wolf方程**：

互强度的传播由Wolf方程描述：
$J(P_1, P_2) = \iint K^*(P_1, Q_1)K(P_2, Q_2)J_0(Q_1, Q_2) dQ_1 dQ_2$

其中 $K$ 是从源平面到观察平面的传播核（Green函数）。

**交叉谱密度**：

时域互强度的傅里叶变换给出交叉谱密度：
$W(\mathbf{r}_1, \mathbf{r}_2, \omega) = \int J(\mathbf{r}_1, \mathbf{r}_2, \tau) \exp(i\omega\tau) d\tau$

它描述了不同频率分量的空间相干性，满足：
$J(\mathbf{r}_1, \mathbf{r}_2, \tau) = \frac{1}{2\pi} \int W(\mathbf{r}_1, \mathbf{r}_2, \omega) \exp(-i\omega\tau) d\omega$

**相干模式分解**：

互强度可分解为相干模式的叠加：
$J(\mathbf{r}_1, \mathbf{r}_2) = \sum_{n} \lambda_n \psi_n^*(\mathbf{r}_1)\psi_n(\mathbf{r}_2)$

其中 $\lambda_n$ 是特征值，$\psi_n(\mathbf{r})$ 是正交模式，满足：
$\int J(\mathbf{r}, \mathbf{r}')\psi_n(\mathbf{r}') d\mathbf{r}' = \lambda_n \psi_n(\mathbf{r})$

这提供了部分相干光的模式表示。

### 18.1.3 二阶统计量

**强度相关函数**：
$G^{(2)}(\mathbf{r}_1, t_1; \mathbf{r}_2, t_2) = \langle I(\mathbf{r}_1, t_1)I(\mathbf{r}_2, t_2) \rangle$

二阶相关函数描述强度涨落的相关性，在量子光学和经典统计光学中都有重要意义。

**归一化二阶相关函数**：
$g^{(2)}(\mathbf{r}_1, t_1; \mathbf{r}_2, t_2) = \frac{G^{(2)}(\mathbf{r}_1, t_1; \mathbf{r}_2, t_2)}{\langle I(\mathbf{r}_1, t_1) \rangle \langle I(\mathbf{r}_2, t_2) \rangle}$

对于平稳过程：
$g^{(2)}(\mathbf{r}_1, \mathbf{r}_2, \tau) = \frac{\langle I(\mathbf{r}_1, t)I(\mathbf{r}_2, t + \tau) \rangle}{\langle I(\mathbf{r}_1) \rangle \langle I(\mathbf{r}_2) \rangle}$

**Siegert关系**：

对于高斯随机过程，存在重要关系：
$g^{(2)}(\mathbf{r}_1, t_1; \mathbf{r}_2, t_2) = 1 + |g^{(1)}(\mathbf{r}_1, t_1; \mathbf{r}_2, t_2)|^2$

这将二阶统计与一阶统计联系起来，是高斯场的特征。对于非高斯场，此关系不成立，偏差量化了非高斯性。

**高斯矩定理**：

对于零均值复高斯随机过程，所有高阶矩可由二阶矩完全确定：
$\langle U_1U_2...U_nU^*_{n+1}...U^*_{2n} \rangle = \sum_{\text{all pairings}} \prod \langle U_iU^*_j \rangle$

奇数阶矩为零：
$\langle U_1U_2...U_n \rangle = 0 \quad (\text{n is odd})$

这极大简化了统计计算。例如，四阶矩：
$\langle U_1U_2U_3^*U_4^* \rangle = \langle U_1U_3^* \rangle \langle U_2U_4^* \rangle + \langle U_1U_4^* \rangle \langle U_2U_3^* \rangle$

**强度涨落的物理含义**：

二阶相关函数描述强度涨落：
- $g^{(2)}(0) = 1$：泊松统计（相干态）
- $g^{(2)}(0) > 1$：超泊松统计（热光、混沌光）
- $g^{(2)}(0) < 1$：亚泊松统计（需要量子描述）

强度方差与二阶相关的关系：
$\text{Var}(I) = \langle I^2 \rangle - \langle I \rangle^2 = \langle I \rangle^2[g^{(2)}(0) - 1]$

**Hanbury Brown-Twiss效应**：

强度相关可用于测量恒星角直径。对于扩展热光源：
$g^{(2)}(d) = 1 + |\gamma^{(1)}(d)|^2$

其中 $d$ 是探测器间距。当 $d$ 增大时，$|\gamma^{(1)}(d)|$ 下降，通过测量 $g^{(2)}$ 的空间分布可推断源的角尺寸：
$\theta \approx \lambda/d_0$

其中 $d_0$ 是 $g^{(2)}$ 降到 1.5 时的间距。

**聚束与反聚束**：

- **聚束（Bunching）**：$g^{(2)}(0) > g^{(2)}(\infty)$
  - 光子倾向于成群到达
  - 热光的典型特征
  
- **反聚束（Antibunching）**：$g^{(2)}(0) < g^{(2)}(\infty)$
  - 光子倾向于等间隔到达
  - 量子光源的标志

**高阶相关函数**：

n阶强度相关：
$g^{(n)}(\tau_1, ..., \tau_{n-1}) = \frac{\langle I(t)I(t+\tau_1)...I(t+\tau_{n-1}) \rangle}{\langle I \rangle^n}$

对于高斯光：
$g^{(n)} = n! \quad (\text{完全相干})$
$g^{(n)} = 1 \quad (\text{完全非相干})$

**相关时间和相关面积**：

强度相关时间：
$\tau_i = \int_{0}^{\infty} [g^{(2)}(\tau) - 1] d\tau$

强度相关面积：
$A_i = \iint [g^{(2)}(\mathbf{r}_1, \mathbf{r}_2) - 1] d^2\mathbf{r}_2$

这些量表征了强度涨落的时空尺度。

### 18.1.4 空间相干函数

定义互相干函数：
$\Gamma(\mathbf{r}_1, \mathbf{r}_2, \tau) = \langle U^*(\mathbf{r}_1, t)U(\mathbf{r}_2, t + \tau) \rangle$

这是描述光场空间相干性的基本函数，满足Hermite性质：
$\Gamma(\mathbf{r}_1, \mathbf{r}_2, \tau) = \Gamma^*(\mathbf{r}_2, \mathbf{r}_1, -\tau)$

复相干度：
$\gamma(\mathbf{r}_1, \mathbf{r}_2, \tau) = \frac{\Gamma(\mathbf{r}_1, \mathbf{r}_2, \tau)}{\sqrt{\Gamma(\mathbf{r}_1, \mathbf{r}_1, 0)\Gamma(\mathbf{r}_2, \mathbf{r}_2, 0)}}$

其模值满足：$0 \le |\gamma(\mathbf{r}_1, \mathbf{r}_2, \tau)| \le 1$

**相干面积和体积**：

相干面积定义为：
$A_c = \iint |\gamma(\mathbf{r}_1, \mathbf{r}_2, 0)|^2 d^2\mathbf{r}_2$

它量化了保持相干性的空间范围。对于均匀场：
$A_c \approx (\lambda z/D)^2$

其中 $D$ 是源的特征尺寸。

**van Cittert-Zernike定理**：

非相干扩展源产生的场，其空间相干性由源的强度分布的傅里叶变换决定：

$\gamma_{12} = \frac{\iint I(\xi, \eta) \exp[ik(\xi x_{12} + \eta y_{12})/z] d\xi d\eta}{\iint I(\xi, \eta) d\xi d\eta}$

其中 $(\xi, \eta)$ 是源坐标，$x_{12} = x_2 - x_1$， $y_{12} = y_2 - y_1$。

这个定理的重要推论：
- 均匀圆形源：$\gamma(\rho) = \frac{2J_1(kD\rho/2z)}{kD\rho/2z}$
- 均匀矩形源：$\gamma(x, y) = \text{sinc}(kD_x x/z) \times \text{sinc}(kD_y y/z)$

**相干长度的估计**：

横向相干长度（空间相干性）：
$lc_{\perp} \approx \lambda z/D$

这是第一个零点的位置，表示场保持相干的横向距离。

纵向相干长度（时间相干性）：
$lc_{\parallel} = c\tau_c = \lambda^2/\Delta\lambda$

其中 $\tau_c$ 是相干时间，$\Delta\lambda$ 是光源带宽。

对于高斯谱线：
$lc_{\parallel} = (2\ln2/\pi) \times (\lambda^2/\Delta\lambda) \approx 0.44\lambda^2/\Delta\lambda$

**相干体积**：

三维相干体积定义为：
$V_c = \iiint |\gamma(\mathbf{r}, \mathbf{r}', 0)|^2 d^3\mathbf{r}'$

对于各向同性场：
$V_c \approx lc_{\perp}^2 \times lc_{\parallel} \approx (\lambda z/D)^2 \times (\lambda^2/\Delta\lambda)$

它表征了光场中相干单元的大小，对理解散斑统计至关重要。

**传播中的相干性演化**：

自由空间传播中，相干函数通过广义Huygens-Fresnel原理演化：
$\Gamma(P_1, P_2) = \iint K^*(P_1, Q_1)K(P_2, Q_2)\Gamma_0(Q_1, Q_2) dQ_1 dQ_2$

其中 $K$ 是传播核（Green函数）：
$K(P, Q) = \frac{\exp(ikr)}{i\lambda r} \times \frac{1 + \cos \theta}{2}$

对于傍轴近似：
$K(P, Q) \approx \frac{\exp(ikz)}{i\lambda z} \exp\left[ik|\mathbf{r}_p - \mathbf{r}_q|^2/2z\right]$

**相干性的增强与降低**：

1. **空间滤波增强相干性**：
   - 针孔滤波：选择单一相干模式
   - 模式选择：保留低阶模式
   
2. **部分相干性的产生**：
   - 旋转毛玻璃：时间平均降低相干性
   - 多模光纤：模式混合
   - 扩展源：空间非相干叠加

### 18.1.5 高斯随机过程模型

许多光学系统中，光场可以近似为圆复高斯随机过程。其特点是：

1. 实部和虚部独立且具有相同的高斯分布
2. 概率密度函数：
   $p(U) = \frac{1}{\pi\sigma^2} \exp(-|U|^2/\sigma^2)$
   
3. 所有高阶矩可由二阶矩确定（高斯矩定理）

**中心极限定理的应用**：

当大量独立随机贡献叠加时，根据中心极限定理：
$U = \sum_{n} U_n \to \text{高斯分布} \quad (N \to \infty)$

这解释了为什么散斑场通常是高斯的。

**联合概率分布**：

对于两个位置的场，联合高斯分布为：
$p(U_1, U_2) = \frac{1}{\pi^2\text{det}[\mathbf{C}]} \exp(-\mathbf{U}^{\dagger}\mathbf{C}^{-1}\mathbf{U})$

其中协方差矩阵：
$\mathbf{C} = \begin{bmatrix} \langle|U_1|^2\rangle & \langle U_1U_2^* \rangle \\ \langle U_1^*U_2 \rangle & \langle|U_2|^2\rangle \end{bmatrix}$

**条件概率和预测**：

给定 $U_1$， $U_2$ 的条件期望值：
$\langle U_2|U_1 \rangle = \left(\frac{\langle U_2U_1^* \rangle}{\langle|U_1|^2\rangle}\right) U_1$

这在自适应光学和相位恢复中有重要应用。

### 18.1.6 功率谱密度

通过Wiener-Khintchine定理，相干函数与功率谱密度通过傅里叶变换相联系：

$S(\mathbf{k}, \omega) = \iint \Gamma(\mathbf{r}, \tau) \exp[i(\mathbf{k}\cdot\mathbf{r} - \omega\tau)] d^3\mathbf{r} d\tau$

其中 $\mathbf{k}$ 是空间频率，$\omega$ 是时间频率。

**相干长度和相干时间**：
$lc = \int |\gamma(\mathbf{r}, 0)|^2 d\mathbf{r}$
$\tau_c = \int |\gamma(0, \tau)|^2 d\tau$

**谱相干度**：

定义谱相干度函数：
$W(\mathbf{r}_1, \mathbf{r}_2, \omega) = \int \Gamma(\mathbf{r}_1, \mathbf{r}_2, \tau) \exp(i\omega\tau) d\tau$

它描述了不同频率分量的空间相干性。

**准单色近似**：

当 $\Delta\omega/\omega_0 \ll 1$ 时：
$S(\mathbf{k}, \omega) \approx S(\mathbf{k}) \times S(\omega)$

空间和时间特性可分离处理。

**功率谱的物理意义**：

- $S(\mathbf{k})$ 描述空间结构的尺度分布
- $S(\omega)$ 描述时间涨落的频率成分
- 带宽 $\Delta\omega$ 决定相干时间：$\tau_c \approx 2\pi/\Delta\omega$

**测量方法**：

1. 直接法：傅里叶变换测量的相关函数
2. 间接法：通过干涉仪扫描获得
3. 光谱分析仪：直接测量 $S(\omega)$

## 18.2 散斑的形成与统计

散斑是相干光照射粗糙表面或通过随机介质后形成的随机强度图案。理解散斑的统计特性对于光学系统设计和成像质量分析至关重要。

### 18.2.1 散斑现象的物理起源

当相干光照射到粗糙度与波长相当或更大的表面时，不同位置的散射光具有随机相位差。在观察平面上，这些散射光相干叠加形成散斑图案。

考虑 $N$ 个散射点，观察点的复振幅为：

$U = \sum_{k=1}^{N} A_k \exp(i\varphi_k)$

其中 $A_k$ 和 $\varphi_k$ 分别是第 $k$ 个散射光的振幅和相位。

对于完全发展的散斑（$N \to \infty$），根据中心极限定理，$U$ 的实部和虚部趋于高斯分布。

**相位随机性的来源**：

1. **表面高度变化**：
   $\varphi_k = 2k h_k \cos \theta_i + \mathbf{k}\cdot\mathbf{r}_k$
   
   其中 $h_k$ 是表面高度，$\theta_i$ 是入射角，$\mathbf{r}_k$ 是散射点位置。

2. **路径长度差异**：
   对于体散射介质：
   $\varphi_k = k \int n(s) ds$
   
   其中积分沿第 $k$ 条散射路径。

**完全发展散斑的条件**：

1. 表面粗糙度：$\sigma_h > \lambda/4$
2. 散射点数：$N \gg 1$
3. 相位均匀分布：$\varphi_k \in [0, 2\pi]$
4. 振幅统计独立性

**散斑的空间尺度**：

散斑的典型尺寸由系统的点扩散函数决定：
- 自由空间：$\delta_s \approx \lambda z/D$
- 成像系统：$\delta_s \approx \lambda f\#$
- 其中 $f\# = f/D$ 是光圈数

### 18.2.2 振幅和强度的一阶统计

**振幅分布**：

设 $U = U_r + iU_i$，其中 $U_r$ 和 $U_i$ 是独立的零均值高斯随机变量，方差为 $\sigma^2$。

振幅 $A = |U| = \sqrt{U_r^2 + U_i^2}$ 的概率密度函数为：

$p(A) = \frac{A}{\sigma^2} \exp(-A^2/2\sigma^2), \quad A \ge 0$

这是Rayleigh分布。

**强度分布**：

强度 $I = |U|^2 = A^2$ 的概率密度函数为：

$p(I) = \frac{1}{\langle I \rangle} \exp(-I/\langle I \rangle), \quad I \ge 0$

其中 $\langle I \rangle = 2\sigma^2$ 是平均强度。这是负指数分布。

**统计特性**：

振幅的矩：
- 平均值：$\langle A \rangle = \sigma\sqrt{\pi/2} \approx 1.253\sigma$
- 方差：$\text{Var}(A) = (2 - \pi/2)\sigma^2 \approx 0.429\sigma^2$
- 最概然值：$A_{mp} = \sigma$

强度的矩：
- 平均值：$\langle I \rangle = 2\sigma^2$
- 方差：$\text{Var}(I) = \langle I \rangle^2$
- 标准差：$\sigma_i = \langle I \rangle$
- n阶矩：$\langle I^n \rangle = n! \langle I \rangle^n$

**概率分布的验证**：

实验上可通过直方图验证：
1. 测量大量独立散斑点的强度
2. 归一化：$I' = I/\langle I \rangle$
3. 验证：$p(I') = \exp(-I')$

**极值统计**：

$N$个独立散斑中的最大强度：
$P(I_{\text{max}} < I) = [1 - \exp(-I/\langle I \rangle)]^N$

大$N$时，$I_{\text{max}} \approx \langle I \rangle \ln N$

### 18.2.3 部分发展散斑：Rice分布

当存在确定性分量（如镜面反射）时，复振幅为：

$U = U_s + U_n$

其中 $U_s$ 是确定性分量，$U_n$ 是随机分量。

振幅的概率密度函数变为Rice分布：

$p(A) = \frac{A}{\sigma^2} \exp\left(-\frac{A^2 + |U_s|^2}{2\sigma^2}\right) I_0\left(\frac{A|U_s|}{\sigma^2}\right)$

其中 $I_0$ 是零阶修正贝塞尔函数。

**Rice参数**：
$K = |U_s|^2/(2\sigma^2)$

- $K = 0$：完全发展散斑（Rayleigh分布）
- $K \gg 1$：接近高斯分布

### 18.2.4 散斑的二阶统计

**强度相关函数**：

对于平稳散斑场，归一化强度相关函数为：

$g^{(2)}(\Delta\mathbf{r}) = \frac{\langle I(\mathbf{r})I(\mathbf{r} + \Delta\mathbf{r}) \rangle}{\langle I \rangle^2} = 1 + |\mu(\Delta\mathbf{r})|^2$

其中 $\mu(\Delta\mathbf{r})$ 是复相干因子。

**散斑尺寸**：

平均散斑尺寸可通过相关函数的半高宽定义：

$\delta_s = \int |\mu(\Delta\mathbf{r})|^2 d^2(\Delta\mathbf{r})$

对于远场散斑：
$\delta_s \approx \lambda z/D$

其中 $\lambda$ 是波长，$z$ 是观察距离，$D$ 是照明区域尺寸。

### 18.2.5 散斑对比度

散斑对比度定义为：

$C = \sigma_i/\langle I \rangle$

其中 $\sigma_i$ 是强度标准差。

对于不同类型的散斑：

1. **完全发展散斑**：$C = 1$
2. **部分发展散斑**：$C = \sqrt{1/(1 + K)}$
3. **多重独立散斑叠加**：$C = 1/\sqrt{N}$

### 18.2.6 偏振散斑

当考虑偏振时，需要使用Stokes参数描述：

$\mathbf{S}_i = [S_0, S_1, S_2, S_3]^T$

其中：
- $S_0 = \langle I_x \rangle + \langle I_y \rangle$（总强度）
- $S_1 = \langle I_x \rangle - \langle I_y \rangle$（线偏振）
- $S_2 = \langle 2\text{Re}(U_xU_y^*) \rangle$（45°线偏振）
- $S_3 = \langle 2\text{Im}(U_xU_y^*) \rangle$（圆偏振）

偏振度：
$P = \sqrt{S_1^2 + S_2^2 + S_3^2}/S_0$

## 18.3 粗糙表面散射

粗糙表面的散射是产生散斑的主要机制之一。本节将建立粗糙表面的统计模型，并分析其散射特性。

### 18.3.1 表面粗糙度的统计描述

粗糙表面可用高度函数 $h(x, y)$ 描述，其统计特性包括：

**高度分布**：

假设高度服从高斯分布：
$p(h) = \frac{1}{\sigma_h\sqrt{2\pi}} \exp(-h^2/2\sigma_h^2)$

其中 $\sigma_h$ 是均方根粗糙度。

**自相关函数**：

表面高度的自相关函数：
$C(\xi, \eta) = \langle h(x, y)h(x + \xi, y + \eta) \rangle$

常用模型：
- 高斯相关：$C(\rho) = \sigma_h^2 \exp(-\rho^2/l_c^2)$
- 指数相关：$C(\rho) = \sigma_h^2 \exp(-\rho/l_c)$

其中 $l_c$ 是相关长度，$\rho = \sqrt{\xi^2 + \eta^2}$。

**功率谱密度**：

表面的功率谱密度是自相关函数的傅里叶变换：
$S(k_x, k_y) = \iint C(\xi, \eta) \exp[-i(k_x\xi + k_y\eta)] d\xi d\eta$

### 18.3.2 Kirchhoff近似

对于缓变粗糙表面（$\sigma_h \ll lc$），可使用Kirchhoff近似（也称物理光学近似）。

**边界条件**：

在表面上，场满足：
$U^{+} - U^{-} = 0$
$\partial U^{+}/\partial n - \partial U^{-}/\partial n = 0$

其中 $n$ 是表面法向量。

**散射振幅**：

入射平面波 $\exp(i\mathbf{k}_i\cdot\mathbf{r})$ 的散射场：

$U_s(\mathbf{r}) = \frac{ik}{2\pi} \iint \frac{\exp(ikR)}{R} [1 + \cos \theta] \exp[i(\mathbf{k}_i - \mathbf{k}_s)\cdot\mathbf{r}'] dS'$

其中 $R = |\mathbf{r} - \mathbf{r}'|$，$\theta$ 是散射角。

### 18.3.3 相位扰动模型

对于小粗糙度（$\sigma_h \ll \lambda$），可使用相位扰动近似：

**反射系数**：

$r(x, y) \approx r_0 \exp[i2k_i \cos \theta_i h(x, y)]$

其中 $r_0$ 是光滑表面的反射系数，$\theta_i$ 是入射角。

**散射截面**：

微分散射截面：
$d\sigma/d\Omega = A \cos^2 \theta_i \cos \theta_s |r_0|^2 S[(\mathbf{k}_s - \mathbf{k}_i)_{\parallel}]$

其中 $(\mathbf{k}_s - \mathbf{k}_i)_{\parallel}$ 是动量转移的平行分量。

### 18.3.4 散射系数的统计特性

**相干散射**（镜面方向）：

$\langle U_s \rangle = r_0 \exp(-g)$

其中 $g = (2k \cos \theta_i \sigma_h)^2$ 是粗糙度因子。

**非相干散射**：

$\langle|U_s|^2\rangle - |\langle U_s \rangle|^2 = |r_0|^2 [\exp(2g \text{Re}[\rho]) - \exp(-2g)]$

其中 $\rho$ 是归一化相关函数。

**角度分布**：

散射强度的角度分布取决于表面相关函数：

$I(\theta_s) \propto \iint C(\xi, \eta) \exp[ik(\sin \theta_s - \sin \theta_i)\xi] d\xi d\eta$

### 18.3.5 多重散射效应

对于强粗糙表面，需考虑多重散射：

**增强后向散射**：

在精确后向方向（$\theta_s = -\theta_i$），由于相干效应，散射强度增强：

$I(180^\circ)/I(\theta) \approx 2$

**阴影效应**：

掠入射时的阴影函数：
$S(\theta_i, \theta_s) = [1 + \Lambda(\theta_i)]^{-1}[1 + \Lambda(\theta_s)]^{-1}$

其中 $\Lambda(\theta) = (\sqrt{2} \sigma_h/lc) \cot \theta$。

### 18.3.6 与体积渲染的联系

粗糙表面散射可纳入体积渲染框架：

**等效体积表示**：

将表面散射转换为薄层体积散射：

$\sigma_s(z) = \sigma_0 \delta(z - h(x, y))$

其中 $\sigma_0$ 包含表面反射特性。

**BSDF与相位函数**：

表面BSDF可表示为等效相位函数：

$p(\mathbf{\omega}_i \to \mathbf{\omega}_s) = \frac{4\pi}{\sigma_0} f_s(\mathbf{\omega}_i, \mathbf{\omega}_s) \delta(z)$

这建立了表面散射与体积渲染方程的联系。

## 18.4 相干成像中的散斑

相干成像系统中的散斑特性与成像系统的点扩散函数密切相关。本节分析成像过程中散斑的形成和传播。

### 18.4.1 成像系统的统计描述

**成像方程**：

相干成像系统的输出场：

$U_i(\mathbf{r}) = \iint h(\mathbf{r} - \mathbf{r}') U_0(\mathbf{r}') d^2\mathbf{r}'$

其中 $h(\mathbf{r})$ 是相干点扩散函数，$U_0(\mathbf{r})$ 是物平面的场分布。

**统计传递**：

输出场的统计特性：

$\langle U_i(\mathbf{r}_1)U_i^*(\mathbf{r}_2) \rangle = \iiiint h(\mathbf{r}_1 - \mathbf{r}')h^*(\mathbf{r}_2 - \mathbf{r}'') \times \langle U_0(\mathbf{r}')U_0^*(\mathbf{r}'') \rangle d^2\mathbf{r}' d^2\mathbf{r}''$

### 18.4.2 图像散斑与物体散斑

**物体散斑**：

当物体本身产生散斑场时：
- 散斑尺寸由照明条件决定
- 成像系统对散斑进行空间滤波
- 输出散斑受系统分辨率限制

**图像散斑**：

成像系统引入的散斑：
- 由系统孔径的衍射效应产生
- 散斑尺寸：$\delta_i \approx 1.22\lambda f/D$
- 其中 $f$ 是焦距，$D$ 是孔径直径

### 18.4.3 传递函数的统计特性

**相干传递函数**：

$\text{CTF} = \mathcal{F}\{h(\mathbf{r})\} = P(\mathbf{k})$

其中 $P(\mathbf{k})$ 是光瞳函数。

**散斑传递函数**：

对于散斑输入，有效传递函数变为：

$\text{STF}(\mathbf{k}) = |\text{CTF}(\mathbf{k})|^2 = |P(\mathbf{k})|^2$

这相当于非相干成像的光学传递函数(OTF)。

### 18.4.4 分辨率与散斑的关系

**Rayleigh判据的修正**：

在散斑场中，两点的可分辨性取决于：

1. 确定性分辨率：$\Delta r > 1.22\lambda/\text{NA}$
2. 统计分辨率：需要足够的散斑平均

**散斑平均数**：

分辨单元内的独立散斑数：

$N = (\text{物体尺寸}/\text{散斑尺寸}) \times (\text{成像孔径}/\text{衍射极限})$

**信噪比**：

散斑场中的信噪比：
$\text{SNR} = \sqrt{N} \times (\text{信号对比度})$

### 18.4.5 动态散斑

**时间相关性**：

动态散斑的时间相关函数：

$g^{(1)}(\tau) = \frac{\langle U^*(t)U(t + \tau) \rangle}{\langle|U|^2\rangle}$

**速度测量**：

通过散斑的时间相关性可测量物体运动：

$\tau_c = \lambda/(2v \sin(\theta/2))$

其中 $v$ 是横向速度，$\theta$ 是散射角。

**散斑跟踪**：

横向位移 $\Delta x$ 导致的相位变化：
$\Delta\varphi = 2\pi(\Delta x/\lambda) \sin \theta$

### 18.4.6 部分相干成像

**van Cittert-Zernike定理的应用**：

扩展光源的空间相干性：

$\mu_{12} = \frac{\iint I(\xi, \eta) \exp[ik(\xi x_{12} + \eta y_{12})/z] d\xi d\eta}{\iint I(\xi, \eta) d\xi d\eta}$

**散斑对比度降低**：

部分相干照明下的散斑对比度：

$C = C_0 \times |\mu_{12}|$

其中 $C_0$ 是完全相干时的对比度。

### 18.4.7 计算成像中的散斑模拟

**数值生成**：

散斑场的数值模拟：

$U(x, y) = \sum_{n} A_n \exp[i(k_{xn}x + k_{yn}y + \varphi_n)]$

其中 $\mathbf{f}_n$ 是空间频率，$\varphi_n$ 是随机相位。

**与体积渲染的结合**：

在体积渲染中加入散斑效应：

$L(x, y) = \int \sigma_t(s) \exp\left[-\int \sigma_t(s') ds'\right] \times [L_d(s) + S(s, x, y)] ds$

其中 $S(s, x, y)$ 是散斑调制项。

## 18.5 散斑的应用与抑制

散斑既是相干成像的限制因素，也是强大的测量工具。本节探讨散斑的实际应用和抑制技术。

### 18.5.1 散斑干涉测量

**电子散斑干涉（ESPI）**：

物体变形前后的散斑场：
$U_1 = U_r + U_0$
$U_2 = U_r + U_0 \exp(i\Delta\varphi)$

其中 $U_r$ 是参考光，$\Delta\varphi$ 是变形引起的相位变化。

干涉条纹对应于：
$\Delta\varphi = 2\pi(2\mathbf{d}\cdot\mathbf{n})/\lambda = 2\pi m$

其中 $\mathbf{d}$ 是位移矢量，$\mathbf{n}$ 是观察方向。

**剪切散斑干涉**：

测量位移梯度：
$I(x, y) = |U(x, y) + U(x + \delta x, y)|^2$

条纹对应于：
$\partial u/\partial x = m\lambda/(2\delta x)$

### 18.5.2 散斑相关技术

**数字图像相关（DIC）**：

互相关函数：
$C(\Delta x, \Delta y) = \iint I_1(x, y)I_2(x + \Delta x, y + \Delta y) dx dy$

峰值位置给出位移量。

**激光散斑测速（LSV）**：

速度与相关时间的关系：
$v = \delta_s/\tau_c$

其中 $\delta_s$ 是散斑尺寸，$\tau_c$ 是相关时间。

**散斑相关光谱**：

通过散斑的光谱分析获得动态信息：
$G(\tau) = \langle I(t)I(t + \tau) \rangle = \langle I \rangle^2[1 + g_1^2(\tau)]$

### 18.5.3 多样性方法抑制散斑

**角度多样性**：

使用 $N$ 个独立角度照明：
$C = C_0/\sqrt{N}$

实现方法：旋转散射体或多角度照明。

**频率多样性**：

使用带宽 $\Delta\lambda$ 的光源：
相干长度：$lc = \lambda^2/\Delta\lambda$
独立散斑数：$N \approx L/lc$

其中 $L$ 是光程差范围。

**偏振多样性**：

组合正交偏振态：
$I = I_x + I_y$

散斑对比度降低因子：$1/\sqrt{2}$

### 18.5.4 计算机图形学中的散斑模拟

**程序化散斑纹理**：

散斑强度分布：
$I(x, y) = \left|\sum_{n} \exp[i(2\pi\mathbf{f}_n\cdot\mathbf{r} + \varphi_n)]\right|^2$

其中 $\mathbf{f}_n$ 是空间频率，$\varphi_n$ 是随机相位。

**基于物理的散斑渲染**：

1. 表面微结构建模
2. 相干光传播计算
3. 统计特性保持

**实时散斑效果**：

使用预计算的散斑图案：
$I_{\text{final}} = I_{\text{base}} \times (1 + m \times S(u, v))$

其中 $m$ 是调制深度，$S(u, v)$ 是散斑纹理。

### 18.5.5 与体积渲染的统一

**散斑作为体积现象**：

将散斑纳入体积渲染方程：

$L(\mathbf{x}, \mathbf{\omega}) = \int_{0}^{L} \sigma_t(s) \exp\left[-\int_{0}^{s} \sigma_t(t) dt\right] \times [L_e(s) + L_s(s) \times M(s, \mathbf{x})] ds$

其中 $M(s, \mathbf{x})$ 是散斑调制函数。

**相干体积渲染**：

考虑相位的体积渲染：
$U(\mathbf{x}) = \int_{0}^{L} \sigma_s(s) \exp[ikn(s)s] \exp\left[-\int_{0}^{s} \sigma_t(t) dt\right] U_s(s) ds$

**蒙特卡洛散斑**：

路径积分中加入相位：
$\langle L \rangle = \int L(\mathbf{x}) \left|\sum_{\text{paths}} \exp(i\varphi_{\text{path}})\right|^2 d\mathbf{x}$

### 18.5.6 先进散斑控制技术

**自适应光学补偿**：

使用空间光调制器补偿散斑相位：
$U_{\text{corrected}} = U_{\text{speckle}} \times \exp(-i\varphi_{\text{speckle}})$

**计算去散斑**：

基于统计模型的去噪：
$I_{\text{denoised}} = \arg \min_I [||I - I_{\text{observed}}||^2 + \lambda R(I)]$

其中 $R(I)$ 是正则化项。

**深度学习方法**：

训练网络学习散斑到清晰图像的映射：
$I_{\text{clear}} = \text{CNN}(I_{\text{speckle}}, \theta)$

### 18.5.7 未来展望

**量子散斑**：

利用量子纠缠光源的散斑特性：
- 亚散粒噪声成像
- 量子光学相干断层扫描

**计算散斑全息**：

结合散斑与全息技术：
- 散斑场的相位恢复
- 三维散斑重建

## 本章小结

本章建立了统计光学的基本框架，深入分析了散斑现象的物理机理和统计特性。关键要点：

1. **统计描述**：光场的随机性通过一阶和二阶相关函数完整描述，高斯随机过程模型适用于大多数情况
2. **散斑统计**：完全发展散斑的强度服从负指数分布，振幅服从Rayleigh分布；部分发展散斑服从Rice分布
3. **粗糙表面散射**：Kirchhoff近似和相位扰动模型提供了散射场的统计特性，可转化为体积渲染框架
4. **成像系统**：散斑通过成像系统传播时受点扩散函数调制，影响分辨率和信噪比
5. **应用技术**：散斑既可用于精密测量（ESPI、DIC），也可通过多样性方法有效抑制
6. **统一框架**：散斑现象可纳入体积渲染方程，通过相位调制实现物理准确的渲染

## 练习题

### 基础题

**练习18.1**：推导高斯随机光场的强度概率分布
已知复振幅 $U = U_r + iU_i$，其中 $U_r$ 和 $U_i$ 是独立的零均值高斯随机变量，方差均为 $\sigma^2$。证明强度 $I = |U|^2$ 服从参数为 $\langle I \rangle = 2\sigma^2$ 的负指数分布。

*提示*：使用变量变换从 $(U_r, U_i)$ 到 $(I, \varphi)$，其中 $\varphi$ 是相位。

<details>
<summary>答案</summary>

联合概率密度：
$p(U_r, U_i) = \frac{1}{2\pi\sigma^2} \exp\left(-\frac{U_r^2 + U_i^2}{2\sigma^2}\right)$

变换到极坐标：$U = \sqrt{I} \exp(i\varphi)$
雅可比行列式：$|J| = 1/2$

$p(I, \varphi) = p(U_r, U_i)|J| = \frac{1}{4\pi\sigma^2} \exp(-I/2\sigma^2)$

对 $\varphi$ 积分（0 到 2π）：
$p(I) = \frac{1}{2\sigma^2} \exp(-I/2\sigma^2) = \frac{1}{\langle I \rangle} \exp(-I/\langle I \rangle)$

其中 $\langle I \rangle = 2\sigma^2$。
</details>

**练习18.2**：计算散斑对比度
$N$ 个独立的完全发展散斑场叠加，每个场的平均强度为 $I_0$。求合成散斑场的对比度。

*提示*：利用独立随机变量和的方差性质。

<details>
<summary>答案</summary>

每个散斑场：$\langle I_i \rangle = I_0$，$\text{Var}(I_i) = I_0^2$（负指数分布）

总强度：$I = \sum_{i=1}^{N} I_i$
平均值：$\langle I \rangle = N I_0$
方差：$\text{Var}(I) = N I_0^2$（独立性）

对比度：
$C = \sigma_I/\langle I \rangle = \sqrt{N I_0^2}/(N I_0) = 1/\sqrt{N}$
</details>

**练习18.3**：Kirchhoff近似下的镜面反射
使用Kirchhoff近似，计算正入射平面波从均方根粗糙度为 $\sigma_h$ 的高斯粗糙表面的镜面反射系数。

*提示*：计算相干散射分量 $\langle U_s \rangle$。

<details>
<summary>答案</summary>

相位变化：$\Delta\varphi = 2k h(x,y)$（正入射）

反射系数：$r = r_0 \exp(i2kh)$

平均反射系数：
$\langle r \rangle = r_0\langle\exp(i2kh)\rangle = r_0 \exp(-2k^2\langle h^2 \rangle)$
    $= r_0 \exp(-2k^2\sigma_h^2) = r_0 \exp(-g/2)$

其中 $g = (2k\sigma_h)^2$ 是粗糙度因子。
</details>

### 挑战题

**练习18.4**：部分相干光的散斑统计
推导部分相干光（相干度 $|\gamma_{12}| < 1$）产生的散斑对比度。考虑两个部分相干的点源。

*提示*：使用 van Cittert-Zernike 定理和二阶相关函数。

<details>
<summary>答案</summary>

两点源的场：$U = U_1 + U_2$
互相干函数：$\langle U_1^*U_2 \rangle = \sqrt{I_1I_2} \gamma_{12}$

强度：$I = |U_1|^2 + |U_2|^2 + 2\text{Re}(U_1^*U_2)$

平均强度：$\langle I \rangle = I_1 + I_2$
强度方差：$\text{Var}(I) = I_1^2 + I_2^2 + 4I_1I_2|\gamma_{12}|^2$

对于 $I_1 = I_2 = I_0/2$：
$C^2 = \text{Var}(I)/\langle I \rangle^2 = (1 + |\gamma_{12}|^2)/2$

因此：$C = \sqrt{(1 + |\gamma_{12}|^2)/2}$
</details>

**练习18.5**：动态散斑的速度测量极限
分析使用散斑相关技术测量横向速度的理论分辨率极限。给定激光波长 $\lambda$、散斑尺寸 $\delta_s$ 和探测器积分时间 $T$。

*提示*：考虑相关函数衰减到 1/e 的时间。

<details>
<summary>答案</summary>

散斑移动距离 $\delta_s$ 时相关性降到 1/e：
$\tau_c = \delta_s/v$

探测条件：$T < \tau_c$（避免运动模糊）

速度分辨率由相关峰宽度决定：
$\Delta\tau \approx \lambda/(\pi\delta_s) \times (\delta_s/v) = \lambda/(\pi v)$

相对速度分辨率：
$\Delta v/v \approx \lambda/(\pi v T) = \lambda/(\pi\delta_s) \quad (\text{当 } T = \tau_c)$

最小可测速度：$v_{\text{min}} \approx \delta_s/T$
最大可测速度：$v_{\text{max}} \approx \delta_s \times (\text{采样率})$
</details>

**练习18.6**：散斑场的相位恢复
给定散斑强度分布 $I(x,y)$ 和部分相位信息（如边界条件），讨论相位恢复的可能性和唯一性。

*提示*：考虑 Gerchberg-Saxton 算法的收敛性。

<details>
<summary>答案</summary>

相位恢复问题：已知 $|U|^2 = I$，求 $\varphi$ 使 $U = \sqrt{I} \exp(i\varphi)$

约束条件：
1. 强度约束：$|U|^2 = I$（已知）
2. 支撑域约束：$U$ 在某区域外为零
3. 平滑性约束：$\varphi$ 连续可微

迭代算法：
1. 傅里叶域：施加孔径约束
2. 空间域：施加强度约束
3. 重复直到收敛

唯一性：
- 无额外约束时有平移和共轭模糊性
- 过采样因子 > 2 时通常可唯一恢复
- 散斑的随机性有助于打破对称性

收敛性依赖于：
- 初始猜测质量
- 噪声水平
- 约束强度
</details>

## 常见陷阱与错误

1. **混淆散斑尺寸定义**
   - 错误：使用强度半高宽作为散斑尺寸
   - 正确：使用自相关函数积分或第一零点

2. **忽略部分相干效应**
   - 错误：假设所有激光都产生 $C = 1$ 的散斑
   - 正确：考虑光源相干性，多模激光器 $C < 1$

3. **错误的统计假设**
   - 错误：对非完全发展散斑使用 Rayleigh 分布
   - 正确：使用 Rice 分布或实验确定分布

4. **相位扰动模型的误用**
   - 错误：对大粗糙度表面使用相位扰动近似
   - 正确：检查 $k\sigma_h \ll 1$ 的条件

5. **忽略偏振效应**
   - 错误：标量理论处理所有散斑问题
   - 正确：金属表面等需要考虑偏振

## 最佳实践检查清单

### 散斑分析设计审查

- [ ] **统计模型选择**
  - 确认光源相干性（时间和空间）
  - 验证散斑类型（完全/部分发展）
  - 选择合适的概率分布

- [ ] **测量系统设计**
  - 散斑尺寸与探测器分辨率匹配
  - 采样满足 Nyquist 准则
  - 动态范围覆盖强度分布

- [ ] **数值计算考虑**
  - 相关函数的有效计算（FFT方法）
  - 统计量的样本数充足性
  - 边界效应的处理

- [ ] **抑制方法评估**
  - 多样性方法的独立性验证
  - 对比度降低与分辨率损失权衡
  - 实时性要求的满足

- [ ] **应用特定要求**
  - 干涉测量：相位稳定性
  - 速度测量：时间分辨率
  - 成像系统：信噪比优化

- [ ] **与渲染系统集成**
  - 物理准确性 vs 计算效率
  - 相干效应的选择性包含
  - 艺术控制参数的暴露

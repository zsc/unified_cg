# 第2章：光场与近场光学

本章建立了光场的数学框架，从几何光学的4D光场表示扩展到包含相干性和衍射的波动光学描述。我们将探讨光场的各种表示方法，特别是Wigner分布函数如何统一几何光学和波动光学。通过学习菲涅尔和夫琅禾费衍射理论，以及空间和时间相干性概念，我们为后续章节中的高级渲染技术奠定基础。

## 学习目标

完成本章后，您将能够：
1. 理解并推导4D和5D光场的数学表示
2. 掌握全光函数的概念及其在渲染中的应用
3. 使用Wigner分布函数分析光学系统
4. 从第一性原理推导菲涅尔和夫琅禾费衍射公式
5. 量化分析光源的空间和时间相干性
6. 将光场理论应用于计算机图形学问题

## 2.1 4D和5D光场表示

### 2.1.1 光场的几何定义

光场 $L(x, y, \theta, \varphi)$ 描述了空间中每个点 $(x, y)$ 沿每个方向 $(\theta, \varphi)$ 传播的光强度。在自由空间中，这个4D函数完全描述了光的分布。

对于计算机图形学，我们通常使用两平面参数化：
$$L(u, v, s, t)$$

其中 $(u, v)$ 是第一个平面上的坐标，$(s, t)$ 是第二个平面上的坐标。光线由连接这两点的直线定义。

#### 参数化之间的转换

从光线参数 $(x, y, \theta, \varphi)$ 到两平面参数化的转换关系为：
$$u = x - d_1 \tan \theta \cos \varphi$$
$$v = y - d_1 \tan \theta \sin \varphi$$
$$s = x + d_2 \tan \theta \cos \varphi$$
$$t = y + d_2 \tan \theta \sin \varphi$$

其中 $d_1$ 和 $d_2$ 是到两个参考平面的距离。

#### 光场的积分表示

光场可以通过场景的辐射度函数积分得到：
$$L(u, v, s, t) = \int_{\mathcal{S}} \rho(\mathbf{p}) V(\mathbf{p}, \mathbf{r}_{u,v,s,t}) G(\mathbf{p}, \mathbf{r}_{u,v,s,t}) dA(\mathbf{p})$$

其中：
- $\rho(\mathbf{p})$ 是表面点p的反射率
- $V$ 是可见性函数（0或1）
- $G$ 是几何项：$G = \frac{\cos \theta_i \cos \theta_o}{||\mathbf{p} - \mathbf{r}||^2}$

### 2.1.2 5D光场：包含时间和波长

完整的5D光场表示为：
$$L(x, y, \theta, \varphi, t, \lambda)$$

或在两平面参数化下：
$$L(u, v, s, t, \tau, \lambda)$$

其中 $\tau$ 是时间延迟，$\lambda$ 是波长。这种表示对于分析时变场景和色散效应至关重要。

#### 时间维度的物理意义

时间维度允许我们描述：
1. **运动模糊**：快门开启期间的积分
   $$L_{blur}(u, v, s, t) = \frac{1}{T} \int_0^T L(u, v, s, t, \tau) d\tau$$

2. **频闪效应**：周期性照明下的采样
   $$L_{strobe}(u, v, s, t, n) = L(u, v, s, t, nT_s)$$

3. **光脉冲传播**：超快成像中的时间分辨
   $$L_{pulse}(u, v, s, t, \tau) = L_0(u, v, s, t) h(\tau - d/c)$$
   
   其中 $h$ 是脉冲形状函数，$d$ 是传播距离。

#### 光谱维度的重要性

波长维度捕获：
1. **色散**：不同波长的折射率差异
   $$n(\lambda) = A + \frac{B}{\lambda^2} + \frac{C}{\lambda^4} + ...$$ (Cauchy公式)

2. **光谱BRDF**：材质的波长依赖反射
   $$f_r(\mathbf{\omega}_i, \mathbf{\omega}_o, \lambda)$$

3. **荧光效应**：波长转换
   $$L_o(\lambda_o) = \int f_{fluorescence}(\lambda_i, \lambda_o) L_i(\lambda_i) d\lambda_i$$

### 2.1.3 光场的守恒性质

在无损介质中，光场满足亮度守恒：
$$\frac{\partial L}{\partial s} = 0$$

这意味着沿光线的辐射亮度保持不变（在几何光学近似下）。

#### 刘维尔定理的应用

更一般地，光场在相空间中满足刘维尔定理：
$$\frac{dL}{dt} = \frac{\partial L}{\partial t} + \{\mathcal{H}, L\} = 0$$

其中 $\{\cdot,\cdot\}$ 是泊松括号，$\mathcal{H}$ 是哈密顿量。

泊松括号的展开形式：
$$\{f, g\} = \sum_i \left(\frac{\partial f}{\partial q_i}\frac{\partial g}{\partial p_i} - \frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q_i}\right)$$

对于自由传播：
$$\mathcal{H} = c||\mathbf{k}|| = c\sqrt{k_x^2 + k_y^2 + k_z^2}$$

运动方程通过哈密顿正则方程给出：
$$\frac{dq_i}{dt} = \frac{\partial \mathcal{H}}{\partial p_i}, \quad \frac{dp_i}{dt} = -\frac{\partial \mathcal{H}}{\partial q_i}$$

在傍轴近似下（$k_z \gg k_x, k_y$）：
$$\mathcal{H} \approx ck_z + \frac{c}{2k_z}(k_x^2 + k_y^2)$$

这导出了光场传播的基本方程：
$$L(u', v', s', t') = L(u, v, s, t)$$

当光线从 $(u,v)$ 传播到 $(s',t')$ 时。

#### 相空间体积守恒

在哈密顿光学框架下，光场在相空间中的演化保持体积不变：
$$\iiint\int L(x, y, p_x, p_y) dx dy dp_x dp_y = \text{const}$$

其中 $(p_x, p_y) = (n\sin\theta_x, n\sin\theta_y)$ 是光学动量。

这一守恒律导致了重要的光学不变量：
$$n^2 A \Omega = \text{const}$$

其中 $A$ 是光束截面积，$\Omega$ 是立体角，$n$ 是折射率。这就是著名的étendue（光学扩展量）守恒。

推导过程：考虑光束通过光学系统，入射和出射参数满足：
$$n_1^2 A_1 \Omega_1 = n_2^2 A_2 \Omega_2$$

对于小立体角：$\Omega = \pi \sin^2\theta_{max}$

这导出了数值孔径的不变性：
$$n_1 A_1 \text{NA}_1^2 = n_2 A_2 \text{NA}_2^2$$

其中 $\text{NA} = n \sin \theta$ 是数值孔径。

#### 光场采样定理

Chai等人证明了光场的采样要求与场景深度相关：
$$\Delta u \cdot \Delta s \geq \frac{\lambda(z_{max} - z_{min})}{z_{min}}$$

这是光场相机设计的基础约束。

更精确的分析表明，光场的频谱支撑区域是一个双锥：
$$\Omega_{LF} = \{(f_u, f_v, f_s, f_t) : |f_s| \leq \frac{z_{max}}{\lambda}|f_u|, |f_t| \leq \frac{z_{max}}{\lambda}|f_v|\}$$

这导出了最优采样策略：
- 角度分辨率：$\Delta\theta < \frac{\lambda}{D_{aperture}}$
- 空间分辨率：$\Delta x < \frac{\lambda z_{min}}{D_{scene}}$

### 2.1.4 光场的变换与投影

#### 光场的仿射变换

考虑相机的运动，光场经历仿射变换：
$$L'(u', v', s', t') = L(Au + Bs + E, Cv + Dt + F, s, t)$$

其中变换矩阵编码了相机的平移和旋转。

完整的4D变换矩阵形式：
$$\begin{pmatrix} u' \\ v' \\ s' \\ t' \end{pmatrix} = \mathbf{T} \begin{pmatrix} u \\ v \\ s \\ t \end{pmatrix} + \mathbf{d}$$

对于纯平移 $(t_x, t_y, t_z)$：
$$\begin{pmatrix} u' \\ v' \\ s' \\ t' \end{pmatrix} = \begin{pmatrix} 1 & 0 & -t_z/z_0 & 0 \\ 0 & 1 & 0 & -t_z/z_0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} u \\ v \\ s \\ t \end{pmatrix} + \begin{pmatrix} t_x \\ t_y \\ 0 \\ 0 \end{pmatrix}$$

对于旋转角度 $\theta$（绕 z 轴）：
$$\mathbf{T}_{rot} = \begin{pmatrix} 
\cos\theta & \sin\theta & 0 & 0 \\
-\sin\theta & \cos\theta & 0 & 0 \\
0 & 0 & \cos\theta & \sin\theta \\
0 & 0 & -\sin\theta & \cos\theta
\end{pmatrix}$$

#### 光场的投影操作

从4D光场生成2D图像的过程是一个投影操作：

1. **针孔相机**：
   $$I(x, y) = L(x, y, x, y)$$
   选择共线的 $(u, v) = (s, t)$

2. **有限孔径相机**：
   $$I(x, y) = \int\int_{aperture} L(u, v, x, y) A(u, v) du dv$$
   其中 $A(u, v)$ 是孔径函数
   
   对于圆形孔径：$A(u, v) = \text{circ}(\sqrt{u^2 + v^2}/R)$
   
   景深公式：
   $$\text{DOF} = \frac{2Nc(z_f^2 - z_n^2)}{f^2}$$
   其中 $N$ 是 f 数，$c$ 是混淆圆直径，$z_f, z_n$ 是远近焦平面。

3. **光场相机**（聚焦后）：
   $$I(x, y) = \int\int L(u, v, u + \alpha(x-u), v + \alpha(y-v)) du dv$$
   其中 $\alpha$ 控制焦平面深度
   
   焦平面深度与 $\alpha$ 的关系：
   $$z_{focus} = \frac{z_{uv} \cdot z_{st}}{z_{st} - \alpha(z_{st} - z_{uv})}$$

4. **积分成像**（多视点）：
   $$I_n(x, y) = L(u_n, v_n, x, y)$$
   其中 $(u_n, v_n)$ 是第 n 个视点位置

#### 光场的剪切操作

聚焦对应于光场的4D剪切：
$$L_{focused}(u, v, s, t) = L(u, v, s + \alpha u, t + \alpha v)$$

剪切矩阵表示：
$$\mathbf{S}_\alpha = \begin{pmatrix} 
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
\alpha & 0 & 1 & 0 \\
0 & \alpha & 0 & 1
\end{pmatrix}$$

在频域中，这是一个相位调制：
$$\tilde{L}_{focused}(f_u, f_v, f_s, f_t) = \tilde{L}(f_u - \alpha f_s, f_v - \alpha f_t, f_s, f_t)$$

剪切定理的一般形式：
$$\mathcal{F}\{f(x - \alpha y, y)\} = e^{-2\pi i \alpha f_x f_y} \tilde{f}(f_x, f_y)$$

### 2.1.5 光场的稀疏表示

#### 基于表面的光场参数化

对于大多数真实场景，光场具有固有的2D结构（表面约束）：
$$L(u, v, s, t) = \sum_i R_i(u, v, s, t) \delta(d_i(u, v, s, t))$$

其中 $d_i$ 是到第i个表面的符号距离函数，$R_i$ 是表面的反射特性。

#### 光场的低秩分解

利用场景的低维结构，光场可以分解为：
$$L(u, v, s, t) \approx \sum_{k=1}^r \sigma_k U_k(u, v) V_k(s, t)$$

其中 $r$ 是有效秩，通常 $r \ll \min(UV, ST)$。

这种分解的误差界为：
$$||L - L_r||_F \leq \sqrt{\sum_{k=r+1}^{\min(UV,ST)} \sigma_k^2}$$

#### 光场的稀疏傅里叶表示

在频域中，光场能量集中在低频区域。定义稀疏度：
$$\text{Sparsity}(L) = \frac{||\tilde{L}||_0}{||\tilde{L}||_2}$$

对于典型场景，超过99%的能量集中在不到1%的频率分量中。

## 2.2 全光函数

### 2.2.1 七维全光函数

Adelson和Bergen提出的全光函数 $P(x, y, z, \theta, \varphi, \lambda, t)$ 描述了观察者在位置 $(x, y, z)$、时间 $t$、观察方向 $(\theta, \varphi)$ 和波长 $\lambda$ 下感知的光强度。

这个函数可以分解为：
$$P(x, y, z, \theta, \varphi, \lambda, t) = \int L(x', y', z', \theta', \varphi', \lambda, t') \delta(\text{光线约束}) dx'dy'dz'd\theta'd\varphi'$$

光线约束的具体形式：
$$\delta(\text{光线约束}) = \delta(\mathbf{x} - \mathbf{x}' - s\mathbf{\omega}) \delta(\mathbf{\omega} - \mathbf{\omega}')$$

其中 $s$ 是从观察点到场景点的距离。

#### 全光函数的物理解释

全光函数捕获了视觉世界的完整信息：
1. **空间位置** $(x, y, z)$：观察者的3D位置
2. **观察方向** $(\theta, \varphi)$：视线的球面角度
3. **波长** $\lambda$：电磁辐射的光谱成分（380-780 nm 可见光）
4. **时间** $t$：动态场景的时间演化

从信息论角度，全光函数包含了产生任何可能视觉体验所需的所有信息。

维度分析：
- 空间：3D连续流形
- 方向：2D球面 S²
- 光谱：1D正实轴
- 时间：1D实轴
- 总维度：7D流形 $\mathbb{R}^3 \times S^2 \times \mathbb{R}^+ \times \mathbb{R}$

#### 全光函数的积分形式

考虑场景中的光传播，全光函数可以表示为：
$$P(\mathbf{x}, \mathbf{\omega}, \lambda, t) = \int_0^\infty L(\mathbf{x} - s\mathbf{\omega}, \mathbf{\omega}, \lambda, t - s/c) e^{-\int_0^s \sigma_t(\mathbf{x} - s'\mathbf{\omega}, \lambda) ds'} ds$$

其中：
- $\sigma_t$ 是介质的消光系数（吸收 + 散射）
- 积分沿着观察方向 $\mathbf{\omega}$ 进行
- 指数项表示介质中的衰减（Beer-Lambert定律）

对于散射介质，需要加入内散射项：
$$P(\mathbf{x}, \mathbf{\omega}, \lambda, t) = \int_0^\infty \left[L_e + \int_{4\pi} f_p(\mathbf{\omega}', \mathbf{\omega}) L_{in}(\mathbf{\omega}') d\mathbf{\omega}'\right] T(s) ds$$

其中：
- $L_e$ 是自发光项
- $f_p$ 是相位函数
- $T(s) = e^{-\int_0^s \sigma_t ds'}$ 是透射率

### 2.2.2 降维与实际应用

在实际应用中，我们通过各种假设降低全光函数的维度：

#### 常见的降维策略

1. **静态场景** (6D)：消除 $t$
   $$P_{static}(x, y, z, \theta, \varphi, \lambda) = P(x, y, z, \theta, \varphi, \lambda, t_0)$$

2. **单色光** (6D)：消除 $\lambda$  
   $$P_{mono}(x, y, z, \theta, \varphi, t) = \int S(\lambda) P(x, y, z, \theta, \varphi, \lambda, t) d\lambda$$
   
   其中 $S(\lambda)$ 是光谱灵敏度函数。

3. **自由空间** (4D)：降至4D光场
   $$L(u, v, s, t) = P(x(u,s), y(v,t), z, \theta(u,s), \varphi(v,t))$$
   
   在无遮挡空间中，光线参数化简化了表示。

4. **平面约束** (3D)：相机在平面上移动
   $$P_{planar}(x, y, \theta, \varphi) = P(x, y, z_0, \theta, \varphi)$$

#### 渲染技术与全光函数维度

不同渲染技术对应不同的全光函数采样：

| 技术 | 采样维度 | 全光函数约束 |
|------|----------|--------------|
| 光线追踪 | 2D (像素) | 固定相机位置 |
| 光场渲染 | 4D | 自由空间假设 |
| 体积渲染 | 5D | 包含空间位置 |
| 时空渲染 | 6D | 包含时间维度 |

### 2.2.3 与渲染方程的关系

全光函数与渲染方程通过以下关系连接：
$$L_o(\mathbf{x}, \mathbf{\omega}_o) = L_e(\mathbf{x}, \mathbf{\omega}_o) + \int_\Omega f_r(\mathbf{x}, \mathbf{\omega}_i, \mathbf{\omega}_o) L_i(\mathbf{x}, \mathbf{\omega}_i) (\mathbf{\omega}_i \cdot \mathbf{n}) d\mathbf{\omega}_i$$

其中光场 $L$ 对应于入射和出射辐射亮度。

#### 从全光函数推导渲染方程

给定全光函数 $P$，表面点 $\mathbf{x}$ 的入射辐射度为：
$$L_i(\mathbf{x}, \mathbf{\omega}_i) = P(\mathbf{x}, -\mathbf{\omega}_i, \lambda, t)$$

出射辐射度通过BRDF积分得到：
$$L_o(\mathbf{x}, \mathbf{\omega}_o) = L_e(\mathbf{x}, \mathbf{\omega}_o) + \int_\Omega f_r(\mathbf{x}, \mathbf{\omega}_i, \mathbf{\omega}_o) P(\mathbf{x}, -\mathbf{\omega}_i, \lambda, t) (\mathbf{\omega}_i \cdot \mathbf{n}) d\mathbf{\omega}_i$$

这建立了局部着色模型与全局光传输的联系。

#### 路径积分表述

将渲染方程递归展开，得到路径积分：
$$L(\mathbf{x}_0 \to \mathbf{x}_1) = \sum_{n=1}^\infty \int_{\mathcal{P}_n} L_e(\mathbf{x}_n \to \mathbf{x}_{n-1}) \prod_{i=1}^{n-1} f_r(\mathbf{x}_{i+1} \to \mathbf{x}_i \to \mathbf{x}_{i-1}) G(\mathbf{x}_i \leftrightarrow \mathbf{x}_{i+1}) d\mathcal{P}_n$$

其中 $\mathcal{P}_n$ 是长度为 n 的路径空间。

#### 测量方程

相机传感器的响应是全光函数的加权积分：
$$I_{pixel} = \int_A \int_\Omega \int_T \int_\Lambda W(\mathbf{x}, \mathbf{\omega}, t, \lambda) P(\mathbf{x}, \mathbf{\omega}, \lambda, t) d\lambda dt d\mathbf{\omega} dA$$

其中 $W$ 是像素的响应函数，包含：
- 空间滤波器（像素形状）
- 角度滤波器（孔径）
- 时间滤波器（快门）
- 光谱滤波器（色彩滤镜）

### 2.2.4 全光函数的现代应用

#### 神经辐射场 (NeRF)

NeRF本质上是学习场景的5D全光函数近似：
$$F_\theta: (x, y, z, \theta, \varphi) \to (r, g, b, \sigma)$$

其中神经网络 $F_\theta$ 隐式编码了全光函数。

训练过程可以表述为最小化重建误差：
$$\mathcal{L} = \sum_{rays} ||C(r) - \hat{C}(r)||^2$$

其中 $\hat{C}(r)$ 是沿光线 $r$ 的体积渲染：
$$\hat{C}(r) = \int_0^\infty T(t) \sigma(r(t)) c(r(t), d) dt$$

透射率 $T(t) = \exp(-\int_0^t \sigma(r(s)) ds)$ 编码了遮挡关系。

#### 光场相机

光场相机直接采样4D光场：
$$L_{captured}(u, v, s, t) = \int W_{microlens}(u, v, s, t) L_{scene}(u, v, s, t) dudvdsdt$$

微透镜阵列实现了空间-角度的联合采样。

设计参数的权衡：
- 微透镜直径 $d$：角度分辨率 $\propto 1/d$
- 微透镜焦距 $f$：空间分辨率 $\propto f/d$
- 传感器像素尺寸 $p$：总分辨率受限于 min(空间, 角度)

光场相机的调制传递函数（MTF）：
$$\text{MTF}(f_x, f_\theta) = \text{sinc}(d f_x) \cdot \text{sinc}(p f_\theta/f)$$

#### 计算摄影

许多计算摄影技术可以理解为全光函数的特殊采样：
- **HDR成像**：扩展动态范围的$\lambda$维采样
  $$P_{HDR}(x, y, \theta, \varphi) = \sum_i w_i(L) P(x, y, \theta, \varphi, t_i)$$
  其中 $w_i$ 是基于亮度的权重函数

- **光场显微镜**：高分辨率的4D光场采样
  $$P_{micro}(x, y, z, \theta, \varphi) = \int \text{PSF}(x', y', z') P(x-x', y-y', z-z', \theta, \varphi) dx'dy'dz'$$
  
- **飞行时间成像**：利用时间维度测量深度
  $$d(x, y) = \frac{c}{2} \arg\max_\tau \{P(x, y, \tau) \star h(\tau)\}$$
  其中 $h(\tau)$ 是发射脉冲形状

### 2.2.5 全光函数的信息理论分析

#### 信息容量

全光函数的信息容量受物理约束限制：
$$I_{max} = \frac{V_{scene} \cdot \Omega_{view} \cdot T \cdot B}{\lambda^3 \cdot c}$$

其中：
- $V_{scene}$：场景体积
- $\Omega_{view}$：观察立体角
- $T$：时间窗口
- $B$：光谱带宽

对于典型室内场景（10m³，4$\pi$立体角，1秒，可见光谱）：
$$I_{max} \approx 10^{25} \text{ bits}$$

#### 压缩与冗余

真实场景的全光函数具有大量冗余：
1. **空间相干性**：相邻光线相似
2. **时间连续性**：运动平滑
3. **光谱相关性**：颜色通道相关
4. **几何约束**：满足极线几何

有效压缩比可达 $10^6:1$ 而不显著损失感知质量。

#### 采样理论的信息论视角

根据率失真理论，给定失真容限 $D$，所需采样率：
$$R(D) \geq H(P) - H(D)$$

其中 $H(P)$ 是全光函数的熵，$H(D)$ 是容许失真的熵。

这导出了自适应采样策略：
- 高纹理区域：密集采样
- 平滑区域：稀疏采样
- 遮挡边界：超采样

## 2.3 光学中的Wigner分布函数

### 2.3.1 Wigner分布的定义

对于光场 $U(x)$，其Wigner分布函数定义为：
$$W(x, k) = \int U^*(x - \xi/2) U(x + \xi/2) e^{-ik\xi} d\xi$$

其中 $x$ 是位置，$k$ 是空间频率（与传播角度相关）。

等价的频域表示：
$$W(x, k) = \frac{1}{2\pi} \int \tilde{U}^*(k - \kappa/2) \tilde{U}(k + \kappa/2) e^{ix\kappa} d\kappa$$

对于2D光场，Wigner分布推广为：
$$W(x, y, k_x, k_y) = \int\int U^*(x - \xi_x/2, y - \xi_y/2) U(x + \xi_x/2, y + \xi_y/2) e^{-i(k_x \xi_x + k_y \xi_y)} d\xi_x d\xi_y$$

#### 从量子力学到光学的类比

Wigner分布最初在量子力学中引入，用于相空间中的准概率分布。在光学中：
- 位置 $x \leftrightarrow$ 横向坐标
- 动量 $p \leftrightarrow$ 横向波矢 $k = k_0\sin \theta$（$k_0 = 2\pi/\lambda$）
- 波函数 $\psi \leftrightarrow$ 光场振幅 $U$
- 普朗克常数 $\hbar \leftrightarrow$ 波长 $\lambda/2\pi$

对应关系：
$$W_{quantum}(x, p) \leftrightarrow W_{optics}(x, k)$$
$[x, p] = i\hbar \leftrightarrow \Delta x \Delta k \geq 1/2$

#### Wigner分布的性质

1. **实值性**：$W(x, k) \in \mathbb{R}$，即使对复光场
   证明：$W^*(x, k) = W(x, k)$

2. **边缘分布**：
   $$\int W(x, k) dk = |U(x)|^2$$ （强度分布）
   $$\int W(x, k) dx = |\tilde{U}(k)|^2$$ （角谱分布）
   
3. **非正定性**：$W(x, k)$ 可以为负，因此是准概率分布
   负值区域表示量子/波动干涉效应

4. **归一化**：$$\int\int W(x, k) dx dk = \int |U(x)|^2 dx$$ （总功率）

5. **平移不变性**：
   若 $U'(x) = U(x - x_0)$，则 $W'(x, k) = W(x - x_0, k)$

6. **调制性质**：
   若 $U'(x) = U(x)e^{ik_0 x}$，则 $W'(x, k) = W(x, k - k_0)$

7. **尺度变换**：
   若 $U'(x) = \sqrt{a}U(ax)$，则 $W'(x, k) = W(ax, k/a)$

### 2.3.2 相空间表示

Wigner分布提供了光场在位置-角度相空间中的表示。对于4D光场，我们有：
$$W(x, y, k_x, k_y) = \int\int L(x - \xi_x/2, y - \xi_y/2, x + \xi_x/2, y + \xi_y/2) e^{-i(k_x \xi_x + k_y \xi_y)} d\xi_x d\xi_y$$

#### 相空间的几何解释

在相空间 $(x, k)$ 中：
- **点**：表示具有确定位置和传播方向的光线
- **水平线**：平面波（确定的 $k$，所有 $x$）
- **垂直线**：点源（确定的 $x$，所有 $k$）
- **倾斜线**：会聚/发散球面波

#### 相空间体积与信息容量

根据相空间不确定性原理：
$$\Delta x \cdot \Delta k \geq \frac{1}{2}$$

这限制了光学系统的信息容量。对于有限孔径系统：
$$N_{DOF} \approx \frac{A_{object} \cdot \Omega_{NA}}{\lambda^2}$$

其中 $A_{object}$ 是物体面积，$\Omega_{NA}$ 是数值孔径对应的立体角。

### 2.3.3 传播与变换

通过自由空间传播距离 $z$ 后，Wigner分布变换为：
$$W'(x, k) = W(x - zk/k_0, k)$$

这是一个相空间中的剪切变换，其中 $k_0 = 2\pi/\lambda$。

#### ABCD矩阵形式

对于一般的傍轴光学系统，Wigner分布通过线性正则变换：
$$\begin{pmatrix} x' \\ k' \end{pmatrix} = \begin{pmatrix} A & B \\ C & D \end{pmatrix} \begin{pmatrix} x \\ k \end{pmatrix}$$

其中 ABCD 是光学系统的传输矩阵：
$$\begin{pmatrix} x' \\ k' \end{pmatrix} = \begin{pmatrix} A & B \\ C & D \end{pmatrix} \begin{pmatrix} x \\ k \end{pmatrix}$$

常见光学元件的ABCD矩阵：
- 自由传播 $z$：$\begin{pmatrix} 1 & z/k_0 \\ 0 & 1 \end{pmatrix}$
- 薄透镜 $f$：$\begin{pmatrix} 1 & 0 \\ -k_0/f & 1 \end{pmatrix}$
- 傅里叶变换：$\begin{pmatrix} 0 & 1/k_0 \\ -k_0 & 0 \end{pmatrix}$

#### 分数傅里叶变换

Wigner分布在分数傅里叶变换下旋转角度 $\alpha$：
$$W_\alpha(x_\alpha, k_\alpha) = W(x\cos \alpha - k\sin \alpha/k_0, k_0(x\sin \alpha + k\cos \alpha/k_0))$$

这提供了相空间中的统一描述：
- $\alpha = 0$：恒等变换
- $\alpha = \pi/2$：标准傅里叶变换
- $\alpha = \pi$：反演
- $0 < \alpha < \pi/2$：分数域

### 2.3.4 与光线追踪的联系

在几何光学极限下，Wigner分布退化为光线的相空间表示：
$$W(x, k) \to \sum_i I_i \delta(x - x_i) \delta(k - k_i)$$

#### 从波动到几何的过渡

考虑高斯光束的Wigner分布：
$$W(x, k) = \frac{1}{\pi} \exp\left[-\frac{x^2}{w^2} - \frac{k^2w^2}{k_0^2}\right]$$

其中 $w$ 是束腰半径。当 $w \gg \lambda$ 时：
- 相空间分布变得高度局域化
- 可以用光线束近似
- 衍射效应可忽略

#### 光线光学中的相空间

在光线光学中，Wigner分布简化为：
$$W_{ray}(x, \theta) = \sum_i I_i(s) \delta(x - x_i(s)) \delta(\theta - \theta_i(s))$$

其中 $s$ 是沿光线的参数。这导出了：
- **光程函数**：$S(x) = \int n(s) ds$
- **程函方程**：$|\nabla S|^2 = n^2$
- **光线方程**：$\frac{d}{ds}\left(n\frac{d\mathbf{r}}{ds}\right) = \nabla n$

### 2.3.5 Wigner分布的计算方法

#### 直接计算

从光场$U(x)$计算Wigner分布：
```
1. 对每个位置x和频率k：
   - 计算U*(x - ξ/2)U(x + ξ/2)
   - 乘以exp(-ikξ)
   - 对ξ积分
2. 处理边界效应（零填充或周期延拓）
```

#### 基于FFT的快速算法

利用Wigner分布与模糊函数的关系：
$$W(x, k) = \mathcal{F}_\xi\{U^*(x - \xi/2)U(x + \xi/2)\}$$

计算复杂度：O(N²log N)，其中N是采样点数。

#### 相空间层析

从多个投影重建Wigner分布：
$$R_\theta(t) = \int W(t\cos \theta - s\sin \theta, t\sin \theta + s\cos \theta) ds$$

这是Radon变换，可用于从强度测量重建相位信息。

## 2.4 菲涅尔和夫琅禾费衍射

### 2.4.1 惠更斯-菲涅尔原理

从波动方程出发，光场传播满足：
$$U(P) = \frac{1}{i\lambda} \int\int_\Sigma U(Q) \frac{e^{ikr}}{r} \cos(\mathbf{n}, \mathbf{r}) d\Sigma$$

其中 $P$ 是观察点，$Q$ 是孔径上的点，$r = |P - Q|$ 是它们之间的距离。

惠更斯原理的数学表述：
- 每个波前元素都是次级球面波源
- 后续波前是所有次级波的包络
- 菲涅尔贡献：引入了相位和倾斜因子

倾斜因子的物理意义：
$$K(\chi) = \frac{1 + \cos \chi}{2}$$
其中 $\chi$ 是法向与观察方向的夹角。

#### 基尔霍夫衍射理论

从标量波动方程出发：
$$\nabla^2 U + k^2 U = 0$$

亥姆霍兹方程的格林函数解：
$$G(\mathbf{r}, \mathbf{r}') = \frac{e^{ik|\mathbf{r} - \mathbf{r}'|}}{4\pi|\mathbf{r} - \mathbf{r}'|}$$

应用格林定理，得到基尔霍夫积分：
$$U(P) = \frac{1}{4\pi} \int\int_\Sigma \left[U\frac{\partial G}{\partial n} - G\frac{\partial U}{\partial n}\right] d\Sigma$$

展开后：
$$U(P) = \frac{1}{i\lambda} \int\int_\Sigma U(Q) \frac{e^{ikr}}{r} \frac{1 + \cos \chi}{2} d\Sigma$$

基尔霍夫边界条件（在孔径$\Sigma$上）：
- $U = U_{incident}$（孔径内）
- $U = 0$（屏幕上）
- $\partial U/\partial n = \partial U_{incident}/\partial n$（孔径内）
- $\partial U/\partial n = 0$（屏幕上）

#### Rayleigh-Sommerfeld衍射公式

第一类Rayleigh-Sommerfeld公式（更精确）：
$$U(P) = \frac{1}{i\lambda} \int\int_\Sigma U(Q) \frac{e^{ikr}}{r} \cos(\mathbf{n}, \mathbf{r}) d\Sigma$$

推导使用了更合理的边界条件：
- 仅指定U或$\partial U/\partial n$之一
- 避免了基尔霍夫理论的数学不一致性

第二类公式：
$$U(P) = -\frac{1}{i\lambda} \int\int_\Sigma \frac{\partial U}{\partial n} \frac{e^{ikr}}{r} d\Sigma$$

两种公式的选择：
- 第一类：已知孔径平面的场分布U
- 第二类：已知孔径平面的场梯度$\partial U/\partial n$

### 2.4.2 菲涅尔衍射

在菲涅尔近似下（近场），传播核为：
$$h(x, y; x', y') = \frac{e^{ikz}}{i\lambda z} \exp\left[\frac{ik}{2z}[(x-x')^2 + (y-y')^2]\right]$$

衍射场为：
$$U(x, y) = \int\int U_0(x', y') h(x, y; x', y') dx' dy'$$

#### 菲涅尔近似的条件

菲涅尔近似成立的条件：
$$z^3 \gg \frac{\pi}{4\lambda}[(x-x')^2 + (y-y')^2]_{max}^2$$

这确保了相位误差小于 $\pi/8$。

#### 例1：圆孔的菲涅尔衍射

考虑半径为 $a$ 的圆孔，入射平面波。轴上场分布：
$$U(0, 0, z) = U_0 \left[1 - e^{ika^2/2z}\right]$$

强度分布：
$$I(0, 0, z) = 4I_0 \sin^2\left(\frac{ka^2}{4z}\right)$$

菲涅尔区数：$N = \frac{a^2}{\lambda z}$

当 $N$ 为奇数时，中心为亮点；$N$ 为偶数时，中心为暗点。

#### 例2：直边的菲涅尔衍射

半平面屏的菲涅尔衍射，使用菲涅尔积分：
$$U(x) = \frac{1+i}{2} \left[C(v) + iS(v)\right]$$

其中菲涅尔参数：$v = x\sqrt{\frac{2}{\lambda z}}$

菲涅尔积分：
$$C(v) = \int_0^v \cos\left(\frac{\pi t^2}{2}\right) dt$$
$$S(v) = \int_0^v \sin\left(\frac{\pi t^2}{2}\right) dt$$

#### 例3：菲涅尔波带片

菲涅尔波带片的透过率函数：
$$t(r) = \frac{1}{2}\left[1 + \cos\left(\frac{\pi r^2}{\lambda f}\right)\right]$$

这相当于一个衍射透镜，焦距为 $f$。主焦点强度：
$$I_{focus} = \left(\frac{\pi a^2}{\lambda f}\right)^2 I_0$$

### 2.4.3 夫琅禾费衍射

在远场（夫琅禾费区域），衍射简化为傅里叶变换：
$$U(x, y) = \frac{e^{ikz}}{i\lambda z} e^{ik(x^2+y^2)/2z} \mathcal{F}\{U_0\}\left(\frac{x}{\lambda z}, \frac{y}{\lambda z}\right)$$

其中 $\mathcal{F}$ 表示傅里叶变换。

#### 夫琅禾费近似条件

夫琅禾费近似要求：
$$z \gg \frac{k(x'^2 + y'^2)_{max}}{2} = \frac{\pi a^2}{\lambda}$$

这确保孔径上的二次相位项可忽略。

#### 例4：矩形孔径的夫琅禾费衍射

矩形孔径 (宽度 $a \times b$) 的衍射图样：
$$U(x, y) = U_0 ab \frac{e^{ikz}}{i\lambda z} \text{sinc}\left(\frac{ax}{\lambda z}\right) \text{sinc}\left(\frac{by}{\lambda z}\right)$$

其中 $\text{sinc}(x) = \frac{\sin(\pi x)}{\pi x}$。

强度分布：
$$I(x, y) = I_0 \left(\frac{ab}{\lambda z}\right)^2 \text{sinc}^2\left(\frac{ax}{\lambda z}\right) \text{sinc}^2\left(\frac{by}{\lambda z}\right)$$

第一暗环位置：$x = \pm\frac{\lambda z}{a}$，$y = \pm\frac{\lambda z}{b}$

#### 例5：圆孔的夫琅禾费衍射（艾里斑）

圆孔（半径 $a$）的衍射：
$$U(r, \theta) = U_0 \frac{\pi a^2}{i\lambda z} e^{ikz} \frac{2J_1(kar/z)}{kar/z}$$

其中 $J_1$ 是一阶贝塞尔函数。

艾里斑特征：
- 中心最大值：$I_0 = \left(\frac{\pi a^2}{\lambda z}\right)^2 I_{incident}$
- 第一暗环：$r = 1.22\frac{\lambda z}{2a}$ (数值孔径 $\text{NA} = a/z$)
- 中心亮斑包含总功率的 83.8%

#### 例6：光栅的夫琅禾费衍射

N 条缝的光栅（缝宽 $a$，周期 $d$）：
$$I(\theta) = I_0 \left(\frac{\sin(N\pi d\sin \theta/\lambda)}{N\sin(\pi d\sin \theta/\lambda)}\right)^2 \left(\frac{\sin(\pi a\sin \theta/\lambda)}{\pi a\sin \theta/\lambda}\right)^2$$

主极大位置（光栅方程）：$d\sin \theta_m = m\lambda$

分辨本领：$R = mN$（$m$ 是衍射级次）

### 2.4.4 数值计算考虑

菲涅尔数 $F = a^2/(\lambda z)$ 决定了使用哪种近似：
- $F \gg 1$：几何光学
- $F \approx 1$：菲涅尔衍射
- $F \ll 1$：夫琅禾费衍射

#### 衍射积分的数值方法

1. **直接积分法**：
   - 适用于任意孔径形状
   - 计算复杂度：O(N_output $\times$ N_aperture)

2. **FFT方法**（夫琅禾费）：
   - 利用衍射是傅里叶变换
   - 计算复杂度：O(N log N)

3. **角谱法**：
   $$U(x, y, z) = \mathcal{F}^{-1}\left\{\mathcal{F}\{U_0\} \exp\left[ikz\sqrt{1-(\lambda f_x)^2-(\lambda f_y)^2}\right]\right\}$$
   - 精确，无近似
   - 适用于近场传播

4. **菲涅尔卷积法**：
   - 使用卷积定理
   - 需要适当的采样以避免混叠

#### 采样要求

空间采样间隔：
- 菲涅尔衍射：$\Delta x < \sqrt{\frac{\lambda z}{N}}$（N是采样点数）
- 夫琅禾费衍射：$\Delta x < \frac{\lambda z}{L}$（L是观察区域尺寸）

避免混叠的条件：
$$\Delta x_{aperture} \times \Delta x_{observation} \leq \frac{\lambda z}{N}$$

## 2.5 空间和时间相干性

### 2.5.1 互相干函数

两点间的互相干函数定义为：
$$\Gamma(\mathbf{r}_1, \mathbf{r}_2, \tau) = \langle U^*(\mathbf{r}_1, t) U(\mathbf{r}_2, t + \tau) \rangle$$

归一化后得到复相干度：
$$\gamma(\mathbf{r}_1, \mathbf{r}_2, \tau) = \frac{\Gamma(\mathbf{r}_1, \mathbf{r}_2, \tau)}{\sqrt{I(\mathbf{r}_1)I(\mathbf{r}_2)}}$$

### 2.5.2 时间相干性

时间相干性由相干时间 $\tau_c$ 和相干长度 $l_c = c\tau_c$ 表征。对于光谱宽度 $\Delta\nu$ 的光源：
$$\tau_c \approx \frac{1}{\Delta\nu}$$

时间相干性影响干涉条纹的可见度：
$$V = |\gamma(\tau)| = |\mathcal{F}\{S(\nu)\}|$$

其中 $S(\nu)$ 是归一化光谱密度。

### 2.5.3 空间相干性

空间相干性由相干面积表征。根据van Cittert-Zernike定理，扩展非相干源产生的场的空间相干性为：
$$\gamma(\mathbf{r}_1, \mathbf{r}_2) = \frac{\int\int I(\xi, \eta) e^{ik[(\mathbf{r}_1-\mathbf{r}_2)\cdot(\xi,\eta)]/z} d\xi d\eta}{\int\int I(\xi, \eta) d\xi d\eta}$$

对于均匀圆形源，相干半径为：
$$\rho_c \approx 1.22 \frac{\lambda z}{D}$$

其中 $D$ 是源的直径，$z$ 是观察距离。

### 2.5.4 相干性与渲染

在计算机图形学中，相干性影响：
- 散斑图案的形成
- 薄膜干涉的可见度
- 全息显示的质量
- 激光投影系统的设计

## 本章小结

本章建立了光场理论的数学基础，主要概念包括：

1. **4D/5D光场表示**：$L(u, v, s, t, \lambda)$ 完整描述光的空间、角度和光谱分布
2. **全光函数**：七维函数 $P(x, y, z, \theta, \varphi, \lambda, t)$ 描述所有可能的视觉信息
3. **Wigner分布**：$W(x, k)$ 统一了几何光学和波动光学的相空间表示
4. **衍射理论**：菲涅尔数 $F = a^2/(\lambda z)$ 决定使用几何、菲涅尔或夫琅禾费近似
5. **相干性**：时间相干长度 $l_c \approx c/\Delta\nu$，空间相干半径 $\rho_c \approx 1.22\lambda z/D$

这些概念为理解现代渲染技术（如光场相机、全息显示）提供了理论基础。

## 练习题

### 基础题

**练习 2.1**：证明两平面参数化的4D光场在自由空间传播时满足亮度守恒。  
*提示*：使用光线的参数方程和雅可比行列式。

<details>
<summary>答案</summary>

设光线从平面1 $(u, v)$ 传播到平面2 $(s, t)$，距离为 $d$。光线方程为：
$$s = u + d \cdot \tan \theta_x, \quad t = v + d \cdot \tan \theta_y$$

雅可比行列式：
$$J = \left|\frac{\partial(s,t)}{\partial(u,v)}\right| = 1$$

因此 $L(u, v, s, t) = L(u', v', s', t')$，证明了亮度守恒。
</details>

**练习 2.2**：对于直径 $D = 1 \text{ mm}$ 的圆形非相干光源，波长 $\lambda = 500 \text{ nm}$，计算在距离 $z = 1 \text{ m}$ 处的空间相干半径。  
*提示*：使用van Cittert-Zernike定理的结果。

<details>
<summary>答案</summary>

使用公式 $\rho_c \approx 1.22\lambda z/D$：
$$\rho_c = 1.22 \times \frac{500 \times 10^{-9} \times 1}{10^{-3}} = 0.61 \text{ mm}$$

这意味着相距超过0.61 mm的两点基本不相干。
</details>

**练习 2.3**：推导菲涅尔数 $F = 1$ 时的衍射场强度分布。  
*提示*：考虑圆孔的菲涅尔衍射积分。

<details>
<summary>答案</summary>

对于半径 $a$ 的圆孔，$F = a^2/(\lambda z) = 1$ 意味着 $z = a^2/\lambda$。

轴上强度：
$$I(0) = I_0 |1 - e^{ika^2/2z}|^2 = I_0 |1 - e^{i\pi}|^2 = 4I_0$$

这是几何阴影强度的4倍，展示了菲涅尔衍射的聚焦效应。
</details>

### 挑战题

**练习 2.4**：推导Wigner分布函数通过薄透镜的变换规律。透镜焦距为 $f$。  
*提示*：薄透镜相位变换为 $\exp[-ik(x^2 + y^2)/2f]$。

<details>
<summary>答案</summary>

薄透镜引入二次相位：
$$U'(x) = U(x) \exp\left[-\frac{ik x^2}{2f}\right]$$

Wigner分布变换为：
$$W'(x, k) = W(x, k + \frac{x}{f})$$

这是相空间中的另一种剪切变换，与自由传播的剪切方向正交。组合得到ABCD矩阵光学。
</details>

**练习 2.5**：证明部分相干光的4D光场表示需要考虑互相干函数 $\Gamma(\mathbf{r}_1, \mathbf{r}_2)$。推导相应的广义光场方程。  
*提示*：从统计光学出发，考虑场的二阶相关。

<details>
<summary>答案</summary>

部分相干光场不能用确定性函数$L(u, v, s, t)$描述。需要引入互强度：
$$J(u_1, v_1, u_2, v_2) = \langle L^*(u_1, v_1, s, t) L(u_2, v_2, s, t) \rangle$$

传播方程：
$$J(s_1, t_1, s_2, t_2) = \int\int\int\int h^*(s_1-u_1, t_1-v_1) h(s_2-u_2, t_2-v_2) J(u_1, v_1, u_2, v_2) du_1dv_1du_2dv_2$$

其中$h$是传播核。这导致光场维度翻倍：4D $\to$ 8D。
</details>

**练习 2.6**：设计一个算法，从多视角图像重建4D光场。分析采样要求和重建误差。  
*提示*：考虑极线几何约束和频域分析。

<details>
<summary>答案</summary>

算法框架：
1. 极线约束：对应点满足 $\mathbf{x}'^\top\mathbf{F}\mathbf{x} = 0$
2. 视差与深度：$d = f \cdot B/Z$
3. 光场重建：$L(u, v, s, t) = I(u + d(s-u)/z, v + d(t-v)/z)$

采样要求（Nyquist）：
- 角度采样：$\Delta\theta < \lambda/D$（D是场景尺度）
- 空间采样：$\Delta x < \lambda z/D$（z是深度范围）

重建误差主要来自：
- 遮挡区域：O(遮挡面积/总面积)
- 深度误差：$\delta L/L \approx \delta Z/Z$
- 混叠：当违反Nyquist条件时
</details>

### 开放性思考题

**练习 2.7**：讨论如何将传统光线追踪算法扩展以处理波动光学效应。考虑计算复杂度和精度权衡。  
*提示*：思考混合几何-波动方法。

<details>
<summary>答案</summary>

可能的方法：
1. **局部波动修正**：在需要时（如边缘、小孔）切换到波动计算
2. **相位光线追踪**：追踪复振幅而非仅强度
3. **Wigner函数光线追踪**：在相空间追踪，自然包含衍射

复杂度分析：
- 纯几何：O(N log N)（N是光线数）
- 局部波动：O(N log N + M$\cdot$K²)（M是衍射区域，K是采样点）
- 全波动：O(N³)（体积离散化）

建议：根据菲涅尔数自适应选择方法。
</details>

**练习 2.8**：探讨机器学习如何用于光场压缩和重建。设计一个网络架构并分析其理论基础。  
*提示*：考虑光场的低秩结构和稀疏性。

<details>
<summary>答案</summary>

网络架构建议：
1. **编码器**：4D卷积提取光场特征
2. **低秩分解**：类似TensoRF的向量-矩阵分解
3. **解码器**：上采样重建完整光场

理论基础：
- 光场在无遮挡区域是低秩的（秩 $\leq$ 场景复杂度）
- 傅里叶域稀疏性（大部分能量集中在低频）
- 极线一致性提供额外约束

损失函数：
$$\mathcal{L} = ||L - L_{gt}||_2 + \lambda_1||\nabla L||_1 + \lambda_2 \mathcal{L}_{epi}$$

其中 $\mathcal{L}_{epi}$ 是极线一致性损失。
</details>

## 常见陷阱与错误

### 陷阱 1：混淆4D光场的不同参数化
**错误**：直接在 $(x, y, \theta, \varphi)$ 和 $(u, v, s, t)$ 之间转换  
**正确**：需要考虑坐标系和采样结构的差异

### 陷阱 2：忽略相干性在干涉中的作用
**错误**：假设所有光源都能产生稳定干涉  
**正确**：只有相干长度内的光程差才能产生干涉

### 陷阱 3：错误使用菲涅尔近似
**错误**：在近场直接使用夫琅禾费公式  
**正确**：检查菲涅尔数F，选择合适的近似

### 陷阱 4：Wigner分布的物理解释
**错误**：将$W(x, k)$解释为位置x处沿方向k的强度  
**正确**：W可以为负，是准概率分布

### 陷阱 5：光场采样不足
**错误**：使用稀疏相机阵列重建光场  
**正确**：满足角度和空间Nyquist条件

## 最佳实践检查清单

### 光场表示选择
- [ ] 根据应用选择合适的参数化（光线角度 vs 两平面）
- [ ] 考虑是否需要时间和光谱维度
- [ ] 评估存储和计算需求

### 衍射计算
- [ ] 计算菲涅尔数确定适用的近似
- [ ] 验证采样满足Nyquist条件
- [ ] 考虑边界效应和零填充

### 相干性分析
- [ ] 识别光源类型（激光/LED/热光源）
- [ ] 计算相干长度和相干面积
- [ ] 在设计光学系统时考虑相干性限制

### 数值实现
- [ ] 使用FFT加速卷积计算
- [ ] 注意相位展开和2$\pi$模糊
- [ ] 验证能量守恒

### 算法优化
- [ ] 利用光场的稀疏性和低秩性
- [ ] 考虑空间-角度分辨率权衡
- [ ] 实现自适应采样策略

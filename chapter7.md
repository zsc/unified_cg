# 第7章：动态神经辐射场

现实世界本质上是动态的。虽然静态神经辐射场（NeRF）彻底改变了3D场景表示，但将其扩展以处理时间变化为建模从人类表演到流体动力学的一切事物开辟了可能性。本章探讨了如何将时间作为神经辐射场中的一个额外维度，转换体渲染方程以处理动态场景，同时保持我们统一框架的数学优雅。

动态场景带来了独特的挑战：在没有混叠的情况下捕捉快速运动，保持时间连贯性，处理拓扑变化，以及在将时间作为第四维度添加时管理数据复杂性的指数增长。我们将看到，仔细的数学公式，结合物理学和信号处理的见解，如何带来实用且高效的解决方案。

## 学习目标

完成本章后，您将能够：

1.  **公式化**动态场景的时间扩展体渲染方程
2.  **设计**在规范空间和观测空间之间映射的变形场
3.  **从**体渲染框架中**推导**光学流和场景流约束
4.  **分析**不同时间表示的计算和内存权衡
5.  **实现**时间连贯重建的正则化策略
6.  **评估**动态神经辐射场的收敛特性

## 先决条件

本章建立在以下内容之上：
- 第6章：神经辐射场（静态公式）
- 向量微积分的理解（雅可比矩阵，散度定理）
- 熟悉优化和正则化技术
- 连续介质力学基础知识（有帮助但非必需）

## 7.1 时变辐射场表示

动态辐射场要求我们建模几何（密度 $\sigma$）和外观（颜色 $\mathbf{c}$）如何随时间演变。关键的见解是，虽然数学框架保持优雅，但时间表示的选择深刻影响表达能力和计算效率。

### 7.1.1 扩展体渲染方程

第6章中的静态体渲染方程：

$$C(\mathbf{r}) = \int_{t_n}^{t_f} T(t)\sigma(\mathbf{r}(t))\mathbf{c}(\mathbf{r}(t), \mathbf{d}) dt$$

其中 $T(t) = \exp\left(-\int_{t_n}^{t} \sigma(\mathbf{r}(s))ds\right)$

对于动态场景，我们引入时间 $\tau$ 作为附加参数（注意：我们使用 $\tau$ 表示场景时间以区别于射线参数 $t$）：

$$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T(t, \tau)\sigma(\mathbf{r}(t), \tau)\mathbf{c}(\mathbf{r}(t), \mathbf{d}, \tau) dt$$

透射率也变得与时间相关：

$$T(t, \tau) = \exp\left(-\int_{t_n}^{t} \sigma(\mathbf{r}(s), \tau)ds\right)$$

这个看似简单的扩展具有深远的影响：

1.  **维度增加**：辐射场现在是 $f: \mathbb{R}^3 \times \mathbb{S}^2 \times \mathbb{R} \rightarrow \mathbb{R}^4$
2.  **时间连贯性**：相邻时间帧应产生相似的辐射值
3.  **运动模糊**：曝光时间内的快速运动需要对 $\tau$ 进行积分

**时间扩展场的数学性质：**

时间扩展辐射场必须满足几个连续性条件以实现物理合理性：

1.  **时间连续性**：对于平滑运动，我们要求：
    $$\lim_{\Delta\tau \to 0} \|\sigma(\mathbf{x}, \tau + \Delta\tau) - \sigma(\mathbf{x}, \tau)\| = 0$$

2.  **时空耦合**：梯度张量：
    $$\nabla_{(\mathbf{x},\tau)} f = \begin{bmatrix} \nabla_\mathbf{x} f \\ \frac{\partial f}{\partial \tau} \end{bmatrix}$$
    
    捕捉了空间和时间变化如何相互作用。

3.  **守恒定律**：对于物理场景，密度变化必须遵守质量守恒：
    $$\frac{\partial \rho}{\partial \tau} + \nabla \cdot (\rho \mathbf{v}) = 0$$
    
    其中 $\rho = \sigma$（将密度视为质量密度）且 $\mathbf{v}$ 是速度场。

**运动模糊公式：**

为了实现真实的渲染，我们必须考虑有限的曝光时间。观察到的辐射变为时间积分：

$$C_{\text{blur}}(\mathbf{r}) = \frac{1}{\Delta\tau} \int_{\tau_0}^{\tau_0+\Delta\tau} C(\mathbf{r}, \tau) d\tau$$

用体渲染方程展开：

$$C_{\text{blur}}(\mathbf{r}) = \frac{1}{\Delta\tau} \int_{\tau_0}^{\tau_0+\Delta\tau} \int_{t_n}^{t_f} T(t, \tau)\sigma(\mathbf{r}(t), \tau)\mathbf{c}(\mathbf{r}(t), \mathbf{d}, \tau) dt d\tau$$

根据富比尼定理（在适当的可积性条件下），我们可以交换积分顺序：

$$C_{\text{blur}}(\mathbf{r}) = \int_{t_n}^{t_f} \left[\frac{1}{\Delta\tau} \int_{\tau_0}^{\tau_0+\Delta\tau} T(t, \tau)\sigma(\mathbf{r}(t), \tau)\mathbf{c}(\mathbf{r}(t), \mathbf{d}, \tau) d\tau\right] dt$$

这个公式揭示了运动模糊可以被解释为使用时间平均密度和颜色场进行渲染。

**快门函数：**

真实相机没有随时间均匀的曝光。我们用快门函数 $s(\tau)$ 进行泛化：

$$C_{\text{blur}}(\mathbf{r}) = \int_{\tau_0}^{\tau_0+\Delta\tau} s(\tau) C(\mathbf{r}, \tau) d\tau$$

其中 $\int s(\tau) d\tau = 1$。常见的快门函数包括：

1.  **盒式快门**（理想）：$s(\tau) = \frac{1}{\Delta\tau}\mathbb{1}_{[\tau_0, \tau_0+\Delta\tau]}(\tau)$
2.  **三角快门**：$s(\tau) = \frac{2}{\Delta\tau}\left(1 - \frac{2|\tau - \tau_0 - \Delta\tau/2|}{\Delta\tau}\right)$
3.  **高斯快门**：$s(\tau) = \frac{1}{\sqrt{2\pi}\sigma_s}\exp\left(-\frac{(\tau - \tau_0 - \Delta\tau/2)^2}{2\sigma_s^2}\right)$

快门函数的选择影响运动模糊特性，可以从数据中学习或为艺术效果而设计。

### 7.1.2 离散时间与连续时间建模

离散时间表示和连续时间表示之间的选择涉及表达能力、内存效率和优化复杂性之间的基本权衡。

**离散时间表示：**
对于 $N$ 个时间步的序列 $\{\tau_i\}_{i=1}^N$，我们可以将辐射场表示为：

$$\mathcal{F}_{\text{discrete}} = \{f_{\theta_i}: \mathbb{R}^3 \times \mathbb{S}^2 \rightarrow \mathbb{R}^4\}_{i=1}^N$$

每个 $f_{\theta_i}$ 是时间 $\tau_i$ 的完整神经辐射场。时间演变通过参数集序列 $\{\theta_i\}$ 捕获。

优点：
- 在采样时间点完美重建
- 在离散点没有时间混叠
- 每帧独立优化
- 可以处理任意时间不连续性

缺点：
- 内存复杂度：$O(N \cdot |\theta|)$
- 没有插值就无法生成中间帧
- 帧之间可能存在时间不连续性
- 静态场景元素的冗余学习

**离散场之间的插值：**

对于中间时间 $\tau \in [\tau_i, \tau_{i+1}]$，我们需要插值策略：

1.  **参数插值**：
    $$\theta(\tau) = (1-\alpha)\theta_i + \alpha\theta_{i+1}, \quad \alpha = \frac{\tau - \tau_i}{\tau_{i+1} - \tau_i}$$
    
    由于神经网络参数空间的非凸性，这存在问题。

2.  **输出插值**：
    $$f(\mathbf{x}, \mathbf{d}, \tau) = (1-\alpha)f_{\theta_i}(\mathbf{x}, \mathbf{d}) + \alpha f_{\theta_{i+1}}(\mathbf{x}, \mathbf{d})$$
    
    更稳定但可能产生重影伪影。

3.  **潜在插值**：
    学习共享的编码器-解码器架构，其中只有潜在代码随时间变化。

**连续时间表示：**
一个将时间作为输入的单一网络：

$$f_\theta: \mathbb{R}^3 \times \mathbb{S}^2 \times \mathbb{R} \rightarrow \mathbb{R}^4$$

时间维度通常使用位置编码：
$$\gamma(\tau) = [\sin(2^0\pi\tau), \cos(2^0\pi\tau), ..., \sin(2^{L-1}\pi\tau), \cos(2^{L-1}\pi\tau)]$$

**位置编码的频率分析：**

$L$ 的选择决定了可表示的最大时间频率：
$$f_{\text{max}} = 2^{L-1} \text{ cycles per unit time}$$

对于持续时间为 $T$ 的序列，这转化为 $2^{L-1}T$ 个可分辨的时间特征。

优点：
- 内存复杂度：$O(|\theta|)$（与序列长度无关）
- 自然时间插值
- 平滑运动轨迹
- 紧凑表示

缺点：
- 有限的时间分辨率（受网络容量限制）
- 高频运动可能出现时间混叠
- 更难拟合快速变化或不连续性
- 所有时间瞬间的全局耦合

**混合方法：**

1.  **时间条件超网络：**
    $$f_\theta(\mathbf{x}, \mathbf{d}, \tau) = g_{\phi(\tau)}(\mathbf{x}, \mathbf{d})$$
    
    其中 $\phi: \mathbb{R} \rightarrow \mathbb{R}^{|\theta_g|}$ 为主网络 $g$ 生成参数。

2.  **专家混合：**
    $$f(\mathbf{x}, \mathbf{d}, \tau) = \sum_{k=1}^K w_k(\tau) f_{\theta_k}(\mathbf{x}, \mathbf{d})$$
    
    其中 $w_k(\tau)$ 是时间相关的门控函数，且 $\sum_k w_k(\tau) = 1$。

3.  **分层时间建模：**
    $$f(\mathbf{x}, \mathbf{d}, \tau) = f_{\text{coarse}}(\mathbf{x}, \mathbf{d}, \lfloor\tau\rfloor) + f_{\text{fine}}(\mathbf{x}, \mathbf{d}, \tau)$$
    
    粗略网络处理低频变化，精细网络捕获细节。

**内存-计算权衡：**

设 $M$ 为内存预算，$C$ 为每次查询的计算量，$Q$ 为重建质量：

- 离散：$M = O(N|\theta|)$，$C = O(|\theta|)$，$Q = $ 在 $\{\tau_i\}$ 处完美
- 连续：$M = O(|\theta|)$，$C = O(|\theta|)$，$Q = $ 受网络容量限制
- 混合：$M = O(K|\theta|)$，$C = O(K|\theta|)$，$Q = $ 自适应

最佳选择取决于：
- 序列长度 $N$
- 运动复杂性
- 可用内存
- 实时要求

### 7.1.3 时间基函数

为了平衡表达能力和效率，我们可以使用基函数分解时变场：

$$\sigma(\mathbf{x}, \tau) = \sum_{k=1}^K \sigma_k(\mathbf{x}) \cdot \phi_k(\tau)$$

$$\mathbf{c}(\mathbf{x}, \mathbf{d}, \tau) = \sum_{k=1}^K \mathbf{c}_k(\mathbf{x}, \mathbf{d}) \cdot \phi_k(\tau)$$

这种可分离表示将问题简化为学习 $K$ 个空间场和时间基函数。

**数学基础：**

该分解可以看作是张量分解。设 $\mathcal{T} \in \mathbb{R}^{X \times Y \times Z \times T}$ 为 4D 辐射张量。我们近似：

$$\mathcal{T} \approx \sum_{k=1}^K \mathbf{S}_k \otimes \boldsymbol{\phi}_k$$

其中 $\mathbf{S}_k \in \mathbb{R}^{X \times Y \times Z}$ 是空间分量，$\boldsymbol{\phi}_k \in \mathbb{R}^T$ 是时间分量。

**常见基函数选择：**

1.  **傅里叶基**（用于周期性运动）：
    $$\phi_k(\tau) = \begin{cases}
    \sin(2\pi k\tau/T) & k \text{ 为奇数} \\
    \cos(2\pi k\tau/T) & k \text{ 为偶数}
    \end{cases}$$
    
    属性：
    - 全局支持（影响整个时间线）
    - 正交基：$\langle\phi_i, \phi_j\rangle = \delta_{ij}$
    - 适用于循环运动
    - 通过 FFT 进行频谱解释
    - 帕塞瓦尔定理：$\|\mathbf{f}\|^2 = \sum_k |c_k|^2$
    
    **截断误差分析：**
    对于有界变差 $V(f) < \infty$ 的函数，截断误差：
    $$\|f - f_K\| \leq \frac{V(f)}{\pi K}$$

2.  **B样条基**（用于局部控制）：
    $$\phi_k(\tau) = B_n\left(\frac{\tau - \tau_k}{\Delta\tau}\right)$$
    
    其中 $B_n$ 是 $n$ 阶 B样条核，递归定义为：
    $$B_0(t) = \begin{cases} 1 & |t| < 0.5 \\ 0 & \text{否则} \end{cases}$$
    $$B_n(t) = \int_{-0.5}^{0.5} B_{n-1}(t-s) ds$$
    
    属性：
    - 紧支集：$\phi_k(\tau) = 0$ 对于 $|\tau - \tau_k| > n\Delta\tau/2$
    - $C^{n-1}$ 连续性
    - 单位分解：$\sum_k \phi_k(\tau) = 1$
    - 局部编辑能力
    - 通过 Cox-de Boor 公式高效评估

3.  **小波基**（多分辨率）：
    $$\phi_{j,k}(\tau) = 2^{j/2}\psi(2^j\tau - k)$$
    
    其中 $\psi$ 是母小波（例如，Daubechies，Meyer）。
    
    属性：
    - 多分辨率分析
    - 在时间和频率上都局部化
    - 瞬态特征的高效压缩
    - 快速小波变换：$O(T)$

4.  **学习基**（数据驱动）：
    $$\phi_k(\tau) = \text{softmax}(\text{MLP}_k(\gamma(\tau)))_k$$
    
    属性：
    - 适应数据分布
    - 可以捕获复杂的时间模式
    - 需要正则化以实现平滑性：
      $$\mathcal{L}_{\text{smooth}} = \sum_k \int \left\|\frac{d^2\phi_k}{d\tau^2}\right\|^2 d\tau$$

**最佳基选择：**

选择取决于预期的运动特性：
- 周期性场景 → 傅里叶基
- 局部变化 → B样条或小波
- 多尺度动态 → 小波
- 复杂、非周期性 → 学习基

**基截断分析：**

截断参数 $K$ 控制时间带宽。根据采样定理：
$$K \geq 2 \cdot f_{\text{max}} \cdot T$$

其中 $f_{\text{max}}$ 是最大运动频率，$T$ 是序列持续时间。

**压缩比：**
$$\text{CR} = \frac{\text{原始大小}}{\text{压缩大小}} = \frac{XYZ \cdot T}{K(XYZ + T)}$$

对于 $T \gg K$，我们实现了大约 $T/K$ 的压缩。

**自适应基选择：**

我们可以根据重建误差自适应地选择 $K$：
$$K^* = \arg\min_K \left\{\|f - f_K\|^2 + \lambda K\right\}$$

其中 $\lambda$ 控制稀疏性-准确性权衡。

### 7.1.4 时间变化的频率分析

理解动态场景的频率内容对于避免混叠和确保足够的时间分辨率至关重要。

**奈奎斯特-香农采样定理：**
无混叠重建的基本要求：

$$f_{\text{sample}} \geq 2 f_{\text{max}}$$

其中 $f_{\text{max}}$ 是场景变化的最大频率。

**时间频率谱：**

对于时变辐射场 $f(\mathbf{x}, \tau)$，在每个空间位置的时间傅里叶变换：

$$\hat{f}(\mathbf{x}, \omega) = \int_{-\infty}^{\infty} f(\mathbf{x}, \tau) e^{-i\omega\tau} d\tau$$

功率谱密度：
$$S(\mathbf{x}, \omega) = |\hat{f}(\mathbf{x}, \omega)|^2$$

揭示了每个点的运动频率内容。

**运动特定分析：**

1.  **周期性运动**（例如，旋转物体）：
    对于角频率为 $\omega$ 的运动：
    $$\mathbf{x}(\tau) = \mathbf{R}(\omega\tau)\mathbf{x}_0$$
    
    辐射场变为：
    $$f(\mathbf{x}, \tau) = f_0(\mathbf{R}^T(\omega\tau)\mathbf{x})$$
    
    频率内容：在 $\omega/2\pi$ Hz 处有单个尖峰
    
    每周期最小采样数：
    $$N_{\text{min}} = \frac{4\pi}{\omega \Delta\tau} \geq 4$$

2.  **线性运动**（恒定速度 $\mathbf{v}$）：
    $$\mathbf{x}(\tau) = \mathbf{x}_0 + \mathbf{v}\tau$$
    
    时空梯度：
    $$\frac{\partial f}{\partial \tau} = -\mathbf{v} \cdot \nabla f$$
    
    频率内容取决于空间频率：
    $$\omega_{\text{max}} = \|\mathbf{v}\| \cdot k_{\text{max}}$$
    
    其中 $k_{\text{max}} = 2\pi/\lambda_{\text{min}}$ 是最大空间频率。

3.  **加速运动**：
    $$\mathbf{x}(\tau) = \mathbf{x}_0 + \mathbf{v}_0\tau + \frac{1}{2}\mathbf{a}\tau^2$$
    
    瞬时频率（通过驻相近似）：
    $$\omega(\tau) = (\mathbf{v}_0 + \mathbf{a}\tau) \cdot \mathbf{k}$$
    
    类似啁啾的行为需要时频分析。

4.  **可变形运动**：
    对于一般变形 $\mathbf{W}(\mathbf{x}, \tau)$：
    $$\omega_{\text{local}}(\mathbf{x}) = \left\|\frac{\partial \mathbf{W}}{\partial \tau}\right\| \cdot \|\nabla f\|$$

**混叠分析：**

当 $f_{\text{motion}} > f_{\text{Nyquist}} = f_{\text{sample}}/2$ 时发生混叠。

**混叠频率映射：**
$$f_{\text{alias}} = |f_{\text{true}} - n \cdot f_{\text{sample}}|$$

其中 $n = \text{round}(f_{\text{true}}/f_{\text{sample}})$。

**预混叠滤波器设计：**

为防止混叠，在采样前应用低通滤波器：
$$H(\omega) = \begin{cases}
1 & |\omega| < \omega_c \\
\cos^2\left(\frac{\pi(|\omega| - \omega_c)}{2(\omega_s - \omega_c)}\right) & \omega_c \leq |\omega| \leq \omega_s \\
0 & |\omega| > \omega_s
\end{cases}$$

其中 $\omega_c = 0.8 \cdot \pi f_{\text{sample}}$ 且 $\omega_s = \pi f_{\text{sample}}$。

**实用采样策略：**

给定捕获率 $f_{\text{capture}}$，可表示的运动带宽为：

$$f_{\text{motion}} < \frac{f_{\text{capture}}}{2}$$

对于具有位置编码级别 $L$ 的神经表示：
$$f_{\text{neural}} \leq 2^{L-1}$$

有效时间分辨率为：
$$f_{\text{effective}} = \min(f_{\text{motion}}, f_{\text{neural}})$$

**抗混叠策略：**

1.  **时间超采样**：
    $$f'(\tau) = \frac{1}{K}\sum_{k=0}^{K-1} f\left(\tau + \frac{k}{K}\Delta\tau\right)$$

2.  **运动自适应采样**：
    $$\Delta\tau_{\text{local}} = \frac{C}{\|\partial \mathbf{W}/\partial \tau\|_{\text{max}}}$$
    
    其中 $C$ 是用户定义的常数。

3.  **频域滤波**：
    在重建前在傅里叶域应用带限。

**时间分辨率要求：**
好的，这是您要求的中文翻译，并已将数学公式转换为 LaTeX 格式：

| 运动类型 | 典型频率 | 所需帧率 |
|------------|------------------|--------------|
| 人类行走 | 1-2 Hz | 4-8 |
| 跑步 | 3-4 Hz | 8-16 |
| 手势 | 5-10 Hz | 20-40 |
| 快速运动 | 10-20 Hz | 40-80 |
| 振动 | 20-100 Hz | 80-400 |

## 7.2 变形场建模

我们不直接学习时变辐射场，而是将问题分解为学习一个静态的规范表示和一个时变变形场。这种方法利用了许多动态场景表现出可以通过形变捕获的强时空相关性的洞察。

### 7.2.1 前向变形场

核心思想是建模一个变形场，它将点从规范（参考）空间映射到每个时刻的观测空间：

$$\mathbf{x}_{\text{obs}} = \mathbf{W}(\mathbf{x}_{\text{can}}, \tau)$$

其中：
- $\mathbf{x}_{\text{can}} \in \mathbb{R}^3$: 规范空间中的点
- $\mathbf{x}_{\text{obs}} \in \mathbb{R}^3$: 在时间 $\tau$ 对应的观测空间中的点
- $\mathbf{W}: \mathbb{R}^3 \times \mathbb{R} \rightarrow \mathbb{R}^3$: 变形场

**带变形的体渲染：**

关键挑战是正确地转换体渲染积分。从以下公式开始：

$$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T(t, \tau)\sigma(\mathbf{r}(t), \tau)\mathbf{c}(\mathbf{r}(t), \mathbf{d}, \tau) dt$$

我们代入形变后的坐标，并且必须考虑变量的变化：

$$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T(t)\sigma_{\text{can}}(\mathbf{W}^{-1}(\mathbf{r}(t), \tau))\mathbf{c}_{\text{can}}(\mathbf{W}^{-1}(\mathbf{r}(t), \tau), \mathbf{d}_{\text{can}}) \left|\det J_{\mathbf{W}^{-1}}\right| dt$$

其中：
- $J_{\mathbf{W}^{-1}} = \frac{\partial \mathbf{W}^{-1}}{\partial \mathbf{x}}$ 是逆变形的雅可比矩阵
- $\mathbf{d}_{\text{can}}$ 是变换后的视角方向

**密度变换：**
密度根据以下公式变换：
$$\sigma_{\text{obs}}(\mathbf{x}, \tau) = \sigma_{\text{can}}(\mathbf{W}^{-1}(\mathbf{x}, \tau)) \left|\det J_{\mathbf{W}^{-1}}(\mathbf{x}, \tau)\right|$$

这确保了质量守恒：在任何区域上积分密度都会在规范空间和观测空间中给出相同的总质量。

### 7.2.2 体积保持约束

许多物理变形在局部保持体积（不可压缩性），这为正则化提供了强大的约束。

**不可压缩条件：**
对于体积保持变形：

$$\det(J_{\mathbf{W}}) = 1$$

其中 $J_{\mathbf{W}} = \frac{\partial \mathbf{W}}{\partial \mathbf{x}_{\text{can}}}$ 是变形梯度。

**速度场公式：**
将速度场定义为变形的时间导数：

$$\mathbf{v}(\mathbf{x}, \tau) = \frac{\partial \mathbf{W}(\mathbf{x}, \tau)}{\partial \tau}$$

欧拉形式的不可压缩约束：

$$\nabla \cdot \mathbf{v} = 0$$

**实际实施：**

1. **软约束**（惩罚方法）：
   $$\mathcal{L}_{\text{volume}} = \lambda \int_{\Omega} (\det(J_{\mathbf{W}}) - 1)^2 d\mathbf{x}$$

2. **投影方法**：
   将速度分解为不可压缩分量和势分量：
   $$\mathbf{v} = \mathbf{v}_{\text{incomp}} + \nabla \phi$$
   
   然后求解泊松方程：
   $$\nabla^2 \phi = \nabla \cdot \mathbf{v}_{\text{initial}}$$
   
   并投影：$\mathbf{v}_{\text{incomp}} = \mathbf{v}_{\text{initial}} - \nabla \phi$

3. **流函数**（2D）或**矢量势**（3D）：
   对于 2D：$\mathbf{v} = \nabla^\perp \psi = (\frac{\partial \psi}{\partial y}, -\frac{\partial \psi}{\partial x})$
   
   自动满足 $\nabla \cdot \mathbf{v} = 0$

**近似体积保持：**
对于小变形，将约束线性化：
$$\det(J_{\mathbf{W}}) \approx 1 + \text{tr}(J_{\mathbf{W}} - I) = 1 + \nabla \cdot (\mathbf{W} - \mathbf{x})$$

导致更简单的约束：
$$\nabla \cdot \mathbf{u} = 0$$

其中 $\mathbf{u} = \mathbf{W} - \mathbf{x}$ 是位移场。

### 7.2.3 正则化策略

变形场具有许多自由度，需要仔细正则化以确保物理上合理和稳定的解决方案。

**1. 空间平滑度正则化：**
惩罚高频空间变化：

$$\mathcal{L}_{\text{smooth}} = \int_{\Omega} \int_{\mathcal{T}} \|\nabla^2 \mathbf{W}(\mathbf{x}, \tau)\|_F^2 d\mathbf{x} d\tau$$

替代的一阶平滑度：
$$\mathcal{L}_{\text{grad}} = \int_{\Omega} \int_{\mathcal{T}} \|\nabla \mathbf{W} - I\|_F^2 d\mathbf{x} d\tau$$

**2. 刚性正则化：**
鼓励局部刚性变换（旋转 + 平移）：

$$\mathcal{L}_{\text{rigid}} = \int_{\Omega} \|J_{\mathbf{W}}^T J_{\mathbf{W}} - I\|_F^2 d\mathbf{x}$$

当 $J_{\mathbf{W}}$ 是旋转矩阵时，此项为零。为了计算效率，使用极分解：
$$J_{\mathbf{W}} = \mathbf{R}\mathbf{S}$$

其中 $\mathbf{R}$ 是旋转，$\mathbf{S}$ 是对称矩阵。然后惩罚：
$$\mathcal{L}_{\text{rigid}} = \|\mathbf{S} - I\|_F^2$$

**3. 时间相干性：**
确保随时间平滑运动：

$$\mathcal{L}_{\text{temporal}} = \int_{\Omega} \int_{\mathcal{T}} \left\|\frac{\partial^2 \mathbf{W}}{\partial \tau^2}\right\|^2 d\mathbf{x} d\tau$$

对于离散时间步长：
$$\mathcal{L}_{\text{temporal}} = \sum_{i} \|\mathbf{W}(\mathbf{x}, \tau_{i+1}) - 2\mathbf{W}(\mathbf{x}, \tau_i) + \mathbf{W}(\mathbf{x}, \tau_{i-1})\|^2$$

**4. 弹性势能：**
基于连续介质力学，惩罚变形能量：

$$\mathcal{L}_{\text{elastic}} = \int_{\Omega} \frac{\lambda}{2}(\text{tr}(\mathbf{E}))^2 + \mu \text{tr}(\mathbf{E}^2) d\mathbf{x}$$

其中 $\mathbf{E} = \frac{1}{2}(J_{\mathbf{W}}^T J_{\mathbf{W}} - I)$ 是格林应变张量，$\lambda, \mu$ 是拉梅参数。

**5. 尽可能刚性（ARAP）：**
$$\mathcal{L}_{\text{ARAP}} = \sum_{i} \sum_{j \in \mathcal{N}(i)} w_{ij} \|(\mathbf{W}(\mathbf{x}_i) - \mathbf{W}(\mathbf{x}_j)) - \mathbf{R}_i(\mathbf{x}_i - \mathbf{x}_j)\|^2$$

其中 $\mathbf{R}_i$ 是点 $i$ 处的最佳旋转。

**组合正则化：**
$$\mathcal{L}_{\text{reg}} = \lambda_1 \mathcal{L}_{\text{smooth}} + \lambda_2 \mathcal{L}_{\text{rigid}} + \lambda_3 \mathcal{L}_{\text{temporal}} + \lambda_4 \mathcal{L}_{\text{volume}}$$

权重 $\lambda_i$ 控制相对重要性，并取决于场景类型：
- 流体场景：高 $\lambda_1$，低 $\lambda_2$
- 关节对象：高 $\lambda_2$，中等 $\lambda_3$
- 软体：平衡权重

### 7.2.4 双射映射保证

确保变形场是可逆的（双射）对于稳定的优化和物理上有意义的结果至关重要。

**1. 残差公式：**
将变形参数化为恒等式加上小位移：

$$\mathbf{W}(\mathbf{x}, \tau) = \mathbf{x} + \epsilon \cdot \mathbf{u}(\mathbf{x}, \tau)$$

**可逆性条件：**
根据巴拿赫不动点定理，如果：
$$\|\nabla \mathbf{u}\|_{\text{op}} < \frac{1}{\epsilon}$$

那么 $\mathbf{W}$ 是一个微分同胚（光滑双射）。

**证明草图：**
逆满足：$\mathbf{y} = \mathbf{W}^{-1}(\mathbf{y}) + \epsilon \mathbf{u}(\mathbf{W}^{-1}(\mathbf{y}))$

定义算子 $T(\mathbf{x}) = \mathbf{y} - \epsilon \mathbf{u}(\mathbf{x})$。如果 $\epsilon\|\nabla \mathbf{u}\|_{\text{op}} < 1$，则 $T$ 是一个收缩映射，保证存在唯一不动点。

**2. 神经 ODE 公式：**
通过连续流建模变形：

$$\frac{d\mathbf{x}}{d\tau} = \mathbf{v}(\mathbf{x}, \tau), \quad \mathbf{x}(0) = \mathbf{x}_0$$

**优点：**
- 可逆性由 ODE 理论保证（如果 $\mathbf{v}$ 是 Lipschitz 连续的）
- 逆通过向后积分计算：$\frac{d\mathbf{x}}{d\tau} = -\mathbf{v}(\mathbf{x}, -\tau)$
- 自然的时间连续性

**实现：**
$$\mathbf{W}(\mathbf{x}_0, \tau) = \mathbf{x}_0 + \int_0^\tau \mathbf{v}(\mathbf{x}(s), s) ds$$

使用数值 ODE 求解器（例如，龙格-库塔法）和自适应步长。

**3. 微分同胚配准：**
通过正则化速度场确保光滑双射：

$$\mathcal{L}_{\text{diffeo}} = \int_0^T \int_{\Omega} \|\mathbf{L}\mathbf{v}(\mathbf{x}, t)\|^2 d\mathbf{x} dt$$

其中 $\mathbf{L}$ 是微分算子（例如，$\mathbf{L} = \alpha\nabla^2 + \beta I$）。

**4. 障碍方法：**
添加一个惩罚项，当雅可比行列式趋近于零时，该惩罚项趋近于无穷大：

$$\mathcal{L}_{\text{barrier}} = \int_{\Omega} \phi(\det(J_{\mathbf{W}})) d\mathbf{x}$$

其中 $\phi(s) = -\log(s)$ 对于 $s > 0$ 或 $\phi(s) = \frac{1}{s} - 1$ 对于 $s > 0$。

**实用指南：**
- 初始化接近恒等式：$\mathbf{W}(\mathbf{x}, 0) = \mathbf{x}$
- 训练期间监控 $\min_{\mathbf{x}} \det(J_{\mathbf{W}})$
- 使用渐进优化：从强正则化开始，逐渐放松
- 对于神经网络，使用谱归一化来控制 Lipschitz 常数

## 7.3 光流和场景流约束

### 7.3.1 3D 场景流定义

场景流 $\mathbf{s}: \mathbb{R}^3 \times \mathbb{R} \rightarrow \mathbb{R}^3$ 描述 3D 运动场：

$$\mathbf{s}(\mathbf{x}, \tau) = \frac{\partial \mathbf{W}(\mathbf{x}, \tau)}{\partial \tau}$$

### 7.3.2 从体渲染中获取 2D 光流

图像空间中的 2D 光流可以通过渲染方程从场景流中导出：

$$\mathbf{u}_{\text{2D}}(\mathbf{p}) = \frac{\int_{t_n}^{t_f} T(t)\sigma(\mathbf{r}(t)) w(\mathbf{r}(t)) \Pi(\mathbf{s}(\mathbf{r}(t))) dt}{\int_{t_n}^{t_f} T(t)\sigma(\mathbf{r}(t)) w(\mathbf{r}(t)) dt}$$

其中：
- $\Pi$ 是从 3D 到 2D 的投影算子
- $w(\mathbf{x}) = \frac{\partial C}{\partial \mathbf{x}}$ 是贡献权重

### 7.3.3 流一致性方程

**亮度恒定：**
$$\frac{\partial I}{\partial \tau} + \nabla I \cdot \mathbf{u}_{\text{2D}} = 0$$

**体渲染一致性：**
$$\frac{\partial C(\mathbf{r}, \tau)}{\partial \tau} + \int_{t_n}^{t_f} T(t)\sigma(\mathbf{r}(t), \tau) \nabla_{\mathbf{x}} \mathbf{c} \cdot \mathbf{s} dt = 0$$

### 7.3.4 多视角流约束

对于具有中心 $\{\mathbf{o}_i\}$ 和方向 $\{\mathbf{d}_i\}$ 的多个相机：

**流的对极约束：**
$$(\mathbf{u}_{\text{2D}}^{(i)} - \mathbf{u}_{\text{2D}}^{(j)})^T \mathbf{F}_{ij} \mathbf{p} = 0$$

其中 $\mathbf{F}_{ij}$ 是视图 $i$ 和 $j$ 之间的基本矩阵。

**深度-流一致性：**
$$z_i \mathbf{u}_{\text{2D}}^{(i)} = \mathbf{K}_i \mathbf{R}_i \mathbf{s} + \frac{\partial z_i}{\partial \tau} \mathbf{K}_i \mathbf{R}_i \mathbf{d}_i$$

## 7.4 基于形变的运动表示

### 7.4.1 时空形变函数

通用形变函数：
$$\Psi: \mathbb{R}^3 \times \mathbb{R} \rightarrow \mathbb{R}^3 \times \mathbb{R}$$
$$(\mathbf{x}', \tau') = \Psi(\mathbf{x}, \tau)$$

**可分离形变：**
$$\Psi(\mathbf{x}, \tau) = (\mathbf{W}(\mathbf{x}, \tau), \tau)$$

**不可分离形变：**
$$\Psi(\mathbf{x}, \tau) = (\mathbf{W}(\mathbf{x}, \tau), T(\tau))$$

### 7.4.2 分层运动分解

$$\mathbf{W}(\mathbf{x}, \tau) = \mathbf{W}_{\text{global}}(\tau) \circ \mathbf{W}_{\text{local}}(\mathbf{x}, \tau)$$

其中：
- $\mathbf{W}_{\text{global}}$: 刚体运动（6 自由度）
- $\mathbf{W}_{\text{local}}$: 非刚性变形

**全局运动的李群表示：**
$$\mathbf{W}_{\text{global}}(\tau) = \exp(\sum_{i=1}^{6} \theta_i(\tau) \mathbf{G}_i)$$

其中 $\{\mathbf{G}_i\}$ 是 SE(3) 的生成元。

### 7.4.3 时间插值方案

**线性插值：**
$$\mathbf{W}(\mathbf{x}, \tau) = (1-\alpha)\mathbf{W}(\mathbf{x}, \tau_i) + \alpha\mathbf{W}(\mathbf{x}, \tau_{i+1})$$

**球面线性插值（用于旋转）：**
$$\mathbf{R}(\tau) = \mathbf{R}_i (\mathbf{R}_i^T \mathbf{R}_{i+1})^{\alpha}$$

**三次 Hermite 样条：**
$$\mathbf{W}(\mathbf{x}, \tau) = \sum_{j=0}^{3} h_j(\alpha) \mathbf{c}_j(\mathbf{x})$$

其中 $h_j$ 是 Hermite 基函数。

### 7.4.4 逆形变计算

**不动点迭代：**
$$\mathbf{x}_{n+1} = \mathbf{x}_0 - \mathbf{W}(\mathbf{x}_n, \tau) + \mathbf{x}_{\text{target}}$$

收敛保证当：
$$\|\nabla_{\mathbf{x}} \mathbf{W}\|_{\text{op}} < 1$$

**神经逆网络：**
$$\mathbf{W}^{-1} = f_{\phi}(\mathbf{x}, \tau)$$

具有循环一致性损失：
$$\mathcal{L}_{\text{cycle}} = \|\mathbf{W}(\mathbf{W}^{-1}(\mathbf{x}, \tau), \tau) - \mathbf{x}\|^2$$

## 7.5 规范空间方法

### 7.5.1 规范空间定义

规范空间表示场景的参考配置：

$$\mathcal{C} = \{(\sigma_{\text{can}}(\mathbf{x}), \mathbf{c}_{\text{can}}(\mathbf{x}, \mathbf{d})) : \mathbf{x} \in \Omega_{\text{can}}\}$$

常见选择：
- **静止姿态**：对于关节对象
- **平均形状**：$\mathbf{x}_{\text{can}} = \frac{1}{|\mathcal{T}|}\int_{\mathcal{T}} \mathbf{W}^{-1}(\mathbf{x}, \tau) d\tau$
- **学习的规范空间**：与变形联合优化

### 7.5.2 变形到观测空间

规范空间中的渲染方程：

$$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T_{\text{can}}(t')\sigma_{\text{can}}(\mathbf{r}_{\text{can}}(t'))\mathbf{c}_{\text{can}}(\mathbf{r}_{\text{can}}(t'), \mathbf{d}_{\text{can}}) \left|\frac{dt'}{dt}\right| dt$$

其中：
- $t' = t'(t, \tau)$ 是形变后的射线参数
- $\mathbf{r}_{\text{can}}(t') = \mathbf{W}^{-1}(\mathbf{r}(t), \tau)$
- $\mathbf{d}_{\text{can}} = \frac{d\mathbf{r}_{\text{can}}}{dt'} / \|\frac{d\mathbf{r}_{\text{can}}}{dt'}\|$

### 7.5.3 模板学习策略

**联合优化：**
$$\min_{\theta, \phi} \sum_{i,j} \|C_{ij} - \hat{C}_{ij}(\theta, \phi)\|^2 + \lambda \mathcal{R}(\phi)$$

其中：
- $\theta$: 规范辐射场参数
- $\phi$: 变形场参数
- $\mathcal{R}$: 正则化项

**渐进训练：**
1. 在 $\tau_0$ 处用静态重建初始化
2. 逐渐添加时间样本
3. 联合细化规范空间和变形

### 7.5.4 处理拓扑变化

对于拓扑变化的场景（例如，流体），我们扩展规范空间：

$$\sigma_{\text{can}}(\mathbf{x}) = \sigma_{\text{base}}(\mathbf{x}) + \sum_{k} \alpha_k(\tau) \sigma_{\text{residual}}^{(k)}(\mathbf{x})$$

其中 $\alpha_k(\tau) \in [0,1]$ 控制组件的出现/消失。

### 7.5.5 规范空间中的正则化

**空间平滑度：**
$$\mathcal{L}_{\text{can-smooth}} = \int_{\Omega_{\text{can}}} \|\nabla \sigma_{\text{can}}\|^2 + \|\nabla \mathbf{c}_{\text{can}}\|^2 d\mathbf{x}$$

**稀疏性先验：**
$$\mathcal{L}_{\text{sparse}} = \int_{\Omega_{\text{can}}} \log(1 + \sigma_{\text{can}}^2/\epsilon^2) d\mathbf{x}$$

**变形正则化：**
$$\mathcal{L}_{\text{def-reg}} = \int_{\mathcal{T}} \int_{\Omega} \|\mathbf{W}(\mathbf{x}, \tau) - \mathbf{x}\|^2 \rho(\mathbf{x}) d\mathbf{x} d\tau$$

其中 $\rho(\mathbf{x})$ 权重区域（例如，刚性部分更高）。

## 本章小结

动态神经辐射场通过将时间作为附加维度来扩展静态公式。关键的数学框架包括：

1. **时间扩展体渲染方程：**
   $$C(\mathbf{r}, \tau) = \int_{t_n}^{t_f} T(t, \tau)\sigma(\mathbf{r}(t), \tau)\mathbf{c}(\mathbf{r}(t), \mathbf{d}, \tau) dt$$

2. **变形场公式：**
   $$\mathbf{x}_{\text{obs}} = \mathbf{W}(\mathbf{x}_{\text{can}}, \tau)$$
   
   具有体积保持：$\det(J_{\mathbf{W}}) = 1$

3. **场景流一致性：**
   $$\mathbf{s}(\mathbf{x}, \tau) = \frac{\partial \mathbf{W}(\mathbf{x}, \tau)}{\partial \tau}$$

4. **规范空间渲染：**
   $$C(\mathbf{r}, \tau) = \int T_{\text{can}}\sigma_{\text{can}}(\mathbf{W}^{-1}(\mathbf{r}(t), \tau))\mathbf{c}_{\text{can}} \left|\det J_{\mathbf{W}^{-1}}\right| dt$$

5. **正则化框架：**
   - 空间平滑度：$\mathcal{L}_{\text{smooth}}$
   - 时间相干性：$\mathcal{L}_{\text{temporal}}$
   - 刚性约束：$\mathcal{L}_{\text{rigid}}$

这些公式提供了一种建模动态场景的原则性方法，同时保持了体渲染框架的优雅。

## 练习题

### 基础题

**习题 7.1**: 时间采样分析
考虑一个以频率 $f = 2$ Hz 旋转的物体。如果我们以 30 FPS 采样，是否满足 Nyquist 条件？如果物体同时包含 5 Hz 的振动，需要什么采样率？

*提示*: Nyquist 频率 $f_s \geq 2f_{\max}$

<details>
<summary>答案</summary>

对于 2 Hz 旋转：
- Nyquist 要求: $f_s \geq 2 \times 2 = 4$ Hz
- 30 FPS 采样率满足条件 (30 > 4)

对于 5 Hz 振动：
- Nyquist 要求: $f_s \geq 2 \times 5 = 10$ Hz  
- 需要至少 10 FPS

组合运动：$f_{\max} = 5$ Hz，需要至少 10 FPS 采样率。

</details>

**习题 7.2**: 体积保持约束
给定 2D 变形场 $\mathbf{W}(x,y) = (x + 0.1\sin(y), y + 0.1\cos(x))$，验证这是否满足体积保持约束。
*提示*: 计算 Jacobian 行列式

<details>
<summary>答案</summary>

Jacobian 矩阵：
$$J = \begin{bmatrix}
\frac{\partial W_x}{\partial x} & \frac{\partial W_x}{\partial y} \\
\frac{\partial W_y}{\partial x} & \frac{\partial W_y}{\partial y}
\end{bmatrix} = \begin{bmatrix}
1 & 0.1\cos(y) \\
-0.1\sin(x) & 1
\end{bmatrix}$$

行列式：
$$\det(J) = 1 \times 1 - (0.1\cos(y))(-0.1\sin(x)) = 1 + 0.01\sin(x)\cos(y)$$

由于 $\det(J) \neq 1$，不满足体积保持约束。最大偏差为 0.01。

</details>

**习题 7.3**: 光流计算
给定场景流 $\mathbf{s}(x,y,z) = (0, 0, -v)$ (沿 z 轴的匀速运动) 和相机内参矩阵 $\mathbf{K} = \text{diag}(f, f, 1)$，推导 2D 光流。

*提示*: 使用投影关系 $u = fx/z$, $v = fy/z$

<details>
<summary>答案</summary>

3D 点投影：$(u, v) = (fx/z, fy/z)$

时间导数：
$$\frac{du}{dt} = f \frac{d}{dt}\left(\frac{x}{z}\right) = f \frac{\dot{x}z - x\dot{z}}{z^2}$$

代入 $\dot{x} = 0$, $\dot{z} = -v$：
$$\frac{du}{dt} = f \frac{xv}{z^2} = \frac{uv}{z}$$

类似地：
$$\frac{dv}{dt} = \frac{vv}{z}$$

2D 光流：$\mathbf{u}_{2D} = \frac{v}{z}(u, v)$ - 径向扩张模式。

</details>

### 挑战题

**习题 7.4**: 变形场设计
设计一个变形场 $\mathbf{W}(\mathbf{x}, \tau)$ 来表示心跳运动，满足：
1. 周期为 $T = 1$ 秒
2. 体积在 $\pm 10\%$ 范围内变化
3. 运动主要沿径向

*提示*: 考虑球坐标系和时变缩放

<details>
<summary>答案</summary>

球坐标表示：$\mathbf{x} = (r, \theta, \phi)$

径向缩放变形：
$$\mathbf{W}(r, \theta, \phi, \tau) = s(\tau) \cdot r \cdot \hat{\mathbf{r}}$$

其中缩放因子：
$$s(\tau) = 1 + 0.1\sin(2\pi \tau/T)$$

笛卡尔坐标：
$$\mathbf{W}(\mathbf{x}, \tau) = s(\tau) \cdot \mathbf{x}$$

体积变化：
$$\det(J_{\mathbf{W}}) = s(\tau)^3 \approx 1 + 0.3\sin(2\pi\tau/T)$$

满足 $\pm 10\%$ 体积变化（一阶近似）。

</details>

**习题 7.5**: 正则化分析
证明对于小变形 $\mathbf{W}(\mathbf{x}) = \mathbf{x} + \epsilon \mathbf{u}(\mathbf{x})$，如果 $\|\nabla \mathbf{u}\|_{op} < 1/\epsilon$，则变形是可逆的。

*提示*: 使用 Banach 不动点定理

<details>
<summary>答案</summary>

求逆映射需要解：$\mathbf{y} = \mathbf{x} + \epsilon \mathbf{u}(\mathbf{x})$

定义算子：$T(\mathbf{x}) = \mathbf{y} - \epsilon \mathbf{u}(\mathbf{x})$

对于两点 $\mathbf{x}_1, \mathbf{x}_2$：
$$\|T(\mathbf{x}_1) - T(\mathbf{x}_2)\| = \epsilon\|\mathbf{u}(\mathbf{x}_1) - \mathbf{u}(\mathbf{x}_2)\|$$

由中值定理：
$$\|\mathbf{u}(\mathbf{x}_1) - \mathbf{u}(\mathbf{x}_2)\| \leq \|\nabla \mathbf{u}\|_{op} \|\mathbf{x}_1 - \mathbf{x}_2\|$$

因此：
$$\|T(\mathbf{x}_1) - T(\mathbf{x}_2)\| \leq \epsilon \|\nabla \mathbf{u}\|_{op} \|\mathbf{x}_1 - \mathbf{x}_2\|$$

当 $\epsilon \|\nabla \mathbf{u}\|_{op} < 1$ 时，$T$ 是压缩映射，由 Banach 定理存在唯一不动点，即逆映射存在。

</details>

**习题 7.6**: 多视图一致性
推导两个视图之间的场景流一致性约束。给定两个相机的投影矩阵 $\mathbf{P}_1, \mathbf{P}_2$，以及 3D 场景流 $\mathbf{s}$。

*提示*: 利用极线几何和流的投影关系

<details>
<summary>答案</summary>

3D 点 $\mathbf{X}$ 在两个视图的投影：
- 视图 1: $\mathbf{x}_1 = \mathbf{P}_1 \mathbf{X}$
- 视图 2: $\mathbf{x}_2 = \mathbf{P}_2 \mathbf{X}$

场景流后：$\mathbf{X}' = \mathbf{X} + \mathbf{s}$

投影变化：
- $\mathbf{x}'_1 = \mathbf{P}_1(\mathbf{X} + \mathbf{s})$
- $\mathbf{x}'_2 = \mathbf{P}_2(\mathbf{X} + \mathbf{s})$

2D 流：
- $\mathbf{u}_1 = \mathbf{x}'_1 - \mathbf{x}_1 = \mathbf{P}_1 \mathbf{s}$
- $\mathbf{u}_2 = \mathbf{x}'_2 - \mathbf{x}_2 = \mathbf{P}_2 \mathbf{s}$

消除 $\mathbf{s}$：如果 $\mathbf{P}_1^+$ 是 $\mathbf{P}_1$ 的伪逆，则：
$$\mathbf{u}_2 = \mathbf{P}_2 \mathbf{P}_1^+ \mathbf{u}_1$$

这建立了两视图光流的线性关系。

</details>

**习题 7.7**: 时间插值优化
给定离散时间点 $\{t_i\}$ 的变形场样本 $\{\mathbf{W}_i\}$，设计一个保持 $C^2$ 连续性的插值方案，并分析计算复杂度。

*提示*: 考虑三次样条插值

<details>
<summary>答案</summary>

三次 Hermite 样条插值：

对于 $t \in [t_i, t_{i+1}]$，令 $s = (t - t_i)/(t_{i+1} - t_i)$：

$$\mathbf{W}(t) = h_{00}(s)\mathbf{W}_i + h_{10}(s)\mathbf{m}_i\Delta t + h_{01}(s)\mathbf{W}_{i+1} + h_{11}(s)\mathbf{m}_{i+1}\Delta t$$

其中 Hermite 基函数：
- $h_{00}(s) = 2s^3 - 3s^2 + 1$
- $h_{10}(s) = s^3 - 2s^2 + s$
- $h_{01}(s) = -2s^3 + 3s^2$
- $h_{11}(s) = s^3 - s^2$

切线 $\mathbf{m}_i$ 通过解三对角系统获得（保证 $C^2$ 连续）：

$$\mathbf{m}_{i-1} + 4\mathbf{m}_i + \mathbf{m}_{i+1} = \frac{3}{\Delta t}(\mathbf{W}_{i+1} - \mathbf{W}_{i-1})$$

复杂度：
- 预计算：$O(n)$ 解三对角系统
- 查询：$O(\log n)$ 二分查找 + $O(1)$ 插值
- 空间：$O(n)$ 存储切线

</details>

**习题 7.8**: 拓扑变化建模
设计一个数学框架来处理场景中物体的出现和消失（例如：倒水时水的体积变化）。要求保持时间连续性。

*提示*: 考虑使用软掩码和连续变化的密度场

<details>
<summary>答案</summary>

时变密度场分解：
$$\sigma(\mathbf{x}, t) = \sigma_{\text{static}}(\mathbf{x}) + \sum_k \alpha_k(t) \sigma_k(\mathbf{x})$$

软出现/消失函数：
$$\alpha_k(t) = \begin{cases}
0 & t < t_k^{\text{start}} \\
\frac{1}{2}(1 - \cos(\pi \frac{t - t_k^{\text{start}}}{\Delta t})) & t \in [t_k^{\text{start}}, t_k^{\text{start}} + \Delta t] \\
1 & t \in [t_k^{\text{start}} + \Delta t, t_k^{\text{end}}] \\
\frac{1}{2}(1 + \cos(\pi \frac{t - t_k^{\text{end}}}{\Delta t})) & t \in [t_k^{\text{end}}, t_k^{\text{end}} + \Delta t] \\
0 & t > t_k^{\text{end}} + \Delta t
\end{cases}$$

体积守恒（对于流体）：
$$\int_{\Omega} \sigma(\mathbf{x}, t) d\mathbf{x} = V_0 + \int_0^t Q(\tau) d\tau$$

其中 $Q(t)$ 是流入/流出率。

渲染方程保持不变，自然处理拓扑变化：
$$C(\mathbf{r}, t) = \int T(s,t) \sigma(\mathbf{r}(s), t) \mathbf{c}(\mathbf{r}(s), \mathbf{d}, t) ds$$

</details>

## 常见陷阱与错误 (Gotchas)

### 1. 时间混叠问题

**问题**: 快速运动导致时间欠采样，产生伪影
```
症状：闪烁、运动模糊、鬼影
原因：违反 Nyquist 采样定理
```

**解决方案**:
- 增加时间采样密度
- 使用运动模糊建模
- 时间超分辨率技术

### 2. 非双射变形

**问题**: 变形场产生自相交或折叠
```
症状：渲染伪影、优化不收敛
检测：det(J) < 0 或 det(J) >> 1
```

**解决方案**:
- 限制变形幅度：$\|\mathbf{W} - \mathbf{x}\| < \epsilon$
- 使用残差参数化
- 添加 Jacobian 正则化

### 3. 规范空间歧义

**问题**: 多个规范空间-变形对产生相同观察
```
例子：膨胀的小球 vs 静止的大球
影响：优化陷入局部最优
```

**解决方案**:
- 强先验约束（如最小变形）
- 多阶段优化策略
- 利用多视图约束

### 4. 计算复杂度爆炸

**问题**: 4D 表示导致内存和计算需求剧增
```
静态 NeRF: O(XYZ)
动态 NeRF: O(XYZT)
```

**解决方案**:
- 分解表示（如 K-Planes）
- 自适应采样
- 增量式训练

### 5. 时间插值伪影

**问题**: 离散时间采样之间的插值不自然
```
线性插值：运动不流畅
高阶插值：可能过拟合
```

**解决方案**:
- 物理约束的插值
- 学习式时间编码
- 运动先验正则化

### 6. 光流-几何不一致

**问题**: 2D 光流与 3D 场景流投影不匹配
```
原因：遮挡、透明度变化
表现：优化震荡
```

**解决方案**:
- 遮挡感知的流计算
- 鲁棒损失函数
- 分层优化策略

## 最佳实践检查清单

### 设计阶段
- [ ] **时间建模选择**
  - [ ] 分析场景运动特性（周期性？局部性？）
  - [ ] 估算最大运动频率
  - [ ] 选择合适的时间表示（连续/离散/混合）

- [ ] **变形场设计**
  - [ ] 确定是否需要拓扑变化支持
  - [ ] 选择参数化方式（前向/后向/双向）
  - [ ] 设计正则化策略

- [ ] **内存预算规划**
  - [ ] 估算 4D 表示的内存需求
  - [ ] 考虑分解或压缩方案
  - [ ] 规划增量训练策略

### 实现阶段
- [ ] **数值稳定性**
  - [ ] 检查 Jacobian 计算的数值稳定性
  - [ ] 实现梯度裁剪
  - [ ] 添加 epsilon 防止除零

- [ ] **采样策略**
  - [ ] 实现自适应时间采样
  - [ ] 平衡空间和时间分辨率
  - [ ] 考虑重要性采样

- [ ] **优化流程**
  - [ ] 从粗到细的训练策略
  - [ ] 交替优化规范空间和变形
  - [ ] 监控各项损失的平衡

### 验证阶段
- [ ] **质量指标**
  - [ ] 时间一致性度量
  - [ ] 运动平滑度评估
  - [ ] 多视图一致性检查

- [ ] **性能分析**
  - [ ] 渲染时间 vs 质量权衡
  - [ ] 内存使用监控
  - [ ] 训练收敛速度

- [ ] **鲁棒性测试**
  - [ ] 快速运动场景
  - [ ] 大变形情况
  - [ ] 稀疏视图条件

### 调试技巧
- [ ] **可视化工具**
  - [ ] 变形场可视化（箭头图/网格）
  - [ ] Jacobian 热力图
  - [ ] 时间切片对比

- [ ] **诊断输出**
  - [ ] 逐帧光流检查
  - [ ] 规范空间稳定性
  - [ ] 正则化项贡献分析

- [ ] **消融实验**
  - [ ] 移除时间维度验证静态情况
  - [ ] 固定变形场测试规范空间
  - [ ] 简化场景逐步增加复杂度


# 第16章：相干理论

本章介绍了光学相干性的基本概念，弥合了确定性波动光学与光场统计描述之间的鸿沟。我们建立了描述部分相干性的数学框架，确立了控制相干性如何在光学系统中转换和在空间中传播的关键定理。这些概念对于理解涉及真实光源的干涉和衍射的高级渲染效果至关重要。

## 16.1 时间相干性与光谱

时间相干性描述了光波自身在不同时间延迟下的相关性。对于标量光场 $E(t)$，我们定义时间相干函数：

$$
\Gamma(\tau) = \langle E^*(t)E(t+\tau) \rangle
$$

其中 $\langle \cdot \rangle$ 表示时间平均，而 $*$ 表示复共轭。

对于平稳随机过程，这可以更明确地写为：

$$
\Gamma(\tau) = \lim_{T\to\infty} \frac{1}{T} \int_{-T/2}^{T/2} E^*(t)E(t+\tau) dt
$$

物理解释很简单：$\Gamma(\tau)$ 衡量了光场与其自身时间延迟版本之间的相似程度。对于 $\tau = 0$，我们有 $\Gamma(0) = \langle |E(t)|^2 \rangle$，这是时间平均强度。

### 统计基础

光场 $E(t)$ 被视为一个复值随机过程。对于遍历过程，时间平均等于系综平均：

$$
\Gamma(\tau) = E[E^*(t)E(t+\tau)]
$$

其中 $E[\cdot]$ 表示期望值。自相关函数满足：

$$
\Gamma(-\tau) = \Gamma^*(\tau)
$$

（厄米对称性）

这源于平稳性：$\langle E^*(t)E(t-\tau) \rangle = \langle E^*(t+\tau)E(t) \rangle = \langle E(t)E^*(t+\tau) \rangle^*$。

### 与干涉测量的联系

在路径差为 $\Delta = c\tau$ 的迈克尔逊干涉仪中，输出处的强度为：

$$
I_{out} = I_1 + I_2 + 2\text{Re}[\sqrt{I_1I_2}\gamma(\tau)\exp(i\omega_0\tau)]
$$

其中 $I_1, I_2$ 是来自两个臂的强度，$\gamma(\tau) = \Gamma(\tau)/\Gamma(0)$ 是归一化相干函数。条纹可见度直接衡量 $| \gamma(\tau) |$：

$$
V(\tau) = \frac{2\sqrt{I_1I_2}|\gamma(\tau)|}{I_1 + I_2}
$$

对于平衡臂 ($I_1 = I_2$)，我们有 $V(\tau) = |\gamma(\tau)|$，这提供了一种直接测量时间相干性的实验方法。

### 16.1.1 相干时间与相干长度

归一化时间相干度为：

$$
\gamma(\tau) = \frac{\Gamma(\tau)}{\Gamma(0)}
$$

这个复值函数满足几个重要性质：
- $\gamma(0) = 1$ (零延迟时完美自相关)
- $|\gamma(\tau)| \le 1$ (施瓦茨不等式)
- $\gamma(-\tau) = \gamma^*(\tau)$ (厄米对称性)
- 对于有限带宽光源，当 $\tau \to \infty$ 时，$\gamma(\tau) \to 0$

相干时间 $\tau_c$ 可以通过几种等效方式定义：

1.  **$1/e$ 定义**：$|\gamma(\tau)| = 1/e$ 时的时间
2.  **FWHM 定义**：$|\gamma(\tau)|^2$ 的半高全宽
3.  **积分定义**：$\tau_c = \int_0^\infty |\gamma(\tau)|^2 d\tau$

积分定义是最基本的，因为它表示光场保持相关性的有效持续时间。相干长度则为：

$$
l_c = c \cdot \tau_c
$$

其中 $c$ 是介质中的光速（对于折射率 $n$， $c = c_0/n$）。

### 16.1.2 光谱宽度关系

根据维纳-辛钦定理（详见第16.3节），时间相干函数是功率谱密度 $S(\omega)$ 的傅里叶变换：

$$
\Gamma(\tau) = \int_{-\infty}^\infty S(\omega)e^{i\omega\tau} d\omega
$$

$$
S(\omega) = \frac{1}{2\pi} \int_{-\infty}^\infty \Gamma(\tau)e^{-i\omega\tau} d\tau
$$

这种傅里叶变换关系意味着一个基本的不确定性原理：

$$
\Delta\omega \cdot \Delta\tau \ge K
$$

其中 $K$ 是一个量级为一的常数，取决于宽度的定义方式。

#### 不确定性原理的推导

从傅里叶变换的施瓦茨不等式开始：

$$
\int|f(t)|^2 dt \cdot \int|F(\omega)|^2 d\omega \ge \frac{1}{4\pi}\left|\int f(t) dt\right|^2
$$

将其应用于 $f(t) = t\Gamma(t)$ 并使用分部积分：

$$
\langle\tau^2\rangle\langle\omega^2\rangle \ge \frac{1}{4}
$$

其中 $\langle\tau^2\rangle = \frac{\int\tau^2|\gamma(\tau)|^2 d\tau}{\int|\gamma(\tau)|^2 d\tau}$ 且 $\langle\omega^2\rangle = \frac{\int(\omega-\omega_0)^2S(\omega) d\omega}{\int S(\omega) d\omega}$。

这给出了最小不确定性乘积：

$$
\Delta\tau_{rms} \cdot \Delta\omega_{rms} \ge \frac{1}{2}
$$

对于其他定义：
- FWHM：$\Delta\tau_{FWHM} \cdot \Delta\omega_{FWHM} \ge 4\ln(2)/\pi \approx 0.88$
- $1/e$ 宽度：$\Delta\tau_{1/e} \cdot \Delta\omega_{1/e} \ge 1$

对于特定的光谱形状：

**1. 高斯光谱：**
$$
S(\omega) = S_0 \exp\left[-\frac{(\omega-\omega_0)^2}{2\sigma^2}\right]
$$

其中 $\sigma = 2\pi\Delta\nu/(2\sqrt{2\ln2})$ 对于 FWHM $\Delta\nu$，得到：

$$
\gamma(\tau) = \exp(i\omega_0\tau) \exp(-\sigma^2\tau^2/2)
$$

$$
\tau_c = \int_0^\infty \exp(-\sigma^2\tau^2) d\tau = \frac{\sqrt{\pi}}{2\sigma} \approx 0.44/\Delta\nu
$$

**2. 洛伦兹光谱：**
$$
S(\omega) = \frac{S_0\Gamma_0}{\pi[(\omega-\omega_0)^2 + \Gamma_0^2]}
$$

这给出：
$$
\gamma(\tau) = \exp(i\omega_0\tau) \exp(-\Gamma_0|\tau|)
$$

$$
\tau_c = \frac{1}{2\Gamma_0} = \frac{1}{2\pi\Delta\nu}
$$

**3. 矩形光谱：**
$$
S(\omega) = S_0 \quad \text{for } |\omega-\omega_0| < \pi\Delta\nu, \quad 0 \quad \text{otherwise}
$$

$$
\gamma(\tau) = \exp(i\omega_0\tau) \text{sinc}(\pi\Delta\nu\tau)
$$

$$
\tau_c \approx 1/\Delta\nu
$$

高斯光谱给出最小的时间-带宽积，而矩形光谱在给定带宽下产生最长的相干时间。

### 16.1.3 相干时间示例

不同的光源表现出截然不同的相干特性：

**1. 稳频氦氖激光器（单模）：**
- $\Delta\nu \sim 1 \text{ kHz} - 1 \text{ MHz}$
- $\tau_c \sim 10^{-3} - 1 \text{ s}$
- $l_c \sim 300 \text{ km} - 300,000 \text{ km}$
- 应用：干涉测量、全息术
- 光谱轮廓：由于腔体动力学，接近洛伦兹线形

**2. 半导体激光二极管：**
- $\Delta\nu \sim 1-10 \text{ MHz}$ (单模)
- $\tau_c \sim 10^{-7} - 10^{-8} \text{ s}$
- $l_c \sim 30 - 300 \text{ m}$
- 应用：光纤通信、CD/DVD 读取器
- 温度依赖性：$\Delta\nu \propto T^{3/2}$ (载流子散射)

**3. 发光二极管 (LED)：**
- $\Delta\nu \sim 10 \text{ THz}$ ($\Delta\lambda \sim 30 \text{ nm}$ 对于可见光)
- $\tau_c \sim 10^{-14} \text{ s}$
- $l_c \sim 3 \text{ μm}$
- 应用：显示技术、照明
- 光谱形状：带边发射近似高斯形

**4. 过滤后的阳光（1 nm 滤光片）：**
- $\Delta\nu \sim 1 \text{ THz}$
- $\tau_c \sim 10^{-12} \text{ s}$
- $l_c \sim 0.3 \text{ mm}$
- 自然照明参考
- 黑体光谱：$S(\omega) \propto \omega^3/(\exp(\hbar\omega/k_BT) - 1)$

**5. 白光（未过滤）：**
- $\Delta\nu \sim 300 \text{ THz}$ (400-700 nm 范围)
- $\tau_c \sim 10^{-15} \text{ s}$
- $l_c \sim 0.3 \text{ μm}$
- 对于大多数应用来说基本不相干
- 相干函数：$\gamma(\tau) \approx \text{rect}(\tau/\tau_c)\exp(i\omega_c\tau)$

**6. 钠 D 线（低压灯）：**
- $\Delta\nu \sim 500 \text{ MHz}$ (多普勒展宽)
- $\tau_c \sim 10^{-9} \text{ s}$
- $l_c \sim 0.3 \text{ m}$
- 经典原子光谱线
- 精细结构：D₁ (589.6 nm) 和 D₂ (589.0 nm) 双线

**7. 同步辐射：**
- $\Delta\nu/\nu \sim 10^{-3}$ (典型波荡器)
- $\tau_c \sim 10^{-12} \text{ s}$ 在 $\lambda = 1 \text{ nm}$
- $l_c \sim 0.3 \text{ mm}$
- 高方向性、部分相干光束

**8. 自由电子激光器 (FEL)：**
- $\Delta\nu/\nu \sim 10^{-4} - 10^{-3}$
- $\tau_c \sim 10^{-11} - 10^{-12} \text{ s}$ (X 射线范围)
- $l_c \sim 3-30 \text{ mm}$
- SASE 过程产生部分时间相干性

### 16.1.4 相干时间的测量

相干时间可以通过几种方法测量：

**1. 迈克尔逊干涉测量：**
测量条纹可见度 $V$ 作为路径差 $\Delta$ 的函数：
$$
V(\Delta) = |\gamma(\Delta/c)|
$$

相干长度 $l_c$ 是 $V$ 下降到 $1/e$ 或 $1/2$ 的位置。

对于准单色光：
$$
I(\Delta) = I_0[1 + V(\Delta)\cos(k_0\Delta + \phi)]
$$

可见度包络 $V(\Delta)$ 直接描绘了 $|\gamma(\tau)|$。关键考虑因素：
- 机械稳定性：$\Delta$ 稳定到 $\lambda/20$
- 等臂平衡：补偿色散
- 探测器响应：必须能分辨条纹

**2. 光谱分析：**
使用光谱仪直接测量功率谱 $S(\omega)$，然后计算：
$$
\tau_c = 2\pi/\Delta\omega_{eff}
$$

其中 $\Delta\omega_{eff}$ 是有效光谱宽度。

有效宽度定义：
- RMS 宽度：$\Delta\omega_{rms} = \sqrt{\langle\omega^2\rangle - \langle\omega\rangle^2}$
- FWHM：半高全宽
- 等效宽度：$\Delta\omega_{eq} = \int S(\omega)d\omega/S_{max}$

分辨率要求：
- 光谱仪分辨率 $\delta\omega \ll \Delta\omega$
- 自由光谱范围 > 全光谱宽度

**3. 强度相关：**
对于热光，测量强度相关函数：
$$
g^{(2)}(\tau) = \frac{\langle I(t)I(t+\tau)\rangle}{\langle I(t)\rangle^2}
$$

与场相关性相关：
$$
g^{(2)}(\tau) = 1 + |\gamma(\tau)|^2
$$

这就是汉伯里-布朗-特维斯效应。实现：
- 快速探测器：响应时间 $\ll \tau_c$
- 相关器：数字或模拟
- 光子计数：用于弱信号

**4. 傅里叶变换光谱学：**
在大的 $\Delta$ 范围内扫描迈克尔逊干涉仪：
$$
I(\Delta) = \int S(\omega)[1 + \cos(\omega\Delta/c)]d\omega
$$

傅里叶变换给出 $S(\omega)$：
$$
S(\omega) \propto \mathcal{F}\{I(\Delta) - I_\infty\}
$$

优点：
- 多路复用优势 (Fellgett)
- 高集光能力 (Jacquinot)
- 精确波长校准

## 16.2 空间相干性与杨氏实验

空间相干性描述了光场在同一时间不同空间点之间的相关性。杨氏双缝实验提供了理解空间相干性的经典框架。

空间相干性的基本问题是：给定空间中的两个点 $P_1$ 和 $P_2$，这些点的光场相关性如何？这种相关性决定了当来自这两个点的光结合时是否能观察到干涉条纹。

### 历史背景与现代应用

杨氏1801年的实验不仅证明了光的波动性，而且确立了相干性决定干涉可见度的原理。现代应用包括：

- **恒星干涉测量**：使用基线相关性测量恒星直径
- **光学相干断层扫描**：使用相干门控进行深度成像
- **光刻**：图案转移中的部分相干效应
- **量子光学**：EPR 相关性和贝尔不等式检验

### 16.2.1 互强度与可见度

考虑由光源照射的两个点 $P_1$ 和 $P_2$。互强度为：

$$
J_{12} = \langle E^*(P_1,t)E(P_2,t) \rangle
$$

这个复数量包含相关性的幅度和相位信息。

#### 数学结构

互强度形成一个厄米矩阵：

$$
J = \begin{bmatrix} J_{11} & J_{12} \\ J_{21} & J_{22} \end{bmatrix}
$$

具有以下性质：
- $J_{11}, J_{22} \ge 0$ (强度)
- $J_{21} = J_{12}^*$ (厄米性)
- $\det(J) \ge 0$ (半正定)
- $|J_{12}|^2 \le J_{11}J_{22}$ (施瓦茨不等式)

#### 杨氏双缝分析

在杨氏实验中，缝位于 $r_1$ 和 $r_2$ 处，观察点 $P$ 处的光场为：

$$
E(P,t) = K_1E(r_1,t-t_1) + K_2E(r_2,t-t_2)
$$

其中 $K_i$ 是传播因子，$t_i = |P-r_i|/c$ 是传播时间。

强度变为：

$$
I(P) = |K_1|^2I_1 + |K_2|^2I_2 + 2\text{Re}[K_1^*K_2J_{12} \exp(ik\Delta)]
$$

其中 $\Delta = |P-r_2| - |P-r_1|$ 是路径差。

对于等传播因子的近轴情况：

$$
I(P) = I_1 + I_2 + 2\sqrt{I_1I_2}|J_{12}|\cos(k\Delta + \phi_{12})
$$

其中 $\phi_{12} = \text{arg}(J_{12})$ 是互强度的相位。

#### 条纹可见度

强度在以下范围变化：
- $I_{max} = I_1 + I_2 + 2\sqrt{I_1I_2}|\gamma_{12}|$
- $I_{min} = I_1 + I_2 - 2\sqrt{I_1I_2}|\gamma_{12}|$

可见度定义为：

$$
V = \frac{I_{max} - I_{min}}{I_{max} + I_{min}} = \frac{2\sqrt{I_1I_2}|\gamma_{12}|}{I_1 + I_2}
$$

特殊情况：
1.  **等强度** ($I_1 = I_2 = I_0$)：$V = |\gamma_{12}|$
2.  **完全相干光**：$|\gamma_{12}| = 1$, $V = \frac{2\sqrt{I_1I_2}}{I_1 + I_2} \le 1$
3.  **非相干光**：$\gamma_{12} = 0$, $V = 0$ (无条纹)

### 16.2.2 复相干度

复相干度定义为：

$$
\gamma_{12} = \frac{J_{12}}{\sqrt{J_{11}J_{22}}} = \frac{J_{12}}{\sqrt{I_1I_2}}
$$

#### 数学性质

1.  **施瓦茨不等式**：$|\gamma_{12}| \le 1$
    证明：根据柯西-施瓦茨不等式，$|\langle E_1^*E_2\rangle|^2 \le \langle|E_1|^2\rangle\langle|E_2|^2\rangle$

2.  **厄米对称性**：$\gamma_{12} = \gamma_{21}^*$
    因为 $J_{12} = \langle E_1^*E_2\rangle = \langle E_2^*E_1\rangle^* = J_{21}^*$

3.  **自相干性**：$\gamma_{11} = \gamma_{22} = 1$
    每个点都与自身完美相干

4.  **相位信息**：$\gamma_{12} = |\gamma_{12}|\exp(i\phi_{12})$
    相位 $\phi_{12}$ 影响条纹位置，不影响可见度

#### 物理解释

- **$|\gamma_{12}| = 1$**：完全相干 - 完美相关
- **$0 < |\gamma_{12}| < 1$**：部分相干 - 条纹可见度降低
- **$|\gamma_{12}| = 0$**：非相干 - 无相关性，无条纹

#### 与经典相关性的联系

相干度类似于统计学中的相关系数：

$$
\gamma_{12} = \frac{\text{Cov}(E_1,E_2)}{\sigma_1\sigma_2}
$$

其中 Cov 是协方差，$\sigma_i$ 是标准差。

### 16.2.3 扩展光源

对于扩展光源，光源上的每个点独立贡献（对于非相干光源），总互强度是贡献的总和。

#### 一般公式

对于具有强度分布 $I_s(r_s)$ 的扩展非相干光源，观察点 $r_1, r_2$ 处的互强度为：

$$
J_{12} = \iiint I_s(r_s) K^*(r_1,r_s)K(r_2,r_s) d^3r_s
$$

其中 $K(r,r_s)$ 是从光源点 $r_s$ 到观察点 $r$ 的传播核。

在菲涅尔近似中：
$$
K(r,r_s) = \left(\frac{\exp(ikR)}{R}\right) \times \exp\left[ik\frac{(r-r_s)^2}{2R}\right]
$$

其中 $R = |r_z - r_{s,z}|$ 是轴向距离。

#### 圆形光源分析

考虑一个直径为 $D$ 的均匀照明圆形非相干光源，距离观察平面 $z$。对于相距 $d$ 的两个观察点 $P_1$ 和 $P_2$：

$$
J_{12} = \iint_{source} I_s(\xi,\eta) \exp\left[ik\frac{(r_2-r_1)\cdot r_s}{z}\right] d\xi d\eta
$$

在近轴近似下，令 $r_2-r_1 = (d,0,0)$：

$$
\gamma_{12} = \frac{2}{\pi a^2} \iint_{|\rho|<a} \exp\left(ik\frac{d\xi}{z}\right) d\xi d\eta \quad \text{where } a = D/2
$$

转换为极坐标并积分：

$$
|\gamma_{12}| = \left|\frac{2J_1(\pi Dd/\lambda z)}{\pi Dd/\lambda z}\right|
$$

其中 $J_1$ 是第一类一阶贝塞尔函数。

#### 详细推导

从光源平面中的极坐标 $(\rho,\theta)$ 开始：

$$
\gamma_{12} = \frac{1}{\pi a^2} \int_0^{2\pi} \int_0^a \exp\left(ik\frac{d\rho\cos(\theta)}{z}\right) \rho d\rho d\theta
$$

角度积分给出：
$$
\int_0^{2\pi} \exp\left(ik\frac{d\rho\cos(\theta)}{z}\right) d\theta = 2\pi J_0\left(\frac{kd\rho}{z}\right)
$$

其中 $J_0$ 是零阶贝塞尔函数。径向积分变为：

$$
\gamma_{12} = \frac{2}{a^2} \int_0^a J_0\left(\frac{kd\rho}{z}\right) \rho d\rho
$$

使用恒等式 $\int_0^x tJ_0(t)dt = xJ_1(x)$：

$$
\gamma_{12} = \frac{2J_1(kda/z)}{kda/z} = \frac{2J_1(\pi Dd/\lambda z)}{\pi Dd/\lambda z}
$$

这就是著名的相干艾里斑。

好的，以下是逐字翻译为中文并转换数学公式为LaTeX的版本：

#### 相干半径

横向相干半径 $\rho_c$ 定义为 $|γ_{12}|$ 首次达到零的位置：

**$\pi D d_{coh}/\lambda z = 3.832$** (J₁ 的第一个零点)

**$d_{coh} = \rho_c \approx 1.22\lambda z/D = 0.61\lambda z/a$**

这与衍射中的艾里斑半径相同！

#### 其他光源几何形状

**1. 矩形光源 ($L_x \times L_y$):**
**$\gamma_{12} = \text{sinc}(\pi L_x \Delta x/\lambda z) \times \text{sinc}(\pi L_y \Delta y/\lambda z)$**

相干长度：$\rho_x = \lambda z/L_x$, $\rho_y = \lambda z/L_y$

**2. 双星 (两个点光源，间隔 $\theta$):**
**$\gamma_{12} = [\exp(i\pi\theta d/\lambda) + \exp(-i\pi\theta d/\lambda)]/2 = \cos(\pi\theta d/\lambda)$**

第一个零点在 $d = \lambda/(2\theta)$ - 迈克尔逊恒星干涉仪的基础

**3. 高斯光源 (1/e 半径 $w_s$):**
**$\gamma_{12} = \exp(-\pi^2 w_s^2 d^2/\lambda^2 z^2)$**

相干半径 (1/e): $\rho_c = \lambda z/(\pi w_s)$

#### 角直径与相干性

对于张角为 $\theta_s = D/z$ 的光源：

**$\rho_c \approx \lambda/\theta_s$**

这个基本关系表明：
- 较小的角光源 → 较大的相干区域
- 恒星 (微角秒) 产生米级相干性
- 太阳 ($\theta_s \approx 0.5^\circ$) 在 500 nm 处产生 $\rho_c \approx 0.07$ mm 的相干性

## 16.3 互相关函数与维纳-辛钦定理

互相关函数提供了部分相干场完整的二阶统计描述。本节将建立连接时域相干性与频域光谱特性的数学框架。

### 16.3.1 一般定义

对于标量场 $U(\mathbf{r},t)$，互相关函数为：

**$\Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau) = \langle U^*(\mathbf{r}_1,t)U(\mathbf{r}_2,t+\tau) \rangle$**

这个四点函数 (两个空间，一个时间，一个系综) 捕获了所有二阶相干特性。

#### 互相关函数的性质

1.  **厄米对称性**: $\Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau) = \Gamma^*(\mathbf{r}_2,\mathbf{r}_1,-\tau)$

2.  **正半定性**: 对于任意函数 $f_1(\mathbf{r})$, $f_2(\mathbf{r})$:
    $\iint \iint f_1^*(\mathbf{r}_1)\Gamma(\mathbf{r}_1,\mathbf{r}_2,0)f_2(\mathbf{r}_2) d\mathbf{r}_1 d\mathbf{r}_2 \ge 0$

3.  **边界值**:
    - $\Gamma(\mathbf{r},\mathbf{r},0) = I(\mathbf{r})$ (强度)
    - $|\Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau)| \le \sqrt{[\Gamma(\mathbf{r}_1,\mathbf{r}_1,0)\Gamma(\mathbf{r}_2,\mathbf{r}_2,0)]}$

4.  **平稳性**: 对于平稳场：
    $\Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau)$ 仅依赖于 $\tau$，而不依赖于绝对时间

等时互相关函数 (互强度) 为：

**$J(\mathbf{r}_1,\mathbf{r}_2) = \Gamma(\mathbf{r}_1,\mathbf{r}_2,0) = \langle U^*(\mathbf{r}_1,t)U(\mathbf{r}_2,t) \rangle$**

归一化互相关函数：

**$\gamma(\mathbf{r}_1,\mathbf{r}_2,\tau) = \Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau)/\sqrt{[\Gamma(\mathbf{r}_1,\mathbf{r}_1,0)\Gamma(\mathbf{r}_2,\mathbf{r}_2,0)]}$**

### 16.3.2 交叉谱密度

交叉谱密度是互相关函数关于时间延迟的傅里叶变换：

**$W(\mathbf{r}_1,\mathbf{r}_2,\omega) = \int_{-\infty}^\infty \Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau)e^{-i\omega\tau} d\tau$**

**$\Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau) = (1/2\pi) \int_{-\infty}^\infty W(\mathbf{r}_1,\mathbf{r}_2,\omega)e^{i\omega\tau} d\omega$**

#### 物理意义

$W(\mathbf{r}_1,\mathbf{r}_2,\omega)$ 表示在频率 $\omega$ 处，点 $\mathbf{r}_1$ 和 $\mathbf{r}_2$ 之间谱分量的相关性。它可以写成：

**$W(\mathbf{r}_1,\mathbf{r}_2,\omega) = \langle \hat{U}^*(\mathbf{r}_1,\omega)\hat{U}(\mathbf{r}_2,\omega) \rangle$**

其中 $\hat{U}(\mathbf{r},\omega)$ 是 $U(\mathbf{r},t)$ 的傅里叶变换。

#### 交叉谱密度的性质

1.  **厄米性**: $W(\mathbf{r}_1,\mathbf{r}_2,\omega) = W^*(\mathbf{r}_2,\mathbf{r}_1,\omega)$

2.  **对角线非负**: $W(\mathbf{r},\mathbf{r},\omega) \ge 0$ (功率谱)

3.  **谱相干度**:
    **$\mu(\mathbf{r}_1,\mathbf{r}_2,\omega) = W(\mathbf{r}_1,\mathbf{r}_2,\omega)/\sqrt{[W(\mathbf{r}_1,\mathbf{r}_1,\omega)W(\mathbf{r}_2,\mathbf{r}_2,\omega)]}$**

    且 $|\mu(\mathbf{r}_1,\mathbf{r}_2,\omega)| \le 1$

4.  **总强度**:
    **$I(\mathbf{r}) = (1/2\pi) \int_{-\infty}^\infty W(\mathbf{r},\mathbf{r},\omega) d\omega$**

### 16.3.3 维纳-辛钦定理

维纳-辛钦定理建立了自相关与功率谱之间的基本联系。

#### 定理陈述

对于位置 $\mathbf{r}$ 处的平稳场：

**$\Gamma(\mathbf{r},\mathbf{r},\tau) = \langle U^*(\mathbf{r},t)U(\mathbf{r},t+\tau) \rangle = \int_{-\infty}^\infty S(\mathbf{r},\omega)e^{i\omega\tau} d\omega$**

**$S(\mathbf{r},\omega) = (1/2\pi) \int_{-\infty}^\infty \Gamma(\mathbf{r},\mathbf{r},\tau)e^{-i\omega\tau} d\tau$**

其中 $S(\mathbf{r},\omega) = W(\mathbf{r},\mathbf{r},\omega)/(2\pi)$ 是功率谱密度。

#### 含义

1.  **能量守恒**:
    **$I(\mathbf{r}) = \langle|U(\mathbf{r},t)|^2\rangle = \Gamma(\mathbf{r},\mathbf{r},0) = \int_{-\infty}^\infty S(\mathbf{r},\omega) d\omega$**

2.  **相干-带宽积**:
    该定理意味着 $\Delta\omega \cdot \tau_c \sim 2\pi$，一种不确定性原理的形式

3.  **白光极限**:
    对于 $S(\omega) = S_0$ (常数)，$\Gamma(\tau) = 2\pi S_0 \delta(\tau)$ - 完美非相干

#### 广义维纳-辛钦定理

对于非平稳场，我们使用维格纳分布：

**$W_U(t,\omega) = \int U^*(t-\tau/2)U(t+\tau/2)e^{-i\omega\tau} d\tau$**

系综平均给出：
**$\langle W_U(t,\omega) \rangle = \int \Gamma(t-\tau/2,t+\tau/2,\tau)e^{-i\omega\tau} d\tau$**

### 16.3.4 准单色近似

许多实际光源的光谱宽度相对于中心频率较窄：$\Delta\omega/\omega_0 \ll 1$。

#### 数学公式

对于中心频率为 $\omega_0$ 的窄带光：

**$U(\mathbf{r},t) = A(\mathbf{r},t)e^{-i\omega_0 t}$**

其中 $A(\mathbf{r},t)$ 是一个缓慢变化的复振幅。

互相关函数变为：

**$\Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau) = \langle A^*(\mathbf{r}_1,t)A(\mathbf{r}_2,t+\tau) \rangle e^{i\omega_0\tau}$**

**$= J(\mathbf{r}_1,\mathbf{r}_2)e^{i\omega_0\tau}g(\tau)$**

其中：
- $J(\mathbf{r}_1,\mathbf{r}_2) = \langle A^*(\mathbf{r}_1,t)A(\mathbf{r}_2,t) \rangle$ 是互强度
- $g(\tau)$ 是一个缓慢变化的包络，且 $g(0) = 1$

#### 有效条件

1.  **光谱条件**: $\Delta\omega/\omega_0 \ll 1$
2.  **时间条件**: $A(t)$ 的变化时间尺度 $\gg 1/\omega_0$
3.  **传播条件**: $\Delta k \cdot L \ll 2\pi$ ($L$ 是传播距离)

#### 交叉谱密度

在准单色近似下：

**$W(\mathbf{r}_1,\mathbf{r}_2,\omega) \approx J(\mathbf{r}_1,\mathbf{r}_2)G(\omega-\omega_0)$**

其中 $G(\omega)$ 是 $g(\tau)$ 的傅里叶变换，中心在 $\omega = 0$。

#### 应用

1.  **干涉**: 相位差 $k\Delta \approx k_0\Delta$ 使用中心波数
2.  **衍射**: 单频计算即可
3.  **相干传播**: 简化为单色情况
4.  **光谱学**: 原子跃迁的自然线宽

## 16.4 范西特-泽尼克定理

范西特-泽尼克定理是相干理论中最重要的结果之一，它阐明了空间相干性如何从非相干光源通过传播而产生。

### 16.4.1 定理陈述

对于源平面中强度分布为 $I_s(\xi,\eta)$ 的平面非相干光源，在距离 $z$ 处的观测平面中的复相干度为：

**$\gamma_{12} = \frac{\iint I_s(\xi,\eta) \exp[ik(\mathbf{r}_1-\mathbf{r}_2)\cdot(\xi,\eta)/z] d\xi d\eta}{\iint I_s(\xi,\eta) d\xi d\eta}$**

在傍轴近似下，这变为：

**$\gamma_{12} = \frac{\mathcal{F}\{I_s(\xi,\eta)\}|_{(u,v)=(x_1-x_2)/\lambda z, (y_1-y_2)/\lambda z}}{\iint I_s(\xi,\eta) d\xi d\eta}$**

其中 $\mathcal{F}$ 表示二维傅里叶变换。

#### 严格推导

从非相干光源的互强度开始：

**$J(\mathbf{r}_1,\mathbf{r}_2) = \iint I_s(\mathbf{r}_s) G^*(\mathbf{r}_1,\mathbf{r}_s)G(\mathbf{r}_2,\mathbf{r}_s) d^2\mathbf{r}_s$**

其中 $G$ 是格林函数。在自由空间中：

**$G(\mathbf{r},\mathbf{r}_s) = \exp(ik|\mathbf{r}-\mathbf{r}_s|)/4\pi|\mathbf{r}-\mathbf{r}_s|$**

对于 $z \gg |\mathbf{r}_\perp|$, $|\mathbf{r}_s|$ 的傍轴传播：

**$|\mathbf{r}-\mathbf{r}_s| \approx z + (\mathbf{r}_\perp - \mathbf{r}_s)^2/2z$**

代入并简化：

**$J(\mathbf{r}_1,\mathbf{r}_2) = (k/2\pi z)^2 \exp[ik(|\mathbf{r}_1|^2 - |\mathbf{r}_2|^2)/2z] \times \iint I_s(\mathbf{r}_s) \exp[ik(\mathbf{r}_2-\mathbf{r}_1)\cdot\mathbf{r}_s/z] d^2\mathbf{r}_s$**

归一化相干度变为：

**$\gamma_{12} = \exp[ik(|\mathbf{r}_1|^2 - |\mathbf{r}_2|^2)/2z] \times \mathcal{F}\{I_s\}(k(\mathbf{r}_2-\mathbf{r}_1)/2\pi z) / I_{total}$**

对于与轴距离相等的点，二次相位项抵消，得到经典结果。

### 16.4.2 物理诠释

该定理揭示：
1.  相干度是光源强度分布的归一化傅里叶变换
2.  较大的光源产生较小的相干区域
3.  相干函数继承了光源的对称性

### 16.4.3 示例与应用

**圆形光源:**
对于半径为 $a$ 的均匀圆形光源：

**$\gamma_{12} = \frac{2J_1(2\pi a|\mathbf{r}_1-\mathbf{r}_2|/\lambda z)}{(2\pi a|\mathbf{r}_1-\mathbf{r}_2|/\lambda z)}$**

相干半径 (第一个零点) 为：
**$\rho_c = 0.61\lambda z/a$**

**矩形光源:**
对于尺寸为 $L_x \times L_y$ 的矩形光源：

**$\gamma_{12} = \text{sinc}(\pi L_x(x_1-x_2)/\lambda z) \times \text{sinc}(\pi L_y(y_1-y_2)/\lambda z)$**

**双星:**
对于间隔为 $\theta$ 的两个点光源：

**$\gamma_{12} = \cos(\pi\theta|\mathbf{r}_1-\mathbf{r}_2|/\lambda)$**

这是恒星干涉测量的基础。

### 16.4.4 广义范西特-泽尼克定理

对于非平面几何和任意传播距离，广义形式使用格林函数：

**$J(\mathbf{r}_1,\mathbf{r}_2) = \iint I_s(\mathbf{r}_s)G^*(\mathbf{r}_1,\mathbf{r}_s)G(\mathbf{r}_2,\mathbf{r}_s) d^2\mathbf{r}_s$**

其中 $G$ 是适用于该几何的格林函数。

## 16.5 部分相干光的传播

理解相干特性在传播过程中如何变化对于光学系统的精确建模至关重要。

### 16.5.1 互相关函数的传播定律

互相关函数满足一对波动方程：

**$\nabla_1^2 J(\mathbf{r}_1,\mathbf{r}_2) + k^2 J(\mathbf{r}_1,\mathbf{r}_2) = 0$**
**$\nabla_2^2 J(\mathbf{r}_1,\mathbf{r}_2) + k^2 J(\mathbf{r}_1,\mathbf{r}_2) = 0$**

其中 $\nabla_i^2$ 作用于坐标 $\mathbf{r}_i$。这些被称为沃尔夫方程。

### 16.5.2 自由空间中的传播

对于从平面 $z=0$ 到平面 $z$ 的传播，使用菲涅尔近似：

**$J(x_1,y_1,x_2,y_2;z) = (k/2\pi z)^2 \iiint\int J_0(\xi_1,\eta_1,\xi_2,\eta_2) \times \exp[ik/2z((x_1-\xi_1)^2 + (y_1-\eta_1)^2 - (x_2-\xi_2)^2 - (y_2-\eta_2)^2)] d\xi_1 d\eta_1 d\xi_2 d\eta_2$**

### 16.5.3 相干模式表示

任何部分相干场都可以分解为相干模式：

**$J(\mathbf{r}_1,\mathbf{r}_2) = \sum_n \lambda_n \phi_n^*(\mathbf{r}_1)\phi_n(\mathbf{r}_2)$**

其中 $\lambda_n$ 是特征值，$\phi_n$ 是满足以下条件的正交特征函数：

**$\int J(\mathbf{r}_1,\mathbf{r}_2)\phi_n(\mathbf{r}_2) d^2\mathbf{r}_2 = \lambda_n\phi_n(\mathbf{r}_1)$**

这被称为相干模式表示，类似于主成分分析。

### 16.5.4 薛尔模型光源

一类重要的光源满足：

**$J(\mathbf{r}_1,\mathbf{r}_2) = \sqrt{[I(\mathbf{r}_1)I(\mathbf{r}_2)]} \mu(\mathbf{r}_1-\mathbf{r}_2)$**

其中 $\mu$ 是仅依赖于 $\mathbf{r}_1-\mathbf{r}_2$ 的相干函数。这些薛尔模型光源在传播过程中保持这种形式：

**$J_z(\mathbf{r}_1,\mathbf{r}_2) = \sqrt{[I_z(\mathbf{r}_1)I_z(\mathbf{r}_2)]} \mu_z(\mathbf{r}_1-\mathbf{r}_2)$**

### 16.5.5 与体渲染的联系

在部分相干照明的体渲染背景下，每个点的辐射度变为：

**$L(\mathbf{x},\omega) = \iint W(\mathbf{x},\mathbf{x}',\omega)\sigma_s(\mathbf{x}')p(\mathbf{x}',\omega'\to\omega) d\omega' d^3\mathbf{x}'$**

其中 $W(\mathbf{x},\mathbf{x}',\omega)$ 是照明的交叉谱密度。这概括了标准体渲染方程，以包含相干效应，这对于以下方面的精确建模至关重要：
- 激光扫描显微镜
- 光学相干断层扫描
- 全息显示
- 干涉成像

## 章总结

本章建立了光学相干性的数学框架，将光场的统计特性与可观测的干涉现象联系起来。关键概念包括：

1.  **时间相干性**: 通过傅里叶变换关系与光谱带宽相关。相干时间 $\tau_c \approx 1/\Delta\nu$ 决定了可以发生干涉的最大光程差。

2.  **空间相干性**: 由复相干度 $\gamma_{12}$ 量化，可通过干涉实验中的条纹可见度直接观测。

3.  **互相关函数**: $\Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau)$ 提供了部分相干场的完整二阶统计描述，交叉谱密度 $W$ 是其频域对应物。

4.  **范西特-泽尼克定理**: 确立了远场空间相干性是光源强度分布的傅里叶变换，对于理解扩展光源的相干性至关重要。

5.  **传播定律**: 沃尔夫方程控制着相干特性在传播过程中的演变，在包含统计效应的同时保持了波动性质。

体渲染方程概括为：
**$L(\mathbf{x},\omega) = \iint W(\mathbf{x},\mathbf{x}',\omega)\sigma_s(\mathbf{x}')p(\mathbf{x}',\omega'\to\omega) d\omega' d^3\mathbf{x}'$**

通过交叉谱密度纳入了相干效应。

## 练习

### 练习 16.1: 相干时间计算
氦氖激光器的谱线宽度为 1 GHz。计算：
a) 相干时间 $\tau_c$
b) 相干长度 $l_c$
c) 干涉条纹可见度 > 0.5 的最大光程差

*提示: 对于高斯光谱，使用关系 $\tau_c \approx 0.44/\Delta\nu$。*

<details>
<summary>解答</summary>

a) $\tau_c = 0.44/\Delta\nu = 0.44/(10^9 \text{ Hz}) = 4.4 \times 10^{-10} \text{ s} = 0.44 \text{ ns}$

b) $l_c = c \cdot \tau_c = (3 \times 10^8 \text{ m/s})(4.4 \times 10^{-10} \text{ s}) = 0.132 \text{ m} = 13.2 \text{ cm}$

c) 对于高斯光谱，可见度 $V = \exp(-\pi^2\tau^2\Delta\nu^2/2\ln2)$
   设 $V = 0.5$: $\tau = \sqrt{2\ln2 \cdot \ln2}/\pi\Delta\nu \approx 0.37/\Delta\nu$
   最大光程差 $= c \cdot \tau = 0.37c/\Delta\nu = 11.1 \text{ cm}$
</details>

### 练习 16.2: 扩展光源的杨氏双缝实验
钠灯 ($\lambda = 589 \text{ nm}$) 具有直径为 2 mm 的圆形孔径，在 1 m 距离处照射间距为 0.5 mm 的杨氏双缝。计算：
a) 缝处的相干度
b) 条纹可见度
c) 可见度降至零的狭缝间距

*提示: 对圆形光源应用范西特-泽尼克定理。*

<details>
<summary>解答</summary>

a) 使用圆形光源的范西特-泽尼克定理：
   $\gamma_{12} = \frac{2J_1(\pi Dd/\lambda z)}{(\pi Dd/\lambda z)}$
   其中 $D = 2 \text{ mm}$, $d = 0.5 \text{ mm}$, $z = 1 \text{ m}$, $\lambda = 589 \text{ nm}$

   $\pi Dd/\lambda z = \pi(2\times10^{-3})(0.5\times10^{-3})/(589\times10^{-9})(1) = 5.33$
   $\gamma_{12} = 2J_1(5.33)/5.33 \approx 2(-0.327)/5.33 = -0.123$

b) 对于等强度狭缝：$V = |\gamma_{12}| = 0.123$

c) 当 $\pi Dd/\lambda z = 3.83$ (J₁ 的第一个零点) 时，第一个零点出现
   $d = 3.83\lambda z/\pi D = 3.83(589\times10^{-9})(1)/\pi(2\times10^{-3}) = 0.36 \text{ mm}$
</details>

### 练习 16.3: 维纳-辛钦应用
光源具有洛伦兹谱：$S(\omega) = S_0\Gamma^2/[(\omega-\omega_0)^2 + \Gamma^2]$
推导：
a) 时间相干函数 $\Gamma(\tau)$
b) 相干时间
c) 与相同半高全宽的高斯光谱进行比较

*提示: 使用傅里叶变换的围道积分。*

<details>
<summary>解答</summary>

a) $\Gamma(\tau) = \int S(\omega)e^{i\omega\tau} d\omega = S_0\Gamma^2 \int \frac{e^{i\omega\tau}}{[(\omega-\omega_0)^2 + \Gamma^2]} d\omega$

   对于 $\tau > 0$，使用留数定理，极点在 $\omega = \omega_0 + i\Gamma$:
   $\Gamma(\tau) = 2\pi i S_0\Gamma^2 \cdot \frac{e^{i(\omega_0+i\Gamma)\tau}}{2i\Gamma} = \pi S_0\Gamma e^{i\omega_0\tau}e^{-\Gamma\tau}$

   归一化：$\gamma(\tau) = e^{i\omega_0\tau}e^{-\Gamma|\tau|}$

b) $\tau_c = \int_0^\infty |\gamma(\tau)|^2 d\tau = \int_0^\infty e^{-2\Gamma\tau} d\tau = 1/(2\Gamma)$

   对于 FWHM $= 2\Gamma$: $\tau_c = 1/\text{FWHM}$

c) 高斯：$\tau_c = 0.44/\Delta\nu$
   洛伦兹：$\tau_c = 1/(2\pi\Delta\nu) \approx 0.16/\Delta\nu$
   洛伦兹由于其延伸的翼部，具有更短的相干时间。
</details>
好的，我将逐字翻译您提供的文本为中文，并将所有数学公式转换为 LaTeX 格式。

</details>

### 练习 16.4：相干模式分解
对于具有高斯强度 $I(x) = I_0\exp(-x^2/w_0^2)$ 和高斯相干性 $\mu(\Delta x) = \exp(-\Delta x^2/2\sigma_c^2)$ 的 Schell 模型光源，找出前三个相干模式。

*提示：使用厄米-高斯函数作为基底。*

<details>
<summary>解答</summary>

特征值方程：$\int J(x,x')\varphi_n(x') dx' = \lambda_n\varphi_n(x)$

对于此高斯 Schell 模型，特征函数是厄米-高斯函数：
$\varphi_n(x) = (2^n n!\sqrt{\pi} \sigma)^{-1/2} H_n(x/\sigma) \exp(-x^2/2\sigma^2)$

其中 $\sigma^4 = w_0^2\sigma_c^2/2$

特征值：$\lambda_n = \lambda_0(\sigma_c^2/(\sigma_c^2 + w_0^2))^n$
其中 $\lambda_0 = I_0\sqrt{2\pi\sigma_c^2w_0^2/(\sigma_c^2 + w_0^2)}$

前三个模式：
- $n=0$: $\varphi_0(x) = (\pi\sigma^2)^{-1/4} \exp(-x^2/2\sigma^2)$, $\lambda_0$
- $n=1$: $\varphi_1(x) = (\pi\sigma^2)^{-1/4} \sqrt{2}x/\sigma \exp(-x^2/2\sigma^2)$, $\lambda_1$
- $n=2$: $\varphi_2(x) = (\pi\sigma^2)^{-1/4} (2x^2/\sigma^2 - 1)/\sqrt{2} \exp(-x^2/2\sigma^2)$, $\lambda_2$
</details>

### 练习 16.5：相干性传播（挑战）
一个部分相干光束的初始互强度为 $J_0(x_1,x_2) = \exp(-|x_1-x_2|^2/2\sigma_0^2)\exp(-(x_1^2+x_2^2)/4w_0^2)$。求传播距离 $z$ 后的 $J(x_1,x_2,z)$。

*提示：使用互强度的菲涅尔传播积分。*

<details>
<summary>解答</summary>

使用菲涅尔传播：
$J(x_1,x_2,z) = (k/2\pi z)^2 \iint J_0(\xi_1,\xi_2) \exp[ik(x_1-\xi_1)^2/2z] \exp[-ik(x_2-\xi_2)^2/2z] d\xi_1d\xi_2$

代入高斯形式并配方：
$J(x_1,x_2,z) = A(z) \exp(-|x_1-x_2|^2/2\sigma^2(z)) \exp(-(x_1^2+x_2^2)/4w^2(z))$

其中：
- $w^2(z) = w_0^2(1 + z^2/z_R^2)$, $z_R = \pi w_0^2/\lambda$
- $\sigma^2(z) = \sigma_0^2 + \lambda^2z^2/(4\pi^2\sigma_0^2)$
- $A(z)$ 包含归一化因子

光束保持高斯-Schell 形式，参数随传播演变。
</details>

### 练习 16.6：恒星干涉测量中的 Van Cittert-Zernike 定理（挑战）
两台望远镜基线分离 $B$，观测一颗角分离为 $\theta$、强度比为 $R$ 的双星。推导可见度曲线 $V(B)$ 并展示如何提取 $\theta$ 和 $R$。

*提示：建模为两个非相干点源。*

<details>
<summary>解答</summary>

对于位于 $\pm\theta/2$ 角度、强度为 $I_1, I_2$ 的两颗恒星：
$\gamma_{12} = [I_1\exp(ik\theta B/2) + I_2\exp(-ik\theta B/2)]/(I_1 + I_2)$

令 $R = I_2/I_1$，则：
$\gamma_{12} = [\exp(ik\theta B/2) + R\cdot\exp(-ik\theta B/2)]/(1 + R)$
$\quad = [(1+R)\cos(k\theta B/2) + i(1-R)\sin(k\theta B/2)]/(1 + R)$

可见度：$V(B) = |\gamma_{12}| = \sqrt{\cos^2(\pi\theta B/\lambda) + ((1-R)/(1+R))^2\sin^2(\pi\theta B/\lambda)}$

分析：
- 在 $B = 0$ 时：$V = 1$（完全相干）
- 第一个最小值在 $\pi\theta B/\lambda = \pi/2$ 处，给出 $\theta = \lambda/2B_{min}$
- 最小值处的可见度：$V_{min} = |1-R|/(1+R)$ 给出 $R$
- 对于等强度恒星 ($R=1$)：$V = |\cos(\pi\theta B/\lambda)|$
</details>

### 练习 16.7：体渲染中的相干性
推导激光扫描显微镜的修正体渲染方程，其中照明具有半径为 $\rho_c$ 的高斯空间相干性。

*提示：从相干体渲染方程开始，并对相干函数进行平均。*

<details>
<summary>解答</summary>

从点 $x'$ 处的相干照明开始：
$L_{coh}(x,\omega) = \int E(x')\sigma_s(x')p(x',\omega'\to\omega)G(x',x) d^3x'$

对于具有互强度 $J(x_1,x_2)$ 的部分相干照明：
$L(x,\omega) = \iint J(x_1,x_2)\sigma_s(x_1)\sigma_s^*(x_2)p(x_1,\omega'\to\omega)p^*(x_2,\omega'\to\omega)G(x_1,x)G^*(x_2,x) d^3x_1d^3x_2$

对于高斯相干性：$J(x_1,x_2) = I(x_1)\exp(-|x_1-x_2|^2/2\rho_c^2)$

在 $\rho_c \to 0$（非相干）的极限下：
$L(x,\omega) = \int I(x')|\sigma_s(x')|^2|p(x',\omega'\to\omega)|^2|G(x',x)|^2 d^3x'$

对于有限的 $\rho_c$，相干体积效应发生在半径 $\rho_c$ 内，导致在距离 $z$ 处产生特征尺寸约为 $\sim\lambda z/\rho_c$ 的散斑图案。
</details>

### 练习 16.8：交叉谱纯度（开放式）
研究部分相干场可以进行谱分解的条件，使得每个频率分量都是完全相干的。哪些物理光源满足此条件？

*提示：考虑当 $W(r_1,r_2,\omega)$ 可以分解为 $U^*(r_1,\omega)U(r_2,\omega)$ 时。*

<details>
<summary>解答</summary>

如果一个场满足以下条件，则它是交叉谱纯的：
$W(r_1,r_2,\omega) = U^*(r_1,\omega)U(r_2,\omega)S(\omega)$

这要求谱相干度：
$\mu(r_1,r_2,\omega) = W(r_1,r_2,\omega)/\sqrt{W(r_1,r_1,\omega)W(r_2,r_2,\omega)} = 1$

对于所有 $S(\omega) \neq 0$ 的 $\omega$。

物理示例：
1.  **滤波后的热光**：白光通过窄带滤波器
2.  **锁模激光器**：每个纵模都是相干的
3.  **稳态 Schell 模型光源**：在特定传播条件下

反例：
- 移动光源（多普勒展宽破坏谱纯度）
- 非线性过程（频率混合）
- 时变介质

对渲染的影响：交叉谱纯光源允许逐频率计算相干效应，大大简化了计算。
</details>

## 常见陷阱和错误

1.  **混淆相干时间和相干长度**
    - 错误：在需要 $l_c$ 的地方使用 $\tau_c$
    - 修正：记住 $l_c = c\cdot\tau_c$，相干长度的单位是距离

2.  **误用 van Cittert-Zernike 定理**
    - 错误：将其用于相干或部分相干光源
    - 修正：该定理仅适用于非相干光源；首先检查光源特性

3.  **相干度归一化不正确**
    - 错误：忘记除以 $\sqrt{I_1I_2}$ 进行归一化
    - 修正：始终计算 $\gamma_{12} = J_{12}/\sqrt{J_{11}J_{22}}$

4.  **复相干性中的符号错误**
    - 错误：在干涉计算中忽略 $\gamma_{12}$ 的相位
    - 修正：注意其复数性质；相位影响条纹位置

5.  **假设相干性得以保留**
    - 错误：将相干性视为传播过程中不变的量
    - 修正：使用 Wolf 方程或适当的传播定律

6.  **混淆时间相干性和空间相干性**
    - 错误：将时间相干性公式用于空间问题
    - 修正：确定是处理时间延迟（时间相干性）还是空间分离（空间相干性）

7.  **准单色近似的误用**
    - 错误：应用于宽带光源
    - 修正：在使用近似之前检查 $\Delta\nu/\nu_0 \ll 1$

8.  **忘记统计平均**
    - 错误：使用瞬时场而不是系综平均
    - 修正：所有相干函数都涉及时间或系综平均 $\langle\cdot\rangle$

## 最佳实践清单

### 设计评审
- [ ] 确定是处理时间相干性、空间相干性还是两者兼有
- [ ] 验证光源相干特性（相干/部分相干/非相干）
- [ ] 选择适当的相干性度量（$\Gamma$, $J$, $W$, $\gamma$）
- [ ] 检查近似的有效性（傍轴、准单色）

### 数学验证
- [ ] 相干函数正确归一化
- [ ] 复共轭位置正确
- [ ] 傅里叶变换对一致（时间 $\leftrightarrow$ 频率，空间 $\leftrightarrow$ 空间频率）
- [ ] 整个过程中的单位量纲正确

### 物理约束
- [ ] $| \gamma_{12} | \le 1$ 处处成立
- [ ] 相干函数是厄米共轭的：$\Gamma(r_1,r_2,\tau) = \Gamma^*(r_2,r_1,-\tau)$
- [ ] 功率谱密度非负：$S(\omega) \ge 0$
- [ ] 传播中的能量守恒

### 计算考量
- [ ] 相干尺度采样充足（$\rho_c$ 的奈奎斯特采样）
- [ ] 振荡积分的数值积分稳定
- [ ] 相干模式展开的模式截断误差估计
- [ ] 蒙特卡洛方法的统计收敛性验证

### 实验验证
- [ ] 相干长度与光谱宽度测量结果一致
- [ ] 可见度测量结果与理论 $\gamma_{12}$ 一致
- [ ] 传播效应与 Wolf 方程预测一致
- [ ] Van Cittert-Zernike 定理对扩展光源的验证


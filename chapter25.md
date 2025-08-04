# 第25章：超材料与变换光学

超材料开启了超越自然材料限制的光学设计新纪元。通过精心设计的亚波长结构，我们能够实现负折射、隐身斗篷等奇异现象。本章将超材料理论与计算机图形学的体积渲染框架相结合，探讨如何在渲染系统中建模和利用这些非凡的光学特性。

## 学习目标

完成本章后，您将能够：
1. 理解超材料的基本原理和负折射率的物理基础
2. 推导变换光学的数学框架并设计隐身器件
3. 分析梯度折射率介质中的光线传播
4. 设计超表面以实现任意相位分布
5. 将超材料概念集成到图形学渲染管线中

## 1. 超材料基础

超材料是由亚波长尺度的人工结构单元周期性排列而成的复合材料，其电磁响应由结构而非成分决定。关键概念是有效介质理论：

### 1.1 有效介电常数和磁导率

对于特征尺寸 $a \ll \lambda$ 的周期结构，可以定义有效参数：

$$\varepsilon_{\text{eff}}(\omega) = \varepsilon_0 \left(1 + \chi_e^{\text{eff}}(\omega)\right)$$

$$\mu_{\text{eff}}(\omega) = \mu_0 \left(1 + \chi_m^{\text{eff}}(\omega)\right)$$

其中 $\chi_e^{\text{eff}}$ 和 $\chi_m^{\text{eff}}$ 是通过均匀化理论计算的有效电磁极化率。

### 1.2 色散关系

在均匀各向同性超材料中，平面波的色散关系为：

$$k^2 = \frac{\omega^2}{c^2} n^2 = \frac{\omega^2}{c^2} \varepsilon_r(\omega) \mu_r(\omega)$$

其中折射率 $n = \pm\sqrt{\varepsilon_r \mu_r}$，符号选择遵循因果性原则。

### 1.3 与体积渲染的联系

在体积渲染框架中，超材料可以表示为空间变化的复折射率场：

$$L(\mathbf{x}, \omega) = \int_0^s T(\mathbf{x}, \mathbf{x}') \sigma_s(\mathbf{x}') L_s(\mathbf{x}', \omega) ds' + T(\mathbf{x}, \mathbf{x}_s) L_0(\mathbf{x}_s, \omega)$$

其中传输函数 $T$ 需要考虑复折射率：

$$T(\mathbf{x}, \mathbf{x}') = \exp\left(-\int_{\mathbf{x}}^{\mathbf{x}'} \sigma_t(\mathbf{s}) + ik_0 \Delta n(\mathbf{s}) ds\right)$$

这里 $\Delta n(\mathbf{s}) = n(\mathbf{s}) - n_0$ 是相对于背景的折射率变化。

## 2. 负折射率材料

### 2.1 理论基础

负折射率材料（NIM）同时具有负的介电常数和磁导率。Veselago在1968年首次预言了这类材料的存在，直到2000年才通过超材料技术实现。

#### 2.1.1 构成关系

对于各向同性NIM，本构关系为：

$$\mathbf{D} = \varepsilon \mathbf{E}, \quad \mathbf{B} = \mu \mathbf{H}$$

其中 $\varepsilon < 0$ 和 $\mu < 0$。折射率定义为：

$$n = -\sqrt{\varepsilon_r \mu_r} < 0$$

负号确保能量沿正确方向传播（Poynting矢量方向）。

#### 2.1.2 逆向波传播

在NIM中，相速度和群速度方向相反：

$$\mathbf{k} \cdot \mathbf{S} < 0$$

其中 $\mathbf{S} = \mathbf{E} \times \mathbf{H}$ 是Poynting矢量。这导致：
- 逆向切伦科夫辐射
- 反常多普勒效应
- 负相速度

### 2.2 负折射定律

在正折射率介质（PIM）和NIM界面上，Snell定律修正为：

$$n_1 \sin\theta_1 = -n_2 \sin\theta_2$$

其中 $n_2 < 0$ 是NIM的折射率。这导致入射光和折射光在法线同侧。

### 2.3 Veselago透镜

平板NIM可以实现完美成像。对于折射率 $n = -1$ 的平板，传递函数为：

$$G(\mathbf{k}_\perp, z) = \begin{cases}
e^{ik_z z} & \text{空气中} \\
e^{-ik_z (z-d)} & \text{NIM中}
\end{cases}$$

其中 $k_z = \sqrt{k_0^2 - |\mathbf{k}_\perp|^2}$。关键是倏逝波（$|\mathbf{k}_\perp| > k_0$）在NIM中被放大而非衰减：

$$|G(\mathbf{k}_\perp, d)| = e^{\sqrt{|\mathbf{k}_\perp|^2 - k_0^2} \cdot d}$$

### 2.4 实现方法

#### 2.4.1 分裂环谐振器（SRR）

磁响应通过SRR阵列实现，其有效磁导率：

$$\mu_{\text{eff}} = 1 - \frac{F\omega^2}{\omega^2 - \omega_0^2 + i\gamma\omega}$$

其中 $F$ 是填充因子，$\omega_0$ 是谐振频率，$\gamma$ 是损耗系数。

#### 2.4.2 金属线阵列

电响应通过亚波长金属线实现，产生等离子体型响应：

$$\varepsilon_{\text{eff}} = 1 - \frac{\omega_p^2}{\omega^2 + i\gamma_e\omega}$$

### 2.5 体积渲染中的NIM

在包含NIM区域的场景中，光线追踪需要修改：

```
1. 界面处理：检测PIM-NIM界面
2. 折射计算：应用修正的Snell定律
3. 相位累积：考虑负相速度
4. 能量守恒：确保Poynting矢量连续
```

传输方程变为：

$$\frac{dL}{ds} = -\sigma_t L + \sigma_s \int_{4\pi} p(\omega, \omega') L(\omega') d\omega' + Q$$

其中路径积分需要考虑相位共轭效应。

## 3. 隐身斗篷原理

### 3.1 变换光学理论

变换光学基于这样的原理：麦克斯韦方程在坐标变换下保持形式不变，但材料参数会相应变化。

#### 3.1.1 坐标变换

考虑从虚拟空间 $(x', y', z')$ 到物理空间 $(x, y, z)$ 的变换：

$$x = x(x', y', z'), \quad y = y(x', y', z'), \quad z = z(x', y', z')$$

Jacobian矩阵：

$$\Lambda_{ij} = \frac{\partial x_i}{\partial x'_j}$$

#### 3.1.2 材料参数变换

在新坐标系中，材料张量变换为：

$$\varepsilon^{ij} = \frac{\Lambda^i_k \Lambda^j_l \varepsilon'^{kl}}{\det(\Lambda)}$$

$$\mu^{ij} = \frac{\Lambda^i_k \Lambda^j_l \mu'^{kl}}{\det(\Lambda)}$$

对于各向同性介质 $\varepsilon' = \varepsilon_0$, $\mu' = \mu_0$，变换后通常得到各向异性张量。

### 3.2 圆柱形隐身斗篷

最简单的隐身斗篷将内半径 $a$ 的圆柱区域映射到外半径 $b$ 的环形区域。

#### 3.2.1 坐标映射

圆柱坐标下的径向压缩变换：

$$r = a + \frac{r'(b-a)}{b}, \quad \theta = \theta', \quad z = z'$$

其中 $r' \in [0, b]$ 映射到 $r \in [a, b]$。

#### 3.2.2 斗篷参数

通过变换光学计算得到的材料参数：

$$\varepsilon_r = \mu_r = \frac{r-a}{r}$$

$$\varepsilon_\theta = \mu_\theta = \frac{r}{r-a}$$

$$\varepsilon_z = \mu_z = \left(\frac{b}{b-a}\right)^2 \frac{r-a}{r}$$

注意在内边界 $r = a$ 处，$\varepsilon_r = \mu_r = 0$ 且 $\varepsilon_\theta = \mu_\theta = \infty$。

### 3.3 完美隐身条件

#### 3.3.1 阻抗匹配

外边界处必须满足阻抗匹配以避免反射：

$$Z(b) = \sqrt{\frac{\mu_r(b)}{\varepsilon_r(b)}} = Z_0$$

#### 3.3.2 相位匹配

光程长度必须保持不变：

$$\int_{\text{斗篷}} n(s) ds = \int_{\text{自由空间}} n_0 ds$$

### 3.4 实际限制

#### 3.4.1 带宽限制

色散关系限制了隐身效果的带宽：

$$n(\omega) = n_0 + \frac{\partial n}{\partial \omega}\bigg|_{\omega_0} (\omega - \omega_0) + O((\omega - \omega_0)^2)$$

#### 3.4.2 损耗效应

实际材料的损耗导致：

$$\varepsilon = \varepsilon' + i\varepsilon''$$

其中 $\varepsilon'' > 0$ 造成能量吸收。

### 3.5 简化设计

#### 3.5.1 地毯式隐身斗篷

通过准共形映射实现近似隐身：

$$\varepsilon = \mu = n^2(x,y) = \left|\frac{\partial w}{\partial z}\right|^2$$

其中 $w = u + iv$ 是复变函数。

#### 3.5.2 多层近似

使用分层各向同性材料近似连续分布：

$$n_i = n(r_i), \quad i = 1, 2, ..., N$$

误差估计：

$$\Delta \phi \sim \frac{(b-a)^2}{N} \max\left|\frac{\partial^2 n}{\partial r^2}\right|$$

## 4. 梯度折射率光学

梯度折射率（GRIN）材料的折射率随空间连续变化，能够实现光线的平滑弯曲而非突变折射。

### 4.1 光线方程

在GRIN介质中，光线路径由光程泛函的极值决定：

$$\delta \int_A^B n(\mathbf{r}) ds = 0$$

导出光线方程：

$$\frac{d}{ds}\left(n\frac{d\mathbf{r}}{ds}\right) = \nabla n$$

### 4.2 光线追踪

#### 4.2.1 Hamilton形式

定义光学Hamilton量：

$$H = \frac{1}{2n^2}|\mathbf{p}|^2$$

其中 $\mathbf{p} = n\frac{d\mathbf{r}}{ds}$ 是光学动量。Hamilton方程：

$$\frac{d\mathbf{r}}{ds} = \frac{\partial H}{\partial \mathbf{p}} = \frac{\mathbf{p}}{n^2}$$

$$\frac{d\mathbf{p}}{ds} = -\frac{\partial H}{\partial \mathbf{r}} = \frac{|\mathbf{p}|^2}{n^3}\nabla n$$

#### 4.2.2 数值积分

使用4阶Runge-Kutta方法：

$$\mathbf{r}_{n+1} = \mathbf{r}_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

其中 $k_i$ 是中间步骤的导数估计。

### 4.3 经典GRIN器件

#### 4.3.1 Luneburg透镜

球对称折射率分布：

$$n(r) = \sqrt{2 - \left(\frac{r}{R}\right)^2}$$

其中 $R$ 是透镜半径。特性：
- 将平行光聚焦到球面上的点
- 无球差
- 广角性能优异

#### 4.3.2 Maxwell鱼眼透镜

折射率分布：

$$n(r) = \frac{n_0}{1 + (r/R)^2}$$

特性：
- 将球面上任意点成像到对径点
- 完美成像（包括倏逝波）
- 共形不变性

#### 4.3.3 自聚焦光纤

抛物线折射率分布：

$$n(r) = n_0\sqrt{1 - 2\Delta(r/a)^2}$$

其中 $\Delta = (n_0^2 - n_1^2)/(2n_0^2)$ 是相对折射率差。

光线轨迹为正弦曲线：

$$r(z) = r_0\cos(gz) + \frac{r'_0}{g}\sin(gz)$$

其中 $g = \sqrt{2\Delta}/a$ 是聚焦参数。

### 4.4 GRIN的体积渲染

在体积渲染中，GRIN效应通过修改光线路径实现：

$$L(\mathbf{x}, \omega) = \int_{\Gamma(\mathbf{x},\omega)} \sigma_s(\mathbf{s}) L_s(\mathbf{s}, \omega(\mathbf{s})) e^{-\tau(\mathbf{x},\mathbf{s})} ds$$

其中 $\Gamma(\mathbf{x},\omega)$ 是弯曲光线路径，$\omega(\mathbf{s})$ 是沿路径变化的方向。

### 4.5 逆向设计

#### 4.5.1 目标函数

给定期望的光线映射 $\mathbf{r}_{\text{out}} = F(\mathbf{r}_{\text{in}})$，寻找折射率分布：

$$J[n] = \int_{\Omega} |F(\mathbf{r}) - F_n(\mathbf{r})|^2 d\mathbf{r}$$

#### 4.5.2 伴随方法

使用拉格朗日乘子法处理约束：

$$\mathcal{L} = J[n] + \int \lambda \left(\frac{d}{ds}\left(n\frac{d\mathbf{r}}{ds}\right) - \nabla n\right) ds$$

梯度计算：

$$\frac{\delta J}{\delta n} = -\int_{\Gamma} \lambda \cdot \nabla\left(\frac{d\mathbf{r}}{ds}\right) ds$$

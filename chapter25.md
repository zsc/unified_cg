# 第25章：超材料与变换光学

超材料开启了超越自然材料限制的光学设计新纪元。通过精心设计的亚波长结构，我们能够实现负折射、隐身斗篷等奇异现象。本章将超材料理论与计算机图形学的体积渲染框架相结合，探讨如何在渲染系统中建模和利用这些非凡的光学特性。超材料的核心在于通过结构而非化学成分来控制电磁响应，这种范式转变为光学设计带来了前所未有的自由度。

## 学习目标

完成本章后，您将能够：
1. 理解超材料的基本原理和负折射率的物理基础
2. 推导变换光学的数学框架并设计隐身器件
3. 分析梯度折射率介质中的光线传播
4. 设计超表面以实现任意相位分布
5. 将超材料概念集成到图形学渲染管线中
6. 评估超材料设计的性能指标和实际限制

## 1. 超材料基础

超材料是由亚波长尺度的人工结构单元（meta-atoms）周期性排列而成的复合材料，其电磁响应由结构几何而非材料成分决定。这种设计理念突破了自然材料的Kramers-Kronig关系限制，实现了前所未有的光学特性。

### 1.1 有效介质理论

当结构特征尺寸远小于工作波长时，电磁波"看不到"单个结构单元，而是感受到均匀化的有效响应。这种均匀化过程的数学基础是多尺度渐近展开。

#### 1.1.1 有效介电常数和磁导率

对于特征尺寸 $a \ll \lambda$ 的周期结构，通过求解单元胞的边值问题可得有效参数：

$$\varepsilon_{\text{eff}}(\omega) = \varepsilon_0 \left(1 + \chi_e^{\text{eff}}(\omega)\right)$$

$$\mu_{\text{eff}}(\omega) = \mu_0 \left(1 + \chi_m^{\text{eff}}(\omega)\right)$$

其中有效极化率通过体积平均获得：

$$\chi_e^{\text{eff}} = \frac{1}{V_{\text{cell}}} \int_{V_{\text{cell}}} \chi_e^{\text{local}}(\mathbf{r}, \omega) d^3\mathbf{r}$$

#### 1.1.2 非局域效应与空间色散

当 $ka \sim 0.1-0.3$ 时，需要考虑空间色散：

$$\mathbf{D}(\mathbf{r}, \omega) = \int \varepsilon(\mathbf{r}-\mathbf{r}', \omega) \mathbf{E}(\mathbf{r}', \omega) d^3\mathbf{r}'$$

在动量空间中：

$$\mathbf{D}(\mathbf{k}, \omega) = \varepsilon(\mathbf{k}, \omega) \cdot \mathbf{E}(\mathbf{k}, \omega)$$

这导致有效介电常数依赖于波矢：$\varepsilon_{\text{eff}} = \varepsilon_{\text{eff}}(\omega, \mathbf{k})$。

### 1.2 色散关系与能带结构

超材料的色散特性决定了其独特的光学响应。

#### 1.2.1 一般色散关系

在均匀各向同性超材料中，平面波解满足：

$$k^2 = \frac{\omega^2}{c^2} n^2(\omega) = \frac{\omega^2}{c^2} \varepsilon_r(\omega) \mu_r(\omega)$$

折射率的符号选择基于因果性和能量守恒：

$$n = \text{sgn}(\varepsilon_r' + \varepsilon_r'') \times \text{sgn}(\mu_r' + \mu_r'') \times \sqrt{|\varepsilon_r \mu_r|}$$

#### 1.2.2 Bloch波与能带

对于周期性超材料，场解具有Bloch形式：

$$\mathbf{E}(\mathbf{r}) = \mathbf{u}_{\mathbf{k}}(\mathbf{r}) e^{i\mathbf{k} \cdot \mathbf{r}}$$

其中 $\mathbf{u}_{\mathbf{k}}(\mathbf{r})$ 具有晶格周期性。能带结构通过求解特征值问题获得：

$$\mathcal{L}_{\mathbf{k}} \mathbf{u}_{\mathbf{k}} = \omega^2(\mathbf{k}) \mathbf{u}_{\mathbf{k}}$$

### 1.3 与体积渲染的深度集成

将超材料响应纳入体积渲染需要扩展传统框架以包含复折射率和非局域效应。

#### 1.3.1 广义传输方程

考虑复折射率场 $n(\mathbf{r}, \omega) = n'(\mathbf{r}, \omega) + in''(\mathbf{r}, \omega)$，辐射传输方程修正为：

$$\frac{\partial L}{\partial s} + ik_0[n(\mathbf{r})-1]L = -\sigma_t L + \int_{4\pi} \sigma_s p(\omega \to \omega') L(\omega') d\omega' + Q$$

其中相位项 $ik_0[n(\mathbf{r})-1]$ 描述了折射引起的相位积累。

#### 1.3.2 传输函数的解析形式

沿光线路径的传输函数包含吸收和相位两部分：

$$T(\mathbf{x}, \mathbf{x}') = \exp\left(-\int_{\mathbf{x}}^{\mathbf{x}'} [\sigma_t(\mathbf{s}) + 2k_0 n''(\mathbf{s})] ds - i\int_{\mathbf{x}}^{\mathbf{x}'} k_0[n'(\mathbf{s})-1] ds\right)$$

振幅衰减因子：
$$A = \exp\left(-\int_{\mathbf{x}}^{\mathbf{x}'} [\sigma_t + 2k_0 n''] ds\right)$$

相位因子：
$$\Phi = \exp\left(-i\int_{\mathbf{x}}^{\mathbf{x}'} k_0[n'-1] ds\right)$$

#### 1.3.3 相干效应的处理

当考虑部分相干光时，需要使用互相干函数：

$$\Gamma(\mathbf{r}_1, \mathbf{r}_2, \tau) = \langle E^*(\mathbf{r}_1, t) E(\mathbf{r}_2, t+\tau) \rangle$$

在超材料中的传播由广义Van Cittert-Zernike定理描述。

### 1.4 超材料的分类与实现

#### 1.4.1 按维度分类

- **3D超材料**：体积型结构，如分裂环谐振器阵列
- **2D超表面**：平面型结构，厚度远小于波长
- **1D超材料**：多层膜结构，仅在一个方向周期性

#### 1.4.2 按功能分类

- **电磁超材料**：控制介电常数 $\varepsilon$
- **磁性超材料**：控制磁导率 $\mu$
- **双各向异性超材料**：电场和磁场耦合
- **手性超材料**：对左右圆偏振光响应不同

### 1.5 设计原理与优化

#### 1.5.1 逆向设计方法

给定目标响应函数 $f_{\text{target}}(\omega)$，通过优化找到结构参数：

$$\min_{\mathbf{p}} \int |\varepsilon_{\text{eff}}(\omega; \mathbf{p}) - \varepsilon_{\text{target}}(\omega)|^2 d\omega + \lambda R(\mathbf{p})$$

其中 $R(\mathbf{p})$ 是正则化项，确保可制造性。

#### 1.5.2 拓扑优化

使用密度方法进行材料分布优化：

$$\rho(\mathbf{r}) \in [0, 1]$$

材料插值使用SIMP (Solid Isotropic Material with Penalization)：

$$\varepsilon(\mathbf{r}) = \varepsilon_{\text{air}} + \rho^p(\mathbf{r})[\varepsilon_{\text{material}} - \varepsilon_{\text{air}}]$$

其中 $p > 1$ 惩罚中间密度值。

## 2. 负折射率材料

负折射率材料代表了超材料最引人注目的成就之一。这类材料同时具有负的介电常数和磁导率，导致电磁波传播的根本性改变。

### 2.1 理论基础

负折射率材料（NIM）的理论基础建立在Veselago 1968年的开创性工作之上。他预言了这类"左手材料"的存在，其中电场、磁场和波矢构成左手坐标系。

#### 2.1.1 本构关系与因果性

对于各向同性NIM，频域本构关系为：

$$\mathbf{D}(\omega) = \varepsilon(\omega) \mathbf{E}(\omega), \quad \mathbf{B}(\omega) = \mu(\omega) \mathbf{H}(\omega)$$

其中 $\text{Re}[\varepsilon(\omega)] < 0$ 和 $\text{Re}[\mu(\omega)] < 0$。折射率的正确定义需要考虑因果性：

$$n(\omega) = \pm\sqrt{\varepsilon_r(\omega) \mu_r(\omega)}$$

符号选择遵循以下原则：
- 对于无损介质：$n = -\sqrt{\varepsilon_r \mu_r}$ 当 $\varepsilon_r < 0, \mu_r < 0$
- 对于有损介质：确保 $\text{Im}[n] > 0$ 以满足能量衰减

#### 2.1.2 逆向波传播

在NIM中，相速度 $\mathbf{v}_p$ 和群速度 $\mathbf{v}_g$ 方向相反：

$$\mathbf{v}_p = \frac{\omega}{k}\hat{\mathbf{k}}, \quad \mathbf{v}_g = \frac{\partial \omega}{\partial \mathbf{k}}$$

能量流动方向由Poynting矢量决定：

$$\mathbf{S} = \frac{1}{2}\text{Re}[\mathbf{E} \times \mathbf{H}^*]$$

在NIM中：$\mathbf{k} \cdot \mathbf{S} < 0$，表明相位传播与能量流动相反。

#### 2.1.3 奇异物理现象

负折射导致多种反常现象：

- **逆向切伦科夫辐射**：辐射锥指向粒子运动的反方向
- **反常多普勒效应**：源接近时频率降低
- **逆向Goos-Hänchen位移**：全反射时的横向位移反向

### 2.2 界面现象与负折射定律

#### 2.2.1 广义Snell定律

在PIM-NIM界面，边界条件要求切向分量连续：

$$\mathbf{k}_1^{\parallel} = \mathbf{k}_2^{\parallel}$$

考虑到NIM中的负相速度，Snell定律修正为：

$$n_1 \sin\theta_1 = -|n_2| \sin\theta_2$$

折射角与入射角在法线同侧，违反了传统的折射直觉。

#### 2.2.2 Fresnel系数

PIM-NIM界面的反射和透射系数需要仔细推导。对于TE偏振：

$$r_{TE} = \frac{n_1\cos\theta_1 + |n_2|\cos\theta_2}{n_1\cos\theta_1 - |n_2|\cos\theta_2}$$

$$t_{TE} = \frac{2n_1\cos\theta_1}{n_1\cos\theta_1 - |n_2|\cos\theta_2}$$

注意系数中出现的符号变化确保了能量守恒。

#### 2.2.3 临界角与全反射

NIM中的临界角条件：

$$\sin\theta_c = \frac{|n_2|}{n_1}$$

当 $|n_2| > n_1$ 时，从PIM到NIM不存在全反射，这与传统情况相反。

### 2.3 Veselago平板透镜

最引人注目的NIM应用是Veselago透镜——一个简单的平行平板能够实现完美成像。

#### 2.3.1 成像原理

对于厚度为 $d$ 的 $n = -1$ 平板，点源的场在频域表示为：

$$\tilde{E}(\mathbf{k}_\perp, z) = \tilde{E}_0(\mathbf{k}_\perp) G(\mathbf{k}_\perp, z)$$

其中Green函数：

$$G(\mathbf{k}_\perp, z) = \begin{cases}
e^{ik_z z} & 0 < z < z_s \\
e^{-ik_z (z-z_s)} e^{ik_z (z-z_1)} & z_s < z < z_s + d \\
e^{ik_z z_s} e^{-ik_z d} e^{ik_z (z-z_s-d)} & z > z_s + d
\end{cases}$$

#### 2.3.2 倏逝波放大

关键创新在于倏逝波的处理。对于 $|\mathbf{k}_\perp| > k_0$：

$$k_z = i\kappa, \quad \kappa = \sqrt{|\mathbf{k}_\perp|^2 - k_0^2}$$

在普通材料中倏逝波指数衰减，但在NIM中：

$$|G(\mathbf{k}_\perp, d)|_{NIM} = e^{\kappa d}$$

这种指数增长补偿了源到透镜的衰减，实现亚波长分辨率。

#### 2.3.3 完美成像条件

完美成像需要满足：
1. 阻抗匹配：$Z_{NIM} = Z_0$
2. 厚度匹配：$d = z_s + z_i$（源距+像距）
3. 无损耗：$\text{Im}[\varepsilon] = \text{Im}[\mu] = 0$

### 2.4 实现方法与结构设计

#### 2.4.1 分裂环谐振器（SRR）

SRR提供磁响应，等效电路模型给出：

$$\mu_{\text{eff}}(\omega) = 1 - \frac{F\omega^2}{\omega^2 - \omega_0^2 + i\gamma\omega}$$

其中：
- $\omega_0 = 1/\sqrt{LC}$ 是LC谐振频率
- $F = \pi r^2/a^2$ 是几何填充因子
- $\gamma = R/L$ 是损耗率

设计参数：
- 环半径 $r \sim \lambda/10$
- 缝隙宽度 $g \sim \lambda/100$
- 金属厚度 $t \sim 3\delta_s$（趋肤深度）

#### 2.4.2 金属线阵列

亚波长金属线产生等离子体型电响应：

$$\varepsilon_{\text{eff}}(\omega) = 1 - \frac{\omega_p^2}{\omega(\omega + i\gamma_e)}$$

等离子体频率：

$$\omega_p = \sqrt{\frac{2\pi c^2}{a^2 \ln(a/r)}}$$

其中 $a$ 是晶格常数，$r$ 是线半径。

#### 2.4.3 复合结构优化

同时实现负 $\varepsilon$ 和负 $\mu$ 需要仔细的参数匹配：

$$\text{FOM} = \frac{-\text{Re}[n]}{|\text{Im}[n]|}$$

优化目标是最大化品质因子（Figure of Merit, FOM）。

### 2.5 体积渲染中的NIM建模

将NIM集成到渲染管线需要根本性的算法修改。

#### 2.5.1 修正的光线追踪

光线-NIM相交需要特殊处理：

```
function traceRayInNIM(ray, nimRegion):
    // 1. 计算入射点和入射角
    hitPoint = intersect(ray, nimRegion.boundary)
    normal = getNormal(nimRegion.boundary, hitPoint)
    cosTheta1 = dot(ray.direction, normal)
    
    // 2. 应用负折射定律
    n_ratio = ray.medium.n / (-abs(nimRegion.n))
    sinTheta2 = n_ratio * sin(acos(cosTheta1))
    
    // 3. 计算折射方向（注意同侧折射）
    refractDir = computeNegativeRefraction(ray.direction, normal, n_ratio)
    
    // 4. 追踪内部路径，考虑逆向相位
    internalRay = Ray(hitPoint, refractDir)
    internalRay.phase_velocity = -1  // 标记负相速度
    
    return propagateInNIM(internalRay, nimRegion)
```

#### 2.5.2 相位累积计算

在NIM中传播时，相位累积方向与传播方向相反：

$$\phi_{total} = \int_{path} k(s) \cdot ds = \int_{PIM} k_0 n_{PIM} ds - \int_{NIM} k_0 |n_{NIM}| ds$$

#### 2.5.3 相干叠加

考虑多路径干涉时，必须正确处理相位：

$$L_{total} = \sum_i A_i \exp(i\phi_i)$$

其中 $\phi_i$ 包含了NIM区域的负相位贡献。

### 2.6 高级效应与应用

#### 2.6.1 非线性NIM

考虑Kerr非线性：

$$n = n_0 + n_2 |E|^2$$

在NIM中，$n_0 < 0$ 可能导致孤子的反常传播。

#### 2.6.2 时变NIM

时间调制的NIM参数：

$$\varepsilon(t) = \varepsilon_0(1 + m\cos(\Omega t))$$

产生频率转换和参量放大效应。

#### 2.6.3 量子效应

NIM中的量子真空涨落导致反常Casimir力，为负的：

$$F_{Casimir} = -\frac{\pi^2 \hbar c}{240 d^4} f(n_1, n_2)$$

其中 $f(n_1, n_2)$ 在NIM情况下改变符号。

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

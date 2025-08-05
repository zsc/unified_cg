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

## 5. 超表面与相位调控

超表面作为超材料的2D版本，通过亚波长结构单元的精确排列实现对电磁波的全方位调控。与体超材料相比，超表面具有低损耗、易集成、设计灵活等优势，已成为平面光学器件的核心技术。

### 5.1 超表面的理论基础

#### 5.1.1 表面等效原理

超表面可视为具有特定边界条件的不连续界面。根据表面等效原理，任意电磁场分布可由等效的表面电流和磁流产生：

$$\mathbf{J}_s = \hat{\mathbf{n}} \times (\mathbf{H}_2 - \mathbf{H}_1)$$

$$\mathbf{M}_s = -\hat{\mathbf{n}} \times (\mathbf{E}_2 - \mathbf{E}_1)$$

其中 $\hat{\mathbf{n}}$ 是从介质1指向介质2的法向量。

#### 5.1.2 表面阻抗模型

对于电薄超表面（厚度 $t \ll \lambda$），可用表面阻抗张量描述：

$$\mathbf{Z}_s = \begin{pmatrix}
Z_{xx} & Z_{xy} \\
Z_{yx} & Z_{yy}
\end{pmatrix}$$

边界条件变为：

$$\mathbf{E}_\parallel^+ - \mathbf{E}_\parallel^- = \mathbf{Z}_s \cdot (\hat{\mathbf{n}} \times \mathbf{H}_{avg})$$

#### 5.1.3 极化率张量

单个超原子的电磁响应由极化率张量描述：

$$\mathbf{p} = \overleftrightarrow{\alpha}_e \cdot \mathbf{E}_{loc}$$

$$\mathbf{m} = \overleftrightarrow{\alpha}_m \cdot \mathbf{H}_{loc}$$

考虑晶格耦合后，有效极化率：

$$\overleftrightarrow{\alpha}_{eff} = \frac{\overleftrightarrow{\alpha}}{1 - \overleftrightarrow{\alpha} \cdot \overleftrightarrow{G}}$$

其中 $\overleftrightarrow{G}$ 是格林函数张量，描述晶格相互作用。

### 5.2 广义Snell定律

超表面通过引入相位梯度扩展了经典折射定律：

$$n_t\sin\theta_t - n_i\sin\theta_i = \frac{\lambda_0}{2\pi}\frac{d\Phi}{dx}$$

其中 $d\Phi/dx$ 是界面上的相位梯度。

#### 5.2.1 动量匹配条件

从动量守恒角度理解，超表面提供额外动量：

$$k_t \sin\theta_t = k_i \sin\theta_i + \frac{d\Phi}{dx}$$

这可以推广到2D情况：

$$\mathbf{k}_t^\parallel = \mathbf{k}_i^\parallel + \nabla_s\Phi$$

其中 $\nabla_s$ 是表面梯度算子。

#### 5.2.2 临界角与异常折射

当相位梯度足够大时，可能出现异常现象：

1. **无临界角全透射**：当 $\frac{d\Phi}{dx} > k_i$ 时
2. **负折射**：适当设计可实现 $\theta_t < 0$
3. **表面波激发**：当 $|\mathbf{k}_t^\parallel| > k_t$ 时

### 5.3 相位分布设计

超表面的核心是实现任意相位分布 $\Phi(x,y)$。主要方法包括：

#### 5.3.1 几何相位（Pancharatnam-Berry相位）

通过旋转各向异性结构单元获得几何相位：

$$\Phi_{PB} = 2\sigma\theta$$

其中 $\sigma = \pm 1$ 是入射圆偏振光的手性，$\theta$ 是结构的方位角。

**Jones矩阵描述**：

旋转角为 $\theta$ 的半波片：

$$\mathbf{J}(\theta) = \begin{pmatrix}
\cos(2\theta) & \sin(2\theta) \\
\sin(2\theta) & -\cos(2\theta)
\end{pmatrix}$$

对于圆偏振基：

$$|L\rangle = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ i \end{pmatrix}, \quad |R\rangle = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -i \end{pmatrix}$$

作用结果：

$$\mathbf{J}(\theta)|L\rangle = e^{-i2\theta}|R\rangle$$
$$\mathbf{J}(\theta)|R\rangle = e^{i2\theta}|L\rangle$$

**优势与限制**：
- 优势：色散小，相位范围完整 $[0, 2\pi]$
- 限制：需要圆偏振光，存在手性转换

#### 5.3.2 传播相位（动态相位）

通过改变结构的有效折射率或高度实现：

$$\Phi_{prop} = k_0 \int_0^h n_{eff}(z) dz \approx k_0 n_{eff} h$$

**谐振结构的相位响应**：

对于Lorentz谐振器：

$$\Phi(\omega) = \arg\left[\frac{1}{\omega_0^2 - \omega^2 - i\gamma\omega}\right]$$

在谐振频率附近：

$$\Phi(\omega) \approx \pi - \arctan\left(\frac{2(\omega-\omega_0)}{\gamma}\right)$$

提供快速相位变化。

#### 5.3.3 混合相位方法

结合几何和传播相位：

$$\Phi_{total} = \Phi_{prop} + \Phi_{PB} = k_0 n_{eff} h + 2\sigma\theta$$

这允许独立控制不同偏振态的相位。

### 5.4 超表面的功能实现

#### 5.4.1 平面透镜

实现聚焦需要的相位分布：

$$\Phi(x,y) = -k_0\left(\sqrt{x^2 + y^2 + f^2} - f\right)$$

其中 $f$ 是焦距。离散化设计：

$$\Phi_{discrete}(x_i, y_j) = \Phi(x_i, y_j) \mod 2\pi$$

**数值孔径与分辨率**：

$$NA = \sin\left(\arctan\frac{D}{2f}\right) \approx \frac{D}{2f}$$

衍射极限分辨率：

$$\Delta x = \frac{0.61\lambda}{NA}$$

#### 5.4.2 光束偏折器

线性相位梯度实现光束偏折：

$$\Phi(x,y) = k_x x + k_y y$$

偏折角：

$$\theta_x = \arcsin\left(\frac{k_x}{k_0 n}\right), \quad \theta_y = \arcsin\left(\frac{k_y}{k_0 n}\right)$$

#### 5.4.3 涡旋光束生成器

轨道角动量光束的相位分布：

$$\Phi(r,\phi) = \ell\phi$$

其中 $\ell$ 是拓扑荷数。总角动量：

$$\mathbf{J} = \mathbf{L} + \mathbf{S} = \ell\hbar\hat{\mathbf{z}} + \sigma\hbar\hat{\mathbf{z}}$$

#### 5.4.4 全息超表面

计算全息的相位分布通过干涉计算：

$$\Phi(x,y) = \arg\left[\sum_j A_j e^{ik|\mathbf{r}-\mathbf{r}_j|}\right]$$

其中 $A_j$ 和 $\mathbf{r}_j$ 是目标图像的振幅和位置。

### 5.5 高效率超表面设计

#### 5.5.1 局域相位匹配

为实现高效率，需要满足局域阻抗匹配：

$$\eta = \frac{|t|^2 \text{Re}[Z_2/Z_1]}{|1 + \Gamma|^2}$$

其中 $\Gamma$ 是反射系数。对于无反射条件：

$$Z_{eff} = \sqrt{Z_1 Z_2}$$

#### 5.5.2 互易性与时间反演对称性

对于无损互易超表面，散射矩阵满足：

$$\mathbf{S}^T = \mathbf{S}, \quad \mathbf{S}^\dagger\mathbf{S} = \mathbf{I}$$

这限制了可实现的功能。打破互易性需要：
- 磁光材料：$\varepsilon_{ij} \neq \varepsilon_{ji}$
- 时变调制：$\varepsilon(t) = \varepsilon_0(1 + m\cos(\Omega t))$
- 非线性效应：$\mathbf{P} = \varepsilon_0(\chi^{(1)}\mathbf{E} + \chi^{(2)}\mathbf{E}\mathbf{E} + ...)$

#### 5.5.3 宽带和宽角设计

色散工程通过多谐振耦合：

$$\alpha(\omega) = \sum_n \frac{f_n}{\omega_n^2 - \omega^2 - i\gamma_n\omega}$$

优化目标函数：

$$F = \int_{\omega_1}^{\omega_2} \int_{\theta_1}^{\theta_2} |\eta(\omega,\theta) - \eta_{target}|^2 d\omega d\theta$$

### 5.6 超表面与体积渲染的集成

将超表面效应纳入体积渲染需要修改边界条件和传输模型。

#### 5.6.1 修正的边界条件

在包含超表面的界面处，辐射度不连续：

$$L(\mathbf{x}^+, \omega_{out}) = \int_{2\pi} f_{ms}(\mathbf{x}, \omega_{in} \to \omega_{out}) L(\mathbf{x}^-, \omega_{in}) |\cos\theta_{in}| d\omega_{in}$$

其中 $f_{ms}$ 是超表面的双向散射分布函数(BSDF)，包含相位调制：

$$f_{ms}(\omega_{in} \to \omega_{out}) = |t(\omega_{in})|^2 \delta(\omega_{out} - \omega_{refract}) e^{i\Phi(\mathbf{x})}$$

#### 5.6.2 相干体积渲染

考虑相位积累的体积渲染方程：

$$\tilde{L}(\mathbf{x}, \omega) = \int_0^s T(t) \sigma_s(\mathbf{x}_t) \tilde{L}_s(\mathbf{x}_t, \omega) e^{i\phi(t)} dt + T(s) \tilde{L}_0 e^{i\phi(s)}$$

其中复振幅 $\tilde{L} = |L|e^{i\phi}$，相位积累：

$$\phi(t) = \int_0^t k(\mathbf{x}_s) ds + \sum_i \Phi_{ms,i}$$

第二项是经过超表面的相位跳变。

#### 5.6.3 多层超表面系统

对于级联的超表面，使用传输矩阵方法：

$$\begin{pmatrix} E^+_N \\ E^-_N \end{pmatrix} = \prod_{i=1}^{N} \mathbf{T}_i \cdot \mathbf{P}_i \begin{pmatrix} E^+_0 \\ E^-_0 \end{pmatrix}$$

其中 $\mathbf{T}_i$ 是第i层超表面的传输矩阵：

$$\mathbf{T}_i = \frac{1}{t_i} \begin{pmatrix} 1 & r_i \\ r_i & 1 \end{pmatrix}$$

$\mathbf{P}_i$ 是传播矩阵：

$$\mathbf{P}_i = \begin{pmatrix} e^{ik_z d_i} & 0 \\ 0 & e^{-ik_z d_i} \end{pmatrix}$$

### 5.7 先进超表面概念

#### 5.7.1 动态可调超表面

通过外场调控实现动态功能：

1. **电调控**：利用液晶、石墨烯或二维材料
   $$n_{eff}(V) = n_0 + \Delta n \cdot f(V)$$

2. **光调控**：利用非线性或相变材料
   $$\varepsilon(I) = \varepsilon_0 + \chi^{(3)} |E|^2$$

3. **热调控**：利用相变材料如VO₂
   $$\varepsilon(T) = \begin{cases}
   \varepsilon_{insulator} & T < T_c \\
   \varepsilon_{metal} & T > T_c
   \end{cases}$$

#### 5.7.2 非局域超表面

考虑空间色散效应：

$$\mathbf{D}(\mathbf{r}, \omega) = \int \overleftrightarrow{\varepsilon}(\mathbf{r}-\mathbf{r}', \omega) \cdot \mathbf{E}(\mathbf{r}', \omega) d^3\mathbf{r}'$$

傅里叶空间表示：

$$\mathbf{D}(\mathbf{k}, \omega) = \overleftrightarrow{\varepsilon}(\mathbf{k}, \omega) \cdot \mathbf{E}(\mathbf{k}, \omega)$$

这允许实现角度选择性和高Q谐振。

#### 5.7.3 拓扑超表面

利用拓扑保护实现鲁棒的光传输：

$$H = \sum_{\mathbf{k}} \mathbf{c}^\dagger_{\mathbf{k}} \mathcal{H}(\mathbf{k}) \mathbf{c}_{\mathbf{k}}$$

其中 $\mathcal{H}(\mathbf{k})$ 是动量空间哈密顿量。拓扑不变量：

$$C = \frac{1}{2\pi} \oint_{BZ} \mathbf{A}(\mathbf{k}) \cdot d\mathbf{k}$$

保证边缘态的存在。

## 6. 超材料在图形学中的应用

超材料概念为计算机图形学带来了新的渲染可能性和设计工具。

### 6.1 超材料BRDF模型

#### 6.1.1 频率依赖BRDF

考虑超材料的色散特性：

$$f_r(\omega_i, \omega_o, \nu) = \frac{|r(\nu)|^2}{4\pi} \frac{\partial \omega_o}{\partial \omega_i} \delta(\omega_o - g(\omega_i, \nu))$$

其中 $g(\omega_i, \nu)$ 描述方向映射，可能包含负折射。

#### 6.1.2 各向异性超材料BRDF

对于各向异性超材料：

$$f_r(\omega_i, \omega_o) = \sum_{m,n} a_{mn} Y_m(\omega_i) Y_n^*(\omega_o)$$

其中 $Y_m$ 是球谐函数，系数 $a_{mn}$ 由材料张量决定：

$$a_{mn} = \int \int \overleftrightarrow{\varepsilon}(\theta, \phi) Y_m(\theta_i, \phi_i) Y_n^*(\theta_o, \phi_o) d\Omega_i d\Omega_o$$

#### 6.1.3 非线性BRDF

强光下的非线性效应：

$$f_r(\omega_i, \omega_o, L_i) = f_r^{(1)} + f_r^{(2)} L_i + f_r^{(3)} L_i^2 + ...$$

这可产生谐波生成等效应。

### 6.2 超材料纹理与着色

#### 6.2.1 程序化超材料纹理

使用噪声函数生成超材料参数分布：

$$\varepsilon(\mathbf{x}) = \varepsilon_0 + \Delta\varepsilon \cdot \text{noise}(\mathbf{x}/\lambda_{pattern})$$

$$\mu(\mathbf{x}) = \mu_0 + \Delta\mu \cdot \text{noise}(\mathbf{x}/\lambda_{pattern} + \mathbf{offset})$$

#### 6.2.2 物理驱动的图案生成

基于物理约束优化超材料图案：

$$\min_{\rho} \int_\Omega |\mathbf{E}_{target} - \mathbf{E}(\rho)|^2 d\Omega + \lambda R(\rho)$$

其中 $\rho(\mathbf{x})$ 是材料密度分布。

#### 6.2.3 多尺度渲染

使用mipmap思想处理多尺度超材料：

$$\varepsilon_{LOD}(\mathbf{x}) = \begin{cases}
\varepsilon_{micro}(\mathbf{x}) & d < d_0 \\
\langle\varepsilon\rangle_{eff} & d > d_1 \\
\text{blend} & d_0 < d < d_1
\end{cases}$$

### 6.3 逆向渲染与超材料设计

#### 6.3.1 目标驱动设计

给定期望的渲染效果，优化超材料参数：

$$\mathbf{p}^* = \arg\min_{\mathbf{p}} \|I_{render}(\mathbf{p}) - I_{target}\|^2 + \lambda \|\mathbf{p}\|^2$$

使用自动微分计算梯度：

$$\frac{\partial I}{\partial \mathbf{p}} = \frac{\partial I}{\partial f_r} \frac{\partial f_r}{\partial \varepsilon} \frac{\partial \varepsilon}{\partial \mathbf{p}}$$

#### 6.3.2 可制造性约束

加入物理可实现性约束：

1. **最小特征尺寸**：$\|\nabla \rho\| < \beta_{max}$
2. **连通性**：确保结构连通
3. **材料限制**：$\varepsilon_{min} < \varepsilon < \varepsilon_{max}$

#### 6.3.3 多目标优化

同时优化多个性能指标：

$$F = w_1 F_{optical} + w_2 F_{mechanical} + w_3 F_{thermal}$$

使用Pareto前沿分析权衡。

### 6.4 实时渲染考虑

#### 6.4.1 GPU加速的超材料着色器

GLSL实现示例结构：

```glsl
vec3 metamaterialBRDF(vec3 wi, vec3 wo, vec2 uv) {
    // 查找超材料参数
    vec2 epsilon = texture(epsilonMap, uv).xy;
    vec2 mu = texture(muMap, uv).xy;
    
    // 计算折射率
    float n = sqrt(epsilon.x * mu.x);
    float k = sqrt(epsilon.y * mu.y);
    
    // 修正的Fresnel项
    vec3 fresnel = complexFresnel(wi, wo, n, k);
    
    // 各向异性响应
    mat3 epsilonTensor = getEpsilonTensor(uv);
    vec3 response = epsilonTensor * fresnel;
    
    return response;
}
```

#### 6.4.2 预计算与缓存

对于复杂的超材料响应，使用查找表：

$$f_r(\theta_i, \phi_i, \theta_o, \phi_o) \approx \sum_{k=1}^K \sigma_k u_k(\theta_i, \phi_i) v_k(\theta_o, \phi_o)$$

SVD分解降低存储需求。

#### 6.4.3 层次细节(LOD)策略

根据观察距离切换模型：

- 近距离：完整超材料模型
- 中距离：有效介质近似
- 远距离：简化BRDF

### 6.5 新兴应用方向

#### 6.5.1 计算成像集成

将超材料光学元件与计算重建结合：

$$I_{sensor} = \mathcal{A}\{I_{scene}\}$$

其中 $\mathcal{A}$ 是包含超材料调制的成像算子。重建：

$$I_{scene} = \arg\min_I \|\mathcal{A}\{I\} - I_{sensor}\|^2 + \lambda R(I)$$

#### 6.5.2 光场调控

使用超表面阵列生成复杂光场：

$$E(\mathbf{r}) = \sum_j A_j(\mathbf{r}) \exp(i\Phi_j(\mathbf{r}))$$

每个超表面贡献特定的振幅和相位分布。

#### 6.5.3 量子渲染集成

考虑超材料的量子光学效应：

$$\hat{\rho}_{out} = \text{Tr}_{env}[\hat{U}(\hat{\rho}_{in} \otimes \hat{\rho}_{env})\hat{U}^\dagger]$$

其中 $\hat{U}$ 包含超材料的量子响应。

## 本章小结

本章深入探讨了超材料和变换光学的基本原理及其在计算机图形学中的应用。我们从有效介质理论出发，建立了超材料的数学框架，详细分析了负折射率材料的奇异性质，包括逆向波传播、完美透镜等现象。通过变换光学方法，我们展示了如何设计隐身斗篷等非凡器件。梯度折射率光学提供了另一种操控光线的途径，而超表面作为2D超材料，通过相位调控实现了平面光学元件的革命。

关键概念总结：
- **有效介质理论**：$\varepsilon_{eff} = \varepsilon_0(1 + \chi_e^{eff})$，当结构尺寸远小于波长时的均匀化描述
- **负折射定律**：$n_1\sin\theta_1 = -|n_2|\sin\theta_2$，折射光线与入射光线在法线同侧
- **变换光学**：$\varepsilon^{ij} = \frac{\Lambda^i_k \Lambda^j_l \varepsilon'^{kl}}{\det(\Lambda)}$，坐标变换导致材料参数变换
- **广义Snell定律**：$n_t\sin\theta_t - n_i\sin\theta_i = \frac{\lambda_0}{2\pi}\frac{d\Phi}{dx}$，超表面引入的相位梯度
- **体积渲染集成**：$T(\mathbf{x}, \mathbf{x}') = \exp(-\int[\sigma_t + 2k_0 n''] ds - i\int k_0[n'-1] ds)$

这些概念不仅拓展了我们对光与物质相互作用的理解，也为图形学渲染提供了新的物理基础和算法工具。

## 练习题

### 基础题

**1. 有效介质计算**
考虑周期为 $a = \lambda/20$ 的二维方形晶格，其中圆形金属柱占据面积分数 $f = 0.3$。假设金属的介电常数为 $\varepsilon_m = -10 + 0.5i$，背景为空气。

a) 使用Maxwell-Garnett公式估算TE偏振的有效介电常数
b) 判断该结构在给定频率下是否表现为等离子体行为
c) 计算有效折射率的实部和虚部

*提示：Maxwell-Garnett公式为 $\frac{\varepsilon_{eff}-\varepsilon_h}{\varepsilon_{eff}+\varepsilon_h} = f\frac{\varepsilon_i-\varepsilon_h}{\varepsilon_i+\varepsilon_h}$*

<details>
<summary>答案</summary>

a) 使用Maxwell-Garnett公式：
$$\frac{\varepsilon_{eff}-1}{\varepsilon_{eff}+1} = 0.3 \times \frac{-10+0.5i-1}{-10+0.5i+1} = 0.3 \times \frac{-11+0.5i}{-9+0.5i}$$

计算复数除法：
$$\frac{-11+0.5i}{-9+0.5i} = \frac{(-11+0.5i)(-9-0.5i)}{81+0.25} = \frac{99+5.5i-4.5i+0.25}{81.25} = \frac{99.25+i}{81.25}$$

因此：
$$\frac{\varepsilon_{eff}-1}{\varepsilon_{eff}+1} = 0.3 \times 1.222 + 0.0123i = 0.367 + 0.0037i$$

解得：
$$\varepsilon_{eff} = \frac{1+0.367+0.0037i}{1-0.367-0.0037i} = 2.15 + 0.012i$$

b) 由于 $\text{Re}[\varepsilon_{eff}] > 0$，该结构不表现为等离子体行为（等离子体需要 $\text{Re}[\varepsilon] < 0$）

c) 有效折射率：
$$n_{eff} = \sqrt{\varepsilon_{eff}} = \sqrt{2.15 + 0.012i} \approx 1.466 + 0.004i$$

实部：$n' = 1.466$，虚部：$n'' = 0.004$
</details>

**2. 负折射界面分析**
光从空气（$n_1 = 1$）以 $45°$ 入射到负折射率材料（$n_2 = -1.5$）。

a) 计算折射角
b) 画出入射、反射和折射光线的示意图
c) 计算TE和TM偏振的功率反射率

*提示：注意负折射时折射角的定义*

<details>
<summary>答案</summary>

a) 使用负折射的Snell定律：
$$n_1\sin\theta_1 = -|n_2|\sin\theta_2$$
$$1 \times \sin(45°) = -(-1.5)\sin\theta_2$$
$$\frac{1}{\sqrt{2}} = 1.5\sin\theta_2$$
$$\sin\theta_2 = \frac{1}{1.5\sqrt{2}} = 0.471$$
$$\theta_2 = 28.1°$$

折射光线在界面法线的入射侧，与入射光线同侧。

b) [示意图应显示入射光线、反射光线在常规位置，但折射光线在法线的同一侧]

c) Fresnel系数（修正版）：
TE偏振：
$$r_{TE} = \frac{n_1\cos\theta_1 + |n_2|\cos\theta_2}{n_1\cos\theta_1 - |n_2|\cos\theta_2} = \frac{0.707 + 1.5 \times 0.882}{0.707 - 1.323} = \frac{2.03}{-0.616} = -3.29$$

功率反射率：$R_{TE} = |r_{TE}|^2 = 10.8$ （大于1表示存在增益，需要考虑能量守恒的完整分析）

实际上，对于无损介质，正确的分析应考虑能流方向，功率反射率应为：
$$R_{TE} = \left|\frac{0.707 - 1.323}{0.707 + 1.323}\right|^2 = 0.092 = 9.2\%$$
</details>

**3. 超表面相位设计**
设计一个工作在 $\lambda = 600$ nm的超表面透镜，焦距 $f = 100\mu m$，直径 $D = 50\mu m$。

a) 推导所需的相位分布函数 $\Phi(r)$
b) 如果使用8个离散相位级（0, π/4, π/2, ..., 7π/4），计算离散化误差
c) 估算该透镜的数值孔径和理论分辨率

*提示：考虑从透镜上任意点到焦点的光程*

<details>
<summary>答案</summary>

a) 相位分布来自等光程条件：
$$\Phi(r) = k_0[\sqrt{r^2 + f^2} - f] = \frac{2\pi}{\lambda}[\sqrt{r^2 + f^2} - f]$$

代入数值：
$$\Phi(r) = \frac{2\pi}{600 \times 10^{-9}}[\sqrt{r^2 + (100 \times 10^{-6})^2} - 100 \times 10^{-6}]$$

在边缘 $r = 25\mu m$：
$$\Phi(25\mu m) = 10.47 \times 10^6 \times [103.08 - 100] \times 10^{-6} = 32.2 \text{ rad}$$

b) 8级离散化，相位量化步长 $\Delta\Phi = \pi/4$。最大量化误差：
$$|\Delta\Phi|_{max} = \pi/8 = 0.393 \text{ rad}$$

相对相位误差：
$$\frac{\Delta\Phi}{\Phi_{max}} = \frac{0.393}{32.2} = 1.2\%$$

c) 数值孔径：
$$NA = \sin\left(\arctan\frac{D/2}{f}\right) = \sin\left(\arctan\frac{25}{100}\right) = \sin(14°) = 0.242$$

理论分辨率（Rayleigh判据）：
$$\Delta x = \frac{0.61\lambda}{NA} = \frac{0.61 \times 600 \times 10^{-9}}{0.242} = 1.51\mu m$$
</details>

### 挑战题

**4. 变换光学隐身斗篷优化**
设计一个工作在可见光波段的简化隐身斗篷，内径 $a = 1\mu m$，外径 $b = 3\mu m$。

a) 推导理想斗篷所需的材料参数张量
b) 如果限制材料参数在 $0.1 < \varepsilon, \mu < 10$ 范围内，设计一个准共形映射实现近似隐身
c) 分析带宽限制和入射角度依赖性

*提示：考虑使用准共形映射 $w = f(z)$ 其中 $|f'(z)|$ 变化较小*

<details>
<summary>答案</summary>

a) 理想圆柱斗篷的坐标变换：
$$r = a + \frac{r'(b-a)}{b}, \quad \theta = \theta', \quad z = z'$$

材料参数：
$$\varepsilon_r = \mu_r = \frac{r-a}{r} = \frac{r-1}{r}$$
$$\varepsilon_\theta = \mu_\theta = \frac{r}{r-a} = \frac{r}{r-1}$$
$$\varepsilon_z = \mu_z = \left(\frac{b}{b-a}\right)^2 \frac{r-a}{r} = \frac{9}{4} \times \frac{r-1}{r}$$

在内边界 $r = 1\mu m$：$\varepsilon_r = 0$，$\varepsilon_\theta = \infty$ （不可实现）
在外边界 $r = 3\mu m$：$\varepsilon_r = 2/3$，$\varepsilon_\theta = 3/2$，$\varepsilon_z = 3/2$

b) 准共形映射设计：
使用映射 $w = z + \frac{c}{z}$，选择 $c$ 使得 $|w'| = |1 - c/z^2|$ 变化较小。

对于 $1 < |z| < 3$，选择 $c = 0.5$：
$$n(r) = |w'| = |1 - 0.5/r^2|$$

在 $r = 1$：$n = 0.5$
在 $r = 3$：$n = 0.944$

调整到允许范围：$n_{adjusted}(r) = 0.5 + 0.444(r-1)/2$

c) 带宽限制：
色散关系导致相位失配：
$$\Delta\phi = \int_a^b [k(\omega) - k(\omega_0)]dr \approx \frac{\partial k}{\partial \omega}\Delta\omega \times 2\mu m$$

对于 $10\%$ 带宽，相位误差约 $\pi/2$，导致隐身效果下降。

角度依赖性：
非垂直入射时，有效路径长度增加：
$$L_{eff} = L_0/\cos\theta$$

当 $\theta > 30°$ 时，隐身效果显著降低。
</details>

**5. 超材料体积渲染算法**
开发一个考虑超材料色散和空间变化的体积渲染算法。

a) 推导包含复折射率的辐射传输方程
b) 设计一个自适应步长的光线行进算法
c) 分析算法的计算复杂度和收敛性

*提示：考虑相位相干性对采样要求的影响*

<details>
<summary>答案</summary>

a) 复折射率辐射传输方程：
$$\frac{dL}{ds} + (\sigma_t + 2k_0n'')L + ik_0(n'-1)L = \sigma_s \int p(\omega' \to \omega)L(\omega')d\omega' + Q$$

分离实部和虚部：
$$\frac{d|L|}{ds} = -(\sigma_t + 2k_0n'')|L| + \text{Re}[\sigma_s \int p L^* d\omega'] + Q$$
$$\frac{d\phi}{ds} = -k_0(n'-1) + \frac{1}{|L|}\text{Im}[\sigma_s \int p L^* d\omega']$$

b) 自适应步长算法：
```
function adaptiveRayMarch(ray, volume):
    L_total = 0
    phase_total = 0
    s = 0
    
    while s < s_max:
        n_complex = volume.getRefractiveIndex(ray.position)
        grad_n = volume.getGradient(ray.position)
        
        // 步长基于相位变化率和折射率梯度
        ds = min(
            lambda / (8 * |n_complex|),  // 相位采样
            0.1 / |grad_n|,              // 梯度采样
            (s_max - s)                  // 剩余距离
        )
        
        // 更新光线方向（GRIN效应）
        ray.direction += ds * grad_n / n_complex
        ray.direction.normalize()
        
        // 积累辐射度和相位
        L_local = volume.emission(ray.position)
        T = exp(-integral_absorption)
        L_total += T * L_local * exp(i * phase_total) * ds
        
        phase_total += k0 * real(n_complex) * ds
        s += ds
    
    return L_total
```

c) 复杂度分析：
- 空间复杂度：$O(N)$，其中 $N$ 是体素数
- 时间复杂度：$O(M \times S)$，其中 $M$ 是光线数，$S$ 是平均步数

收敛性：
误差主要来自相位离散化：
$$\epsilon_{phase} \sim \frac{k_0 \max|n'|\Delta s^2}{2}$$

要求 $\epsilon_{phase} < \pi/4$ 得到步长约束：
$$\Delta s < \sqrt{\frac{\pi}{2k_0\max|n'|}} \approx \frac{\lambda}{4\sqrt{\max|n'|}}$$
</details>

**6. 非线性超表面全息图**
设计一个利用非线性超表面实现的动态全息显示系统。

a) 推导二次谐波生成(SHG)的相位匹配条件
b) 计算实现任意二维图案所需的超表面相位分布
c) 分析泵浦功率对图像质量的影响

*提示：考虑非线性过程的相位匹配和能量守恒*

<details>
<summary>答案</summary>

a) SHG相位匹配条件：
基本过程：$\omega + \omega \to 2\omega$

动量守恒（相位匹配）：
$$\mathbf{k}_{2\omega} = 2\mathbf{k}_\omega + \mathbf{G}$$

其中 $\mathbf{G}$ 是超表面提供的倒格矢。展开：
$$n_{2\omega}(\theta_{out}) = 2n_\omega\cos(\theta_{in}) + \frac{\lambda_\omega}{2\pi}\frac{d\Phi}{dx}$$

对于垂直入射和垂直出射：
$$n_{2\omega} = 2n_\omega + \frac{\lambda_\omega}{2\pi}\frac{d\Phi}{dx}$$

由于通常 $n_{2\omega} > n_\omega$（正常色散），需要：
$$\frac{d\Phi}{dx} = \frac{2\pi}{\lambda_\omega}(n_{2\omega} - 2n_\omega) > 0$$

b) 目标图案的相位分布：
设目标场分布为 $E_{target}(x,y)$，使用Gerchberg-Saxton算法：

1. 初始化：$\Phi_0(x,y) = \text{random}$
2. 迭代：
   - 正向传播：$E_{far} = \mathcal{F}\{A_{in}e^{i\Phi_n}\}$
   - 施加约束：$E'_{far} = |E_{target}|e^{i\arg(E_{far})}$
   - 反向传播：$E_{near} = \mathcal{F}^{-1}\{E'_{far}\}$
   - 更新相位：$\Phi_{n+1} = \arg(E_{near})$

考虑非线性效率：
$$E_{2\omega} \propto \chi^{(2)} E_\omega^2 e^{i\Phi_{NL}(x,y)}$$

因此需要预补偿：
$$\Phi_{design} = \Phi_{target} - 2\Phi_{pump}$$

c) 泵浦功率影响：
SHG效率：
$$\eta_{SHG} = \frac{P_{2\omega}}{P_\omega^2} = \frac{8\pi^2 d_{eff}^2 L^2}{n_\omega^2 n_{2\omega} \epsilon_0 c \lambda_\omega^2} \times \text{sinc}^2(\Delta k L/2)$$

图像质量指标：
- 信噪比：$SNR \propto \sqrt{P_{pump}}$（散粒噪声限制）
- 对比度：随功率增加而饱和
- 空间分辨率：受热效应限制

功率阈值：
- 损伤阈值：$I_{damage} \sim 1 \text{ GW/cm}^2$（脉冲）
- 热透镜效应：$P_{thermal} \sim 100 \text{ mW}$（连续）
- 最优工作点：$P_{opt} \approx 0.1 P_{damage}$
</details>

## 常见陷阱与错误

1. **因果性违反**
   - 错误：在NIM中使用 $n = +\sqrt{\varepsilon\mu}$
   - 正确：考虑损耗确保 $\text{Im}[n] > 0$
   - 调试：检查能流方向与相位传播方向

2. **边界条件处理**
   - 错误：直接应用PIM的Fresnel公式
   - 正确：修正折射角定义和系数符号
   - 调试：验证能量守恒和动量匹配

3. **数值稳定性**
   - 错误：在共振频率直接计算
   - 正确：使用解析延拓或正则化
   - 调试：监控条件数和数值误差

4. **物理可实现性**
   - 错误：设计要求 $\varepsilon = 0$ 或 $\mu = \infty$
   - 正确：使用有界参数的近似设计
   - 调试：检查Kramers-Kronig关系

5. **多尺度问题**
   - 错误：统一采样率处理所有尺度
   - 正确：自适应采样和多尺度方法
   - 调试：分析不同尺度的误差贡献

## 最佳实践检查清单

### 设计阶段
- [ ] 验证设计参数的物理可实现性
- [ ] 检查工作频率下的材料色散
- [ ] 考虑制造公差和误差容限
- [ ] 评估热效应和非线性阈值
- [ ] 确保满足互易性或明确打破机制

### 仿真阶段
- [ ] 选择合适的数值方法（FDTD、FEM、BEM）
- [ ] 验证网格收敛性
- [ ] 正确设置PML或吸收边界
- [ ] 包含材料色散模型
- [ ] 验证能量守恒

### 优化阶段
- [ ] 定义明确的目标函数
- [ ] 包含制造约束
- [ ] 使用适当的正则化
- [ ] 考虑多目标优化
- [ ] 验证局部vs全局最优

### 集成阶段
- [ ] 与现有渲染管线兼容性
- [ ] 性能基准测试
- [ ] 内存使用分析
- [ ] 精度vs速度权衡
- [ ] 用户参数暴露

### 验证阶段
- [ ] 与解析解对比（如适用）
- [ ] 物理量守恒检查
- [ ] 极限情况测试
- [ ] 与实验数据对比
- [ ] 误差敏感性分析

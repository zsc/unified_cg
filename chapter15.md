# 第15章：标量波动光学基础

在本章中，我们将从几何光学和体渲染过渡到波动光学，为理解光作为电磁波建立数学基础。我们将从麦克斯韦方程组推导出标量波近似，并探讨当波长尺度效应变得至关重要时，变得关键的基本衍射现象。这种从基于射线到基于波的描述的转变丰富了我们对光传输的理解，并为计算机图形学中的高级光学现象奠定了基础。

从射线到波的转变从根本上改变了我们对光传输的建模方式。几何光学将光视为遵循直线路径的无穷小射线，而波动光学揭示了光会扩散、绕过边缘衍射并自身干涉。当以下情况发生时，这些效应变得至关重要：
- 特征尺寸接近光的波长（约400-700纳米）
- 存在相干照明（激光、某些LED）
- 需要高保真渲染光学现象
- 微尺度表面结构产生视觉效果

我们将看到体渲染方程如何通过格林函数形式自然地扩展以包含波动现象，提供了一个包含射线和波两种机制的统一框架。

## 数学基础与背景

在深入研究波动方程之前，让我们建立数学背景。从几何光学到波动光学的转变代表了我们描述光传播方式的根本性转变：

**几何光学**：光强度 $I(\mathbf{x},\omega)$ 沿射线遵循：
$\frac{dI}{ds} = -\sigma_t I$ 沿由 $s$ 参数化的射线

**波动光学**：复场振幅 $u(\mathbf{x},t)$ 满足波动方程：
$\nabla^2 u - \frac{1}{c^2}\frac{\partial^2 u}{\partial t^2} = 0$

这种联系通过**程函近似**出现。对于高度振荡的场 $u = A \exp(ikS)$，其中 $k \gg 1$：
- 振幅 $A$ 变化缓慢
- 相位 $S$ 满足程函方程：$|\nabla S|^2 = n^2$
- 射线与等相位面正交

本章探讨当 $k$ 有限时会发生什么，揭示衍射、干涉和光的波动性质。

## 15.1 从麦克斯韦方程组到亥姆霍兹方程

### 矢量波动方程

我们从无源、均匀介质中的麦克斯韦方程组开始：

$\nabla \times \mathbf{E} = -\frac{\partial \mathbf{B}}{\partial t}$ (法拉第定律)
$\nabla \times \mathbf{H} = \frac{\partial \mathbf{D}}{\partial t}$ (安培-麦克斯韦定律)
$\nabla \cdot \mathbf{D} = 0$ (无自由电荷)
$\nabla \cdot \mathbf{B} = 0$ (无磁单极子)

对于线性、各向同性介质：$\mathbf{D} = \varepsilon \mathbf{E}$ 和 $\mathbf{B} = \mu \mathbf{H}$，其中 $\varepsilon = \varepsilon_0 \varepsilon_r$ 和 $\mu = \mu_0 \mu_r$。

这些方程体现了基本的电磁原理：
- **法拉第定律**：时变磁场产生电场
- **安培-麦克斯韦定律**：时变电场（位移电流）和传导电流产生磁场
- **高斯定律**：电场散度与电荷密度相关
- **无单极子**：磁力线总是闭合回路

对法拉第定律取旋度：
$\nabla \times (\nabla \times \mathbf{E}) = -\nabla \times (\frac{\partial \mathbf{B}}{\partial t}) = -\frac{\partial}{\partial t}(\nabla \times \mathbf{B}) = -\mu\frac{\partial}{\partial t}(\nabla \times \mathbf{H})$

使用矢量恒等式 $\nabla \times (\nabla \times \mathbf{E}) = \nabla(\nabla \cdot \mathbf{E}) - \nabla^2 \mathbf{E}$ 并注意到在无源区域 $\nabla \cdot \mathbf{E} = 0$：

$-\nabla^2 \mathbf{E} = -\mu\frac{\partial}{\partial t}(\nabla \times \mathbf{H}) = -\mu\frac{\partial}{\partial t}(\frac{\partial \mathbf{D}}{\partial t}) = -\mu\varepsilon\frac{\partial^2 \mathbf{E}}{\partial t^2}$

这产生了矢量波动方程：

$\nabla^2 \mathbf{E} - \mu\varepsilon(\frac{\partial^2 \mathbf{E}}{\partial t^2}) = 0$

或更紧凑的形式：
$\nabla^2 \mathbf{E} - \frac{1}{v^2}(\frac{\partial^2 \mathbf{E}}{\partial t^2}) = 0$

波速为 $v = \frac{1}{\sqrt{\mu\varepsilon}} = \frac{c}{n}$，其中：
- $c = \frac{1}{\sqrt{\mu_0\varepsilon_0}} \approx 2.998 \times 10^8 \text{ m/s}$ 是真空中的光速
- $n = \sqrt{\varepsilon_r\mu_r} \approx \sqrt{\varepsilon_r}$ 是折射率（因为对于大多数光学材料 $\mu_r \approx 1$）

磁场 $\mathbf{H}$ 也有一个相同的方程。这些矢量方程通过以下方式耦合场的三个空间分量：
1. 材料界面处的**边界条件**
2. **横向约束**：$\nabla \cdot \mathbf{E} = 0$ 意味着对于平面波 $\mathbf{k} \cdot \mathbf{E} = 0$
3. **阻抗关系**：$Z = \sqrt{\mu/\varepsilon}$ 关联 $\mathbf{E}$ 和 $\mathbf{H}$ 的大小

### 数学结构

矢量波动方程表现出几个关键的数学性质：

**线性**：解可以叠加
如果 $\mathbf{E}_1$ 和 $\mathbf{E}_2$ 是解，那么 $\alpha\mathbf{E}_1 + \beta\mathbf{E}_2$ 也是解

**时间反演对称性**：将 $t \rightarrow -t$ 替换会产生有效解
向前和向后传播的波同样有效

**规范不变性**：在真空中，我们可以选择 $\nabla \cdot \mathbf{A} = 0$（库仑规范）或
$\frac{\partial\Phi}{\partial t} + \nabla \cdot \mathbf{A} = 0$（洛伦兹规范），其中 $\mathbf{E} = -\nabla\Phi - \frac{\partial\mathbf{A}}{\partial t}$

**能量守恒**：坡印亭矢量 $\mathbf{S} = \mathbf{E} \times \mathbf{H}$ 满足：
$\frac{\partial u}{\partial t} + \nabla \cdot \mathbf{S} = 0$
其中 $u = \frac{1}{2}(\varepsilon|\mathbf{E}|^2 + \mu|\mathbf{H}|^2)$ 是电磁能量密度

### 标量波近似

对于许多光学现象，我们可以用标量场 $U(\mathbf{r},t)$ 来近似矢量场。这种近似在以下情况下有效：
- 介质在波长尺度上是均匀的（$\nabla n \cdot \lambda \ll n$）
- 偏振效应可以忽略不计（非偏振或固定偏振）
- 场相对于波长变化缓慢（近轴近似）
- 我们远离材料界面，其中边界条件耦合分量

#### 严格推导

为了系统地推导标量近似，我们从每个分量的矢量亥姆霍兹方程开始。考虑一个主要沿 $z$ 传播的波，其电场为：

$\mathbf{E} = E_x \mathbf{\hat{x}} + E_y \mathbf{\hat{y}} + E_z \mathbf{\hat{z}}$

根据麦克斯韦方程组，横向条件 $\nabla \cdot \mathbf{E} = 0$ 给出：
$\frac{\partial E_x}{\partial x} + \frac{\partial E_y}{\partial y} + \frac{\partial E_z}{\partial z} = 0$

对于近轴波，其中 $\frac{\partial}{\partial z} \sim ik$（快速相位变化）但横向导数很小：
$E_z \approx -\frac{1}{ik}(\frac{\partial E_x}{\partial x} + \frac{\partial E_y}{\partial y})$

这表明对于近轴传播，$E_z \ll E_x, E_y$，从而证明了对横向分量的关注是合理的。

每个横向分量满足：
$\nabla^2 E_\perp - \frac{1}{v^2}(\frac{\partial^2 E_\perp}{\partial t^2}) = 0$

对于角频率为 $\omega$ 的单色场：
$E_\perp(\mathbf{r},t) = \text{Re}[e_\perp(\mathbf{r})e^{-i\omega t}]$

复振幅表示将时间和空间分离：
- 时间：$e^{-i\omega t}$，其中 $\frac{\partial^2}{\partial t^2} \rightarrow -\omega^2$
- 空间：$e_\perp(\mathbf{r})$ 包含所有空间变化

代入：
$\nabla^2 e_\perp + (\frac{\omega^2}{v^2})e_\perp = 0$

定义波数 $k = \omega/v = 2\pi n/\lambda$，我们得到：

$\nabla^2 u + k^2 u = 0$

这就是**亥姆霍兹方程**，其中 $u$ 代表场的任何标量分量。

#### 有效性限制

标量近似在以下情况下失效：

1. **高数值孔径**（NA > 0.6）：
   - 矢量效应：纵向场变得显著
   - 紧密聚焦中的偏振耦合
   - 使用矢量衍射理论（Richards-Wolf）

2. **靠近材料界面**：
   - 边界条件耦合场分量
   - 菲涅耳系数取决于偏振
   - 表面等离子体和导模

3. **亚波长结构**：
   - 近场增强
   - 倏逝波占主导
   - 需要完整的矢量处理

4. **双折射介质**：
   - 正交偏振的不同传播
   - 耦合波动方程
   - 需要琼斯或穆勒微积分

### 物理诠释

亥姆霍兹方程描述了单色波传播，其中：
- $k = 2\pi n/\lambda$ 表示空间频率（弧度/米）
- 该方程平衡了空间曲率（$\nabla^2 u$）与相位累积（$k^2 u$）
- 解构成了任意场的完整基

#### 基本解

1. **平面波**：$u = A \exp(i\mathbf{k}\cdot\mathbf{r})$
   - 波矢量 $\mathbf{k} = k(\sin \theta \cos \phi \mathbf{\hat{x}} + \sin \theta \sin \phi \mathbf{\hat{y}} + \cos \theta \mathbf{\hat{z}})$
   - $|\mathbf{k}| = k = 2\pi n/\lambda$
   - 恒定振幅面垂直于 $\mathbf{k}$
   - 能量流沿 $\mathbf{k}$ 方向
   - 角谱表示的基础

   验证：$\nabla^2[\exp(i\mathbf{k}\cdot\mathbf{r})] = -k^2\exp(i\mathbf{k}\cdot\mathbf{r})$ ✓

2. **球面波**：$u = (A/r)\exp(\pm ikr)$
   - 原点处的点源/汇
   - $\pm$ 用于出射/入射波
   - 振幅 $\propto 1/r$（能量守恒）
   - 强度 $\propto 1/r^2$（平方反比定律）
   - 相位面是同心球

   在具有径向对称性的球坐标中：
   $\nabla^2 u = \frac{1}{r^2}\frac{d}{dr}(r^2\frac{du}{dr}) = \frac{A}{r^2}\frac{d}{dr}[r^2\frac{d}{dr}(\frac{1}{r}e^{ikr})]$
   计算后：$\nabla^2 u = -k^2 u$ ✓

3. **高斯光束**：$u = (A_0/q(z))\exp[ikz + ikr^2/2q(z)]$

   复光束参数：$q(z) = z - iz_0$，其中 $z_0 = \pi w_0^2/\lambda$

   光束特性：
   - 光束宽度：$w(z) = w_0\sqrt{1 + (z/z_0)^2}$
   - 波前半径：$R(z) = z(1 + (z_0/z)^2)$
   - 戈伊相移：$\zeta(z) = \arctan(z/z_0)$
   - 瑞利范围：$z_0$（光束面积加倍）

   近轴（$r \ll w$）：满足近轴波动方程
   $\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + 2ik\frac{\partial u}{\partial z} = 0$

4. **贝塞尔光束**：$u = J_0(k_\perp r)\exp(ik_z z)$
   - 非衍射解
   - $k_\perp^2 + k_z^2 = k^2$
   - 无限能量（物理上不可实现）
   - 通过有限孔径近似

5. **厄米-高斯模式**：$u_{mn} = H_m(\sqrt{2}x/w)H_n(\sqrt{2}y/w)\exp(-r^2/w^2)\times[\text{高斯光束因子}]$
   - $H_m$：厄米多项式
   - 正交模式基
   - 矩形对称性
   - 对激光腔很重要

6. **拉盖尔-高斯模式**：$u_{pl} = (r/w)^{|l|} L_p^{|l|}(2r^2/w^2)\exp(-r^2/w^2)\exp(il\phi)\times[\text{高斯光束因子}]$
   - $L_p^{|l|}$：广义拉盖尔多项式
   - 轨道角动量：每个光子 $l\hbar$
   - 圆柱对称性
   - 对于 $l \neq 0$ 的光学涡旋

### 与体渲染的联系

亥姆霍兹方程通过格林函数形式自然地与我们的体渲染框架联系起来。这种联系揭示了波动光学如何从辐射传输方程中产生并扩展。

#### 格林函数公式

考虑相干照明的频域渲染方程：

$L(\mathbf{x},\omega) = L_0(\mathbf{x},\omega) + \int \sigma_s(\mathbf{x}')p(\mathbf{x}',\omega'\rightarrow\omega)G(\mathbf{x},\mathbf{x}')L(\mathbf{x}',\omega')dV'$

格林函数 $G(\mathbf{x},\mathbf{x}')$ 表示从 $\mathbf{x}'$ 到 $\mathbf{x}$ 的相干传播，并满足：

$(\nabla^2 + k^2)G(\mathbf{x},\mathbf{x}') = -\delta(\mathbf{x}-\mathbf{x}')$

这个非齐次亥姆霍兹方程有基本解：
$G(\mathbf{x},\mathbf{x}') = \frac{\exp(ik|\mathbf{x}-\mathbf{x}'|)}{4\pi|\mathbf{x}-\mathbf{x}'|}$

物理诠释：
- 从 $\mathbf{x}'$ 处的点源发出的出射球面波
- 相位累积：$k|\mathbf{x}-\mathbf{x}'|$
- 振幅衰减：$1/|\mathbf{x}-\mathbf{x}'|$
- 满足索末菲辐射条件

#### 尺度相关机制

参数 $kr = k|\mathbf{x}-\mathbf{x}'|$ 决定了传播机制：

1. **几何光学极限**（$kr \gg 1$）：
   - 快速相位振荡：$\exp(ikr)$
   - 驻相近似适用
   - 只有 $\nabla\phi = 0$ 的路径有贡献
   - 格林函数 $\rightarrow$ 沿射线的 $\delta$-函数

   $G(\mathbf{x},\mathbf{x}') \approx \frac{\delta(s - |\mathbf{x}-\mathbf{x}'|)}{|\partial s/\partial \mathbf{x}'|}$

   其中 $s$ 参数化从 $\mathbf{x}'$ 到 $\mathbf{x}$ 的射线

2. **波机制**（$kr \sim 1$）：
   - 相位和振幅相当
   - 衍射效应显著
   - 多条路径之间的干涉
   - 需要完整的格林函数

3. **近场**（$kr \ll 1$）：
   - 静态场近似
   - $G(\mathbf{x},\mathbf{x}') \approx \frac{1}{4\pi|\mathbf{x}-\mathbf{x}'|}$（类库仑）
   - 倏逝波占主导
   - 非传播近场耦合

#### 相干与非相干渲染

相干和非相干渲染之间的转换取决于光源相干性：

**相干光源**（激光、单模光纤）：
- 场叠加：$\mathbf{E}_{\text{total}} = \mathbf{E}_1 + \mathbf{E}_2 + \dots$
- 交叉项：$|\mathbf{E}_{\text{total}}|^2 = |\mathbf{E}_1|^2 + |\mathbf{E}_2|^2 + 2\text{Re}(\mathbf{E}_1^*\mathbf{E}_2) + \dots$
- 干涉条纹，可见度 $V = |\gamma_{12}|$

**部分相干光源**（LED、热源）：
- 互相关函数：$\Gamma(\mathbf{x}_1,\mathbf{x}_2) = \langle \mathbf{E}^*(\mathbf{x}_1)\mathbf{E}(\mathbf{x}_2) \rangle$
- 范西特-泽尼克定理将光源尺寸与相干性关联起来
- 相干长度：$l_c = \lambda^2/(\Delta\lambda)$
- 相干面积：$A_c = \lambda^2 R^2/A_{\text{source}}$

**非相干极限**（大多数渲染）：
- 强度叠加：$I_{\text{total}} = I_1 + I_2 + \dots$
- 无干涉项
- 射线光学足够
- 标准渲染方程适用

#### 扩展体渲染方程

对于部分相干性，渲染方程推广为：

$L(\mathbf{x}_1,\mathbf{x}_2,\omega) = L_0(\mathbf{x}_1,\mathbf{x}_2,\omega) + \iint \sigma_s(\mathbf{x}'_1)\sigma_s(\mathbf{x}'_2)p_1 p_2 G^*(\mathbf{x}_1,\mathbf{x}'_1)G(\mathbf{x}_2,\mathbf{x}'_2)L(\mathbf{x}'_1,\mathbf{x}'_2,\omega)d\mathbf{x}'_1 d\mathbf{x}'_2$

这个6D方程简化为：
- 标准渲染（$\mathbf{x}_1 = \mathbf{x}_2$，仅对角项）
- 相干渲染（可分解的 $L$）
- 散斑/干涉（非对角项）

#### 实际意义

1. **多尺度渲染**：
   - 几何光学：$\lambda \ll$ 特征尺寸
   - 波修正：$\lambda \sim$ 特征尺寸
   - 全波解：$\lambda \gg$ 特征尺寸

2. **统一算法**：
   - 带有相位跟踪的路径追踪
   - 光束传播方法
   - 混合射线-波技术

3. **新现象**：
   - 边缘/孔径衍射
   - 薄膜干涉
   - 粗糙表面散斑
   - 超越几何极限的聚焦和焦散

## 15.2 惠更斯-菲涅耳原理

### 历史发展

克里斯蒂安·惠更斯（1678年）提出波前上的每个点都可视为次级球面波的波源。奥古斯丁-让·菲涅耳（1815年）补充了干涉原理，通过这些波的相干叠加解释了衍射图样。

### 数学公式

考虑在时间 $t$ 的波前 $\Sigma$。在时间 $t + \Delta t$ 处点 $P$ 的场为：

$u(P) = \frac{1}{i\lambda} \iint_\Sigma u(Q) \frac{e^{ikr}}{r} K(\chi) dS$

其中：
- $Q$ 是波前 $\Sigma$ 上的一个点
- $r = |P - Q|$ 是距离
- $K(\chi)$ 是倾斜因子
- $\chi$ 是法线与 $P-Q$ 方向之间的角度

### 倾斜因子

菲涅耳最初提出 $K(\chi) = (1 + \cos \chi)/2$，它：
- 对于向前传播（$\chi = 0$）等于1
- 对于向后传播（$\chi = \pi$）等于0
- 提供平滑的角度依赖性

物理意义：
- $\cos \chi$ 项：波小波在观察方向上的投影
- 常数项：各向同性贡献
- 共同作用：心形辐射图样

这个倾斜因子确保：
1. 没有向后传播的波（因果关系）
2. 向前方向贡献最大
3. 平滑变化防止不连续性
4. 远场能量守恒

### 基尔霍夫严格公式

古斯塔夫·基尔霍夫（1882年）使用格林定理从亥姆霍兹方程推导了惠更斯-菲涅耳原理：

$u(P) = \frac{1}{4\pi} \iint_\Sigma [\frac{e^{ikr}}{r} \frac{\partial u}{\partial n} - u \frac{\partial}{\partial n}(\frac{e^{ikr}}{r})] dS$

对于不透明屏幕中的孔径，入射场为 $u_{\text{inc}}$：
- 在孔径上：$u = u_{\text{inc}}$，$\frac{\partial u}{\partial n} = \frac{\partial u_{\text{inc}}}{\partial n}$
- 在屏幕上：$u = 0$，$\frac{\partial u}{\partial n} = 0$

这产生了**基尔霍夫衍射公式**：

$u(P) = \frac{1}{i\lambda} \iint_{\text{aperture}} u_{\text{inc}}(Q) \frac{e^{ikr}}{r} \frac{(1 + \cos \chi)}{2} dS$

### 与渲染的联系

惠更斯-菲涅耳原理与渲染中的重要性采样并行：

| 波动光学 | 渲染 |
|-------------|----------|
| 次级波源 | 采样点 |
| 波小波叠加 | 蒙特卡洛积分 |
| 倾斜因子 $K(\chi)$ | 余弦加权 ($\mathbf{N}\cdot\mathbf{L}$) |
| 相干叠加 | 复相量和 |
| 强度 = $|\sum \text{fields}|^2$ | 辐射度累积 |

主要区别：
- 波动光学：带有相位的复振幅
- 渲染：实值强度
- 相干性引入了非相干渲染中不存在的干涉效应

这表明渲染的扩展：
1. 相干光源的复值路径追踪
2. 相位感知的重要性采样
3. 材料模型中的干涉效应

## 15.3 菲涅耳衍射积分

### 近场几何

考虑 $z=0$ 平面上的一个平面孔径，由场 $u_0(x_0,y_0)$ 照亮。观察点 $P(x,y,z)$ 处的场由基尔霍夫积分给出。对于近场衍射，我们必须仔细展开出现在振幅和相位项中的距离 $r$。

设 $\mathbf{r} = (x,y,z)$ 为观察点，$\mathbf{r}_0 = (x_0,y_0,0)$ 为孔径中的一个点，则：

$r = |\mathbf{r} - \mathbf{r}_0| = \sqrt{(x-x_0)^2 + (y-y_0)^2 + z^2}$

关键的见解是相位变化比振幅变化快得多：
- 相位变化：$kr \sim 10^6 \text{ rad/m}$ 对于可见光
- 振幅变化：$1/r$ 在波长尺度上变化缓慢

这允许对相位和振幅项采用不同的近似阶数。
### 菲涅尔近似

对于 $z \gg (x-x_0), (y-y_0)$，我们使用二项式展开：

$r = z\sqrt{1 + \frac{(x-x_0)^2 + (y-y_0)^2}{z^2}}$

令 $\rho^2 = (x-x_0)^2 + (y-y_0)^2$。使用 $(1+\varepsilon)^{1/2} \approx 1 + \varepsilon/2 - \varepsilon^2/8 + \dots$ 对于 $\varepsilon \ll 1$：

$r \approx z\left[1 + \frac{\rho^2}{2z^2} - \frac{\rho^4}{8z^4} + \dots\right]$

对于相位项 $kr$，我们保留导致相位误差小于 $\pi/2$ 的项：
- 一阶项：$kz$
- 二阶项：$k\rho^2/2z$ (菲涅尔项)
- 三阶项：$-k\rho^4/8z^3$ (通常忽略)

对于振幅项 $1/r$，我们只保留主导项：
$1/r \approx 1/z$

这得到：
$\frac{e^{ikr}}{r} \approx \left(\frac{e^{ikz}}{z}\right) \exp\left[\frac{ik\rho^2}{2z}\right] = \left(\frac{e^{ikz}}{z}\right) \exp\left[\frac{ik}{2z}((x-x_0)^2 + (y-y_0)^2)\right]$

### 菲涅尔衍射公式

代入基尔霍夫积分：

$u(x,y,z) = \left(\frac{e^{ikz}}{i\lambda z}\right) \iint_{\text{aperture}} u_0(x_0,y_0) \exp\left[\frac{ik}{2z}((x-x_0)^2 + (y-y_0)^2)\right] dx_0dy_0$

展开二次项：

$u(x,y,z) = \left(\frac{e^{ikz}}{i\lambda z}\right) \exp\left[\frac{ik}{2z}(x^2 + y^2)\right] \times$
           $\iint u_0(x_0,y_0) \exp\left[\frac{ik}{2z}(x_0^2 + y_0^2)\right] \exp\left[-\frac{ik}{z}(xx_0 + yy_0)\right] dx_0dy_0$

### 有效性条件

菲涅尔近似的有效性取决于被忽略项引起的相位误差。四次项贡献的相位为：

$\Phi_4 = -\frac{k\rho^4}{8z^3}$

要求 $|\Phi_4|_{\text{max}} < \pi/2$：

$\frac{k\rho^4_{\text{max}}}{8z^3} < \frac{\pi}{2}$

代入 $k = 2\pi/\lambda$ 并求解：

$z^3 > \frac{\rho^4_{\text{max}}}{4\lambda} = \frac{[(x-x_0)^2 + (y-y_0)^2]^2_{\text{max}}}{4\lambda}$

定义**菲涅尔数**：

$F = \frac{a^2}{\lambda z}$

其中 $a$ 是特征孔径尺寸。近似的区域：

1. **$F \gg 1$**：几何阴影 (几何光学)
2. **$F \sim 1$**：菲涅尔衍射 (近场)
3. **$F \ll 1$**：夫琅和费衍射 (远场)

物理诠释：
- $F$ 比较孔径面积 ($a^2$) 与衍射面积 ($\lambda z$)
- 大 $F$：可见多个菲涅尔带，接近几何极限
- 小 $F$：单个菲涅尔带，纯衍射

### 计算方法

#### 1. 直接积分
菲涅尔积分的数值求积：

$u(x,y,z) = \left(\frac{e^{ikz}}{i\lambda z}\right) \iint u_0(x_0,y_0) \exp\left[\frac{ik}{2z}((x-x_0)^2 + (y-y_0)^2)\right] dx_0dy_0$

- 复杂度：对于 $N \times N$ 网格为 $O(N^4)$
- 精确但计算成本高昂
- 适用于不规则孔径或稀疏采样

#### 2. FFT 卷积法
将菲涅尔积分重写为卷积形式：

$u(x,y,z) = C \times [u_0(x,y)\exp(ik(x^2+y^2)/2z)] \otimes \exp(ik(x^2+y^2)/2z)$

其中 $C = \exp(ikz)/(i\lambda z)$ 且 $\otimes$ 表示卷积。

实现：
1. 将输入乘以二次相位 (啁啾)
2. FFT 到频域
3. 乘以传递函数
4. 逆 FFT
5. 乘以输出啁啾

- 复杂度：$O(N^2 \log N)$
- 需要仔细采样以避免混叠
- 需要零填充以提高精度

#### 3. 角谱法
在空间频域中传播场：

$u(x,y,z) = \mathcal{F}^{-1}\{\mathcal{F}\{u_0(x_0,y_0)\} \times H(f_x,f_y,z)\}$

其中传递函数：
$H(f_x,f_y,z) = \exp\left[ikz\sqrt{1-(\lambda f_x)^2-(\lambda f_y)^2}\right]$

对于 $(\lambda f_x)^2 + (\lambda f_y)^2 < 1$：传播波
对于 $(\lambda f_x)^2 + (\lambda f_y)^2 > 1$：倏逝波 (指数衰减)

优点：
- 最有效：$O(N^2 \log N)$
- 在采样限制内精确
- 处理任意传播距离
- 倏逝波的自然处理

采样要求：
$\Delta x < \lambda z/(2X)$ 其中 $X$ 是场范围

## 15.4 夫琅和费衍射与傅里叶光学

### 远场近似

在夫琅和费 (远场) 区域，我们通过假设观察距离 $z$ 足够大来进一步近似菲涅尔积分：

$z \gg k(x_0^2 + y_0^2)_{\text{max}}/2$

这允许我们将二次相位项移到积分号外：

$u(x,y,z) = \left(\frac{e^{ikz}}{i\lambda z}\right) \exp\left[\frac{ik}{2z}(x^2 + y^2)\right] \times$
           $\iint u_0(x_0,y_0) \exp\left[-\frac{ik}{z}(xx_0 + yy_0)\right] dx_0dy_0$

### 傅里叶变换关系

该积分现在是孔径场的一个二维傅里叶变换：

$u(x,y,z) = \left(\frac{e^{ikz}}{i\lambda z}\right) \exp\left[\frac{ik}{2z}(x^2 + y^2)\right] \times \mathcal{F}\{u_0(x_0,y_0)\}|_{f_x=x/\lambda z, f_y=y/\lambda z}$

对于入射到孔径上的平面波 ($u_0 = A(x_0,y_0)$，其中 $A$ 是孔径函数)：

$u(x,y,z) \propto \mathcal{F}\{A(x_0,y_0)\}$

**关键见解**：远场衍射图样是孔径的傅里叶变换。

### 示例

1. **矩形孔径** $A(x_0,y_0) = \text{rect}(x_0/a)\text{rect}(y_0/b)$：
   $u(x,y) \propto \text{sinc}(ax/\lambda z)\text{sinc}(by/\lambda z)$

2. **圆形孔径**，半径为 $a$：
   $u(r,\theta) \propto \frac{2J_1(kar/z)}{kar/z}$
   其中 $J_1$ 是第一类贝塞尔函数。

3. **双缝**，间距为 $d$：
   $u(x) \propto \text{sinc}(ax/\lambda z)\cos(\pi dx/\lambda z)$

### 角谱表示

任何场都可以分解为沿不同方向传播的平面波：

$u(x,y,z) = \iint A(k_x,k_y) \exp[i(k_xx + k_yy + k_rz)] dk_xdk_y$

其中波矢的 $z$ 分量：
$k_r = \sqrt{k^2 - k_x^2 - k_y^2}$

出现两种情况：

1. **传播波** ($k_x^2 + k_y^2 < k^2$)：
   - $k_r$ 是实数
   - 平面波无衰减传播
   - 方向余弦：$(\alpha,\beta,\gamma) = (k_x/k, k_y/k, k_r/k)$
   - 物理角度：$\theta_x = \arcsin(k_x/k)$, $\theta_y = \arcsin(k_y/k)$

2. **倏逝波** ($k_x^2 + k_y^2 > k^2$)：
   - $k_r = i\kappa$ 其中 $\kappa = \sqrt{k_x^2 + k_y^2 - k^2}$
   - 指数衰减：$\exp(-\kappa z)$
   - 局限于近场 ($z \sim 1/\kappa \sim \lambda$)
   - 携带亚波长信息

$z = 0$ 处的角谱：
$A(k_x,k_y) = \left(\frac{1}{2\pi}\right)^2 \iint u(x,y,0) \exp[-i(k_xx + k_yy)] dxdy = \mathcal{F}\{u(x,y,0)\}$

这种表示提供了：
- 场的完整描述
- 自然传播：乘以 $\exp(ik_rz)$
- 与傅里叶光学的直接联系
- 理解分辨率限制的基础

### 与渲染的联系

傅里叶光学框架为渲染提供了深刻的见解：

#### 1. 材料的频率分析
BRDF 在角频率空间中充当传递函数：

- **空间 BRDF**：$\rho(\mathbf{x},\omega_0,\omega_i)$
- **角谱**：$\tilde{\rho}(\mathbf{k},\omega_0,\omega_i) = \mathcal{F}_x\{\rho(\mathbf{x},\omega_0,\omega_i)\}$
- **带宽**：决定所需的采样率

镜面：$\tilde{\rho} \sim \delta(\mathbf{k})$ (所有频率)
漫反射：$\tilde{\rho} \sim \text{sinc}(\mathbf{k})$ (低通)
光泽：中等带宽

#### 2. 采样理论应用

渲染上下文中的**奈奎斯特-香农定理**：
- 空间：$\Delta x < 1/(2f_{\text{max}})$，其中 $f_{\text{max}}$ 是最高空间频率
- 角度：$\Delta \omega < \pi/k_{\text{max}}$ 用于 BRDF 采样
- 时间：$\Delta t < 1/(2f_{\text{motion}})$ 用于运动模糊

**实际意义**：
- 纹理过滤：基于频率内容的 mipmap 级别
- 阴影贴图分辨率：由光频率决定
- 重要性采样：在 $|\tilde{\rho}|$ 较大的地方集中采样

#### 3. 抗锯齿作为滤波

频域中的渲染管线：

1. **场景谱**：$\tilde{S}(\mathbf{k}) = \mathcal{F}\{\text{scene geometry/materials}\}$
2. **采样**：与梳状函数相乘
3. **重建**：与滤波器核卷积
4. **显示**：受像素网格的带宽限制

最佳抗锯齿：
- 预滤波以去除高于奈奎斯特频率的频率
- 常见滤波器：Box (sinc), Gaussian (高斯), Lanczos (窗函数 sinc)
- 权衡：锐度与混叠

#### 4. 光场分析

4D 光场 $L(x,y,u,v)$ 具有 4D 傅里叶变换：
$\tilde{L}(k_x,k_y,k_u,k_v)$

关键见解：
- 朗伯表面：能量集中在 $k_u = k_v = 0$
- 镜面：能量沿 $k_x = \lambda k_u, k_y = \lambda k_v$
- 深度在频域中产生剪切
- 实现最佳采样策略

#### 5. 相干渲染效果

扩展渲染方程以实现相干性：

$L(\mathbf{x},\omega) = L_0(\mathbf{x},\omega) + \int \rho(\mathbf{x},\omega'\to\omega)L(\mathbf{x},\omega')V(\mathbf{x},\mathbf{x}')G(\mathbf{x},\mathbf{x}')d\mathbf{x}'$

其中 $V(\mathbf{x},\mathbf{x}')$ 是互相关函数：
$V(\mathbf{x},\mathbf{x}') = \frac{\langle E^*(\mathbf{x})E(\mathbf{x}')\rangle}{\sqrt{I(\mathbf{x})I(\mathbf{x}')}}$

这使得：
- 激光散斑模拟
- 全息显示
- 薄膜干涉
- 相干次表面散射

建模方法：
- 具有相关长度 $\xi$ 的高度场 $h(x,y)$
- 相位变化：$\phi = 2kh \cos\theta$
- 散斑尺寸：$\Delta x \sim \lambda z/\xi$
- 实现为法线贴图衍射

**虹彩**：
- 薄膜干涉
- 周期性纳米结构产生的结构色
- 波长相关的反射
- 需要基于波的 BRDF 模型

## 15.5 衍射受限成像系统

### 点扩散函数 (PSF)

理想的成像系统将每个物体点映射到唯一的图像点。然而，衍射限制了这种理想行为。点光源的图像是**点扩散函数 (PSF)**。

对于直径为 $D$ 焦距为 $f$ 的圆形孔径，图像平面中的 PSF 为：

$h(r) = |\mathcal{F}\{P(x,y)\}|^2 = \left[\frac{2J_1(\pi Dr/\lambda f)}{\pi Dr/\lambda f}\right]^2$

其中 $P(x,y)$ 是瞳孔函数 (孔径内为 1，孔径外为 0)。

### 艾里斑与分辨率

圆形孔径的 PSF 形成**艾里图样**：

$h(r) = \left[\frac{2J_1(\pi Dr/\lambda f)}{\pi Dr/\lambda f}\right]^2$

特点：
- 中心亮斑 (艾里斑) 包含总能量的 83.8%
- 第一个暗环在 $r_1 = 1.22\lambda f/D$ 处
- 第一个亮环：7.2% 的能量
- 第二个亮环：2.8% 的能量
- 环半径：$r_n \approx (n + 0.22)\lambda f/D$ 对于 $n \ge 1$

艾里斑半径 (第一个零点)：

$r_0 = 1.22\lambda f/D = 1.22\lambda F\#$

其中 $F\# = f/D$ 是光圈数。

**能量分布**：
- 在 $r_0$ 内：83.8%
- 在 $2r_0$ 内：91.0%
- 在 $3r_0$ 内：93.8%

能量集中在中心圆盘中是艾里斑半径作为分辨率实用度量的原因。

### 瑞利判据

当一个艾里斑的最大值落在另一个艾里斑的第一个最小值上时，两个点光源“刚好分辨”：

$\theta_{\text{min}} = 1.22\lambda/D$

这个角分辨率极限是所有成像系统的基本限制。

### 相干与非相干成像

**非相干成像** (自然光典型)：
- 强度相加：$I_{\text{total}} = I_1 + I_2$
- 图像强度 = $|Object|^2 \otimes |PSF|^2$
- 强度线性

**相干成像** (激光照明)：
- 场相加：$U_{\text{total}} = U_1 + U_2$
- 图像场 = $Object \otimes PSF$
- 复振幅线性
- 可表现出干涉效应

### 传递函数

非相干成像的**光学传递函数 (OTF)**：
$\text{OTF}(f) = \mathcal{F}\{|PSF|^2\}$

**调制传递函数 (MTF)**：
$\text{MTF}(f) = |\text{OTF}(f)|$

对于圆形孔径：
$\text{MTF}(\nu) = \frac{2}{\pi}\left[\arccos(\nu) - \nu\sqrt{1-\nu^2}\right]$ 对于 $\nu \le 1$
$\text{MTF}(\nu) = 0$ 对于 $\nu > 1$

其中 $\nu = \lambda f \cdot f_{\text{spatial}}/D$ 是归一化空间频率。

### 对计算机图形学的影响

#### 1. 景深与衍射极限

弥散圆 (CoC) 有两个贡献：
- **几何**：$C_{\text{geom}} = D|z - z_f|/z_f$ (散焦)
- **衍射**：$C_{\text{diff}} = 2.44\lambda F\#$ (艾里斑直径)

总 CoC：$C_{\text{total}} = \sqrt{C_{\text{geom}}^2 + C_{\text{diff}}^2}$

结果：
- 最佳光圈处的最小 CoC：$F\# = \sqrt{|z - z_f|\lambda/(2.44z_f)}$
- 在可见光下，当 $F\# > 8-11$ 时受衍射限制
- 超焦距：$H = f^2/(F\#c) + f$，其中 $c$ 包含衍射

#### 2. 基于物理的焦外成像 (Bokeh)

焦外成像形状取决于：

**几何极限** ($F\# < 5.6$)：
- 形状与光圈几何形状匹配
- 光圈叶片产生锐利边缘
- 均匀的强度分布

**过渡区域** ($F\# \sim 5.6-11$)：
- 衍射使边缘柔化
- 亮度变化：中心更亮
- 卷积：Bokeh = Aperture $\otimes$ Airy

**衍射极限** ($F\# > 11$)：
- 无论光圈形状如何，都是圆形
- 艾里图样占主导
- 高对比度下可见光环

实现方法：
1. 计算几何焦外成像核
2. 与波长相关的艾里函数卷积
3. 对可见光谱求和以获得颜色效果

#### 3. 波光学材料效应

**闪光和闪烁**：
- 由粗糙表面的相干反射引起
- 每个微面产生衍射图样
- 附近微面之间的干涉
- 统计：$I = |E_1 + E_2 + \dots|^2$ 遵循散斑统计

建模方法：
- 具有相关长度 $\xi$ 的高度场 $h(x,y)$
- 相位变化：$\phi = 2kh \cos\theta$
- 散斑尺寸：$\Delta x \sim \lambda z/\xi$
- 实现为法线贴图衍射

**虹彩**：
- 薄膜干涉
- 周期性纳米结构产生的结构色
- 波长相关的反射
- 需要基于波的 BRDF 模型

#### 4. 高级相机模型

**超越薄透镜**：
1. **波前像差**：$\Phi(x,y) = \sum Z_n(x,y)$
   - 泽尼克多项式 $Z_n$ 描述像差
   - $\text{PSF} = |\mathcal{F}\{P(x,y)\exp(ik\Phi(x,y))\}|^2$
   - 空间变化的模糊核

2. **色差效应**：
   - 纵向：焦距 $f(\lambda)$
   - 横向：放大率 $m(\lambda)$
   - PSF 随波长变化
   - 自然色差

3. **偏振**：
   - 菲涅尔系数取决于偏振
   - 透镜系统中的偏振滤光片
   - 带有偏振的天空模型

4. **相干效应**：
   - 扩展光源的部分相干性
   - 相干区域：$A_c \sim \lambda^2 R^2/A_s$
- 影响对比度和分辨率

## 总结

本章建立了波动光学的数学基础，从麦克斯韦方程组过渡到实用的衍射公式。关键概念包括：

1.  **亥姆霍兹方程（Helmholtz Equation）**：$\nabla^2u + k^2u = 0$ - 单色波传播的基本方程
2.  **惠更斯-菲涅尔原理（Huygens-Fresnel Principle）**：每个波前点都可视为次级波源
3.  **菲涅尔衍射（Fresnel Diffraction）**：近场，采用二次相位近似
4.  **夫琅禾费衍射（Fraunhofer Diffraction）**：远场，简化为傅里叶变换
5.  **分辨率极限（Resolution Limits）**：衍射从根本上限制了成像分辨率

光的波动性引入了几何光学之外的现象：
- 干涉和衍射图样
- 基本分辨率极限（瑞利判据）
- 成像中的相干性效应
- 光学系统的频域分析

这些概念通过以下方式与计算机图形学联系起来：
- 带有衍射的物理相机模型
- 基于波的材质外观（闪光、虹彩）
- 渲染算法的傅里叶分析
- 通过格林函数与体渲染的联系

## 练习

### 基本理解（3个问题）

**练习 15.1**：亥姆霍兹方程解
证明 $u(r) = (A/r)\exp(ikr)$ 是球坐标系中亥姆霍兹方程的解。这代表了哪种物理波？

*提示*：对于球对称函数，使用球坐标拉普拉斯算子：$\nabla^2u = (1/r^2)d/dr(r^2du/dr)$。

<details>
<summary>解答</summary>

对于 $u(r) = (A/r)\exp(ikr)$：

$du/dr = A[(-1/r^2)\exp(ikr) + (ik/r)\exp(ikr)] = (A/r)\exp(ikr)[ik - 1/r]$

$r^2du/dr = A \cdot r \cdot \exp(ikr)[ik - 1/r] = A \cdot \exp(ikr)[ikr - 1]$

$d/dr(r^2du/dr) = A[ik \cdot \exp(ikr) \cdot [ikr - 1] + \exp(ikr) \cdot ik]$
                $= A \cdot \exp(ikr)[-k^2r + 2ik]$

$\nabla^2u = (A/r^2)\exp(ikr)[-k^2r + 2ik] = (A/r)\exp(ikr)[-k^2]$

因此：$\nabla^2u + k^2u = (A/r)\exp(ikr)[-k^2 + k^2] = 0 \checkmark$

这代表了从点源发出的向外传播的球面波。
</details>

**练习 15.2**：菲涅尔数
一束平面波（$\lambda = 500\text{nm}$）照射一个半径为 $a = 1\text{mm}$ 的圆形孔径。在什么距离 $z$ 处，菲涅尔数 $F = a^2/\lambda z$ 等于 1？这属于哪种近似区域？

*提示*：菲涅尔近似在 $F \gtrsim 1$ 时有效，而夫琅禾费近似要求 $F \ll 1$。

<details>
<summary>解答</summary>

已知：$\lambda = 500 \times 10^{-9}\text{ m}$，$a = 1 \times 10^{-3}\text{ m}$

$F = a^2/\lambda z = 1$

解 $z$：
$z = a^2/\lambda = (10^{-3})^2 / (500 \times 10^{-9}) = 10^{-6} / (5 \times 10^{-7}) = 2\text{ m}$

在 $z = 2\text{m}$ 处，我们处于菲涅尔（近场）和夫琅禾费（远场）区域的过渡点。对于 $z < 2\text{m}$，使用菲涅尔衍射；对于 $z \gg 2\text{m}$，夫琅禾费近似有效。
</details>

**练习 15.3**：艾里斑大小
一个相机镜头焦距 $f = 50\text{mm}$，光圈直径 $D = 25\text{mm}$（f/2）。计算绿光（$\lambda = 550\text{nm}$）的艾里斑半径。这与典型的像素尺寸相比如何？

*提示*：艾里斑半径为 $r_0 = 1.22\lambda f/D$。

<details>
<summary>解答</summary>

已知：$f = 50\text{mm}$，$D = 25\text{mm}$，$\lambda = 550\text{nm} = 550 \times 10^{-9}\text{ m}$

$r_0 = 1.22\lambda f/D = 1.22 \times (550 \times 10^{-9}) \times (50 \times 10^{-3}) / (25 \times 10^{-3})$
   $= 1.22 \times 550 \times 10^{-9} \times 2$
   $= 1.342 \times 10^{-6}\text{ m} = 1.34\text{ μm}$

直径 $= 2r_0 = 2.68\text{ μm}$

现代相机传感器的像素尺寸为 $1-5\text{ μm}$，因此艾里斑大约跨越 $1-3$ 个像素。这表明许多相机接近衍射极限，尤其是在小光圈下。
</details>

### 高级问题（3个问题）

**练习 15.4**：傅里叶光学与渲染
证明渲染方程在傅里叶域中变为卷积。从以下方程开始：
$L_o(\mathbf{x},\omega_o) = \int \rho(\mathbf{x},\omega_o,\omega_i)L(\mathbf{x},\omega_i)(\omega_o \cdot \mathbf{n})d\omega_i$

*提示*：进行二维空间傅里叶变换并使用卷积定理。

<details>
<summary>解答</summary>

对 $\mathbf{x}$ 进行二维傅里叶变换：

$\mathcal{F}\{L_o(\mathbf{x},\omega_o)\} = \mathcal{F}\{\int \rho(\mathbf{x},\omega_o,\omega_i)L(\mathbf{x},\omega_i)(\omega_o \cdot \mathbf{n})d\omega_i\}$

对于空间不变的 BRDF $\rho(\mathbf{x},\omega_o,\omega_i) = \rho(\omega_o,\omega_i)$：

$\tilde{L}_o(\mathbf{k},\omega_o) = \int \rho(\omega_o,\omega_i)\mathcal{F}\{L(\mathbf{x},\omega_i)\}(\omega_o \cdot \mathbf{n})d\omega_i$
          $= \int \rho(\omega_o,\omega_i)\tilde{L}(\mathbf{k},\omega_i)(\omega_o \cdot \mathbf{n})d\omega_i$

对于纹理表面，其中 $\rho$ 随 $\mathbf{x}$ 变化：

$\tilde{L}_o(\mathbf{k},\omega_o) = \int [\tilde{\rho}(\mathbf{k},\omega_o,\omega_i) \otimes \tilde{L}(\mathbf{k},\omega_i)]  (\omega_o \cdot \mathbf{n})d\omega_i$

这表明空间纹理变化会导致频域卷积，如果采样不当，会导致模糊和混叠。
</details>

**练习 15.5**：基尔霍夫边界条件
使用格林定理从亥姆霍兹方程推导基尔霍夫衍射公式。说明不透明屏幕上的边界条件为何存在问题。

*提示*：使用格林函数 $G = \exp(ikr)/r$ 和格林定理：$\iiint_V (\psi\nabla^2\phi - \phi\nabla^2\psi)dV = \iint_S (\psi\partial\phi/\partial n - \phi\partial\psi/\partial n)dS$

<details>
<summary>解答</summary>

令 $u$ 满足 $(\nabla^2 + k^2)u = 0$，且 $G = \exp(ikr)/r$ 满足 $(\nabla^2 + k^2)G = -4\pi\delta(\mathbf{r})$。

将 $\psi = G$ 和 $\phi = u$ 应用于格林定理：

$\iiint_V [G\nabla^2u - u\nabla^2G]dV = \iint_S [G\partial u/\partial n - u\partial G/\partial n]dS$

由于 $\nabla^2u = -k^2u$ 且 $\nabla^2G = -k^2G - 4\pi\delta(\mathbf{r}-\mathbf{r}_0)$：

$-4\pi u(\mathbf{r}_0) = \iint_S [G\partial u/\partial n - u\partial G/\partial n]dS$

$u(P) = (1/4\pi) \iint_S [\exp(ikr)/r \partial u/\partial n - u \partial/\partial n(\exp(ikr)/r)]dS$

基尔霍夫边界条件假设：
- 在孔径上：$u = u_{\text{incident}}$，$\partial u/\partial n = \partial u_{\text{incident}}/\partial n$
- 在屏幕上：$u = 0$，$\partial u/\partial n = 0$

问题：这些条件在孔径边缘不一致，因为 $u$ 必须从 $u_{\text{incident}}$ 不连续地跳变到 $0$，这违反了波动方程。这就是“基尔霍夫悖论”——尽管存在理论上的不一致，但该近似在实践中效果良好。
</details>

**练习 15.6**：体渲染连接
说明带有散射的体渲染方程在适当的极限下如何简化为惠更斯-菲涅尔原理。考虑：
$L(\mathbf{x},\omega) = \int \sigma_s(\mathbf{x}')p(\mathbf{x}',\omega'\to\omega)G(\mathbf{x},\mathbf{x}')L(\mathbf{x}',\omega')d\mathbf{x}'$

*提示*：考虑一个薄散射层和亥姆霍兹方程的格林函数。

<details>
<summary>解答</summary>

对于单色光，格林函数满足：
$(\nabla^2 + k^2)G(\mathbf{x},\mathbf{x}') = -\delta(\mathbf{x}-\mathbf{x}')$

在自由空间中：$G(\mathbf{x},\mathbf{x}') = \exp(ik|\mathbf{x}-\mathbf{x}'|)/(4\pi|\mathbf{x}-\mathbf{x}'|)$

对于 $z = 0$ 处的薄散射层，其中 $\sigma_s(\mathbf{x}') = \sigma_0\delta(z')A(x',y')$：

$L(x,y,z) = \iint \sigma_0A(x',y')p(\theta)G(\mathbf{x},\mathbf{x}')L_0(x',y')dx'dy'$

对于前向散射 $p(\theta) \approx (1 + \cos \theta)/2$ 和入射场 $L_0$：

$L(x,y,z) = \sigma_0/(4\pi) \iint A(x',y')L_0(x',y') \times$
           $[\exp(ikr)/r][(1 + \cos \chi)/2]dx'dy'$

设置 $\sigma_0/(4\pi) = 1/(i\lambda)$ 即可恢复惠更斯-菲涅尔公式：

$u(P) = (1/i\lambda) \iint u_0(Q)[\exp(ikr)/r]K(\chi)dS$

这表明惠更斯-菲涅尔原理在具有适当散射特性的薄层极限下从体散射中产生。
</details>

### 挑战问题（2个问题）

**练习 15.7**：计算复杂度
比较三种计算菲涅尔衍射图样的方法的计算复杂度：
1.  直接数值积分
2.  基于 FFT 的卷积
3.  角谱法

对于 $N \times N$ 采样网格，推导复杂度并讨论权衡。

*提示*：考虑计算成本和内存需求。

<details>
<summary>解答</summary>

1.  **直接积分**：
    - 对于每个输出点（总共 $N^2$ 个），对 $N^2$ 个输入点进行积分
    - 复杂度：$O(N^4)$
    - 内存：$O(N^2)$
    - 精确但对于大 $N$ 来说速度极慢

2.  **FFT 卷积**：
    - 菲涅尔积分作为与啁啾函数的卷积
    - 步骤：FFT 输入 ($O(N^2\log N)$)，乘法 ($O(N^2)$)，逆 FFT ($O(N^2\log N)$)
    - 复杂度：$O(N^2\log N)$
    - 内存：$O(N^2)$
    - 需要仔细采样以避免混叠

3.  **角谱法**：
    - 在傅里叶域中传播：$H(f_x,f_y) = \exp[ikz\sqrt{1-(\lambda f_x)^2-(\lambda f_y)^2}]$
    - 步骤：FFT ($O(N^2\log N)$)，乘以 $H$ ($O(N^2)$)，逆 FFT ($O(N^2\log N)$)
    - 复杂度：$O(N^2\log N)$
    - 内存：$O(N^2)$
    - 最有效，正确处理倏逝波

权衡：
- 直接：最灵活（任意几何形状）但最慢
- FFT 卷积：快速但可能存在二次相位引起的混叠问题
- 角谱法：对于平面几何形状最快、最精确

对于典型的 $N = 1024$：直接方法需要约 $10^{12}$ 次操作，而 FFT 方法需要约 $10^7$ 次操作。
</details>

**练习 15.8**：统一框架
开发一个统一的数学框架，包含几何光线追踪和波动光学。说明如何根据菲涅尔数在不同区域之间平滑过渡。

*提示*：考虑驻相近似以及光线和波前之间的关系。

<details>
<summary>解答</summary>

**统一框架**：维格纳分布函数（WDF）

WDF $W(\mathbf{x},\mathbf{k})$ 结合了位置和动量（方向）信息：

$W(\mathbf{x},\mathbf{k},z) = \int u^*(\mathbf{x} - \xi/2,z)u(\mathbf{x} + \xi/2,z)\exp(-i\mathbf{k}\cdot\xi)d\xi$

属性：
- 边缘分布给出强度和角谱：$\int W d\mathbf{k} = |u(\mathbf{x})|^2$，$\int W d\mathbf{x} = |\tilde{u}(\mathbf{k})|^2$
- 演化：$\partial W/\partial z + (\mathbf{k}/k_0)\cdot\nabla W = 0$（自由空间）
- 在几何极限下简化为光线密度

**区域过渡**：

定义归一化尺度参数：$\varepsilon = \lambda z/a^2 = 1/F$

1.  **几何光学**（$\varepsilon \to 0$，$F \to \infty$）：
    - WDF $\to$ 光线相空间密度
    - $W(\mathbf{x},\mathbf{k}) = \sum_i \delta(\mathbf{x} - \mathbf{x}_i(z))\delta(\mathbf{k} - \mathbf{k}_i)$
    - 光线追踪有效

2.  **菲涅尔区域**（$\varepsilon \sim 1$，$F \sim 1$）：
    - 二次相位近似
    - W 在相空间中扩散
    - 使用菲涅尔积分

3.  **夫琅禾费区域**（$\varepsilon \gg 1$，$F \ll 1$）：
    - 位置-动量不确定性最大化
    - $W(\mathbf{x},\mathbf{k}) \approx W_0(\mathbf{x})\tilde{W}_0(\mathbf{k})$
    - 傅里叶光学适用

**平滑过渡**：

传播算子：$P(z) = \exp[iz(k^2/2k_0 + \Phi(\mathbf{x},\mathbf{k},\varepsilon))]$

其中 $\Phi$ 插值：
- $\Phi \to 0$ 当 $\varepsilon \to 0$（几何）
- $\Phi \to$ 高阶项 当 $\varepsilon$ 增加

该框架统一了：
- 光线追踪（$\varepsilon \to 0$）
- 高斯光束传播（中间 $\varepsilon$）
- 全波动光学（任意 $\varepsilon$）

WDF 提供了一种相空间表示，可以平滑地在粒子状光线和波状衍射之间过渡，由菲涅尔数控制。
</details>

## 常见陷阱和错误

### 近似有效性
1.  **标量近似**：以下情况无效：
    - 强聚焦（NA > 0.6）
    - 亚波长特征的近场
    - 偏振相关效应

2.  **菲涅尔与夫琅禾费**：
    - 菲涅尔：$F \gtrsim 1$（近场）
    - 夫琅禾费：$F \ll 1$（远场）
    - 过渡区域需要仔细处理

3.  **傍轴近似**：以下情况失效：
    - 大角度（> 15-20°）
    - 宽孔径系统
    - 离轴点

### 数值考虑

1.  **采样要求**：
    - 菲涅尔积分中的二次相位需要密集采样
    - 奈奎斯特准则：$\Delta x < \lambda z/(2X)$，其中 $X$ 是场范围
    - 混叠导致人工条纹

2.  **FFT 伪影**：
    - 周期性边界条件产生环绕效应
    - 精确卷积需要零填充
    - 窗函数减少边缘效应

3.  **相位解缠**：
    - 计算相位限制在 $[-\pi, \pi]$
    - 连续相位需要解缠算法
    - 对噪声和欠采样敏感

4.  **数值精度**：
    - 大 $k$ 值导致精度损失
    - $\exp(ikr)$ 对于大 $r$ 快速振荡
    - 长距离传播使用微分传播

## 最佳实践清单

### 设计审查点

✓ **物理有效性**
- [ ] 指定波长范围
- [ ] 定义相干性属性
- [ ] 考虑偏振效应
- [ ] 需要时包含材料色散

✓ **近似选择**
- [ ] 计算菲涅尔数
- [ ] 选择适当的区域
- [ ] 估计误差范围
- [ ] 识别边缘情况

✓ **数值实现**
- [ ] 采样率满足奈奎斯特准则
- [ ] FFT 大小包含填充
- [ ] 正确处理边界条件
- [ ] 相位计算精度足够

✓ **性能优化**
- [ ] 分析算法复杂度
- [ ] 估计内存需求
- [ ] 考虑 GPU 加速
- [ ] 评估多尺度方法

✓ **验证策略**
- [ ] 验证分析测试用例
- [ ] 检查能量守恒
- [ ] 保持互易性
- [ ] 与几何极限比较

✓ **图形集成**
- [ ] 渲染管线兼容性
- [ ] 评估实时约束
- [ ] 定义细节层次策略
- [ ] 评估感知重要性

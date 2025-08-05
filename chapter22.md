# 第22章：偏振光学基础

本章介绍偏振光学的数学基础，为后续章节中的偏振渲染提供理论框架。我们将系统地探讨偏振态的各种描述方法，包括琼斯矢量、斯托克斯参数和庞加莱球表示，并建立偏振器件的矩阵描述体系。

## 学习目标

完成本章后，您将能够：
1. 使用琼斯矢量和斯托克斯参数描述任意偏振态
2. 运用矩阵方法分析偏振器件的作用
3. 在庞加莱球上可视化偏振态的演化
4. 推导偏振度与相干性之间的关系
5. 设计复合偏振器件系统

## 22.1 偏振态的描述

### 22.1.1 电磁波的偏振现象

光作为横波，其电场矢量 $\mathbf{E}$ 在垂直于传播方向的平面内振动。考虑沿 z 轴传播的单色平面波，电场可表示为：

$$ \mathbf{E}(z,t) = \mathbf{E}_0 \exp[\mathrm{i}(kz - \omega t)] $$

其中 $\mathbf{E}_0 = (Ex_0, Ey_0, 0)$ 是复振幅矢量。偏振描述的是电场矢量端点在 xy 平面内的运动轨迹。

实电场为：
$$ \mathbf{E}^\mathrm{R}(z,t) = \mathrm{Re}[\mathbf{E}(z,t)] = \mathrm{Re}[Ex_0 \exp[\mathrm{i}(kz - \omega t)]]\mathbf{e}_x + \mathrm{Re}[Ey_0 \exp[\mathrm{i}(kz - \omega t)]]\mathbf{e}_y $$

写成分量形式：
$$ Ex(z,t) = |Ex_0| \cos(kz - \omega t + \varphi_x) $$
$$ Ey(z,t) = |Ey_0| \cos(kz - \omega t + \varphi_y) $$

其中 $\varphi_x, \varphi_y$ 是各分量的初相位。相位差 $\delta = \varphi_y - \varphi_x$ 决定了偏振态的性质。

### 22.1.2 线偏振、圆偏振与椭圆偏振

通过消去时间参数，可得到电场矢量端点的轨迹方程：

$$ \left(\frac{Ex}{|Ex_0|}\right)^2 + \left(\frac{Ey}{|Ey_0|}\right)^2 - 2\left(\frac{Ex}{|Ex_0|}\right)\left(\frac{Ey}{|Ey_0|}\right)\cos \delta = \sin^2 \delta $$

这是一般的椭圆方程。根据振幅比和相位差的不同值，可得到不同的偏振态：

**线偏振**（$\delta = 0$ 或 $\pi$）：
- 电场沿固定方向振动
- 偏振方向与 x 轴夹角：$\theta = \arctan(|Ey_0|/|Ex_0|)$
- 轨迹为直线

**圆偏振**（$|Ex_0| = |Ey_0| = E_0$，$\delta = \pm\pi/2$）：
- 右旋圆偏振（RCP）：$\delta = -\pi/2$
- 左旋圆偏振（LCP）：$\delta = +\pi/2$
- 电场矢量端点描绘半径为 $E_0$ 的圆

**椭圆偏振**（一般情况）：
- 电场矢量端点描绘椭圆
- 椭圆长轴与 x 轴夹角：$\psi = (1/2)\arctan\left[\frac{2|Ex_0||Ey_0|\cos \delta}{|Ex_0|^2 - |Ey_0|^2}\right]$
- 椭圆率：$e = \tan \chi$，其中 $\tan 2\chi = \frac{2|Ex_0||Ey_0|\sin \delta}{|Ex_0|^2 - |Ey_0|^2}$

### 22.1.3 偏振态的数学表示

完全偏振光的状态可用以下参数唯一确定：

**几何参数表示**：
- 方位角 $\psi \in [-\pi/2, \pi/2]$：椭圆长轴方向
- 椭圆率角 $\chi \in [-\pi/4, \pi/4]$：描述椭圆的"扁平"程度
- 强度 $I = |Ex_0|^2 + |Ey_0|^2$

**复振幅表示**：
- $Ex_0 = |Ex_0| \exp(\mathrm{i}\varphi_x)$
- $Ey_0 = |Ey_0| \exp(\mathrm{i}\varphi_y)$

两种表示之间的转换关系：
$$ |Ex_0| = \sqrt{I} \cos \psi \cos \chi - \mathrm{i}\sqrt{I} \sin \psi \sin \chi $$
$$ |Ey_0| = \sqrt{I} \sin \psi \cos \chi + \mathrm{i}\sqrt{I} \cos \psi \sin \chi $$

这些关系构成了后续琼斯矢量和斯托克斯参数描述的基础。

### 22.1.4 偏振椭圆的几何构造

**从瞬时电场到偏振椭圆**：

在固定位置 $z = 0$，电场矢量随时间的轨迹形成偏振椭圆。设：
$$ Ex(t) = a_x \cos(\omega t + \varphi_x) $$
$$ Ey(t) = a_y \cos(\omega t + \varphi_y) $$

引入无量纲坐标 $\xi = Ex/a_x$，$\eta = Ey/a_y$，消去时间得椭圆方程：
$$ \xi^2 + \eta^2 - 2\xi\eta \cos \delta = \sin^2 \delta $$

其中 $\delta = \varphi_y - \varphi_x$。

**椭圆参数的推导**：

将椭圆方程转换到主轴坐标系 $(\xi', \eta')$：
$$ \begin{pmatrix} \xi' \\ \eta' \end{pmatrix} = \begin{pmatrix} \cos \psi & \sin \psi \\ -\sin \psi & \cos \psi \end{pmatrix} \begin{pmatrix} \xi \\ \eta \end{pmatrix} $$

主轴方向 $\psi$ 由下式确定：
$$ \tan 2\psi = \frac{2a_x a_y \cos \delta}{a_x^2 - a_y^2} $$

椭圆半长轴 $a$ 和半短轴 $b$：
$$ a^2 + b^2 = a_x^2 + a_y^2 $$
$$ a^2 b^2 = (a_x a_y \sin \delta)^2 $$

椭圆率（ellipticity）：
$$ \varepsilon = b/a = |\tan \chi| $$

其中椭圆率角 $\chi$ 满足：
$$ \sin 2\chi = \frac{2a_x a_y \sin \delta}{a_x^2 + a_y^2} $$

**手性判定**：

偏振光的手性（旋向）由电场矢量的旋转方向决定：
- 右旋（$\delta < 0$）：从光源看，电场顺时针旋转
- 左旋（$\delta > 0$）：从光源看，电场逆时针旋转

手性也可通过角动量密度判定：
$$ \mathbf{L} = \varepsilon_0 \mathbf{E} \times \mathbf{E}/2\mathrm{i}\omega $$

对于圆偏振光，$|\mathbf{L}| = \varepsilon_0 E_0^2/2\omega$。

### 22.1.5 偏振态的能量考虑

**能量密度与Poynting矢量**：

电磁场能量密度：
$$ u = (1/2)(\varepsilon_0|\mathbf{E}|^2 + \mu_0|\mathbf{H}|^2) = \varepsilon_0|\mathbf{E}|^2 $$

时间平均能量密度：
$$ \langle u \rangle = (\varepsilon_0/2)(|Ex_0|^2 + |Ey_0|^2) = (\varepsilon_0/2)I $$

Poynting矢量：
$$ \mathbf{S} = \mathbf{E} \times \mathbf{H} = (1/\mu_0 c)\mathbf{E} \times \mathbf{B} $$

时间平均Poynting矢量：
$$ \langle \mathbf{S} \rangle = (c\varepsilon_0/2)I \mathbf{\hat{z}} $$

注意能量流与偏振态无关，但自旋角动量密度依赖于偏振态。

**偏振态的正交性与完备性**：

两个偏振态 $\mathbf{E}_1$ 和 $\mathbf{E}_2$ 的正交性定义为：
$$ \langle \mathbf{E}_1|\mathbf{E}_2 \rangle = Ex_1^*Ex_2 + Ey_1^*Ey_2 = 0 $$

任意偏振态可在正交基下展开：
$$ \mathbf{E} = \alpha_1 \mathbf{E}_1 + \alpha_2 \mathbf{E}_2 $$

其中展开系数：
$$ \alpha_1 = \langle \mathbf{E}_1|\mathbf{E} \rangle/\langle \mathbf{E}_1|\mathbf{E}_1 \rangle $$
$$ \alpha_2 = \langle \mathbf{E}_2|\mathbf{E} \rangle/\langle \mathbf{E}_2|\mathbf{E}_2 \rangle $$

能量守恒要求：
$$ |\alpha_1|^2\langle \mathbf{E}_1|\mathbf{E}_1 \rangle + |\alpha_2|^2\langle \mathbf{E}_2|\mathbf{E}_2 \rangle = \langle \mathbf{E}|\mathbf{E} \rangle $$

## 22.2 琼斯矢量与琼斯矩阵

### 22.2.1 琼斯矢量的定义

对于完全偏振的单色光，琼斯矢量定义为电场复振幅的归一化列矢量：

$$ \mathbf{J} = \begin{pmatrix} Ex_0 \\ Ey_0 \end{pmatrix} $$

通常采用归一化形式，使得 $|Ex_0|^2 + |Ey_0|^2 = 1$。琼斯矢量完全描述了偏振态，但不包含强度信息。

琼斯矢量的一般形式可写为：
$$ \mathbf{J} = \begin{pmatrix} \cos \theta \exp(\mathrm{i}\varphi_x) \\ \sin \theta \exp(\mathrm{i}\varphi_y) \end{pmatrix} $$

其中 $\theta$ 决定振幅比，$\varphi_x$ 和 $\varphi_y$ 是相位。利用整体相位的任意性，可选择 $\varphi_x = 0$：

$$ \mathbf{J} = \begin{pmatrix} \cos \theta \\ \sin \theta \exp(\mathrm{i}\delta) \end{pmatrix} $$

这里 $\delta = \varphi_y - \varphi_x$ 是相位差。

### 22.2.2 基本偏振态的琼斯表示

**水平线偏振（H）**：
$$ \mathbf{J}_\mathrm{H} = \begin{pmatrix} 1 \\ 0 \end{pmatrix} $$

**垂直线偏振（V）**：
$$ \mathbf{J}_\mathrm{V} = \begin{pmatrix} 0 \\ 1 \end{pmatrix} $$

**45° 线偏振（D）**：
$$ \mathbf{J}_\mathrm{D} = (1/\sqrt{2})\begin{pmatrix} 1 \\ 1 \end{pmatrix} $$

**-45° 线偏振（A）**：
$$ \mathbf{J}_\mathrm{A} = (1/\sqrt{2})\begin{pmatrix} 1 \\ -1 \end{pmatrix} $$

**右旋圆偏振（R）**：
$$ \mathbf{J}_\mathrm{R} = (1/\sqrt{2})\begin{pmatrix} 1 \\ -\mathrm{i} \end{pmatrix} $$

**左旋圆偏振（L）**：
$$ \mathbf{J}_\mathrm{L} = (1/\sqrt{2})\begin{pmatrix} 1 \\ \mathrm{i} \end{pmatrix} $$

这六个状态构成三组正交基：$\{\mathrm{H},\mathrm{V}\}$、$\{\mathrm{D},\mathrm{A}\}$、$\{\mathrm{R},\mathrm{L}\}$。任意偏振态可在任一组基下展开。

### 22.2.3 琼斯矩阵方法

偏振器件对琼斯矢量的作用可用 $2\times2$ 复矩阵表示：

$$ \mathbf{J}_\mathrm{out} = \mathbf{M} \mathbf{J}_\mathrm{in} $$

其中 $\mathbf{M}$ 是琼斯矩阵。对于无损器件，矩阵应满足幺正性：$\mathbf{M}^\dagger\mathbf{M} = \mathbf{I}$。

**基本琼斯矩阵**：

线偏振器（透振方向与 x 轴夹角 $\theta$）：
$$ \mathbf{M}_\mathrm{pol}(\theta) = \begin{pmatrix} \cos^2\theta & \cos \theta \sin \theta \\ \cos \theta \sin \theta & \sin^2\theta \end{pmatrix} $$

相位延迟器（快轴沿 x，相位延迟 $\delta$）：
$$ \mathbf{M}_\mathrm{ret}(\delta) = \begin{pmatrix} 1 & 0 \\ 0 & \exp(\mathrm{i}\delta) \end{pmatrix} $$

旋转器（旋转角 $\theta$）：
$$ \mathbf{M}_\mathrm{rot}(\theta) = \begin{pmatrix} \cos \theta & -\sin \theta \\ \sin \theta & \cos \theta \end{pmatrix} $$

### 22.2.4 级联系统的矩阵运算

多个器件级联时，总的琼斯矩阵为各矩阵的乘积（注意顺序）：

$$ \mathbf{M}_\mathrm{total} = \mathbf{M}_n \mathbf{M}_{n-1} \dots \mathbf{M}_2 \mathbf{M}_1 $$

**例：四分之一波片 + 线偏振器**

设四分之一波片快轴沿 x（$\delta = \pi/2$），后接水平偏振器：

$$ \mathbf{M}_\mathrm{total} = \mathbf{M}_\mathrm{pol}(0) \mathbf{M}_\mathrm{ret}(\pi/2) = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & \mathrm{i} \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} $$

输入右旋圆偏振光：
$$ \mathbf{J}_\mathrm{out} = \mathbf{M}_\mathrm{total} \mathbf{J}_\mathrm{R} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}(1/\sqrt{2})\begin{pmatrix} 1 \\ -\mathrm{i} \end{pmatrix} = (1/\sqrt{2})\begin{pmatrix} 1 \\ 0 \end{pmatrix} $$

输出为水平线偏振光，强度减半。

**任意取向相位延迟器的矩阵**：

快轴与 x 轴夹角为 $\theta$ 的相位延迟器：
$$ \mathbf{M}(\theta,\delta) = \mathbf{R}(-\theta) \mathbf{M}_\mathrm{ret}(\delta) \mathbf{R}(\theta) $$

展开得：
$$ \mathbf{M}(\theta,\delta) = \begin{pmatrix} \cos^2\theta + \sin^2\theta \exp(\mathrm{i}\delta) & \sin \theta \cos \theta(1 - \exp(\mathrm{i}\delta)) \\ \sin \theta \cos \theta(1 - \exp(\mathrm{i}\delta)) & \sin^2\theta + \cos^2\theta \exp(\mathrm{i}\delta) \end{pmatrix} $$

这个公式在分析复杂偏振系统时非常有用。

### 22.2.5 琼斯演算的应用实例

**偏振态转换设计**：

问题：设计一个系统将水平线偏振转换为左旋圆偏振。

解：需要找到矩阵 $\mathbf{M}$ 使得：
$$ \mathbf{M}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = (1/\sqrt{2})\begin{pmatrix} 1 \\ \mathrm{i} \end{pmatrix} $$

一种解决方案：45°取向的四分之一波片
$$ \mathbf{M} = \mathbf{R}(-\pi/4)\mathbf{M}_\mathrm{ret}(\pi/2)\mathbf{R}(\pi/4) = (1/\sqrt{2})\begin{pmatrix} 1+\mathrm{i} & 1-\mathrm{i} \\ 1-\mathrm{i} & 1+\mathrm{i} \end{pmatrix} $$

验证：$\mathbf{M}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = (1/\sqrt{2})\begin{pmatrix} 1+\mathrm{i} \\ 1-\mathrm{i} \end{pmatrix} = (1/\sqrt{2})\begin{pmatrix} 1 \\ \mathrm{i} \end{pmatrix} \exp(\mathrm{i}\pi/4)$

相位因子不影响偏振态。

**偏振模式转换器**：

光纤通信中常需要在TE和TM模式间转换。使用半波片：

TE→TM：
$$ \mathbf{M}_\mathrm{HWP}(\pi/4) = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} $$

这实现了x和y分量的交换。

**消偏振器设计**：

理想消偏振器将任意输入偏振态转换为非偏振光。但琼斯矩阵无法描述这一过程，需要使用Mueller矩阵。

伪消偏振器（Lyot消偏振器）：
- 两个延迟量不同的波片级联
- $\delta_1 = 2\pi$，$\delta_2 = \pi$
- 对窄带光源产生快速变化的输出偏振态

### 22.2.6 琼斯矩阵的数学性质

**特征值与特征矢量**：

任意琼斯矩阵 $\mathbf{M}$ 的特征方程：
$$ \det(\mathbf{M} - \lambda\mathbf{I}) = 0 $$

特征矢量代表不变偏振态（eigenpolarizations）。

例：四分之一波片（快轴0°）
特征值：$\lambda_1 = 1, \lambda_2 = \mathrm{i}$
特征矢量：$\mathbf{v}_1 = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$（水平偏振），$\mathbf{v}_2 = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$（垂直偏振）

**矩阵分解**：

任意琼斯矩阵可分解为：
$$ \mathbf{M} = \mathbf{U}\Lambda\mathbf{U}^\dagger $$

其中 $\mathbf{U}$ 是幺正矩阵，$\Lambda$ 是对角矩阵。

物理意义：
1. $\mathbf{U}^\dagger$：转换到特征偏振基
2. $\Lambda$：在特征基下的相位延迟
3. $\mathbf{U}$：转换回原始基

**群论性质**：

无损偏振器件的琼斯矩阵构成SU(2)群：
- 行列式为1（能量守恒）
- 幺正性（可逆过程）
- 群运算：矩阵乘法

SU(2)的生成元（Pauli矩阵）：
$$ \sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \quad \sigma_2 = \begin{pmatrix} 0 & -\mathrm{i} \\ \mathrm{i} & 0 \end{pmatrix} \quad \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix} $$

任意SU(2)元素：
$$ \mathbf{M} = \exp(\mathrm{i}\mathbf{n}\cdot\mathbf{\sigma}\theta/2) = \cos(\theta/2)\mathbf{I} + \mathrm{i} \sin(\theta/2)\mathbf{n}\cdot\mathbf{\sigma} $$

这与量子力学中的旋转算符同构。

## 22.3 斯托克斯参数

### 22.3.1 斯托克斯参数的定义

斯托克斯参数提供了描述任意偏振态（包括部分偏振光）的完整方法。定义四个实参数：

$$ S_0 = \langle|Ex|^2\rangle + \langle|Ey|^2\rangle = I \text{（总强度）} $$
$$ S_1 = \langle|Ex|^2\rangle - \langle|Ey|^2\rangle $$
$$ S_2 = 2\langle Ex^*Ey \rangle^\mathrm{R} = 2\langle|Ex||Ey|\cos \delta\rangle $$
$$ S_3 = 2\langle Ex^*Ey \rangle^\mathrm{I} = 2\langle|Ex||Ey|\sin \delta\rangle $$

其中 $\langle\cdot\rangle$ 表示时间平均或统计平均。斯托克斯矢量定义为：

$$ \mathbf{S} = \begin{pmatrix} S_0 \\ S_1 \\ S_2 \\ S_3 \end{pmatrix} $$

对于完全偏振光，有约束条件：
$$ S_0^2 = S_1^2 + S_2^2 + S_3^2 $$

### 22.3.2 部分偏振光的描述

部分偏振光可分解为完全偏振部分和完全非偏振部分：

$$ \mathbf{S} = \mathbf{S}_\mathrm{pol} + \mathbf{S}_\mathrm{unpol} $$

其中：
$$ \mathbf{S}_\mathrm{unpol} = \begin{pmatrix} S_0 - \sqrt{S_1^2 + S_2^2 + S_3^2} \\ 0 \\ 0 \\ 0 \end{pmatrix} $$

$$ \mathbf{S}_\mathrm{pol} = \begin{pmatrix} \sqrt{S_1^2 + S_2^2 + S_3^2} \\ S_1 \\ S_2 \\ S_3 \end{pmatrix} $$

偏振度定义为：
$$ P = \frac{\sqrt{S_1^2 + S_2^2 + S_3^2}}{S_0} $$

满足 $0 \le P \le 1$，其中 $P = 0$ 表示完全非偏振光，$P = 1$ 表示完全偏振光。

### 22.3.3 偏振度的计算

**与相干性的关系**：

对于准单色光，偏振度与相干度密切相关：

$$ P^2 = \frac{S_1^2 + S_2^2 + S_3^2}{S_0^2} = \frac{4|\langle Ex^*Ey \rangle|^2}{(\langle|Ex|^2\rangle + \langle|Ey|^2\rangle)^2} $$

当 $Ex$ 和 $Ey$ 完全相干时，$P = 1$；完全不相干时，$P = 0$。

**偏振态参数的提取**：

从斯托克斯参数可提取偏振椭圆的参数：

方位角：$\psi = (1/2)\arctan(S_2/S_1)$
椭圆率角：$\chi = (1/2)\arcsin\left(S_3/\sqrt{S_1^2 + S_2^2 + S_3^2}\right)$

对于部分偏振光，这些参数描述其偏振部分的性质。

### 22.3.4 测量方法

斯托克斯参数可通过一系列强度测量获得：

**基本测量方案**：

$S_0 = I_\mathrm{H} + I_\mathrm{V}$（水平 + 垂直偏振器）
$S_1 = I_\mathrm{H} - I_\mathrm{V}$
$S_2 = I_\mathrm{D} - I_\mathrm{A}$（+45° 和 -45° 偏振器）
$S_3 = I_\mathrm{R} - I_\mathrm{L}$（右旋和左旋圆偏振分析器）

其中圆偏振分析器由四分之一波片加线偏振器组成。

**矩阵表示**：

测量可表示为：
$$ I = (1/2)\mathbf{a}^\mathrm{T}\mathbf{S} $$

其中 $\mathbf{a}$ 是分析器矢量。例如：
- 水平偏振器：$\mathbf{a}_\mathrm{H} = [1, 1, 0, 0]^\mathrm{T}$
- 右旋圆偏振分析器：$\mathbf{a}_\mathrm{R} = [1, 0, 0, 1]^\mathrm{T}$

**Mueller 矩阵**：

斯托克斯矢量通过光学元件的变换用 $4\times4$ Mueller 矩阵描述：

$$ \mathbf{S}_\mathrm{out} = \mathbf{M} \mathbf{S}_\mathrm{in} $$

Mueller 矩阵包含了器件的全部偏振特性，包括偏振相关损耗和去偏振效应。

常见器件的 Mueller 矩阵：

线偏振器（水平）：
$$ \mathbf{M}_\mathrm{pol} = (1/2)\begin{pmatrix} 1 & 1 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix} $$

四分之一波片（快轴水平）：
$$ \mathbf{M}_\mathrm{QWP} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & -1 & 0 \end{pmatrix} $$

### 22.3.5 Mueller矩阵的深入分析

**Mueller矩阵的约束条件**：

物理可实现的Mueller矩阵必须满足：

1. **正定性约束**：
   对任意输入Stokes矢量 $\mathbf{S}_\mathrm{in}$（满足$S_0^2 \ge S_1^2 + S_2^2 + S_3^2$），
   输出 $\mathbf{S}_\mathrm{out} = \mathbf{M}\mathbf{S}_\mathrm{in}$ 也必须满足相同约束。

2. **特征值约束**：
   相干矩阵 $\mathbf{C} = \mathbf{M}^\mathrm{T}\mathbf{M}$ 的特征值必须非负。

3. **Cloude分解**：
   任意Mueller矩阵可分解为：
   $$ \mathbf{M} = \sum_i p_i\mathbf{M}_i $$
   其中 $p_i \ge 0$，$\sum p_i = 1$，$\mathbf{M}_i$ 是非去偏振Mueller矩阵。

**去偏振度的量化**：

Mueller矩阵的去偏振能力可用多个指标描述：

1. **去偏振指数**：
   $$ \mathrm{DI} = \sqrt{\sum_{i,j\ne00} M_{ij}^2}/\sqrt{3M_{00}^2} $$

2. **偏振纯度指数**：
   $P_1 = |\mathbf{m}_1|/m_{00}$（第一列）
   $P_2 = |\mathbf{m}^\mathrm{T}|/m_{00}$（第一行）
   
   其中 $\mathbf{m}_1 = [M_{10}, M_{20}, M_{30}]^\mathrm{T}$

3. **平均去偏振度**：
   $$ \Delta = 1 - (P_1 + P_2)/2 $$

**Mueller-Jones对应**：

对于非去偏振系统，Mueller矩阵与Jones矩阵的关系：

$$ \mathbf{M} = \mathbf{A}(\mathbf{J} \otimes \mathbf{J}^*)\mathbf{A}^{-1} $$

其中 $\otimes$ 是Kronecker积，$\mathbf{A}$ 是变换矩阵：

$$ \mathbf{A} = \begin{pmatrix} 1 & 0 & 0 & 1 \\ 1 & 0 & 0 & -1 \\ 0 & 1 & 1 & 0 \\ 0 & \mathrm{i} & -\mathrm{i} & 0 \end{pmatrix} $$

逆变换：
$$ \mathbf{A}^{-1} = (1/2)\begin{pmatrix} 1 & 1 & 0 & 0 \\ 0 & 0 & 1 & -\mathrm{i} \\ 0 & 0 & 1 & \mathrm{i} \\ 1 & -1 & 0 & 0 \end{pmatrix} $$

### 22.3.6 偏振态的统计描述

**相干矩阵表示**：

对于部分偏振光，引入$2\times2$相干矩阵：

$$ \mathbf{\Phi} = \langle\mathbf{E}\mathbf{E}^\dagger\rangle = \begin{pmatrix} \langle|Ex|^2\rangle & \langle Ex^*Ey \rangle \\ \langle ExEy^* \rangle & \langle|Ey|^2\rangle \end{pmatrix} $$

与Stokes参数的关系：
$$ \mathbf{\Phi} = (1/2)\begin{pmatrix} S_0 + S_1 & S_2 - \mathrm{i}S_3 \\ S_2 + \mathrm{i}S_3 & S_0 - S_1 \end{pmatrix} $$

相干矩阵性质：
- Hermitian：$\mathbf{\Phi} = \mathbf{\Phi}^\dagger$
- 半正定：特征值 $\lambda_1, \lambda_2 \ge 0$
- 迹：$\mathrm{Tr}(\mathbf{\Phi}) = S_0$（总强度）

**偏振熵**：

类比信息论，定义偏振熵：

$$ H = -\sum_i (\lambda_i/S_0)\log_2(\lambda_i/S_0) $$

其中 $\lambda_i$ 是相干矩阵的特征值。

- $H = 0$：完全偏振光（纯态）
- $H = 1$：完全非偏振光（最大混合态）

**偏振态的概率解释**：

部分偏振光可视为不同偏振态的统计混合：

$$ \mathbf{\Phi} = \sum_i p_i|\mathbf{E}_i\rangle\langle\mathbf{E}_i| $$

其中 $p_i$ 是第i个偏振态的概率权重。

这与量子力学的密度矩阵形式完全类似。

## 22.4 庞加莱球表示

### 22.4.1 球面几何与偏振态

庞加莱球（Poincaré sphere）提供了偏振态的优雅几何表示，将所有可能的偏振态映射到单位球面上。这种表示不仅直观，而且揭示了偏振态之间的深层几何关系。

**球面坐标与斯托克斯参数**：

对于完全偏振光，归一化斯托克斯参数满足：
$$ s_1^2 + s_2^2 + s_3^2 = 1 $$

其中 $s_i = S_i/S_0$（$i = 1,2,3$）。这定义了三维空间中的单位球面。

球面上每点的笛卡尔坐标直接对应归一化斯托克斯参数：
- $x = s_1 = \cos 2\chi \cos 2\psi$
- $y = s_2 = \cos 2\chi \sin 2\psi$  
- $z = s_3 = \sin 2\chi$

其中 $\psi$ 是偏振椭圆的方位角，$\chi$ 是椭圆率角。

**球面坐标表示**：

使用球坐标 $(\theta, \varphi)$：
- 极角 $\theta \in [0, \pi]$：$\theta = \pi/2 - 2\chi$
- 方位角 $\varphi \in [0, 2\pi]$：$\varphi = 2\psi$

斯托克斯参数的球坐标形式：
$s_1 = \sin \theta \cos \varphi$
$s_2 = \sin \theta \sin \varphi$
$s_3 = \cos \theta$

**特殊偏振态的位置**：

- 赤道（$\theta = \pi/2$）：所有线偏振态
  - H态：$(1,0,0)$，$\varphi = 0$
  - V态：$(-1,0,0)$，$\varphi = \pi$
  - D态：$(0,1,0)$，$\varphi = \pi/2$
  - A态：$(0,-1,0)$，$\varphi = 3\pi/2$

- 极点：圆偏振态
  - 北极 L态：$(0,0,1)$，左旋圆偏振
  - 南极 R态：$(0,0,-1)$，右旋圆偏振

- 一般位置：椭圆偏振态
  - 北半球：左旋椭圆偏振
  - 南半球：右旋椭圆偏振

### 22.4.2 偏振态的球面映射

**从琼斯矢量到庞加莱球**：

给定归一化琼斯矢量 $\mathbf{J} = \begin{pmatrix} a \\ b \end{pmatrix}$，其中 $|a|^2 + |b|^2 = 1$，对应的球面坐标为：

$s_1 = |a|^2 - |b|^2 = 2|a||b|\cos(\arg(b) - \arg(a))$
$s_2 = 2\mathrm{Re}(a^*b) = 2|a||b|\cos(\arg(b/a))$
$s_3 = 2\mathrm{Im}(a^*b) = 2|a||b|\sin(\arg(b/a))$

**立体投影表示**：

庞加莱球可通过立体投影映射到复平面。从南极向北投影：

$$ \zeta = \frac{s_1 + \mathrm{i}s_2}{1 + s_3} = \tan(\theta/2)\exp(\mathrm{i}\varphi) $$

这个复数 $\zeta$ 完全确定了偏振态。对于琼斯矢量 $\begin{pmatrix} a \\ b \end{pmatrix}$，有：
$$ \zeta = b/a $$

这建立了琼斯矢量复数比与庞加莱球的直接联系。

**对径点的物理意义**：

球面上的对径点表示正交偏振态：
- 若 $\mathbf{S}$ 和 $\mathbf{S}'$ 是对径点，则对应的偏振态正交
- $\mathbf{S}' = -\mathbf{S}$
- 正交性：$\langle\mathbf{J}|\mathbf{J}'\rangle = 0$

这个性质在偏振分析和量子信息中非常重要。

**大圆的意义**：

通过球心的平面与球面的交线形成大圆。任意大圆上的点构成一组特殊的偏振态家族：
- 赤道大圆：所有线偏振态
- 经线大圆：固定椭圆率，变化方位角
- 一般大圆：表示可由两个正交态线性叠加得到的所有态

### 22.4.3 偏振态演化的几何描述

**相位延迟器的作用**：

相位延迟器在庞加莱球上的作用是绕某轴的旋转。设延迟器快轴方向为 $\alpha$，相位延迟 $\delta$：

旋转轴：$\mathbf{n} = [\cos 2\alpha, \sin 2\alpha, 0]^\mathrm{T}$
旋转角：$\delta$

旋转矩阵（Rodrigues公式）：
$$ \mathbf{R}(\delta,\mathbf{n}) = \mathbf{I} + \sin \delta [\mathbf{n}]_\times + (1-\cos \delta)[\mathbf{n}]_\times^2 $$

其中 $[\mathbf{n}]_\times$ 是反对称矩阵：
$$ [\mathbf{n}]_\times = \begin{pmatrix} 0 & 0 & -\sin 2\alpha \\ 0 & 0 & \cos 2\alpha \\ \sin 2\alpha & -\cos 2\alpha & 0 \end{pmatrix} $$

**特殊情况**：

1. 半波片（$\delta = \pi$）：绕快轴方向旋转 180°
   - 保持与快轴平行的线偏振态不变
   - 反转垂直于快轴的偏振分量

2. 四分之一波片（$\delta = \pi/2$）：绕快轴方向旋转 90°
   - 将与快轴成 45° 的线偏振转为圆偏振
   - 将圆偏振转为线偏振

**偏振态序列的几何轨迹**：

连续变化的偏振器件产生的偏振态在庞加莱球上描绘连续轨迹：

1. 旋转线偏振器：轨迹为赤道上的点
2. 可变延迟器：轨迹为通过固定轴的大圆弧
3. 旋光器：轨迹为绕 z 轴的纬线圆

**偏振模色散（PMD）的几何描述**：

在光纤中，偏振态演化可表示为庞加莱球上的随机游走：
$$ \mathrm{d}\mathbf{s}/\mathrm{d}z = \mathbf{\Omega}(z) \times \mathbf{s} $$

其中 $\mathbf{\Omega}(z)$ 是局部双折射矢量。这导致：
- 短距离：确定性旋转
- 长距离：扩散行为
- 相关长度：定义偏振保持能力

### 22.4.4 与量子比特的类比

庞加莱球与量子力学中的布洛赫球（Bloch sphere）在数学上同构，这揭示了经典偏振光学与量子信息的深刻联系。

**量子态的对应**：

量子比特态：$|\psi\rangle = \alpha|0\rangle + \beta|1\rangle$，其中 $|\alpha|^2 + |\beta|^2 = 1$

对应的布洛赫矢量：
$$ \mathbf{r} = [2\mathrm{Re}(\alpha^*\beta), 2\mathrm{Im}(\alpha^*\beta), |\alpha|^2 - |\beta|^2]^\mathrm{T} $$

与偏振态的对应：
- $|0\rangle \leftrightarrow$ 水平偏振 $|\mathrm{H}\rangle$
- $|1\rangle \leftrightarrow$ 垂直偏振 $|\mathrm{V}\rangle$
- $(|0\rangle + |1\rangle)/\sqrt{2} \leftrightarrow$ 45°偏振 $|\mathrm{D}\rangle$
- $(|0\rangle + \mathrm{i}|1\rangle)/\sqrt{2} \leftrightarrow$ 左旋圆偏振 $|\mathrm{L}\rangle$

**泡利矩阵与斯托克斯参数**：

密度矩阵表示：
$$ \rho = (1/2)(\sigma_0 + s_1\sigma_1 + s_2\sigma_2 + s_3\sigma_3) $$

其中 $\sigma_i$ 是泡利矩阵：
$$ \sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \quad \sigma_2 = \begin{pmatrix} 0 & -\mathrm{i} \\ \mathrm{i} & 0 \end{pmatrix} \quad \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix} $$

这建立了斯托克斯参数与量子力学密度矩阵的直接联系。

**幺正变换与偏振变换**：

量子门操作 $U$ 对应偏振变换：
$|\psi'\rangle = U|\psi\rangle \leftrightarrow \mathbf{J}' = \mathbf{U}\mathbf{J}$

在布洛赫球上，这表现为旋转：
$$ \mathbf{r}' = \mathbf{R}\mathbf{r} $$

其中旋转矩阵 $\mathbf{R}$ 与幺正矩阵 $U$ 通过如下关系联系：
$$ U = \exp(-\mathrm{i}\theta\mathbf{n}\cdot\mathbf{\sigma}/2) \leftrightarrow \mathbf{R} = \exp(-\theta[\mathbf{n}]_\times) $$

**应用实例**：

1. **量子密钥分发（QKD）**：
   - BB84协议使用两组正交偏振基
   - 在庞加莱球上对应两组对径点

2. **量子态层析**：
   - 通过偏振测量重构量子态
   - 对应庞加莱球上的态重构

3. **几何相位（Berry相位）**：
   - 偏振态在参数空间闭合路径的相位
   - 等于庞加莱球上对应路径所围立体角的一半

## 22.5 偏振器件的矩阵描述

### 22.5.1 线偏振器

线偏振器是最基本的偏振器件，只允许特定方向的线偏振分量通过。理想线偏振器可用琼斯矩阵和Mueller矩阵完整描述。

**理想线偏振器的琼斯矩阵**：

透振方向与x轴夹角为$\theta$的线偏振器：

$$ \mathbf{P}(\theta) = \begin{pmatrix} \cos^2\theta & \cos \theta \sin \theta \\ \cos \theta \theta \sin \theta & \sin^2\theta \end{pmatrix} $$

这可以分解为：
$$ \mathbf{P}(\theta) = \mathbf{\hat{p}}\mathbf{\hat{p}}^\dagger $$

其中 $\mathbf{\hat{p}} = [\cos \theta, \sin \theta]^\mathrm{T}$ 是透振方向的单位矢量。

**特殊角度的偏振器**：

水平偏振器（$\theta = 0$）：
$$ \mathbf{P}_\mathrm{H} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} $$

垂直偏振器（$\theta = \pi/2$）：
$$ \mathbf{P}_\mathrm{V} = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix} $$

45°偏振器（$\theta = \pi/4$）：
$$ \mathbf{P}_\mathrm{D} = (1/2)\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} $$

**消光比和非理想性**：

实际偏振器存在有限消光比 $\varepsilon$：

$$ \mathbf{P}_\mathrm{real}(\theta) = \begin{pmatrix} \cos^2\theta + \varepsilon \sin^2\theta & (1-\varepsilon)\cos \theta \sin \theta \\ (1-\varepsilon)\cos \theta \sin \theta & \sin^2\theta + \varepsilon \cos^2\theta \end{pmatrix} $$

消光比定义：$\mathrm{ER} = 10 \log_{10}(1/\varepsilon) \text{ dB}$

典型值：
- 薄膜偏振器：20-30 dB
- 晶体偏振器：40-60 dB
- 纳米线栅偏振器：>40 dB

**Mueller矩阵表示**：

理想线偏振器的Mueller矩阵：

$$ \mathbf{M}_\mathrm{pol}(\theta) = (1/2)\begin{pmatrix} 1 + \cos 2\theta & \sin 2\theta & 0 & 0 \\ \cos 2\theta & \cos^2 2\theta & \sin 2\theta \cos 2\theta & 0 \\ \sin 2\theta & \sin 2\theta \cos 2\theta & \sin^2 2\theta & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix} $$

注意最后一行全零，表示完全去偏振化。

**偏振器的物理实现**：

1. **二向色性材料**（如Polaroid）：
   - 基于各向异性吸收
   - 宽带工作，但损耗较大
   
2. **双折射晶体**（如方解石）：
   - 利用o光和e光分离
   - 高消光比，窄带工作
   
3. **布儒斯特角反射**：
   - p偏振无反射
   - 消光比依赖入射角精度

4. **金属线栅**：
   - 亚波长周期结构
   - 宽带、高功率承受能力

### 22.5.2 相位延迟器

相位延迟器（波片）在两个正交偏振分量间引入相位差，是偏振控制的核心器件。

**一般相位延迟器的琼斯矩阵**：

快轴沿x轴，相位延迟$\delta$：

$$ \mathbf{R}_0(\delta) = \begin{pmatrix} 1 & 0 \\ 0 & \exp(\mathrm{i}\delta) \end{pmatrix} $$

快轴与x轴夹角$\theta$的一般延迟器：

$$ \mathbf{R}(\theta,\delta) = \mathbf{T}(-\theta)\mathbf{R}_0(\delta)\mathbf{T}(\theta) $$

展开得：
$$ \mathbf{R}(\theta,\delta) = \begin{pmatrix} \cos^2\theta + \sin^2\theta \mathrm{e}^{\mathrm{i}\delta} & \sin \theta \cos \theta(1 - \mathrm{e}^{\mathrm{i}\delta}) \\ \sin \theta \cos \theta(1 - \mathrm{e}^{\mathrm{i}\delta}) & \sin^2\theta + \cos^2\theta \mathrm{e}^{\mathrm{i}\delta} \end{pmatrix} $$

其中旋转矩阵：
$$ \mathbf{T}(\theta) = \begin{pmatrix} \cos \theta & -\sin \theta \\ \sin \theta & \cos \theta \end{pmatrix} $$

**常用波片**：

1. **四分之一波片（QWP）**，$\delta = \pi/2$：
   
   $$ \mathbf{R}_\mathrm{QWP}(\theta) = \begin{pmatrix} \cos^2\theta + \mathrm{i} \sin^2\theta & (1-\mathrm{i})\sin \theta \cos \theta \\ (1-\mathrm{i})\sin \theta \cos \theta & \sin^2\theta + \mathrm{i} \cos^2\theta \end{pmatrix} $$
   
   特殊取向：
   - $\theta = 0$：$\mathbf{R}_\mathrm{QWP} = \mathrm{diag}(1, \mathrm{i})$
   - $\theta = 45^\circ$：线偏振$\leftrightarrow$圆偏振转换器

2. **半波片（HWP）**，$\delta = \pi$：
   
   $$ \mathbf{R}_\mathrm{HWP}(\theta) = \begin{pmatrix} \cos 2\theta & \sin 2\theta \\ \sin 2\theta & -\cos 2\theta \end{pmatrix} $$
   
   作用：将偏振方向旋转$2\theta$

3. **全波片（FWP）**，$\delta = 2\pi$：
   
   $\mathbf{R}_\mathrm{FWP} = \mathbf{I}$（恒等变换）
   
   用于补偿和相位调节

**色散与温度效应**：

实际波片的相位延迟依赖波长和温度：

$$ \delta(\lambda,T) = 2\pi(n_\mathrm{e} - n_\mathrm{o})d/\lambda + \alpha(T - T_0) $$

其中：
- $(n_\mathrm{e} - n_\mathrm{o})$：双折射率差
- $d$：晶体厚度
- $\alpha$：温度系数

**零级、多级和复合波片**：

1. **零级波片**：$\delta = \delta_0$（单片，很薄）
   - 优点：色散小
   - 缺点：机械脆弱

2. **多级波片**：$\delta = 2\pi m + \delta_0$
   - 优点：机械强度好
   - 缺点：色散大，温度敏感

3. **复合波片**：两片相减
   $\delta = \delta_1 - \delta_2 = \delta_0$
   - 平衡了色散和强度

**Mueller矩阵表示**：

相位延迟器的Mueller矩阵：

$$ \mathbf{M}_\mathrm{ret}(\theta,\delta) = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos^22\theta+\sin^22\theta \cos \delta & \sin 2\theta \cos 2\theta(1-\cos \delta) & -\sin 2\theta \sin \delta \\ 0 & \sin 2\theta \cos 2\theta(1-\cos \delta) & \sin^22\theta+\cos^22\theta \cos \delta & \cos 2\theta \sin \delta \\ 0 & \sin 2\theta \sin \delta & -\cos 2\theta \sin \delta & \cos \delta \end{pmatrix} $$

### 22.5.3 旋光器

旋光器旋转线偏振光的偏振方向，而不改变偏振态的其他性质。这种效应源于圆双折射。

**旋光的物理机制**：

1. **自然旋光性**（石英、糖溶液）：
   - 分子手性导致
   - 旋转角：$\alpha = \pi\Delta n \cdot d/\lambda$
   - $\Delta n = n_\mathrm{L} - n_\mathrm{R}$（左右圆偏振折射率差）

2. **法拉第旋光**（磁光效应）：
   - 外加磁场诱导
   - 旋转角：$\alpha = V \cdot B \cdot d$
   - $V$：Verdet常数

**琼斯矩阵表示**：

旋光角为$\alpha$的旋光器：

$$ \mathbf{O}(\alpha) = \begin{pmatrix} \cos \alpha & -\sin \alpha \\ \sin \alpha & \cos \alpha \end{pmatrix} $$

这与坐标旋转矩阵相同，但物理意义不同：
- 坐标旋转：观察者视角改变
- 旋光：偏振态实际旋转

**旋光器的重要性质**：

1. **非互易性**（法拉第旋光）：
   - 正向传播：旋转$+\alpha$
   - 反向传播：再旋转$+\alpha$（而非$-\alpha$）
   - 总旋转：$2\alpha$

2. **互易性**（自然旋光）：
   - 正向：$+\alpha$
   - 反向：$-\alpha$
   - 往返后无净旋转

**旋光器与波片的组合**：

旋光器可由圆偏振基下的相位延迟实现：

$$ \mathbf{O}(\alpha) = \mathbf{T}_\mathrm{cir} \mathbf{R}_0(2\alpha) \mathbf{T}_\mathrm{cir}^\dagger $$

其中 $\mathbf{T}_\mathrm{cir}$ 是从线偏振基到圆偏振基的变换：

$$ \mathbf{T}_\mathrm{cir} = (1/\sqrt{2})\begin{pmatrix} 1 & 1 \\ \mathrm{i} & -\mathrm{i} \end{pmatrix} $$

**Mueller矩阵**：

$$ \mathbf{M}_\mathrm{rot}(\alpha) = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos 2\alpha & -\sin 2\alpha & 0 \\ 0 & \sin 2\alpha & \cos 2\alpha & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} $$

在庞加莱球上表现为绕$S_3$轴（z轴）旋转$2\alpha$。

### 22.5.4 复合器件系统

实际偏振系统由多个器件级联构成。理解其综合效应需要矩阵运算和系统分析。

**级联系统的矩阵描述**：

N个器件级联：
- 琼斯矩阵：$\mathbf{J}_\mathrm{total} = \mathbf{J}_N\cdots\mathbf{J}_2\mathbf{J}_1$
- Mueller矩阵：$\mathbf{M}_\mathrm{total} = \mathbf{M}_N\cdots\mathbf{M}_2\mathbf{M}_1$

注意顺序：光先通过器件1，最后通过器件N。

**偏振补偿器设计**：

**Babinet-Soleil补偿器**：
- 两个楔形波片，相对滑动
- 连续可调相位延迟：$\delta = k(d_1 - d_2)$
- 矩阵：$\mathbf{R}(\theta, \delta(d))$

**Berek补偿器**：
- 倾斜波片
- 延迟量：$\delta = \delta_0/\cos \varphi$（$\varphi$为倾斜角）
- 同时改变延迟量和快轴方向

**偏振控制器**：

标准三片式控制器（QWP-HWP-QWP）：

$$ \mathbf{P}_\mathrm{total} = \mathbf{R}_\mathrm{QWP}(\theta_3)\mathbf{R}_\mathrm{HWP}(\theta_2)\mathbf{R}_\mathrm{QWP}(\theta_1) $$

可实现任意偏振态变换（SU(2)群的任意元素）。

参数与庞加莱球旋转的关系：
- $\theta_1$：绕$S_1$轴旋转
- $\theta_2$：绕变换后的$S_3$轴旋转
- $\theta_3$：绕最终$S_1$轴旋转

**偏振分析器设计**：

**完整Stokes偏振计**：

需要至少4次测量。典型配置：
1. $I_{0^\circ}$：水平偏振器
2. $I_{90^\circ}$：垂直偏振器
3. $I_{45^\circ}$：45°偏振器
4. $I_\mathrm{RCP}$：QWP + 水平偏振器

Stokes参数：
$S_0 = I_{0^\circ} + I_{90^\circ}$
$S_1 = I_{0^\circ} - I_{90^\circ}$
$S_2 = 2I_{45^\circ} - S_0$
$S_3 = 2I_\mathrm{RCP} - S_0$

**旋转波片偏振计**：

连续旋转QWP，固定偏振器：
$$ I(\omega_1 t) = (1/2)[S_0 + S_1\cos(4\omega_1 t) + S_2\sin(4\omega_1 t)\cos \delta - S_3\sin(4\omega_1 t)\sin \delta] $$

傅里叶分析提取Stokes参数。

**矩阵条件数与误差传播**：

器件组合的数值稳定性：
$$ \kappa(\mathbf{M}) = \|\mathbf{M}\| \cdot \|\mathbf{M}^{-1}\| $$

条件数大表示对输入扰动敏感。优化设计应最小化条件数。

误差传播：
$$ \delta\mathbf{S}_\mathrm{out} \approx \mathbf{M} \delta\mathbf{S}_\mathrm{in} + \delta\mathbf{M} \mathbf{S}_\mathrm{in} $$

需要考虑：
- 器件制造公差
- 对准误差
- 波长依赖性
- 温度漂移

### 22.5.5 实用偏振测量与应用

偏振测量技术在科学研究和工业应用中扮演着重要角色。本节探讨实际偏振测量系统的设计原理及其在各领域的应用。

**高精度偏振测量**：

**双旋转延迟器偏振计**：

最精确的偏振测量方法之一，使用两个同步旋转的延迟器：

探测强度：
$$ I(t) = a_0 + \sum_{n=1}^4 [a_n\cos(n\omega_1 t) + b_n\sin(n\omega_1 t)] + \sum_{m=1}^4 [c_m\cos(m\omega_2 t) + d_m\sin(m\omega_2 t)] $$

其中系数与Stokes参数的关系通过傅里叶分析确定。

优势：
- 自校准能力
- 消除系统误差
- 宽光谱范围
- 精度可达 $\Delta P/P \sim 10^{-5}$

**分振幅偏振计**：

使用偏振分束器同时测量两个正交分量：

设计考虑：
1. 分束器消光比 > $10^3:1$
2. 探测器匹配：响应度差异 < 0.1%
3. 同步采集消除时间变化

校准矩阵：
$$ \mathbf{A} = \begin{pmatrix} a_{11} & a_{12} & a_{13} & a_{14} \\ a_{21} & a_{22} & a_{23} & a_{24} \\ a_{31} & a_{32} & a_{33} & a_{34} \\ a_{41} & a_{42} & a_{43} & a_{44} \end{pmatrix} $$

Stokes参数反演：
$$ \mathbf{S} = \mathbf{A}^{-1}\mathbf{I} $$

其中 $\mathbf{I} = [I_1, I_2, I_3, I_4]^\mathrm{T}$ 是四个通道的强度。

**光谱偏振测量**：

**通道光谱偏振计（Channeled Spectropolarimeter）**：

利用厚延迟器产生的光谱调制：

$$ I(\sigma) = (1/2)[S_0(\sigma) + S_1(\sigma)\cos(2\pi\delta_1\sigma) + S_2(\sigma)\sin(2\pi\delta_1\sigma)\cos \theta_1 - S_3(\sigma)\sin(2\pi\delta_1\sigma)\sin \theta_1] $$

其中 $\sigma = 1/\lambda$ 是波数，$\delta_1$ 是延迟器厚度。

通过傅里叶变换恢复Stokes光谱：
$$ S_0(\sigma) = 2\int I(\sigma')\mathrm{d}\sigma' $$
$$ S_1(\sigma) = 4\int I(\sigma')\cos(2\pi\delta_1\sigma')\mathrm{d}\sigma' $$

分辨率与延迟器厚度的关系：
$$ \Delta\sigma = 1/(2\delta_1) $$

**成像偏振测量**：

**分焦平面偏振相机**：

微偏振器阵列集成在探测器上：
- $2\times2$超像素：0°, 45°, 90°, 135°
- 瞬时Stokes参数计算
- 空间分辨率降低4倍

插值算法恢复全分辨率：
$S_0(x,y) = I_{0^\circ} + I_{90^\circ} + 2I_\mathrm{interp}$
$S_1(x,y) = I_{0^\circ} - I_{90^\circ}$
$S_2(x,y) = I_{45^\circ} - I_{135^\circ}$

其中 $I_\mathrm{interp}$ 是邻近像素插值。

**液晶可调谐偏振成像**：

使用液晶可变延迟器（LCVR）：
$$ \delta(V) = 2\pi\Delta n(V)d/\lambda $$

时序控制获取完整Stokes图像：
1. $V_1 \to \delta = 0$（参考）
2. $V_2 \to \delta = \pi/2$（QWP）
3. $V_3 \to \delta = \pi$（HWP）
4. 不同快轴角度

优势：无机械运动，快速切换（~10ms）。

**椭偏测量技术**：

**光谱椭偏仪原理**：

测量反射光偏振态变化确定材料光学常数：

反射系数比：
$$ \rho = r_\mathrm{p}/r_\mathrm{s} = \tan \psi \exp(\mathrm{i}\Delta) $$

其中 $\psi$ 和 $\Delta$ 是椭偏参数。

与材料参数的关系（各向同性材料）：
$$ \rho = \frac{\sqrt{n_2^2 - n_1^2\sin^2\theta_1} - n_1\cos \theta_1}{\sqrt{n_2^2 - n_1^2\sin^2\theta_1} + n_1\cos \theta_1} \times \frac{n_2^2\cos \theta_1 - n_1\sqrt{n_2^2 - n_1^2\sin^2\theta_1}}{n_2^2\cos \theta_1 + n_1\sqrt{n_2^2 - n_1^2\sin^2\theta_1}} $$

多层结构需要传输矩阵方法：
$$ \mathbf{M} = \prod_i \mathbf{M}_i $$

模型拟合提取：
- 层厚度
- 折射率 $n(\lambda)$
- 消光系数 $k(\lambda)$
- 表面粗糙度

**穆勒矩阵椭偏测量**：

完整16元素Mueller矩阵测量：

归一化Mueller矩阵：
$$ \mathbf{m} = \mathbf{M}/M_{00} $$

去偏振度分析：
- 偏振纯度：$P_\mathrm{u} = |\mathbf{m}_1|$（第一列）
- 双衰减：$D = 1 - \min(\text{eigenvalues})$

各向异性参数提取：
- 线性双折射：$\mathrm{LB} = |m_{12} - m_{21}|$
- 圆双折射：$\mathrm{CB} = |m_{13} - m_{31}|$
- 线性二向色性：$\mathrm{LD}' = |m_{01} - m_{10}|$

**偏振在各领域的应用**：

**遥感与大气光学**：

**气溶胶特性反演**：

Mie散射的偏振特征：
$$ P(\theta) = \frac{|S_2(\theta)|^2 - |S_1(\theta)|^2}{|S_2(\theta)|^2 + |S_1(\theta)|^2} $$

其中 $S_1, S_2$ 是散射振幅函数。

多角度偏振测量反演：
- 粒径分布 $n(r)$
- 复折射率 $m = n + \mathrm{i}k$
- 形状因子

反演算法：
$$ \min \sum_i [P_\mathrm{obs}(\theta_i) - P_\mathrm{model}(\theta_i,\mathbf{p})]^2/\sigma_i^2 $$

其中 $\mathbf{p}$ 是待反演参数。

**偏振雷达**：

降水粒子识别：
- 差分反射率：$\mathrm{ZDR} = 10\log(Z_\mathrm{hh}/Z_\mathrm{vv})$
- 相关系数：$\rho_\mathrm{hv} = |\langle S_\mathrm{hh}^*S_\mathrm{vv}\rangle|/\sqrt{\langle|S_\mathrm{hh}|^2\rangle\langle|S_\mathrm{vv}|^2\rangle}$
- 比差分相位：$\mathrm{KDP} = (\Phi_\mathrm{DP}(r_2) - \Phi_\mathrm{DP}(r_1))/(2(r_2-r_1))$

水凝物分类：
- 雨滴：$\mathrm{ZDR} > 0.5 \text{ dB}$, $\rho_\mathrm{hv} > 0.97$
- 冰雹：$\mathrm{ZDR} \sim 0 \text{ dB}$, $\rho_\mathrm{hv} < 0.95$
- 融化层：ZDR峰值，$\rho_\mathrm{hv}$下降

**生物医学成像**：

**偏振敏感OCT（PS-OCT）**：

组织双折射测量：
$$ \delta(z) = 2k_0\int_0^z \Delta n(z')\mathrm{d}z' $$

其中 $\Delta n$ 是局部双折射。

从Jones矩阵提取：
$$ \mathbf{J}(z) = \begin{pmatrix} \sqrt{R_{11}} & \sqrt{R_{12}} \exp(\mathrm{i}\varphi_{12}) \\ \sqrt{R_{21}} \exp(\mathrm{i}\varphi_{21}) & \sqrt{R_{22}} \end{pmatrix} $$

相位延迟：
$$ \delta = \arg(J_{11}J_{22}) - \arg(J_{12}J_{21}) $$

快轴方向：
$$ \theta = (1/2)\arctan\left[\frac{2\mathrm{Re}(J_{12}J_{21})}{|J_{11}|^2 - |J_{22}|^2}\right] $$

临床应用：
- 视网膜神经纤维层（双折射 $\sim0.5^\circ/\mu\mathrm{m}$）
- 动脉粥样硬化斑块表征
- 烧伤深度评估

**穆勒矩阵显微镜**：

组织病理诊断参数：
1. 线性延迟：$\delta = \arccos[(m_{22} + m_{33})/2]$
2. 圆二向色性：$\psi = \arctan(m_{34}/m_{24})$
3. 去偏振度：$\Delta = 1 - |\det(\mathbf{m}_{3\times3})|^{1/3}$

其中 $\mathbf{m}_{3\times3}$ 是$3\times3$子矩阵。

癌变组织特征：
- 去偏振度增加（散射增强）
- 双折射模式改变（胶原重组）
- 各向异性增强

**量子通信与密码学**：

**BB84协议实现**：

四态编码：
$|0\rangle_\mathrm{H} = |\mathrm{H}\rangle$, $|1\rangle_\mathrm{H} = |\mathrm{V}\rangle$（直线基）
$|0\rangle_\mathrm{D} = |\mathrm{D}\rangle$, $|1\rangle_\mathrm{D} = |\mathrm{A}\rangle$（对角基）

安全性分析：
- 量子比特错误率（QBER）：$e = N_\mathrm{error}/N_\mathrm{total}$
- 安全阈值：$e < 11\%$（理想情况）
- 密钥率：$R = 1 - h(e) - h(e)$（信息论界限）

其中 $h(x) = -x \log_2 x - (1-x)\log_2(1-x)$ 是二元熵函数。

**连续变量QKD**：

使用偏振调制的相干态：
$|\alpha\rangle_\mathrm{H} + |\beta\rangle_\mathrm{V}$

协方差矩阵：
$$ \mathbf{V} = \begin{pmatrix} \langle\Delta X^2\rangle & \langle\Delta X\Delta P\rangle \\ \langle\Delta P\Delta X\rangle & \langle\Delta P^2\rangle \end{pmatrix} $$

安全条件（反向协调）：
$I(A:B) > \chi(B:E)$

其中 $I$ 是互信息，$\chi$ 是Holevo信息。

**工业检测应用**：

**光弹性应力分析**：

应力诱导双折射：
$$ \Delta n = C(\sigma_1 - \sigma_2) $$

其中 $C$ 是光弹性系数，$\sigma_i$ 是主应力。

等色线（恒定相位延迟）：
$$ N = \delta/(2\pi) = (t/\lambda)C(\sigma_1 - \sigma_2) $$

等倾线（恒定主应力方向）：
通过圆偏振光照明消除

全场应力分离：
- 六步相移法
- RGB光弹性法
- 数字图像相关

**表面检测**：

偏振散射测量缺陷：

BRDF的偏振分量：
$$ f_{pq} = \mathrm{d}L_p/(E_q \mathrm{d}\Omega \cos \theta_i) $$

其中 $p,q \in \{\mathrm{s},\mathrm{p}\}$ 表示出射和入射偏振。

表面粗糙度参数提取：
- RMS高度：$\sigma_\mathrm{h}$
- 相关长度：$l_\mathrm{c}$
- 功率谱密度：$\mathrm{PSD}(k)$

与散射的关系（Beckmann模型）：
$$ f \propto \exp[-(4\pi\sigma_\mathrm{h} \cos \theta_i/\lambda)^2] $$

缺陷分类：
- 划痕：强各向异性散射
- 颗粒污染：去偏振增强
- 膜层缺陷：偏振对比度异常

**先进偏振技术展望**：

**超表面偏振器件**：

亚波长结构实现任意偏振变换：
- 几何相位（Pancharatnam-Berry相位）
- 传播相位调控
- 局部偏振转换

设计原理：
$$ \mathbf{J}_\mathrm{out} = \mathbf{T}(x,y)\mathbf{J}_\mathrm{in} $$

其中 $\mathbf{T}(x,y)$ 是空间变化的Jones矩阵。

应用：
- 偏振全息
- 矢量光束生成
- 偏振复用通信

**量子偏振纠缠**：

双光子偏振纠缠态：
$$ |\Psi\rangle = (1/\sqrt{2})(|\mathrm{H}\rangle_\mathrm{A}|\mathrm{V}\rangle_\mathrm{B} - |\mathrm{V}\rangle_\mathrm{A}|\mathrm{H}\rangle_\mathrm{B}) $$

Bell不等式测试：
$$ S = |E(a,b) - E(a,b') + E(a',b) + E(a',b')| \le 2 \text{（经典）} $$
$S_\mathrm{quantum} = 2\sqrt{2}$（量子力学预测）

其中 $E(a,b)$ 是关联函数。

应用前景：
- 量子计算
- 量子成像
- 量子传感

## 本章小结

本章系统介绍了偏振光学的数学基础，建立了从基本概念到实际应用的完整理论体系：

1. **偏振态描述**：深入探讨了电磁波的偏振现象，建立了线偏振、圆偏振和椭圆偏振的数学描述。通过偏振椭圆的几何构造，揭示了不同偏振态之间的内在联系。

2. **琼斯矢量方法**：介绍了完全偏振光的复振幅表示，建立了琼斯矩阵演算体系。通过矩阵运算描述偏振器件的作用，为偏振系统设计提供了强大工具。重要公式包括：
   - 琼斯矢量：$\mathbf{J} = [Ex_0, Ey_0]^\mathrm{T}$
   - 器件作用：$\mathbf{J}_\mathrm{out} = \mathbf{M}\mathbf{J}_\mathrm{in}$
   - 级联系统：$\mathbf{M}_\mathrm{total} = \mathbf{M}_n\cdots\mathbf{M}_2\mathbf{M}_1$

3. **斯托克斯参数**：扩展到部分偏振光的描述，定义了四个实参数 $S_0, S_1, S_2, S_3$。建立了偏振度概念 $P = \sqrt{S_1^2 + S_2^2 + S_3^2}/S_0$，连接了偏振与相干性。Mueller矩阵方法允许描述去偏振效应。

4. **庞加莱球表示**：将所有偏振态映射到单位球面，提供了偏振态演化的直观几何图像。揭示了与量子比特的深刻类比，建立了经典偏振光学与量子信息的桥梁。关键关系：
   - 球面坐标：$(s_1, s_2, s_3) = (\cos 2\chi \cos 2\psi, \cos 2\chi \sin 2\psi, \sin 2\chi)$
   - 偏振变换：球面旋转
   - 正交态：对径点

5. **偏振器件矩阵**：详细推导了线偏振器、相位延迟器、旋光器的琼斯矩阵和Mueller矩阵。探讨了实际器件的非理想性，如有限消光比、色散效应等。复合系统设计原则强调了矩阵运算顺序的重要性。

6. **实用测量与应用**：从高精度偏振测量技术到各领域应用，展示了偏振光学的广泛影响。涵盖遥感、生物医学成像、量子通信、工业检测等前沿应用，突出了偏振作为信息载体的独特优势。

本章建立的数学框架为后续章节的偏振渲染奠定了坚实基础，特别是Mueller矩阵方法将直接应用于描述材料的偏振散射特性。

## 练习题

### 基础题

**练习22.1** 证明任意偏振态可以分解为两个圆偏振态的叠加。
*提示*：使用圆偏振基 $\{|\mathrm{R}\rangle, |\mathrm{L}\rangle\}$，其中 $|\mathrm{R}\rangle = (1/\sqrt{2})\begin{pmatrix} 1 \\ -\mathrm{i} \end{pmatrix}$，$|\mathrm{L}\rangle = (1/\sqrt{2})\begin{pmatrix} 1 \\ \mathrm{i} \end{pmatrix}$。

<details>
<summary>答案</summary>

任意归一化琼斯矢量 $\mathbf{J} = \begin{pmatrix} a \\ b \end{pmatrix}$，其中 $|a|^2 + |b|^2 = 1$。

定义圆偏振基：
$|\mathrm{R}\rangle = (1/\sqrt{2})\begin{pmatrix} 1 \\ -\mathrm{i} \end{pmatrix}$，$|\mathrm{L}\rangle = (1/\sqrt{2})\begin{pmatrix} 1 \\ \mathrm{i} \end{pmatrix}$

验证正交性：$\langle\mathrm{R}|\mathrm{L}\rangle = (1/2)(1 - \mathrm{i}^2) = 1 \checkmark$

展开系数：
$\alpha_\mathrm{R} = \langle\mathrm{R}|\mathbf{J}\rangle = (1/\sqrt{2})(a^* + \mathrm{i}b^*)$
$\alpha_\mathrm{L} = \langle\mathrm{L}|\mathbf{J}\rangle = (1/\sqrt{2})(a^* - \mathrm{i}b^*)$

验证：
$$ \alpha_\mathrm{R}|\mathrm{R}\rangle + \alpha_\mathrm{L}|\mathrm{L}\rangle = (1/2)\left[(a^* + \mathrm{i}b^*)\begin{pmatrix} 1 \\ -\mathrm{i} \end{pmatrix} + (a^* - \mathrm{i}b^*)\begin{pmatrix} 1 \\ \mathrm{i} \end{pmatrix}\right] $$
$$ = (1/2)\begin{pmatrix} 2a^* \\ -\mathrm{i}a^* - \mathrm{i}^2b^* + \mathrm{i}a^* - \mathrm{i}^2b^* \end{pmatrix} = (1/2)\begin{pmatrix} 2a^* \\ -\mathrm{i}a^* + b^* + \mathrm{i}a^* + b^* \end{pmatrix} = \begin{pmatrix} a \\ b \end{pmatrix} = \mathbf{J} \checkmark $$

物理意义：任意偏振光可视为左旋和右旋圆偏振光的相干叠加。
</details>

**练习22.2** 一束线偏振光通过两个偏振器，第一个偏振器透振方向与入射光偏振方向夹角30°，第二个与第一个夹角60°。计算最终透射光强度（初始光强为$I_0$）。
*提示*：使用Malus定律或琼斯矩阵方法。

<details>
<summary>答案</summary>

方法1：Malus定律
第一个偏振器后：$I_1 = I_0\cos^2(30^\circ) = I_0(\sqrt{3}/2)^2 = 3I_0/4$
第二个偏振器后：$I_2 = I_1\cos^2(60^\circ) = (3I_0/4)(1/2)^2 = 3I_0/16$

方法2：琼斯矩阵
设入射光沿x方向偏振：$\mathbf{J}_\mathrm{in} = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$

第一个偏振器（30°）：
$$ \mathbf{P}_1 = \begin{pmatrix} \cos^230^\circ & \cos30^\circ\sin30^\circ \\ \cos30^\circ\sin30^\circ & \sin^230^\circ \end{pmatrix} = \begin{pmatrix} 3/4 & \sqrt{3}/4 \\ \sqrt{3}/4 & 1/4 \end{pmatrix} $$

第二个偏振器（与第一个夹角60°，即与x轴夹角90°）：
$$ \mathbf{P}_2 = \begin{pmatrix} \cos^290^\circ & \cos90^\circ\sin90^\circ \\ \cos90^\circ\sin90^\circ & \sin^290^\circ \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix} $$

总矩阵：$\mathbf{M} = \mathbf{P}_2\mathbf{P}_1 = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 3/4 & \sqrt{3}/4 \\ \sqrt{3}/4 & 1/4 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ \sqrt{3}/4 & 1/4 \end{pmatrix}$

输出：$\mathbf{J}_\mathrm{out} = \mathbf{M}\begin{pmatrix} 1 \\ 0 \end{pmatrix} = \begin{pmatrix} 0 \\ \sqrt{3}/4 \end{pmatrix}$

强度：$I = |J_x|^2 + |J_y|^2 = 0^2 + (\sqrt{3}/4)^2 = 3/16 I_0 \checkmark$
</details>

**练习22.3** 推导斯托克斯参数与偏振椭圆参数（方位角$\psi$，椭圆率角$\chi$）之间的关系。
*提示*：从电场分量的参数表示开始。

<details>
<summary>答案</summary>

设归一化电场：
$Ex = \cos \theta$
$Ey = \sin \theta \exp(\mathrm{i}\delta)$

其中 $\theta$ 决定振幅比，$\delta$ 是相位差。

偏振椭圆参数与$(\theta,\delta)$的关系：
$\tan 2\psi = \tan 2\theta \cos \delta$
$\sin 2\chi = \sin 2\theta \sin \delta$

计算斯托克斯参数：
$S_0 = |Ex|^2 + |Ey|^2 = \cos^2\theta + \sin^2\theta = 1$（归一化）
$S_1 = |Ex|^2 - |Ey|^2 = \cos^2\theta - \sin^2\theta = \cos 2\theta$
$S_2 = 2\mathrm{Re}(Ex^*Ey) = 2\cos \theta \sin \theta \cos \delta = \sin 2\theta \cos \delta$
$S_3 = 2\mathrm{Im}(Ex^*Ey) = 2\cos \theta \sin \theta \sin \delta = \sin 2\theta \sin \delta$

从上述关系解出：
$\cos 2\theta = S_1$
$\sin 2\theta \cos \delta = S_2$
$\sin 2\theta \sin \delta = S_3$

因此：
$\tan 2\psi = S_2/S_1$
$\sin 2\chi = S_3/\sqrt{S_1^2 + S_2^2 + S_3^2} = S_3$（完全偏振光）

反解：
$S_1 = \cos 2\chi \cos 2\psi$
$S_2 = \cos 2\chi \sin 2\psi$
$S_3 = \sin 2\chi$
</details>

### 挑战题

**练习22.4** 设计一个偏振系统，将任意输入偏振态转换为预定的输出偏振态。给出所需器件的最少数量和参数。
*提示*：考虑SU(2)群的生成元数量。

<details>
<summary>答案</summary>

理论基础：任意偏振变换属于SU(2)群，需要3个独立参数。

最小系统：两个可旋转的延迟器
1. 第一个延迟器：相位延迟$\delta_1$，快轴角度$\theta_1$
2. 第二个延迟器：相位延迟$\delta_2$，快轴角度$\theta_2$

总共4个参数，但由于整体相位任意，实际3个自由度。

通用解：QWP-HWP-QWP组合
$$ \mathbf{M} = \mathbf{R}_\mathrm{QWP}(\theta_3)\mathbf{R}_\mathrm{HWP}(\theta_2)\mathbf{R}_\mathrm{QWP}(\theta_1) $$

参数确定步骤：
1. 将输入态$\mathbf{J}_\mathrm{in}$转到庞加莱球坐标$(\theta_\mathrm{in}, \varphi_\mathrm{in})$
2. 将输出态$\mathbf{J}_\mathrm{out}$转到球坐标$(\theta_\mathrm{out}, \varphi_\mathrm{out})$
3. 计算所需旋转：
   - 先绕x轴转$\theta_1$角
   - 再绕新z轴转$\theta_2$角
   - 最后绕新x轴转$\theta_3$角

数学表达：
$\theta_1 = \arctan[\sin(\theta_\mathrm{out} - \theta_\mathrm{in})/(1 - \cos(\theta_\mathrm{out} - \theta_\mathrm{in})\cos(\varphi_\mathrm{out} - \varphi_\mathrm{in}))]$
$\theta_2 = \varphi_\mathrm{out} - \varphi_\mathrm{in}$
$\theta_3 = \arctan[\sin(\theta_\mathrm{out} - \theta_\mathrm{in})/(\cos(\theta_\mathrm{out} - \theta_\mathrm{in}) - \cos(\varphi_\mathrm{out} - \varphi_\mathrm{in}))]$

特殊情况优化：
- 线偏振$\to$线偏振：单个HWP
- 线偏振$\to$圆偏振：单个QWP
- 相位调节：单个延迟器
</details>

**练习22.5** 证明部分偏振光的偏振度P与相干矩阵的特征值$\lambda_1$、$\lambda_2$的关系为：$P = |\lambda_1 - \lambda_2|/(\lambda_1 + \lambda_2)$。
*提示*：使用相干矩阵的迹和行列式。

<details>
<summary>答案</summary>

相干矩阵：
$$ \mathbf{\Phi} = \begin{pmatrix} \langle|Ex|^2\rangle & \langle Ex^*Ey \rangle \\ \langle ExEy^* \rangle & \langle|Ey|^2\rangle \end{pmatrix} $$

与斯托克斯参数关系：
$$ \mathbf{\Phi} = (1/2)\begin{pmatrix} S_0 + S_1 & S_2 - \mathrm{i}S_3 \\ S_2 + \mathrm{i}S_3 & S_0 - S_1 \end{pmatrix} $$

计算不变量：
$\mathrm{Tr}(\mathbf{\Phi}) = \lambda_1 + \lambda_2 = S_0$
$\det(\mathbf{\Phi}) = \lambda_1\lambda_2 = (1/4)[S_0^2 - S_1^2 - S_2^2 - S_3^2]$

特征方程：
$$ \lambda^2 - S_0\lambda + (1/4)[S_0^2 - S_1^2 - S_2^2 - S_3^2] = 0 $$

解得：
$$ \lambda_{1,2} = (S_0/2) \pm (1/2)\sqrt{S_1^2 + S_2^2 + S_3^2} $$

因此：
$\lambda_1 - \lambda_2 = \sqrt{S_1^2 + S_2^2 + S_3^2}$
$\lambda_1 + \lambda_2 = S_0$

偏振度：
$$ P = \frac{\sqrt{S_1^2 + S_2^2 + S_3^2}}{S_0} = \frac{|\lambda_1 - \lambda_2|}{\lambda_1 + \lambda_2} \checkmark $$

物理意义：
- $\lambda_1 = \lambda_2$：完全非偏振（$P = 0$）
- $\lambda_2 = 0$：完全偏振（$P = 1$）
- 一般情况：$P$反映两特征值的相对差异
</details>

**练习22.6** 推导通过散射介质后偏振度退化的表达式。假设介质引入随机相位延迟$\delta$，其概率分布为$p(\delta)$。
*提示*：计算期望值$\langle\exp(\mathrm{i}\delta)\rangle$。

<details>
<summary>答案</summary>

初始偏振态：$\mathbf{J}_0 = \begin{pmatrix} a \\ b \end{pmatrix}$

通过随机延迟器后：
$$ \mathbf{J}(\delta) = \begin{pmatrix} a \\ b \exp(\mathrm{i}\delta) \end{pmatrix} $$

计算相干矩阵元素：
$\Phi_{11} = \langle|a|^2\rangle = |a|^2$（不变）
$\Phi_{22} = \langle|b|^2\rangle = |b|^2$（不变）
$\Phi_{12} = \langle ab^* \exp(\mathrm{i}\delta)\rangle = ab^*\langle\exp(\mathrm{i}\delta)\rangle$
$\Phi_{21} = \langle a^*b \exp(-\mathrm{i}\delta)\rangle = a^*b\langle\exp(-\mathrm{i}\delta)\rangle$

关键：计算特征函数
$$ \chi(1) = \langle\exp(\mathrm{i}\delta)\rangle = \int p(\delta)\exp(\mathrm{i}\delta)\mathrm{d}\delta $$

常见分布：
1. 均匀分布 $p(\delta) = 1/(2\pi)$，$\delta \in [0, 2\pi]$：
   $\chi(1) = 0 \to$ 完全去偏振

2. 高斯分布 $p(\delta) = (1/\sqrt{2\pi\sigma^2})\exp(-\delta^2/2\sigma^2)$：
   $\chi(1) = \exp(-\sigma^2/2)$

3. Von Mises分布 $p(\delta) = \exp(\kappa\cos \delta)/2\pi I_0(\kappa)$：
   $\chi(1) = I_1(\kappa)/I_0(\kappa)$

初始和最终斯托克斯参数：
初始：$S_1^{(0)} = |a|^2 - |b|^2$，$S_2^{(0)} = 2\mathrm{Re}(ab^*)$，$S_3^{(0)} = 2\mathrm{Im}(ab^*)$

最终：
$S_1 = S_1^{(0)}$（不变）
$S_2 = S_2^{(0)}\mathrm{Re}[\chi(1)]$
$S_3 = S_3^{(0)}\mathrm{Re}[\chi(1)]$

偏振度变化：
$$ P_\mathrm{final}/P_\mathrm{initial} = \frac{\sqrt{S_1^2 + (S_2^{(0)})^2|\chi(1)|^2 + (S_3^{(0)})^2|\chi(1)|^2}}{\sqrt{S_1^2 + (S_2^{(0)})^2 + (S_3^{(0)})^2}} $$

特殊情况：
- 圆偏振入射（$S_1 = 0$）：$P_\mathrm{final} = P_\mathrm{initial}|\chi(1)|$
- 线偏振入射（$S_3 = 0$）：部分保持
</details>

**练习22.7** 分析Pancharatnam-Berry几何相位。当偏振态在庞加莱球上沿闭合路径C演化时，计算获得的几何相位。
*提示*：几何相位等于路径所围立体角的一半。

<details>
<summary>答案</summary>

偏振态参数化：
$|\psi(t)\rangle = \cos(\theta/2)|\mathrm{H}\rangle + \sin(\theta/2)\exp(\mathrm{i}\varphi)|\mathrm{V}\rangle$

对应庞加莱球坐标：$(\theta, \varphi)$

沿闭路C演化后的几何相位：
$$ \gamma = (1/2)\oint_C \mathbf{A}\cdot\mathrm{d}\mathbf{l} $$

其中$\mathbf{A}$是"矢势"：
$A_\theta = 0$
$A_\varphi = \cos \theta - 1$

因此：
$$ \gamma = (1/2)\oint_C (\cos \theta - 1)\mathrm{d}\varphi = -(1/2)\Omega $$

其中$\Omega$是路径所围立体角。

具体例子：

1. 赤道大圆（所有线偏振态）：
   路径：$\theta = \pi/2$，$\varphi \in [0, 2\pi]$
   立体角：$\Omega = 2\pi$
   几何相位：$\gamma = -\pi$

2. 纬线圆（固定椭圆率）：
   路径：$\theta = \theta_0$，$\varphi \in [0, 2\pi]$
   立体角：$\Omega = 2\pi(1 - \cos \theta_0)$
   几何相位：$\gamma = -\pi(1 - \cos \theta_0)$

3. 测地三角形：
   顶点：北极、$(\theta=\pi/2,\varphi=0)$、$(\theta=\pi/2,\varphi=\pi/2)$
   立体角：$\Omega = \pi/2$
   几何相位：$\gamma = -\pi/4$

实验验证：
使用三个波片序列实现闭合路径
相位测量：干涉法或偏振分析

应用：
- 几何相位调制器
- 鲁棒量子计算
- 拓扑光子学器件
</details>

**练习22.8** 设计一个Mueller矩阵椭偏测量系统，要求能够完整测量16个矩阵元素。分析测量精度与系统参数的关系。
*提示*：需要至少16个独立测量，考虑条件数优化。

<details>
<summary>答案</summary>

测量方程：
$$ I_k = \mathbf{a}_k^\mathrm{T} \mathbf{M} \mathbf{s}_k $$

其中$\mathbf{s}_k$是入射Stokes矢量，$\mathbf{a}_k$是分析器矢量。

最优设计原则：
1. 输入态集合$\{\mathbf{s}_k\}$应均匀分布在庞加莱球上
2. 分析器集合$\{\mathbf{a}_k\}$也应均匀分布
3. 最小化条件数$\kappa(\mathbf{W})$，其中$\mathbf{W}$是测量矩阵

具体配置：

输入态（4个）：
- H：$[1, 1, 0, 0]^\mathrm{T}$
- V：$[1, -1, 0, 0]^\mathrm{T}$  
- D：$[1, 0, 1, 0]^\mathrm{T}$
- R：$[1, 0, 0, 1]^\mathrm{T}$

分析器（4个）：
- 同上配置

测量矩阵（$16\times16$）：
$W_{ij} = (a_i \otimes s_j)$

Mueller矩阵反演：
$$ \mathrm{vec}(\mathbf{M}) = \mathbf{W}^{-1} \mathbf{I} $$

其中$\mathbf{I} = [I_1, I_2, \dots, I_{16}]^\mathrm{T}$ 是四个通道的强度。

误差分析：

1. 强度噪声：$\sigma_I$
   Mueller元素误差：$\sigma_M = \kappa(\mathbf{W})\sigma_I/\sqrt{N}$

2. 器件误差：
   - 延迟器误差$\delta\varphi \to \Delta M \propto \sin(2\theta)\delta\varphi$
   - 方位角误差$\delta\theta \to \Delta M \propto \cos(2\theta)\delta\theta$

3. 系统校准：
   使用已知样品（空气、标准延迟器）
   最小二乘拟合修正系统误差

优化方案：

1. 双旋转延迟器：
   - 连续测量，过定系统
   - 傅里叶分析提取$M_{ij}$
   - 自校准能力

2. 分振幅配置：
   - 4个通道同时测量
   - 减少测量时间
   - 需要精确通道校准

3. 光谱域测量：
   - 利用色散同时获得多个延迟
   - 单次测量获得完整矩阵
   - 光谱分辨率限制

精度极限：
- 散粒噪声：$\Delta M/M \sim 1/\sqrt{N_\mathrm{photons}}$
- 系统稳定性：$\sim0.1\%$（温控）
- 最终精度：$\sim10^{-3} - 10^{-4}$
</details>

## 常见陷阱与错误

### 概念混淆

1. **偏振与偏振化**
   - 错误：认为非偏振光通过偏振器后"变成"偏振光
   - 正确：偏振器选择性透过某一偏振分量，不改变光的偏振性质

2. **相位差的符号约定**
   - 错误：混淆$\delta = \varphi_y - \varphi_x$还是$\varphi_x - \varphi_y$
   - 正确：保持一致的约定，注意与旋向定义的关系

3. **Jones矢量的归一化**
   - 错误：忘记Jones矢量可以包含强度信息
   - 正确：明确是否使用归一化形式，注意物理含义

### 计算错误

4. **矩阵乘法顺序**
   - 错误：$\mathbf{M}_1\mathbf{M}_2$与$\mathbf{M}_2\mathbf{M}_1$混淆
   - 正确：光先通过的器件矩阵在右边

5. **坐标系旋转**
   - 错误：混淆器件旋转与坐标系旋转
   - 正确：$\mathbf{M}' = \mathbf{R}(-\theta)\mathbf{M}\mathbf{R}(\theta)$为器件旋转

6. **庞加莱球角度因子**
   - 错误：忘记庞加莱球上使用$2\psi$和$2\chi$
   - 正确：球面角度是偏振椭圆角度的两倍

### 实验误差

7. **消光比假设**
   - 错误：假设偏振器是理想的
   - 正确：考虑有限消光比对测量的影响

8. **波长依赖性**
   - 错误：忽略波片的色散
   - 正确：校准应在工作波长进行

9. **温度效应**
   - 错误：忽略温度对双折射的影响
   - 正确：监控温度或使用温度补偿设计

### 物理理解

10. **部分偏振光的处理**
    - 错误：试图用Jones矢量描述部分偏振光
    - 正确：使用Stokes参数或相干矩阵

11. **去偏振机制**
    - 错误：认为去偏振是偏振态的破坏
    - 正确：理解为偏振态的统计平均

12. **几何相位的理解**
    - 错误：认为几何相位是动力学效应
    - 正确：纯粹的几何/拓扑效应

## 最佳实践检查清单

### 理论分析
- [ ] 明确偏振态描述方法的适用范围
- [ ] 保持符号约定的一致性
- [ ] 验证能量守恒和幺正性条件
- [ ] 检查极限情况的物理合理性

### 系统设计
- [ ] 选择合适的偏振元件组合
- [ ] 考虑器件的非理想特性
- [ ] 优化系统的条件数
- [ ] 包含校准和对准程序

### 实验测量
- [ ] 使用适当的偏振态采样
- [ ] 实施系统误差校正
- [ ] 监控环境参数（温度、振动）
- [ ] 验证测量的自洽性

### 数据处理
- [ ] 正确实施矩阵反演算法
- [ ] 评估测量不确定度
- [ ] 检查物理约束条件
- [ ] 使用适当的拟合和优化方法

### 应用考虑
- [ ] 匹配理论模型与实际系统
- [ ] 考虑带宽和时间响应
- [ ] 评估环境稳定性要求
- [ ] 制定维护和校准计划

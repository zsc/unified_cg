# 第21章：半经典光-物质相互作用

本章介绍半经典光学理论，其中电磁场用经典方式描述，而物质系统用量子力学处理。这种方法成功解释了大多数激光-原子相互作用现象，包括吸收、受激发射、拉比振荡和饱和效应。我们将建立从经典渲染到量子光学的桥梁，为理解现代光学现象提供必要的理论基础。

## 学习目标

完成本章后，您将能够：
1. 推导并求解光学Bloch方程
2. 分析二能级系统在激光场中的动力学
3. 计算拉比频率和饱和强度
4. 理解功率展宽和光谱烧孔现象
5. 识别半经典理论的适用范围和局限性
6. 将半经典概念与体积渲染框架联系起来

## 21.1 经典场与量子化物质基础

### 21.1.1 半经典近似的物理动机

在许多光学现象中，光场强度足够大，使得光子数不确定度相对较小。此时，可以将光场视为经典电磁波，而保持原子系统的量子描述。这种半经典近似在以下条件下有效：

1. **大光子数条件**：$\langle\hat{n}\rangle \gg 1$，其中 $\hat{n}$ 是光子数算符
2. **相干态近似**：激光场接近相干态 $|\alpha\rangle$
3. **忽略量子涨落**：$\Delta n/\langle n\rangle \ll 1$

定量分析：对于典型激光器（功率 P = 1 mW，波长 $\lambda = 632.8 \text{ nm}$），单模光子数：

$\langle n\rangle = P/(\hbar\omega\gamma) \approx 10^{12} \gg 1$

其中 $\gamma$ 是腔线宽。相对涨落：

$\Delta n/\langle n\rangle = 1/\sqrt{\langle n\rangle} \approx 10^{-6}$

因此半经典近似极其精确。

**相干态的经典对应**：
相干态 $|\alpha\rangle$ 是湮灭算符的本征态：$\hat{a}|\alpha\rangle = \alpha|\alpha\rangle$，其中 $\alpha = |\alpha|e^{i\varphi}$ 是复振幅。电场期望值：

$\langle\alpha|\hat{E}(\mathbf{r},t)|\alpha\rangle = E_0\cos(\mathbf{k}\cdot\mathbf{r} - \omega t + \varphi)$

其中 $E_0 = 2|\alpha|\sqrt{\hbar\omega/2\varepsilon_0V}$，V 是量子化体积。这正是经典电磁场的形式，证明了半经典处理的合理性。

**半经典近似的适用边界**：
1. **单光子物理失效**：当 $\langle n\rangle \sim 1$ 时，必须考虑光子统计
2. **强耦合腔QED**：当 $g\sqrt{n} > \kappa,\gamma$（g是耦合强度，$\kappa$是腔损耗）
3. **非线性量子光学**：光子-光子相互作用变得重要

**从路径积分视角理解**：
光场的路径积分在大光子数极限下由鞍点主导：

$Z = \int \mathcal{D}A \exp[iS[A]/\hbar] \approx \exp[iS[A_{cl}]/\hbar]$

其中 $A_{cl}$ 是经典场构型。量子修正项 $\sim O(\hbar/\langle n\rangle E)$，对于大光子数可忽略。

**半经典近似的层次结构**：
不同物理过程需要不同层次的量子描述：

1. **零阶半经典**：经典场 + 经典物质（几何光学）
2. **一阶半经典**：经典场 + 量子物质（本章内容）
3. **二阶半经典**：量子化场的经典极限 + 量子物质
4. **全量子**：量子场 + 量子物质（第27-28章）

**与其他近似方法的关系**：
- **WKB近似**：处理空间缓变场，$\mathbf{k}\cdot\delta\mathbf{r} \gg 1$
- **绝热近似**：时间缓变场，$|\dot{\omega}/\omega^2| \ll 1$
- **Born-Oppenheimer近似**：分离电子与核运动
- **平均场近似**：多体系统的有效单粒子描述

**实验判据**：
判断是否需要超越半经典的实验特征：

1. **光子反聚束**：$g^{(2)}(0) < 1$，单光子源特征
2. **压缩光**：$\Delta X_1\Delta X_2 < 1/4$，低于散粒噪声极限
3. **纠缠**：违反Bell不等式
4. **Wigner函数负值**：相空间的量子特征

**计算复杂度考虑**：
- 半经典方法：$O(N)$，N是原子数
- 全量子方法：$O(N^{2m})$，m是光子数截断
- 这解释了为什么半经典方法在实际计算中如此重要

### 21.1.2 经典电磁场描述

考虑单色平面波：

$\mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0 \cos(\mathbf{k}\cdot\mathbf{r} - \omega t + \varphi) = \frac{1}{2}\mathbf{E}_0[e^{i(\mathbf{k}\cdot\mathbf{r} - \omega t + \varphi)} + c.c.]$

其中：
- $\mathbf{E}_0$ 是电场振幅矢量
- $\mathbf{k} = 2\pi/\lambda \hat{\mathbf{n}}$ 是波矢（$\hat{\mathbf{n}}$ 是传播方向）
- $\omega = 2\pi c/\lambda$ 是角频率
- $\varphi$ 是初相位

在偶极近似下（原子尺度 $a_0 \ll$ 波长 $\lambda$），可以忽略空间变化：

$\mathbf{E}(t) = \mathbf{E}_0 \cos(\omega t - \varphi)$

偶极近似的有效性：对于可见光 $\lambda \approx 500 \text{ nm}$ 和 Bohr 半径 $a_0 \approx 0.05 \text{ nm}$：

$\mathbf{k}\cdot a_0 \approx 2\pi(a_0/\lambda) \approx 6 \times 10^{-4} \ll 1$

因此 $e^{i\mathbf{k}\cdot\mathbf{r}} \approx 1 + i\mathbf{k}\cdot\mathbf{r} \approx 1$ 在原子尺度上。

**电磁场的完整数学描述**：
从Maxwell方程组出发，在无源空间中：

$\nabla \times \mathbf{E} = -\partial\mathbf{B}/\partial t$
$\nabla \times \mathbf{B} = \mu_0\varepsilon_0\partial\mathbf{E}/\partial t$
$\nabla \cdot \mathbf{E} = 0$
$\nabla \cdot \mathbf{B} = 0$

波动方程：$\nabla^2\mathbf{E} - (1/c^2)\partial^2\mathbf{E}/\partial t^2 = 0$

平面波解的一般形式：
$\mathbf{E}(\mathbf{r},t) = \text{Re}[\mathbf{E}_0 \exp(i\mathbf{k}\cdot\mathbf{r} - i\omega t)]$

其中色散关系：$\omega^2 = c^2k^2$

**偏振态的完整描述**：
任意偏振态可分解为两个正交分量：

$\mathbf{E} = E_x \hat{\mathbf{e}}_x + E_y \hat{\mathbf{e}}_y$

其中复振幅：
$E_x = |E_x|\exp(i\varphi_x)$
$E_y = |E_y|\exp(i\varphi_y)$

偏振椭圆的参数：
- 椭圆度：$\tan(\chi) = \pm|E_y|/|E_x|$（当$\Delta\varphi = \pm\pi/2$）
- 方位角：$\tan(2\psi) = 2|E_x||E_y|\cos(\Delta\varphi)/(|E_x|^2 - |E_y|^2)$
- 相位差：$\Delta\varphi = \varphi_y - \varphi_x$

特殊情况：
1. 线偏振：$\Delta\varphi = 0$ 或 $\pi$
2. 圆偏振：$|E_x| = |E_y|$，$\Delta\varphi = \pm\pi/2$
3. 椭圆偏振：一般情况

**超越偶极近似**：
对于某些情况需要考虑高阶项：

$\mathbf{E}(\mathbf{r}, t) = \mathbf{E}(\mathbf{r}_0, t)[1 + i\mathbf{k}\cdot(\mathbf{r}-\mathbf{r}_0) - \frac{1}{2}(\mathbf{k}\cdot(\mathbf{r}-\mathbf{r}_0))^2 + ...]$

1. **电四极跃迁**：当偶极矩禁戒时，保留到 $\mathbf{k}\cdot\mathbf{r}$ 项
2. **磁偶极跃迁**：来自磁场 $\mathbf{B} = \mathbf{k} \times \mathbf{E}/\omega$ 的贡献
3. **多光子过程**：需要保留空间梯度以满足动量守恒

**光场的规范选择**：
在Coulomb规范下（$\nabla\cdot\mathbf{A} = 0$），电磁场可表示为：

$\mathbf{E} = -\partial\mathbf{A}/\partial t$，$\mathbf{B} = \nabla \times \mathbf{A}$

矢势的量子化形式：

$\mathbf{A}(\mathbf{r},t) = \sum_{\mathbf{k},\lambda} \sqrt{\hbar/2\varepsilon_0\omega_{\mathbf{k}}V} \hat{\mathbf{e}}_{\mathbf{k},\lambda} [\hat{a}_{\mathbf{k},\lambda} e^{i(\mathbf{k}\cdot\mathbf{r}-\omega_{\mathbf{k}}t)} + h.c.]$

其中 $\hat{\mathbf{e}}_{\mathbf{k},\lambda}$ 是偏振矢量，满足 $\mathbf{k}\cdot\hat{\mathbf{e}}_{\mathbf{k},\lambda} = 0$。

**脉冲光场的描述**：
对于有限持续时间的激光脉冲：

$\mathbf{E}(\mathbf{r},t) = \mathbf{E}_0\varepsilon(t)\cos[\mathbf{k}\cdot\mathbf{r} - \int_0^t \omega(t')dt' + \varphi(t)]$

其中：
- $\varepsilon(t)$ 是包络函数（如高斯型：$\varepsilon(t) = \exp(-t^2/2\tau^2)$）
- $\omega(t)$ 允许频率啁啾
- $\varphi(t)$ 包含相位调制

**光强与Poynting矢量**：
瞬时光强：$I(t) = \varepsilon_0c|\mathbf{E}(t)|^2$
时间平均光强：$\langle I\rangle = \varepsilon_0c|\mathbf{E}_0|^2/2$
能流密度：$\mathbf{S} = \mathbf{E} \times \mathbf{H} = \varepsilon_0c^2\mathbf{E} \times \mathbf{B}$

**复杂光场的模式分解**：
任意光场可分解为模式叠加：

$\mathbf{E}(\mathbf{r},t) = \sum_n c_n\mathbf{E}_n(\mathbf{r})e^{-i\omega_n t}$

常见模式基：
1. **平面波基**：$\exp(i\mathbf{k}\cdot\mathbf{r})$，适合自由空间
2. **球面波基**：$j_l(kr)Y_{lm}(\theta,\varphi)$，适合散射问题
3. **Hermite-Gaussian模式**：$\text{HG}_{nm}$，适合激光腔
4. **Laguerre-Gaussian模式**：$\text{LG}_{pl}$，携带轨道角动量

**结构光场**：
现代光学中的结构光场具有复杂的空间分布：

1. **涡旋光束**：$\mathbf{E} \propto \exp(il\varphi)$，轨道角动量 $l\hbar$
2. **贝塞尔光束**：$\mathbf{E} \propto J_n(k_\perp r)\exp(ik_z z)$，无衍射传播
3. **艾里光束**：$\mathbf{E} \propto \text{Ai}(x)$，自加速传播
4. **矢量光束**：空间变化的偏振态

**光场的相干性描述**：
部分相干光需要用相关函数描述：

$\Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau) = \langle E^*(\mathbf{r}_1,t)E(\mathbf{r}_2,t+\tau)\rangle$

相干度：$\gamma_{12}(\tau) = \Gamma(\mathbf{r}_1,\mathbf{r}_2,\tau)/\sqrt{[\Gamma(\mathbf{r}_1,\mathbf{r}_1,0)\Gamma(\mathbf{r}_2,\mathbf{r}_2,0)]}$

**光场的统计性质**：
不同光源的统计特性：

1. **相干态（激光）**：泊松统计，$\langle\Delta n^2\rangle = \langle n\rangle$
2. **热光**：超泊松统计，$\langle\Delta n^2\rangle = \langle n\rangle(1 + \langle n\rangle)$
3. **压缩光**：亚泊松统计，$\langle\Delta n^2\rangle < \langle n\rangle$
4. **单光子态**：Fock态，$\langle\Delta n^2\rangle = 0$

**非线性光学效应的唯象描述**：
强场下需考虑非线性极化：

$\mathbf{P} = \varepsilon_0[\chi^{(1)}\mathbf{E} + \chi^{(2)}\mathbf{E}\mathbf{E} + \chi^{(3)}\mathbf{E}\mathbf{E}\mathbf{E} + ...]$

导致：
- 二次谐波产生（SHG）
- 和频/差频产生（SFG/DFG）
- 四波混频（FWM）
- 光学Kerr效应

**光场的时空耦合**：
超短脉冲中时间和空间自由度耦合：

1. **脉冲前沿倾斜**：$\partial t_{pulse}/\partial x \neq 0$
2. **空间啁啾**：$\omega = \omega(x)$
3. **角色散**：$\mathbf{k} = \mathbf{k}(\omega)$
4. **时空涡旋**：螺旋型时空相位

### 21.1.3 原子的量子力学描述

考虑具有离散能级的原子系统。系统状态由波函数描述：

$|\psi(t)\rangle = \sum_n c_n(t)e^{-iE_n t/\hbar}|n\rangle$

其中 $|n\rangle$ 是能量本征态，$E_n$ 是对应能量。展开系数满足归一化：

$\sum_n |c_n(t)|^2 = 1$

密度算符形式（纯态）：

$\hat{\rho}(t) = |\psi(t)\rangle\langle\psi(t)|$

更一般的混合态：

$\hat{\rho}(t) = \sum_i p_i|\psi_i(t)\rangle\langle\psi_i(t)|$

其中 $p_i$ 是系综权重，$\sum_i p_i = 1$。

密度矩阵元素：

$\rho_{mn}(t) = \langle m|\hat{\rho}(t)|n\rangle = c_m(t)c_n^*(t)e^{-i(E_m-E_n)t/\hbar}$

物理意义：
- 对角元 $\rho_{nn}$：能级 $|n\rangle$ 的布居数概率
- 非对角元 $\rho_{mn} (m\neq n)$：量子相干性

**量子态的几何表示**：
在N维Hilbert空间中，纯态构成复射影空间$\text{CP}^{N-1}$。对于二能级系统：

$|\psi\rangle = \alpha|g\rangle + \beta|e\rangle$， $|\alpha|^2 + |\beta|^2 = 1$

可参数化为：
$\alpha = \cos(\theta/2)e^{i\varphi_1}$
$\beta = \sin(\theta/2)e^{i\varphi_2}$

整体相位不影响物理，取$\varphi_1 = 0$，得到Bloch球表示：
$|\psi\rangle = \cos(\theta/2)|g\rangle + e^{i\varphi}\sin(\theta/2)|e\rangle$

**密度矩阵的谱分解**：
任意密度矩阵可对角化：

$\hat{\rho} = \sum_i \lambda_i|\lambda_i\rangle\langle\lambda_i|$

其中$\lambda_i$是本征值（$0 \le \lambda_i \le 1$），$|\lambda_i\rangle$是对应本征态。

纯度参数：$P = \text{Tr}(\hat{\rho}^2) = \sum_i \lambda_i^2$
- $P = 1$：纯态（仅一个$\lambda_i = 1$）
- $P < 1$：混合态
- $P = 1/N$：最大混合态（所有$\lambda_i = 1/N$）

**量子态的纠缠测度**：
对于复合系统$\hat{\rho}_{AB}$，部分迹给出约化密度矩阵：

$\hat{\rho}_A = \text{Tr}_B(\hat{\rho}_{AB})$

纠缠熵：$S(\hat{\rho}_A) = -\text{Tr}(\hat{\rho}_A \ln \hat{\rho}_A)$
- $S = 0$：可分离态
- $S > 0$：纠缠态
- $S = \ln(d_A)$：最大纠缠（$d_A$是子系统A的维度）

**密度矩阵的性质与约束**：
1. **厄米性**：$\hat{\rho}^\dagger = \hat{\rho} \implies \rho_{mn} = \rho_{nm}^*$
2. **迹归一**：$\text{Tr}(\hat{\rho}) = \sum_n \rho_{nn} = 1$
3. **正定性**：$\langle\psi|\hat{\rho}|\psi\rangle \ge 0$ 对所有 $|\psi\rangle$
4. **纯度**：$\text{Tr}(\hat{\rho}^2) \le 1$，等号成立当且仅当纯态

**von Neumann熵**：
量子系统的熵定义为：

$S = -k_B \text{Tr}(\hat{\rho} \ln \hat{\rho}) = -k_B \sum_n \lambda_n \ln \lambda_n$

其中 $\lambda_n$ 是密度矩阵的本征值。
- 纯态：$S = 0$
- 最大混合态：$S = k_B \ln N$（N是维度）

**相干性的量化**：
$l_1$-范数相干性：$C_{l_1}(\rho) = \sum_{m\neq n} |\rho_{mn}|$
相对熵相干性：$C_r(\rho) = S(\rho_{diag}) - S(\rho)$

其中 $\rho_{diag}$ 是去除非对角元后的密度矩阵。

**开放系统动力学**：
考虑系统与环境耦合，总密度矩阵演化：

$\hat{\rho}_{total}(t) = \hat{U}(t)\hat{\rho}_{total}(0)\hat{U}^\dagger(t)$

对环境求迹得约化密度矩阵：

$\hat{\rho}_S(t) = \text{Tr}_E[\hat{\rho}_{total}(t)]$

这导致非幺正演化和退相干。

**主方程形式**：
在马尔可夫近似下，系统演化遵循Lindblad主方程：

$d\hat{\rho}/dt = -i[\hat{H},\hat{\rho}]/\hbar + \sum_k \gamma_k(\hat{L}_k\hat{\rho}\hat{L}_k^\dagger - \frac{1}{2}\{\hat{L}_k^\dagger\hat{L}_k,\hat{\rho}\})$

其中 $\hat{L}_k$ 是Lindblad算符，描述不同的耗散通道。

**原子能级结构的精细细节**：
实际原子能级包含多重结构：

1. **精细结构**：自旋-轨道耦合
   - 能级分裂：$\Delta E_{fs} \sim \alpha^2E_n$（$\alpha$是精细结构常数）
   - 导致J量子数：$\mathbf{J} = \mathbf{L} + \mathbf{S}$

2. **超精细结构**：核自旋耦合
   - 能级分裂：$\Delta E_{hfs} \sim (m_e/M_p)\Delta E_{fs}$
   - 导致F量子数：$\mathbf{F} = \mathbf{J} + \mathbf{I}$

3. **Zeeman效应**：外磁场下的分裂
   - 线性Zeeman：$\Delta E = \mu_B g_J m_J B$
   - 二次Zeeman：$\Delta E \propto B^2$（强场）

4. **Stark效应**：外电场下的分裂
   - 线性Stark：存在于有永久偶极矩的态
   - 二次Stark：$\Delta E \propto E^2$（普遍存在）

**多电子原子的处理**：
对于多电子系统，需要考虑：

1. **中心场近似**：
   $V_{eff}(\mathbf{r}) = V_{nuc}(\mathbf{r}) + V_{screen}(\mathbf{r})$
   
2. **组态相互作用**：
   $|\psi\rangle = \sum_i c_i|\Phi_i\rangle$（$\Phi_i$是Slater行列式）

3. **LS耦合与jj耦合**：
   - 轻原子：LS耦合，$\mathbf{L} = \sum_i\mathbf{l}_i$，$\mathbf{S} = \sum_i\mathbf{s}_i$
   - 重原子：jj耦合，$\mathbf{j}_i = \mathbf{l}_i + \mathbf{s}_i$

**量子缺陷理论**：
对于Rydberg态，能级可表示为：

$E_n = -R_\infty/(n - \delta_l)^2$

其中$\delta_l$是量子缺陷，反映核心电子的屏蔽效应。

**时间依赖的量子系统**：
对于含时哈密顿量$\hat{H}(t)$，演化算符：

$\hat{U}(t,t_0) = \mathcal{T} \exp[-i\int_{t_0}^t \hat{H}(t')dt'/\hbar]$

其中$\mathcal{T}$是时间排序算符。

**绝热定理与Berry相位**：
缓慢变化的哈密顿量下，系统保持在瞬时本征态，但获得几何相位：

$\gamma_n = i\oint\langle n(\mathbf{R})|\nabla_{\mathbf{R}}|n(\mathbf{R})\rangle\cdot d\mathbf{R}$

这是Berry相位，纯几何起源，与动力学相位不同。

**退相干机制**：
主要退相干源包括：

1. **自发辐射**：$T_1$过程，能量弛豫
2. **弹性碰撞**：纯失相，保持能量
3. **非弹性碰撞**：同时影响$T_1$和$T_2$
4. **磁场涨落**：导致相位扩散
5. **温度涨落**：热声子散射

### 21.1.4 相互作用哈密顿量

在电偶极近似下，光-物质相互作用哈密顿量为：

$\hat{H}_{int} = -\hat{\mathbf{d}}\cdot\mathbf{E}(t)$

其中 $\hat{\mathbf{d}} = -e\hat{\mathbf{r}}$ 是电偶极矩算符（$e > 0$ 为基本电荷）。

总哈密顿量：

$\hat{H} = \hat{H}_0 + \hat{H}_{int} = \hat{H}_0 - \hat{\mathbf{d}}\cdot\mathbf{E}(t)$

偶极矩矩阵元素：

$\mathbf{d}_{mn} = \langle m|\hat{\mathbf{d}}|n\rangle = -e\langle m|\hat{\mathbf{r}}|n\rangle$

选择定则：
- 宇称选择定则：$\Delta l = \pm 1$（电偶极跃迁）
- 角动量选择定则：$\Delta J = 0, \pm 1$（但 $J = 0 \leftrightarrow J = 0$ 禁戒）
- 自旋选择定则：$\Delta S = 0$（LS耦合下）

相互作用能量尺度：

$E_{int} \approx |\mathbf{d}||\mathbf{E}| \approx ea_0E_0$

对于 $E_0 = 10^6 \text{ V/m}$（典型激光场），$E_{int} \approx 10^{-23} \text{ J} \approx 10^{-4} \text{ eV}$。

**精确的选择定则推导**：
从宇称算符$\hat{P}$的性质出发：$\hat{P}|nlm\rangle = (-1)^l|nlm\rangle$

偶极算符具有奇宇称：$\hat{P}\hat{\mathbf{r}}\hat{P}^{-1} = -\hat{\mathbf{r}}$

因此跃迁矩阵元：
$\langle n'l'm'|\hat{\mathbf{r}}|nlm\rangle \propto \langle l'm'|(-1)^{l'+l+1}|lm\rangle$

非零条件：$(-1)^{l'+l+1} = 1 \implies l' + l = \text{奇数} \implies \Delta l = \pm 1, \pm 3, ...$

但球谐函数的正交性进一步限制：$\Delta l = \pm 1$

**相互作用能量的数量级估计**：
不同尺度的比较：
- 原子能级间隔：$\Delta E \sim \text{eV} - \text{keV}$
- 热能：$kT \sim 0.025 \text{ eV}$ (室温)
- 相互作用能：$E_{int} \sim 10^{-4} \text{ eV}$ ($1 \text{ MW/cm}^2$)
- 真空场涨落：$\delta E \sim \hbar\omega/V^{1/3} \sim 10^{-12} \text{ eV}$

强场判据：
- 微扰区：$E_{int} \ll \Delta E$
- 非微扰区：$E_{int} \sim \Delta E$
- 超强场：$E_{int} > \Delta E$（电离）

**偶极矩的对称性分析**：
利用Wigner-Eckart定理，偶极矩矩阵元素可写为：

$\langle n'l'm'|\hat{d}_q|nlm\rangle = \langle n'l'||\hat{d}||nl\rangle\langle l'm'|lm;1q\rangle$

其中：
- $\langle n'l'||\hat{d}||nl\rangle$ 是约化矩阵元素
- $\langle l'm'|lm;1q\rangle$ 是Clebsch-Gordan系数
- $q = 0, \pm 1$ 对应球坐标分量

**高阶多极矩展开**：
完整的相互作用哈密顿量包含：

$\hat{H}_{int} = -\hat{\mathbf{d}}\cdot\mathbf{E} - \hat{\mathbf{Q}}:\nabla\mathbf{E} - \hat{\mathbf{m}}\cdot\mathbf{B} + ...$

其中：
- $\hat{Q}_{ij} = e \sum_k \hat{r}_i^{(k)}\hat{r}_j^{(k)}$ 是电四极矩张量
- $\hat{\mathbf{m}} = (e/2m)\hat{\mathbf{L}} + g_s(e/2m)\hat{\mathbf{S}}$ 是磁偶极矩

**规范不变性**：
在不同规范下，相互作用形式不同：
1. **长度规范**：$\hat{H}_{int} = -\hat{\mathbf{d}}\cdot\mathbf{E}(\mathbf{r}_0,t)$
2. **速度规范**：$\hat{H}_{int} = -(e/m)\hat{\mathbf{p}}\cdot\mathbf{A}(\mathbf{r}_0,t)$
3. **加速度规范**：用于强场物理

规范变换：$\hat{U} = \exp[ie \mathbf{A}(\mathbf{r},t)\cdot\mathbf{r}/\hbar]$

**跃迁强度与振子强度**：
振子强度定义：

$f_{mn} = (2m/\hbar^2)\omega_{mn}|\langle m|\hat{\mathbf{r}}|n\rangle|^2$

满足Thomas-Reich-Kuhn求和规则：

$\sum_n f_{n0} = N$（电子数）

**相互作用图像**：
在相互作用图像中，态矢和算符演化分离：

$|\psi_I(t)\rangle = e^{i\hat{H}_0 t/\hbar}|\psi_S(t)\rangle$
$\hat{O}_I(t) = e^{i\hat{H}_0 t/\hbar}\hat{O}_S e^{-i\hat{H}_0 t/\hbar}$

演化方程：
$i\hbar\partial|\psi_I\rangle/\partial t = \hat{H}_{int,I}(t)|\psi_I\rangle$

### 21.1.5 旋转波近似

在近共振条件下（$|\omega - \omega_0| \ll \omega_0$），可以采用旋转波近似（RWA），忽略快速振荡项。

详细推导：将电场写成正负频率分量：

$\mathbf{E}(t) = \mathbf{E}^+(t) + \mathbf{E}^-(t) = \frac{1}{2}\mathbf{E}_0e^{-i\omega t} + \frac{1}{2}\mathbf{E}_0^*e^{i\omega t}$

相互作用哈密顿量：

$\hat{H}_{int} = -\hat{\mathbf{d}}\cdot\mathbf{E}(t) = -(\hat{d}_+ + \hat{d}_-)\cdot(\mathbf{E}^+ + \mathbf{E}^-)$

其中 $\hat{d}_+ = |e\rangle\langle g|\mathbf{d}_{eg}$，$\hat{d}_- = |g\rangle\langle e|\mathbf{d}_{ge}$。

展开得四项，其中两项以频率 $\omega + \omega_0$ 快速振荡（反旋项），在时间尺度 $1/|\omega - \omega_0|$ 上平均为零。保留慢变项：

$\hat{H}_{int}^{RWA} = -\hbar\Omega/2(\hat{\sigma}_+e^{-i\omega t} + \hat{\sigma}_-e^{i\omega t})$

其中：
- $\Omega = \mathbf{d}_{eg}\cdot\mathbf{E}_0/\hbar$ 是拉比频率（实数，取相位使其为正）
- $\hat{\sigma}_+ = |e\rangle\langle g|$, $\hat{\sigma}_- = |g\rangle\langle e|$ 是升降算符

RWA有效条件：
1. 近共振：$|\omega - \omega_0| \ll \omega_0$
2. 弱场：$\Omega \ll \omega_0$
3. 时间分辨率：$\Delta t \gg 1/\omega_0$

## 21.2 二能级系统与Bloch方程

### 21.2.1 二能级系统模型

考虑最简单的量子系统：二能级原子，基态 $|g\rangle$ 和激发态 $|e\rangle$：

- 能量：$E_g = 0, E_e = \hbar\omega_0$
- 跃迁频率：$\omega_0 = (E_e - E_g)/\hbar$
- 偶极矩：$\mathbf{d} = \langle e|\hat{\mathbf{d}}|g\rangle$

二能级模型的普适性：
1. **原子系统**：精细结构跃迁（如 Na D线）
2. **分子系统**：振动或电子跃迁
3. **量子点**：导带-价带跃迁
4. **NV色心**：自旋态跃迁
5. **超导量子比特**：约瑟夫森结的两个最低能级

完整哈密顿量（旋转坐标系中）：

$\hat{H} = \hbar\omega_0/2 \hat{\sigma}_z - \hbar\Omega/2(\hat{\sigma}_+e^{-i\omega t} + \hat{\sigma}_-e^{i\omega t})$

其中$\hat{\sigma}_z = |e\rangle\langle e| - |g\rangle\langle g|$ 是Pauli-z算符。

**二能级系统的数学结构**：
Hilbert空间 $\mathcal{H} = \text{span}\{|g\rangle, |e\rangle\} \cong \mathbb{C}^2$

基矢的完备性和正交性：
$|g\rangle\langle g| + |e\rangle\langle e| = \hat{I}$
$\langle g|e\rangle = \langle e|g\rangle = 0$
$\langle g|g\rangle = \langle e|e\rangle = 1$

任意算符的展开：
$\hat{O} = O_{gg}|g\rangle\langle g| + O_{ee}|e\rangle\langle e| + O_{ge}|g\rangle\langle e| + O_{eg}|e\rangle\langle g|$

**偶极矩的详细结构**：
复偶极矩矢量：$\mathbf{d} = |\mathbf{d}|e^{i\varphi_d}\hat{\mathbf{e}}_d$

其中：
- $|\mathbf{d}| = |\langle e|e\hat{\mathbf{r}}|g\rangle|$ 是偶极矩大小
- $\hat{\mathbf{e}}_d$ 是偶极矩方向
- $\varphi_d$ 是偶极矩相位

实偶极矩的对称性：
$\mathbf{d}_{eg} = \langle e|\hat{\mathbf{d}}|g\rangle = \langle g|\hat{\mathbf{d}}|e\rangle^* = \mathbf{d}_{ge}^*$

对角元为零（宇称禁戒）：
$\langle g|\hat{\mathbf{d}}|g\rangle = \langle e|\hat{\mathbf{d}}|e\rangle = 0$

**Pauli算符的完整代数**：
定义Pauli算符：
- $\hat{\sigma}_x = \hat{\sigma}_+ + \hat{\sigma}_- = |e\rangle\langle g| + |g\rangle\langle e|$
- $\hat{\sigma}_y = -i(\hat{\sigma}_+ - \hat{\sigma}_-) = -i|e\rangle\langle g| + i|g\rangle\langle e|$
- $\hat{\sigma}_z = |e\rangle\langle e| - |g\rangle\langle g|$

满足对易关系：
$[\hat{\sigma}_i, \hat{\sigma}_j] = 2i\varepsilon_{ijk}\hat{\sigma}_k$

反对易关系：
$\{\hat{\sigma}_i, \hat{\sigma}_j\} = 2\delta_{ij}\hat{I}$

**二能级系统的SU(2)对称性**：
任意二能级态可参数化为：

$|\psi\rangle = \cos(\theta/2)|g\rangle + e^{i\varphi}\sin(\theta/2)|e\rangle$

这对应于Bloch球上的点 $(\theta, \varphi)$。幺正演化对应SO(3)旋转。

**有效二能级系统的实现**：
在多能级系统中，当满足以下条件时可约化为二能级：

1. **大失谐条件**：$|\Delta_{ij}| \gg \Omega_{ij}$（其他跃迁失谐远大于拉比频率）
2. **选择规则**：其他跃迁被选择规则禁戒
3. **频率选择**：激光带宽小于能级间隔

**绝热消除示例**：
三能级系统 $|g\rangle, |e\rangle, |r\rangle$，当 $|\Delta_r| \gg \Omega_r$ 时：

$\hat{H}_{eff} \approx \hbar\delta|e\rangle\langle e| - \hbar\Omega_{eff}/2(|e\rangle\langle g| + h.c.)$

其中 $\delta = \Delta + |\Omega_r|^2/4\Delta_r$ 是AC Stark位移。

### 21.2.2 密度矩阵演化

密度矩阵形式：

$\hat{\rho} = \begin{pmatrix} \rho_{ee} & \rho_{eg} \\ \rho_{ge} & \rho_{gg} \end{pmatrix}$

满足约束：
- $\text{Tr}(\hat{\rho}) = \rho_{ee} + \rho_{gg} = 1$（归一化）
- $\rho_{ge} = \rho_{eg}^*$（厄米性）
- $0 \le \rho_{ee}, \rho_{gg} \le 1$（正定性）
- $\det(\hat{\rho}) \ge 0, \text{Tr}(\hat{\rho}^2) \le 1$（物理密度矩阵条件）

物理意义的深化：
- **纯态**：$\text{Tr}(\hat{\rho}^2) = 1$，如 $|\psi\rangle = \alpha|g\rangle + \beta|e\rangle$
- **混合态**：$\text{Tr}(\hat{\rho}^2) < 1$，表示统计系综或退相干
- **相干性**：$|\rho_{eg}|^2 \le \rho_{ee}\rho_{gg}$（Cauchy-Schwarz不等式）

可观测量的期望值：

$\langle\hat{O}\rangle = \text{Tr}(\hat{\rho}\hat{O}) = \sum_{ij} \rho_{ij}\hat{O}_{ji}$

### 21.2.3 光学Bloch方程推导

从von Neumann方程出发：

$i\hbar\partial\hat{\rho}/\partial t = [\hat{H}, \hat{\rho}]$

详细推导过程：

1. **计算对易子**：
   $[\hat{H}, \hat{\rho}]_{-e} = -\hbar\omega_0\rho_{eg}/2 - \hbar\Omega e^{-i\omega t}(\rho_{gg} - \rho_{ee})/2$
   $[\hat{H}, \hat{\rho}]_{e-} = -\hbar\omega_0\rho_{eg}/2 - \hbar\Omega e^{i\omega t}(\rho_{ee} - \rho_{gg})/2$
   $[\hat{H}, \hat{\rho}]_{--} = -[\hat{H}, \hat{\rho}]_{ee} = -\hbar\Omega(\rho_{eg}e^{i\omega t} - \rho_{ge}e^{-i\omega t})/2$

2. **引入旋转坐标系**：
   $\tilde{\rho}_{eg} = \rho_{eg} e^{i\omega t}$
   $\tilde{\rho}_{ge} = \rho_{ge} e^{-i\omega t}$
   
   这消除了快速振荡的$e^{\pm i\omega_0 t}$因子。

3. **加入弛豫项**（唯象处理）：
   - 纵向弛豫：描述能级布居数衰减
   - 横向弛豫：描述相干性丧失

最终得到光学Bloch方程：

$d\rho_{ee}/dt = i\Omega/2(\tilde{\rho}_{eg} - \tilde{\rho}_{ge}) - \Gamma\rho_{ee}$
$d\tilde{\rho}_{eg}/dt = i\Delta\tilde{\rho}_{eg} + i\Omega/2(\rho_{ee} - \rho_{gg}) - \gamma\tilde{\rho}_{eg}$

其中：
- $\Delta = \omega - \omega_0$ 是失谐
- $\Gamma = 1/T_1$ 是纵向弛豫率（布居数衰减）
- $\gamma = 1/T_2$ 是横向弛豫率（相干性衰减）

弛豫时间关系：
- **纯退相干**：$T_2 = 2T_1$（只有自发辐射）
- **实际系统**：$T_2 < 2T_1$（额外退相干机制）

**完整的矩阵形式推导**：
将密度矩阵写成2×2形式：

$\hat{\rho} = \begin{pmatrix} \rho_{gg} & \rho_{ge} \\ \rho_{eg} & \rho_{ee} \end{pmatrix}$

哈密顿量矩阵：

$\hat{H} = \hbar/2 \begin{pmatrix} -\omega_0 & -\Omega e^{-i\omega t} \\ -\Omega e^{i\omega t} & \omega_0 \end{pmatrix}$

对易子$[\hat{H}, \hat{\rho}]$的矩阵元素：

$[\hat{H}, \hat{\rho}]_{gg} = -\hbar\Omega/2(\rho_{eg}e^{-i\omega t} - \rho_{ge}e^{i\omega t})$
$[\hat{H}, \hat{\rho}]_{ee} = \hbar\Omega/2(\rho_{eg}e^{-i\omega t} - \rho_{ge}e^{i\omega t})$
$[\hat{H}, \hat{\rho}]_{ge} = \hbar\omega_0\rho_{ge} + \hbar\Omega/2(\rho_{gg} - \rho_{ee})e^{i\omega t}$
$[\hat{H}, \hat{\rho}]_{eg} = -\hbar\omega_0\rho_{eg} - \hbar\Omega/2(\rho_{gg} - \rho_{ee})e^{-i\omega t}$

**弛豫项的微观起源**：
从系统-环境耦合的微观模型出发：

$\hat{H}_{total} = \hat{H}_S + \hat{H}_E + \hat{H}_{SE}$

其中$\hat{H}_{SE} = \sum_k g_k(\hat{\sigma}_+\hat{a}_k + \hat{\sigma}_-\hat{a}_k^\dagger)$描述系统与环境模式的耦合。

在Born-Markov近似下，导出主方程：

$d\hat{\rho}_S/dt = -i[\hat{H}_S, \hat{\rho}_S]/\hbar + \mathcal{L}[\hat{\rho}_S]$

其中Lindblad超算符：

$\mathcal{L}[\hat{\rho}] = \Gamma(\bar{n} + 1)(\hat{\sigma}_-\hat{\rho}\hat{\sigma}_+ - \frac{1}{2}\{\hat{\sigma}_+\hat{\sigma}_-,\hat{\rho}\})$
      $+ \Gamma\bar{n}(\hat{\sigma}_+\hat{\rho}\hat{\sigma}_- - \frac{1}{2}\{\hat{\sigma}_-\hat{\sigma}_+,\hat{\rho}\})$

这里$\bar{n}$是热平均光子数，$\Gamma$是自发辐射率。

**纯失相过程**：
额外的失相机制可表示为：

$\mathcal{L}_{deph}[\hat{\rho}] = \gamma_\varphi/2(\hat{\sigma}_z\hat{\rho}\hat{\sigma}_z - \hat{\rho})$

导致：$d\rho_{eg}/dt|_{deph} = -\gamma_\varphi\rho_{eg}/2$

总的横向弛豫率：$\gamma = \Gamma/2 + \gamma_\varphi$

### 21.2.4 Bloch矢量表示

定义Bloch矢量分量：

$u = \tilde{\rho}_{eg} + \tilde{\rho}_{ge} = 2\text{Re}(\tilde{\rho}_{eg})$
$v = i(\tilde{\rho}_{eg} - \tilde{\rho}_{ge}) = -2\text{Im}(\tilde{\rho}_{eg})$
$w = \rho_{ee} - \rho_{gg}$

物理意义：
- **u**：同相分量（与驱动场同相的偶极矩振荡）
- **v**：正交分量（与驱动场90°相位差的偶极矩振荡）
- **w**：反转布居数（w > 0 表示布居数反转）

Bloch方程的矢量形式：

$du/dt = \Delta v - u/T_2$
$dv/dt = -\Delta u + \Omega w - v/T_2$
$dw/dt = -\Omega v - (w - w_0)/T_1$

其中 $w_0 = -1$ 是热平衡时的布居数差。

紧凑形式：

$d\mathbf{R}/dt = \mathbf{\Omega}_{eff} \times \mathbf{R} - \mathbf{\Gamma}\cdot(\mathbf{R} - \mathbf{R}_0)$

其中：
- $\mathbf{\Omega}_{eff} = (\Omega, 0, \Delta)$ 是有效磁场
- $\mathbf{R}_0 = (0, 0, -1)$ 是平衡态
- $\mathbf{\Gamma} = \text{diag}(1/T_2, 1/T_2, 1/T_1)$ 是弛豫张量

### 21.2.5 Bloch球几何

Bloch矢量 $\mathbf{R} = (u, v, w)$ 在单位球内运动：

$|\mathbf{R}|^2 = u^2 + v^2 + w^2 \le 1$

- 纯态：$|\mathbf{R}| = 1$（球面上）
- 混合态：$|\mathbf{R}| < 1$（球内部）
- 最大混合态：$\mathbf{R} = 0$（球心）

**特殊点与态的对应**：
- 北极 (0,0,1)：激发态 $|e\rangle$
- 南极 (0,0,-1)：基态 $|g\rangle$
- 赤道 (cosφ,sinφ,0)：相干叠加态 $(|g\rangle+e^{i\varphi}|e\rangle)/\sqrt{2}$

**动力学解释**：
1. **无弛豫情况**：Bloch矢量绕有效磁场 $\mathbf{\Omega}_{eff}$ 进动
2. **共振情况**（$\Delta=0$）：绕x轴旋转（拉比振荡）
3. **大失谐情况**（$|\Delta|\gg\Omega$）：绕z轴快速进动
4. **弛豫情况**：螺旋向平衡点靠拢

**几何计算示例**：
任意纯态可写为：
$|\psi\rangle = \cos(\theta/2)|g\rangle + e^{i\varphi}\sin(\theta/2)|e\rangle$

对应Bloch矢量：
$\mathbf{R} = (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta)$

其中$\theta, \varphi$是球坐标角度。

## 21.3 拉比振荡与光学章动

### 21.3.1 共振情况（$\Delta = 0$）

在精确共振时，Bloch方程简化为：

$du/dt = -u/T_2$
$dv/dt = \Omega w - v/T_2$
$dw/dt = -\Omega v - (w - w_0)/T_1$

忽略弛豫（$T_1, T_2 \to \infty$），得到：

$d^2w/dt^2 + \Omega^2w = 0$

解为：

$w(t) = \cos(\Omega t)$
$v(t) = -\sin(\Omega t)$
$u(t) = 0$

这描述了布居数的拉比振荡。

### 21.3.2 失谐情况（$\Delta \neq 0$）

考虑失谐但无弛豫的情况，定义广义拉比频率：

$\Omega' = \sqrt{\Omega^2 + \Delta^2}$

布居数演化：

$w(t) = (\Delta^2/\Omega'^2) + (\Omega^2/\Omega'^2)\cos(\Omega't)$

最大激发概率：

$P_{max} = \Omega^2/(\Omega^2 + \Delta^2)$

### 21.3.3 脉冲面积定理

对于时变电场 $\mathbf{E}(t)$，定义脉冲面积：

$\theta = \int_{-\infty}^{\infty} \Omega(t)dt = (1/\hbar)\int_{-\infty}^{\infty} \mathbf{d}\cdot\mathbf{E}(t)dt$

特殊脉冲：
- $\pi$脉冲：$\theta = \pi$，完全反转布居数
- $2\pi$脉冲：$\theta = 2\pi$，返回初态
- $\pi/2$脉冲：$\theta = \pi/2$，创建相干叠加态

**脉冲传播方程**（McCall-Hahn）：
对于在二能级介质中传播的短脉冲：

$\partial\theta/\partial z + (1/v_g)\partial\theta/\partial t = \alpha_0\sin\theta$

其中：
- $v_g$ 是群速度
- $\alpha_0 = 2\pi N|\mathbf{d}|^2/(\hbar c)$ 是吸收系数
- N 是原子数密度

**自感应透明**（SIT）：
- $2\pi$脉冲无损传播
- 形成孤子解：$\theta(z,t) = 4\arctan[\exp(\pm(z-vt)/L)]$
- 脉冲宽度与传播距离成反比

**实验应用**：
1. **量子信息**：量子门操作
2. **光存储**：相干态准备与读出
3. **超快光学**：相干控制

### 21.3.4 章动现象

在旋转坐标系中，Bloch矢量绕有效磁场 $\mathbf{\Omega}_{eff} = (\Omega, 0, \Delta)$ 进动：

$d\mathbf{R}/dt = \mathbf{\Omega}_{eff} \times \mathbf{R} - \mathbf{\Gamma}\cdot\mathbf{R}$

其中 $\mathbf{\Gamma}$ 是弛豫张量。

**几何理解**：
1. **无弛豫情况**：
   - 进动角频率：$\Omega' = |\mathbf{\Omega}_{eff}| = \sqrt{\Omega^2 + \Delta^2}$
   - 进动轴：$\hat{\mathbf{e}} = \mathbf{\Omega}_{eff}/\Omega'$
   - 进动角：$\theta = \arctan(\Omega/\Delta)$

2. **有弛豫情况**：
   - 螺旋轨迹向平衡点收敛
   - 特征时间：$\tau = 1/\sqrt{\Omega'^2 + 1/T_1T_2}$

**章动频率测量**：
通过测量偶极矩振荡的拍频，可以推断失谐：

$\Delta\nu_{beat} = |\Omega' - \omega_0| = |\Delta|$

### 21.3.5 绝热跟随

当场参数缓慢变化时（$|d\Omega/dt| \ll \Omega^2$），系统绝热跟随瞬时本征态：

$w_{ad}(t) \approx -\Delta(t)/\Omega'(t)$

这是绝热快速通道（ARP）和STIRAP技术的基础。

**绝热条件**（Landau-Zener）：

绝热参数：$Q = \Omega^2/(|d\Delta/dt|) \gg 1$

非绝热跃迁概率：

$P_{na} = \exp(-\pi\Omega^2/2|d\Delta/dt|)$

**STIRAP**（受激拉曼绝热通道）：
三能级系统中的相干布居数转移：
1. 使用两个光场：pump和Stokes
2. 反直觉脉冲顺序：Stokes先于pump
3. 通过“dark state”实现100%转移效率
4. 对失谐和场强波动鲁棒

**应用领域**：
- 量子信息处理
- 同位素分离
- 光学冷却
- 相干控制化学

## 21.4 饱和与功率展宽

### 21.4.1 稳态解

在连续波激发下，设 $d/dt = 0$，得到稳态解：

$u_0 = -2\Omega\Delta T_2^2/(1 + \Delta^2 T_2^2 + \Omega^2 T_1 T_2)$
$v_0 = -2\Omega T_2/(1 + \Delta^2 T_2^2 + \Omega^2 T_1 T_2)$
$w_0 = -(1 + \Delta^2 T_2^2)/(1 + \Delta^2 T_2^2 + \Omega^2 T_1 T_2)$

### 21.4.2 饱和强度

定义饱和参数：

$s = \Omega^2 T_1 T_2/(1 + \Delta^2 T_2^2) = I/I_{sat}$

其中饱和强度：

$I_{sat} = \hbar^2/(2|\mathbf{d}|^2 T_1 T_2) \times (1 + \Delta^2 T_2^2)$

稳态激发态布居数：

$\rho_{ee} = s/(2(1 + s))$

### 21.4.3 功率展宽

吸收线型：

$L(\omega) = (\gamma'/2\pi)/[(\omega - \omega_0)^2 + (\gamma'/2)^2]$

其中展宽的线宽：

$\gamma' = \gamma\sqrt{1 + s} = \gamma\sqrt{1 + I/I_{sat}}$

这就是功率展宽效应。

### 21.4.4 均匀与非均匀展宽

**均匀展宽**：所有原子具有相同跃迁频率
- 自然展宽：$\gamma_n = \Gamma/2$
- 碰撞展宽：$\gamma_c \propto p/\sqrt{T}$

**非均匀展宽**：原子跃迁频率分布
- 多普勒展宽：$\Delta\omega_D = \omega_0\sqrt{2kT/mc^2}$
- 晶格无序：场的空间不均匀性

### 21.4.5 烧孔效应

在非均匀展宽介质中，强激光选择性激发特定速度类原子，在吸收谱中"烧"出一个孔：

孔宽度：$\Delta\omega_{hole} \approx \gamma'$
孔深度：$\propto s/(1 + s)$

## 21.5 从经典到量子的过渡

### 21.5.1 半经典理论的局限性

半经典理论无法解释：
1. **自发辐射**：需要量子化光场
2. **光子统计**：亚泊松分布、压缩态
3. **量子关联**：纠缠、非经典关联
4. **真空涨落**：Casimir效应、Lamb位移

### 21.5.2 自发辐射的唯象引入

在半经典框架中，通过唯象地引入自发辐射率 $A_{21}$：

$d\rho_{ee}/dt|_{sp} = -A_{21}\rho_{ee}$

Einstein A系数与偶极矩的关系：

$A_{21} = \omega_0^3|\mathbf{d}|^2/(3\pi\varepsilon_0\hbar c^3)$

修正的Bloch方程：

$dw/dt = -\Omega v - (w + 1)/T_1 - A_{21}(w + 1)/2$

### 21.5.3 速率方程近似

当失谐 $|\Delta| \gg \gamma$ 或强场 $\Omega \gg \gamma$ 时，相干项快速衰减，可采用速率方程：

$d\rho_{ee}/dt = R_{12}\rho_{gg} - R_{21}\rho_{ee} - A_{21}\rho_{ee}$
$d\rho_{gg}/dt = -R_{12}\rho_{gg} + R_{21}\rho_{ee} + A_{21}\rho_{ee}$

其中受激跃迁率：

$R_{12} = R_{21} = (\Omega^2/4)\cdot(\gamma/2\pi)/[(\omega - \omega_0)^2 + (\gamma/2)^2]$

### 21.5.4 与体积渲染方程的联系

将原子介质视为具有频率依赖吸收/发射的体积：

$dL(\omega)/ds = -\sigma_a(\omega)n[1 - \rho_{ee}(\omega)]L(\omega) + \sigma_e(\omega)n\rho_{ee}(\omega)L_e(\omega)$

其中：
- $\sigma_a(\omega), \sigma_e(\omega)$ 是吸收/发射截面
- n 是原子数密度
- $\rho_{ee}(\omega)$ 是频率$\omega$处的激发态布居

这建立了微观原子物理与宏观渲染方程的联系。

### 21.5.5 量子修正与展望

完整量子理论预测的修正：

1. **真空Rabi分裂**：强耦合腔QED
2. **共振荧光谱**：Mollow三峰
3. **光子反聚束**：$g^{(2)}(0) < 1$
4. **纠缠光子对**：参量下转换

这些效应需要第27-28章的完整量子光学处理。

## 本章小结

本章建立了半经典光-物质相互作用理论：

1. **基本框架**：经典场 + 量子化原子
2. **核心方程**：光学Bloch方程描述二能级动力学
3. **重要现象**：
   - 拉比振荡：相干布居转移
   - 饱和效应：强场下的非线性响应
   - 功率展宽：线宽的强度依赖
4. **理论界限**：自发辐射需要量子场论
5. **渲染联系**：微观原子响应→宏观光学性质

关键公式汇总：

- 拉比频率：$\Omega = \mathbf{d}\cdot\mathbf{E}_0/\hbar$
- Bloch方程：$d\mathbf{R}/dt = \mathbf{\Omega}_{eff} \times \mathbf{R} - \mathbf{\Gamma}\cdot\mathbf{R}$
- 饱和强度：$I_{sat} = \hbar^2/(2|\mathbf{d}|^2 T_1 T_2)$
- 功率展宽：$\gamma' = \gamma\sqrt{1 + I/I_{sat}}$

## 练习题

### 基础题

**21.1** 推导二能级系统在共振$\pi$脉冲作用下的演化。证明初始处于基态的原子在脉冲后完全转移到激发态。

<details>
<summary>提示</summary>

使用Bloch方程，设$\Delta = 0$且忽略弛豫。对于矩形脉冲，$\Omega(t) = \Omega_0$（$0 < t < \tau$），脉冲面积$\theta = \Omega_0\tau = \pi$。

</details>

<details>
<summary>答案</summary>

在共振条件下（$\Delta = 0$）且忽略弛豫，Bloch方程为：

$dw/dt = -\Omega v$
$dv/dt = \Omega w$

初始条件：$w(0) = -1, v(0) = 0$（原子在基态）

解得：
$w(t) = -\cos(\Omega t)$
$v(t) = -\sin(\Omega t)$

对于$\pi$脉冲，$\Omega\tau = \pi$，故：
$w(\tau) = -\cos(\pi) = 1$
$v(\tau) = -\sin(\pi) = 0$

因此$\rho_{ee} = (1 + w)/2 = 1$，原子完全在激发态。

</details>

**21.2** 计算氢原子1S-2P跃迁的饱和强度。已知跃迁波长$\lambda = 121.6 \text{ nm}$，偶极矩$|\mathbf{d}| = 2.5 \times 10^{-29} \text{ C}\cdot\text{m}$，$T_1 = T_2 = 1.6 \text{ ns}$。

<details>
<summary>提示</summary>

使用公式 $I_{sat} = \hbar^2/(2|\mathbf{d}|^2 T_1 T_2)$，注意单位换算。

</details>

<details>
<summary>答案</summary>

饱和强度：
$I_{sat} = \hbar^2/(2|\mathbf{d}|^2 T_1 T_2)$

代入数值：
- $\hbar = 1.055 \times 10^{-34} \text{ J}\cdot\text{s}$
- $|\mathbf{d}| = 2.5 \times 10^{-29} \text{ C}\cdot\text{m}$
- $T_1 = T_2 = 1.6 \times 10^{-9} \text{ s}$

$I_{sat} = (1.055 \times 10^{-34})^2/[2 \times (2.5 \times 10^{-29})^2 \times (1.6 \times 10^{-9})^2]$
    $= 1.11 \times 10^{-68}/(3.2 \times 10^{-67})$
    $= 0.35 \text{ W/m}^2$

</details>

**21.3** 证明在大失谐极限（$|\Delta| \gg \Omega, \gamma$）下，二能级系统的AC Stark位移为 $\delta E = \hbar\Omega^2/(4\Delta)$。

<details>
<summary>提示</summary>

考虑有效哈密顿量，使用二阶微扰理论或绝热消除激发态。

</details>

<details>
<summary>答案</summary>

在大失谐下，激发态布居很小，可绝热消除。从Bloch方程：

$d\tilde{\rho}_{eg}/dt \approx 0 \implies \tilde{\rho}_{eg} \approx -\Omega/(2\Delta)\cdot\rho_{gg}$

有效哈密顿量的对角元素给出能级移动：

$\delta E_g = -|\Omega|^2/(4\Delta)$ （基态下移）
$\delta E_e = +|\Omega|^2/(4\Delta)$ （激发态上移）

总的AC Stark位移：$\delta E = \hbar|\Omega|^2/(4\Delta)$

</details>

### 挑战题

**21.4** 分析双色激光场下的二能级系统动力学。考虑两个频率$\omega_1$和$\omega_2$的激光，推导有效拉比频率和共振条件。

<details>
<summary>提示</summary>

使用Floquet理论或旋转波近似的推广。考虑和频与差频过程。

</details>

<details>
<summary>答案</summary>

双色场：$\mathbf{E}(t) = \mathbf{E}_1\cos(\omega_1 t) + \mathbf{E}_2\cos(\omega_2 t)$

在适当的旋转坐标系中，当满足双光子共振条件 $\omega_1 + \omega_2 = 2\omega_0$ 时，有效拉比频率：

$\Omega_{eff} = \Omega_1\Omega_2/(2\delta)$

其中$\delta = \omega_1 - \omega_0$是单光子失谐。

这导致双光子拉比振荡，频率为$\Omega_{eff}$。中间态的虚激发产生AC Stark位移。

</details>

**21.5** 推导包含自发辐射的光学Bloch方程的稳态荧光谱（Mollow谱）。说明为什么半经典理论无法完全解释三峰结构。

<details>
<summary>提示</summary>

计算稳态偶极矩的频谱，考虑强场下的dressed states。

</details>

<details>
<summary>答案</summary>

稳态荧光强度谱：

$S(\omega) \propto \text{Re}[\int_0^\infty \langle\hat{d}^\dagger(t)\hat{d}(0)\rangle e^{i\omega t}dt]$

半经典理论预测中心峰和由于拉比振荡的边带，但无法解释：
1. 三峰的相对强度（需要量子回归定理）
2. 边带的非对称性（需要考虑量子涨落）
3. 光子统计性质（需要二阶相干函数）

完整解释需要量子化光场和主方程方法。

</details>

**21.6** 设计一个绝热快速通道（ARP）脉冲序列，实现99%的布居反转效率。给出脉冲形状和参数选择标准。

<details>
<summary>提示</summary>

使用啁啾脉冲，满足绝热条件 $|d\Delta/dt| \ll \Omega^2$。

</details>

<details>
<summary>答案</summary>

线性啁啾脉冲：
$\Delta(t) = \alpha t$ （$-T/2 < t < T/2$）
$\Omega(t) = \Omega_0\text{sech}(t/\tau)$

绝热条件：$\alpha \ll \Omega_0^2/\tau$

对于99%效率，需要：
1. 扫描范围：$|\alpha T/2| > 5\Omega_0$
2. 脉冲面积：$\int\Omega(t)dt > 5\pi$
3. 绝热参数：$\alpha T^3/(\Omega_0\tau^2) < 0.1$

典型参数：$\Omega_0 = 10\gamma, \tau = 1/\gamma, \alpha T = 20\gamma$

</details>

### 开放性思考题

**21.7** 讨论如何将半经典光-物质相互作用理论应用于计算机图形学中的荧光材料渲染。考虑多能级系统和能量转移过程。

<details>
<summary>思考方向</summary>

- 荧光的Stokes位移
- 量子产率和辐射寿命
- 浓度猝灭效应
- 时间分辨渲染
- 与PBR材质模型的集成

</details>

**21.8** 探讨超快激光脉冲（飞秒量级）与物质相互作用时，半经典理论的修正。这对渲染超快现象有何启示？

<details>
<summary>思考方向</summary>

- 脉冲宽度 $< T_2$时的相干效应
- 非马尔可夫动力学
- 载波包络相位效应
- 高次谐波产生
- 时间分辨光谱的可视化

</details>

## 常见陷阱与错误

1. **旋转波近似的误用**
   - 错误：在强场或大失谐时仍使用RWA
   - 正确：检查条件 $\Omega, |\Delta| \ll \omega_0$

2. **忽略逆过程**
   - 错误：只考虑吸收，忽略受激发射
   - 正确：完整的Bloch方程自动包含两个过程

3. **弛豫时间的混淆**
   - 错误：认为$T_1 = T_2$
   - 正确：$T_2 \le 2T_1$（纯失相导致$T_2 < 2T_1$）

4. **稳态近似的不当使用**
   - 错误：对脉冲激发使用稳态解
   - 正确：检查脉冲宽度 $\gg T_1, T_2$

5. **绝热条件的违反**
   - 错误：快速改变参数时假设绝热跟随
   - 正确：验证 $|\dot{\omega}/\omega^2| \ll 1$

## 最佳实践检查清单

设计半经典光-物质相互作用系统时，确保：

- [ ] 验证半经典近似的有效性（大光子数）
- [ ] 正确选择旋转坐标系（减少数值刚性）
- [ ] 包含所有相关的弛豫过程
- [ ] 检查能量和概率守恒
- [ ] 考虑功率展宽和饱和效应
- [ ] 评估量子修正的必要性
- [ ] 实现数值稳定的积分方案
- [ ] 验证极限情况（弱场、强场、共振、大失谐）

---

*继续到[第22章：偏振光学基础](chapter22.md) →*


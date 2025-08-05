# 第19章：衍射理论与计算方法

本章深入探讨电磁波与边缘、拐角及其他几何不连续性相互作用产生的衍射现象。我们将从几何光学的局限性出发，逐步建立能够准确描述衍射效应的数学框架，并探讨其在计算机图形学中的应用。通过学习几何绕射理论（GTD）、物理绕射理论（PTD/UTD）和Ufimtsev的边缘波理论，您将掌握分析和计算复杂衍射场的方法，并理解如何将这些理论应用于真实感渲染。

## 学习目标

完成本章后，您将能够：
1. 理解几何光学的局限性和衍射的物理本质
2. 掌握GTD的基本原理和射线追踪扩展
3. 应用PTD/UTD计算边缘和拐角的衍射场
4. 理解Ufimtsev理论中的物理光学和边缘波分离
5. 计算复射线和衍射系数
6. 在图形学渲染中实现衍射效应

## 19.1 几何绕射理论（GTD）

几何绕射理论由Joseph Keller在1962年提出，是对传统几何光学的重要扩展。GTD通过引入衍射射线的概念，使得射线光学方法能够处理边缘、拐角和其他几何不连续性产生的衍射效应。

### 19.1.1 几何光学的局限性

传统几何光学基于费马原理和程函方程：

$$(\nabla S)^2 = n^2(\mathbf{r})$$

其中 $S$ 是相位函数（程函），$n(\mathbf{r})$ 是折射率分布。几何光学场可以表示为：

$$u(\mathbf{r}) = A(\mathbf{r}) e^{ik_0 S(\mathbf{r})}$$

在高频极限 $k_0 \to \infty$ 下，这种近似在以下情况失效：

1. **阴影边界**：几何光学预测场的不连续跳变
2. **焦散区域**：振幅 $A(\mathbf{r})$ 变为无穷大
3. **边缘和拐角**：场的导数不连续

实际物理场必须满足：
- 连续性：$u$ 在空间中连续
- 有限性：$|u| < \infty$ 处处成立
- 因果性：满足辐射条件

### 19.1.2 Keller的射线理论扩展

Keller的关键洞察是：衍射可以通过引入新的射线类型来描述。除了传统的入射、反射和折射射线，GTD引入了：

1. **衍射射线**：从边缘或拐角发出
2. **爬行波射线**：沿曲面传播
3. **表面衍射射线**：从表面不连续性发出

衍射射线遵循广义费马原理：

$$\delta \int_A^B n(\mathbf{r}) ds = 0$$

但路径必须通过衍射点 $Q_D$：

$$\delta \left[ \int_A^{Q_D} n ds + \int_{Q_D}^B n ds \right] = 0$$

### 19.1.3 衍射射线的构造

对于直边缘衍射，考虑入射射线 $\hat{\mathbf{s}}^i$ 打到边缘上的点 $Q_D$。衍射射线方向 $\hat{\mathbf{s}}^d$ 满足：

1. **Keller锥条件**：
   $$\hat{\mathbf{s}}^i \cdot \hat{\mathbf{t}} = \hat{\mathbf{s}}^d \cdot \hat{\mathbf{t}} = \cos \beta_0$$
   
   其中 $\hat{\mathbf{t}}$ 是边缘切向，$\beta_0$ 是与边缘的夹角。

2. **衍射锥**：衍射射线形成以边缘为轴、半角为 $\beta_0$ 的圆锥。

3. **局部坐标系**：在衍射点建立：
   - $\hat{\mathbf{t}}$：边缘切向
   - $\hat{\boldsymbol{\phi}}^i$：入射平面法向
   - $\hat{\mathbf{n}}^i = \hat{\boldsymbol{\phi}}^i \times \hat{\mathbf{t}}$

衍射射线方向由方位角 $\phi$ 参数化：

$$\hat{\mathbf{s}}^d = \sin\beta_0(\cos\phi \hat{\mathbf{n}}^i + \sin\phi \hat{\boldsymbol{\phi}}^i) + \cos\beta_0 \hat{\mathbf{t}}$$

#### 边缘局部坐标系的严格定义

为了精确计算衍射，必须建立明确的坐标系统。对于楔形边缘：

1. **楔形几何定义**：
   - 面0：$\phi = 0$ 的表面
   - 面n：$\phi = n\pi$ 的表面（$n = 2\pi/\alpha - 1$）
   - 楔形内部：$0 < \phi < n\pi$

2. **标准坐标基**：
   $$\begin{aligned}
   \hat{\mathbf{t}} &= \text{边缘单位切向量} \\
   \hat{\boldsymbol{\beta}}_0 &= \text{面0的外法向} \\
   \hat{\boldsymbol{\beta}}_n &= \text{面n的外法向}
   \end{aligned}$$

3. **入射和衍射平面**：
   - 入射平面：包含 $\hat{\mathbf{s}}^i$ 和 $\hat{\mathbf{t}}$ 的平面
   - 衍射平面：包含 $\hat{\mathbf{s}}^d$ 和 $\hat{\mathbf{t}}$ 的平面
   - 方位角关系：$\phi^d - \phi^i = \Delta\phi$

#### 广义费马原理的应用

衍射射线路径满足稳相条件：

$$\frac{\partial}{\partial Q_D}\left[\Phi^i(S, Q_D) + k_0|Q_D - P|\right] = 0$$

其中 $\Phi^i$ 是入射场在 $Q_D$ 点的相位。这导致：

1. **路径长度极值**：衍射路径是从源点经边缘到观察点的极值路径
2. **斯涅尔定律的推广**：入射角等于衍射角（相对于边缘）
3. **因果性**：衍射射线沿最小时间路径传播

### 19.1.4 衍射系数的渐近形式

GTD中的衍射场表示为：

$$u^d(P) = u^i(Q_D) \cdot D \cdot A^d(P) \cdot e^{ik_0 s}$$

其中：
- $u^i(Q_D)$：衍射点的入射场
- $D$：二阶衍射系数张量
- $A^d(P)$：衍射射线的扩散因子
- $s$：从 $Q_D$ 到观察点 $P$ 的距离

对于完美导体楔形（楔角 $\alpha$），标量衍射系数为：

$$D(\phi, \phi'; \beta_0; \alpha) = -\frac{e^{i\pi/4}}{2\sqrt{2\pi k_0}} \frac{1}{\sin\beta_0} \left[ \frac{1}{\cos\frac{\pi(\phi-\phi')}{\alpha} - \cos\frac{\pi}{\alpha}} + \frac{1}{\cos\frac{\pi(\phi+\phi')}{\alpha} - \cos\frac{\pi}{\alpha}} \right]$$

其中 $\phi'$ 是入射方位角。

扩散因子考虑射线管的展开：

$$A^d(P) = \sqrt{\frac{\rho}{s(\rho + s)}}$$

其中 $\rho$ 是边缘的曲率半径（直边缘时 $\rho \to \infty$）。

#### 衍射系数的物理解释

衍射系数 $D$ 的结构揭示了重要的物理含义：

1. **频率依赖性**：$D \propto k_0^{-1/2}$，表明衍射强度随频率增加而减弱
2. **几何奇点**：当 $\phi = \pm\phi' + 2n\pi/N$ 时出现奇点，对应于几何光学射线方向
3. **对称性**：$D(\phi, \phi') = D(-\phi, -\phi')$，反映了互易原理

衍射系数可以分解为软边界和硬边界贡献：

$$D = D_s \hat{\mathbf{e}}_s \otimes \hat{\mathbf{e}}_s + D_h \hat{\mathbf{e}}_h \otimes \hat{\mathbf{e}}_h$$

其中 $\hat{\mathbf{e}}_s$ 和 $\hat{\mathbf{e}}_h$ 分别是软（TE）和硬（TM）极化的单位矢量。

### 19.1.5 GTD的高阶修正

标准GTD是 $k_0 \to \infty$ 的首阶渐近近似。高阶修正提供更精确的结果：

$$u^d = u_0^d \left[1 + \frac{1}{ik_0}\mathcal{L}_1 + \frac{1}{(ik_0)^2}\mathcal{L}_2 + \cdots \right]$$

其中 $\mathcal{L}_n$ 是 $n$ 阶微分算子，作用于振幅和相位函数。

#### 曲率修正

对于弯曲边缘，需要考虑边缘曲率 $\kappa(l)$ 和挠率 $\tau(l)$ 的影响：

$$D_{curved} = D_{straight} \cdot \exp\left[i\int_0^s \kappa(l') \cos\beta(l') dl'\right]$$

这导致额外的相位积累和振幅调制。

#### 高阶衍射过程

除了单次衍射，GTD框架还包括：

1. **顶点衍射**：三维拐角产生的球面波
   $$u^v \sim \frac{e^{ik_0 r}}{r} V(\hat{\mathbf{r}}, \hat{\mathbf{s}}^i)$$
   
   其中 $V$ 是顶点衍射系数，依赖于观察和入射方向。

2. **爬行波**：沿凸曲面传播的表面波
   $$u^{creeping} = \sum_m A_m(s) e^{i\nu_m s} e^{-\alpha_m s}$$
   
   其中 $\nu_m$ 是复传播常数，$\alpha_m$ 是衰减系数。

3. **斜率衍射**：表面曲率不连续处的贡献
   $$D_{slope} \sim \frac{1}{k_0} \frac{\partial n}{\partial s}\bigg|_{edge}$$

### 19.1.6 GTD的适用范围与局限

#### 有效性条件

GTD的准确性依赖于以下条件：

1. **高频条件**：$k_0 a \gg 1$，其中 $a$ 是最小几何特征尺寸
2. **局部性条件**：衍射点附近的几何形状决定衍射场
3. **远场条件**：观察距离 $r \gg \lambda$
4. **避开奇点**：远离阴影边界和焦散区

#### 失效区域

GTD在以下区域失效，需要特殊处理：

1. **阴影边界（SB）和反射边界（RB）**：
   - 衍射系数出现奇点
   - 场的不连续性
   - 需要UTD过渡函数修正

2. **焦散区域**：
   - 射线汇聚导致振幅无穷大
   - 需要复射线或波动方法

3. **边缘附近**：
   - $r < \lambda$ 时射线近似失效
   - 需要完整波动方程求解

4. **多重衍射区域**：
   - 高阶衍射贡献可能显著
   - 计算复杂度急剧增加

#### 与其他方法的衔接

为了克服GTD的局限性，常与其他方法结合：

$$u_{total} = \begin{cases}
u_{MoM} & r < \lambda \\
u_{PO} + u_{PTD} & \lambda < r < 10\lambda \\
u_{GTD/UTD} & r > 10\lambda
\end{cases}$$

其中：
- MoM：矩量法（Method of Moments）
- PO：物理光学（Physical Optics）
- PTD：物理衍射理论（Physical Theory of Diffraction）

过渡区域需要特殊的匹配技术确保场的连续性。

GTD在阴影边界和反射边界附近失效，因为衍射系数包含奇点。一致性衍射理论（UTD）通过引入过渡函数修正了这些问题，使得场在所有区域都保持连续。

### 19.2.1 过渡区域的修正

在阴影边界附近，定义距离参数：

$$L = s \sin^2 \beta_0$$

其中 $s$ 是从衍射点到观察点的距离。过渡区域的特征长度为：

$$\delta \sim \sqrt{\frac{\lambda L}{2\pi}}$$

引入Fresnel积分参数：

$$\xi = \sqrt{\frac{2k_0 L}{\pi}} \sin\left(\frac{\phi - \phi'}{2}\right)$$

过渡函数定义为：

$$F(\xi) = 2i\sqrt{\pi} \xi e^{i\xi^2} \int_\xi^\infty e^{-it^2} dt$$

这确保了：
- $|\xi| \gg 1$：退化为GTD结果
- $\xi = 0$：场连续通过边界
- $|\xi| \ll 1$：平滑过渡

### 19.2.2 一致性衍射理论（UTD）

UTD衍射系数包含过渡函数：

$$D_{UTD} = -\frac{e^{i\pi/4}}{2\sqrt{2\pi k_0}} \frac{1}{\sin\beta_0} \sum_{n=1,2} \cot\left(\frac{\pi + (-1)^n(\phi - \phi')}{2N}\right) F(k_0 L a_n^{\pm})$$

其中：
- $N = \alpha/\pi$（楔形参数）
- $a_n^{\pm}$ 是距离参数，依赖于几何配置

对于完美导体，软边界条件（TE极化）和硬边界条件（TM极化）的系数不同：

$$D_s = D_{UTD} \cdot \text{sgn}(\text{入射角})$$
$$D_h = -D_{UTD} \cdot \text{sgn}(\text{反射角})$$

### 19.2.3 斜入射的衍射公式

当入射射线不垂直于边缘时（$\beta_0 \neq \pi/2$），需要考虑三维效应。定义：

1. **边缘固定坐标系**：
   - $\hat{\mathbf{e}}_1$：边缘切向
   - $\hat{\mathbf{e}}_2$：边缘法向（指向楔内）
   - $\hat{\mathbf{e}}_3 = \hat{\mathbf{e}}_1 \times \hat{\mathbf{e}}_2$

2. **入射和衍射角**：
   $$\hat{\mathbf{s}}^i = \sin\beta_0^i(\cos\phi^i \hat{\mathbf{e}}_2 + \sin\phi^i \hat{\mathbf{e}}_3) + \cos\beta_0^i \hat{\mathbf{e}}_1$$
   $$\hat{\mathbf{s}}^d = \sin\beta_0^d(\cos\phi^d \hat{\mathbf{e}}_2 + \sin\phi^d \hat{\mathbf{e}}_3) + \cos\beta_0^d \hat{\mathbf{e}}_1$$

3. **斜入射衍射系数**：
   $$D_{3D} = D_{2D} \cdot \sqrt{\frac{\sin\beta_0^i}{\sin\beta_0^d}} \cdot P(\beta_0^i, \beta_0^d)$$
   
   其中 $P$ 是极化修正因子。

### 19.2.4 多重衍射的处理

实际场景中常遇到多次衍射：

1. **级联衍射**：射线依次经过多个边缘
   $$u^{(n)} = u^{(0)} \prod_{j=1}^n D_j \cdot A_j \cdot e^{ik_0 s_j}$$

2. **耦合效应**：考虑边缘间的相互作用
   - 近场耦合：边缘距离 $< \lambda$
   - 远场近似：独立衍射叠加

3. **高阶衍射**：
   - 拐角衍射：需要特殊的顶点衍射系数
   - 曲面-边缘过渡：爬行波与衍射波的转换

收敛性分析表明，对于大多数几何配置，考虑到第二阶衍射即可达到 $O(k_0^{-2})$ 的精度。

#### 双重衍射的UTD公式

对于两个边缘的连续衍射，总场为：

$$u^{dd}(P) = u^i(Q_1) \cdot D_1 \cdot T_{12} \cdot D_2 \cdot A_{tot} \cdot e^{ik_0(s_1 + s_2)}$$

其中：
- $Q_1, Q_2$：第一和第二衍射点
- $T_{12}$：边缘间的传输因子
- $A_{tot}$：总的扩散因子

传输因子考虑了中间路径的相位和振幅变化：

$$T_{12} = \sqrt{\frac{s_1}{s_1 + s_2}} \cdot \exp\left[i\Phi_{12}\right]$$

#### 多边缘系统的矩阵方法

对于复杂的多边缘系统，可以使用散射矩阵方法：

$$\mathbf{u}^{out} = \mathbf{S} \cdot \mathbf{u}^{in}$$

其中散射矩阵元素：

$$S_{ij} = \begin{cases}
D_{ii} & i = j \text{（自衍射）} \\
T_{ij} D_{jj} & i \neq j \text{（交叉衍射）}
\end{cases}$$

总场通过迭代求解：

$$\mathbf{u} = (\mathbf{I} - \mathbf{S})^{-1} \mathbf{u}^{inc}$$

### 19.2.5 UTD的数值实现考虑

#### 过渡函数的高效计算

UTD的计算瓶颈通常在过渡函数 $F(\xi)$ 的评估。常用方法：

1. **级数展开**（小参数 $|\xi| < 2$）：
   $$F(\xi) = \sqrt{\pi} \xi \sum_{n=0}^{\infty} \frac{(i\xi^2)^n}{n!(2n+1)}$$

2. **渐近展开**（大参数 $|\xi| > 2$）：
   $$F(\xi) \approx 1 + \sum_{n=1}^{\infty} \frac{i^n (2n-1)!!}{2^n n! \xi^{2n}}$$

3. **有理函数近似**：
   $$F(\xi) \approx \frac{P_m(\xi)}{Q_n(\xi)} + \epsilon$$
   
   其中 $\epsilon < 10^{-6}$ 对于适当选择的多项式阶数。

#### 射线追踪的优化策略

1. **空间索引**：
   - 边缘的层次包围盒（BVH）
   - 八叉树加速结构
   - 距离场预计算

2. **重要性采样**：
   $$p(\phi) \propto |D(\phi)|^2 \cdot \cos(\phi - \phi_{obs})$$
   
   集中计算资源在对观察点贡献大的方向。

3. **自适应细分**：
   - 基于场变化率的边缘细分
   - 视角相关的LOD系统
   - 误差驱动的递归细化

#### 边界处理的鲁棒性

在阴影和反射边界附近，数值稳定性至关重要：

1. **连续性测试**：
   $$\lim_{\phi \to \phi_{SB}} u_{lit} = \lim_{\phi \to \phi_{SB}} u_{shadow}$$

2. **奇点消除**：
   使用 L'Hôpital 法则处理 $0/0$ 型不定式：
   $$\lim_{\xi \to 0} \frac{F(\xi)}{\xi} = \sqrt{\pi}$$

3. **相位展开**：
   避免 $2\pi$ 跳变导致的不连续性，维护连续相位历史。

## 19.3 Ufimtsev的边缘波理论

Pyotr Ufimtsev在1960年代发展的物理衍射理论（PTD）提供了一种不同于GTD的方法。PTD基于物理光学（PO）近似，并明确分离出边缘贡献，这一方法后来成为隐身技术发展的理论基础。

### 19.3.1 物理光学近似

物理光学近似假设表面电流仅在照明区存在：

$$\mathbf{J}_s = \begin{cases}
2\hat{\mathbf{n}} \times \mathbf{H}^i & \text{照明区} \\
0 & \text{阴影区}
\end{cases}$$

对于平面波入射 $\mathbf{E}^i = E_0 e^{ik_0 \hat{\mathbf{k}}^i \cdot \mathbf{r}}$，散射场为：

$$\mathbf{E}^{PO}(\mathbf{r}) = \frac{ik_0 Z_0}{4\pi} \int_S \mathbf{J}_s \times \hat{\mathbf{R}} \frac{e^{ik_0 R}}{R} dS'$$

其中 $\mathbf{R} = \mathbf{r} - \mathbf{r}'$，$Z_0$ 是自由空间阻抗。

PO的主要缺陷：
1. 在阴影边界产生不连续
2. 忽略了边缘的贡献
3. 不满足边界条件的连续性

### 19.3.2 边缘波的分离

Ufimtsev的关键洞察是将总场分解为：

$$\mathbf{E}^{total} = \mathbf{E}^{PO} + \mathbf{E}^{edge}$$

边缘波 $\mathbf{E}^{edge}$ 补偿PO近似的不足。对于直边缘，边缘电流可以表示为：

$$\mathbf{J}^{edge} = I(l) \delta(n) \hat{\mathbf{t}}$$

其中：
- $I(l)$：沿边缘的电流分布
- $\delta(n)$：垂直于边缘的δ函数
- $\hat{\mathbf{t}}$：边缘切向

边缘散射场：

$$\mathbf{E}^{edge} = \int_{edge} I(l') G(\mathbf{r}, l') dl'$$

其中 $G$ 是并矢格林函数。

### 19.3.3 物理衍射理论（PTD）

PTD通过匹配边界条件确定边缘电流。对于完美导体楔形，非均匀电流为：

$$I^{PTD} = I^{exact} - I^{PO}$$

其中：
- $I^{exact}$：精确解的边缘电流（从Sommerfeld解导出）
- $I^{PO}$：物理光学电流的等效边缘贡献

PTD衍射系数：

$$D^{PTD} = D^{exact} - D^{PO}$$

对于远场，这给出：

$$D^{PTD}(\phi, \phi'; \alpha) = \frac{e^{i\pi/4}}{\sqrt{2\pi k_0}} \left[ f(\phi - \phi') - f(\phi + \phi') \right]$$

其中：

$$f(\psi) = \frac{\sin(\pi/N)}{N[\cos(\pi/N) - \cos(\psi/N)]}$$

$N = 2\pi/\alpha$ 是楔形参数。

### 19.3.4 与GTD/UTD的关系

PTD和GTD/UTD在高频极限下是等价的，但有重要区别：

1. **物理解释**：
   - GTD：纯射线光学扩展
   - PTD：物理光学 + 边缘修正

2. **适用范围**：
   - GTD：$k_0 \to \infty$ 的渐近展开
   - PTD：中高频段（$k_0 a > 1$，$a$ 是特征尺寸）

3. **计算复杂度**：
   - GTD/UTD：$O(N_{edges})$
   - PTD：$O(N_{surface} + N_{edges})$

4. **精度比较**：
   $$|E^{GTD} - E^{exact}| \sim O(k_0^{-3/2})$$
   $$|E^{PTD} - E^{exact}| \sim O(k_0^{-1})$$

转换关系：

$$D^{UTD} = D^{PTD} \cdot T(\xi)$$

其中 $T(\xi)$ 是过渡函数，确保在阴影边界的连续性。

## 19.4 复射线与衍射系数

复射线理论将实射线推广到复空间，为处理焦散、爬行波和隧穿等现象提供了统一框架。这种方法在处理阴影区域的指数衰减场和过渡区域的振荡行为时特别有效。

### 19.4.1 复射线的概念

传统射线理论中，射线路径由实空间中的轨迹描述。复射线理论允许射线在复空间中传播，其位置和方向可以具有虚部：

$$\mathbf{r}_c = \mathbf{r}_r + i\mathbf{r}_i$$
$$\hat{\mathbf{s}}_c = \hat{\mathbf{s}}_r + i\hat{\mathbf{s}}_i$$

复射线场的一般形式：

$$u_c(\mathbf{r}) = A_c(\mathbf{r}) \exp[ik_0 S_c(\mathbf{r})]$$

其中复程函 $S_c = S_r + iS_i$ 满足复程函方程：

$$(\nabla S_c)^2 = n^2(\mathbf{r}_c)$$

#### 高斯束表示

最常用的复射线是高斯束，其在垂直于传播方向的平面上具有高斯分布：

$$u_{GB}(\mathbf{r}) = \frac{A_0}{\sqrt{q(s)}} \exp\left[ik_0\left(s + \frac{(\mathbf{r}_\perp)^2}{2q(s)}\right)\right]$$

其中：
- $s$：沿射线的弧长
- $\mathbf{r}_\perp$：垂直于射线的位移
- $q(s)$：复曲率参数

复曲率参数满足Riccati方程：

$$\frac{dq}{ds} = 1 + \frac{q^2}{R^2(s)}$$

初始条件：$q(0) = q_0 = w_0^2/(ik_0)$，其中 $w_0$ 是束腰半径。

#### 复射线的物理意义

1. **实部**：描述能量集中的区域
2. **虚部**：描述场的指数衰减或增长
3. **相位**：包含几何相位和Gouy相位

在阴影区域，复射线可以描述隧穿效应：

$$u_{shadow} \sim \exp[-k_0 \text{Im}(S_c)]$$

### 19.4.2 衍射系数的解析延拓

衍射系数在复角度平面上是解析函数，这允许我们通过解析延拓处理各种特殊情况。

#### Sommerfeld积分表示

对于楔形衍射，精确解可以表示为Sommerfeld积分：

$$u(\mathbf{r}, \phi) = \frac{1}{2\pi i} \int_{\Gamma} \frac{e^{ik_0 r \cos(\theta - \phi)}}{\sin(\theta \pi/\alpha)} d\theta$$

其中积分路径 $\Gamma$ 是Sommerfeld轮廓。

#### 衍射系数的解析结构

衍射系数 $D(\phi, \phi'; \alpha)$ 在复 $\phi$ 平面上具有：

1. **极点**：对应于几何光学射线
   - 位置：$\phi = \pm \phi' + 2n\alpha$（$n \in \mathbb{Z}$）
   - 留数：给出反射/透射系数

2. **分支点**：对应于掠射
   - 位置：$\phi = 0, \pm \alpha$
   - 分支切割：定义场的物理片

3. **鞍点**：对应于衍射射线
   - 通过最陡下降法确定
   - 给出GTD射线方向

#### Watson变换

将Sommerfeld积分转换为级数表示：

$$u = \sum_n \text{Res}[f(\theta)] + \int_{\text{SDP}} f(\theta) d\theta$$

其中：
- 第一项：几何光学贡献（极点留数）
- 第二项：衍射贡献（鞍点积分）

### 19.4.3 因果性与解析性

物理场必须满足因果性，这对衍射系数施加了解析性约束。

#### Kramers-Kronig关系

对于频率相关的衍射系数 $D(\omega)$：

$$\text{Re}[D(\omega)] = \frac{1}{\pi} \mathcal{P} \int_{-\infty}^{\infty} \frac{\text{Im}[D(\omega')]}{\omega' - \omega} d\omega'$$

$$\text{Im}[D(\omega)] = -\frac{1}{\pi} \mathcal{P} \int_{-\infty}^{\infty} \frac{\text{Re}[D(\omega')]}{\omega' - \omega} d\omega'$$

其中 $\mathcal{P}$ 表示主值积分。

#### 时域因果性

脉冲响应 $h(t)$ 必须满足：

$$h(t) = 0, \quad t < 0$$

这要求频域传递函数 $H(\omega)$ 在上半复平面解析。

对于衍射问题：

$$D(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} D(\omega) e^{-i\omega t} d\omega$$

因果性要求积分路径可以在上半平面闭合。

### 19.4.4 数值计算方法

实际计算衍射系数需要高效的数值方法，特别是在过渡区域。

#### 特殊函数计算

许多衍射问题涉及特殊函数：

1. **Fresnel积分**：
   $$C(\xi) + iS(\xi) = \int_0^\xi e^{i\pi t^2/2} dt$$
   
   快速计算：使用有理函数近似或级数展开

2. **过渡函数**：
   $$F(\xi) = 2i\sqrt{\pi} \xi e^{i\xi^2} \text{erfc}(-i\xi)$$
   
   其中 $\text{erfc}$ 是互补误差函数

3. **Fock函数**：用于爬行波
   $$v(z) = \frac{1}{\sqrt{\pi}} \int_{\infty e^{-i\pi/3}}^{\infty} e^{zt - t^3/3} dt$$

#### 快速算法

1. **查表法**：预计算常用配置的衍射系数
   - 内存需求：$O(N_\phi \times N_{\phi'} \times N_\alpha)$
   - 精度：通过插值达到 $10^{-4}$

2. **渐近展开**：
   - 大参数：$|k_0 L| \gg 1$
   $$D \approx D_0 + \frac{D_1}{k_0 L} + \frac{D_2}{(k_0 L)^2} + \cdots$$
   
   - 小参数：$|k_0 L| \ll 1$
   $$D \approx a_0 + a_1(k_0 L) + a_2(k_0 L)^2 + \cdots$$

3. **数值积分**：
   - 自适应Gauss-Kronrod积分
   - 复平面上的轮廓积分
   - FFT加速卷积计算

#### 矢量衍射系数

电磁场的矢量性质需要二阶张量衍射系数：

$$\mathbf{E}^d = \overleftrightarrow{D} \cdot \mathbf{E}^i$$

对于边缘衍射：

$$\overleftrightarrow{D} = D_\parallel \hat{\mathbf{e}}_\parallel \hat{\mathbf{e}}_\parallel + D_\perp \hat{\mathbf{e}}_\perp \hat{\mathbf{e}}_\perp$$

其中：
- $\hat{\mathbf{e}}_\parallel$：平行于入射面的单位矢量
- $\hat{\mathbf{e}}_\perp$：垂直于入射面的单位矢量

软硬边界条件下：

$$D_\perp^{soft} = -D_{scalar}$$
$$D_\parallel^{hard} = D_{scalar}$$

#### 误差估计

数值计算的误差来源：

1. **截断误差**：级数或积分截断
   $$\epsilon_{trunc} \sim O(N^{-p})$$
   
2. **离散化误差**：数值积分
   $$\epsilon_{disc} \sim O(h^q)$$
   
3. **舍入误差**：浮点运算
   $$\epsilon_{round} \sim O(\epsilon_{machine})$$

总误差：
$$\epsilon_{total} \leq C_1 N^{-p} + C_2 h^q + C_3 \epsilon_{machine}$$

优化策略：平衡各项误差，使 $N^{-p} \approx h^q \approx \sqrt{\epsilon_{machine}}$。

## 19.5 计算机图形学中的衍射效应

衍射效应虽然在日常场景中通常较微弱，但在某些情况下对视觉真实感至关重要。本节探讨如何将衍射理论集成到现代渲染系统中，平衡物理准确性和计算效率。

### 19.5.1 衍射在渲染中的重要性

#### 视觉显著的衍射场景

1. **锐利边缘和刀口**：
   - 剃刀刀片边缘的彩虹色条纹
   - 建筑物边缘在逆光下的光晕
   - 特征尺度：$d \sim 1-10\mu m$

2. **狭缝和小孔**：
   - 百叶窗的衍射图案
   - 针孔相机的艾里斑
   - 特征尺度：$d \sim 0.1-1mm$

3. **光学仪器**：
   - 望远镜的衍射极限
   - 相机光圈的星芒效果
   - 显微镜的分辨率限制

4. **特殊材料**：
   - CD/DVD的彩虹反射
   - 衍射光栅和全息图
   - 蝴蝶翅膀等生物结构

#### 衍射的尺度分析

判断是否需要考虑衍射的准则：

$$\text{菲涅尔数} = \frac{a^2}{L\lambda}$$

其中：
- $a$：特征尺寸（边缘到观察点的横向距离）
- $L$：传播距离
- $\lambda$：波长

分类：
- $F \gg 1$：几何光学区域，衍射可忽略
- $F \sim 1$：菲涅尔衍射区域，需要考虑
- $F \ll 1$：夫琅禾费衍射区域，远场近似

### 19.5.2 基于GTD的射线追踪扩展

#### 衍射射线的生成

扩展传统射线追踪器以支持衍射：

1. **边缘检测与参数化**：
   ```
   对每个几何体：
     提取锐利边缘（二面角 < 阈值）
     参数化边缘曲线 r(t)
     存储边缘切向 t(t) 和法向信息
   ```

2. **衍射点采样**：
   - 均匀采样：$t_i = i\Delta t$
   - 自适应采样：基于观察角度和距离
   - 重要性采样：基于入射场强度

3. **衍射锥构造**：
   对每个衍射点 $Q_D$：
   ```
   计算入射角 β₀ = arccos(s^i · t)
   构造Keller锥（半角 = β₀）
   在锥上采样衍射方向
   ```

#### 射线树的扩展

传统射线树：`入射 → 反射/折射 → ...`

GTD射线树：`入射 → 反射/折射/衍射 → ...`

递归深度控制：
- 几何光学射线：深度 $\leq D_{geo}$
- 衍射射线：深度 $\leq D_{diff}$（通常 $D_{diff} < D_{geo}$）

#### 衍射贡献的计算

对于点 $P$ 的衍射场：

$$L_d(P) = \sum_{edges} \int_{edge} L_i(Q) \cdot D(Q,P) \cdot V(Q,P) \cdot G(Q,P) dl$$

其中：
- $L_i(Q)$：边缘点 $Q$ 的入射辐射度
- $D(Q,P)$：衍射系数（GTD/UTD）
- $V(Q,P)$：可见性函数
- $G(Q,P)$：几何衰减因子

离散化：

$$L_d(P) \approx \sum_{i} L_i(Q_i) \cdot D_i \cdot V_i \cdot G_i \cdot \Delta l_i$$

### 19.5.3 实时衍射近似

实时渲染需要高效的近似方法，牺牲一定精度换取性能。

#### 屏幕空间衍射

1. **边缘检测**（后处理）：
   - 深度不连续：$|\nabla z| > \epsilon_z$
   - 法向不连续：$|\nabla \mathbf{n}| > \epsilon_n$
   - 材质边界：ID buffer

2. **衍射核卷积**：
   ```glsl
   vec3 diffraction = vec3(0);
   for(int i = -N; i <= N; i++) {
     vec2 uv = texCoord + i * edgeNormal * spacing;
     float phase = k * i * spacing;
     vec3 color = texture(colorBuffer, uv).rgb;
     diffraction += color * sinc(phase) * exp(-abs(i)/sigma);
   }
   ```

3. **色散近似**：
   - 使用RGB通道的不同相位
   - 简化的波长依赖：$k_R < k_G < k_B$

#### 预计算衍射纹理

对于已知几何形状，预计算衍射图案：

1. **刀口衍射**：
   $$I(x) = I_0 \left| \frac{1 + i}{2} - (C(\xi) + iS(\xi)) \right|^2$$
   
   存储为1D查找表

2. **狭缝衍射**：
   $$I(\theta) = I_0 \left( \frac{\sin(\beta)}{\beta} \right)^2, \quad \beta = \frac{\pi a \sin\theta}{\lambda}$$
   
   参数化存储：$(a/\lambda, \theta) \to I$

3. **圆孔衍射**（艾里图案）：
   $$I(r) = I_0 \left( \frac{2J_1(x)}{x} \right)^2, \quad x = \frac{2\pi a r}{\lambda f}$$

#### GPU加速实现

利用GPU并行性加速衍射计算：

1. **边缘缓冲区**：
   ```hlsl
   struct Edge {
     float3 start, end;    // 端点
     float3 tangent;       // 切向
     float2 angles;        // 楔角
     uint materialID;      // 材质
   };
   StructuredBuffer<Edge> edgeBuffer;
   ```

2. **并行衍射计算**：
   ```hlsl
   [numthreads(32,32,1)]
   void DiffractionCS(uint3 id : SV_DispatchThreadID) {
     float3 P = GetWorldPos(id.xy);
     float3 diffraction = 0;
     
     for(uint i = 0; i < numEdges; i++) {
       Edge e = edgeBuffer[i];
       float3 contrib = ComputeDiffraction(P, e);
       diffraction += contrib;
     }
     
     outputTexture[id.xy] += float4(diffraction, 1);
   }
   ```

3. **层次化加速**：
   - 边缘的空间划分（BVH/octree）
   - 屏幕空间瓦片化
   - 重要性驱动的LOD

### 19.5.4 与波动光学方法的比较

#### 方法对比

| 方法 | 精度 | 速度 | 内存 | 适用场景 |
|------|------|------|------|----------|
| 完整波动方程 | 最高 | 最慢 | 最大 | 科学计算 |
| BPM/FDTD | 高 | 慢 | 大 | 光学设计 |
| 物理光学(PO) | 中 | 中 | 中 | 雷达散射 |
| GTD/UTD | 中 | 快 | 小 | 实时渲染 |
| 屏幕空间近似 | 低 | 最快 | 最小 | 游戏/VR |

#### 混合方法

结合不同方法的优势：

1. **近场-远场分离**：
   - 近场（$r < 10\lambda$）：波动方程
   - 过渡区：物理光学
   - 远场（$r > 100\lambda$）：GTD/UTD

2. **频率分解**：
   - 低频（$ka < 1$）：矩量法(MoM)
   - 中频（$1 < ka < 100$）：物理光学
   - 高频（$ka > 100$）：射线方法

3. **自适应精度**：
   ```
   if (visualImportance > threshold1) {
     使用完整波动方法
   } else if (visualImportance > threshold2) {
     使用GTD/UTD
   } else {
     使用屏幕空间近似
   }
   ```

#### 验证与校准

确保近似方法的准确性：

1. **解析解对比**：
   - 半平面衍射（Sommerfeld解）
   - 楔形衍射（精确级数解）
   - 圆柱衍射（贝塞尔函数解）

2. **数值收敛性**：
   $$\text{误差} = \|u_{approx} - u_{exact}\|_2 / \|u_{exact}\|_2$$
   
   要求：误差 $< 5\%$ 在视觉重要区域

3. **感知度量**：
   - SSIM（结构相似性）
   - 色差 $\Delta E$
   - 用户研究验证

### 19.5.5 物理效应与应用

衍射不仅是一个理论概念，它在许多实际渲染场景中产生重要的视觉效果。理解这些效应有助于决定何时以及如何在渲染中包含衍射计算。

#### 大气光学中的衍射

大气中的衍射效应经常被观察到，特别是在特定的气象条件下：

1. **光晕（Corona）现象**：
   - 成因：云层中均匀水滴的衍射
   - 特征：围绕光源的彩色同心圆环
   - 理论：米氏散射理论的小粒子极限
   
   强度分布：
   $$I(\theta) = I_0 \left| \sum_{n} a_n J_n(kr\sin\theta) \right|^2$$
   
   其中 $J_n$ 是贝塞尔函数，$a_n$ 是米氏系数。

2. **光荣（Glory）现象**：
   - 成因：背散射的复杂干涉
   - 观察：飞机影子周围的彩虹圆环
   - 机制：表面波、内反射和衍射的组合
   
   近似模型：
   $$I_{glory}(\theta) \propto \sum_{p=0}^{\infty} |A_p|^2 \exp[-(\theta/\theta_p)^2]$$
   
   其中 $p$ 是内反射次数。

3. **衍射辉光（Diffraction Spikes）**：
   - 成因：大气湍流造成的相位屏
   - 特征：星光的闪烁和径向辉光
   - 建模：随机相位屏方法
   
   相位扰动：
   $$\phi(\mathbf{r}) = \int C_n^2(\mathbf{k}) |\mathbf{k}|^{-11/3} e^{i\mathbf{k}\cdot\mathbf{r}} d^2\mathbf{k}$$

#### 光学系统中的衍射极限

理解衍射极限对模拟真实光学系统至关重要：

1. **瑞利判据**：
   两个点源刚好可分辨的条件：
   $$\theta_{min} = 1.22 \frac{\lambda}{D}$$
   
   其中 $D$ 是孔径直径。

2. **点扩散函数（PSF）**：
   理想圆孔径的PSF（艾里函数）：
   $$\text{PSF}(r) = \left( \frac{2J_1(\pi D r/\lambda f)}{\pi D r/\lambda f} \right)^2$$
   
   其中 $f$ 是焦距。

3. **光学传递函数（OTF）**：
   $$\text{OTF}(\nu) = \mathcal{F}\{\text{PSF}\}(\nu)$$
   
   截止频率：$\nu_{cutoff} = D/\lambda f$

#### 纳米结构的衍射效应

现代材料科学创造了许多具有特殊衍射性质的纳米结构：

1. **光子晶体**：
   - 周期性介电结构
   - 布拉格衍射条件：$2d\sin\theta = m\lambda$
   - 带隙效应：某些频率完全反射
   
   反射率：
   $$R(\omega) = \tanh^2\left(\frac{L}{\xi(\omega)}\right)$$
   
   其中 $\xi$ 是衰减长度。

2. **超表面（Metasurfaces）**：
   - 亚波长结构阵列
   - 相位调控：$\phi(x,y) = \arg[t(x,y)]$
   - 广义斯涅尔定律：
   $$n_t\sin\theta_t - n_i\sin\theta_i = \frac{\lambda_0}{2\pi}\frac{d\phi}{dx}$$

3. **等离激元结构**：
   - 金属纳米粒子的共振
   - 场增强因子：$|E/E_0|^2 \sim 10^3-10^6$
   - 颜色选择性：依赖于粒子形状和尺寸

#### 生物光学中的衍射

自然界中许多生物结构产生衍射效应：

1. **蝴蝶翅膀**：
   - 多层薄膜干涉
   - 光子晶体结构
   - 角度相关的虹彩色
   
   反射光谱：
   $$R(\lambda, \theta) = \left| \sum_n r_n \exp(2ik_z d_n) \right|^2$$

2. **鸟类羽毛**：
   - 准有序纳米结构
   - 相干散射产生结构色
   - 偏振相关效应

3. **昆虫复眼**：
   - 微透镜阵列
   - 抗反射纳米结构
   - 宽视角衍射效应

#### 工程应用中的衍射设计

利用衍射原理的工程设计：

1. **衍射光学元件（DOE）**：
   - 菲涅尔透镜：台阶化相位轮廓
   - 光束整形器：定制强度分布
   - 分束器：多级衍射光栅
   
   设计方程（IFTA算法）：
   $$\phi_{n+1} = \arg\left\{\mathcal{F}^{-1}\{A_{target} \exp(i\arg[\mathcal{F}\{A_{input}e^{i\phi_n}\}])\}\right\}$$

2. **全息光学元件（HOE）**：
   - 体全息光栅
   - 波长选择性
   - 角度多路复用
   
   耦合波方程：
   $$\frac{dR}{dz} = -i\kappa S$$
   $$\frac{dS}{dz} = -i\kappa^* R + i\Delta\beta S$$

3. **增强现实显示**：
   - 波导耦合光栅
   - 衍射扩瞳技术
   - 彩虹效应补偿

#### 衍射伪影与校正

识别和处理渲染中的衍射伪影：

1. **数值衍射伪影**：
   - 采样不足：摩尔纹
   - 边界截断：吉布斯现象
   - 相位展开错误：跳变伪影
   
   抗锯齿策略：
   $$u_{filtered} = \int u(\mathbf{x}) h(\mathbf{x} - \mathbf{x}_0) d\mathbf{x}$$

2. **物理伪影**：
   - 鬼像：多重反射衍射
   - 杂散光：边缘衍射
   - 色差：波长相关衍射
   
   校正方法：
   - 预补偿设计
   - 自适应光学
   - 计算后处理

3. **感知优化**：
   - 利用人眼MTF
   - 考虑观察条件
   - 实时质量调整

#### 未来展望

衍射理论在计算机图形学中的发展方向：

1. **机器学习加速**：
   - 神经网络预测衍射图案
   - 数据驱动的简化模型
   - 实时逆向设计

2. **量子渲染**：
   - 考虑光子统计
   - 纠缠光子效应
   - 量子相干性

3. **多物理场耦合**：
   - 热-光耦合（热透镜效应）
   - 力-光耦合（光镊模拟）
   - 电磁-机械耦合（MEMS）

## 本章小结

本章深入探讨了衍射理论及其在计算机图形学中的应用。我们从几何光学的局限性出发，逐步建立了能够准确描述衍射现象的数学框架。

### 关键概念

1. **几何绕射理论（GTD）**：
   - 将射线光学扩展到包含衍射射线
   - 基于广义费马原理和Keller锥
   - 提供高频渐近解

2. **一致性衍射理论（UTD）**：
   - 修正GTD在阴影边界的奇点
   - 引入过渡函数确保场连续性
   - 适用于工程计算

3. **物理衍射理论（PTD）**：
   - 分离物理光学和边缘贡献
   - Ufimtsev的边缘波理论
   - 中频段的有效方法

4. **复射线理论**：
   - 处理焦散和阴影区域
   - 高斯束方法
   - 解析延拓技术

### 核心公式

**GTD衍射系数**：
$$D(\phi, \phi'; \beta_0; \alpha) = -\frac{e^{i\pi/4}}{2\sqrt{2\pi k_0}} \frac{1}{\sin\beta_0} \sum_{i=1,2} \cot\left(\frac{\pi \pm (\phi-\phi')}{2N}\right)$$

**UTD过渡函数**：
$$F(\xi) = 2i\sqrt{\pi} \xi e^{i\xi^2} \int_\xi^\infty e^{-it^2} dt$$

**衍射场**：
$$u^d(P) = u^i(Q_D) \cdot D \cdot \sqrt{\frac{\rho}{s(\rho + s)}} \cdot e^{ik_0 s}$$

**菲涅尔数判据**：
$$F = \frac{a^2}{L\lambda}$$

### 计算考虑

1. **精度与效率权衡**：
   - 完整波动方程：最精确但计算密集
   - GTD/UTD：高频近似，计算高效
   - 屏幕空间方法：实时但精度有限

2. **数值稳定性**：
   - 过渡区域的特殊处理
   - 避免奇点和不连续性
   - 相位展开的正确处理

3. **适用范围**：
   - 频率限制：$k_0 a \gg 1$
   - 几何限制：局部凸性
   - 观察距离：远场条件

### 实际应用指南

1. **何时考虑衍射**：
   - 锐利边缘和小孔径
   - 波长尺度的特征
   - 光学系统模拟
   - 特殊材料效果

2. **实现策略**：
   - 射线追踪扩展：添加衍射射线
   - 后处理近似：屏幕空间卷积
   - 预计算方法：查找表和纹理
   - GPU加速：并行化和层次结构

3. **验证方法**：
   - 与解析解比较
   - 数值收敛性测试
   - 感知度量评估
   - 用户研究验证

## 练习题

### 基础题

**练习19.1** 考虑一个半无限大完美导体平面（$y = 0, x > 0$）。平面波以角度 $\phi^i = 60°$ 入射，波长 $\lambda = 500$ nm。
a) 使用GTD计算 $\phi^d = 120°$ 方向上的衍射系数。
b) 确定阴影边界和反射边界的位置。
c) 解释为什么标准GTD在这些边界失效。

*提示*：使用半平面的衍射系数公式，注意 $\alpha = 2\pi$。

<details>
<summary>答案</summary>

a) 半平面衍射系数（$\alpha = 2\pi$，$N = 1$）：
$$D = -\frac{e^{i\pi/4}}{2\sqrt{2\pi k_0}} \frac{1}{\sin\beta_0} \left[ \frac{1}{\cos\frac{\phi-\phi'}{2}} + \frac{1}{\cos\frac{\phi+\phi'}{2}} \right]$$

代入 $\phi^i = 60°$，$\phi^d = 120°$，$\beta_0 = 90°$（垂直入射）：
$$D = -\frac{e^{i\pi/4}}{2\sqrt{2\pi k_0}} \left[ \frac{1}{\cos 30°} + \frac{1}{\cos 90°} \right]$$

第二项发散，表明需要UTD修正。

b) 阴影边界：$\phi_{SB} = 180° - \phi^i = 120°$
   反射边界：$\phi_{RB} = 180° + \phi^i = 240°$

c) GTD在边界失效因为衍射系数包含 $1/\cos 90° = \infty$ 的奇点。物理上，这对应于衍射射线与几何光学射线重合的情况。

</details>

**练习19.2** 推导楔形边缘（楔角 $\alpha$）的Keller锥条件。证明入射射线和衍射射线与边缘的夹角相等。

*提示*：从广义费马原理出发，考虑通过边缘点的路径。

<details>
<summary>答案</summary>

设源点 $S$，边缘点 $Q$，观察点 $P$。路径长度：
$$L = |SQ| + |QP| = \sqrt{(x_s-x_q)^2 + (y_s-y_q)^2 + (z_s-z_q)^2} + \sqrt{(x_p-x_q)^2 + (y_p-y_q)^2 + (z_p-z_q)^2}$$

对边缘参数 $t$ 求导并令其为零：
$$\frac{\partial L}{\partial t} = -\hat{\mathbf{s}}^i \cdot \hat{\mathbf{t}} + \hat{\mathbf{s}}^d \cdot \hat{\mathbf{t}} = 0$$

因此：$\hat{\mathbf{s}}^i \cdot \hat{\mathbf{t}} = \hat{\mathbf{s}}^d \cdot \hat{\mathbf{t}} = \cos\beta_0$

这证明了入射角等于衍射角（相对于边缘）。

</details>

**练习19.3** 比较UTD和PTD方法计算直边缘衍射的计算复杂度。考虑 $N$ 个观察点和 $M$ 个边缘的情况。

*提示*：分析每种方法需要的基本运算。

<details>
<summary>答案</summary>

UTD复杂度：
- 对每个观察点和边缘对：$O(1)$ 计算衍射系数
- 总复杂度：$O(NM)$

PTD复杂度：
- 物理光学积分：$O(N \cdot S)$，其中 $S$ 是表面采样点数
- 边缘修正：$O(NM)$
- 总复杂度：$O(N(S + M))$

当 $S \gg M$ 时，PTD计算量显著大于UTD。

</details>

### 挑战题

**练习19.4** 设计一个衍射光栅，在正入射下将白光（400-700 nm）分散到 ±30° 范围内。
a) 确定光栅周期。
b) 计算各波长的衍射效率。
c) 讨论高阶衍射的影响。

*提示*：使用光栅方程 $d(\sin\theta_m - \sin\theta_i) = m\lambda$。

<details>
<summary>答案</summary>

a) 光栅方程（正入射，$\theta_i = 0$）：
$$d \sin\theta_m = m\lambda$$

第一级（$m = 1$），最大偏转角 30°，最长波长 700 nm：
$$d = \frac{\lambda_{max}}{\sin 30°} = \frac{700 \text{ nm}}{0.5} = 1400 \text{ nm}$$

b) 衍射效率（标量理论）：
$$\eta_m = \left| \text{sinc}\left(\frac{\pi a}{\lambda}\sin\theta_m\right) \right|^2$$

其中 $a$ 是光栅线宽。对于 $a = d/2$：
- 400 nm：$\theta_1 = 16.6°$，$\eta_1 \approx 0.41$
- 550 nm：$\theta_1 = 23.1°$，$\eta_1 \approx 0.41$
- 700 nm：$\theta_1 = 30.0°$，$\eta_1 \approx 0.41$

c) 二阶衍射（$m = 2$）：
- 仅对 $\lambda < 700$ nm 存在
- 造成光谱重叠
- 可通过闪耀光栅优化抑制

</details>

**练习19.5** 推导复射线（高斯束）通过楔形边缘的衍射。假设束腰位于边缘，初始束腰半径 $w_0$。

*提示*：使用复源点方法和解析延拓。

<details>
<summary>答案</summary>

高斯束可视为从复源点发出的球面波：
$$\mathbf{r}_s = \mathbf{r}_0 + i\frac{w_0^2}{2k_0}\hat{\mathbf{s}}^i$$

衍射场通过解析延拓获得：
$$u^d_{GB} = D(\phi_c, \phi_c'; \beta_{0c}) \frac{e^{ik_0 R_c}}{\sqrt{R_c}}$$

其中复角度：
$$\cos\beta_{0c} = \frac{\hat{\mathbf{s}}^i_c \cdot \hat{\mathbf{t}}}{|\hat{\mathbf{s}}^i_c|}$$

复距离：
$$R_c = |\mathbf{r} - \mathbf{r}_s|$$

结果是衍射高斯束，具有修正的传播参数和额外的Gouy相位。

</details>

**练习19.6** 分析计算机图形学中相机光圈叶片产生的星芒效应。考虑 $N$ 片叶片的多边形光圈。

*提示*：将问题分解为多个直边的衍射，考虑对称性。

<details>
<summary>答案</summary>

$N$ 边正多边形光圈的衍射图案：

1. 每条边产生垂直于边的衍射条纹
2. 总场是所有边贡献的相干叠加
3. 对称性导致 $2N$ 条射线（$N$ 为奇数）或 $N$ 条射线（$N$ 为偶数）

远场强度分布：
$$I(\theta) \propto \left| \sum_{n=0}^{N-1} \text{sinc}\left(\frac{ka}{2}\sin(\theta - 2\pi n/N)\right) e^{i\phi_n} \right|^2$$

其中 $a$ 是边长，$\phi_n$ 是相对相位。

星芒亮度随孔径减小而增强，因为衍射效应 $\propto 1/a$。

</details>

### 思考题

**练习19.7** 讨论为什么隐身飞机的设计需要考虑PTD理论。哪些几何特征最容易产生强衍射回波？

*提示*：考虑雷达散射截面（RCS）和边缘贡献。

<details>
<summary>答案</summary>

PTD理论对隐身设计至关重要因为：

1. **边缘是主要散射源**：
   - 机翼前缘和后缘
   - 舱门和面板接缝
   - 进气口边缘

2. **设计原则**：
   - 避免垂直边缘对准威胁方向
   - 使用锯齿边缘分散衍射能量
   - 边缘处理（雷达吸波材料）

3. **PTD预测**：
   - 精确计算边缘散射贡献
   - 优化边缘几何降低RCS
   - 评估涂层效果

4. **关键几何特征**：
   - 直角和锐角产生强后向散射
   - 平行边缘造成相干增强
   - 腔体开口形成谐振增强

</details>

**练习19.8** 设计一个实验来验证GTD预测的衍射场。需要什么设备？如何区分衍射贡献和其他散射机制？

*提示*：考虑测量精度和隔离不同贡献的方法。

<details>
<summary>答案</summary>

实验设计：

1. **设备需求**：
   - 相干光源（激光器）
   - 精密定位系统（分辨率 < λ/10）
   - 高灵敏度探测器
   - 消声/暗室环境

2. **测试几何**：
   - 标准楔形或半平面
   - 已知材料参数
   - 可变楔角验证理论

3. **测量方案**：
   - 角度扫描：固定距离，变化观察角
   - 距离扫描：固定角度，变化距离
   - 频率扫描：验证 $k_0^{-1/2}$ 依赖性

4. **数据分析**：
   - 背景扣除（暗场测量）
   - 几何光学贡献分离（遮挡法）
   - 与理论预测对比（最小二乘拟合）

5. **误差源**：
   - 有限光束尺寸
   - 多重散射
   - 测量噪声
   - 对准误差

</details>
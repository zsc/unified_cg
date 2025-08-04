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

## 19.2 物理绕射理论（PTD/UTD）

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

### 19.4.1 复射线的概念

### 19.4.2 衍射系数的解析延拓

### 19.4.3 因果性与解析性

### 19.4.4 数值计算方法

## 19.5 计算机图形学中的衍射效应

### 19.5.1 衍射在渲染中的重要性

### 19.5.2 基于GTD的射线追踪扩展

### 19.5.3 实时衍射近似

### 19.5.4 与波动光学方法的比较

## 本章小结

## 练习题

## 常见陷阱与错误

## 最佳实践检查清单
# 第17章：波前整形与自适应光学

自适应光学和波前整形代表了计算控制与物理光学相结合的前沿领域。在本章中，我们将探讨如何通过对光学波前的受控操纵，实现穿透散射介质的聚焦、像差校正以及光传输的优化。我们将使用泽尼克多项式建立数学基础，研究空间光调制器，并将这些概念与计算机图形学中的自适应采样策略联系起来。

从像差波前到衍射极限成像的历程，与从朴素蒙特卡洛采样到渲染中复杂自适应技术的发展异曲同工。正如自适应光学实时校正大气畸变一样，现代渲染算法动态调整采样模式以最小化方差。这种深层联系超越了类比——波前优化和重要性采样的数学框架共享植根于信号重建和优化理论的基本原理。

## 17.1 泽尼克多项式与波前描述

### 17.1.1 数学基础

泽尼克多项式在单位圆盘上形成一个完整的正交基，非常适合描述圆形瞳孔中的波前像差。这些多项式以弗里茨·泽尼克命名，他因发明相衬显微镜而获得1953年诺贝尔奖，它们具有独特的性质，使其在光学中不可或缺。这些多项式定义为：

$$Z_n^m(\rho, \theta) = R_n^{|m|}(\rho) \times \begin{cases}
\cos(m\theta) & m \geq 0 \\
\sin(|m|\theta) & m < 0
\end{cases}$$

其中径向多项式为：

$$R_n^m(\rho) = \sum_{k=0}^{(n-m)/2} \frac{(-1)^k (n-k)!}{k! \left(\frac{n+m}{2}-k\right)! \left(\frac{n-m}{2}-k\right)!} \rho^{n-2k}$$

索引遵循以下约束：
- n ≥ 0 (径向阶数)
- -n ≤ m ≤ n (方位角频率)
- n - |m| 为偶数

最后一个约束确保 $R_n^m(\rho)$ 只包含 $\rho$ 的偶数或奇数次幂，从而保持多项式的奇偶性。径向多项式也可以用雅可比多项式表示：

$$R_n^m(\rho) = (-1)^{(n-m)/2} \rho^m P_{(n-m)/2}^{(m,0)}(2\rho^2 - 1)$$

其中 $P_k^{(\alpha,\beta)}$ 是雅可比多项式。这种联系揭示了更深层的数学结构，并能够使用递推关系进行高效计算。

**递推关系：**
雅可比多项式表示允许通过三项递推进行高效计算：

$$P_k^{(\alpha,\beta)}(x) = \frac{(2k+\alpha+\beta-1)[(2k+\alpha+\beta)(2k+\alpha+\beta-2)x + \alpha^2 - \beta^2]}{2k(k+\alpha+\beta)(2k+\alpha+\beta-2)} P_{k-1}^{(\alpha,\beta)}(x) - \frac{2(k+\alpha-1)(k+\beta-1)(2k+\alpha+\beta)}{2k(k+\alpha+\beta)(2k+\alpha+\beta-2)} P_{k-2}^{(\alpha,\beta)}(x)$$

对于 $\alpha = m, \beta = 0$ 的泽尼克径向多项式：

$$R_n^m(\rho) = \frac{4(n-1)\rho^2 - 2(n+m-2)}{n-m} R_{n-2}^m(\rho) - \frac{(n+m-2)}{n-m} R_{n-4}^m(\rho)$$

这种递推关系将高阶多项式的计算复杂度从 $O(n^2)$ 大幅降低到 $O(n)$。

**生成函数：**
泽尼克径向多项式的生成函数提供了另一种计算途径：

$$\sum_{n=m}^{\infty} R_n^m(\rho) t^{(n-m)/2} = \frac{\rho^m}{(1-t)^{m+1}} \cdot _2F_1\left(\frac{m+1}{2}, \frac{m+2}{2}; m+1; \frac{4\rho^2 t}{(1-t)^2}\right)$$

其中 $_2F_1$ 是超几何函数。这种表示与量子力学角动量本征函数和球谐函数相关联。

为了归一化，我们包含因子：
$$N_n^m = \sqrt{\frac{2(n+1)}{1 + \delta_{m0}}}$$

使得归一化多项式为：
$$\tilde{Z}_n^m(\rho, \theta) = N_n^m Z_n^m(\rho, \theta)$$

归一化确保了在单位圆盘上积分时每个模式的单位方差，这对于比较不同阶数的像差强度至关重要。

**与其他多项式系统的联系：**
泽尼克多项式与其他正交多项式族相关：
- **勒让德多项式**：对于 $m = 0$，有 $R_n^0(\rho) = P_{n/2}(2\rho^2 - 1)$
- **切比雪夫多项式**：通过雅可比多项式的特例连接
- **贝塞尔函数**：大 $n$ 时的渐近行为遵循贝塞尔函数振荡

这些联系使得可以利用现有的数学工具进行分析和计算。

### 17.1.2 正交性

在单位圆盘上的正交性关系构成了波前分析的数学基础：

$$\int_0^{2\pi} \int_0^1 Z_n^m(\rho, \theta) Z_{n'}^{m'}(\rho, \theta) \rho d\rho d\theta = \frac{\pi}{2n+2} \delta_{nn'} \delta_{mm'}$$

这种正交性源于两个独立的积分，由于多项式结构而解耦：

**角度正交性：**
$$\int_0^{2\pi} \cos(m\theta)\cos(m'\theta) d\theta = \pi\delta_{mm'} \quad (m, m' > 0)$$
$$\int_0^{2\pi} \sin(m\theta)\sin(m'\theta) d\theta = \pi\delta_{mm'} \quad (m, m' > 0)$$
$$\int_0^{2\pi} \cos(m\theta)\sin(m'\theta) d\theta = 0$$

这些关系遵循三角函数在完整周期上的正交性。交叉项由于被积函数的奇对称性而消失。

**径向正交性：**
$$\int_0^1 R_n^m(\rho) R_{n'}^m(\rho) \rho d\rho = \frac{1}{2(n+1)} \delta_{nn'}$$

积分中的权重函数 $\rho$ 源于极坐标中的雅可比行列式，并确保在圆形域上的适当正交性。这个权重至关重要——没有它，径向多项式将不正交。$R_n^m$ 的特定形式经过精心构造，以实现相对于此测度的正交性。

**径向正交性证明：**
径向正交性可以使用雅可比多项式的罗德里格斯公式证明：

$$P_n^{(\alpha,\beta)}(x) = \frac{(-1)^n}{2^n n!}(1-x)^{-\alpha}(1+x)^{-\beta} \frac{d^n}{dx^n}[(1-x)^{\alpha+n}(1+x)^{\beta+n}]$$

转换为泽尼克径向多项式域，其中 $x = 2\rho^2 - 1$：

$$\int_0^1 R_n^m(\rho) R_{n'}^m(\rho) \rho d\rho = \frac{1}{2}\int_{-1}^1 P_{(n-m)/2}^{(m,0)}(x) P_{(n'-m)/2}^{(m,0)}(x) (1+x)^{m/2} dx$$

雅可比多项式与权重 $(1-x)^\alpha(1+x)^\beta$ 的正交性产生了所需的结果。

**广义正交性关系：**
对于内半径为 $\epsilon$ 的环形瞳孔：

$$\int_0^{2\pi} \int_\epsilon^1 Z_n^m(\rho, \theta) Z_{n'}^{m'}(\rho, \theta) \rho d\rho d\theta = h_{nm}(\epsilon) \delta_{nn'} \delta_{mm'}$$

其中归一化因子 $h_{nm}(\epsilon)$ 取决于遮挡比。这种修改对于具有中心遮挡的望远镜系统至关重要。

**离散正交性：**
对于 $N$ 个点 $(\rho_i, \theta_i)$ 的离散网格上的数值实现：

$$\sum_{i=1}^N Z_n^m(\rho_i, \theta_i) Z_{n'}^{m'}(\rho_i, \theta_i) w_i \approx \frac{\pi}{2n+2} \delta_{nn'} \delta_{mm'}$$

其中 $w_i$ 是求积权重。圆盘上的求积规则提供了最佳点分布：
- 高斯-雅可比径向节点：$\rho_j$ 是 $P_n^{(1,0)}(2\rho^2-1)$ 的根
- 等间距角度节点：$\theta_k = 2\pi k/M$

**完备性与帕塞瓦尔恒等式：**
单位圆盘上的任何平方可积函数都可以用泽尼克多项式展开，帕塞瓦尔恒等式保证了能量守恒：

$$\int_0^{2\pi} \int_0^1 |W(\rho, \theta)|^2 \rho d\rho d\theta = \sum_{n=0}^{\infty} \sum_{m=-n}^n |a_n^m|^2 \frac{\pi}{2n+2}$$

这种关系使得可以直接从泽尼克系数计算波前方差，这对于光学设计中的像差预算至关重要。

**贝塞尔不等式与收敛性：**
对于 $N$ 阶的任何有限截断：

$$\sum_{n=0}^{N} \sum_{m=-n}^n |a_n^m|^2 \frac{\pi}{2n+2} \leq \int_0^{2\pi} \int_0^1 |W(\rho, \theta)|^2 \rho d\rho d\theta$$

等式仅在 $N \to \infty$ 的极限下成立，为截断展开提供了收敛准则。

### 17.1.3 波前展开

任何波前 $W(\rho, \theta)$ 都可以展开为：

$$W(\rho, \theta) = \sum_{n=0}^{\infty} \sum_{m=-n}^n a_n^m Z_n^m(\rho, \theta)$$

系数通过内积计算：

$$a_n^m = \frac{2n+2}{\pi \epsilon_m} \int_0^{2\pi} \int_0^1 W(\rho, \theta) Z_n^m(\rho, \theta) \rho d\rho d\theta$$

其中 $\epsilon_m = 2$ 对于 $m = 0$，否则 $\epsilon_m = 1$。这个因子解释了 $m = 0$（旋转对称）项的不同归一化。

**矩阵公式：**
对于在点 $(\rho_i, \theta_i)$ 采样的波前，展开变为线性系统：

$$\mathbf{w} = \mathbf{Z} \mathbf{a} + \boldsymbol{\epsilon}$$

其中：
- $\mathbf{w} \in \mathbb{R}^M$：波前测量值
- $\mathbf{Z} \in \mathbb{R}^{M \times J}$：泽尼克多项式矩阵，其中 $Z_{ij} = Z_j(\rho_i, \theta_i)$
- $\mathbf{a} \in \mathbb{R}^J$：系数向量
- $\boldsymbol{\epsilon}$：测量噪声

**高效系数计算：**
对于离散网格上的测量波前数据，系数可以使用最小二乘法计算：

$$\mathbf{a} = (\mathbf{Z}^T \mathbf{Z})^{-1} \mathbf{Z}^T \mathbf{w}$$

其中 $\mathbf{Z}$ 是测量点处泽尼克多项式值的矩阵，$\mathbf{w}$ 是波前测量值的向量。对于具有适当采样的规则网格，$\mathbf{Z}^T \mathbf{Z}$ 变得近似对角，从而简化了求逆。

**加权最小二乘法：**
当测量不确定性在空间上变化时，使用加权最小二乘法：

$$\mathbf{a} = (\mathbf{Z}^T \mathbf{W} \mathbf{Z})^{-1} \mathbf{Z}^T \mathbf{W} \mathbf{w}$$

其中 $\mathbf{W} = \text{diag}(1/\sigma_i^2)$ 包含逆方差。这种方法自然处理：
- 瞳孔上可变的信噪比
- 缺失数据（将权重设置为零）
- 非均匀瞳孔照明

**奇异值分解方法：**
对于病态问题，使用 SVD 正则化：

$$\mathbf{Z} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T$$

$$\mathbf{a} = \sum_{i=1}^{r} \frac{\mathbf{u}_i^T \mathbf{w}}{\sigma_i} \mathbf{v}_i$$

其中 $r$ 是由以下确定的有效秩：
$$\sigma_r / \sigma_1 > \epsilon_{machine}$$

这可以防止在采样不良的模式中噪声放大。

在实践中，我们将展开截断到某个最大阶数 $N$：
$$W(\rho, \theta) \approx \sum_{n=0}^{N} \sum_{m=-n}^n a_n^m Z_n^m(\rho, \theta)$$

直到 $N$ 阶的项数为：
$$J = \frac{(N+1)(N+2)}{2}$$

例如：
- $N = 3$：$J = 10$ 项（直到彗差）
- $N = 4$：$J = 15$ 项（包括球差）
- $N = 8$：$J = 45$ 项（高阶像差）
- $N = 20$：$J = 231$ 项（极端自适应光学）

**最佳截断阶数：**
最佳截断平衡了拟合误差与噪声放大：

$$N_{opt} = \arg\min_N \left[ \sigma_{fit}^2(N) + \sigma_{noise}^2(N) \right]$$

其中：
- $\sigma_{fit}^2(N) \propto N^{-\alpha}$（拟合误差）
- $\sigma_{noise}^2(N) \propto N$（噪声传播）

对于具有测量噪声的科尔莫哥洛夫湍流：
$$N_{opt} \approx 2.2(D/r_0)^{0.6}(\text{SNR})^{0.3}$$

**截断误差分析：**
在 $N$ 阶截断的 RMS 误差遵循：

$$\sigma_{truncation}^2 = \sum_{n=N+1}^{\infty} \sum_{m=-n}^n |a_n^m|^2$$

对于大气湍流，此误差随 $N^{-\alpha/2}$ 减小，其中 $\alpha \approx 11/3$，为根据所需精度选择截断阶数提供了指导。

**模态方差标度：**
在科尔莫哥洛夫湍流中，每个模式的预期方差：

$$\langle |a_n^m|^2 \rangle = A_n (D/r_0)^{5/3}$$

其中对于 $n \gg 1$，有 $A_n \propto (n+1)^{-11/3}$。这种幂律使得可以预测高阶贡献。

### 17.1.4 索引转换

泽尼克多项式存在几种索引方案，每种方案在不同应用中都有其特定优势：

**Noll 符号（单索引 j）：**
$$j = \frac{n(n+1)}{2} + |m| + \begin{cases}
0 & \text{if } m \leq 0 \text{ and } n \bmod 4 \in \{0,1\} \\
0 & \text{if } m > 0 \text{ and } n \bmod 4 \in \{2,3\} \\
1 & \text{otherwise}
\end{cases}$$

这种排序在正弦项和余弦项之间交替，以保持像差类型的逻辑进展。前几个 Noll 索引对应于：
- $j = 1$：活塞（$Z_0^0$）
- $j = 2,3$：倾斜（$Z_1^1, Z_1^{-1}$）
- $j = 4,5,6$：散光和离焦（$Z_2^{-2}, Z_2^0, Z_2^2$）

**OSA/ANSI 标准：**
$$j = \frac{n(n+2) + m}{2}$$

这个更简单的公式创建了一个更容易计算的单调排序，但对于像差分析来说不太直观。

**Wyant 排序：**
按多项式阶数分组，然后按方位角频率分组。这种排序有助于：
- 多项式复杂度的系统性增加
- 径向阶数的清晰分离
- 简化的误差传播分析

**Fringe/亚利桑那大学符号：**
按空间频率递增排序，对于干涉分析很有用，其中较高频率对应于更精细的条纹图案。

排序的选择影响系数解释，但不影响数学性质。排序之间的转换矩阵是稀疏的，可以预先计算以提高效率。

### 17.1.5 常见像差

前几个泽尼克项对应于常见的光学像差，每种像差都有独特的物理起源和视觉效果：

**低阶项（n ≤ 2）：**
- $Z_0^0 = 1$：**活塞**（恒定相位偏移）
  - 对图像质量无影响，仅影响绝对相位
  - 对于干涉测量和相干成像很重要
  
- $Z_1^{-1} = 2\rho\sin\theta$：**垂直倾斜**（y-倾斜）
- $Z_1^1 = 2\rho\cos\theta$：**水平倾斜**（x-倾斜）
  - 导致图像位移而不模糊
  - 源于光学元件中的楔形或未对准
  
- $Z_2^{-2} = \sqrt{6}\rho^2\sin(2\theta)$：**斜散光**（45°）
- $Z_2^2 = \sqrt{6}\rho^2\cos(2\theta)$：**垂直散光**（0°/90°）
  - 在不同距离处创建正交线焦点
  - 常见于柱面光学元件和眼睛像差
  
- $Z_2^0 = \sqrt{3}(2\rho^2 - 1)$：**离焦**
  - 对称模糊，改变最佳焦点位置
  - 等效于像平面的纵向移动

**三阶项（n = 3）：**
- $Z_3^{-3} = \sqrt{8}\rho^3\sin(3\theta)$：**三叶像差**
- $Z_3^3 = \sqrt{8}\rho^3\cos(3\theta)$：**三叶像差**
  - 三重对称畸变
  - 通常来自受力光学支架
  
- $Z_3^{-1} = \sqrt{8}(3\rho^3 - 2\rho)\sin\theta$：**垂直彗差**
- $Z_3^1 = \sqrt{8}(3\rho^3 - 2\rho)\cos\theta$：**水平彗差**
  - 彗星状模糊，随视场角增大
  - 离轴成像的基本限制

**四阶项（n = 4）：**
- $Z_4^{-4} = \sqrt{10}\rho^4\sin(4\theta)$：**四叶像差**
- $Z_4^4 = \sqrt{10}\rho^4\cos(4\theta)$：**四叶像差**
  - 四重对称，方形畸变
  
- $Z_4^{-2} = \sqrt{10}(4\rho^4 - 3\rho^2)\sin(2\theta)$：**次级散光**
- $Z_4^2 = \sqrt{10}(4\rho^4 - 3\rho^2)\cos(2\theta)$：**次级散光**
  - 高阶散光效应
  
- $Z_4^0 = \sqrt{5}(6\rho^4 - 6\rho^2 + 1)$：**初级球差**
  - 径向对称模糊，随孔径增大
  - 球面表面的基本限制

**视觉影响：**
每种像差都会产生特征性的点扩散函数（PSF）畸变：
- 彗差：彗星状 PSF，尾部径向指向
- 散光：椭圆形 PSF，随焦点位置旋转
- 球差：圆形光晕，核心明亮
- 三叶/四叶像差：多瓣 PSF 图案

### 17.1.6 物理解释与赛德尔连接

泽尼克多项式通过坐标变换与经典赛德尔像差相关。对于视场角 $\alpha$ 和瞳孔坐标 $(\rho, \theta)$ 处的点：

**赛德尔到泽尼克映射：**
- 球差：$W_{040} = a_4^0 Z_4^0$
- 彗差：$W_{131} = a_3^1 Z_3^1 \cos\alpha + a_3^{-1} Z_3^{-1} \sin\alpha$
- 散光：$W_{222} = a_2^2 Z_2^2 \cos(2\alpha) + a_2^{-2} Z_2^{-2} \sin(2\alpha)$
- 场曲：$W_{220} = a_2^0 Z_2^0$
- 畸变：$W_{311} = a_1^1 Z_1^1 \cos\alpha + a_1^{-1} Z_1^{-1} \sin\alpha$

经典赛德尔像差系数来自光线追迹，而泽尼克系数来自波前拟合。它们之间的变换涉及：

$$W_{Seidel}(h, \rho, \theta, \phi) = \sum_{i,j,k,l} S_{ijkl} h^i \rho^j \cos^k(\theta-\phi)$$

其中 $h$ 是归一化视场高度，索引满足 $i + j = 2(k + l)$。这种幂级数展开可以通过三角恒等式改写为泽尼克形式。

**视场依赖性：**
与仅依赖瞳孔的泽尼克多项式不同，全视场像差需要额外的视场坐标：
$$W(h_x, h_y, \rho, \theta) = \sum_{p,q,n,m} b_{pqnm} h_x^p h_y^q Z_n^m(\rho, \theta)$$

这种视场相关的泽尼克展开使得：
- 宽视场系统的节点像差理论
- 测量像差视场的有效存储
- 视场平均性能的直接优化

每个像差的波前方差贡献：
$$\sigma_n^2 = \sum_{m=-n}^n (a_n^m)^2$$

这种分解揭示了哪些阶数在像差预算中占主导地位，从而指导校正策略。

### 17.1.7 RMS 波前误差

均方根波前误差在泽尼克基中优雅地表达：

$$\sigma_{RMS} = \sqrt{\sum_{n,m} (a_n^m)^2}$$

这种正交性使得泽尼克系数成为优化算法的理想选择。方差可以按阶数分解：

$$\sigma_{RMS}^2 = \sum_{n=0}^{\infty} \sigma_n^2 = \sum_{n=0}^{\infty} \sum_{m=-n}^n (a_n^m)^2$$

**物理意义：**
RMS 波前误差直接与光学性能指标相关：
- 斯特雷尔比（峰值强度）：对于 $\sigma_{RMS} < 1$ 弧度，有 $S \approx \exp(-\sigma_{RMS}^2)$
- 包围能量：给定半径内的光分数随 $\sigma_{RMS}$ 减小
- 分辨率：有效 PSF 宽度近似随 $(1 + \sigma_{RMS}^2)^{1/2}$ 增加

对于遵循科尔莫哥洛夫统计的大气湍流：
$$\langle (a_n^m)^2 \rangle \propto (n + 1)^{-\alpha}$$

其中对于充分发展的湍流，$\alpha \approx 11/3$。这种幂律衰减证明了在有限阶数处截断展开的合理性。

**弗里德参数连接：**
视场受限分辨率由弗里德参数 $r_0$ 表征：
$$\sigma_{RMS}^2 = 1.03(D/r_0)^{5/3}$$

其中 $D$ 是望远镜直径。对于 $D \gg r_0$，波前误差主要由低阶模式主导，特别是倾斜模式，它包含总方差的约 87%。

**时间演化：**
大气像差随特征时间尺度演化：
$$\tau_n \approx \frac{r_0}{v_{wind}} \cdot n^{-3/5}$$

高阶像差变化更快，需要更快的校正循环。这种时间谱指导自适应光学控制系统设计。

### 17.1.8 斯特雷尔比与性能指标

对于小像差，斯特雷尔比（峰值强度比）近似为：

$$S \approx \exp(-\sigma_{RMS}^2) \approx 1 - \sigma_{RMS}^2$$

其中 $\sigma_{RMS}$ 以弧度表示。这提供了波前质量与成像性能之间的直接联系。

**扩展马雷夏尔近似：**
$$S \approx \exp\left[-\sigma_{RMS}^2 + \frac{1}{2}\sigma_{RMS}^4(\kappa_4 - 1)\right]$$

其中 $\kappa_4$ 是相位分布的归一化四阶矩。对于高斯统计，$\kappa_4 = 3$。非高斯相位屏（例如，来自强闪烁）需要此校正。

**高阶统计：**
除了 RMS，相位结构函数提供了额外的洞察：
$$D_\phi(\mathbf{r}) = \langle[W(\mathbf{x} + \mathbf{r}) - W(\mathbf{x})]^2\rangle$$

对于科尔莫哥洛夫湍流：
$$D_\phi(r) = 6.88(r/r_0)^{5/3}$$

此结构函数确定像差的相关长度，并指导可变形镜中的执行器间距。

**衍射极限判据：**
- 瑞利判据：λ/4 峰谷值 → $S \approx 0.80$
- 马雷夏尔判据：λ/14 RMS → $S \approx 0.80$
- 实际限制：$\sigma_{RMS} < \lambda/20$ → $S > 0.94$

具有像差的点扩散函数（PSF）：
$$\text{PSF}(x,y) = \left|\mathcal{F}\left\{P(\xi,\eta)\exp\left[i\frac{2\pi}{\lambda}W(\xi,\eta)\right]\right\}\right|^2$$

其中 $P$ 是瞳孔函数，$W$ 是波前误差。

**部分校正效应：**
当仅校正到泽尼克阶数 $N$ 时：
$$S_{corrected} = S_{uncorrected} \cdot \exp\left(\sum_{n=0}^N \sum_{m=-n}^n |a_n^m|^2\right)$$

这表明每个校正模式都会带来指数级的改进，从而促使分层校正方案。

### 17.1.9 自适应光学性能

校正后的残余波前误差分解为独立的误差源：
$$\sigma_{residual}^2 = \sigma_{fitting}^2 + \sigma_{temporal}^2 + \sigma_{measurement}^2 + \sigma_{calibration}^2 + \sigma_{anisoplanatic}^2$$

每个术语都代表自适应光学系统中的基本限制：

**拟合误差（有限执行器）：**
$$\sigma_{fitting}^2 \approx 0.28\left(\frac{d}{r_0}\right)^{5/3}$$

其中 $d$ 是执行器间距，$r_0$ 是弗里德参数。此误差源于可变形镜无法再现高空间频率。系数 0.28 假设科尔莫哥洛夫湍流和连续面片镜。

**时间误差（有限带宽）：**
$$\sigma_{temporal}^2 \approx \left(\frac{f_G}{f_0}\right)^{5/3}$$

其中格林伍德频率为：
$$f_G = 0.43 \frac{v_{wind}}{r_0}$$

$f_0$ 是控制环路带宽。当大气变化速度超过校正系统时，此误差会累积。

**测量噪声传播：**
$$\sigma_{measurement}^2 = \left(\frac{\partial W}{\partial s}\right)^2 \sigma_s^2$$

其中 $s$ 代表传感器测量值。对于夏克-哈特曼传感器：
$$\sigma_{measurement}^2 \propto \frac{1}{N_{photons}} + \frac{\sigma_{read}^2}{N_{photons}^2}$$

第一项是光子噪声，第二项是读出噪声贡献。

**不等视场误差：**
$$\sigma_{anisoplanatic}^2 = \left(\frac{\theta}{\theta_0}\right)^{5/3}$$

其中 $\theta$ 是与导星的角分离，$\theta_0$ 是等视场角：
$$\theta_0 = 0.314 \frac{r_0}{H_{eff}}$$

其中 $H_{eff}$ 是有效湍流高度。

**误差预算优化：**
总误差最小化需要平衡：
- 更多执行器减少拟合误差但增加成本
- 更高带宽减少时间误差但放大噪声
- 更亮的导星减少测量误差但限制天空覆盖

此优化问题与蒙特卡洛渲染中的方差减少并行，其中样本分配必须平衡不同的误差源。

## 17.2 空间光调制器（SLM）原理

### 17.2.1 相位调制机制

空间光调制器能够实现光学相位、振幅或偏振的像素级控制。存在多种技术，每种技术都具有独特的特性：

**硅基液晶（LCoS）：**
最常见的 SLM 技术使用电控双折射。对于向列液晶，相位调制深度 $\delta$ 取决于：

$$\delta = \frac{2\pi}{\lambda} \Delta n \cdot d$$

其中 $\Delta n$ 是双折射变化，$d$ 是单元厚度。电压相关的折射率遵循：

$$n(V) = n_o + \frac{n_e - n_o}{1 + (V_{th}/V)^2}$$

其中 $n_o$ 和 $n_e$ 是寻常和非寻常折射率，$V_{th}$ 是阈值电压。

**指向矢取向模型：**
液晶指向矢角度 $\theta(z)$ 穿过单元厚度满足：
$$K\frac{d^2\theta}{dz^2} = \epsilon_0\Delta\epsilon E^2\sin\theta\cos\theta$$

其中 $K$ 是弹性常数。产生的相移：
$$\phi = \frac{2\pi}{\lambda}\int_0^d n_{eff}(\theta(z))dz$$

有效折射率为：
$$n_{eff}(\theta) = \frac{n_o n_e}{\sqrt{n_o^2\sin^2\theta + n_e^2\cos^2\theta}}$$

**响应动力学：**
响应时间尺度为：
$$\tau_{rise} = \frac{\gamma d^2}{K\epsilon_0\Delta\epsilon(V^2 - V_{th}^2)}$$
$$\tau_{fall} = \frac{\gamma d^2}{K\pi^2}$$

其中 $\gamma$ 是旋转粘度。注意不对称性：弛豫通常比激活慢。

**温度依赖性：**
双折射和响应时间都随温度变化：
$$\Delta n(T) = \Delta n_0(1 - T/T_c)^\beta$$
$$\gamma(T) = \gamma_0\exp(E_a/k_BT)$$

其中 $T_c$ 是清亮点，$\beta \approx 0.2$，$E_a$ 是活化能。

**数字微镜器件（DMD）：**
通过倾斜微镜进行二元振幅调制：
- 倾斜角度：通常 ±12°
- 切换时间：约 10 μs
- 衍射效率：约 88%（进入所需阶次）
- 对比度：>5000:1

**可变形镜：**
用于相位控制的连续表面变形：
- 执行器类型：压电、静电、磁性
- 行程：典型 1-10 μm
- 带宽：1-10 kHz
- 通过影响函数实现执行器间耦合

### 17.2.2 复振幅调制

纯相位 SLM 可以通过几种编码方案实现复振幅调制：

**1. 离轴计算机生成全息图（Lee 方法）：**
将复场 $U = A \exp(i\psi)$ 编码为：
$$\phi(x,y) = \arg[A(x,y)e^{i\psi(x,y)}] + 2\pi f_c x$$

其中 $f_c$ 是载波频率。第一衍射阶近似所需的复场。

载波频率必须满足：
$$f_c > f_{max} + \frac{W}{2\lambda z}$$

其中 $f_{max}$ 是 $U$ 中的最大空间频率，$W$ 是光束宽度。这确保了衍射阶的分离。

**效率分析：**
- 一阶效率：$\eta_1 \approx |\langle U \rangle|^2/\langle |U|^2 \rangle$
- 零阶泄漏：$\eta_0 = |\langle \exp(i\phi) \rangle|^2$
- 信噪比：$\text{SNR} \propto A^2/(1-A^2)$

**2. 双相位振幅编码：**
使用两个相距 $z$ 的相位掩模 $\phi_1$ 和 $\phi_2$：

$$U_{out} = \mathcal{F}^{-1}\{\mathcal{F}\{e^{i\phi_1}\} \cdot H(z) \cdot \mathcal{F}\{e^{i\phi_2}\}\}$$

其中 $H(z)$ 是菲涅尔传播核：
$$H(k_x, k_y; z) = \exp\left[iz\sqrt{k^2 - k_x^2 - k_y^2}\right]$$

相位掩模源自：
$$\phi_1 = \arg[U] + \arg[\mathcal{F}^{-1}\{|U|^{1/2}\}]$$
$$\phi_2 = -\arg[\mathcal{F}^{-1}\{|U|^{1/2}\}]$$

**3. 迭代傅里叶变换算法：**

**Gerchberg-Saxton (GS)：**
$$\phi_{n+1} = \arg[\mathcal{F}^{-1}\{|A_{target}|e^{i\arg[\mathcal{F}\{|A_{input}|e^{i\phi_n}\}]}\}]$$

**带反馈的加权 GS：**
$$\phi_{n+1} = \phi_n + \alpha \arg[\mathcal{F}^{-1}\{(A_{target} - A_n)e^{i\psi_n}\}]$$

其中 $\alpha \in (0,2)$ 是反馈强度。

**误差减少指标：**
$$\epsilon = \frac{\sum|I_{target} - I_{measured}|^2}{\sum I_{target}^2}$$

典型收敛：在 20-50 次迭代内 $\epsilon < 0.01$。

**4. 超像素方法：**
将 $N \times N$ 像素分组以编码振幅和相位：
- 中心像素：编码相位
- 边界像素：通过闪耀光栅控制振幅

通过衍射效率控制振幅：
$$A_{effective} = \eta(\theta_{blaze}) = \text{sinc}^2\left(\frac{N\pi\sin\theta_{blaze}}{\lambda/p}\right)$$

其中 $p$ 是像素间距。

### 17.2.3 衍射效率与限制

相位光栅的衍射效率 $\eta$：

$$\eta_m = \left|\frac{\sin(m\pi\Delta\phi/2\pi)}{m\pi\Delta\phi/2\pi}\right|^2$$

对于二元相位掩模（0, π）：
- $\eta_1 = 4/\pi^2 \approx 40.5\%$（一阶）
- $\eta_0 = 0$（零阶被抑制）

像素化效应引入 sinc 包络：

$$E_{pixel}(u,v) = \text{sinc}(au)\text{sinc}(bv)$$

其中 $a, b$ 是像素尺寸。填充因子 $F$ 降低效率：

$$\eta_{effective} = F^2 \cdot \eta_{ideal}$$

## 17.3 相位共轭与时间反演

### 17.3.1 光学相位共轭理论

相位共轭产生一个沿原始路径向后传播的波，从而消除畸变。对于前向波：

$$E_{forward}(\mathbf{r}, t) = A(\mathbf{r})e^{i[\mathbf{k}\cdot\mathbf{r} - \omega t + \phi(\mathbf{r})]}$$

相位共轭波为：

$$E_{conjugate}(\mathbf{r}, t) = A^*(\mathbf{r})e^{i[-\mathbf{k}\cdot\mathbf{r} - \omega t - \phi(\mathbf{r})]}$$

在四波混频中，共轭场产生于：

$$E_4 = \chi^{(3)} E_1 E_2 E_3^*$$

其中 $\chi^{(3)}$ 是三阶非线性磁化率。相位匹配条件：

$$\mathbf{k}_4 = \mathbf{k}_1 + \mathbf{k}_2 - \mathbf{k}_3$$

确保动量守恒。

### 17.3.2 数字相位共轭

通过全息记录和回放的数字实现：

**记录阶段：**
$$H(\mathbf{r}) = |E_{signal} + E_{reference}|^2 = |E_s|^2 + |E_r|^2 + E_s E_r^* + E_s^* E_r$$

**用共轭参考回放：**
$$E_{out} = H \cdot E_r^* = E_s|E_r|^2 + \text{other terms}$$

第一项代表相位共轭波。对于 SLM 实现：

$$\phi_{SLM}(x,y) = -\phi_{measured}(x,y) + \phi_{carrier}$$

### 17.3.3 时间反演对称性

波动方程的时间反演不变性使得光能够穿透无序介质聚焦：

$$\nabla^2 E - \frac{n^2(\mathbf{r})}{c^2}\frac{\partial^2 E}{\partial t^2} = 0$$

将 $t \to -t$ 代入，方程保持不变。对于单色场：

$$[\nabla^2 + k^2 n^2(\mathbf{r})]E = 0$$

如果 $E(\mathbf{r})$ 是一个解，那么 $E^*(\mathbf{r})$ 也是一个反向传播的解。这个原理是以下的基础：

1. **互易性**：对于传输矩阵 $T_{ij} = T_{ji}$
2. **聚焦不变性**：如果光通过无序介质从 A 聚焦到 B，则共轭将 B 聚焦到 A
3. **散射抵消**：多重散射路径在原始光源处相干叠加

## 17.4 穿透散射介质聚焦

### 17.4.1 传输矩阵形式

传输矩阵 $T$ 将通过散射介质的输入和输出光学模式关联起来：

$$\mathbf{E}_{out} = \mathbf{T} \cdot \mathbf{E}_{in}$$

对于 $N$ 个输入模式和 $M$ 个输出模式，$T$ 是一个 $M \times N$ 复数矩阵。每个元素：

$$T_{mn} = |T_{mn}|e^{i\phi_{mn}}$$

表示模式之间的振幅和相位耦合。奇异值分解：

$$\mathbf{T} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^†$$

揭示了最佳传输通道。传输特征值 $\tau_i = \sigma_i^2$ 遵循马尔琴科-帕斯图尔分布：

$$\rho(\tau) = \frac{1}{2\pi\gamma\tau}\sqrt{\frac{(\tau_+ - \tau)(\tau - \tau_-)}{1 + \tau}}$$

其中 $\gamma = M/N$ 且 $\tau_\pm = (1 \pm \sqrt{\gamma})^2$。

### 17.4.2 波前优化算法

**1. 迭代优化：**
通过相位控制最大化目标位置的强度：

$$I_{target} = \left|\sum_{n=1}^N T_{mn} A_n e^{i\phi_n}\right|^2$$

梯度上升更新：
$$\phi_n^{(k+1)} = \phi_n^{(k)} + \alpha \frac{\partial I_{target}}{\partial \phi_n}$$

其中：
$$\frac{\partial I_{target}}{\partial \phi_n} = 2\text{Im}[E_{target}^* T_{mn} A_n e^{i\phi_n}]$$

**2. 遗传算法：**
相位模式群体通过以下方式演化：
- 选择：基于适应度 $I_{target}$ 的轮盘赌
- 交叉：相位模式的均匀或算术混合
- 变异：高斯扰动 $\sigma \sim \pi/10$

**3. 哈达玛基测量：**
使用正交模式系统地测量 $T$：

$$\phi_{Hadamard}^{(k)} = H_{kn} \cdot \pi$$

其中 $H$ 是哈达玛矩阵。需要 4N 次测量才能获得完整的复数 $T$。

### 17.4.3 记忆效应与相关性

**角度记忆效应：**
散斑图案随输入角度变化而平移：

$$C(\Delta\theta) = \langle I(\theta)I(\theta + \Delta\theta) \rangle / \langle I \rangle^2$$

相关宽度：
$$\Delta\theta_{memory} \approx \frac{\lambda}{2\pi L}$$

其中 $L$ 是介质厚度。在此范围内：

$$T_{m,n+\Delta n} \approx T_{mn} e^{i\phi_{shift}(\Delta n)}$$

**光谱记忆效应：**
频率相关函数：

$$C(\Delta\omega) = \exp\left[-\left(\frac{\Delta\omega}{\Delta\omega_c}\right)^2\right]$$

相关带宽为：
$$\Delta\omega_c = \frac{2\pi c}{L_{path}\sqrt{\ln 2}}$$

其中 $L_{path}$ 是通过介质的平均路径长度。

**对成像的影响：**
1. 单次优化可在记忆效应范围内进行扫描
2. 在光谱相关范围内可实现宽带聚焦
3. 时间门控检测隔离弹道光子：$\tau_{gate} \sim L/c$

## 17.5 图形学中的自适应采样与优化

### 17.5.1 与重要性采样的联系

波前优化与蒙特卡洛渲染中的重要性采样并行。渲染方程：

$$L_o(\mathbf{x}, \omega_o) = \int_{\Omega} f_r(\mathbf{x}, \omega_i, \omega_o) L_i(\mathbf{x}, \omega_i) |\cos\theta_i| d\omega_i$$

重要性采样选择与被积函数成比例的方向 $\omega_i$。类似地，波前整形优化：

$$I_{focus} = \left|\int_A \psi_{in}(\mathbf{r}) G(\mathbf{r}, \mathbf{r}_{target}) dA\right|^2$$

其中 $G$ 是通过介质的格林函数。最佳输入场：

$$\psi_{in}^{opt}(\mathbf{r}) \propto G^*(\mathbf{r}, \mathbf{r}_{target})$$

这是传输核的相位共轭——类似于按 BSDF × 入射辐射度采样的比例。

### 17.5.2 自适应路径引导

**空间-方向树：**
根据辐射度变化自适应地划分 $(x, \omega)$ 空间：

$$\text{Split criterion} = \frac{\text{Var}[L(\mathbf{x}, \omega)]}{\text{Mean}[L(\mathbf{x}, \omega)]^2}$$

引导采样 PDF：

$$p_{guide}(\omega) = \alpha p_{learned}(\omega) + (1-\alpha) p_{BSDF}(\omega)$$

其中 $\alpha$ 根据学习置信度进行调整。

**高斯混合模型：**
将入射辐射场表示为：

$$L_i(\mathbf{x}, \omega) \approx \sum_{k=1}^K w_k \mathcal{N}(\omega; \mu_k, \Sigma_k)$$

通过 EM 算法更新：
- E 步：计算责任 $\gamma_k$
- M 步：更新 $\mu_k, \Sigma_k, w_k$

**路径空间马尔可夫链：**
在保持详细平衡的同时变异路径：

$$\frac{p(\mathbf{x} \to \mathbf{y})}{p(\mathbf{y} \to \mathbf{x})} = \frac{f(\mathbf{y})}{f(\mathbf{x})}$$

扰动由学习到的辐射度梯度引导：

$$\Delta\mathbf{x} = \epsilon \nabla_{\mathbf{x}} \log L(\mathbf{x})$$

### 17.5.3 学习光传输

**神经辐射缓存：**
训练网络预测入射辐射度：

$$L_i^{predicted}(\mathbf{x}, \omega) = \mathcal{N}_\theta(\mathbf{x}, \omega, \text{features})$$

损失函数结合 L2 和相对误差：

$$\mathcal{L} = \|L_i - L_i^{predicted}\|^2 + \lambda \frac{\|L_i - L_i^{predicted}\|^2}{\|L_i\|^2 + \epsilon}$$

**学习透射函数：**
对于参与介质，预测光学深度：

$$\tau_{predicted}(t) = \int_0^t \sigma_t^{predicted}(s) ds$$

网络架构利用对称性：
- 平移不变性 → 卷积层
- 旋转等变性 → 球谐函数特征
- 尺度不变性 → 多分辨率编码

**可微分渲染循环：**
通过梯度流连接到波前优化：

$$\frac{\partial L_{pixel}}{\partial \phi_{SLM}} = \sum_{paths} \frac{\partial L_{pixel}}{\partial \mathbf{x}_i} \frac{\partial \mathbf{x}_i}{\partial \phi_{SLM}}$$

这使得可以联合优化：
1. 采样分布
2. 去噪网络
3. 光传输近似

## 章总结

本章探讨了波前整形和自适应光学的数学基础和实际应用，将物理光学与计算图形学联系起来：

**关键概念：**
- **泽尼克多项式**提供了波前描述的正交基，实现了系统像差分析和校正
- **空间光调制器**通过液晶技术实现可编程相位控制，效率受像素化和填充因子限制
- **相位共轭**利用时间反演对称性消除传播畸变，实现穿透散射介质聚焦
- **传输矩阵**形式将光通过复杂介质的传播描述为模式之间的线性变换
- **记忆效应**允许单次优化在有限角度和光谱范围内工作
- 图形学中的**自适应采样**与光学中的波前优化共享数学原理

**关键方程：**
- 泽尼克展开：$W(\rho, \theta) = \sum_{n,m} a_n^m Z_n^m(\rho, \theta)$
- 相位共轭波：$E_{conjugate} = A^* e^{i[-\mathbf{k}\cdot\mathbf{r} - \omega t - \phi]}$
- 传输矩阵：$\mathbf{E}_{out} = \mathbf{T} \cdot \mathbf{E}_{in}$
- 记忆效应范围：$\Delta\theta_{memory} \approx \lambda/(2\pi L)$
- 最佳重要性采样：$\psi_{in}^{opt} \propto G^*(\mathbf{r}, \mathbf{r}_{target})$

## 练习

### 基本练习（建立直觉）

**练习 17.1：** 泽尼克多项式正交性
证明 $Z_2^0 = \sqrt{3}(2\rho^2 - 1)$（离焦）和 $Z_2^2 = \sqrt{6}\rho^2\cos(2\theta)$（散光）在单位圆盘上是正交的。
<details>
<summary>提示</summary>
使用正交性关系并记住 $\int_0^{2\pi} \cos(2\theta)d\theta = 0$。
</details>

<details>
<summary>答案</summary>

正交积分：
$$\int_0^{2\pi}\int_0^1 Z_2^0 Z_2^2 \rho d\rho d\theta = \sqrt{18}\int_0^{2\pi}\cos(2\theta)d\theta \int_0^1 (2\rho^2-1)\rho^3 d\rho$$

角度积分：$\int_0^{2\pi}\cos(2\theta)d\theta = 0$

因此，无论径向积分如何，乘积都为零，证明了正交性。
</details>

**练习 17.2：** SLM 相位缠绕
SLM 提供 0-2π 相位调制。为了在孔径上实现 4π 的总相移，必须添加什么空间频率光栅？进入一阶的衍射效率是多少？

<details>
<summary>提示</summary>
考虑闪耀光栅图案并使用相位光栅的衍射效率公式。
</details>

<details>
<summary>答案</summary>

添加线性相位斜坡：$\phi(x) = 4\pi x/D \bmod 2\pi$，创建一个周期为 D/2 的二元光栅。

对于二元（0,π）光栅，一阶效率：
$$\eta_1 = \left|\frac{\sin(\pi/2)}{\pi/2}\right|^2 = \frac{4}{\pi^2} \approx 0.405$$

40.5% 的光衍射到所需阶次。
</details>

**练习 17.3：** 记忆效应范围
一个厚度为 1mm，传输平均自由程 $\ell^*$ = 50μm 的散射介质在 $\lambda$ = 632.8nm 处被照亮。计算角度记忆效应范围以及 10° 视场中的独立散斑图案数量。

<details>
<summary>提示</summary>
使用 $\Delta\theta \approx \lambda/(2\pi L)$ 并考虑有多少个不重叠的范围适合总视场。
</details>

<details>
<summary>答案</summary>

记忆效应范围：
$$\Delta\theta = \frac{632.8 \times 10^{-9}}{2\pi \times 10^{-3}} = 1.01 \times 10^{-4} \text{ rad} = 0.0058°$$

独立图案数量：
$$N = \left(\frac{10°}{0.0058°}\right)^2 \approx 2.98 \times 10^6$$

视场中存在近 300 万个独立的优化区域。
</details>

### 挑战练习（扩展概念）

**练习 17.4：** 传输矩阵特征值边界
对于具有 $M$ 个输出模式和 $N$ 个输入模式（$M > N$）的传输矩阵 $T$，证明最大传输特征值满足 $\tau_{max} \leq (1 + \sqrt{M/N})^2$。

<details>
<summary>提示</summary>
使用马尔琴科-帕斯图尔分布并考虑能量守恒约束。
</details>

<details>
<summary>答案</summary>

根据马尔琴科-帕斯图尔分布，支持区间为 $[\tau_-, \tau_+]$，其中：
$$\tau_\pm = (1 \pm \sqrt{M/N})^2$$

能量守恒要求 $\text{Tr}(T^\dagger T) = N$（平均传输 = 1）。最大特征值出现在上边缘：
$$\tau_{max} = \tau_+ = (1 + \sqrt{M/N})^2$$

取传输特征值的平方根得到奇异值：
$$\sigma_{max} = 1 + \sqrt{M/N}$$

这个边界表明过采样（$M > N$）可以实现与 $\sqrt{M/N}$ 成比例的聚焦增强。
</details>

**练习 17.5：** 相位共轭保真度
相位共轭系统测量场 $E_{measured} = E_{true} + E_{noise}$，其中 $E_{noise}$ 是高斯噪声。推导保真度 $F = |⟨E_{true}|E_{conjugate}⟩|^2/|E_{true}|^2$ 作为 SNR 的函数。

<details>
<summary>提示</summary>
用测量场表示共轭场并计算重叠积分。
</details>

<details>
<summary>答案</summary>

共轭场：$E_{conjugate} = E_{measured}^* = E_{true}^* + E_{noise}^*$

保真度计算：
$$F = \frac{|⟨E_{true}|E_{true}^* + E_{noise}^*⟩|^2}{|E_{true}|^2}$$

$$F = \frac{||E_{true}|^2 + ⟨E_{true}|E_{noise}^*⟩|^2}{|E_{true}|^2}$$

对于方差为 $\sigma^2$ 的不相关高斯噪声：
$$F \approx \frac{|E_{true}|^4}{|E_{true}|^2(|E_{true}|^2 + \sigma^2)} = \frac{1}{1 + \sigma^2/|E_{true}|^2}$$

$$F = \frac{\text{SNR}}{1 + \text{SNR}}$$

高保真度需要 SNR >> 1。
</details>

**练习 17.6：** 自适应采样收敛
在路径引导中，样本从 $p_{guide}(\omega) = \alpha p_{learned}(\omega) + (1-\alpha)p_{BSDF}(\omega)$ 中抽取。推导最小化给定学习质量度量 $Q \in [0,1]$ 的方差的最佳 $\alpha$。

<details>
<summary>提示</summary>
最小化重要性采样估计量的方差，相对于 $\alpha$。
</details>

<details>
<summary>答案</summary>

重要性采样的方差：
$$\text{Var}[I] = \int \frac{f^2(\omega)}{p_{guide}(\omega)} d\omega - I^2$$

其中 $f(\omega)$ 是被积函数。定义学习质量：
$$Q = \frac{\int p_{learned} \cdot f}{\int p_{BSDF} \cdot f}$$

最小化相对于 $\alpha$ 的方差得到：
$$\alpha_{opt} = \frac{Q}{1 + Q}$$

当 $Q = 0$（学习效果差）时，$\alpha = 0$（仅使用 BSDF）
当 $Q = 1$（学习完美）时，$\alpha = 0.5$（平衡组合）
当 $Q \to \infty$（学习分布好得多）时，$\alpha \to 1$
</details>

**练习 17.7：** 开放问题 - 计算成像的波前编码
设计一个相位掩模 $\phi(x,y)$，同时满足：
1. 将景深扩展 10 倍
2. 在奈奎斯特频率处保持 MTF > 0.5
3. 实现单次深度估计

描述您的方法并推导点扩散函数。

<details>
<summary>提示</summary>
考虑三次相位掩模或优化的二元掩模。思考相位如何影响 PSF 的深度依赖性。
</details>

<details>
<summary>答案</summary>

这是一个开放的研究问题。一种潜在的方法：

**三次相位掩模：** $\phi(x,y) = \alpha(x^3 + y^3)$

PSF 变为：
$$h(x,y;z) = \left|\mathcal{F}\{P(u,v)\exp[i\phi(u,v)]\exp[i\frac{z}{2k}(u^2+v^2)]\}\right|^2$$

关键见解：
1. 三次相位使 PSF 在范围 $\Delta z \propto \alpha^{2/3}$ 内几乎与深度无关
2. 需要反卷积来恢复图像清晰度
3. PSF 不对称性使得可以从离焦中获取深度
4. 优化框架：

$$\min_\phi \int_{z_{min}}^{z_{max}} \text{Var}[h(x,y;z)] dz$$

受限于：$\text{MTF}(\nu_{Nyquist}) > 0.5$

现代方法使用与重建网络端到端优化的学习相位掩模。
</details>

## 常见陷阱与调试技巧

### 波前整形中的陷阱

1. **相位缠绕模糊**
   - 问题：相位测量是模 2π
   - 解决方案：使用时间解缠绕或多波长
   - 调试：检查相位梯度是否存在突然跳变

2. **SLM 校准误差**
   - 问题：非线性电压-相位响应
   - 解决方案：测量完整的校准曲线
   - 调试：使用已知相位图案（光栅）进行测试

3. **相干性要求**
   - 问题：部分相干性降低对比度
   - 解决方案：确保相干长度 > 光程差
   - 调试：测量干涉条纹的可见度

4. **采样混叠**
   - 问题：SLM 像素对高频相位欠采样
   - 解决方案：应用抗混叠滤波器或增加采样
   - 调试：检查傅里叶谱是否存在混叠伪影

5. **偏振敏感性**
   - 问题：SLM 依赖于偏振
   - 解决方案：使用适当的输入偏振态
   - 调试：使用偏振器/分析器对进行测试

### 性能优化技巧

1. **记忆效应利用**
   - 测量一次，在记忆效应范围内扫描
   - 对附近目标使用分组优化

2. **计算捷径**
   - 适用时使用 FFT 进行菲涅尔传播
   - 非线性相位响应的查找表
   - 矩阵运算的 GPU 加速

3. **降噪**
   - 平均多次测量
   - 对弱信号使用锁相检测
   - 在优化中实现异常值剔除

## 最佳实践清单

### 系统设计审查

- [ ] **相干性分析**
  - 相干长度 > 最大路径差？
  - 时间稳定性足以满足测量时间？
  - 空间相干性与光学系统匹配？

- [ ] **采样考虑**
  - 相位图案满足奈奎斯特准则？
  - SLM 分辨率足以满足所需 NA？
  - 相机动态范围足够？

- [ ] **校准协议**
  - 相位-电压校准已测量？
  - 像素串扰已表征？
  - 系统像差已补偿？

### 算法实现

- [ ] **优化策略**
  - 适用于 N 的算法（迭代式 vs. 矩阵式）？
  - 收敛准则已定义？
  - 局部最小值避免策略？

- [ ] **噪声处理**
  - SNR 估计已实现？
  - 针对噪声测量值的鲁棒优化？
  - 病态问题的正则化？

- [ ] **验证方法**
  - 是否有地面真值比较？
  - 性能指标已定义？
  - 边缘情况测试已完成？

### 与图形学集成

- [ ] **数学一致性**
  - 单位和约定与图形管线匹配？
  - 坐标系已正确转换？
  - 能量守恒已保持？

- [ ] **性能扩展**
  - 复杂度分析已执行？
  - 实时约束可实现？
  - 内存要求可接受？

- [ ] **泛化**
  - 方法是否扩展到相关的图形问题？
  - 假设已清晰记录？
  - 故障模式已理解？

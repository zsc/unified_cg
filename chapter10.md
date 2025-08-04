# 第10章：3D高斯溅射

3D高斯溅射（3D Gaussian Splatting）代表了神经渲染领域的一个重要突破，它将经典的高斯混合模型与现代可微渲染技术相结合，实现了高质量的实时新视图合成。本章将从统一体积渲染方程的角度深入探讨3D高斯溅射的数学基础，建立其与前述渲染技术的联系，并分析其在计算效率上的优势。我们将看到，通过巧妙地利用高斯函数的数学性质，3D高斯溅射在保持表达能力的同时，实现了比隐式神经表示快几个数量级的渲染速度。

## 学习目标

完成本章后，您将能够：

1. **推导3D高斯溅射的体积渲染方程形式**，理解其与统一框架的关系
2. **分析各向异性高斯核的数学性质**，掌握协方差矩阵的参数化方法
3. **理解可微光栅化的数学原理**，推导梯度的解析形式
4. **设计自适应密度控制算法**，分析其收敛性质
5. **评估实时渲染的计算复杂度**，理解并行化策略

## 10.1 场景的高斯混合模型

### 10.1.1 从离散点到连续表示

考虑一个3D场景由N个带属性的点组成。在第3章中，我们将点云表示为δ函数的集合：

$$\rho(\mathbf{x}) = \sum_{i=1}^{N} \alpha_i \delta(\mathbf{x} - \mathbf{x}_i)$$

其中$\alpha_i$是第i个点的不透明度，$\mathbf{x}_i$是其位置。然而，δ函数表示存在若干限制：
- **不连续性**：δ函数在空间中产生不连续的密度场
- **采样困难**：需要精确命中点位置才能获得非零值
- **缺乏空间支撑**：每个点的影响范围为零测度

3D高斯溅射通过用连续的高斯函数替换δ函数来解决这些问题：

$$\rho(\mathbf{x}) = \sum_{i=1}^{N} \alpha_i G_i(\mathbf{x})$$

其中$G_i(\mathbf{x})$是中心在$\mathbf{x}_i$的3D高斯函数。这种表示的优势包括：
- **平滑性**：$C^\infty$连续，便于梯度计算
- **有限支撑**：虽然理论上无限，但实际影响范围有限
- **解析性质**：高斯函数具有丰富的数学性质

**从点云到体积的过渡**：考虑点云的核密度估计（KDE）：
$$\hat{\rho}(\mathbf{x}) = \frac{1}{N} \sum_{i=1}^{N} K_h(\mathbf{x} - \mathbf{x}_i)$$

其中$K_h$是带宽为$h$的核函数。当选择高斯核时：
$$K_h(\mathbf{x}) = \frac{1}{(2\pi h^2)^{3/2}} \exp\left(-\frac{\|\mathbf{x}\|^2}{2h^2}\right)$$

3D高斯溅射可视为各向异性KDE的推广，其中每个点具有独立的带宽矩阵。

### 10.1.2 高斯基函数的数学性质

3D高斯函数定义为：

$$G_i(\mathbf{x}) = \frac{1}{(2\pi)^{3/2}|\mathbf{\Sigma}_i|^{1/2}} \exp\left(-\frac{1}{2}(\mathbf{x} - \mathbf{x}_i)^T \mathbf{\Sigma}_i^{-1} (\mathbf{x} - \mathbf{x}_i)\right)$$

其中$\mathbf{\Sigma}_i$是3×3协方差矩阵，决定了高斯的形状和方向。

**关键性质**：

1. **归一化**：$\int_{\mathbb{R}^3} G_i(\mathbf{x}) d\mathbf{x} = 1$
   
   证明：使用变量替换$\mathbf{z} = \mathbf{\Sigma}_i^{-1/2}(\mathbf{x} - \mathbf{x}_i)$

2. **仿射变换下的协变性**：若$\mathbf{y} = \mathbf{A}\mathbf{x} + \mathbf{b}$，则：
   $$G'(\mathbf{y}) = \frac{1}{|\det(\mathbf{A})|} G(\mathbf{A}^{-1}(\mathbf{y} - \mathbf{b}))$$
   
   新高斯的参数：$\mathbf{\mu}' = \mathbf{A}\mathbf{\mu} + \mathbf{b}$，$\mathbf{\Sigma}' = \mathbf{A}\mathbf{\Sigma}\mathbf{A}^T$

3. **卷积性质**：两个高斯的卷积仍是高斯：
   $$G_1 * G_2 = G_{12}, \quad \mathbf{\Sigma}_{12} = \mathbf{\Sigma}_1 + \mathbf{\Sigma}_2$$
   
   均值：$\mathbf{\mu}_{12} = \mathbf{\mu}_1 + \mathbf{\mu}_2$

4. **傅里叶变换**：高斯的傅里叶变换仍是高斯：
   $$\mathcal{F}\{G(\mathbf{x}; \mathbf{0}, \mathbf{\Sigma})\} = \exp\left(-\frac{1}{2}\mathbf{k}^T\mathbf{\Sigma}\mathbf{k}\right)$$

5. **矩生成函数**：
   $$M(\mathbf{t}) = \mathbb{E}[e^{\mathbf{t}^T\mathbf{x}}] = \exp\left(\mathbf{t}^T\mathbf{\mu} + \frac{1}{2}\mathbf{t}^T\mathbf{\Sigma}\mathbf{t}\right)$$

6. **条件分布**：高斯的条件分布仍是高斯
   
   若$\mathbf{x} = [\mathbf{x}_1^T, \mathbf{x}_2^T]^T \sim \mathcal{N}(\mathbf{\mu}, \mathbf{\Sigma})$，则：
   $$\mathbf{x}_1|\mathbf{x}_2 \sim \mathcal{N}(\mathbf{\mu}_{1|2}, \mathbf{\Sigma}_{1|2})$$
   
   其中$\mathbf{\mu}_{1|2} = \mathbf{\mu}_1 + \mathbf{\Sigma}_{12}\mathbf{\Sigma}_{22}^{-1}(\mathbf{x}_2 - \mathbf{\mu}_2)$

7. **信息形式**：使用精度矩阵$\mathbf{\Lambda} = \mathbf{\Sigma}^{-1}$：
   $$G(\mathbf{x}) \propto \exp\left(-\frac{1}{2}\mathbf{x}^T\mathbf{\Lambda}\mathbf{x} + \mathbf{\eta}^T\mathbf{x}\right)$$
   
   其中$\mathbf{\eta} = \mathbf{\Lambda}\mathbf{\mu}$是信息向量

**数值考虑**：
- 当$\|\mathbf{x} - \mathbf{\mu}\| > 3\sqrt{\lambda_{max}(\mathbf{\Sigma})}$时，高斯值小于$0.01 \times$峰值
- 可以安全地将高斯截断在$3\sigma$范围内，误差$< 0.3\%$

### 10.1.3 高斯混合模型的概率解释

从概率论角度，场景可视为一个高斯混合模型（GMM）：

$$p(\mathbf{x}) = \sum_{i=1}^{N} \pi_i \mathcal{N}(\mathbf{x}; \mathbf{\mu}_i, \mathbf{\Sigma}_i)$$

其中$\pi_i$是混合权重，满足$\sum_i \pi_i = 1$。在渲染中，我们将$\pi_i$与不透明度$\alpha_i$关联。

**从概率到渲染的映射**：
- **混合权重**：$\pi_i \propto \alpha_i$（未归一化）
- **密度函数**：$\rho(\mathbf{x}) = \sum_i \alpha_i \mathcal{N}(\mathbf{x}; \mathbf{\mu}_i, \mathbf{\Sigma}_i)$
- **颜色分布**：每个高斯关联颜色$\mathbf{c}_i$

**期望最大化（EM）视角**：
给定观测图像集$\{I_j\}$，可以通过EM算法优化GMM参数：

1. **E步**：计算后验概率
   $$\gamma_{ji} = \frac{\pi_i \mathcal{N}(\mathbf{x}_j; \mathbf{\mu}_i, \mathbf{\Sigma}_i)}{\sum_k \pi_k \mathcal{N}(\mathbf{x}_j; \mathbf{\mu}_k, \mathbf{\Sigma}_k)}$$

2. **M步**：更新参数
   $$\mathbf{\mu}_i = \frac{\sum_j \gamma_{ji} \mathbf{x}_j}{\sum_j \gamma_{ji}}, \quad \mathbf{\Sigma}_i = \frac{\sum_j \gamma_{ji} (\mathbf{x}_j - \mathbf{\mu}_i)(\mathbf{x}_j - \mathbf{\mu}_i)^T}{\sum_j \gamma_{ji}}$$

但在实际的3D高斯溅射中，我们采用梯度下降直接优化渲染损失。

**信息论解释**：
高斯混合模型的熵：
$$H[p] = -\int p(\mathbf{x}) \log p(\mathbf{x}) d\mathbf{x}$$

对于单个高斯：$H[\mathcal{N}(\mathbf{\mu}, \mathbf{\Sigma})] = \frac{3}{2}\log(2\pi e) + \frac{1}{2}\log|\mathbf{\Sigma}|$

GMM的熵没有解析形式，但可以通过上下界估计：
$$\max_i H[\mathcal{N}(\mathbf{\mu}_i, \mathbf{\Sigma}_i)] \leq H[p] \leq \sum_i \pi_i H[\mathcal{N}(\mathbf{\mu}_i, \mathbf{\Sigma}_i)] - \sum_i \pi_i \log \pi_i$$

### 10.1.4 与统一体积渲染方程的联系

回顾第3章的统一体积渲染方程：

$$L(\mathbf{r}) = \int_0^{\infty} T(t) \sigma(t) c(t) dt$$

其中$T(t) = \exp\left(-\int_0^t \sigma(s) ds\right)$是透射率。

对于高斯混合表示，密度场变为：

$$\sigma(\mathbf{r}(t)) = \sum_{i=1}^{N} \alpha_i G_i(\mathbf{r}(t))$$

代入体积渲染方程：

$$L(\mathbf{r}) = \int_0^{\infty} \exp\left(-\int_0^t \sum_{j=1}^{N} \alpha_j G_j(\mathbf{r}(s)) ds\right) \sum_{i=1}^{N} \alpha_i G_i(\mathbf{r}(t)) \mathbf{c}_i dt$$

这个积分一般没有解析解，但通过适当的近似，我们可以得到高效的计算方法。

**离散化近似**：
将射线离散化为$M$个采样点：
$$L(\mathbf{r}) \approx \sum_{k=1}^{M} T_k \sigma_k c_k \Delta t$$

其中：
- $\sigma_k = \sum_{i=1}^{N} \alpha_i G_i(\mathbf{r}(t_k))$
- $T_k = \exp\left(-\sum_{j=1}^{k-1} \sigma_j \Delta t\right)$
- $c_k = \frac{\sum_{i=1}^{N} \alpha_i G_i(\mathbf{r}(t_k)) \mathbf{c}_i}{\sigma_k}$

**局部线性化近似**：
在每个高斯的支撑域内，假设射线近似为直线：
$$\mathbf{r}(t) \approx \mathbf{r}_0 + t\mathbf{d}$$

这使得我们可以将沿射线的积分转化为对高斯的求和。

**多分辨率表示**：
为了加速计算，可以构建高斯的层次结构：
$$\sigma_{coarse}(\mathbf{x}) = \sum_{j} \beta_j G_j^{coarse}(\mathbf{x})$$

其中粗糍级别的高斯是精细级别的聚合。

**与NeRF的对比**：
- NeRF：$\sigma(\mathbf{x}) = MLP(\gamma(\mathbf{x}))$（隐式表示）
- 3D GS：$\sigma(\mathbf{x}) = \sum_i \alpha_i G_i(\mathbf{x})$（显式表示）

显式表示的优势：
- 可以预计算和索引高斯
- 避免了MLP的重复计算
- 更容易并行化

## 10.2 各向异性核与协方差

### 10.2.1 3D高斯的参数化

每个3D高斯由以下参数完全确定：
- 位置：$\mathbf{\mu}_i \in \mathbb{R}^3$
- 协方差矩阵：$\mathbf{\Sigma}_i \in \mathbb{R}^{3 \times 3}$（对称正定）
- 不透明度：$\alpha_i \in [0, 1]$
- 颜色：$\mathbf{c}_i$（可以是RGB或球谐系数）

协方差矩阵必须满足对称正定条件：
$$\mathbf{\Sigma}_i = \mathbf{\Sigma}_i^T, \quad \mathbf{x}^T\mathbf{\Sigma}_i\mathbf{x} > 0, \forall \mathbf{x} \neq \mathbf{0}$$

**参数计数**：
- 位置：3个参数
- 协方差：6个独立参数（对称矩阵）
- 不透明度：1个参数
- 颜色：3个（RGB）或更多（球谐系数）

总计：每个高斯至少13个参数

**颜色的视角依赖性**：
使用球谐函数表示颜色：
$$\mathbf{c}(\mathbf{d}) = \sum_{l=0}^{L} \sum_{m=-l}^{l} c_{lm} Y_l^m(\mathbf{d})$$

其中$\mathbf{d}$是视角方向，$Y_l^m$是球谐基函数。

常用阶数：
- $L=0$：常数颜色（1个系数）
- $L=1$：线性变化（4个系数）
- $L=2$：二次变化（9个系数）
- $L=3$：三次变化（16个系数）

### 10.2.2 协方差矩阵的几何意义

协方差矩阵的特征分解揭示了高斯的几何结构：

$$\mathbf{\Sigma} = \mathbf{V}\mathbf{\Lambda}\mathbf{V}^T$$

其中：
- $\mathbf{V} = [\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3]$：特征向量矩阵（正交）
- $\mathbf{\Lambda} = \text{diag}(\lambda_1, \lambda_2, \lambda_3)$：特征值对角矩阵

几何解释：
- 特征向量$\mathbf{v}_i$定义了椭球的主轴方向
- 特征值$\lambda_i$决定了沿各主轴的方差
- 等概率密度面是椭球：$(\mathbf{x} - \mathbf{\mu})^T\mathbf{\Sigma}^{-1}(\mathbf{x} - \mathbf{\mu}) = c$

**椭球半轴长度**：
对于给定的置信水平$p$，椭球半轴长度为：
$$a_i = \sqrt{\lambda_i \chi^2_3(p)}$$

其中$\chi^2_3(p)$是自由度为3的卡方分布的$p$分位数。

常用值：
- $p = 0.68$：$\chi^2_3(0.68) \approx 3.51$ （1-$\sigma$椭球）
- $p = 0.95$：$\chi^2_3(0.95) \approx 7.81$ （2-$\sigma$椭球）
- $p = 0.997$：$\chi^2_3(0.997) \approx 14.16$ （3-$\sigma$椭球）

**体积与表面积**：
- 椭球体积：$V = \frac{4}{3}\pi \sqrt{\det(\mathbf{\Sigma})} \cdot (\chi^2_3(p))^{3/2}$
- 有效半径：$r_{eff} = (\det(\mathbf{\Sigma}))^{1/6}$

**各向异性度量**：
定义各向异性指标：
$$\text{anisotropy} = 1 - \frac{\lambda_{min}}{\lambda_{max}}$$

值域：
- 0：完全各向同性（球形）
- 接近1：高度各向异性（细长椭球）

**条件数与数值稳定性**：
$$\kappa(\mathbf{\Sigma}) = \frac{\lambda_{max}}{\lambda_{min}}$$

当$\kappa$很大时，矩阵接近奇异，可能导致数值问题。

### 10.2.3 旋转与缩放的分解

为了保证正定性并简化优化，我们采用如下参数化：

$$\mathbf{\Sigma} = \mathbf{R}\mathbf{S}\mathbf{S}^T\mathbf{R}^T$$

其中：
- $\mathbf{R} \in SO(3)$：旋转矩阵
- $\mathbf{S} = \text{diag}(s_1, s_2, s_3)$：缩放矩阵，$s_i > 0$

这种参数化的优点：
1. **自动保证正定性**：只要$s_i > 0$
2. **直观的几何意义**：旋转+缩放
3. **参数数量合理**：3（旋转）+ 3（缩放）= 6个自由度

旋转矩阵可用四元数$\mathbf{q} = (q_w, q_x, q_y, q_z)$表示：

$$\mathbf{R}(\mathbf{q}) = \begin{bmatrix}
1-2(q_y^2+q_z^2) & 2(q_x q_y - q_w q_z) & 2(q_x q_z + q_w q_y) \\
2(q_x q_y + q_w q_z) & 1-2(q_x^2+q_z^2) & 2(q_y q_z - q_w q_x) \\
2(q_x q_z - q_w q_y) & 2(q_y q_z + q_w q_x) & 1-2(q_x^2+q_y^2)
\end{bmatrix}$$

**四元数的优势**：
- 避免万向锁
- 插值平滑（SLERP）
- 只需归一化约束：$\|\mathbf{q}\| = 1$
- 存储4个参数，优化时可使用3个自由度

**从四元数到旋转矩阵的梯度**：
$$\frac{\partial \mathbf{R}}{\partial q_w} = 2\begin{bmatrix}
0 & -q_z & q_y \\
q_z & 0 & -q_x \\
-q_y & q_x & 0
\end{bmatrix}$$

类似地可以计算$\frac{\partial \mathbf{R}}{\partial q_x}$, $\frac{\partial \mathbf{R}}{\partial q_y}$, $\frac{\partial \mathbf{R}}{\partial q_z}$。

**参数化的等价性**：
注意到$(\mathbf{R}, \mathbf{S})$和$(-\mathbf{R}, \mathbf{S})$产生相同的$\mathbf{\Sigma}$。这在优化时可能导致不唯一性。

**替代参数化方案**：
1. **欧拉角**：简单但存在万向锁
2. **轴角表示**：$\mathbf{R} = \exp([\mathbf{\omega}]_\times)$
3. **Cayley变换**：$\mathbf{R} = (\mathbf{I} - \mathbf{K})(\mathbf{I} + \mathbf{K})^{-1}$

### 10.2.4 保证正定性的参数化技巧

在优化过程中，我们需要确保协方差矩阵始终保持正定。常用技巧：

1. **对数尺度参数化**：
   $$s_i = \exp(\tilde{s}_i)$$
   优化$\tilde{s}_i \in \mathbb{R}$而非$s_i$
   
   梯度传播：$\frac{\partial \mathcal{L}}{\partial \tilde{s}_i} = \frac{\partial \mathcal{L}}{\partial s_i} \cdot s_i$

2. **Cholesky分解**：
   $$\mathbf{\Sigma} = \mathbf{L}\mathbf{L}^T$$
   其中$\mathbf{L}$是下三角矩阵，对角元素为正
   
   参数化：
   $$\mathbf{L} = \begin{bmatrix}
   \exp(l_{11}) & 0 & 0 \\
   l_{21} & \exp(l_{22}) & 0 \\
   l_{31} & l_{32} & \exp(l_{33})
   \end{bmatrix}$$

3. **正则化**：添加小的对角项
   $$\mathbf{\Sigma}_{reg} = \mathbf{\Sigma} + \epsilon\mathbf{I}$$
   
   典型值：$\epsilon = 10^{-4}$到$10^{-6}$

4. **投影到正定锥**：
   若$\mathbf{\Sigma}$变得非正定，投影到最近的正定矩阵：
   $$\mathbf{\Sigma}_{proj} = \arg\min_{\mathbf{\Sigma}' \succ 0} \|\mathbf{\Sigma}' - \mathbf{\Sigma}\|_F$$
   
   解法：特征分解后将负特征值替换为$\epsilon$

**数值稳定性考虑**：
- 避免过小的特征值（$\lambda_i > \epsilon$）
- 条件数控制：$\kappa(\mathbf{\Sigma}) = \lambda_{max}/\lambda_{min} < \kappa_{max}$
- 梯度裁剪防止数值爆炸

**梯度裁剪策略**：
$$\nabla_{clipped} = \begin{cases}
\nabla & \text{if } \|\nabla\| \leq \tau \\
\tau \frac{\nabla}{\|\nabla\|} & \text{otherwise}
\end{cases}$$

其中$\tau$是裁剪阈值，典型值为1.0到5.0。

**自适应步长**：
根据条件数调整学习率：
$$\eta_{adaptive} = \frac{\eta_{base}}{1 + \beta \cdot (\kappa - \kappa_{target})}$$

其中$\kappa_{target} \approx 100$是目标条件数。

## 10.3 可微光栅化

### 10.3.1 从3D到2D的投影

给定相机参数（内参$\mathbf{K}$，外参$[\mathbf{R}_c|\mathbf{t}_c]$），3D高斯投影到图像平面的过程包括：

1. **世界坐标到相机坐标**：
   $$\mathbf{x}_c = \mathbf{R}_c(\mathbf{x}_w - \mathbf{t}_c)$$

2. **透视投影**：
   $$\mathbf{x}_{ndc} = \begin{bmatrix} x_c/z_c \\ y_c/z_c \\ 1 \end{bmatrix}$$

3. **归一化设备坐标到像素坐标**：
   $$\mathbf{x}_{pixel} = \mathbf{K}\mathbf{x}_{ndc}$$

关键问题：3D高斯投影后还是高斯吗？

**内参矩阵的结构**：
$$\mathbf{K} = \begin{bmatrix}
f_x & 0 & c_x \\
0 & f_y & c_y \\
0 & 0 & 1
\end{bmatrix}$$

其中$(f_x, f_y)$是焦距，$(c_x, c_y)$是主点。

**完整的投影变换**：
$$\begin{bmatrix} u \\ v \\ 1 \end{bmatrix} = \frac{1}{z_c} \mathbf{K} \mathbf{R}_c (\mathbf{x}_w - \mathbf{t}_c)$$

注意到这个变换包含了非线性的除法操作（$1/z_c$），因此不是仿射变换。

### 10.3.2 高斯的仿射变换性质

**定理**：3D高斯在仿射变换下保持高斯性。

证明：设3D高斯$G(\mathbf{x}; \mathbf{\mu}, \mathbf{\Sigma})$，仿射变换$\mathbf{y} = \mathbf{A}\mathbf{x} + \mathbf{b}$，则：

$$G'(\mathbf{y}) = G(\mathbf{y}; \mathbf{A}\mathbf{\mu} + \mathbf{b}, \mathbf{A}\mathbf{\Sigma}\mathbf{A}^T)$$

但透视投影不是仿射变换！为了保持高斯性，我们采用局部仿射近似。

**局部仿射近似**：在高斯中心$\mathbf{\mu}$处，用雅可比矩阵近似投影：

$$\mathbf{J} = \frac{\partial \mathbf{x}_{2D}}{\partial \mathbf{x}_{3D}}\bigg|_{\mathbf{x}=\mathbf{\mu}}$$

对于透视投影：
$$\mathbf{J} = \begin{bmatrix}
f_x/z & 0 & -f_x x/z^2 \\
0 & f_y/z & -f_y y/z^2
\end{bmatrix}$$

其中$(f_x, f_y)$是焦距，$(x, y, z)$是相机坐标。

投影后的2D协方差：
$$\mathbf{\Sigma}_{2D} = \mathbf{J}\mathbf{\Sigma}_{3D}\mathbf{J}^T$$

**近似误差分析**：
局部仿射近似的误差为$\mathcal{O}(\|\mathbf{x} - \mathbf{\mu}\|^2)$。对于典型的高斯（$3\sigma$范围），相对误差通常小于5%。

**投影雅可比的详细推导**：
设相机坐标为$\mathbf{x}_c = [x, y, z]^T$，投影后的像素坐标为：
$$u = f_x \frac{x}{z} + c_x, \quad v = f_y \frac{y}{z} + c_y$$

计算偏导数：
$$\frac{\partial u}{\partial x} = \frac{f_x}{z}, \quad \frac{\partial u}{\partial y} = 0, \quad \frac{\partial u}{\partial z} = -\frac{f_x x}{z^2}$$
$$\frac{\partial v}{\partial x} = 0, \quad \frac{\partial v}{\partial y} = \frac{f_y}{z}, \quad \frac{\partial v}{\partial z} = -\frac{f_y y}{z^2}$$

考虑到世界到相机的变换，完整的雅可比为：
$$\mathbf{J}_{full} = \mathbf{J} \cdot \mathbf{R}_c$$

**2D高斯的有效计算**：
给定2D协方差$\mathbf{\Sigma}_{2D}$，其逆矩阵为：
$$\mathbf{\Sigma}_{2D}^{-1} = \frac{1}{\det(\mathbf{\Sigma}_{2D})} \begin{bmatrix}
\sigma_{22} & -\sigma_{12} \\
-\sigma_{12} & \sigma_{11}
\end{bmatrix}$$

其中$\det(\mathbf{\Sigma}_{2D}) = \sigma_{11}\sigma_{22} - \sigma_{12}^2$。

### 10.3.3 α-blending与排序

给定像素处的多个高斯贡献，最终颜色通过α-blending计算：

$$C = \sum_{i=1}^{N} c_i \alpha_i \prod_{j=1}^{i-1} (1 - \alpha_j)$$

其中排序按深度从前到后。这等价于：

$$C = \sum_{i=1}^{N} c_i \alpha_i T_i, \quad T_i = \prod_{j=1}^{i-1} (1 - \alpha_j)$$

每个高斯在像素$(u,v)$处的贡献：
$$\alpha_i(u,v) = \alpha_i^{base} \cdot G_{2D}(u,v; \mathbf{\mu}_{2D,i}, \mathbf{\Sigma}_{2D,i})$$

**深度排序算法**：
1. **快速排序**：$\mathcal{O}(N \log N)$，但在GPU上效率不高
2. **基数排序**：$\mathcal{O}(N)$，GPU友好，适合固定精度深度
3. **层次深度缓冲**：将深度范围划分为桶，在桶内排序

**透明度累积的数值稳定性**：
由于连续乘积$(1-\alpha_j)$，当$N$很大时可能出现数值下溢。解决方案：

1. **提前终止**：当$T_i < \epsilon$时停止累积
2. **对数空间计算**：
   $$\log T_i = \sum_{j=1}^{i-1} \log(1 - \alpha_j)$$

**像素级别的高斯筛选**：
为了减少计算量，只处理对像素有显著贡献的高斯：
1. **边界框测试**：计算2D高斯的边界框
2. **阈值测试**：$\alpha_i \cdot G_{2D}(u,v) > \tau_{min}$
3. **视锥体剪裁**：在投影前剔除不可见的高斯

### 10.3.4 梯度计算与反向传播

为了优化高斯参数，我们需要计算损失函数对所有参数的梯度。设损失函数$\mathcal{L} = \|C_{rendered} - C_{gt}\|^2$。

**链式法则**：
$$\frac{\partial \mathcal{L}}{\partial \theta} = \frac{\partial \mathcal{L}}{\partial C} \cdot \frac{\partial C}{\partial \alpha} \cdot \frac{\partial \alpha}{\partial G_{2D}} \cdot \frac{\partial G_{2D}}{\partial \mathbf{\Sigma}_{2D}} \cdot \frac{\partial \mathbf{\Sigma}_{2D}}{\partial \theta}$$

其中$\theta$可以是位置$\mathbf{\mu}$、旋转$\mathbf{q}$、缩放$\mathbf{s}$等参数。

**关键梯度**：

1. **颜色对不透明度**：
   $$\frac{\partial C}{\partial \alpha_i} = c_i T_i - \sum_{j=i+1}^{N} c_j \alpha_j T_j / (1-\alpha_i)$$

2. **2D高斯对协方差**：
   $$\frac{\partial G_{2D}}{\partial \mathbf{\Sigma}_{2D}} = \frac{1}{2}G_{2D} \left[\mathbf{\Sigma}_{2D}^{-1}(\mathbf{x}-\mathbf{\mu})(\mathbf{x}-\mathbf{\mu})^T\mathbf{\Sigma}_{2D}^{-1} - \mathbf{\Sigma}_{2D}^{-1}\right]$$

3. **2D协方差对3D参数**：通过雅可比矩阵的导数计算

**数值稳定性**：
- 使用对数空间计算避免下溢
- 梯度裁剪防止爆炸
- 自适应学习率调整

**完整的梯度计算流程**：

1. **损失对颜色的梯度**：
   $$\frac{\partial \mathcal{L}}{\partial C} = 2(C_{rendered} - C_{gt})$$

2. **颜色对高斯属性的梯度**：
   - 对基础不透明度：$\frac{\partial C}{\partial \alpha_i^{base}} = G_{2D}(u,v) \cdot T_i \cdot (c_i - C_{behind})$
   - 对颜色：$\frac{\partial C}{\partial c_i} = \alpha_i T_i$
   - 对位置：$\frac{\partial C}{\partial \mathbf{\mu}_i} = \alpha_i T_i (c_i - C_{behind}) \cdot \frac{\partial G_{2D}}{\partial \mathbf{\mu}_{2D}} \cdot \mathbf{J}$

其中$C_{behind} = \frac{\sum_{j=i+1}^{N} c_j \alpha_j T_j}{T_{i+1}}$是i之后所有高斯的加权颜色。

3. **坐标变换的梯度传播**：
   $$\frac{\partial \mathbf{\mu}_{2D}}{\partial \mathbf{\mu}_{3D}} = \mathbf{J}$$
   $$\frac{\partial \mathbf{\Sigma}_{2D}}{\partial \mathbf{\Sigma}_{3D}} = \mathbf{J} \otimes \mathbf{J}$$

其中$\otimes$表示Kronecker积。

4. **参数化的梯度**：
   - 对旋转四元数：$\frac{\partial \mathbf{\Sigma}}{\partial q_k} = \frac{\partial \mathbf{R}}{\partial q_k} \mathbf{S}\mathbf{S}^T\mathbf{R}^T + \mathbf{R}\mathbf{S}\mathbf{S}^T\frac{\partial \mathbf{R}^T}{\partial q_k}$
   - 对缩放：$\frac{\partial \mathbf{\Sigma}}{\partial s_k} = 2\mathbf{R}\text{diag}(0,...,s_k,...,0)\mathbf{R}^T$

**优化策略**：
1. **Adam优化器**：适应性学习率，动量项
2. **学习率调度**：
   - 位置：$\eta_{\mu} = 1.6 \times 10^{-4}$
   - 旋转：$\eta_q = 1.0 \times 10^{-3}$
   - 缩放：$\eta_s = 5.0 \times 10^{-3}$
   - 不透明度：$\eta_\alpha = 5.0 \times 10^{-2}$
3. **梯度裁剪阈值**：通常设置为1.0到5.0

## 10.4 自适应密度控制

### 10.4.1 高斯的分裂与克隆

在优化过程中，固定数量的高斯可能无法充分表示复杂场景。自适应密度控制通过动态调整高斯数量来提高重建质量。

**分裂策略**：当高斯满足以下条件时进行分裂：
1. **梯度阈值**：$\|\nabla_{\mathbf{\mu}} \mathcal{L}\| > \tau_{grad}$
2. **尺寸阈值**：$\max(s_1, s_2, s_3) > \tau_{size}$

分裂操作：
$$\begin{aligned}
\mathbf{\mu}_{new,1} &= \mathbf{\mu} + \epsilon \mathbf{v}_1 \\
\mathbf{\mu}_{new,2} &= \mathbf{\mu} - \epsilon \mathbf{v}_1 \\
\mathbf{s}_{new} &= \mathbf{s} / \phi
\end{aligned}$$

其中$\mathbf{v}_1$是最大特征值对应的特征向量，$\phi > 1$是缩放因子。

**克隆策略**：对于梯度大但尺寸小的高斯，进行克隆：
$$\mathbf{\mu}_{clone} = \mathbf{\mu} + \epsilon \cdot \text{randn}(3)$$

### 10.4.2 修剪策略

为了控制内存使用和计算成本，需要移除贡献小的高斯：

1. **不透明度阈值**：$\alpha < \tau_{\alpha}$
2. **视锥体剔除**：高斯完全在视锥体外
3. **屏幕空间贡献**：投影面积 < $\tau_{area}$

**贡献度量**：
$$\text{contribution}_i = \alpha_i \cdot \text{area}_{2D,i} \cdot \exp(-d_i^2/2\sigma_d^2)$$

其中$d_i$是到相机的距离，$\sigma_d$控制距离衰减。

### 10.4.3 密度的自适应调整

**目标函数**：平衡重建质量和高斯数量：
$$\mathcal{L}_{total} = \mathcal{L}_{render} + \lambda_{sparse} \|\boldsymbol{\alpha}\|_1 + \lambda_{compact} \sum_i \text{vol}(\mathbf{\Sigma}_i)$$

其中：
- $\mathcal{L}_{render}$：渲染损失（L1或L2）
- $\|\boldsymbol{\alpha}\|_1$：稀疏性正则化
- $\text{vol}(\mathbf{\Sigma}_i) = \sqrt{\det(\mathbf{\Sigma}_i)}$：体积正则化

**自适应调整算法**：
```
每K次迭代：
1. 计算所有高斯的贡献度
2. 分裂高梯度大尺寸的高斯
3. 克隆高梯度小尺寸的高斯
4. 修剪低贡献的高斯
5. 更新正则化权重
```

### 10.4.4 收敛性分析

**定理**：在适当的条件下，自适应密度控制算法收敛到局部最优。

**假设**：
1. 损失函数$\mathcal{L}$是Lipschitz连续的
2. 分裂/克隆操作保持有界性
3. 学习率满足Robbins-Monro条件

**收敛速率**：设$N_t$是第$t$次迭代的高斯数量，则：
$$\mathbb{E}[\mathcal{L}_t - \mathcal{L}^*] \leq \mathcal{O}\left(\frac{\log N_t}{\sqrt{t}}\right)$$

**实践考虑**：
- 初始化影响收敛速度
- 自适应操作的频率影响稳定性
- 正则化权重需要动态调整

**数值实验表明**：
- 典型场景需要$10^4 - 10^6$个高斯
- 收敛通常在$10^4 - 10^5$次迭代内
- 自适应策略比固定数量提升20-50%质量

## 10.5 实时考虑

### 10.5.1 计算复杂度分析

3D高斯溅射的渲染复杂度取决于多个因素：

**前向渲染复杂度**：
- 投影变换：$\mathcal{O}(N)$
- 深度排序：$\mathcal{O}(N \log N)$或$\mathcal{O}(N)$（基数排序）
- 光栅化：$\mathcal{O}(N \cdot P)$，其中$P$是每个高斯覆盖的平均像素数

总复杂度：$\mathcal{O}(N \log N + N \cdot P)$

**内存复杂度**：
- 高斯参数存储：$\mathcal{O}(N)$
- 排序缓冲区：$\mathcal{O}(N)$
- 帧缓冲区：$\mathcal{O}(W \times H)$

**与其他方法比较**：
| 方法 | 时间复杂度 | 内存复杂度 |
|------|------------|------------|
| NeRF | $\mathcal{O}(W \times H \times S)$ | $\mathcal{O}(M)$ |
| 3D GS | $\mathcal{O}(N \log N + N \cdot P)$ | $\mathcal{O}(N)$ |
| Mesh | $\mathcal{O}(T)$ | $\mathcal{O}(V + T)$ |

其中$S$是采样点数，$M$是网络参数数，$T$是三角形数，$V$是顶点数。

### 10.5.2 并行化策略

**GPU并行化**：

1. **高斯投影并行**：
   ```
   并行对每个高斯：
     计算相机坐标
     计算2D投影参数
     计算边界框
   ```

2. **分块光栅化**：
   ```
   将图像分成tiles（如16×16）
   并行对每个tile：
     收集相交的高斯
     局部深度排序
     α-blending
   ```

3. **排序优化**：
   - 使用GPU基数排序
   - 分层深度缓冲
   - 近似排序（对远处物体）

**SIMD优化**：
- 向量化高斯计算
- 批量矩阵运算
- 融合操作减少内存访问

### 10.5.3 内存带宽优化

内存带宽是主要瓶颈，优化策略：

1. **数据布局**：
   - Structure of Arrays (SoA) vs Array of Structures (AoS)
   - 对齐访问模式
   - 紧凑表示（如压缩球谐系数）

2. **缓存优化**：
   - 空间局部性：分块处理
   - 时间局部性：融合kernel
   - 纹理缓存利用

3. **带宽计算**：
   ```
   每帧带宽 = N × (参数读取 + 投影写入) + 
              像素数 × (深度测试 + 颜色累积)
   ```

   典型值：$10^6$个高斯，1080p分辨率
   - 参数读取：~200MB
   - 中间结果：~100MB
   - 帧缓冲：~8MB
   总计：~300MB/帧 → 9GB/s @ 30fps

### 10.5.4 与传统光栅化的比较

**优势**：
1. **连续表示**：无需拓扑信息
2. **软边界**：自然的抗锯齿
3. **透明度**：原生支持
4. **可微性**：易于优化

**劣势**：
1. **排序开销**：需要深度排序
2. **过度绘制**：重叠高斯
3. **内存占用**：每个高斯多个参数

**混合方案**：
- 近处使用高斯溅射（细节）
- 远处使用传统网格（效率）
- LOD系统动态切换

**性能指标**（典型硬件）：
- 渲染速度：100-200 FPS @ 1080p
- 训练速度：5-30分钟（取决于场景）
- 内存使用：200MB-2GB（场景相关）

**优化建议**：
1. 使用视锥体剔除减少处理的高斯数
2. 实现层次化数据结构（八叉树）
3. 动态调整质量参数
4. 利用时间连贯性（运动场景）

## 本章小结

3D高斯溅射通过将场景表示为高斯混合模型，成功地在表达能力和计算效率之间取得了平衡。关键创新包括：

1. **统一体积渲染视角**：3D高斯溅射可以视为体积渲染方程的一种特殊实现，其中密度场由高斯混合表示
2. **各向异性表示**：通过协方差矩阵的巧妙参数化，实现了灵活的形状表示同时保证数值稳定性
3. **可微光栅化**：利用高斯的数学性质，实现了高效的梯度计算和端到端优化
4. **自适应密度控制**：动态调整高斯数量，自动适应场景复杂度
5. **实时性能**：通过并行化和优化，达到了实时渲染速度

重要公式回顾：
- 高斯混合密度场：$\sigma(\mathbf{x}) = \sum_{i=1}^{N} \alpha_i G_i(\mathbf{x})$
- 协方差参数化：$\mathbf{\Sigma} = \mathbf{R}\mathbf{S}\mathbf{S}^T\mathbf{R}^T$
- 2D投影：$\mathbf{\Sigma}_{2D} = \mathbf{J}\mathbf{\Sigma}_{3D}\mathbf{J}^T$
- α-blending：$C = \sum_{i=1}^{N} c_i \alpha_i \prod_{j=1}^{i-1} (1 - \alpha_j)$

## 练习题

### 基础题

**练习10.1**：证明3D高斯函数在仿射变换下保持高斯性。具体地，若$\mathbf{y} = \mathbf{A}\mathbf{x} + \mathbf{b}$，证明变换后的函数仍是高斯。

*提示*：使用变量替换和雅可比行列式。

<details>
<summary>答案</summary>

设原始高斯为：
$$G(\mathbf{x}) = \frac{1}{(2\pi)^{3/2}|\mathbf{\Sigma}|^{1/2}} \exp\left(-\frac{1}{2}(\mathbf{x} - \mathbf{\mu})^T \mathbf{\Sigma}^{-1} (\mathbf{x} - \mathbf{\mu})\right)$$

进行变换$\mathbf{y} = \mathbf{A}\mathbf{x} + \mathbf{b}$，则$\mathbf{x} = \mathbf{A}^{-1}(\mathbf{y} - \mathbf{b})$。

雅可比行列式：$|\frac{\partial \mathbf{x}}{\partial \mathbf{y}}| = |\mathbf{A}^{-1}| = \frac{1}{|\mathbf{A}|}$

代入得：
$$G'(\mathbf{y}) = \frac{1}{|\mathbf{A}|} G(\mathbf{A}^{-1}(\mathbf{y} - \mathbf{b}))$$

展开指数项：
$$(\mathbf{A}^{-1}(\mathbf{y} - \mathbf{b}) - \mathbf{\mu})^T \mathbf{\Sigma}^{-1} (\mathbf{A}^{-1}(\mathbf{y} - \mathbf{b}) - \mathbf{\mu})$$
$$= (\mathbf{y} - (\mathbf{A}\mathbf{\mu} + \mathbf{b}))^T (\mathbf{A}^{-1})^T \mathbf{\Sigma}^{-1} \mathbf{A}^{-1} (\mathbf{y} - (\mathbf{A}\mathbf{\mu} + \mathbf{b}))$$

因此$G'(\mathbf{y})$是均值为$\mathbf{A}\mathbf{\mu} + \mathbf{b}$，协方差为$\mathbf{A}\mathbf{\Sigma}\mathbf{A}^T$的高斯。

</details>

**练习10.2**：推导2D高斯在像素$(u,v)$处的值对3D位置$\mathbf{\mu}$的梯度。

*提示*：使用链式法则，考虑投影变换的雅可比矩阵。

<details>
<summary>答案</summary>

设2D高斯为：
$$G_{2D}(u,v) = \exp\left(-\frac{1}{2}\mathbf{d}^T\mathbf{\Sigma}_{2D}^{-1}\mathbf{d}\right)$$

其中$\mathbf{d} = [u,v]^T - \mathbf{\mu}_{2D}$。

梯度计算：
$$\frac{\partial G_{2D}}{\partial \mathbf{\mu}} = \frac{\partial G_{2D}}{\partial \mathbf{\mu}_{2D}} \cdot \frac{\partial \mathbf{\mu}_{2D}}{\partial \mathbf{\mu}}$$

第一项：
$$\frac{\partial G_{2D}}{\partial \mathbf{\mu}_{2D}} = G_{2D} \cdot \mathbf{\Sigma}_{2D}^{-1}\mathbf{d}$$

第二项是投影的雅可比矩阵$\mathbf{J}$的转置。

最终：
$$\frac{\partial G_{2D}}{\partial \mathbf{\mu}} = G_{2D} \cdot \mathbf{J}^T \mathbf{\Sigma}_{2D}^{-1}\mathbf{d}$$

</details>

**练习10.3**：给定N个高斯，推导α-blending公式中透射率$T_i$的递推关系。

*提示*：考虑$T_{i+1}$与$T_i$的关系。

<details>
<summary>答案</summary>

透射率定义为：
$$T_i = \prod_{j=1}^{i-1} (1 - \alpha_j)$$

递推关系：
- 初始：$T_1 = 1$（没有前面的高斯）
- 递推：$T_{i+1} = T_i \cdot (1 - \alpha_i)$

这个递推关系允许高效的顺序计算，避免重复乘积。

验证：
$$T_{i+1} = \prod_{j=1}^{i} (1 - \alpha_j) = \prod_{j=1}^{i-1} (1 - \alpha_j) \cdot (1 - \alpha_i) = T_i \cdot (1 - \alpha_i)$$

</details>

### 挑战题

**练习10.4**：分析3D高斯溅射中深度排序的必要性。考虑如果不进行深度排序，而是使用无序α-blending，会产生什么误差？给出误差的数学界限。

*提示*：考虑两个高斯的情况，比较正确排序和错误排序的结果差异。

<details>
<summary>答案</summary>

考虑两个高斯，正确顺序（1在前，2在后）：
$$C_{correct} = c_1\alpha_1 + c_2\alpha_2(1-\alpha_1)$$

错误顺序：
$$C_{wrong} = c_2\alpha_2 + c_1\alpha_1(1-\alpha_2)$$

误差：
$$|C_{correct} - C_{wrong}| = |\alpha_1\alpha_2||c_1 - c_2|$$

对于N个高斯的一般情况，最坏情况误差界限：
$$\mathcal{E}_{max} \leq \sum_{i<j} \alpha_i\alpha_j|c_i - c_j|$$

当所有$\alpha_i$都很小时，误差近似为$\mathcal{O}(\alpha^2)$，这解释了为什么半透明物体的排序更重要。

</details>

**练习10.5**：设计一个自适应采样策略，在保持渲染质量的同时减少需要处理的高斯数量。给出理论分析和复杂度界限。

*提示*：考虑基于视角的重要性采样和层次化数据结构。

<details>
<summary>答案</summary>

自适应采样策略：

1. **重要性度量**：
   $$I_i = \alpha_i \cdot \frac{\text{area}_{2D,i}}{d_i^2} \cdot \cos\theta_i$$
   
   其中$d_i$是距离，$\theta_i$是视角。

2. **层次化结构**：使用八叉树，每个节点存储：
   - 边界框
   - 子高斯的聚合属性
   - 重要性上界

3. **剔除算法**：
   ```
   function cull(node, threshold):
     if node.importance_bound < threshold:
       return []
     if node.is_leaf:
       return filter(node.gaussians, I > threshold)
     else:
       return concat([cull(child, threshold) for child in node.children])
   ```

复杂度分析：
- 构建：$\mathcal{O}(N \log N)$
- 查询：$\mathcal{O}(K \log N)$，其中K是输出数量

误差界限：设剔除阈值为$\tau$，则渲染误差：
$$\mathcal{E} \leq \tau \cdot \text{pixel_count}$$

</details>

**练习10.6**：推导3D高斯溅射的信息论解释。将场景编码为高斯混合模型，分析其率失真性能。

*提示*：使用KL散度衡量重建误差，考虑高斯数量与精度的权衡。

<details>
<summary>答案</summary>

将场景视为概率分布$p(\mathbf{x})$，高斯混合近似为$q(\mathbf{x})$。

KL散度：
$$D_{KL}(p||q) = \int p(\mathbf{x}) \log \frac{p(\mathbf{x})}{q(\mathbf{x})} d\mathbf{x}$$

率失真函数：
$$R(D) = \min_{q: D_{KL}(p||q) \leq D} I(X;Y)$$

其中$I(X;Y)$是互信息，表示编码率。

对于N个高斯的混合模型：
- 编码率：$R = N \cdot (d_{pos} + d_{cov} + d_{color})$ bits
- 失真：$D \approx \mathcal{O}(1/N^{2/d})$（d维空间）

最优分配：使用Lloyd算法的变体，最小化：
$$\mathcal{L} = D_{KL}(p||q) + \lambda R$$

这给出了高斯数量与重建质量的理论权衡。

</details>

## 常见陷阱与错误 (Gotchas)

1. **数值不稳定**
   - 问题：协方差矩阵可能变成非正定
   - 解决：使用稳定的参数化（如Cholesky分解）
   - 调试：检查特征值，添加正则化项

2. **排序错误**
   - 问题：深度排序不正确导致渲染artifacts
   - 解决：使用稳定排序算法，处理深度相等情况
   - 调试：可视化深度图，检查排序一致性

3. **梯度爆炸**
   - 问题：小高斯的梯度可能非常大
   - 解决：梯度裁剪，自适应学习率
   - 调试：监控梯度范数，使用梯度累积

4. **内存溢出**
   - 问题：自适应分裂导致高斯数量爆炸
   - 解决：设置数量上限，aggressive修剪
   - 调试：跟踪高斯数量变化，profile内存使用

5. **性能瓶颈**
   - 问题：排序成为主要瓶颈
   - 解决：使用GPU友好的排序算法，考虑近似方法
   - 调试：性能profiling，识别热点

## 最佳实践检查清单

### 实现设计
- [ ] 选择合适的协方差参数化方式
- [ ] 实现数值稳定的投影变换
- [ ] 设计高效的数据结构（SoA vs AoS）
- [ ] 考虑GPU内存访问模式

### 优化策略
- [ ] 实现视锥体剔除
- [ ] 使用层次化数据结构
- [ ] 批量处理减少kernel调用
- [ ] 利用shared memory和texture cache

### 质量控制
- [ ] 监控高斯数量和分布
- [ ] 检查协方差矩阵条件数
- [ ] 验证梯度计算正确性
- [ ] 测试不同场景的泛化性

### 调试技巧
- [ ] 可视化单个高斯的贡献
- [ ] 检查深度排序结果
- [ ] 分析梯度流
- [ ] Profile性能瓶颈

### 扩展考虑
- [ ] 支持动态场景
- [ ] 集成环境光照
- [ ] 处理大规模场景
- [ ] 与其他表示混合使用
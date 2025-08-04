# 第9章：显式神经表示

在前面的章节中，我们探讨了将辐射场表示为连续函数的隐式神经表示（NeRF）。本章转向显式表示方法，这些方法直接在离散空间结构中存储场景信息。我们将展示如何将这些显式表示统一到体积渲染框架中，并分析它们相对于隐式方法的计算和存储权衡。

## 9.1 体素网格与八叉树

### 9.1.1 经典体素表示

最直接的显式表示是规则体素网格。对于分辨率为 $N^3$ 的体素网格，我们将空间离散化为：

$$\mathbf{x}_{ijk} = \mathbf{x}_{\min} + \left(\frac{i}{N}, \frac{j}{N}, \frac{k}{N}\right) \odot (\mathbf{x}_{\max} - \mathbf{x}_{\min})$$

其中 $i,j,k \in \{0,1,...,N-1\}$，$\odot$ 表示逐元素乘积。

每个体素存储局部体积属性：
- 密度：$\sigma_{ijk}$
- 辐射：$\mathbf{c}_{ijk}$ 或更一般的球谐系数

体积渲染方程在离散情况下变为：

$$C(\mathbf{r}) = \sum_{n=1}^{N_{\text{samples}}} T_n \alpha_n \mathbf{c}_n$$

其中：
- $T_n = \exp\left(-\sum_{m=1}^{n-1} \sigma_m \delta_m\right)$ 是透射率
- $\alpha_n = 1 - \exp(-\sigma_n \delta_n)$ 是不透明度
- $\delta_n$ 是样本间距

### 9.1.2 三线性插值

为避免块状伪影，我们在体素中心之间进行三线性插值：

$$f(\mathbf{x}) = \sum_{i=0}^{1}\sum_{j=0}^{1}\sum_{k=0}^{1} f_{ijk} \prod_{d \in \{x,y,z\}} (1-|x_d - x_{d,ijk}|)$$

这保证了 $C^0$ 连续性但在体素边界处导数不连续。

### 9.1.3 稀疏八叉树

对于大多数场景，空间是稀疏的。八叉树通过自适应细分提供了高效的表示：

```
OctreeNode {
    如果是叶子节点：存储 (σ, c)
    否则：指向8个子节点的指针
}
```

八叉树遍历的期望复杂度为 $O(\log N)$，其中 $N$ 是叶节点数。

存储复杂度从密集网格的 $O(N^3)$ 降低到 $O(N_{\text{occupied}})$，其中 $N_{\text{occupied}}$ 是非空体素的数量。

## 9.2 Plenoxels：球谐函数体素

### 9.2.1 球谐函数回顾

球谐函数 $Y_{\ell}^m(\theta, \phi)$ 形成了球面上平方可积函数的正交基：

$$Y_{\ell}^m(\theta, \phi) = \sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}} P_{\ell}^m(\cos\theta) e^{im\phi}$$

其中 $P_{\ell}^m$ 是关联勒让德多项式。

对于实值函数，我们使用实球谐函数：
- $Y_{\ell m} = \sqrt{2} \Re(Y_{\ell}^m)$ 当 $m > 0$
- $Y_{\ell 0} = Y_{\ell}^0$
- $Y_{\ell m} = \sqrt{2} \Im(Y_{\ell}^{|m|})$ 当 $m < 0$

### 9.2.2 Plenoxels 表示

Plenoxels 在每个体素存储：
1. 密度 $\sigma \in \mathbb{R}^+$
2. 球谐系数 $\mathbf{k} \in \mathbb{R}^{(L+1)^2 \times 3}$ 用于 RGB 颜色

视角相关的颜色通过球谐展开计算：

$$\mathbf{c}(\mathbf{d}) = \text{sigmoid}\left(\sum_{\ell=0}^{L} \sum_{m=-\ell}^{\ell} \mathbf{k}_{\ell m} Y_{\ell m}(\mathbf{d})\right)$$

其中 $L$ 是最大阶数（通常 $L=2$ 对应9个系数）。

### 9.2.3 直接优化

与NeRF不同，Plenoxels直接优化体素参数，无需神经网络：

$$\mathcal{L} = \sum_{\mathbf{r} \in \mathcal{R}} \|\hat{C}(\mathbf{r}) - C(\mathbf{r})\|^2 + \lambda_{\text{TV}} \mathcal{L}_{\text{TV}} + \lambda_{\text{sparse}} \|\sigma\|_1$$

其中：
- $\mathcal{L}_{\text{TV}}$ 是总变差正则化：$\sum_{i,j,k} \|\nabla \sigma_{ijk}\|_2$
- 稀疏性项促进空体素

优化使用标准梯度下降，梯度通过可微体积渲染反向传播。

## 9.3 TensoRF：张量分解辐射场

### 9.3.1 张量分解基础

TensoRF 将 4D 辐射场张量 $\mathcal{T} \in \mathbb{R}^{X \times Y \times Z \times C}$ 分解为低秩分量。

**CP分解（CANDECOMP/PARAFAC）**：
$$\mathcal{T} = \sum_{r=1}^{R} \mathbf{u}_r \otimes \mathbf{v}_r \otimes \mathbf{w}_r \otimes \mathbf{s}_r$$

其中 $\otimes$ 是外积，$R$ 是秩。

**向量-矩阵（VM）分解**：
$$\mathcal{T} = \sum_{r=1}^{R} (\mathbf{M}_{r,1} \circ \mathbf{v}_{r,1}) \otimes \mathbf{b}_{r,1} + (\mathbf{M}_{r,2} \circ \mathbf{v}_{r,2}) \otimes \mathbf{b}_{r,2} + (\mathbf{M}_{r,3} \circ \mathbf{v}_{r,3}) \otimes \mathbf{b}_{r,3}$$

其中 $\mathbf{M}_{r,i} \in \mathbb{R}^{D_1 \times D_2}$ 是矩阵，$\mathbf{v}_{r,i} \in \mathbb{R}^{D_3}$ 是向量，$\circ$ 表示模式积。

### 9.3.2 密度和外观分解

TensoRF 分别建模密度和外观：

**密度场**：
$$\sigma(\mathbf{x}) = \sum_{r=1}^{R_\sigma} \prod_{d=1}^{3} f_d^{\sigma,r}(x_d)$$

**外观场**：
$$\mathbf{c}(\mathbf{x}, \mathbf{d}) = \mathcal{F}_\theta\left(\sum_{r=1}^{R_c} \prod_{d=1}^{3} f_d^{c,r}(x_d), \mathbf{d}\right)$$

其中 $\mathcal{F}_\theta$ 是小型MLP解码器。

### 9.3.3 存储和计算复杂度

存储需求从 $O(N^3)$ 降低到 $O(RN)$：
- CP分解：$3RN$ 参数
- VM分解：$R(2N^2 + N)$ 参数

计算复杂度：
- 查询：$O(R)$ 而非 $O(1)$
- 但 $R \ll N$，因此实践中更高效

## 9.4 低秩与稀疏表示

### 9.4.1 数学基础

辐射场的低秩结构源于场景的内在规律性。考虑奇异值分解（SVD）：

$$\mathbf{A} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^T$$

对于秩为 $r$ 的矩阵，最优 $k$ 秩近似（在Frobenius范数下）是：
$$\mathbf{A}_k = \sum_{i=1}^{k} \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$

误差界：
$$\|\mathbf{A} - \mathbf{A}_k\|_F = \sqrt{\sum_{i=k+1}^{r} \sigma_i^2}$$

### 9.4.2 稀疏性诱导

许多场景在变换域中是稀疏的。考虑小波变换 $\mathcal{W}$：

$$\mathbf{f} = \mathcal{W}^{-1}(\mathbf{w})$$

其中 $\mathbf{w}$ 是稀疏系数向量。

通过 $\ell_1$ 正则化促进稀疏性：
$$\min_{\mathbf{w}} \|\mathcal{Y} - \mathcal{A}\mathcal{W}^{-1}(\mathbf{w})\|_2^2 + \lambda\|\mathbf{w}\|_1$$

这是基追踪去噪问题，可通过迭代收缩阈值算法（ISTA）求解。

### 9.4.3 压缩感知视角

当测量矩阵 $\mathcal{A}$ 满足限制等距性质（RIP）时，我们可以从 $O(s\log N)$ 个测量中恢复 $s$-稀疏信号。

对于辐射场，"测量"是渲染的像素，稀疏性在适当选择的基中实现。

## 9.5 显式与隐式权衡

### 9.5.1 内存占用

**隐式（NeRF）**：
- 固定大小：$O(W \times L)$，其中 $W$ 是宽度，$L$ 是层数
- 与场景复杂度无关
- 典型值：~5MB

**显式（体素）**：
- 密集：$O(N^3 \times C)$
- 稀疏：$O(N_{\text{occupied}} \times C)$
- 低秩：$O(R \times N \times C)$
- 典型值：100MB-10GB

### 9.5.2 渲染速度

**查询复杂度**：
- 隐式：$O(L)$ 神经网络前向传播
- 显式：$O(1)$ 数组查找或 $O(\log N)$ 树遍历

**缓存效率**：
- 显式方法具有更好的空间局部性
- 可以利用GPU纹理硬件

### 9.5.3 训练特性

**收敛速度**：
- 显式：直接优化，收敛更快
- 隐式：通过梯度下降，需要更多迭代

**表示能力**：
- 隐式：平滑先验，更好的插值
- 显式：可以表示高频细节，但可能过拟合

### 9.5.4 统一视角

两种方法都解决相同的体积渲染方程：

$$C(\mathbf{r}) = \int_{t_n}^{t_f} T(t)\sigma(t)\mathbf{c}(t,\mathbf{d})dt$$

区别在于如何表示 $\sigma$ 和 $\mathbf{c}$：
- 隐式：$\sigma = f_\theta(\gamma(\mathbf{x}))$
- 显式：$\sigma = \text{interp}(\{\sigma_{ijk}\})$

选择取决于具体应用需求。

## 本章小结

本章探讨了神经辐射场的显式表示方法，展示了如何将离散空间结构统一到体积渲染框架中：

**核心概念**：
1. **体素网格**：最直接的显式表示，通过三线性插值实现连续性
2. **Plenoxels**：使用球谐函数编码视角相关效果，直接优化无需神经网络
3. **TensoRF**：通过张量分解实现紧凑表示，将存储从 $O(N^3)$ 降至 $O(RN)$
4. **低秩结构**：利用场景的内在规律性进行压缩

**关键公式**：
- 离散体积渲染：$C(\mathbf{r}) = \sum_{n=1}^{N} T_n \alpha_n \mathbf{c}_n$
- 球谐展开：$\mathbf{c}(\mathbf{d}) = \sum_{\ell,m} \mathbf{k}_{\ell m} Y_{\ell m}(\mathbf{d})$
- 张量分解：$\mathcal{T} = \sum_{r=1}^{R} \mathbf{u}_r \otimes \mathbf{v}_r \otimes \mathbf{w}_r \otimes \mathbf{s}_r$

**实践要点**：
- 显式方法提供更快的查询速度但需要更多内存
- 低秩和稀疏表示在保持质量的同时大幅减少存储
- 选择显式或隐式取决于速度、内存和质量的权衡

## 练习题

### 基础题

**练习 9.1：八叉树遍历复杂度**
证明在平衡八叉树中，从根到叶的平均遍历深度是 $O(\log N)$，其中 $N$ 是叶节点数。

<details>
<summary>提示</summary>
考虑深度为 $d$ 的完全八叉树最多包含多少叶节点。
</details>

<details>
<summary>答案</summary>

深度为 $d$ 的完全八叉树最多有 $8^d$ 个叶节点。

如果有 $N$ 个叶节点，则最小深度 $d$ 满足：
$$8^{d-1} < N \leq 8^d$$

取对数：
$$d-1 < \log_8 N \leq d$$

因此 $d = \lceil \log_8 N \rceil = O(\log N)$。

平均情况下，假设叶节点均匀分布，平均深度约为 $\log_8 N$。
</details>

**练习 9.2：球谐函数正交性**
证明前两阶球谐函数在单位球面上正交：
$$\int_{S^2} Y_{\ell m}(\mathbf{d}) Y_{\ell' m'}(\mathbf{d}) d\mathbf{d} = \delta_{\ell\ell'}\delta_{mm'}$$

<details>
<summary>提示</summary>
使用球坐标 $(\theta, \phi)$ 并利用关联勒让德多项式的正交性。
</details>

<details>
<summary>答案</summary>

在球坐标中：
$$\int_{S^2} Y_{\ell m} Y_{\ell' m'} d\mathbf{d} = \int_0^{2\pi} \int_0^{\pi} Y_{\ell m}(\theta,\phi) Y_{\ell' m'}(\theta,\phi) \sin\theta d\theta d\phi$$

代入球谐函数定义：
$$= \sqrt{\frac{(2\ell+1)(2\ell'+1)}{16\pi^2} \frac{(\ell-m)!(\ell'-m')!}{(\ell+m)!(\ell'+m')!}} \times$$
$$\int_0^{2\pi} e^{i(m-m')\phi} d\phi \int_0^{\pi} P_{\ell}^m(\cos\theta) P_{\ell'}^{m'}(\cos\theta) \sin\theta d\theta$$

方位角积分给出 $2\pi\delta_{mm'}$。

对于极角积分，令 $x = \cos\theta$：
$$\int_{-1}^{1} P_{\ell}^m(x) P_{\ell'}^{m}(x) dx = \frac{2}{2\ell+1} \frac{(\ell+m)!}{(\ell-m)!} \delta_{\ell\ell'}$$

组合所有项得到 $\delta_{\ell\ell'}\delta_{mm'}$。
</details>

**练习 9.3：三线性插值连续性**
证明三线性插值在体素边界处是 $C^0$ 连续但不是 $C^1$ 连续。

<details>
<summary>提示</summary>
检查跨越体素边界时函数值和导数的行为。
</details>

<details>
<summary>答案</summary>

考虑1D情况的线性插值：
$$f(x) = f_0(1-x) + f_1 x, \quad x \in [0,1]$$

在边界 $x=1$ 处，从当前体素：$f(1^-) = f_1$
从相邻体素（新坐标系）：$f(0^+) = f_1$

因此函数值连续（$C^0$）。

但导数：
- 当前体素：$f'(1^-) = f_1 - f_0$
- 相邻体素：$f'(0^+) = f_2 - f_1$

一般情况下 $f'(1^-) \neq f'(0^+)$，所以不是 $C^1$ 连续。

三线性情况类似，在每个维度独立应用。
</details>

### 挑战题

**练习 9.4：张量分解误差界**
给定秩为 $r$ 的 3 阶张量 $\mathcal{T} \in \mathbb{R}^{n \times n \times n}$，推导 CP 分解的最优 $k$ 项近似误差界。

<details>
<summary>提示</summary>
类比矩阵 SVD 的 Eckart-Young 定理，但注意张量情况更复杂。
</details>

<details>
<summary>答案</summary>

与矩阵不同，张量的最优低秩近似不一定由截断的CP分解给出。但我们可以建立界限。

设 $\mathcal{T}_k$ 是最优 $k$ 秩近似，$\tilde{\mathcal{T}}_k$ 是 $k$ 项CP分解。则：

$$\|\mathcal{T} - \tilde{\mathcal{T}}_k\|_F \leq \sqrt{1 + \epsilon_k} \|\mathcal{T} - \mathcal{T}_k\|_F$$

其中 $\epsilon_k$ 依赖于张量的条件数。

对于"良态"张量（如低秩加噪声），CP分解给出近最优近似：
$$\|\mathcal{T} - \tilde{\mathcal{T}}_k\|_F \leq C\sqrt{\sum_{i=k+1}^r \lambda_i^2}$$

其中 $\lambda_i$ 是张量的"奇异值"（通过高阶SVD定义），$C$ 是常数。
</details>

**练习 9.5：稀疏八叉树的期望存储**
考虑在单位立方体中均匀分布的 $M$ 个点，每个点占据半径 $\epsilon$ 的球。推导存储这些点所需的八叉树节点数的期望值。

<details>
<summary>提示</summary>
计算不同深度级别的期望占用体素数。
</details>

<details>
<summary>答案</summary>

在深度 $d$ 处，体素大小为 $2^{-d}$。

一个半径 $\epsilon$ 的球与体素相交的条件是球心到体素的距离 $\leq \epsilon + \frac{\sqrt{3}}{2} \cdot 2^{-d}$。

在深度 $d$ 的期望占用体素数：
$$E[N_d] = \min\left(8^d, M \cdot V_{\text{intersection}}(d) \cdot 8^d\right)$$

其中 $V_{\text{intersection}}(d)$ 是相交体积。

当 $2^{-d} \gg \epsilon$：一个球占据 $O(1)$ 个体素
当 $2^{-d} \ll \epsilon$：一个球占据 $O((\epsilon \cdot 2^d)^3)$ 个体素

最优深度约为 $d^* = \log_2(1/\epsilon)$，总存储：
$$E[N_{\text{total}}] = O(M) + O(M(\epsilon \cdot 2^{d^*})^3) = O(M)$$

因此稀疏八叉树存储与点数成线性关系。
</details>

**练习 9.6：Plenoxels vs NeRF 收敛分析**
分析为什么直接优化体素参数（Plenoxels）比通过神经网络（NeRF）收敛更快。考虑优化landscape和梯度流。

<details>
<summary>提示</summary>
比较参数到输出的映射复杂度和梯度传播路径。
</details>

<details>
<summary>答案</summary>

**Plenoxels 优化景观**：
- 参数直接映射到输出：$\theta_{ijk} \rightarrow C(\mathbf{r})$
- 凸优化问题（对每个体素独立）
- 梯度：$\frac{\partial \mathcal{L}}{\partial \theta_{ijk}} = \sum_{\mathbf{r} \text{ through } ijk} \frac{\partial \mathcal{L}}{\partial C(\mathbf{r})} \frac{\partial C(\mathbf{r})}{\partial \theta_{ijk}}$

**NeRF 优化景观**：
- 非凸，多个局部最小值
- 梯度通过多层网络：$\frac{\partial \mathcal{L}}{\partial \mathbf{W}} = \frac{\partial \mathcal{L}}{\partial C} \prod_{l} \frac{\partial f_l}{\partial f_{l-1}} \frac{\partial f_1}{\partial \mathbf{W}}$
- 梯度消失/爆炸问题

**收敛速度差异**：
1. **条件数**：Plenoxels 的 Hessian 对角占优，条件数更好
2. **局部性**：每个体素参数主要影响局部区域
3. **并行性**：体素更新可以并行进行

理论上，Plenoxels 的收敛速度为 $O(1/t)$（凸优化），而 NeRF 为 $O(1/\sqrt{t})$（非凸）。
</details>

## 常见陷阱与错误

1. **内存爆炸**：密集体素网格的 $O(N^3)$ 增长
   - 调试技巧：始终从粗分辨率开始，监控内存使用

2. **混叠伪影**：采样不足导致的摩尔纹
   - 调试技巧：使用抗混叠滤波器，增加超采样

3. **数值不稳定**：张量分解中的退化情况
   - 调试技巧：添加小的正则化项 $\epsilon I$，使用稳定的初始化

4. **边界伪影**：体素边界的不连续性
   - 调试技巧：扩展边界体素，使用更高阶插值

5. **稀疏性与质量权衡**：过度稀疏化导致细节丢失
   - 调试技巧：使用自适应阈值，保留重要区域

6. **球谐截断误差**：低阶球谐不能表示高频
   - 调试技巧：分析频谱内容，自适应选择阶数

## 最佳实践检查清单

设计显式神经表示时，确保：

- [ ] **分辨率选择**：基于场景尺度和可用内存选择合适的体素分辨率
- [ ] **稀疏性策略**：实现高效的稀疏数据结构（八叉树、哈希表）
- [ ] **插值方法**：选择适当的插值（三线性、三次样条）平衡质量和速度
- [ ] **压缩技术**：应用低秩分解或量化减少内存占用
- [ ] **正则化设计**：包含空间平滑性和稀疏性约束
- [ ] **初始化策略**：使用场景先验（如视觉外壳）初始化占用
- [ ] **自适应细化**：实现由误差度量引导的渐进式细化
- [ ] **缓存优化**：组织数据布局以最大化GPU缓存命中率
- [ ] **数值稳定性**：在所有计算中处理退化情况
- [ ] **评估指标**：测量渲染质量、速度和内存使用的权衡
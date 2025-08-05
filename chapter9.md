# 第9章：显式神经表示

在前面的章节中，我们探讨了将辐射场表示为连续函数的隐式神经表示（NeRF）。本章转向显式表示方法，这些方法直接在离散空间结构中存储场景信息。显式表示的核心优势在于查询效率和直观的空间组织，而主要挑战在于内存需求和离散化带来的伪影。我们将展示如何将这些显式表示统一到体积渲染框架中，并分析它们相对于隐式方法的计算和存储权衡。

## 9.1 体素网格与八叉树

### 9.1.1 经典体素表示

最直接的显式表示是规则体素网格。对于分辨率为 $N^3$ 的体素网格，我们将空间离散化为：

$$\mathbf{x}_{ijk} = \mathbf{x}_{\min} + \left(\frac{i}{N}, \frac{j}{N}, \frac{k}{N}\right) \odot (\mathbf{x}_{\max} - \mathbf{x}_{\min})$$

其中 $i,j,k \in \{0,1,...,N-1\}$，$\odot$ 表示逐元素乘积。

每个体素存储局部体积属性：
- 密度：$\sigma_{ijk} \in \mathbb{R}^+$
- 辐射：$\mathbf{c}_{ijk} \in \mathbb{R}^3$ 或更一般的球谐系数 $\mathbf{k}_{ijk} \in \mathbb{R}^{(L+1)^2 \times 3}$
- 可选：法线 $\mathbf{n}_{ijk}$、材质属性等

体积渲染方程在离散情况下变为黎曼和近似：

$$C(\mathbf{r}) = \sum_{n=1}^{N_{\text{samples}}} T_n \alpha_n \mathbf{c}_n$$

其中：
- $T_n = \exp\left(-\sum_{m=1}^{n-1} \sigma_m \delta_m\right)$ 是透射率
- $\alpha_n = 1 - \exp(-\sigma_n \delta_n)$ 是不透明度  
- $\delta_n$ 是样本间距

这个离散化引入了 $O(\delta^2)$ 的误差，根据中点法则：
$$\left|\int_a^b f(t)dt - \sum_{n} f(t_n)\delta\right| \leq \frac{(b-a)\delta^2}{24}\max_{t \in [a,b]}|f''(t)|$$

**采样策略与误差分析**

为了减少离散化误差，可以采用更高阶的数值积分方法：

1. **梯形法则**：使用端点值的平均，误差降至 $O(\delta^2)$
   $$\int_{t_i}^{t_{i+1}} f(t)dt \approx \frac{\delta}{2}[f(t_i) + f(t_{i+1})]$$

2. **Simpson法则**：使用二次插值，误差为 $O(\delta^4)$
   $$\int_{t_i}^{t_{i+2}} f(t)dt \approx \frac{\delta}{3}[f(t_i) + 4f(t_{i+1}) + f(t_{i+2})]$$

3. **自适应采样**：在高频区域（边界、细节）增加采样密度
   - 基于梯度的细化：$\delta_{\text{local}} = \delta_0 / (1 + \|\nabla\sigma\|)$
   - 基于曲率的细化：考虑二阶导数 $\nabla^2\sigma$

**体素遍历算法**

光线与体素网格的相交使用3D DDA（Digital Differential Analyzer）算法：

1. **初始化**：计算光线进入边界框的参数 $t_{\text{min}}$
2. **步进计算**：
   $$t_{\text{next},d} = t_{\text{current}} + \frac{\text{sign}(\mathbf{d}_d)}{\mathbf{d}_d} \cdot \text{voxelSize}_d$$
   其中 $d \in \{x,y,z\}$
3. **选择最近边界**：$t_{\text{next}} = \min(t_{\text{next},x}, t_{\text{next},y}, t_{\text{next},z})$
4. **更新体素索引**：根据最小 $t$ 对应的维度递增

该算法的复杂度为 $O(N)$，其中 $N$ 是光线穿过的体素数，最坏情况下为 $O(N_{\text{grid}})$。

### 9.1.2 三线性插值

为避免块状伪影，我们在体素中心之间进行三线性插值。给定查询点 $\mathbf{x}$，首先找到包含它的体素 $(i,j,k)$，然后计算局部坐标：

$$\mathbf{u} = \left(\frac{x - x_{i}}{x_{i+1} - x_i}, \frac{y - y_{j}}{y_{j+1} - y_j}, \frac{z - z_{k}}{z_{k+1} - z_k}\right)$$

三线性插值公式为：

$$f(\mathbf{x}) = \sum_{i'=0}^{1}\sum_{j'=0}^{1}\sum_{k'=0}^{1} f_{i+i',j+j',k+k'} \prod_{d \in \{x,y,z\}} \begin{cases}
1-u_d & \text{if } d'=0 \\
u_d & \text{if } d'=1
\end{cases}$$

这等价于三次线性插值的嵌套：
$$f(\mathbf{x}) = \text{lerp}_z\left(\text{lerp}_y\left(\text{lerp}_x(f_{000}, f_{100}, u_x), \text{lerp}_x(f_{010}, f_{110}, u_x), u_y\right), \ldots, u_z\right)$$

这保证了 $C^0$ 连续性但在体素边界处导数不连续，导致法线计算时出现阶跃。

**插值权重的几何解释**

三线性插值的权重具有直观的几何意义——每个顶点的贡献与查询点到对角顶点形成的子体积成正比：

$$w_{ijk} = \prod_{d \in \{x,y,z\}} \begin{cases}
(1-u_d) & \text{if bit}_d(ijk) = 0 \\
u_d & \text{if bit}_d(ijk) = 1
\end{cases}$$

其中 $\text{bit}_d(ijk)$ 提取索引的第 $d$ 位。

**梯度计算**

尽管函数值连续，梯度在体素边界不连续。体素内部的梯度为：

$$\nabla f(\mathbf{x}) = \begin{pmatrix}
\frac{\partial f}{\partial x} \\
\frac{\partial f}{\partial y} \\
\frac{\partial f}{\partial z}
\end{pmatrix} = \begin{pmatrix}
\frac{1}{\Delta x} \sum_{j',k'} (f_{1,j',k'} - f_{0,j',k'}) w_{j',k'}^{yz} \\
\frac{1}{\Delta y} \sum_{i',k'} (f_{i',1,k'} - f_{i',0,k'}) w_{i',k'}^{xz} \\
\frac{1}{\Delta z} \sum_{i',j'} (f_{i',j',1} - f_{i',j',0}) w_{i',j'}^{xy}
\end{pmatrix}$$

其中 $w^{yz}$ 表示仅在 $y,z$ 维度的权重积。

**高阶插值方案**

为获得更平滑的结果，可以使用高阶插值：

1. **三次样条插值**：保证 $C^2$ 连续性
   $$f(\mathbf{x}) = \sum_{i=-1}^{2}\sum_{j=-1}^{2}\sum_{k=-1}^{2} f_{i+i_0,j+j_0,k+k_0} \cdot B(u_x-i)B(u_y-j)B(u_z-k)$$
   其中 $B(t)$ 是三次B样条基函数

2. **Catmull-Rom插值**：通过4个点的三次插值
   $$f(u) = \frac{1}{2}\begin{pmatrix} u^3 & u^2 & u & 1 \end{pmatrix} \begin{pmatrix}
   -1 & 3 & -3 & 1 \\
   2 & -5 & 4 & -1 \\
   -1 & 0 & 1 & 0 \\
   0 & 2 & 0 & 0
   \end{pmatrix} \begin{pmatrix} f_{-1} \\ f_0 \\ f_1 \\ f_2 \end{pmatrix}$$

3. **Kaiser-Bessel窗插值**：最小化频谱泄漏
   $$w(r) = \begin{cases}
   \frac{I_0(\beta\sqrt{1-r^2})}{I_0(\beta)} & |r| \leq 1 \\
   0 & |r| > 1
   \end{cases}$$
   其中 $I_0$ 是修正贝塞尔函数，$\beta$ 控制窗口形状

### 9.1.3 稀疏八叉树

对于大多数场景，空间是稀疏的——大部分体积为空。八叉树通过自适应细分提供了高效的表示。每个节点代表一个立方体区域，可以是：
- **叶节点**：存储实际数据 $(\sigma, \mathbf{c})$
- **内部节点**：包含8个子节点的指针

形式化定义：
```
struct OctreeNode {
    Vector3 center;     // 节点中心
    float halfWidth;    // 半边长
    union {
        struct { float sigma; Vector3 color; } leaf;
        OctreeNode* children[8];
    };
    bool isLeaf;
}
```

子节点索引通过位操作确定：
$$\text{childIndex} = (x > x_c) + 2(y > y_c) + 4(z > z_c)$$

八叉树遍历算法：
1. 从根节点开始
2. 计算射线与当前节点的交点
3. 如果是叶节点，采样并累积
4. 否则递归访问相交的子节点

遍历的期望复杂度为 $O(\log N)$，其中 $N$ 是叶节点数。证明基于平衡树的高度为 $\lceil \log_8 N \rceil$。

存储复杂度从密集网格的 $O(N^3)$ 降低到 $O(N_{\text{occupied}})$，其中 $N_{\text{occupied}}$ 是非空体素的数量。实际压缩率取决于场景稀疏性。

**高效光线-八叉树相交**

优化的遍历算法利用空间相干性：

1. **参数化相交测试**：对于光线 $\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$，计算与节点边界的进入/退出参数：
   $$t_{\text{near}} = \max_{i \in \{x,y,z\}} \frac{(\mathbf{c}_i - h_i - \mathbf{o}_i)}{\mathbf{d}_i}, \quad t_{\text{far}} = \min_{i \in \{x,y,z\}} \frac{(\mathbf{c}_i + h_i - \mathbf{o}_i)}{\mathbf{d}_i}$$
   其中 $\mathbf{c}$ 是节点中心，$h$ 是半边长

2. **子节点访问顺序**：基于光线方向确定遍历顺序，优先访问更近的子节点：
   ```
   childOrder = ComputeTraversalOrder(rayDirection)
   for i in childOrder:
       if RayIntersectsChild(ray, children[i]):
           TraverseChild(ray, children[i])
   ```

3. **早期终止**：当累积透明度低于阈值时停止遍历：
   $$T < \epsilon \Rightarrow \text{terminate}$$

**自适应细分准则**

决定何时细分节点的准则包括：

1. **基于密度的细分**：当节点内密度变化超过阈值
   $$\text{var}(\sigma) > \tau_{\sigma} \Rightarrow \text{subdivide}$$

2. **基于梯度的细分**：在边界区域增加分辨率
   $$\|\nabla\sigma\| > \tau_{\nabla} \Rightarrow \text{subdivide}$$

3. **视点相关细分**：基于投影大小动态调整
   $$\frac{\text{nodeSize}}{\|\mathbf{x} - \mathbf{o}_{\text{camera}}\|} > \tau_{\text{pixel}} \Rightarrow \text{subdivide}$$

**内存优化技术**

1. **线性八叉树**：使用Morton编码将树结构线性化
   $$\text{MortonCode}(x,y,z) = \sum_{i=0}^{L-1} 2^{3i}[(x_i) + 2(y_i) + 4(z_i)]$$
   其中 $x_i, y_i, z_i$ 是坐标的第 $i$ 位

2. **节点池分配**：预分配大块内存，减少碎片
   ```
   nodePool = AllocateBlock(MAX_NODES * sizeof(OctreeNode))
   freeList = InitializeFreeList(nodePool)
   ```

3. **LOD（细节层次）**：存储多分辨率版本
   - 粗糙级别用于远处物体
   - 精细级别用于近处细节
   - 级别间的平滑过渡

**并行构建算法**

利用GPU并行性加速八叉树构建：

1. **自顶向下构建**：
   - 并行计算每个节点的占用状态
   - 使用原子操作分配子节点
   - 复杂度：$O(\log N)$ 并行步骤

2. **自底向上构建**：
   - 首先对所有点进行空间排序（Morton编码）
   - 并行识别节点边界
   - 自底向上合并构建树
   - 复杂度：$O(N \log N)$ 工作量，$O(\log N)$ 深度

## 9.2 Plenoxels：球谐函数体素

### 9.2.1 球谐函数回顾

球谐函数 $Y_{\ell}^m(\theta, \phi)$ 形成了球面 $S^2$ 上平方可积函数的完备正交基。它们是拉普拉斯-贝尔特拉米算子在球面上的本征函数：

$$\nabla^2_{S^2} Y_{\ell}^m = -\ell(\ell+1) Y_{\ell}^m$$

复球谐函数定义为：

$$Y_{\ell}^m(\theta, \phi) = \sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}} P_{\ell}^m(\cos\theta) e^{im\phi}$$

其中 $P_{\ell}^m$ 是关联勒让德多项式，满足：
$$P_{\ell}^m(x) = (-1)^m (1-x^2)^{m/2} \frac{d^m}{dx^m} P_{\ell}(x)$$

对于计算机图形学中的实值函数，我们使用实球谐函数基：
- $Y_{\ell m} = \sqrt{2} \Re(Y_{\ell}^m) = \sqrt{2} N_{\ell}^m P_{\ell}^m(\cos\theta) \cos(m\phi)$ 当 $m > 0$
- $Y_{\ell 0} = N_{\ell}^0 P_{\ell}(\cos\theta)$ 当 $m = 0$
- $Y_{\ell m} = \sqrt{2} \Im(Y_{\ell}^{|m|}) = \sqrt{2} N_{\ell}^{|m|} P_{\ell}^{|m|}(\cos\theta) \sin(|m|\phi)$ 当 $m < 0$

前几阶的具体形式：
- $Y_{00} = \frac{1}{2\sqrt{\pi}}$ （常数项）
- $Y_{1,-1} = \sqrt{\frac{3}{4\pi}} y$，$Y_{10} = \sqrt{\frac{3}{4\pi}} z$，$Y_{11} = \sqrt{\frac{3}{4\pi}} x$
- $Y_{2,-2} = \frac{1}{2}\sqrt{\frac{15}{\pi}} xy$，等等

### 9.2.2 Plenoxels 表示

Plenoxels（"Plenoptic Voxels"的缩写）是一种纯显式优化方法，完全避免了神经网络。在每个体素中存储：
1. **体积密度** $\sigma \in \mathbb{R}^+$：控制不透明度
2. **球谐系数** $\mathbf{k} \in \mathbb{R}^{(L+1)^2 \times 3}$：编码视角相关的RGB颜色

视角相关的辐射通过球谐展开计算：

$$\mathbf{c}(\mathbf{x}, \mathbf{d}) = \text{sigmoid}\left(\sum_{\ell=0}^{L} \sum_{m=-\ell}^{\ell} \mathbf{k}_{\ell m}(\mathbf{x}) Y_{\ell m}(\mathbf{d})\right)$$

其中：
- $L$ 是最大阶数（$L=0$：1个系数，$L=1$：4个系数，$L=2$：9个系数）
- sigmoid函数确保颜色值在 $[0,1]$ 范围内
- 高阶项捕获高频视角变化（镜面反射等）

### 9.2.3 直接优化

与NeRF不同，Plenoxels直接优化体素参数，无需神经网络。这导致了一个大规模但结构化的优化问题。损失函数包含数据项和正则化项：

$$\mathcal{L} = \mathcal{L}_{\text{data}} + \lambda_{\text{TV}} \mathcal{L}_{\text{TV}} + \lambda_{\text{sparse}} \mathcal{L}_{\text{sparse}}$$

**数据项**衡量渲染图像与真实图像的差异：
$$\mathcal{L}_{\text{data}} = \sum_{\mathbf{r} \in \mathcal{R}} \|\hat{C}(\mathbf{r}) - C(\mathbf{r})\|^2$$

**总变差（TV）正则化**促进空间平滑性：
$$\mathcal{L}_{\text{TV}} = \sum_{i,j,k} \left(\|\sigma_{i+1,j,k} - \sigma_{i,j,k}\|^2 + \|\sigma_{i,j+1,k} - \sigma_{i,j,k}\|^2 + \|\sigma_{i,j,k+1} - \sigma_{i,j,k}\|^2\right)$$

这等价于各向异性TV范数，防止密度场出现噪声。

**稀疏性正则化**鼓励空体素：
$$\mathcal{L}_{\text{sparse}} = \sum_{i,j,k} \log(1 + \sigma_{ijk}^2/\epsilon^2)$$

这是 $\ell_1$ 范数的平滑近似，避免了不可微点。

优化使用Adam优化器，学习率调度遵循：
$$\eta(t) = \eta_0 \cdot 0.1^{t/T}$$

梯度通过可微体积渲染反向传播：
$$\frac{\partial \mathcal{L}}{\partial \sigma_{ijk}} = \sum_{\mathbf{r} \text{ through } (i,j,k)} \frac{\partial \mathcal{L}}{\partial C(\mathbf{r})} \frac{\partial C(\mathbf{r})}{\partial \sigma_{ijk}}$$

## 9.3 TensoRF：张量分解辐射场

### 9.3.1 张量分解基础

TensoRF的核心洞察是辐射场具有低秩结构——场景中的规律性（平面、对称性、重复纹理）导致信息冗余。通过张量分解，我们可以用紧凑的因子表示高维数据。

考虑4D辐射场张量 $\mathcal{T} \in \mathbb{R}^{X \times Y \times Z \times C}$，其中前三维是空间坐标，第四维是特征通道。

**CP分解（CANDECOMP/PARAFAC）**将张量表示为秩-1张量的和：
$$\mathcal{T} = \sum_{r=1}^{R} \lambda_r \cdot \mathbf{u}_r \otimes \mathbf{v}_r \otimes \mathbf{w}_r \otimes \mathbf{s}_r$$

其中：
- $\lambda_r$ 是权重
- $\mathbf{u}_r \in \mathbb{R}^X$，$\mathbf{v}_r \in \mathbb{R}^Y$，$\mathbf{w}_r \in \mathbb{R}^Z$，$\mathbf{s}_r \in \mathbb{R}^C$ 是因子向量
- $\otimes$ 表示外积

元素形式：
$$\mathcal{T}_{ijkc} = \sum_{r=1}^{R} \lambda_r u_{ri} v_{rj} w_{rk} s_{rc}$$

**向量-矩阵（VM）分解**是CP分解的推广，将某些模式组合成矩阵：
$$\mathcal{T} = \sum_{r=1}^{R_1} \mathbf{M}_{r,xy} \otimes \mathbf{v}_{r,z} \otimes \mathbf{b}_{r,1} + \sum_{r=1}^{R_2} \mathbf{M}_{r,xz} \otimes \mathbf{v}_{r,y} \otimes \mathbf{b}_{r,2} + \sum_{r=1}^{R_3} \mathbf{M}_{r,yz} \otimes \mathbf{v}_{r,x} \otimes \mathbf{b}_{r,3}$$

其中 $\mathbf{M}_{r,\cdot} \in \mathbb{R}^{D_1 \times D_2}$ 是矩阵因子。这种分解更灵活，可以捕获平面结构。

### 9.3.2 密度和外观分解

TensoRF 分别建模密度和外观，利用它们的不同特性：

**密度场**使用VM分解：
$$\sigma(\mathbf{x}) = \text{ReLU}\left(\sum_{r=1}^{R_\sigma} \left\langle \mathbf{A}_r^{\sigma}, \mathbf{x} \right\rangle + b_\sigma\right)$$

其中 $\mathbf{A}_r^{\sigma}$ 是通过VM分解得到的：
$$\mathbf{A}_r^{\sigma} = \sum_{i=1}^{3} \mathbf{M}_{r,i}^{\sigma} \otimes \mathbf{v}_{r,i}^{\sigma}$$

这种分解特别适合表示平面结构（墙壁、地板）和轴对齐的几何。

**外观场**结合张量分解和小型神经网络：
$$\mathbf{c}(\mathbf{x}, \mathbf{d}) = \mathcal{F}_\theta\left(\mathbf{f}_{\text{app}}(\mathbf{x}), \mathbf{d}\right)$$

其中外观特征 $\mathbf{f}_{\text{app}}(\mathbf{x})$ 通过VM分解计算：
$$\mathbf{f}_{\text{app}}(\mathbf{x}) = \sum_{r=1}^{R_c} \left\langle \mathbf{A}_r^{c}, \mathbf{x} \right\rangle$$

解码器 $\mathcal{F}_\theta$ 是一个2层MLP，将空间特征和视角方向映射到RGB颜色：
$$\mathcal{F}_\theta: \mathbb{R}^{R_c} \times S^2 \rightarrow [0,1]^3$$

这种混合方法平衡了表达能力和计算效率。

### 9.3.3 存储和计算复杂度

**存储分析**：

对于分辨率 $N^3$ 的密集网格，原始存储需求为 $O(N^3C)$，其中 $C$ 是每个体素的通道数。

TensoRF的存储需求：
- **CP分解**：每个秩-1分量需要 $3N$ 个参数，总共 $3RN$ 参数
- **VM分解**：每个平面 $N^2$ 参数，每个向量 $N$ 参数，总共 $R(N^2 + N) \approx RN^2$ 参数

压缩率：
$$\text{压缩率} = \frac{N^3C}{RN^2} = \frac{NC}{R}$$

典型情况下，$N=300$，$C=4$，$R=16$，压缩率约为 75×。

**计算复杂度**：

查询单个点的值：
- 密集网格：$O(1)$ 数组访问
- CP分解：$O(R)$ 向量内积计算
- VM分解：$O(R)$ 矩阵-向量乘法

虽然单点查询更慢，但：
1. $R \ll N$（通常 $R \sim 16-48$）
2. 向量运算可以高效并行化
3. 缓存局部性更好（因子向量可以驻留在缓存中）

## 9.4 低秩与稀疏表示

### 9.4.1 数学基础

辐射场的低秩结构源于场景的内在规律性——平坦表面、重复纹理、对称性等。这些规律性在适当的基下表现为数据矩阵的低秩性。

**奇异值分解（SVD）**提供了最优低秩近似。对于矩阵 $\mathbf{A} \in \mathbb{R}^{m \times n}$：

$$\mathbf{A} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^T = \sum_{i=1}^{r} \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$

其中：
- $\mathbf{U} \in \mathbb{R}^{m \times r}$：左奇异向量
- $\mathbf{\Sigma} = \text{diag}(\sigma_1, ..., \sigma_r)$：奇异值，$\sigma_1 \geq \sigma_2 \geq ... \geq \sigma_r > 0$
- $\mathbf{V} \in \mathbb{R}^{n \times r}$：右奇异向量

**Eckart-Young定理**：最优 $k$ 秩近似（在任何酉不变范数下）是：
$$\mathbf{A}_k = \sum_{i=1}^{k} \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$

近似误差：
$$\|\mathbf{A} - \mathbf{A}_k\|_F = \sqrt{\sum_{i=k+1}^{r} \sigma_i^2}$$
$$\|\mathbf{A} - \mathbf{A}_k\|_2 = \sigma_{k+1}$$

**有效秩**衡量矩阵的"近似秩"：
$$r_{\text{eff}}(\mathbf{A}) = \frac{\|\mathbf{A}\|_F^2}{\|\mathbf{A}\|_2^2} = \frac{\sum_{i=1}^r \sigma_i^2}{\sigma_1^2}$$

当奇异值快速衰减时，有效秩远小于真实秩，表明数据可压缩。

### 9.4.2 稀疏性诱导

许多场景在适当的变换域中是稀疏的。关键是找到使信号稀疏的基。

**稀疏表示问题**：给定过完备字典 $\mathbf{D} \in \mathbb{R}^{n \times p}$（$p > n$），寻找稀疏系数 $\mathbf{w}$：
$$\mathbf{f} = \mathbf{D}\mathbf{w}, \quad \|\mathbf{w}\|_0 \ll p$$

其中 $\|\mathbf{w}\|_0$ 计算非零元素个数。

由于 $\ell_0$ 优化是NP难的，我们使用 $\ell_1$ 松弛：
$$\min_{\mathbf{w}} \frac{1}{2}\|\mathbf{y} - \mathbf{D}\mathbf{w}\|_2^2 + \lambda\|\mathbf{w}\|_1$$

这是LASSO问题，可通过多种算法求解：

**迭代收缩阈值算法（ISTA）**：
$$\mathbf{w}^{(k+1)} = \mathcal{S}_{\lambda/L}\left(\mathbf{w}^{(k)} + \frac{1}{L}\mathbf{D}^T(\mathbf{y} - \mathbf{D}\mathbf{w}^{(k)})\right)$$

其中 $\mathcal{S}_\tau(x) = \text{sign}(x)\max(|x|-\tau, 0)$ 是软阈值算子，$L$ 是 $\mathbf{D}^T\mathbf{D}$ 的最大特征值。

**快速ISTA（FISTA）**通过动量加速收敛：
$$\mathbf{z}^{(k+1)} = \mathbf{w}^{(k)} + \frac{k-1}{k+2}(\mathbf{w}^{(k)} - \mathbf{w}^{(k-1)})$$

收敛率从 $O(1/k)$ 提升到 $O(1/k^2)$。

### 9.4.3 压缩感知视角

压缩感知理论表明，稀疏信号可以从远少于奈奎斯特速率的测量中精确重建。

**限制等距性质（RIP）**：矩阵 $\mathbf{A}$ 满足阶为 $s$ 的RIP，如果存在 $\delta_s \in (0,1)$ 使得对所有 $s$-稀疏向量 $\mathbf{x}$：
$$(1-\delta_s)\|\mathbf{x}\|_2^2 \leq \|\mathbf{A}\mathbf{x}\|_2^2 \leq (1+\delta_s)\|\mathbf{x}\|_2^2$$

**重建保证**：如果 $\mathbf{A}$ 满足 $\delta_{2s} < \sqrt{2} - 1$，则 $\ell_1$ 最小化的解 $\mathbf{x}^*$ 满足：
$$\|\mathbf{x}^* - \mathbf{x}\|_2 \leq \frac{C}{\sqrt{s}}\|\mathbf{x} - \mathbf{x}_s\|_1 + C'\epsilon$$

其中 $\mathbf{x}_s$ 是 $\mathbf{x}$ 的最佳 $s$ 项近似，$\epsilon$ 是测量噪声。

对于辐射场重建：
- **测量**：多视角图像（每个像素是一个线性测量）
- **信号**：体素化的场景表示
- **稀疏基**：小波、DCT或学习的字典

所需测量数：$M = O(s\log(N/s))$，其中 $s$ 是稀疏度，$N$ 是信号维度。

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
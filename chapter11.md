# 第11章：逆向渲染的数学基础

逆向渲染将计算机图形学的前向问题颠倒过来：给定观测图像，我们寻求恢复场景的底层参数——几何、材质和光照。这一逆问题在计算机视觉、机器学习和物理模拟的交叉点上，形成了可微渲染和神经渲染等现代技术的数学基础。本章建立了解决这些逆问题的严格数学框架，从基本的不适定性分析到实用的优化策略。

逆向渲染的核心挑战在于其本质的不适定性：多种场景配置可能产生相同的图像。例如，明亮的材质在弱光下与暗材质在强光下可能看起来相同。这种歧义性要求我们引入额外的约束和先验知识。同时，渲染过程的非线性和不连续性（如阴影边界）使得优化变得复杂。本章将系统地分析这些挑战，并提供数学工具来应对它们。

从应用角度看，逆向渲染支撑着许多重要技术：从传统的形状重建和材质估计，到现代的神经辐射场优化和可微路径追踪。理解其数学基础不仅有助于改进现有方法，更能指导新算法的设计。

## 学习目标

完成本章后，您将能够：

1. 将逆向渲染形式化为数学逆问题，理解其固有的不适定性
2. 应用正则化技术使不适定问题变得良态
3. 推导并实现计算渲染梯度的高效方法
4. 在贝叶斯框架下构建逆向渲染问题，结合先验知识
5. 分析优化景观并选择合适的求解算法
6. 识别并解决逆向渲染中的常见数值挑战

## 11.1 渲染方程的逆问题

### 前向渲染回顾

标准渲染方程描述了前向问题：

$$L_o(\mathbf{x}, \boldsymbol{\omega}_o) = L_e(\mathbf{x}, \boldsymbol{\omega}_o) + \int_{\Omega} f_r(\mathbf{x}, \boldsymbol{\omega}_i, \boldsymbol{\omega}_o) L_i(\mathbf{x}, \boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n}) \, d\boldsymbol{\omega}_i$$

其中：
- $L_o$：出射辐射度
- $L_e$：自发光
- $f_r$：BRDF
- $L_i$：入射辐射度
- $\Omega$：半球立体角

这个积分方程可以通过Neumann级数展开为路径积分形式：

$$L_o = L_e + \mathcal{T}L_e + \mathcal{T}^2L_e + \cdots = \sum_{k=0}^{\infty} \mathcal{T}^k L_e$$

其中传输算子 $\mathcal{T}$ 定义为：

$$(\mathcal{T}L)(\mathbf{x}, \boldsymbol{\omega}_o) = \int_{\Omega} f_r(\mathbf{x}, \boldsymbol{\omega}_i, \boldsymbol{\omega}_o) L(\mathbf{x}_i, -\boldsymbol{\omega}_i) V(\mathbf{x}, \mathbf{x}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n}) \, d\boldsymbol{\omega}_i$$

传输算子的谱半径 $\rho(\mathcal{T}) < 1$ 保证了Neumann级数的收敛性。这个条件等价于场景中的能量守恒：没有表面能反射超过100%的入射光。

对于体积渲染，我们有：

$$L(\mathbf{r}, \boldsymbol{\omega}) = \int_0^t T(s) \sigma_s(\mathbf{r}(s)) L_s(\mathbf{r}(s), \boldsymbol{\omega}) \, ds + T(t) L_{\text{bg}}$$

其中透射率 $T(s) = \exp\left(-\int_0^s \sigma_t(\mathbf{r}(u)) \, du\right)$。

体积渲染方程同样可以写成算子形式：

$$L = \mathcal{V}[L_s, \sigma_s, \sigma_t]$$

这种算子视角为理解逆问题提供了数学框架。

**路径空间表述**

从路径积分角度，像素值可表示为：

$$I_j = \int_{\mathcal{P}} f_j(\bar{\mathbf{x}}) \, d\mu(\bar{\mathbf{x}})$$

其中 $\mathcal{P}$ 是所有可能光路的空间，$\bar{\mathbf{x}} = (\mathbf{x}_0, \mathbf{x}_1, ..., \mathbf{x}_k)$ 是一条光路，$f_j$ 是路径贡献函数：

$$f_j(\bar{\mathbf{x}}) = W_j(\mathbf{x}_0 \leftarrow \mathbf{x}_1) G(\mathbf{x}_0 \leftrightarrow \mathbf{x}_1) \prod_{i=1}^{k-1} f_r(\mathbf{x}_{i-1} \leftarrow \mathbf{x}_i \leftarrow \mathbf{x}_{i+1}) G(\mathbf{x}_i \leftrightarrow \mathbf{x}_{i+1}) L_e(\mathbf{x}_k \rightarrow \mathbf{x}_{k-1})$$

这里 $W_j$ 是像素 $j$ 的重要性函数，$G$ 是几何项。路径空间视角揭示了逆向渲染的核心挑战：从有限的积分值（像素）推断被积函数（场景参数）。

**测量方程的完整形式**

实际的图像形成过程包含更多细节：

$$I_j = \int_{\lambda} \int_t \int_A \int_{\Omega} R_j(\lambda) S(t) W(\mathbf{x}, \boldsymbol{\omega}) L(\mathbf{x}, \boldsymbol{\omega}, t, \lambda) \cos\theta \, d\boldsymbol{\omega} \, dA \, dt \, d\lambda + n_j$$

其中：
- $R_j(\lambda)$：传感器光谱响应
- $S(t)$：快门函数
- $W(\mathbf{x}, \boldsymbol{\omega})$：像素滤波器
- $n_j$：传感器噪声

这个完整模型对于处理实际捕获的图像至关重要。

### 逆问题形式化

设 $\boldsymbol{\theta}$ 表示所有场景参数（几何、材质、光照），$\mathbf{I}_{\text{obs}}$ 为观测图像。逆向渲染寻求：

$$\boldsymbol{\theta}^* = \arg\min_{\boldsymbol{\theta}} \mathcal{L}(\mathcal{R}(\boldsymbol{\theta}), \mathbf{I}_{\text{obs}}) + \mathcal{R}_{\text{reg}}(\boldsymbol{\theta})$$

其中：
- $\mathcal{R}$：渲染算子
- $\mathcal{L}$：损失函数（如 $L^2$ 或感知损失）
- $\mathcal{R}_{\text{reg}}$：正则化项

从函数分析角度，渲染算子 $\mathcal{R}: \Theta \rightarrow \mathcal{I}$ 映射参数空间到图像空间。逆问题寻求逆映射 $\mathcal{R}^{-1}$，但这个逆映射通常不存在或不唯一。

**算子的数学性质**

渲染算子 $\mathcal{R}$ 具有以下关键性质：

1. **非线性性**：$\mathcal{R}(\alpha\boldsymbol{\theta}_1 + \beta\boldsymbol{\theta}_2) \neq \alpha\mathcal{R}(\boldsymbol{\theta}_1) + \beta\mathcal{R}(\boldsymbol{\theta}_2)$
   
2. **非单射性**：存在 $\boldsymbol{\theta}_1 \neq \boldsymbol{\theta}_2$ 使得 $\mathcal{R}(\boldsymbol{\theta}_1) = \mathcal{R}(\boldsymbol{\theta}_2)$

3. **紧算子性质**：在适当的函数空间中，$\mathcal{R}$ 常表现为紧算子，将有界集映射到预紧集

4. **Fréchet可微性**：在大多数实际情况下，$\mathcal{R}$ 是Fréchet可微的，即存在线性算子 $D\mathcal{R}[\boldsymbol{\theta}]$ 使得：
   $$\lim_{\|\mathbf{h}\| \to 0} \frac{\|\mathcal{R}(\boldsymbol{\theta} + \mathbf{h}) - \mathcal{R}(\boldsymbol{\theta}) - D\mathcal{R}[\boldsymbol{\theta}]\mathbf{h}\|}{\|\mathbf{h}\|} = 0$$

**损失函数的选择**

不同的损失函数对应不同的统计假设和感知特性：

1. **$L^2$ 损失**：
   $$\mathcal{L}_{L^2} = \sum_j |\mathcal{R}(\boldsymbol{\theta})_j - I_{\text{obs},j}|^2$$
   对应高斯噪声假设，但可能过度惩罚异常值。

2. **$L^1$ 损失**：
   $$\mathcal{L}_{L^1} = \sum_j |\mathcal{R}(\boldsymbol{\theta})_j - I_{\text{obs},j}|$$
   对异常值更鲁棒，对应拉普拉斯噪声。

3. **感知损失**：
   $$\mathcal{L}_{\text{perceptual}} = \sum_l \lambda_l \|\phi_l(\mathcal{R}(\boldsymbol{\theta})) - \phi_l(\mathbf{I}_{\text{obs}})\|^2$$
   其中 $\phi_l$ 是预训练网络的第 $l$ 层特征。

4. **对抗损失**：
   $$\mathcal{L}_{\text{adv}} = -\log D(\mathcal{R}(\boldsymbol{\theta}))$$
   其中 $D$ 是判别器网络。

5. **结构相似度(SSIM)**：
   $$\mathcal{L}_{\text{SSIM}} = 1 - \text{SSIM}(\mathcal{R}(\boldsymbol{\theta}), \mathbf{I}_{\text{obs}})$$

### 参数空间分解

我们通常将 $\boldsymbol{\theta}$ 分解为：

$$\boldsymbol{\theta} = (\boldsymbol{\theta}_g, \boldsymbol{\theta}_m, \boldsymbol{\theta}_l)$$

- **几何参数** $\boldsymbol{\theta}_g$：顶点位置、法线、SDF值
- **材质参数** $\boldsymbol{\theta}_m$：反照率、粗糙度、折射率
- **光照参数** $\boldsymbol{\theta}_l$：环境贴图、点光源位置/强度

每个参数空间有其特定的流形结构：
- 几何空间：可能受限于可制造性或物理约束
- 材质空间：受物理定律约束（如能量守恒）
- 光照空间：通常为正值函数空间

### 观测模型

实际观测包含噪声：

$$\mathbf{I}_{\text{obs}} = \mathcal{R}(\boldsymbol{\theta}_{\text{true}}) + \boldsymbol{\epsilon}$$

其中 $\boldsymbol{\epsilon}$ 表示传感器噪声、量化误差等。

噪声模型的选择影响优化目标：
- 高斯噪声 $\boldsymbol{\epsilon} \sim \mathcal{N}(0, \sigma^2\mathbf{I})$ 导致 $L^2$ 损失
- 泊松噪声（光子计数）需要不同的似然函数
- 混合噪声模型更贴近实际传感器

### 不唯一性示例

逆向渲染的根本挑战在于解的不唯一性：

**示例1：反照率-光照歧义**
给定观测亮度 $B$，无法区分：
- 暗材质 + 强光照：$\rho_1 L_1 = B$，其中 $\rho_1 = 0.2, L_1 = 5$
- 亮材质 + 弱光照：$\rho_2 L_2 = B$，其中 $\rho_2 = 1.0, L_2 = 1$

数学上，这对应于零空间的存在：
$$\mathcal{N}(\mathcal{R}) = \{\Delta\boldsymbol{\theta} : \mathcal{R}(\boldsymbol{\theta} + \Delta\boldsymbol{\theta}) = \mathcal{R}(\boldsymbol{\theta})\}$$

零空间的维度可通过分析渲染算子的雅可比矩阵秩来估计。对于典型场景：
$$\text{dim}(\mathcal{N}(\mathcal{R})) \geq \text{dim}(\Theta) - \text{dim}(\mathcal{I})$$

**示例2：凹凸贴图vs几何细节**
表面法线扰动可由以下产生：
- 真实几何位移：$\mathbf{p}' = \mathbf{p} + h(\mathbf{p})\mathbf{n}$
- 法线贴图：$\mathbf{n}' = \text{normalize}(\mathbf{n} + \nabla h)$
- 位移贴图

这些在单视角下可能产生相同的着色。数学上，对于小位移 $|h| \ll 1$：
$$\mathbf{n}' \approx \mathbf{n} + (\mathbf{I} - \mathbf{n}\mathbf{n}^T)\nabla h$$

两种方法在一阶近似下等价，导致局部不可区分性。

**示例3：形状-BRDF歧义**
镜面球和漫反射椭球在特定视角和光照下可能产生相同图像。设球面参数化为 $\mathbf{x}(\theta, \phi)$，BRDF为 $f_r$，则存在椭球参数化 $\mathbf{y}(\theta, \phi)$ 和 BRDF $g_r$ 使得：

$$\int_{\Omega} f_r(\mathbf{x}, \boldsymbol{\omega}_i, \boldsymbol{\omega}_o) L_i \cos\theta_i \, d\boldsymbol{\omega}_i = \int_{\Omega} g_r(\mathbf{y}, \boldsymbol{\omega}_i, \boldsymbol{\omega}_o) L_i \cos\theta_i' \, d\boldsymbol{\omega}_i$$

这种歧义可通过Helmholtz互易性部分解决：观察位置和光源位置交换应产生相同结果。

**示例4：全局光照中的路径等价**
多次反射可通过不同路径达到相同结果：
- 直接光照 + 强环境光
- 弱直接光 + 多次反射

这在渲染方程的Neumann级数展开中表现为不同截断阶数的等价性：
$$L = \sum_{k=0}^{n} \mathcal{T}^k L_e + \mathcal{T}^{n+1} L_{\text{env}}$$

当 $\|\mathcal{T}^{n+1}\|$ 很小时，高阶项可被环境光近似。

**示例5：频率混叠**
高频纹理细节在有限分辨率下产生混叠：
$$I_{\text{observed}} = \int_{\Omega} \text{sinc}(f_s(\mathbf{x} - \mathbf{x}_0)) T(\mathbf{x}) \, d\mathbf{x}$$

其中 $f_s$ 是采样频率，$T$ 是纹理函数。当 $T$ 包含频率高于 $f_s/2$ 的成分时，不同的高频纹理可能产生相同的采样结果。

**示例6：体积密度-发射歧义**
在体积渲染中，沿射线的不同密度和发射分布可产生相同积分：
$$\int_0^L T_1(t)\sigma_{s,1}(t)c_1(t) \, dt = \int_0^L T_2(t)\sigma_{s,2}(t)c_2(t) \, dt$$

即使 $(\sigma_{s,1}, c_1) \neq (\sigma_{s,2}, c_2)$。这种歧义在单视角下无法解决。

## 11.2 不适定性与正则化

### Hadamard不适定性

问题称为良态（well-posed）需满足：
1. 解存在（存在性）
2. 解唯一（唯一性）
3. 解连续依赖于数据（稳定性）

逆向渲染通常违反条件2和3，使其成为不适定问题。

形式化地，对于算子方程 $\mathcal{R}(\boldsymbol{\theta}) = \mathbf{I}$：
- **存在性**：$\mathbf{I} \in \text{Range}(\mathcal{R})$ 可能不满足（如物理不可实现的图像）
- **唯一性**：$\text{dim}(\mathcal{N}(\mathcal{R})) > 0$ （零空间非平凡）
- **稳定性**：小扰动 $\|\Delta\mathbf{I}\| = \epsilon$ 可能导致大变化 $\|\Delta\boldsymbol{\theta}\| \gg \epsilon$

### 条件数分析

考虑线性化渲染算子 $\mathbf{J} = \partial \mathcal{R} / \partial \boldsymbol{\theta}$。条件数：

$$\kappa(\mathbf{J}) = \frac{\sigma_{\max}(\mathbf{J})}{\sigma_{\min}(\mathbf{J})}$$

大条件数表示数值不稳定性。对于典型渲染问题：
- 几何参数：$\kappa \sim 10^2 - 10^3$（中等病态）
- 材质参数：$\kappa \sim 10^3 - 10^5$（严重病态）
- 光照参数：$\kappa \sim 10^4 - 10^6$（极度病态）

**奇异值分解视角**：
设 $\mathbf{J} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T$，则：
$$\Delta\boldsymbol{\theta} = \sum_{i=1}^r \frac{\mathbf{u}_i^T \Delta\mathbf{I}}{\sigma_i} \mathbf{v}_i$$

小奇异值 $\sigma_i$ 放大噪声，导致解的不稳定。

### Tikhonov正则化

最基本的正则化方法添加 $L^2$ 惩罚：

$$\boldsymbol{\theta}^* = \arg\min_{\boldsymbol{\theta}} \|\mathcal{R}(\boldsymbol{\theta}) - \mathbf{I}_{\text{obs}}\|^2 + \lambda \|\boldsymbol{\theta} - \boldsymbol{\theta}_0\|^2$$

其中 $\boldsymbol{\theta}_0$ 是先验估计，$\lambda$ 控制正则化强度。

**正则化的几何解释**：
在参数空间中，正则化项定义了一个"信任区域"：
$$\{\boldsymbol{\theta} : \|\boldsymbol{\theta} - \boldsymbol{\theta}_0\|^2 \leq \delta\}$$

解是数据拟合和正则化约束的平衡点。

**广义Tikhonov正则化**：
$$\boldsymbol{\theta}^* = \arg\min_{\boldsymbol{\theta}} \|\mathcal{R}(\boldsymbol{\theta}) - \mathbf{I}_{\text{obs}}\|^2 + \lambda \|\mathbf{L}(\boldsymbol{\theta} - \boldsymbol{\theta}_0)\|^2$$

其中 $\mathbf{L}$ 是正则化算子，如：
- $\mathbf{L} = \mathbf{I}$：标准Tikhonov
- $\mathbf{L} = \nabla$：梯度正则化（平滑性）
- $\mathbf{L} = \nabla^2$：曲率正则化

**正则化参数选择**

选择最优 $\lambda$ 是关键挑战。常用方法包括：

1. **L曲线方法**：绘制 $\log\|\mathcal{R}(\boldsymbol{\theta}_\lambda) - \mathbf{I}_{\text{obs}}\|$ vs $\log\|\boldsymbol{\theta}_\lambda - \boldsymbol{\theta}_0\|$，选择曲线拐点

2. **广义交叉验证(GCV)**：
   $$\lambda_{\text{GCV}} = \arg\min_\lambda \frac{\|\mathcal{R}(\boldsymbol{\theta}_\lambda) - \mathbf{I}_{\text{obs}}\|^2}{[1 - \text{tr}(\mathbf{A}_\lambda)/n]^2}$$
   其中 $\mathbf{A}_\lambda$ 是影响矩阵

3. **Morozov不符原则**：选择 $\lambda$ 使得残差等于噪声水平：
   $$\|\mathcal{R}(\boldsymbol{\theta}_\lambda) - \mathbf{I}_{\text{obs}}\| = \delta$$
   其中 $\delta$ 是噪声估计

**分数阶Tikhonov正则化**

对于具有幂律衰减特性的问题：
$$\boldsymbol{\theta}^* = \arg\min_{\boldsymbol{\theta}} \|\mathcal{R}(\boldsymbol{\theta}) - \mathbf{I}_{\text{obs}}\|^2 + \lambda \|(-\Delta)^{\alpha/2}(\boldsymbol{\theta} - \boldsymbol{\theta}_0)\|^2$$

其中 $(-\Delta)^{\alpha/2}$ 是分数阶拉普拉斯算子，$0 < \alpha < 2$。在傅里叶域：
$$\widehat{(-\Delta)^{\alpha/2}f}(\boldsymbol{k}) = |\boldsymbol{k}|^\alpha \hat{f}(\boldsymbol{k})$$

### 稀疏性先验

**L1正则化**促进稀疏解：

$$\mathcal{R}_{\text{L1}}(\boldsymbol{\theta}) = \lambda_1 \|\boldsymbol{\theta}\|_1$$

稀疏性在渲染中的应用：
- 环境光的球谐系数
- 纹理的小波系数
- 几何的特征点

**总变分（TV）正则化**保持边缘：

$$\mathcal{R}_{\text{TV}}(\mathbf{u}) = \lambda_{\text{TV}} \int_{\Omega} |\nabla \mathbf{u}| \, d\mathbf{x}$$

对于离散图像：

$$\mathcal{R}_{\text{TV}}(\mathbf{u}) = \lambda_{\text{TV}} \sum_{i,j} \sqrt{(u_{i+1,j} - u_{i,j})^2 + (u_{i,j+1} - u_{i,j})^2 + \epsilon}$$

**各向异性TV**：
$$\mathcal{R}_{\text{ATV}}(\mathbf{u}) = \lambda_{\text{TV}} \sum_{i,j} |u_{i+1,j} - u_{i,j}| + |u_{i,j+1} - u_{i,j}|$$

计算更简单但可能产生阶梯效应。

### 物理约束

**能量守恒**：
$$\int_{\Omega} f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_i \cdot \mathbf{n}) \, d\boldsymbol{\omega}_i \leq 1$$

实际实现中，通过反照率上界约束：
$$\rho(\lambda) \leq 1, \quad \forall \lambda \in [380\text{nm}, 780\text{nm}]$$

对于参数化BRDF模型，能量守恒可表示为参数空间的约束：
$$\mathcal{C}_{\text{energy}} = \{\boldsymbol{\alpha} : \rho_{\text{hd}}(\boldsymbol{\omega}_o; \boldsymbol{\alpha}) \leq 1, \forall \boldsymbol{\omega}_o\}$$

其中半球方向反射率：
$$\rho_{\text{hd}}(\boldsymbol{\omega}_o; \boldsymbol{\alpha}) = \int_{\Omega^+} f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o; \boldsymbol{\alpha}) (\boldsymbol{\omega}_i \cdot \mathbf{n}) \, d\boldsymbol{\omega}_i$$

**互易性**：
$$f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = f_r(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)$$

对于参数化BRDF，这要求参数满足特定对称性。例如，对于微表面模型：
$$D(\mathbf{h}) G(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o, \mathbf{h}) = D(\mathbf{h}) G(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i, \mathbf{h})$$

在优化中，可通过对称化强制互易性：
$$f_r^{\text{sym}}(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \frac{1}{2}[f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) + f_r(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)]$$

**正定性**：
$$f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) \geq 0$$

这些可作为硬约束或软惩罚项：
- 硬约束：投影到可行集 $\mathcal{C} = \{\boldsymbol{\theta} : g(\boldsymbol{\theta}) \leq 0\}$
- 软约束：惩罚项 $\lambda \sum_i \max(0, g_i(\boldsymbol{\theta}))^2$
- 障碍方法：$-\mu \sum_i \log(-g_i(\boldsymbol{\theta}))$，其中 $\mu \to 0$

**几何约束**：
- 表面法线单位化：$\|\mathbf{n}\|^2 = 1$
- 凸性约束：$\nabla^2 f(\mathbf{x}) \succeq 0$（Hessian半正定）
- 体积保持：$\int_V dV = V_0$
- 表面积正则化：$\mathcal{R}_{\text{area}} = \int_S \sqrt{1 + |\nabla f|^2} \, dS$

**材质物理约束**：
- Fresnel约束：反射率随入射角变化遵循Fresnel方程
- 色散关系：折射率的实部和虚部通过Kramers-Kronig关系相联
- 微表面法线分布归一化：$\int_{\Omega^+} D(\mathbf{h}) (\mathbf{h} \cdot \mathbf{n}) \, d\mathbf{h} = 1$

**光照约束**：
- 非负性：$L(\mathbf{x}, \boldsymbol{\omega}) \geq 0$
- 远场衰减：$L(\mathbf{x}, \boldsymbol{\omega}) = O(1/\|\mathbf{x}\|^2)$ 当 $\|\mathbf{x}\| \to \infty$
- 球谐系数衰减：$|L_{lm}| \leq C \cdot l^{-\alpha}$，典型 $\alpha > 2$

## 11.3 梯度计算：伴随法与自动微分

高效的梯度计算是逆向渲染的核心。对于复杂的渲染系统，手动推导梯度既繁琐又容易出错。本节介绍三种主要方法：解析推导、伴随状态方法和自动微分。

### 解析梯度推导

对于简单情况，可直接推导梯度。考虑Lambert着色：

$$I = \mathbf{n} \cdot \mathbf{l} \cdot \rho$$

对反照率的梯度：
$$\frac{\partial I}{\partial \rho} = \mathbf{n} \cdot \mathbf{l}$$

对法线的梯度（需要投影到切平面）：
$$\frac{\partial I}{\partial \mathbf{n}} = (\mathbf{I} - \mathbf{n}\mathbf{n}^T)(\mathbf{l} \cdot \rho)$$

投影矩阵 $\mathbf{P} = \mathbf{I} - \mathbf{n}\mathbf{n}^T$ 确保梯度在切平面内，保持 $\|\mathbf{n}\| = 1$。

**Phong模型的梯度**：
$$I = k_d(\mathbf{n} \cdot \mathbf{l}) + k_s(\mathbf{r} \cdot \mathbf{v})^n$$

其中反射向量 $\mathbf{r} = 2(\mathbf{n} \cdot \mathbf{l})\mathbf{n} - \mathbf{l}$。

对镜面指数的梯度：
$$\frac{\partial I}{\partial n} = k_s (\mathbf{r} \cdot \mathbf{v})^n \ln(\mathbf{r} \cdot \mathbf{v})$$

注意当 $\mathbf{r} \cdot \mathbf{v} \approx 0$ 时的数值稳定性问题。

### 伴随状态方法

对于复杂系统，伴随法高效计算梯度。给定约束优化：

$$\min_{\boldsymbol{\theta}} J(\mathbf{u}, \boldsymbol{\theta}) \quad \text{s.t.} \quad F(\mathbf{u}, \boldsymbol{\theta}) = 0$$

其中 $\mathbf{u}$ 是状态变量（如辐射度场），$F$ 是约束（如渲染方程）。

拉格朗日函数：
$$\mathcal{L} = J(\mathbf{u}, \boldsymbol{\theta}) + \boldsymbol{\lambda}^T F(\mathbf{u}, \boldsymbol{\theta})$$

伴随方程：
$$\left(\frac{\partial F}{\partial \mathbf{u}}\right)^T \boldsymbol{\lambda} = -\left(\frac{\partial J}{\partial \mathbf{u}}\right)^T$$

参数梯度：
$$\frac{dJ}{d\boldsymbol{\theta}} = \frac{\partial J}{\partial \boldsymbol{\theta}} + \boldsymbol{\lambda}^T \frac{\partial F}{\partial \boldsymbol{\theta}}$$

**渲染方程的伴随**：
对于渲染方程 $L = L_e + \mathcal{T}L$，伴随辐射度 $L^*$ 满足：
$$L^* = \frac{\partial J}{\partial L} + \mathcal{T}^* L^*$$

其中 $\mathcal{T}^*$ 是伴随传输算子。

**计算优势**：
- 前向求解：$O(N_{\text{params}} \times N_{\text{pixels}})$
- 伴随方法：$O(N_{\text{pixels}})$，与参数数量无关

对于高维参数空间（如纹理、环境光），伴随法显著提高效率。

### 自动微分原理

**前向模式**计算雅可比-向量积：
$$\dot{y} = \frac{\partial f}{\partial \mathbf{x}} \dot{\mathbf{x}}$$

复杂度：$O(n)$ 对于 $n$ 个输入。

**反向模式**计算向量-雅可比积：
$$\bar{\mathbf{x}} = \left(\frac{\partial f}{\partial \mathbf{x}}\right)^T \bar{y}$$

复杂度：$O(m)$ 对于 $m$ 个输出。

对于损失函数（$m=1$），反向模式更高效。

### 计算图优化

**检查点技术**权衡内存与计算：
- 前向传播时仅存储关键中间结果
- 反向传播时重新计算需要的值

**梯度累积**处理大批量：
$$\nabla_{\boldsymbol{\theta}} \mathcal{L} = \frac{1}{N} \sum_{i=1}^{N} \nabla_{\boldsymbol{\theta}} \mathcal{L}_i$$

可分批计算并累积。

## 11.4 贝叶斯推断框架

### 贝叶斯公式

后验分布：
$$p(\boldsymbol{\theta} | \mathbf{I}_{\text{obs}}) = \frac{p(\mathbf{I}_{\text{obs}} | \boldsymbol{\theta}) p(\boldsymbol{\theta})}{p(\mathbf{I}_{\text{obs}})}$$

- $p(\mathbf{I}_{\text{obs}} | \boldsymbol{\theta})$：似然函数
- $p(\boldsymbol{\theta})$：先验分布
- $p(\mathbf{I}_{\text{obs}})$：证据（归一化常数）

### 似然模型

假设高斯噪声：
$$p(\mathbf{I}_{\text{obs}} | \boldsymbol{\theta}) = \mathcal{N}(\mathcal{R}(\boldsymbol{\theta}), \sigma^2 \mathbf{I})$$

对数似然：
$$\log p(\mathbf{I}_{\text{obs}} | \boldsymbol{\theta}) = -\frac{1}{2\sigma^2} \|\mathcal{R}(\boldsymbol{\theta}) - \mathbf{I}_{\text{obs}}\|^2 + \text{const}$$

### 先验设计

**材质先验**：
- 反照率：Beta分布约束到[0,1]
- 粗糙度：截断高斯分布
- 法线：von Mises-Fisher分布在单位球上

**几何先验**：
- 平滑性：高斯过程
- 稀疏性：Laplace分布

**光照先验**：
- 自然光照：球谐系数的能量衰减

### 最大后验(MAP)估计

$$\boldsymbol{\theta}_{\text{MAP}} = \arg\max_{\boldsymbol{\theta}} p(\boldsymbol{\theta} | \mathbf{I}_{\text{obs}}) = \arg\max_{\boldsymbol{\theta}} \log p(\mathbf{I}_{\text{obs}} | \boldsymbol{\theta}) + \log p(\boldsymbol{\theta})$$

这等价于带正则化的优化问题。

### 不确定性量化

**拉普拉斯近似**：
$$p(\boldsymbol{\theta} | \mathbf{I}_{\text{obs}}) \approx \mathcal{N}(\boldsymbol{\theta}_{\text{MAP}}, \mathbf{H}^{-1})$$

其中Hessian矩阵：
$$\mathbf{H} = -\nabla^2 \log p(\boldsymbol{\theta} | \mathbf{I}_{\text{obs}})|_{\boldsymbol{\theta}_{\text{MAP}}}$$

**马尔可夫链蒙特卡洛(MCMC)**：
- Metropolis-Hastings算法
- 哈密顿蒙特卡洛(HMC)利用梯度信息

## 11.5 凸优化与非凸优化

### 凸性分析

函数 $f$ 为凸当且仅当：
$$f(\lambda \mathbf{x} + (1-\lambda)\mathbf{y}) \leq \lambda f(\mathbf{x}) + (1-\lambda)f(\mathbf{y})$$

对于二阶可微函数，等价于 $\nabla^2 f \succeq 0$（半正定）。

### 渲染中的非凸性来源

1. **可见性函数**：阶跃不连续
2. **BRDF**：镜面波瓣产生多峰
3. **全局光照**：多次反射的递归性
4. **参数化**：欧拉角的周期性

### 凸松弛技术

**凸包络**：最大的凸下界
$$f_{\text{env}}(\mathbf{x}) = \sup\{g(\mathbf{x}) : g \text{ 凸且 } g \leq f\}$$

**半定规划(SDP)松弛**：
将 $\mathbf{x}\mathbf{x}^T$ 松弛为半正定矩阵 $\mathbf{X} \succeq 0$。

### 优化算法

**一阶方法**：
- 梯度下降：$\boldsymbol{\theta}_{k+1} = \boldsymbol{\theta}_k - \alpha \nabla f(\boldsymbol{\theta}_k)$
- 动量法：$\mathbf{v}_{k+1} = \beta \mathbf{v}_k + \nabla f(\boldsymbol{\theta}_k)$
- Adam：自适应学习率

**二阶方法**：
- 牛顿法：$\boldsymbol{\theta}_{k+1} = \boldsymbol{\theta}_k - [\nabla^2 f(\boldsymbol{\theta}_k)]^{-1} \nabla f(\boldsymbol{\theta}_k)$
- L-BFGS：有限内存拟牛顿

**交替方向乘子法(ADMM)**：
将问题分解为子问题：
$$\min_{\mathbf{x}, \mathbf{z}} f(\mathbf{x}) + g(\mathbf{z}) \quad \text{s.t.} \quad \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{z} = \mathbf{c}$$

迭代更新：
1. $\mathbf{x}_{k+1} = \arg\min_{\mathbf{x}} f(\mathbf{x}) + \frac{\rho}{2}\|\mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{z}_k - \mathbf{c} + \mathbf{u}_k\|^2$
2. $\mathbf{z}_{k+1} = \arg\min_{\mathbf{z}} g(\mathbf{z}) + \frac{\rho}{2}\|\mathbf{A}\mathbf{x}_{k+1} + \mathbf{B}\mathbf{z} - \mathbf{c} + \mathbf{u}_k\|^2$
3. $\mathbf{u}_{k+1} = \mathbf{u}_k + \mathbf{A}\mathbf{x}_{k+1} + \mathbf{B}\mathbf{z}_{k+1} - \mathbf{c}$

### 收敛性保证

**Lipschitz连续性**：
$$\|\nabla f(\mathbf{x}) - \nabla f(\mathbf{y})\| \leq L \|\mathbf{x} - \mathbf{y}\|$$

梯度下降收敛率（步长 $\alpha = 1/L$）：
$$f(\mathbf{x}_k) - f(\mathbf{x}^*) \leq \frac{2L\|\mathbf{x}_0 - \mathbf{x}^*\|^2}{k}$$

**强凸性**：
$$f(\mathbf{y}) \geq f(\mathbf{x}) + \nabla f(\mathbf{x})^T(\mathbf{y} - \mathbf{x}) + \frac{\mu}{2}\|\mathbf{y} - \mathbf{x}\|^2$$

线性收敛率：
$$\|\mathbf{x}_k - \mathbf{x}^*\| \leq \left(1 - \frac{\mu}{L}\right)^k \|\mathbf{x}_0 - \mathbf{x}^*\|$$

## 本章小结

本章建立了逆向渲染的数学基础，从基本的逆问题形式化到实用的优化技术。关键要点包括：

1. **逆问题本质**：逆向渲染是不适定的，存在解的不唯一性和对噪声的敏感性
2. **正则化必要性**：通过Tikhonov、稀疏性和物理约束使问题良态化
3. **梯度计算**：伴随法和自动微分提供了高效的梯度计算方法
4. **贝叶斯框架**：将先验知识系统地融入优化过程
5. **优化策略**：理解凸性有助于算法选择和收敛性分析

核心数学工具：
- 变分法和伴随状态方法
- 正则化理论
- 自动微分
- 贝叶斯推断
- 凸优化理论

## 练习题

### 基础题

**练习11.1**：证明渲染方程关于BRDF参数的雅可比矩阵
考虑简化的直接光照：
$$L = \int_{\Omega} f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o; \boldsymbol{\alpha}) L_i(\boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n}) \, d\boldsymbol{\omega}_i$$

其中 $\boldsymbol{\alpha}$ 是BRDF参数。推导 $\partial L / \partial \boldsymbol{\alpha}$。

<details>
<summary>提示</summary>
利用Leibniz积分法则，将导数移入积分内。注意BRDF的参数化形式。
</details>

<details>
<summary>答案</summary>

应用Leibniz法则：
$$\frac{\partial L}{\partial \boldsymbol{\alpha}} = \int_{\Omega} \frac{\partial f_r}{\partial \boldsymbol{\alpha}}(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o; \boldsymbol{\alpha}) L_i(\boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n}) \, d\boldsymbol{\omega}_i$$

对于具体BRDF模型，如Phong模型 $f_r = k_d/\pi + k_s \frac{n+2}{2\pi}(\mathbf{r} \cdot \boldsymbol{\omega}_o)^n$：
- $\partial f_r / \partial k_d = 1/\pi$
- $\partial f_r / \partial k_s = \frac{n+2}{2\pi}(\mathbf{r} \cdot \boldsymbol{\omega}_o)^n$
- $\partial f_r / \partial n = k_s \frac{1}{2\pi}[(\mathbf{r} \cdot \boldsymbol{\omega}_o)^n + (n+2)(\mathbf{r} \cdot \boldsymbol{\omega}_o)^n \ln(\mathbf{r} \cdot \boldsymbol{\omega}_o)]$
</details>

**练习11.2**：分析Tikhonov正则化的偏差-方差权衡
给定线性逆问题 $\mathbf{A}\mathbf{x} = \mathbf{b} + \boldsymbol{\epsilon}$，其中 $\boldsymbol{\epsilon} \sim \mathcal{N}(0, \sigma^2\mathbf{I})$。
Tikhonov解为：
$$\mathbf{x}_{\lambda} = (\mathbf{A}^T\mathbf{A} + \lambda\mathbf{I})^{-1}\mathbf{A}^T\mathbf{b}$$

推导偏差和方差关于 $\lambda$ 的表达式。

<details>
<summary>提示</summary>
使用SVD分解 $\mathbf{A} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T$，并分析每个奇异值的贡献。
</details>

<details>
<summary>答案</summary>

设 $\mathbf{A} = \sum_i \sigma_i \mathbf{u}_i \mathbf{v}_i^T$，真实解 $\mathbf{x}_{\text{true}}$。

偏差：
$$\text{Bias}[\mathbf{x}_{\lambda}] = \mathbb{E}[\mathbf{x}_{\lambda}] - \mathbf{x}_{\text{true}} = -\lambda(\mathbf{A}^T\mathbf{A} + \lambda\mathbf{I})^{-1}\mathbf{x}_{\text{true}}$$

在SVD基下：
$$\|\text{Bias}[\mathbf{x}_{\lambda}]\|^2 = \sum_i \frac{\lambda^2}{(\sigma_i^2 + \lambda)^2} |\mathbf{v}_i^T \mathbf{x}_{\text{true}}|^2$$

方差：
$$\text{Var}[\mathbf{x}_{\lambda}] = \sigma^2 \text{tr}[(\mathbf{A}^T\mathbf{A} + \lambda\mathbf{I})^{-1}\mathbf{A}^T\mathbf{A}(\mathbf{A}^T\mathbf{A} + \lambda\mathbf{I})^{-1}]$$

$$= \sigma^2 \sum_i \frac{\sigma_i^2}{(\sigma_i^2 + \lambda)^2}$$

总误差（MSE）：
$$\text{MSE} = \sum_i \left[\frac{\lambda^2}{(\sigma_i^2 + \lambda)^2} |\mathbf{v}_i^T \mathbf{x}_{\text{true}}|^2 + \frac{\sigma^2\sigma_i^2}{(\sigma_i^2 + \lambda)^2}\right]$$

最优 $\lambda$ 平衡两项。
</details>

**练习11.3**：实现伴随法计算体积渲染梯度
对于简化的体积渲染：
$$I = \int_0^L T(t) \sigma(t) c(t) \, dt$$
其中 $T(t) = \exp(-\int_0^t \sigma(s) \, ds)$。

推导 $\partial I / \partial \sigma(t)$ 使用伴随方法。

<details>
<summary>提示</summary>
定义伴随变量 $\lambda(t)$ 满足适当的终端条件，利用分部积分。
</details>

<details>
<summary>答案</summary>

定义哈密顿量：
$$H = T(t)\sigma(t)c(t) + \lambda(t)[-T(t)\sigma(t)]$$

伴随方程（从 $t=L$ 反向）：
$$\frac{d\lambda}{dt} = \frac{\partial H}{\partial T} = \sigma(t)[c(t) - \lambda(t)]$$

边界条件：$\lambda(L) = 0$

梯度：
$$\frac{\partial I}{\partial \sigma(t)} = T(t)[c(t) - \int_t^L \sigma(s)c(s)T(s)/T(t) \, ds]$$

这给出了每个点对最终图像的贡献，考虑了遮挡效应。
</details>

### 挑战题

**练习11.4**：证明渲染逆问题的不适定性
考虑单一像素观测下的球面反照率恢复。设球面由Lambert材质组成，在均匀环境光下。证明存在无穷多反照率分布产生相同观测。

<details>
<summary>提示</summary>
利用球谐函数展开，分析零空间维度。
</details>

<details>
<summary>答案</summary>

球面上的反照率函数：$\rho(\theta, \phi)$
观测积分：
$$I = \int_{S^2} \rho(\theta, \phi) \max(0, \cos\theta) \, d\Omega$$

将 $\rho$ 展开为球谐函数：
$$\rho(\theta, \phi) = \sum_{l,m} a_{lm} Y_l^m(\theta, \phi)$$

观测变为：
$$I = \sum_{l,m} a_{lm} \int_{S^2} Y_l^m(\theta, \phi) \max(0, \cos\theta) \, d\Omega$$

由于 $\max(0, \cos\theta)$ 只有有限的球谐系数非零（主要是 $l=0,1$），高阶项 $(l \geq 2)$ 的贡献可能为零。

零空间包含所有满足：
$$\int_{S^2} Y_l^m(\theta, \phi) \max(0, \cos\theta) \, d\Omega = 0$$
的球谐基函数。

因此，任何形如 $\rho + \sum_{(l,m) \in \text{null}} b_{lm} Y_l^m$ 的反照率都产生相同观测，证明了不唯一性。
</details>

**练习11.5**：设计自适应正则化策略
给定多尺度问题，不同参数需要不同正则化强度。设计一个基于局部曲率的自适应正则化方案。

<details>
<summary>提示</summary>
考虑Hessian矩阵的谱分解，根据条件数调整正则化。
</details>

<details>
<summary>答案</summary>

局部Hessian近似：
$$\mathbf{H}_i = \nabla^2 f(\boldsymbol{\theta})_i$$

谱分解：
$$\mathbf{H}_i = \mathbf{Q}_i \boldsymbol{\Lambda}_i \mathbf{Q}_i^T$$

自适应正则化矩阵：
$$\mathbf{R}_i = \mathbf{Q}_i \text{diag}(\lambda_1/\sigma_1, \ldots, \lambda_n/\sigma_n) \mathbf{Q}_i^T$$

其中 $\lambda_j = \lambda_0 \cdot \max(1, \sigma_j/\sigma_{\text{threshold}})$

更新规则：
$$\boldsymbol{\theta}_{k+1} = \boldsymbol{\theta}_k - (\mathbf{H} + \mathbf{R})^{-1} \nabla f(\boldsymbol{\theta}_k)$$

这在病态方向（小特征值）施加更强正则化，在良态方向保持原始收敛速度。
</details>

**练习11.6**：分析非凸渲染优化的鞍点
考虑带镜面反射的场景，目标函数包含多个局部最优。分析鞍点的性质并设计逃逸策略。

<details>
<summary>提示</summary>
研究Hessian的负特征值，考虑随机扰动或二阶方法。
</details>

<details>
<summary>答案</summary>

鞍点处Hessian有正负特征值混合。设：
$$\mathbf{H} = \sum_{i: \lambda_i > 0} \lambda_i \mathbf{v}_i \mathbf{v}_i^T + \sum_{j: \lambda_j < 0} \lambda_j \mathbf{v}_j \mathbf{v}_j^T$$

逃逸策略：

1. **特征向量方法**：
   沿最负特征向量方向移动：
   $$\boldsymbol{\theta}_{k+1} = \boldsymbol{\theta}_k - \epsilon \mathbf{v}_{\min}$$

2. **扰动梯度下降**：
   $$\boldsymbol{\theta}_{k+1} = \boldsymbol{\theta}_k - \alpha \nabla f + \sqrt{2\alpha\beta} \boldsymbol{\xi}$$
   其中 $\boldsymbol{\xi} \sim \mathcal{N}(0, \mathbf{I})$

3. **信赖域修正**：
   解修正的牛顿系统：
   $$(\mathbf{H} + \mu\mathbf{I})\mathbf{d} = -\nabla f$$
   选择 $\mu > |\lambda_{\min}|$ 确保正定性。

对于镜面反射，鞍点常出现在高光边缘，这些策略帮助优化器探索参数空间。
</details>

**练习11.7**：推导变分贝叶斯逆向渲染
使用平均场近似，推导材质和光照联合估计的变分更新方程。

<details>
<summary>提示</summary>
分解后验 $q(\boldsymbol{\theta}_m, \boldsymbol{\theta}_l) = q(\boldsymbol{\theta}_m)q(\boldsymbol{\theta}_l)$，最小化KL散度。
</details>

<details>
<summary>答案</summary>

变分目标（ELBO）：
$$\mathcal{L} = \mathbb{E}_q[\log p(\mathbf{I}|\boldsymbol{\theta}_m, \boldsymbol{\theta}_l)] - \text{KL}[q(\boldsymbol{\theta}_m)||p(\boldsymbol{\theta}_m)] - \text{KL}[q(\boldsymbol{\theta}_l)||p(\boldsymbol{\theta}_l)]$$

坐标上升更新：

材质更新：
$$\log q(\boldsymbol{\theta}_m) = \mathbb{E}_{q(\boldsymbol{\theta}_l)}[\log p(\mathbf{I}|\boldsymbol{\theta}_m, \boldsymbol{\theta}_l)] + \log p(\boldsymbol{\theta}_m) + \text{const}$$

光照更新：
$$\log q(\boldsymbol{\theta}_l) = \mathbb{E}_{q(\boldsymbol{\theta}_m)}[\log p(\mathbf{I}|\boldsymbol{\theta}_m, \boldsymbol{\theta}_l)] + \log p(\boldsymbol{\theta}_l) + \text{const}$$

对于高斯近似：
$$q(\boldsymbol{\theta}_m) = \mathcal{N}(\boldsymbol{\mu}_m, \boldsymbol{\Sigma}_m)$$
$$q(\boldsymbol{\theta}_l) = \mathcal{N}(\boldsymbol{\mu}_l, \boldsymbol{\Sigma}_l)$$

均值更新类似MAP，协方差捕获不确定性：
$$\boldsymbol{\Sigma}_m^{-1} = -\mathbb{E}_{q(\boldsymbol{\theta}_l)}[\nabla^2_{\boldsymbol{\theta}_m} \log p(\mathbf{I}|\boldsymbol{\theta}_m, \boldsymbol{\theta}_l)] - \nabla^2 \log p(\boldsymbol{\theta}_m)$$
</details>

**练习11.8**：开放问题 - 神经网络正则化
现代逆向渲染常用神经网络表示场景。讨论如何将传统正则化概念（平滑性、稀疏性）转化为网络架构和训练策略。考虑：
- 架构偏置（如卷积、注意力）
- 隐式正则化（如dropout、谱归一化）
- 与经典方法的联系

<details>
<summary>思考方向</summary>
- 卷积 = 平移不变性先验
- Skip连接 = 多尺度正则化
- 权重衰减 = Tikhonov正则化
- Dropout = 贝叶斯模型平均
- 谱归一化 = Lipschitz约束
</details>

## 常见陷阱与错误

1. **数值精度问题**
   - 错误：直接计算 $T(t) = \exp(-\int_0^t \sigma \, ds)$ 可能下溢
   - 正确：使用对数空间计算或截断小值

2. **梯度消失/爆炸**
   - 错误：忽视深度网络中的梯度范数
   - 正确：梯度裁剪、归一化、仔细的初始化

3. **局部最优陷阱**
   - 错误：单一初始化，陷入次优解
   - 正确：多重启策略，渐进优化（由粗到细）

4. **正则化参数选择**
   - 错误：固定 $\lambda$ 值
   - 正确：交叉验证、L曲线方法、贝叶斯优化

5. **离散化误差**
   - 错误：忽视积分离散化带来的偏差
   - 正确：高阶积分方法、自适应采样

6. **不可微操作**
   - 错误：硬阈值、排序等破坏梯度流
   - 正确：软化近似、重参数化技巧

## 最佳实践检查清单

### 问题设置
- [ ] 明确定义参数空间和观测模型
- [ ] 分析问题的不适定程度
- [ ] 选择合适的损失函数（L2、感知、对抗）
- [ ] 设计物理一致的参数化

### 正则化设计
- [ ] 根据先验知识选择正则化类型
- [ ] 使用验证集调整正则化强度
- [ ] 考虑多尺度/自适应正则化
- [ ] 保持物理约束（能量守恒等）

### 优化策略
- [ ] 分析目标函数的凸性
- [ ] 选择合适的优化算法
- [ ] 实现梯度检查和数值稳定性测试
- [ ] 使用渐进策略（分辨率、复杂度）

### 实现考虑
- [ ] 利用自动微分框架
- [ ] 优化内存使用（检查点、梯度累积）
- [ ] 并行化可能的计算
- [ ] 监控收敛指标和中间结果

### 验证与调试
- [ ] 合成数据验证算法正确性
- [ ] 分析失败案例和模式
- [ ] 量化不确定性
- [ ] 与基准方法比较

# 第14章：神经逆向渲染

神经逆向渲染将深度学习与物理渲染原理相结合，从观测图像中推断场景的三维结构、材质属性和光照条件。本章探讨如何利用神经表示解决渲染方程的逆问题，包括分解表示、生成先验的引入、少样本学习以及实时系统的设计。我们将看到，通过巧妙地参数化和约束，可以将不适定的逆向问题转化为可优化的目标函数。

## 学习目标

完成本章后，您将能够：
1. 设计神经场的逆优化框架，理解梯度流与收敛性
2. 构建分解表示，分离几何、材质和光照的相互作用
3. 利用生成模型作为先验知识，正则化逆向重建
4. 实现少样本重建算法，从有限观测中恢复完整场景
5. 分析实时逆向渲染系统的设计权衡

## 14.1 神经场的逆优化

### 14.1.1 从正向渲染到逆向推断

在第6-9章中，我们学习了神经辐射场的正向渲染：给定场景表示 $\Theta$，生成图像 $I = \mathcal{R}(\Theta)$。逆向渲染则是求解逆问题：

$$\Theta^* = \arg\min_\Theta \mathcal{L}(\mathcal{R}(\Theta), I_{\text{obs}}) + \lambda \mathcal{R}_{\text{reg}}(\Theta)$$

其中 $I_{\text{obs}}$ 是观测图像，$\mathcal{L}$ 是重建损失，$\mathcal{R}_{\text{reg}}$ 是正则化项。

对于神经场表示，参数 $\Theta$ 包含：
- 网络权重 $\mathbf{W} = \{W_i, b_i\}_{i=1}^L$
- 位置编码参数 $\gamma$
- 特征网格或其他显式表示

### 14.1.2 神经场参数化的梯度计算

考虑标准的体积渲染方程：

$$C(\mathbf{r}) = \int_{t_n}^{t_f} T(t) \sigma(\mathbf{r}(t)) \mathbf{c}(\mathbf{r}(t), \mathbf{d}) dt$$

其中透射率 $T(t) = \exp\left(-\int_{t_n}^t \sigma(\mathbf{r}(s)) ds\right)$。

对网络参数 $\mathbf{W}$ 的梯度通过链式法则计算：

$$\frac{\partial \mathcal{L}}{\partial \mathbf{W}} = \sum_{\mathbf{r}} \frac{\partial \mathcal{L}}{\partial C(\mathbf{r})} \frac{\partial C(\mathbf{r})}{\partial \mathbf{W}}$$

关键在于计算 $\frac{\partial C(\mathbf{r})}{\partial \mathbf{W}}$：

$$\frac{\partial C(\mathbf{r})}{\partial \mathbf{W}} = \int_{t_n}^{t_f} \left[ \frac{\partial T(t)}{\partial \mathbf{W}} \sigma(\mathbf{r}(t)) \mathbf{c}(\mathbf{r}(t), \mathbf{d}) + T(t) \frac{\partial}{\partial \mathbf{W}}[\sigma(\mathbf{r}(t)) \mathbf{c}(\mathbf{r}(t), \mathbf{d})] \right] dt$$

其中透射率的梯度涉及路径积分：

$$\frac{\partial T(t)}{\partial \mathbf{W}} = -T(t) \int_{t_n}^t \frac{\partial \sigma(\mathbf{r}(s))}{\partial \mathbf{W}} ds$$

### 14.1.3 多视图一致性约束

给定 $N$ 个视图 $\{I_i\}_{i=1}^N$ 及其相机参数 $\{\mathbf{P}_i\}_{i=1}^N$，多视图重建损失为：

$$\mathcal{L}_{\text{multi}} = \sum_{i=1}^N \sum_{\mathbf{p} \in I_i} \rho\left( \|\mathcal{R}(\Theta, \mathbf{P}_i, \mathbf{p}) - I_i(\mathbf{p})\|_2 \right)$$

其中 $\rho$ 是鲁棒损失函数（如 Huber loss）。

**几何一致性**通过深度图约束实现：

$$\mathcal{L}_{\text{depth}} = \sum_{i,j} \sum_{\mathbf{p}} w_{ij}(\mathbf{p}) \|D_i(\mathbf{p}) - \Pi_{ij}(D_j(\Pi_{ji}(\mathbf{p})))\|$$

其中 $\Pi_{ij}$ 是从视图 $j$ 到视图 $i$ 的投影变换，$w_{ij}$ 是可见性权重。

### 14.1.4 正则化策略与收敛性分析

**密度正则化**：促进紧凑的几何表示

$$\mathcal{R}_{\text{density}} = \int_{\mathcal{V}} \sigma(\mathbf{x}) \log \sigma(\mathbf{x}) d\mathbf{x}$$

这是信息熵的负值，鼓励稀疏的密度分布。

**平滑性正则化**：

$$\mathcal{R}_{\text{smooth}} = \int_{\mathcal{V}} \|\nabla \sigma(\mathbf{x})\|^2 + \|\nabla \mathbf{c}(\mathbf{x})\|^2 d\mathbf{x}$$

**收敛性分析**：考虑梯度下降 $\Theta_{k+1} = \Theta_k - \eta_k \nabla_\Theta \mathcal{L}$

在适当的条件下（Lipschitz 连续梯度，凸性），收敛速率为：

$$\mathcal{L}(\Theta_k) - \mathcal{L}(\Theta^*) \leq \frac{\|\Theta_0 - \Theta^*\|^2}{2\sum_{i=0}^{k-1} \eta_i}$$

对于非凸神经场，通常只能保证收敛到局部极小值。

## 14.2 分解表示：几何、材质、光照

### 14.2.1 内在图像分解理论

场景的外观可以分解为内在成分的乘积：

$$I(\mathbf{x}) = A(\mathbf{x}) \cdot S(\mathbf{x}) \cdot \text{vis}(\mathbf{x})$$

其中：
- $A(\mathbf{x})$：反照率（材质颜色）
- $S(\mathbf{x})$：着色（光照效果）
- $\text{vis}(\mathbf{x})$：可见性/阴影

在神经渲染框架下，我们将这种分解扩展到体积表示。

### 14.2.2 神经场中的分解架构

**分解的神经辐射场**表示为多个专门网络：

1. **几何网络** $f_\sigma: \mathbb{R}^3 \rightarrow \mathbb{R}^+$
   $$\sigma(\mathbf{x}) = f_\sigma(\gamma_{\text{geo}}(\mathbf{x}); \Theta_\sigma)$$

2. **材质网络** $f_A: \mathbb{R}^3 \rightarrow \mathbb{R}^3$
   $$\mathbf{a}(\mathbf{x}) = f_A(\gamma_{\text{mat}}(\mathbf{x}); \Theta_A)$$

3. **光照网络** $f_L: \mathbb{R}^3 \times \mathbb{S}^2 \rightarrow \mathbb{R}^3$
   $$\mathbf{l}(\mathbf{x}, \boldsymbol{\omega}) = f_L(\mathbf{x}, \boldsymbol{\omega}; \Theta_L)$$

最终颜色通过渲染方程计算：

$$\mathbf{c}(\mathbf{x}, \boldsymbol{\omega}_o) = \mathbf{a}(\mathbf{x}) \odot \int_{\mathbb{S}^2} f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o, \mathbf{n}) \mathbf{l}(\mathbf{x}, \boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n})^+ d\boldsymbol{\omega}_i$$

其中 $f_r$ 是 BRDF，$\mathbf{n}$ 是表面法线（从密度梯度计算）。

### 14.2.3 物理约束与解耦损失

**反照率一致性**：相同材质在不同光照下应保持恒定

$$\mathcal{L}_{\text{albedo}} = \sum_{i,j} \sum_{\mathbf{x} \in \mathcal{M}} \|\mathbf{a}_i(\mathbf{x}) - \mathbf{a}_j(\mathbf{x})\|^2$$

其中 $\mathcal{M}$ 是跨视图对应的表面点集。

**光照平滑性**：自然光照通常是低频的

$$\mathcal{L}_{\text{light}} = \int_{\mathbb{S}^2} \|\nabla_{\boldsymbol{\omega}} \mathbf{l}(\mathbf{x}, \boldsymbol{\omega})\|^2 d\boldsymbol{\omega}$$

**白平衡约束**：场景平均反照率接近灰色

$$\mathcal{L}_{\text{white}} = \left\| \frac{1}{|\mathcal{V}|} \int_{\mathcal{V}} \mathbf{a}(\mathbf{x}) d\mathbf{x} - \mathbf{a}_{\text{ref}} \right\|^2$$

### 14.2.4 歧义性与规范化

分解问题存在固有歧义性：
- **尺度歧义**：$(k\mathbf{a}, \mathbf{l}/k)$ 产生相同外观
- **颜色偏移**：全局颜色变换的不确定性

**规范化策略**：

1. **固定尺度**：约束 $\|\mathbf{a}\|_\infty \leq 1$
2. **参考锚点**：指定某些点的已知反照率
3. **统计先验**：利用自然图像统计

$$\mathcal{L}_{\text{prior}} = D_{\text{KL}}(p(\mathbf{a}) \| p_{\text{natural}}(\mathbf{a}))$$

## 14.3 生成模型作为先验

### 14.3.1 变分自编码器（VAE）先验

将场景表示编码到潜在空间 $\mathbf{z} \in \mathbb{R}^d$：

$$q(\mathbf{z}|\Theta) = \mathcal{N}(\mu_\phi(\Theta), \Sigma_\phi(\Theta))$$

解码器生成场景参数：

$$p(\Theta|\mathbf{z}) = \mathcal{N}(\mu_\psi(\mathbf{z}), \Sigma_\psi(\mathbf{z}))$$

VAE 损失包含重建项和 KL 散度：

$$\mathcal{L}_{\text{VAE}} = \mathbb{E}_{q(\mathbf{z}|\Theta)}[\log p(I_{\text{obs}}|\Theta)] - D_{\text{KL}}(q(\mathbf{z}|\Theta) \| p(\mathbf{z}))$$

其中先验 $p(\mathbf{z}) = \mathcal{N}(0, \mathbf{I})$。

### 14.3.2 扩散模型在逆向渲染中的应用

扩散模型通过逐步去噪过程生成数据：

$$\mathbf{x}_t = \sqrt{\alpha_t} \mathbf{x}_0 + \sqrt{1-\alpha_t} \boldsymbol{\epsilon}$$

在逆向渲染中，我们使用条件扩散模型：

$$p_\theta(\mathbf{x}_{t-1}|\mathbf{x}_t, I_{\text{obs}}) = \mathcal{N}(\mu_\theta(\mathbf{x}_t, t, I_{\text{obs}}), \Sigma_\theta(t))$$

**分数匹配目标**：

$$\mathcal{L}_{\text{diffusion}} = \mathbb{E}_{t,\mathbf{x}_0,\boldsymbol{\epsilon}}\left[\|\boldsymbol{\epsilon} - \boldsymbol{\epsilon}_\theta(\mathbf{x}_t, t, I_{\text{obs}})\|^2\right]$$

### 14.3.3 条件生成与后验采样

给定观测 $I_{\text{obs}}$，我们需要采样后验分布：

$$p(\Theta|I_{\text{obs}}) \propto p(I_{\text{obs}}|\Theta) p(\Theta)$$

使用 **Langevin 动力学**：

$$\Theta_{k+1} = \Theta_k + \frac{\eta}{2} \nabla_\Theta \log p(\Theta_k|I_{\text{obs}}) + \sqrt{\eta} \boldsymbol{\xi}_k$$

其中 $\boldsymbol{\xi}_k \sim \mathcal{N}(0, \mathbf{I})$。

梯度可以分解为：

$$\nabla_\Theta \log p(\Theta|I_{\text{obs}}) = \nabla_\Theta \log p(I_{\text{obs}}|\Theta) + \nabla_\Theta \log p(\Theta)$$

### 14.3.4 先验强度与重建保真度的权衡

总体目标函数：

$$\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{recon}} + \lambda_{\text{prior}} \mathcal{L}_{\text{prior}} + \lambda_{\text{reg}} \mathcal{L}_{\text{reg}}$$

**自适应权重策略**：

$$\lambda_{\text{prior}}(k) = \lambda_0 \exp(-k/\tau)$$

随着优化进行，逐渐减少先验影响，允许更精确的重建。

**不确定性加权**：

$$w(\mathbf{x}) = \frac{1}{\sigma^2_{\text{recon}}(\mathbf{x}) + \sigma^2_{\text{prior}}(\mathbf{x})}$$

在不确定区域更依赖先验。

## 14.4 少样本重建

### 14.4.1 元学习框架

少样本重建的目标是从极少的观测（通常 1-5 张图像）恢复完整的三维场景。元学习提供了一个强大的框架：

**Model-Agnostic Meta-Learning (MAML)** 应用于神经场：

初始化参数 $\Theta_0$ 通过多任务学习获得：

$$\Theta_0 = \arg\min_\Theta \sum_{i=1}^{N_{\text{tasks}}} \mathcal{L}_i(\Theta - \alpha \nabla_\Theta \mathcal{L}_i^{\text{support}}(\Theta))$$

其中 $\mathcal{L}_i^{\text{support}}$ 是任务 $i$ 在支持集上的损失。

对于新场景，快速适应：

$$\Theta_{\text{adapted}} = \Theta_0 - \alpha \sum_{k=1}^{K} \nabla_\Theta \mathcal{L}^{\text{query}}(\Theta)$$

**条件神经过程 (CNP)** 框架：

编码器聚合上下文信息：

$$\mathbf{r} = h\left(\frac{1}{N}\sum_{i=1}^N g(\mathbf{x}_i, I(\mathbf{x}_i))\right)$$

解码器基于全局表示预测：

$$p(I(\mathbf{x}_*)|I_{1:N}) = \mathcal{N}(\mu_\theta(\mathbf{x}_*, \mathbf{r}), \sigma^2_\theta(\mathbf{x}_*, \mathbf{r}))$$

### 14.4.2 神经表示的快速适应

**超网络方法**：使用超网络 $H$ 生成场景特定的权重：

$$\mathbf{W}_{\text{scene}} = H(\mathbf{z}_{\text{context}}; \Phi)$$

其中 $\mathbf{z}_{\text{context}}$ 是从观测图像提取的上下文向量。

**模块化适应**：将网络分为共享模块和适应模块：

$$f(\mathbf{x}) = f_{\text{adapt}}(f_{\text{shared}}(\mathbf{x}; \Theta_{\text{shared}}); \Theta_{\text{adapt}})$$

只更新 $\Theta_{\text{adapt}}$，保持 $\Theta_{\text{shared}}$ 固定。

**低秩适应 (LoRA)**：

$$\mathbf{W}_{\text{adapted}} = \mathbf{W}_0 + \Delta \mathbf{W} = \mathbf{W}_0 + \mathbf{B}\mathbf{A}$$

其中 $\mathbf{B} \in \mathbb{R}^{d \times r}$，$\mathbf{A} \in \mathbb{R}^{r \times k}$，$r \ll \min(d, k)$。

### 14.4.3 稀疏视图综合

**几何先验**：利用深度先验网络：

$$D_{\text{prior}}(\mathbf{x}) = f_{\text{depth}}(I_{\text{obs}}, \mathbf{x}; \Theta_{\text{depth}})$$

融合到密度预测：

$$\sigma(\mathbf{x}) = \sigma_{\text{neural}}(\mathbf{x}) + \lambda \cdot \delta(d(\mathbf{x}) - D_{\text{prior}}(\mathbf{x}))$$

**多尺度特征聚合**：

$$\mathbf{F}(\mathbf{x}) = \sum_{l=1}^L w_l \cdot \text{Interp}(\mathbf{F}_l, \mathbf{x})$$

其中 $\mathbf{F}_l$ 是第 $l$ 层的特征图。

**跨视图注意力机制**：

$$\mathbf{f}_{\text{agg}}(\mathbf{x}) = \sum_{i=1}^N \text{softmax}\left(\frac{q(\mathbf{x}) \cdot k_i(\Pi_i(\mathbf{x}))}{\sqrt{d}}\right) v_i(\Pi_i(\mathbf{x}))$$

其中 $q$、$k$、$v$ 是查询、键、值变换。

### 14.4.4 不确定性量化

**认知不确定性**：通过集成方法估计：

$$\mu(\mathbf{x}) = \frac{1}{M}\sum_{m=1}^M f_m(\mathbf{x})$$

$$\sigma^2_{\text{epistemic}}(\mathbf{x}) = \frac{1}{M}\sum_{m=1}^M (f_m(\mathbf{x}) - \mu(\mathbf{x}))^2$$

**偶然不确定性**：直接预测：

$$f(\mathbf{x}) = (\mu(\mathbf{x}), \sigma^2_{\text{aleatoric}}(\mathbf{x}))$$

**总不确定性**：

$$\sigma^2_{\text{total}}(\mathbf{x}) = \sigma^2_{\text{epistemic}}(\mathbf{x}) + \sigma^2_{\text{aleatoric}}(\mathbf{x})$$

**主动视图选择**：基于信息增益选择下一个最佳视图：

$$\mathbf{v}_{\text{next}} = \arg\max_{\mathbf{v}} \mathbb{E}[H(\Theta) - H(\Theta|I_{\mathbf{v}})]$$

其中 $H$ 是熵函数。

## 14.5 实时逆向渲染系统

### 14.5.1 轻量级神经表示

**哈希编码**：使用多分辨率哈希表存储特征：

$$\mathbf{f}(\mathbf{x}) = \bigoplus_{l=1}^L \text{HashTable}_l[\text{hash}(\lfloor \mathbf{x} \cdot 2^l \rfloor)]$$

内存复杂度：$O(L \cdot T \cdot F)$，其中 $T$ 是表大小，$F$ 是特征维度。

**量化与剪枝**：

1. **权重量化**：$w_q = \text{round}(w/s) \cdot s$
2. **结构化剪枝**：移除重要性低的通道
3. **知识蒸馏**：$\mathcal{L}_{\text{distill}} = \|f_{\text{student}}(\mathbf{x}) - f_{\text{teacher}}(\mathbf{x})\|^2$

**张量分解**：使用 CP 或 Tucker 分解：

$$\mathcal{T} \approx \sum_{r=1}^R \mathbf{a}_r \otimes \mathbf{b}_r \otimes \mathbf{c}_r$$

### 14.5.2 增量优化策略

**关键帧策略**：维护关键帧集合 $\mathcal{K}$：

$$\mathcal{K}_{t+1} = \begin{cases}
\mathcal{K}_t \cup \{I_t\} & \text{if } d(I_t, \mathcal{K}_t) > \tau \\
\mathcal{K}_t & \text{otherwise}
\end{cases}$$

**局部更新**：只更新观测区域的参数：

$$\Theta_{t+1}(\mathbf{x}) = \begin{cases}
\Theta_t(\mathbf{x}) - \eta \nabla \mathcal{L} & \text{if } \mathbf{x} \in \mathcal{V}_{\text{observed}} \\
\Theta_t(\mathbf{x}) & \text{otherwise}
\end{cases}$$

**滑动窗口优化**：

$$\mathcal{L}_{\text{window}} = \sum_{i=t-W}^t w_{t-i} \mathcal{L}_i$$

其中 $w_i$ 是时间衰减权重。

### 14.5.3 硬件加速架构

**GPU 并行化策略**：

1. **射线并行**：每个 CUDA 线程处理一条射线
2. **体素并行**：并行评估空间网格
3. **特征并行**：在特征维度上分割计算

**专用硬件设计**：

$$\text{Throughput} = \frac{N_{\text{rays}} \times N_{\text{samples}}}{\text{Latency}_{\text{kernel}}}$$

优化目标：
- 最大化内存带宽利用率
- 最小化分支发散
- 利用张量核心加速

**混合精度训练**：

$$\mathcal{L}_{\text{fp16}} = \text{scale} \times \mathcal{L}_{\text{original}}$$

使用动态损失缩放防止梯度下溢。

### 14.5.4 交互式编辑与反馈

**实时预览**：使用低分辨率代理：

$$I_{\text{preview}} = \text{Upsample}(\mathcal{R}_{\text{low}}(\Theta))$$

**增量细化**：

$$\Theta_{k+1} = \Theta_k + \alpha \cdot \text{UserEdit}(k)$$

**约束传播**：用户编辑作为硬约束：

$$\mathcal{L}_{\text{edit}} = \sum_{\mathbf{x} \in \mathcal{E}} \|\Theta(\mathbf{x}) - \Theta_{\text{target}}(\mathbf{x})\|^2$$

**性能指标**：
- 延迟：$< 16.7$ ms（60 FPS）
- 吞吐量：$> 10^6$ 射线/秒
- 内存：$< 4$ GB VRAM

## 本章小结

神经逆向渲染通过结合深度学习和物理渲染原理，实现了从图像到三维场景的逆向推断。关键贡献包括：

1. **逆优化框架**：将渲染方程的逆问题转化为可微分优化，通过梯度下降求解场景参数
2. **分解表示**：通过物理约束和解耦损失，分离几何、材质和光照的内在成分
3. **生成先验**：利用 VAE 和扩散模型编码场景统计，正则化不适定的逆问题
4. **少样本学习**：通过元学习和条件神经过程，从极少观测恢复完整场景
5. **实时系统**：通过轻量级表示和硬件优化，实现交互式逆向渲染

核心公式回顾：
- 逆向渲染目标：$\Theta^* = \arg\min_\Theta \mathcal{L}(\mathcal{R}(\Theta), I_{\text{obs}}) + \lambda \mathcal{R}_{\text{reg}}(\Theta)$
- 分解表示：$\mathbf{c}(\mathbf{x}, \boldsymbol{\omega}_o) = \mathbf{a}(\mathbf{x}) \odot \int_{\mathbb{S}^2} f_r \mathbf{l}(\mathbf{x}, \boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n})^+ d\boldsymbol{\omega}_i$
- 后验采样：$p(\Theta|I_{\text{obs}}) \propto p(I_{\text{obs}}|\Theta) p(\Theta)$

## 练习题

### 基础题

**练习 14.1**：推导体积渲染中透射率 $T(t)$ 对网络参数 $\mathbf{W}$ 的梯度表达式。

*提示*：使用链式法则和路径积分的导数。

<details>
<summary>答案</summary>

从透射率定义：
$$T(t) = \exp\left(-\int_{t_n}^t \sigma(\mathbf{r}(s)) ds\right)$$

对 $\mathbf{W}$ 求导：
$$\frac{\partial T(t)}{\partial \mathbf{W}} = T(t) \cdot \frac{\partial}{\partial \mathbf{W}}\left[-\int_{t_n}^t \sigma(\mathbf{r}(s)) ds\right]$$

$$= -T(t) \int_{t_n}^t \frac{\partial \sigma(\mathbf{r}(s))}{\partial \mathbf{W}} ds$$

这表明透射率梯度与路径上所有点的密度梯度相关。
</details>

**练习 14.2**：证明在分解表示中，尺度歧义 $(k\mathbf{a}, \mathbf{l}/k)$ 保持渲染结果不变。

*提示*：将缩放后的值代入渲染方程。

<details>
<summary>答案</summary>

原始渲染：
$$\mathbf{c} = \mathbf{a} \odot \mathbf{l}$$

缩放后：
$$\mathbf{c}' = (k\mathbf{a}) \odot (\mathbf{l}/k) = k \cdot \frac{1}{k} \cdot \mathbf{a} \odot \mathbf{l} = \mathbf{a} \odot \mathbf{l} = \mathbf{c}$$

因此渲染结果保持不变，这就是尺度歧义的来源。
</details>

**练习 14.3**：计算哈希编码的内存需求，给定 $L=16$ 层，每层表大小 $T=2^{19}$，特征维度 $F=2$。

*提示*：总内存 = 层数 × 表大小 × 特征维度 × 每个浮点数字节。

<details>
<summary>答案</summary>

内存需求：
$$\text{Memory} = L \times T \times F \times 4 \text{ bytes}$$
$$= 16 \times 2^{19} \times 2 \times 4$$
$$= 16 \times 524288 \times 8$$
$$= 67,108,864 \text{ bytes} \approx 64 \text{ MB}$$

这是相对紧凑的表示，适合实时应用。
</details>

### 挑战题

**练习 14.4**：设计一个联合优化框架，同时进行场景重建和相机标定。写出目标函数和优化策略。

*提示*：考虑相机参数 $\mathbf{P}$ 和场景参数 $\Theta$ 的联合优化。

<details>
<summary>答案</summary>

联合优化目标：
$$(\Theta^*, \mathbf{P}^*) = \arg\min_{\Theta, \mathbf{P}} \sum_{i=1}^N \mathcal{L}_{\text{photo}}(\mathcal{R}(\Theta, \mathbf{P}_i), I_i) + \lambda_1 \mathcal{R}_{\text{scene}}(\Theta) + \lambda_2 \mathcal{R}_{\text{camera}}(\mathbf{P})$$

其中相机正则化：
$$\mathcal{R}_{\text{camera}}(\mathbf{P}) = \sum_{i,j} \|\mathbf{P}_i - \mathbf{P}_j\|^2_{\text{smooth}} + \sum_i \|\mathbf{P}_i - \mathbf{P}_i^{\text{init}}\|^2$$

交替优化策略：
1. 固定 $\mathbf{P}$，优化 $\Theta$：$\Theta_{k+1} = \Theta_k - \eta_\Theta \nabla_\Theta \mathcal{L}$
2. 固定 $\Theta$，优化 $\mathbf{P}$：$\mathbf{P}_{k+1} = \mathbf{P}_k - \eta_P \nabla_P \mathcal{L}$
3. 重复直到收敛

关键是平衡两者的学习率，防止退化解。
</details>

**练习 14.5**：推导条件扩散模型在逆向渲染中的分数函数 $\nabla_{\mathbf{x}_t} \log p(\mathbf{x}_t|I_{\text{obs}})$。

*提示*：使用贝叶斯定理和分数匹配。

<details>
<summary>答案</summary>

从贝叶斯定理：
$$p(\mathbf{x}_t|I_{\text{obs}}) = \frac{p(I_{\text{obs}}|\mathbf{x}_t) p(\mathbf{x}_t)}{p(I_{\text{obs}})}$$

取对数梯度：
$$\nabla_{\mathbf{x}_t} \log p(\mathbf{x}_t|I_{\text{obs}}) = \nabla_{\mathbf{x}_t} \log p(I_{\text{obs}}|\mathbf{x}_t) + \nabla_{\mathbf{x}_t} \log p(\mathbf{x}_t)$$

第二项是无条件分数，由预训练模型提供：
$$\nabla_{\mathbf{x}_t} \log p(\mathbf{x}_t) = -\frac{\boldsymbol{\epsilon}_\theta(\mathbf{x}_t, t)}{\sqrt{1-\alpha_t}}$$

第一项需要通过渲染计算：
$$\nabla_{\mathbf{x}_t} \log p(I_{\text{obs}}|\mathbf{x}_t) = -\frac{1}{\sigma^2} \nabla_{\mathbf{x}_t} \|\mathcal{R}(\mathbf{x}_t) - I_{\text{obs}}\|^2$$

这需要可微渲染器 $\mathcal{R}$。
</details>

**练习 14.6**：分析元学习 MAML 在神经场中的计算复杂度，包括前向传播和反向传播。

*提示*：考虑内外循环的梯度计算。

<details>
<summary>答案</summary>

设网络有 $P$ 个参数，批大小为 $B$，内循环步数为 $K$。

**内循环（适应）**：
- 前向：$O(KBP)$
- 反向：$O(KBP)$

**外循环（元优化）**：
- 需要计算二阶导数
- 前向：$O(BP)$
- 反向：$O(KBP^2)$（由于需要通过内循环梯度）

**总复杂度**：$O(N_{\text{tasks}} \cdot K \cdot B \cdot P^2)$

内存需求：需要存储计算图，$O(KP)$

优化：使用一阶近似（FOMAML）降低到 $O(KBP)$。
</details>

**练习 14.7**（开放题）：设计一个结合物理仿真的神经逆向渲染系统，用于推断动态场景的物理参数（如质量、弹性系数）。

*提示*：考虑如何将物理约束集成到优化框架中。

<details>
<summary>答案</summary>

系统设计：

1. **表示**：
   - 静态几何：$\Theta_{\text{geo}}$
   - 物理参数：$\Theta_{\text{phys}} = \{m_i, k_i, c_i\}$（质量、刚度、阻尼）
   - 动态状态：$\mathbf{x}_i(t), \mathbf{v}_i(t)$

2. **物理约束**：
   $$\mathbf{F} = m\mathbf{a} = -k\mathbf{x} - c\mathbf{v}$$

3. **优化目标**：
   $$\mathcal{L} = \mathcal{L}_{\text{render}} + \lambda_1 \mathcal{L}_{\text{physics}} + \lambda_2 \mathcal{L}_{\text{temporal}}$$

   其中：
   - $\mathcal{L}_{\text{physics}} = \sum_t \|m\ddot{\mathbf{x}} + k\mathbf{x} + c\dot{\mathbf{x}}\|^2$
   - $\mathcal{L}_{\text{temporal}} = \sum_t \|\mathbf{x}_{t+1} - \text{Integrate}(\mathbf{x}_t, \mathbf{v}_t, \Theta_{\text{phys}})\|^2$

4. **可微分仿真**：使用隐式积分器保证稳定性

5. **应用**：材质识别、碰撞预测、机器人抓取规划
</details>

**练习 14.8**：证明在少样本设置下，使用 $N$ 个视图时重建误差的期望界限。

*提示*：使用 PAC 学习理论。

<details>
<summary>答案</summary>

设场景空间复杂度为 $\mathcal{H}$，使用 $N$ 个独立同分布的视图。

根据 PAC 界限，以概率至少 $1-\delta$：

$$\mathbb{E}[\mathcal{L}_{\text{test}}] \leq \mathcal{L}_{\text{train}} + \sqrt{\frac{2\log(|\mathcal{H}|/\delta)}{N}} + \mathcal{B}_{\text{approx}}$$

其中：
- $\mathcal{L}_{\text{train}}$：训练误差
- $|\mathcal{H}|$：假设空间大小（与网络容量相关）
- $\mathcal{B}_{\text{approx}}$：函数逼近误差

对于神经场，$\log|\mathcal{H}| \approx O(P \log P)$，其中 $P$ 是参数数量。

因此：
$$\mathbb{E}[\mathcal{L}_{\text{test}}] = O\left(\mathcal{L}_{\text{train}} + \sqrt{\frac{P \log P}{N}}\right)$$

这表明需要 $N = \Omega(P \log P)$ 个视图才能获得良好的泛化。
</details>

## 常见陷阱与错误

1. **梯度消失/爆炸**
   - 问题：深层网络中梯度不稳定
   - 解决：使用梯度裁剪、层归一化、残差连接

2. **局部极小值**
   - 问题：非凸优化陷入次优解
   - 解决：多次随机初始化、课程学习、模拟退火

3. **分解歧义**
   - 问题：材质-光照分离不唯一
   - 解决：强物理先验、多视图约束、参考标定

4. **过拟合**
   - 问题：少样本设置下记忆训练视图
   - 解决：正则化、数据增强、早停

5. **计算效率**
   - 问题：实时要求与精度冲突
   - 解决：自适应采样、层级表示、混合精度

6. **数值稳定性**
   - 问题：体积渲染积分的数值误差
   - 解决：对数空间计算、稳定的积分方案

## 最佳实践检查清单

### 设计阶段
- [ ] 明确逆向渲染目标（几何/材质/光照）
- [ ] 选择合适的神经表示（隐式/显式/混合）
- [ ] 设计物理合理的分解架构
- [ ] 确定先验知识的引入方式
- [ ] 规划计算资源和实时性要求

### 实现阶段
- [ ] 实现稳定的梯度计算和反向传播
- [ ] 添加必要的正则化项
- [ ] 设置合理的初始化策略
- [ ] 实现高效的采样和积分方案
- [ ] 优化内存使用和计算并行性

### 评估阶段
- [ ] 定量评估重建精度（PSNR, SSIM, LPIPS）
- [ ] 测试泛化能力（新视角、新光照）
- [ ] 分析计算性能（FPS, 内存占用）
- [ ] 验证物理合理性（能量守恒、互易性）
- [ ] 用户研究（如适用）

### 部署阶段
- [ ] 模型压缩和量化
- [ ] 硬件适配优化
- [ ] 鲁棒性测试（噪声、遮挡）
- [ ] 增量更新机制
- [ ] 用户交互接口

---

*下一章：[第15章：标量波动光学基础](chapter15.md)*
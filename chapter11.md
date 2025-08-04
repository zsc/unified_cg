# 第11章：逆向渲染的数学基础

逆向渲染将计算机图形学的前向问题颠倒过来：给定观测图像，我们寻求恢复场景的底层参数——几何、材质和光照。这一逆问题在计算机视觉、机器学习和物理模拟的交叉点上，形成了可微渲染和神经渲染等现代技术的数学基础。本章建立了解决这些逆问题的严格数学框架，从基本的不适定性分析到实用的优化策略。

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

对于体积渲染，我们有：

$$L(\mathbf{r}, \boldsymbol{\omega}) = \int_0^t T(s) \sigma_s(\mathbf{r}(s)) L_s(\mathbf{r}(s), \boldsymbol{\omega}) \, ds + T(t) L_{\text{bg}}$$

其中透射率 $T(s) = \exp\left(-\int_0^s \sigma_t(\mathbf{r}(u)) \, du\right)$。

### 逆问题形式化

设 $\boldsymbol{\theta}$ 表示所有场景参数（几何、材质、光照），$\mathbf{I}_{\text{obs}}$ 为观测图像。逆向渲染寻求：

$$\boldsymbol{\theta}^* = \arg\min_{\boldsymbol{\theta}} \mathcal{L}(\mathcal{R}(\boldsymbol{\theta}), \mathbf{I}_{\text{obs}}) + \mathcal{R}_{\text{reg}}(\boldsymbol{\theta})$$

其中：
- $\mathcal{R}$：渲染算子
- $\mathcal{L}$：损失函数（如 $L^2$ 或感知损失）
- $\mathcal{R}_{\text{reg}}$：正则化项

### 参数空间分解

我们通常将 $\boldsymbol{\theta}$ 分解为：

$$\boldsymbol{\theta} = (\boldsymbol{\theta}_g, \boldsymbol{\theta}_m, \boldsymbol{\theta}_l)$$

- **几何参数** $\boldsymbol{\theta}_g$：顶点位置、法线、SDF值
- **材质参数** $\boldsymbol{\theta}_m$：反照率、粗糙度、折射率
- **光照参数** $\boldsymbol{\theta}_l$：环境贴图、点光源位置/强度

### 观测模型

实际观测包含噪声：

$$\mathbf{I}_{\text{obs}} = \mathcal{R}(\boldsymbol{\theta}_{\text{true}}) + \boldsymbol{\epsilon}$$

其中 $\boldsymbol{\epsilon}$ 表示传感器噪声、量化误差等。

### 不唯一性示例

逆向渲染的根本挑战在于解的不唯一性：

**示例1：反照率-光照歧义**
给定观测亮度 $B$，无法区分：
- 暗材质 + 强光照：$\rho_1 L_1 = B$，其中 $\rho_1 = 0.2, L_1 = 5$
- 亮材质 + 弱光照：$\rho_2 L_2 = B$，其中 $\rho_2 = 1.0, L_2 = 1$

**示例2：凹凸贴图vs几何细节**
表面法线扰动可由以下产生：
- 真实几何位移
- 法线贴图
- 位移贴图

这些在单视角下可能产生相同的着色。

## 11.2 不适定性与正则化

### Hadamard不适定性

问题称为良态（well-posed）需满足：
1. 解存在
2. 解唯一
3. 解连续依赖于数据

逆向渲染通常违反条件2和3，使其成为不适定问题。

### 条件数分析

考虑线性化渲染算子 $\mathbf{J} = \partial \mathcal{R} / \partial \boldsymbol{\theta}$。条件数：

$$\kappa(\mathbf{J}) = \frac{\sigma_{\max}(\mathbf{J})}{\sigma_{\min}(\mathbf{J})}$$

大条件数表示数值不稳定性。

### Tikhonov正则化

最基本的正则化方法添加 $L^2$ 惩罚：

$$\boldsymbol{\theta}^* = \arg\min_{\boldsymbol{\theta}} \|\mathcal{R}(\boldsymbol{\theta}) - \mathbf{I}_{\text{obs}}\|^2 + \lambda \|\boldsymbol{\theta} - \boldsymbol{\theta}_0\|^2$$

其中 $\boldsymbol{\theta}_0$ 是先验估计，$\lambda$ 控制正则化强度。

### 稀疏性先验

**L1正则化**促进稀疏解：

$$\mathcal{R}_{\text{L1}}(\boldsymbol{\theta}) = \lambda_1 \|\boldsymbol{\theta}\|_1$$

**总变分（TV）正则化**保持边缘：

$$\mathcal{R}_{\text{TV}}(\mathbf{u}) = \lambda_{\text{TV}} \int_{\Omega} |\nabla \mathbf{u}| \, d\mathbf{x}$$

对于离散图像：

$$\mathcal{R}_{\text{TV}}(\mathbf{u}) = \lambda_{\text{TV}} \sum_{i,j} \sqrt{(u_{i+1,j} - u_{i,j})^2 + (u_{i,j+1} - u_{i,j})^2 + \epsilon}$$

### 物理约束

**能量守恒**：
$$\int_{\Omega} f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_i \cdot \mathbf{n}) \, d\boldsymbol{\omega}_i \leq 1$$

**互易性**：
$$f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = f_r(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)$$

**正定性**：
$$f_r(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) \geq 0$$

这些可作为硬约束或软惩罚项。

## 11.3 梯度计算：伴随法与自动微分

### 解析梯度推导

对于简单情况，可直接推导梯度。考虑Lambert着色：

$$I = \mathbf{n} \cdot \mathbf{l} \cdot \rho$$

对反照率的梯度：
$$\frac{\partial I}{\partial \rho} = \mathbf{n} \cdot \mathbf{l}$$

对法线的梯度（需要投影到切平面）：
$$\frac{\partial I}{\partial \mathbf{n}} = (\mathbf{I} - \mathbf{n}\mathbf{n}^T)(\mathbf{l} \cdot \rho)$$

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

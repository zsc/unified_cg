# 第13章：材质与几何重建

本章探讨逆向渲染的核心问题：从观察图像中恢复场景的材质属性和几何形状。我们将建立在前面章节的可微渲染基础上，展示如何通过优化框架同时重建材质和几何。这些技术在计算机视觉、虚拟现实和数字资产创建中有着广泛应用。

## 学习目标

完成本章后，您将能够：
1. 从多视图图像估计BRDF和BSSRDF参数
2. 理解并推导shape-from-shading的变分公式
3. 分析多视图立体重建中的优化景观
4. 设计联合材质-几何优化的能量函数
5. 应用物理约束提高重建质量
6. 评估不同先验对重建结果的影响

## 13.1 BRDF/BSSRDF估计

### 13.1.1 逆向渲染方程

给定观察图像 $I(\mathbf{x})$，我们寻求恢复表面反射属性。渲染方程的逆问题可表述为：

$$\min_{\rho} \sum_{\mathbf{x}} \left\| I(\mathbf{x}) - \int_{\Omega} L(\mathbf{x}, \boldsymbol{\omega}_o) d\boldsymbol{\omega}_o \right\|^2$$

其中出射辐射度由BRDF $\rho$ 决定：

$$L(\mathbf{x}, \boldsymbol{\omega}_o) = \int_{\Omega} \rho(\mathbf{x}, \boldsymbol{\omega}_i, \boldsymbol{\omega}_o) L_i(\mathbf{x}, \boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n}) d\boldsymbol{\omega}_i$$

这是一个典型的不适定逆问题，因为：
1. **非唯一性**：不同的BRDF和光照组合可能产生相同图像
2. **不稳定性**：观察中的小扰动可能导致解的大变化
3. **不完备性**：有限视角无法观察到所有反射方向

为了使问题适定，我们引入正则化项：

$$\mathcal{L}[\rho] = \sum_{\mathbf{x}} \left\| I(\mathbf{x}) - \mathcal{R}[\rho](\mathbf{x}) \right\|^2 + \lambda \mathcal{R}_{prior}[\rho]$$

其中 $\mathcal{R}[\rho]$ 表示渲染算子，$\mathcal{R}_{prior}$ 为先验正则化。

### 13.1.2 参数化BRDF模型

实际应用中，我们通常采用参数化BRDF模型：

**微表面模型**：
$$\rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \frac{F(\boldsymbol{\omega}_i, \mathbf{h}) G(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) D(\mathbf{h})}{4(\boldsymbol{\omega}_i \cdot \mathbf{n})(\boldsymbol{\omega}_o \cdot \mathbf{n})}$$

其中：
- $F$: Fresnel项，参数为折射率 $\eta$
- $G$: 几何衰减，参数为粗糙度 $\alpha$
- $D$: 法线分布，如GGX分布
- $\mathbf{h} = (\boldsymbol{\omega}_i + \boldsymbol{\omega}_o)/\|\boldsymbol{\omega}_i + \boldsymbol{\omega}_o\|$

**Fresnel项**（Schlick近似）：
$$F(\boldsymbol{\omega}_i, \mathbf{h}) = F_0 + (1 - F_0)(1 - \boldsymbol{\omega}_i \cdot \mathbf{h})^5$$

其中 $F_0 = \left(\frac{\eta - 1}{\eta + 1}\right)^2$ 为垂直入射时的反射率。

**精确Fresnel方程**：
对于非偏振光：
$$F = \frac{1}{2}(F_s + F_p)$$

其中：
$$F_s = \left|\frac{n_1\cos\theta_i - n_2\cos\theta_t}{n_1\cos\theta_i + n_2\cos\theta_t}\right|^2$$
$$F_p = \left|\frac{n_2\cos\theta_i - n_1\cos\theta_t}{n_2\cos\theta_i + n_1\cos\theta_t}\right|^2$$

利用Snell定律 $n_1\sin\theta_i = n_2\sin\theta_t$ 可得：
$$\cos\theta_t = \sqrt{1 - \left(\frac{n_1}{n_2}\right)^2\sin^2\theta_i}$$

**GGX法线分布**：
$$D_{GGX}(\mathbf{h}) = \frac{\alpha^2}{\pi((\mathbf{n} \cdot \mathbf{h})^2(\alpha^2 - 1) + 1)^2}$$

**各向异性GGX**：
$$D_{aniso}(\mathbf{h}) = \frac{1}{\pi\alpha_x\alpha_y} \frac{1}{\left(\frac{(\mathbf{h}\cdot\mathbf{t})^2}{\alpha_x^2} + \frac{(\mathbf{h}\cdot\mathbf{b})^2}{\alpha_y^2} + (\mathbf{h}\cdot\mathbf{n})^2\right)^2}$$

其中 $\mathbf{t}$ 和 $\mathbf{b}$ 为切线和副切线方向。

**Smith遮蔽函数**：
$$G(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = G_1(\boldsymbol{\omega}_i) G_1(\boldsymbol{\omega}_o)$$

其中：
$$G_1(\boldsymbol{\omega}) = \frac{2(\mathbf{n} \cdot \boldsymbol{\omega})}{(\mathbf{n} \cdot \boldsymbol{\omega}) + \sqrt{\alpha^2 + (1 - \alpha^2)(\mathbf{n} \cdot \boldsymbol{\omega})^2}}$$

**Height-correlated Smith函数**：
考虑入射和出射方向的相关性：
$$G_{corr}(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \frac{1}{1 + \Lambda(\boldsymbol{\omega}_i) + \Lambda(\boldsymbol{\omega}_o)}$$

其中：
$$\Lambda(\boldsymbol{\omega}) = \frac{\sqrt{1 + \alpha^2\tan^2\theta} - 1}{2}$$

**完整参数向量**：
$$\boldsymbol{\theta} = [\alpha, \eta, \mathbf{k}_d, \mathbf{k}_s, \mathbf{k}_e]$$

包含粗糙度、折射率、漫反射率、镜面反射率和自发光。

**Disney BRDF参数化**：
更直观的艺术家友好参数：
- baseColor: 基础颜色
- metallic: 金属度 $\in [0,1]$
- roughness: 粗糙度 $\in [0,1]$
- specular: 镜面反射强度
- specularTint: 镜面反射染色
- anisotropic: 各向异性 $\in [0,1]$
- sheen: 织物光泽
- sheenTint: 光泽染色
- clearcoat: 清漆层
- clearcoatGloss: 清漆光泽度

参数映射：
$$\alpha = roughness^2$$
$$F_0 = lerp(0.08 \cdot specular, baseColor, metallic)$$

**分层材质模型**：
$$\rho_{total} = \rho_{coat} + (1 - F_{coat}) \cdot \rho_{base}$$

其中 $F_{coat}$ 为清漆层的Fresnel反射率。

### 13.1.3 梯度计算

对BRDF参数 $\boldsymbol{\theta} = [\alpha, \eta, \mathbf{k}_d, \mathbf{k}_s]$ 的梯度：

$$\frac{\partial \mathcal{L}}{\partial \boldsymbol{\theta}} = -2\sum_{\mathbf{x}} (I(\mathbf{x}) - \hat{I}(\mathbf{x})) \frac{\partial \hat{I}(\mathbf{x})}{\partial \boldsymbol{\theta}}$$

其中渲染图像对参数的导数通过链式法则计算：

$$\frac{\partial \hat{I}}{\partial \boldsymbol{\theta}} = \int_{\Omega} \frac{\partial \rho}{\partial \boldsymbol{\theta}} L_i (\boldsymbol{\omega}_i \cdot \mathbf{n}) d\boldsymbol{\omega}_i$$

**对粗糙度的导数**：
$$\frac{\partial \rho}{\partial \alpha} = \rho \left[ \frac{\partial \ln D}{\partial \alpha} + \frac{\partial \ln G}{\partial \alpha} \right]$$

其中：
$$\frac{\partial \ln D_{GGX}}{\partial \alpha} = \frac{2}{\alpha} - \frac{4\alpha((\mathbf{n} \cdot \mathbf{h})^2 - 1)}{(\mathbf{n} \cdot \mathbf{h})^2(\alpha^2 - 1) + 1}$$

**Smith函数对粗糙度的导数**：
$$\frac{\partial G_1}{\partial \alpha} = -\frac{G_1^2(\boldsymbol{\omega})}{2(\mathbf{n} \cdot \boldsymbol{\omega})} \cdot \frac{\alpha}{\sqrt{\alpha^2 + (1 - \alpha^2)(\mathbf{n} \cdot \boldsymbol{\omega})^2}}$$

**对折射率的导数**：
$$\frac{\partial \rho}{\partial \eta} = \frac{G D}{4(\boldsymbol{\omega}_i \cdot \mathbf{n})(\boldsymbol{\omega}_o \cdot \mathbf{n})} \frac{\partial F}{\partial \eta}$$

$$\frac{\partial F}{\partial \eta} = \frac{\partial F_0}{\partial \eta}(1 - (1 - \boldsymbol{\omega}_i \cdot \mathbf{h})^5)$$

其中 $\frac{\partial F_0}{\partial \eta} = \frac{4}{(\eta + 1)^3}$。

**对漫反射系数的导数**：
对于包含漫反射项的完整BRDF：
$$\rho_{total} = \frac{\mathbf{k}_d}{\pi}(1 - F) + \rho_{spec}$$

梯度为：
$$\frac{\partial \rho_{total}}{\partial \mathbf{k}_d} = \frac{1 - F}{\pi}$$

**蒙特卡洛估计**：
实际计算中使用重要性采样：
$$\frac{\partial \hat{I}}{\partial \boldsymbol{\theta}} \approx \frac{1}{N} \sum_{i=1}^N \frac{\partial \rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o)}{\partial \boldsymbol{\theta}} \frac{L_i(\boldsymbol{\omega}_i)(\boldsymbol{\omega}_i \cdot \mathbf{n})}{p(\boldsymbol{\omega}_i)}$$

**方差减少技术**：
1. **多重重要性采样**（MIS）：
   $$\frac{\partial \hat{I}}{\partial \boldsymbol{\theta}} = \sum_{k} \frac{1}{N_k} \sum_{i=1}^{N_k} w_k(\boldsymbol{\omega}_i) \frac{\partial \rho}{\partial \boldsymbol{\theta}} \frac{L_i(\boldsymbol{\omega}_i \cdot \mathbf{n})}{p_k(\boldsymbol{\omega}_i)}$$
   
   其中权重函数：
   $$w_k(\boldsymbol{\omega}) = \frac{N_k p_k(\boldsymbol{\omega})}{\sum_j N_j p_j(\boldsymbol{\omega})}$$

2. **控制变量**：
   使用已知期望的辅助函数减少方差：
   $$\frac{\partial \hat{I}}{\partial \boldsymbol{\theta}} = \frac{\partial \hat{I}_{MC}}{\partial \boldsymbol{\theta}} - c(\frac{\partial \hat{I}_{CV}}{\partial \boldsymbol{\theta}} - E[\frac{\partial I_{CV}}{\partial \boldsymbol{\theta}}])$$

**Hessian计算**：
二阶导数对于优化算法的收敛性分析至关重要：
$$\frac{\partial^2 \mathcal{L}}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} = 2\sum_{\mathbf{x}} \left[ \frac{\partial \hat{I}}{\partial \boldsymbol{\theta}_i} \frac{\partial \hat{I}}{\partial \boldsymbol{\theta}_j} - (I - \hat{I}) \frac{\partial^2 \hat{I}}{\partial \boldsymbol{\theta}_i \partial \boldsymbol{\theta}_j} \right]$$

**自动微分**：
使用反向模式自动微分（backpropagation）高效计算梯度：
1. 前向传播：计算渲染值
2. 反向传播：从损失函数反向计算梯度

对于路径追踪，需要处理不连续性：
- 使用边缘采样处理可见性变化
- 重参数化技巧处理采样相关的导数

### 13.1.4 BSSRDF估计

对于半透明材质，需要考虑次表面散射：

$$L_o(\mathbf{x}_o, \boldsymbol{\omega}_o) = \int_A \int_{\Omega} S(\mathbf{x}_o, \boldsymbol{\omega}_o, \mathbf{x}_i, \boldsymbol{\omega}_i) L_i(\mathbf{x}_i, \boldsymbol{\omega}_i) (\boldsymbol{\omega}_i \cdot \mathbf{n}_i) d\boldsymbol{\omega}_i dA$$

**完整BSSRDF分解**：
$$S(\mathbf{x}_o, \boldsymbol{\omega}_o, \mathbf{x}_i, \boldsymbol{\omega}_i) = S^{(1)}(\mathbf{x}_o, \boldsymbol{\omega}_o, \mathbf{x}_i, \boldsymbol{\omega}_i) + S^{d}(\mathbf{x}_o, \boldsymbol{\omega}_o, \mathbf{x}_i, \boldsymbol{\omega}_i)$$

其中：
- $S^{(1)}$: 单次散射项
- $S^{d}$: 多次散射（扩散）项

偶极子近似下的BSSRDF：
$$S(\mathbf{x}_o, \mathbf{x}_i) = \frac{1}{\pi} F_t(\eta) R_d(\|\mathbf{x}_o - \mathbf{x}_i\|) F_t(\eta)$$

其中 $R_d(r)$ 为扩散剖面，依赖于吸收系数 $\sigma_a$ 和散射系数 $\sigma_s$。

**扩散剖面**：
$$R_d(r) = \frac{\alpha'}{4\pi} \left[ \frac{z_r(1 + \sigma_{tr}d_r)e^{-\sigma_{tr}d_r}}{d_r^3} + \frac{z_v(1 + \sigma_{tr}d_v)e^{-\sigma_{tr}d_v}}{d_v^3} \right]$$

其中：
- $\sigma_{tr} = \sqrt{3\sigma_a(\sigma_a + \sigma_s')}$ 为有效传输系数
- $\sigma_s' = \sigma_s(1 - g)$ 为约化散射系数
- $\alpha' = \sigma_s'/\sigma_{tr}$ 为反照率
- $z_r = 1/\sigma_{tr}$，$z_v = z_r + 4AD$ 为偶极子深度
- $d_r = \sqrt{r^2 + z_r^2}$，$d_v = \sqrt{r^2 + z_v^2}$ 为偶极子距离
- $A = \frac{1 + F_{dr}}{1 - F_{dr}}$，$D = 1/(3\sigma_{tr})$ 为边界条件参数

**改进的扩散模型**：

1. **量化扩散**（Quantized Diffusion）：
   处理多层材质的离散化方法：
   $$R_d(r) = \sum_{n=0}^{\infty} R_n \exp(-\sigma_n r)$$
   
   其中 $R_n$ 和 $\sigma_n$ 通过求解扩散方程的特征值问题得到。

2. **光子束扩散**（Photon Beam Diffusion）：
   更准确地处理各向异性和几何边界：
   $$R_{PBD}(r, \theta) = \int_0^{\infty} J_0(sr) \tilde{R}(s) e^{-\mu_0 s \cos\theta} s ds$$
   
   其中 $J_0$ 为零阶贝塞尔函数，$\tilde{R}(s)$ 为频域反射率。

**参数化优化**：
我们优化材质参数 $\boldsymbol{\phi} = [\sigma_a, \sigma_s, g, \eta]$：

$$\min_{\boldsymbol{\phi}} \sum_{\mathbf{x}} \left\| I(\mathbf{x}) - \int_A S[\boldsymbol{\phi}](\mathbf{x}, \mathbf{x}') L(\mathbf{x}') dA' \right\|^2$$

**梯度计算**：
对于扩散参数的导数：
$$\frac{\partial R_d}{\partial \sigma_a} = R_d \left[ \frac{\partial \ln \alpha'}{\partial \sigma_a} + \sum_{i \in \{r,v\}} \frac{\partial \ln T_i}{\partial \sigma_a} \right]$$

其中传输项：
$$T_i = \frac{z_i(1 + \sigma_{tr}d_i)e^{-\sigma_{tr}d_i}}{d_i^3}$$

**采样策略**：
1. **重要性采样**：
   根据扩散剖面采样入射点：
   $$p(\mathbf{x}_i | \mathbf{x}_o) \propto R_d(\|\mathbf{x}_o - \mathbf{x}_i\|)$$
   
2. **分层采样**：
   将采样域分为近场（$r < 3l_t$）和远场（$r \geq 3l_t$）

**分层表示**：
对于多层材质，使用分层BSSRDF：
$$S_{total} = \sum_{k=1}^K w_k S_k(\sigma_{a,k}, \sigma_{s,k})$$

其中 $w_k$ 为各层权重，满足 $\sum_k w_k = 1$。

**快速近似**：
1. **可分离近似**：
   $$S(\mathbf{x}_o, \mathbf{x}_i) \approx S_r(\|\mathbf{x}_o - \mathbf{x}_i\|) S_{\omega}(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)$$
   
2. **高斯和近似**：
   $$R_d(r) \approx \sum_{i=1}^n w_i \frac{e^{-r^2/2v_i}}{2\pi v_i}$$
   
   其中 $w_i$ 和 $v_i$ 通过最小二乘拟合得到。

**非均匀介质**：
对于空间变化的散射参数：
$$\nabla \cdot (D(\mathbf{x})\nabla\phi(\mathbf{x})) - \sigma_a(\mathbf{x})\phi(\mathbf{x}) + Q(\mathbf{x}) = 0$$

其中 $\phi$ 为光通量密度，$Q$ 为源项。使用有限元方法（FEM）求解。

## 13.2 形状从明暗恢复

### 13.2.1 经典形状从明暗

对于朗伯表面在正交投影下，亮度方程为：

$$I(x,y) = \rho \mathbf{n}(x,y) \cdot \mathbf{l}$$

其中表面法线与深度函数 $z(x,y)$ 的关系：

$$\mathbf{n} = \frac{(-\partial z/\partial x, -\partial z/\partial y, 1)}{\sqrt{1 + (\partial z/\partial x)^2 + (\partial z/\partial y)^2}}$$

引入梯度符号 $\mathbf{p} = -\partial z/\partial x$，$\mathbf{q} = -\partial z/\partial y$，则：
$$\mathbf{n} = \frac{(p, q, 1)}{\sqrt{1 + p^2 + q^2}}$$

**反射图方程**：
当光源方向 $\mathbf{l} = (l_x, l_y, l_z)$ 已知时，亮度方程变为：
$$I(x,y) = \frac{\rho(p l_x + q l_y + l_z)}{\sqrt{1 + p^2 + q^2}}$$

这是一个关于 $(p, q)$ 的非线性偏微分方程，称为 **Eikonal方程**。

**特征条带方法**：
沿着特征曲线，方程可以简化为常微分方程。特征曲线满足：
$$\frac{dx}{dt} = \frac{\partial H}{\partial p}, \quad \frac{dy}{dt} = \frac{\partial H}{\partial q}$$

其中 $H(p, q) = \sqrt{1 + p^2 + q^2} I(x,y) - \rho(p l_x + q l_y + l_z)$ 为Hamiltonian。

### 13.2.2 变分公式

能量泛函包含数据项和正则化项：

$$E[z] = \int_{\Omega} (I(x,y) - \rho \mathbf{n}[z] \cdot \mathbf{l})^2 dxdy + \lambda \int_{\Omega} |\nabla z|^2 dxdy$$

**泛函变分推导**：
设 $z + \epsilon h$ 为深度的扰动，其中 $h$ 为测试函数，$\epsilon \to 0$。

首先计算法线的变分：
$$\mathbf{n}[z + \epsilon h] = \frac{(-\nabla(z + \epsilon h), 1)}{\sqrt{1 + |\nabla(z + \epsilon h)|^2}}$$

对 $\epsilon$ 求导并在 $\epsilon = 0$ 处计算：
$$\frac{d\mathbf{n}}{d\epsilon}\bigg|_{\epsilon=0} = \frac{(-\nabla h, 0)}{\sqrt{1 + |\nabla z|^2}} - \frac{(-\nabla z, 1)(\nabla z \cdot \nabla h)}{(1 + |\nabla z|^2)^{3/2}}$$

欧拉-拉格朗日方程：

$$\frac{\partial}{\partial z}\left(\frac{(I - \rho \mathbf{n} \cdot \mathbf{l})^2}{\sqrt{1 + |\nabla z|^2}}\right) - \lambda \nabla^2 z = 0$$

**完整形式的欧拉-拉格朗日方程**：
对数据项进行变分：
$$\frac{\delta E_{data}}{\delta z} = -2(I - \rho \mathbf{n} \cdot \mathbf{l}) \rho \frac{\partial (\mathbf{n} \cdot \mathbf{l})}{\partial z}$$

其中：
$$\frac{\partial (\mathbf{n} \cdot \mathbf{l})}{\partial z} = \frac{\mathbf{l} \cdot \nabla^2 z \cdot \mathbf{n} - (\mathbf{l} \cdot \mathbf{n})(\mathbf{n} \cdot \nabla^2 z \cdot \mathbf{n})}{(1 + |\nabla z|^2)^{3/2}}$$

**高阶正则化**：
考虑曲率正则化以获得更平滑的解：
$$E_{reg}[z] = \int_{\Omega} \left(\kappa_1^2 + \kappa_2^2\right) dxdy$$

其中主曲率：
$$\kappa_1, \kappa_2 = \frac{z_{xx} + z_{yy} \pm \sqrt{(z_{xx} - z_{yy})^2 + 4z_{xy}^2}}{2(1 + |\nabla z|^2)^{3/2}}$$

**数值解法**：
1. **固定点迭代**：
   $$z^{(k+1)} = z^{(k)} - \tau \left[ \frac{\delta E_{data}}{\delta z} + \lambda \nabla^2 z^{(k)} \right]$$
   
   其中步长 $\tau$ 通过线搜索确定。

2. **多重网格方法**：
   从粗网格到细网格逐步求解，避免局部最小值
   - 限制算子：$I_h^{2h}: \Omega_h \to \Omega_{2h}$
   - 延拓算子：$I_{2h}^h: \Omega_{2h} \to \Omega_h$
   - V-循环或W-循环迭代

3. **线性化方法**：
   在每次迭代中将非线性项线性化：
   $$I \approx \rho(\mathbf{n}^{(k)} \cdot \mathbf{l}) + \rho \frac{\partial(\mathbf{n} \cdot \mathbf{l})}{\partial z}\bigg|_{z^{(k)}} (z - z^{(k)})$$

4. **半隐式方案**：
   处理非线性和刚性问题：
   $$\frac{z^{(k+1)} - z^{(k)}}{\Delta t} = -\frac{\delta E}{\delta z}[z^{(k+1)}]_{linear} - \frac{\delta E}{\delta z}[z^{(k)}]_{nonlinear}$$

**边界条件处理**：
1. **Dirichlet条件**：$z|_{\partial\Omega} = z_0$
2. **Neumann条件**：$\frac{\partial z}{\partial n}|_{\partial\Omega} = 0$
3. **混合条件**：结合已知深度和法线信息

**稳定性分析**：
线性化算子的谱半径决定收敛性：
$$\rho(L) = \max_i |\lambda_i(L)| < 1$$

其中 $L$ 为迭代矩阵。CFL条件：
$$\tau < \frac{2}{\lambda_{max}(H)}$$

其中 $H$ 为Hessian矩阵。

### 13.2.3 多光源形状从明暗

给定 $m$ 个光照条件下的图像 $\{I_j\}_{j=1}^m$：

$$E[z] = \sum_{j=1}^m \int_{\Omega} (I_j - \rho \mathbf{n}[z] \cdot \mathbf{l}_j)^2 dxdy$$

最优性条件产生非线性PDE系统。

**约束优化公式**：
引入辅助变量 $\mathbf{N} = (n_x, n_y, n_z)$ 表示法线场：

$$\min_{z, \mathbf{N}} \sum_{j=1}^m \int_{\Omega} (I_j - \rho \mathbf{N} \cdot \mathbf{l}_j)^2 dxdy$$

s.t. $\mathbf{N} = \frac{(-z_x, -z_y, 1)}{\sqrt{1 + z_x^2 + z_y^2}}$, $\|\mathbf{N}\| = 1$

**交替方向乘子法**（ADMM）：
$$\mathcal{L} = \sum_j \|I_j - \rho \mathbf{N} \cdot \mathbf{l}_j\|^2 + \frac{\mu}{2}\|\mathbf{N} - \mathbf{n}[z] + \mathbf{u}\|^2$$

迭代步骤：
1. **N-子问题**：$\mathbf{N}^{(k+1)} = \arg\min_{\|\mathbf{N}\|=1} \mathcal{L}(\mathbf{N}, z^{(k)}, \mathbf{u}^{(k)})$
2. **z-子问题**：$z^{(k+1)} = \arg\min_z \|\mathbf{N}^{(k+1)} - \mathbf{n}[z] + \mathbf{u}^{(k)}\|^2$
3. **对偶更新**：$\mathbf{u}^{(k+1)} = \mathbf{u}^{(k)} + \mathbf{N}^{(k+1)} - \mathbf{n}[z^{(k+1)}]$

**鲁棒估计**：
使用Huber损失函数处理异常值：
$$\rho_{Huber}(r) = \begin{cases} 
\frac{1}{2}r^2 & |r| \leq \delta \\
\delta(|r| - \frac{\delta}{2}) & |r| > \delta
\end{cases}$$

### 13.2.4 光度立体

当 $m \geq 3$ 时，可以线性求解法线：

$$\begin{bmatrix} I_1 \\ I_2 \\ \vdots \\ I_m \end{bmatrix} = \rho \begin{bmatrix} \mathbf{l}_1^T \\ \mathbf{l}_2^T \\ \vdots \\ \mathbf{l}_m^T \end{bmatrix} \mathbf{n}$$

通过最小二乘：$\mathbf{n} = \frac{1}{\rho}(\mathbf{L}^T\mathbf{L})^{-1}\mathbf{L}^T\mathbf{I}$

**鲁棒光度立体**：
处理阴影和高光的鲁棒估计：
1. **RANSAC方法**：
   随机选择3个光源，估计法线，计算内点
   
2. **稀疏回归**：
   $$\min_{\mathbf{n}, \mathbf{e}} \|\mathbf{I} - \rho\mathbf{L}\mathbf{n} - \mathbf{e}\|^2 + \lambda\|\mathbf{e}\|_1$$
   
   其中 $\mathbf{e}$ 为稀疏误差项（阴影/高光）。

3. **矩阵秩最小化**：
   $$\min_{\mathbf{N}, \mathbf{E}} rank(\mathbf{N}) + \lambda\|\mathbf{E}\|_0$$
   s.t. $\mathbf{I} = \mathbf{N} + \mathbf{E}$

**法线积分**：
从法线场恢复深度需要满足可积性条件：
$$\frac{\partial n_x/n_z}{\partial y} = \frac{\partial n_y/n_z}{\partial x}$$

使用泊松方程求解：
$$\nabla^2 z = \nabla \cdot \left( \frac{n_x}{n_z}, \frac{n_y}{n_z} \right)$$

**快速傅里叶变换求解**：
在频域中：
$$\hat{z}(u,v) = \frac{i2\pi(u\hat{n}_x + v\hat{n}_y)/\hat{n}_z}{-4\pi^2(u^2 + v^2)}$$

其中 $\hat{·}$ 表示傅里叶变换。

**最小二乘积分**：
考虑噪声的影响，最小化：
$$E[z] = \int_{\Omega} \left[ \left(\frac{\partial z}{\partial x} + \frac{n_x}{n_z}\right)^2 + \left(\frac{\partial z}{\partial y} + \frac{n_y}{n_z}\right)^2 \right] dxdy$$

**非朗伯表面的光度立体**：
对于一般的BRDF $\rho(\mathbf{n}, \mathbf{l}, \mathbf{v})$：
$$I_j = \rho(\mathbf{n}, \mathbf{l}_j, \mathbf{v})$$

**球谐函数方法**：
使用球谐函数展开：
$$\rho(\mathbf{n}, \mathbf{l}, \mathbf{v}) \approx \sum_{l=0}^{L} \sum_{m=-l}^{l} a_{lm}(\mathbf{v}) Y_{lm}(\mathbf{n}) Y_{lm}(\mathbf{l})$$

这将非线性问题转化为线性系统。

对于各向同性BRDF，使用Zernike多项式：
$$\rho(\theta_h) = \sum_{n=0}^{N} c_n P_n(\cos\theta_h)$$

其中 $\theta_h$ 为半角，$P_n$ 为Legendre多项式。

**双向纹理函数**（BTF）：
考虑空间变化的BRDF：
$$I_j(\mathbf{x}) = \rho(\mathbf{x}, \mathbf{n}(\mathbf{x}), \mathbf{l}_j, \mathbf{v})$$

使用张量分解：
$$\rho(\mathbf{x}, \boldsymbol{\omega}_i, \boldsymbol{\omega}_o) \approx \sum_{k=1}^{K} a_k(\mathbf{x}) b_k(\boldsymbol{\omega}_i) c_k(\boldsymbol{\omega}_o)$$

**校准光度立体**：
使用已知形状的参考物体校准光源方向：
$$\mathbf{l}_j = \arg\min_{\|\mathbf{l}\|=1} \sum_{\mathbf{x} \in \mathcal{R}} (I_j(\mathbf{x}) - \rho_{ref} \mathbf{n}_{ref}(\mathbf{x}) \cdot \mathbf{l})^2$$

**自校准方法**：
利用可积性约束同时估计光源和法线：
$$\min_{\mathbf{L}, \mathbf{N}} \|\mathbf{I} - \mathbf{L}\mathbf{N}\|_F^2 + \lambda \int_{\Omega} \left\|\frac{\partial \mathbf{p}}{\partial y} - \frac{\partial \mathbf{q}}{\partial x}\right\|^2 dxdy$$

其中 $\mathbf{p} = n_x/n_z$，$\mathbf{q} = n_y/n_z$。

**深度学习方法**：
使用CNN直接从图像预测法线：
$$\mathbf{n} = f_{\theta}(\{I_j\}_{j=1}^m)$$

其中 $f_{\theta}$ 为神经网络，$\theta$ 为学习参数。

## 13.3 多视图立体重建中的优化

### 13.3.1 体积雕刻能量

定义占用场 $\phi: \mathbb{R}^3 \rightarrow \{0,1\}$，光度一致性能量：

$$E[\phi] = \sum_{i,j} \int_{\partial\Omega[\phi]} \|I_i(\pi_i(\mathbf{x})) - I_j(\pi_j(\mathbf{x}))\|^2 dS$$

其中 $\pi_i$ 为相机 $i$ 的投影函数。

**投影函数**：
对于透视相机：
$$\pi_i(\mathbf{x}) = \mathbf{K}_i \mathbf{R}_i (\mathbf{x} - \mathbf{c}_i)$$

其中：
- $\mathbf{K}_i$ 为内参矩阵
- $\mathbf{R}_i$ 为旋转矩阵
- $\mathbf{c}_i$ 为相机中心

**可见性约束**：
仅在点 $\mathbf{x}$ 对两个相机都可见时计算光度一致性：
$$E[\phi] = \sum_{i,j} \int_{\partial\Omega[\phi]} V_i(\mathbf{x}) V_j(\mathbf{x}) \|I_i - I_j\|^2 dS$$

其中 $V_i(\mathbf{x}) = 1$ 当 $\mathbf{x}$ 对相机 $i$ 可见。

**空间雕刻算法**：
1. 初始化为包含所有点的体积
2. 迭代删除不一致的体素：
   $$\phi^{(k+1)}(\mathbf{x}) = \begin{cases}
   0 & \text{if } \sum_{i,j} V_i V_j \|I_i - I_j\|^2 > \tau \\
   \phi^{(k)}(\mathbf{x}) & \text{otherwise}
   \end{cases}$$

### 13.3.2 变分水平集方法

使用隐式表面 $\{\mathbf{x}: \phi(\mathbf{x}) = 0\}$：

$$E[\phi] = \sum_{i,j} \int_{\mathbb{R}^3} \|I_i - I_j\|^2 \delta(\phi) |\nabla\phi| d\mathbf{x}$$

演化方程：
$$\frac{\partial \phi}{\partial t} = \delta(\phi) \left[ \sum_{i,j} (I_i - I_j)(\nabla I_i - \nabla I_j) \cdot \frac{\nabla\phi}{|\nabla\phi|} + \lambda \kappa \right]$$

**完整的能量泛函**：
$$E[\phi] = E_{photo}[\phi] + \lambda_1 E_{smooth}[\phi] + \lambda_2 E_{silhouette}[\phi]$$

其中：
- **光度一致性**：$E_{photo} = \sum_{i,j} \int \rho(I_i - I_j) \delta(\phi) |\nabla\phi| d\mathbf{x}$
- **平滑性**：$E_{smooth} = \int \delta(\phi) |\nabla\phi| d\mathbf{x}$ （最小表面积）
- **轮廓约束**：$E_{silhouette} = \sum_i \int (S_i - \hat{S}_i[\phi])^2 d\mathbf{x}$

**数值实现**：
1. **窄带水平集**：仅在 $|\phi| < \epsilon$ 的窄带内更新
2. **重初始化**：周期性将 $\phi$ 重初始化为符号距离函数
3. **CFL条件**：$\Delta t < \frac{\Delta x}{\max|F|}$，其中 $F$ 为速度场

**拓扑变化**：
水平集方法自然处理拓扑变化（合并、分裂）：
$$\phi_{merge} = \min(\phi_1, \phi_2)$$
$$\phi_{split} = \max(\phi_1, -\phi_2)$$

### 13.3.3 PatchMatch立体

离散优化公式，为每个像素分配深度 $d_p$：

$$E(\{d_p\}) = \sum_p C(p, d_p) + \sum_{(p,q) \in \mathcal{N}} V(d_p, d_q)$$

其中：
- $C(p,d)$: 匹配代价
- $V(d_p,d_q)$: 平滑项

**匹配代价函数**：
基于局部窗口的归一化互相关（NCC）：
$$C_{NCC}(p, d) = 1 - \frac{\sum_{q \in W(p)} (I_{ref}(q) - \bar{I}_{ref})(I_{src}(\pi(q,d)) - \bar{I}_{src})}{\sqrt{\sum_q (I_{ref}(q) - \bar{I}_{ref})^2} \sqrt{\sum_q (I_{src} - \bar{I}_{src})^2}}$$

其中 $W(p)$ 为以 $p$ 为中心的窗口，$\pi(q,d)$ 为深度 $d$ 处的重投影。

**PatchMatch传播**：
核心思想是利用空间连贯性：
1. **随机初始化**：$d_p \sim U[d_{min}, d_{max}]$
2. **传播**：从邻域传播好的深度值
   $$d_p^{(k+1)} = \arg\min_{d \in \{d_p^{(k)}, d_{p-1}^{(k)}, d_{p-w}^{(k)}\}} C(p, d)$$
3. **随机搜索**：在当前值附近随机采样
   $$d_{test} = d_p + \Delta d \cdot \alpha^i, \quad \Delta d \sim U[-1, 1]$$
   其中 $\alpha < 1$ 为搜索半径缩减系数。

**平面拟合扩展**：
为每个像素估计平面参数 $(a_p, b_p, c_p)$：
$$d(x,y) = a_p x + b_p y + c_p$$

这提供了亚像素精度和更好的传播。

### 13.3.4 平面扫描体积

构建代价体积 $\mathcal{C}(x,y,d)$：

$$\mathcal{C}(x,y,d) = \frac{1}{|\mathcal{V}|} \sum_{i \in \mathcal{V}} \|I_{ref}(x,y) - I_i(\mathbf{H}_i(x,y,d))\|$$

其中 $\mathbf{H}_i$ 为深度 $d$ 处的单应变换。

**单应变换推导**：
对于参考视图中的点 $(x,y)$ 在深度 $d$ 处：
1. 3D点：$\mathbf{X} = d \mathbf{K}_{ref}^{-1}[x, y, 1]^T$
2. 投影到源视图：$[x', y', 1]^T \sim \mathbf{K}_i \mathbf{R}_i (\mathbf{X} - \mathbf{c}_i)$

单应矩阵：
$$\mathbf{H}_i(d) = \mathbf{K}_i (\mathbf{R}_i - \frac{\mathbf{t}_i \mathbf{n}^T}{d}) \mathbf{K}_{ref}^{-1}$$

其中 $\mathbf{n}$ 为参考视图中的平面法线（通常为 $[0,0,1]^T$）。

**代价体积正则化**：
使用加权中值滤波或双边滤波：
$$\tilde{\mathcal{C}}(x,y,d) = \frac{\sum_{(x',y',d') \in \mathcal{N}} w(x',y',d') \mathcal{C}(x',y',d')}{\sum_{(x',y',d') \in \mathcal{N}} w(x',y',d')}$$

权重函数：
$$w(x',y',d') = \exp\left(-\frac{\|(x-x',y-y')\|^2}{2\sigma_s^2} - \frac{(d-d')^2}{2\sigma_d^2}\right)$$

**深度图提取**：
1. **Winner-takes-all**：$d^*(x,y) = \arg\min_d \mathcal{C}(x,y,d)$
2. **亚像素精度**：使用抛物线拟合
   $$d_{sub} = d^* - \frac{\mathcal{C}(d^*+1) - \mathcal{C}(d^*-1)}{2(\mathcal{C}(d^*+1) + \mathcal{C}(d^*-1) - 2\mathcal{C}(d^*))}$$

## 13.4 联合材质-几何优化

### 13.4.1 耦合优化框架

同时优化几何 $\mathcal{G}$ 和材质 $\mathcal{M}$：

$$E[\mathcal{G}, \mathcal{M}] = E_{photo}[\mathcal{G}, \mathcal{M}] + \lambda_g E_{geom}[\mathcal{G}] + \lambda_m E_{mat}[\mathcal{M}]$$

**详细的能量项**：

1. **光度一致性**：
   $$E_{photo} = \sum_{i} \int_{\Omega} V_i(\mathbf{x}) \|I_i(\pi_i(\mathbf{x})) - \mathcal{R}[\mathcal{G}, \mathcal{M}](\mathbf{x}, \boldsymbol{\omega}_i)\|^2 d\mathbf{x}$$
   
   其中 $\mathcal{R}$ 为渲染算子，考虑几何和材质。

2. **几何正则化**：
   $$E_{geom} = \alpha_1 \int_{\mathcal{S}} H^2 dS + \alpha_2 \int_{\mathcal{S}} K^2 dS$$
   
   其中 $H$ 为平均曲率，$K$ 为高斯曲率。

3. **材质正则化**：
   $$E_{mat} = \beta_1 \int_{\mathcal{S}} \|\nabla_{\mathcal{S}} \mathcal{M}\|^2 dS + \beta_2 D_{KL}(p(\mathcal{M}) \| p_{prior}(\mathcal{M}))$$
   
   其中 $\nabla_{\mathcal{S}}$ 为表面梯度，$D_{KL}$ 为KL散度。

### 13.4.2 交替优化

**几何步骤**（固定材质）：
$$\mathcal{G}^{(k+1)} = \arg\min_{\mathcal{G}} E[\mathcal{G}, \mathcal{M}^{(k)}]$$

**材质步骤**（固定几何）：
$$\mathcal{M}^{(k+1)} = \arg\min_{\mathcal{M}} E[\mathcal{G}^{(k+1)}, \mathcal{M}]$$

**收敛性分析**：
定义增广拉格朗日函数：
$$\mathcal{L}(\mathcal{G}, \mathcal{M}, \mathbf{z}) = E[\mathcal{G}, \mathcal{M}] + \frac{\rho}{2}\|\mathcal{F}(\mathcal{G}) - \mathcal{M} + \mathbf{z}\|^2$$

其中 $\mathcal{F}$ 为从几何到材质的映射（如法线到反照率）。

**收敛条件**：
如果 $E$ 是下有界的，且每个子问题都有唯一解，则：
$$E[\mathcal{G}^{(k+1)}, \mathcal{M}^{(k+1)}] \leq E[\mathcal{G}^{(k)}, \mathcal{M}^{(k)}]$$

且序列收敛到驻点。

**加速策略**：
1. **动量方法**：
   $$\mathcal{G}^{(k+1)} = \mathcal{G}^{(k)} - \alpha \nabla_{\mathcal{G}} E + \beta(\mathcal{G}^{(k)} - \mathcal{G}^{(k-1)})$$

2. **多尺度优化**：
   从粗糙到精细的几何表示

3. **预条件**：
   使用Hessian的对角近似

### 13.4.3 联合梯度下降

计算完整梯度：
$$\begin{bmatrix} \mathcal{G}^{(k+1)} \\ \mathcal{M}^{(k+1)} \end{bmatrix} = \begin{bmatrix} \mathcal{G}^{(k)} \\ \mathcal{M}^{(k)} \end{bmatrix} - \alpha \begin{bmatrix} \nabla_{\mathcal{G}} E \\ \nabla_{\mathcal{M}} E \end{bmatrix}$$

需要考虑Hessian的条件数。

**Hessian矩阵结构**：
$$\mathbf{H} = \begin{bmatrix} 
\nabla_{\mathcal{G}\mathcal{G}}^2 E & \nabla_{\mathcal{G}\mathcal{M}}^2 E \\
\nabla_{\mathcal{M}\mathcal{G}}^2 E & \nabla_{\mathcal{M}\mathcal{M}}^2 E
\end{bmatrix}$$

其中交叉项 $\nabla_{\mathcal{G}\mathcal{M}}^2 E$ 反映了几何-材质耦合。

**预条件器设计**：
使用块对角预条件器：
$$\mathbf{P} = \begin{bmatrix} 
(\nabla_{\mathcal{G}\mathcal{G}}^2 E)^{-1} & 0 \\
0 & (\nabla_{\mathcal{M}\mathcal{M}}^2 E)^{-1}
\end{bmatrix}$$

**自适应步长**：
使用Armijo-Goldstein线搜索：
$$\alpha^* = \arg\max_{\alpha} \{\alpha : E(\mathbf{x} - \alpha\nabla E) \leq E(\mathbf{x}) - c\alpha\|\nabla E\|^2\}$$

其中 $c \in (0, 0.5)$ 为常数。

**二阶方法**：
牛顿法更新：
$$\begin{bmatrix} \Delta\mathcal{G} \\ \Delta\mathcal{M} \end{bmatrix} = -\mathbf{H}^{-1} \begin{bmatrix} \nabla_{\mathcal{G}} E \\ \nabla_{\mathcal{M}} E \end{bmatrix}$$

使用拟Newton方法（L-BFGS）避免存储完整Hessian。

### 13.4.4 分解歧义性

材质-几何-光照分解存在固有歧义：
- 尺度歧义：$(\alpha\mathcal{M}, \mathcal{L}/\alpha)$ 产生相同图像
- 低频歧义：光照与反照率的低频分量难以区分

## 13.5 物理约束与先验

### 13.5.1 材质物理约束

**能量守恒**：
$$\int_{\Omega} \rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o \leq 1$$

**Helmholtz互易性**：
$$\rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \rho(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)$$

**单调性约束**（对于粗糙度）：
$$\alpha_1 < \alpha_2 \Rightarrow D_{\alpha_1}(\mathbf{h}) > D_{\alpha_2}(\mathbf{h}), \forall \mathbf{h} \neq \mathbf{n}$$

### 13.5.2 几何先验

**最小表面积**：
$$E_{area}[\mathcal{S}] = \int_{\mathcal{S}} dS$$

**平均曲率流**：
$$E_{smooth}[\mathcal{S}] = \int_{\mathcal{S}} H^2 dS$$

其中 $H = (\kappa_1 + \kappa_2)/2$ 为平均曲率。

**体积保持**：
$$\int_{\Omega} \phi d\mathbf{x} = V_0$$

### 13.5.3 统计先验

**材质分布先验**（对数正态）：
$$p(\alpha) = \frac{1}{\alpha\sigma\sqrt{2\pi}} \exp\left(-\frac{(\ln\alpha - \mu)^2}{2\sigma^2}\right)$$

**空间相关性**（MRF）：
$$p(\mathcal{M}) \propto \exp\left(-\sum_{(i,j) \in \mathcal{N}} \psi(\mathcal{M}_i, \mathcal{M}_j)\right)$$

### 13.5.4 学习先验

使用神经网络编码先验：
$$E_{prior}[\mathcal{M}] = \|\mathcal{M} - G_\theta(z)\|^2$$

其中 $G_\theta$ 为预训练的材质生成器。

## 本章小结

本章介绍了材质与几何重建的核心技术：

**关键概念**：
- BRDF/BSSRDF参数估计的优化框架
- 形状从明暗的变分公式和数值解法
- 多视图立体的连续和离散优化方法
- 联合材质-几何优化的耦合问题
- 物理约束和统计先验的作用

**核心方程**：
1. **逆向渲染**：$\min_{\rho} \|I - \mathcal{R}[\rho]\|^2$
2. **形状从明暗**：$(I - \rho\mathbf{n}\cdot\mathbf{l})^2 + \lambda|\nabla z|^2$
3. **多视图立体**：$\sum_{i,j}\|I_i(\pi_i(\mathbf{x})) - I_j(\pi_j(\mathbf{x}))\|^2$
4. **联合优化**：$E_{photo} + \lambda_g E_{geom} + \lambda_m E_{mat}$

## 常见陷阱与错误

1. **局部最小值**：非凸优化容易陷入局部最优，需要良好初始化
2. **尺度歧义**：忽略归一化导致材质-光照分解不稳定
3. **过拟合**：参数过多时容易过拟合观察数据
4. **数值稳定性**：法线计算中分母接近零需要正则化
5. **边界处理**：形状边界的可见性变化需要特殊处理

## 最佳实践检查清单

- [ ] 选择合适的BRDF参数化以平衡表达能力和稳定性
- [ ] 使用多尺度优化策略避免局部最小值
- [ ] 加入物理约束确保结果合理性
- [ ] 考虑光照-材质-几何的耦合关系
- [ ] 验证优化算法的收敛性
- [ ] 评估重建结果的不确定性

## 练习题

### 练习 13.1：微表面BRDF的能量守恒
证明GGX微表面模型满足能量守恒约束。具体地，证明：
$$\int_{\Omega} \rho_{GGX}(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o \leq 1$$

**提示**：利用微表面理论中的遮蔽-阴影函数 $G$ 的性质。

<details>
<summary>解答</summary>

对于微表面BRDF：
$$\rho = \frac{F(\boldsymbol{\omega}_i, \mathbf{h}) G(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) D(\mathbf{h})}{4(\boldsymbol{\omega}_i \cdot \mathbf{n})(\boldsymbol{\omega}_o \cdot \mathbf{n})}$$

积分可以变换到半角向量空间：
$$\int_{\Omega} \rho (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o = \int_{\Omega_h} F G D \frac{(\mathbf{h} \cdot \boldsymbol{\omega}_i)}{(\boldsymbol{\omega}_i \cdot \mathbf{n})} d\mathbf{h}$$

利用Smith遮蔽函数的性质：
$$\int_{\Omega} G(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) \frac{(\mathbf{h} \cdot \boldsymbol{\omega}_i)}{(\boldsymbol{\omega}_i \cdot \mathbf{n})} D(\mathbf{h}) d\mathbf{h} = 1$$

由于 $F \leq 1$，因此：
$$\int_{\Omega} \rho (\boldsymbol{\omega}_o \cdot \mathbf{n}) d\boldsymbol{\omega}_o \leq 1$$
</details>

### 练习 13.2：形状从明暗的唯一性
考虑朗伯表面的形状从明暗问题。证明在以下条件下解不唯一：
1. 单一光源
2. 无边界条件

构造两个不同的表面产生相同的图像。

**提示**：考虑凸/凹歧义。

<details>
<summary>解答</summary>

给定亮度方程：$I(x,y) = \mathbf{n}(x,y) \cdot \mathbf{l}$

设光源方向 $\mathbf{l} = (0, 0, 1)$（顶光照明）。则：
$$I(x,y) = \frac{1}{\sqrt{1 + p^2 + q^2}}$$

其中 $p = \partial z/\partial x$, $q = \partial z/\partial y$。

两个解：
1. 凸解：$z_1(x,y) = \sqrt{R^2 - x^2 - y^2}$（上半球）
2. 凹解：$z_2(x,y) = -\sqrt{R^2 - x^2 - y^2}$（下半球）

两者产生相同的图像，因为法线的 $z$ 分量相同：
$$n_z = \frac{1}{\sqrt{1 + (x/z)^2 + (y/z)^2}} = \frac{|z|}{\sqrt{x^2 + y^2 + z^2}} = \frac{|z|}{R}$$
</details>

### 练习 13.3：光度立体的最优光源配置
给定 $n = 3$ 个光源，求使光度立体重建最稳定的光源方向配置。稳定性用条件数 $\kappa(\mathbf{L}^T\mathbf{L})$ 衡量。

**提示**：最小化条件数等价于最大化最小特征值。

<details>
<summary>解答</summary>

光源矩阵：
$$\mathbf{L} = \begin{bmatrix} \mathbf{l}_1^T \\ \mathbf{l}_2^T \\ \mathbf{l}_3^T \end{bmatrix}$$

条件数：$\kappa = \lambda_{max}/\lambda_{min}$

对于 $n = 3$，最优配置使 $\mathbf{L}^T\mathbf{L} = I$，即光源正交：
$$\mathbf{l}_i \cdot \mathbf{l}_j = \delta_{ij}$$

一个最优配置：
$$\mathbf{l}_1 = (1, 0, 0), \quad \mathbf{l}_2 = (0, 1, 0), \quad \mathbf{l}_3 = (0, 0, 1)$$

此时 $\kappa = 1$（最小可能值）。

实际中考虑阴影，通常选择：
$$\mathbf{l}_i = (\sin\theta\cos(2\pi i/3), \sin\theta\sin(2\pi i/3), \cos\theta)$$
其中 $\theta \approx 45°$。
</details>

### 练习 13.4：BSSRDF的扩散近似误差
评估偶极子近似对薄材质的误差。设材质厚度为 $d$，推导当 $d \ll l_t$（输运平均自由程）时的相对误差。

**提示**：比较精确解（平板几何）与偶极子近似。

<details>
<summary>解答</summary>

平板几何的精确反射率：
$$R_{exact} = \frac{1 - \exp(-2\tau)}{1 + \exp(-2\tau)}$$
其中 $\tau = d/l_t$ 为光学厚度。

偶极子近似：
$$R_{dipole} \approx \frac{A}{1 + A}$$
其中 $A = (1 + F_{dr})/(1 - F_{dr})$，$F_{dr}$ 为内部反射系数。

泰勒展开对小 $\tau$：
$$R_{exact} \approx \tau - \frac{2\tau^3}{3} + O(\tau^5)$$
$$R_{dipole} \approx \tau - \tau^2 + O(\tau^3)$$

相对误差：
$$\epsilon = \frac{|R_{dipole} - R_{exact}|}{R_{exact}} \approx \tau = \frac{d}{l_t}$$

因此误差与厚度成正比，薄材质时偶极子近似失效。
</details>

### 练习 13.5：联合优化的收敛性分析
考虑简化的联合材质-几何优化：
$$E(\alpha, z) = (I - \alpha h(z))^2 + \lambda z^2$$
其中 $h(z)$ 是 $z$ 的非线性函数。分析交替优化的收敛条件。

**提示**：构造Lyapunov函数。

<details>
<summary>解答</summary>

交替优化步骤：
1. 固定 $z$：$\alpha^{(k+1)} = I h(z^{(k)})/h^2(z^{(k)})$
2. 固定 $\alpha$：$z^{(k+1)} = \arg\min_z (I - \alpha^{(k+1)} h(z))^2 + \lambda z^2$

定义增广能量：
$$\tilde{E}(\alpha, z, \beta) = (I - \alpha h(z))^2 + \lambda z^2 + \mu(\alpha - \beta)^2$$

可以证明：
$$E(\alpha^{(k+1)}, z^{(k+1)}) \leq E(\alpha^{(k)}, z^{(k)})$$

收敛条件：
1. $h(z)$ Lipschitz连续：$|h(z_1) - h(z_2)| \leq L|z_1 - z_2|$
2. $\lambda > L^2 I^2/(4h_{min}^2)$

其中 $h_{min} = \min_z |h(z)|$。
</details>

### 练习 13.6：多视图立体的基线选择
推导多视图立体中深度估计误差与相机基线的关系。设深度 $z$，基线 $b$，视差误差 $\Delta d$。

**提示**：使用三角测量的误差传播。

<details>
<summary>解答</summary>

三角测量关系：
$$z = \frac{fb}{d}$$

其中 $f$ 为焦距，$d$ 为视差。

深度误差：
$$\Delta z = \frac{\partial z}{\partial d} \Delta d = -\frac{fb}{d^2} \Delta d = -\frac{z^2}{fb} \Delta d$$

相对误差：
$$\frac{\Delta z}{z} = \frac{z}{fb} \Delta d$$

结论：
1. 深度误差与深度平方成正比
2. 增大基线 $b$ 减小误差
3. 但基线过大导致遮挡和匹配困难

最优基线选择需平衡精度和匹配可靠性：
$$b_{opt} \approx 0.1 \cdot z_{avg}$$
</details>

### 练习 13.7：材质先验的信息论分析
使用最大熵原理推导粗糙度参数 $\alpha \in [0,1]$ 的先验分布。已知约束：$E[\alpha] = \mu$，$E[\log\alpha] = \nu$。

**提示**：拉格朗日乘数法求解约束最大熵问题。

<details>
<summary>解答</summary>

最大熵问题：
$$\max_{p(\alpha)} H[p] = -\int_0^1 p(\alpha) \log p(\alpha) d\alpha$$

约束条件：
1. $\int_0^1 p(\alpha) d\alpha = 1$
2. $\int_0^1 \alpha p(\alpha) d\alpha = \mu$
3. $\int_0^1 \log\alpha \cdot p(\alpha) d\alpha = \nu$

拉格朗日函数：
$$\mathcal{L} = H[p] + \lambda_0(1 - \int p) + \lambda_1(\mu - \int \alpha p) + \lambda_2(\nu - \int \log\alpha \cdot p)$$

变分得：
$$p(\alpha) = \frac{1}{Z} \exp(-\lambda_1 \alpha - \lambda_2 \log\alpha) = \frac{1}{Z} \alpha^{-\lambda_2} \exp(-\lambda_1 \alpha)$$

这是广义逆高斯分布的特例。当 $\lambda_2 = 1$ 时，退化为指数分布。
</details>

### 练习 13.8：物理约束的投影算法
设计投影算子 $\Pi$ 将任意函数投影到满足互易性的BRDF空间。定义合适的度量并证明投影的最优性。

**提示**：考虑 $L^2$ 度量下的投影。

<details>
<summary>解答</summary>

互易性约束：$\rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \rho(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)$

对任意函数 $f$，定义投影：
$$\Pi[f] = \arg\min_{\rho \in \mathcal{R}} \int_{\Omega^2} |f(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) - \rho(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o)|^2 d\boldsymbol{\omega}_i d\boldsymbol{\omega}_o$$

其中 $\mathcal{R}$ 为互易BRDF集合。

解为对称化：
$$\Pi[f](\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) = \frac{1}{2}[f(\boldsymbol{\omega}_i, \boldsymbol{\omega}_o) + f(\boldsymbol{\omega}_o, \boldsymbol{\omega}_i)]$$

证明最优性：设 $\rho \in \mathcal{R}$，则：
$$\|f - \rho\|^2 = \|\frac{f + f^T}{2} - \rho\|^2 + \|\frac{f - f^T}{2}\|^2 \geq \|f - \Pi[f]\|^2$$

因为 $\rho$ 与 $(f - f^T)/2$ 正交。
</details>
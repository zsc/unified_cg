# 第4章：基于图像的渲染

基于图像的渲染（Image-Based Rendering, IBR）代表了计算机图形学中的一个范式转变：从基于几何的渲染转向基于采样的渲染。本章将探讨如何使用预先捕获的图像来合成新视角，并将这些技术统一在体积渲染框架下。我们将深入研究光场的数学表示、视图合成的几何约束，以及采样理论在IBR中的应用。

IBR的核心思想是将渲染问题转化为信号重建问题。与传统的基于几何的方法不同，IBR直接对光的分布进行采样和插值，避免了复杂的光照计算和材质建模。这种方法特别适合于捕获和重现真实世界的复杂场景，包括那些难以用传统方法建模的现象，如次表面散射、复杂反射和半透明材质。

## 学习目标

完成本章后，您将能够：
1. 将光场渲染表述为4D函数的插值问题
2. 推导视图合成中的极线约束及其在体积渲染中的应用
3. 分析光图与表面光场的数学关系
4. 应用傅里叶切片定理优化光场重建
5. 计算IBR系统的最小采样要求和计算复杂度
6. 将IBR技术统一在体积渲染方程框架下

## 4.1 光场渲染作为4D插值

### 4.1.1 光场参数化

回顾第2章的全光函数（Plenoptic Function），我们将其简化为4D光场：

$$L(s,t,u,v) = \int_{\lambda} L_{\lambda}(s,t,u,v) d\lambda$$

其中 $(s,t)$ 和 $(u,v)$ 分别表示两个平行平面上的坐标。这种双平面参数化（Two-Plane Parameterization）避免了球面参数化的奇点问题。

**参数化的几何意义**：考虑两个平行平面 $\Pi_1$ 和 $\Pi_2$，相距 $d$。任何不平行于这两个平面的光线都可以唯一地由其与两平面的交点确定：
- $(s,t) \in \Pi_1$：光线与第一个平面的交点
- $(u,v) \in \Pi_2$：光线与第二个平面的交点

光线的参数方程为：
$$\mathbf{r}(\tau) = (1-\tau)\begin{pmatrix}s\\t\\0\end{pmatrix} + \tau\begin{pmatrix}u\\v\\d\end{pmatrix}, \quad \tau \in [0,1]$$

**从光线参数到方向向量的映射**：
光线方向向量为：
$$\mathbf{d} = \begin{pmatrix}u-s\\v-t\\d\end{pmatrix}$$

归一化方向：
$$\omega = \frac{\mathbf{d}}{\|\mathbf{d}\|} = \frac{1}{\sqrt{(u-s)^2 + (v-t)^2 + d^2}}\begin{pmatrix}u-s\\v-t\\d\end{pmatrix}$$

**参数化的选择准则**：
1. **完备性**：能表示所有相关光线
2. **唯一性**：每条光线有唯一表示
3. **连续性**：相邻光线在参数空间中相邻
4. **计算效率**：易于采样和插值
5. **均匀性**：参数空间的度量与物理空间的度量相关

**参数化的雅可比行列式**：
从双平面参数到光线空间 $(\mathbf{x}, \omega)$ 的映射的雅可比行列式为：
$$J = \frac{\partial(\mathbf{x}, \omega)}{\partial(s,t,u,v)} = \frac{d}{(\omega \cdot \mathbf{n})^3}$$

其中 $\mathbf{n}$ 是平面法向量。这个雅可比在积分变换中至关重要。

除了双平面参数化，常见的还有：
- **球面参数化**：$L(\mathbf{x}, \theta, \phi)$，适合全景捕获
- **平面+方向参数化**：$L(s,t,\theta,\phi)$，适合相机阵列
- **表面参数化**：$L(u,v,\theta,\phi)$，适合已知几何
- **光线空间参数化**：使用Plücker坐标 $(\mathbf{l}, \mathbf{m})$，其中 $\mathbf{l}$ 是方向，$\mathbf{m} = \mathbf{p} \times \mathbf{l}$ 是矩

### 4.1.2 离散采样与插值

给定一组采样光线 $\{L_{ijkl}\}$，新视角的光线通过4D插值获得：

$$\hat{L}(s,t,u,v) = \sum_{i,j,k,l} L_{ijkl} \cdot K(s-s_i, t-t_j, u-u_k, v-v_l)$$

其中 $K$ 是4D插值核。常见选择包括：
- **最近邻**：$K = \delta$，计算最快但质量最差
- **三线性**：$K = \prod_{d} (1-|x_d|)_+$，平衡了质量和效率
- **Lanczos**：$K = \prod_{d} \text{sinc}(x_d)\text{sinc}(x_d/a)$，高质量但计算密集
- **高斯核**：$K = \prod_{d} \exp(-x_d^2/2\sigma^2)$，平滑但有模糊
- **Mitchell-Netravali**：$K = \prod_{d} B(x_d)$，其中 $B$ 是双三次B样条

**插值核的数学性质**：
理想的插值核应满足以下条件：
1. **分离条件**：$K(0,0,0,0) = 1$，$K(i,j,k,l) = 0$ 对所有非零整数
2. **归一化**：$\sum_{ijkl} K(s-s_i, t-t_j, u-u_k, v-v_l) = 1$
3. **紧支撑**：$K(x) = 0$ 当 $|x| > r$ 对某个有限 $r$
4. **平滑性**：$K \in C^k$ 对某个 $k \geq 0$
5. **正定性**：$K(x) \geq 0$ 对所有 $x$（某些核如Lanczos可能违反）

**广义插值理论**：
将插值视为在再生核Hilbert空间（RKHS）中的投影：
$$\hat{L} = \arg\min_{f \in \mathcal{H}} \sum_{ijkl} (f(s_i,t_j,u_k,v_l) - L_{ijkl})^2 + \lambda\|f\|_{\mathcal{H}}^2$$

解具有形式：
$$\hat{L}(s,t,u,v) = \sum_{ijkl} \alpha_{ijkl} K((s,t,u,v), (s_i,t_j,u_k,v_l))$$

其中系数 $\alpha$ 通过求解线性系统 $(\mathbf{K} + \lambda\mathbf{I})\alpha = \mathbf{L}$ 获得。

**频域分析**：
插值核的频率响应决定了重建质量。理想低通滤波器：
$$\hat{K}(\omega) = \begin{cases}
1, & |\omega| < \omega_c \\
0, & |\omega| \geq \omega_c
\end{cases}$$

实际核的频率响应：
- 最近邻：$\hat{K}(\omega) = \text{sinc}^4(\omega/2\pi)$，有严重的频谱泄漏
- 三线性：$\hat{K}(\omega) = \text{sinc}^8(\omega/2\pi)$，改进但仍有混叠
- Lanczos：接近理想低通，但有Gibbs现象
- 高斯：$\hat{K}(\omega) = \exp(-\|\omega\|^2\sigma^2/2)$，无Gibbs但过度平滑

**预滤波理论**：
为避免混叠，采样前应用抗混叠滤波器：
$$L_{filtered}(s,t,u,v) = (L * \phi)(s,t,u,v)$$

其中 $\phi$ 是预滤波核，满足 $\hat{\phi}(\omega) = 0$ 当 $|\omega| > \pi/\Delta$。

### 4.1.3 体积渲染方程形式

将光场渲染统一到体积渲染框架，定义隐式光场体：

$$\sigma(\mathbf{x}) = \sum_{ijkl} \sigma_{ijkl} \cdot \delta(\mathbf{x} - \mathbf{x}_{ijkl})$$

$$c(\mathbf{x},\omega) = \sum_{ijkl} c_{ijkl} \cdot \delta(\mathbf{x} - \mathbf{x}_{ijkl}) \cdot \delta(\omega - \omega_{ijkl})$$

其中 $\mathbf{x}_{ijkl}$ 和 $\omega_{ijkl}$ 是从光线参数 $(s_i,t_j,u_k,v_l)$ 导出的3D位置和方向。

**从光线参数到3D坐标的映射**：
给定光线参数 $(s,t,u,v)$，沿光线的3D点为：
$$\mathbf{x}(\lambda) = \begin{pmatrix}
s + \lambda(u-s)/d \\
t + \lambda(v-t)/d \\
\lambda
\end{pmatrix}$$

**光场体的数学性质**：
定义光场测度：
$$dL = L(s,t,u,v) ds\,dt\,du\,dv$$

与体积渲染测度的关系：
$$dV = \sigma(\mathbf{x})c(\mathbf{x},\omega) d\mathbf{x}\,d\omega$$

通过Radon变换建立联系：
$$\mathcal{R}[\sigma c](s,t,u,v) = \int_{\mathbb{R}} \sigma(\mathbf{x}(\lambda))c(\mathbf{x}(\lambda),\omega) d\lambda = L(s,t,u,v)$$

体积渲染方程变为：

$$L(\mathbf{r}) = \int_0^{\infty} T(\tau) \cdot \sigma(\mathbf{r}(\tau)) \cdot c(\mathbf{r}(\tau), -\mathbf{r}'(\tau)) d\tau$$

其中透射率：
$$T(\tau) = \exp\left(-\int_0^{\tau} \sigma(\mathbf{r}(t)) dt\right)$$

**离散化与求积公式**：
使用分段常数近似：
$$L \approx \sum_{n=0}^{N-1} T_n \cdot \alpha_n \cdot c_n$$

其中：
- $\alpha_n = 1 - \exp(-\sigma_n \Delta_n)$ 是不透明度
- $T_n = \prod_{k=0}^{n-1}(1-\alpha_k)$ 是累积透射率
- $\Delta_n = t_{n+1} - t_n$ 是步长

**误差分析**：
离散化误差满足：
$$|L - L_{discrete}| \leq C \cdot \max_n \Delta_n \cdot \sup_t |\sigma''(t)|$$

使用Richardson外推可获得高阶精度：
$$L_{extrap} = \frac{4L(\Delta/2) - L(\Delta)}{3} + O(\Delta^4)$$

这种形式揭示了IBR与体积渲染的深层联系：光场的每个采样可视为空间中的一个方向性点光源，而连续光场则对应于体积密度场。

### 4.1.4 插值误差分析

插值误差可以通过Taylor展开分析：

$$\|L - \hat{L}\|_2 \leq C \cdot h^{p+1} \cdot \|D^{p+1}L\|_{\infty}$$

其中 $h$ 是采样间隔，$p$ 是插值核的阶数。

**误差的概率表征**：
假设光场是随机过程，具有协方差函数 $C_L(\Delta s, \Delta t, \Delta u, \Delta v)$。插值误差的方差为：

$$\text{Var}[L - \hat{L}] = C_L(0,0,0,0) - \sum_{ijkl} K^2(s-s_i, t-t_j, u-u_k, v-v_l) C_L(s-s_i, t-t_j, u-u_k, v-v_l)$$

对于平稳过程，这简化为：
$$\text{Var}[L - \hat{L}] = C_L(0,0,0,0)(1 - \int_{\mathbb{R}^4} |\hat{K}(\omega)|^2 S_L(\omega) d\omega)$$

其中 $S_L$ 是功率谱密度。

**误差的频域表征**：
对于带限信号，频域分析给出更紧的界：

$$\mathcal{E}(\omega) = |1 - \hat{K}(\omega)| \cdot |\tilde{L}(\omega)|$$

总误差能量：
$$E_{error} = \int_{\mathbb{R}^4} |\mathcal{E}(\omega)|^2 d\omega = \int_{\mathbb{R}^4} |1 - \hat{K}(\omega)|^2 |\tilde{L}(\omega)|^2 d\omega$$

**局部误差估计**：
使用多元Taylor展开，在点 $(s_0,t_0,u_0,v_0)$ 附近：

$$L(s,t,u,v) = \sum_{|\alpha| \leq p} \frac{D^{\alpha}L(s_0,t_0,u_0,v_0)}{\alpha!}(s-s_0)^{\alpha_1}(t-t_0)^{\alpha_2}(u-u_0)^{\alpha_3}(v-v_0)^{\alpha_4} + R_p$$

余项 $R_p$ 满足：
$$|R_p| \leq \frac{M_{p+1}}{(p+1)!} \|(s-s_0, t-t_0, u-u_0, v-v_0)\|^{p+1}$$

其中 $M_{p+1} = \sup |D^{p+1}L|$。

**Bramble-Hilbert引理应用**：
对于Sobolev空间 $W^{k,p}$ 中的函数，插值误差满足：
$$\|L - \Pi_h L\|_{W^{m,p}} \leq C h^{k-m} |L|_{W^{k,p}}$$

其中 $\Pi_h$ 是插值算子，$|\cdot|_{W^{k,p}}$ 是半范数。

**自适应采样策略**：
基于局部误差估计，可设计自适应采样密度：
$$\rho(s,t,u,v) \propto \left(\sum_{|\alpha|=2} |D^{\alpha}L(s,t,u,v)|^2\right)^{1/2}$$

**最优采样理论**：
给定总采样数 $N$，最优采样密度最小化积分误差：
$$\rho^* = \arg\min_{\rho} \int_{\Omega} \frac{\|\nabla L\|^2}{\rho} d\Omega \quad \text{s.t.} \quad \int_{\Omega} \rho d\Omega = N$$

使用Lagrange乘数法，解为：
$$\rho^*(s,t,u,v) = N \frac{\|\nabla L(s,t,u,v)\|^{2/3}}{\int_{\Omega} \|\nabla L\|^{2/3} d\Omega}$$

这确保高频区域获得更密集的采样。

## 4.2 视图合成与极线约束

### 4.2.1 极线几何基础

给定两个相机矩阵 $P_1 = K_1[R_1|t_1]$ 和 $P_2 = K_2[R_2|t_2]$，基础矩阵 $F$ 满足：

$$\mathbf{x}_2^T F \mathbf{x}_1 = 0$$

其中 $\mathbf{x}_1$, $\mathbf{x}_2$ 是对应点的齐次坐标。

**基础矩阵的导出**：
从3D点 $\mathbf{X}$ 投影到两个视图：
$$\mathbf{x}_1 = P_1\mathbf{X}, \quad \mathbf{x}_2 = P_2\mathbf{X}$$

消去 $\mathbf{X}$ 得到极线约束。基础矩阵可分解为：

$$F = K_2^{-T} [t_{21}]_{\times} R_{21} K_1^{-1}$$

这里 $[t]_{\times}$ 是反对称矩阵：
$$[t]_{\times} = \begin{pmatrix}
0 & -t_z & t_y \\
t_z & 0 & -t_x \\
-t_y & t_x & 0
\end{pmatrix}$$

$R_{21} = R_2 R_1^T$，$t_{21} = t_2 - R_{21}t_1$。

**几何解释**：
极线约束的几何意义是：点 $\mathbf{x}_1$ 对应的极线 $l_2 = F\mathbf{x}_1$ 是所有可能的对应点 $\mathbf{x}_2$ 的轨迹。这条极线是光心 $C_1$、点 $\mathbf{x}_1$ 和极点 $\mathbf{e}_2$ 构成的平面与图像平面 $\Pi_2$ 的交线。

**基础矩阵的性质**：
1. **秩为2**：$\text{rank}(F) = 2$，因为 $\det([t]_{\times}) = 0$
2. **极点关系**：$F^T\mathbf{e}_2 = 0$，$F\mathbf{e}_1 = 0$，其中 $\mathbf{e}_1, \mathbf{e}_2$ 是极点
3. **对应关系**：点 $\mathbf{x}_1$ 在第二幅图像中的极线为 $l_2 = F\mathbf{x}_1$
4. **对偶性**：若 $F$ 是从视图1到视图2的基础矩阵，则 $F^T$ 是从视图2到视图1的基础矩阵

**SVD分解**：
基础矩阵的SVD分解：
$$F = U\Sigma V^T = U\begin{pmatrix}\sigma_1 & 0 & 0\\0 & \sigma_2 & 0\\0 & 0 & 0\end{pmatrix}V^T$$

其中 $\sigma_1 \geq \sigma_2 > 0$。极点为：
- $\mathbf{e}_1 = V(:,3)$（右零空间）
- $\mathbf{e}_2 = U(:,3)$（左零空间）

**本质矩阵**：
对于校准相机，本质矩阵 $E = [t_{21}]_{\times} R_{21}$ 满足：
$$\hat{\mathbf{x}}_2^T E \hat{\mathbf{x}}_1 = 0$$

其中 $\hat{\mathbf{x}}$ 是归一化坐标。本质矩阵的特殊性质：
- 两个非零奇异值相等：$\sigma_1 = \sigma_2$
- 满足约束：$2E^TE\cdot E - \text{tr}(E^TE)\cdot E = 0$
- 5个自由度（3个旋转 + 2个平移方向）

**从本质矩阵恢复运动**：
给定SVD分解 $E = U\text{diag}(1,1,0)V^T$，可能的解为：
$$R = U\text{diag}(1,1,\pm 1)V^T$$
$$[t]_{\times} = U\text{diag}(1,1,0)W^{\pm}U^T$$

其中 $W^{\pm} = \begin{pmatrix}0 & \mp 1 & 0\\\pm 1 & 0 & 0\\0 & 0 & 1\end{pmatrix}$。有四种可能的解，需要通过深度约束消歧义。

### 4.2.2 深度引导的视图合成

给定参考视图的深度图 $D(\mathbf{x}_1)$，目标视图的像素可通过以下变换获得：

$$\mathbf{x}_2 = K_2 R_{21} K_1^{-1} \mathbf{x}_1 + \frac{K_2 \mathbf{t}_{21}}{D(\mathbf{x}_1)}$$

**推导过程**：
设3D点 $\mathbf{X}$ 在相机1坐标系中为：
$$\mathbf{X}_1 = D(\mathbf{x}_1) K_1^{-1} \mathbf{x}_1$$

转换到相机2坐标系：
$$\mathbf{X}_2 = R_{21}\mathbf{X}_1 + \mathbf{t}_{21}$$

投影到图像平面：
$$\mathbf{x}_2 = K_2\frac{\mathbf{X}_2}{Z_2} = K_2\left(R_{21}K_1^{-1}\mathbf{x}_1 + \frac{\mathbf{t}_{21}}{D(\mathbf{x}_1)}\right)$$

**单应矩阵形式**：
定义单应矩阵：
$$H = K_2(R_{21} + \mathbf{t}_{21}\mathbf{n}^T/d)K_1^{-1}$$

其中 $\mathbf{n}$ 是平面法向量，$d$ 是平面到相机1的距离。对于平面上的点：
$$\mathbf{x}_2 = H\mathbf{x}_1$$

**视差形式**：
对于平行相机配置（纯水平平移），公式简化为：

$$x_2 - x_1 = d(\mathbf{x}_1) = \frac{bf}{D(\mathbf{x}_1)}$$

其中 $b$ 是基线长度，$f$ 是焦距，$d$ 是视差。

**逆深度参数化**：
使用逆深度 $\rho = 1/D$ 可以线性化变换：
$$\mathbf{x}_2 = K_2 R_{21} K_1^{-1} \mathbf{x}_1 + \rho(\mathbf{x}_1) K_2 \mathbf{t}_{21}$$

这在优化中更稳定，因为逆深度的分布更均匀。

**深度不确定性传播**：
深度误差对重投影的影响：
$$\frac{\partial \mathbf{x}_2}{\partial D} = -\frac{K_2 \mathbf{t}_{21}}{D^2}$$

给定深度标准差 $\sigma_D$，重投影误差协方差：
$$\Sigma_{\mathbf{x}_2} = \left(\frac{\partial \mathbf{x}_2}{\partial D}\right)\sigma_D^2\left(\frac{\partial \mathbf{x}_2}{\partial D}\right)^T$$

**一阶误差传播**：
考虑完整的误差源（深度误差、标定误差、姿态误差）：
$$\delta\mathbf{x}_2 = \frac{\partial \mathbf{x}_2}{\partial D}\delta D + \frac{\partial \mathbf{x}_2}{\partial \mathbf{K}_1}\delta\mathbf{K}_1 + \frac{\partial \mathbf{x}_2}{\partial \mathbf{R}}\delta\mathbf{R} + \frac{\partial \mathbf{x}_2}{\partial \mathbf{t}}\delta\mathbf{t}$$

总误差协方差：
$$\Sigma_{total} = J_D\Sigma_D J_D^T + J_K\Sigma_K J_K^T + J_R\Sigma_R J_R^T + J_t\Sigma_t J_t^T$$

其中 $J_*$ 是对应的雅可比矩阵。

### 4.2.3 遮挡处理与体积积分

遮挡可以通过体积渲染自然处理。定义遮挡感知的传输函数：

$$T(\mathbf{x}_1 \to \mathbf{x}_2) = \exp\left(-\int_0^1 \sigma(\gamma(t)) dt\right)$$

其中 $\gamma(t)$ 是连接对应3D点的路径。

**前向映射与后向映射**：
- **前向映射**：从源视图到目标视图
  $$I_2^{forward}(\mathbf{x}_2) = \sum_{\mathbf{x}_1 \to \mathbf{x}_2} T(\mathbf{x}_1 \to \mathbf{x}_2) \cdot I_1(\mathbf{x}_1)$$
  
- **后向映射**：从目标视图到源视图
  $$I_2^{backward}(\mathbf{x}_2) = \sum_{\mathbf{x}_1} W(\mathbf{x}_1, \mathbf{x}_2) \cdot I_1(\mathbf{x}_1)$$

**贝叶斯遮挡推理**：
将遮挡建模为二元随机变量 $O \in \{0,1\}$：
$$P(I_2|I_1,D_1,D_2) = \sum_{O} P(I_2|I_1,D_1,D_2,O)P(O|D_1,D_2)$$

其中：
- $P(O=1|D_1,D_2) = \exp(-\lambda|D_2 - D_{proj}|)$（遮挡概率）
- $P(I_2|I_1,O=1) = \mathcal{N}(I_1, \sigma_I^2)$（可见时的强度分布）
- $P(I_2|I_1,O=0) = \mathcal{U}(0,255)$（遮挡时的均匀分布）

**遮挡检测**：
使用深度一致性检查：
$$\mathcal{O}(\mathbf{x}_1, \mathbf{x}_2) = \begin{cases}
1, & |D_2(\mathbf{x}_2) - D_{proj}(\mathbf{x}_1)| < \epsilon \\
0, & \text{otherwise}
\end{cases}$$

其中 $D_{proj}$ 是从视图1投影到视图2的深度。

**多视图遮挡累积**：
对于 $N$ 个视图，可见性测度：
$$V(\mathbf{x}) = \prod_{i=1}^N V_i(\mathbf{x}) = \prod_{i=1}^N \exp\left(-\int_{C_i}^{\mathbf{x}} \sigma(\mathbf{r}) d\mathbf{r}\right)$$

其中 $C_i$ 是第 $i$ 个相机中心。

**软遮挡模型**：
使用sigmoid函数建模遮挡边界：
$$W_{occ}(\mathbf{x}) = \frac{1}{1 + \exp(-k(D_{front} - D_{back}))}$$

**光线空间遮挡处理**：
在光场空间 $(s,t,u,v)$ 中，遮挡表现为不连续性。定义遮挡函数：
$$\mathcal{O}(s,t,u,v) = H(d(s,t,u,v) - d_{threshold})$$

其中 $H$ 是Heaviside函数，$d$ 是深度函数。

**变分遮挡处理**：
通过最小化能量泛函处理遮挡：
$$E[u,v] = \int_{\Omega} (I_2 - I_1 \circ \phi)^2 v\,dx + \lambda_1\int_{\Omega} |\nabla u|\,dx + \lambda_2\int_{\Omega} |\nabla v|\,dx$$

其中 $u$ 是深度场，$v$ 是遮挡指示函数，$\phi$ 是变形场。

合成的像素值为：

$$I_2(\mathbf{x}_2) = \sum_{\mathbf{x}_1} T(\mathbf{x}_1 \to \mathbf{x}_2) \cdot I_1(\mathbf{x}_1) \cdot K(\mathbf{x}_2 - \pi_2(\pi_1^{-1}(\mathbf{x}_1)))$$

### 4.2.4 多视图约束

对于 $N$ 个视图，约束系统变为：

$$\mathcal{M}\mathbf{x} = \begin{bmatrix}
F_{12} & -I & 0 & \cdots \\
F_{13} & 0 & -I & \cdots \\
\vdots & & \ddots & \\
F_{1N} & 0 & \cdots & -I
\end{bmatrix} \begin{bmatrix}
\mathbf{x}_1 \\ \mathbf{x}_2 \\ \vdots \\ \mathbf{x}_N
\end{bmatrix} = 0$$

**约束系统的秩分析**：
系统矩阵 $\mathcal{M}$ 的秩为：
$$\text{rank}(\mathcal{M}) = \min(3(N-1), 2N-\text{rank}(\mathcal{C}))$$

其中 $\mathcal{C}$ 是相机中心的配置矩阵。当所有相机中心共线时，系统退化。

**三视图张量**：
对于三个视图，存在三焦张量 $\mathcal{T}_{ijk}$ 满足：
$$[\mathbf{x}_2]_{\times}^i [\mathbf{x}_3]_{\times}^j \mathcal{T}_{ijk} \mathbf{x}_1^k = 0$$

三焦张量的分解：
$$\mathcal{T}_{ijk} = \sum_{\alpha=1}^3 a_i^{\alpha} b_j^{\alpha} c_k^{\alpha}$$

其中：
- $a_i^{\alpha} = [P_2^{\alpha}]_i$（第2个相机矩阵的第$\alpha$列）
- $b_j^{\alpha} = [P_3^{\alpha}]_j$（第3个相机矩阵的第$\alpha$列）  
- $c_k^{\alpha} = (-1)^{\alpha+1}[P_1^{\bar{\alpha}}]_k$（第1个相机矩阵去掉第$\alpha$列）

三焦张量编码了三视图间的所有几何关系，包含27个元素但只有18个自由度。

**从三焦张量提取基础矩阵**：
$$F_{21} = [\mathbf{e}_2]_{\times}[\mathcal{T}_{1jk}\mathbf{e}_3^j\mathbf{e}_3^k]$$
$$F_{31} = [\mathbf{e}_3]_{\times}[\mathcal{T}_{i1k}\mathbf{e}_2^i\mathbf{e}_2^k]$$
$$F_{32} = [\mathbf{e}_3]_{\times}[\mathcal{T}_{ij1}\mathbf{e}_2^i\mathbf{e}_2^j]$$

**鲁棒估计**：
使用RANSAC算法估计多视图几何：
1. 随机选择最小样本集（F矩院7点，三焦张量6点）
2. 计算几何模型（F矩阵或三焦张量）
3. 计算内点数：$|\{\mathbf{x} : d(\mathbf{x}, \mathcal{M}) < \tau\}|$
4. 迭代直到找到足够好的模型

**迭代次数理论界限**：
$$N_{iter} = \frac{\log(1-p)}{\log(1-w^m)}$$

其中 $p$ 是置信度，$w$ 是内点比例，$m$ 是最小样本数。

**Bundle Adjustment**：
联合优化所有相机参数和3D点：
$$\min_{\{P_i\}, \{\mathbf{X}_j\}} \sum_{i,j} \rho\left(\|\mathbf{x}_{ij} - \pi(P_i\mathbf{X}_j)\|^2\right)$$

其中 $\rho$ 是鲁棒核函数（如Huber函数）：
$$\rho_{Huber}(x) = \begin{cases}
\frac{1}{2}x^2, & |x| \leq \delta \\
\delta(|x| - \frac{1}{2}\delta), & |x| > \delta
\end{cases}$$

**Schur补优化**：
Bundle Adjustment的Hessian矩阵具有特殊结构：
$$H = \begin{bmatrix}
U & W \\
W^T & V
\end{bmatrix}$$

其中 $U$ 对应相机参数，$V$ 对应3D点，$W$ 是交叉项。使用Schur补：
$$S = U - WV^{-1}W^T$$

可以高效求解法方程。

## 4.3 光图与表面光场

### 4.3.1 几何代理的引入

光图（Lumigraph）通过引入近似几何 $\mathcal{S}$ 改进了纯光场方法：

$$L(s,t,u,v) = L_{\mathcal{S}}(x(s,t,u,v), \omega(s,t,u,v))$$

其中 $(x,\omega)$ 是光线与表面 $\mathcal{S}$ 的交点和方向。

### 4.3.2 表面参数化

对于表面点 $x \in \mathcal{S}$，定义局部参数化：

$$x = x(u,v), \quad \omega = \omega(\theta, \phi)$$

表面光场变为：

$$L_{\mathcal{S}}(u,v,\theta,\phi) = L(x(u,v), \omega(\theta,\phi))$$

这将6D函数降为4D，显著减少存储需求。

### 4.3.3 双层表示

考虑带有透明度的双层模型：

$$L_{out} = \alpha_1 L_1 + (1-\alpha_1)\alpha_2 L_2 + (1-\alpha_1)(1-\alpha_2)L_{bg}$$

这可以推广到 $N$ 层：

$$L_{out} = \sum_{i=1}^N L_i \prod_{j=1}^{i-1}(1-\alpha_j) + L_{bg}\prod_{j=1}^N(1-\alpha_j)$$

### 4.3.4 与体积渲染的统一

表面光场可视为体积密度在表面的 $\delta$ 函数：

$$\sigma(x) = \sigma_s(x) \cdot \delta(d(x))$$

其中 $d(x)$ 是到表面的符号距离。体积渲染方程简化为：

$$L(r) = \sum_{i} T_i \cdot c_i(r(t_i), -r'(t_i))$$

其中求和遍历所有光线-表面交点。

## 4.4 傅里叶切片定理应用

### 4.4.1 光场的频域表示

4D光场的傅里叶变换：

$$\tilde{L}(\omega_s, \omega_t, \omega_u, \omega_v) = \mathcal{F}\{L(s,t,u,v)\}$$

对于Lambertian场景，频谱集中在低频：

$$|\tilde{L}(\omega)| \propto \frac{1}{|\omega|^2}$$

### 4.4.2 切片-投影定理

2D图像是4D光场的切片：

$$I(u,v) = L(s_0, t_0, u, v)$$

其傅里叶变换满足：

$$\tilde{I}(\omega_u, \omega_v) = \tilde{L}(0, 0, \omega_u, \omega_v)$$

这是经典投影-切片定理的推广。

### 4.4.3 频域重建

利用切片定理，可从多个2D切片重建4D光场：

$$\tilde{L}(\omega) = \sum_k W_k(\omega) \cdot \tilde{I}_k(\Pi_k(\omega))$$

其中 $\Pi_k$ 是到第 $k$ 个切片的投影算子，$W_k$ 是重建权重。

### 4.4.4 带宽分析

场景深度范围 $[D_{min}, D_{max}]$ 决定了最大频率：

$$\omega_{max} = \frac{2\pi B}{D_{min}}$$

其中 $B$ 是相机基线。这给出了最小采样率：

$$\Delta s = \Delta u < \frac{\pi}{\omega_{max}} = \frac{D_{min}}{2B}$$

## 4.5 IBR的采样理论与计算复杂度

### 4.5.1 Nyquist采样率推导

考虑场景中点 $P$ 在深度 $D$ 处，其在相邻相机间的视差为：

$$\delta = \frac{B \cdot f}{D}$$

为避免混叠，采样间隔必须满足：

$$\Delta_{cam} \leq \frac{\lambda D}{2f}$$

其中 $\lambda$ 是像素间距。对于典型参数（$\lambda = 10\mu m$, $f = 50mm$），在 $D = 1m$ 处：

$$\Delta_{cam} \leq 0.1mm$$

### 4.5.2 最小采样密度

考虑 $M \times N$ 的相机阵列覆盖 $A \times B$ 的区域，总采样数为：

$$N_{total} = MN \cdot W \cdot H$$

其中 $W \times H$ 是图像分辨率。存储需求为：

$$S = N_{total} \cdot b \cdot c$$

其中 $b$ 是每通道位深，$c$ 是通道数。对于4K RGB图像的100×100阵列：

$$S = 10^4 \cdot 4096 \cdot 2160 \cdot 3 \cdot 8 \approx 2.1 \text{ PB}$$

### 4.5.3 计算复杂度分析

**渲染复杂度**：
- 暴力搜索：$O(MN)$ per pixel
- 使用加速结构：$O(\log(MN))$ per pixel
- 总复杂度：$O(WH\log(MN))$

**预处理复杂度**：
- 深度估计：$O(MN \cdot WH \cdot D)$，其中 $D$ 是视差搜索范围
- 几何重建：$O((MN)^2 \cdot WH)$ 最坏情况

### 4.5.4 渐进式采样策略

定义重要性函数：

$$I(s,t) = \|\nabla L(s,t,\cdot,\cdot)\|_2$$

自适应采样通过求解优化问题：

$$\min_{\{(s_i,t_i)\}} \int\int |L(s,t) - \hat{L}(s,t)|^2 \, ds \, dt$$

使用贪心算法：
1. 初始化均匀采样
2. 计算重建误差 $\epsilon(s,t)$
3. 在 $\arg\max \epsilon(s,t)$ 处添加样本
4. 重复直到 $\|\epsilon\|_{\infty} < \tau$

### 4.5.5 压缩与流式传输

**频域压缩**：
利用光场的频谱稀疏性：

$$\tilde{L}_{compressed} = \mathcal{T}_{\epsilon}(\tilde{L})$$

其中 $\mathcal{T}_{\epsilon}$ 是阈值算子，保留大于 $\epsilon$ 的系数。

**视图预测编码**：
$$L_i = L_{ref} + \Delta L_i$$

其中 $\Delta L_i$ 是残差，通常具有更好的压缩率。

压缩比可达：
$$CR = \frac{S_{original}}{S_{compressed}} \approx 100:1 \text{ 到 } 1000:1$$

### 4.5.6 实时渲染优化

**GPU并行化**：
```
每个线程处理一个输出像素：
1. 计算光线参数 (s,t,u,v)
2. 定位最近的4×4×4×4采样点
3. 执行4D插值
4. 输出颜色值
```

复杂度降为 $O(1)$ per pixel，实现60+ FPS的4K渲染。

**层次化表示**：
构建光场金字塔：

$$L_l = \mathcal{D}_2(L_{l-1})$$

其中 $\mathcal{D}_2$ 是2倍下采样。根据需要的精度选择层级：

$$l^* = \max\{l : \text{res}(L_l) \geq \text{res}_{target}\}$$

## 本章小结

本章将基于图像的渲染技术统一在体积渲染框架下，主要贡献包括：

1. **4D光场插值**：展示了光场渲染本质上是4D函数的重建问题，通过适当的插值核可以实现高质量的视图合成。

2. **极线几何约束**：推导了多视图间的几何关系，展示了如何利用这些约束改进视图合成质量。

3. **表面光场统一**：通过引入 $\delta$ 函数形式的体积密度，将表面光场纳入体积渲染方程。

4. **频域分析**：应用傅里叶切片定理分析了光场的频谱特性，为采样率选择提供理论依据。

5. **复杂度界限**：建立了IBR系统的存储需求 $O(MN \cdot WH)$ 和计算复杂度 $O(WH\log(MN))$ 的精确界限。

关键公式汇总：
- 4D插值：$\hat{L}(s,t,u,v) = \sum_{ijkl} L_{ijkl} \cdot K(s-s_i, t-t_j, u-u_k, v-v_l)$
- 极线约束：$x_2^T F x_1 = 0$
- 最小采样率：$\Delta_{cam} \leq \frac{\lambda D}{2f}$
- 频带限制：$\omega_{max} = \frac{2\pi B}{D_{min}}$

## 练习题

### 基础题

**习题 4.1**：推导双平面参数化下，光线方向 $\omega$ 与参数 $(s,t,u,v)$ 的关系。

*提示*：考虑两平面间的距离为 $d$，使用相似三角形。

<details>
<summary>答案</summary>

光线方向为：
$$\omega = \frac{1}{\sqrt{(u-s)^2 + (v-t)^2 + d^2}} \begin{pmatrix} u-s \\ v-t \\ d \end{pmatrix}$$

归一化确保 $|\omega| = 1$。
</details>

**习题 4.2**：证明对于Lambertian表面，4D光场可以用2D纹理完全表示。

*提示*：Lambertian表面的辐射与观察方向无关。

<details>
<summary>答案</summary>

对于Lambertian表面，$L(x, \omega) = \rho(x) \cdot L_i(x)$，与 $\omega$ 无关。因此：
$$L(s,t,u,v) = \rho(x(s,t,u,v)) \cdot L_i(x(s,t,u,v))$$
只需存储每个表面点的反射率 $\rho(x)$，这是2D参数化。
</details>

**习题 4.3**：计算存储一个 $10 \times 10 \times 10 \times 10$ 的光场所需的内存，假设每个样本是1024×1024的RGB图像。

*提示*：考虑每个通道8位精度。

<details>
<summary>答案</summary>

内存需求：
$$M = 10^4 \times 1024^2 \times 3 \times 1 \text{ byte} = 3.15 \text{ GB}$$

若使用16位精度，需求翻倍至6.3 GB。
</details>

### 挑战题

**习题 4.4**：推导非均匀采样下的光场重建误差界，假设采样密度函数为 $\rho(s,t)$。

*提示*：使用Voronoi单元分析局部误差。

<details>
<summary>答案</summary>

定义Voronoi单元 $V_i$ 对应样本 $i$，局部误差：
$$\epsilon_i = \sup_{(s,t) \in V_i} |L(s,t) - L(s_i,t_i)|$$

使用Lipschitz连续性：
$$\epsilon_i \leq K \cdot \text{diam}(V_i) \leq \frac{K}{\sqrt{\rho(s_i,t_i)}}$$

总误差：
$$\|\epsilon\|_2 = \sqrt{\sum_i |V_i| \epsilon_i^2} \leq K \sqrt{\int\int \frac{1}{\rho(s,t)} ds dt}$$
</details>

**习题 4.5**：设计一个自适应IBR系统，根据场景内容动态分配采样密度。推导最优采样分布。

*提示*：将其形式化为变分问题。

<details>
<summary>答案</summary>

定义能量泛函：
$$E[\rho] = \int\int |L - \hat{L}_\rho|^2 ds dt + \lambda \int\int \rho(s,t) ds dt$$

使用变分法，最优密度满足：
$$\rho^*(s,t) \propto \|\nabla L(s,t)\|^{2/3}$$

这给出了基于局部频率内容的自适应采样。
</details>

**习题 4.6**：证明对于带限光场，存在有限采样率使得重建误差为零。

*提示*：使用Shannon采样定理的4D推广。

<details>
<summary>答案</summary>

若 $\tilde{L}(\omega) = 0$ 对所有 $|\omega| > \Omega$，则Shannon定理给出：
$$L(s,t,u,v) = \sum_{ijkl} L_{ijkl} \cdot \text{sinc}(\pi(s-s_i)/\Delta) \cdots$$

当采样率 $1/\Delta > 2\Omega/\pi$ 时，重建是精确的。对于实际场景，深度范围限制了最大频率，因此存在有限采样率。
</details>

### 计算实现题

**习题 4.7**：实现4D光场的频域压缩算法，分析压缩率与重建质量的关系。

*提示*：使用4D DCT或小波变换。

<details>
<summary>答案</summary>

算法框架：
1. 计算4D DCT：$\tilde{L} = \text{DCT4D}(L)$
2. 阈值处理：保留最大的 $k$ 个系数
3. 逆变换：$\hat{L} = \text{IDCT4D}(\tilde{L}_{sparse})$

压缩率：$CR = N^4 / k$
质量度量：$\text{PSNR} = 20\log_{10}(\text{MAX} / \text{RMSE})$

典型结果：CR = 100时，PSNR > 35dB。
</details>

**习题 4.8**：推导并实现基于深度的光场渲染加速算法。

*提示*：利用深度信息减少4D搜索到2D。

<details>
<summary>答案</summary>

给定深度 $D(u,v)$，光线参数约束为：
$$s = u - \frac{(u-u_0)d}{D(u,v)}, \quad t = v - \frac{(v-v_0)d}{D(u,v)}$$

这将4D查找降为2D，复杂度从 $O(N^4)$ 降至 $O(N^2)$。实现时使用深度缓冲加速遮挡测试。
</details>

## 常见陷阱与错误

1. **混淆光场参数化**：$(s,t,u,v)$ 不是空间坐标，而是光线参数。常见错误是直接将其当作4D空间点。

2. **忽略极线约束的退化情况**：当相机光心共线时，基础矩阵秩降为2，标准算法失效。需要特殊处理。

3. **采样不足导致的混叠**：未考虑场景最小深度，导致采样率过低。症状是视角变化时出现跳变。

4. **插值核选择不当**：使用高阶插值核但采样稀疏，导致振铃效应。应根据采样密度选择合适的核。

5. **频域分析的周期性假设**：DFT假设信号周期延拓，导致边界伪影。使用窗函数或镜像填充缓解。

6. **深度不连续处理**：在遮挡边界简单插值产生"橡皮片"效应。需要显式的遮挡处理。

## 最佳实践检查清单

### 系统设计阶段
- [ ] 分析场景深度范围，计算所需最小采样率
- [ ] 选择合适的参数化方式（双平面、球面、表面）
- [ ] 评估存储需求，设计压缩策略
- [ ] 确定实时性要求，选择相应的加速结构

### 数据采集阶段
- [ ] 校准所有相机的内外参数
- [ ] 确保采样密度满足Nyquist准则
- [ ] 记录光照条件，考虑时变光照的影响
- [ ] 采集几何代理（深度图、点云）以改进重建

### 算法实现阶段
- [ ] 实现鲁棒的极线几何估计
- [ ] 选择适合数据特性的插值核
- [ ] 实现高效的4D数据结构（如4D树）
- [ ] 添加遮挡处理和边界混合

### 质量评估阶段
- [ ] 使用保留视图验证重建质量
- [ ] 分析频谱确认无混叠
- [ ] 测试极端视角下的表现
- [ ] 评估压缩对质量的影响

### 优化阶段
- [ ] 实现GPU并行渲染管线
- [ ] 使用层次化表示加速
- [ ] 应用视锥体裁剪减少计算
- [ ] 实现自适应质量控制
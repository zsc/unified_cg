# 第4章：基于图像的渲染

基于图像的渲染（Image-Based Rendering, IBR）代表了计算机图形学中的一个范式转变：从基于几何的渲染转向基于采样的渲染。本章将探讨如何使用预先捕获的图像来合成新视角，并将这些技术统一在体积渲染框架下。我们将深入研究光场的数学表示、视图合成的几何约束，以及采样理论在IBR中的应用。

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

### 4.1.2 离散采样与插值

给定一组采样光线 $\{L_{ijkl}\}$，新视角的光线通过4D插值获得：

$$\hat{L}(s,t,u,v) = \sum_{i,j,k,l} L_{ijkl} \cdot K(s-s_i, t-t_j, u-u_k, v-v_l)$$

其中 $K$ 是4D插值核。常见选择包括：
- 最近邻：$K = \delta$
- 三线性：$K = \prod_{d} (1-|x_d|)_+$
- Lanczos：$K = \prod_{d} \text{sinc}(x_d)\text{sinc}(x_d/a)$

### 4.1.3 体积渲染方程形式

将光场渲染统一到体积渲染框架，定义隐式光场体：

$$\sigma(x) = \sum_{ijkl} \sigma_{ijkl} \cdot \delta(x - x_{ijkl})$$

$$c(x,\omega) = \sum_{ijkl} c_{ijkl} \cdot \delta(x - x_{ijkl}) \cdot \delta(\omega - \omega_{ijkl})$$

其中 $x_{ijkl}$ 和 $\omega_{ijkl}$ 是从光线参数 $(s_i,t_j,u_k,v_l)$ 导出的3D位置和方向。

体积渲染方程变为：

$$L(r) = \int_0^{\infty} T(t) \cdot \sigma(r(t)) \cdot c(r(t), -r'(t)) dt$$

这种形式揭示了IBR与体积渲染的深层联系。

### 4.1.4 插值误差分析

插值误差可以通过Taylor展开分析：

$$\|L - \hat{L}\|_2 \leq C \cdot h^{p+1} \cdot \|D^{p+1}L\|_{\infty}$$

其中 $h$ 是采样间隔，$p$ 是插值核的阶数。对于带限信号，频域分析给出更紧的界：

$$\mathcal{E}(\omega) = |1 - \hat{K}(\omega)| \cdot |\tilde{L}(\omega)|$$

## 4.2 视图合成与极线约束

### 4.2.1 极线几何基础

给定两个相机矩阵 $P_1 = K_1[R_1|t_1]$ 和 $P_2 = K_2[R_2|t_2]$，基础矩阵 $F$ 满足：

$$x_2^T F x_1 = 0$$

其中 $x_1$, $x_2$ 是对应点的齐次坐标。基础矩阵可分解为：

$$F = K_2^{-T} [t_{21}]_{\times} R_{21} K_1^{-1}$$

这里 $[t]_{\times}$ 是反对称矩阵，$R_{21} = R_2 R_1^T$，$t_{21} = t_2 - R_{21}t_1$。

### 4.2.2 深度引导的视图合成

给定参考视图的深度图 $D(x_1)$，目标视图的像素可通过以下变换获得：

$$x_2 = K_2 R_{21} K_1^{-1} x_1 + \frac{K_2 t_{21}}{D(x_1)}$$

这个公式可以改写为视差（disparity）形式：

$$x_2 - x_1 = d(x_1) \cdot v_{21}$$

其中 $d(x_1) = \frac{b}{D(x_1)}$ 是视差，$b$ 是基线长度。

### 4.2.3 遮挡处理与体积积分

遮挡可以通过体积渲染自然处理。定义遮挡感知的传输函数：

$$T(x_1 \to x_2) = \exp\left(-\int_0^1 \sigma(\gamma(t)) dt\right)$$

其中 $\gamma(t)$ 是连接 $x_1$ 和 $x_2$ 的3D路径。

合成的像素值为：

$$I_2(x_2) = \sum_{x_1} T(x_1 \to x_2) \cdot I_1(x_1) \cdot K(x_2 - \pi_2(\pi_1^{-1}(x_1)))$$

### 4.2.4 多视图约束

对于 $N$ 个视图，约束系统变为：

$$\begin{bmatrix}
F_{12} & 0 & \cdots \\
F_{13} & 0 & \cdots \\
\vdots & \ddots & \\
0 & \cdots & F_{(N-1)N}
\end{bmatrix} \begin{bmatrix}
x_1 \\ x_2 \\ \vdots \\ x_N
\end{bmatrix} = 0$$

这个超定系统可通过最小二乘求解，提供鲁棒的对应关系。

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
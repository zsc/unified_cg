# 第24章：全息显示与计算全息

## 章节概要

本章探讨全息术在计算机图形学中的应用，将传统的全息记录原理与现代计算方法相结合。我们将从物理全息的基本原理出发，发展到计算机生成全息图（CGH）的算法实现，并讨论如何将全息技术集成到现代渲染管线中。通过将全息术表述为体积渲染方程的特殊形式，我们建立了与前述章节的理论联系。

### 学习目标

完成本章后，您将能够：
1. 理解全息记录与重建的物理原理及其数学表述
2. 推导计算机生成全息图的各种算法
3. 分析空间光调制器的工作原理与限制
4. 实现相位恢复算法解决全息重建问题
5. 设计完整的全息渲染管线并分析其计算复杂度

## 24.1 全息记录与重建原理

### 24.1.1 全息的物理基础

全息术基于光波的干涉与衍射原理。与传统成像只记录光的强度不同，全息同时记录光波的振幅和相位信息。考虑物光波 $U_o(\mathbf{x})$ 和参考光波 $U_r(\mathbf{x})$，在记录平面上的总光场为：

$$U_{total}(\mathbf{x}) = U_o(\mathbf{x}) + U_r(\mathbf{x})$$

记录介质响应光强度 $I(\mathbf{x}) = |U_{total}(\mathbf{x})|^2$：

$$I(\mathbf{x}) = |U_o|^2 + |U_r|^2 + U_o U_r^* + U_o^* U_r$$

其中后两项包含了物光波的相位信息，这是全息记录的关键。

将复振幅写成极坐标形式 $U = |U|e^{i\phi}$，干涉项展开为：

$$U_o U_r^* + U_o^* U_r = 2|U_o||U_r|\cos(\phi_o - \phi_r)$$

这表明干涉条纹的对比度由振幅乘积决定，条纹间距由相位差梯度决定：

$$\Lambda = \frac{2\pi}{|\nabla(\phi_o - \phi_r)|}$$

对于典型的离轴全息配置，参考光与物光夹角为 $\theta$，条纹间距约为：

$$\Lambda \approx \frac{\lambda}{2\sin(\theta/2)}$$

### 24.1.2 菲涅尔全息图

对于菲涅尔全息，参考光为球面波。设参考点源位于 $\mathbf{r}_r$，则：

$$U_r(\mathbf{x}) = \frac{A_r}{|\mathbf{x} - \mathbf{r}_r|} \exp\left(ik|\mathbf{x} - \mathbf{r}_r|\right)$$

在近轴近似下，设全息平面位于 $z = 0$，参考源位于 $(x_r, y_r, z_r)$，则：

$$|\mathbf{x} - \mathbf{r}_r| \approx z_r + \frac{(x - x_r)^2 + (y - y_r)^2}{2z_r}$$

参考光相位在全息平面上的分布为：

$$\phi_r(x, y) = k\left[z_r + \frac{(x - x_r)^2 + (y - y_r)^2}{2z_r}\right]$$

记录的全息图模式为：

$$H(\mathbf{x}) = |U_o|^2 + |U_r|^2 + 2|U_o||U_r|\cos[\phi_o(\mathbf{x}) - \phi_r(\mathbf{x})]$$

当参考光远强于物光时（$|U_r| >> |U_o|$），可简化为：

$$H(\mathbf{x}) \approx |U_r|^2[1 + 2\frac{|U_o|}{|U_r|}\cos(\phi_o - \phi_r)]$$

这种线性记录条件下，全息图的调制深度正比于物光振幅。

### 24.1.3 重建过程

用相同的参考光照明全息图，透射光场为：

$$U_{trans}(\mathbf{x}) = H(\mathbf{x}) \cdot U_r(\mathbf{x})$$

展开后得到四项：

$$U_{trans} = |U_o|^2 U_r + |U_r|^2 U_r + U_o |U_r|^2 + U_o^* U_r^2$$

各项的物理意义：
1. **零级衍射**：$(|U_o|^2 + |U_r|^2)U_r$ - 直透光
2. **+1级衍射**：$U_o |U_r|^2$ - 虚像（原始物光波）
3. **-1级衍射**：$U_o^* U_r^2$ - 实像（共轭波）

虚像项准确重现原始物光波：

$$U_{virtual}(\mathbf{x}) = U_o(\mathbf{x}) |U_r(\mathbf{x})|^2$$

实像项产生共轭波：

$$U_{real}(\mathbf{x}) = U_o^*(\mathbf{x}) U_r^2(\mathbf{x})$$

在空间上，如果物体位于 $z = z_o$，则：
- 虚像位于原物体位置 $z = z_o$
- 实像位于 $z = -z_o + 2z_r$（相对于参考源的镜像位置）

离轴配置的优势在于这三个衍射级在空间上分离，便于观察单一重建像。

### 24.1.4 体积全息与布拉格条件

对于厚全息图，需考虑体积内的布拉格衍射。记录时在介质内形成三维折射率光栅：

$$n(\mathbf{r}) = n_0 + n_1 \cos(\mathbf{K} \cdot \mathbf{r})$$

其中光栅矢量 $\mathbf{K} = \mathbf{k}_o - \mathbf{k}_r$ 决定了光栅的周期和方向。

布拉格条件要求入射光满足动量匹配：

$$\mathbf{k}_{in} + \mathbf{K} = \mathbf{k}_{out}$$

在标量形式下，对于对称几何（入射角等于衍射角），布拉格角为：

$$2d\sin\theta_B = m\lambda/n_0$$

其中 $d = 2\pi/|\mathbf{K}|$ 是光栅周期。

衍射效率由耦合波理论（Kogelnik理论）给出。对于透射型相位光栅：

$$\eta = \sin^2\left(\frac{\pi n_1 d}{\lambda \cos\theta_B}\right)$$

对于反射型相位光栅：

$$\eta = \tanh^2\left(\frac{\pi n_1 d}{\lambda \cos\theta_B}\right)$$

其中 $n_1$ 是折射率调制深度，$d$ 是全息图厚度。

角度选择性（半高全宽）为：

$$\Delta\theta = \frac{\lambda}{d\sin(2\theta_B)}$$

波长选择性为：

$$\Delta\lambda = \frac{\lambda^2}{d|\cos\theta_B|}$$

这种高选择性使体积全息图可用于波长复用和角度复用存储。

## 24.2 计算机生成全息图（CGH）

### 24.2.1 从物理到计算

计算机生成全息图通过数值计算模拟物光波的传播和干涉过程。CGH的核心优势在于：
1. 可以生成物理上不存在的物体的全息图
2. 精确控制光波的振幅和相位分布
3. 无需物理干涉系统，避免了振动和相干性要求

对于离散化的3D场景，物光波可表示为点源叠加：

$$U_o(\mathbf{x}) = \sum_{j=1}^N A_j \frac{\exp(ik|\mathbf{x} - \mathbf{r}_j|)}{|\mathbf{x} - \mathbf{r}_j|}$$

其中 $\mathbf{r}_j$ 是第 $j$ 个点源的位置，$A_j$ 是其复振幅。

在实际计算中，需要考虑：
- **采样定理**：全息平面的采样间隔 $\Delta x$ 必须满足：
  $$\Delta x < \frac{\lambda z_{min}}{L}$$
  其中 $L$ 是物体横向尺寸，$z_{min}$ 是最近物点距离。

- **数值孔径限制**：可记录的最大空间频率：
  $$f_{max} = \frac{2NA}{\lambda}$$
  
- **量化效应**：数字表示的有限精度导致量化噪声。

### 24.2.2 点源法（Point Source Method）

最直接的CGH算法是点源叠加法。每个物点被视为球面波源，在全息平面上产生的光场为：

$$U_j(\mathbf{x}) = A_j \frac{\exp(ik r_j)}{r_j} \cdot \text{obliquity}(\theta_j)$$

其中 $r_j = |\mathbf{x} - \mathbf{r}_j|$，倾斜因子 $\text{obliquity}(\theta) = \frac{1 + \cos\theta}{2}$ 考虑了大角度传播的修正。

对于全息平面上的每个采样点 $\mathbf{x}_i$：

$$H(\mathbf{x}_i) = \left|\sum_{j=1}^N A_j \frac{\exp(ik|\mathbf{x}_i - \mathbf{r}_j|)}{|\mathbf{x}_i - \mathbf{r}_j|} + U_r(\mathbf{x}_i)\right|^2$$

计算复杂度为 $O(MN)$，其中 $M$ 是全息图像素数，$N$ 是场景点数。

**优化策略**：
1. **查找表加速**：预计算 $\exp(ikr)/r$ 对于量化的 $r$ 值
2. **并行计算**：每个全息像素独立计算，适合GPU加速
3. **自适应采样**：根据物点分布密度调整计算精度

**误差分析**：
主要误差源包括：
- 离散采样误差：$\epsilon_{sampling} \propto 1/\sqrt{N}$
- 有限孔径截断：$\epsilon_{aperture} \propto \lambda z/D$
- 数值精度：浮点运算引入的舍入误差

### 24.2.3 多边形法（Polygon-based Method）

对于多边形物体，可通过解析积分提高效率。对于三角形面片 $T$，其贡献为：

$$U_T(\mathbf{x}) = \iint_T \frac{A(\mathbf{r}')\exp(ik|\mathbf{x} - \mathbf{r}'|)}{|\mathbf{x} - \mathbf{r}'|} d\mathbf{r}'$$

对于平面三角形，可使用解析方法。设三角形顶点为 $\mathbf{v}_1, \mathbf{v}_2, \mathbf{v}_3$，使用重心坐标：

$$\mathbf{r}'(u,v) = (1-u-v)\mathbf{v}_1 + u\mathbf{v}_2 + v\mathbf{v}_3$$

积分变换为：

$$U_T(\mathbf{x}) = 2A_{triangle} \int_0^1 \int_0^{1-u} \frac{A(u,v)\exp(ikr(u,v))}{r(u,v)} dv du$$

**Babinet原理优化**：
对于不透明多边形，可利用Babinet原理：

$$U_{polygon} = U_{aperture} - U_{background}$$

这将复杂形状的计算转化为简单孔径的计算。

**近似方法**：
1. **恒定相位近似**：当多边形远小于到观察点的距离时
   $$U_T \approx \frac{A_{avg} \cdot S_T \exp(ikr_c)}{r_c}$$
   其中 $r_c$ 是到多边形质心的距离，$S_T$ 是面积。

2. **Fresnel近似**：在近轴条件下使用二次相位展开
   $$r \approx z + \frac{(x-x')^2 + (y-y')^2}{2z}$$

### 24.2.4 波前记录平面法（Wavefront Recording Plane）

WRP方法通过引入中间虚拟平面，将3D问题分解为多个2D传播问题，显著减少计算量。

**基本原理**：
1. 在物体附近放置多个WRP
2. 计算物点到最近WRP的短距离传播
3. 计算WRP到全息平面的长距离传播

**算法步骤**：

1. **物体到WRP的传播**（使用Rayleigh-Sommerfeld积分）：
   $$U_{WRP}(\mathbf{u}) = \sum_{j} A_j \frac{\exp(ik\rho_j)}{\rho_j}$$
   其中 $\rho_j = |\mathbf{u} - \mathbf{r}_j|$ 是物点到WRP的距离。

2. **WRP到全息平面的传播**（使用角谱方法）：
   $$U_h(\mathbf{x}) = \mathcal{F}^{-1}\{\mathcal{F}\{U_{WRP}\} \cdot H(f_x, f_y)\}$$
   
   传播传递函数：
   $$H(f_x, f_y) = \exp\left(ikd\sqrt{1 - \lambda^2(f_x^2 + f_y^2)}\right)$$
   
   对于大传播距离，可使用Fresnel近似：
   $$H(f_x, f_y) \approx \exp(ikd)\exp\left(-i\pi\lambda d(f_x^2 + f_y^2)\right)$$

**优化考虑**：
- WRP位置选择：通常放置在物体的包围盒表面
- WRP分辨率：由物体细节和传播距离决定
- 多WRP策略：对复杂物体使用多个WRP，每个负责一部分物点

**计算复杂度**：
从 $O(MN)$ 降低到 $O(M\log M + KN)$，其中 $K$ 是WRP像素数，通常 $K << M$。

### 24.2.5 层析法（Layer-based Method）

层析法将3D场景沿深度方向切片，特别适合体积数据和半透明物体的全息计算。

**基本原理**：
将3D场景分解为多个深度层，每层独立计算后叠加：

$$U_o(\mathbf{x}) = \sum_{l=1}^L U_l(\mathbf{x}) * h_{z_l}(\mathbf{x})$$

其中 $h_{z_l}$ 是从深度 $z_l$ 到全息平面的传播核。

**Fresnel传播核**：
$$h_z(\mathbf{x}) = \frac{\exp(ikz)}{i\lambda z}\exp\left(\frac{ik|\mathbf{x}|^2}{2z}\right)$$

**频域实现**（更高效）：
$$U_o = \sum_{l=1}^L \mathcal{F}^{-1}\{\mathcal{F}\{U_l\} \cdot H_l\}$$

其中 $H_l(f_x, f_y) = \exp(ikz_l)\exp(-i\pi\lambda z_l(f_x^2 + f_y^2))$

**层间距选择**：
根据采样定理，层间距应满足：
$$\Delta z \leq \frac{\lambda}{2(NA)^2}$$

其中 NA 是系统数值孔径。这确保了轴向分辨率。

**优化策略**：
1. **非均匀层分布**：在物体密集区域使用更多层
2. **自适应层数**：根据场景复杂度动态调整
3. **层间插值**：使用三线性插值减少所需层数

**遮挡处理**：
对于不透明物体，需要考虑遮挡：
$$U_l(\mathbf{x}) = A_l(\mathbf{x}) \cdot V_l(\mathbf{x})$$

其中 $V_l(\mathbf{x})$ 是可见性函数，可通过深度缓冲或光线投射计算。

**计算复杂度**：
$O(LM\log M)$，其中 $L$ 是层数，利用FFT加速每层的传播计算。

## 24.3 空间光调制器显示技术

### 24.3.1 SLM的工作原理

空间光调制器（SLM）是实现动态全息显示的关键器件。主要类型包括：

1. **液晶SLM（LC-SLM）**：通过电场控制液晶分子取向改变折射率
2. **数字微镜器件（DMD）**：通过微镜阵列的机械偏转调制光
3. **硅基液晶（LCoS）**：结合液晶和CMOS技术

对于相位型SLM，其传递函数为：

$$t_{SLM}(\mathbf{x}) = \exp[i\phi_{SLM}(\mathbf{x})]$$

其中 $\phi_{SLM} \in [0, 2\pi]$ 是可控相位延迟。

### 24.3.2 像素化效应与衍射

SLM的像素结构导致额外的衍射效应。对于像素间距 $p$，衍射角由光栅方程给出：

$$\sin\theta_m = m\frac{\lambda}{p}$$

有效视场角（FOV）受限于：

$$\text{FOV} = 2\arcsin\left(\frac{\lambda}{2p}\right)$$

### 24.3.3 振幅与相位调制

纯相位SLM无法直接实现复数调制。常用编码方法包括：

1. **双相位编码**：
   $$A\exp(i\phi) = \frac{1}{2}[\exp(i\phi_1) + \exp(i\phi_2)]$$
   其中 $\phi_1 = \phi + \arccos(A)$，$\phi_2 = \phi - \arccos(A)$

2. **误差扩散法**：
   将复数值量化到最近的可实现相位值，并将误差扩散到邻近像素

### 24.3.4 时分复用与空分复用

提高显示质量的复用技术：

1. **时分复用**：快速切换多个全息图
   $$H_{avg} = \frac{1}{T}\sum_{t=1}^T H_t$$

2. **空分复用**：将SLM分割为多个子区域
   $$H_{total}(\mathbf{x}) = \sum_{k} W_k(\mathbf{x}) H_k(\mathbf{x})$$
   
其中 $W_k$ 是窗函数。

## 24.4 相位恢复算法

### 24.4.1 相位恢复问题的数学表述

相位恢复是从强度测量中重建复数光场的逆问题。给定目标强度分布 $I_{target}(\mathbf{x}) = |U_{target}(\mathbf{x})|^2$，求解相位 $\phi(\mathbf{x})$ 使得：

$$U(\mathbf{x}) = \sqrt{I_{target}(\mathbf{x})} \exp[i\phi(\mathbf{x})]$$

这是一个非凸优化问题，存在多个局部最优解。

### 24.4.2 Gerchberg-Saxton算法

最经典的迭代相位恢复算法：

1. 初始化随机相位：$\phi_0(\mathbf{x}) = \text{random}[0, 2\pi]$
2. 迭代过程：
   - 前向传播：$U_{far}^{(k)} = \mathcal{F}\{A_{near}\exp(i\phi_{near}^{(k)})\}$
   - 施加远场约束：$U_{far}^{(k+1)} = A_{far}\exp(i\arg[U_{far}^{(k)}])$
   - 逆向传播：$U_{near}^{(k+1)} = \mathcal{F}^{-1}\{U_{far}^{(k+1)}\}$
   - 施加近场约束：$\phi_{near}^{(k+1)} = \arg[U_{near}^{(k+1)}]$

收敛条件：$\|A_{far} - |U_{far}^{(k)}|\|_2 < \epsilon$

### 24.4.3 加权Gerchberg-Saxton算法

引入权重因子改善收敛性：

$$U_{far}^{(k+1)} = w \cdot A_{far}\exp(i\arg[U_{far}^{(k)}]) + (1-w) \cdot U_{far}^{(k)}$$

其中 $w \in [0,1]$ 是混合权重。

### 24.4.4 梯度下降法

定义损失函数：

$$L = \int ||\mathcal{F}\{A_{near}\exp(i\phi)\}|^2 - I_{far}|^2 d\mathbf{x}$$

梯度更新：

$$\phi^{(k+1)} = \phi^{(k)} - \alpha \nabla_\phi L$$

其中梯度通过自动微分或解析推导获得。

### 24.4.5 相位多样性方法

使用多个测量约束提高重建质量：

$$L = \sum_{j=1}^J \||\mathcal{P}_j\{U\}|^2 - I_j\|^2$$

其中 $\mathcal{P}_j$ 是不同的传播算子（如不同距离或波长）。

## 24.5 全息渲染管线

### 24.5.1 从传统渲染到全息渲染

将传统图形管线扩展到全息领域需要考虑波动性质。全息渲染管线的主要阶段：

1. **几何处理**：与传统管线相同
2. **光波计算**：将几何转换为复数光场
3. **传播模拟**：计算光波传播
4. **全息编码**：生成SLM驱动信号

### 24.5.2 体积渲染方程的全息形式

将体积渲染方程扩展到复数域：

$$U(\mathbf{x}) = \int_V \sigma(\mathbf{r}) A(\mathbf{r}) \frac{\exp(ik|\mathbf{x} - \mathbf{r}|)}{|\mathbf{x} - \mathbf{r}|} d\mathbf{r}$$

其中 $\sigma(\mathbf{r})$ 是体密度，$A(\mathbf{r})$ 是复振幅。这与第3章的统一体积渲染方程形式一致，但在复数域工作。

### 24.5.3 加速结构与优化

1. **八叉树加速**：
   $$U(\mathbf{x}) = \sum_{node} U_{node}(\mathbf{x}) \cdot \mathbb{1}_{visible}(node)$$

2. **层次细节（LOD）**：
   根据观察距离选择不同分辨率：
   $$N_{samples} = \min\left(\frac{c}{\Delta\theta \cdot d}, N_{max}\right)$$
   
   其中 $\Delta\theta$ 是角分辨率，$d$ 是距离。

3. **GPU并行化**：
   利用FFT的并行性和光波传播的独立性

### 24.5.4 实时全息渲染

实现实时性能的关键技术：

1. **查找表方法**：
   预计算传播核：
   $$H_{LUT}[i,j] = \frac{\exp(ikr_{ij})}{r_{ij}}$$

2. **稀疏表示**：
   只计算显著贡献的点：
   $$U(\mathbf{x}) \approx \sum_{j \in S} A_j H_{LUT}[\mathbf{x}, \mathbf{r}_j]$$
   
   其中 $S = \{j : |A_j| > \epsilon\}$

3. **时间相干性利用**：
   $$U_t(\mathbf{x}) = \alpha U_{t-1}(\mathbf{x}) + (1-\alpha)\Delta U_t(\mathbf{x})$$

### 24.5.5 质量评估指标

全息重建质量的定量评估：

1. **信噪比（SNR）**：
   $$\text{SNR} = 10\log_{10}\frac{\sum|U_{target}|^2}{\sum|U_{recon} - U_{target}|^2}$$

2. **结构相似性（SSIM）**：
   应用于强度和相位分布

3. **斑点对比度**：
   $$C = \frac{\sigma_I}{\langle I \rangle}$$
   
   衡量相干噪声水平。

## 本章小结

本章建立了从物理全息到计算全息的完整理论框架：

1. **全息原理**：通过干涉记录振幅和相位，通过衍射重建原始光场
2. **CGH算法**：点源法、多边形法、WRP法和层析法，各有不同的效率-质量权衡
3. **SLM技术**：理解像素化、调制限制和复用技术对显示质量的影响
4. **相位恢复**：从强度约束反演相位的迭代算法
5. **渲染集成**：将全息计算纳入统一的体积渲染框架

关键数学工具：
- 菲涅尔-基尔霍夫衍射积分
- 快速傅里叶变换（FFT）
- 非凸优化与相位恢复
- 复数域的体积渲染方程

## 练习题

### 基础题

**24.1** 推导菲涅尔全息图的重建过程，说明为什么会产生孪生像。

<details>
<summary>提示</summary>
考虑 $H \cdot U_r$ 展开后的四项，分析每项的物理意义。
</details>

<details>
<summary>答案</summary>
展开 $H \cdot U_r = (|U_o|^2 + |U_r|^2)U_r + U_o|U_r|^2 + U_o^*U_r^2$。第三项重现原始物光波，第四项产生共轭波（孪生像），位于参考光源的镜像位置。
</details>

**24.2** 给定SLM像素间距 $p = 8\mu m$，波长 $\lambda = 633nm$，计算最大衍射角和视场角。

<details>
<summary>提示</summary>
使用光栅方程和FOV公式。
</details>

<details>
<summary>答案</summary>
最大衍射角：$\theta_{max} = \arcsin(\lambda/p) = \arcsin(633×10^{-9}/8×10^{-6}) = 4.54°$。
视场角：$FOV = 2\theta_{max} = 9.08°$。
</details>

**24.3** 证明Gerchberg-Saxton算法每次迭代不会增加误差。

<details>
<summary>提示</summary>
考虑投影算子的性质。
</details>

<details>
<summary>答案</summary>
G-S算法在近场和远场约束集之间交替投影。投影算子是非扩张的：$\|P(x) - P(y)\| \leq \|x - y\|$。因此误差单调递减。
</details>

### 挑战题

**24.4** 设计一个自适应采样算法for CGH，根据局部相位梯度调整采样密度。

<details>
<summary>提示</summary>
相位变化率与所需采样率相关。考虑奈奎斯特采样定理。
</details>

<details>
<summary>答案</summary>
局部空间频率 $f_{local} = \frac{1}{2\pi}|\nabla\phi|$。采样间隔应满足 $\Delta x < \frac{1}{2f_{local}}$。实现四叉树结构，当 $|\nabla\phi| > \frac{\pi}{\Delta x}$ 时细分。
</details>

**24.5** 推导使用两个正交偏振态同时编码两个独立全息图的方法。

<details>
<summary>提示</summary>
利用偏振的正交性和琼斯矢量表示。
</details>

<details>
<summary>答案</summary>
设两个全息图为 $H_1$, $H_2$，编码为：
$\mathbf{E} = H_1\hat{\mathbf{x}} + H_2\hat{\mathbf{y}}$。
使用偏振分束器分离：$I_x = |\hat{\mathbf{x}} \cdot \mathbf{E}|^2 = |H_1|^2$，$I_y = |\hat{\mathbf{y}} \cdot \mathbf{E}|^2 = |H_2|^2$。
</details>

**24.6** 分析层析CGH方法中层数 $L$ 与重建质量的关系，给出最优层数的估计。

<details>
<summary>提示</summary>
考虑深度分辨率、计算复杂度和衍射效应。
</details>

<details>
<summary>答案</summary>
层间距应小于景深：$\Delta z < \frac{\lambda}{NA^2}$。对于深度范围 $D$，最优层数 $L_{opt} = \frac{D \cdot NA^2}{\lambda}$。过多层数增加计算但不改善质量，因为衍射模糊了细节。
</details>

## 常见陷阱与错误

1. **混淆强度全息与相位全息**
   - 错误：认为SLM可以直接显示强度全息图
   - 正确：需要编码方法将强度信息转换为相位调制

2. **忽略采样要求**
   - 错误：使用不足的采样率导致混叠
   - 正确：确保采样满足奈奎斯特准则，特别是在高NA系统中

3. **相位恢复的初值选择**
   - 错误：总是使用随机初始相位
   - 正确：利用先验知识（如光滑性）选择更好的初值

4. **忽略SLM的物理限制**
   - 错误：假设理想的相位调制范围和分辨率
   - 正确：考虑量化误差、死区和串扰效应

5. **计算效率问题**
   - 错误：直接计算所有点对的贡献
   - 正确：使用FFT、查找表和自适应采样

## 最佳实践检查清单

### 算法选择
- [ ] 根据场景特性选择CGH算法（稀疏/密集）
- [ ] 评估实时性要求vs质量要求
- [ ] 考虑可用的硬件加速（GPU/FPGA）

### 系统设计
- [ ] 匹配SLM分辨率与目标应用
- [ ] 优化照明光源的相干性
- [ ] 设计合适的光学系统（f数、视场）

### 质量优化
- [ ] 实施迭代相位恢复提高图像质量
- [ ] 使用多重约束减少斑点噪声
- [ ] 应用预补偿校正系统像差

### 性能优化
- [ ] 利用对称性减少计算
- [ ] 实现多分辨率/LOD策略
- [ ] 使用时间相干性加速动态场景

### 验证测试
- [ ] 定量评估重建质量（SNR, SSIM）
- [ ] 测试不同观察角度和距离
- [ ] 验证实时性能指标

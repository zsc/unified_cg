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

### 24.1.2 菲涅尔全息图

对于菲涅尔全息，参考光为球面波。设参考点源位于 $\mathbf{r}_r$，则：

$$U_r(\mathbf{x}) = \frac{A_r}{|\mathbf{x} - \mathbf{r}_r|} \exp\left(ik|\mathbf{x} - \mathbf{r}_r|\right)$$

在近轴近似下，记录的全息图模式为：

$$H(\mathbf{x}) = |U_o|^2 + |U_r|^2 + 2|U_o||U_r|\cos[\phi_o(\mathbf{x}) - \phi_r(\mathbf{x})]$$

### 24.1.3 重建过程

用相同的参考光照明全息图，透射光场为：

$$U_{trans}(\mathbf{x}) = H(\mathbf{x}) \cdot U_r(\mathbf{x})$$

展开后得到四项，其中第三项重现原始物光波：

$$U_{recon}(\mathbf{x}) \propto U_o(\mathbf{x}) |U_r(\mathbf{x})|^2$$

### 24.1.4 体积全息与布拉格条件

对于厚全息图，需考虑体积内的布拉格衍射。光栅矢量 $\mathbf{K} = \mathbf{k}_o - \mathbf{k}_r$，布拉格条件为：

$$\mathbf{k}_{in} + \mathbf{K} = \mathbf{k}_{out}$$

衍射效率由耦合波理论给出：

$$\eta = \sin^2\left(\frac{\pi n_1 d}{\lambda \cos\theta_B}\right)$$

其中 $n_1$ 是折射率调制深度，$d$ 是全息图厚度。

## 24.2 计算机生成全息图（CGH）

### 24.2.1 从物理到计算

计算机生成全息图通过数值计算模拟物光波的传播和干涉过程。对于离散化的3D场景，物光波可表示为：

$$U_o(\mathbf{x}) = \sum_{j=1}^N A_j \frac{\exp(ik|\mathbf{x} - \mathbf{r}_j|)}{|\mathbf{x} - \mathbf{r}_j|}$$

其中 $\mathbf{r}_j$ 是第 $j$ 个点源的位置，$A_j$ 是其复振幅。

### 24.2.2 点源法（Point Source Method）

最直接的CGH算法是点源叠加法。对于全息平面上的每个采样点 $\mathbf{x}_i$：

$$H(\mathbf{x}_i) = \left|\sum_{j=1}^N A_j \frac{\exp(ik|\mathbf{x}_i - \mathbf{r}_j|)}{|\mathbf{x}_i - \mathbf{r}_j|} + U_r(\mathbf{x}_i)\right|^2$$

计算复杂度为 $O(MN)$，其中 $M$ 是全息图像素数，$N$ 是场景点数。

### 24.2.3 多边形法（Polygon-based Method）

对于多边形物体，可通过解析积分提高效率。对于三角形面片 $T$，其贡献为：

$$U_T(\mathbf{x}) = \iint_T \frac{A(\mathbf{r}')\exp(ik|\mathbf{x} - \mathbf{r}'|)}{|\mathbf{x} - \mathbf{r}'|} d\mathbf{r}'$$

使用解析近似或数值积分方法计算此积分。

### 24.2.4 波前记录平面法（Wavefront Recording Plane）

引入中间波前记录平面（WRP）减少计算量：

1. 计算物体到WRP的传播：
   $$U_{WRP}(\mathbf{u}) = \mathcal{F}\{U_{obj}\} \exp(ikz\sqrt{1 - \lambda^2(f_x^2 + f_y^2)})$$

2. 从WRP到全息平面的传播：
   $$U_h(\mathbf{x}) = \mathcal{F}^{-1}\{\mathcal{F}\{U_{WRP}\} \exp(ikd\sqrt{1 - \lambda^2(f_x^2 + f_y^2)})\}$$

### 24.2.5 层析法（Layer-based Method）

将3D场景分解为多个深度层，每层独立计算后叠加：

$$U_o(\mathbf{x}) = \sum_{l=1}^L U_l(\mathbf{x}) * h_{z_l}(\mathbf{x})$$

其中 $h_{z_l}$ 是从深度 $z_l$ 到全息平面的传播核：

$$h_z(\mathbf{x}) = \frac{\exp(ikz)}{i\lambda z}\exp\left(\frac{ik|\mathbf{x}|^2}{2z}\right)$$

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

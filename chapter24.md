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

**干涉条纹的形成机理**：
当物光波和参考光波在记录介质上相遇时，空间各点的相位差决定了干涉强度。对于平面参考波 $U_r = A_r \exp(i\mathbf{k}_r \cdot \mathbf{x})$ 和来自点源的球面物波：

$$U_o = \frac{A_o}{r} \exp(ikr)$$

相位差为：
$$\Delta\phi = kr - \mathbf{k}_r \cdot \mathbf{x}$$

等相位面（亮条纹位置）满足 $\Delta\phi = 2\pi m$，形成三维驻波模式。

**记录介质的响应特性**：
全息记录材料的响应可用透过率调制度表征：

$$T(I) = T_0 + \beta I$$

其中 $\beta$ 是材料的响应系数。对于线性记录：

$$T(\mathbf{x}) = T_0 + \beta[|U_o|^2 + |U_r|^2 + 2|U_o||U_r|\cos(\phi_o - \phi_r)]$$

非线性响应会产生高阶衍射项，影响重建质量。常见记录材料包括：
- 卤化银乳胶：高灵敏度，分辨率>5000线/mm
- 光致聚合物：实时记录，动态范围大
- 光折变晶体：可擦写，适合动态全息

**空间频率分析**：
全息图可视为物光波的空间频谱载波调制。设物光波的频谱为 $\tilde{U}_o(f_x, f_y)$，参考光引入载波频率 $(f_{xr}, f_{yr})$：

$$\tilde{I}(f_x, f_y) = \tilde{U}_o \otimes \tilde{U}_o^* + \delta(f_x, f_y) + \tilde{U}_o(f_x - f_{xr}, f_y - f_{yr}) + \tilde{U}_o^*(-(f_x + f_{xr}), -(f_y + f_{yr}))$$

这解释了为什么离轴配置可以空间分离不同衍射级。载波频率必须满足：

$$f_{xr} > 3f_{x,max}, \quad f_{yr} > 3f_{y,max}$$

以避免频谱重叠。

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

将全息图透过率函数 $H(\mathbf{x}) = T_0 + \beta I(\mathbf{x})$ 代入：

$$U_{trans} = [T_0 + \beta(|U_o|^2 + |U_r|^2 + U_o U_r^* + U_o^* U_r)]U_r$$

展开后得到多项：

$$U_{trans} = T_0 U_r + \beta(|U_o|^2 + |U_r|^2)U_r + \beta U_o |U_r|^2 + \beta U_o^* U_r^2$$

各项的物理意义和空间分布：
1. **零级衍射**：$[T_0 + \beta(|U_o|^2 + |U_r|^2)]U_r$ - 直透光，沿参考光方向传播
2. **+1级衍射**：$\beta U_o |U_r|^2$ - 虚像（原始物光波的准确重现）
3. **-1级衍射**：$\beta U_o^* U_r^2$ - 实像（共轭波，产生赝像）

**虚像的形成机理**：
虚像项 $U_{virtual}(\mathbf{x}) = \beta U_o(\mathbf{x}) |U_r(\mathbf{x})|^2$ 准确重现原始物光波。当 $|U_r|$ 为常数（平面参考波）时：

$$U_{virtual} \propto U_o(\mathbf{x})$$

这解释了为什么虚像保持原始物体的所有光学特性，包括深度信息和视差。

**实像的共轭性质**：
实像项 $U_{real}(\mathbf{x}) = \beta U_o^*(\mathbf{x}) U_r^2(\mathbf{x})$ 产生相位共轭波。对于球面波物光：

$$U_o = \frac{A_o}{r_o} \exp(ikr_o) \Rightarrow U_o^* = \frac{A_o^*}{r_o} \exp(-ikr_o)$$

共轭波向内汇聚而非向外发散，形成实像。空间位置关系：
- 虚像位于原物体位置 $z = z_o$（保持原始深度）
- 实像位于 $z = -z_o + 2z_r$（相对于参考源的镜像位置）

**角度分离条件**：
离轴配置使不同衍射级在角度上分离。设参考光入射角为 $\theta_r$，物光平均出射角为 $\theta_o$，则：
- 零级：沿 $\theta_r$ 方向
- +1级：沿 $\theta_o$ 方向
- -1级：沿 $2\theta_r - \theta_o$ 方向

为避免重叠，需要：
$$|\theta_o - \theta_r| > \Delta\theta_{obj} + \Delta\theta_{ref}$$

其中 $\Delta\theta$ 是光束发散角。

**衍射效率分析**：
对于振幅全息图，一级衍射效率定义为：

$$\eta = \frac{|U_{+1}|^2}{|U_{incident}|^2} = \frac{|\beta U_o |U_r|^2|^2}{|U_r|^2}$$

在线性记录条件下（$|U_o| << |U_r|$），且透过率调制度 $m = \beta |U_o||U_r|/T_0$：

$$\eta = \left(\frac{m T_0}{2}\right)^2$$

对于理想的正弦光栅（$m = 1$），当 $T_0 = 0.5$ 时：
$$\eta_{max} = \left(\frac{1 \times 0.5}{2}\right)^2 = 0.0625 = 6.25\%$$

**厚全息图的耦合波分析**：
对于体积全息图，使用Kogelnik的耦合波理论。定义耦合参数：
$$\nu = \frac{\pi n_1 d}{\lambda \cos\theta_B}$$

透射型相位光栅的衍射效率：
$$\eta = \sin^2(\nu) = \sin^2\left(\frac{\pi n_1 d}{\lambda \cos\theta_B}\right)$$

当 $\nu = \pi/2$ 时达到100%效率，需要：
$$n_1 d = \frac{\lambda \cos\theta_B}{2}$$

反射型相位光栅：
$$\eta = \tanh^2(\nu) = \tanh^2\left(\frac{\pi n_1 d}{\lambda \cos\theta_B}\right)$$

大耦合强度下趋向100%。这解释了为什么相位全息图比振幅全息图效率更高。

**像质评估**：
重建像的质量受多种因素影响，可通过系统分析量化：

1. **记录介质的MTF（调制传递函数）**：
   $$MTF(f) = \frac{M_{out}(f)}{M_{in}(f)} = \frac{|H(f)|}{|H(0)|}$$
   
   其中 $M$ 是调制度，$H(f)$ 是系统传递函数。对于典型的卤化银乳胶：
   $$MTF(f) \approx \exp\left(-\frac{f^2}{2f_c^2}\right)$$
   
   截止频率 $f_c \approx 3000-5000$ 线/mm。高频衰减导致细节损失，表现为边缘模糊。

2. **参考光的相干性要求**：
   
   **时间相干性**：相干长度必须大于最大光程差
   $$l_c = c\tau_c = \frac{\lambda^2}{\Delta\lambda} > \Delta_{max}$$
   
   其中 $\Delta_{max} = |r_{o,max} - r_{o,min}| + |r_{r,max} - r_{r,min}|$。
   
   对于He-Ne激光器（$\lambda = 633nm$, $\Delta\lambda \approx 0.002nm$）：
   $$l_c \approx \frac{(633 \times 10^{-9})^2}{0.002 \times 10^{-9}} \approx 20cm$$
   
   **空间相干性**：相干宽度决定干涉条纹可见度
   $$w_c = \frac{\lambda D}{s}$$
   
   其中 $D$ 是到光源的距离，$s$ 是光源尺寸。条纹可见度：
   $$V = \frac{I_{max} - I_{min}}{I_{max} + I_{min}} = |\gamma_{12}|$$
   
   其中 $\gamma_{12}$ 是复相干度。对于扩展光源：
   $$\gamma_{12} = \frac{\sin(\pi s \sin\theta/\lambda)}{\pi s \sin\theta/\lambda}$$

3. **记录几何的像差分析**：
   球面参考波引入的波前像差可展开为Zernike多项式：
   $$W(\rho, \theta) = \sum_{n=0}^{\infty} \sum_{m=-n}^{n} a_{nm} Z_n^m(\rho, \theta)$$
   
   主要像差项：
   - 离焦（$n=2, m=0$）：$a_{20} = \frac{\Delta z}{8(f/\#)^2}$
   - 球差（$n=4, m=0$）：$a_{40} = \frac{1}{48(f/\#)^4}$
   - 彗差（$n=3, m=1$）：$a_{31} = \frac{\theta_{off}}{8(f/\#)^3}$
   
   其中 $f/\#$ 是F数，$\theta_{off}$ 是离轴角。
   
   **Strehl比**评估像质：
   $$S = \exp\left(-\left(\frac{2\pi\sigma_W}{\lambda}\right)^2\right)$$
   
   其中 $\sigma_W$ 是波前误差的RMS值。$S > 0.8$ 认为是衍射极限成像。

**波长选择性与角度选择性**：
重建时使用不同波长 $\lambda'$ 或不同角度 $\theta'$ 会导致像质变化：

1. **波长变化的影响**：
   
   **横向放大率**：根据光栅方程
   $$\sin\theta_d = \sin\theta_r + \frac{\lambda'}{\Lambda}$$
   
   其中 $\Lambda$ 是干涉条纹间距。横向放大率：
   $$M_x = \frac{\lambda'}{\lambda} \cdot \frac{\cos\theta_r}{\cos\theta_d}$$
   
   **轴向位置变化**：
   $$z_i' = z_i \cdot \frac{\lambda}{\lambda'}$$
   
   **色差**：不同波长成像在不同深度，产生色散：
   $$\Delta z_{chromatic} = z_i \left(1 - \frac{\lambda_{min}}{\lambda_{max}}\right)$$

2. **角度变化的影响**：
   
   **布拉格失配**：偏离布拉格角 $\Delta\theta$ 导致效率下降
   $$\eta(\Delta\theta) = \eta_0 \cdot \text{sinc}^2\left(\frac{\pi d \Delta\theta}{\Lambda}\right)$$
   
   **像差引入**：主要是球差和彗差
   $$W_{spherical} = \frac{(\Delta\theta)^2 \rho^4}{32(f/\#)^3}$$
   $$W_{coma} = \frac{\Delta\theta \cdot \rho^3 \cos\phi}{8(f/\#)^2}$$
   
   其中 $(\rho, \phi)$ 是归一化极坐标。

3. **补偿方法**：
   
   **预畸变**：记录时引入相反的畸变
   $$\phi_{comp}(x,y) = -\frac{2\pi}{\lambda}\left(\frac{\lambda'}{\lambda} - 1\right)\sqrt{x^2 + y^2 + z^2}$$
   
   **动态补偿**：使用SLM实时校正波前误差

**数字重建方法**：
除了光学重建，可通过数值计算模拟重建过程，这在数字全息显微镜中尤其重要：

$$U_{recon}(\xi, \eta) = \iint H(x,y) U_r(x,y) h(x,y;\xi,\eta) dx dy$$

其中 $h$ 是从全息平面到重建平面的传播核。根据Rayleigh-Sommerfeld衍射理论：

$$h(x,y;\xi,\eta) = \frac{1}{i\lambda} \frac{\exp(ikr)}{r} \cos\theta$$

其中 $r = \sqrt{(\xi-x)^2 + (\eta-y)^2 + z^2}$，$\cos\theta = z/r$。

**频域快速算法**：
利用卷积定理和FFT加速计算：

$$U_{recon} = \mathcal{F}^{-1}\{\mathcal{F}\{H \cdot U_r\} \cdot \mathcal{F}\{h\}\}$$

对于Fresnel近似，传播传递函数简化为：
$$H_z(f_x, f_y) = \exp(ikz)\exp\left[-i\pi\lambda z(f_x^2 + f_y^2)\right]$$

**数字聚焦技术**：
通过改变传播距离 $z$ 实现后聚焦：
$$U(x,y,z) = \mathcal{F}^{-1}\{\mathcal{F}\{U(x,y,0)\} \cdot H_z(f_x,f_y)\}$$

搜索最佳聚焦平面的判据：
- 强度梯度最大化：$\max_z \sum|\nabla I(x,y,z)|^2$
- 频谱能量集中度：$\max_z \frac{\sum f^2|\tilde{I}(f)|^2}{\sum|\tilde{I}(f)|^2}$
- Tamura系数：$\max_z \frac{\sigma_I^2}{\mu_I}$

**数字像差补偿**：
测量或计算系统像差后，引入补偿相位：
$$U_{corrected} = U_{recon} \cdot \exp[-i\phi_{aberration}]$$

常见补偿项：
- 倾斜：$\phi_{tilt} = k(x\sin\theta_x + y\sin\theta_y)$
- 离焦：$\phi_{defocus} = \frac{k(x^2+y^2)}{2R}$
- 球差：$\phi_{spherical} = \frac{k(x^2+y^2)^2}{8R^3}$

### 24.1.4 体积全息与布拉格条件

**三维光栅的形成**：
对于厚全息图（厚度 $d >> \Lambda$，其中 $\Lambda$ 是条纹间距），需考虑体积内的布拉格衍射。记录时，干涉图样在整个体积内形成三维强度分布：

$$I(\mathbf{r}) = |U_o(\mathbf{r}) + U_r(\mathbf{r})|^2$$

对于线性记录材料，折射率调制正比于曝光强度：

$$n(\mathbf{r}) = n_0 + n_1 \cos(\mathbf{K} \cdot \mathbf{r} + \phi_0)$$

其中光栅矢量 $\mathbf{K} = \mathbf{k}_o - \mathbf{k}_r$ 决定了光栅的周期和方向。光栅周期：

$$\Lambda = \frac{2\pi}{|\mathbf{K}|} = \frac{\lambda}{2\sin(\theta/2)}$$

其中 $\theta$ 是物光和参考光的夹角。

**体积光栅的类型**：
1. **透射型光栅**：$\mathbf{K} \perp$ 表面，用于透射几何
2. **反射型光栅**：$\mathbf{K} \parallel$ 表面，用于反射几何
3. **倾斜光栅**：$\mathbf{K}$ 与表面成任意角，混合特性

布拉格条件要求入射光满足动量匹配：

$$\mathbf{k}_{in} + \mathbf{K} = \mathbf{k}_{out}$$

在标量形式下，对于对称几何（入射角等于衍射角），布拉格角为：

$$2d\sin\theta_B = m\lambda/n_0$$

其中 $d = 2\pi/|\mathbf{K}|$ 是光栅周期。

**耦合波理论分析**：
Kogelnik的耦合波理论提供了体积光栅衍射效率的解析解。定义耦合常数：

$$\kappa = \frac{\pi n_1}{\lambda \cos\theta_B}$$

和失谐参数（偏离布拉格条件）：

$$\xi = \frac{\Delta\theta \cdot K d}{2\cos\theta_B}$$

**透射型相位光栅**：
衍射效率为：
$$\eta = \frac{\sin^2(\sqrt{\nu^2 + \xi^2})}{1 + \xi^2/\nu^2}$$

其中 $\nu = \kappa d = \frac{\pi n_1 d}{\lambda \cos\theta_B}$。布拉格条件下（$\xi = 0$）：

$$\eta = \sin^2(\nu) = \sin^2\left(\frac{\pi n_1 d}{\lambda \cos\theta_B}\right)$$

最大效率100%出现在 $\nu = \pi/2$。

**反射型相位光栅**：
衍射效率为：
$$\eta = \frac{\sinh^2(\sqrt{s^2 - \xi^2})}{1 + s^2/\xi^2}$$

其中 $s = \kappa d$。布拉格条件下：

$$\eta = \tanh^2(s) = \tanh^2\left(\frac{\pi n_1 d}{\lambda \cos\theta_B}\right)$$

大耦合强度下趋向100%。

**Q参数与光栅分类**：
Klein-Cook参数区分薄光栅和厚光栅：

$$Q = \frac{2\pi\lambda d}{n_0\Lambda^2}$$

- $Q < 1$：Raman-Nath衍射（薄光栅），多级衍射
- $Q > 10$：Bragg衍射（厚光栅），单级衍射
- $1 < Q < 10$：过渡区域

**选择性分析**：

**角度选择性**（半高全宽）：
对于透射光栅：
$$\Delta\theta = \frac{\lambda\cos\theta_B}{\pi d\sin\theta_B} \approx \frac{\Lambda}{d}$$

对于反射光栅：
$$\Delta\theta = \frac{\lambda}{2d\sin\theta_B}$$

角度选择性与光栅厚度成反比，厚光栅具有更高的角度选择性。

**波长选择性**：
布拉格波长偏移的容忍度：
$$\Delta\lambda = \frac{\lambda^2}{2nd\sin\theta_B}$$

对于垂直入射的反射光栅：
$$\Delta\lambda = \frac{\lambda^2}{2nd}$$

**温度和应力效应**：
温度变化引起的波长偏移：
$$\Delta\lambda_T = \lambda \left(\alpha + \frac{1}{n}\frac{dn}{dT}\right) \Delta T$$

其中 $\alpha$ 是热膨胀系数。应力引起的双折射：
$$\Delta n = C_{ij} \sigma_{ij}$$

其中 $C_{ij}$ 是光弹系数张量。

**多重全息与复用技术**：
利用选择性可在同一体积内记录多个全息图：

1. **角度复用**：不同角度记录 $N$ 个全息图
   $$N_{angle} \approx \frac{\theta_{range}}{\Delta\theta} = \frac{\theta_{range} \cdot d}{\Lambda}$$

2. **波长复用**：不同波长记录
   $$N_{wavelength} \approx \frac{\Delta\lambda_{source}}{\Delta\lambda} = \frac{\Delta\lambda_{source} \cdot 2nd}{\lambda^2}$$

3. **空间复用**：分区记录
4. **相位编码复用**：使用正交相位码

总存储容量：
$$C = N_{angle} \times N_{wavelength} \times N_{spatial} \times \frac{A}{(\Delta x)^2}$$

其中 $A$ 是全息图面积，$\Delta x$ 是空间分辨率。

## 24.2 计算机生成全息图（CGH）

### 24.2.1 从物理到计算

计算机生成全息图通过数值计算模拟物光波的传播和干涉过程。CGH的核心优势在于：
1. 可以生成物理上不存在的物体的全息图
2. 精确控制光波的振幅和相位分布
3. 无需物理干涉系统，避免了振动和相干性要求
4. 可以引入计算优化和预补偿

**从连续到离散的表示**：
对于连续3D场景，物光波为：
$$U_o(\mathbf{x}) = \iiint_{V} \rho(\mathbf{r}') G(\mathbf{x}, \mathbf{r}') d\mathbf{r}'$$

其中 $\rho(\mathbf{r}')$ 是物体的复振幅分布，$G$ 是格林函数：
$$G(\mathbf{x}, \mathbf{r}') = \frac{\exp(ik|\mathbf{x} - \mathbf{r}'|)}{4\pi|\mathbf{x} - \mathbf{r}'|}$$

离散化后，使用点采样：
$$U_o(\mathbf{x}) \approx \sum_{j=1}^N A_j \frac{\exp(ik|\mathbf{x} - \mathbf{r}_j|)}{|\mathbf{x} - \mathbf{r}_j|}$$

其中 $A_j = \rho(\mathbf{r}_j) \Delta V$ 包含了体积元素。

**光波传播的数学框架**：
CGH计算基于标量衍射理论。根据不同近似程度：

1. **Rayleigh-Sommerfeld衍射**（最精确）：
   $$U(\mathbf{x}) = \frac{1}{i\lambda} \iint_{\Sigma} U_0(\mathbf{x}') \frac{\exp(ikr)}{r} \cos(\mathbf{n}, \mathbf{r}) d\mathbf{x}'$$

2. **Fresnel衍射**（近轴近似）：
   $$U(x,y,z) = \frac{\exp(ikz)}{i\lambda z} \iint U_0(x',y') \exp\left[\frac{ik}{2z}[(x-x')^2 + (y-y')^2]\right] dx'dy'$$

3. **Fraunhofer衍射**（远场近似）：
   $$U(x,y) = \frac{\exp(ikz)\exp\left[\frac{ik(x^2+y^2)}{2z}\right]}{i\lambda z} \mathcal{F}\{U_0\}\left(\frac{x}{\lambda z}, \frac{y}{\lambda z}\right)$$

在实际计算中，需要考虑多个约束条件：

**采样要求**：
1. **横向采样**（Nyquist准则）：
   $$\Delta x < \frac{\lambda z_{min}}{2L}$$
   其中 $L$ 是物体横向尺寸，$z_{min}$ 是最近物点距离。

2. **角谱采样**：为避免混叠
   $$\Delta u < \frac{1}{2L}, \quad \Delta v < \frac{1}{2L}$$
   其中 $(u,v) = (\sin\theta_x/\lambda, \sin\theta_y/\lambda)$ 是方向余弦。

**数值孔径与带宽限制**：
系统的空间带宽受限于：
$$B_{spatial} = \frac{2NA}{\lambda} = \frac{2\sin\theta_{max}}{\lambda}$$

对应的最小可分辨特征：
$$\delta_{min} = \frac{\lambda}{2NA}$$

**量化与噪声分析**：
1. **相位量化**：$B$位量化产生量化噪声
   $$\sigma_\phi^2 = \frac{(2\pi)^2}{12 \cdot 2^{2B}}$$
   
2. **振幅量化**：影响动态范围
   $$SNR = 20\log_{10}\left(\frac{2^B - 1}{\sigma_n}\right) \text{ dB}$$

3. **计算精度**：浮点运算误差
   $$\epsilon_{total} = \epsilon_{round} + N \cdot \epsilon_{accum}$$
   其中 $N$ 是累加次数。

**Rayleigh-Sommerfeld衍射理论基础**：
CGH计算基于标量衍射理论。第一类Rayleigh-Sommerfeld积分给出：

$$U(\mathbf{x}) = \frac{1}{i\lambda} \iint_{\Sigma} U_0(\mathbf{x}') \frac{\exp(ikr)}{r} \cos(\mathbf{n}, \mathbf{r}) d\mathbf{x}'$$

其中 $r = |\mathbf{x} - \mathbf{x}'|$，$\cos(\mathbf{n}, \mathbf{r})$ 是倾斜因子。这是CGH算法的理论基础。

**带限信号的完美重建条件**：
根据Whittaker-Shannon采样定理，要完美重建带限信号，采样率必须满足：

$$f_s > 2B$$

其中 $B$ 是信号带宽。对于全息图，空间带宽由物体的角谱范围决定：

$$B_x = \frac{NA_x}{\lambda}, \quad B_y = \frac{NA_y}{\lambda}$$

因此全息图的采样间隔应满足：

$$\Delta x < \frac{\lambda}{2NA_x}, \quad \Delta y < \frac{\lambda}{2NA_y}$$

**复数编码策略**：
由于大多数空间光调制器只能调制振幅或相位，需要编码方法表示复数值：

1. **Kinoform（纯相位编码）**：
   忽略振幅，只编码相位：
   $$H_{kino} = \exp(i\arg[U_o])$$
   
   效率损失：$\eta = |\langle U_o, H_{kino} \rangle|^2 / \|U_o\|^2$
   
2. **迭代傅里叶变换算法（IFTA）**：
   交替投影在全息平面和重建平面约束间：
   - 全息平面：$\phi_{holo}^{(k+1)} = \arg[U_{holo}^{(k)}]$
   - 重建平面：$U_{image}^{(k+1)} = A_{target} \exp(i\arg[U_{image}^{(k)}])$
   
   收敛速度：$O(1/\sqrt{k})$，其中 $k$ 是迭代次数。

3. **双相位分解**：
   $$A e^{i\phi} = \frac{1}{2}[e^{i\phi_1} + e^{i\phi_2}]$$
   
   求解约束：
   $$\cos(\phi_1 - \phi) = \cos(\phi - \phi_2) = A$$
   
   得到：$\phi_1 = \phi + \arccos(A)$，$\phi_2 = \phi - \arccos(A)$
   
   限制：要求 $0 \leq A \leq 1$。

4. **误差扩散编码**：
   将量化误差扩散到邻近像素：
   $$e_{i,j} = U_{target} - U_{quantized}$$
   $$U_{i+1,j} \leftarrow U_{i+1,j} + \alpha e_{i,j}$$
   
   其中 $\alpha \approx 7/16$ 为Floyd-Steinberg权重。

5. **张量分解方法**：
   将复振幅矩阵分解为低秩近似：
   $$\mathbf{U} \approx \sum_{r=1}^R \sigma_r \mathbf{u}_r \mathbf{v}_r^T$$
   
   每个分量可以独立编码。

**计算精度与数值稳定性**：
浮点运算的有限精度影响CGH质量，需要特殊处理：

1. **相位卷绕处理**：
   大传播距离时相位迅速增长：
   $$\phi = kr = \frac{2\pi}{\lambda}\sqrt{(x-x')^2 + (y-y')^2 + (z-z')^2}$$
   
   使用双精度计算或相位展开：
   $$\exp(ikr) = \exp(i\phi_{wrapped}) \cdot \exp(i2\pi m)$$
   其中 $\phi_{wrapped} = \mod(kr, 2\pi)$，$m = \lfloor kr/2\pi \rfloor$。

2. **近场奇异性处理**：
   当 $r \to 0$ 时，使用正则化：
   $$G_{reg}(r) = \begin{cases}
   \frac{\exp(ikr) - 1}{ikr} + 1/\epsilon & r < \epsilon \\
   \frac{\exp(ikr)}{r} & r \geq \epsilon
   \end{cases}$$
   
   其中 $\epsilon = \lambda/100$ 避免除零。

3. **FFT数值误差分析**：
   
   **舆入误差**：FFT的误差传播
   $$\|\text{FFT}(x + \delta x) - \text{FFT}(x)\| \leq \sqrt{N} \|\delta x\|$$
   
   **累积误差**：$N$点FFT的总误差
   $$\epsilon_{total} \approx \sqrt{N\log_2 N} \cdot \epsilon_{machine}$$
   
   对于 $N = 4096^2$：
   - 单精度：$\epsilon \approx 10^{-4}$
   - 双精度：$\epsilon \approx 10^{-13}$

4. **数值稳定性优化**：
   
   **预条件化**：归一化坐标和振幅
   $$\tilde{x} = x/L, \quad \tilde{U} = U/U_{max}$$
   
   **Kahan求和**：减少累加误差
   ```
   sum = 0; c = 0
   for j in 1:N
       y = a[j] - c
       t = sum + y
       c = (t - sum) - y
       sum = t
   ```
   
   **分块计算**：避免过大数组

### 24.2.2 点源法（Point Source Method）

**基本原理**：
点源法是最直观的CGH算法，基于惠更斯-菲涅尔原理。每个物点被视为次级球面波源：

$$U_j(\mathbf{x}) = A_j \frac{\exp(ik r_j)}{r_j} \cdot K(\theta_j)$$

其中：
- $r_j = |\mathbf{x} - \mathbf{r}_j| = \sqrt{(x-x_j)^2 + (y-y_j)^2 + (z-z_j)^2}$
- $K(\theta)$ 是倾斜因子，有多种选择：

**倾斜因子的选择**：
1. **Kirchhoff倾斜因子**：
   $$K_{Kirchhoff}(\theta) = \frac{1 + \cos\theta}{2}$$

2. **Rayleigh-Sommerfeld倾斜因子**：
   $$K_{RS}(\theta) = \cos\theta = \frac{z}{r}$$

3. **无倾斜因子**（近轴近似）：
   $$K_{paraxial}(\theta) = 1$$

对于大多数应用，Rayleigh-Sommerfeld因子给出最准确的结果。

**全息图计算**：
对于全息平面上的每个采样点 $\mathbf{x}_i$，计算总光场：

$$U_{total}(\mathbf{x}_i) = \sum_{j=1}^N A_j \frac{\exp(ik r_{ij})}{r_{ij}} K(\theta_{ij}) + U_r(\mathbf{x}_i)$$

其中 $r_{ij} = |\mathbf{x}_i - \mathbf{r}_j|$。

**干涉图样的记录**：
1. **强度全息图**：
   $$H_I(\mathbf{x}_i) = |U_{total}(\mathbf{x}_i)|^2$$

2. **相位全息图**（Kinoform）：
   $$H_\phi(\mathbf{x}_i) = \arg[U_{total}(\mathbf{x}_i)]$$

3. **复振幅全息图**：
   $$H_C(\mathbf{x}_i) = U_{total}(\mathbf{x}_i)$$

**计算复杂度分析**：
- 时间复杂度：$O(MN)$
  - $M = M_x \times M_y$ 是全息图像素数
  - $N$ 是场景点数
- 空间复杂度：$O(M + N)$
- 每像素计算：
  - 距离计算：9次乘法，6次加法，1次平方根
  - 相位计算：1次乘法
  - 复数运算：2次三角函数

**优化策略**：

1. **查找表（LUT）加速**：
   预计算常用函数值：
   ```
   对于 r ∈ [r_min, r_max]，步长 Δr
   LUT[i] = exp(ik*r[i])/r[i]
   ```
   
   内存需求：$N_{LUT} = (r_{max} - r_{min})/\Delta r$
   精度权衡：$\Delta r < \lambda/10$ 保证相位精度

2. **GPU并行加速**：
   ```
   每个线程计算一个全息像素
   thread(i,j):
       U = 0
       for each object point k:
           U += A[k] * G(x[i,j], r[k])
       H[i,j] = |U + U_ref|^2
   ```
   
   加速比：$S = N_{cores} \times \eta_{occupancy}$
   典型值：$S \approx 100-1000$ 对于现代GPU

3. **分层距离计算**：
   根据距离分组，使用不同精度：
   - 近场（$r < 10\lambda$）：完整计算
   - 中场（$10\lambda < r < 1000\lambda$）：Fresnel近似
   - 远场（$r > 1000\lambda$）：Fraunhofer近似

4. **空间划分优化**：
   使用八叉树或k-d树：
   - 只计算影响显著的点
   - 剪枝条件：$|A_j|/r_j^2 < \epsilon$
   - 复杂度降低到：$O(M\log N)$

5. **自适应精度控制**：
   根据局部相位梯度调整：
   $$N_{local} = \max\left(N_{min}, \frac{|\nabla\phi|}{2\pi} \cdot N_{base}\right)$$

**误差分析**：

1. **离散采样误差**：
   由于有限采样点表示连续物体：
   $$\epsilon_{sampling} = \frac{\sigma_U}{\sqrt{N}} \approx \frac{\langle|A|\rangle}{\sqrt{N}}$$
   
   其中 $\sigma_U$ 是光场标准差。

2. **有限孔径误差**：
   全息图有限尺寸导致频谱截断：
   $$\epsilon_{aperture} = \int_{|f| > f_{max}} |\tilde{U}(f)|^2 df$$
   
   对于点源，近似为：
   $$\epsilon_{aperture} \approx \frac{\lambda z}{\pi D} \cdot \frac{|A|}{r}$$
   
   其中 $D$ 是全息图孔径。

3. **量化误差**：
   数字表示的有限位数：
   $$\epsilon_{quant} = \frac{\Delta}{2\sqrt{3}}$$
   
   其中 $\Delta$ 是量化步长。

4. **总误差估计**：
   假设各误差源独立：
   $$\epsilon_{total} = \sqrt{\epsilon_{sampling}^2 + \epsilon_{aperture}^2 + \epsilon_{quant}^2 + \epsilon_{numeric}^2}$$

5. **信噪比分析**：
   $$SNR = 10\log_{10}\left(\frac{\sum|U_{signal}|^2}{\sum|U_{noise}|^2}\right)$$
   
   典型值：
   - 单精度计算：SNR ≈ 40-50 dB
   - 双精度计算：SNR ≈ 80-90 dB
   - 8位量化：SNR ≈ 48 dB

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

**液晶SLM的物理机制**：
向列型液晶的双折射特性使其折射率随分子取向变化：

$$\Delta n(\theta) = n_e \cos^2\theta + n_o \sin^2\theta - n_o$$

其中 $n_e$、$n_o$ 分别是非常光和寻常光折射率，$\theta$ 是液晶分子倾角。相位延迟为：

$$\phi = \frac{2\pi d \Delta n(\theta)}{\lambda}$$

其中 $d$ 是液晶层厚度。通过施加电压 $V$ 控制倾角：

$$\theta(V) = \begin{cases}
0 & V < V_{th} \\
\arcsin\sqrt{\frac{V^2 - V_{th}^2}{V^2}} & V > V_{th}
\end{cases}$$

**DMD的工作特性**：
DMD由微镜阵列组成，每个微镜可在 ±12° 间切换。其调制特性：

$$t_{DMD}(x,y) = \begin{cases}
1 & \text{镜子处于 ON 状态} \\
0 & \text{镜子处于 OFF 状态}
\end{cases}$$

通过脉宽调制（PWM）实现灰度：

$$t_{avg} = \frac{t_{ON}}{t_{ON} + t_{OFF}}$$

刷新率可达 22kHz，适合时分复用全息显示。

**响应时间与带宽**：
- 液晶SLM：响应时间 ~10ms，刷新率 60-120Hz
- DMD：切换时间 ~10μs，二进制帧率 >20kHz
- 铁电液晶：响应时间 ~100μs，适合中速应用

### 24.3.2 像素化效应与衍射

SLM的像素结构导致额外的衍射效应。对于像素间距 $p$，衍射角由光栅方程给出：

$$\sin\theta_m = m\frac{\lambda}{p}$$

有效视场角（FOV）受限于：

$$\text{FOV} = 2\arcsin\left(\frac{\lambda}{2p}\right)$$

**像素化的傅里叶分析**：
SLM可建模为理想光场与像素函数的卷积：

$$t_{pixelated}(x,y) = t_{ideal}(x,y) \otimes \text{rect}\left(\frac{x}{w}, \frac{y}{w}\right) * \text{comb}\left(\frac{x}{p}, \frac{y}{p}\right)$$

其中 $w$ 是像素宽度（填充因子 $FF = w/p$）。频域表示：

$$\tilde{t}_{pixelated}(f_x, f_y) = \tilde{t}_{ideal}(f_x, f_y) \cdot \text{sinc}(wf_x, wf_y) \otimes \text{comb}(pf_x, pf_y)$$

这产生了多个衍射级，强度分布为：

$$I_m = I_0 \cdot \text{sinc}^2\left(\frac{m\lambda w}{p}\right)$$

**串扰与对比度**：
相邻像素间的串扰影响调制质量：

$$\text{Crosstalk} = \frac{I_{neighbor}}{I_{pixel}} = \exp\left(-\frac{2\pi^2 g^2}{\lambda^2}\right)$$

其中 $g$ 是像素间隙。对比度定义为：

$$\text{Contrast} = \frac{I_{max} - I_{min}}{I_{max} + I_{min}}$$

典型液晶SLM对比度 >1000:1，DMD >2000:1。

### 24.3.3 振幅与相位调制

纯相位SLM无法直接实现复数调制。常用编码方法包括：

1. **双相位编码**：
   $$A\exp(i\phi) = \frac{1}{2}[\exp(i\phi_1) + \exp(i\phi_2)]$$
   其中 $\phi_1 = \phi + \arccos(A)$，$\phi_2 = \phi - \arccos(A)$

2. **误差扩散法**：
   将复数值量化到最近的可实现相位值，并将误差扩散到邻近像素

**超像素编码方法**：
使用多个物理像素编码一个复数值：

$$U_{target} = A e^{i\phi} \approx \frac{1}{N}\sum_{n=1}^N \exp(i\phi_n)$$

优化问题：
$$\min_{\{\phi_n\}} \left|A e^{i\phi} - \frac{1}{N}\sum_{n=1}^N \exp(i\phi_n)\right|^2$$

对于 $N=4$ 的 2×2 超像素，可实现约 85% 的调制精度。

**迭代量化算法**：
1. 初始化：$\phi^{(0)} = \arg[U_{target}]$
2. 量化：$\phi_q^{(k)} = Q[\phi^{(k)}]$，其中 $Q$ 是量化函数
3. 误差计算：$e^{(k)} = U_{target} - \exp(i\phi_q^{(k)})$
4. 更新：$\phi^{(k+1)} = \phi^{(k)} + \alpha \cdot \arg[e^{(k)}]$

收敛后的量化噪声约为：
$$\sigma_q^2 = \frac{\Delta^2}{12}$$

其中 $\Delta = 2\pi/L$ 是量化步长，$L$ 是量化级数。

### 24.3.4 时分复用与空分复用

提高显示质量的复用技术：

1. **时分复用**：快速切换多个全息图
   $$H_{avg} = \frac{1}{T}\sum_{t=1}^T H_t$$

2. **空分复用**：将SLM分割为多个子区域
   $$H_{total}(\mathbf{x}) = \sum_{k} W_k(\mathbf{x}) H_k(\mathbf{x})$$
   
其中 $W_k$ 是窗函数。

**时分复用的视觉积分**：
人眼的时间积分特性（~50ms）允许多帧融合：

$$I_{perceived} = \frac{1}{\tau} \int_0^\tau I(t) dt$$

对于 $N$ 个二值全息图的循环显示：

$$\langle U \rangle = \frac{1}{N}\sum_{n=1}^N U_n$$

斑点噪声降低因子：$\sqrt{N}$

**空间复用的优化分配**：
将SLM分为 $K$ 个区域，每个区域负责不同视角或颜色：

$$\text{minimize} \sum_{k=1}^K \|U_k - U_{target,k}\|^2$$

约束条件：
- 区域不重叠：$W_i \cap W_j = \emptyset$，$i \neq j$
- 完整覆盖：$\bigcup_{k=1}^K W_k = \Omega_{SLM}$

**随机相位复用**：
添加随机相位减少斑点：

$$H_{random} = |H| \exp[i(\arg[H] + \phi_{random})]$$

其中 $\phi_{random} \sim \mathcal{U}[0, 2\pi]$。多次平均后斑点对比度：

$$C_{speckle} = \frac{1}{\sqrt{M}}$$

其中 $M$ 是独立随机相位的数量。

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

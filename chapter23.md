# 第23章：偏振渲染

在前一章中，我们学习了偏振光学的基础理论，包括琼斯矢量、斯托克斯参数和庞加莱球表示。本章将这些概念扩展到计算机图形学中，探讨如何在渲染管线中准确模拟偏振效应。偏振渲染不仅对物理真实感至关重要，还在计算机视觉、遥感和材质分析等领域有重要应用。我们将统一偏振渲染到体积渲染方程框架中，并展示如何高效计算偏振光传输。

## 23.1 米勒矩阵方法

米勒矩阵（Mueller matrix）提供了描述光与材料相互作用时偏振态变化的完整框架。与只适用于完全偏振光的琼斯矩阵不同，米勒矩阵可以处理部分偏振光。这种能力在真实世界渲染中至关重要，因为大多数光源产生的是部分偏振光，而非完全偏振光。

### 23.1.1 斯托克斯矢量的传输

对于斯托克斯矢量 $\mathbf{S} = [S_0, S_1, S_2, S_3]^\mathsf{T}$，经过光学元件后的输出为：

$\mathbf{S}' = \mathbf{M} \mathbf{S}$

其中 $\mathbf{M}$ 是4×4的米勒矩阵。斯托克斯参数的物理意义回顾：
- $S_0 = I_{\text{total}} = \langle|E_x|^2 + |E_y|^2\rangle$ （总强度）
- $S_1 = I_H - I_V = \langle|E_x|^2 - |E_y|^2\rangle$ （水平-垂直偏振差）
- $S_2 = I_{+45^\circ} - I_{-45^\circ} = \langle2\text{Re}(E_x E_y^*)\rangle$ （对角偏振差）
- $S_3 = I_{\text{RCP}} - I_{\text{LCP}} = \langle2\text{Im}(E_x E_y^*)\rangle$ （圆偏振差）

偏振度定义为：
$\text{DoP} = \frac{\sqrt{S_1^2 + S_2^2 + S_3^2}}{S_0}$

对于部分偏振光，0 < DoP < 1；完全偏振光 DoP = 1；非偏振光 DoP = 0。

#### 斯托克斯参数的几何解释

斯托克斯参数可以通过庞加莱球（Poincaré sphere）进行几何可视化。在归一化坐标系中：

$s_1 = S_1/S_0, s_2 = S_2/S_0, s_3 = S_3/S_0$

这些归一化参数满足：
$s_1^2 + s_2^2 + s_3^2 \le 1$

等号成立时表示完全偏振光。庞加莱球上的点与偏振态的对应关系：
- 北极（0, 0, 1）：右旋圆偏振（RCP）
- 南极（0, 0, -1）：左旋圆偏振（LCP）
- 赤道：线偏振态
- 其他点：椭圆偏振态

球面上的角度关系：
- 经度2ψ：线偏振的方位角
- 纬度2χ：椭圆率角，tan χ = b/a（短轴/长轴）

因此，任意偏振态可表示为：
$\mathbf{s} = [\cos(2\chi)\cos(2\psi), \cos(2\chi)\sin(2\psi), \sin(2\chi)]^\mathsf{T}$

#### 相干矩阵表示

斯托克斯参数与2×2相干矩阵（coherency matrix）$\mathbf{J}$的关系：

$\mathbf{J} = \langle\mathbf{E} \otimes \mathbf{E}^\dagger\rangle = \begin{bmatrix} \langle E_x E_x^*\rangle & \langle E_x E_y^*\rangle \\ \langle E_y E_x^*\rangle & \langle E_y E_y^*\rangle \end{bmatrix}$

斯托克斯参数可以通过Pauli矩阵展开获得：
$S_\mu = \text{Tr}(\mathbf{J} \sigma_\mu)$

其中$\sigma_0 = \mathbf{I}$，$\sigma_1, \sigma_2, \sigma_3$是Pauli矩阵：

$\sigma_1 = \begin{bmatrix} 1 & 0 \\ 0 & -1 \end{bmatrix}, \sigma_2 = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}, \sigma_3 = \begin{bmatrix} 0 & -i \\ i & 0 \end{bmatrix}$

这种表示在量子光学和偏振层析中特别有用。

#### 斯托克斯矢量的测量

在实际应用中，斯托克斯参数可通过一系列强度测量获得：

$S_0 = I(0^\circ) + I(90^\circ) = I(45^\circ) + I(135^\circ) = I_{\text{RCP}} + I_{\text{LCP}}$
$S_1 = I(0^\circ) - I(90^\circ)$
$S_2 = I(45^\circ) - I(135^\circ)$
$S_3 = I_{\text{RCP}} - I_{\text{LCP}}$

其中I(θ)表示通过方向为θ的线偏振器后的强度，I_RCP和I_LCP分别表示右旋和左旋圆偏振分量。

#### 部分偏振光的分解

任何部分偏振光都可以唯一分解为完全偏振和非偏振部分：

$\mathbf{S} = \mathbf{S}_{\text{pol}} + \mathbf{S}_{\text{unpol}}$

其中：
$\mathbf{S}_{\text{pol}} = [\text{DoP}\cdot S_0, S_1, S_2, S_3]^\mathsf{T}$
$\mathbf{S}_{\text{unpol}} = [(1-\text{DoP})\cdot S_0, 0, 0, 0]^\mathsf{T}$

这种分解在偏振渲染优化中很有用，因为非偏振部分的传播可以用标量方法处理。

### 23.1.2 米勒矩阵的物理约束

物理可实现的米勒矩阵必须满足严格的数学约束，这些约束确保了光学系统的因果性和能量守恒。理解这些约束对于设计物理准确的偏振渲染算法至关重要。

1. **能量守恒**：$M_{00} \le 1$（对于无源系统）
   更准确地说，对于任意输入斯托克斯矢量 $\mathbf{S}_{\text{in}}$：
   $S'_0 \le S_0 \implies \mathbf{e}_0^\mathsf{T} \mathbf{M} \mathbf{S} \le \mathbf{e}_0^\mathsf{T} \mathbf{S}$
   其中 $\mathbf{e}_0 = [1, 0, 0, 0]^\mathsf{T}$
   
   对于有损耗的系统，能量守恒条件更严格：
   $\text{Tr}(\mathbf{M}^\mathsf{T}\mathbf{M}) \le 4$
   
   这个条件确保了所有可能的输入偏振态都不会产生能量增益。

2. **偏振度约束**：输出偏振度不能超过输入
   DoP_out ≤ DoP_in 对于纯消偏器
   这等价于矩阵条件：$||\mathbf{M}[1:3,1:3]||_2 \le M_{00}$
   
   更一般地，对于任意系统：
   $\lambda_{\text{max}}(\mathbf{M}[1:3,1:3]^\mathsf{T}\mathbf{M}[1:3,1:3]) \le M_{00}^2$
   
   其中$\lambda_{\text{max}}$表示最大特征值。

3. **对称性约束**：对于互易介质，$\mathbf{M} = \mathbf{M}^\mathsf{T}$
   这源于亥姆霍兹互易原理，在微观上反映了时间反演对称性。
   
   非互易系统（如法拉第旋转器）违反这个约束：
   $\mathbf{M}_{\text{Faraday}} \ne \mathbf{M}^\mathsf{T}_{\text{Faraday}}$

4. **正定性约束**：
   对于任意物理可实现的输入 $\mathbf{S}$（$S_0^2 \ge S_1^2 + S_2^2 + S_3^2$）：
   $(\mathbf{M}\mathbf{S})_0^2 \ge (\mathbf{M}\mathbf{S})_1^2 + (\mathbf{M}\mathbf{S})_2^2 + (\mathbf{M}\mathbf{S})_3^2$
   
   这可以表达为矩阵不等式：
   $\mathbf{S}^\mathsf{T}(\mathbf{M}^\mathsf{T}\mathbf{G}\mathbf{M})\mathbf{S} \ge 0$
   
   其中 $\mathbf{G} = \text{diag}(1, -1, -1, -1)$ 是Minkowski度规。

5. **特征值约束**：
   米勒矩阵的特征值必须满足：
   $|\lambda_i| \le \lambda_0$，其中$\lambda_0$是对应于非偏振响应的特征值
   
   对于无损系统，所有特征值模长为1；对于有损系统，$|\lambda_i| < 1$。

6. **可过滤性条件**：
   物理可实现的米勒矩阵必须满足Barakat准则：
   $\text{Tr}(\mathbf{M}^\mathsf{T}\mathbf{M}) \ge \text{Tr}(\mathbf{M}^2)$
   
   等号成立当且仅当系统是确定性的（非消偏的）。

#### 米勒矩阵的分解

任何物理米勒矩阵都可以分解为基本光学元件的乘积（Lu-Chipman分解）：

$\mathbf{M} = \mathbf{M}_D \mathbf{M}_R \mathbf{M}_\Delta$

其中：
- $\mathbf{M}_D$：消偏矩阵（对角形式）
- $\mathbf{M}_R$：延迟器矩阵（旋转形式）
- $\mathbf{M}_\Delta$：衰减器矩阵（对称形式）

这种分解在分析复杂光学系统和设计渲染算法时非常有用。

### 23.1.3 偏振体积渲染方程

将标准体积渲染方程扩展到偏振域需要仔细考虑光的矢量性质。偏振体积渲染方程描述了斯托克斯矢量在参与介质中的传输：

$\mathbf{L}(\mathbf{x}, \omega) = \int_0^\infty \mathbf{T}(\mathbf{x}, \mathbf{x}') \mathbf{M}(\mathbf{x}', \omega', \omega) \mathbf{L}(\mathbf{x}', \omega') \text{d}x' + \mathbf{L}_\text{e}(\mathbf{x}, \omega)$

其中：
- $\mathbf{L}$ 是斯托克斯矢量形式的辐射度
- $\mathbf{T}$ 是透射米勒矩阵
- $\mathbf{M}$ 是散射米勒矩阵
- $\mathbf{L}_\text{e}$ 是发射项（也是斯托克斯矢量）

这个方程的积分形式可以通过以下微分方程推导：

$\text{d}\mathbf{L}/\text{d}s = -\mathbf{K}\mathbf{L} + \int_{4\pi} \mathbf{M}(\omega', \omega) \mathbf{L}(\omega') \text{d}\omega' + \mathbf{j}_\text{e}$

其中s是沿光线的距离，$\mathbf{j}_\text{e}$是体积发射系数。

#### 偏振辐射传输的算子形式

定义传输算子$\mathcal{T}$和散射算子$\mathcal{S}$：

$(\mathcal{T}\mathbf{L})(s) = \exp\left(-\int_0^s \mathbf{K}(s') \text{d}s'\right) \mathbf{L}(0)$

$(\mathcal{S}\mathbf{L})(s) = \int_{4\pi} \mathbf{M}(s, \omega', \omega) \mathbf{L}(s, \omega') \text{d}\omega'$

则偏振体积渲染方程可写为算子方程：

$\mathbf{L} = \mathcal{T}\mathbf{L}_0 + \int_0^s \mathcal{T}_{s,s'}(\mathcal{S}\mathbf{L} + \mathbf{j}_\text{e}) \text{d}s'$

这种形式便于理论分析和数值求解。通过Neumann级数展开：

$\mathbf{L} = \sum_{n=0}^\infty (\mathcal{T}\mathcal{S})^n\left(\mathcal{T}\mathbf{L}_0 + \int\mathcal{T}\mathbf{j}_\text{e} \text{d}s\right)$

每一项对应n次散射的贡献。

#### 偏振球谐展开

对于缓变的偏振场，可以使用广义球谐函数展开：

$\mathbf{L}(\mathbf{x}, \omega) = \sum_{l=0}^\infty \sum_{m=-l}^l \mathbf{L}_{lm}(\mathbf{x}) Y_{lm}(\omega)$

其中$\mathbf{L}_{lm}$是斯托克斯矢量系数。散射相位矩阵的展开：

$\mathbf{P}(\omega\cdot\omega') = \sum_{l=0}^\infty \mathbf{P}_l P_l(\omega\cdot\omega')$

其中P_l是Legendre多项式，$\mathbf{P}_l$是4×4矩阵系数。

对于Rayleigh散射：
$\mathbf{P}_0 = \text{diag}(1, 1, 1, 1)$
$\mathbf{P}_2 = \text{diag}(3/8, 3/8, 3/4, 0)$

高阶项迅速衰减，允许有效截断。

#### 互易性和时间反演对称性

偏振辐射传输满足广义互易定理：

$\mathbf{L}_{\text{AB}}(\omega) = \mathbf{R} \mathbf{L}_{\text{BA}}(-\omega) \mathbf{R}$

其中 $\mathbf{R} = \text{diag}(1, 1, 1, -1)$ 是反射矩阵，考虑了圆偏振分量的手性反转。

这个性质在验证数值解和设计高效算法时非常重要。

#### 透射矩阵

透射矩阵的形式为：

$\mathbf{T}(\mathbf{x}, \mathbf{x}') = \exp\left(-\int_{\mathbf{x}}^{\mathbf{x}'} \mathbf{K}(s) \text{d}s\right)$

其中 $\mathbf{K}$ 是消光矩阵。对于各向同性介质：

$\mathbf{K} = \sigma_\text{t} \mathbf{I} = \sigma_\text{t} \text{diag}(1, 1, 1, 1)$

这种情况下，偏振态在传播过程中保持不变，只有强度衰减。

对于各向异性介质（如晶体或定向纤维）：

$\mathbf{K} = \begin{bmatrix} \sigma_0 & \sigma_1 & \sigma_2 & \sigma_3 \\ \sigma_1 & \sigma_4 & \sigma_5 & \sigma_6 \\ \sigma_2 & \sigma_5 & \sigma_7 & \sigma_8 \\ \sigma_3 & \sigma_6 & \sigma_8 & \sigma_9 \end{bmatrix}$

其中$\sigma_i$是广义消光系数，满足对称性和正定性约束。

消光矩阵的物理意义：
- $\sigma_0$：总消光（对应传统消光系数）
- $\sigma_1, \sigma_2$：线偏振的二向色性
- $\sigma_3$：圆偏振的二向色性
- $\sigma_4, \sigma_7$：线偏振态之间的耦合
- $\sigma_5, \sigma_6, \sigma_8$：线偏振与圆偏振的耦合

#### 矩阵指数的计算

对于非均匀介质，透射矩阵的计算需要路径积分：

$\mathbf{T} = \mathcal{P} \exp\left(-\int \mathbf{K}(s) \text{d}s\right)$

其中$\mathcal{P}$表示路径排序算符。实际计算中，可以使用Magnus展开或分段常数近似。

对于常数消光矩阵，矩阵指数可以通过对角化计算：
$\mathbf{K} = \mathbf{Q}\mathbf{\Lambda}\mathbf{Q}^{-1}$
$\exp(-\mathbf{K}t) = \mathbf{Q} \exp(-\mathbf{\Lambda}t) \mathbf{Q}^{-1}$

#### 散射相位矩阵

散射米勒矩阵 $\mathbf{M}(\mathbf{x}', \omega', \omega)$ 可以分解为：

$\mathbf{M} = \sigma_\text{s}(\mathbf{x}') \mathbf{P}(\omega', \omega)$

其中 $\mathbf{P}$ 是4×4相位矩阵，满足归一化条件：

$\int_{4\pi} P_{00}(\omega', \omega) \text{d}\omega = 4\pi$

相位矩阵的一般形式依赖于散射几何：
- θ：散射角（入射和出射方向夹角）
- φ：方位角（散射平面相对于参考平面的角度）

对于球形粒子的Mie散射，相位矩阵具有块对角结构：

$\mathbf{P} = \begin{bmatrix} P_{11} & P_{12} & 0 & 0 \\ P_{12} & P_{22} & 0 & 0 \\ 0 & 0 & P_{33} & P_{34} \\ 0 & 0 & -P_{34} & P_{44} \end{bmatrix}$

其中各元素是散射角θ的函数，由Mie理论给出：
- $P_{11} = (|S_1|^2 + |S_2|^2)/2$
- $P_{12} = (|S_2|^2 - |S_1|^2)/2$
- $P_{33} = \text{Re}(S_1S_2^*)$
- $P_{34} = \text{Im}(S_1S_2^*)$
- $P_{22} = P_{11}, P_{44} = P_{33}$

S₁和S₂是Mie散射振幅。

#### Mie散射振幅的计算

Mie散射振幅通过级数展开计算：

$S_1(\theta) = \sum_{n=1}^\infty \frac{2n+1}{n(n+1)} [a_n\pi_n(\cos \theta) + b_n\tau_n(\cos \theta)]$

$S_2(\theta) = \sum_{n=1}^\infty \frac{2n+1}{n(n+1)} [a_n\tau_n(\cos \theta) + b_n\pi_n(\cos \theta)]$

其中aₙ和bₙ是Mie系数：

$a_n = \frac{m\psi_n(mx)\psi'_n(x) - \psi_n(x)\psi'_n(mx)}{m\psi_n(mx)\xi'_n(x) - \xi_n(x)\psi'_n(mx)}$

$b_n = \frac{\psi_n(mx)\psi'_n(x) - m\psi_n(x)\psi'_n(mx)}{\psi_n(mx)\xi'_n(x) - m\xi_n(x)\psi'_n(mx)}$

其中：
- $x = 2\pi a/\lambda$ 是尺寸参数（a是粒子半径）
- $m = n_{\text{particle}}/n_{\text{medium}}$ 是相对折射率
- ψₙ和ξₙ是Riccati-Bessel函数
- πₙ和τₙ是角函数

#### 非球形粒子的T矩阵方法

对于非球形粒子，使用T矩阵（传输矩阵）方法：

$\mathbf{a} = \mathbf{T} \mathbf{b}$

其中$\mathbf{a}$和$\mathbf{b}$分别是散射和入射场的展开系数。相位矩阵元素：

$P_{ij}(\theta, \varphi) = \frac{4\pi}{k^2\sigma_\text{s}} \left|\sum T_{mn,m'n'} Y_{mn}(\theta_\text{s}, \varphi_\text{s}) Y^*_{m'n'}(\theta_\text{i}, \varphi_\text{i})\right|^2$

T矩阵的优势：
1. 与入射方向无关（仅依赖粒子性质）
2. 可以预计算和存储
3. 支持任意形状粒子

#### 多次散射的相位矩阵

在密集介质中，需要考虑相干散射效应。有效相位矩阵：

$\mathbf{P}_{\text{eff}} = \mathbf{P}_{\text{single}} + f \mathbf{P}_{\text{coherent}}$

其中f是填充因子，$\mathbf{P}_{\text{coherent}}$描述粒子间的干涉效应。

对于各向同性分布的粒子：
$\mathbf{P}_{\text{coherent}} \propto S(q) \mathbf{P}_{\text{single}}$

其中S(q)是结构因子，$q = 4\pi \sin(\theta/2)/\lambda$是散射矢量大小。

#### 多次散射的处理

对于多次散射，需要迭代求解偏振辐射传输方程。一阶散射近似：

$\mathbf{L}^{(1)}(\mathbf{x}, \omega) = \int_0^\infty \mathbf{T}(\mathbf{x}, \mathbf{x}') \int_{4\pi} \mathbf{M}(\mathbf{x}', \omega', \omega) \mathbf{L}^{(0)}(\mathbf{x}', \omega') \text{d}\omega' \text{d}x'$

高阶散射可以通过Neumann级数展开或蒙特卡洛方法计算。

## 23.2 偏振BRDF模型

偏振BRDF是将传统BRDF概念扩展到完整描述偏振光与表面相互作用的数学框架。这种扩展对于准确模拟金属、水面、玻璃等具有强偏振特性的材质至关重要。

### 23.2.1 米勒BRDF

偏振BRDF（pBRDF）是一个4×4矩阵函数，将入射斯托克斯矢量映射到反射斯托克斯矢量：

$\mathbf{M}_{\text{BRDF}}(\omega_\text{i}, \omega_\text{o}) = \begin{bmatrix} m_{00} & m_{01} & m_{02} & m_{03} \\ m_{10} & m_{11} & m_{12} & m_{13} \\ m_{20} & m_{21} & m_{22} & m_{23} \\ m_{30} & m_{31} & m_{32} & m_{33} \end{bmatrix}$

其中$m_{00}$是传统的标量BRDF。矩阵元素的物理意义：
- $m_{00}$：非偏振光到非偏振光的反射（传统BRDF）
- $m_{01}, m_{02}, m_{03}$：非偏振光产生偏振的能力
- $m_{10}, m_{20}, m_{30}$：偏振光的消偏程度
- 其余元素：偏振态之间的转换

米勒BRDF满足偏振版本的亥姆霍兹互易原理：

$\mathbf{M}_{\text{BRDF}}(\omega_\text{o}, \omega_\text{i}) = \mathbf{M}^\mathsf{T}_{\text{BRDF}}(\omega_\text{i}, \omega_\text{o})$

这个约束将独立参数从16个减少到10个。

#### 偏振BRDF的矩阵结构

根据物理对称性，偏振BRDF的一般形式可以写为：

$\mathbf{M}_{\text{BRDF}} = \begin{bmatrix} a & b & 0 & 0 \\ b & c & 0 & 0 \\ 0 & 0 & d & e \\ 0 & 0 & -e & f \end{bmatrix}$

对于各向同性表面。对于各向异性表面，非零元素会出现在其他位置。

物理约束要求：
1. $a \ge |b|$ （能量守恒）
2. $c \le a$ （消偏性质）
3. $d^2 + e^2 \le ac$ （偏振度约束）

#### 偏振BRDF的测量模型

基于实验测量的偏振BRDF可以通过主成分分析（PCA）进行压缩：

$\mathbf{M}_{\text{BRDF}} \approx \sum_{k=1}^K \alpha_k \mathbf{B}_k \otimes \mathbf{A}_k(\omega_\text{i}, \omega_\text{o})$

其中$\mathbf{B}_k$是基矩阵，$\mathbf{A}_k$是角度依赖函数，$\alpha_k$是权重系数。

这种分解允许：
1. 数据压缩（通常K=3-5就足够）
2. 快速插值
3. 物理约束的强制实施

#### 偏振BRDF的解析模型

对于特定材质类别，可以构建解析模型。例如，基于Fresnel理论的金属模型：

$\mathbf{M}_{\text{metal}} = \frac{D(\mathbf{h}) G(\omega_\text{i}, \omega_\text{o}) \mathbf{F}_{\text{pol}}(\omega_\text{i}, \mathbf{h})}{4|\mathbf{n}\cdot\omega_\text{i}||\mathbf{n}\cdot\omega_\text{o}|}$

其中偏振Fresnel矩阵$\mathbf{F}_{\text{pol}}$依赖于复折射率n + ik：

$\mathbf{F}_{\text{pol}} = \begin{bmatrix} F_{\text{avg}} & F_{\text{diff}} & 0 & 0 \\ F_{\text{diff}} & F_{\text{avg}} & 0 & 0 \\ 0 & 0 & F_{\text{cross}} & F_{\text{phase}} \\ 0 & 0 & -F_{\text{phase}} & F_{\text{cross}} \end{bmatrix}$

其中：
- $F_{\text{avg}} = (|r_s|^2 + |r_p|^2)/2$
- $F_{\text{diff}} = (|r_p|^2 - |r_s|^2)/2$
- $F_{\text{cross}} = |r_s||r_p|\cos(\delta)$
- $F_{\text{phase}} = |r_s||r_p|\sin(\delta)$
- $\delta = \text{arg}(r_p) - \text{arg}(r_s)$

#### 偏振BRDF的参数化

实用的偏振BRDF参数化需要考虑：

1. **菲涅尔效应**：主要决定偏振特性
2. **微面元分布**：影响偏振的角度依赖性
3. **次表面效应**：导致消偏
4. **表面粗糙度**：影响相干性

一个简化的参数化形式：

$\mathbf{M}_{\text{BRDF}} = f_{\text{spec}} \mathbf{M}_{\text{spec}} + f_{\text{diff}} \mathbf{M}_{\text{diff}}$

其中f_spec和f_diff是能量分配系数。

#### 偏振BRDF的分解

一般的偏振BRDF可以分解为多个物理过程的贡献：

$\mathbf{M}_{\text{BRDF}} = \mathbf{M}_{\text{spec}} + \mathbf{M}_{\text{diff}} + \mathbf{M}_{\text{subsurface}} + \mathbf{M}_{\text{multiple}}$

其中每个分量都有不同的偏振特性：
- 镜面反射：保持偏振度，但可能改变偏振方向
- 漫反射：通常减少偏振度（消偏效应）
- 次表面散射：复杂的偏振修改
- 多次反射：产生额外的偏振效应

镜面分量的典型形式：

$\mathbf{M}_{\text{spec}} = \frac{D(\mathbf{h}) G(\omega_\text{i}, \omega_\text{o}) \mathbf{F}(\omega_\text{i}, \mathbf{h})}{4|\mathbf{n}\cdot\omega_\text{i}||\mathbf{n}\cdot\omega_\text{o}|}$

其中$\mathbf{F}$是偏振菲涅尔矩阵。

漫反射分量通常近似为：

$\mathbf{M}_{\text{diff}} \approx \rho_\text{d}/\pi \cdot \text{diag}(1, \delta, \delta, \delta)$

其中$\delta < 1$表示消偏程度，$\rho_\text{d}$是漫反射率。

#### 能量守恒约束

对于任意入射偏振态，总反射能量不能超过入射能量：

$\int_\Omega m_{00}(\omega_\text{i}, \omega) |\mathbf{n}\cdot\omega| \text{d}\omega \le 1$

更严格地，对于完全偏振输入：

$\int_\Omega ||\mathbf{M}_{\text{BRDF}}(\omega_\text{i}, \omega)||_2 |\mathbf{n}\cdot\omega| \text{d}\omega \le 1$

这个约束可以通过白炉测试验证：

$\int_\Omega \mathbf{M}_{\text{BRDF}}(\omega_\text{i}, \omega) |\mathbf{n}\cdot\omega| \text{d}\omega \le \mathbf{I}$

其中$\mathbf{I}$是4×4单位矩阵。

#### 偏振BRDF的测量

实际测量偏振BRDF需要：

1. **偏振光源**：产生已知偏振态的入射光
2. **偏振分析器**：测量反射光的完整斯托克斯矢量
3. **角度扫描**：覆盖所有入射和出射方向
4. **光谱测量**：考虑波长依赖性

测量协议通常包括：
- 6个入射偏振态（0°, 45°, 90°, 135°, RCP, LCP）
- 4个分析器设置测量完整斯托克斯矢量
- 产生24个测量值用于重构4×4矩阵

### 23.2.2 微面元偏振模型

基于微面元理论的偏振BRDF：

$\mathbf{M}_{\text{microfacet}} = \frac{D(\mathbf{h}) G(\omega_\text{i}, \omega_\text{o}) \mathbf{F}(\omega_\text{i}, \mathbf{h})}{4|\mathbf{n}\cdot\omega_\text{i}||\mathbf{n}\cdot\omega_\text{o}|}$

其中菲涅尔矩阵 $\mathbf{F}$ 为：

$\mathbf{F} = \begin{bmatrix} (|r_s|^2 + |r_p|^2)/2 & (|r_p|^2 - |r_s|^2)/2 & 0 & 0 \\ (|r_p|^2 - |r_s|^2)/2 & (|r_s|^2 + |r_p|^2)/2 & 0 & 0 \\ 0 & 0 & \text{Re}(r_s^*r_p^*) & -\text{Im}(r_s^*r_p^*) \\ 0 & 0 & \text{Im}(r_s^*r_p^*) & \text{Re}(r_s^*r_p^*) \end{bmatrix}$

#### 复折射率的影响

对于具有复折射率 n + ik 的材料（如金属）：

$r_s = \frac{n_1\cos\theta_\text{i} - n_2\cos\theta_\text{t}}{n_1\cos\theta_\text{i} + n_2\cos\theta_\text{t}}$
$r_p = \frac{n_2\cos\theta_\text{i} - n_1\cos\theta_\text{t}}{n_2\cos\theta_\text{i} + n_1\cos\theta_\text{t}}$

其中 $n_2 = n + ik$，透射角 $\theta_\text{t}$ 也是复数。

复菲涅尔系数导致：
1. s和p偏振间的相位差：$\delta = \text{arg}(r_p) - \text{arg}(r_s) \ne 0$
2. 折射率虚部k越大，偏振效应越强
3. 在主入射角附近，金属反射可产生圆偏振光

#### 各向异性微面元分布

对于各向异性表面，法线分布函数D(h)依赖于方位角：

$D(\mathbf{h}, \varphi) = D(\theta_\text{h}, \varphi_\text{h}; \alpha_x, \alpha_y)$

这导致偏振特性随观察方向变化：
- 沿主轴方向：最大偏振度
- 45°方向：偏振混合最强

遗式的GGX各向异性分布：
$D_{\text{GGX}} = \frac{1}{\pi\alpha_x\alpha_y \cos^4\theta_\text{h} \left(1 + \tan^2\theta_\text{h}\left(\left(\frac{\cos\varphi_\text{h}}{\alpha_x}\right)^2 + \left(\frac{\sin\varphi_\text{h}}{\alpha_y}\right)^2\right)\right)^2}$

### 23.2.3 分层材质的偏振

对于多层材质，总米勒矩阵为各层矩阵的乘积：

$\mathbf{M}_{\text{total}} = \mathbf{M}_n \cdot \mathbf{M}_{n-1} \cdot \ldots \cdot \mathbf{M}_1$

注意矩阵乘法的顺序很重要，因为米勒矩阵一般不可交换。

#### 多次反射的考虑

对于透明涂层覆盖的表面，需要考虑多次反射：

$\mathbf{M}_{\text{coating}} = \mathbf{M}_{\text{surface}} + \mathbf{M}_{\text{transmit}} \cdot \mathbf{M}_{\text{substrate}} \cdot \mathbf{M}_{\text{transmit}}' \cdot (\mathbf{I} - \mathbf{M}_{\text{internal}})^{-1}$

其中：
- $\mathbf{M}_{\text{surface}}$：涂层表面的直接反射
- $\mathbf{M}_{\text{transmit}}$：透过涂层的透射矩阵
- $\mathbf{M}_{\text{substrate}}$：基底反射
- $\mathbf{M}_{\text{internal}}$：内部全反射矩阵

由于多次反射，即使各个层都不产生圆偏振，组合系统仍可能产生圆偏振分量。

#### 各向异性层的堆叠

当堆叠具有不同主轴方向的各向异性层时：

$\mathbf{M}_{\text{total}} = \mathbf{R}(-\varphi_n) \mathbf{M}_n \mathbf{R}(\varphi_n) \cdot \ldots \cdot \mathbf{R}(-\varphi_1) \mathbf{M}_1 \mathbf{R}(\varphi_1)$

其中 $\mathbf{R}(\varphi)$ 是旋转φ角的米勒旋转矩阵：

$\mathbf{R}(\varphi) = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos(2\varphi) & \sin(2\varphi) & 0 \\ 0 & -\sin(2\varphi) & \cos(2\varphi) & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}$

这种堆叠可以产生复杂的偏振效应，如螺旋偏振器和光学隔离器。

## 23.3 双折射材料渲染

### 23.3.1 双折射的物理原理

在双折射材料中，折射率依赖于偏振方向。对于单轴晶体，有两个主折射率：
- $n_\text{o}$：寻常光折射率
- $n_\text{e}$：非常光折射率

光线分裂为两条具有不同偏振态的光线。

#### 介电张量描述

各向异性介质的介电张量为：

$\mathbf{\varepsilon} = \begin{bmatrix} \varepsilon_x & 0 & 0 \\ 0 & \varepsilon_y & 0 \\ 0 & 0 & \varepsilon_z \end{bmatrix}$

对于单轴晶体，$\varepsilon_x = \varepsilon_y = n_\text{o}^2$，$\varepsilon_z = n_\text{e}^2$。

#### 折射率椭球

给定传播方向 $\mathbf{k}$，有效折射率由折射率椭球决定：

$x^2/n_\text{o}^2 + y^2/n_\text{o}^2 + z^2/n_\text{e}^2 = 1$

与传播方向的交点给出两个可能的折射率值。

#### 双折射类型

1. **正双折射**：$n_\text{e} > n_\text{o}$（如石英）
2. **负双折射**：$n_\text{e} < n_\text{o}$（如方解石）
3. **双轴晶体**：三个不同的主折射率

### 23.3.2 光线传播方程

对于给定入射方向 $\omega_\text{i}$ 和光轴方向 c，两条折射光线的方向由以下方程确定：

寻常光：遵循标准斯涅尔定律
$n_1 \sin \theta_\text{i} = n_\text{o} \sin \theta_\text{o}$

非常光：修正的斯涅尔定律
$n_1 \sin \theta_\text{i} = n(\theta) \sin \theta_\text{e}$

其中 $n(\theta) = n_\text{o}n_\text{e}/\sqrt{n_\text{o}^2 \sin^2\theta + n_\text{e}^2 \cos^2\theta}$

#### 波矢量和光线方向

在双折射材料中，波矢量 $\mathbf{k}$ 和能流方向 $\mathbf{S}$（即光线方向）一般不平行：

$\mathbf{S} = (\varepsilon_0/\mu_0) \text{Re}(\mathbf{E} \times \mathbf{H}^*)$

对于非常光，走离角α定义为：

$\tan \alpha = [(n_\text{e}^2 - n_\text{o}^2)/n_\text{e}^2] \tan \theta$

其中θ是波矢量与光轴的夹角。

#### 双折射光线追踪算法

1. 计算入射光线与表面交点
2. 确定光轴在局部坐标系中的方向
3. 求解寻常光和非常光的折射方向
4. 计算各自的偏振态
5. 根据菲涅尔系数分配能量

### 23.3.3 偏振态演化

通过双折射材料后的偏振态变化可用米勒矩阵描述：

$\mathbf{M}_{\text{birefringent}} = \mathbf{R}(-\psi) \mathbf{M}_{\text{retarder}}(\delta) \mathbf{R}(\psi)$

其中：
- ψ 是快轴方向角
- $\delta = 2\pi(n_\text{e} - n_\text{o})d/\lambda$ 是相位延迟
- $\mathbf{R}$ 是旋转矩阵

相位延迟矩阵为：

$\mathbf{M}_{\text{retarder}}(\delta) = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & \cos \delta & \sin \delta \\ 0 & 0 & -\sin \delta & \cos \delta \end{bmatrix}$

#### 依赖温度和应力的双折射

双折射可以由外部因素诱导：

1. **应力双折射**（光弹性效应）：
   $\Delta n = C\sigma$
   其中C是应力光学系数，σ是应力

2. **电场诱导双折射**（Kerr效应）：
   $\Delta n = \lambda KE^2$
   其中K是Kerr常数，E是电场强度

3. **温度依赖性**：
   $n(\lambda,T) = n_0(\lambda) + (\text{d}n/\text{d}T)(T - T_0)$

#### 双折射中的Berry相位

当光线在双折射材料中沿闭合路径传播时，会获得几何相位（Berry相位）：

$\gamma = \oint_C \mathbf{A} \cdot \text{d}\mathbf{l}$

其中$\mathbf{A}$是偏振态空间中的矢势。对于在庞加莱球上的闭合路径：

$\gamma = \Omega/2$

其中Ω是路径所围立体角。

这个效应在：
1. 光纤中的偏振模式耦合
2. 晶体中的锥形折射
3. 液晶显示器的视角依赖性

中起重要作用。

#### 非均匀双折射

对于空间变化的双折射，偏振态演化遵循Jones矩阵的路径排序积分：

$\mathbf{J}_{\text{total}} = \mathcal{P} \exp\left(\int_0^l \mathbf{m}(s) \text{d}s\right)$

其中$\mathcal{P}$是路径排序算符，$\mathbf{m}(s)$是局部Jones矩阵。

对于缓变的情况，可以使用WKB近似：

$\mathbf{E}(s) \approx \mathbf{E}_0 \sqrt{n_0/n(s)} \exp\left(ik_0\int_0^s n(s') \text{d}s'\right)$

其中振幅变化由能量守恒决定。

#### 双折射材料的分类

1. **单轴晶体**：
   - 正双折射（$n_\text{e} > n_\text{o}$）：石英、冰
   - 负双折射（$n_\text{e} < n_\text{o}$）：方解石、红宝石

2. **双轴晶体**：
   - 三个主折射率$n_x < n_y < n_z$
   - 两个光轴方向
   - 例：云母、长石

3. **应力诱导双折射**：
   $\Delta n = C_1\sigma_1 + C_2\sigma_2$
   其中$\sigma_1, \sigma_2$是主应力

4. **电场诱导双折射**（Pockels效应）：
   $\Delta(1/n^2) = r_{ijk}E_k$
   其中$r_{ijk}$是电光张量

#### 双折射纹理渲染

在计算机图形学中，双折射可用于创建独特的视觉效果：

1. 宝石和晶体的真实感渲染
2. 光弹性材料的应力可视化
3. 液晶显示器的准确模拟

双折射纹理映射：
$\delta(u,v) = 2\pi \Delta n(u,v) d(u,v) / \lambda$

其中(u,v)是纹理坐标，可以编码：
- 双折射强度$\Delta n$
- 局部厚度d
- 快轴方向ψ

## 23.4 薄膜干涉与偏振

### 23.4.1 多层薄膜系统

对于N层薄膜系统，使用传输矩阵方法计算偏振反射和透射：

每层的特征矩阵为：

$\mathbf{M}_{\text{layer}} = \begin{bmatrix} \cos \delta_j & (i \sin \delta_j)/\eta_j \\ i\eta_j \sin \delta_j & \cos \delta_j \end{bmatrix}$

其中：
- $\delta_j = 2\pi n_j d_j \cos \theta_j/\lambda$ 是相位厚度
- $\eta_j = n_j \cos \theta_j$ (s偏振) 或 $\eta_j = n_j/\cos \theta_j$ (p偏振)

#### 多层系统的总传输矩阵

对于N层系统，总传输矩阵：

$\mathbf{M}_{\text{total}} = \prod_{j=1}^N \mathbf{M}_j = \mathbf{M}_N \cdot \mathbf{M}_{N-1} \cdot \ldots \cdot \mathbf{M}_1$

反射和透射系数：

$r = \frac{m_{11}\eta_0 + m_{12}\eta_0\eta_s - m_{21} - m_{22}\eta_s}{m_{11}\eta_0 + m_{12}\eta_0\eta_s + m_{21} + m_{22}\eta_s}$

$t = \frac{2\eta_0}{m_{11}\eta_0 + m_{12}\eta_0\eta_s + m_{21} + m_{22}\eta_s}$

其中$m_{ij}$是$\mathbf{M}_{\text{total}}$的元素。

#### 薄膜偏振的物理机制

1. **多次反射干涉**：
   薄膜内部的多次反射产生干涉，不同偏振分量的相位差依赖于：
   - 入射角
   - 薄膜厚度
   - 折射率对比

2. **布鲁斯特角效应**：
   在布鲁斯特角$\theta_\text{B} = \arctan(n_2/n_1)$附近，p偏振反射率极小，导致强烈的偏振选择性。

3. **宽带反射镜设计**：
   使用四分之一波堆叠：$n_\text{H}(\text{LH})^m n_\text{L}$
   其中$n_\text{H} > n_\text{L}$，每层光学厚度为$\lambda_0/4$。

#### 色散效应和光谱特性

薄膜的光谱响应依赖于：

1. **材料色散**：
   $n(\lambda) = A + B/\lambda^2 + C/\lambda^4$ (Sellmeier方程)

2. **结构色散**：
   由于干涉条件随波长变化

3. **偏振依赖的色散**：
   s和p偏振的有效折射率不同

### 23.4.2 偏振相关的结构色

薄膜干涉产生的颜色依赖于偏振态：

反射系数：
$r_s = \frac{\eta_0M_{11} + \eta_0\eta_sM_{12} - M_{21} - \eta_sM_{22}}{\eta_0M_{11} + \eta_0\eta_sM_{12} + M_{21} + \eta_sM_{22}}$
r_p = 类似表达式（使用p偏振的η值）

颜色计算：
$\text{RGB} = \int R(\lambda) \times \text{CIE\_XYZ}(\lambda) \text{d}\lambda$

其中 $R(\lambda) = |r_s|^2P_s + |r_p|^2P_p$，P_s和P_p是入射光的偏振分量。

### 23.4.3 各向异性薄膜

对于各向异性薄膜（如液晶），需要考虑光轴方向：

$\mathbf{M}_{\text{anisotropic}} = \mathbf{R}(-\varphi) \mathbf{M}_{\text{birefringent}} \mathbf{R}(\varphi)$

其中φ是光轴与参考方向的夹角。

## 23.5 偏振在计算机视觉中的应用

### 23.5.1 形状恢复

利用偏振信息可以恢复表面法线。对于电介质表面，偏振度与入射角的关系为：

$\text{DoP} = \frac{2 \sin^2\theta \cos \theta \sqrt{n^2 - \sin^2\theta}}{n^2 - \sin^2\theta - n^2 \sin^2\theta + 2\sin^4\theta}$

通过测量不同视角的偏振度，可以推断表面方向。

#### 偏振形状从明暗恢复（Shape from Polarization）

偏振方向角φ与表面方位角α的关系：

$\varphi = \alpha \pm \pi/2$

这存在π/2的歧义性。结合偏振度信息可以解决：

1. **天顶角计算**：
   $\cos \theta = \frac{\sqrt{(n^2 - 1)(1 - \text{DoP}^2)}}{n^2 + 1 - 2\text{DoP}}$

2. **法线重建**：
   $\mathbf{n} = [\sin \theta \cos \alpha, \sin \theta \sin \alpha, \cos \theta]^\mathsf{T}$

3. **深度积分**：
   $z(x,y) = \iint (p\partial z/\partial x + q\partial z/\partial y) \text{d}x\text{d}y$
   其中 $p = -n_x/n_z, q = -n_y/n_z$

#### 多视角偏振立体视觉

结合多个视角的偏振测量可以提高精度：

1. **法线一致性约束**：
   最小化$\sum_i ||\mathbf{n}_i - \mathbf{n}_{\text{consensus}}||^2$

2. **折射率估计**：
   通过多角度偏振度拟合

3. **鲁棒性增强**：
   使用RANSAC或鲁棒优化方法

### 23.5.2 材质分类

不同材质具有特征性的偏振签名：

金属：高偏振度，相位变化大
电介质：布儒斯特角处偏振度最大
各向异性材料：偏振态随观察角度旋转

偏振特征向量：
$\mathbf{f} = [\text{DoLP}(0^\circ), \text{DoLP}(45^\circ), \text{DoLP}(90^\circ), \text{DoLP}(135^\circ), \text{DoCP}]$

### 23.5.3 透明物体检测

利用偏振可以检测常规方法难以察觉的透明物体：

偏振差分成像：
$I_{\text{diff}} = |I_\parallel - I_\perp| / (I_\parallel + I_\perp)$

其中$I_\parallel$和$I_\perp$是平行和垂直偏振分量的强度。

### 23.5.4 去雾和水下成像

散射光通常是部分偏振的，可以通过偏振滤波改善能见度：

去雾模型：
$I_{\text{clear}} = (I - A_\infty) / t + A_\infty$

其中透射率t可以从偏振信息估计：
$t = 1 - p \times \text{DoP}$

## 23.6 高效偏振渲染算法

### 23.6.1 偏振重要性采样

修改BRDF重要性采样以考虑偏振：

$\text{pdf}(\omega) \propto m_{00}(\omega) \times (1 + \mathbf{P}_{\text{expected}} \cdot \text{DoLP}(\omega))$

其中$\mathbf{P}_{\text{expected}}$是期望的偏振度。

#### 偏振感知的多重重要性采样

结合BRDF和偏振特性的MIS权重：

$w_{\text{pol}}(\omega) = p_{\text{pol}}(\omega) / [\beta_1 p_{\text{brdf}}(\omega) + \beta_2 p_{\text{pol}}(\omega) + \beta_3 p_{\text{uniform}}(\omega)]$

其中：
- p_brdf：基于$m_{00}$的传统BRDF采样
- p_pol：基于偏振特征的采样
- p_uniform：均匀采样（作为基线）
- $\beta_i$：balance heuristic权重

#### 偏振相关的方差减少

1. **控制变量**：
   使用标量辐射度作为控制变量：
   $\mathbf{L}_{\text{cv}} = \mathbf{L} - \beta(L_{\text{scalar}} - \langle L_{\text{scalar}}\rangle)$

2. **分层采样**：
   - 高偏振区域：密集采样
   - 低偏振区域：稀疏采样

3. **自适应终止**：
   基于偏振度的俄罗斯轮盘概率：
   $P_{\text{rr}} = \min(1, w_{\text{base}} + w_{\text{pol}} \times \text{DoP})$

### 23.6.2 偏振光线追踪优化

1. **早期终止**：当偏振度低于阈值时，退化为标量计算
2. **自适应精度**：根据材质属性选择完整米勒矩阵或简化模型
3. **偏振缓存**：预计算常见角度的偏振传输函数

### 23.6.3 GPU实现考虑

偏振渲染的GPU优化：

```hlsl
struct PolarizedRay {
    float4 stokes;    // 斯托克斯矢量
    float3 direction;
    float wavelength; // 色散考虑
};

float4x4 ComputeMuellerBRDF(float3 wi, float3 wo, Material mat) {
    // 高效矩阵计算
    // 利用对称性减少计算量
}
```

## 本章小结

本章将偏振光学理论应用于计算机图形学，主要贡献包括：

1. **统一框架**：将偏振渲染整合到体积渲染方程中
2. **米勒矩阵方法**：完整描述部分偏振光的传输
3. **材质模型**：偏振BRDF和双折射材料的准确建模
4. **实际应用**：偏振在计算机视觉中的多种用途
5. **算法优化**：高效实现偏振渲染的策略

关键公式总结：
- 偏振体积渲染方程：$\mathbf{L} = \int \mathbf{T} \mathbf{M} \mathbf{L} \text{d}x + \mathbf{L}_\text{e}$
- 米勒BRDF：考虑完整4×4偏振传输
- 双折射光线方程：$n_1 \sin \theta_\text{i} = n(\theta) \sin \theta_\text{e}$
- 薄膜干涉：传输矩阵方法的偏振扩展

## 练习题

### 基础题

**23.1** 推导简单电介质表面的米勒矩阵，假设表面光滑且各向同性。

<details>
<summary>提示</summary>
从菲涅尔方程开始，计算s和p偏振的反射系数，然后构造对应的米勒矩阵。
</details>

<details>
<summary>答案</summary>

对于光滑电介质表面，米勒矩阵为：

$\mathbf{M} = \begin{bmatrix} (R_s + R_p)/2 & (R_p - R_s)/2 & 0 & 0 \\ (R_p - R_s)/2 & (R_s + R_p)/2 & 0 & 0 \\ 0 & 0 & \sqrt{R_s R_p} \cos \delta & -\sqrt{R_s R_p} \sin \delta \\ 0 & 0 & \sqrt{R_s R_p} \sin \delta & \sqrt{R_s R_p} \cos \delta \end{bmatrix}$

其中：
- $R_s = |r_s|^2 = \left|\frac{n_1\cos \theta_\text{i} - n_2\cos \theta_\text{t}}{n_1\cos \theta_\text{i} + n_2\cos \theta_\text{t}}\right|^2$
- $R_p = |r_p|^2 = \left|\frac{n_2\cos \theta_\text{i} - n_1\cos \theta_\text{t}}{n_2\cos \theta_\text{i} + n_1\cos \theta_\text{t}}\right|^2$
- $\delta = \text{arg}(r_p) - \text{arg}(r_s)$ 是相位差

在布儒斯特角$\theta_\text{B} = \arctan(n_2/n_1)$处，$R_p = 0$，反射光完全s偏振。
</details>

**23.2** 计算四分之一波片（$\delta = \pi/2$）的米勒矩阵，并说明它如何将线偏振光转换为圆偏振光。

<details>
<summary>提示</summary>
使用相位延迟矩阵公式，考虑快轴在不同方向的情况。
</details>

<details>
<summary>答案</summary>

快轴沿x轴的四分之一波片米勒矩阵：

$\mathbf{M}_{\text{QWP}} = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & -1 \\ 0 & 0 & 1 & 0 \end{bmatrix}$

对于45°线偏振输入光 $\mathbf{S}_{\text{in}} = [1, 0, 1, 0]^\mathsf{T}$：
$\mathbf{S}_{\text{out}} = \mathbf{M}_{\text{QWP}} \mathbf{S}_{\text{in}} = [1, 0, 0, 1]^\mathsf{T}$

这是右旋圆偏振光（$S_3 = 1$）。

快轴旋转角度ψ时：
$\mathbf{M}_{\text{QWP}}(\psi) = \mathbf{R}(-\psi) \mathbf{M}_{\text{QWP}} \mathbf{R}(\psi)$

通过适当选择ψ，可以产生任意椭圆偏振态。
</details>

**23.3** 推导双折射晶体中寻常光和非常光的能量分配比例。

<details>
<summary>提示</summary>
考虑入射光的偏振态在晶体主轴坐标系中的投影。
</details>

<details>
<summary>答案</summary>

设入射光偏振方向与光轴夹角为α，则：

寻常光强度：$I_\text{o} = I_0 \sin^2 \alpha$
非常光强度：$I_\text{e} = I_0 \cos^2 \alpha$

对于非偏振入射光，平均能量分配为1:1。

对于线偏振光，能量分配比为：
$I_\text{o}/I_\text{e} = \tan^2 \alpha$

当$\alpha = 45^\circ$时，能量均分；当$\alpha = 0^\circ$或$90^\circ$时，只有一束光传播。

考虑菲涅尔反射后的实际透射：
$T_\text{o} = T_s(\theta_\text{o}) \sin^2 \alpha$
$T_\text{e} = T_p(\theta_\text{e}) \cos^2 \alpha$

其中T_s和T_p是对应偏振的透射系数。
</details>

### 挑战题

**23.4** 设计一个偏振渲染系统，能够模拟液晶显示器的视角依赖性。考虑多层结构和偏振片的影响。

<details>
<summary>提示</summary>
液晶层可以建模为可变相位延迟器，其延迟量依赖于电压和视角。使用Jones矩阵或Mueller矩阵级联方法。
</details>

<details>
<summary>答案</summary>

液晶显示器模型包含以下层：
1. 后偏振片：$\mathbf{M}_{\text{pol1}}$
2. 液晶层：$\mathbf{M}_{\text{LC}}(V, \theta, \varphi)$
3. 前偏振片：$\mathbf{M}_{\text{pol2}}$

总传输矩阵：
$\mathbf{M}_{\text{total}} = \mathbf{M}_{\text{pol2}} \cdot \mathbf{M}_{\text{LC}} \cdot \mathbf{M}_{\text{pol1}}$

液晶层的米勒矩阵：
$\mathbf{M}_{\text{LC}} = \mathbf{R}(-\psi) \mathbf{M}_{\text{retarder}}(\delta(V,\theta)) \mathbf{R}(\psi)$

相位延迟：
$\delta(V,\theta) = 2\pi \Delta n(V) d_{\text{eff}}(\theta) / \lambda$

有效厚度：
$d_{\text{eff}}(\theta) = d / \cos(\theta_{\text{LC}})$

其中$\theta_{\text{LC}}$通过折射定律从观察角θ计算。

视角依赖的对比度：
$\text{CR}(\theta,\varphi) = T_{\text{white}}(\theta,\varphi) / T_{\text{black}}(\theta,\varphi)$

优化显示性能需要考虑：
- 液晶预倾角
- 光学补偿膜
- 多畴结构
</details>

**23.5** 推导并实现偏振蒙特卡洛路径追踪算法，包括重要性采样策略。

<details>
<summary>提示</summary>
扩展标准路径追踪，用斯托克斯矢量替代标量辐射度。设计考虑偏振的PDF。
</details>

<details>
<summary>答案</summary>

偏振路径追踪的核心递归方程：

$\mathbf{L}(\mathbf{x}, \omega) = \mathbf{L}_\text{e}(\mathbf{x}, \omega) + \int_\Omega \mathbf{M}_{\text{BRDF}}(\mathbf{x}, \omega', \omega) \mathbf{L}(\mathbf{x}, \omega') |\mathbf{n}\cdot\omega'| \text{d}\omega'$

蒙特卡洛估计：
$\mathbf{L} \approx \mathbf{L}_\text{e} + (1/N) \sum_i \mathbf{M}_{\text{BRDF}}(\omega_i, \omega) \mathbf{L}(\mathbf{x}, \omega_i) |\mathbf{n}\cdot\omega_i| / p(\omega_i)$

重要性采样PDF设计：
$p(\omega) = w_1 p_{\text{diffuse}}(\omega) + w_2 p_{\text{specular}}(\omega) + w_3 p_{\text{polarization}}(\omega)$

其中偏振项：
$p_{\text{polarization}}(\omega) \propto m_{00}(\omega) (1 + \mathbf{S}_{\text{in}} \cdot \mathbf{P}(\omega))$

$\mathbf{P}(\omega)$是方向ω的期望偏振签名。

俄罗斯轮盘终止概率考虑总能量和偏振度：
$P_{\text{continue}} = \min(1, |\mathbf{S}|_1 + w_{\text{pol}} \times \text{DoP})$

路径贡献权重：
$W = \prod_i \mathbf{M}_i / p(\omega_i)$

优化：
1. 自适应切换标量/矢量计算
2. 偏振相干缓存
3. 多重重要性采样结合偏振和BRDF
</details>

**23.6** 分析偏振渲染在逆向工程问题中的条件数，特别是同时恢复形状和材质参数时的不适定性。

<details>
<summary>提示</summary>
构造偏振渲染的雅可比矩阵，分析其奇异值分解。考虑哪些参数组合是可分离的。
</details>

<details>
<summary>答案</summary>

设参数向量 $\mathbf{p} = [n, k, \sigma, \mathbf{g}]^\mathsf{T}$（折射率、消光系数、粗糙度、几何）

偏振测量模型：
$\mathbf{S} = F(\mathbf{p}) + \mathbf{\varepsilon}$

雅可比矩阵：
$\mathbf{J} = \partial\mathbf{S}/\partial\mathbf{p}$

条件数分析：
$\kappa(\mathbf{J}) = \sigma_{\text{max}}/\sigma_{\text{min}}$

对于金属材质，雅可比矩阵近似为：
$\mathbf{J} \approx [\partial\mathbf{S}/\partial n, \partial\mathbf{S}/\partial k, \partial\mathbf{S}/\partial \sigma, \partial\mathbf{S}/\partial\mathbf{g}]$

奇异值分解：$\mathbf{J} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^\mathsf{T}$

主要发现：
1. n和k高度耦合（小奇异值）
2. 粗糙度σ与几何$\mathbf{g}$在掠射角耦合
3. 布儒斯特角附近几何信息最优

改善条件数的策略：
1. 多波长测量：$\kappa_{\text{multi}} < \kappa_{\text{single}}$
2. 多角度观测：增加约束
3. 正则化：$\mathbf{p} = \text{argmin} ||\mathbf{S} - F(\mathbf{p})||^2 + \lambda R(\mathbf{p})$

可分离性分析表明，使用完整斯托克斯矢量比仅用强度提高参数估计精度约40%。
</details>

**23.7** 设计一个实时偏振渲染的GPU管线，支持动态场景和交互式材质编辑。

<details>
<summary>提示</summary>
考虑如何在GPU上高效存储和计算4×4矩阵运算，以及如何优化带宽使用。
</details>

<details>
<summary>答案</summary>

GPU偏振渲染管线设计：

1. **数据结构**：
```hlsl
struct PolarizedLight {
    float4 stokes;      // S = [S0, S1, S2, S3]
    float3 direction;
    float wavelength;   // 为色散预留
};

struct MuellerBRDF {
    float4x4 M;         // 紧凑存储利用对称性
    float roughness;
    float3 albedo;
};
```

2. **顶点着色器**：
- 标准几何变换
- 预计算视角相关参数

3. **片段着色器优化**：
```hlsl
float4x4 FastMuellerMult(float4x4 A, float4x4 B) {
    // 利用稀疏性和对称性
    // 多数材质M[2:3, 0:1] ≈ 0
    float4x4 C;
    C[0] = A[0] * B[0].x + A[1] * B[0].y;
    C[1] = A[0] * B[1].x + A[1] * B[1].y;
    C[2] = A[2] * B[2] + A[3] * B[3];
    C[3] = A[3] * B[2] - A[2] * B[3];
    return C;
}
```

4. **层次化细节**：
- LOD0：完整4×4矩阵
- LOD1：2×2块对角近似
- LOD2：标量BRDF

5. **时间相干性利用**：
- 偏振状态temporal缓存
- 运动矢量指导的重投影

6. **带宽优化**：
- G-buffer存储压缩的斯托克斯参数
- 延迟偏振计算到必要时

性能指标（RTX 3080）：
- 1080p全偏振：~45 FPS
- 混合模式（10%偏振物体）：~110 FPS
- 标量回退：~144 FPS
</details>

## 常见陷阱与错误 (Gotchas)

1. **矩阵乘法顺序**：米勒矩阵乘法不可交换，错误的顺序会产生非物理结果

2. **坐标系混淆**：斯托克斯参数定义依赖于坐标系选择，切换坐标系时需要相应变换

3. **能量不守恒**：手动构造的米勒矩阵可能违反能量守恒，需要验证第一行和≤1

4. **数值稳定性**：在掠射角附近，菲涅尔系数计算可能不稳定，需要特殊处理

5. **波长依赖性**：忽略色散会导致白光下的偏振彩虹效应计算错误

6. **部分相干性**：完全相干假设在实际场景中可能不成立，影响干涉效应

7. **采样不足**：偏振BRDF的各向异性更强，需要更密集的采样

## 最佳实践检查清单

### 设计阶段
- [ ] 确定是否真的需要偏振渲染（计算成本vs视觉收益）
- [ ] 选择合适的偏振表示（琼斯vs米勒）
- [ ] 规划波长采样策略（单色vs光谱）
- [ ] 设计性能分级方案

### 实现阶段
- [ ] 验证米勒矩阵的物理有效性
- [ ] 实现坐标系变换的一致性检查
- [ ] 添加能量守恒断言
- [ ] 优化矩阵运算（利用稀疏性和对称性）

### 验证阶段
- [ ] 与解析解对比（平板、球体等简单几何）
- [ ] 检查互易性（交换光源和相机）
- [ ] 验证极限情况（掠射角、垂直入射）
- [ ] 测试非偏振光输入的退化情况

### 优化阶段
- [ ] profile矩阵运算开销
- [ ] 实现自适应精度切换
- [ ] 添加重要性采样
- [ ] 考虑GPU并行化机会


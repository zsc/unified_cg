# 第27章：量子光学基础

## 章节大纲

1. **开篇与学习目标**
2. **光的量子化**
   - 从经典电磁场到量子化场
   - 光子的产生与湮灭算符
   - Fock态与光子数态
   - 真空态与零点能
3. **相干态与压缩态**
   - 相干态的定义与性质
   - 位移算符与相干态生成
   - 压缩态与压缩算符
   - 最小不确定性态
4. **光子统计**
   - 光子计数分布
   - 泊松分布、超泊松与亚泊松光
   - Mandel Q参数
   - 光子聚束与反聚束
5. **二阶相干函数g^(2)**
   - 强度关联函数
   - Hanbury Brown-Twiss实验
   - g^(2)的物理意义
   - 经典与量子光的g^(2)界限
6. **量子噪声与散粒噪声**
   - 散粒噪声的量子起源
   - 标准量子极限
   - 压缩态降噪
   - 噪声谱密度
7. **本章小结**
8. **练习题**
9. **常见陷阱与错误**
10. **最佳实践检查清单**

---

## 开篇段落

本章介绍量子光学的基础概念，为理解现代光学现象和量子成像技术奠定基础。我们将从光场的量子化开始，探讨相干态、压缩态等量子光态的数学描述，深入研究光子统计特性和量子关联，最后讨论量子噪声的物理起源及其在精密测量中的影响。这些概念不仅对理解量子光学实验至关重要，也为下一章的量子成像与计算提供必要的理论工具。

### 学习目标

完成本章后，您将能够：
1. 推导光场的量子化过程，理解光子的粒子性描述
2. 计算相干态和压缩态的量子特性，包括不确定性关系
3. 分析不同光源的光子统计分布，区分经典光与非经典光
4. 使用二阶相干函数g^(2)表征光场的量子特性
5. 评估量子噪声对测量精度的影响，理解标准量子极限

---

## 27.1 光的量子化

### 27.1.1 从经典到量子

#### 经典电磁场的能量

经典电磁场的能量密度包含电场和磁场贡献：
$$u_{em} = \frac{1}{2}\left[\epsilon_0 E^2(r,t) + \frac{1}{\mu_0}B^2(r,t)\right]$$

总能量通过空间积分得到：
$$H_{classical} = \int d^3r \, u_{em} = \frac{1}{2}\int d^3r \left[\epsilon_0 E^2(r,t) + \frac{1}{\mu_0}B^2(r,t)\right]$$

这个能量是连续的，可以取任意值。然而，Planck的黑体辐射理论暗示电磁场能量应该是量子化的。

#### 模式展开

在有限体积$V$的谐振腔中，满足边界条件的电磁场可以展开为正交模式的叠加。对于立方腔（边长$L$），完美导体边界条件要求：
$$\mathbf{E}_{\parallel}|_{boundary} = 0, \quad \mathbf{B}_{\perp}|_{boundary} = 0$$

这导致离散的本征模式，由波矢$\mathbf{k} = \frac{2\pi}{L}(n_x, n_y, n_z)$标记，其中$n_i$为正整数。每个$\mathbf{k}$对应频率$\omega_k = c|\mathbf{k}|$。

#### 矢势的选择

在Coulomb规范下（$\nabla \cdot \mathbf{A} = 0$），标势$\Phi = 0$，电磁场完全由矢势描述。矢势展开式：
$$\mathbf{A}(r,t) = \sum_{k,s} \sqrt{\frac{\hbar}{2\epsilon_0\omega_k V}} \epsilon_{k,s} \left[a_{k,s}(t)e^{i\mathbf{k} \cdot \mathbf{r}} + a_{k,s}^*(t)e^{-i\mathbf{k} \cdot \mathbf{r}}\right]$$

其中：
- $s=1,2$标记两个正交偏振方向
- $\epsilon_{k,s}$是单位偏振矢量，满足$\epsilon_{k,s} \cdot \mathbf{k} = 0$（横波条件）
- $\epsilon_{k,1} \cdot \epsilon_{k,2} = 0$（正交性）
- $\epsilon_{k,1} \times \epsilon_{k,2} = \mathbf{k}/|\mathbf{k}|$（右手系）

归一化因子$\sqrt{\frac{\hbar}{2\epsilon_0\omega_k V}}$的选择将在量子化后变得清晰。

#### 电场和磁场

电场和磁场由矢势导出：
$$\mathbf{E} = -\frac{\partial \mathbf{A}}{\partial t}, \quad \mathbf{B} = \nabla \times \mathbf{A}$$

代入矢势展开式，并假设时间依赖为$a_{k,s}(t) = a_{k,s}e^{-i\omega_k t}$：
$$\mathbf{E}(r,t) = i\sum_{k,s} \sqrt{\frac{\hbar\omega_k}{2\epsilon_0 V}} \epsilon_{k,s} \left[a_{k,s}e^{i(\mathbf{k} \cdot \mathbf{r} - \omega_k t)} - a_{k,s}^*e^{-i(\mathbf{k} \cdot \mathbf{r} - \omega_k t)}\right]$$

$$\mathbf{B}(r,t) = i\sum_{k,s} \sqrt{\frac{\hbar}{2\epsilon_0\omega_k V}} (\mathbf{k} \times \epsilon_{k,s}) \left[a_{k,s}e^{i(\mathbf{k} \cdot \mathbf{r} - \omega_k t)} - a_{k,s}^*e^{-i(\mathbf{k} \cdot \mathbf{r} - \omega_k t)}\right]$$

#### 经典哈密顿量的模式表示

将场的模式展开代入经典哈密顿量，利用模式的正交性：
$$\int_V d^3r \, e^{i(\mathbf{k} - \mathbf{k}') \cdot \mathbf{r}} = V\delta_{\mathbf{k},\mathbf{k}'}$$

$$\epsilon_{k,s} \cdot \epsilon_{k,s'} = \delta_{s,s'}, \quad (\mathbf{k} \times \epsilon_{k,s}) \cdot (\mathbf{k} \times \epsilon_{k,s'}) = k^2\delta_{s,s'}$$

经过计算得到：
$$H_{classical} = \sum_{k,s} \hbar\omega_k |a_{k,s}|^2$$

这正是无穷多个谐振子的能量之和，每个模式$(k,s)$对应一个频率为$\omega_k$的谐振子。

### 27.1.2 量子化过程

#### 正则变量的识别

正则量子化要求识别共轭变量对。从经典谐振子类比，对于每个模式$(k,s)$，我们定义：
- 广义坐标：$q_{k,s} = \frac{1}{\sqrt{2}}(a_{k,s} + a_{k,s}^*)$
- 广义动量：$p_{k,s} = \frac{i}{\sqrt{2\omega_k}}(a_{k,s}^* - a_{k,s})$

这些变量是实数，满足经典泊松括号：
$$\{q_{k,s}, p_{k',s'}\} = \delta_{k,k'}\delta_{s,s'}$$

可以验证经典哈密顿量表示为：
$$H_{classical} = \sum_{k,s} \frac{1}{2}[\omega_k^2 q_{k,s}^2 + p_{k,s}^2]$$

这确认了每个模式确实是一个谐振子。

#### 正则量子化

量子化时，将泊松括号替换为对易子：
$$[\hat{q}_{k,s}, \hat{p}_{k',s'}] = i\hbar\delta_{k,k'}\delta_{s,s'}$$

所有其他对易子为零：
$$[\hat{q}_{k,s}, \hat{q}_{k',s'}] = [\hat{p}_{k,s}, \hat{p}_{k',s'}] = 0$$

#### 产生湮灭算符

引入算符：
$$\hat{a}_{k,s} = \sqrt{\frac{\omega_k}{2\hbar}}\hat{q}_{k,s} + \frac{i}{\sqrt{2\hbar\omega_k}}\hat{p}_{k,s}$$
$$\hat{a}_{k,s}^\dagger = \sqrt{\frac{\omega_k}{2\hbar}}\hat{q}_{k,s} - \frac{i}{\sqrt{2\hbar\omega_k}}\hat{p}_{k,s}$$

反过来：
$$\hat{q}_{k,s} = \sqrt{\frac{\hbar}{2\omega_k}}(\hat{a}_{k,s} + \hat{a}_{k,s}^\dagger)$$
$$\hat{p}_{k,s} = i\sqrt{\frac{\hbar\omega_k}{2}}(\hat{a}_{k,s}^\dagger - \hat{a}_{k,s})$$

#### 对易关系

利用$[\hat{q}, \hat{p}] = i\hbar$，可以导出：
$$[\hat{a}_{k,s}, \hat{a}_{k',s'}^\dagger] = \delta_{k,k'}\delta_{s,s'}$$
$$[\hat{a}_{k,s}, \hat{a}_{k',s'}] = [\hat{a}_{k,s}^\dagger, \hat{a}_{k',s'}^\dagger] = 0$$

这些是玻色子的标准对易关系。

#### 量子哈密顿量

将正则变量的算符表示代入哈密顿量：
$$\hat{H} = \sum_{k,s} \frac{1}{2}[\omega_k^2 \hat{q}_{k,s}^2 + \hat{p}_{k,s}^2]$$

使用产生湮灭算符表示：
$$\hat{q}_{k,s}^2 = \frac{\hbar}{2\omega_k}(\hat{a}_{k,s} + \hat{a}_{k,s}^\dagger)^2$$
$$\hat{p}_{k,s}^2 = -\frac{\hbar\omega_k}{2}(\hat{a}_{k,s} - \hat{a}_{k,s}^\dagger)^2$$

展开并利用对易关系$[\hat{a}, \hat{a}^\dagger] = 1$：
$$\hat{H}_{k,s} = \frac{\hbar\omega_k}{4}[(\hat{a}_{k,s} + \hat{a}_{k,s}^\dagger)^2 - (\hat{a}_{k,s} - \hat{a}_{k,s}^\dagger)^2]$$
$$= \frac{\hbar\omega_k}{2}[\hat{a}_{k,s}\hat{a}_{k,s}^\dagger + \hat{a}_{k,s}^\dagger\hat{a}_{k,s}]$$
$$= \hbar\omega_k\left(\hat{a}_{k,s}^\dagger\hat{a}_{k,s} + \frac{1}{2}\right)$$

总哈密顿量：
$$\hat{H} = \sum_{k,s} \hbar\omega_k\left(\hat{a}_{k,s}^\dagger\hat{a}_{k,s} + \frac{1}{2}\right)$$

#### 光子数算符

定义光子数算符：
$$\hat{n}_{k,s} = \hat{a}_{k,s}^\dagger\hat{a}_{k,s}$$

哈密顿量简化为：
$$\hat{H} = \sum_{k,s} \hbar\omega_k\left(\hat{n}_{k,s} + \frac{1}{2}\right)$$

每个模式的能量量子化为$\hbar\omega_k$的整数倍，加上零点能$\hbar\omega_k/2$。

### 27.1.3 Fock态

#### 真空态

定义真空态$|0\rangle$为所有模式的基态：
$$\hat{a}_{k,s}|0\rangle = 0, \quad \forall k,s$$

真空态是所有模式都处于最低能量状态的态。它满足：
$$\hat{n}_{k,s}|0\rangle = 0, \quad \forall k,s$$

但真空能量不为零：
$$E_{vacuum} = \langle 0|\hat{H}|0\rangle = \sum_{k,s} \frac{\hbar\omega_k}{2}$$

#### 单模Fock态

对于单个模式$(k,s)$，Fock态（光子数态）通过反复作用产生算符构造：
$$|n\rangle_{k,s} = \frac{(\hat{a}_{k,s}^\dagger)^n}{\sqrt{n!}}|0\rangle$$

归一化因子$1/\sqrt{n!}$确保$\langle n|n\rangle = 1$。可以通过数学归纳法证明：
$$\hat{a}_{k,s}^\dagger|n\rangle_{k,s} = \sqrt{n+1}|n+1\rangle_{k,s}$$
$$\hat{a}_{k,s}|n\rangle_{k,s} = \sqrt{n}|n-1\rangle_{k,s}$$

这些关系可记忆为：
- $\hat{a}^\dagger$产生一个光子，系数$\sqrt{n+1}$反映了玻色增强
- $\hat{a}$湮灭一个光子，系数$\sqrt{n}$确保$\hat{a}|0\rangle = 0$

#### 多模Fock态

一般的多模Fock态：
$$|\{n_{k,s}\}\rangle = \prod_{k,s} \frac{(\hat{a}_{k,s}^\dagger)^{n_{k,s}}}{\sqrt{n_{k,s}!}}|0\rangle$$

简记为$|n_1, n_2, ...\rangle$，其中每个$n_i$表示模式$i$中的光子数。

#### 算符作用规则

光子数算符的本征方程：
$$\hat{n}_{k,s}|n_{k,s}\rangle = n_{k,s}|n_{k,s}\rangle$$

利用$\hat{n} = \hat{a}^\dagger\hat{a}$和对易关系，可以导出：
$$\hat{n}\hat{a}^\dagger|n\rangle = \hat{a}^\dagger(\hat{n}+1)|n\rangle = (n+1)\hat{a}^\dagger|n\rangle$$

这确认了$\hat{a}^\dagger|n\rangle \propto |n+1\rangle$。

#### 能量本征值

Fock态是哈密顿量的本征态：
$$\hat{H}|\{n_{k,s}\}\rangle = E_{\{n_{k,s}\}}|\{n_{k,s}\}\rangle$$

其中能量本征值：
$$E_{\{n_{k,s}\}} = \sum_{k,s} \hbar\omega_k\left(n_{k,s} + \frac{1}{2}\right)$$

能级间隔：
- 单光子激发：$\Delta E = \hbar\omega_k$
- 多光子激发：$\Delta E = \sum_{k,s} m_{k,s}\hbar\omega_k$

#### Fock态的正交完备性

正交性：
$$\langle\{n_{k,s}\}|\{n'_{k,s}\}\rangle = \prod_{k,s} \delta_{n_{k,s},n'_{k,s}}$$

完备性：
$$\sum_{\{n_{k,s}\}} |\{n_{k,s}\}\rangle\langle\{n_{k,s}\}| = \mathbb{I}$$

任意态可展开：
$$|\psi\rangle = \sum_{\{n_{k,s}\}} c_{\{n_{k,s}\}}|\{n_{k,s}\}\rangle$$

其中$c_{\{n_{k,s}\}} = \langle\{n_{k,s}\}|\psi\rangle$是概率幅。

### 27.1.4 真空涨落

#### 零点能问题

真空态虽然没有光子，但具有非零能量：
$$E_0 = \langle 0|\hat{H}|0\rangle = \sum_{k,s} \frac{\hbar\omega_k}{2}$$

这个和发散，因为模式数无穷。处理方法：
1. **正规序**：重定义哈密顿量$:\hat{H}: = \sum_{k,s} \hbar\omega_k\hat{a}_{k,s}^\dagger\hat{a}_{k,s}$
2. **物理截断**：认识到实际系统有最高频率$\omega_{max}$
3. **重整化**：只有能量差有物理意义

#### 场的真空期望值

电场算符的真空期望值为零：
$$\langle 0|\hat{\mathbf{E}}(r,t)|0\rangle = 0$$

这是因为$\hat{\mathbf{E}} \propto (\hat{a} - \hat{a}^\dagger)$，而$\langle 0|\hat{a}|0\rangle = \langle 0|\hat{a}^\dagger|0\rangle = 0$。

但均方涨落非零。对于单模：
$$\langle 0|\hat{E}_{k,s}^2|0\rangle = \frac{\hbar\omega_k}{2\epsilon_0 V}$$

总的真空涨落：
$$\langle 0|\hat{\mathbf{E}}^2(r,t)|0\rangle = \sum_{k,s} \frac{\hbar\omega_k}{2\epsilon_0 V}|\epsilon_{k,s}|^2$$

#### 真空涨落的物理图像

真空涨落可理解为：
- **时间-能量不确定性**：$\Delta E \Delta t \geq \hbar/2$允许短时间内能量涨落
- **虚粒子**：光子可以短暂出现又消失，只要满足不确定性关系
- **量子零点运动**：类似谐振子的零点振动

#### 可观测的真空效应

1. **Casimir效应**
   
   两平行导体板（间距$d$）改变了允许的模式：
   $$k_z = \frac{n\pi}{d}, \quad n = 1,2,3,...$$
   
   真空能依赖于$d$：
   $$E_{Casimir}(d) = -\frac{\pi^2\hbar c}{720d^3}A$$
   
   产生吸引力：
   $$F = -\frac{\partial E}{\partial d} = -\frac{\pi^2\hbar c}{240d^4}A$$

2. **Lamb位移**
   
   原子能级因与真空场耦合而移动。对于氢原子2S₁/₂和2P₁/₂能级：
   $$\Delta E_{Lamb} \approx 1057 \text{ MHz}$$
   
   主要贡献来自电子位置算符与真空电场的二阶微扰。

3. **自发辐射**
   
   激发态原子的衰减率（Einstein A系数）：
   $$A_{if} = \frac{\omega_{if}^3}{3\pi\epsilon_0\hbar c^3}|\langle f|\hat{\mathbf{d}}|i\rangle|^2$$
   
   其中$\hat{\mathbf{d}}$是电偶极矩算符。这可以理解为真空涨落诱导的跃迁。

#### 真空涨落的实验验证

1. **Casimir力的测量**：使用原子力显微镜可以精确测量
2. **动态Casimir效应**：快速移动的镜子可以从真空中产生真实光子
3. **量子电动力学的精密测试**：Lamb位移的精确测量验证了QED理论

### 27.1.5 场的正交分量

#### 正交振幅算符

定义单模场的正交振幅算符：
$$\hat{X} = \frac{1}{\sqrt{2}}(\hat{a} + \hat{a}^\dagger), \quad \hat{P} = \frac{i}{\sqrt{2}}(\hat{a}^\dagger - \hat{a})$$

这些算符的物理意义：
- $\hat{X}$：场的"位置"分量，对应于$\cos$相位
- $\hat{P}$：场的"动量"分量，对应于$\sin$相位

反演关系：
$$\hat{a} = \frac{1}{\sqrt{2}}(\hat{X} - i\hat{P}), \quad \hat{a}^\dagger = \frac{1}{\sqrt{2}}(\hat{X} + i\hat{P})$$

#### 对易关系和不确定性

这些算符满足：
$$[\hat{X}, \hat{P}] = i$$

导出过程：
$$[\hat{X}, \hat{P}] = \frac{i}{2}[(\hat{a} + \hat{a}^\dagger), (\hat{a}^\dagger - \hat{a})]$$
$$= \frac{i}{2}([\hat{a}, \hat{a}^\dagger] - [\hat{a}^\dagger, \hat{a}]) = \frac{i}{2} \cdot 2 = i$$

对应的不确定性关系：
$$\Delta X \Delta P \geq \frac{1}{2}$$

#### 真空态的正交分量

在真空态中：
$$\langle 0|\hat{X}|0\rangle = \langle 0|\hat{P}|0\rangle = 0$$

方差：
$$\langle 0|\hat{X}^2|0\rangle = \frac{1}{2}\langle 0|(\hat{a} + \hat{a}^\dagger)^2|0\rangle = \frac{1}{2}\langle 0|\hat{a}\hat{a}^\dagger + \hat{a}^\dagger\hat{a}|0\rangle = \frac{1}{2}$$

类似地：
$$\langle 0|\hat{P}^2|0\rangle = \frac{1}{2}$$

因此：$\Delta X = \Delta P = \frac{1}{\sqrt{2}}$，乘积$\Delta X \Delta P = \frac{1}{2}$达到最小不确定性。

#### 相位空间表示

在$(X,P)$相位空间中：
- 真空态：以原点为中心的圆形高斯分布
- 不确定性圆：半径$\sim 1/\sqrt{2}$
- 面积：$\pi \Delta X \Delta P = \pi/2$（最小相空间面积）

#### 与经典场的联系

经典相干场可写为：
$$E(t) = E_0\cos(\omega t + \phi) = E_0[\cos\phi\cos(\omega t) - \sin\phi\sin(\omega t)]$$

定义：
- $X_{cl} = E_0\cos\phi$（同相分量）
- $P_{cl} = E_0\sin\phi$（正交分量）

量子对应：
$$\hat{E}(t) \propto \hat{X}\cos(\omega t) - \hat{P}\sin(\omega t)$$

#### 旋转的正交分量

更一般地，可定义任意相位的正交分量：
$$\hat{X}_\theta = \frac{1}{\sqrt{2}}(\hat{a}e^{-i\theta} + \hat{a}^\dagger e^{i\theta})$$
$$\hat{P}_\theta = \frac{i}{\sqrt{2}}(\hat{a}^\dagger e^{i\theta} - \hat{a}e^{-i\theta})$$

这相当于在相位空间中旋转角度$\theta$：
$$\begin{pmatrix} \hat{X}_\theta \\ \hat{P}_\theta \end{pmatrix} = \begin{pmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{pmatrix} \begin{pmatrix} \hat{X} \\ \hat{P} \end{pmatrix}$$

对易关系保持不变：$[\hat{X}_\theta, \hat{P}_\theta] = i$。

---

## 27.2 相干态与压缩态

### 27.2.1 相干态定义

相干态代表了量子光学中最接近经典光的量子态。它们是湮灭算符的本征态：
$$\hat{a}|\alpha\rangle = \alpha|\alpha\rangle$$

其中$\alpha = |\alpha|e^{i\phi}$是复数本征值，$|\alpha|$对应经典振幅，$\phi$对应相位。

#### 历史背景

相干态的概念：
- 1926年：Schrödinger首次研究谐振子的"最小不确定性波包"
- 1963年：Glauber引入光场的相干态，奠定量子光学基础
- 1964年：Sudarshan证明了任意密度算符可用相干态表示（P表示）

#### Fock基展开

在光子数基下，相干态表示为：
$$|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^{\infty}\frac{\alpha^n}{\sqrt{n!}}|n\rangle$$

推导过程：设$|\alpha\rangle = \sum_n c_n|n\rangle$，代入本征方程：
$$\hat{a}|\alpha\rangle = \sum_n c_n\sqrt{n}|n-1\rangle = \alpha\sum_n c_n|n\rangle$$

比较$|n\rangle$的系数：
$$c_n\sqrt{n+1} = \alpha c_{n+1}$$

递推关系：$c_{n+1} = \frac{\alpha}{\sqrt{n+1}}c_n = \frac{\alpha^{n+1}}{\sqrt{(n+1)!}}c_0$

归一化条件：
$$1 = \langle\alpha|\alpha\rangle = |c_0|^2\sum_n \frac{|\alpha|^{2n}}{n!} = |c_0|^2 e^{|\alpha|^2}$$

因此$c_0 = e^{-|\alpha|^2/2}$（选择相位为0）。

#### 相干态的非正交性

两个相干态的内积：
$$\langle\alpha|\beta\rangle = e^{-\frac{|\alpha|^2+|\beta|^2}{2}}\sum_n \frac{(\alpha^*)^n\beta^n}{n!} = e^{-\frac{|\alpha|^2+|\beta|^2}{2}}e^{\alpha^*\beta}$$

简化为：
$$\langle\alpha|\beta\rangle = e^{-\frac{1}{2}|\alpha-\beta|^2}e^{i\text{Im}(\alpha^*\beta)}$$

物理意义：
- 重叠随距离$|\alpha-\beta|$呈高斯衰减
- 当$|\alpha-\beta| > 3$时，$|\langle\alpha|\beta\rangle|^2 < 10^{-4}$（近似正交）
- 相位因子$e^{i\text{Im}(\alpha^*\beta)}$反映了量子相位

#### 过完备性

相干态构成过完备基：
$$\frac{1}{\pi}\int d^2\alpha |\alpha\rangle\langle\alpha| = \mathbb{I}$$

证明：利用Fock基展开
$$\frac{1}{\pi}\int d^2\alpha |\alpha\rangle\langle\alpha| = \frac{1}{\pi}\int d^2\alpha e^{-|\alpha|^2}\sum_{n,m}\frac{\alpha^n(\alpha^*)^m}{\sqrt{n!m!}}|n\rangle\langle m|$$

极坐标下$\alpha = re^{i\theta}$：
$$\int_0^{2\pi}d\theta e^{i(n-m)\theta} = 2\pi\delta_{nm}$$

径向积分：
$$\int_0^\infty r dr \, r^{2n}e^{-r^2} = \frac{n!}{2}$$

因此：
$$\frac{1}{\pi}\int d^2\alpha |\alpha\rangle\langle\alpha| = \sum_n |n\rangle\langle n| = \mathbb{I}$$

#### 相干态的生成

实验上产生相干态的方法：
1. **激光输出**：理想单模激光产生相干态
2. **强衰减相干光**：$|\alpha| \ll 1$时近似单光子源
3. **位移真空态**：$|\alpha\rangle = \hat{D}(\alpha)|0\rangle$

### 27.2.2 相干态性质

#### 光子统计

1. **平均光子数**：
   $$\langle\hat{n}\rangle = \langle\alpha|\hat{a}^\dagger\hat{a}|\alpha\rangle = |\alpha|^2$$
   
   推导：利用$\hat{a}|\alpha\rangle = \alpha|\alpha\rangle$
   $$\langle\hat{n}\rangle = \langle\alpha|\hat{a}^\dagger\hat{a}|\alpha\rangle = \alpha^*\alpha\langle\alpha|\alpha\rangle = |\alpha|^2$$

2. **高阶矩**：
   $$\langle\hat{n}^k\rangle = \langle\alpha|(\hat{a}^\dagger\hat{a})^k|\alpha\rangle$$
   
   利用正规序和$\hat{a}|\alpha\rangle = \alpha|\alpha\rangle$：
   $$\langle:\hat{n}^k:\rangle = |\alpha|^{2k}$$
   
   但实际的$k$阶矩包含了正规序修正。

3. **光子数方差**：
   $$\langle\hat{n}^2\rangle = \langle\alpha|\hat{n}(\hat{n}-1)+\hat{n}|\alpha\rangle = |\alpha|^4 + |\alpha|^2$$
   
   因此：
   $$(\Delta n)^2 = \langle\hat{n}^2\rangle - \langle\hat{n}\rangle^2 = |\alpha|^2 = \langle\hat{n}\rangle$$
   
   这是泊松统计的特征：方差等于平均值。

4. **光子数分布**（泊松分布）：
   $$P(n) = |\langle n|\alpha\rangle|^2 = \frac{|\alpha|^{2n}}{n!}e^{-|\alpha|^2} = \frac{\bar{n}^n}{n!}e^{-\bar{n}}$$
   
   其中$\bar{n} = |\alpha|^2$是平均光子数。

#### 场的期望值

电场算符（单模）：
$$\hat{E}(r,t) = i\sqrt{\frac{\hbar\omega}{2\epsilon_0 V}}\epsilon\left[\hat{a}e^{i(k \cdot r - \omega t)} - \hat{a}^\dagger e^{-i(k \cdot r - \omega t)}\right]$$

期望值：
$$\langle\alpha|\hat{E}(r,t)|\alpha\rangle = i\sqrt{\frac{\hbar\omega}{2\epsilon_0 V}}\epsilon\left[\alpha e^{i(k \cdot r - \omega t)} - \alpha^* e^{-i(k \cdot r - \omega t)}\right]$$

设$\alpha = |\alpha|e^{i\phi}$：
$$\langle\hat{E}(r,t)\rangle = 2\sqrt{\frac{\hbar\omega}{2\epsilon_0 V}}|\alpha|\epsilon\sin(k \cdot r - \omega t + \phi)$$

这正是经典相干波，振幅$E_0 = 2|\alpha|\sqrt{\frac{\hbar\omega}{2\epsilon_0 V}}$。

#### 相位空间表示

正交分量的期望值：
$$\langle\hat{X}\rangle = \frac{1}{\sqrt{2}}\langle\alpha|\hat{a} + \hat{a}^\dagger|\alpha\rangle = \frac{1}{\sqrt{2}}(\alpha + \alpha^*) = \sqrt{2}\text{Re}(\alpha)$$
$$\langle\hat{P}\rangle = \frac{i}{\sqrt{2}}\langle\alpha|\hat{a}^\dagger - \hat{a}|\alpha\rangle = \frac{i}{\sqrt{2}}(\alpha^* - \alpha) = \sqrt{2}\text{Im}(\alpha)$$

方差：
$$(\Delta X)^2 = \langle\hat{X}^2\rangle - \langle\hat{X}\rangle^2 = \frac{1}{2}$$
$$(\Delta P)^2 = \langle\hat{P}^2\rangle - \langle\hat{P}\rangle^2 = \frac{1}{2}$$

相干态特征：
- 最小不确定性：$\Delta X \Delta P = \frac{1}{2}$
- 各向同性：$\Delta X = \Delta P$
- 高斯分布：Wigner函数为高斯形

#### 时间演化

在自由演化下：
$$|\alpha(t)\rangle = e^{-i\hat{H}t/\hbar}|\alpha(0)\rangle = e^{-i\omega t/2}|\alpha(0)e^{-i\omega t}\rangle$$

相干态保持相干态，只是复振幅旋转：
$$\alpha(t) = \alpha(0)e^{-i\omega t}$$

这对应经典谐振子的运动。

#### 相干态的准经典性

相干态被称为"准经典态"的原因：
1. **最小不确定性**：达到海森堡极限
2. **泊松光子统计**：类似经典随机过程
3. **场的经典行为**：期望值满足经典方程
4. **相位空间局域化**：Wigner函数为正定高斯
5. **动力学对应**：演化遵循经典轨迹

### 27.2.3 位移算符

位移算符定义：
$$\hat{D}(\alpha) = \exp(\alpha\hat{a}^\dagger - \alpha^*\hat{a})$$

#### 物理意义

位移算符在相位空间中平移量子态：
- 将真空态变为相干态
- 实现相位空间的平移变换
- 对应经典场的叠加

#### Baker-Campbell-Hausdorff公式

对于两个算符$\hat{A}$和$\hat{B}$，若$[\hat{A}, \hat{B}] = c$（c数），则：
$$e^{\hat{A}+\hat{B}} = e^{-c/2}e^{\hat{A}}e^{\hat{B}} = e^{c/2}e^{\hat{B}}e^{\hat{A}}$$

应用于位移算符，设$\hat{A} = \alpha\hat{a}^\dagger$，$\hat{B} = -\alpha^*\hat{a}$：
$$[\hat{A}, \hat{B}] = \alpha(-\alpha^*)[\hat{a}^\dagger, \hat{a}] = |\alpha|^2$$

因此：
$$\hat{D}(\alpha) = e^{-|\alpha|^2/2}e^{\alpha\hat{a}^\dagger}e^{-\alpha^*\hat{a}} = e^{|\alpha|^2/2}e^{-\alpha^*\hat{a}}e^{\alpha\hat{a}^\dagger}$$

#### 位移算符的基本性质

1. **幺正性**：
   $$\hat{D}^\dagger(\alpha) = \exp(\alpha^*\hat{a} - \alpha\hat{a}^\dagger) = \hat{D}(-\alpha)$$
   
   因此：$\hat{D}^\dagger(\alpha)\hat{D}(\alpha) = \hat{D}(-\alpha)\hat{D}(\alpha) = \mathbb{I}$

2. **群性质**：
   $$\hat{D}(\alpha)\hat{D}(\beta) = e^{i\text{Im}(\alpha^*\beta)}\hat{D}(\alpha+\beta)$$
   
   证明：使用BCH公式和$[\alpha\hat{a}^\dagger - \alpha^*\hat{a}, \beta\hat{a}^\dagger - \beta^*\hat{a}] = 2i\text{Im}(\alpha^*\beta)$

3. **变换性质**：
   
   利用$e^{\hat{A}}\hat{B}e^{-\hat{A}} = \hat{B} + [\hat{A}, \hat{B}] + \frac{1}{2!}[\hat{A}, [\hat{A}, \hat{B}]] + ...$
   
   对于$\hat{A} = \alpha\hat{a}^\dagger - \alpha^*\hat{a}$：
   $$[\hat{A}, \hat{a}] = -\alpha, \quad [\hat{A}, \hat{a}^\dagger] = \alpha^*$$
   
   高阶对易子为零，因此：
   $$\hat{D}^\dagger(\alpha)\hat{a}\hat{D}(\alpha) = \hat{a} + \alpha$$
   $$\hat{D}^\dagger(\alpha)\hat{a}^\dagger\hat{D}(\alpha) = \hat{a}^\dagger + \alpha^*$$

4. **生成相干态**：
   $$|\alpha\rangle = \hat{D}(\alpha)|0\rangle$$
   
   验证：
   $$\hat{a}|\alpha\rangle = \hat{a}\hat{D}(\alpha)|0\rangle = \hat{D}(\alpha)[\hat{D}^\dagger(\alpha)\hat{a}\hat{D}(\alpha)]|0\rangle$$
   $$= \hat{D}(\alpha)(\hat{a} + \alpha)|0\rangle = \alpha\hat{D}(\alpha)|0\rangle = \alpha|\alpha\rangle$$

#### 位移算符的表示

1. **Fock基表示**：
   $$\langle n|\hat{D}(\alpha)|m\rangle = \sqrt{\frac{m!}{n!}}e^{-|\alpha|^2/2}\alpha^{n-m}L_m^{n-m}(|\alpha|^2)$$
   
   其中$L_m^{n-m}$是关联Laguerre多项式（$n \geq m$时）。

2. **相干态表示**：
   $$\hat{D}(\alpha) = \int \frac{d^2\beta}{\pi}|\beta+\alpha\rangle\langle\beta|$$

#### 实验实现

位移操作的实验方法：
1. **经典场注入**：将弱相干态与强局域振荡器混合
2. **参量放大器**：使用简并参量下转换
3. **电光调制**：通过Pockels效应实现相位空间位移

### 27.2.4 压缩态

压缩算符定义：
$$\hat{S}(\xi) = \exp\left[\frac{1}{2}(\xi^*\hat{a}^2 - \xi\hat{a}^{\dagger 2})\right]$$

其中$\xi = re^{i\theta}$，$r \geq 0$是压缩强度，$\theta$是压缩方向。

#### 物理起源

压缩态的产生机制：
1. **参量下转换**：二阶非线性过程$\chi^{(2)}$
2. **四波混频**：三阶非线性过程$\chi^{(3)}$
3. **原子系综**：集体自旋压缩
4. **光机械系统**：辐射压力耦合

#### 压缩算符的李代数结构

定义$SU(1,1)$生成元：
$$\hat{K}_+ = \frac{1}{2}\hat{a}^{\dagger 2}, \quad \hat{K}_- = \frac{1}{2}\hat{a}^2, \quad \hat{K}_0 = \frac{1}{2}(\hat{a}^\dagger\hat{a} + \frac{1}{2})$$

对易关系：
$$[\hat{K}_-, \hat{K}_+] = 2\hat{K}_0, \quad [\hat{K}_0, \hat{K}_\pm] = \pm\hat{K}_\pm$$

压缩算符表示为：
$$\hat{S}(\xi) = \exp(\xi^*\hat{K}_- - \xi\hat{K}_+)$$

#### 压缩算符的分解

利用$SU(1,1)$的BCH公式：
$$\hat{S}(\xi) = \exp\left[\frac{\tanh r}{2}e^{-i\theta}\hat{a}^{\dagger 2}\right]\exp\left[-\ln(\cosh r)(\hat{a}^\dagger\hat{a} + \frac{1}{2})\right]\exp\left[-\frac{\tanh r}{2}e^{i\theta}\hat{a}^2\right]$$

这个分解对计算矩阵元很有用。

#### 压缩真空态

压缩真空态：
$$|\xi\rangle = \hat{S}(\xi)|0\rangle$$

Fock基展开：
$$|\xi\rangle = \frac{1}{\sqrt{\cosh r}}\sum_{n=0}^{\infty}\frac{\sqrt{(2n)!}}{2^n n!}(-e^{i\theta}\tanh r)^n|2n\rangle$$

推导要点：
- $\hat{a}^2|0\rangle = 0$使得$\exp[-\frac{\tanh r}{2}e^{i\theta}\hat{a}^2]|0\rangle = |0\rangle$
- 中间项贡献因子$(\cosh r)^{-1/2}$
- $\hat{a}^{\dagger 2}$只产生偶数光子态

物理特征：
1. **双光子关联**：只有偶数光子态
2. **平均光子数**：$\langle\hat{n}\rangle = \sinh^2 r$
3. **光子数方差**：$(\Delta n)^2 = 2\sinh^2 r \cosh 2r$（超泊松）

#### 压缩相干态

一般的压缩相干态：
$$|\alpha,\xi\rangle = \hat{D}(\alpha)\hat{S}(\xi)|0\rangle$$

注意算符顺序。另一种定义：
$$|\alpha,\xi\rangle' = \hat{S}(\xi)\hat{D}(\beta)|0\rangle$$

其中$\beta$与$\alpha$的关系由Bogoliubov变换确定。

#### 压缩变换下的算符

Bogoliubov变换：
$$\hat{S}^\dagger(\xi)\hat{a}\hat{S}(\xi) = \mu\hat{a} - \nu\hat{a}^\dagger$$
$$\hat{S}^\dagger(\xi)\hat{a}^\dagger\hat{S}(\xi) = \mu^*\hat{a}^\dagger - \nu^*\hat{a}$$

其中：
$$\mu = \cosh r, \quad \nu = e^{i\theta}\sinh r$$

满足$|\mu|^2 - |\nu|^2 = 1$（保持对易关系）。

#### 双模压缩

双模压缩算符：
$$\hat{S}_{12}(\xi) = \exp[\xi^*\hat{a}_1\hat{a}_2 - \xi\hat{a}_1^\dagger\hat{a}_2^\dagger]$$

双模压缩真空态（EPR态）：
$$|\xi\rangle_{12} = \frac{1}{\cosh r}\sum_{n=0}^{\infty}(-e^{i\theta}\tanh r)^n|n\rangle_1|n\rangle_2$$

特性：
- 完美的光子数关联：$\langle\hat{n}_1\hat{n}_2\rangle - \langle\hat{n}_1\rangle\langle\hat{n}_2\rangle = \sinh^2 r \cosh^2 r$
- 连续变量纠缠：正交分量的EPR关联
- 量子通信资源：用于量子隐形传态和密集编码

### 27.2.5 压缩态的不确定性

#### 正交分量的变换

在压缩变换下：
$$\hat{S}^\dagger(\xi)\hat{a}\hat{S}(\xi) = \cosh r \cdot \hat{a} - e^{i\theta}\sinh r \cdot \hat{a}^\dagger$$

定义旋转的正交分量：
$$\hat{X}_\phi = \frac{1}{\sqrt{2}}(\hat{a}e^{-i\phi} + \hat{a}^\dagger e^{i\phi})$$
$$\hat{P}_\phi = \frac{i}{\sqrt{2}}(\hat{a}^\dagger e^{i\phi} - \hat{a}e^{-i\phi})$$

对于压缩真空态，当$\phi = \theta/2$时：
$$(\Delta X_{\theta/2})^2 = \frac{1}{2}e^{-2r}, \quad (\Delta P_{\theta/2})^2 = \frac{1}{2}e^{2r}$$

不确定性乘积：$\Delta X_{\theta/2} \Delta P_{\theta/2} = \frac{1}{2}$（最小不确定性）

#### 压缩的物理意义

1. **噪声重分配**：压缩不减少总噪声，而是将噪声从一个正交分量转移到另一个
2. **量子优势**：在压缩方向上，噪声低于真空涨落（散粒噪声极限）
3. **应用**：
   - 引力波探测（LIGO使用压缩光提高灵敏度）
   - 量子密钥分发（连续变量QKD）
   - 超分辨成像

#### 双模压缩

双模压缩算符：
$$\hat{S}_{12}(\xi) = \exp[\xi^*\hat{a}_1\hat{a}_2 - \xi\hat{a}_1^\dagger\hat{a}_2^\dagger]$$

产生纠缠的双模压缩真空态：
$$|\xi\rangle_{12} = \frac{1}{\cosh r}\sum_{n=0}^{\infty}(-e^{i\theta}\tanh r)^n|n\rangle_1|n\rangle_2$$

这是Einstein-Podolsky-Rosen (EPR)态的一种实现。

---

## 27.3 光子统计

### 27.3.1 光子计数分布

光子统计描述了在给定时间窗口内探测到特定数目光子的概率。对于一般的量子态$\hat{\rho}$，光子数分布为：
$$P(n) = \langle n|\hat{\rho}|n\rangle = \text{Tr}(\hat{\rho}|n\rangle\langle n|)$$

#### 矩和累积量

光子数的各阶矩：
$$\langle n^k\rangle = \text{Tr}(\hat{\rho}\hat{n}^k) = \sum_{n=0}^{\infty} n^k P(n)$$

特别重要的是前两阶矩：
- 一阶矩（平均值）：$\langle n\rangle = \text{Tr}(\hat{\rho}\hat{a}^\dagger\hat{a})$
- 二阶矩：$\langle n^2\rangle = \text{Tr}(\hat{\rho}\hat{a}^\dagger\hat{a}\hat{a}^\dagger\hat{a})$

利用对易关系$[\hat{a}, \hat{a}^\dagger] = 1$：
$$\langle n^2\rangle = \langle\hat{a}^\dagger\hat{a}^\dagger\hat{a}\hat{a}\rangle + \langle\hat{n}\rangle$$

#### 阶乘矩

阶乘矩在光子统计中特别有用：
$$\langle n^{(k)}\rangle = \langle n(n-1)...(n-k+1)\rangle = \langle\hat{a}^{\dagger k}\hat{a}^k\rangle$$

例如：
- $\langle n^{(1)}\rangle = \langle n\rangle$
- $\langle n^{(2)}\rangle = \langle n(n-1)\rangle = \langle\hat{a}^\dagger\hat{a}^\dagger\hat{a}\hat{a}\rangle$

### 27.3.2 典型光源的统计分布

#### 相干光（泊松分布）

对于相干态$|\alpha\rangle$：
$$P(n) = \frac{|\alpha|^{2n}}{n!}e^{-|\alpha|^2} = \frac{\bar{n}^n}{n!}e^{-\bar{n}}$$

其中$\bar{n} = |\alpha|^2$。

统计特性：
- 平均值：$\langle n\rangle = \bar{n}$
- 方差：$(\Delta n)^2 = \bar{n}$
- 所有阶乘矩：$\langle n^{(k)}\rangle = \bar{n}^k$

#### 热光（玻色-爱因斯坦分布）

热平衡光场的光子数分布：
$$P(n) = \frac{\bar{n}^n}{(1+\bar{n})^{n+1}}$$

这是几何分布，源于黑体辐射的量子统计。

统计特性：
- 平均值：$\langle n\rangle = \bar{n}$
- 方差：$(\Delta n)^2 = \bar{n} + \bar{n}^2$
- 二阶阶乘矩：$\langle n^{(2)}\rangle = 2\bar{n}^2$

#### Fock态

纯光子数态$|m\rangle$的分布：
$$P(n) = \delta_{n,m}$$

统计特性：
- 平均值：$\langle n\rangle = m$
- 方差：$(\Delta n)^2 = 0$（无涨落）
- 这是唯一具有确定光子数的态

#### 压缩态

压缩真空态的光子数分布：
$$P(n) = \begin{cases}
\frac{(2m)!}{2^{2m}(m!)^2}\frac{(\tanh r)^{2m}}{\cosh r} & n = 2m \\
0 & n = 2m+1
\end{cases}$$

统计特性：
- 平均值：$\langle n\rangle = \sinh^2 r$
- 方差：$(\Delta n)^2 = 2\sinh^2 r(\cosh 2r)$
- 表现出超泊松统计

### 27.3.3 统计参数

#### Mandel Q参数

定义：
$$Q = \frac{(\Delta n)^2 - \langle n\rangle}{\langle n\rangle} = \frac{\langle n^2\rangle - \langle n\rangle^2 - \langle n\rangle}{\langle n\rangle}$$

使用阶乘矩表示：
$$Q = \frac{\langle n^{(2)}\rangle}{\langle n\rangle} - 1$$

物理分类：
- $Q = 0$：泊松统计（相干光）
- $Q > 0$：超泊松统计（聚束光，如热光$Q = \bar{n}$）
- $Q < 0$：亚泊松统计（反聚束光，纯量子效应）
- $Q = -1$：Fock态（最小值）

#### Fano因子

定义：
$$F = \frac{(\Delta n)^2}{\langle n\rangle} = Q + 1$$

物理意义：
- $F = 1$：泊松噪声水平
- $F > 1$：超泊松（噪声增强）
- $F < 1$：亚泊松（噪声抑制）

#### 相对涨落

$$\frac{\Delta n}{\langle n\rangle} = \frac{1}{\sqrt{\langle n\rangle}}\sqrt{F}$$

对于大光子数，相干光的相对涨落按$1/\sqrt{\langle n\rangle}$减小。

### 27.3.4 光子聚束与反聚束

#### 物理图像

**聚束（Bunching）**：
- 光子倾向于成群到达探测器
- 经典波动性的体现
- 源于光场振幅的经典涨落

**反聚束（Antibunching）**：
- 光子倾向于分开到达
- 纯量子效应，无经典对应
- 体现光的粒子性

#### 条件概率解释

设在时刻$t$探测到一个光子，则在$t+\tau$探测到另一个光子的条件概率：

对于聚束光：
$$P(t+\tau|t) > P_{random}$$（高于随机情况）

对于反聚束光：
$$P(t+\tau|t) < P_{random}$$（低于随机情况）

#### 与二阶相干函数的关系

聚束/反聚束可通过$g^{(2)}(\tau)$定量描述：
- 聚束：$g^{(2)}(0) > 1$
- 随机（相干）：$g^{(2)}(0) = 1$
- 反聚束：$g^{(2)}(0) < 1$

### 27.3.5 光子统计的测量

#### 直接光子计数

使用单光子探测器（如APD、PMT）直接计数：
1. 设定计数时间窗口$T$
2. 重复测量获得分布$P(n)$
3. 计算统计参数

#### 探测器影响

实际探测器的影响：
- **量子效率**$\eta < 1$：测得的分布是真实分布与二项分布的卷积
- **暗计数**：增加泊松背景
- **死时间**：限制最大计数率

修正后的分布：
$$P_{measured}(m) = \sum_{n=m}^{\infty} P_{true}(n)\binom{n}{m}\eta^m(1-\eta)^{n-m}$$

#### 间接方法

1. **强度关联测量**：通过HBT实验测量$g^{(2)}(0)$推断Q参数
2. **平衡零拍探测**：测量场的正交分量涨落
3. **光子数分辨探测**：使用超导纳米线等技术直接分辨光子数

---

## 27.4 二阶相干函数g^(2)

### 27.4.1 定义

归一化二阶相干函数：
$$g^{(2)}(\tau) = \frac{\langle\hat{a}^\dagger(t)\hat{a}^\dagger(t+\tau)\hat{a}(t+\tau)\hat{a}(t)\rangle}{\langle\hat{a}^\dagger(t)\hat{a}(t)\rangle^2}$$

对于平稳过程：
$$g^{(2)}(\tau) = \frac{\langle\hat{a}^\dagger\hat{a}^\dagger(\tau)\hat{a}(\tau)\hat{a}\rangle}{\langle\hat{n}\rangle^2}$$

### 27.4.2 零延迟值g^(2)(0)

$$g^{(2)}(0) = \frac{\langle\hat{n}(\hat{n}-1)\rangle}{\langle\hat{n}\rangle^2} = \frac{\langle\hat{n}^2\rangle - \langle\hat{n}\rangle}{\langle\hat{n}\rangle^2}$$

与Mandel Q参数的关系：
$$g^{(2)}(0) = 1 + \frac{Q}{\langle n\rangle}$$

### 27.4.3 不同光源的g^(2)(0)

1. **相干光**：
   $$g^{(2)}(0) = 1$$

2. **热光**（混沌光）：
   $$g^{(2)}(0) = 2$$

3. **Fock态**$|n\rangle$：
   $$g^{(2)}(0) = \frac{n(n-1)}{n^2} = 1 - \frac{1}{n}$$
   
   特别地，单光子态：$g^{(2)}(0) = 0$

### 27.4.4 Hanbury Brown-Twiss实验

实验装置测量强度关联：
$$G^{(2)}(\tau) = \langle I(t)I(t+\tau)\rangle$$

归一化：
$$g^{(2)}(\tau) = \frac{G^{(2)}(\tau)}{\langle I\rangle^2}$$

### 27.4.5 经典与量子界限

经典光场的Cauchy-Schwarz不等式：
$$g^{(2)}(0) \geq 1$$

量子光可以违反此界限：
$$0 \leq g^{(2)}(0) < 1$$ （反聚束，纯量子效应）

---

## 27.5 量子噪声与散粒噪声

### 27.5.1 散粒噪声的起源

光电探测中的电流：
$$I(t) = e\sum_{i}\delta(t-t_i)$$

其中$t_i$是光电子到达时间。

平均电流：
$$\langle I\rangle = e\langle\dot{N}\rangle = e\eta P/\hbar\omega$$

其中$\eta$是量子效率，$P$是光功率。

### 27.5.2 噪声功率谱

散粒噪声的功率谱密度（白噪声）：
$$S_I(f) = 2e\langle I\rangle$$

对于相干光，光子数涨落导致的电流噪声：
$$\langle\Delta I^2\rangle = e^2\langle\Delta n^2\rangle/T^2 = e\langle I\rangle/T$$

### 27.5.3 信噪比

光电探测的信噪比：
$$\text{SNR} = \frac{\langle I\rangle}{\sqrt{\langle\Delta I^2\rangle}} = \sqrt{\frac{\eta P T}{\hbar\omega}}$$

这定义了散粒噪声极限。

### 27.5.4 标准量子极限

相位测量的不确定性：
$$\Delta\phi \geq \frac{1}{2\sqrt{\langle n\rangle}}$$

这是使用相干光的标准量子极限（SQL）。

### 27.5.5 压缩态降噪

使用压缩光可以突破标准量子极限：

1. **振幅压缩**：降低强度噪声
   $$(\Delta I)_{squeezed} = e^{-r}(\Delta I)_{coherent}$$

2. **相位压缩**：提高相位测量精度
   $$(\Delta\phi)_{squeezed} = e^{-r}(\Delta\phi)_{coherent}$$

### 27.5.6 量子噪声在成像中的影响

成像系统的量子噪声限制：
- 每个像素的光子数涨落：$\Delta n_{pixel} = \sqrt{n_{pixel}}$
- 图像信噪比：$\text{SNR} = \sqrt{n_{pixel}}$
- 需要的总光子数：$N_{total} = N_{pixels} \times \text{SNR}^2$

---

## 本章小结

本章介绍了量子光学的核心概念：

1. **光的量子化**：从经典电磁场过渡到量子场论描述，引入产生湮灭算符和Fock态
2. **相干态**：最接近经典光的量子态，具有泊松光子统计和最小不确定性
3. **压缩态**：通过重新分配量子涨落，可在某个正交分量上突破标准量子极限
4. **光子统计**：Mandel Q参数区分经典（Q≥0）和非经典光（Q<0）
5. **二阶相干函数**：g^(2)(0)<1标志着量子反聚束效应
6. **量子噪声**：散粒噪声源于光的粒子性，定义了测量的基本极限

这些概念为理解量子成像、量子计算和未来光学技术奠定了基础。

---

## 练习题

### 基础题

**27.1** 证明相干态不正交：计算$\langle\alpha|\beta\rangle$并说明其物理意义。

*提示*：使用相干态的Fock基展开式。

<details>
<summary>答案</summary>

$$\langle\alpha|\beta\rangle = e^{-\frac{1}{2}(|\alpha|^2+|\beta|^2-2\alpha^*\beta)}$$

当$|\alpha-\beta|^2 \gg 1$时，两态近似正交。这反映了相干态的准经典特性。
</details>

**27.2** 对于热光场，证明$g^{(2)}(0) = 2$。假设热光服从玻色-爱因斯坦分布。

*提示*：计算$\langle n^2\rangle$和$\langle n\rangle$的关系。

<details>
<summary>答案</summary>

对于热光：$P(n) = \frac{\bar{n}^n}{(1+\bar{n})^{n+1}}$

计算得：$\langle n^2\rangle = 2\bar{n}^2 + \bar{n}$

因此：$g^{(2)}(0) = \frac{\langle n^2\rangle - \langle n\rangle}{\langle n\rangle^2} = 2$
</details>

**27.3** 计算压缩真空态的光子数分布$P(n)$。

*提示*：只有偶数光子数态有非零概率。

<details>
<summary>答案</summary>

$$P(2m) = \frac{(2m)!}{2^{2m}(m!)^2}\frac{(\tanh r)^{2m}}{\cosh r}$$
$$P(2m+1) = 0$$

平均光子数：$\langle n\rangle = \sinh^2 r$
</details>

### 挑战题

**27.4** 推导压缩相干态$|\alpha,\xi\rangle = \hat{D}(\alpha)\hat{S}(\xi)|0\rangle$的g^(2)(0)。

*提示*：先计算$\langle\hat{n}\rangle$和$\langle\hat{n}^2\rangle$。

<details>
<summary>答案</summary>

经过复杂计算：
$$g^{(2)}(0) = 1 + \frac{2\sinh^2 r}{(|\alpha|^2\cosh 2r + \sinh^2 r)^2}[\cosh 2r - \text{Re}(\alpha^2e^{-i\theta}/|\alpha|^2)]$$

其中$\xi = re^{i\theta}$。当$r=0$（无压缩）时，回到$g^{(2)}(0)=1$。
</details>

**27.5** 考虑双模压缩真空态（参量下转换产生）。证明两个模式间存在完美关联。

*提示*：计算联合光子数分布$P(n_1,n_2)$。

<details>
<summary>答案</summary>

双模压缩态：$|\psi\rangle = \frac{1}{\cosh r}\sum_{n=0}^{\infty}(\tanh r)^n|n\rangle_1|n\rangle_2$

联合分布：$P(n_1,n_2) = \delta_{n_1,n_2}\frac{(\tanh r)^{2n_1}}{\cosh^2 r}$

完美关联：测量模式1得到n个光子，模式2必定也是n个光子。
</details>

**27.6** 设计一个实验方案，区分单光子源和衰减的相干光源。两者平均光子数都很小（$\langle n\rangle \ll 1$）。

*提示*：利用g^(2)(0)的差异。

<details>
<summary>答案</summary>

使用Hanbury Brown-Twiss装置测量g^(2)(0)：

1. 单光子源：g^(2)(0) = 0（理想情况）
2. 衰减相干光：g^(2)(0) = 1

即使$\langle n\rangle$相同，通过测量符合计数率可以明确区分。实际单光子源的$g^{(2)}(0) < 0.5$即可认为是量子光源。
</details>

### 开放性思考题

**27.7** 在量子密钥分发(QKD)中，为什么单光子源比衰减激光更安全？从光子统计角度分析。

**27.8** 讨论如何将量子光学概念应用于计算机图形学的全局照明算法。考虑光子映射中的"光子"与量子光学中光子的本质区别。

---

## 常见陷阱与错误

1. **混淆经典相干性与量子相干性**
   - 错误：认为激光是"量子光"
   - 正确：激光是相干态，最接近经典光，g^(2)(0)=1

2. **误解光子数态**
   - 错误：认为"n个光子"的Fock态容易制备
   - 正确：Fock态极难制备，需要特殊的非线性过程

3. **压缩态的误用**
   - 错误：认为压缩可以同时减小所有噪声
   - 正确：压缩只是重新分配噪声，总不确定性不变

4. **g^(2)函数的测量**
   - 错误：直接测量光强关联
   - 正确：需要考虑探测器响应时间、死时间等因素

5. **量子效率的忽视**
   - 错误：假设探测器完美
   - 正确：实际探测器η<1会改变测量的光子统计

---

## 最佳实践检查清单

### 理论分析
- [ ] 明确区分经典场与量子场描述
- [ ] 正确使用产生湮灭算符的对易关系
- [ ] 检查态的归一化条件
- [ ] 验证物理量的期望值合理性

### 实验设计
- [ ] 考虑探测器的量子效率和暗计数
- [ ] 评估需要的积分时间达到统计显著性
- [ ] 选择合适的光源功率避免饱和
- [ ] 考虑背景光和杂散光的影响

### 数值计算
- [ ] 使用适当的Hilbert空间截断
- [ ] 检查数值精度对高阶矩的影响
- [ ] 验证概率分布的归一化
- [ ] 考虑数值不稳定性（如大光子数）

### 应用考虑
- [ ] 评估量子优势的实际可达性
- [ ] 考虑退相干和损耗的影响
- [ ] 比较量子方案与经典方案的资源需求
- [ ] 分析系统的可扩展性
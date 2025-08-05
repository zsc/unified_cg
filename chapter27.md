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

经典电磁场的能量密度包含电场和磁场贡献：
$$u_{em} = \frac{1}{2}\left[\epsilon_0 E^2(r,t) + \frac{1}{\mu_0}B^2(r,t)\right]$$

总能量通过空间积分得到：
$$H_{classical} = \int d^3r \, u_{em} = \frac{1}{2}\int d^3r \left[\epsilon_0 E^2(r,t) + \frac{1}{\mu_0}B^2(r,t)\right]$$

在有限体积$V$的谐振腔中，满足边界条件的电磁场可以展开为正交模式的叠加。对于立方腔（边长$L$），模式由波矢$\mathbf{k} = \frac{2\pi}{L}(n_x, n_y, n_z)$标记，其中$n_i$为整数。

矢势展开式：
$$\mathbf{A}(r,t) = \sum_{k,s} \sqrt{\frac{\hbar}{2\epsilon_0\omega_k V}} \epsilon_{k,s} \left[a_{k,s}(t)e^{i\mathbf{k} \cdot \mathbf{r}} + a_{k,s}^*(t)e^{-i\mathbf{k} \cdot \mathbf{r}}\right]$$

其中$s=1,2$标记两个正交偏振方向，$\epsilon_{k,s}$是单位偏振矢量，满足$\epsilon_{k,s} \cdot \mathbf{k} = 0$（横波条件）和$\epsilon_{k,1} \cdot \epsilon_{k,2} = 0$（正交性）。

电场和磁场由矢势导出：
$$\mathbf{E} = -\frac{\partial \mathbf{A}}{\partial t}, \quad \mathbf{B} = \nabla \times \mathbf{A}$$

代入矢势展开式：
$$\mathbf{E}(r,t) = i\sum_{k,s} \sqrt{\frac{\hbar\omega_k}{2\epsilon_0 V}} \epsilon_{k,s} \left[a_{k,s}(t)e^{i\mathbf{k} \cdot \mathbf{r}} - a_{k,s}^*(t)e^{-i\mathbf{k} \cdot \mathbf{r}}\right]$$

### 27.1.2 量子化过程

正则量子化要求识别共轭变量。对于电磁场，我们定义：
- 广义坐标：$q_{k,s} = \frac{1}{\sqrt{2}}(a_{k,s} + a_{k,s}^*)$
- 广义动量：$p_{k,s} = \frac{i}{\sqrt{2\omega_k}}(a_{k,s}^* - a_{k,s})$

满足泊松括号：$\{q_{k,s}, p_{k',s'}\} = \delta_{k,k'}\delta_{s,s'}$

量子化时，将泊松括号替换为对易子：
$$[q_{k,s}, p_{k',s'}] = i\hbar\delta_{k,k'}\delta_{s,s'}$$

引入算符：
$$\hat{a}_{k,s} = \sqrt{\frac{\omega_k}{2\hbar}}\hat{q}_{k,s} + \frac{i}{\sqrt{2\hbar\omega_k}}\hat{p}_{k,s}$$
$$\hat{a}_{k,s}^\dagger = \sqrt{\frac{\omega_k}{2\hbar}}\hat{q}_{k,s} - \frac{i}{\sqrt{2\hbar\omega_k}}\hat{p}_{k,s}$$

这些算符满足玻色对易关系：
$$[\hat{a}_{k,s}, \hat{a}_{k',s'}^\dagger] = \delta_{k,k'}\delta_{s,s'}$$
$$[\hat{a}_{k,s}, \hat{a}_{k',s'}] = [\hat{a}_{k,s}^\dagger, \hat{a}_{k',s'}^\dagger] = 0$$

单模哈密顿量成为：
$$\hat{H}_{k,s} = \hbar\omega_k\left(\hat{a}_{k,s}^\dagger\hat{a}_{k,s} + \frac{1}{2}\right)$$

总哈密顿量：
$$\hat{H} = \sum_{k,s} \hbar\omega_k\left(\hat{a}_{k,s}^\dagger\hat{a}_{k,s} + \frac{1}{2}\right)$$

### 27.1.3 Fock态

定义真空态$|0\rangle$为所有模式的基态：
$$\hat{a}_{k,s}|0\rangle = 0, \quad \forall k,s$$

Fock态（光子数态）通过反复作用产生算符构造：
$$|n_{k,s}\rangle = \frac{(\hat{a}_{k,s}^\dagger)^n}{\sqrt{n!}}|0\rangle$$

多模Fock态：
$$|\{n_{k,s}\}\rangle = \prod_{k,s} \frac{(\hat{a}_{k,s}^\dagger)^{n_{k,s}}}{\sqrt{n_{k,s}!}}|0\rangle$$

算符作用规则：
- $\hat{a}_{k,s}^\dagger|n_{k,s}\rangle = \sqrt{n_{k,s}+1}|n_{k,s}+1\rangle$ （产生一个光子）
- $\hat{a}_{k,s}|n_{k,s}\rangle = \sqrt{n_{k,s}}|n_{k,s}-1\rangle$ （湮灭一个光子）
- $\hat{n}_{k,s}|n_{k,s}\rangle = n_{k,s}|n_{k,s}\rangle$ （光子数本征态）

能量本征值：
$$\hat{H}|\{n_{k,s}\}\rangle = \left(\sum_{k,s} \hbar\omega_k\left(n_{k,s} + \frac{1}{2}\right)\right)|\{n_{k,s}\}\rangle$$

### 27.1.4 真空涨落

真空态虽然没有光子，但具有非零能量：
$$E_0 = \langle 0|\hat{H}|0\rangle = \sum_{k,s} \frac{\hbar\omega_k}{2}$$

这个无穷大可通过正规序（normal ordering）处理，定义$:\hat{H}: = \sum_{k,s} \hbar\omega_k\hat{a}_{k,s}^\dagger\hat{a}_{k,s}$。

更重要的是场算符的真空期望值。虽然：
$$\langle 0|\hat{\mathbf{E}}(r,t)|0\rangle = 0$$

但均方涨落非零：
$$\langle 0|\hat{\mathbf{E}}^2(r,t)|0\rangle = \sum_{k,s} \frac{\hbar\omega_k}{2\epsilon_0 V}$$

单模电场涨落：
$$(\Delta E)_{vacuum} = \sqrt{\frac{\hbar\omega}{2\epsilon_0 V}}$$

这个真空涨落具有可观测效应：
1. **Casimir效应**：两平行导体板间的真空能依赖于板间距
2. **Lamb位移**：原子能级因与真空场耦合而移动
3. **自发辐射**：激发态原子在真空涨落驱动下跃迁

### 27.1.5 场的正交分量

定义单模场的正交振幅算符：
$$\hat{X} = \frac{1}{\sqrt{2}}(\hat{a} + \hat{a}^\dagger), \quad \hat{P} = \frac{i}{\sqrt{2}}(\hat{a}^\dagger - \hat{a})$$

这些算符满足：
$$[\hat{X}, \hat{P}] = i$$

对应的不确定性关系：
$$\Delta X \Delta P \geq \frac{1}{2}$$

在真空态中：
$$\langle 0|\hat{X}^2|0\rangle = \langle 0|\hat{P}^2|0\rangle = \frac{1}{2}$$

因此：$\Delta X = \Delta P = \frac{1}{\sqrt{2}}$，达到最小不确定性。

---

## 27.2 相干态与压缩态

### 27.2.1 相干态定义

相干态代表了量子光学中最接近经典光的量子态。它们是湮灭算符的本征态：
$$\hat{a}|\alpha\rangle = \alpha|\alpha\rangle$$

其中$\alpha = |\alpha|e^{i\phi}$是复数本征值，$|\alpha|$对应经典振幅，$\phi$对应相位。

#### Fock基展开

在光子数基下，相干态表示为：
$$|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^{\infty}\frac{\alpha^n}{\sqrt{n!}}|n\rangle$$

这个展开可以通过求解本征方程得到。设$|\alpha\rangle = \sum_n c_n|n\rangle$，则：
$$\hat{a}|\alpha\rangle = \sum_n c_n\sqrt{n}|n-1\rangle = \alpha\sum_n c_n|n\rangle$$

比较系数得递推关系：$c_{n+1} = \frac{\alpha}{\sqrt{n+1}}c_n$

归一化条件$\langle\alpha|\alpha\rangle = 1$确定$c_0 = e^{-|\alpha|^2/2}$。

#### 相干态的非正交性

两个相干态的内积：
$$\langle\alpha|\beta\rangle = e^{-\frac{1}{2}(|\alpha|^2+|\beta|^2-2\alpha^*\beta)} = e^{-\frac{1}{2}|\alpha-\beta|^2}e^{i\text{Im}(\alpha^*\beta)}$$

这表明相干态不正交，但当$|\alpha-\beta| \gg 1$时近似正交。

#### 过完备性

相干态构成过完备基：
$$\frac{1}{\pi}\int d^2\alpha |\alpha\rangle\langle\alpha| = \mathbb{I}$$

其中$d^2\alpha = d(\text{Re}\alpha)d(\text{Im}\alpha)$。这允许任意态展开为相干态的叠加。

### 27.2.2 相干态性质

#### 光子统计

1. **平均光子数**：
   $$\langle\hat{n}\rangle = \langle\alpha|\hat{a}^\dagger\hat{a}|\alpha\rangle = |\alpha|^2$$

2. **光子数方差**：
   $$\langle\hat{n}^2\rangle = \langle\alpha|\hat{a}^\dagger\hat{a}\hat{a}^\dagger\hat{a}|\alpha\rangle = |\alpha|^4 + |\alpha|^2$$
   
   因此：$(\Delta n)^2 = \langle\hat{n}^2\rangle - \langle\hat{n}\rangle^2 = |\alpha|^2 = \langle\hat{n}\rangle$

3. **光子数分布**（泊松分布）：
   $$P(n) = |\langle n|\alpha\rangle|^2 = \frac{|\alpha|^{2n}}{n!}e^{-|\alpha|^2} = \frac{\langle n\rangle^n}{n!}e^{-\langle n\rangle}$$

#### 场的期望值

电场算符的期望值：
$$\langle\alpha|\hat{E}(r,t)|\alpha\rangle = 2\sqrt{\frac{\hbar\omega}{2\epsilon_0 V}}|\alpha|\cos(k \cdot r - \omega t + \phi)$$

这正是经典相干波的形式，振幅正比于$|\alpha|$。

#### 相位空间表示

在$(X,P)$相位空间中，相干态表现为：
- 中心位置：$\langle X\rangle = \sqrt{2}\text{Re}(\alpha)$，$\langle P\rangle = \sqrt{2}\text{Im}(\alpha)$
- 不确定性：$\Delta X = \Delta P = \frac{1}{\sqrt{2}}$（各向同性高斯分布）

### 27.2.3 位移算符

位移算符定义：
$$\hat{D}(\alpha) = \exp(\alpha\hat{a}^\dagger - \alpha^*\hat{a})$$

#### Baker-Campbell-Hausdorff公式

利用$[\hat{a}, \hat{a}^\dagger] = 1$，可得：
$$\hat{D}(\alpha) = e^{-|\alpha|^2/2}e^{\alpha\hat{a}^\dagger}e^{-\alpha^*\hat{a}} = e^{|\alpha|^2/2}e^{-\alpha^*\hat{a}}e^{\alpha\hat{a}^\dagger}$$

#### 位移算符性质

1. **幺正性**：$\hat{D}^\dagger(\alpha) = \hat{D}^{-1}(\alpha) = \hat{D}(-\alpha)$

2. **乘法规则**：
   $$\hat{D}(\alpha)\hat{D}(\beta) = e^{i\text{Im}(\alpha^*\beta)}\hat{D}(\alpha+\beta)$$

3. **变换性质**：
   $$\hat{D}^\dagger(\alpha)\hat{a}\hat{D}(\alpha) = \hat{a} + \alpha$$
   $$\hat{D}^\dagger(\alpha)\hat{a}^\dagger\hat{D}(\alpha) = \hat{a}^\dagger + \alpha^*$$

4. **生成相干态**：
   $$|\alpha\rangle = \hat{D}(\alpha)|0\rangle$$

### 27.2.4 压缩态

压缩算符定义：
$$\hat{S}(\xi) = \exp\left[\frac{1}{2}(\xi^*\hat{a}^2 - \xi\hat{a}^{\dagger 2})\right]$$

其中$\xi = re^{i\theta}$，$r \geq 0$是压缩强度，$\theta$是压缩方向。

#### 压缩算符的分解

利用$SU(1,1)$李代数，可将压缩算符分解为：
$$\hat{S}(\xi) = \exp\left[\frac{\tanh r}{2}e^{-i\theta}\hat{a}^{\dagger 2}\right]\exp\left[-\ln(\cosh r)\hat{a}^\dagger\hat{a}\right]\exp\left[-\frac{\tanh r}{2}e^{i\theta}\hat{a}^2\right]$$

#### 压缩真空态

压缩真空态的Fock基展开：
$$|\xi\rangle = \hat{S}(\xi)|0\rangle = \frac{1}{\sqrt{\cosh r}}\sum_{n=0}^{\infty}\frac{\sqrt{(2n)!}}{2^n n!}(-e^{i\theta}\tanh r)^n|2n\rangle$$

注意只有偶数光子态出现，这是双光子过程的体现。

#### 压缩相干态

一般的压缩相干态：
$$|\alpha,\xi\rangle = \hat{D}(\alpha)\hat{S}(\xi)|0\rangle$$

顺序很重要，因为$[\hat{D}(\alpha), \hat{S}(\xi)] \neq 0$。

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
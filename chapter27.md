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

经典电磁场的哈密顿量为：
$$H_{classical} = \frac{1}{2}\int d^3r \left[\epsilon_0 E^2(r,t) + \frac{1}{\mu_0}B^2(r,t)\right]$$

在腔内，电磁场可展开为模式：
$$E(r,t) = \sum_k \sqrt{\frac{\hbar\omega_k}{2\epsilon_0 V}} \epsilon_k \left[a_k(t)e^{ik \cdot r} + a_k^*(t)e^{-ik \cdot r}\right]$$

其中$V$是量子化体积，$\epsilon_k$是偏振矢量。

### 27.1.2 量子化过程

将经典场振幅$a_k$和$a_k^*$提升为算符$\hat{a}_k$和$\hat{a}_k^\dagger$，满足对易关系：
$$[\hat{a}_k, \hat{a}_{k'}^\dagger] = \delta_{k,k'}$$
$$[\hat{a}_k, \hat{a}_{k'}] = [\hat{a}_k^\dagger, \hat{a}_{k'}^\dagger] = 0$$

单模哈密顿量变为：
$$\hat{H}_k = \hbar\omega_k\left(\hat{a}_k^\dagger\hat{a}_k + \frac{1}{2}\right)$$

### 27.1.3 Fock态

光子数态（Fock态）定义为：
$$|n\rangle = \frac{(\hat{a}^\dagger)^n}{\sqrt{n!}}|0\rangle$$

满足：
- $\hat{a}^\dagger|n\rangle = \sqrt{n+1}|n+1\rangle$ （产生算符）
- $\hat{a}|n\rangle = \sqrt{n}|n-1\rangle$ （湮灭算符）
- $\hat{n}|n\rangle = \hat{a}^\dagger\hat{a}|n\rangle = n|n\rangle$ （数算符）

### 27.1.4 真空涨落

真空态$|0\rangle$的能量：
$$\langle 0|\hat{H}|0\rangle = \frac{1}{2}\hbar\omega$$

电场的真空涨落：
$$\langle 0|\hat{E}^2|0\rangle = \frac{\hbar\omega}{2\epsilon_0 V}$$

这导致了零点能和真空涨落现象。

---

## 27.2 相干态与压缩态

### 27.2.1 相干态定义

相干态$|\alpha\rangle$是湮灭算符的本征态：
$$\hat{a}|\alpha\rangle = \alpha|\alpha\rangle$$

其中$\alpha = |\alpha|e^{i\phi}$是复数。

在Fock基下展开：
$$|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^{\infty}\frac{\alpha^n}{\sqrt{n!}}|n\rangle$$

### 27.2.2 相干态性质

1. **平均光子数**：
   $$\langle\hat{n}\rangle = \langle\alpha|\hat{a}^\dagger\hat{a}|\alpha\rangle = |\alpha|^2$$

2. **光子数分布**（泊松分布）：
   $$P(n) = |\langle n|\alpha\rangle|^2 = \frac{|\alpha|^{2n}}{n!}e^{-|\alpha|^2}$$

3. **最小不确定性**：
   $$\Delta X \Delta P = \frac{\hbar}{2}$$

### 27.2.3 位移算符

相干态可由真空态通过位移算符生成：
$$|\alpha\rangle = \hat{D}(\alpha)|0\rangle$$

其中：
$$\hat{D}(\alpha) = \exp(\alpha\hat{a}^\dagger - \alpha^*\hat{a})$$

### 27.2.4 压缩态

压缩算符定义为：
$$\hat{S}(\xi) = \exp\left[\frac{1}{2}(\xi^*\hat{a}^2 - \xi\hat{a}^{\dagger 2})\right]$$

其中$\xi = re^{i\theta}$是压缩参数。

压缩真空态：
$$|\xi\rangle = \hat{S}(\xi)|0\rangle$$

### 27.2.5 压缩态的不确定性

定义正交振幅：
$$\hat{X} = \frac{1}{\sqrt{2}}(\hat{a} + \hat{a}^\dagger), \quad \hat{P} = \frac{i}{\sqrt{2}}(\hat{a}^\dagger - \hat{a})$$

对于压缩态：
$$\Delta X = \frac{e^{-r}}{2}, \quad \Delta P = \frac{e^r}{2}$$

仍满足最小不确定性：$\Delta X \Delta P = \frac{1}{4}$（设$\hbar=1$）。

---

## 27.3 光子统计

### 27.3.1 光子计数分布

一般光态$\hat{\rho}$的光子数分布：
$$P(n) = \langle n|\hat{\rho}|n\rangle$$

一阶矩和二阶矩：
$$\langle n\rangle = \text{Tr}(\hat{\rho}\hat{n})$$
$$\langle n^2\rangle = \text{Tr}(\hat{\rho}\hat{n}^2)$$

### 27.3.2 光子数方差

方差定义：
$$(\Delta n)^2 = \langle n^2\rangle - \langle n\rangle^2$$

对于不同光源：
- **相干光**（泊松统计）：$(\Delta n)^2 = \langle n\rangle$
- **热光**（玻色-爱因斯坦统计）：$(\Delta n)^2 = \langle n\rangle + \langle n\rangle^2$
- **亚泊松光**：$(\Delta n)^2 < \langle n\rangle$

### 27.3.3 Mandel Q参数

定义：
$$Q = \frac{(\Delta n)^2 - \langle n\rangle}{\langle n\rangle} = \frac{\langle n^2\rangle - \langle n\rangle^2 - \langle n\rangle}{\langle n\rangle}$$

物理意义：
- $Q = 0$：泊松统计（相干光）
- $Q > 0$：超泊松统计（聚束光）
- $Q < 0$：亚泊松统计（反聚束光，非经典）

### 27.3.4 Fano因子

另一个常用参数：
$$F = \frac{(\Delta n)^2}{\langle n\rangle}$$

关系：$F = Q + 1$

### 27.3.5 光子聚束与反聚束

**聚束**（bunching）：光子倾向于成群到达
- 热光表现出聚束效应
- $g^{(2)}(0) > 1$

**反聚束**（antibunching）：光子倾向于分开到达
- 单光子源表现出反聚束
- $g^{(2)}(0) < 1$（量子效应）

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
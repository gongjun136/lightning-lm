# BCH 公式

Baker-Campbell-Hausdorff 公式简称 BCH 公式。它回答的问题是：**两个李群元素相乘以后，结果在李代数中应该如何表示**。

对普通向量，有

$$
\mathbf{a}+\mathbf{b}=\mathbf{b}+\mathbf{a}
$$

但旋转矩阵乘法一般不可交换，因此通常没有

$$
\operatorname{Exp}(\boldsymbol{\phi}_1)\operatorname{Exp}(\boldsymbol{\phi}_2)
=
\operatorname{Exp}(\boldsymbol{\phi}_1+\boldsymbol{\phi}_2)
$$

BCH 公式正是用来描述

$$
\operatorname{Log}\!\left(
\operatorname{Exp}(\boldsymbol{\phi}_1)
\operatorname{Exp}(\boldsymbol{\phi}_2)
\right)
$$

和 $\boldsymbol{\phi}_1,\boldsymbol{\phi}_2$ 之间关系的工具。ESKF、优化和预积分推导中经常要把李群上的乘法扰动转换成李代数上的加法扰动，因此 BCH 的一阶近似非常常用。

## 基本符号

本文主要讨论 $\operatorname{SO}(3)$。旋转向量写作 $\boldsymbol{\phi}\in\mathbb{R}^3$，对应的李代数矩阵为 $[\boldsymbol{\phi}]_\times\in\mathfrak{so}(3)$：

$$
\operatorname{Exp}(\boldsymbol{\phi})
=
\exp\!\left([\boldsymbol{\phi}]_\times\right),
\qquad
\operatorname{Log}(\mathbf{R})
=
\log(\mathbf{R})^\vee
$$

李括号定义为：

$$
[\mathbf{A},\mathbf{B}]
=
\mathbf{A}\mathbf{B}-\mathbf{B}\mathbf{A}
$$

在 $\mathfrak{so}(3)$ 中，若 $\mathbf{A}=[\boldsymbol{\phi}_1]_\times$、$\mathbf{B}=[\boldsymbol{\phi}_2]_\times$，则

$$
[\mathbf{A},\mathbf{B}]
=
\left[\boldsymbol{\phi}_1\times\boldsymbol{\phi}_2\right]_\times
$$

这个式子右侧的 $\left[\boldsymbol{\phi}_1\times\boldsymbol{\phi}_2\right]_\times$ 表示：先计算向量叉乘 $\boldsymbol{\phi}_1\times\boldsymbol{\phi}_2$，再把得到的向量转换成反对称矩阵。也就是说，**矩阵李括号的结果仍然是一个 $\mathfrak{so}(3)$ 中的反对称矩阵**。在只写旋转向量时，李括号也常简写成：

$$
[\boldsymbol{\phi}_1,\boldsymbol{\phi}_2]
=
\boldsymbol{\phi}_1\times\boldsymbol{\phi}_2
$$

## BCH 完整形式

对两个李代数元素 $\mathbf{A},\mathbf{B}$，BCH 公式为：

$$
\begin{aligned}
\log\!\left(\exp(\mathbf{A})\exp(\mathbf{B})\right)
&=
\mathbf{A}+\mathbf{B}
+\frac{1}{2}[\mathbf{A},\mathbf{B}] \\
&\quad
+\frac{1}{12}[\mathbf{A},[\mathbf{A},\mathbf{B}]]
+\frac{1}{12}[\mathbf{B},[\mathbf{B},\mathbf{A}]]
+\cdots
\end{aligned}
$$

等价地，在 $\operatorname{SO}(3)$ 的旋转向量记法下：

$$
\begin{aligned}
&\operatorname{Log}\!\left(
\operatorname{Exp}(\boldsymbol{\phi}_1)
\operatorname{Exp}(\boldsymbol{\phi}_2)
\right) \\
&=
\boldsymbol{\phi}_1+\boldsymbol{\phi}_2
+\frac{1}{2}(\boldsymbol{\phi}_1\times\boldsymbol{\phi}_2)
+\frac{1}{12}\boldsymbol{\phi}_1\times
(\boldsymbol{\phi}_1\times\boldsymbol{\phi}_2)
+\frac{1}{12}\boldsymbol{\phi}_2\times
(\boldsymbol{\phi}_2\times\boldsymbol{\phi}_1)
+\cdots
\end{aligned}
$$

如果 $\boldsymbol{\phi}_1,\boldsymbol{\phi}_2$ 都是小量，保留到二阶项可得：

$$
\operatorname{Log}\!\left(
\operatorname{Exp}(\boldsymbol{\phi}_1)
\operatorname{Exp}(\boldsymbol{\phi}_2)
\right)
\approx
\boldsymbol{\phi}_1+\boldsymbol{\phi}_2
+\frac{1}{2}\boldsymbol{\phi}_1\times\boldsymbol{\phi}_2
$$

若只保留一阶项，则进一步退化为：

$$
\operatorname{Log}\!\left(
\operatorname{Exp}(\boldsymbol{\phi}_1)
\operatorname{Exp}(\boldsymbol{\phi}_2)
\right)
\approx
\boldsymbol{\phi}_1+\boldsymbol{\phi}_2
$$

**一阶近似只适用于两个旋转向量都足够小的情况；如果其中一个旋转不是小量，就需要引入左雅可比或右雅可比。**

## 左右雅可比

设

$$
\boldsymbol{\phi}=\theta\mathbf{n},
\qquad
\theta=\|\boldsymbol{\phi}\|,
\qquad
\|\mathbf{n}\|=1
$$

$\operatorname{SO}(3)$ 的左雅可比为：

$$
\begin{aligned}
\mathbf{J}_l(\boldsymbol{\phi})
&=
\frac{\sin\theta}{\theta}\mathbf{I}
+\left(1-\frac{\sin\theta}{\theta}\right)\mathbf{n}\mathbf{n}^\mathsf{T}
+\frac{1-\cos\theta}{\theta}[\mathbf{n}]_\times \\
&=
\mathbf{I}
+\frac{1-\cos\theta}{\theta^2}[\boldsymbol{\phi}]_\times
+\frac{\theta-\sin\theta}{\theta^3}[\boldsymbol{\phi}]_\times^2
\end{aligned}
$$

左雅可比逆为：

$$
\mathbf{J}_l^{-1}(\boldsymbol{\phi})
=
\mathbf{I}
-\frac{1}{2}[\boldsymbol{\phi}]_\times
+\left(
\frac{1}{\theta^2}
-\frac{1+\cos\theta}{2\theta\sin\theta}
\right)
[\boldsymbol{\phi}]_\times^2
$$

> **左雅可比逆的验证：**
>
> 记
>
> $$
> \mathbf{A}=[\boldsymbol{\phi}]_\times,
> \qquad
> \theta=\|\boldsymbol{\phi}\|
> $$
>
> 叉乘矩阵满足
>
> $$
> \mathbf{A}^3=-\theta^2\mathbf{A},
> \qquad
> \mathbf{A}^4=-\theta^2\mathbf{A}^2
> $$
>
> 左雅可比可以写成
>
> $$
> \mathbf{J}_l(\boldsymbol{\phi})
> =
> \mathbf{I}+a\mathbf{A}+b\mathbf{A}^2
> $$
>
> 其中
>
> $$
> a=\frac{1-\cos\theta}{\theta^2},
> \qquad
> b=\frac{\theta-\sin\theta}{\theta^3}
> $$
>
> 假设它的逆也可以写成同样的二次多项式形式：
>
> $$
> \mathbf{J}_l^{-1}(\boldsymbol{\phi})
> =
> \mathbf{I}+c\mathbf{A}+d\mathbf{A}^2
> $$
>
> 两者相乘：
>
> $$
> \begin{aligned}
> \mathbf{J}_l\mathbf{J}_l^{-1}
> &=
> (\mathbf{I}+a\mathbf{A}+b\mathbf{A}^2)
> (\mathbf{I}+c\mathbf{A}+d\mathbf{A}^2) \\
> &=
> \mathbf{I}
> +\left(a+c-\theta^2(ad+bc)\right)\mathbf{A} \\
> &\quad
> +\left(b+d+ac-\theta^2bd\right)\mathbf{A}^2
> \end{aligned}
> $$
>
> 若它确实为逆矩阵，则 $\mathbf{A}$ 和 $\mathbf{A}^2$ 的系数都应为 0。代入
>
> $$
> c=-\frac{1}{2},
> \qquad
> d=
> \frac{1}{\theta^2}
> -\frac{1+\cos\theta}{2\theta\sin\theta}
> $$
>
> 其中 $\mathbf{A}$ 的系数可以化简为
>
> $$
> \begin{aligned}
> a+c-\theta^2(ad+bc)
> &=
> \left(a-\frac{1}{2}\right)
> -\theta^2
> \left[
> \frac{1-\cos\theta}{\theta^4}
> -\frac{1}{2\theta^2}
> \right] \\
> &=
> \left(a-\frac{1}{2}\right)
> -\left(a-\frac{1}{2}\right)
> =
> 0
> \end{aligned}
> $$
>
> $\mathbf{A}^2$ 的系数可以化简为
>
> $$
> \begin{aligned}
> b+d+ac-\theta^2bd
> &=
> \frac{\theta-\sin\theta}{\theta^3}
> +\frac{1+\cos\theta}{2\theta^2}
> -\frac{1+\cos\theta}{2\theta\sin\theta} \\
> &\quad
> -\frac{\theta-\sin\theta}{\theta^3}
> +\frac{(\theta-\sin\theta)(1+\cos\theta)}
> {2\theta^2\sin\theta} \\
> &=
> \frac{1+\cos\theta}{2}
> \left[
> \frac{1}{\theta^2}
> -\frac{1}{\theta\sin\theta}
> +\frac{\theta-\sin\theta}{\theta^2\sin\theta}
> \right] \\
> &=
> 0
> \end{aligned}
> $$
>
> 因此
>
> $$
> \mathbf{J}_l(\boldsymbol{\phi})
> \mathbf{J}_l^{-1}(\boldsymbol{\phi})
> =
> \mathbf{I}
> $$
>

右雅可比可以由左雅可比得到：

$$
\mathbf{J}_r(\boldsymbol{\phi})
=
\mathbf{J}_l(-\boldsymbol{\phi})=\mathbf{J}_l(\boldsymbol{\phi})^\mathsf{T}
$$

> 由叉乘矩阵的反对称性还可以得到一个常用关系。令
>
> $$
> \mathbf{A}=[\boldsymbol{\phi}]_\times,
> \qquad
> \mathbf{A}^\mathsf{T}=-\mathbf{A},
> \qquad
> (\mathbf{A}^2)^\mathsf{T}=\mathbf{A}^2
> $$
>
> 左雅可比的闭式为 $\mathbf{J}_l(\boldsymbol{\phi})=\mathbf{I}+a\mathbf{A}+b\mathbf{A}^2$，其中 $a,b$ 只依赖
> $\theta=\|\boldsymbol{\phi}\|$。因此
> $$
> \begin{aligned}
> \mathbf{J}_l(\boldsymbol{\phi})^\mathsf{T}
> &=
> \mathbf{I}-a\mathbf{A}+b\mathbf{A}^2 \\
> &=
> \mathbf{J}_l(-\boldsymbol{\phi})
> \end{aligned}
> $$
>

即：
$$
\mathbf{J}_r(\boldsymbol{\phi})
=
\mathbf{I}
-\frac{1-\cos\theta}{\theta^2}[\boldsymbol{\phi}]_\times
+\frac{\theta-\sin\theta}{\theta^3}[\boldsymbol{\phi}]_\times^2
$$

右雅可比逆为：

$$
\mathbf{J}_r^{-1}(\boldsymbol{\phi})
=
\mathbf{J}_l^{-1}(-\boldsymbol{\phi})
=
\mathbf{I}
+\frac{1}{2}[\boldsymbol{\phi}]_\times
+\left(
\frac{1}{\theta^2}
-\frac{1+\cos\theta}{2\theta\sin\theta}
\right)
[\boldsymbol{\phi}]_\times^2
$$

当 $\boldsymbol{\phi}$ 是小量时：
$$
\begin{aligned}
\mathbf{J}_l(\boldsymbol{\phi})
&\approx
\mathbf{I}+\frac{1}{2}[\boldsymbol{\phi}]_\times
+\frac{1}{6}[\boldsymbol{\phi}]_\times^2 \\
\mathbf{J}_r(\boldsymbol{\phi})
&\approx
\mathbf{I}-\frac{1}{2}[\boldsymbol{\phi}]_\times
+\frac{1}{6}[\boldsymbol{\phi}]_\times^2
\end{aligned}
$$

对应的逆在小量时为：

$$
\begin{aligned}
\mathbf{J}_l^{-1}(\boldsymbol{\phi})
&\approx
\mathbf{I}-\frac{1}{2}[\boldsymbol{\phi}]_\times
+\frac{1}{12}[\boldsymbol{\phi}]_\times^2 \\
\mathbf{J}_r^{-1}(\boldsymbol{\phi})
&\approx
\mathbf{I}+\frac{1}{2}[\boldsymbol{\phi}]_\times
+\frac{1}{12}[\boldsymbol{\phi}]_\times^2
\end{aligned}
$$

> **小量近似的证明：**
>
> 当 $\theta$ 很小时，三角函数的泰勒展开为
>
> $$
> \sin\theta
> =
> \theta-\frac{\theta^3}{6}+O(\theta^5),
> \qquad
> \cos\theta
> =
> 1-\frac{\theta^2}{2}+\frac{\theta^4}{24}+O(\theta^6)
> $$
>
> 因此
>
> $$
> \frac{1-\cos\theta}{\theta^2}
> =
> \frac{1}{2}-\frac{\theta^2}{24}+O(\theta^4),
> \qquad
> \frac{\theta-\sin\theta}{\theta^3}
> =
> \frac{1}{6}-\frac{\theta^2}{120}+O(\theta^4)
> $$
>
> 代入
>
> $$
> \mathbf{J}_l(\boldsymbol{\phi})
> =
> \mathbf{I}
> +\frac{1-\cos\theta}{\theta^2}[\boldsymbol{\phi}]_\times
> +\frac{\theta-\sin\theta}{\theta^3}[\boldsymbol{\phi}]_\times^2
> $$
>
> 并保留到 $[\boldsymbol{\phi}]_\times^2$，得到
>
> $$
> \mathbf{J}_l(\boldsymbol{\phi})
> \approx
> \mathbf{I}
> +\frac{1}{2}[\boldsymbol{\phi}]_\times
> +\frac{1}{6}[\boldsymbol{\phi}]_\times^2
> $$
>
> 又因为 $\mathbf{J}_r(\boldsymbol{\phi})=\mathbf{J}_l(-\boldsymbol{\phi})$，所以一次项变号、二次项不变：
>
> $$
> \mathbf{J}_r(\boldsymbol{\phi})
> \approx
> \mathbf{I}
> -\frac{1}{2}[\boldsymbol{\phi}]_\times
> +\frac{1}{6}[\boldsymbol{\phi}]_\times^2
> $$
>
> 对逆矩阵中的系数
>
> $$
> \frac{1}{\theta^2}
> -\frac{1+\cos\theta}{2\theta\sin\theta}
> $$
>
> 做同样的泰勒展开，可得
>
> $$
> \frac{1}{\theta^2}
> -\frac{1+\cos\theta}{2\theta\sin\theta}
> =
> \frac{1}{12}+O(\theta^2)
> $$
>
> 代入 $\mathbf{J}_l^{-1}$ 和 $\mathbf{J}_r^{-1}$，即可得到上面的逆矩阵小量近似。

一阶使用时，也经常直接取 $\mathbf{J}_l(\boldsymbol{\phi})\approx\mathbf{I}$、$\mathbf{J}_r(\boldsymbol{\phi})\approx\mathbf{I}$。

## 李群乘扰动转为李代数加扰动

### 小量左乘

若 $\delta\boldsymbol{\phi}$ 是小量，则：

$$
\operatorname{Exp}(\delta\boldsymbol{\phi})
\operatorname{Exp}(\boldsymbol{\phi})
\approx
\operatorname{Exp}\!\left(
\boldsymbol{\phi}
+\mathbf{J}_l^{-1}(\boldsymbol{\phi})\delta\boldsymbol{\phi}
\right)
$$

也就是：

$$
\operatorname{Log}\!\left(
\operatorname{Exp}(\delta\boldsymbol{\phi})
\operatorname{Exp}(\boldsymbol{\phi})
\right)
\approx
\boldsymbol{\phi}
+\mathbf{J}_l^{-1}(\boldsymbol{\phi})\delta\boldsymbol{\phi}
$$

左乘的小量用左雅可比逆搬到李代数加法中。

### 小量右乘

若 $\delta\boldsymbol{\phi}$ 是小量，则：

$$
\operatorname{Exp}(\boldsymbol{\phi})
\operatorname{Exp}(\delta\boldsymbol{\phi})
\approx
\operatorname{Exp}\!\left(
\boldsymbol{\phi}
+\mathbf{J}_r^{-1}(\boldsymbol{\phi})\delta\boldsymbol{\phi}
\right)
$$

也就是：

$$
\operatorname{Log}\!\left(
\operatorname{Exp}(\boldsymbol{\phi})
\operatorname{Exp}(\delta\boldsymbol{\phi})
\right)
\approx
\boldsymbol{\phi}
+\mathbf{J}_r^{-1}(\boldsymbol{\phi})\delta\boldsymbol{\phi}
$$

右乘的小量用右雅可比逆搬到李代数加法中。

> **和 BCH 一阶截断的关系：**
>
> 若 $\boldsymbol{\phi}$ 本身也很小，则
>
> $$
> \mathbf{J}_l^{-1}(\boldsymbol{\phi})\approx\mathbf{I},
> \qquad
> \mathbf{J}_r^{-1}(\boldsymbol{\phi})\approx\mathbf{I}
> $$
>
> 于是左乘扰动和右乘扰动分别退化为
>
> $$
> \begin{aligned}
> \operatorname{Log}\!\left(
> \operatorname{Exp}(\delta\boldsymbol{\phi})
> \operatorname{Exp}(\boldsymbol{\phi})
> \right)
> &\approx
> \boldsymbol{\phi}+\delta\boldsymbol{\phi} \\
> \operatorname{Log}\!\left(
> \operatorname{Exp}(\boldsymbol{\phi})
> \operatorname{Exp}(\delta\boldsymbol{\phi})
> \right)
> &\approx
> \boldsymbol{\phi}+\delta\boldsymbol{\phi}
> \end{aligned}
> $$
>
> 如果需要保留二阶精度，再加入 $\frac{1}{2}\boldsymbol{\phi}\times\delta\boldsymbol{\phi}$ 这类李括号项。

## 李代数加扰动转为李群乘扰动

上一节处理的是“李群上乘一个小扰动，等价于李代数上加什么”。反过来，如果李代数中出现 $\boldsymbol{\phi}+\delta\boldsymbol{\phi}$，也可以把它改写成李群上的左乘或右乘扰动。

### 加扰动转为左乘

若 $\delta\boldsymbol{\phi}$ 是小量，则：

$$
\operatorname{Exp}(\boldsymbol{\phi}+\delta\boldsymbol{\phi})
\approx
\operatorname{Exp}\!\left(
\mathbf{J}_l(\boldsymbol{\phi})\delta\boldsymbol{\phi}
\right)
\operatorname{Exp}(\boldsymbol{\phi})
$$

### 加扰动转为右乘

若 $\delta\boldsymbol{\phi}$ 是小量，则：

$$
\operatorname{Exp}(\boldsymbol{\phi}+\delta\boldsymbol{\phi})
\approx
\operatorname{Exp}(\boldsymbol{\phi})
\operatorname{Exp}\!\left(
\mathbf{J}_r(\boldsymbol{\phi})\delta\boldsymbol{\phi}
\right)
$$

这两式和上一节互为逆过程：从乘法扰动变成加法扰动时用雅可比逆；从加法扰动变成乘法扰动时用雅可比本身。

## 小结

上面的近似式成立需要注意以下几点：

1. 被线性化的小扰动必须足够小，例如 $\delta\boldsymbol{\phi}$、$\delta\boldsymbol{\theta}$、噪声项或零偏误差。
2. 大量 $\boldsymbol{\phi}$ 不一定要很小，但应位于 $\operatorname{SO}(3)$ 对数映射的良好区域，避免接近 $\theta=\pi$ 时出现分支不连续或雅可比病态。
3. 如果两个量都很小，可以直接使用 BCH 二阶截断；如果一个量大、另一个量小，应优先使用左/右雅可比形式。
4. 左乘扰动和右乘扰动不能混用。小量在李群左侧相乘时对应 $\mathbf{J}_l$，在右侧相乘时对应 $\mathbf{J}_r$。
5. 符号里的 $\operatorname{Exp}(\boldsymbol{\phi})$ 表示向量形式指数映射；$\exp([\boldsymbol{\phi}]_\times)$ 表示矩阵指数映射。二者在 $\operatorname{SO}(3)$ 中等价，但写法含义不同。

常用关系可以总结为：

$$
\begin{aligned}
\operatorname{Exp}(\delta\boldsymbol{\phi})\operatorname{Exp}(\boldsymbol{\phi})
&\approx
\operatorname{Exp}\!\left(
\boldsymbol{\phi}
+\mathbf{J}_l^{-1}(\boldsymbol{\phi})\delta\boldsymbol{\phi}
\right) \\
\operatorname{Exp}(\boldsymbol{\phi})\operatorname{Exp}(\delta\boldsymbol{\phi})
&\approx
\operatorname{Exp}\!\left(
\boldsymbol{\phi}
+\mathbf{J}_r^{-1}(\boldsymbol{\phi})\delta\boldsymbol{\phi}
\right) \\
\operatorname{Exp}(\boldsymbol{\phi}+\delta\boldsymbol{\phi})
&\approx
\operatorname{Exp}\!\left(
\mathbf{J}_l(\boldsymbol{\phi})\delta\boldsymbol{\phi}
\right)\operatorname{Exp}(\boldsymbol{\phi}) \\
\operatorname{Exp}(\boldsymbol{\phi}+\delta\boldsymbol{\phi})
&\approx
\operatorname{Exp}(\boldsymbol{\phi})
\operatorname{Exp}\!\left(
\mathbf{J}_r(\boldsymbol{\phi})\delta\boldsymbol{\phi}
\right)
\end{aligned}
$$

一句话记忆：**李群乘法转李代数加法，用左/右雅可比的逆；李代数加法转李群乘法，用左/右雅可比本身。小量左乘用左雅可比，小量右乘用右雅可比。**

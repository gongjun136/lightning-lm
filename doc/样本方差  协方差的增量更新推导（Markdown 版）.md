# 样本方差 / 协方差的增量更新推导（Markdown 版）

本文推导在「已有 \(m\) 个样本的均值和方差/协方差已知」的前提下，  
当新增一个样本点时，如何**直接从定义**得到方差/协方差的增量更新公式。

---

## 1 一维样本方差的增量更新推导

假设已经有 \(m\) 个标量样本 \(x_1,\dots,x_m\)，其**样本均值**和**样本方差**（带贝塞尔修正）为：

$$
\mu_m = \frac{1}{m}\sum_{i=1}^m x_i,
\qquad
s_m^2 = \frac{1}{m-1}\sum_{i=1}^m (x_i-\mu_m)^2.
$$

现在来了一个新样本 \(p\)，总样本数变为 \(m+1\)。  
记更新后的样本均值、方差为 \(\mu_{m+1}, s_{m+1}^2\)。

---

### 1.1 均值的增量更新

新均值可以直接写成：

$$
\mu_{m+1} 
= \frac{1}{m+1}\left(\sum_{i=1}^m x_i + p\right)
= \frac{m\mu_m + p}{m+1}.
$$

定义“增量”：

$$
\delta = p - \mu_m,
$$

则有常用的增量形式：

$$
\boxed{
\mu_{m+1} = \mu_m + \frac{1}{m+1}\,\delta
}
$$

---

### 1.2 方差的直接增量推导

样本方差是由**平方偏差和**得到的。先定义：

$$
S_m = \sum_{i=1}^m (x_i-\mu_m)^2
\quad\Longrightarrow\quad
s_m^2 = \frac{S_m}{m-1}.
$$

加入新点 \(p\) 后，新的平方偏差和为：

$$
S_{m+1} = \sum_{i=1}^{m+1} (x_i - \mu_{m+1})^2
= \underbrace{\sum_{i=1}^{m} (x_i - \mu_{m+1})^2}_{\text{旧 }m\text{ 个点相对新均值的偏差}}
+ (p - \mu_{m+1})^2.
$$

把「旧点相对新均值」写成「旧点相对旧均值 + 均值平移」：

$$
x_i - \mu_{m+1}
= (x_i - \mu_m) + (\mu_m - \mu_{m+1}).
$$

记

$$
a_i = x_i - \mu_m,\qquad
u = \mu_m - \mu_{m+1},
$$

于是：

$$
\begin{aligned}
\sum_{i=1}^{m} (x_i - \mu_{m+1})^2
&= \sum_{i=1}^{m} (a_i + u)^2 \\
&= \sum_{i=1}^{m}\left(a_i^2 + 2a_i u + u^2\right) \\
&= \sum_{i=1}^{m} a_i^2 
+ 2u\sum_{i=1}^{m} a_i
+ m u^2.
\end{aligned}
$$

注意到：

$$
\sum_{i=1}^m a_i
= \sum_{i=1}^m (x_i-\mu_m)
= \left(\sum_{i=1}^m x_i\right) - m\mu_m
= 0,
$$

因此 \(\sum 2a_i u = 2u\sum a_i = 0\)，上式化简为：

$$
\sum_{i=1}^{m} (x_i - \mu_{m+1})^2
= \sum_{i=1}^{m} a_i^2 + m u^2
= S_m + m u^2.
$$

---

再看新点那一项。根据均值增量公式：

$$
\mu_{m+1} = \mu_m + \frac{1}{m+1}\delta
\quad\Longrightarrow\quad
\mu_m - \mu_{m+1} = -\,\frac{1}{m+1}\delta = u.
$$

于是：

$$
\begin{aligned}
p - \mu_{m+1}
&= (p - \mu_m) + (\mu_m - \mu_{m+1}) \\
&= \delta + u
= \delta - \frac{1}{m+1}\delta
= \frac{m}{m+1}\delta,
\end{aligned}
$$

因此：

$$
(p - \mu_{m+1})^2
= \left(\frac{m}{m+1}\right)^2 \delta^2.
$$

---

将两部分合并，得到：

$$
\begin{aligned}
S_{m+1}
&= \underbrace{\bigl(S_m + m u^2\bigr)}_{\sum_{i=1}^{m} (x_i - \mu_{m+1})^2}
+ \underbrace{\left(\frac{m}{m+1}\right)^2 \delta^2}_{(p - \mu_{m+1})^2} \\
&= S_m + m\left(\frac{\delta}{m+1}\right)^2
+ \frac{m^2}{(m+1)^2}\delta^2 \\
&= S_m + \left(\frac{m}{(m+1)^2} + \frac{m^2}{(m+1)^2}\right)\delta^2 \\
&= S_m + \frac{m(m+1)}{(m+1)^2}\delta^2 \\
&= S_m + \frac{m}{m+1}\,\delta^2.
\end{aligned}
$$

更新后的样本方差为：

$$
\begin{aligned}
s_{m+1}^2
&= \frac{S_{m+1}}{m}
= \frac{1}{m}\left(S_m + \frac{m}{m+1}\delta^2\right) \\
&= \frac{S_m}{m}
+ \frac{1}{m+1}\delta^2
= \frac{m-1}{m}\,s_m^2
+ \frac{1}{m+1}\delta^2.
\end{aligned}
$$

因此得到一维方差的增量更新公式：

$$
\boxed{
s_{m+1}^2
=
\frac{m-1}{m}\,s_m^2
+ \frac{1}{m+1}\bigl(p-\mu_m\bigr)^2
}
$$

这里的 \(s_m^2\) 和 \(s_{m+1}^2\) 都是带贝塞尔修正的**样本方差**。

---

## 2 多维协方差矩阵的直接增量推导

多维情况只是把标量平方变成「外积矩阵」，推导步骤完全类似。

已有 \(m\) 个向量样本 \(\mathbf{x}_1,\dots,\mathbf{x}_m \in \mathbb{R}^d\)，  
样本均值与样本协方差定义为：

$$
\boldsymbol{\mu}_m = \frac{1}{m}\sum_{i=1}^m \mathbf{x}_i,
\qquad
\boldsymbol{\Sigma}_m = \frac{1}{m-1}\sum_{i=1}^{m}
(\mathbf{x}_i-\boldsymbol{\mu}_m)(\mathbf{x}_i-\boldsymbol{\mu}_m)^{\mathrm T}.
$$

同样定义“平方偏差和矩阵”：

$$
\mathbf{S}_m
= \sum_{i=1}^{m}
(\mathbf{x}_i-\boldsymbol{\mu}_m)(\mathbf{x}_i-\boldsymbol{\mu}_m)^{\mathrm T}
\quad\Longrightarrow\quad
\boldsymbol{\Sigma}_m = \frac{1}{m-1}\mathbf{S}_m.
$$

来了一个新样本 \(\mathbf{p}\)，记：

$$
\delta = \mathbf{p} - \boldsymbol{\mu}_m.
$$

---

### 2.1 均值的增量更新

和一维完全一致：

$$
\boxed{
\boldsymbol{\mu}_{m+1}
= \boldsymbol{\mu}_m + \frac{1}{m+1}\,\delta
}
$$

---

### 2.2 平方偏差和矩阵的增量更新

新样本加入后：

$$
\mathbf{S}_{m+1}
= \sum_{i=1}^{m+1}
(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})^{\mathrm T}.
$$

将其拆成两部分：

$$
\mathbf{S}_{m+1}
= \underbrace{
\sum_{i=1}^{m}(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})
(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})^{\mathrm T}
}_{\text{旧 }m\text{ 个点相对新均值}}
+
(\mathbf{p}-\boldsymbol{\mu}_{m+1})
(\mathbf{p}-\boldsymbol{\mu}_{m+1})^{\mathrm T}.
$$

设

$$
\mathbf{a}_i = \mathbf{x}_i-\boldsymbol{\mu}_m,\qquad
\mathbf{u} = \boldsymbol{\mu}_m - \boldsymbol{\mu}_{m+1},
$$

则有

$$
\mathbf{x}_i-\boldsymbol{\mu}_{m+1} = \mathbf{a}_i + \mathbf{u}.
$$

于是第一部分为：

$$
\begin{aligned}
\sum_{i=1}^{m}(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})
(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})^{\mathrm T}
&= \sum_{i=1}^{m}
(\mathbf{a}_i+\mathbf{u})(\mathbf{a}_i+\mathbf{u})^{\mathrm T} \\
&= \sum_{i=1}^{m}
\Bigl(
\mathbf{a}_i \mathbf{a}_i^{\mathrm T}
+ \mathbf{a}_i \mathbf{u}^{\mathrm T}
+ \mathbf{u}\mathbf{a}_i^{\mathrm T}
+ \mathbf{u}\mathbf{u}^{\mathrm T}
\Bigr).
\end{aligned}
$$

注意到：

$$
\sum_{i=1}^m \mathbf{a}_i
= \sum_{i=1}^m (\mathbf{x}_i-\boldsymbol{\mu}_m)
= \left(\sum_{i=1}^m \mathbf{x}_i\right) - m\boldsymbol{\mu}_m
= \mathbf{0},
$$

因此两项线性项抵消，得到：

$$
\sum_{i=1}^{m}(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})
(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})^{\mathrm T}
= \mathbf{S}_m + m \mathbf{u}\mathbf{u}^{\mathrm T}.
$$

---

根据均值增量公式：

$$
\boldsymbol{\mu}_{m+1}
= \boldsymbol{\mu}_m + \frac{1}{m+1}\delta
\quad\Longrightarrow\quad
\mathbf{u}
= \boldsymbol{\mu}_m - \boldsymbol{\mu}_{m+1}
= -\,\frac{1}{m+1}\delta.
$$

新点相对新均值的偏差为：

$$
\begin{aligned}
\mathbf{p}-\boldsymbol{\mu}_{m+1}
&= (\mathbf{p}-\boldsymbol{\mu}_m)
+ (\boldsymbol{\mu}_m - \boldsymbol{\mu}_{m+1}) \\
&= \delta + \mathbf{u}
= \delta - \frac{1}{m+1}\delta
= \frac{m}{m+1}\,\delta.
\end{aligned}
$$

因此第二部分为：

$$
(\mathbf{p}-\boldsymbol{\mu}_{m+1})
(\mathbf{p}-\boldsymbol{\mu}_{m+1})^{\mathrm T}
= \left(\frac{m}{m+1}\right)^2 \delta\delta^{\mathrm T}.
$$

---

将两部分合并，得到：

$$
\begin{aligned}
\mathbf{S}_{m+1}
&= \bigl(\mathbf{S}_m + m \mathbf{u}\mathbf{u}^{\mathrm T}\bigr)
+ \left(\frac{m}{m+1}\right)^2 \delta\delta^{\mathrm T} \\
&= \mathbf{S}_m
+ m\left(\frac{\delta}{m+1}\right)\left(\frac{\delta}{m+1}\right)^{\mathrm T}
+ \frac{m^2}{(m+1)^2}\delta\delta^{\mathrm T} \\
&= \mathbf{S}_m
+ \left(
\frac{m}{(m+1)^2}
+ \frac{m^2}{(m+1)^2}
\right)\delta\delta^{\mathrm T} \\
&= \mathbf{S}_m
+ \frac{m(m+1)}{(m+1)^2}\delta\delta^{\mathrm T} \\
&= \mathbf{S}_m
+ \frac{m}{m+1}\,\delta\delta^{\mathrm T}.
\end{aligned}
$$

---

### 2.3 协方差矩阵的增量更新公式

更新后的样本协方差为：

$$
\boldsymbol{\Sigma}_{m+1}
= \frac{1}{m}\mathbf{S}_{m+1}
= \frac{1}{m}\left(
\mathbf{S}_m + \frac{m}{m+1}\delta\delta^{\mathrm T}
\right).
$$

利用

$$
\mathbf{S}_m = (m-1)\boldsymbol{\Sigma}_m,
$$

得到：

$$
\begin{aligned}
\boldsymbol{\Sigma}_{m+1}
&= \frac{(m-1)}{m}\boldsymbol{\Sigma}_m
+ \frac{1}{m+1}\delta\delta^{\mathrm T}.
\end{aligned}
$$

最终得到多维样本协方差的增量更新公式：

$$
\boxed{
\boldsymbol{\Sigma}_{m+1}
=
\frac{m-1}{m}\,\boldsymbol{\Sigma}_m
+
\frac{1}{m+1}
\bigl(\mathbf{p}-\boldsymbol{\mu}_m\bigr)
\bigl(\mathbf{p}-\boldsymbol{\mu}_m\bigr)^{\mathrm T}
}
$$

> 如果只关心每一维的方差（对角线），可以将  
> \(\delta\delta^{\mathrm T}\) 换成逐元素乘积 \(\delta \odot \delta\)，  
> 即在实现中使用 `delta.cwiseProduct(delta)`，只保留对角项。

---





# 样本方差 / 协方差的增量更新推导

本文推导在「已有 \(m\) 个样本的均值和方差/协方差已知」的前提下，  
当新增一个样本点时，如何**直接从定义**得到方差/协方差的增量更新公式。

记号约定：

- **向量**：用黑体加粗小写字母表示，例如：\(\mathbf{x}, \mathbf{p}, \mathbf{a}, \mathbf{u}, \boldsymbol{\mu}, \boldsymbol{\delta}\)。
- **矩阵**：用黑体加粗大写字母表示，例如：\(\mathbf{S}, \boldsymbol{\Sigma}\)（大写拉丁或希腊字母）。

---

## 1 一维样本方差的增量更新推导

一维情形下，样本是标量，不使用向量/矩阵记号。

假设已经有 \(m\) 个标量样本 \(x_1,\dots,x_m\)，其**样本均值**和**样本方差**（带贝塞尔修正）为：

$$
\mu_m = \frac{1}{m}\sum_{i=1}^m x_i,
\qquad
s_m^2 = \frac{1}{m-1}\sum_{i=1}^m (x_i-\mu_m)^2.
$$

现在来了一个新样本 \(p\)，总样本数变为 \(m+1\)。  
记更新后的样本均值、方差为 \(\mu_{m+1}, s_{m+1}^2\)。

---

### 1.1 均值的增量更新

新均值可以直接写成：

$$
\mu_{m+1} 
= \frac{1}{m+1}\left(\sum_{i=1}^m x_i + p\right)
= \frac{m\mu_m + p}{m+1}.
$$

定义“增量”：

$$
\delta = p - \mu_m,
$$

则有常用的增量形式：

$$
\boxed{
\mu_{m+1} = \mu_m + \frac{1}{m+1}\,\delta
}
$$

---

### 1.2 方差的直接增量推导

样本方差是由**平方偏差和**得到的。先定义：

$$
S_m = \sum_{i=1}^m (x_i-\mu_m)^2
\quad\Longrightarrow\quad
s_m^2 = \frac{S_m}{m-1}.
$$

加入新点 \(p\) 后，新的平方偏差和为：

$$
S_{m+1} = \sum_{i=1}^{m+1} (x_i - \mu_{m+1})^2
= \underbrace{\sum_{i=1}^{m} (x_i - \mu_{m+1})^2}_{\text{旧 }m\text{ 个点相对新均值的偏差}}
+ (p - \mu_{m+1})^2.
$$

把「旧点相对新均值」写成「旧点相对旧均值 + 均值平移」：

$$
x_i - \mu_{m+1}
= (x_i - \mu_m) + (\mu_m - \mu_{m+1}).
$$

记

$$
a_i = x_i - \mu_m,\qquad
u = \mu_m - \mu_{m+1},
$$

于是：

$$
\begin{aligned}
\sum_{i=1}^{m} (x_i - \mu_{m+1})^2
&= \sum_{i=1}^{m} (a_i + u)^2 \\
&= \sum_{i=1}^{m}\left(a_i^2 + 2a_i u + u^2\right) \\
&= \sum_{i=1}^{m} a_i^2 
+ 2u\sum_{i=1}^{m} a_i
+ m u^2.
\end{aligned}
$$

注意到：

$$
\sum_{i=1}^m a_i
= \sum_{i=1}^m (x_i-\mu_m)
= \left(\sum_{i=1}^m x_i\right) - m\mu_m
= 0,
$$

因此 \(\sum 2a_i u = 2u\sum a_i = 0\)，上式化简为：

$$
\sum_{i=1}^{m} (x_i - \mu_{m+1})^2
= \sum_{i=1}^{m} a_i^2 + m u^2
= S_m + m u^2.
$$

---

再看新点那一项。根据均值增量公式：

$$
\mu_{m+1} = \mu_m + \frac{1}{m+1}\delta
\quad\Longrightarrow\quad
\mu_m - \mu_{m+1} = -\,\frac{1}{m+1}\delta = u.
$$

于是：

$$
\begin{aligned}
p - \mu_{m+1}
&= (p - \mu_m) + (\mu_m - \mu_{m+1}) \\
&= \delta + u
= \delta - \frac{1}{m+1}\delta
= \frac{m}{m+1}\delta,
\end{aligned}
$$

因此：

$$
(p - \mu_{m+1})^2
= \left(\frac{m}{m+1}\right)^2 \delta^2.
$$

---

将两部分合并，得到：

$$
\begin{aligned}
S_{m+1}
&= \underbrace{\bigl(S_m + m u^2\bigr)}_{\sum_{i=1}^{m} (x_i - \mu_{m+1})^2}
+ \underbrace{\left(\frac{m}{m+1}\right)^2 \delta^2}_{(p - \mu_{m+1})^2} \\
&= S_m + m\left(\frac{\delta}{m+1}\right)^2
+ \frac{m^2}{(m+1)^2}\delta^2 \\
&= S_m + \left(\frac{m}{(m+1)^2} + \frac{m^2}{(m+1)^2}\right)\delta^2 \\
&= S_m + \frac{m(m+1)}{(m+1)^2}\delta^2 \\
&= S_m + \frac{m}{m+1}\,\delta^2.
\end{aligned}
$$

更新后的样本方差为：

$$
\begin{aligned}
s_{m+1}^2
&= \frac{S_{m+1}}{m}
= \frac{1}{m}\left(S_m + \frac{m}{m+1}\delta^2\right) \\
&= \frac{S_m}{m}
+ \frac{1}{m+1}\delta^2
= \frac{m-1}{m}\,s_m^2
+ \frac{1}{m+1}\delta^2.
\end{aligned}
$$

因此得到一维方差的增量更新公式：

$$
\boxed{
s_{m+1}^2
=
\frac{m-1}{m}\,s_m^2
+ \frac{1}{m+1}\bigl(p-\mu_m\bigr)^2
}
$$

这里的 \(s_m^2\) 和 \(s_{m+1}^2\) 都是带贝塞尔修正的**样本方差**。

---

## 2 多维协方差矩阵的直接增量推导

多维情况下，样本为向量，方差推广为协方差矩阵。

已有 \(m\) 个向量样本 \(\mathbf{x}_1,\dots,\mathbf{x}_m \in \mathbb{R}^d\)，  
样本均值与样本协方差定义为：

$$
\boldsymbol{\mu}_m = \frac{1}{m}\sum_{i=1}^m \mathbf{x}_i,
\qquad
\boldsymbol{\Sigma}_m = \frac{1}{m-1}\sum_{i=1}^{m}
(\mathbf{x}_i-\boldsymbol{\mu}_m)(\mathbf{x}_i-\boldsymbol{\mu}_m)^{\mathrm T}.
$$

同样定义“平方偏差和矩阵”：

$$
\mathbf{S}_m
= \sum_{i=1}^{m}
(\mathbf{x}_i-\boldsymbol{\mu}_m)(\mathbf{x}_i-\boldsymbol{\mu}_m)^{\mathrm T}
\quad\Longrightarrow\quad
\boldsymbol{\Sigma}_m = \frac{1}{m-1}\mathbf{S}_m.
$$

来了一个新样本向量 \(\mathbf{p}\)，记：

$$
\boldsymbol{\delta} = \mathbf{p} - \boldsymbol{\mu}_m.
$$

---

### 2.1 均值的增量更新

与一维完全一致，只是改为向量形式：

$$
\boxed{
\boldsymbol{\mu}_{m+1}
= \boldsymbol{\mu}_m + \frac{1}{m+1}\,\boldsymbol{\delta}
}
$$

---

### 2.2 平方偏差和矩阵的增量更新

新样本加入后：

$$
\mathbf{S}_{m+1}
= \sum_{i=1}^{m+1}
(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})^{\mathrm T}.
$$

将其拆成两部分：

$$
\mathbf{S}_{m+1}
= \underbrace{
\sum_{i=1}^{m}(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})
(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})^{\mathrm T}
}_{\text{旧 }m\text{ 个点相对新均值}}
+
(\mathbf{p}-\boldsymbol{\mu}_{m+1})
(\mathbf{p}-\boldsymbol{\mu}_{m+1})^{\mathrm T}.
$$

设

$$
\mathbf{a}_i = \mathbf{x}_i-\boldsymbol{\mu}_m,\qquad
\mathbf{u} = \boldsymbol{\mu}_m - \boldsymbol{\mu}_{m+1},
$$

则有

$$
\mathbf{x}_i-\boldsymbol{\mu}_{m+1} = \mathbf{a}_i + \mathbf{u}.
$$

于是第一部分为：

$$
\begin{aligned}
\sum_{i=1}^{m}(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})
(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})^{\mathrm T}
&= \sum_{i=1}^{m}
(\mathbf{a}_i+\mathbf{u})(\mathbf{a}_i+\mathbf{u})^{\mathrm T} \\
&= \sum_{i=1}^{m}
\Bigl(
\mathbf{a}_i \mathbf{a}_i^{\mathrm T}
+ \mathbf{a}_i \mathbf{u}^{\mathrm T}
+ \mathbf{u}\mathbf{a}_i^{\mathrm T}
+ \mathbf{u}\mathbf{u}^{\mathrm T}
\Bigr).
\end{aligned}
$$

注意到：

$$
\sum_{i=1}^m \mathbf{a}_i
= \sum_{i=1}^m (\mathbf{x}_i-\boldsymbol{\mu}_m)
= \left(\sum_{i=1}^m \mathbf{x}_i\right) - m\boldsymbol{\mu}_m
= \mathbf{0},
$$

因此两项线性项抵消，得到：

$$
\sum_{i=1}^{m}(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})
(\mathbf{x}_i-\boldsymbol{\mu}_{m+1})^{\mathrm T}
= \mathbf{S}_m + m \mathbf{u}\mathbf{u}^{\mathrm T}.
$$

---

根据均值增量公式：

$$
\boldsymbol{\mu}_{m+1}
= \boldsymbol{\mu}_m + \frac{1}{m+1}\boldsymbol{\delta}
\quad\Longrightarrow\quad
\mathbf{u}
= \boldsymbol{\mu}_m - \boldsymbol{\mu}_{m+1}
= -\,\frac{1}{m+1}\boldsymbol{\delta}.
$$

新点相对新均值的偏差为：

$$
\begin{aligned}
\mathbf{p}-\boldsymbol{\mu}_{m+1}
&= (\mathbf{p}-\boldsymbol{\mu}_m)
+ (\boldsymbol{\mu}_m - \boldsymbol{\mu}_{m+1}) \\
&= \boldsymbol{\delta} + \mathbf{u}
= \boldsymbol{\delta} - \frac{1}{m+1}\boldsymbol{\delta}
= \frac{m}{m+1}\,\boldsymbol{\delta}.
\end{aligned}
$$

因此第二部分为：

$$
(\mathbf{p}-\boldsymbol{\mu}_{m+1})
(\mathbf{p}-\boldsymbol{\mu}_{m+1})^{\mathrm T}
= \left(\frac{m}{m+1}\right)^2 \boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T}.
$$

---

将两部分合并，得到：

$$
\begin{aligned}
\mathbf{S}_{m+1}
&= \bigl(\mathbf{S}_m + m \mathbf{u}\mathbf{u}^{\mathrm T}\bigr)
+ \left(\frac{m}{m+1}\right)^2 \boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T} \\
&= \mathbf{S}_m
+ m\left(\frac{\boldsymbol{\delta}}{m+1}\right)\left(\frac{\boldsymbol{\delta}}{m+1}\right)^{\mathrm T}
+ \frac{m^2}{(m+1)^2}\boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T} \\
&= \mathbf{S}_m
+ \left(
\frac{m}{(m+1)^2}
+ \frac{m^2}{(m+1)^2}
\right)\boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T} \\
&= \mathbf{S}_m
+ \frac{m(m+1)}{(m+1)^2}\boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T} \\
&= \mathbf{S}_m
+ \frac{m}{m+1}\,\boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T}.
\end{aligned}
$$

---

### 2.3 协方差矩阵的增量更新公式

更新后的样本协方差为：

$$
\boldsymbol{\Sigma}_{m+1}
= \frac{1}{m}\mathbf{S}_{m+1}
= \frac{1}{m}\left(
\mathbf{S}_m + \frac{m}{m+1}\boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T}
\right).
$$

利用

$$
\mathbf{S}_m = (m-1)\boldsymbol{\Sigma}_m,
$$

得到：

$$
\begin{aligned}
\boldsymbol{\Sigma}_{m+1}
&= \frac{(m-1)}{m}\boldsymbol{\Sigma}_m
+ \frac{1}{m+1}\boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T}.
\end{aligned}
$$

最终得到多维样本协方差的增量更新公式：

$$
\boxed{
\boldsymbol{\Sigma}_{m+1}
=
\frac{m-1}{m}\,\boldsymbol{\Sigma}_m
+
\frac{1}{m+1}
\bigl(\mathbf{p}-\boldsymbol{\mu}_m\bigr)
\bigl(\mathbf{p}-\boldsymbol{\mu}_m\bigr)^{\mathrm T}
}
$$

> 如果只关心每一维的方差（对角线），可以将  
> \(\boldsymbol{\delta}\boldsymbol{\delta}^{\mathrm T}\) 换成逐元素乘积 \(\boldsymbol{\delta} \odot \boldsymbol{\delta}\)，  
> 即在实现中使用 `delta.cwiseProduct(delta)`，只保留对角项。

---









### n = 1 时的增量更新公式（向量粗体、矩阵粗体）

当只新增 1 个样本点 \(\mathbf{p}\) 时，更新前后记号为：

- 旧样本数：\(m\)
- 旧均值：\(\boldsymbol{\mu}_m \in \mathbb{R}^d\)
- 旧协方差矩阵（样本协方差）：\(\boldsymbol{\Sigma}_m \in \mathbb{R}^{d\times d}\)
- 新点：\(\mathbf{p} \in \mathbb{R}^d\)

---

**1. 新均值**

\[
$$
\boldsymbol{\mu}_{m+1}
=
\frac{m\boldsymbol{\mu}_m + \mathbf{p}}{m+1}
=
\boldsymbol{\mu}_m + \frac{1}{m+1}\bigl(\mathbf{p}-\boldsymbol{\mu}_m\bigr)
$$
\]

---

**2. 新协方差矩阵（样本协方差，带贝塞尔修正）**

先记增量向量：

\[
$$
\boldsymbol{\delta} = \mathbf{p} - \boldsymbol{\mu}_m
$$
\]

则协方差更新为：

\[
$$
\boldsymbol{\Sigma}_{m+1}
=
\frac{m-1}{m}\,\boldsymbol{\Sigma}_m
+
\frac{1}{m+1}\,
\boldsymbol{\delta}\,\boldsymbol{\delta}^{\mathrm T}
$$
\]

也可以直接写成：

\[
$$
\boxed{
\boldsymbol{\Sigma}_{m+1}
=
\frac{m-1}{m}\,\boldsymbol{\Sigma}_m
+
\frac{1}{m+1}
\bigl(\mathbf{p}-\boldsymbol{\mu}_m\bigr)
\bigl(\mathbf{p}-\boldsymbol{\mu}_m\bigr)^{\mathrm T}
}
$$
\]



当只新增 1 个样本点 \(\mathbf{p}\) 时，令
\(\boldsymbol{\delta} = \mathbf{p} - \boldsymbol{\mu}_m\)，则协方差的单点增量更新公式为：

\[
$$
\boldsymbol{\Sigma}_{m+1}
=
\frac{m-1}{m}\,\boldsymbol{\Sigma}_m
+
\frac{1}{m+1}\,
\boldsymbol{\delta}\,\boldsymbol{\delta}^{\mathrm T},
\quad
\boldsymbol{\delta} = \mathbf{p} - \boldsymbol{\mu}_m.
$$
\]
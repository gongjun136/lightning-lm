# ESKF

本文从 [运动学](../运动学/运动学.md) 的离散模型出发，把 IMU 测量模型代入刚体运动学，得到 IMU 驱动的名义状态传播；推导误差状态的线性传播、协方差预测和观测更新。这里默认 IMU 坐标系为 $\mathcal{F}_i$，世界系为 $\mathcal{F}_w$，旋转矩阵 $\mathbf{R}_{wi}$ 表示从 IMU 系到世界系的旋转，满足 ${}^{w}\mathbf{v}=\mathbf{R}_{wi}\,{}^{i}\mathbf{v}$。

## 1 IMU 离散运动模型

### 1.1 运动学

运动学文档给出的刚体**连续模型**为：

$$
\begin{aligned}
{}^{w}\dot{\mathbf{p}} &= {}^{w}\mathbf{v} \\
{}^{w}\dot{\mathbf{v}} &= {}^{w}\mathbf{a} \\
\dot{\mathbf{R}}_{wi} &= \mathbf{R}_{wi}\left[{}^{i}\boldsymbol{\omega}\right]_{\times}
\end{aligned}
$$

其中 ${}^{w}\mathbf{p}$ 是 IMU 原点在世界系下的位置坐标，${}^{w}\mathbf{v}$ 是速度，${}^{w}\mathbf{a}$ 是世界系下的线加速度，${}^{i}\boldsymbol{\omega}$ 是 IMU 系下的角速度。后续要把 IMU 测量模型代入离散运动学，所以先把运动学**离散模型**改写成第 $k$ 个时间步到第 $k+1$ 个时间步的形式：

$$
\begin{aligned}
{}^{w}\mathbf{p}_{k+1}
&=
{}^{w}\mathbf{p}_{k}
+{}^{w}\mathbf{v}_{k}\Delta t
+\frac{1}{2}{}^{w}\mathbf{a}_{k}\Delta t^2 \\
{}^{w}\mathbf{v}_{k+1}
&=
{}^{w}\mathbf{v}_{k}
+{}^{w}\mathbf{a}_{k}\Delta t \\
\mathbf{R}_{wi,k+1}
&=
\mathbf{R}_{wi,k}
\operatorname{Exp}\!\left({}^{i}\boldsymbol{\omega}_{k}\Delta t\right) \\
\mathbf{q}_{wi,k+1}
&=
\mathbf{q}_{wi,k}\otimes
\begin{bmatrix}
\cos\frac{\left\|{}^{i}\boldsymbol{\omega}_{k}\right\|\Delta t}{2} \\
\frac{{}^{i}\boldsymbol{\omega}_{k}}{\left\|{}^{i}\boldsymbol{\omega}_{k}\right\|}
\sin\frac{\left\|{}^{i}\boldsymbol{\omega}_{k}\right\|\Delta t}{2}
\end{bmatrix}
\approx
\mathbf{q}_{wi,k}\otimes
\begin{bmatrix}
1 \\
\frac{1}{2}{}^{i}\boldsymbol{\omega}_{k}\Delta t
\end{bmatrix}
\end{aligned}
$$

本文后续主要使用旋转矩阵形式，四元数形式只用于说明同一个旋转增量也可以用单位四元数表达。若 $\left\|{}^{i}\boldsymbol{\omega}_{k}\right\|$ 很小，上式中的四元数向量部按极限处理。

上面的位置更新写成常加速度离散形式。后面推导误差状态传播矩阵时，为了让主线更清楚，会先采用一阶欧拉形式 ${}^{w}\mathbf{p}_{k+1}={}^{w}\mathbf{p}_{k}+{}^{w}\mathbf{v}_{k}\Delta t$；如果实际实现使用常加速度位置更新，只需在位置误差雅可比中补上对应的 $\frac{1}{2}\Delta t^2$ 项。

### 1.2 IMU 测量模型

IMU 不直接给出无噪声的 ${}^{i}\boldsymbol{\omega}$ 和 ${}^{i}\mathbf{a}$，而是给出带零偏和白噪声的测量值：

$$
\begin{aligned}
{}^{i}\overline{\boldsymbol{\omega}}
&= {}^{i}\boldsymbol{\omega}+{}^{i}\mathbf{b}_g+{}^{i}\boldsymbol{\eta}_g \\
{}^{i}\overline{\mathbf{a}}
&= {}^{i}\mathbf{a}+{}^{i}\mathbf{b}_a+{}^{i}\boldsymbol{\eta}_a
\end{aligned}
$$

因此：

$$
\begin{aligned}
{}^{i}\boldsymbol{\omega}
&= {}^{i}\overline{\boldsymbol{\omega}}-{}^{i}\mathbf{b}_g-{}^{i}\boldsymbol{\eta}_g \\
{}^{i}\mathbf{a}
&= {}^{i}\overline{\mathbf{a}}-{}^{i}\mathbf{b}_a-{}^{i}\boldsymbol{\eta}_a
\end{aligned}
$$

加速度计测得的是比力。代入运动学方程，世界系加速度为：

> 这里的比力可以理解为“单位质量受到的非重力力”，也就是除重力以外的支撑力、推力、接触力等带来的加速度。一个直观例子是：IMU 静止放在桌面上时，真实运动加速度为零，但加速度计仍会测到约 $9.8\,\mathrm{m/s^2}$ 的读数；这个读数来自桌面对 IMU 的支撑力，而不是 IMU 真的在向上加速。因此加速度计读数旋转到世界系后，还需要加回重力项，才能得到运动学方程中的世界系加速度。

$$
{}^{w}\mathbf{a}
= \mathbf{R}_{wi}\left({}^{i}\overline{\mathbf{a}}-{}^{i}\mathbf{b}_a-{}^{i}\boldsymbol{\eta}_a\right)+{}^{w}\mathbf{g}
$$

零偏通常建模为随机游走：

$$
\begin{aligned}
{}^{i}\mathbf{b}_{g,k+1}
&= {}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{bg,k}\Delta t \\
{}^{i}\mathbf{b}_{a,k+1}
&= {}^{i}\mathbf{b}_{a,k}+{}^{i}\boldsymbol{\eta}_{ba,k}\Delta t
\end{aligned}
$$

### 1.3 真值状态传播

令真值状态为：

$$
\mathbf{x}_k=
\left(
\mathbf{R}_{wi,k},
{}^{w}\mathbf{p}_k,
{}^{w}\mathbf{v}_k,
{}^{i}\mathbf{b}_{g,k},
{}^{i}\mathbf{b}_{a,k},
{}^{w}\mathbf{g}_k
\right)
$$

在 $\Delta t$ 内认为 IMU 测量、零偏和噪声为常值，则一阶离散传播为：

$$
\begin{aligned}
\mathbf{R}_{wi,k+1}
&= \mathbf{R}_{wi,k}
\operatorname{Exp}\!\left(
\left({}^{i}\overline{\boldsymbol{\omega}}_k-{}^{i}\mathbf{b}_{g,k}-{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
\right) \\
{}^{w}\mathbf{v}_{k+1}
&= {}^{w}\mathbf{v}_k+
\left[
\mathbf{R}_{wi,k}
\left({}^{i}\overline{\mathbf{a}}_k-{}^{i}\mathbf{b}_{a,k}-{}^{i}\boldsymbol{\eta}_{a,k}\right)
+{}^{w}\mathbf{g}_k
\right]\Delta t \\
{}^{w}\mathbf{p}_{k+1}
&= {}^{w}\mathbf{p}_k+{}^{w}\mathbf{v}_k\Delta t \\
{}^{i}\mathbf{b}_{g,k+1}
&= {}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{bg,k}\Delta t \\
{}^{i}\mathbf{b}_{a,k+1}
&= {}^{i}\mathbf{b}_{a,k}+{}^{i}\boldsymbol{\eta}_{ba,k}\Delta t \\
{}^{w}\mathbf{g}_{k+1}
&= {}^{w}\mathbf{g}_k
\end{aligned}
$$

上式的**位置项采用一阶欧拉积分**，便于后续误差传播保持简洁。工程实现中常用

$$
{}^{w}\mathbf{p}_{k+1}
= {}^{w}\mathbf{p}_k+{}^{w}\mathbf{v}_k\Delta t
+\frac{1}{2}
\left[
\mathbf{R}_{wi,k}
\left({}^{i}\overline{\mathbf{a}}_k-{}^{i}\mathbf{b}_{a,k}-{}^{i}\boldsymbol{\eta}_{a,k}\right)
+{}^{w}\mathbf{g}_k
\right]\Delta t^2
$$

这会让位置误差雅可比额外出现与姿态、加速度零偏、重力和加速度噪声相关的 $\frac{1}{2}\Delta t^2$ 项。

### 1.4 名义状态传播

滤波器实际维护的是名义状态。**名义状态传播使用同一组 IMU 测量，但不代入噪声**：

$$
\begin{aligned}
{}^{i}\hat{\boldsymbol{\omega}}_k
&= {}^{i}\overline{\boldsymbol{\omega}}_k-{}^{i}\hat{\mathbf{b}}_{g,k}^{+} \\
{}^{i}\hat{\mathbf{a}}_k
&= {}^{i}\overline{\mathbf{a}}_k-{}^{i}\hat{\mathbf{b}}_{a,k}^{+}
\end{aligned}
$$

$$
\begin{aligned}
\hat{\mathbf{R}}_{wi,k+1}^{-}
&= \hat{\mathbf{R}}_{wi,k}^{+}
\operatorname{Exp}\!\left({}^{i}\hat{\boldsymbol{\omega}}_k\Delta t\right) \\
&= \hat{\mathbf{R}}_{wi,k}^{+}
\operatorname{Exp}\!\left(
\left({}^{i}\overline{\boldsymbol{\omega}}_k-{}^{i}\hat{\mathbf{b}}_{g,k}^{+}\right)\Delta t
\right) \\
{}^{w}\hat{\mathbf{v}}_{k+1}^{-}
&= {}^{w}\hat{\mathbf{v}}_k^{+}
+\left(\hat{\mathbf{R}}_{wi,k}^{+}\,{}^{i}\hat{\mathbf{a}}_k+{}^{w}\hat{\mathbf{g}}_k^{+}\right)\Delta t \\
&= {}^{w}\hat{\mathbf{v}}_k^{+}
+\left[
\hat{\mathbf{R}}_{wi,k}^{+}
\left({}^{i}\overline{\mathbf{a}}_k-{}^{i}\hat{\mathbf{b}}_{a,k}^{+}\right)
+{}^{w}\hat{\mathbf{g}}_k^{+}
\right]\Delta t \\
{}^{w}\hat{\mathbf{p}}_{k+1}^{-}
&= {}^{w}\hat{\mathbf{p}}_k^{+}
+{}^{w}\hat{\mathbf{v}}_k^{+}\Delta t \\
{}^{i}\hat{\mathbf{b}}_{g,k+1}^{-}
&= {}^{i}\hat{\mathbf{b}}_{g,k}^{+} \\
{}^{i}\hat{\mathbf{b}}_{a,k+1}^{-}
&= {}^{i}\hat{\mathbf{b}}_{a,k}^{+} \\
{}^{w}\hat{\mathbf{g}}_{k+1}^{-}
&= {}^{w}\hat{\mathbf{g}}_k^{+}
\end{aligned}
$$

后续推导误差状态传播时，直接把真值传播式与这里展开后的名义传播式相减即可。这样可以清楚看到误差项来自哪里：姿态误差来自 $\mathbf{R}_{wi,k}$ 与 $\hat{\mathbf{R}}_{wi,k}^{+}$ 的差，速度误差还会额外来自加速度零偏误差、加速度噪声和重力误差。为了避免公式过于拥挤，后面的误差传播推导中省略名义状态在 $k$ 时刻的上标 $+$ 和在 $k+1$ 时刻的上标 $-$；也就是说，$\hat{\mathbf{x}}_k$ 默认表示 $\hat{\mathbf{x}}_k^{+}$，$\hat{\mathbf{x}}_{k+1}$ 默认表示 $\hat{\mathbf{x}}_{k+1}^{-}$。

若名义位置使用匀加速积分，只需把第三行替换为：

$$
\begin{aligned}
{}^{w}\hat{\mathbf{p}}_{k+1}^{-}
&= {}^{w}\hat{\mathbf{p}}_k^{+}
+{}^{w}\hat{\mathbf{v}}_k^{+}\Delta t
+\frac{1}{2}
\left(\hat{\mathbf{R}}_{wi,k}^{+}\,{}^{i}\hat{\mathbf{a}}_k+{}^{w}\hat{\mathbf{g}}_k^{+}\right)\Delta t^2 \\
&= {}^{w}\hat{\mathbf{p}}_k^{+}
+{}^{w}\hat{\mathbf{v}}_k^{+}\Delta t
+\frac{1}{2}
\left[
\hat{\mathbf{R}}_{wi,k}^{+}
\left({}^{i}\overline{\mathbf{a}}_k-{}^{i}\hat{\mathbf{b}}_{a,k}^{+}\right)
+{}^{w}\hat{\mathbf{g}}_k^{+}
\right]\Delta t^2
\end{aligned}
$$

## 2 ESKF 误差状态传播

### 2.1 误差状态定义

误差状态采用右扰动：

$$
\mathbf{x}_k=\hat{\mathbf{x}}_k\boxplus\tilde{\mathbf{x}}_k
$$

展开为：

$$
\tilde{\mathbf{x}}_k=
\begin{bmatrix}
\delta\boldsymbol{\theta}_k \\
\delta{}^{w}\mathbf{p}_k \\
\delta{}^{w}\mathbf{v}_k \\
\delta{}^{i}\mathbf{b}_{g,k} \\
\delta{}^{i}\mathbf{b}_{a,k} \\
\delta{}^{w}\mathbf{g}_k
\end{bmatrix}
=
\begin{bmatrix}
\operatorname{Log}\!\left((\hat{\mathbf{R}}_{wi,k})^\mathsf{T}\mathbf{R}_{wi,k}\right) \\
{}^{w}\mathbf{p}_k-{}^{w}\hat{\mathbf{p}}_k \\
{}^{w}\mathbf{v}_k-{}^{w}\hat{\mathbf{v}}_k \\
{}^{i}\mathbf{b}_{g,k}-{}^{i}\hat{\mathbf{b}}_{g,k} \\
{}^{i}\mathbf{b}_{a,k}-{}^{i}\hat{\mathbf{b}}_{a,k} \\
{}^{w}\mathbf{g}_k-{}^{w}\hat{\mathbf{g}}_k
\end{bmatrix}
\in\mathbb{R}^{18}
$$

姿态误差满足：

$$
\mathbf{R}_{wi,k}
= \hat{\mathbf{R}}_{wi,k}\operatorname{Exp}(\delta\boldsymbol{\theta}_k)
$$

其中 $\delta\boldsymbol{\theta}_k$ 是小量，表达在名义 IMU 坐标系附近。

### 2.2 姿态误差传播

记

$$
\boldsymbol{\Omega}_k
= {}^{i}\hat{\boldsymbol{\omega}}_k\Delta t
= \left({}^{i}\overline{\boldsymbol{\omega}}_k-{}^{i}\hat{\mathbf{b}}_{g,k}\right)\Delta t
$$

由真值姿态传播项和已展开的名义姿态传播项可得：

$$
\begin{aligned}
\delta\boldsymbol{\theta}_{k+1}
&=
\operatorname{Log}\!\left(
(\hat{\mathbf{R}}_{wi,k+1})^\mathsf{T}\mathbf{R}_{wi,k+1}
\right) \\
&=
\operatorname{Log}\!\left(
\left[
\hat{\mathbf{R}}_{wi,k}
\operatorname{Exp}(\boldsymbol{\Omega}_k)
\right]^\mathsf{T}
\mathbf{R}_{wi,k}
\operatorname{Exp}\!\left(
\left({}^{i}\overline{\boldsymbol{\omega}}_k-{}^{i}\mathbf{b}_{g,k}-{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
\right)
\right) \\
&=
\operatorname{Log}\!\left(
\operatorname{Exp}(-\boldsymbol{\Omega}_k)
\operatorname{Exp}(\delta\boldsymbol{\theta}_k)
\operatorname{Exp}\!\left(
\boldsymbol{\Omega}_k
-\left(\delta{}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
\right)
\right)
\end{aligned}
$$

> **从展开式到三段指数：**
>
> 在上式第二行中代入姿态误差和陀螺零偏误差：
>
> $$
> \begin{aligned}
> \mathbf{R}_{wi,k}
> &=
> \hat{\mathbf{R}}_{wi,k}\operatorname{Exp}(\delta\boldsymbol{\theta}_k) \\
> {}^{i}\mathbf{b}_{g,k}
> &=
> {}^{i}\hat{\mathbf{b}}_{g,k}
> +\delta{}^{i}\mathbf{b}_{g,k}
> \end{aligned}
> $$
>
> 真值姿态增量就可以写成
>
> $$
> \left({}^{i}\overline{\boldsymbol{\omega}}_k-{}^{i}\mathbf{b}_{g,k}-{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
> =
> \boldsymbol{\Omega}_k
> -\left(\delta{}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
> $$
>
> 因此
>
> $$
> \begin{aligned}
> &\left[
> \hat{\mathbf{R}}_{wi,k}
> \operatorname{Exp}(\boldsymbol{\Omega}_k)
> \right]^\mathsf{T}
> \mathbf{R}_{wi,k}
> \operatorname{Exp}\!\left(
> \left({}^{i}\overline{\boldsymbol{\omega}}_k-{}^{i}\mathbf{b}_{g,k}-{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
> \right) \\
> &=
> \operatorname{Exp}(-\boldsymbol{\Omega}_k)
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \operatorname{Exp}\!\left(
> \boldsymbol{\Omega}_k
> -\left(\delta{}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
> \right)
> \end{aligned}
> $$

对小量 $\delta\boldsymbol{\theta}_k$、$\delta{}^{i}\mathbf{b}_{g,k}$ 和 ${}^{i}\boldsymbol{\eta}_{g,k}$ 做一阶近似：

$$
\delta\boldsymbol{\theta}_{k+1}
\approx
\operatorname{Exp}(-\boldsymbol{\Omega}_k)\,\delta\boldsymbol{\theta}_k
-\mathbf{J}_r(\boldsymbol{\Omega}_k)
\left(\delta{}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
$$

> **一阶近似过程：**
>
> 令
>
> $$
> \boldsymbol{\epsilon}_k
> =
> \left(\delta{}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{g,k}\right)\Delta t
> $$
>
> 上一式中的最后一项为 $\operatorname{Exp}(\boldsymbol{\Omega}_k-\boldsymbol{\epsilon}_k)$。当 $\boldsymbol{\epsilon}_k$ 是小量时，有右扰动近似
>
> $$
> \operatorname{Exp}(\boldsymbol{\Omega}_k-\boldsymbol{\epsilon}_k)
> \approx
> \operatorname{Exp}(\boldsymbol{\Omega}_k)
> \operatorname{Exp}\!\left(-\mathbf{J}_r(\boldsymbol{\Omega}_k)\boldsymbol{\epsilon}_k\right)
> $$
>
> 因此
>
> $$
> \begin{aligned}
> &\operatorname{Exp}(-\boldsymbol{\Omega}_k)
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \operatorname{Exp}(\boldsymbol{\Omega}_k-\boldsymbol{\epsilon}_k) \\
> &\approx
> \operatorname{Exp}(-\boldsymbol{\Omega}_k)
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \operatorname{Exp}(\boldsymbol{\Omega}_k)
> \operatorname{Exp}\!\left(-\mathbf{J}_r(\boldsymbol{\Omega}_k)\boldsymbol{\epsilon}_k\right)
> \end{aligned}
> $$
>
> 对前三项使用伴随关系
>
> $$
> \mathbf{R}\operatorname{Exp}(\boldsymbol{\phi})\mathbf{R}^\mathsf{T}
> =
> \operatorname{Exp}(\mathbf{R}\boldsymbol{\phi})
> $$
>
> 取 $\mathbf{R}=\operatorname{Exp}(-\boldsymbol{\Omega}_k)$，得到
>
> $$
> \operatorname{Exp}(-\boldsymbol{\Omega}_k)
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \operatorname{Exp}(\boldsymbol{\Omega}_k)
> =
> \operatorname{Exp}\!\left(
> \operatorname{Exp}(-\boldsymbol{\Omega}_k)\delta\boldsymbol{\theta}_k
> \right)
> $$
>
> 再对两个小量的指数乘积使用 BCH 一阶近似
> $\operatorname{Log}(\operatorname{Exp}(\mathbf{a})\operatorname{Exp}(\mathbf{b}))\approx\mathbf{a}+\mathbf{b}$，得到
>
> $$
> \delta\boldsymbol{\theta}_{k+1}
> \approx
> \operatorname{Exp}(-\boldsymbol{\Omega}_k)\delta\boldsymbol{\theta}_k
> -\mathbf{J}_r(\boldsymbol{\Omega}_k)\boldsymbol{\epsilon}_k
> $$
>
> 将 $\boldsymbol{\epsilon}_k$ 展开，就得到正文中的姿态误差传播式。

这里 $\operatorname{Exp}(-\boldsymbol{\Omega}_k)$ 是旋转矩阵，对向量左乘表示把旧姿态误差搬运到新名义 IMU 系附近；$\mathbf{J}_r(\cdot)$ 是 SO(3) 右雅可比。

### 2.3 位置误差传播

在一阶欧拉位置传播下，直接用真值位置传播项减去已展开的名义位置传播项：

$$
\begin{aligned}
\delta{}^{w}\mathbf{p}_{k+1}
&= {}^{w}\mathbf{p}_{k+1}-{}^{w}\hat{\mathbf{p}}_{k+1} \\
&=
\left({}^{w}\mathbf{p}_k+{}^{w}\mathbf{v}_k\Delta t\right)
-\left({}^{w}\hat{\mathbf{p}}_k+{}^{w}\hat{\mathbf{v}}_k\Delta t\right) \\
&= \delta{}^{w}\mathbf{p}_k+\delta{}^{w}\mathbf{v}_k\Delta t
\end{aligned}
$$

若位置采用匀加速离散，这一行还应加入：
$$
\frac{1}{2}
\left[
-\hat{\mathbf{R}}_{wi,k}\left[{}^{i}\hat{\mathbf{a}}_k\right]_\times\delta\boldsymbol{\theta}_k
-\hat{\mathbf{R}}_{wi,k}\delta{}^{i}\mathbf{b}_{a,k}
+\delta{}^{w}\mathbf{g}_k
-\hat{\mathbf{R}}_{wi,k}{}^{i}\boldsymbol{\eta}_{a,k}
\right]\Delta t^2
$$

> **匀加速位置项说明：**
>
> 若位置也使用常加速度积分，则位置真值项和名义项分别包含
>
> $$
> \frac{1}{2}
> \left[
> \mathbf{R}_{wi,k}
> \left({}^{i}\overline{\mathbf{a}}_k-{}^{i}\mathbf{b}_{a,k}-{}^{i}\boldsymbol{\eta}_{a,k}\right)
> +{}^{w}\mathbf{g}_k
> \right]\Delta t^2
> $$
>
> 和
>
> $$
> \frac{1}{2}
> \left[
> \hat{\mathbf{R}}_{wi,k}
> \left({}^{i}\overline{\mathbf{a}}_k-{}^{i}\hat{\mathbf{b}}_{a,k}\right)
> +{}^{w}\hat{\mathbf{g}}_k
> \right]\Delta t^2
> $$
>
> 两者相减后，括号里的加速度误差项与**速度误差传播**中的加速度误差项相同，因此得到正文中的 $\frac{1}{2}\Delta t^2$ 补偿项。

### 2.4 速度误差传播

由速度真值传播减去名义传播：

$$
\begin{aligned}
\delta{}^{w}\mathbf{v}_{k+1}
&= \delta{}^{w}\mathbf{v}_k
+\left[
\mathbf{R}_{wi,k}
\left({}^{i}\overline{\mathbf{a}}_k-{}^{i}\mathbf{b}_{a,k}-{}^{i}\boldsymbol{\eta}_{a,k}\right)
-\hat{\mathbf{R}}_{wi,k}\,\left({}^{i}\overline{\mathbf{a}}_k
-{}^{i}\hat{\mathbf{b}}_{a,k}\right)
+\delta{}^{w}\mathbf{g}_k
\right]\Delta t
\end{aligned}
$$

将姿态误差和加速度计零偏误差代入：

$$
\mathbf{R}_{wi,k}
=\hat{\mathbf{R}}_{wi,k}\operatorname{Exp}(\delta\boldsymbol{\theta}_k),
\qquad
{}^{i}\mathbf{b}_{a,k}
= {}^{i}\hat{\mathbf{b}}_{a,k}+\delta{}^{i}\mathbf{b}_{a,k}
$$

并忽略二阶小量，得到：

$$
\delta{}^{w}\mathbf{v}_{k+1}
\approx
\delta{}^{w}\mathbf{v}_k
-\hat{\mathbf{R}}_{wi,k}\left[{}^{i}\hat{\mathbf{a}}_k\right]_\times
\delta\boldsymbol{\theta}_k\Delta t
-\hat{\mathbf{R}}_{wi,k}
\delta{}^{i}\mathbf{b}_{a,k}\Delta t
-\hat{\mathbf{R}}_{wi,k}{}^{i}\boldsymbol{\eta}_{a,k}\Delta t
+\delta{}^{w}\mathbf{g}_k\Delta t
$$

> **速度误差的一阶展开：**
>
> 从正文已经代入名义状态传播式后的方括号项开始，记
>
> $$
> \mathbf{s}_k
> =
> \mathbf{R}_{wi,k}
> \left(
> {}^{i}\overline{\mathbf{a}}_k
> -{}^{i}\mathbf{b}_{a,k}
> -{}^{i}\boldsymbol{\eta}_{a,k}
> \right)
> -\hat{\mathbf{R}}_{wi,k}
> \left(
> {}^{i}\overline{\mathbf{a}}_k
> -{}^{i}\hat{\mathbf{b}}_{a,k}
> \right)
> $$
>
> 再代入
>
> $$
> \mathbf{R}_{wi,k}
> =
> \hat{\mathbf{R}}_{wi,k}
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k),
> \qquad
> {}^{i}\mathbf{b}_{a,k}
> =
> {}^{i}\hat{\mathbf{b}}_{a,k}
> +\delta{}^{i}\mathbf{b}_{a,k}
> $$
>
> 得到
>
> $$
> \begin{aligned}
> \mathbf{s}_k
> &=
> \hat{\mathbf{R}}_{wi,k}
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \left(
> {}^{i}\overline{\mathbf{a}}_k
> -{}^{i}\hat{\mathbf{b}}_{a,k}
> -\delta{}^{i}\mathbf{b}_{a,k}
> -{}^{i}\boldsymbol{\eta}_{a,k}
> \right)
> -\hat{\mathbf{R}}_{wi,k}
> \left(
> {}^{i}\overline{\mathbf{a}}_k
> -{}^{i}\hat{\mathbf{b}}_{a,k}
> \right) \\
> &=
> \hat{\mathbf{R}}_{wi,k}
> \left[
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \left(
> {}^{i}\overline{\mathbf{a}}_k
> -{}^{i}\hat{\mathbf{b}}_{a,k}
> \right)
> -\left(
> {}^{i}\overline{\mathbf{a}}_k
> -{}^{i}\hat{\mathbf{b}}_{a,k}
> \right)
> \right] \\
> &\quad
> -\hat{\mathbf{R}}_{wi,k}
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \left(
> \delta{}^{i}\mathbf{b}_{a,k}
> +{}^{i}\boldsymbol{\eta}_{a,k}
> \right)
> \end{aligned}
> $$
>
> 为了简化后续书写，再把名义去零偏加速度记为
>
> $$
> {}^{i}\hat{\mathbf{a}}_k
> =
> {}^{i}\overline{\mathbf{a}}_k
> -{}^{i}\hat{\mathbf{b}}_{a,k}
> $$
>
> 对姿态小量做一阶近似：
>
> $$
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \approx
> \mathbf{I}+[\delta\boldsymbol{\theta}_k]_\times
> $$
>
> 因此
>
> $$
> \begin{aligned}
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k){}^{i}\hat{\mathbf{a}}_k
> &\approx
> {}^{i}\hat{\mathbf{a}}_k
> +\delta\boldsymbol{\theta}_k\times{}^{i}\hat{\mathbf{a}}_k \\
> &=
> {}^{i}\hat{\mathbf{a}}_k
> -\left[{}^{i}\hat{\mathbf{a}}_k\right]_\times
> \delta\boldsymbol{\theta}_k
> \end{aligned}
> $$
>
> 同时，
>
> $$
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k)
> \left(
> \delta{}^{i}\mathbf{b}_{a,k}
> +{}^{i}\boldsymbol{\eta}_{a,k}
> \right)
> \approx
> \delta{}^{i}\mathbf{b}_{a,k}
> +{}^{i}\boldsymbol{\eta}_{a,k}
> $$
>
> 因为展开后多出来的
>
> $$
> [\delta\boldsymbol{\theta}_k]_\times
> \left(
> \delta{}^{i}\mathbf{b}_{a,k}
> +{}^{i}\boldsymbol{\eta}_{a,k}
> \right)
> $$
>
> 是姿态误差和零偏/噪声误差的乘积，属于二阶小量，在一阶误差传播中忽略。
>
> 综上：
>
> $$
> \begin{aligned}
> \mathbf{s}_k
> &\approx
> -\hat{\mathbf{R}}_{wi,k}
> \left[{}^{i}\hat{\mathbf{a}}_k\right]_\times
> \delta\boldsymbol{\theta}_k
> -\hat{\mathbf{R}}_{wi,k}\delta{}^{i}\mathbf{b}_{a,k}
> -\hat{\mathbf{R}}_{wi,k}{}^{i}\boldsymbol{\eta}_{a,k}
> \end{aligned}
> $$
>
> 再加上 $\delta{}^{w}\mathbf{v}_k$ 和 $\delta{}^{w}\mathbf{g}_k\Delta t$，就得到正文中的速度误差传播式。



### 2.5 零偏和重力误差传播

零偏随机游走和重力常值模型同样按“真值项 - 已展开名义项”得到：

$$
\begin{aligned}
\delta{}^{i}\mathbf{b}_{g,k+1}
&=
\left({}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{bg,k}\Delta t\right)
-{}^{i}\hat{\mathbf{b}}_{g,k+1} \\
&=
\left({}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{bg,k}\Delta t\right)
-{}^{i}\hat{\mathbf{b}}_{g,k} \\
&= \delta{}^{i}\mathbf{b}_{g,k}+{}^{i}\boldsymbol{\eta}_{bg,k}\Delta t \\
\delta{}^{i}\mathbf{b}_{a,k+1}
&=
\left({}^{i}\mathbf{b}_{a,k}+{}^{i}\boldsymbol{\eta}_{ba,k}\Delta t\right)
-{}^{i}\hat{\mathbf{b}}_{a,k+1} \\
&=
\left({}^{i}\mathbf{b}_{a,k}+{}^{i}\boldsymbol{\eta}_{ba,k}\Delta t\right)
-{}^{i}\hat{\mathbf{b}}_{a,k} \\
&= \delta{}^{i}\mathbf{b}_{a,k}+{}^{i}\boldsymbol{\eta}_{ba,k}\Delta t \\
\delta{}^{w}\mathbf{g}_{k+1}
&=
{}^{w}\mathbf{g}_{k+1}-{}^{w}\hat{\mathbf{g}}_{k+1} \\
&=
{}^{w}\mathbf{g}_k-{}^{w}\hat{\mathbf{g}}_k \\
&= \delta{}^{w}\mathbf{g}_k
\end{aligned}
$$

> **零偏和重力误差代入说明：**
>
> 名义零偏传播不代入随机游走噪声：
>
> $$
> {}^{i}\hat{\mathbf{b}}_{g,k+1}
> =
> {}^{i}\hat{\mathbf{b}}_{g,k},
> \qquad
> {}^{i}\hat{\mathbf{b}}_{a,k+1}
> =
> {}^{i}\hat{\mathbf{b}}_{a,k}
> $$
>
> 真值零偏传播包含随机游走噪声：
>
> $$
> {}^{i}\mathbf{b}_{g,k+1}
> =
> {}^{i}\mathbf{b}_{g,k}
> +{}^{i}\boldsymbol{\eta}_{bg,k}\Delta t,
> \qquad
> {}^{i}\mathbf{b}_{a,k+1}
> =
> {}^{i}\mathbf{b}_{a,k}
> +{}^{i}\boldsymbol{\eta}_{ba,k}\Delta t
> $$
>
> 两者相减后，零偏误差会继承真值随机游走噪声。重力在真值和名义模型中都按常值传播，因此重力误差保持不变。

### 2.6 误差状态矩阵形式

把误差传播写成线性形式：

$$
\tilde{\mathbf{x}}_{k+1}^{-}
=
\mathbf{F}_{\tilde{x},k}\tilde{\mathbf{x}}_k^{+}
+\mathbf{F}_{\eta,k}\boldsymbol{\eta}_k
$$

其中

$$
\boldsymbol{\eta}_k=
\begin{bmatrix}
{}^{i}\boldsymbol{\eta}_{g,k} \\
{}^{i}\boldsymbol{\eta}_{a,k} \\
{}^{i}\boldsymbol{\eta}_{bg,k} \\
{}^{i}\boldsymbol{\eta}_{ba,k}
\end{bmatrix}
$$

一阶欧拉位置传播对应的状态雅可比为：

$$
\mathbf{F}_{\tilde{x},k}
=
\begin{bmatrix}
\operatorname{Exp}(-\boldsymbol{\Omega}_k) & \mathbf{0} & \mathbf{0} & -\mathbf{J}_r(\boldsymbol{\Omega}_k)\Delta t & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{I} & \mathbf{I}\Delta t & \mathbf{0} & \mathbf{0} & \mathbf{0} \\
-\hat{\mathbf{R}}_{wi,k}\left[{}^{i}\hat{\mathbf{a}}_k\right]_\times\Delta t & \mathbf{0} & \mathbf{I} & \mathbf{0} & -\hat{\mathbf{R}}_{wi,k}\Delta t & \mathbf{I}\Delta t \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{I} & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{I} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{I}
\end{bmatrix}
$$

其中$\boldsymbol{\Omega}_k
= {}^{i}\hat{\boldsymbol{\omega}}_k\Delta t$，噪声雅可比为：
$$
\mathbf{F}_{\eta,k}
=
\begin{bmatrix}
-\mathbf{J}_r(\boldsymbol{\Omega}_k)\Delta t & \mathbf{0} & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & -\hat{\mathbf{R}}_{wi,k}\Delta t & \mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{I}\Delta t & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{I}\Delta t \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & \mathbf{0}
\end{bmatrix}
$$

若位置采用匀加速离散，应在 $\mathbf{F}_{\tilde{x},k}$ 的第二行加入：

$$
\begin{aligned}
\mathbf{F}_{p\theta} &= -\frac{1}{2}\hat{\mathbf{R}}_{wi,k}\left[{}^{i}\hat{\mathbf{a}}_k\right]_\times\Delta t^2 \\
\mathbf{F}_{p b_a} &= -\frac{1}{2}\hat{\mathbf{R}}_{wi,k}\Delta t^2 \\
\mathbf{F}_{p g} &= \frac{1}{2}\mathbf{I}\Delta t^2
\end{aligned}
$$

并在 $\mathbf{F}_{\eta,k}$ 的位置行加入：

$$
\mathbf{F}_{p\eta_a}=-\frac{1}{2}\hat{\mathbf{R}}_{wi,k}\Delta t^2
$$

协方差预测为：

$$
\tilde{\mathbf{P}}_{k+1}^{-}
=
\mathbf{F}_{\tilde{x},k}
\tilde{\mathbf{P}}_{k}^{+}
\mathbf{F}_{\tilde{x},k}^{\mathsf{T}}
+\mathbf{F}_{\eta,k}
\mathbf{W}_k
\mathbf{F}_{\eta,k}^{\mathsf{T}}
$$

> **为什么这里传播的是误差状态协方差？**
>
> 上面的协方差公式来自前面得到的误差状态线性传播：
>
> $$
> \tilde{\mathbf{x}}_{k+1}^{-}
> =
> \mathbf{F}_{\tilde{x},k}\tilde{\mathbf{x}}_k^{+}
> +
> \mathbf{F}_{\eta,k}\boldsymbol{\eta}_k.
> $$
>
> 也就是说，$\mathbf{F}_{\tilde{x},k}$ 是误差状态对误差状态的雅可比，$\mathbf{F}_{\eta,k}$ 是误差状态对 IMU 噪声的雅可比。因此这里传播的是误差状态协方差 $\tilde{\mathbf{P}}$，不是名义状态 $\hat{\mathbf{x}}$ 的协方差。
>
> 这和普通 KF 的协方差传播并不矛盾。KF 中的
>
> $$
> \mathbf{P}_{k}^{-}
> =
> \mathbf{F}_{k}\mathbf{P}_{k-1}^{+}\mathbf{F}_{k}^{\mathsf{T}}
> +\mathbf{W}_{k}
> $$
>
> 本质上也是从误差传播
>
> $$
> \tilde{\mathbf{x}}_{k}^{-}
> =
> \mathbf{F}_{k}\tilde{\mathbf{x}}_{k-1}^{+}
> +\mathbf{w}_{k}
> $$
>
> 推出来的。只是线性 KF 的状态空间本身就是欧式向量空间，状态和状态误差在同一个线性空间里，所以常简称为“状态协方差”。ESKF 中姿态等变量不适合直接做普通加减，因此必须明确维护和传播切空间中的误差状态协方差。
>
> 这里 $\mathbf{W}_k$ 是 $\boldsymbol{\eta}_k$ 的协方差。**ESKF 在每次状态注入后会把误差状态均值重置为零，因此预测阶段最重要的是协方差如何被过程模型和 IMU 噪声传播。**

## 3 普通观测更新

设传感器观测模型为：

$$
\mathbf{z}_k
=
\mathbf{h}(\mathbf{x}_k)+\mathbf{n}_k,
\qquad
\mathbf{n}_k\sim\mathcal{N}(\mathbf{0},\mathbf{N}_k)
$$

把真值状态写成名义状态加误差状态：

$$
\mathbf{x}_k
=
\hat{\mathbf{x}}_k^{-}\boxplus\tilde{\mathbf{x}}_k^{-}
$$

在 $\tilde{\mathbf{x}}_k^{-}=\mathbf{0}$ 处线性化：

$$
\begin{aligned}
\mathbf{z}_k-\mathbf{h}(\hat{\mathbf{x}}_k^{-})
&\approx
\mathbf{H}_k\tilde{\mathbf{x}}_k^{-}+\mathbf{n}_k \\
\mathbf{H}_k
&=
\left.
\frac{\partial \mathbf{h}(\hat{\mathbf{x}}_k^{-}\boxplus\tilde{\mathbf{x}})}
{\partial\tilde{\mathbf{x}}}
\right|_{\tilde{\mathbf{x}}=\mathbf{0}}
\end{aligned}
$$

这一定义比直接写 $\frac{\partial\mathbf{h}}{\partial\mathbf{x}}\frac{\partial\mathbf{x}}{\partial\tilde{\mathbf{x}}}$ 更稳妥，因为状态中包含旋转，普通向量空间的导数容易混淆左扰动和右扰动。若某个观测只依赖位置、速度、零偏或重力，相关块就是普通向量导数；若观测依赖姿态，则必须按上面的 $\boxplus$ 定义对 $\delta\boldsymbol{\theta}$ 求导。

卡尔曼更新为：

$$
\begin{aligned}
\mathbf{K}_k
&=
\tilde{\mathbf{P}}_k^{-}\mathbf{H}_k^\mathsf{T}
\left(
\mathbf{H}_k\tilde{\mathbf{P}}_k^{-}\mathbf{H}_k^\mathsf{T}
+\mathbf{N}_k
\right)^{-1} \\
\tilde{\mathbf{x}}_k^{+}
&=
\mathbf{K}_k
\left(
\mathbf{z}_k-\mathbf{h}(\hat{\mathbf{x}}_k^{-})
\right) \\
\tilde{\mathbf{P}}_k^{+}
&=
\left(\mathbf{I}-\mathbf{K}_k\mathbf{H}_k\right)
\tilde{\mathbf{P}}_k^{-}
\left(\mathbf{I}-\mathbf{K}_k\mathbf{H}_k\right)^\mathsf{T}
+\mathbf{K}_k\mathbf{N}_k\mathbf{K}_k^\mathsf{T}
\end{aligned}
$$

最后一行采用 Joseph 形式，数值上比简单写成 $(\mathbf{I}-\mathbf{K}_k\mathbf{H}_k)\tilde{\mathbf{P}}_k^{-}$ 更不容易破坏协方差的对称半正定性。

## 4 迭代观测更新

上面的观测更新只写了一次线性化。实际 LiDAR 紧耦合系统中，ICP、NDT 等点云配准观测通常是非线性的，数据关联也可能随位姿变化。若只在预测位姿 $\hat{\mathbf{x}}_k^{-}$ 处线性化一次，残差较大时容易得到不稳定的更新。因此工程上通常在同一个时刻 $k$ 内做若干次内循环线性化，写成**迭代 ESKF** 或 **iterated ESKF**。

需要注意：迭代观测更新并不是每次都重新做 IMU 预测。IMU 预测给出的先验

$$
\tilde{\mathbf{x}}_k^{-}\sim
\mathcal{N}\!\left(\mathbf{0},\tilde{\mathbf{P}}_k^{-}\right)
$$

在整个观测内循环中保持不变。变化的是观测线性化点 $\hat{\mathbf{x}}_k^\kappa$，以及先验分布被表达在哪一个误差坐标系中。

### 4.1 先验映射与误差状态观测方程

初始化时令

$$
\hat{\mathbf{x}}_k^0=\hat{\mathbf{x}}_k^{-}.
$$

第 $\kappa$ 次迭代时，真值状态仍然写成当前线性化点加误差状态：

$$
\tilde{\mathbf{x}}_k^\kappa
=
\mathbf{x}_k \boxminus \hat{\mathbf{x}}_k^\kappa,
\qquad
\mathbf{x}_k
=
\hat{\mathbf{x}}_k^\kappa
\boxplus
\tilde{\mathbf{x}}_k^\kappa.
$$

这里 $\tilde{\mathbf{x}}_k^\kappa$ 是**相对于第 $\kappa$ 次迭代点的局部误差状态**。特别地，因为 $\hat{\mathbf{x}}_k^0=\hat{\mathbf{x}}_k^{-}$，所以 $\tilde{\mathbf{x}}_k^0$ 就是相对于预测名义状态的误差状态，和前面写的先验误差 $\tilde{\mathbf{x}}_k^{-}$ 表达在同一个误差坐标系中。

但是先验项需要衡量真值状态相对预测名义状态 $\hat{\mathbf{x}}_k^{-}$ 的误差。为避免和局部误差 $\tilde{\mathbf{x}}_k^\kappa$ 混淆，把当前迭代点相对预测点的累计偏移记为

$$
\boldsymbol{\xi}_k^\kappa
=
\hat{\mathbf{x}}_k^\kappa \boxminus \hat{\mathbf{x}}_k^{-}
=
\begin{bmatrix}
\boldsymbol{\xi}_{\theta,k}^\kappa \\
\boldsymbol{\xi}_{p,k}^\kappa \\
\boldsymbol{\xi}_{v,k}^\kappa \\
\boldsymbol{\xi}_{bg,k}^\kappa \\
\boldsymbol{\xi}_{ba,k}^\kappa \\
\boldsymbol{\xi}_{g,k}^\kappa
\end{bmatrix},
\qquad
\boldsymbol{\xi}_{\theta,k}^\kappa
=
\operatorname{Log}\!\left(
\left(\hat{\mathbf{R}}_{wi,k}^{-}\right)^\mathsf{T}
\hat{\mathbf{R}}_{wi,k}^{\kappa}
\right),
$$

其余向量块就是当前迭代名义值减去预测名义值。这里的 $\boxminus$ 是按状态各分量定义的广义减法：旋转用 $\operatorname{Log}$，普通向量用直接相减。

由于先验协方差 $\tilde{\mathbf{P}}_k^{-}$ 表达的是预测点误差坐标中的不确定性，要把当前局部误差 $\tilde{\mathbf{x}}_k^\kappa$ 映射回预测点的误差坐标：

$$
\tilde{\mathbf{x}}_k^0
=
\mathbf{x}_k \boxminus \hat{\mathbf{x}}_k^{-}
\approx
\boldsymbol{\xi}_k^\kappa
+
\mathbf{J}_\kappa
\tilde{\mathbf{x}}_k^\kappa.
$$

对右扰动姿态误差，有

$$
\operatorname{Log}\!\left(
\operatorname{Exp}(\boldsymbol{\xi}_{\theta,k}^\kappa)
\operatorname{Exp}(\delta\boldsymbol{\theta}_k^\kappa)
\right)
\approx
\boldsymbol{\xi}_{\theta,k}^\kappa
+
\mathbf{J}_r^{-1}(\boldsymbol{\xi}_{\theta,k}^\kappa)
\delta\boldsymbol{\theta}_k^\kappa,
$$

所以

$$
\mathbf{J}_\kappa
=
\operatorname{diag}\!\left(
\mathbf{J}_r^{-1}(\boldsymbol{\xi}_{\theta,k}^\kappa),
\mathbf{I},
\mathbf{I},
\mathbf{I},
\mathbf{I},
\mathbf{I}
\right).
$$

当 $\kappa=0$ 时，$\boldsymbol{\xi}_k^0=\mathbf{0}$ 且 $\mathbf{J}_0=\mathbf{I}$，因此 $\tilde{\mathbf{x}}_k^0$ 就是预测点误差坐标中的先验误差。

有了先验在当前局部误差坐标中的表达，接下来构造同一坐标系下的观测方程。

为了统一写法，把观测模型改写成“残差期望为零”的形式：

$$
\mathbf{r}_k(\mathbf{x}_k)+\mathbf{n}_k=\mathbf{0},
\qquad
\mathbf{n}_k\sim\mathcal{N}(\mathbf{0},\mathbf{N}_k).
$$

对于标准观测 $\mathbf{z}_k=\mathbf{h}(\mathbf{x}_k)+\mathbf{n}_k$，可以取 $\mathbf{r}_k(\mathbf{x}_k)=\mathbf{h}(\mathbf{x}_k)-\mathbf{z}_k$；若取 $\mathbf{z}_k-\mathbf{h}(\mathbf{x}_k)$，只会改变残差和雅可比的符号，最终更新等价。点云配准更自然地直接定义几何残差，例如点到面距离、点到线距离或 NDT 的白化残差。

在当前迭代点 $\hat{\mathbf{x}}_k^\kappa$ 处线性化：

$$
\begin{aligned}
\mathbf{r}_k\!\left(
\hat{\mathbf{x}}_k^\kappa
\boxplus
\tilde{\mathbf{x}}_k^\kappa
\right)
&\approx
\mathbf{r}_k^\kappa
+
\mathbf{H}_k^\kappa
\tilde{\mathbf{x}}_k^\kappa \\
\mathbf{r}_k^\kappa
&=
\mathbf{r}_k(\hat{\mathbf{x}}_k^\kappa) \\
\mathbf{H}_k^\kappa
&=
\left.
\frac{\partial
\mathbf{r}_k(\hat{\mathbf{x}}_k^\kappa\boxplus\boldsymbol{\epsilon})}
{\partial\boldsymbol{\epsilon}}
\right|_{\boldsymbol{\epsilon}=\mathbf{0}}
\end{aligned}
$$

因此一次迭代中的线性观测方程是

$$
\mathbf{0}
\approx
\mathbf{r}_k^\kappa
+
\mathbf{H}_k^\kappa
\tilde{\mathbf{x}}_k^\kappa
+
\mathbf{n}_k.
$$

这里的未知量是当前迭代点附近的局部误差 $\tilde{\mathbf{x}}_k^\kappa$，不是预测点处的误差 $\tilde{\mathbf{x}}_k^0$。先验通过上面的 $\tilde{\mathbf{x}}_k^0\approx\boldsymbol{\xi}_k^\kappa+\mathbf{J}_\kappa\tilde{\mathbf{x}}_k^\kappa$ 进入方程。

### 4.2 点云残差、堆叠与单次迭代更新

上面的观测方程是抽象形式。对 LiDAR 点云配准来说，$\mathbf{r}_k^\kappa$ 和 $\mathbf{H}_k^\kappa$ 来自每个有效匹配点的几何残差与雅可比。

设 LiDAR 外参为 $\mathbf{T}_{il}=(\mathbf{R}_{il},{}^{i}\mathbf{t}_{il})$，即从 LiDAR 坐标系 $\mathcal{F}_l$ 到 IMU 坐标系 $\mathcal{F}_i$。LiDAR 点 $j$ 在 $\mathcal{F}_l$ 下的坐标完整写作 ${}^{l}\mathbf{p}_{lj}$，右下标 $lj$ 表示从 $o_l$ 指向点 $j$。经过外参变换后，该点在 IMU 坐标系下的坐标记为 ${}^{i}\mathbf{p}_{ij}$：

$$
{}^{i}\mathbf{p}_{ij}
=
\mathbf{R}_{il}\,{}^{l}\mathbf{p}_{lj}
+
{}^{i}\mathbf{t}_{il}.
$$

在第 $\kappa$ 次迭代名义位姿下，该点的世界坐标为

$$
{}^{w}\hat{\mathbf{p}}_{wj}^{\kappa}
=
\hat{\mathbf{R}}_{wi,k}^{\kappa}
{}^{i}\mathbf{p}_{ij}
+
{}^{w}\hat{\mathbf{p}}_k^\kappa.
$$

对当前迭代点加入局部误差

$$
\tilde{\mathbf{x}}_k^\kappa
=
\begin{bmatrix}
\delta\boldsymbol{\theta}_k^\kappa \\
\delta{}^{w}\mathbf{p}_k^\kappa \\
\delta{}^{w}\mathbf{v}_k^\kappa \\
\delta{}^{i}\mathbf{b}_{g,k}^\kappa \\
\delta{}^{i}\mathbf{b}_{a,k}^\kappa \\
\delta{}^{w}\mathbf{g}_k^\kappa
\end{bmatrix},
$$

则

$$
\begin{aligned}
{}^{w}\mathbf{p}_{wj}
\!\left(\tilde{\mathbf{x}}_k^\kappa\right)
&=
\hat{\mathbf{R}}_{wi,k}^{\kappa}
\operatorname{Exp}(\delta\boldsymbol{\theta}_k^\kappa)
{}^{i}\mathbf{p}_{ij}
+
{}^{w}\hat{\mathbf{p}}_k^\kappa
+
\delta{}^{w}\mathbf{p}_k^\kappa \\
&\approx
{}^{w}\hat{\mathbf{p}}_{wj}^{\kappa}
-
\hat{\mathbf{R}}_{wi,k}^{\kappa}
\left[{}^{i}\mathbf{p}_{ij}\right]_\times
\delta\boldsymbol{\theta}_k^\kappa
+
\delta{}^{w}\mathbf{p}_k^\kappa .
\end{aligned}
$$

不同配准方法只是选择了不同的残差函数：

- 点到点 ICP：残差通常是当前点到匹配点的 $3$ 维差；
- 点到线 ICP：残差通常取垂直于匹配线方向的 $2$ 维分量；
- 点到面 ICP：残差通常是 $1$ 维点面距离；
- NDT：残差通常来自点落入对应高斯单元后的白化距离。

以点到面 ICP 为例，若局部平面法向为 ${}^{w}\mathbf{n}_j$，平面上一点为 $q_j$，其世界坐标完整写作 ${}^{w}\mathbf{q}_{wj}$，则

$$
r_j^\kappa
=
{}^{w}\mathbf{n}_j^{\mathsf{T}}
\left(
{}^{w}\hat{\mathbf{p}}_{wj}^{\kappa}
-
{}^{w}\mathbf{q}_{wj}
\right).
$$

对应雅可比行是

$$
\mathbf{H}_j^\kappa
=
\begin{bmatrix}
-{}^{w}\mathbf{n}_j^{\mathsf{T}}
\hat{\mathbf{R}}_{wi,k}^{\kappa}
\left[{}^{i}\mathbf{p}_{ij}\right]_\times
&
{}^{w}\mathbf{n}_j^{\mathsf{T}}
&
\mathbf{0}
&
\mathbf{0}
&
\mathbf{0}
&
\mathbf{0}
\end{bmatrix}.
$$

列块顺序与本文误差状态

$$
\tilde{\mathbf{x}}_k=
\begin{bmatrix}
\delta\boldsymbol{\theta}_k &
\delta{}^{w}\mathbf{p}_k &
\delta{}^{w}\mathbf{v}_k &
\delta{}^{i}\mathbf{b}_{g,k} &
\delta{}^{i}\mathbf{b}_{a,k} &
\delta{}^{w}\mathbf{g}_k
\end{bmatrix}^{\mathsf{T}}
$$

一致。点云残差只直接依赖当前位姿，所以速度、零偏、重力列为零；这些状态仍会通过先验协方差中的相关项被间接修正。

上面得到的 $r_j^\kappa$ 和 $\mathbf{H}_j^\kappa$ 只是第 $j$ 个有效匹配点在当前迭代点处的线性化结果。第 $\kappa$ 次迭代中，需要先用 $\hat{\mathbf{x}}_k^\kappa$ 把当前扫描投到地图中，完成数据关联、几何一致性筛选和鲁棒权重计算；保留下来的有效约束再统一进入观测更新。严格来说，有效约束数会随迭代点变化，可写成 $J_k^\kappa$，下面为简洁仍记为 $J_k$。

对第 $j$ 个有效约束，统一写成

$$
\mathbf{0}
\approx
\mathbf{r}_j^\kappa
+
\mathbf{H}_j^\kappa
\tilde{\mathbf{x}}_k^\kappa
+
\mathbf{n}_j,
\qquad
\mathbf{n}_j\sim\mathcal{N}(\mathbf{0},\mathbf{N}_j),
$$

其中

$$
\mathbf{r}_j^\kappa\in\mathbb{R}^{d_j},
\qquad
\mathbf{H}_j^\kappa\in\mathbb{R}^{d_j\times n_x},
\qquad
\mathbf{N}_j\in\mathbb{R}^{d_j\times d_j}.
$$

这里 $n_x$ 是误差状态维度，$d_j$ 是单个约束的残差维度。点到面残差 $d_j=1$，点到线残差通常 $d_j=2$，点到点残差通常 $d_j=3$。若使用鲁棒核或匹配质量权重，通常在堆叠前把同一个白化/加权因子乘到 $\mathbf{r}_j^\kappa$ 和 $\mathbf{H}_j^\kappa$ 上，等价于调整该约束对应的 $\mathbf{N}_j$。

把所有有效约束按行堆叠，得到整帧观测残差、雅可比和噪声：

$$
\mathbf{r}_k^\kappa
=
\begin{bmatrix}
\mathbf{r}_1^\kappa\\
\mathbf{r}_2^\kappa\\
\vdots\\
\mathbf{r}_{J_k}^\kappa
\end{bmatrix},
\qquad
\mathbf{H}_k^\kappa
=
\begin{bmatrix}
\mathbf{H}_1^\kappa\\
\mathbf{H}_2^\kappa\\
\vdots\\
\mathbf{H}_{J_k}^\kappa
\end{bmatrix},
\qquad
\mathbf{n}_k
=
\begin{bmatrix}
\mathbf{n}_1\\
\mathbf{n}_2\\
\vdots\\
\mathbf{n}_{J_k}
\end{bmatrix}.
$$

若各点残差噪声相互独立，则

$$
\mathbf{N}_k
=
\operatorname{blkdiag}\!\left(
\mathbf{N}_1,
\mathbf{N}_2,
\ldots,
\mathbf{N}_{J_k}
\right).
$$

因此整帧线性观测方程就是

$$
\mathbf{0}
\approx
\mathbf{r}_k^\kappa
+
\mathbf{H}_k^\kappa
\tilde{\mathbf{x}}_k^\kappa
+
\mathbf{n}_k,
\qquad
\mathbf{n}_k\sim\mathcal{N}(\mathbf{0},\mathbf{N}_k).
$$

所以 $\mathbf{H}_k^\kappa$ 不是由单个点直接得到的，而是所有有效点约束的雅可比按行堆叠后的结果；$\mathbf{r}_k^\kappa$ 和 $\mathbf{N}_k$ 同理。

第 $\kappa$ 次迭代时，前面已经把预测点误差坐标和当前迭代点误差坐标联系起来：

$$
\tilde{\mathbf{x}}_k^0
\approx
\boldsymbol{\xi}_k^\kappa
+
\mathbf{J}_\kappa
\tilde{\mathbf{x}}_k^\kappa
$$

由于预测先验为 $\tilde{\mathbf{x}}_k^0\sim\mathcal{N}(\mathbf{0},\tilde{\mathbf{P}}_k^{-})$，当前迭代点附近的局部误差状态服从近似先验分布：

$$
\tilde{\mathbf{x}}_k^\kappa
\sim
\mathcal{N}\!\left(
\tilde{\boldsymbol{\mu}}_k^{\kappa-},
\tilde{\mathbf{P}}_k^{\kappa-}
\right),
\qquad
\tilde{\boldsymbol{\mu}}_k^{\kappa-}
=
-\mathbf{J}_\kappa^{-1}\boldsymbol{\xi}_k^\kappa,
\qquad
\tilde{\mathbf{P}}_k^{\kappa-}
=
\mathbf{J}_\kappa^{-1}
\tilde{\mathbf{P}}_k^{-}
\mathbf{J}_\kappa^{-\mathsf{T}}.
$$

用上面堆叠得到的整帧残差和雅可比，当前迭代点处的线性观测方程为

$$
\mathbf{0}
\approx
\mathbf{r}_k^\kappa
+
\mathbf{H}_k^\kappa
\tilde{\mathbf{x}}_k^\kappa
+
\mathbf{n}_k,
\qquad
\mathbf{n}_k\sim\mathcal{N}(\mathbf{0},\mathbf{N}_k).
$$

把它改写成标准线性观测形式：

$$
\mathbf{z}_k^\kappa
\approx
\mathbf{H}_k^\kappa
\tilde{\mathbf{x}}_k^\kappa
+
\mathbf{n}_k,
\qquad
\mathbf{z}_k^\kappa=-\mathbf{r}_k^\kappa.
$$

于是这一轮可以直接按卡尔曼更新计算。卡尔曼增益为

$$
\mathbf{K}_\kappa
=
\tilde{\mathbf{P}}_k^{\kappa-}
\left(\mathbf{H}_k^\kappa\right)^{\mathsf{T}}
\left(
\mathbf{H}_k^\kappa
\tilde{\mathbf{P}}_k^{\kappa-}
\left(\mathbf{H}_k^\kappa\right)^{\mathsf{T}}
+
\mathbf{N}_k
\right)^{-1}.
$$

局部误差状态后验均值为

$$
\tilde{\boldsymbol{\mu}}_k^{\kappa+}
=
\tilde{\boldsymbol{\mu}}_k^{\kappa-}
+
\mathbf{K}_\kappa
\left(
-\mathbf{r}_k^\kappa
-
\mathbf{H}_k^\kappa
\tilde{\boldsymbol{\mu}}_k^{\kappa-}
\right).
$$

这里的均值是局部误差状态 $\tilde{\mathbf{x}}_k^\kappa$ 的均值，所以写作 $\tilde{\boldsymbol{\mu}}$；名义状态 $\hat{\mathbf{x}}$ 不在这个线性高斯空间里直接做均值加法，而是通过后面的注入式更新。

局部误差状态后验协方差可以用 Joseph 形式更新：

$$
\begin{aligned}
\tilde{\mathbf{P}}_k^{\kappa+}
&=
\left(
\mathbf{I}
-
\mathbf{K}_\kappa
\mathbf{H}_k^\kappa
\right)
\tilde{\mathbf{P}}_k^{\kappa-}
\left(
\mathbf{I}
-
\mathbf{K}_\kappa
\mathbf{H}_k^\kappa
\right)^{\mathsf{T}}
+
\mathbf{K}_\kappa
\mathbf{N}_k
\mathbf{K}_\kappa^{\mathsf{T}}.
\end{aligned}
$$

得到 $\tilde{\boldsymbol{\mu}}_k^{\kappa+}$ 后，把它注入当前迭代点：

$$
\hat{\mathbf{x}}_k^{\kappa+1}
=
\hat{\mathbf{x}}_k^\kappa
\boxplus
\tilde{\boldsymbol{\mu}}_k^{\kappa+}.
$$

然后重新计算累计偏移

$$
\boldsymbol{\xi}_k^{\kappa+1}
=
\hat{\mathbf{x}}_k^{\kappa+1}
\boxminus
\hat{\mathbf{x}}_k^{-},
$$

并根据新的位姿重新做数据关联、残差和雅可比。若 $\|\tilde{\boldsymbol{\mu}}_k^{\kappa+}\|$ 足够小，或达到设定的最大迭代次数，则停止内循环。设停止时的迭代索引为 $\kappa_\star$，令

$$
\hat{\mathbf{x}}_k^{+}
=
\hat{\mathbf{x}}_k^{\kappa_\star}.
$$

若只做一次迭代，且 $\hat{\mathbf{x}}_k^0=\hat{\mathbf{x}}_k^{-}$，则 $\boldsymbol{\xi}_k^0=\mathbf{0}$、$\mathbf{J}_0=\mathbf{I}$，上式退化为普通 ESKF 观测更新。

### 4.3 高维观测等效与优化视角

**高维等效。**

点云残差维度通常远大于误差状态维度。若当前扫描有 $J_k$ 个点，点到面残差每个点给 $1$ 维约束，观测维度就是 $J_k$；点到点残差每个点给 $3$ 维约束，观测维度就是 $3J_k$。直接使用卡尔曼增益

$$
\mathbf{K}
=
\tilde{\mathbf{P}}^{-}\mathbf{H}^{\mathsf{T}}
\left(
\mathbf{H}\tilde{\mathbf{P}}^{-}\mathbf{H}^{\mathsf{T}}
+
\mathbf{N}
\right)^{-1}
$$

需要对观测维度的矩阵求逆，点数很多时计算代价很高。

这个高维求逆可以等效改写成信息形式。对当前局部误差状态，卡尔曼后验均值也满足

$$
\left[
\left(\tilde{\mathbf{P}}_k^{\kappa-}\right)^{-1}
+
\left(\mathbf{H}_k^\kappa\right)^{\mathsf{T}}
\mathbf{N}_k^{-1}
\mathbf{H}_k^\kappa
\right]
\tilde{\boldsymbol{\mu}}_k^{\kappa+}
=
\left(\tilde{\mathbf{P}}_k^{\kappa-}\right)^{-1}
\tilde{\boldsymbol{\mu}}_k^{\kappa-}
-
\left(\mathbf{H}_k^\kappa\right)^{\mathsf{T}}
\mathbf{N}_k^{-1}
\mathbf{r}_k^\kappa.
$$

代入 $\tilde{\boldsymbol{\mu}}_k^{\kappa-}=-\mathbf{J}_\kappa^{-1}\boldsymbol{\xi}_k^\kappa$ 和 $\tilde{\mathbf{P}}_k^{\kappa-}=\mathbf{J}_\kappa^{-1}\tilde{\mathbf{P}}_k^{-}\mathbf{J}_\kappa^{-\mathsf{T}}$，得到不显式构造 $\mathbf{J}_\kappa^{-1}$ 的等价线性系统：

$$
\mathbf{A}_\kappa
\tilde{\boldsymbol{\mu}}_k^{\kappa+}
=
\mathbf{b}_\kappa,
\qquad
\begin{aligned}
\mathbf{A}_\kappa
&=
\mathbf{J}_\kappa^{\mathsf{T}}
\left(\tilde{\mathbf{P}}_k^{-}\right)^{-1}
\mathbf{J}_\kappa
+
\left(\mathbf{H}_k^\kappa\right)^{\mathsf{T}}
\mathbf{N}_k^{-1}
\mathbf{H}_k^\kappa \\
\mathbf{b}_\kappa
&=
-
\mathbf{J}_\kappa^{\mathsf{T}}
\left(\tilde{\mathbf{P}}_k^{-}\right)^{-1}
\boldsymbol{\xi}_k^\kappa
-
\left(\mathbf{H}_k^\kappa\right)^{\mathsf{T}}
\mathbf{N}_k^{-1}
\mathbf{r}_k^\kappa .
\end{aligned}
$$

这和上面的卡尔曼增益更新等价，但不需要在观测空间中求逆。若各点残差噪声相互独立，$\mathbf{N}_k$ 是分块对角矩阵，还可以逐点累积观测信息：

$$
\begin{aligned}
\mathbf{A}_\kappa
&=
\mathbf{J}_\kappa^{\mathsf{T}}
\left(\tilde{\mathbf{P}}_k^{-}\right)^{-1}
\mathbf{J}_\kappa
+
\sum_{j=1}^{J_k}
\left(\mathbf{H}_j^\kappa\right)^{\mathsf{T}}
\mathbf{N}_j^{-1}
\mathbf{H}_j^\kappa \\
\mathbf{b}_\kappa
&=
-
\mathbf{J}_\kappa^{\mathsf{T}}
\left(\tilde{\mathbf{P}}_k^{-}\right)^{-1}
\boldsymbol{\xi}_k^\kappa
-
\sum_{j=1}^{J_k}
\left(\mathbf{H}_j^\kappa\right)^{\mathsf{T}}
\mathbf{N}_j^{-1}
\mathbf{r}_j^\kappa .
\end{aligned}
$$

对标量点到面残差，若第 $j$ 个残差方差为 $\sigma_j^2$，鲁棒核或匹配质量给出的权重为 $w_j$，则累积项就是

$$
\mathbf{A}_\kappa
\leftarrow
\mathbf{A}_\kappa
+
\frac{w_j}{\sigma_j^2}
\left(\mathbf{H}_j^\kappa\right)^{\mathsf{T}}
\mathbf{H}_j^\kappa,
\qquad
\mathbf{b}_\kappa
\leftarrow
\mathbf{b}_\kappa
-
\frac{w_j}{\sigma_j^2}
\left(\mathbf{H}_j^\kappa\right)^{\mathsf{T}}
r_j^\kappa.
$$

这就是高维观测的等效处理：所有点的约束被压缩成状态维度上的信息矩阵和信息向量。它不是丢弃观测，而是利用矩阵求逆引理把卡尔曼增益中的“观测空间求逆”改写成“状态空间求解”。对 ESKF 来说，误差状态通常只有 $18$ 维，状态空间求解比在几千维观测空间中求逆稳定得多。

迭代结束后的协方差也可以在最终线性化点处用信息矩阵近似写为

$$
\tilde{\mathbf{P}}_k^{+}
\approx
\left[
\mathbf{J}_{\kappa_\star}^{\mathsf{T}}
\left(\tilde{\mathbf{P}}_k^{-}\right)^{-1}
\mathbf{J}_{\kappa_\star}
+
\left(\mathbf{H}_k^{\kappa_\star}\right)^{\mathsf{T}}
\mathbf{N}_k^{-1}
\mathbf{H}_k^{\kappa_\star}
\right]^{-1}.
$$

如果几何退化导致 $\mathbf{A}_\kappa$ 条件数很差，可以加入 Levenberg-Marquardt 阻尼或对退化方向做特征值约束。这属于数值求解层面的处理，不改变上面的滤波模型。

**优化视角。**

迭代 ESKF 观测更新和非线性最小二乘中的 Gauss-Newton 本质上是同一个 MAP 问题。区别只在叙述角度：优化视角说“带先验的 Gauss-Newton”，滤波视角说“迭代 ESKF 观测更新”。

如果把待估状态写成当前迭代点加局部误差，目标函数就是

$$
J\!\left(\tilde{\mathbf{x}}_k^\kappa\right)
=
\left\|
\boldsymbol{\xi}_k^\kappa
+
\mathbf{J}_\kappa
\tilde{\mathbf{x}}_k^\kappa
\right\|_{\tilde{\mathbf{P}}_k^{-}}^2
+
\left\|
\mathbf{r}_k^\kappa
+
\mathbf{H}_k^\kappa
\tilde{\mathbf{x}}_k^\kappa
\right\|_{\mathbf{N}_k}^2.
$$

对这个目标函数做 Gauss-Newton，一阶正规方程正是高维观测等效处理中得到的信息形式方程：

$$
\mathbf{A}_\kappa
\tilde{\mathbf{x}}_k^\kappa
=
\mathbf{b}_\kappa.
$$

因此该方程的解就是卡尔曼信息形式中的后验均值 $\tilde{\boldsymbol{\mu}}_k^{\kappa+}$。从这个角度看，高维观测等效处理既可以理解为卡尔曼更新的代数等价形式，也可以理解为带 IMU 先验的 Gauss-Newton 正规方程。

当先验协方差很大时，$\left(\tilde{\mathbf{P}}_k^{-}\right)^{-1}$ 很小，问题接近普通 ICP/NDT 的无先验位姿优化；当先验协方差较小时，IMU 预测会更强地约束点云配准结果。速度、零偏、重力虽然不直接出现在点云残差中，但只要它们与位姿在 $\tilde{\mathbf{P}}_k^{-}$ 中存在相关性，观测更新就会通过先验项把修正量传递到这些状态块。

## 5 状态注入与误差重置

普通观测更新得到的是误差状态估计 $\tilde{\mathbf{x}}_k^{+}$，需要把它注入名义状态。迭代观测更新中，内循环已经把每次卡尔曼更新得到的局部误差后验均值 $\tilde{\boldsymbol{\mu}}_k^{\kappa+}$ 临时注入到了 $\hat{\mathbf{x}}_k^\kappa$，收敛后可以直接把最终迭代点作为 $\hat{\mathbf{x}}_k^{+}$；若实现选择最后再统一注入，形式仍然是：

$$
\hat{\mathbf{x}}_k^{+}
=
\hat{\mathbf{x}}_k^{-}\boxplus\tilde{\mathbf{x}}_k^{+}
$$

如果采用高维观测等效处理中的信息形式，并且 $\tilde{\mathbf{P}}_k^{+}$ 已经按最终线性化点 $\hat{\mathbf{x}}_k^{\kappa_\star}$ 计算，则可以把它看成已经表达在最终名义状态附近的误差坐标中；下面的重置雅可比主要对应“先得到一次非零误差更新，再把它注入名义状态”的写法。

展开为：

$$
\begin{aligned}
\hat{\mathbf{R}}_{wi,k}^{+}
&=
\hat{\mathbf{R}}_{wi,k}^{-}
\operatorname{Exp}(\delta\boldsymbol{\theta}_k^{+}) \\
{}^{w}\hat{\mathbf{p}}_k^{+}
&=
{}^{w}\hat{\mathbf{p}}_k^{-}
+\delta{}^{w}\mathbf{p}_k^{+} \\
{}^{w}\hat{\mathbf{v}}_k^{+}
&=
{}^{w}\hat{\mathbf{v}}_k^{-}
+\delta{}^{w}\mathbf{v}_k^{+} \\
{}^{i}\hat{\mathbf{b}}_{g,k}^{+}
&=
{}^{i}\hat{\mathbf{b}}_{g,k}^{-}
+\delta{}^{i}\mathbf{b}_{g,k}^{+} \\
{}^{i}\hat{\mathbf{b}}_{a,k}^{+}
&=
{}^{i}\hat{\mathbf{b}}_{a,k}^{-}
+\delta{}^{i}\mathbf{b}_{a,k}^{+} \\
{}^{w}\hat{\mathbf{g}}_k^{+}
&=
{}^{w}\hat{\mathbf{g}}_k^{-}
+\delta{}^{w}\mathbf{g}_k^{+}
\end{aligned}
$$

注入后把误差状态均值重置为零：

$$
\tilde{\mathbf{x}}_k^{+}\leftarrow\mathbf{0}
$$

严格来说，重置会改变误差坐标系，因此协方差还应乘以重置雅可比：

$$
\tilde{\mathbf{P}}_k^{+}
\leftarrow
\mathbf{G}_k\tilde{\mathbf{P}}_k^{+}\mathbf{G}_k^\mathsf{T}
$$

对右扰动姿态误差，主要影响在姿态块：

$$
\mathbf{G}_{\theta\theta}
\approx
\mathbf{I}-\frac{1}{2}\left[\delta\boldsymbol{\theta}_k^{+}\right]_\times
$$

> 证明：
>
> 根据 [SO3扰动雅可比.md](../SO3扰动雅可比.md) 中的右扰动约定，注入前后同一个真值姿态满足
>
> $$
> \operatorname{Exp}(\delta\boldsymbol{\theta}^{\mathrm{reset}})
> =
> \operatorname{Exp}(-\delta\boldsymbol{\theta}_k^{+})
> \operatorname{Exp}(\delta\boldsymbol{\theta}_k^{+}+\boldsymbol{\eta}),
> $$
>
> 其中 $\boldsymbol{\eta}$ 是注入均值附近的小扰动。由 [BCH.md](../BCH/BCH.md) 中的结论
> $\operatorname{Exp}(\boldsymbol{\phi}+\boldsymbol{\eta})
> \approx
> \operatorname{Exp}(\boldsymbol{\phi})\operatorname{Exp}(\mathbf{J}_r(\boldsymbol{\phi})\boldsymbol{\eta})$，可得
>
> $$
> \delta\boldsymbol{\theta}^{\mathrm{reset}}
> \approx
> \mathbf{J}_r(\delta\boldsymbol{\theta}_k^{+})\boldsymbol{\eta}.
> $$
>
> 因此 $\mathbf{G}_{\theta\theta}\approx\mathbf{J}_r(\delta\boldsymbol{\theta}_k^{+})$。再由 [BCH.md](../BCH/BCH.md) 的右雅可比一阶展开
> $\mathbf{J}_r(\boldsymbol{\phi})\approx\mathbf{I}-\frac{1}{2}[\boldsymbol{\phi}]_\times$，得到
>
> $$
> \mathbf{G}_{\theta\theta}
> \approx
> \mathbf{I}-\frac{1}{2}\left[\delta\boldsymbol{\theta}_k^{+}\right]_\times.
> $$

其他普通向量块通常为单位阵。若更新量很小，工程实现中有时近似取 $\mathbf{G}_k=\mathbf{I}$；但在推导文档中应保留这个重置步骤，避免把“误差状态清零”和“协方差不变”混为一谈。

## 6 小结

ESKF 的主线可以概括为：

1. 用 IMU 测量和名义零偏传播名义状态 $\hat{\mathbf{x}}$。
2. 对真值传播和名义传播做差，得到误差状态线性模型 $\tilde{\mathbf{x}}_{k+1}^{-}=\mathbf{F}_{\tilde{x},k}\tilde{\mathbf{x}}_k^{+}+\mathbf{F}_{\eta,k}\boldsymbol{\eta}_k$。
3. 用线性模型预测误差协方差 $\tilde{\mathbf{P}}$。
4. 观测更新只估计误差状态，再把误差注入名义状态；点云紧耦合更新通常把所有点的 ICP/NDT 残差堆叠后做迭代 ESKF，本质上等价于带 IMU 先验的 Gauss-Newton。
5. 注入后重置误差状态均值，并用重置雅可比修正协方差。

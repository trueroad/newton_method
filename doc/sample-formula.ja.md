<!-- -*- coding: utf-8 -*- -->

---
title: サンプルプログラムの計算式について
author: 細田 真道
---

\newcommand\Mat[1]{\boldsymbol{#1}}
\newcommand\Vect[1]{\boldsymbol{#1}}

# はじめに

[本ライブラリ](https://github.com/trueroad/newton_method)
の三種類のサンプルプログラム

* 最小二乗法を使わないサンプル
    + [
    sample.cc
    ](https://github.com/trueroad/newton_method/blob/master/sample.cc)
* 最小二乗法を使ったサンプル（重み無し）
    + [
    sample-non_weighted.cc
    ](https://github.com/trueroad/newton_method/blob/master/sample-non_weighted.cc)
* 最小二乗法を使ったサンプル（重み付き）
    + [
    sample-weighted.cc
    ](https://github.com/trueroad/newton_method/blob/master/sample-weighted.cc)

は、いずれもGPS衛星の座標と距離からGPS受信機の位置を計算する、
というものになっています。
基本的な考え方は文献[1]に拠っていますが、使用している計算式は異なります。
ここでは、サンプルプログラムで使用している計算式と、
その計算式を使ったサンプルプログラムの処理について説明します。
サンプルプログラムで使用している座標、距離、重み等のデータについては
[
サンプルプログラムのデータについて
](https://github.com/trueroad/newton_method/blob/master/doc/sample-data.ja.md)
をご覧ください。

# 計算式

## 衛星の座標と距離

GPS は 4 つ以上の GPS 衛星の位置を既知、
それぞれの衛星から GPS 受信機までの距離を既知、
として未知である GPS 受信機の位置を求めるものです。

### 座標

$n$ 個の衛星の座標（既知）を
$S_i \left( x_i, y_i, z_i \right)$,
$\left( i = 1, 2, ..., n \right)$
として、受信機の座標（未知）を
$P \left( x_{\mathrm{p}}, y_{\mathrm{p}}, z_{\mathrm{p}} \right)$
とします。

### 距離

各衛星から受信機までの距離は、
誤差を含んだ仮の距離（疑似距離）として計測されます（既知）。
これを
$R_i$
とします。
一方、衛星と受信機の真の距離は

$$
\sqrt{ \left( x_i - x_{\mathrm{p}} \right)^2
+ \left( y_i - y_{\mathrm{p}} \right)^2
+ \left( z_i - z_{\mathrm{p}} \right)^2 }
$$

となります。
この、真の距離と疑似距離との差は
受信機の時計オフセット（未知）によるものとなります。
これを $\Delta S$ で表すと、

$$
R_i =
\sqrt{ \left( x_i - x_{\mathrm{p}} \right)^2
+ \left( y_i - y_{\mathrm{p}} \right)^2
+ \left( z_i - z_{\mathrm{p}} \right)^2 }
+ \Delta S
$$ {#eq:distance_relation}

という関係があることになります。

### 方程式

未知数は
$x_{\mathrm{p}}$, $y_{\mathrm{p}}$, $z_{\mathrm{p}}$, $\Delta S$
の 4 個ですので、
衛星が4個あれば、式 ([@eq:distance_relation]) より、

$$
\left\{
\begin{array}{l}
R_1 =
\sqrt{ \left( x_1 - x_{\mathrm{p}} \right)^2
+ \left( y_1 - y_{\mathrm{p}} \right)^2
+ \left( z_1 - z_{\mathrm{p}} \right)^2 }
+ \Delta S \\
R_2 =
\sqrt{ \left( x_2 - x_{\mathrm{p}} \right)^2
+ \left( y_2 - y_{\mathrm{p}} \right)^2
+ \left( z_2 - z_{\mathrm{p}} \right)^2 }
+ \Delta S \\
R_3 =
\sqrt{ \left( x_3 - x_{\mathrm{p}} \right)^2
+ \left( y_3 - y_{\mathrm{p}} \right)^2
+ \left( z_3 - z_{\mathrm{p}} \right)^2 }
+ \Delta S \\
R_4 =
\sqrt{ \left( x_4 - x_{\mathrm{p}} \right)^2
+ \left( y_4 - y_{\mathrm{p}} \right)^2
+ \left( z_4 - z_{\mathrm{p}} \right)^2 }
+ \Delta S
\end{array}
\right.
$$

という連立方程式を解くことによって受信機の座標
$P \left( x_{\mathrm{p}}, y_{\mathrm{p}}, z_{\mathrm{p}} \right)$
および真の距離と疑似距離の誤差
$\Delta S$
を求めることができます。
しかし、この方程式は非線形であるため、解析的に解を得ることは困難です。
そこでニュートン法で近似解を得る方法をとります。

## 関数の準備

本ライブラリでは、解きたい連立方程式について、ベクトル関数の形で表した
$\Vect{F}$
と、そのヤコビ行列 $\Mat{J}$
を用意する必要があります。

### F

まず、式 ([@eq:distance_relation]) を以下のように関数の形に書き換えます。

$$
f_i \left( x, y, z, \Delta S \right) =
\sqrt{ \left( x_i - x \right)^2
+ \left( y_i - y \right)^2
+ \left( z_i - z \right)^2 }
+ \Delta S - R_i
$$ {#eq:function_i}

そうすると、解くべき方程式は

$$
f_i \left( x_{\mathrm{p}}, y_{\mathrm{p}}, z_{\mathrm{p}}, \Delta S \right) = 0
$$

となります。これをベクトルの形で書き直すことにします。

$$
\Vect{F} \left( \Vect{X} \right) =
\left(
\begin{array}{c}
f_1 \left( x, y, z, \Delta S \right) \\
f_2 \left( x, y, z, \Delta S \right) \\
\vdots \\
f_n \left( x, y, z, \Delta S \right)
\end{array}
\right)
$$ {#eq:vector_F}

ただし、

$$
\Vect{X} =
\left(
\begin{array}{c}
x \\
y \\
z \\
\Delta S
\end{array}
\right)
,
\Vect {0} =
\left(
\begin{array}{c}
0 \\
0 \\
\vdots \\
0
\end{array}
\right)
$$

とします。
すると、解くべき方程式は、

$$
\Vect{F} \left( \Vect{X} \right) = \Vect{0}
$$ {#eq:vector_equation}

になります。

### J

次に、ここからヤコビ行列 $\Mat{J}$ を導きます。
ヤコビ行列は、

$$
\Mat{J} = \frac{\partial \Vect{F}}{\partial \Vect{X}}
$$

なので、成分で表すと、

$$
\Mat{J} \left( \Vect{X} \right) =
\left(
\begin{array}{cccc}
\frac{\partial f_1}{\partial x} & \frac{\partial f_1}{\partial y} &
\frac{\partial f_1}{\partial z} & \frac{\partial f_1}{\partial \Delta S} \\
\frac{\partial f_2}{\partial x} & \frac{\partial f_2}{\partial y} &
\frac{\partial f_2}{\partial z} & \frac{\partial f_2}{\partial \Delta S} \\
\vdots & \vdots & \vdots & \vdots \\
\frac{\partial f_n}{\partial x} & \frac{\partial f_n}{\partial y} &
\frac{\partial f_n}{\partial z} & \frac{\partial f_n}{\partial \Delta S}
\end{array}
\right)
$$ {#eq:J}

となります。各成分は式 ([@eq:function_i]) を偏微分して、

$$
\frac{\partial f_i}{\partial x}
=
\frac{- \left( x_i - x \right) }{\sqrt{
\left( x_i - x \right)^2 +
\left( y_i - y \right)^2 +
\left( z_i - z \right)^2
}}
$$

$$
\frac{\partial f_i}{\partial y}
=
\frac{- \left( y_i - y \right) }{\sqrt{
\left( x_i - x \right)^2 +
\left( y_i - y \right)^2 +
\left( z_i - z \right)^2
}}
$$

$$
\frac{\partial f_i}{\partial z}
=
\frac{- \left( z_i - z \right) }{\sqrt{
\left( x_i - x \right)^2 +
\left( y_i - y \right)^2 +
\left( z_i - z \right)^2
}}
$$

$$
\frac{\partial f_i}{\partial \Delta S} = 1
$$

となります。

## ニュートン法

ここでは、
[
ライブラリ本体
](https://github.com/trueroad/newton_method/tree/master/newton_method)
のニュートン法の処理を簡単に紹介します。
ニュートン法では、まず、
解の候補となる初期値 $\Vect{X}_0$ を適当に決めます。

$$
\Vect{X}_0 =
\left(
\begin{array}{c}
x^{(0)} \\
y^{(0)} \\
z^{(0)} \\
\Delta S^{(0)}
\end{array}
\right)
$$ {#eq:X0}

そして
$\Vect{F} \left( \Vect{X}_0 \right)$
と
$\Vect{J} \left( \Vect{X}_0 \right)$
から、次の解の候補 $\Vect{X}_1$ への修正値 $\Delta \Vect{X}$ を計算します。

$$
\Delta \Vect{X} =
\left(
\begin{array}{c}
\Delta x \\
\Delta y \\
\Delta z \\
\Delta s
\end{array}
\right)
$$

これを一般化して、$k$ 番目の解の候補を $\Vect{X}_k$ とします。

$$
\Vect{X}_k =
\left(
\begin{array}{c}
x^{(k)} \\
y^{(k)} \\
z^{(k)} \\
\Delta S^{(k)}
\end{array}
\right)
$$

そして、$\Vect{F} \left( \Vect{X}_k \right)$
と
$\Mat{J} \left( \Vect{X}_k \right)$
を使って、次の解の候補への修正値 $\Delta \Vect{X}$
を求め、 $k + 1$ 番目の解の候補を求める、ということを繰り返します。
具体的には、
$\Mat{J}$
が正方行列であれば、つまり、衛星が4個であれば、

$$
\Mat{J} \left( \Vect{X}_k \right) \Delta \Vect{X} =
- \Vect{F} \left( \Vect{X}_k \right)
$$ {#eq:deltaX_equation}

を $\Delta \Vect{X}$ について解くことになります。
これは数学的には逆行列を使えば求めることができます。

$$
\Delta \Vect{X} = - \left( \Mat{J} \left( \Vect{X}_k \right) \right)^{-1}
\Vect{F} \left( \Vect{X}_k \right)
$$

しかし、この方法は数値計算的には効率が悪く、また、
演算結果が安定しない等の問題があるとされているため、
本ライブラリではここで逆行列によらず連立一次方程式を解く
アルゴリズムを使っています。
つまり、$\Mat{A} = \Mat{J} \left( \Vect{X}_k \right)$,
$\Vect{x} = \Delta \Vect{X}$,
$\Vect{b} = - \Vect{F} \left( \Vect{X}_k \right)$
と置いて、連立一次方程式

$$
\Mat{A} \Vect{x} = \Vect{b}
$$ {#eq:simultaneous_linear_equation}

を $\Vect{x}$ について解くアルゴリズムを使います。
そして、これにより求めた $\Delta \Vect{X}$ を使って次の近似解
$\Vect{X}_{k+1}$ を求める、

$$
\Vect{X}_{k + 1} = \Vect{X}_k + \Delta \Vect{X}
$$

という操作を繰り返し行い、
解の候補を真の解に順次近づけていきます。
このとき、一定の誤差 $\varepsilon_F$
や $\varepsilon_{\Delta X}$ を決めておき、

$$
\left\| \Vect{F} \left( \Vect{X}_k \right) \right\| < \varepsilon_F
$$ {#eq:F_epsilon}

を満たした場合や

$$
\left\|
\begin{array}{c}
\frac{\Delta x}{x^{(k)}} \\
\frac{\Delta y}{y^{(k)}} \\
\frac{\Delta z}{z^{(k)}} \\
\frac{\Delta s}{\Delta S^{(k)}}
\end{array}
\right\|
< \varepsilon_{\Delta X}
$$ {#eq:deltaX_epsilon}

を満たした場合に計算を打ち切り
$\Vect{X}_k$
を近似解として出力します。

## 最小二乗法

ここでは、
[
ライブラリ本体
](https://github.com/trueroad/newton_method/tree/master/newton_method)
の最小二乗法の処理を簡単に紹介します。
衛星が4個よりも多いとき、ヤコビ行列 $\Mat{J}$ は縦長となり、
正方行列ではなくなるため、式 ([@eq:deltaX_equation]) を満たす
$\Delta \Vect{X}$ が存在しないということになります。
こういうときに、測定値には誤差があることを前提にして、
最小二乗法を使います。

### 重み無し

いずれの測定値も誤差が同じである場合や、
誤差の情報が無い場合には、すべて平等に誤差があるとして「重み無し」の
最小二乗法を使います。
式 ([@eq:simultaneous_linear_equation]) では、
$\Mat{A}$ が縦長になっているとこの式を満たす
$\Vect{x}$ が存在しません。
そこで、右辺に誤差 $\Vect{\varepsilon}$ を加えた、

$$
\Mat{A} \Vect{x} = \Vect{b} + \Vect{\varepsilon}
$$ {#eq:least_square}

という式を考えます。
その上で、誤差 $\Vect{\varepsilon}$ が最小になるように、
$\Vect{\varepsilon}$ の各成分が平等に小さくなるようにして
解 $\Vect{x}$ を求める、というのが最小二乗法になります。
具体的には、式を変形すると、

$$
\Mat{A} \Vect{x} - \Vect{b} = \Vect{\varepsilon}
$$

となるので、
$\Mat{A} \Vect{x} - \Vect{b}$
が最小になるような $\Vect{x}$ を求めることになります。
そのためには二乗して微分したものがゼロになればよくて、

$$
\frac{\mathrm{d}}{\mathrm{d} \Vect{x}}
\left\| \Mat{A} \Vect{x} - \Vect{b} \right\|^2
=0
$$

を満たせばよいということになります。
行列やベクトルの微分は少々ややこしいので飛ばしますが、

$$
\frac{\mathrm{d}}{\mathrm{d} \Vect{x}}
\left\| \Mat{A} \Vect{x} - \Vect{b} \right\|^2
= 2 \left( \Mat{A}^{\top} \Mat{A} \Vect{x} - \Mat{A}^{\top} \Vect{b} \right)
$$

となるので、ここから数学の教科書的に言う正規方程式

$$
\Mat{A}^{\top} \Mat{A} \Vect{x} = \Mat{A}^{\top} \Vect{b}
$$

が得られます。これをあてはめて、
式 ([@eq:deltaX_equation]) の代わりに、正規方程式

$$
\left( \Mat{J}^{\top} \Mat{J} \right) \Delta \Vect{X} =
- \Mat{J}^{\top} \Vect{F}
$$ {#eq:normal_equation}

を使えば、重み無しの最小二乗法が適用できます。
ただし、この方法は数値計算的には結果が悪くなるので、
正規方程式を使わない方がよいとされています [1-4] 。
本ライブラリがデフォルトで使用する連立一次方程式を解くアルゴリズムは
$\Mat{J}$ が縦長の時、
つまり式 ([@eq:least_square]) で $\Mat{A}$ が
縦長であっても自動的に
$\Vect{\varepsilon}$
が最小の（つまり最小二乗法が適用された）解$\Vect{x}$ が得られるもので、
正規方程式を使わずに $\Delta \Vect{X}$ が得られます。

### 重み付き

各衛星からの測定値について、誤差が既知である場合に「重み付き」の
最小二乗法を使います。
まず、重みのない式 ([@eq:least_square]) を成分で考えます。
方程式の数が $n$ 個（つまり衛星の数が $n$ 個）、
未知数が $m$ 個（GPS の場合は $m=4$ ）で $n>m$ （縦長）として、

$$
\Mat{A} =
\left(
\begin{array}{cccc}
a_{11} & a_{12} & \cdots & a_{1m} \\
a_{21} & a_{22} & \cdots & a_{2m} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \cdots & a_{nm}
\end{array}
\right)
,
\Vect{x} =
\left(
\begin{array}{c}
x_1 \\ x_2 \\ \vdots \\ x_m
\end{array}
\right)
,
\Vect{b} =
\left(
\begin{array}{c}
b_1 \\ b_2 \\ \vdots \\ b_n
\end{array}
\right)
,
\Vect{\varepsilon} =
\left(
\begin{array}{c}
\varepsilon_1 \\ \varepsilon_2 \\ \vdots \\ \varepsilon_n
\end{array}
\right)
$$

と置くと、式 ([@eq:least_square]) は

$$
\left(
\begin{array}{cccc}
a_{11} & a_{12} & \cdots & a_{1m} \\
a_{21} & a_{22} & \cdots & a_{2m} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n1} & a_{n2} & \cdots & a_{nm}
\end{array}
\right)
\left(
\begin{array}{c}
x_1 \\ x_2 \\ \vdots \\ x_m
\end{array}
\right)
=
\left(
\begin{array}{c}
b_1 \\ b_2 \\ \vdots \\ b_n
\end{array}
\right)
+
\left(
\begin{array}{c}
\varepsilon_1 \\ \varepsilon_2 \\ \vdots \\ \varepsilon_n
\end{array}
\right)
$$

なので、

$$
\left(
\begin{array}{c}
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1m}x_m \\
a_{21}x_2 + a_{22}x_2 + \cdots + a_{2m}x_m \\
\vdots \\
a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nm}x_m
\end{array}
\right)
=
\left(
\begin{array}{c}
b_1 \\ b_2 \\ \vdots \\ b_n
\end{array}
\right)
+
\left(
\begin{array}{c}
\varepsilon_1 \\ \varepsilon_2 \\ \vdots \\ \varepsilon_n
\end{array}
\right)
$$

となります。これは両辺とも $n$ 次の列ベクトルになっています。
これをさらにバラバラにします。 $i = 1, 2, \dots, n$ とすると、

$$
a_{i1}x_1 + a_{i2}x_2 + \cdots + a_{im}x_m = b_i + \varepsilon_i
$$

となります。
これは $i$ 番目の方程式、つまり $i$ 番目の衛星に関する式になっています。
そこで、ここに $i$ 番目の衛星の測定値の重み $w_i$ を適用します。
両辺に $w_i$ を掛けると、

$$
w_i \left(
a_{i1}x_1 + a_{i2}x_2 + \cdots + a_{im}x_m \right)
= w_i b_i + w_i \varepsilon_i
$$

となります。
これを列ベクトルの形にします。

$$
\left(
\begin{array}{c}
w_1 \left( a_{11}x_1 + a_{12}x_2 + \cdots + a_{1m}x_m \right) \\
w_2 \left( a_{21}x_2 + a_{22}x_2 + \cdots + a_{2m}x_m \right) \\
\vdots \\
w_n \left( a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nm}x_m \right)
\end{array}
\right)
=
\left(
\begin{array}{c}
w_1 b_1 \\ w_2 b_2 \\ \vdots \\ w_n b_n
\end{array}
\right)
+
\left(
\begin{array}{c}
w_1 \varepsilon_1 \\ w_2 \varepsilon_2 \\ \vdots \\ w_n \varepsilon_n
\end{array}
\right)
$$

ここで、

$$
\Mat{W}^{\frac{1}{2}} =
\left(
\begin{array}{cccc}
w_1 & 0 & \cdots & 0 \\
0 & w_2 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & w_n
\end{array}
\right)
$$ {#eq:Whalf}

という重み行列を置くと、列ベクトルの式は、

$$
\Mat{W}^{\frac{1}{2}}
\left(
\begin{array}{c}
a_{11}x_1 + a_{12}x_2 + \cdots + a_{1m}x_m \\
a_{21}x_2 + a_{22}x_2 + \cdots + a_{2m}x_m \\
\vdots \\
a_{n1}x_1 + a_{n2}x_2 + \cdots + a_{nm}x_m
\end{array}
\right)
=
\Mat{W}^{\frac{1}{2}}
\left(
\begin{array}{c}
b_1 \\ b_2 \\ \vdots \\ b_n
\end{array}
\right)
+
\Mat{W}^{\frac{1}{2}}
\left(
\begin{array}{c}
\varepsilon_1 \\ \varepsilon_2 \\ \vdots \\ \varepsilon_n
\end{array}
\right)
$$

となります。これを
$\Mat{A}$, $\Vect{x}$, $\Vect{b}$, $\Vect{\varepsilon}$
で表すと、

$$
\Mat{W}^{\frac{1}{2}} \Mat{A} \Vect{x}
=
\Mat{W}^{\frac{1}{2}} \Vect{b}+
\Mat{W}^{\frac{1}{2}} \Vect{\varepsilon}
$$

となります。
ここで、$\Mat{A}_w = \Mat{W}^{\frac{1}{2}} \Mat{A}$,
$\Vect{b}_w = \Mat{W}^{\frac{1}{2}} \Vect{b}$,
$\Vect{\varepsilon}_w = \Mat{W}^{\frac{1}{2}} \Vect{\varepsilon}$
と置けば、方程式

$$
\Mat{A}_w \Vect{x} = \Vect{b}_w + \Vect{\varepsilon}_w
$$ {#eq:weighted}

ができますので、これを重みのない式 ([@eq:least_square])
と同じように最小二乗法を適用して解きます。
この際、 $\Vect{\varepsilon}_m$ の各成分
$w_i \varepsilon_i$ はそれぞれ重み $w_i$ を含んだ形で
平等に小さくなるように解 $\Vect{x}$ を求めることになります。
つまり、重み $w_i$ が大きいときは $\varepsilon_i$ が大きく評価され、
他よりも $\varepsilon_i$ が小さくなるように、
つまり $i$ 番目の衛星の測定値が重視されるようになります。
逆に重み $w_i$ が小さいときには $\varepsilon_i$ も小さく評価され、
他よりも $\varepsilon_i$ が大きくても許容される、
つまり $i$ 番目の衛星の測定値が軽視されるようになります。
重みには例えば各測定値の標準偏差 $\sigma_i$ の逆数を使って、

$$
\Mat{W}^{\frac{1}{2}} =
\left(
\begin{array}{cccc}
\frac{1}{\sigma_1} & 0 & \cdots & 0 \\
0 & \frac{1}{\sigma_2} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & \frac{1}{\sigma_n}
\end{array}
\right)
$$

のように設定することになります。
さて、式 ([@eq:weighted]) は
式 ([@eq:least_square]) と同じように解くことができるので、
数学の教科書的には正規方程式を使うことになります。
重み無しと同様に、

$$
\frac{\mathrm{d}}{\mathrm{d} \Vect{x}}
\left\| \Mat{A}_w \Vect{x} - \Vect{b}_w \right\|^2
=0
$$

を満たすようにします。これを $\Mat{A}$, $\Vect{b}$ に戻すと、

$$
\frac{\mathrm{d}}{\mathrm{d} \Vect{x}}
\left\| \Mat{W}^{\frac{1}{2}} \Mat{A} \Vect{x}
- \Mat{W}^{\frac{1}{2}} \Vect{b} \right\|^2
=0
$$

となるので、ここから重み付きの正規方程式

$$
\left( \Mat{A}^{\top} \Mat{W} \Mat{A} \right) \Vect{x}
= \Mat{A}^{\top} \Mat{W} \Vect{b}
$$

が得られます。
これを実際に使う場合は
各測定値の分散 $\sigma_i{}^2$ の逆数を使った対角行列

$$
\Mat{W} =
\left(
\begin{array}{cccc}
\frac{1}{\sigma_1{}^2} & 0 & \cdots & 0 \\
0 & \frac{1}{\sigma_2{}^2} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & \frac{1}{\sigma_n{}^2}
\end{array}
\right)
$$

を重み行列として使い、
式 ([@eq:normal_equation]) の代わりに重み付きの正規方程式

$$
\left( \Mat{J}^{\top} \Mat{W} \Mat{J} \right) \Delta \Vect{X} =
- \Mat{J}^{\top} \Mat{W} \Vect{F}
$$

を使えばよいことになります。
もちろん、数値計算的には正規方程式を使わない方が良いとされます [1-4] 。
正規方程式を使わない方法の場合は、式 ([@eq:weighted]) から得られる

$$
\left( \Mat{W}^{\frac{1}{2}} \Mat{J} \right) \Delta \Vect{X} =
- \Mat{W}^{\frac{1}{2}} \Vect{F}
$$ {#eq:Whalf_equation}

を使い、
重みが無いときと同様に行列が縦長でも連立一次方程式が解けるアルゴリズムを使い、
解 $\Delta \Vect{X}$ を得ます [1, 3] 。
本ライブラリではどちらの方法も使えますが、
サンプルプログラムは後者の正規方程式を使わない方法を使っています。

# サンプルプログラム

サンプルプログラムで使用している座標、距離、重み等のデータについては、
ここでは中身を示すだけにとどめ詳しい説明をしません。
これらの詳細は
[
サンプルプログラムのデータについて
](https://github.com/trueroad/newton_method/blob/master/doc/sample-data.ja.md)
をご覧ください。

## 最小二乗法を使わないサンプル

以下、最小二乗法を使わないサンプル
[
sample.cc
](https://github.com/trueroad/newton_method/blob/master/sample.cc)
からの抜粋です。

### 衛星座標

```c++
// Satellite position (known, given from received data)
std::vector<std::vector<double>> Si
{
  {   -11327938.990000,     9886884.330000,    21895433.227000 },
  {     4755496.711000,    19362623.328000,    18112665.323000 },
  {    -7506201.243000,    24076860.073000,     7092793.940000 },
  {   -23085789.286000,    12409399.010000,     4602891.246000 },
};
```

4 衛星の座標 $S_i$, $\left( i = 1, 2, 3, 4\right)$ を設定しています。
$\Vect{F}$ を求める式 ([@eq:vector_F]) や
$\Mat{J}$ を求める式 ([@eq:J]) の計算に使います。

### 距離

```c++
// Observed distance (known, measurement result)
std::vector<double> Ri
{
     20690632.972000,
     23225588.018000,
     21288081.687000,
     21187099.471000,
};
```

それぞれの衛星から受信機までの距離
$R_i$, $\left( i = 1, 2, 3, 4\right)$
を設定しています。
$\Vect{F}$ を求める式 ([@eq:vector_F]) や
$\Mat{J}$ を求める式 ([@eq:J]) の計算に使います。

### 受信機座標（正解データ）

```c++
// Receiver position (actually unknown)
std::vector<double> unknown_P
{
   -3947719.36876915,   3364403.46661849,   3699487.64248845
};
```

本ライブラリの計算結果が正しいか確かめるための正解データとする目的で、
文献[1]「（その3）」のプログラムを重み無しに設定して計算した座標です。
本ライブラリのニュートン法の計算には使いません。

### 計算打ち切り条件

```c++
// Error epsilon (0.1 mm)
double epsilon {0.0001};
```

本ライブラリでニュートン法の計算を打ち切る基準として使用する値です。
式 ([@eq:F_epsilon]) や式 ([@eq:deltaX_epsilon]) の
$\varepsilon_F$, $\varepsilon_{\Delta X}$ として使います。

### 補助関数

```c++
inline double square (double x)
{
  return x * x;
}
```

二乗を計算する関数です。

```c++
inline double distance (const std::vector<double> &p1,
                        const std::vector<double> &p2)
{
  return sqrt ( square ( p1[0] - p2[0])
                + square ( p1[1] - p2[1])
                + square ( p1[2] - p2[2]) );
}
```

距離を計算する関数です。

### F

```c++
// Calculate F (to be solved as equations)
std::vector<double> calc_f (const std::vector<double> &Xarg)
{
  std::vector<double> fv;

  for (int i = 0; i < Si.size (); ++i)
    fv.push_back (distance (Si[i], Xarg) + Xarg[3] - Ri[i]);

  return fv;
}
```

$\Vect{F}$ を求める式 ([@eq:vector_F]) を計算する関数です。
ライブラリからコールバックされます。

### J

```c++
// Calculate J (Jacobian matrix of F)
std::vector<std::vector<double>> calc_j (const std::vector<double> &Xarg)
{
  std::vector<std::vector<double>> jv;

  for (int i = 0; i < Si.size (); ++i)
    {
      std::vector<double> j_row;
      double d = distance (Si[i], Xarg);

      for (int n = 0; n < 3; ++n)
        j_row.push_back (- (Si[i][n] - Xarg[n]) / d);
      j_row.push_back (1.0);

      jv.push_back (j_row);
    }

  return jv;
}
```

$\Mat{J}$ を求める式 ([@eq:J]) を計算する関数です。
ライブラリからコールバックされます。
未知数4で衛星の数も4なので4行4列の正方行列となります。

### 計算

```c++
  std::vector<double> solution;
  try
    {
      newton_method::newton_method nm;

      nm.set_function (calc_f, calc_j);
      nm.set_epsilon_F (epsilon);
      nm.set_epsilon_deltaX (epsilon);
      std::vector<double> initial_value {0.0, 0.0, 0.0, 0.0};
      solution = nm.solve (initial_value);
    }
```

ニュートン法の計算をしています。

* インスタンス `nm` を作成
* `nm.set_function ()` で $\Vect{F}$, $\Mat{J}$ を計算する関数をセット
    + それぞれ式 ([@eq:vector_F])、式 ([@eq:J]) を計算する関数です
* `nm.set_epsilon_F ()` および `nm.set_epsilon_deltaX ()` で
  計算打ち切り条件を設定
    + それぞれ式 ([@eq:F_epsilon])、式 ([@eq:deltaX_epsilon]) の条件を設定します
* `initial_value` に計算の初期値を設定
    + 式 ([@eq:X0]) に相当します
* `nm.solve ()` でニュートン法の計算を実施し、`solution` に解を得る
    + この中のライブラリの処理は、
    `nm.set_function ()` でセットされた
    $\Vect{F}$, $\Mat{J}$
    を計算する関数を繰り返しコールバックすることで
    近似解を真の解に近づけていき、
    式 ([@eq:F_epsilon]) または式 ([@eq:deltaX_epsilon])
    の条件を満たした近似解を返します

という流れになっています。

### 結果表示

```c++
  std::cout << std::endl << "Solution: P (x y z) and delta_S:" << std::endl;
  for (auto i: solution)
    std::cout << " " << i;
  std::cout << std::endl;

  std::cout << "Error distance:" << std::endl;
  double error_distance = distance (solution, unknown_P);
  std::cout << " " << error_distance << std::endl;
```

最終的に得られた解として受信機の座標
$P \left( x_{\mathrm{p}}, y_{\mathrm{p}}, z_{\mathrm{p}} \right)$
と時計オフセットによる距離の誤差
$\Delta S$
を表示します。
次いで、受信機座標の正解データとの距離を計算して表示します。

### 実行結果

サンプルプログラムの実行結果は以下のようになりました（抜粋）。

```
Solution: P (x y z) and delta_S:
 -3947717.825152 3364407.721345 3699485.385124 -14.272990
Error distance:
 5.057781
```

文献[1]のプログラムを重み無しに設定したものと比べて約 5 m 強の誤差がある、
ということになります。
文献[1]のプログラムは 8 衛星を使って最小二乗法で計算したものに対して、
このサンプルプログラムは 4 衛星しか使っていませんから、
この程度の誤差は出るものと思います。

## 最小二乗法を使ったサンプル（重み無し）

以下、最小二乗法を使ったサンプル（重み無し）
[
sample-non_weighted.cc
](https://github.com/trueroad/newton_method/blob/master/sample-non_weighted.cc)
からの抜粋です。

### 衛星座標

```c++
// Satellite position (known, given from received data)
std::vector<std::vector<double>> Si
{
  {   -11327938.990000,     9886884.330000,    21895433.227000 },
  {     4755496.711000,    19362623.328000,    18112665.323000 },
  {    -7506201.243000,    24076860.073000,     7092793.940000 },
  {   -23085789.286000,    12409399.010000,     4602891.246000 },
  {   -21893190.888000,    -2248546.668000,    14796664.928000 },
  {   -24893247.395000,     3827508.606000,    -8794926.751000 },
  {   -12971740.598000,   -10587013.898000,    21061849.442000 },
  {     7069732.127000,    22267387.067000,    12627670.276000 },
};
```

8 衛星の座標 $S_i$, $\left( i = 1, 2, ..., 8\right)$ を設定しています。
$\Vect{F}$ を求める式 ([@eq:vector_F]) や
$\Mat{J}$ を求める式 ([@eq:J]) の計算に使います。

### 距離

```c++
// Observed distance (known, measurement result)
std::vector<double> Ri
{
     20690632.972000,
     23225588.018000,
     21288081.687000,
     21187099.471000,
     21833271.739000,
     24393427.283000,
     24031767.538000,
     23630886.925000,
};
```

それぞれの衛星から受信機までの距離
$R_i$, $\left( i = 1, 2, ..., 8\right)$
を設定しています。
$\Vect{F}$ を求める式 ([@eq:vector_F]) や
$\Mat{J}$ を求める式 ([@eq:J]) の計算に使います。

### その他

その他の部分は最小二乗法を使わないサンプルと同じです。
この場合、衛星の数が8個に増えているため、
ヤコビ行列
$\Mat{J}$
が正方行列ではなくなり、8 行 4 列の縦長となります。
縦長の場合、本ライブラリはデフォルトで最小二乗法が適用された結果を返します。

### 実行結果

サンプルプログラムの実行結果は以下のようになりました（抜粋）。

```
Solution: P (x y z) and delta_S:
 -3947719.368769 3364403.466618 3699487.642488 -15.392384
Error distance:
 0.000000
```

文献[1]のプログラムを重み無しに設定したものと比べ
（表示桁数の範囲で）誤差ゼロとなりました。
本ライブラリやサンプルプログラムは、
文献[1]とは計算式、処理等が異なりますが、
同じデータ、同じ条件で計算した場合の結果は一致した、ということになります。

## 最小二乗法を使ったサンプル（重み付き）

以下、最小二乗法を使ったサンプル（重み付き）
[
sample-weighted.cc
](https://github.com/trueroad/newton_method/blob/master/sample-weighted.cc)
からの抜粋です。

### 重み行列

```c++
// Observed weight (known, measurement result)
std::vector<double> W
{
          0.34246575,
          0.41806020,
          0.44662796,
          0.33967391,
          0.33411293,
          0.29682398,
          0.30759766,
          0.31046259,
};
```

式 ([@eq:Whalf]) の重み行列 $W^{\frac{1}{2}}$ です。
実際には行列なのですが、対角行列なので、
その成分を並べた `std::vector<double>` を用意します。

### 受信機座標（正解データ）

```c++
// Receiver position (actually unknown)
std::vector<double> unknown_P
{
   -3947719.26542369,   3364403.97164603,   3699487.31861822
};
```

本ライブラリの計算結果が正しいか確かめるための正解データとする目的で、
文献[1]「（その3）」のプログラムを重み付きに設定して計算した座標です。
先のサンプルプログラムでは、
重み無し設定で計算した座標を使っていたため数字が異なります。
本ライブラリの計算結果が正しいか確かめるための正解データとする目的で使い、
本ライブラリのニュートン法の計算には使いません。

### 計算

```c++
  std::vector<double> solution;
  try
    {
      newton_method::newton_method nm;

      nm.set_function (calc_f, calc_j);
      nm.set_epsilon_F (epsilon);
      nm.set_epsilon_deltaX (epsilon);
      nm.set_weight (W);
      std::vector<double> initial_value {0.0, 0.0, 0.0, 0.0};
      solution = nm.solve<newton_method::least_square::weighted>
        (initial_value);
    }
```

重みをつけた計算をしています。
追加された処理は、

* `nm.set_weight ()` で重み行列を設定
* `nm.solve ()` で重みをつけた計算をするように設定
    + 設定した重み行列を式 ([@eq:Whalf]) とみなして、
    式 ([@eq:Whalf_equation]) の方法で計算する設定になります

の2つだけです。

### その他

その他の部分は重み付きのサンプルと同じです。

### 実行結果

サンプルプログラムの実行結果は以下のようになりました（抜粋）。

```
Solution: P (x y z) and delta_S:
 -3947719.265424 3364403.971646 3699487.318618 -15.633489
Error distance:
 0.000000
```

重み付きの場合も、
文献[1]のプログラムをと比べて（表示桁数の範囲で）誤差ゼロとなりました。

もう一度書きますが、本ライブラリやサンプルプログラムは、
文献[1]とは計算式、処理等が異なります。
それでも、同じデータ、同じ条件で計算した場合の結果は一致した、
ということになります。

# 参考文献

[1]
福島荘之介.
理解するためのGPS測位計算プログラム入門.
[
http://www.enri.go.jp/~fks442/K_MUSEN/
](http://www.enri.go.jp/~fks442/K_MUSEN/).

[2]
小柳義夫.
｢最小二乗法」事始め.
応用数理, Vol. 26 (2016) No. 1, pp. 39-42.
[
http://doi.org/10.11540/bjsiam.26.1_39
](http://doi.org/10.11540/bjsiam.26.1_39)

[3]
小柳義夫.
最小2乗法における新しい手法.
応用物理, Vol. 46 (1977) No. 1, pp. 55-60.
[
http://doi.org/10.11470/oubutsu1932.46.55
](http://doi.org/10.11470/oubutsu1932.46.55)

[4]
小柳義夫.
最小二乗法の新しいアルゴリズム.
情報処理, Vol. 23 (1982) No. 2, pp. 99-108.
[
http://id.nii.ac.jp/1001/00006446/
](http://id.nii.ac.jp/1001/00006446/)

# 更新履歴

* 2018-08-24
    + 最小二乗法の項を大幅に拡充
    + ライブラリ 2018-08-23.15 版に伴うサンプルプログラム微修正を反映
* 2017-07-16
    + 初版

# License

Copyright (C) 2017, 2018 Masamichi Hosoda. All rights reserved.

License: BSD-2-Clause

[
LICENSE
](https://github.com/trueroad/newton_method/blob/master/LICENSE) をご覧ください。

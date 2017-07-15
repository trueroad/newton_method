<!-- -*- coding: utf-8 -*- -->
# 連立方程式を解くニュートン法ライブラリ

[ 日本語 / [英語 (English)](./README.md) ]

連立方程式を解くことができるニュートン法のライブラリです。

[
https://github.com/trueroad/newton_method
](https://github.com/trueroad/newton_method)

* 優決定系（変数の数よりも方程式の数の方が多い場合）でも
最小二乗法による解を得ることができます。
* C++ で書かれています。
    + C++11 以降が必要です。
* 内部で行列演算のために
[Eigen](https://eigen.tuxfamily.org/)
というライブラリを利用していますが、
外部インタフェースは C++ 標準 STL コンテナにしてあります。
    + ライブラリを利用するソースでは Eigen を include する必要がなく、
    コンパイル時間を短縮できます。

## サンプル

三種類のサンプルを用意しました。

* 最小二乗法を使わないサンプル
    + [sample.cc](./sample.cc)
* 最小二乗法を使ったサンプル（重み無し）
    + [sample-non_weighted.cc](./sample-non_weighted.cc)
* 最小二乗法を使ったサンプル（重み付き）
    + [sample-weighted.cc](./sample-weighted.cc)

いずれのサンプルも、GPS衛星の座標と距離からGPS受信機の位置を計算する、
というものになっています。
サンプルで使用している座標、距離、重み等は文献[1]のものを
利用させていただいております（
[サンプルプログラムのデータについて](./doc/sample-data.ja.md)
）。
本サンプルで使用している計算式および処理は文献[1]とは異なります（
[
サンプルプログラムの計算式について
](https://trueroad.github.io/newton_method/doc/sample-formula.ja.html)
）。

### サンプルをビルドする

サンプル用の
[Makefile](./Makefile)
をつけてありますので、
必要に応じて書き換えてください。

* C++11 に対応したコンパイラが必要です。
    + g++ 5.4.0 等で動作確認しています。
        - 添付の
        [Makefile](./Makefile)
        では、コンパイラオプションに
        `-std=c++11` をつけてあります。
    + 他のコンパイラを使うときには、
    適宜コンパイラオプションを書き換える等して、
    C++11が使えるようにしてください。
* [Eigen](https://eigen.tuxfamily.org/)
が必要になります。
    + バージョン 3.2.5 や 3.3.4 で動作確認しています。
    + eigen3 / libeigen3-dev / eigen3-devel
    等のパッケージがインストールしてあれば
    （pkg-config で eigen3 が見つかるようになっていれば）
    そのまま `make` コマンドでビルドができます。
    + 手動でインストールした場合等は、`EIGEN_CXXFLAGS` を書き換えて、
    Eigen のインストール先が include できるようにしてください。
* デバッグ表示が有効になっています。
    + コンパイルオプションに `-DDEBUG_NEWTON_METHOD` を
    つけてあり、計算の途中経過を表示するようになっています。
    + 途中経過の表示が必要なければ、このオプションを削除してください。

なお、私の環境では Eigen 3.2.5 だと
[newton.cc](./newton.cc)
のコンパイル中に警告が出ます。
Eigen 3.3.4 なら警告がでません。

## 参考文献

[1]
福島荘之介.
理解するためのGPS測位計算プログラム入門.
[
http://www.enri.go.jp/~fks442/K_MUSEN/
](http://www.enri.go.jp/~fks442/K_MUSEN/).

## License

Copyright (C) 2017 Masamichi Hosoda. All rights reserved.

License: BSD-2-Clause

[LICENSE](./LICENSE) をご覧ください。

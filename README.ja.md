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

サンプルをビルドする
[Makefile](./Makefile) とライブラリをビルドする
[newton_method/Makefile](./newton_method/Makefile)
をつけてありますので、
必要に応じて書き換えてください。

* C++11 に対応したコンパイラが必要です。
    + g++ 7.3.0 等で動作確認しています。
        - 添付の
        [Makefile](./Makefile) と
        [newton_method/Makefile](./newton_method/Makefile)
        では、コンパイラオプションに
        `-std=c++11` をつけてあります。
    + 他のコンパイラを使うときには、
    適宜コンパイラオプションを書き換える等して、
    C++11が使えるようにしてください。
* [Eigen](https://eigen.tuxfamily.org/)
が必要になります。
    + バージョン 3.3.4 で動作確認しています。
    + eigen3 / libeigen3-dev / eigen3-devel
    等のパッケージがインストールしてあれば
    （pkg-config で eigen3 が見つかるようになっていれば）
    そのまま `make` コマンドでビルドができます。
    + 手動でインストールした場合等は、`CPPFLAGS_EIGEN` を書き換えて、
    Eigen のインストール先が include できるようにしてください。
* デバッグ表示が有効になっています。
    + コンパイルオプションに `-DDEBUG_NEWTON_METHOD` を
    つけてあり、計算の途中経過を表示するようになっています。
    + 途中経過の表示が必要なければ、このオプションを削除してください。

## 参考文献

[1]
福島荘之介.
理解するためのGPS測位計算プログラム入門.
[
http://www.enri.go.jp/~fks442/K_MUSEN/
](http://www.enri.go.jp/~fks442/K_MUSEN/).

## 更新履歴

* 2018-08-26
    + ライブラリ 2018-08-26.02 版
        - 高速化を図った `solve_fast ()` を追加（実験的実装）
            - ポインタを使いコンストラクタ、デストラクタ、
              メモリ再配置などに伴うオーバーヘッドを回避
                - ポインタは間違えやすく面倒なので使いたくなくて
                  本当は std::span や mdspan あたりを使いたいのですが、
                  使用している環境ではまだ使えず、
                  仕方なくポインタで実装したため
                  一応実験的実装としています。
                - F と J を一つの関数で計算し途中式の使いまわしを可能に
            - `solve_fast ()` を使った高速化サンプル
              [sample-fast.cc](./sample-fast.cc) を追加
                - 本ライブラリのデバッグモードを off にして
                  sample-weight と sample-fast それぞれの `main ()` を
                  1 万回ループさせるベンチマークを行ったところ、
                  手元の環境で 8 秒台から 7 秒台に短縮でき
                  1 割強の高速化が実現できました。
* 2018-08-24
    + [
サンプルプログラムの計算式について
](https://trueroad.github.io/newton_method/doc/sample-formula.ja.html) 更新
          - 最小二乗法の項を大幅に拡充
          - ライブラリ 2018-08-23.15 版に伴うサンプルプログラム微修正を反映
* 2018-08-23
    + ライブラリ 2018-08-23.15 版
        - アルゴリズムや最小二乗法処理の指定方法を変更
        - テンプレート化およびソースファイル分割により
          使用しないアルゴリズムをリンクしないで済むようになり、
          実行ファイルサイズの大幅減を実現
* 2017-07-20
    + ライブラリ 2017-07-20.13 版
        - イテレーション回数超過時に例外発生するか否か指定可に
        - イテレーション完了ステータス取得可に
* 2017-07-16
    + [
サンプルプログラムの計算式について
](https://trueroad.github.io/newton_method/doc/sample-formula.ja.html) 追加
* 2017-07-01
    + [サンプルプログラムのデータについて](./doc/sample-data.ja.md) 追加
* 2017-06-24
    + ライブラリ 2017-06-24.16 版（初版）

## License

Copyright (C) 2017, 2018 Masamichi Hosoda. All rights reserved.

License: BSD-2-Clause

[LICENSE](./LICENSE) をご覧ください。

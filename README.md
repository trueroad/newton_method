<!-- -*- coding: utf-8 -*- -->
# Newton's method library to solve simultaneous equations

[ [Japanese (日本語)](./README.ja.md) / English ]

This is a Newton's method library that can solve simultaneous equations.

* You can get solutions with the least-squares method
even if overdetermined systems
(the number of equations is greater than the number of variables).
* It is written in C++.
    + C++11 or later is required.
* [Eigen](https://eigen.tuxfamily.org/)
is used internally for matrix operation,
however the external interface is the C++ standard STL container.
    + For sources that use libraries,
    you do not need to include Eigen and you can shorten compile time.

## Samples

There are three samples.

* Non least square
    + [sample.cc](./sample.cc)
* Overdetermined system, non weighted
    + [sample-non_weighted.cc](./sample-non_weighted.cc)
* Overdetermined system, weighted
    + [sample-weighted.cc](./sample-weighted.cc)

All samples are position calculation like GPS/GNSS receiver.
The coordinates, distances, weights etc. used in the samples
are based on the document [1]
([Sample program's data](./doc/sample-data.md))
.
The calculation formulas and processing used in these samples
are different from the document [1].

### Build samples

[Makefile](./Makefile) for samples is contained,
so please rewrite if necessary.

* C++11 or later compiler is required.
    + g++ 5.4.0 can compile it.
        - In the [Makefile](./Makefile) ,
        compiler option `-std=c++11` is used.
    + If you would like another compiler,
    please rewrite the compiler option as appropriate for C++11 can be used.
* [Eigen](https://eigen.tuxfamily.org/) is required.
    + Eigen 3.2.5 and 3.3.4 etc. can be used.
    + If the package such as eigen3 / libeigen3-dev / eigen3-devel etc.
    is installed (that is, pkg-config can find eigen3),
    you can build it as it is with the `make` command.
    + If you installed Eigen without using packages,
    please rewrite `EIGEN_CXXFLAGS`
    for your compiler can include Eigen header files.
* Debug mode is enabled.
    + The compile option `-DDEBUG_NEWTON_METHOD` is enabled.
    It shows the progress of the calculation.
    + Please remove this option if you would not like to show the progress.

Note, in my environment,
g++ with Eigen 3.2.5 shows some warnings
while compiling [newton.cc](./newton.cc).
g++ with Eigen 3.3.4 does not show any warnings.

## References

[1]
FUKUSHIMA Sonosuke.
理解するためのGPS測位計算プログラム入門.
[
http://www.enri.go.jp/~fks442/K_MUSEN/
](http://www.enri.go.jp/~fks442/K_MUSEN/)
(in Japanese).

## License

Copyright (C) 2017 Masamichi Hosoda. All rights reserved.

License: BSD-2-Clause

See [LICENSE](./LICENSE).

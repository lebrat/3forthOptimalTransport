# 3/4 discrete optimal transport
______________________
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3474563.svg)](https://doi.org/10.5281/zenodo.3474563)

This codes are the implementation of the following [paper](https://arxiv.org/abs/1806.09537). It allows to compute the exact *L²*  optimal transport between a **polyline** (union of line segments) and a **point cloud**. This toolbox allows the calculation in **2D** and **3D**. We provide hands-on tutorials on the [wiki](https://github.com/lebrat/3forthOptimalTransport/wiki).

The codes are released for Linux platforms (tested for *Mint 18 Cinnamon 64-bit* and *Ubuntu 19.04 Disco Dingo*). 

The back-end computations are coded in `C++` and make use of the computational geometry library [CGAL](www.cgal.org) and the linear algebra library [Eigen3](http://eigen.tuxfamily.org). We provide a `python 3.7` interface by using the wrapper [swig](http://www.swig.org/).



## Authors
This software was developed by:
* Frédéric de Gournay
* Jonas Kahn
* [Léo Lebrat](lebrat.org)

All members of the [Toulouse institute of Mathematics](https://www.math.univ-toulouse.fr/?lang=en) France, this project was granted by the **ANR-17-CE23-0013**

## License

This software is open source distributed under license MIT.

# CppTensority

## Prerequisites

+ [OpenBLAS](https://github.com/xianyi/OpenBLAS)


## How to speed it up

If you don't mind the portability issue, you may try to opt it by working on the `cblas_dgemm()` using:

+ _MKL_ on a Intel machine, or even
+ GPU computing lib, for example, _cuBLAS_
    * `cublasGemmEx` may draw your interests and you would have to convert the types and prepare (set up) the matrix beforehand


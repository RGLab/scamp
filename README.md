[![Travis-CI Build Status](https://travis-ci.org/RGLab/scamp.svg?branch=master)](https://travis-ci.org/RGLab/scamp)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/RGLab/scamp?branch=master&svg=true)](https://ci.appveyor.com/project/RGLab/scamp)

# scamp

The `scamp` package implements the SCAMP algorithm described in [Selective Clustering Annotated using Modes of Projections](https://arxiv.org/abs/1807.10328).

The version implemented on branch "master" is designed to interoperate with the [FAUST](https://github.com/RGLab/FAUST) algorithm.

The version implemented on branch "arxiv" is designed to reproduce experiments described in the arXiv manuscript.

## Installation

Currently `scamp` must be installed from its source.

The most recent version can be installed from [github](https://github.com/RGlab/scamp) using [devtools](https://github.com/r-lib/devtools) in R.

    library(devtools)
    devtools::install_github("RGLab/scamp")

## Citation

If you end up using `scamp` to cluster data in your work,
please consider citing [Selective Clustering Annotated using Modes of Projections](https://arxiv.org/abs/1807.10328).


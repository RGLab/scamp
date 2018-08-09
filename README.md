[![Travis-CI Build Status](https://travis-ci.org/RGLab/scamp.svg?branch=master)](https://travis-ci.org/RGLab/scamp)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/RGLab/scamp?branch=master&svg=true)](https://ci.appveyor.com/project/RGLab/scamp)

# scamp

The `scamp` package implements the SCAMP algorithm described in [Selective Clustering Annotated using Modes of Projections](https://arxiv.org/abs/1807.10328).

## Installation

Currently `scamp` must be installed from its source.

The most recent version can be installed from [github](https://github.com/RGlab/scamp) using [devtools](https://github.com/r-lib/devtools) in R.

    library(devtools)
    devtools::install_github("RGLab/scamp")

To install with vignettes, instead run

    library(devtools)
    devtools::install_github("RGLab/scamp", build_vignettes=T)

This takes longer since the vignettes must be built from source.

Two vignettes are available.
After loading `scamp`, type `vignette('scampIntro')` to read a vignette discussing how to use the `scamp` function in R.
Type `vignette('nDip')` to read a vignette discussing the N-dip facilities that are used in the search for candidate clusters.

## Reproducing results

Scripts to reproduce the main results of section 4 in [Selective Clustering Annotated using Modes of Projections](https://arxiv.org/abs/1807.10328)
are stored in the directory `reproducibility`.

## Citation

If you end up using `scamp` to cluster data in your work,
please consider citing [Selective Clustering Annotated using Modes of Projections](https://arxiv.org/abs/1807.10328).


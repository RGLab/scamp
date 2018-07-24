# scamp

The `scamp` package implements the SCAMP algorithm described in [Selective Clustering Annotation using Modes of Projections](https://arxiv.org).

## Installation

Currently `scamp` must be installed from its source.

The most recent version can be installed from [github](https://github.com/RGlab/scamp) using [devtools](https://github.com/r-lib/devtools) in R.

    library(devtools)
    devtools::install_github("RGLab/scamp")

To install with a vignette, instead run

    library(devtools)
    devtools::install_github("RGLab/scamp", build_vignettes=T)

Then, after loading `scamp`, type `vignette('scampIntro')`.
    
## Citation

If you end up using `scamp` to cluster data in your work,
please consider citing [Selective Clustering Annotated using Modes of Projections](https://arxiv.org).

## Reproducing results

Scripts to reproduce results in section 4 of [Selective Clustering Annotated using Modes of Projections](https://arxiv.org)
are stored in the directory `reproducibility`.
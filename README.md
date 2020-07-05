
<!-- README.md is generated from README.Rmd. Please edit that file -->
dasper
======

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) [![Codecov test coverage](https://codecov.io/gh/dzhang32/dasper/branch/master/graph/badge.svg)](https://codecov.io/gh/dzhang32/dasper?branch=master) [![BioC status](http://www.bioconductor.org/shields/build/release/bioc/dasper.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/dasper) <!-- badges: end -->

The aim of `dasper` is to **d**etect **a**berrant **sp**licing **e**vents from **R**NA-seq data. By comparing patient RNA-seq data to a set of user-defined controls, `dasper` will score each splicing event in the patient with an anomaly score representing the degree to which that splicing event looks abnormal. This package is still in the experimental stage of development currently and subject to major changes in the future.

Installation instructions
-------------------------

Get the latest stable `R` release from [CRAN](http://cran.r-project.org/). Then install `dasper` using from [Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("dasper")
```

And the development version from [GitHub](https://github.com/) with:

``` r
BiocManager::install("dzhang32/dasper")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library("dasper")
## to be added
```

Citation
--------

Below is the citation output from using `citation('dasper')` in R. Please run this yourself to check for any updates on how to cite **dasper**.

``` r
print(citation("dasper"), bibtex = TRUE)
#> 
#> To cite package 'dasper' in publications use:
#> 
#>   David Zhang (2020). dasper: Detecting abberant splicing events from
#>   RNA-sequencing data. R package version 0.99.0.
#>   https://github.com//dzhang32/dasper
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {dasper: Detecting abberant splicing events from RNA-sequencing data},
#>     author = {David Zhang},
#>     year = {2020},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com//dzhang32/dasper},
#>   }
```

Please note that the `dasper` was only made possible thanks to many other R and bioinformatics software authors, which are cited either in the vignettes and/or the paper(s) describing this package.

Code of Conduct
---------------

Please note that the `dasper` project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

Development tools
-----------------

-   Continuous code testing is possible thanks to [GitHub actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/) through *[usethis](https://CRAN.R-project.org/package=usethis)*, *[remotes](https://CRAN.R-project.org/package=remotes)*, *[sysreqs](https://github.com/r-hub/sysreqs)* and *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized to use [Bioconductor's docker containers](https://www.bioconductor.org/help/docker/) and *[BiocCheck](https://bioconductor.org/packages/3.9/BiocCheck)*.
-   Code coverage assessment is possible thanks to [codecov](https://codecov.io/gh) and *[covr](https://CRAN.R-project.org/package=covr)*.
-   The [documentation website](http://dzhang32.github.io/dasper) is automatically updated thanks to *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
-   The code is styled automatically thanks to *[styler](https://CRAN.R-project.org/package=styler)*.
-   The documentation is formatted thanks to *[devtools](https://CRAN.R-project.org/package=devtools)* and *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using *[biocthis](https://github.com/lcolladotor/biocthis)*.

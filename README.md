
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dasper

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R build
status](https://github.com/dzhang32/dasper/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/dzhang32/dasper/actions)
[![Codecov test
coverage](https://codecov.io/gh/dzhang32/dasper/branch/master/graph/badge.svg)](https://codecov.io/gh/dzhang32/dasper?branch=master)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/dasper.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/dasper)
<!-- badges: end -->

The aim of `dasper` is to **d**etect **a**berrant **sp**licing
**e**vents from **R**NA-seq data. By comparing patient RNA-seq data to a
set of controls, `dasper` will score each splicing event in the patient
representing the degree to which that splicing event looks abnormal. For
a detailed guide on the usage of `dasper`, check out the vignette
[here](https://dzhang32.github.io/dasper/articles/dasper.html).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `dasper` using from
[Bioconductor](http://bioconductor.org/) the following code:

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

## Citation

Below is the citation output from using `citation('dasper')` in R.
Please run this yourself to check for any updates on how to cite
**dasper**.

``` r
print(citation("dasper"), bibtex = TRUE)
#> 
#> dzhang32 (2020). _Detecting abberant splicing events from
#> RNA-sequencing data_. doi: 10.18129/B9.bioc.dasper (URL:
#> https://doi.org/10.18129/B9.bioc.dasper),
#> https://github.com/dzhang32/dasper - R package version 0.99.9, <URL:
#> http://www.bioconductor.org/packages/dasper>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {Detecting abberant splicing events from RNA-sequencing data},
#>     author = {{dzhang32}},
#>     year = {2020},
#>     url = {http://www.bioconductor.org/packages/dasper},
#>     note = {https://github.com/dzhang32/dasper - R package version 0.99.9},
#>     doi = {10.18129/B9.bioc.dasper},
#>   }
#> 
#> dzhang32 (2020). "Detecting abberant splicing events from
#> RNA-sequencing data." _bioRxiv_. doi: 10.1101/TODO (URL:
#> https://doi.org/10.1101/TODO), <URL:
#> https://www.biorxiv.org/content/10.1101/TODO>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Detecting abberant splicing events from RNA-sequencing data},
#>     author = {{dzhang32}},
#>     year = {2020},
#>     journal = {bioRxiv},
#>     doi = {10.1101/TODO},
#>     url = {https://www.biorxiv.org/content/10.1101/TODO},
#>   }
```

Please note that the `dasper` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `dasper` project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Development tools

  - Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*,
    *[sysreqs](https://github.com/r-hub/sysreqs)* and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.11/BiocCheck)*.
  - Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
  - The [documentation website](http://dzhang32.github.io/dasper) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
  - The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
  - The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

In particular, I am very grateful to
[Leo](http://lcolladotor.github.io/) for his time and advice throughout
the development of `dasper`. The transition of `dasper`
Bioconductor-friendly package was made possible thanks to his
*[biocthis](https://github.com/lcolladotor/biocthis)* package.


# `gvcR`: Genotypic Variance Components

###### Version : [0.3.0](https://myaseen208.com/gvcR/); License: [GPL-3](https://www.r-project.org/Licenses/)

##### \*Muhammad Yaseen<sup>1,2,3</sup>, Sami Ullah<sup>4</sup>

1.  [Asian Development Bank (ADB), Islamabad,
    Pakistan](https://myaseen208.com/)
2.  [Benazir Income Support Programme (BISP), Islamabad,
    Pakistan](https://myaseen208.com/)
3.  [Department of Mathematics and Statistics, University of Agriculture
    Faisalabad, Pakistan](https://myaseen208.com/)
4.  College of Agriculutre, University of Sargodha, Pakistan

------------------------------------------------------------------------

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-2.10.0-6666ff.svg)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/gvcR)](https://cran.r-project.org/package=gvcR)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/gvcR?color=green)](https://CRAN.R-project.org/package=gvcR)
<!-- [![packageversion](https://img.shields.io/badge/Package%20version-0.2.3.3-orange.svg)](https://github.com/myaseen208/gvcR) -->

[![develVersion](https://img.shields.io/badge/devel%20version-0.3.0-orange.svg)](https://github.com/myaseen208/gvcR)

<!-- [![GitHub Download Count](https://github-basic-badges.herokuapp.com/downloads/myaseen208/gvcR/total.svg)] -->

[![Project Status:
WIP](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--09--14-yellowgreen.svg)](https://github.com/myaseen208/gvcR)

------------------------------------------------------------------------

## Description

Functionalities to compute model based genetic components i.e. genotypic
variance, phenotypic variance and heritability for given traits of
different genotypes from replicated data using methodology explained by
Burton, G. W. & Devane, E. H. (1953)
([doi:10.2134/agronj1953.00021962004500100005x](https://doi.org/10.2134/agronj1953.00021962004500100005x))
and Allard, R.W. (2010, <ISBN:8126524154>).

   

## Installation

The package can be installed from CRAN as follows:

``` r
install.packages("gvcR", dependencies = TRUE)
```

 

The development version can be installed from github as follows:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("myaseen208/gvcR")
```

   

## What’s new

To know whats new in this version type:

``` r
news(package = "gvcR")
```

## Links

[CRAN page](https://cran.r-project.org/package=gvcR)

[Github page](https://github.com/myaseen208/gvcR)

[Documentation website](https://myaseen208.com/gvcR/)

[Companion website](https://myaseen208.com/EDATR/)

## Citing `gvcR`

To cite the R package `gvcR` in publications use:

``` r
citation("gvcR")
Please, support this project by citing it in your publications!

  Yaseen M, Ullah S (2018). _gvcR: Genotypic Variance Components_.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {gvcR: Genotypic Variance Components},
    author = {Muhammad Yaseen and Sami Ullah},
    year = {2018},
    journal = {The Comprehensive R Archive Network (CRAN)},
  }
```

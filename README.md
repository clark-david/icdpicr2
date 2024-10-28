
<!-- README.md is generated from README.Rmd. Please edit that file -->

# icdpicr2

<!-- badges: start -->
<!-- badges: end -->

The goal of icdpicr2 is to …\[The purpose of this package is\] to
provide injury researchers with easy methods for categorization of
injuries described using ICD-10 diagnosis codes and for deriving
measures of injury severity.

ICDPIC (International Classification of Diseases Programs for Injury
Categorization) was originally developed for ICD Version 9 Clinical
Modification (ICD-9-CM) diagnosis codes using Stata, providing an easy
way to convert these codes to standard injury severity scores and
categories. After the introduction of ICD-10-CM to US hospitals in 2015,
an update to accommodate this change was developed using R, and was
eventually published on CRAN as package “icdpicr”.

Package “icdpicr2” is being developed as a successor to “icdpicr”. The
“Implementation Details” document provided with this package describes
the history of injury severity scoring, ICDPIC, and “icdpicr”; it
discusses issues with the current version of “icdpicr” and the aims for
“icdpicr2”.

Package “icdpicr2” consists primarily of the function cat_trauma2. This
function first reads in user data in a specified format, then calculates
AIS, ISS, NISS, mortality predictions, and injury mechanisms. It returns
the original data file with these fields added. Further details about
the options available in cat_trauma2 are provided in the help file for
this function.

Other functions in the package derive an “International Classification
of Diseases-related Injury Severity Score” (ICISS) and a categorization
by injury type and location presented as a “framework” by the US CDC.
Further information is provided in the Implementation Details and the
help files for each function.

If the so-called “vignette” with the implementation details did not
download with the installation of ICDPICR2 from GitHub, perform the
installation again with the commands
\>devtools::install_github(“clark-david/icdpicr2”, build_vignettes=TRUE)
\>require(icdpicr2) \>help(package=”icdpicr2”)

## Installation

You can install the development version of icdpicr2 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("clark-david/icdpicr2")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(icdpicr2)
#> 
#> Attaching package: 'icdpicr2'
#> The following object is masked from 'package:graphics':
#> 
#>     frame
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

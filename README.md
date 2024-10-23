
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ICDPIC

ICDPIC – International Classification of Diseases Programs for Injury
Categorization was originally developed in 2008 using Stata and using
ICD Version 9 Clinical Modification (ICD-9-CM) diagnosis codes,
providing an easy way to convert these codes to standard injury severity
scores and categories. After the introduction of ICD-10-CM to US
hospitals in 2015, an update to accommodate this change was developed
using R, and was eventually published on CRAN as package icdpicr.

Package icdpicr2 is being developed as a successor to icdpicr. The
“Implementation” vignette in this package gives more detail about the
history of injury severity scoring, ICDPIC, and icdpicr; it discusses
issues with the current version of icdpic and the aims for icdpicr2.

Package icdpicr2 consists primarily of a single function, `cat_trauma2`.
This function first reads in user data in a specified format. It then
calculates AIS, ISS, NISS, mortality predictions from the TQP and NIS
models, and injury mechanisms. It returns the original data file with
these fields added. Further details about the options available in
`cat_trauma2` are provided in the help file for this function.

## Installation

You can install the development version of icdpicr2 from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("clark-david/icdpicr2")
```

## Example

``` r
library(icdpicr2)

df_in <- read.table(header = TRUE, text = "
    ident   dx1       dx2       dx3
    31416   S32110A   S3251     NA
    31417   S72141A   T07XXXA   D62 
")
df_out <- cat_trauma2(df_in, "dx", TRUE, FALSE)
#> Loading required package: tidyverse
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
#> ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
#> ✔ purrr     1.0.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
df_out
#>   ident     dx1 sev_1 issbr_1     dx2 sev_2 issbr_2  dx3 sev_3 issbr_3
#> 1 31416 S32110A     1       A   S3251     1       A <NA>    NA    <NA>
#> 2 31417 S72141A     3       E T07XXXA     1       G  D62    NA    <NA>
#>   mxaisbr_F mxaisbr_H mxaisbr_C mxaisbr_A mxaisbr_E mxaisbr_G maxais riss niss
#> 1         0         0         0         1         0         0      1    1    2
#> 2         0         0         0         0         3         1      3   10   10
#>   ecode_1 mechmaj1 mechmin1 intent1 ecode_2 mechmaj2 mechmin2 intent2 ecode_3
#> 1    <NA>     <NA>       NA    <NA>    <NA>     <NA>       NA    <NA>    <NA>
#> 2    <NA>     <NA>       NA    <NA>    <NA>     <NA>       NA    <NA>    <NA>
#>   mechmaj3 mechmin3 intent3 ecode_4 mechmaj4 mechmin4 intent4  PmortTQIP
#> 1     <NA>       NA    <NA>    <NA>     <NA>       NA    <NA> 0.01167072
#> 2     <NA>       NA    <NA>    <NA>     <NA>       NA    <NA> 0.02082725
#>     PmortNIS
#> 1 0.02081440
#> 2 0.01430465
```

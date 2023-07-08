<!-- badges: start -->
[![](https://img.shields.io/badge/lifecycle-stable-yellow.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/github/last-commit/mrbakhsh/MACP.svg)](https://github.com/mrbakhsh/MACP/commits/main)
[![CRAN Version](https://www.r-pkg.org/badges/version/MACP)](https://cran.r-project.org/package=MACP)
[![Downloads from the RStudio CRAN mirror](https://cranlogs.r-pkg.org/badges/MACP)](https://cranlogs.r-pkg.org/badges/MACP)
<!-- badges: end -->


# MACP
Systematic mapping of multiprotein complexes formed by protein-protein interactions (PPIs) can enhance our knowledge and mechanistic basis of how proteins function in the cells. Co-fractionation coupled with mass spectrometry (CF-MS) is gaining momentum as a cost-effective strategy for charting protein assemblies under native conditions using high-resolution chromatography separation techniques (e.g., size-exclusion and  ion-exchange) without the need for antibodies or tagging of individual proteins. To capture high-quality PPIs from CF-MS co-elution profile, we have developed a well standardized and fully automated CF-MS data analysis software toolkit, referred to as MACP (Macromolecular Assemblies from the Co-elution Profile) in an open-source R package, beginning with the processing of raw co-elution data to reconstruction of high-confidence PPI networks via supervised machine-learning and underlying protein complexes using unsupervised approach.

## Installation

You can install the `MACP` from CRAN using:

```r
install.packages('MACP')
```

To install the development version in `R`, run:
  
```r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}
devtools::install_github("BabuLab-UofR/MACP")
```

For a detailed introduction to MACP, see the [vignette](https://cran.r-project.org/web/packages/MACP/vignettes/MACP_tutorial.html).

## Contribute

Check the github page for [source code](https://github.com/BabuLab-UofR/MACP)

## License
This project is licensed under the MIT License - see the LICENSE.md file for more details.

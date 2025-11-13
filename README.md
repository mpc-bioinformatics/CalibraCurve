
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CalibraCurve

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/mpc-bioinformatics/CalibraCurve)](https://github.com/mpc-bioinformatics/CalibraCurve/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/mpc-bioinformatics/CalibraCurve)](https://github.com/mpc-bioinformatics/CalibraCurve/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/CalibraCurve.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/CalibraCurve)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/CalibraCurve.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/CalibraCurve)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/CalibraCurve.svg)](http://bioconductor.org/packages/stats/bioc/CalibraCurve/)
[![Bioc
support](https://bioconductor.org/shields/posts/CalibraCurve.svg)](https://support.bioconductor.org/tag/CalibraCurve)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/CalibraCurve.svg)](https://bioconductor.org/packages/release/bioc/html/CalibraCurve.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/CalibraCurve.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/CalibraCurve/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/CalibraCurve.svg)](https://bioconductor.org/packages/release/bioc/html/CalibraCurve.html#since)
[![check-bioc](https://github.com/mpc-bioinformatics/CalibraCurve/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/mpc-bioinformatics/CalibraCurve/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/mpc-bioinformatics/CalibraCurve/graph/badge.svg)](https://app.codecov.io/gh/mpc-bioinformatics/CalibraCurve)
<!-- badges: end -->

Targeted mass-spectrometry-based techniques allow accurate quantitative
measurements of analytes in complex matrices. They are used in different
fields like proteomics, lipidomics or metabolomics to validate results
from global analyses. Furthermore, they are used to develop robust
assays with the potential to absolutely quantify the analyte of
interest.

During assay development, often experiments using concentration or
dilution series are conducted. Here, samples with known amount of the
analyte of interest are generated and measured by mass spectrometry,
usually in relicates.

This kind of data is used by CalibraCurve to generate calibration
curves, which can be the basis for predicting concentrations from
intensities in new data. However, a linear relationship between
concentration and intensity is often limited to a certain range of
concentrations, the so called linear range. CalibraCurve calculates the
linear range (lower and upper limits of quantification) using a
well-designed algorithm. Quality of the data and the accuracy of the
produced calibration curve is assessed in several computational steps.

Visualizations of the calibration curves and the underlying data basis
are generated and can be customized.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `CalibraCurve` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("CalibraCurve")
```

And the development version from
[GitHub](https://github.com/mpc-bioinformatics/CalibraCurve) with:

``` r
BiocManager::install("mpc-bioinformatics/CalibraCurve")
```

## Implementation

CalibraCurve is written in R and available as an R package. Furthermore,
a nextflow-workflow running this script is available at
<https://github.com/mpc-bioinformatics/CalibraCurve_NF>.

Formerly, CalibraCurve was available as a KNIME-workflow based on an R
script. The last version of this workflow can be found here:
<https://github.com/mpc-bioinformatics/CalibraCurve/releases/tag/v_2_0>
. However, development of the KNIME-workflow was discontinued and
replaced by the R package and nextflow workflow.

## Usage

For details on the usage of CalibraCurve and some examples please see
the vignettes.

## Publication

Kohl M, Stepath M, Bracht T, Megger DA, Sitek B, Marcus K, Eisenacher M.
CalibraCurve: A Tool for Calibration of Targeted MS-Based Measurements.
Proteomics. 2020 Jun;20(11):e1900143. doi: 10.1002/pmic.201900143. Epub
2020 Mar 6. PMID: 32086983.

## Funding

The development and maintanence of CalibraCurve is funded by de.NBI
(<https://www.denbi.de/>) and CUBiMed.RUB
(<https://www.cubimed.ruhr-uni-bochum.de/index.html.en>). We offer also
other cool tools and consulting for statistics, bioninformatics and
machine learning!

<img width="400" alt="deNBI_Logo_cmyk" src="https://github.com/user-attachments/assets/eac6982f-d2ae-455e-94b4-de2b754700f6" />

<img width="400" height="156" alt="CUBiMedRUB-logo-small" src="https://github.com/user-attachments/assets/1be9ab3e-825c-45b8-8bfe-5a41b8cb2f90" />



## Feedback

If you used this de.NBI service please take the time to answer our short survey to help us improve our services:

<https://de.surveymonkey.com/r/denbi-service?sc=bioinfra-prot&tool=calibracurve>

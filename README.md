# CalibraCurve

## Background

Targeted proteomics techniques allow accurate quantitative measurements of analytes in complex matrices with dynamic linear ranges that span up to 4-5 orders of magnitude. Hence, targeted methods are promising for the development of robust protein assays in several sensitive areas, e.g. in health care. However, exploiting the full method potential requires reliable determination of the dynamic range along with related quantification limits for each analyte.
Here, we present a software named CalibraCurve that enables an automated batch-mode determination of dynamic linear ranges and quantification limits for both targeted proteomics and similar assays. The software uses a variety of measures to assess the accuracy of the calibration, namely precision and trueness. Two different kinds of customizable graphs are created (calibration curves and response factor plots). The accuracy measures and the graphs offer an intuitive, detailed and reliable opportunity to assess the quality of the model fit.

## Installation

You can install the development version of CalibraCurve from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mpc-bioinformatics/CalibraCurve")
library(CalibraCurve)
```

New features will be introduced in github branches, which are merged as soon as the feature has been properly tested.
To install a specific branch of the package please use

``` r
# install.packages("devtools")
devtools::install_github("mpc-bioinformatics/CalibraCurve", ref = "<branchname>")
library(CalibraCurve)
```

## Implementation

CalibraCurve is written in R (currently as an R script, the plan is to ). Furthermore, a nextflow-workflow running this script is available. 

CalibraCurve is freely available under the 3-clause BSD license.
The download also comprises a detailed manual (which however, only covers version 2.0 at this point) and example data along with corresponding example result files.

### Version 3.0

Currently, version 3.0 is under development, including a nextflow workflow that generates interactive graphics by using plotly. 

### Version 2.0

The stable version 2.0 can be downloaded as a release https://github.com/mpc-bioinformatics/CalibraCurve/releases/tag/v_2_0 . It includes a KNIME-workflow, that is replaced by a nextflow workflow in newer versions. The KNIME-workflow will not be maintained further.

## Usage

TODO!


## Publication

Kohl M, Stepath M, Bracht T, Megger DA, Sitek B, Marcus K, Eisenacher M. CalibraCurve: A Tool for Calibration of Targeted MS-Based Measurements. Proteomics. 2020 Jun;20(11):e1900143. doi: 10.1002/pmic.201900143. Epub 2020 Mar 6. PMID: 32086983.

## Funding

The development and maintanence of CalibraCurve is funded by de.NBI (https://www.denbi.de/) and CUBiMed.RUB (https://www.cubimed.ruhr-uni-bochum.de/index.html.en).
We offer also other cool tools and consulting for statistics, bioninformatics and machine learning!

## Feedback

Please fill out the following survey to give feedback:

https://de.surveymonkey.com/r/denbi-service?sc=bioinfra-prot&tool=calibracurve




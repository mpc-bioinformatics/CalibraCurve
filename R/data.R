#' Targeted Proteomics data set Albumin
#'
#' Measurement of the y8 fragment ion of the  peptide LVNEVTEFAK of Albumin (Mohammed et al., 2015).
#' TODO
#'
#' @format ## `D_ALB`
#' A data frame with 21 rows and 2 columns:
#' \describe{
#'   \item{Concentration}{Concentration in fmol/µl}
#'   \item{Measurement}{Measurement as ratio between response of stable isotope standard and native peptide}
#' }
#' @source Y. Mohammed, A. J. Percy, A. G. Chambers, C. H. Borchers, J. Proteome Res. 2015, 14, 1137.
"D_ALB"



#' Targeted Proteomics data set Apolipoprotein
#'
#' Measurement of the y7 fragment ion of the  peptide TPAYYPNAGLIK of Apolipoprotein (Mohammed et al., 2015).
#' TODO
#'
#' @format ## `D_Apolipoprotein`
#' A data frame with 21 rows and 2 columns:
#' \describe{
#'   \item{Concentration}{Concentration in fmol/µl}
#'   \item{Measurement}{Measurement as ratio between response of stable isotope standard and native peptide}
#' }
#' @source Y. Mohammed, A. J. Percy, A. G. Chambers, C. H. Borchers, J. Proteome Res. 2015, 14, 1137.
"D_Apolipoprotein"


#' Targeted Proteomics data set MFAP4
#'
#' MRM validation experiment carried out for the y4 fragment ion of peptide WTVFQK of the microfibril associated protein 4 (MFAP4).
#' The measurements were conducted as part of a proteomics discovery study aiming at biomarker identification for hepatic fibrosis (Bracht et al., 2015).
#' The raw files can be obtained from PASSEL using the identifier PASS00653 (TODO: Link).
#'
#' @format ## `D_MFAP4`
#' A data frame with 55 rows and 2 columns:
#' \describe{
#'   \item{Concentration}{Concentration in fmol/µl}
#'   \item{Measurement}{Measurement  as total area}
#' }
#' @source T. Bracht, V. Schweinsberg, M. Trippler, M. Kohl, M. Ahrens, J. Padden,
#'W. Naboulsi, K. Barkovits, D. A. Megger, M. Eisenacher, C. H.
#'Borchers, J. F. Schlaak, A. C. Hoffmann, F. Weber, H. A. Baba, H. E.
#'Meyer, B. Sitek, J. Proteome Res. 2015, 14, 2278.
"D_MFAP4"



#' CalibraCurve results for data set Albumin
#'
#' CalibraCurve results for the data set \code{\link{D_ALB}}. Default parameters for CalibraCurve were used.
#'
#' @format ## `RES_ALB`
#' Results from CalibraCurve calculations.
"RES_ALB"


#' CalibraCurve results for data set Apolipoprotein
#'
#' CalibraCurve results for the data set \code{\link{D_Apolipoprotein}}. Default parameters for CalibraCurve were used.
#'
#' @format ## `RES_Apolipoprotein`
#' Results from CalibraCurve calculations.
"RES_Apolipoprotein"


#' CalibraCurve results for data set MFAP4
#'
#' CalibraCurve results for the data set \code{\link{D_MFAP4}}. Default parameters for CalibraCurve were used.
#'
#' @format ## `RES_MFAP4`
#' Results from CalibraCurve calculations.
"RES_MFAP4"



# As test data we use data from the msqc1 R package. More precisely, we use the
# dilution series stored in the msqc1_dil object.

# We filter this data in the following way:
# - using only QTRAP, TSQVantage and QExactive instruments (PRM or SRM method)
# - use only data on the level of y-ions, no precursors
# - use only the heavy isotope version of the peptide

# We generate a SummarizedExperiment object for each Peptide sequence and store
# them as .rds files.

# Additionally, for the peptide sequence "GGPFSDSYR" we store the data as xlsx
# tables, separated by instrument and ion type.



################################################################################
################################################################################
### SummarizedExperiment

library(msqc1)
library(tidyr)
library(dplyr)
library(openxlsx)
library(SummarizedExperiment)
data(msqc1_dil)

### from supplement of paper:
### Info about the heavy peptide amount (in fmol) in the samples with relative.amount == 1
peptide_amounts <- c(
  "ALIVLAHSER" = 100,
  "AVQQPDGLAVLGIFLK" = 100,
  "EGHLSPDIVAEQK" = 200,
  "ESDTSYVSLK" = 20,
  "FEDENFILK" = 80,
  "FSTVAGESGSADTVR" = 4,
  "GAGAFGYFEVTHDITK" = 200,
  "GGPFSDSYR" = 1000,
  "GYSIFSYATK" = 4,
  "NLSVEDAAR" = 0.8,
  "SADFTNFDPR" = 20,
  "TAENFR" = 20,
  "VLDALQAIK" = 500,
  "VSFELFADK" = 40
)



instruments <- c("QTRAP", "TSQVantage", "QExactive") # unique(msqc1_dil$instrument)
ions <- c("y10", "y11", "y12", "y4", "y5", "y6", "y7", "y8", "y9")

D <- filter(msqc1_dil, instrument %in% instruments,
            Protein.Name != "iRT-C18 Standard Peptides",
            Fragment.Ion %in% ions,
            Isotope.Label.Type == "heavy"
            )
D$Replicate.Name <- droplevels(D$Replicate.Name)
D$File.Name <- droplevels(D$File.Name)
D$Protein.Name <- droplevels(D$Protein.Name)
D$Peptide.Sequence <- droplevels(D$Peptide.Sequence)
D$Fragment.Ion <- droplevels(D$Fragment.Ion)
D$instrument <- D$instrument

replicatename_split <- limma::strsplit2(D$Replicate.Name, "_")
replicate <- replicatename_split[,6]
replicate[replicate == ""] <- replicatename_split[replicate == "",5]
D$replicate <- replicate
### correct a typo in the replicate name:
D$replicate[D$Replicate.Name == "20140818_004_MSQC1_1_40dil_1"] <- "2"


################
## extract relevant data and save them as SummarizedExperiments-object in rds files (one per peptide sequence)

peptides <- unique(D$Peptide.Sequence)

for (i in seq_along(peptides)) {

  peptide <- peptides[i]

  D_tmp <- dplyr::filter(D, Peptide.Sequence == peptide)
  D_tmp$amount <- peptide_amounts[peptide] * D_tmp$relative.amount
  D_tmp$amount_replicate <- paste(D_tmp$amount, D_tmp$replicate, sep = "_")


  D_tmp_wide <- tidyr::pivot_wider(D_tmp,
                            id_cols = c("instrument", "Fragment.Ion", "Isotope.Label.Type"),
                            names_from = amount_replicate,
                            values_from = Area)

  rowData <- D_tmp_wide[, 1:3]
  rowData$Substance <- paste(rowData$instrument, rowData$Fragment.Ion, sep = "_")

  colData <- limma::strsplit2(colnames(D_tmp_wide)[-c(1:3)], "_")
  colnames(colData) <- c("amount_fmol", "replicate")

  D_SE <- SummarizedExperiment(assays=list(Area=as.data.frame(D_tmp_wide[,-c(1:3)])),
                               rowData=rowData, colData=colData)
  SummarizedExperiment::metadata(D_SE) <- list(peptide = peptide)

  saveRDS(D_SE, file = paste0("inst/extdata/MSQC1/msqc1_dil_", peptide, ".rds"))

}



#######################
## for one Peptide, save data as xlsx files, separately for each ion and instrument

D <- filter(msqc1_dil,
            instrument %in% instruments,
            Peptide.Sequence == "GGPFSDSYR",
            Fragment.Ion %in% ions,
            Isotope.Label.Type == "heavy"
)

ions_GGPFSDSYR <- unique(D$Fragment.Ion)

for (inst in instruments) {
  for (ion in ions_GGPFSDSYR) {
    D_tmp <- filter(D, instrument == inst, Fragment.Ion == ion)
    D_tmp$amount <- peptide_amounts["GGPFSDSYR"] * D_tmp$relative.amount
    openxlsx::write.xlsx(D_tmp, paste0("inst/extdata/MSQC1_xlsx/GGPFSDSYR_", inst, "_", ion, ".xlsx"))
  }
}








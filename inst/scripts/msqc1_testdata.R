
library(tidyverse)
library(openxlsx)

instruments <- c("QTRAP", "TSQVantage", "QExactive") # unique(msqc1_dil$instrument)
peptides <- unique(msqc1_dil$Peptide.Sequence)
isotope <- c("light", "heavy")

### remove iRT peptides
D <- filter(msqc1_dil, Protein.Name != "iRT-C18 Standard Peptides")

for (instr in instruments) {
  for (peptide in peptides) {
    D_tmp <- dplyr::filter(D,
                           instrument == instr,
                           Peptide.Sequence == peptide)
    fragments_tmp <- as.character(unique(D_tmp$Fragment.Ion))
    print(fragments_tmp)
    ### remove data from the precursors:
    fragments_tmp <- fragments_tmp[!(fragments_tmp %in% c("precursor", "precursor [M+1]", "precursor [M+2]"))]

    #print(paste(instrument, peptide))
    #print(fragments_tmp)
    for (fragment in fragments_tmp) {
      for (isotope_type in isotope) {
        print(fragment)
        print(isotope_type)
        D_tmp2 <- dplyr::filter(D_tmp,
                                Fragment.Ion == fragment,
                                Isotope.Label.Type == isotope_type)

        if (nrow(D_tmp2) > 0) {
          file_name <- paste0("../example_data/MSQC1/msqc1_dil_", instr, "_", peptide, "_", fragment, "_", isotope_type, ".xlsx")
          write.xlsx(D_tmp2, file = file_name, rowNames = FALSE)
        }
      }
    }


  }
}


instr = "QExactive"
peptide = "ALIVLAHSER"



library(CalibraCurve)
CalibraCurve(data_folder = "../example_data/MSQC1/",
             output_path = "../example_data/MSQC1/results/",
             conc_col = 14,
             meas_col = 12,
             plot_type = "single_plots")



### TODO:response factor plot komisch???
library(CalibraCurve)
CalibraCurve(data_path = "../example_data/MSQC1/msqc1_dil_QExactive_AVQQPDGLAVLGIFLK_y10_light.xlsx",
             output_path = "../example_data/MSQC1/results/",
             conc_col = 14,
             meas_col = 12,
             plot_type = "single_plots")





################################################################################
################################################################################
### SummarizedExperiment

library(msqc1)
library(tidyverse)
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


#D$amount_replicate <- paste(D$relative.amount, D$replicate, sep = "_")


################

peptides <- unique(D$Peptide.Sequence)

for (i in seq_along(peptides)) {

  peptide <- peptides[i]

  D_tmp <- filter(D, Peptide.Sequence == peptide)
  D_tmp$amount <- peptide_amounts[peptide] * D_tmp$relative.amount
  D_tmp$amount_replicate <- paste(D_tmp$amount, D_tmp$replicate, sep = "_")


  D_tmp_wide <- pivot_wider(D_tmp,
                            id_cols = c("instrument", "Fragment.Ion", "Isotope.Label.Type"),
                            names_from = amount_replicate,
                            values_from = Area)

  rowData <- D_tmp_wide[, 1:3]
  rowData$Substance <- paste(rowData$instrument, rowData$Fragment.Ion, sep = "_")

  colData <- limma::strsplit2(colnames(D_tmp_wide)[-c(1:3)], "_")
  colnames(colData) <- c("amount_fmol", "replicate")

  D_SE <- SummarizedExperiment(assays=list(Area=as.data.frame(D_tmp_wide[,-c(1:3)])),
                               rowData=rowData, colData=colData)
  metadata(D_SE) <- list(peptide = peptide)

  saveRDS(D_SE, file = paste0("inst/extdata/MSQC1/msqc1_dil_", peptide, ".rds"))

}



DATA <- readRDS("inst/extdata/MSQC1/msqc1_dil_ALIVLAHSER.rds")

assays(DATA)$Area
rowData(DATA)
colData(DATA)






################################################################################
################################################################################
################################################################################



# ALIVLAHSER: QExactive hat zu hohe CVs, TSQ kleine linear ranges, QTRAP ganz ok
#X AVQQPDGLAVLGIFLK: QExactive teils kleine ranges, TSQ y9 geht nicht, QTRAp kleine ranges
# EGHLSPDIVAEQK: QExactive geht nicht, TSQ ok, QTRAP ok
# ESDTSYVSLK QExactive y9 geht nicht, TSQ teils kleinen range, QTRAP teils kleinen range
# FEDENFILK: sieht alles gut aus, fast perfekte ranges
# FSTVAGESGSADTVR: QExactive geht nicht, sonst alles ok
#X GAGAFGYFEVTHDITK: TSQ ist komisch
# GGPFSDSYR: sieht alles sehr gut aus!
#X GYSIFSYATK: QTRAP teils etwas komisch
#X NLSVEDAAR: QExactive geht kaputt, sonst sehr kleine ranges
#X SADFTNFDPR:
#X TAENFR
# VLDALQAIK: sieht alles ok aus
# VSFELFADK

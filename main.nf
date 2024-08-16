#!/usr/bin/env nextflow
include { runCalibraCurve as cc } from './modules/calibraCurve'
include { runCalibraCurve as ccXlsx } from './modules/calibraCurve'

// Convert given file from Excel Sheet to CSV
process convertFromXlsxToCsv {
    container 'mpc/calibracurve-r:1.0.0'

    input:
    path input_file
    output:
    path "*.csv"

    """
    #!/usr/bin/env Rscript
    require(openxlsx)

    # We read only first sheet. Use closure {-> ..} to fix escaping of special symbols in name
    data <- read.xlsx("${-> input_file}", 1, sep.names = " ")

    # Output the file with the same name but with `.csv` extension
    output_file <- sub(".xlsx", ".csv", "${-> input_file}")
    
    # We use `write.csv2` function to get ";" as column separators   ### TODO: is this smart? Better use international csv standard with comma as separator?
    write.csv2(data, file = output_file, row.names = FALSE, quote = F)
    """
}

workflow {
    println("The input dir is ${params.inputDir}")
    println("The output dir is ${params.outputDir}")

    // Run Calibra Curve for CSV files
    channel.fromPath("${params.inputDir}/*.csv")
            | cc

    // Run Calibra Curve for XLSX files
    channel.fromPath("${params.inputDir}/*.xlsx")
            | convertFromXlsxToCsv
            | ccXlsx
}

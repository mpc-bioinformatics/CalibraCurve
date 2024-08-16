/**
 * Module `Calibra Curve` introduces the process `runCalibraCurveNaive`, which can be included in the parent
 * workflow and be reused. <br>
 *
 * For example: if converting the input from different types to CSV, one can reuse this
 * module and run the same process.
 */
process runCalibraCurve {
    container 'mpc/calibracurve-r:1.0.0'
    publishDir "${params.outputDir}/raw", mode: "copy"
    publishDir "${params.outputDir}/rfplots-json", mode: "copy", pattern: "RFplot_*.json"
    publishDir "${params.outputDir}/lmplots-json", mode: "copy", pattern: "LMplot_*.json"
    publishDir "${params.outputDir}/rfplots", mode: "copy", pattern: "RFplot_*.png"
    publishDir "${params.outputDir}/lmplots", mode: "copy", pattern: "LMplot_*.png"

    input:
    env MPC_CC_INPUT_FILE
    output:
    path "*.html"
    path "*.json"
    path "*.png"
    path "*.txt"

    // using Rscript declared in the `bin` dir of the `nextflow` directory
    """
    CalibraCurve_v3.0.R
    """
}

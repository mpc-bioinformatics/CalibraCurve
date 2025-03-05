# Writing result files: Preparation of result files
preparationResultFiles <- function(x, resF, resRegF){
  write(paste(resFHeader, x, sep = ""),
        file = resF, append = FALSE)
  sink(resRegF, append = FALSE)
  write(paste(resRegFHeader, x, sep = ""),
        file = resRegF, append = FALSE)
  sink()
}

# Writing result files: CV values for the preliminar linear range
writeCV_forPreliminaryLinearRange  <- function(x, resF, cUnit){
  conc <- names(x)
  write("", resF, append=TRUE)
  write(resFCVHeader, resF, append=TRUE)
  for(i in seq_along(x)){
    textCV <- paste(resFCVPrefix,
                    conc[i], cUnit, ": ",round(x[i], digits = numberDecimals),sep="")
    write(textCV, resF, append=TRUE)
  }
}

# Writing result files: Function, which writes percent bias and response factor values for the final linear range
# or for writing of intermediate percent bias calculations into the logging file (CalibraCurve verbose mode)
writePerBiasRfFinalLinearRange <- function(x, perBiasWeighted, perBiasUnweighted, Rf, resF, colDelim, noDec){
  for(i in seq_along(dataFinal)){
    PerBias_no_weight<- perBiasUnweighted[[i]]
    PerBias_weight <- perBiasWeighted[[i]]
    ResponseFactors <- Rf[[i]]
    x[[i]] <- cbind(x[[i]], PerBias_no_weight, PerBias_weight, ResponseFactors)
  }
  xDF <- do.call(rbind, x)
  header <- paste(colnames(xDF), collapse = colDelim)
  write(header, resF, append=TRUE)
  xDF <- format(xDF, digits = noDec, nsmall = noDec)
  write.table(xDF, resF, append=TRUE, row.names=FALSE, col.names=FALSE, sep=colDelim)
}

# Writing result files: Function, which writes overview information (central tendency and dispersion
# measures) for the percent bias values of each concentration level
writePerBiasOverviewFinalLinearRange  <- function(x, resF, colDelim, message, noDec, perBiasT){
  indices <- unique(sort(which(x$avgPerBias >= perBiasT), decreasing = FALSE))
  warnings <- NULL
  for(i in seq_along(x$avgPerBias)){
    warnings <- c(warnings, "")
  }
  warnings[indices] <- message
  x <- cbind(x, warnings)
  header <- paste(colnames(x), collapse = colDelim)
  header <- paste("expectedConc",header, sep = colDelim)
  if(length(indices)>0){
    mess <- paste("Input data file ",currFileNameNoExt[1],": ",genericWarningOnShell2, sep="")
    print(mess)
    write(mess, logFile, append=TRUE)
    write(warnMessAcc, resF, append=TRUE)
  }
  write(header, resF, append=TRUE)
  x <- format(x, digits = noDec, nsmall = noDec)
  write.table(x, resF, append=TRUE, row.names=TRUE, col.names=FALSE, sep=colDelim)
}

# Writing result files: Function, which writes mean response factor values for the
# concentration levels of the final data range
writeMeanRFvalues <- function(x, resF, colDelim, noDec){
  c <- paste(format (as.numeric(names(x)), digits = noDec, nsmall = noDec), collapse = colDelim)
  c <- paste("Expected concentration", c, sep = colDelim)
  write(c, resF, append=TRUE)
  x <- format(x, digits = noDec, nsmall = noDec)
  rf <- paste(x, collapse = colDelim)
  rf <- paste("Mean response factor", rf, sep = colDelim)
  write(rf, resF, append=TRUE)
}

# Writing result files: Function, which writes summary information for the fitted linear model into the second result file
writeLMSummary  <- function(x, weightedLM, resRegF, weightMet){
  # unweighted model
  write("Summary for the unweighted linear model:", resRegF, append = TRUE)
  write("------------------------------------------", resRegF, append = TRUE)
  sink(resRegF, append = TRUE)
  print(summary(x))
  sink()
  # Weighted model
  header <- paste("Summary for the weighted linear model (Weighting method:",weightMet,"):", sep = "")
  write(header, resRegF, append = TRUE)
  write("------------------------------------------", resRegF, append = TRUE)
  sink(resRegF, append = TRUE)
  print(summary(weightedLM))
  sink()
}

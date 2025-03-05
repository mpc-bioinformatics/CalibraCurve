# Plotting
generatePlots <- function(x, titlePart1, sampleName, titlePart2, LLOQ, ULOQ, unit,
                          x_lab, y_lab, lmW, lmNoW, plotOpt, scale,
                          colNumConc, colNumMeas){
  tiff(paste(sampleName, ".tif", sep = ""), width = plotWidth, height = plotHeight, units = "cm", res = plotResolution)
  par(cex.axis = magnificAxis, cex.lab = magnificLabels, cex.main = magnificTitle)
  if(showLegend){
    par(xpd=T, mar=par()$mar+c(6,0,0,0))
  }
  title <- paste(titlePart1, sampleName, '\n',
                 titlePart2, format(as.numeric(LLOQ), digits=numberDecimals, nsmall =numberDecimals),"-",
                 format(as.numeric(ULOQ), digits=numberDecimals, nsmall =numberDecimals), " ",unit, sep = "")
  if(scale == 'log'){
    plot(x[,colNumConc], x[,colNumMeas], main = title,
         xlab = x_lab, ylab = y_lab, log = "xy")
  }else{
    plot(x[,colNumConc], x[,colNumMeas], main = title,
         xlab = x_lab, ylab = y_lab)
  }

  r2NOW=summary(lmNoW)$r.squared
  r2W=summary(lmW)$r.squared
  # Unweighted regression equation
  if (lmNoW$coefficients[1]>0){
    eq <- paste("y=", format(lmNoW$coefficients[2], digits=2, scientific=TRUE),"*x+",
                format(lmNoW$coefficients[1], digits=2, scientific=TRUE),",",
                " R^2=",format(r2NOW, digits=3), sep="")
  }else{
    eq <- paste("y=", format(lmNoW$coefficients[2], digits=2, scientific=TRUE),"*x",
                format(lmNoW$coefficients[1], digits=2, scientific=TRUE),",",
                " R^2=",format(r2NOW, digits=3), sep="")
  }
  # Weighted regression equation
  if (lmW$coefficients[1]>0){
    eq_w <- paste("y=", format(lmW$coefficients[2], digits=2, scientific=TRUE),"*x+",
                  format(lmW$coefficients[1], digits=2, scientific=TRUE),",",
                  " R^2=",format(r2W, digits=3), sep="")
  }else{
    eq_w <- paste("y=", format(lmW$coefficients[2], digits=2, scientific=TRUE),"*x",
                  format(lmW$coefficients[1], digits=2, scientific=TRUE),",",
                  " R^2=",format(r2W, digits=3), sep="")
  }
  if(plotOpt == "both"){ # including both weighted and unweighted regression lines and corresponding equations into the plot
    # Regression lines
    abline(lmNoW, col="blue", untf=TRUE, xpd=FALSE)
    abline(lmW, col="red", untf=TRUE, xpd=FALSE)
    # Equations
    # Unweighted equation: top left
    text(min(x[,colNumConc]), max(x[,colNumMeas]), eq, adj = c( 0, 1 ), col = "blue", cex = magnificEquation)
    if (scale == 'log'){
      # Weighted equation: below the unweighted equation
      text(min(x[,colNumConc]), max(x[,colNumMeas])/2, eq_w, adj = c( 0, 1 ), col = "red", cex = magnificEquation)
    }else{
      text(max(x[,colNumConc])/2, max(x[,colNumMeas]), eq_w, adj = c( 0, 1 ), col = "red", cex = magnificEquation)
    }
    # Legend
    if(showLegend){
      legend("bottom", # location 1: centered at the bottom
             legend=c(legend_weighted, legend_unweighted),
             inset=c(0,-0.45), # location 2: vertical position
             col=c("red","blue"),
             lty = c(1,1))
    }

  }
  if(plotOpt == "regW"){ # including only the weighted regression line into the plot
    abline(lmW, col="red", untf=TRUE, xpd=FALSE)
    # Weighted equation: top left
    text(min(x[,colNumConc]), max(x[,colNumMeas]), eq_w, adj = c( 0, 1 ), col = "red", cex = magnificEquation )
    if(showLegend){
      legend("bottom", # location 1: centered at the bottom
             legend=c(legend_weighted),
             inset=c(0,-0.45), # location 2: vertical position
             col=c("red"),
             lty = c(1))
    }
  }
  if(plotOpt == "regUW"){ # including only unweighted regression line into the plot
    abline(lmNoW, col="blue", untf=TRUE, xpd=FALSE)
    # Unweighted equation: top left
    text(min(x[,colNumConc]), max(x[,colNumMeas]), eq, adj = c( 0, 1 ), col = "blue", cex = magnificEquation )
    if(showLegend){
      legend("bottom", # location 1: centered at the bottom
             legend=c(legend_unweighted),
             inset=c(0,-0.45), # location 2: vertical position
             col=c("blue"),
             lty = c(1))
    }
  }
  dev.off()
}

# Plotting: The function generates a response factor plot (Response factors vs. concentration) and presents
# both response factors for each measurement and mean response factors for each concentration level
# The plot also shows lines indicating user defined thresholds for percent deviation from the mean
# response factor that is calculated from the data of the final linear range.
generateRfPlot <- function(x, RfDataValid,  RfTLowF, RfTUpperF, avgRFDataF, avgRFDataV, restrictRF,
                           adaptYlim, minYlim, maxYlim,
                           scaleRFX, magnAxisRF, magnLabelsRF,
                           magnTitleRF, plotResRF, plotHRF, plotWRF,
                           titleOpenRF, filename, xLRF, yLRF){
  # Calculation of the mean RF from the RF of data final
  allRf <- NULL
  for (i in seq_along(x)){
    currRf <- x[[i]]
    allRf <- c(allRf, currRf)
  }
  meanRfDF <- mean(allRf)
  # Calculating values for the horizontal threshold lines
  hLineLow <- meanRfDF*RfTLowF
  hLineUpper <- meanRfDF*RfTUpperF
  # Prepare data for plotting from the validated data
  plotReadyDataValid <- NULL
  for (i in seq_along(RfDataValid)){
    plotReadyDataValid <- rbind(plotReadyDataValid, calculateRfDf(RfDataValid[i]))
  }
  # Prepare data for plotting from the final data
  plotReadyDataFinal <- NULL
  for (i in seq_along(x)){
    plotReadyDataFinal <- rbind(plotReadyDataFinal, calculateRfDf(x[i]))
  }
  # Calculation of ylim properties
  if(restrictRF == TRUE){
    ylimMin <- min(plotReadyDataFinal$RF)
    ylimMax <- max(plotReadyDataFinal$RF)
  }else{
    ylimMin <- min(plotReadyDataValid$RF, plotReadyDataFinal$RF)
    ylimMax <- max(plotReadyDataValid$RF, plotReadyDataFinal$RF)
  }
  if(hLineLow < ylimMin){
    ylimMin = hLineLow
  }
  if(hLineUpper > ylimMax){
    ylimMax = hLineUpper
  }
  if(adaptYlim){# Enables user defined adaption of the y axis
    ylimMin <- minYlim
    ylimMax <- maxYlim
  }
  rfPlotName <- paste("RFplot_", filename, ".tif", sep="")
  tiff(rfPlotName, width = plotWRF, height = plotHRF, units = 'cm', res = plotResRF)
  par(cex.axis = magnAxisRF, cex.lab = magnLabelsRF, cex.main = magnTitleRF)
  if(showLegendRF){
    par(xpd=T, mar=par()$mar+c(6,0,0,0))
  }
  # Case: Plot response factors for the complete input data
  title <- paste(titleOpenRF, filename, sep="")
  if (restrictRF == FALSE){
    if (scaleRFX == "log"){
      plot(plotReadyDataValid, log = "x", col = 'blue', ylim = c(ylimMin, ylimMax),
           main = title, xlab = xLRF, ylab = yLRF)
    }else{
      plot(plotReadyDataValid, col = 'blue', ylim = c(ylimMin, ylimMax),
           main = title, xlab = xLRF, ylab = yLRF)
    }
  }else{
    # Recalculation of ylim
    if (scaleRFX == "log"){
      plot(plotReadyDataFinal, log = "x", col = 'red', ylim = c(ylimMin, ylimMax),
           main = title, xlab = xLRF, ylab = yLRF)
    }else{
      plot(plotReadyDataFinal, col = 'red', ylim = c(ylimMin, ylimMax),
           main = title, xlab = xLRF, ylab = yLRF)
    }
  }
  if (restrictRF == FALSE){
    points(plotReadyDataFinal$Concentration, plotReadyDataFinal$RF, col="red")
    lines(names(avgRFDataV), avgRFDataV, col="blue",lty=1, type = "o", pch = 19)
  }
  lines(names(avgRFDataF), avgRFDataF, col="red",lty=1, type = "o", pch = 19)
  abline(h = hLineLow, col = 'green', lty = 2, xpd=FALSE)
  abline(h = hLineUpper, col = 'green', lty = 2, xpd=FALSE)
  if(showLegendRF){
    legend("bottom", # location 1: centered at the bottom
           legend=c(legend_RF),
           inset=c(0,-0.45), # location 2: vertical position
           col=c("blue", "red"),
           lty = c(1,1))
  }
  dev.off()
}

# Auxillary function used for response factor plotting:
# The function generates a data frame with a concentration column and a response factor column
calculateRfDf <- function(x){
  Concentration <- NULL
  RF <- NULL
  conc <- as.numeric(names(x))
  for (i in x){
    Concentration <- c(Concentration, conc)
    RF <- c(RF, i)
  }
  mDF <- data.frame(Concentration, RF)
  return(mDF)
}

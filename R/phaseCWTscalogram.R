#' Generate averaged scalograms from Continuous Wavelet Transforms
#'
#' @description
#' This function generates a scalogram averaged over flies. The heatmap plots amplitude across the time-frequency domain. Input for this function must be an output from the trimData() function. The output of this function is a large plotly object. This function requires the packages "plotly" and "WaveletComp".
#'
#' @param input Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param data A choice between "Activity" and "Sleep" data to analyze using the wavelet transform.
#' @param out.bin Define the desired output bin size (in minutes). This defaults to 15.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param n.days Number of days to analyse.
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep or specific bout durations; sleep.def allows users to change the definition of sleep. The default input is a single value vector of value 5. If users wish to analyse sleep only between 5 to 20 mins, the input must be c(5,20). Only relevant if data = "Sleep".
#' @param low.per Choose the lowest period (in hours) for analysis. This defaults to 1.
#' @param high.per Choose the highest period (in hours) for analysis. This defaults to 35.
#' @param edge.rm A logical (TRUE or FALSE). TRUE will plot the scalogram with the edge effects removed. FALSE will plot the scalogram with the regions of the edge effect still retained in the plot. This defaults to FALSE.
#' @param make.pval A logical (TRUE or FALSE). TRUE will generate p-values for the significance of detected periodicities. FALSE will not. Note that generating p-values significantly increases computation time. Defaults to FALSE.
#' @param n.sim If TRUE, the number of simulations based on which p-values will be computed. Defaults to 1000.
#' @param method The method used for generating surrogate time-series against which p-values are generated. Also see help(analyze.wavelet). Available choices are "white.noise": white noise, "shuffle": shuffling the given time-series, "Fourier.rand": time-series with a similar spectrum, "AR": AR(p), "ARIMA": ARIMA(p,0,q). Defaults to "shuffle".
#' @param z.max The maximum color on the scalogram. Defaults to 0.5.
#'
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' @importFrom grDevices rgb colorRamp
#' @importFrom WaveletComp analyze.wavelet
#' 
#' @return A \code{plotly} \code{htmlwidget}.
#'
#' @export phaseCWTscalogram
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 7, bin = 1, t.cycle = 24)
#' scalogram <- phaseCWTscalogram(input = td, n.days = 3, low.per = 1, high.per = 3)

phaseCWTscalogram <- function (input, data = "Activity", out.bin = 15, t.cycle = 24, n.days, sleep.def = 5, low.per = 1, high.per = 35, edge.rm = FALSE, make.pval = FALSE, n.sim = 1000, method = "shuffle", z.max = 0.5) {
  
  requireNamespace("plotly")
  requireNamespace("WaveletComp")
  
  if (requireNamespace("plotly", quietly = T)) {
    s.per.hr = 60/out.bin
    f1 <- list( # define font details
      family = "Arial, sans-serif",
      size = 16,
      color = "black"
    )
    f2 <- list( # define font details
      family = "Arial, sans-serif",
      size = 12,
      color = "black"
    )
    
    if (data == "Activity") {
      df <- binData(data = input, input.bin = 1, output.bin = out.bin, t.cycle = t.cycle)
    } else if (data == "Sleep") {
      df <- sleepData(data = input, sleep.def = sleep.def, bin = out.bin, t.cycle = t.cycle)
    }
    
    cwt.rds <- list()
    
    if (isTRUE(make.pval)) {
      for (i in 1:(length(df[1,])-1)) {
        if (sum(df[,i+1]) == (out.bin * length(df[,1])) || sum(df[,i+1]) == 0) {
          cwt.rds[[i]] <- NA # check if fly-i is dead, and if it is then assign a value of NA to the CWT for fly-i
        } else {
          cwt.rds[[i]] <- WaveletComp::analyze.wavelet(my.data = df, my.series = paste("I", i, sep = ""),
                                                       loess.span = 0, dt = 1, lowerPeriod = low.per * s.per.hr,
                                                       upperPeriod = high.per * s.per.hr, dj = 1/100,
                                                       make.pval = T, verbose = F, n.sim = n.sim, method = method)
        } # if fly-i is alive, compute the CWT and store value.
      } # for loop to compute CWT for each fly
    } else {
      for (i in 1:(length(df[1,])-1)) {
        if (sum(df[,i+1]) == (out.bin * length(df[,1])) || sum(df[,i+1]) == 0) {
          cwt.rds[[i]] <- NA # check if fly-i is dead, and if it is then assign a value of NA to the CWT for fly-i
        } else {
          cwt.rds[[i]] <- WaveletComp::analyze.wavelet(my.data = df, my.series = paste("I", i, sep = ""),
                                                       loess.span = 0, dt = 1, lowerPeriod = low.per * s.per.hr,
                                                       upperPeriod = high.per * s.per.hr, dj = 1/100,
                                                       make.pval = F, verbose = F, n.sim = n.sim, method = method)
        } # if fly-i is alive, compute the CWT and store value.
      } # for loop to compute CWT for each fly
    }
    
    names(cwt.rds) <- paste("Ch-", 1:32, sep = "")
    
    cwt.no.na <- cwt.rds[which(!is.na(cwt.rds))]
    
    power.mat <- list() # create an empty list for storing power matrices
    power.norm <- list() # create an empty list for storing normalized power matrices
    
    if (edge.rm == FALSE) {
      for (i in 1:length(cwt.no.na)) { # NO edge effect removal for fly-i
        pow <- cwt.no.na[[i]]$Power # storing fly-i's power matrix in a new variable
        power.mat[[i]] <- pow # storing the power matrices in a new list
      }
      
      for (i in 1:length(cwt.no.na)) { # normalization routine for fly-i
        avg.power <- mean(as.vector(power.mat[[i]]), na.rm = T) # calculate average power for fly-i
        power.norm[[i]] <- power.mat[[i]]/avg.power # divide power matrix for fly-i by average power and store in a new variable
      }
      
      power.array <- do.call(cbind, power.norm) # convert the list of normalized power to an array
      power.array <- array(power.array, dim = c(dim(power.norm[[1]]), length(power.norm))) # convert the list of normalized power to an array
      power.mean <- apply(power.array, c(1, 2), mean, na.rm = TRUE) # average across flies
      
      p <- plotly::plot_ly(
        # plot scalograms
      )%>%
        add_trace(
          x = 1:dim(power.mean)[2],
          y = cwt.no.na[[1]]$Period,
          z = power.mean,
          type = "heatmap",
          colors = colorRamp(c("blue", "cyan", "green", "yellow", "red")),
          zauto = F,
          zmin = 0,
          zmax = z.max,
          showscale = T
        )%>%
        layout(
          showlegend = F,
          yaxis = list(
            type = "log",
            showgrid = F,
            showline = T,
            titlefont = f1,
            tickfont = f2,
            title = "Period (h)",
            linecolor = "black",
            autotick = FALSE,
            ticks = "outside",
            tick0 = 0,
            dtick = 1 * s.per.hr,
            ticklen = 7,
            tickcolor = "black",
            range = c(log10((low.per) * s.per.hr), log10((high.per) * s.per.hr)),
            ticktext = as.list(
              c(1,2,4,8,12,17,24,35)
            ),
            tickvals = as.list(
              c(1,2,4,8,12,17,24,35) * s.per.hr
            ),
            tickmode = "array"
          ),
          xaxis = list(
            showgrid = F,
            showline = T,
            titlefont = f1,
            tickfont = f2,
            title = "Days",
            linecolor = "black",
            autotick = FALSE,
            ticks = "outside",
            tick0 = 0,
            dtick = (s.per.hr * 24)/2,
            ticklen = 7,
            tickcolor = "black",
            range = c(0, dim(power.mean)[2]),
            ticktext = as.list(
              seq(0.5, n.days, by = 0.5)
            ),
            tickvals = as.list(
              seq((s.per.hr * 24)/2, s.per.hr * 24 * n.days, (s.per.hr * 24)/2)
            ),
            tickmode = "array"
          )
        )
      return(p)
      
    } else if (edge.rm == TRUE) {
      for (i in 1:length(cwt.no.na)) { # edge effect removal for fly-i
        pow <- cwt.no.na[[i]]$Power # storing fly-i's power matrix in a new variable
        for (j in 1:length(pow[,1])) { # edge effect removal for each period-j for fly-i
          idx2rem <- round((cwt.no.na[[i]]$Period[j]) * 1.5) # calculating the number of columns to assign NA due to edge effects
          pow[j,1:idx2rem] <- NA # assigning NA to the start of the time domain
          pow[j,(((length(pow[1,])-idx2rem)+1):length(pow[1,]))] <- NA # assigning NA to the end of the time domain
        }
        power.mat[[i]] <- pow # storing the edge effect removed power matrices in a new list
      }
      
      power.norm <- list() # create an empty list for storing normalized power matrices
      for (i in 1:length(cwt.no.na)) { # normalization routine for fly-i
        avg.power <- mean(as.vector(power.mat[[i]]), na.rm = T) # calculate average power for fly-i
        power.norm[[i]] <- power.mat[[i]]/avg.power # divide power matrix for fly-i by average power and store in a new variable
      }
      
      power.array <- do.call(cbind, power.norm) # convert the list of normalized power to an array
      power.array <- array(power.array, dim = c(dim(power.norm[[1]]), length(power.norm))) # convert the list of normalized power to an array
      power.mean <- apply(power.array, c(1, 2), mean, na.rm = TRUE) # average across flies
      
      p <- plotly::plot_ly(
        # plot scalograms
      )%>%
        add_trace(
          x = 1:dim(power.mean)[2],
          y = cwt.no.na[[1]]$Period,
          z = power.mean,
          type = "heatmap",
          colors = colorRamp(c("blue", "cyan", "green", "yellow", "red")),
          zauto = F,
          zmin = 0,
          zmax = z.max,
          showscale = T
        )%>%
        layout(
          showlegend = F,
          yaxis = list(
            type = "log",
            showgrid = F,
            showline = T,
            titlefont = f1,
            tickfont = f2,
            title = "Period (h)",
            linecolor = "black",
            autotick = FALSE,
            ticks = "outside",
            tick0 = 0,
            dtick = 1 * s.per.hr,
            ticklen = 7,
            tickcolor = "black",
            range = c(log10((low.per) * s.per.hr), log10((high.per) * s.per.hr)),
            ticktext = as.list(
              c(1,2,4,8,12,17,24,35)
            ),
            tickvals = as.list(
              c(1,2,4,8,12,17,24,35) * s.per.hr
            ),
            tickmode = "array"
          ),
          xaxis = list(
            showgrid = F,
            showline = T,
            titlefont = f1,
            tickfont = f2,
            title = "Days",
            linecolor = "black",
            autotick = FALSE,
            ticks = "outside",
            tick0 = 0,
            dtick = (s.per.hr * 24)/2,
            ticklen = 7,
            tickcolor = "black",
            range = c(0, dim(power.mean)[2]),
            ticktext = as.list(
              seq(0.5, n.days, by = 0.5)
            ),
            tickvals = as.list(
              seq((s.per.hr * 24)/2, s.per.hr * 24 * n.days, (s.per.hr * 24)/2)
            ),
            tickmode = "array"
          )
        )
      
      return(p)
    }
  }
}
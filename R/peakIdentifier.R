#' Phase identifier for activity data
#'
#' @description
#' This function generates a list of outputs. The first element is a plot of raw activity along with the smoothed data and objectively estimated phases of peak of activity data. Smoothing of activity is done using a Savitzky-Golay filter. The input of this function must be the output of the trimData function. This function requires the packages "plotly", "pracma" and "signal". As of now, this function works only for 24-hour T-cycles.
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData(). This assumes that the output from trimData has data binned in 1-min intervals.
#' @param filt.order The filter order. This defaults to 3.
#' @param filt.length The length of filter in indices. This defaults to 51.
#' @param min.peak.dist The minimum distance between peaks to be picked (in indices). This defaults to 100.
#' @param peak.ht.scal A scaling factor to identify peaks. This defaults to 0.5. A value of 0.5 will only find peaks that are at least half as tall as that of the tallest peak.
#' @param ZT A vector defining the regions around which the algorithm must looks for peaks. This defaults to c(0,12). The first value in this vector is always assumed to be the morning peak.
#' @param rm.channels All the channels that users want to remove from their averaging. This must be a vector, i.e., channels must be separated by commas. For instance, if users choose to remove channels 1 to 5, 25 and 32, then the input should be either c(1,2,3,4,5,25,32) or c(1:5,25,32). This defaults to an empty vector, meaning no individuals are removed from analysis.
#'
#' @importFrom plotly plot_ly add_trace subplot %>% layout
#' @importFrom pracma findpeaks isempty
#' @importFrom signal sgolayfilt
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd
#' 
#' @return A \code{list} with two items:
#' \describe{
#' \item{Plots}{A \code{plotly} \code{htmlwidget} with all the averaged activity overlayed with the smoothed data and markers to point out identified peaks in a 4-by-8 array.}
#' \item{Data}{A \code{matrix} \code{array} with 32 rows (one for each fly) and 9 columns (Channel/Fly identity, Phases of morning peak onset, maxima and offset, Morning peak height, Phases of evening peak onset, maxima and offset, and Evening peak height (all phases are measured in ZT)).}
#' }
#' 
#'
#' @export peakIdentifier
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 3, bin = 1, t.cycle = 24)
#' pks <- peakIdentifier(data = td)

peakIdentifier <- function(data, filt.order = 3, filt.length = 51, min.peak.dist = 100, peak.ht.scal = 0.5, ZT = c(0,12), rm.channels = c()) {
  
  requireNamespace("plotly")
  requireNamespace("pracma")
  requireNamespace("signal")
  
  if (requireNamespace("plotly", quietly = T)) {
    # library(plotly)
    # library(pracma)
    # library(signal)
    
    bd <- binData(data = data, input.bin = 1, output.bin = 5, t.cycle = 24)
    
    pro <- profilesAct(data = bd, bin = 5, t.cycle = 24, average.type = "Days", rm.channels = c())
    
    pre.dat <- pro$Profiles
    dat <- rbind(subset(pre.dat, pre.dat$ZT > 18), subset(pre.dat, pre.dat$ZT < 18.01))
    
    p <- list()
    
    f1 <- list(
      family = "Arial, sans-serif",
      size = 20,
      color = "black"
    )
    f2 <- list(
      family = "Arial, sans-serif",
      size = 14,
      color = "black"
    )
    ax <- list(
      showgrid = F,
      showline = T,
      titlefont = f1,
      tickfont = f2,
      title = "ZT (h)",
      linecolor = "black",
      # linewidth = 4,
      # mirror = TRUE,
      autotick = FALSE,
      ticks = "inside",
      tick0 = 0,
      dtick = 12*60/5,
      ticklen = 7,
      tickcolor = "black",
      ticktext = as.list(
        c("18", "00", "06", "12", "18")
      ),
      tickvals = as.list(
        c(0, seq(6*60/5, length(dat[,1]), by = 6*60/5))
      ),
      tickmode = "array",
      # tickwidth = 4,
      range = c(0, length(dat[,1])+1)
    )
    ay <- list(
      showgrid = F,
      showline = T,
      titlefont = f1,
      tickfont = f2,
      title = "Activity",
      linecolor = "black",
      # linewidth = 4,
      # mirror = TRUE,
      autotick = TRUE,
      ticks = "inside",
      tick0 = 0,
      # dtick = max(table.period.power[,"Power"], na.rm = T)/6,
      ticklen = 7,
      tickcolor = "black"
      # range = c(0, max(table.period.power[,"Power"]))
      # tickwidth = 4,
    )
    
    for (i in 1:32) {
      ind = i
      if (requireNamespace("signal", quietly = T)) {
        xx = signal::sgolayfilt(x = na.omit(dat[,1+ind]), p = filt.order, n = filt.length)
      }
      pks <- pracma::findpeaks(xx, minpeakdistance = min.peak.dist, minpeakheight = max(xx)*peak.ht.scal)
      
      
      p[[i]] <- plot_ly(
      )%>%
        add_trace(
          x = 1:length(dat[,1]),
          y = dat[,1+ind],
          type = "scatter",
          mode = "lines",
          line = list(
            color = "black",
            dash = "dash",
            width = 1
          )
        )%>%
        add_trace(
          x = 1:length(xx),
          y = xx,
          type = "scatter",
          mode = "lines",
          line = list(
            color = "red",
            dash = "solid",
            width = 2
          )
        )%>%
        add_trace(
          x = pks[,2],
          y = pks[,1]+3,
          type = "scatter",
          mode = "markers",
          marker = list(
            color = "blue",
            symbol = "triangle-down",
            size = 15
          )
        )%>%
        add_trace(
          x = pks[,3],
          y = 2,
          type = "scatter",
          mode = "markers",
          marker = list(
            color = "green",
            symbol = "star-triangle-down",
            size = 10
          )
        )%>%
        add_trace(
          x = pks[,4],
          y = 2,
          type = "scatter",
          mode = "markers",
          marker = list(
            color = "cyan",
            symbol = "star-triangle-down",
            size = 10
          )
        )%>%
        layout(
          xaxis = ax,
          yaxis = ay
        )
    }
    
    sp <- subplot(p, nrows = 4, shareX = T, shareY = T, margin = 0.01)%>%
      layout(
        showlegend = F
      )
    
    
    phase <- matrix(NA, nrow = 32, ncol = 9)
    colnames(phase) <- c("Channel", "M-Peak.onset", "M-Peak", "M-Peak.offset",
                         "M-Peak.height",
                         "E-Peak.onset", "E-Peak", "E-Peak.offset",
                         "E-Peak.height")
    
    phase[1:32,"Channel"] <- 1:32
    
    idx.for.loop <- setdiff(1:32, rm.channels)
    for (i in idx.for.loop) {
      ind = i
      if (requireNamespace("signal", quietly = T)) {
        xx = signal::sgolayfilt(x = na.omit(dat[,1+ind]), p = filt.order, n = filt.length)
      }
      pks <- pracma::findpeaks(xx, minpeakdistance = min.peak.dist, minpeakheight = max(xx)*peak.ht.scal)
      
      morn.peak.ind <- which(abs(dat[pks[,2],1] - ZT[1]) == min(abs(dat[pks[,2],1] - ZT[1])))
      eve.peak.ind <- which(abs(dat[pks[,2],1] - ZT[2]) == min(abs(dat[pks[,2],1] - ZT[2])))
      
      if (is.na(morn.peak.ind) || ((length(pks[,1]) == 1) && (length(morn.peak.ind) == 1) && length(eve.peak.ind == 1) && min(abs(dat[pks[,2],1] - ZT[2])) < min(abs(dat[pks[,2],1] - ZT[1])))) {
        phase[i,"M-Peak.onset"] = NA
        phase[i,"M-Peak"] = NA
        phase[i,"M-Peak.offset"] = NA
        phase[i,"M-Peak.height"] = NA
      } else {
        phase[i,"M-Peak.onset"] = dat[pks[morn.peak.ind,3],1]
        phase[i,"M-Peak"] = dat[pks[morn.peak.ind,2],1]
        phase[i,"M-Peak.offset"] = dat[pks[morn.peak.ind,4],1]
        phase[i,"M-Peak.height"] = dat[pks[morn.peak.ind,1],1]
      }
      
      if (is.na(eve.peak.ind) || ((length(pks[,1]) == 1) && (length(morn.peak.ind) == 1) && length(eve.peak.ind == 1) && min(abs(dat[pks[,2],1] - ZT[1])) < min(abs(dat[pks[,2],1] - ZT[2])))) {
        phase[i,"E-Peak.onset"] = NA
        phase[i,"E-Peak"] = NA
        phase[i,"E-Peak.offset"] = NA
        phase[i,"E-Peak.height"] = NA
      } else {
        phase[i,"E-Peak.onset"] = dat[pks[eve.peak.ind,3],1]
        phase[i,"E-Peak"] = dat[pks[eve.peak.ind,2],1]
        phase[i,"E-Peak.offset"] = dat[pks[eve.peak.ind,4],1]
        phase[i,"E-Peak.height"] = dat[pks[eve.peak.ind,1],1]
      }
    }
    
    output <- list(
      "Plots" = sp,
      "Data" = phase
    )
    
    return(output)
  }
  
}

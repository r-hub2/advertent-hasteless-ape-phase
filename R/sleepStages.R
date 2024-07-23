#' Plot sleep stages in individual flies
#'
#' @description
#' Allows users to generate individual plots of sleep stages in flies. Sleep stages are defined as follows: 5 to 30-min as short sleep (light blue), 30 to 60-min as intermediate sleep (medium blue) and 60 to 720-min as deep sleep (dark blue). Activity is plotted in red. The input for this function must be the output of the trimData() function. The output of this function is a plot. Note: At this moment, this works accurately only for 24-h days.
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param n.days The number of cycles for which sleep stages must be visualized.
#' @param channel The channel number of the fly to be visualized.
#' @param photoperiod This value determines the duration of photo-phase and scoto-phase of the 24-h day. Defaults to 12.
#' 
#' @importFrom graphics abline lines mtext par polygon
#' @importFrom utils head tail
#' 
#' @return A \code{plot} with the day-time and night-time sleep stages of a user defined fly for the number of cycles provided by the user.
#'
#' @export sleepStages
#'
#' @examples
#' \dontrun{
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 2, bin = 1, t.cycle = 24)
#' sleepStages(data = td, n.days = 1)
#' }

sleepStages <- function (data, n.days, channel = 1, photoperiod = 12) {
  
  if (requireNamespace("plotly", quietly = T)) {
    bd <- binData(data = data, input.bin = 1, output.bin = 1, t.cycle = 24)[,-1]
    bd[bd != 0] <- 1
    sd <- list()
    
    sleep.def <- c(5,30,60,720)
    
    for (i in 1:(length(sleep.def) - 1)) {
      sd[[i]] <- sleepData(data = data, sleep.def = c(sleep.def[i], sleep.def[i+1]), bin = 1, t.cycle = 24)[,-1]
    }
    
    day.startindex <- seq(1, (1440 * n.days), by = 1440)
    day.endindex <- seq((photoperiod * 60), (1440 * n.days), by = 1440)
    night.startindex <- seq((photoperiod * 60)+1, (1440 * n.days), by = 1440)
    night.endindex <- seq(1440, (1440 * n.days), by = 1440)
    
    par(mfrow = c(n.days,2), oma = c(2,4,5,2),
        mai = c(0.2,0.2,0.2,0.2), bg = rgb(1,1,1,1))
    col <- c(
      rgb(130/255, 207/255, 251/255, 1),
      rgb(59/255, 130/255, 246/255, 1),
      rgb(54/255, 53/255, 157/255, 1)
    )
    # DAY
    plot(1:(photoperiod * 60), bd[day.startindex[1]:day.endindex[1],channel]+3,
         type = "h", ylim = c(0,4), xaxt = "n", yaxt = "n",
         col = rgb(251/255, 118/255, 89/255, 0), lwd = 2, main = "", col.main = "black", ylab = paste("Cycle-", 1, sep = ""), col.lab = "black", yaxs = "i")
    box(col = "black")
    mtext(text = paste("Cycle-", 1, sep = ""), side = 2, outer = F, col = "black", line = 1.5, cex = 1.2)
    mtext(text = "Day", side = 3, outer = F, col = "black", line = 1.5, cex = 1.1)
    abline(h = c(3,4), col = "gray")
    polygon(x = c(min(1:(photoperiod * 60)), 1:(photoperiod * 60), max(1:(photoperiod * 60))),
            y = c(min(bd[day.startindex[1]:day.endindex[1],channel]+3),
                  bd[day.startindex[1]:day.endindex[1],channel]+3,
                  min(bd[day.startindex[1]:day.endindex[1],channel]+3)),
            col = rgb(251/255, 118/255, 89/255, 1), border = NA)
    to.add <- c(2,1,0)
    mtext(text = paste("Channel-", channel, sep = ""), side = 3, line = 2, outer = T, col = "black", cex = 1.3)
    for (i in 1:length(sd)) {
      lines(1:(photoperiod * 60), sd[[i]][day.startindex[1]:day.endindex[1],channel]+to.add[i], type = "h", col = rgb(0,0,0,0))
      abline(h = to.add[i], col = "gray")
      polygon(x = c(min(1:(photoperiod * 60)), 1:(photoperiod * 60), max(1:(photoperiod * 60))),
              y = c(min(sd[[i]][day.startindex[1]:day.endindex[1],channel]+to.add[i]),
                    sd[[i]][day.startindex[1]:day.endindex[1],channel]+to.add[i],
                    min(sd[[i]][day.startindex[1]:day.endindex[1],channel]+to.add[i])),
              col = col[i], border = NA)
    }
    # NIGHT
    plot(1:(1440 - (photoperiod * 60)), bd[night.startindex[1]:night.endindex[1],channel]+3,
         type = "h", ylim = c(0,4), xaxt = "n", yaxt = "n",
         col = rgb(251/255, 118/255, 89/255, 0), lwd = 2, main = "", col.main = "black", yaxs = "i")
    mtext(text = "Night", side = 3, outer = F, col = "black", line = 1.5, cex = 1.1)
    box(col = "black")
    abline(h = c(3,4), col = "gray")
    polygon(x = c(min(1:(1440 - (photoperiod * 60))), 1:(1440 - (photoperiod * 60)), max(1:(1440 - (photoperiod * 60)))),
            y = c(min(bd[night.startindex[1]:night.endindex[1],channel]+3),
                  bd[night.startindex[1]:night.endindex[1],channel]+3,
                  min(bd[night.startindex[1]:night.endindex[1],channel]+3)),
            col = rgb(251/255, 118/255, 89/255, 1), border = NA)
    to.add <- c(2,1,0)
    for (i in 1:length(sd)) {
      lines(1:(1440 - (photoperiod * 60)), sd[[i]][night.startindex[1]:night.endindex[1],channel]+to.add[i], type = "h", col = rgb(0,0,0,0))
      abline(h = to.add[i], col = "gray")
      polygon(x = c(min(1:(1440 - (photoperiod * 60))), 1:(1440 - (photoperiod * 60)), max(1:(1440 - (photoperiod * 60)))),
              y = c(min(sd[[i]][night.startindex[1]:night.endindex[1],channel]+to.add[i]),
                    sd[[i]][night.startindex[1]:night.endindex[1],channel]+to.add[i],
                    min(sd[[i]][night.startindex[1]:night.endindex[1],channel]+to.add[i])),
              col = col[i], border = NA)
    }
    # DAY
    for (jj in 2:n.days) {
      plot(1:(photoperiod * 60), bd[day.startindex[jj]:day.endindex[jj],channel]+3,
           type = "h", ylim = c(0,4), xaxt = "n", yaxt = "n",
           col = rgb(251/255, 118/255, 89/255, 0), lwd = 2, main = "", col.main = "black", ylab = paste("Cycle-", jj, sep = ""), col.lab = "black", yaxs = "i", cex.main = 1.1)
      box(col = "black")
      mtext(text = paste("Cycle-", jj, sep = ""), side = 2, outer = F, col = "black", line = 1.5, cex = 1.2)
      abline(h = c(3,4), col = "gray")
      polygon(x = c(min(1:(photoperiod * 60)), 1:(photoperiod * 60), max(1:(photoperiod * 60))),
              y = c(min(bd[day.startindex[jj]:day.endindex[jj],channel]+3),
                    bd[day.startindex[jj]:day.endindex[jj],channel]+3,
                    min(bd[day.startindex[jj]:day.endindex[jj],channel]+3)),
              col = rgb(251/255, 118/255, 89/255, 1), border = NA)
      to.add <- c(2,1,0)
      mtext(text = paste("Channel-", channel, sep = ""), side = 3, line = 2, outer = T, col = "black", cex = 1.3)
      for (i in 1:length(sd)) {
        lines(1:(photoperiod * 60), sd[[i]][day.startindex[jj]:day.endindex[jj],channel]+to.add[i], type = "h", col = rgb(0,0,0,0))
        abline(h = to.add[i], col = "gray")
        polygon(x = c(min(1:(photoperiod * 60)), 1:(photoperiod * 60), max(1:(photoperiod * 60))),
                y = c(min(sd[[i]][day.startindex[jj]:day.endindex[jj],channel]+to.add[i]),
                      sd[[i]][day.startindex[jj]:day.endindex[jj],channel]+to.add[i],
                      min(sd[[i]][day.startindex[jj]:day.endindex[jj],channel]+to.add[i])),
                col = col[i], border = NA)
      }
      # NIGHT
      plot(1:(1440 - (photoperiod * 60)), bd[night.startindex[jj]:night.endindex[jj],channel]+3,
           type = "h", ylim = c(0,4), xaxt = "n", yaxt = "n",
           col = rgb(251/255, 118/255, 89/255, 0), lwd = 2, main = "", col.main = "black", yaxs = "i", cex.main = 1.1)
      box(col = "black")
      abline(h = c(3,4), col = "gray")
      polygon(x = c(min(1:(1440 - (photoperiod * 60))), 1:(1440 - (photoperiod * 60)), max(1:(1440 - (photoperiod * 60)))),
              y = c(min(bd[night.startindex[jj]:night.endindex[jj],channel]+3),
                    bd[night.startindex[jj]:night.endindex[jj],channel]+3,
                    min(bd[night.startindex[jj]:night.endindex[jj],channel]+3)),
              col = rgb(251/255, 118/255, 89/255, 1), border = NA)
      to.add <- c(2,1,0)
      for (i in 1:length(sd)) {
        lines(1:(1440 - (photoperiod * 60)), sd[[i]][night.startindex[jj]:night.endindex[jj],channel]+to.add[i], type = "h", col = rgb(0,0,0,0))
        abline(h = to.add[i], col = "gray")
        polygon(x = c(min(1:(1440 - (photoperiod * 60))), 1:(1440 - (photoperiod * 60)), max(1:(1440 - (photoperiod * 60)))),
                y = c(min(sd[[i]][night.startindex[jj]:night.endindex[jj],channel]+to.add[i]),
                      sd[[i]][night.startindex[jj]:night.endindex[jj],channel]+to.add[i],
                      min(sd[[i]][night.startindex[jj]:night.endindex[jj],channel]+to.add[i])),
                col = col[i], border = NA)
      }
    }
  }
}
#' Computes and tabulates day-time and night-time sleep statistics
#'
#' @description
#' This function allows users to estimate day-time and night-time average sleep bout duration, number and latency. Sleep bout latency is defined as the time taken (in minutes) for the occurrence of the first sleep bout since respective transitions. The input for this function must be the output from the trimData() function. The output of this function is a matrix which contains fly-wise (each row) data.
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep or specific bout durations; sleep.def allows users to change the definition of sleep. The default input is a single value vector of value 5. If users wish to analyse sleep only between 5 to 20 mins, the input must be c(5,20).
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param photoperiod Duration (in hours) of what can be considered day-phase. This defaults to 12.
#'
#' @importFrom grDevices rgb
#' @importFrom stats aggregate fitted lm na.omit sd median
#' 
#' @return A \code{matrix} \code{array} matrix with 32 rows (one for each fly) and 9 columns:
#' \describe{
#' \item{Channel}{Fly identity.}
#' \item{Day.BoutNumber}{Number of sleep bouts in the user defined day time.}
#' \item{Day.BoutDuration.Mean}{Mean sleep duration in the user defined day time.}
#' \item{Day.BoutDuration.Median}{Median sleep duration in the user defined day time.}
#' \item{Day.Latency}{Time taken for the first sleep bout to occur in the user defined day time.}
#' \item{Day.Total}{Total minutes of daytime sleep.}
#' \item{Night.BoutNumber}{Number of sleep bouts in the user defined night time.}
#' \item{Night.BoutDuration.Mean}{Mean sleep duration in the user defined night time.}
#' \item{Night.BoutDuration.Median}{Median sleep duration in the user defined night time.}
#' \item{Night.Latency}{Time taken for the first sleep bout to occur in the user defined night time.}
#' \item{Night.Total}{Total minutes of nighttime sleep.}
#' }
#'
#' @export sleepStat
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 10, bin = 1, t.cycle = 24)
#' slp.stat <- sleepStat(data = td)

sleepStat <- function(data, sleep.def = c(5), t.cycle = 24, photoperiod = 12) {
  
  if (length(sleep.def) == 1) {
    
    raw <- data[,-c(1:10)]
    
    for (i in 1:length(raw[1,])) {
      x <- raw[,i]
      y <- rle(x)
      d_y <- as.data.frame(unclass(y))
      d_y$end <- cumsum(d_y$lengths)
      d_y$start <- d_y$end - d_y$lengths + 1
      
      dd_y <- subset(d_y, d_y$values == 0 & d_y$lengths >= sleep.def)
      
      if(length(dd_y[,1]) == 0) {
        x = 0
      } else {
        for (j in 1:length(dd_y[,1])) {
          x[dd_y[j,"start"]:dd_y[j,"end"]] = -1
        }
      }
      
      x[x > -1] = 0
      x[x == -1] = 1
      raw[,i] <- x
    }
    
  } else if (length(sleep.def) == 2) {
    
    raw <- data[,-c(1:10)]
    
    for (i in 1:length(raw[1,])) {
      x <- raw[,i]
      y <- rle(x)
      d_y <- as.data.frame(unclass(y))
      d_y$end <- cumsum(d_y$lengths)
      d_y$start <- d_y$end - d_y$lengths + 1
      
      dd_y <- subset(d_y, d_y$values == 0 & d_y$lengths >= sleep.def[1] & lengths < sleep.def[2])
      
      if(length(dd_y[,1]) == 0) {
        x = 0
      } else {
        for (j in 1:length(dd_y[,1])) {
          x[dd_y[j,"start"]:dd_y[j,"end"]] = -1
        }
      }
      
      x[x > -1] = 0
      x[x == -1] = 1
      raw[,i] <- x
    }
  }
  
  s_per_day <- 60*t.cycle
  n.days <- length(raw[,1])/s_per_day
  
  cyc.wise.split <- list()
  cyc.wise.index <- seq(1, length(raw[,1]), by = 1440)
  for (i in 1:n.days) {
    cyc.wise.split[[i]] <- raw[cyc.wise.index[i]:(cyc.wise.index[i]+(1440-1)),]
  }
  
  output <- matrix(NA, nrow = 32, ncol = 11)
  colnames(output) <- c("Channel", "Day.BoutNumber", "Day.BoutDuration.Mean", "Day.BoutDuration.Median",
                        "Day.Latency", "Day.Total",
                        "Night.BoutNumber", "Night.BoutDuration.Mean", "Night.BoutDuration.Median",
                        "Night.Latency", "Night.Total")
  for.row.nam <- seq(1, 32, by = 1)
  
  output[,"Channel"] <- for.row.nam
  
  day.cyc.boutnum <- matrix(NA, nrow = 32, ncol = n.days)
  day.cyc.boutdur.mean <- matrix(NA, nrow = 32, ncol = n.days)
  day.cyc.boutdur.median <- matrix(NA, nrow = 32, ncol = n.days)
  day.cyc.boutlatency <- matrix(NA, nrow = 32, ncol = n.days)
  day.tot <- matrix(NA, nrow = 32, ncol = n.days)
  
  night.cyc.boutnum <- matrix(NA, nrow = 32, ncol = n.days)
  night.cyc.boutdur.mean <- matrix(NA, nrow = 32, ncol = n.days)
  night.cyc.boutdur.median <- matrix(NA, nrow = 32, ncol = n.days)
  night.cyc.boutlatency <- matrix(NA, nrow = 32, ncol = n.days)
  night.tot <- matrix(NA, nrow = 32, ncol = n.days)
  
  for (i in 1:length(cyc.wise.split)) {
    day <- cyc.wise.split[[i]][(1:(photoperiod*60)),]
    night <- cyc.wise.split[[i]][(((photoperiod*60) + 1):1440),]
    
    day.tot[,i] = colSums(day)
    night.tot[,i] = colSums(night)
    
    for (j in 1:length(day[1,])) {
      d <- day[,j]
      d.rle <- rle(d)
      df.d <- as.data.frame((unclass(d.rle)))
      df.d$end <- cumsum(df.d$lengths)
      df.d$start <- df.d$end - df.d$lengths + 1
      ddf.d <- subset(df.d, df.d$values == 1)
      
      day.cyc.boutnum[j,i] <- length(ddf.d[,1])
      day.cyc.boutdur.mean[j,i] <- mean(ddf.d$lengths)
      day.cyc.boutdur.median[j,i] <- median(ddf.d$lengths)
      day.cyc.boutlatency[j,i] <- ddf.d[1,"start"]
      
      n <- night[,j]
      n.rle <- rle(n)
      df.n <- as.data.frame((unclass(n.rle)))
      df.n$end <- cumsum(df.n$lengths)
      df.n$start <- df.n$end - df.n$lengths + 1
      ddf.n <- subset(df.n, df.n$values == 1)
      
      night.cyc.boutnum[j,i] <- length(ddf.n[,1])
      night.cyc.boutdur.mean[j,i] <- mean(ddf.n$lengths)
      night.cyc.boutdur.median[j,i] <- median(ddf.n$lengths)
      night.cyc.boutlatency[j,i] <- ddf.n[1,"start"]
    }
  }
  
  for (i in 1:length(raw[1,])) {
    output[i,"Day.BoutNumber"] = mean(day.cyc.boutnum[i,], na.rm = T)
    output[i,"Night.BoutNumber"] = mean(night.cyc.boutnum[i,], na.rm = T)
    
    output[i,"Day.BoutDuration.Mean"] = mean(day.cyc.boutdur.mean[i,], na.rm = T)
    output[i,"Night.BoutDuration.Mean"] = mean(night.cyc.boutdur.mean[i,], na.rm = T)
    
    output[i,"Day.BoutDuration.Median"] = mean(day.cyc.boutdur.median[i,], na.rm = T)
    output[i,"Night.BoutDuration.Median"] = mean(night.cyc.boutdur.median[i,], na.rm = T)
    
    output[i,"Day.Latency"] = mean(day.cyc.boutlatency[i,], na.rm = T)
    output[i,"Night.Latency"] = mean(night.cyc.boutlatency[i,], na.rm = T)
    
    output[i,"Day.Total"] = mean(day.tot[i,], na.rm = T)
    output[i,"Night.Total"] = mean(night.tot[i,], na.rm = T)
  }
  
  return(output)
}
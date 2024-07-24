#' Generate data from Continuous Wavelet Transforms
#'
#' @description
#' This function generates CWT results for each fly. Input for this function must be an output from the trimData() function. The output of this function is a large list. This function requires the packages "WaveletComp" and " boot".
#'
#' @param input Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param data A choice between "Activity" and "Sleep" data to analyze using the wavelet transform.
#' @param out.bin Define the desired output bin size (in minutes). This defaults to 15.
#' @param t.cycle Define the period of the environmental cycle or a single day in hours. This defaults to 24.
#' @param sleep.def Definition of sleep. Traditionally, a single bout of sleep is defined as any duration of inactivity that is equal to or greater than 5-minutes. However, sometimes it may be of interest to examine longer bouts of sleep or specific bout durations; sleep.def allows users to change the definition of sleep. The default input is a single value vector of value 5. If users wish to analyse sleep only between 5 to 20 mins, the input must be c(5,20). Only relevant if data = "Sleep".
#' @param low.per Choose the lowest period (in hours) for analysis. This defaults to 1.
#' @param high.per Choose the highest period (in hours) for analysis. This defaults to 35.
#' @param rm.channels All the channels that users want to remove from their averaging. This must be a vector, i.e., channels must be separated by commas. For instance, if users choose to remove channels 1 to 5, 25 and 32, then the input should be either c(1,2,3,4,5,25,32) or c(1:5,25,32). This defaults to an empty vector, meaning no individuals are removed from analysis.
#' @param make.pval A logical (TRUE or FALSE). TRUE will generate p-values for the significance of detected periodicities. FALSE will not. Note that generating p-values significantly increases computation time. Defaults to FALSE.
#' @param n.sim If TRUE, the number of simulations based on which p-values will be computed. Defaults to 1000.
#' @param method The method used for generating surrogate time-series against which p-values are generated. Also see help(analyze.wavelet). Available choices are "white.noise": white noise, "shuffle": shuffling the given time-series, "Fourier.rand": time-series with a similar spectrum, "AR": AR(p), "ARIMA": ARIMA(p,0,q). Defaults to "shuffle".
#' @param boot A logical (TRUE or FALSE). TRUE will bootstrap 95% Confidence Intervals (CI) for the normalized amplitude for each period value. Defaults to TRUE.
#' @param boot.rep Number of replicates for bootstrapping 95% CIs. Defaults to 1000. Only relevant if boot = TRUE.
#'
#' @importFrom WaveletComp analyze.wavelet
#' @importFrom boot boot.ci
#' 
#' @return A \code{list} with two items:
#' \describe{
#' \item{CWT}{A \code{list} with all the details of the CWT computations for each fly, after removing dead or empty channels.}
#' \item{PeriodPower}{A \code{matrix} \code{array} with and 36 columns. The first columns stores period values. Columns 2 through 33 are for each fly. Column 34 is the averaged normalized amplitude (over flies) for each period value. Columns 35 and 36 store the lower and upper bounds of the 95% CI.}
#' }
#'
#' @export phaseCWT
#'
#' @examples
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 3, bin = 1, t.cycle = 24)
#' cwt.dat <- phaseCWT(input = td, low.per = 1, high.per = 3, boot = F)

phaseCWT <- function (input, data = "Activity", out.bin = 15, t.cycle = 24, sleep.def = 5, low.per = 1, high.per = 35, rm.channels = c(), make.pval = FALSE, n.sim = 1000, method = "shuffle", boot = TRUE, boot.rep = 1000) {
  
  requireNamespace("WaveletComp")
  requireNamespace("boot")
  
  s.per.hr = 60/out.bin
  
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
  
  power.mat <- list() # create an empty list for storing edge effects removed power matrices
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
  
  avg.power <- as.data.frame(matrix(NA, nrow = length(cwt.no.na[[1]]$Period), ncol = 33)) # create an empty data frame to store time averaged power for each fly
  colnames(avg.power) <- c("period", paste(rep("Ch", 32), 1:32)) # name column headers
  avg.power$period <- cwt.no.na[[1]]$Period/s.per.hr # fill period values
  
  for (i in 1:length(cwt.no.na)) { # calculate time averaged power for fly-i
    ind.pow <- power.norm[[i]] # store power matrix for fy-i in a separate variable
    avg.power[,i+1] <- rowMeans(ind.pow, na.rm = T) # carry out average across time for fly-i
  }
  
  avg.power[,rm.channels+1] <- NA
  
  avg.power$mean <- rowMeans(avg.power[,2:33], na.rm = T) # compute mean power across flies for each period
  
  if (isTRUE(boot)) {
    # bootstrap and create 95%CI
    for (i in 1:length(avg.power[,1])) { # use loop to bootstrap 95%CI for period-i
      tryCatch({
        mean.fx <- function(x, index) { # define function to calculate mean
          d <- x[index]
          return(mean(d, na.rm = T))
        }
        boot.x <- boot(data = as.numeric(avg.power[i,2:33]), statistic = mean.fx, R = boot.rep) # bootstrap for period-i
        bci <- boot::boot.ci(boot.out = boot.x, conf = 0.95, type = "bca") # estimate CI
        avg.power[i,"lower.ci95"] <- bci$bca[4] # assign lower bound
        avg.power[i,"upper.ci95"] <- bci$bca[5] # assign upper bound
      }, error = function(e){})
    }
  }
  
  output <- list(
    "CWT" = cwt.no.na,
    "PeriodPower" = avg.power
  )
  
  return(output)
}
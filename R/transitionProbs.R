#' Plot and generate data for transition frequencies between sleep stages
#'
#' @description
#' This function generates a sankey plot enabling the visualization of percentage of times transtions from various states of activity and sleep occur. This function also generates a data table to download individual values of these transition percentages. Sleep stages are defined as follows: 5 to 30-min as short sleep (light blue), 30 to 60-min as intermediate sleep (medium blue) and 60 to 720-min as deep sleep (dark blue). Activity is defined as the standard number of beam crosses. The input for this function must be the output of the trimData() function. The output of this function is a list. Note: At this moment, this works accurately only for 24-h days.
#'
#' @param data Input data file. The input for this function must be the output of the function trimData(). See ??trimData().
#' @param n.days The number of cycles for which sleep stages must be visualized.
#' @param photoperiod This value determines the duration of photo-phase and scoto-phase of the 24-h day. Defaults to 12.
#' 
#' @importFrom wesanderson wes_palette
#' @importFrom plotly plot_ly add_trace layout %>% subplot
#' 
#' @return A \code{list} with two items:
#' \describe{
#' \item{Plots}{A \code{plotly} \code{htmlwidget} with the Sankey diagram (representing the average transition percentages).}
#' \item{Data}{A \code{matrix} \code{array} with individual wise transition percentages.}
#' }
#'
#' @export transitionProbs
#'
#' @examples
#' \dontrun{
#' td <- trimData(data = df, start.date = "19 Dec 20", start.time = "21:00",
#' n.days = 2, bin = 1, t.cycle = 24)
#' tp <- transitionProbs(data = td, n.days = 1)
#' }


transitionProbs <- function (data, n.days, photoperiod = 12) {
  requireNamespace("wesanderson")
  requireNamespace("plotly")
  
  if (requireNamespace("plotly", quietly = T)) {
    day.start.ind <- seq(1, 1440*n.days, by = 1440)
    day.end.ind <- seq((photoperiod * 60), 1440*n.days, by = 1440)
    night.start.ind <- seq((photoperiod * 60)+1, 1440*n.days, by = 1440)
    night.end.ind <- seq(1440, 1440*n.days, by = 1440)
    
    pre.day.td <- matrix(NA, nrow = length(data[,1]), ncol = 32)
    pre.night.td <- matrix(NA, nrow = length(data[,1]), ncol = 32)
    
    for (i in 1:n.days) {
      for (j in 1:32) {
        pre.day.td[day.start.ind[i]:day.end.ind[i],j] <- data[day.start.ind[i]:day.end.ind[i],10+j]
        pre.night.td[night.start.ind[i]:night.end.ind[i],j] <- data[night.start.ind[i]:night.end.ind[i],10+j]
      }
    }
    
    day.td <- as.data.frame(na.omit(pre.day.td))
    night.td <- as.data.frame(na.omit(pre.night.td))
    
    # CODE FOR DIFFERENT SLEEP CRITERIA: 1-Wake, 2-Short sleep, 3-Intermediate sleep, 4-Long sleep
    ####################
    ####################
    ####################
    ### For Daytime ####
    ####################
    ####################
    ####################
    out.table.day <- data.frame(matrix(NA, nrow = 12, ncol = 34))
    colnames(out.table.day)[1:2] <- c("From", "To")
    out.table.day[,"From"] <- c(
      "From Wake", "From Wake", "From Wake",
      "From Short", "From Short", "From Short",
      "From Inter", "From Inter", "From Inter",
      "From Long", "From Long", "From Long"
    )
    out.table.day[,"To"] <- c(
      "To Short", "To Inter", "To Long",
      "To Wake", "To Inter", "To Long",
      "To Wake", "To Short", "To Long",
      "To Wake", "To Short", "To Inter"
    )
    
    for (kk in 1:32) {
      x.d <- day.td[,kk]
      
      x.d[x.d == 0] <- -1 # Sleep
      x.d[x.d > 0] <- -4 # Activity
      
      y.d <- rle(x.d)
      
      d_y.d <- as.data.frame(unclass(y.d))
      d_y.d$end <- cumsum(d_y.d$lengths)
      d_y.d$start <- d_y.d$end - d_y.d$lengths + 1
      
      d_y.d$values[d_y.d$lengths > 5 & d_y.d$lengths < 30] <- -2
      d_y.d$values[d_y.d$lengths > 30 & d_y.d$lengths < 60] <- -3
      
      
      s.d <- subset(d_y.d, d_y.d$values == -2)
      ss.d <- subset(d_y.d, d_y.d$values == -3)
      for (i in 1:length(s.d[,1])) {
        if(!is.na(s.d[1,1])) {
          x.d[s.d$start[i]:s.d$end[i]] <- 2 
        }
      }
      
      for (i in 1:length(ss.d[,1])) {
        if(!is.na(ss.d[1,1])) {
          x.d[ss.d$start[i]:ss.d$end[i]] <- 3
        }
      }
      
      x.d[x.d == -4] <- 4
      x.d[x.d == -1] <- 1
      
      x1.d <- factor(paste0(head(x.d,-1), tail(x.d,-1)), levels = c('12','13','14','21','23','24','31','32','34','41','42','43'))
      
      # xdf <- data.frame(x = x[1:length(x) - 1], y = x[2:length(x)])
      
      t.xdf.d <- table(x1.d)
      
      out.table.day[1,2+kk] <- (as.numeric(t.xdf.d["12"])/sum(as.numeric(t.xdf.d["12"]), as.numeric(t.xdf.d["13"]), as.numeric(t.xdf.d["14"]))) * 100
      out.table.day[2,2+kk] <- (as.numeric(t.xdf.d["13"])/sum(as.numeric(t.xdf.d["12"]), as.numeric(t.xdf.d["13"]), as.numeric(t.xdf.d["14"]))) * 100
      out.table.day[3,2+kk] <- (as.numeric(t.xdf.d["14"])/sum(as.numeric(t.xdf.d["12"]), as.numeric(t.xdf.d["13"]), as.numeric(t.xdf.d["14"]))) * 100
      
      out.table.day[4,2+kk] <- (as.numeric(t.xdf.d["21"])/sum(as.numeric(t.xdf.d["21"]), as.numeric(t.xdf.d["23"]), as.numeric(t.xdf.d["24"]))) * 100
      out.table.day[5,2+kk] <- (as.numeric(t.xdf.d["23"])/sum(as.numeric(t.xdf.d["21"]), as.numeric(t.xdf.d["23"]), as.numeric(t.xdf.d["24"]))) * 100
      out.table.day[6,2+kk] <- (as.numeric(t.xdf.d["24"])/sum(as.numeric(t.xdf.d["21"]), as.numeric(t.xdf.d["23"]), as.numeric(t.xdf.d["24"]))) * 100
      
      out.table.day[7,2+kk] <- (as.numeric(t.xdf.d["31"])/sum(as.numeric(t.xdf.d["31"]), as.numeric(t.xdf.d["32"]), as.numeric(t.xdf.d["34"]))) * 100
      out.table.day[8,2+kk] <- (as.numeric(t.xdf.d["32"])/sum(as.numeric(t.xdf.d["31"]), as.numeric(t.xdf.d["32"]), as.numeric(t.xdf.d["34"]))) * 100
      out.table.day[9,2+kk] <- (as.numeric(t.xdf.d["34"])/sum(as.numeric(t.xdf.d["31"]), as.numeric(t.xdf.d["32"]), as.numeric(t.xdf.d["34"]))) * 100
      
      out.table.day[10,2+kk] <- (as.numeric(t.xdf.d["41"])/sum(as.numeric(t.xdf.d["41"]), as.numeric(t.xdf.d["42"]), as.numeric(t.xdf.d["43"]))) * 100
      out.table.day[11,2+kk] <- (as.numeric(t.xdf.d["42"])/sum(as.numeric(t.xdf.d["41"]), as.numeric(t.xdf.d["42"]), as.numeric(t.xdf.d["43"]))) * 100
      out.table.day[12,2+kk] <- (as.numeric(t.xdf.d["43"])/sum(as.numeric(t.xdf.d["41"]), as.numeric(t.xdf.d["42"]), as.numeric(t.xdf.d["43"]))) * 100
      
      # out.table[1,2+kk] <- as.numeric(t.xdf["12"])
      # out.table[2,2+kk] <- as.numeric(t.xdf["13"])
      # out.table[3,2+kk] <- as.numeric(t.xdf["14"])
      # 
      # out.table[4,2+kk] <- as.numeric(t.xdf["21"])
      # out.table[5,2+kk] <- as.numeric(t.xdf["23"])
      # out.table[6,2+kk] <- as.numeric(t.xdf["24"])
      # 
      # out.table[7,2+kk] <- as.numeric(t.xdf["31"])
      # out.table[8,2+kk] <- as.numeric(t.xdf["32"])
      # out.table[9,2+kk] <- as.numeric(t.xdf["34"])
      # 
      # out.table[10,2+kk] <- as.numeric(t.xdf["41"])
      # out.table[11,2+kk] <- as.numeric(t.xdf["42"])
      # out.table[12,2+kk] <- as.numeric(t.xdf["43"])
      
      # out.table[1,2+kk] <- (t.xdf["1","2"]/sum(t.xdf["1",-1], na.rm = T))*100
      # out.table[2,2+kk] <- (t.xdf["1","3"]/sum(t.xdf["1",-1], na.rm = T))*100
      # out.table[3,2+kk] <- (t.xdf["1","4"]/sum(t.xdf["1",-1], na.rm = T))*100
      # 
      # out.table[4,2+kk] <- (t.xdf["2","1"]/sum(t.xdf["2",-2], na.rm = T))*100
      # out.table[5,2+kk] <- (t.xdf["2","3"]/sum(t.xdf["2",-2], na.rm = T))*100
      # out.table[6,2+kk] <- (t.xdf["2","4"]/sum(t.xdf["2",-2], na.rm = T))*100
      # 
      # out.table[7,2+kk] <- (t.xdf["3","1"]/sum(t.xdf["3",-3], na.rm = T))*100
      # out.table[8,2+kk] <- (t.xdf["3","2"]/sum(t.xdf["3",-3], na.rm = T))*100
      # out.table[9,2+kk] <- (t.xdf["3","4"]/sum(t.xdf["3",-3], na.rm = T))*100
      # 
      # out.table[10,2+kk] <- (t.xdf["4","1"]/sum(t.xdf["4",-4], na.rm = T))*100
      # out.table[11,2+kk] <- (t.xdf["4","2"]/sum(t.xdf["4",-4], na.rm = T))*100
      # out.table[12,2+kk] <- (t.xdf["4","3"]/sum(t.xdf["4",-4], na.rm = T))*100
      
    }
    
    out.table.day$mean <- rowMeans(out.table.day[,-c(1:2)], na.rm = T)
    
    nodes.day <- data.frame(
      name=c(as.character(out.table.day$From), 
             as.character(out.table.day$To)) %>% unique()
    )
    
    out.table.day$IDfrom <- match(out.table.day$From, nodes.day$name)-1
    out.table.day$IDto <- match(out.table.day$To, nodes.day$name)-1
    
    q.day <- plot_ly(
      type = "sankey",
      orientation = "h",
      domain = list(
        x = c(0,0.5),
        y = c(0,1)
      ),
      arrangement = "snap",
      node = list(
        label = as.factor(nodes.day$name),
        x = c(0.05,0.05,0.05,0.05,0.95,0.95,0.95,0.95),
        y = c(0.1,0.2,0.3,0.6,-1,-0.6,-0.2,0.5),
        color = c(
          wes_palette("FantasticFox1")[3], wes_palette("FantasticFox1")[2], wes_palette("FantasticFox1")[1], wes_palette("FantasticFox1")[5],
          wes_palette("FantasticFox1")[2], wes_palette("FantasticFox1")[1], wes_palette("FantasticFox1")[5], wes_palette("FantasticFox1")[3]
        ),
        pad = 15,
        thickness = 20,
        line = list(
          color = "black",
          width = 0.5
        )
      ),
      link = list(
        source = out.table.day$IDfrom,
        target = out.table.day$IDto,
        value = out.table.day$mean,
        color = list(
          rgb(70/255,172/255,200/255,0.2),
          rgb(70/255,172/255,200/255,0.2),
          rgb(70/255,172/255,200/255,0.2),
          rgb(226/255,210/255,0,0.2),
          rgb(226/255,210/255,0,0.2),
          rgb(226/255,210/255,0,0.2),
          rgb(221/255,141/255,41/255,0.2),
          rgb(221/255,141/255,41/255,0.2),
          rgb(221/255,141/255,41/255,0.2),
          rgb(180/255,15/255,32/255,0.2),
          rgb(180/255,15/255,32/255,0.2),
          rgb(180/255,15/255,32/255,0.2)
        )
      )
    )
    
    out.table.for.output.day <- out.table.day[,1:35]
    colnames(out.table.for.output.day) <- c(
      "Source", "Target",
      "Ch1", "Ch2", "Ch3", "Ch4", "Ch5", "Ch6", "Ch7", "Ch8",
      "Ch9", "Ch10", "Ch11", "Ch12", "Ch13", "Ch14", "Ch15", "Ch16",
      "Ch17", "Ch18", "Ch19", "Ch20", "Ch21", "Ch22", "Ch23", "Ch24",
      "Ch25", "Ch26", "Ch27", "Ch28", "Ch29", "Ch30", "Ch31", "Ch32",
      "Mean"
    )
    
    ####################
    ####################
    ####################
    ### For Nighttime ####
    ####################
    ####################
    ####################
    out.table.night <- data.frame(matrix(NA, nrow = 12, ncol = 34))
    colnames(out.table.night)[1:2] <- c("From", "To")
    out.table.night[,"From"] <- c(
      "From Wake", "From Wake", "From Wake",
      "From Short", "From Short", "From Short",
      "From Inter", "From Inter", "From Inter",
      "From Long", "From Long", "From Long"
    )
    out.table.night[,"To"] <- c(
      "To Short", "To Inter", "To Long",
      "To Wake", "To Inter", "To Long",
      "To Wake", "To Short", "To Long",
      "To Wake", "To Short", "To Inter"
    )
    
    for (kk in 1:32) {
      x.n <- night.td[,kk]
      
      x.n[x.n == 0] <- -1 # Sleep
      x.n[x.n > 0] <- -4 # Activity
      
      y.n <- rle(x.n)
      
      d_y.n <- as.data.frame(unclass(y.n))
      d_y.n$end <- cumsum(d_y.n$lengths)
      d_y.n$start <- d_y.n$end - d_y.n$lengths + 1
      
      d_y.n$values[d_y.n$lengths > 5 & d_y.n$lengths < 30] <- -2
      d_y.n$values[d_y.n$lengths > 30 & d_y.n$lengths < 60] <- -3
      
      
      s.n <- subset(d_y.n, d_y.n$values == -2)
      ss.n <- subset(d_y.n, d_y.n$values == -3)
      for (i in 1:length(s.n[,1])) {
        if(!is.na(s.n[1,1])) {
          x.n[s.n$start[i]:s.n$end[i]] <- 2 
        }
      }
      
      for (i in 1:length(ss.n[,1])) {
        if(!is.na(ss.n[1,1])) {
          x.n[ss.n$start[i]:ss.n$end[i]] <- 3
        }
      }
      
      x.n[x.n == -4] <- 4
      x.n[x.n == -1] <- 1
      
      x1.n <- factor(paste0(head(x.n,-1), tail(x.n,-1)), levels = c('12','13','14','21','23','24','31','32','34','41','42','43'))
      
      # xdf <- data.frame(x = x[1:length(x) - 1], y = x[2:length(x)])
      
      t.xdf.n <- table(x1.n)
      
      out.table.night[1,2+kk] <- (as.numeric(t.xdf.n["12"])/sum(as.numeric(t.xdf.n["12"]), as.numeric(t.xdf.n["13"]), as.numeric(t.xdf.n["14"]))) * 100
      out.table.night[2,2+kk] <- (as.numeric(t.xdf.n["13"])/sum(as.numeric(t.xdf.n["12"]), as.numeric(t.xdf.n["13"]), as.numeric(t.xdf.n["14"]))) * 100
      out.table.night[3,2+kk] <- (as.numeric(t.xdf.n["14"])/sum(as.numeric(t.xdf.n["12"]), as.numeric(t.xdf.n["13"]), as.numeric(t.xdf.n["14"]))) * 100
      
      out.table.night[4,2+kk] <- (as.numeric(t.xdf.n["21"])/sum(as.numeric(t.xdf.n["21"]), as.numeric(t.xdf.n["23"]), as.numeric(t.xdf.n["24"]))) * 100
      out.table.night[5,2+kk] <- (as.numeric(t.xdf.n["23"])/sum(as.numeric(t.xdf.n["21"]), as.numeric(t.xdf.n["23"]), as.numeric(t.xdf.n["24"]))) * 100
      out.table.night[6,2+kk] <- (as.numeric(t.xdf.n["24"])/sum(as.numeric(t.xdf.n["21"]), as.numeric(t.xdf.n["23"]), as.numeric(t.xdf.n["24"]))) * 100
      
      out.table.night[7,2+kk] <- (as.numeric(t.xdf.n["31"])/sum(as.numeric(t.xdf.n["31"]), as.numeric(t.xdf.n["32"]), as.numeric(t.xdf.n["34"]))) * 100
      out.table.night[8,2+kk] <- (as.numeric(t.xdf.n["32"])/sum(as.numeric(t.xdf.n["31"]), as.numeric(t.xdf.n["32"]), as.numeric(t.xdf.n["34"]))) * 100
      out.table.night[9,2+kk] <- (as.numeric(t.xdf.n["34"])/sum(as.numeric(t.xdf.n["31"]), as.numeric(t.xdf.n["32"]), as.numeric(t.xdf.n["34"]))) * 100
      
      out.table.night[10,2+kk] <- (as.numeric(t.xdf.n["41"])/sum(as.numeric(t.xdf.n["41"]), as.numeric(t.xdf.n["42"]), as.numeric(t.xdf.n["43"]))) * 100
      out.table.night[11,2+kk] <- (as.numeric(t.xdf.n["42"])/sum(as.numeric(t.xdf.n["41"]), as.numeric(t.xdf.n["42"]), as.numeric(t.xdf.n["43"]))) * 100
      out.table.night[12,2+kk] <- (as.numeric(t.xdf.n["43"])/sum(as.numeric(t.xdf.n["41"]), as.numeric(t.xdf.n["42"]), as.numeric(t.xdf.n["43"]))) * 100
      
      # out.table[1,2+kk] <- as.numeric(t.xdf["12"])
      # out.table[2,2+kk] <- as.numeric(t.xdf["13"])
      # out.table[3,2+kk] <- as.numeric(t.xdf["14"])
      # 
      # out.table[4,2+kk] <- as.numeric(t.xdf["21"])
      # out.table[5,2+kk] <- as.numeric(t.xdf["23"])
      # out.table[6,2+kk] <- as.numeric(t.xdf["24"])
      # 
      # out.table[7,2+kk] <- as.numeric(t.xdf["31"])
      # out.table[8,2+kk] <- as.numeric(t.xdf["32"])
      # out.table[9,2+kk] <- as.numeric(t.xdf["34"])
      # 
      # out.table[10,2+kk] <- as.numeric(t.xdf["41"])
      # out.table[11,2+kk] <- as.numeric(t.xdf["42"])
      # out.table[12,2+kk] <- as.numeric(t.xdf["43"])
      
      # out.table[1,2+kk] <- (t.xdf["1","2"]/sum(t.xdf["1",-1], na.rm = T))*100
      # out.table[2,2+kk] <- (t.xdf["1","3"]/sum(t.xdf["1",-1], na.rm = T))*100
      # out.table[3,2+kk] <- (t.xdf["1","4"]/sum(t.xdf["1",-1], na.rm = T))*100
      # 
      # out.table[4,2+kk] <- (t.xdf["2","1"]/sum(t.xdf["2",-2], na.rm = T))*100
      # out.table[5,2+kk] <- (t.xdf["2","3"]/sum(t.xdf["2",-2], na.rm = T))*100
      # out.table[6,2+kk] <- (t.xdf["2","4"]/sum(t.xdf["2",-2], na.rm = T))*100
      # 
      # out.table[7,2+kk] <- (t.xdf["3","1"]/sum(t.xdf["3",-3], na.rm = T))*100
      # out.table[8,2+kk] <- (t.xdf["3","2"]/sum(t.xdf["3",-3], na.rm = T))*100
      # out.table[9,2+kk] <- (t.xdf["3","4"]/sum(t.xdf["3",-3], na.rm = T))*100
      # 
      # out.table[10,2+kk] <- (t.xdf["4","1"]/sum(t.xdf["4",-4], na.rm = T))*100
      # out.table[11,2+kk] <- (t.xdf["4","2"]/sum(t.xdf["4",-4], na.rm = T))*100
      # out.table[12,2+kk] <- (t.xdf["4","3"]/sum(t.xdf["4",-4], na.rm = T))*100
      
    }
    
    out.table.night$mean <- rowMeans(out.table.night[,-c(1:2)], na.rm = T)
    
    nodes.night <- data.frame(
      name=c(as.character(out.table.night$From), 
             as.character(out.table.night$To)) %>% unique()
    )
    
    out.table.night$IDfrom <- match(out.table.night$From, nodes.night$name)-1
    out.table.night$IDto <- match(out.table.night$To, nodes.night$name)-1
    
    q.night <- plot_ly(
      type = "sankey",
      orientation = "h",
      domain = list(
        x = c(0.5,1),
        y = c(0,1)
      ),
      arrangement = "snap",
      node = list(
        label = as.factor(nodes.night$name),
        x = c(0.05,0.05,0.05,0.05,0.95,0.95,0.95,0.95),
        y = c(0.1,0.2,0.3,0.6,-1,-0.6,-0.2,0.5),
        color = c(
          wes_palette("FantasticFox1")[3], wes_palette("FantasticFox1")[2], wes_palette("FantasticFox1")[1], wes_palette("FantasticFox1")[5],
          wes_palette("FantasticFox1")[2], wes_palette("FantasticFox1")[1], wes_palette("FantasticFox1")[5], wes_palette("FantasticFox1")[3]
        ),
        pad = 15,
        thickness = 20,
        line = list(
          color = "black",
          width = 0.5
        )
      ),
      link = list(
        source = out.table.night$IDfrom,
        target = out.table.night$IDto,
        value = out.table.night$mean,
        color = list(
          rgb(70/255,172/255,200/255,0.2),
          rgb(70/255,172/255,200/255,0.2),
          rgb(70/255,172/255,200/255,0.2),
          rgb(226/255,210/255,0,0.2),
          rgb(226/255,210/255,0,0.2),
          rgb(226/255,210/255,0,0.2),
          rgb(221/255,141/255,41/255,0.2),
          rgb(221/255,141/255,41/255,0.2),
          rgb(221/255,141/255,41/255,0.2),
          rgb(180/255,15/255,32/255,0.2),
          rgb(180/255,15/255,32/255,0.2),
          rgb(180/255,15/255,32/255,0.2)
        )
      )
    )
    
    out.table.for.output.night <- out.table.night[,1:35]
    colnames(out.table.for.output.night) <- c(
      "Source", "Target",
      "Ch1", "Ch2", "Ch3", "Ch4", "Ch5", "Ch6", "Ch7", "Ch8",
      "Ch9", "Ch10", "Ch11", "Ch12", "Ch13", "Ch14", "Ch15", "Ch16",
      "Ch17", "Ch18", "Ch19", "Ch20", "Ch21", "Ch22", "Ch23", "Ch24",
      "Ch25", "Ch26", "Ch27", "Ch28", "Ch29", "Ch30", "Ch31", "Ch32",
      "Mean"
    )
    
    sp <- subplot(
      q.day, q.night, nrows = 1
    )%>%
      layout(
        annotations = list(
          list( 
            x = 0.25,  
            y = 1.0,  
            text = "Daytime",  
            xref = "paper",  
            yref = "paper",  
            xanchor = "center",  
            yanchor = "bottom",  
            showarrow = FALSE,
            font = list(color = 'black',
                        family = 'Arial',
                        size = 16)
          ),
          list( 
            x = 0.75,  
            y = 1.0,  
            text = "Nighttime",  
            xref = "paper",  
            yref = "paper",  
            xanchor = "center",  
            yanchor = "bottom",  
            showarrow = FALSE,
            font = list(color = 'black',
                        family = 'Arial',
                        size = 16)
          )
        )
      )
    
    
    output <- list(
      "Plot" = sp,
      "Day Data" = out.table.for.output.day,
      "Night Data" = out.table.for.output.night
    )
    
    return(output)
  }
}
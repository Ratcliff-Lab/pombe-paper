# Analyze data from Darwin: size distribution for pilot experiment
library(extrafont)
rm(list = ls())

setwd("/Users/rpineau3/Dropbox (GaTech)/BioSci-Ratcliff/Rozenn/3.PombeProject/1-Lab/8-size-distribution-Darwin/3-compare-coulter-darwin/darwin-t60-t125-t344/results")
file_names_list <- list.files()

# set variable names
mean_vec_60t <- c()
mean_vec_125t <- c()
mean_vec_344t <- c()
max_vec_60t <- c()
max_vec_125t <- c()
max_vec_344t <- c()
sd_vec_60t <- c()
sd_vec_125t <- c()
sd_vec_344t <- c()


data <- matrix(data=NA,nrow=10000,ncol=15)
# Every line and data point in the same table
for(line in seq(1,length(file_names_list))) {
  mat <- t(as.matrix(read.table(file_names_list[line], sep="\t", header=F)))
  data[1:length(mat),line] <- mat
}

summary(data)
tup <- 65 # Do not keep clusters with diameter > 65
tlow <- 6

# Build dataset for ANOVA
anova_dt <- data.frame(pop = c(), time = c(), size = c())


for(line in 1:5) {
  

  # 60t
  rm_nas <- data[!is.na(data[,line+10]),line+10] # remove NA
  abv_t <- rm_nas[rm_nas<tup & rm_nas>tlow] # above threshold
  mean_vec_60t <- c(mean_vec_60t, mean(abv_t))
  sd_vec_60t <- c(sd_vec_60t, sd(abv_t))
  
  # for ANOVA matrix
  anova_dt_60t <- abv_t
  
  
  dens <- density(abv_t)
  
  # 125t
  rm_nas <- data[!is.na(data[,line]),line] # remove NA
  abv_t <- rm_nas[rm_nas<tup & rm_nas>tlow] # above threshold
  mean_vec_125t <- c(mean_vec_125t, mean(abv_t))
  sd_vec_125t <- c(sd_vec_125t, sd(abv_t))
  
  # for ANOVA matrix
  anova_dt_125t <- abv_t
  
  # 344t
  rm_nas <- data[!is.na(data[,line+5]),line+5] # remove NA
  abv_t <- rm_nas[rm_nas<tup & rm_nas>tlow] # above threshold
  mean_vec_344t <- c(mean_vec_344t, mean(abv_t))
  sd_vec_344t <- c(sd_vec_344t, sd(abv_t))
  
  # for ANOVA matrix
  anova_dt_344t <- abv_t

  # for ANOVA matrix
  sizevec <- c(anova_dt_60t, anova_dt_125t, anova_dt_344t)
  popvec <- rep(line,length(sizevec))
  timevec <- c(rep("60t", length(anova_dt_60t)), rep("125t", length(anova_dt_125t)), rep("344t", length(anova_dt_344t)))
  anova_dt_curr <- data.frame(pop = popvec, time = timevec, size = sizevec)
  anova_dt <- rbind(anova_dt, anova_dt_curr)
  
}


mean_Darwin_60t <- mean_vec_60t
mean_Darwin_125t <- mean_vec_125t
mean_Darwin_344t <- mean_vec_344t

max_Darwin_60t <- max_vec_60t
max_Darwin_125t <- max_vec_125t
max_Darwin_344t <- max_vec_344t

sd_Darwin_60t <- sd_vec_60t
sd_Darwin_125t <- sd_vec_125t
sd_Darwin_344t <- sd_vec_344t



# Output means in a vector
pop_means_darwin <- c(mean_Darwin_60t, mean_Darwin_125t, mean_Darwin_344t)
pop_max_darwin <- c(max_Darwin_60t, max_Darwin_125t, max_Darwin_344t)
pop_sd_darwin <- c(sd_Darwin_60t, sd_Darwin_125t, sd_vec_344t)


################################################################################
#                      ANOVA for all pop at the same time                      # 
################################################################################

# Two-way ANOVA
summary(anova_dt)

# Pop and Time as factors
anova_res <- aov(anova_dt$size ~ anova_dt$pop + anova_dt$time, data = anova_dt)
summary(anova_res)

# Time only
anova_res_2 <- aov(anova_dt$size ~ anova_dt$time, data = anova_dt)
summary(anova_res_2)
TukeyHSD(anova_res_2)




################################################################################
####                           normalize by biomass                        #####
################################################################################

# Goal: to have the proportion not in terms of cluster number, but biomass importance. 

for (line in 1:dim(data)[2]) {
  
  rm_nas <- data[!is.na(data[,line]),line] # remove NAs
  abv_t <- rm_nas[rm_nas<tup & rm_nas>tlow] # above threshold
  
  # From radius, calculate volume
  volume <- 4/3*pi*(abv_t)^3
  
  # Logs volume - to work with long
  
  # Parameters
  nbins <-1000 # Decide a number of bins (high number because volume)
  
  # Calculate total biomass
  total_biomass <- sum(volume)
  
  # Calculate the size of each bin
  binsize <- (max(volume)-min(volume))/nbins
  
  # Make the vector recording bin sizes
  bin_vec <- c()
  
  for (i in 1:c(nbins+1)){ # has the min and max sizes
    bin_vec <- c(bin_vec, min(volume)+(i-1)*binsize)
  }
  
  # Find the bin for each sample
  biomass <- c()
  
  for (b in 1:nbins){
    temp_size <- c()
    
    for (s in 1:length(volume)){
      if (volume[s]>=bin_vec[b] && volume[s]<=bin_vec[b+1]){
        temp_size <- c(temp_size, (volume[s]/total_biomass))}
    }
    biomass <- c(biomass, sum(temp_size))
    
  }
  
  # Find mean bin size (otherwise distribution is skewed to the left)
  mean_bin_size <- c()
  
  for (i in 1:c(length(bin_vec)-1)){
    mean_bin_size <- c(mean_bin_size, (bin_vec[i]+bin_vec[i+1])/2)
  }
  
  # Convert back to radius
  mean_bin_size_radius <- ((mean_bin_size*3)/(4*pi))^(1/3)
  
  # Output
  newname <- paste("line", line, "norm", sep = "_")
  assign(newname, cbind(mean_bin_size_radius, biomass))
  
}


################################################################################
#                             calculate means                                  # 
################################################################################
# !! Attention, here time 60t comes last in the dataset, hence the following order

pop_means_norm <- c(mean(rep(line_11_norm[,1], line_11_norm[,2]*1000)),
                    mean(rep(line_12_norm[,1], line_12_norm[,2]*1000)),
                    mean(rep(line_13_norm[,1], line_13_norm[,2]*1000)),
                    mean(rep(line_14_norm[,1], line_14_norm[,2]*1000)),
                    mean(rep(line_15_norm[,1], line_15_norm[,2]*1000)),
                    mean(rep(line_1_norm[,1], line_1_norm[,2]*1000)),
                    mean(rep(line_2_norm[,1], line_2_norm[,2]*1000)),
                    mean(rep(line_3_norm[,1], line_3_norm[,2]*1000)),
                    mean(rep(line_4_norm[,1], line_4_norm[,2]*1000)),
                    mean(rep(line_5_norm[,1], line_5_norm[,2]*1000)),
                    mean(rep(line_6_norm[,1], line_6_norm[,2]*1000)),
                    mean(rep(line_7_norm[,1], line_7_norm[,2]*1000)),
                    mean(rep(line_8_norm[,1], line_8_norm[,2]*1000)),
                    mean(rep(line_9_norm[,1], line_9_norm[,2]*1000)),
                    mean(rep(line_10_norm[,1], line_10_norm[,2]*1000)))


################################################################################
#                    pop means as a function of time                           # 
################################################################################
# Define colors
green <- rgb(89/255, 191/255, 93/255)
yellow <- rgb(250/255,202/255,67/255)
brown <-  rgb(239/255,133/255,54/255)
purple <- rgb(141/255, 107/255, 184/255)
blue <- rgb(53/255, 134/255, 151/255)
black <-  rgb(2/255, 2/255,2/255)
pink <- "#dd1c77"
pop_cols <- rep(c(pink, yellow, purple, blue, green))

# Upload data from ancestor
anc <- read.table("/Users/rpineau3/Dropbox (GaTech)/BioSci-Ratcliff/Rozenn/3.PombeProject/1-Lab/8-size-distribution-Darwin/4-ancestor/ancestor_size_GOB360.csv", 
                  sep=",", header = TRUE)
summary(anc)
anc_size <- mean(anc$ancestor_radius_um)

time <- c(1,60,125,334)
#pop_means_c <- c(rep(4,5), pop_means_coulter)
pop_means_d <- c(rep(anc_size,5), pop_means_norm)

# Open empty plot 
plot(1, type = "n", xaxt = 'n', yaxt ='n',
     xlab = "", ylab = '', bty = "n",
     main = "", 
     ylim = c(0, 50),
     xlim = c(0, 360))

for (i in 1:5) {
  
  lines(time,
        pop_means_d[c(1,i+5,i+10,i+15)], 
        'b',
        lty = 1,
        pch = 16,
        col = pop_cols[i],
        cex = 1.7)
  
  
}


axis(1, mgp = c(0, 0.2, 0), 
     lwd.ticks = 0, cex=3)
axis(2,mgp = c(0, 0.2, 0), 
     lwd.ticks = 0)
title(ylab = expression(paste("Mean cluster radius (", mu, "m)", sep = "")), 
      mgp = c(1.5, 2, 1), cex.lab = 1.5)   # Add x-axis text
title(xlab = "Time (days)", 
      mgp = c(1.5, 2, 1), cex.lab = 1.5)    # Add y-axis text




# Increase over time
# Day 0 to day 60
mean_increase <- mean(pop_means_d[6:10]/pop_means_d[1:5]) # 4.39
sd_increase <-  sd(pop_means_d[6:10]/pop_means_d[1:5]) # 0.54

# Day 0 to day 344
mean_increase <- mean(pop_means_d[16:20]/pop_means_d[1:5]) # 5.5
sd_increase <-  sd(pop_means_d[16:20]/pop_means_d[1:5]) # 0.3

# Day 60 to day 344
mean_increase <- mean(pop_means_d[16:20]/pop_means_d[6:10]) # 1.3
sd_increase <-  sd(pop_means_d[16:20]/pop_means_d[6:10]) # 0.2

################################################################################
#               pop mean - volume - as a function of time                      # 
################################################################################


pop_means_d <- c(rep(anc_size,5), pop_means_norm)

# Open empty plot 
plot(1, type = "n", xaxt = 'n', yaxt ='n',
     xlab = "", ylab = '', bty = "n",
     main = "", 
     ylim = c(0, 4E5),
     xlim = c(0, 360))

for (i in 1:5) {
  
  lines(time,
        pop_means_d[c(1,i+5,i+10,i+15)]^3 *4/3*pi , 
        'b',
        lty = 1,
        pch = 16,
        col = pop_cols[i],
        cex = 1.7)
  
  
}

yaxis_ticks <- c(0,10E4, 20E4, 30E4, 40E4)
axis(1, mgp = c(0, 0.2, 0),
     lwd.ticks = 0, cex=3)
axis(2,mgp = c(0, 0.2, 0), at=yaxis_ticks, labels = c("0","10E4", "20E4", "30E4", "40E4"),
     lwd.ticks = 0)
title(ylab = expression(paste("Mean cluster volume (", mu, "m" ^3 , ")", sep = "")), 
      mgp = c(1.5, 2, 1), cex.lab = 1.5)   # Add x-axis text
title(xlab = "Time (days)", 
      mgp = c(1.5, 2, 1), cex.lab = 1.5)    # Add y-axis text


# Add legend
# legend(legend=c("line 1", "line 2", "line 3", "line 4", "line 5"),
#        pch=16,
#        col =pop_cols,
#        'bottomright', 
#        bty = 'n',
#        cex = 1.5)













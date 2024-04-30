# Analyze data from Darwin: size distribution for pilot experiment
library("extrafont")

# Define colors
green <- rgb(89/255, 191/255, 93/255)
yellow <- rgb(250/255,202/255,67/255)
brown <-  rgb(239/255,133/255,54/255)
purple <- rgb(141/255, 107/255, 184/255)
blue <- rgb(53/255, 134/255, 151/255)
black <-  rgb(2/255, 2/255,2/255)
pink <- "#dd1c77"
pop_cols <- rep(c(pink, yellow, purple, blue, green))


# Comparison between lines may be interesting as well to look at difference between lines

# each line has its own file : be in the folder to import tham one at a time

file_names_list <- list.files()

tup <- 100 # Do not keep clusters with diameter > 65
tlow <- 4

data <- matrix(data=NA,nrow=10000,ncol=30)
# Every line and data point in the same table
for(line in seq(1,length(file_names_list))) {
  mat <- t(as.matrix(read.table(file_names_list[line], sep="\t", header=F)))
  # Clean data
  abv_t <- mat[mat<tup & mat>tlow] # above threshold
  #mat <- mat*2 # from radius to diameter
  data[1:length(abv_t),line] <- abv_t
}

summary(data)
# Re order dataset to put in order of time and lines
# Line 1, t0 to tf, etc
line1 <- c(1,6,11,21,16)
line2 <- c(2,7,12,22,17)
line3 <- c(3,8,13,23,18)
line4 <- c(4,9,14,24,19)
line5 <- c(5,10,15,25,20)
lines <- rbind(line1,line2, line3, line4, line5)


################################################################################
#                      NUMBER OF CLUSTERS OVER TIME                            # 
################################################################################

dil <- c(10, 15, 40, 100, 200)
num_clusters <- c()
growth_rate <- c()
PDT <- c()

# Empty plot
plot(1, type = "n", xaxt = 'n', yaxt ='n',
     xlab = "", ylab = '',
     main = "", 
     xlim = c(0, 26),
     ylim = c(0, 12E5), bty="n")


# Fill plot
for(line in 1:5) {
    n_cluster_line <- c()
  for(t in 1:5){
    temp <- length(which(is.na(data[,lines[line,t]])==FALSE))
    num_clusters <- c(num_clusters, temp*dil[t])
    n_cluster_line <- c(n_cluster_line, temp*dil[t])

  }
  #points(c(0,3,6,12,24), n_cluster_line, pch = 16, cex = 1.2, col = col_vec[line])
  lines(c(0,3,6,12,24), n_cluster_line, type = "b",
        pch = 16, cex = 2, col = pop_cols[line],
        lty = 1)
  
  # calculate coefficient of regression 
  growth_rate <- c(growth_rate, log(n_cluster_line[5]/n_cluster_line[3])/18) # 24-6=18
  PDT <- c(PDT,  (18)/(3.32*(log(n_cluster_line[5])-log(n_cluster_line[3]))))
  
}

legend(legend=c("line 1", "line 2", "line 3", "line 4", "line 5"),
       pch=16,
       col = pop_cols,
       'topleft', 
       bty = 'n',
       cex = 1.3)

axis(1, mgp = c(0, 0.2, 0),
     lwd.ticks = 0,  cex=3)
axis(2,mgp = c(0, 0.2, 0), 
     lwd.ticks = 0, at=c(0,2E5,6E5,1E6), 
     labels = c(0,expression(paste("2" ^5, sep="")), expression(paste("6.10" ^5, sep="")),
                expression(paste("1.10" ^6, sep=""))))

title(ylab = "Total number of clusters", 
      mgp = c(1.5, 2, 1), cex.lab = 1.5)   # Add x-axis text
title(xlab = "Time (hours)", 
      mgp = c(1.5, 2, 1), cex.lab = 1.5)    # Add y-axis text



################################################################################
#                              GROWTH RATE                                     # 
################################################################################

# Calculate the exponential growth rate
# Use the formula to find the coefficient of the curve : (xb-xa)/(yb-ya)

growth_rate

mean(growth_rate)
sd(growth_rate)

################################################################################
#                              # DOUBLINGS                                     # 
################################################################################
# calculate the PDT from the linear part of the curve using this equation:
# (t2 - t1)/3,32 x (log n2 - log n1) where t is time and n number of cells.
PDT
mean(PDT) # this is the population doubling time - how much time it takes for the pop to double
sd(PDT)
# 2.35 

################################################################################
#                          TOTAL VOLUME OVER TIME                              # 
################################################################################

dil <- c(10, 15, 40, 100, 200)
vol_clusters <- c()

pdf(file="vol_clust_time.pdf", bg = "transparent",
    width=6, height=5, family = "Times New Roman")

# Empty plot
plot(1, type = "n", xaxt = 'n', yaxt ='n',
xlab = "", ylab = '',
main = "", 
xlim = c(0, 26),
ylim = c(0, 2.5E10))

# Fill plot
for(line in 1:5) {
  v_cluster_line <- c() # this is the vector for each line, every time point
  for(t in 1:5){
    temp <- data[which(is.na(data[,lines[line,t]])==FALSE),lines[line,t]]
    # From radius to volume
    tot_vol <- sum(4/3*pi*temp^3)
    vol_clusters <- c(vol_clusters, tot_vol*dil[t])
    v_cluster_line <- c(v_cluster_line, tot_vol*dil[t])
  }
  #points(c(0,3,6,12,24), n_cluster_line, pch = 16, cex = 1.2, col = col_vec[line])
  lines(c(0,3,6,12,24), v_cluster_line, type = "b",
        pch = 16, cex = 1.5, col = pop_cols[line],
        lty = 1.5)
}

legend(legend=c("line 1", "line 2", "line 3", "line 4", "line 5"),
       pch=16,
       col = pop_cols,
       'topleft', 
       bty = 'n',
       cex = 1)

axis(1, mgp = c(0, 0.2, 0),
     lwd.ticks = 0)
axis(2,mgp = c(0, 0.2, 0), 
     lwd.ticks = 0)
title(ylab = expression(paste("Total volume (", mu, "m3)", sep = "")), 
      mgp = c(1.5, 2, 1))   # Add x-axis text
title(xlab = "Time (hours)", 
      mgp = c(1.5, 2, 1))    # Add y-axis text





################################################################################
#                    pop size as a function of time                           # 
################################################################################




################################################################################
#                                Line 3                                        # 
################################################################################

# extract data from each line, one at a time
dataset <- data[,lines[3,]]

# extract data from times 0,3,6,12,24
h_t0 <- hist(dataset[dataset[,1]<40,1], 100)
h_t3 <- hist(dataset[dataset[,2]<40,2], 100)
h_t6 <- hist(dataset[dataset[,3]<40,3], 50)
h_t12 <- hist(dataset[dataset[,4]<40,4], 100)
h_t24 <- hist(dataset[dataset[,5]<40,5], 50)

# Peak max size
M <- 0.03 
colfunc <- colorRampPalette(c("#4d004b", "#bfd3e6"))
#colfunc <- colorRampPalette(c("#737373", "#bfd3e6"))
colornames <- colfunc(10)

#pdf(file="size_24h_line3.pdf", bg = "transparent",
 #   width=4, height=4, family = "Times New Roman")

# Empty plot
par(family="Times New Roman", cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
plot(1, type = "n", xaxt = 'n', yaxt ='n',
     xlab = "", ylab = '',
     main = "",
     ylim = c(-0.085, 0.06),
     xlim = c(0, 40), bty="n")

lines(0:120, rep(0.03, 121),lty=1,col="#d9d9d9")
lines(0:120, rep(0, 121),lty=1,col="#d9d9d9")
lines(0:120, rep(-.03, 121),lty=1,col="#d9d9d9")
lines(0:120, rep(-0.06, 121),lty=1,col="#d9d9d9")
lines(0:120, rep(-0.09, 121),lty=1,col="#d9d9d9")


meanvec <- c()

Z <- M/max(h_t0$density)
lines(h_t0$mids, h_t0$density*Z + 0.03, type = 'l', lwd = 2, col = "#484DAB", xlim=c(5,40))
lines(rep(mean(dataset[dataset[,1]<40,1], na.rm=T),2), c(0.03, 0.059), type = 'l', lwd = 1.2, col = "#ACACB6", xlim=c(5,40))
meanvec <- c(meanvec, mean(dataset[dataset[,1]<40,1], na.rm=T))

Z <- M/max(h_t3$density)
lines(h_t3$mids, h_t3$density*Z , type = 'l', lwd = 2, col = "#484DAB", xlim=c(5,40))
lines(rep(mean(dataset[dataset[,2]<40,1], na.rm=T),2), c(0, 0.028), type = 'l', lwd = 1.2, col = "#ACACB6", xlim=c(5,40))
meanvec <- c(meanvec, mean(dataset[dataset[,2]<40,1], na.rm=T))

Z <- M/max(h_t6$density)
lines(h_t6$mids, h_t6$density*Z -.03 , type = 'l', lwd = 2, col = "#484DAB", xlim=c(5,40))
lines(rep(mean(dataset[dataset[,3]<40,1], na.rm=T),2), c(-0.03, -0.007), type = 'l', lwd = 1.2, col = "#ACACB6", xlim=c(5,40))
meanvec <- c(meanvec, mean(dataset[dataset[,3]<40,1], na.rm=T))

Z <- M/max(h_t12$density)
lines(h_t12$mids, h_t12$density*Z -.06 , type = 'l', lwd = 2, col = "#484DAB", xlim=c(5,40))
lines(rep(mean(dataset[dataset[,4]<40,1], na.rm=T),2), c(-0.06, -0.035), type = 'l', lwd = 1.2, col = "#ACACB6", xlim=c(5,40))
meanvec <- c(meanvec, mean(dataset[dataset[,4]<40,1], na.rm=T))

Z <- M/max(h_t24$density)
lines(h_t24$mids, h_t24$density*Z -.09 , type = 'l', lwd = 2, col = "#484DAB", xlim=c(5,40))
lines(rep(mean(dataset[dataset[,5]<40,1], na.rm=T),2), c(-0.09, -0.063), type = 'l', lwd = 1.2, col = "#ACACB6", xlim=c(5,40))
meanvec <- c(meanvec, mean(dataset[dataset[,5]<40,1], na.rm=T))


axis(2, at = c(0.03,0,-0.03,-0.06,-0.09),
     lwd.ticks = 0 ,
     labels = c("t0", "3h", "6h", "12h", "24h" ), las =1, mgp = c(0, 0.2, 0))
axis(1,lwd.ticks = 0)
title(xlab = expression(paste("Radius (", mu, "m)", sep = "")), cex.lab = 1.5, mgp = c(2.5, 2, 1))  # Add x-axis text
      #main = "Line 3", 




#dev.off()


# Mean and standard deviation of population size
mean(meanvec)
sd(meanvec)



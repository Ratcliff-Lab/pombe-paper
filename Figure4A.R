# This script makes the mean size distribution in Figure 4 panel A.
library("readxl")
library("extrafont")

rm(list = ls()) # clear environment

################################################################################
####                           load data                                   #####
################################################################################

dt <- read_excel("raw_data.xlsx", sheet = "F4A")
summary(dt)

mean(dt$ancestor_2, na.rm = T) # mean radius size for unicellular ancestor is: 2.95 um

################################################################################
####                           normalize                                   #####
################################################################################

tup <- 65 # Do not keep clusters with radius > 65 (this is the tail of the size distribution, anything above might be to clusters too close together, or dirt)
tlow <- 6 # Do not keep clusters with radius < 6 (at least 2 cells)

for (line in 2:dim(dt)[2]) {
  
  rm_nas <- dt[!is.na(dt[,line]),line] # remove NAs
  abv_t <- rm_nas[rm_nas<tup & rm_nas>tlow] # above threshold
  
  # From radius, calculate volume
  volume <- 4/3*pi*(abv_t)^3
  
  # Logs volume - to work with long
  
  # Parameters
  nbins <- 1000 # Decide a number of bins (high number because volume)
  
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
####                             means                                     #####
################################################################################
anc <- mean(dt$ancestor_2, na.rm = T)
t60 <- c(mean(line_2_norm), mean(line_3_norm), mean(line_4_norm), mean(line_5_norm), mean(line_6_norm))
aceko <- mean(c(mean(line_7_norm), mean(line_8_norm)))

means_all <- c(anc, aceko, t60)

################################################################################
####                              plot                                     #####
################################################################################

# pdf(file="F4A.pdf", bg = "transparent",
#     width=10, height=5, family = "Times New Roman")

reps <- c(1,2,2.8,2.9,3,3.1,3.2)

plot(reps, means_all, pch = 16,
     ylab = expression(paste("Radius size (", mu, "m)", sep = "")), xlab = "", 
     xlim=c(0,5), ylim=c(0,25), cex=2.2, xaxt='n', col = "white")
axis(side =1, at = c(1,2,3), cex.axis=1.5, labels = c("uni ancestor (t0)", expression(italic('ace2')), "evolved (t334)"), col="white", col.ticks="white", col.axis="white")
axis(side =2, col="white", col.ticks="white", col.axis="white", col.lab="white",  cex.axis=1.5)


# dev.off()


################################################################################
####                       testing differences                             #####
################################################################################

anc <- dt$ancestor_2[!is.na(dt$ancestor_2)]
l_anc <- length(anc)

aceko <- dt[!is.na(dt[,7]),7]
l_aceko <- length(aceko$ace2KO_rmp2182_1)
mean(aceko$ace2KO_rmp2182_1)

t60_l3 <- dt[!is.na(dt[,4]),4]
l_t60_l3 <- length(t60_l3$t60_line3)
mean(t60_l3$t60_line3)

dt <- data.frame(lines=c(rep("anc", l_anc), rep("ace2KO", l_aceko), rep("60t", l_t60_l3 )), 
                 size=c(anc,aceko$ace2KO_rmp2182_1,t60_l3$t60_line3))

dt$lines <- as.factor(dt$lines)
summary(dt)


# Compute the analysis of variance
res.aov <- aov(size ~ lines, data = dt)
# Summary of the analysis
summary(res.aov)

# tukey's HSD
TukeyHSD(res.aov, conf.level=.95)


################################################################################
####                             volumes                                   #####
################################################################################

(mean(aceko$ace2KO_rmp2182_1)^3*pi*4/3) / (mean(anc)^3*pi*4/3)
# Volume is KO is 25 times volume of ancestor

1 - ((mean(aceko$ace2KO_rmp2182_1)^3*pi*4/3) / (mean(t60_l3$t60_line3)^3*pi*4/3))
# Volume is KO is 7% volume of evolved






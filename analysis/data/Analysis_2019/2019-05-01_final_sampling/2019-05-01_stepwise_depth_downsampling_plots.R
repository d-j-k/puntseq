## Sitewise sampling of 16S reads for PuntSeq MinION data (Plots)
## mrs72@cam.ac.uk
## Last Update - 01/05/2019 ##
##############################

library(stringr)
library(dendextend)

setwd('/Users/ms37/Desktop/PuntSeq/science/src/')


# 1. Load data
##############

load("/Users/ms37/Desktop/PuntSeq/science/results/2019-04-27_final_sampling/2019-05-01_final_sampling_kraken16S.Rdata")
load("/Users/ms37/Desktop/PuntSeq/science/results/2019-04-27_final_sampling/2019-05-01_final_sampling_kraken16S_40000X_sets_stats.Rdata")


# 2. Plotting functions
#######################

# a. Richness plots
richness.plots <- function(x, type, fixed.xmax, fixed.ymax, real.samples, title){
  
  # Define plot ranges
  if(fixed.xmax == ''){
    xlims <- c(0, max(x$April[,"Counts used"],
                      x$June[,"Counts used"],
                      x$August[,"Counts used"], na.rm = T))
  }else{
    xlims <- c(0, fixed.xmax)
  }
  
  if (type == "genus"){
    ind.richness <- 2
    ylabs <- 'Genera'
    
  }else if (type == "family"){
    ind.richness <- 6
    ylabs <- 'Families'
    
  }else if (type == "order"){
    ind.richness <- 10
    ylabs <- 'Orders'
    
  }else if (type == "class"){
    ind.richness <- 14
    ylabs <- 'Classes'
    
  }else if (type == "phylum"){
    ind.richness <- 18
    ylabs <- 'Phyla'
    
  }
  
  if(fixed.ymax == ''){
    ylims <- c(0, max(x$April[,ind.richness],
                      x$June[,ind.richness],
                      x$August[,ind.richness], na.rm = T)) 
  }else{
    ylims <- c(0, fixed.ymax)
  }
  
  # Plot
  plot(x = x$April[,"Counts used"],
       y = x$April[,ind.richness],
       ylab = ylabs, xlab = '', xaxt = 'n',
       pch = 16, type = 'p', col = 'goldenrod1', cex = 0.5,
       xlim = xlims, ylim = ylims,
       main = paste0(title, ": Richness"), cex.main = 3, cex.lab = 1.7, cex.axis = 1.3)
  
  ## June
  points(x = x$June[,"Counts used"],
         y = x$June[,ind.richness],
         pch = 16,
         col = 'darkorange',
         type = 'p', cex = 0.5)
  
  ## August
  points(x = x$August[,"Counts used"],
         y = x$August[,ind.richness],
         pch = 16,
         col = 'red',
         type = 'p', cex = 0.5)
  
  axis(side = 1, at = seq(f = 0, t = xlims[2], length.out = 5), labels = seq(f = 0, t = xlims[2], length.out = 5))
  
  ## add Michaelis-Menten values
  if(type == 'genus'){
    
    if(x$April[c(nrow(x$April)-4),1] > 40000){
      mm.april <- mmfit(x = x$April, xlims = xlims, type = type)[[1]]
    }else{
      mm.april <- c()
    }
    if(x$June[c(nrow(x$June)-4),1] > 40000){
      mm.june <- mmfit(x = x$June, xlims = xlims, type = type)[[1]]
    }else{
      mm.june <- c()
    }
    if(x$August[c(nrow(x$August)-4),1] > 40000){
      mm.august <- mmfit(x = x$August, xlims = xlims, type = type)[[1]]
    }else{
      mm.august <- c()
    }
    
  }else if(type == 'family'){
    
    if(x$April[c(nrow(x$April)-3),1] > 40000){
      mm.april <- mmfit(x = x$April, xlims = xlims, type = type)[[1]]
    }else{
      mm.april <- c()
    }
    if(x$June[c(nrow(x$June)-3),1] > 40000){
      mm.june <- mmfit(x = x$June, xlims = xlims, type = type)[[1]]
    }else{
      mm.june <- c()
    }
    if(x$August[c(nrow(x$August)-3),1] > 40000){
      mm.august <- mmfit(x = x$August, xlims = xlims, type = type)[[1]]
    }else{
      mm.august <- c()
    }
    
  }else if(type == 'order'){
    
    if(x$April[c(nrow(x$April)-2),1] > 40000){
      mm.april <- mmfit(x = x$April, xlims = xlims, type = type)[[1]]
    }else{
      mm.april <- c()
    }
    if(x$June[c(nrow(x$June)-2),1] > 40000){
      mm.june <- mmfit(x = x$June, xlims = xlims, type = type)[[1]]
    }else{
      mm.june <- c()
    }
    if(x$August[c(nrow(x$August)-2),1] > 40000){
      mm.august <- mmfit(x = x$August, xlims = xlims, type = type)[[1]]
    }else{
      mm.august <- c()
    }
    
  }else if(type == 'class'){
    
    if(x$April[c(nrow(x$April)-1),1] > 40000){
      mm.april <- mmfit(x = x$April, xlims = xlims, type = type)[[1]]
    }else{
      mm.april <- c()
    }
    if(x$June[c(nrow(x$June)-1),1] > 40000){
      mm.june <- mmfit(x = x$June, xlims = xlims, type = type)[[1]]
    }else{
      mm.june <- c()
    }
    if(x$August[c(nrow(x$August)-1),1] > 40000){
      mm.august <- mmfit(x = x$August, xlims = xlims, type = type)[[1]]
    }else{
      mm.august <- c()
    }
    
  }else if(type == 'phylum'){
    
    if(x$April[nrow(x$April),1] > 40000){
      mm.april <- mmfit(x = x$April, xlims = xlims, type = type)[[1]]
    }else{
      mm.april <- c()
    }
    if(x$June[nrow(x$June),1] > 40000){
      mm.june <- mmfit(x = x$June, xlims = xlims, type = type)[[1]]
    }else{
      mm.june <- c()
    }
    if(x$August[nrow(x$August),1] > 40000){
      mm.august <- mmfit(x = x$August, xlims = xlims, type = type)[[1]]
    }else{
      mm.august <- c()
    }
    
  }
  
  ## replot actual scores
  if(real.samples != 'n'){
    
    if (type == "genus"){
      highest.april <- nrow(x$April)-4
      highest.june <- nrow(x$June)-4
      highest.august <- nrow(x$August)-4
    }else if (type == "family"){
      highest.april <- nrow(x$April)-3
      highest.june <- nrow(x$June)-3
      highest.august <- nrow(x$August)-3
    }else if (type == "order"){
      highest.april <- nrow(x$April)-2
      highest.june <- nrow(x$June)-2
      highest.august <- nrow(x$August)-2
    }else if (type == "class"){
      highest.april <- nrow(x$April)-1
      highest.june <- nrow(x$June)-1
      highest.august <- nrow(x$August)-1
    }else if (type == "phylum"){
      highest.april <- nrow(x$April)
      highest.june <- nrow(x$June)
      highest.august <- nrow(x$August)
    }
    
    points(x = x$April[highest.april,"Counts used"],
           y = x$April[highest.april,ind.richness],
           pch = 21,
           col = 'black',
           bg = 'goldenrod1',
           type = 'p', cex = 2)
    
    points(x = x$June[highest.june,"Counts used"],
           y = x$June[highest.june,ind.richness],
           pch = 21,
           col = 'black',
           bg = 'darkorange',
           type = 'p', cex = 2)
    
    points(x = x$August[highest.august,"Counts used"],
           y = x$August[highest.august,ind.richness],
           pch = 21,
           col = 'black',
           bg = 'red',
           type = 'p', cex = 2) 
    
  }
  
  legend.text <- c()
  legend.col <- c()
  if (length(mm.april) != 0){
    legend.text <- c(legend.text, paste0('April, Rmax: ', round(mm.april[1], 0)))
    legend.col <- c(legend.col, 'goldenrod1')
  }
  if (length(mm.june) != 0){
    legend.text <- c(legend.text, paste0('June, Rmax: ', round(mm.june[1], 0)))
    legend.col <- c(legend.col, 'darkorange')
  }
  if (length(mm.august) != 0){
    legend.text <- c(legend.text, paste0('August, Rmax: ', round(mm.august[1], 0)))
    legend.col <- c(legend.col, 'red')
  }
  if (length(legend.text) != 0){
    legend("topleft", legend = legend.text, lwd = 2, 
           pch = 16, col = legend.col, bty = 'n', cex = 1.2)  
  }
  
}
richness.plots.spec <- function(x, type, fixed.xmax, fixed.ymax, real.samples, col){
  
  # Define plot ranges
  if(fixed.xmax == ''){
    xlims <- c(0, max(x[,"Counts used"], na.rm = T))
  }else{
    xlims <- c(0, fixed.xmax)
  }
  
  if (type == "genus"){
    ind.richness <- 2
    ylabs <- 'Taxonomic Genera'
    
  }else if (type == "family"){
    ind.richness <- 6
    ylabs <- 'Taxonomic Families'
    
  }else if (type == "order"){
    ind.richness <- 10
    ylabs <- 'Taxonomic Orders'
    
  }else if (type == "class"){
    ind.richness <- 14
    ylabs <- 'Taxonomic Classes'
    
  }else if (type == "phylum"){
    ind.richness <- 18
    ylabs <- 'Taxonomic Phyla'
    
  }
  
  if(fixed.ymax == ''){
    ylims <- c(0, max(x[,ind.richness], na.rm = T)) 
  }else{
    ylims <- c(0, fixed.ymax)
  }
  
  # Plot
  plot(x = x[,"Counts used"],
       y = x[,ind.richness],
       ylab = ylabs, xlab = 'Sampling Depth', xaxt = 'n',
       pch = 16, type = 'p', col = col, cex = 0.5,
       xlim = xlims, ylim = ylims,
       main = "", cex.main = 2, 
       cex.lab = 1.5, cex.axis = 1.2)
  
  axis(side = 1, at = seq(f = 0, t = xlims[2], length.out = 6), labels = seq(f = 0, t = xlims[2], length.out = 6))
  
  ## replot actual scores
  if(real.samples != 'n'){
    
    if (type == "genus"){
      highest.april <- nrow(x$April)-4
      highest.june <- nrow(x$June)-4
      highest.august <- nrow(x$August)-4
    }else if (type == "family"){
      highest.april <- nrow(x$April)-3
      highest.june <- nrow(x$June)-3
      highest.august <- nrow(x$August)-3
    }else if (type == "order"){
      highest.april <- nrow(x$April)-2
      highest.june <- nrow(x$June)-2
      highest.august <- nrow(x$August)-2
    }else if (type == "class"){
      highest.april <- nrow(x$April)-1
      highest.june <- nrow(x$June)-1
      highest.august <- nrow(x$August)-1
    }else if (type == "phylum"){
      highest.april <- nrow(x$April)
      highest.june <- nrow(x$June)
      highest.august <- nrow(x$August)
    }
    
    points(x = x$April[highest.april,"Counts used"],
           y = x$April[highest.april,ind.richness],
           pch = 21,
           col = 'black',
           bg = 'goldenrod1',
           type = 'p', cex = 2)
    
    points(x = x$June[highest.june,"Counts used"],
           y = x$June[highest.june,ind.richness],
           pch = 21,
           col = 'black',
           bg = 'darkorange',
           type = 'p', cex = 2)
    
    points(x = x$August[highest.august,"Counts used"],
           y = x$August[highest.august,ind.richness],
           pch = 21,
           col = 'black',
           bg = 'red',
           type = 'p', cex = 2) 
    
  }
  
}

# b. Shannon evenness plots
shannon.plots <- function(x, type, fixed.xmax, real.samples, title){
  
  # Define plot ranges
  if(fixed.xmax == ''){
    xlims <- c(0, max(x$April[,"Counts used"],
                      x$June[,"Counts used"],
                      x$August[,"Counts used"], na.rm = T))
  }else{
    xlims <- c(0, fixed.xmax)
  }
  
  if (type == "genus"){
    ind.shannon <- 4
    ylabs <- 'Genera'
    
  }else if (type == "family"){
    ind.shannon <- 8
    ylabs <- 'Families'
    
  }else if (type == "order"){
    ind.shannon <- 12
    ylabs <- 'Orders'
    
  }else if (type == "class"){
    ind.shannon <- 16
    ylabs <- 'Classes'
    
  }else if (type == "phylum"){
    ind.shannon <- 20
    ylabs <- 'Phyla'
    
  }
  ylims <- c(0, 1)
  
  # Plot
  plot(x = x$April[,"Counts used"],
       y = x$April[,ind.shannon],
       ylab = '', xlab = '', xaxt = 'n',
       pch = 16, type = 'p', col = 'goldenrod1', cex = 0.5,
       xlim = xlims, ylim = ylims,
       main = paste0(title, ": Shannon Evenness"), cex.main = 3, 
       cex.lab = 1.7, cex.axis = 1.3)
  
  ## June
  points(x = x$June[,"Counts used"],
         y = x$June[,ind.shannon],
         pch = 16,
         col = 'darkorange',
         type = 'p', cex = 0.5)
  
  ## August
  points(x = x$August[,"Counts used"],
         y = x$August[,ind.shannon],
         pch = 16,
         col = 'red',
         type = 'p', cex = 0.5)
  
  axis(side = 1, at = seq(f = 0, t = xlims[2], length.out = 5), labels = seq(f = 0, t = xlims[2], length.out = 5))
  
  ## replot actual scores
  if(real.samples != 'n'){
    
    if (type == "genus"){
      highest.april <- nrow(x$April)-4
      highest.june <- nrow(x$June)-4
      highest.august <- nrow(x$August)-4
    }else if (type == "family"){
      highest.april <- nrow(x$April)-3
      highest.june <- nrow(x$June)-3
      highest.august <- nrow(x$August)-3
    }else if (type == "order"){
      highest.april <- nrow(x$April)-2
      highest.june <- nrow(x$June)-2
      highest.august <- nrow(x$August)-2
    }else if (type == "class"){
      highest.april <- nrow(x$April)-1
      highest.june <- nrow(x$June)-1
      highest.august <- nrow(x$August)-1
    }else if (type == "phylum"){
      highest.april <- nrow(x$April)
      highest.june <- nrow(x$June)
      highest.august <- nrow(x$August)
    }
    
    points(x = x$April[highest.april,"Counts used"],
           y = x$April[highest.april,ind.shannon],
           pch = 21,
           col = 'black',
           bg = 'goldenrod1',
           type = 'p', cex = 2)
    
    points(x = x$June[highest.june,"Counts used"],
           y = x$June[highest.june,ind.shannon],
           pch = 21,
           col = 'black',
           bg = 'darkorange',
           type = 'p', cex = 2)
    
    points(x = x$August[highest.august,"Counts used"],
           y = x$August[highest.august,ind.shannon],
           pch = 21,
           col = 'black',
           bg = 'red',
           type = 'p', cex = 2) 
    
  }
  
#   legend("bottomright", legend = c('April', 'June', 'August'),
#          lwd = 2, pch = 16, col = c('goldenrod1', 'darkorange', 'red'), bty = 'n', cex = 1.5)
  
}

# c. Simpson's D evenness plots
simpson.plots <- function(x, type, fixed.xmax, real.samples, title){
  
  # Define plot ranges
  if(fixed.xmax == ''){
    xlims <- c(0, max(x$April[,"Counts used"],
                      x$June[,"Counts used"],
                      x$August[,"Counts used"], na.rm = T))
  }else{
    xlims <- c(0, fixed.xmax)
  }
  
  if (type == "genus"){
    ind.simpson <- 5
  }else if (type == "family"){
    ind.simpson <- 9
  }else if (type == "order"){
    ind.simpson <- 13
  }else if (type == "class"){
    ind.simpson <- 17
  }else if (type == "phylum"){
    ind.simpson <- 21
  }
  ylims <- c(0, 1)
  
  # Plot
  plot(x = x$April[,"Counts used"],
       y = x$April[,ind.simpson],
       ylab = '', xlab = '', xaxt = 'n',
       pch = 16, type = 'p', col = 'goldenrod1', cex = 0.5,
       xlim = xlims, ylim = ylims,
       main = paste0(title, ": Simpson's D"), cex.main = 3, 
       cex.lab = 1.7, cex.axis = 1.3)
  
  ## June
  points(x = x$June[,"Counts used"],
         y = x$June[,ind.simpson],
         pch = 16,
         col = 'darkorange',
         type = 'p', cex = 0.5)
  
  ## August
  points(x = x$August[,"Counts used"],
         y = x$August[,ind.simpson],
         pch = 16,
         col = 'red',
         type = 'p', cex = 0.5)
  
  axis(side = 1, at = seq(f = 0, t = xlims[2], length.out = 5), labels = seq(f = 0, t = xlims[2], length.out = 5))
  
  ## replot actual scores
  if(real.samples != 'n'){
    
    if (type == "genus"){
      highest.april <- nrow(x$April)-4
      highest.june <- nrow(x$June)-4
      highest.august <- nrow(x$August)-4
    }else if (type == "family"){
      highest.april <- nrow(x$April)-3
      highest.june <- nrow(x$June)-3
      highest.august <- nrow(x$August)-3
    }else if (type == "order"){
      highest.april <- nrow(x$April)-2
      highest.june <- nrow(x$June)-2
      highest.august <- nrow(x$August)-2
    }else if (type == "class"){
      highest.april <- nrow(x$April)-1
      highest.june <- nrow(x$June)-1
      highest.august <- nrow(x$August)-1
    }else if (type == "phylum"){
      highest.april <- nrow(x$April)
      highest.june <- nrow(x$June)
      highest.august <- nrow(x$August)
    }
    
    points(x = x$April[highest.april,"Counts used"],
           y = x$April[highest.april,ind.simpson],
           pch = 21,
           col = 'black',
           bg = 'goldenrod1',
           type = 'p', cex = 2)
    
    points(x = x$June[highest.june,"Counts used"],
           y = x$June[highest.june,ind.simpson],
           pch = 21,
           col = 'black',
           bg = 'darkorange',
           type = 'p', cex = 2)
    
    points(x = x$August[highest.august,"Counts used"],
           y = x$August[highest.august,ind.simpson],
           pch = 21,
           col = 'black',
           bg = 'red',
           type = 'p', cex = 2) 
    
  }
  
#   legend("bottomright", legend = c('April', 'June', 'August'),
#          lwd = 2, pch = 16, col = c('goldenrod1', 'darkorange', 'red'), bty = 'n', cex = 1.5)
  
}


# 3. Plot: Downsampling effects on each taxonomic level
#######################################################

pdf('../results/2019-04-27_final_sampling/plots/2019-05-01_final_sampling_kraken16S_downsampling_richness_evenness.pdf', 
    width = 18, height = 7)
par(mfcol = c(1,2))
for(type in c('genus', 'family', 'order', 'class', 'phylum')){
  
  if(type == 'genus'){
    fixed.ymax = 1200
  }else if(type == 'family'){
    fixed.ymax = 400
  }else if(type == 'order'){
    fixed.ymax = 250
  }else if(type == 'class'){
    fixed.ymax = 100
  }else if(type == 'phylum'){
    fixed.ymax = 50
  }
  
  for (i in 1:12){
    par(mar = c(4, 6, 4, 3))
    richness.plots(x = sampling[[i]], type = type, fixed.xmax = 80000, fixed.ymax = fixed.ymax, real.samples = 'n',
                   title = names(sampling)[i])
    abline(v = 40000, lty = 2, lwd = 0.5)
    shannon.plots(x = sampling[[i]], type = type, fixed.xmax = 80000, real.samples = 'n',
                  title = names(sampling)[i])
    abline(v = 40000, lty = 2, lwd = 0.5)
  }  
}
dev.off()


# 4. Plot: Sequencing depth per sample on each taxonomic level
##############################################################

cols.all <- c(rep('#2157A4', 3), rep('#3694D1', 3), rep('#65C6E8', 3),
                   rep('#9DD7ED', 3), rep('#F7EC73', 3), rep('#FDCB44', 3),
                   rep('#F1861E', 3), rep('#E63A11', 3), rep('#D61015', 3), 
                   rep('#D61015', 3), rep('grey', 3), rep('black', 3))

pdf('../results/2019-04-27_final_sampling/plots/2019-05-01_final_sampling_kraken16S_depth_cutoff_40000X.pdf', 
    width = 8, height = 6)
barplot(summary.genus[,1], las = 2, col = cols.all, border = cols.all, ylim = c(0,1000000), 
        main = 'Original read depth (Genus)', cex.main = 2, cex.names = 0.7)
abline(h = 40000, lty = 2)
barplot(summary.family[,1], las = 2, col = cols.all, border = cols.all, ylim = c(0,1000000), 
        main = 'Original read depth (Family)', cex.main = 2, cex.names = 0.7)
abline(h = 40000, lty = 2)
barplot(summary.order[,1], las = 2, col = cols.all, border = cols.all, ylim = c(0,1000000), 
        main = 'Original read depth (Order)', cex.main = 2, cex.names = 0.7)
abline(h = 40000, lty = 2)
barplot(summary.class[,1], las = 2, col = cols.all, border = cols.all, ylim = c(0,1000000), 
        main = 'Original read depth (Class)', cex.main = 2, cex.names = 0.7)
abline(h = 40000, lty = 2)
barplot(summary.phylum[,1], las = 2, col = cols.all, border = cols.all, ylim = c(0,1000000),
        main = 'Original read depth (Phylum)', cex.main = 2, cex.names = 0.7)
abline(h = 40000, lty = 2)
dev.off()


# 5. Plot: correlation of original sequencing depth vs. richness at 40000X
##########################################################################

pdf('../results/2019-04-27_final_sampling/plots/2019-05-01_final_sampling_kraken16S_richness_at_40000X_vs_sampling_depth.pdf', 
    width = 8, height = 6)

## Genus
plot(x = summary.genus[1:30,1],
     y = summary.genus[1:30,3], 
     pch = 16, 
     xlab = 'Original sampling Depth',
     ylab = 'Genus richness at 40,000X',
     xlim = c(0,1000000),
     ylim = c(0,1200), col = cols.all[1:30], 
     main = 'Genus richness vs. 16S read depth', cex.main = 2)
pearson.cor <- round(cor(summary.genus[1:30,1], summary.genus[1:30,3], use = 'complete.obs'),3)
lm.out <- summary(lm(summary.genus[1:30,3] ~ summary.genus[1:30,1]))
pval <- round(lm.out[[4]][2,4],3)
abline(coef = lm.out[[4]][,1], lty = 2)
legend('topright', c(paste('Pearson correlation: ', pearson.cor),
                     paste('LM p-value: ', pval)))

## Family
plot(x = summary.family[1:30,1],
     y = summary.family[1:30,3], 
     pch = 16, 
     xlab = 'Original sampling Depth',
     ylab = 'Family richness at 40,000X',
     xlim = c(0,1000000),
     ylim = c(0,400), col = cols.all[1:30], 
     main = 'Family richness vs. 16S read depth', cex.main = 2)
pearson.cor <- round(cor(summary.family[1:30,1], summary.family[1:30,3], use = 'complete.obs'),3)
lm.out <- summary(lm(summary.family[1:30,3] ~ summary.family[1:30,1]))
pval <- round(lm.out[[4]][2,4],3)
abline(coef = lm.out[[4]][,1], lty = 2)
legend('topright', c(paste('Pearson correlation: ', pearson.cor),
                     paste('LM p-value: ', pval)))

## Order
plot(x = summary.order[1:30,1],
     y = summary.order[1:30,3], 
     pch = 16, 
     xlab = 'Original sampling Depth',
     ylab = 'Order richness at 40,000X',
     xlim = c(0,1000000),
     ylim = c(0,250), col = cols.all[1:30], 
     main = 'Order richness vs. 16S read depth', cex.main = 2)
pearson.cor <- round(cor(summary.order[1:30,1], summary.order[1:30,3], use = 'complete.obs'),3)
lm.out <- summary(lm(summary.order[1:30,3] ~ summary.order[1:30,1]))
pval <- round(lm.out[[4]][2,4],3)
abline(coef = lm.out[[4]][,1], lty = 2)
legend('topright', c(paste('Pearson correlation: ', pearson.cor),
                     paste('LM p-value: ', pval)))

## Class
plot(x = summary.class[1:30,1],
     y = summary.class[1:30,3], 
     pch = 16, 
     xlab = 'Original sampling Depth',
     ylab = 'Class richness at 40,000X',
     xlim = c(0,1000000),
     ylim = c(0,100), col = cols.all[1:30], 
     main = 'Class richness vs. 16S read depth', cex.main = 2)
pearson.cor <- round(cor(summary.class[1:30,1], summary.class[1:30,3], use = 'complete.obs'),3)
lm.out <- summary(lm(summary.class[1:30,3] ~ summary.class[1:30,1]))
pval <- round(lm.out[[4]][2,4],3)
abline(coef = lm.out[[4]][,1], lty = 2)
legend('topright', c(paste('Pearson correlation: ', pearson.cor),
                     paste('LM p-value: ', pval)))

## Phylum
plot(x = summary.phylum[1:30,1],
     y = summary.phylum[1:30,3], 
     pch = 16, 
     xlab = 'Original sampling Depth',
     ylab = 'Phylum richness at 40,000X',
     xlim = c(0,1000000),
     ylim = c(0,50), col = cols.all[1:30], 
     main = 'Phylum richness vs. 16S read depth', cex.main = 2)
pearson.cor <- round(cor(summary.phylum[1:30,1], summary.phylum[1:30,3], use = 'complete.obs'),3)
lm.out <- summary(lm(summary.phylum[1:30,3] ~ summary.phylum[1:30,1]))
pval <- round(lm.out[[4]][2,4],3)
abline(coef = lm.out[[4]][,1], lty = 2)
legend('topright', c(paste('Pearson correlation: ', pearson.cor),
                     paste('LM p-value: ', pval)))

dev.off()


# 6. Plot: Supplementary Figure 1
#################################

pdf('../results/2019-04-27_final_sampling/plots/2019-05-01_final_sampling_kraken16S_Supplementary_Figure_1.pdf',
    width = 18, height = 7)
par(mfrow = c(1,2))
par(mar = c(7, 6, 5, 2))

## part A
richness.plots.spec(x = sampling[[1]]$April, type = 'family', 
                    fixed.xmax = 100000, fixed.ymax = 250, real.samples = 'n', col = '#2157A4')
abline(v = 40000, lty = 2, lwd = 0.5)

## part B
sampling.cols <- c('#2157A4', '#3694D1', '#65C6E8', '#9DD7ED', '#F7EC73', 
                   '#FDCB44', '#F1861E', '#E63A11', '#D61015', '#D61015')
names(sampling.cols) <- c(1, 2, 3, 4, 5, 6, 7, 8, 9.1, 9.2)

# fetch richnesses at 40000X
tmp.rich <- cbind(family.40000$`1`$April$stats[,2],
                  family.40000$`1`$June$stats[,2],
                  family.40000$`1`$August$stats[,2],
                  family.40000$`2`$April$stats[,2],
                  family.40000$`2`$August$stats[,2], #
                  family.40000$`3`$April$stats[,2],
                  family.40000$`4`$April$stats[,2],
                  family.40000$`4`$August$stats[,2], #
                  family.40000$`5`$April$stats[,2],
                  family.40000$`5`$August$stats[,2],
                  family.40000$`6`$June$stats[,2],
                  family.40000$`6`$August$stats[,2],
                  family.40000$`7`$April$stats[,2],
                  family.40000$`8`$April$stats[,2],
                  family.40000$`9.1`$August$stats[,2],
                  family.40000$`9.2`$August$stats[,2])
colnames(tmp.rich) <- c('1 - April', '1 - June', '1 - August', '2 - April', '2 - August', '3 - April',
                        '4 - April', '4 - August', '5 - April', '5 - August', '6 - June', '6 - August',
                        '7 - April', '8 - April', '9.1 - August', '9.2 - August')

# fetch original richnesses
family.richness.orig <- summary.family[match(colnames(tmp.rich), rownames(summary.family)),2]

# build ratios
for (i in 1:ncol(tmp.rich)){
  tmp.rich[,i] <- 100*c(tmp.rich[,i]/family.richness.orig[i])
}

# match colors in plot
col.match <- str_split_fixed(colnames(tmp.rich), ' - ', 2)[,1]
cols <- as.character(sampling.cols[match(col.match,names(sampling.cols))])

boxplot(tmp.rich,
        main = "", 
        cex.main = 2, ylab = 'Richness at 40,000 X [%]', ylim = c(0, 100), 
        cex.lab = 1.5, cex.axis = 1, 
        pch = 16, cex = 0.3, notch = T, las = 2, yaxt = 'n',
        col = cols, 
        border = cols)
axis(2, at = seq(f = 0, t = 100, by = 20), las = 3, cex.axis = 1.2)
dev.off()


# 7. Clustering heatmap & PCA
#############################

# A. PCA
bacterial.heatmap <- function(x){
  
  # determine cols
  heat.cols <- colorRampPalette(c("white", "darkgreen"))(n = 100)
  
  # pre-cluster
  x.plot <- x
  
  samples <- x.plot %>% t %>% dist %>% hclust %>% as.dendrogram %>%
    set("branches_k_color", value = 'black', k = 1) %>% 
    set("branches_lwd", 0.5)
  
  taxa <- x.plot %>% dist %>% hclust %>% as.dendrogram %>%
    set("branches_k_color", value = 'black', k = 1)
  
  heatmap(x = x.plot, 
          Rowv = taxa, 
          Colv = samples, 
          col = heat.cols,
          na.rm = T,
          margins = c(10, 15))

}

## summarise data
tmp.family <- sort(unique(c(as.character(family.40000$`1`$April$sets$name),
                            as.character(family.40000$`1`$June$sets$name),
                            as.character(family.40000$`1`$August$sets$name),
                            as.character(family.40000$`2`$April$sets$name),
                            as.character(family.40000$`2`$August$sets$name),
                            as.character(family.40000$`3`$April$sets$name),
                            as.character(family.40000$`4`$April$sets$name),
                            as.character(family.40000$`4`$August$sets$name),
                            as.character(family.40000$`5`$April$sets$name),
                            as.character(family.40000$`5`$August$sets$name),
                            as.character(family.40000$`6`$June$sets$name),
                            as.character(family.40000$`6`$August$sets$name),
                            as.character(family.40000$`7`$April$sets$name),
                            as.character(family.40000$`8`$April$sets$name),
                            as.character(family.40000$`9.1`$August$sets$name),
                            as.character(family.40000$`9.2`$August$sets$name),
                            as.character(family.40000$`P`$April$sets$name),
                            as.character(family.40000$`P`$June$sets$name),
                            as.character(family.40000$`P`$August$sets$name))))
family.heatmap <- matrix(NA, ncol = 19, nrow = length(tmp.family))
rownames(family.heatmap) <- tmp.family
colnames(family.heatmap) <- c('1 - April', '1 - June', '1 - August', '2 - April', '2 - August', '3 - April',
                              '4 - April', '4 - August', '5 - April', '5 - August', '6 - June', '6 - August',
                              '7 - April', '8 - April', '9.1 - August', '9.2 - August', 'P - April', 'P - June',
                              'P - August')
family.heatmap[,1] <- family.40000$`1`$April$set[match(rownames(family.heatmap), family.40000$`1`$April$sets$name),3]
family.heatmap[,2] <- family.40000$`1`$June$set[match(rownames(family.heatmap), family.40000$`1`$June$sets$name),3]
family.heatmap[,3] <- family.40000$`1`$August$set[match(rownames(family.heatmap), family.40000$`1`$August$sets$name),3]
family.heatmap[,4] <- family.40000$`2`$April$set[match(rownames(family.heatmap), family.40000$`2`$April$sets$name),3]
family.heatmap[,5] <- family.40000$`2`$August$set[match(rownames(family.heatmap), family.40000$`2`$August$sets$name),3]
family.heatmap[,6] <- family.40000$`3`$April$set[match(rownames(family.heatmap), family.40000$`3`$April$sets$name),3]
family.heatmap[,7] <- family.40000$`4`$April$set[match(rownames(family.heatmap), family.40000$`4`$April$sets$name),3]
family.heatmap[,8] <- family.40000$`4`$August$set[match(rownames(family.heatmap), family.40000$`4`$August$sets$name),3]
family.heatmap[,9] <- family.40000$`5`$April$set[match(rownames(family.heatmap), family.40000$`5`$April$sets$name),3]
family.heatmap[,10] <- family.40000$`5`$August$set[match(rownames(family.heatmap), family.40000$`5`$August$sets$name),3]
family.heatmap[,11] <- family.40000$`6`$June$set[match(rownames(family.heatmap), family.40000$`6`$June$sets$name),3]
family.heatmap[,12] <- family.40000$`6`$August$set[match(rownames(family.heatmap), family.40000$`6`$August$sets$name),3]
family.heatmap[,13] <- family.40000$`7`$April$set[match(rownames(family.heatmap), family.40000$`7`$April$sets$name),3]
family.heatmap[,14] <- family.40000$`8`$April$set[match(rownames(family.heatmap), family.40000$`8`$April$sets$name),3]
family.heatmap[,15] <- family.40000$`9.1`$August$set[match(rownames(family.heatmap), family.40000$`9.1`$August$sets$name),3]
family.heatmap[,16] <- family.40000$`9.2`$August$set[match(rownames(family.heatmap), family.40000$`9.2`$August$sets$name),3]
family.heatmap[,17] <- family.40000$`P`$April$set[match(rownames(family.heatmap), family.40000$`P`$April$sets$name),3]
family.heatmap[,18] <- family.40000$`P`$June$set[match(rownames(family.heatmap), family.40000$`P`$June$sets$name),3]
family.heatmap[,19] <- family.40000$`P`$August$set[match(rownames(family.heatmap), family.40000$`P`$August$sets$name),3]

## log-transform
family.heatmap[is.na(family.heatmap)] <- 0
family.heatmap <- family.heatmap[!apply(family.heatmap, 1, function(x) {all(x == 0)}),]
family.heatmap <- log10(family.heatmap+1)

pdf('../results/2019-04-27_final_sampling/plots/2019-05-01_final_sampling_kraken16S_heatmap_fixed_40000X_family.pdf',
    width=40, height=40)
mar.default <- c(5,5,5,5) + 0.5
par(mar = mar.default)
bacterial.heatmap(family.heatmap)
dev.off()

# B. PCA
bacterial.PCA <- function(x, title){
  
  # determine cols
  sampling.cols <- c('#2157A4', '#3694D1', '#65C6E8', '#9DD7ED', '#F7EC73', 
                     '#FDCB44', '#F1861E', '#E63A11', '#D61015', '#D61015', 'black')
  names(sampling.cols) <- c(1, 2, 3, 4, 5, 6, 7, 8, 9.1, 9.2, 'P')
  
  # calibrate
  x <- x[rowSums(x) != 0, ]
  y.PCA <- prcomp(t(x), 
                  scale=T, 
                  center=T)
  sum.PCA <- summary(y.PCA)
  y.PCA.components.vars <- sum.PCA$importance[2,1:2]*100
  y.PCA.components.vars <- round(y.PCA.components.vars, 1)
  y.PCA.components <- y.PCA$x[,c(1,2)]
  x.axis.lengths <- round(range(y.PCA.components[,1]), 1)
  x.axis.lengths[1] <- x.axis.lengths[1]-1; x.axis.lengths[2] <- x.axis.lengths[2]+1; 
  y.axis.lengths <- round(range(y.PCA.components[,2]), 1)
  y.axis.lengths[1] <- y.axis.lengths[1]-1; y.axis.lengths[2] <- y.axis.lengths[2]+1; 
  
  # plot
  
  # basic PCA plot
  col.match <- str_split_fixed(names(y.PCA.components[,1]), ' - ', 2)[,1]
  
  plot(y.PCA.components[,1], 
       y.PCA.components[,2], 
       col = 'white', 
       pch = 16, cex = 1.5, 
       xlab = paste0("PC1 [", y.PCA.components.vars[1], " % var.]"), 
       ylab = paste0("PC2 [", y.PCA.components.vars[2], " % var.]"),
       main = title, cex.main = 3, cex.lab=3,
       xlim = x.axis.lengths, 
       ylim = y.axis.lengths,
       xaxt = 'n', yaxt='n')
  
  # samples = 'name'
  samples.names.tmp <- colnames(x)[i]
  samples.names.tmp <- str_split_fixed(samples.names.tmp, ' - ', 2)[,1]
  samples.names.tmp.ind <- grep(samples.names.tmp, rownames(y.PCA.components))
  text(x = y.PCA.components[,1], 
       y = y.PCA.components[,2],
       labels = colnames(x), 
       cex = 1, 
       col = sampling.cols[match(col.match, names(sampling.cols))])
}

pdf('../results/2019-04-27_final_sampling/plots/2019-05-01_final_sampling_kraken16S_PCA_fixed_40000X_family.pdf', 
    width=15, height=15)
mar.default <- c(5,4,4,2) + 0.5
par(mar = mar.default + c(0, 2, 0, 0))
bacterial.PCA(x = family.heatmap, title = 'Family PCA at 40,000X')
dev.off()



# 8. Signatures Extractions
###########################

load('2019-02-14_downsampling_centrifuge_30000_reads.Rdata')

# a. Latest sigfit version
library(rstan)
library(sigfit)
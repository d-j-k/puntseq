## Sitewise sampling of 16S reads for PuntSeq MinION data (Processing)
## mrs72@cam.ac.uk
## Last Update - 01/05/2019 ##
##############################

library(stringr)
library(dendextend)

setwd('/Users/ms37/Desktop/PuntSeq/science/src/')


# 1. Load data
##############

# a. 'Final' classifications with Kraken vs. 16S
classifications <- vector(mode = 'list', length = 3)
names(classifications) <- c("April", "June", "August")
classifications$April <- read.table('../data/2019-04_final_classification/April_kraken16S.tsv', 
                                    sep = '\t', check.names = F, header = T, fill = T)
classifications$June <- read.table('../data/2019-04_final_classification/June_kraken16S.tsv', 
                                    sep = '\t', check.names = F, header = T, fill = T)
classifications$August <- read.table('../data/2019-04_final_classification/August_kraken16S.tsv', 
                                    sep = '\t', check.names = F, header = T, fill = T)

# b. Relable data, in accordance with figures
classifications$April <- classifications$April[,c(1:4,30,c(11:9,7:3,1,2,12,8)+4)]
classifications$June <- classifications$June[,c(1:4,28,c(10:8,6:3,1,2,11,7)+4)]
classifications$June <- cbind(classifications$June[,1:8],c(93,rep(NA, c(nrow(classifications$June)-1))),classifications$June[,9:ncol(classifications$June)])
classifications$August <- classifications$August[,c(1:4,30,c(11:9,7:3,1,2,12,8)+4)]
colnames(classifications$April)[6:17] <- colnames(classifications$June)[6:17] <- colnames(classifications$August)[6:17] <- c('1', '2', '3', '4', '5', '6', '7', '8', '9.1', '9.2', 'N', 'P')


# 2. Sampling functions
#######################

# a. Sample
sample_classifications <- function(x, n_rep, read_steps){
  
  # remove NAs (non-existent species/genera/families/phyla)
  if(nrow(x) > 0){
    
    # subsample multinomially for each
    sampling.x <-  do.call(cbind, lapply(read_steps, function(y) rmultinom(n = n_rep, size = y, prob = x[,3])))
    sampling.x.sets <- cbind(x[,1:2],sampling.x)
    sampling.x.sets <- sampling.x.sets[!apply(sampling.x, 1, function(x){all(x == 0)}),]
    
    # calculate sats
    x <- cbind(x, sampling.x)
    summary.x <- t(rbind(mapply(function(ind) {get_richness_shannon(x[,c(1,ind)])}, ind = 4:ncol(x)),
                         mapply(function(ind) {get_simpson(x[,c(1,ind)])}, ind = 4:ncol(x))))
    
    # output
    summary.x <- matrix(unlist(summary.x), nrow = nrow(summary.x), ncol = 4)
    return(list(stats = summary.x, sets = sampling.x.sets))
    
  }else{
    
    # output
    return(matrix(NA, nrow = length(read_steps), ncol = 4))
    
  }
  
}

# b. Calculate sample richness, Shannon diversity index and evenness
get_richness_shannon <- function(x){
  
  # a. remove NAs/0s (if present)
  if(any(x[,2] == 0 | is.na(x[,2]))){
    x <- x[-which(x[,2] == 0 | is.na(x[,2])),]
  }
  
  # b. individuals in pool
  total_N <- sum(as.numeric(as.character(x[,2])))
  
  # c. type richness
  richness <- nrow(x)
  
  # d. calculate shannon index: iterate over each species/genus/order/phylum
  shannon_H <- -sum(apply(x, 1, function(y) {p_i <- as.numeric(y[2])/total_N; shannon_H <- p_i * log(p_i); return(shannon_H)}))
  
  # e. calculate maximum possible Shannon H
  shannon_H_max <- log(richness)
  
  # f. calculate Shannon E (evenness)
  shannon_E <- shannon_H/shannon_H_max
  
  # g. output
  out <- list("Richness" = richness,
              "Shannon H" = shannon_H,
              "Shannon E" = shannon_E)
  return(out)
}

# c. Calculate Simpson's D
get_simpson <- function(x){
  
  # a. remove NAs/0s (if present)
  if(any(x[,2] == 0 | is.na(x[,2]))){
    x <- x[-which(x[,2] == 0 | is.na(x[,2])),]
  }
  
  # b. individuals in pool
  total_N <- sum(as.numeric(as.character(x[,2])))
  
  # c. type richness
  richness <- nrow(x)
  
  # d. Simpson's D
  ## adding a pseudocount to account for species presence of 1 read (!)
  simpson_D <- 1 - sum(apply(x, 1, function(y) {p_i <- c(c(as.numeric(y[2])+1)*as.numeric(y[2]))/c(c(total_N+1)*c(total_N)); return(p_i)}))
  
  # e. output
  out <- list("Simpson D" = simpson_D)
  return(out)
}

# d. Generate Michaelis Menten fit do study downsampling effects
mmfit <- function(x, xlims, type){
  
  # take values, dose = counts, response = richness
  if (type == 'genus'){
    x.counts <- x[-c(nrow(x)-3,nrow(x)-2,nrow(x)-1,nrow(x)),1]
    x.richness <- x[-c(nrow(x)-3,nrow(x)-2,nrow(x)-1,nrow(x)),2]
  }else if (type == 'family'){
    x.counts <- x[-c(nrow(x)-4,nrow(x)-2,nrow(x)-1,nrow(x)),1]
    x.richness <- x[-c(nrow(x)-4,nrow(x)-2,nrow(x)-1,nrow(x)),6]  
  }else if (type == 'order'){
    x.counts <- x[-c(nrow(x)-4,nrow(x)-3,nrow(x)-1,nrow(x)),1]
    x.richness <- x[-c(nrow(x)-4,nrow(x)-3,nrow(x)-1,nrow(x)),10]  
  }else if (type == 'class'){
    x.counts <- x[-c(nrow(x)-4,nrow(x)-3,nrow(x)-2,nrow(x)),1]
    x.richness <- x[-c(nrow(x)-4,nrow(x)-3,nrow(x)-2,nrow(x)),14]  
  }else if (type == 'phylum'){
    x.counts <- x[-c(nrow(x)-4,nrow(x)-3,nrow(x)-2,nrow(x)-1),1]
    x.richness <- x[-c(nrow(x)-4,nrow(x)-3,nrow(x)-2,nrow(x)-1),18]
  }
  
  if (length(x.counts) > 3 & length(unique(x.richness)) != 1){
    
    # MM model
    datas <- data.frame(x.counts, x.richness)
    datas <- datas[!is.na(datas[,2]),]
    
    MMcurve <- formula(x.richness ~ Vmax*x.counts/(Km+x.counts))
    bestfit <- nls(formula = MMcurve, 
                   data = datas, 
                   start = list(Vmax = max(datas[,2]),
                                Km = max(datas[,1])/2))
    SconcRange <- seq(f = 0, t = xlims[2], by = 10)
    theorLine <- predict(bestfit, list(x.counts = SconcRange))
    
    # Output
    out <- list(summary(bestfit)$parameters[,1], xy = cbind(SconcRange, theorLine))
    return(out)
  }
  
}


# 3. Iterative sampling
#######################

## Assumption 1: all counts are "true" alignments
## Assumption 2: all multi-alignments have been well taken care of by classifiers

# a. Parameters

## Number of reads
read_steps <- seq(f = 100, t = 1800000, by = 100)

## Bootstraps per sampling depth
n_rep <- 1

## Minimum number of supporting reads per genus/family/order/class/phylum for richness/shannon/simpson indices
min_reads <- 5

# c. Iterative sampling

## Prepare summary matrix
sampling <- vector(mode = 'list', length = 12)
names(sampling) <- c('1', '2', '3', '4', '5', '6', '7', '8', '9.1', '9.2', 'N', 'P')
for (i in 1:length(sampling)){
  sampling[[i]] <- vector(mode = 'list', length = 3)
  names(sampling[[i]]) <- c("April", "June", "August")
  sampling[[i]]$April <- data.frame(matrix(NA, ncol = 21, nrow = length(read_steps)*n_rep))
  colnames(sampling[[i]]$April) <- c('Counts used',
                                     'Genus richness', 'Genus Shannon-H', 'Genus Shannon-E', 'Genus Simpson-D', 
                                     'Family richness', 'Family Shannon-H', 'Family Shannon-E', 'Family Simpson-D',
                                     'Order richness', 'Order Shannon-H', 'Order Shannon-E', 'Order Simpson-D',
                                     'Class richness', 'Class Shannon-H', 'Class Shannon-E', 'Class Simpson-D',
                                     'Phylum richness', 'Phylum Shannon-H', 'Phylum Shannon-E', 'Phylum Simpson-D')
  sampling[[i]]$August <- sampling[[i]]$June <- sampling[[i]]$April
}

## Iterate over each Cam site
for (barcode in 1:length(sampling)){
  
  # Iterate over each time point (April, June, August) per Cam site
  for (time in 1:length(sampling[[barcode]])){
    
    # Where we at, ma'es ?
    print(paste0('Processing sample ', 
                 names(sampling)[barcode], ' - ',
                 names(classifications)[time]))
    
    sampling[[barcode]][[time]][, 'Counts used'] <- rep(sort(rep(read_steps, n_rep)))
    
    ## fetch data
    tmp.data <- classifications[[time]][,c(1,2,5+barcode)]
    
    ## sample multinomially on genus, order and phylum level, after removing NAs
    tmp.data.genus <- tmp.data[tmp.data[,'taxRank'] == 'G',]
    tmp.data.family <- tmp.data[tmp.data[,'taxRank'] == 'F',]
    tmp.data.order <- tmp.data[tmp.data[,'taxRank'] == 'O',]
    tmp.data.class <- tmp.data[tmp.data[,'taxRank'] == 'C',]
    tmp.data.phylum <- tmp.data[tmp.data[,'taxRank'] == 'P',]
    
    tmp.data.genus <- tmp.data.genus[which(tmp.data.genus[,3] >= min_reads),]
    tmp.data.family <- tmp.data.family[which(tmp.data.family[,3] >= min_reads),]
    tmp.data.order <- tmp.data.order[which(tmp.data.order[,3] >= min_reads),]
    tmp.data.class <- tmp.data.class[which(tmp.data.class[,3] >= min_reads),]
    tmp.data.phylum <- tmp.data.phylum[which(tmp.data.phylum[,3] >= min_reads),]
    
    ## also generate richness & shannon indeces for each sample
    sampling[[barcode]][[time]][,c(2:5)] <- sample_classifications(tmp.data.genus, n_rep, read_steps)$stats
    sampling[[barcode]][[time]][,c(6:9)] <- sample_classifications(tmp.data.family, n_rep, read_steps)$stats
    sampling[[barcode]][[time]][,c(10:13)] <- sample_classifications(tmp.data.order, n_rep, read_steps)$stats
    sampling[[barcode]][[time]][,c(14:17)] <- sample_classifications(tmp.data.class, n_rep, read_steps)$stats
    sampling[[barcode]][[time]][,c(18:21)] <- sample_classifications(tmp.data.phylum, n_rep, read_steps)$stats
    
    ## remove samplings which would lie higher than the actual number of classifications
    ## therefore: tresholds (slightly) different between genus, family, order, class and phylum, as alignment confidence increases towards the latter
    sampling_thresh <- max(sum(tmp.data.genus[,3], na.rm = T),
                           sum(tmp.data.family[,3], na.rm = T),
                           sum(tmp.data.order[,3], na.rm = T),
                           sum(tmp.data.class[,3], na.rm = T),
                           sum(tmp.data.phylum[,3], na.rm = T))
    sampling[[barcode]][[time]] <- sampling[[barcode]][[time]][!sampling[[barcode]][[time]][,1] > sampling_thresh,,drop = F]
    
    ## set samplings which would lie higher above taxonomy-based thresholds to NA
    if(nrow(sampling[[barcode]][[time]]) > 0){
      sampling[[barcode]][[time]][sampling[[barcode]][[time]][,1] > sum(tmp.data.genus[,3], na.rm = T),2:5] <- NA
      sampling[[barcode]][[time]][sampling[[barcode]][[time]][,1] > sum(tmp.data.family[,3], na.rm = T),6:9] <- NA
      sampling[[barcode]][[time]][sampling[[barcode]][[time]][,1] > sum(tmp.data.order[,3], na.rm = T),10:13] <- NA
      sampling[[barcode]][[time]][sampling[[barcode]][[time]][,1] > sum(tmp.data.class[,3], na.rm = T),14:17] <- NA
      sampling[[barcode]][[time]][sampling[[barcode]][[time]][,1] > sum(tmp.data.phylum[,3], na.rm = T),18:21] <- NA
    }

    ## at the end of the summary matrix, add original classification stats for genus, order and phylum
    sampling[[barcode]][[time]] <- rbind(sampling[[barcode]][[time]],
                                         rep(NA, ncol(sampling[[barcode]][[time]])),
                                         rep(NA, ncol(sampling[[barcode]][[time]])),
                                         rep(NA, ncol(sampling[[barcode]][[time]])),
                                         rep(NA, ncol(sampling[[barcode]][[time]])),
                                         rep(NA, ncol(sampling[[barcode]][[time]])))
    colnames(sampling[[barcode]][[time]]) <- c('Counts used',
                                               'Genus richness', 'Genus Shannon-H', 'Genus Shannon-E', 'Genus Simpson-D',
                                               'Family richness', 'Family Shannon-H', 'Family Shannon-E', 'Family Simpson-D',
                                               'Order richness', 'Order Shannon-H', 'Order Shannon-E', 'Order Simpson-D',
                                               'Class richness', 'Class Shannon-H', 'Class Shannon-E', 'Class Simpson-D',
                                               'Phylum richness', 'Phylum Shannon-H', 'Phylum Shannon-E', 'Phylum Simpson-D')
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-4,1] <- sum(tmp.data.genus[,3], na.rm = T)
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-4,2:4] <- as.numeric(get_richness_shannon(tmp.data.genus[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-4,5] <- as.numeric(get_simpson(tmp.data.genus[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-3,1] <- sum(tmp.data.family[,3], na.rm = T)
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-3,6:8] <- as.numeric(get_richness_shannon(tmp.data.family[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-3,9] <- as.numeric(get_simpson(tmp.data.family[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-2,1] <- sum(tmp.data.order[,3], na.rm = T)
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-2,10:12] <- as.numeric(get_richness_shannon(tmp.data.order[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-2,13] <- as.numeric(get_simpson(tmp.data.order[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-1,1] <- sum(tmp.data.class[,3], na.rm = T)
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-1,14:16] <- as.numeric(get_richness_shannon(tmp.data.class[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]])-1,17] <- as.numeric(get_simpson(tmp.data.class[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]]),1] <- sum(tmp.data.phylum[,3], na.rm = T)
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]]),18:20] <- as.numeric(get_richness_shannon(tmp.data.phylum[,c(1,3)]))
    sampling[[barcode]][[time]][nrow(sampling[[barcode]][[time]]),21] <- as.numeric(get_simpson(tmp.data.phylum[,c(1,3)]))
  }
}

## save
save(file = '../results/2019-04-27_final_sampling/2019-05-01_final_sampling_kraken16S.Rdata',
     sampling)


# 4. Summarise information
##########################

load('../results/2019-04-27_final_sampling/2019-05-01_final_sampling_kraken16S.Rdata')

# a. Create summary tables for each taxon-level
summary.genus <- matrix(NA, 3*12, 4)
colnames(summary.genus) <- c('Sampling depth', 'Sampling richness', 
                             'Richness at depth 40000X', 'Richness at depth 40000X (% orig.)')
rownames(summary.genus) <- sort(c(paste(names(sampling), '- 1April'), 
                                  paste(names(sampling), '- 2June'),
                                  paste(names(sampling), '- 3August')))
rownames(summary.genus) <- gsub('1A|3A', 'A',rownames(summary.genus))
rownames(summary.genus) <- gsub('2J', 'J',rownames(summary.genus))
summary.phylum <- summary.class <- summary.order <- summary.family <- summary.genus
count <- 0
for (i in 1:12){
  for (j in 1:3){
    count <- count + 1
    summary.genus[count,1:2] <- as.numeric(sampling[[i]][[j]][nrow(sampling[[i]][[j]])-4,1:2])
    summary.family[count,1:2] <- as.numeric(sampling[[i]][[j]][nrow(sampling[[i]][[j]])-3,c(1,6)])
    summary.order[count,1:2] <- as.numeric(sampling[[i]][[j]][nrow(sampling[[i]][[j]])-2,c(1,10)])
    summary.class[count,1:2] <- as.numeric(sampling[[i]][[j]][nrow(sampling[[i]][[j]])-1,c(1,14)])
    summary.phylum[count,1:2] <- as.numeric(sampling[[i]][[j]][nrow(sampling[[i]][[j]]),c(1,18)])
    
    if(summary.genus[count,1] > 40000){
      out <- mmfit(sampling[[i]][[j]], xlims = c(0,80000), type = 'genus')
      summary.genus[count,3] <- round(out$xy[which(out$xy[,1] == 40000),2],1)
      summary.genus[count,4] <- round(100*summary.genus[count,3]/summary.genus[count,2],1)
    }
    if(summary.family[count,1] > 40000){
      out <- mmfit(sampling[[i]][[j]], xlims = c(0,80000), type = 'family')
      summary.family[count,3] <- round(out$xy[which(out$xy[,1] == 40000),2],1)
      summary.family[count,4] <- round(100*summary.family[count,3]/summary.family[count,2],1)
    }
    if(summary.order[count,1] > 40000){
      out <- mmfit(sampling[[i]][[j]], xlims = c(0,80000), type = 'order')
      summary.order[count,3] <- round(out$xy[which(out$xy[,1] == 40000),2],1)
      summary.order[count,4] <- round(100*summary.order[count,3]/summary.order[count,2],1)
    }
    if(summary.class[count,1] > 40000){
      out <- mmfit(sampling[[i]][[j]], xlims = c(0,80000), type = 'class')
      summary.class[count,3] <- round(out$xy[which(out$xy[,1] == 40000),2],1)
      summary.class[count,4] <- round(100*summary.class[count,3]/summary.class[count,2],1)
    }
    if(summary.phylum[count,1] > 40000){
      out <- mmfit(sampling[[i]][[j]], xlims = c(0,80000), type = 'phylum')
      summary.phylum[count,3] <- round(out$xy[which(out$xy[,1] == 40000),2],1)
      summary.phylum[count,4] <- round(100*summary.phylum[count,3]/summary.phylum[count,2],1)
    }
    
  }
}

# b. Create summary stats and sets at fixed (40k) downsampling depth for each taxon level
sample_classifications_single <- function(x, n_rep, read_steps, type, min_reads){
  
  # a. pre-process
  if(type == 'genus'){
    x <- x[x[,'taxRank'] == 'G',]
  }else if(type == 'family'){
    x <- x[x[,'taxRank'] == 'F',]
  }else if(type == 'order'){
    x <- x[x[,'taxRank'] == 'O',]
  }else if(type == 'class'){
    x <- x[x[,'taxRank'] == 'C',]
  }else if(type == 'phylum'){
    x <- x[x[,'taxRank'] == 'P',]
  }
  x <- x[which(x[,3] >= min_reads),]
  
  # b. run sampling function
  if(nrow(x) > 0){
    out <- sample_classifications(x, n_rep, read_steps)
    out.stats <- cbind(read_steps, out$stats)
    colnames(out.stats) <- c('Depth', "Richness", "Shannon H", "Shannon E", 'Simpson D')
    return(list(stats = out.stats, sets = out$sets)) 
  }else{
    return(NA)
  }
  
}
phylum.40000 <- class.40000 <- order.40000 <- family.40000 <- genus.40000 <- vector(mode = 'list', length = 12)
names(phylum.40000) <- names(class.40000) <- names(order.40000) <- names(family.40000) <- names(genus.40000) <- c('1', '2', '3', '4', '5', '6', '7', '8', '9.1', '9.2', 'N', 'P')
for (i in 1:length(sampling)){
  print(i)
  genus.40000[[i]] <- vector(mode = 'list', length = 3)
  names(genus.40000[[i]]) <- c("April", "June", "August")
  phylum.40000[[i]] <- class.40000[[i]] <- order.40000[[i]] <- family.40000[[i]] <- genus.40000[[i]]
  for (j in 1:3){
    genus.40000[[i]][[j]] <- sample_classifications_single(classifications[[j]][,c(1,2,5+i)], 100, 40000, 'genus', 5)
    family.40000[[i]][[j]] <- sample_classifications_single(classifications[[j]][,c(1,2,5+i)], 100, 40000, 'family', 5)
    order.40000[[i]][[j]] <- sample_classifications_single(classifications[[j]][,c(1,2,5+i)], 100, 40000, 'order', 5)
    class.40000[[i]][[j]] <- sample_classifications_single(classifications[[j]][,c(1,2,5+i)], 100, 40000, 'class', 5)
    phylum.40000[[i]][[j]] <- sample_classifications_single(classifications[[j]][,c(1,2,5+i)], 100, 40000, 'phylum', 5)
  }
}

## save
save(file = '../results/2019-04-27_final_sampling/2019-05-01_final_sampling_kraken16S_40000X_sets_stats.Rdata',
     phylum.40000, class.40000, order.40000, family.40000, genus.40000,
     summary.phylum, summary.class, summary.order, summary.family, summary.genus)
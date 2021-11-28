library(ape)
library(phytools)
library(paleotree)
library(geiger)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(tidyr)
library(ggplot2)
library(xtable)

day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")

# set whether to run analyses on surface-level or in-habitat niches
ss <- TRUE

# Format phylo ------------------------------------------------------------
# this code chunk is copied from analysis_phylo_eco.R, ForamNiches repo

# data(macroperforateForam)
# 
# # 2 species are assigned to a different genus in the phylogeny
# foramAMb$Species_name <- gsub('Globigerinoides_sacculifer', 'Trilobatus_sacculifer', 
#                               foramAMb$Species_name)
# foramAMb$Species_name <- gsub('Globigerinoides_trilobus', 'Trilobatus_trilobus', 
#                               foramAMb$Species_name)
# 
# # from paleotree example: convert Aze et al's suppmat to paleotree-readable format
# createTaxaData <- function(table){
#   #reorder table by first appearance
#   timetable <- table[order(-as.numeric(table[,3])),]
#   ID <- 1:nrow(table)
#   anc <- sapply(table[,2], function(x){
#     if(!is.na(x)){
#       which(x == table[,1])
#     } else { NA }
#   })
#   stillAlive <- as.numeric(table[,4] == 0)
#   ages <- cbind(as.numeric(table[,3]),
#                 as.numeric(table[,4]))
#   res <- cbind(ID,anc,ages,stillAlive,ID)
#   colnames(res) <- c('taxon.id','ancestor.id','orig.time',
#                      'ext.time','still.alive','looks.like')
#   rownames(res) <- table[,1]
#   return(res)
# }
# 
# # warning: slow!
# taxaAMb <- createTaxaData(foramAMb)
# treeAMb <- taxa2phylo(taxaAMb)
# 
# # drop anagenetic zero-length-terminal-edge ancestors
# phyFull <- dropZLB(treeAMb)
# 
# # rename tips from alphanumeric codes to species names
# sppCodes <- phyFull$tip.label
# rowOrdr <- match(sppCodes, foramAMb$Species_code)
# phyFull$tip.label <- foramAMb$Species_name[rowOrdr]
# 
# saveRDS(phyFull, 'Data/Aze-tree-phylo-object.rds')

phyFull <- readRDS('Data/Aze-tree-phylo-object.rds')

# export figure for phylogeny of study species
  # png('Figs/phylo_black.png', 
  #     width = 7.5, height = 7.5, units = 'in', res = 300)
  # plot.phylo(trTrim, main = '')
  # dev.off()

# Format trait data -------------------------------------------------------

if (ss){
  df <- read.csv('Data/niche-sumry-metrics_SJ-ste_SS_2020-11-15.csv')
} else {
  df <- read.csv('Data/niche-sumry-metrics_SJ-ste_hab_2020-11-15.csv')
}

df$sp <- gsub(' ', '_', df$sp)

# restrict to last 2 glacial intervals
maxAge <- 156
inSpan <- df$bin <= maxAge
df <- df[inSpan,]
bins <- sort(unique(df$bin))

# remove 2 species not in phylogeny
toss <- name.check(phyFull, df$sp, data.names = df$sp)
noPhy <- unique(toss$data_not_tree)
paste(paste(noPhy, collapse = ' '), 'not in tree')
# [1] "Globigerinita_glutinata Neogloboquadrina_incompta not in tree"
if (length(noPhy) > 0){
  rows2toss <- df$sp %in% noPhy
  df <- df[!rows2toss,]
}
spp <- unique(df$sp)

# drop tips not sampled in niche data
phyTrim <- keep.tip(phyFull, spp)

# extract mean occupied temperatures from a given time bin (e.g. most recent)
# give data for all species that are available
# print details of what isn't available, to be dealt with separately later
getTips <- function(df, bin, spp, binNm, trtNm, spNm){
  binBool <- df[, binNm] == bin
  dfTips <- df[binBool,]
  nSamp <- nrow(dfTips)
  if (nSamp < length(spp)){
    notSamp <- setdiff(spp, dfTips[, spNm])
    print(paste(bin, notSamp))
  }
  trait <- dfTips[, trtNm]
  names(trait) <- dfTips[, spNm]
  return(trait)
}
# test <- getTips(bin = 52, df = df, spp = spp,
#   binNm = 'bin', trtNm = 'm', spNm = 'sp')

# vector of sp traits per bin, synchronous across spp and in chronological order
tipLall <- sapply(bins, FUN = getTips, df = df, spp = spp,
                  binNm = 'bin', trtNm = 'm', spNm = 'sp')
names(tipLall) <- paste(bins)
# don't use bins missing spp for trait models; need same spp for fair comparisons
# [1] "52 Globigerinoides_conglobatus"
# [1] "156 Beella_digitata"
tipLchron <- Filter( function(x){ 
  length(x) == length(spp) 
  }, 
  tipLall)

# Randomly sample tips ----------------------------------------------------

# for each species, randomly pick a time bin at which to measure its trait
sampleTips <- function(pool, sppVect, binL){
  tipSamp <- sample(pool, length(sppVect), replace = TRUE)
  names(tipSamp) <- sppVect
  getVal <- function(s){
    spBin <- tipSamp[s]
    binL[[spBin]][s]
  }
  valVect <- sapply(sppVect, getVal)
  names(valVect) <- sppVect
  valVect
}

# set how many times to iterate the random sampling of tips:
n <- 10

# names (character) of time bins with available data to serve as tips
# use whole time interval except for B digitata and G conglobatus
tipPoolAll <- names(tipLall)
tipPoolBd <- tipPoolAll[tipPoolAll != 156]
tipPoolGc <- tipPoolAll[tipPoolAll != 52]

mostSpp <- spp[! spp %in% c('Beella_digitata', 'Globigerinoides_conglobatus')]

# customize sampling pools for species unsampled in some time bins, 
# but combine into one sample for analysis
sampleAvlbTips <- function(restSpp, poolBd, poolGc, poolAll, binL){
  sampBd <-   sampleTips(poolBd, 'Beella_digitata', binL)
  sampGc <-   sampleTips(poolGc, 'Globigerinoides_conglobatus', binL)
  sampRest <- sampleTips(poolAll, restSpp, binL)
  c(sampRest, sampBd, sampGc)
}

tipLrdm <- replicate(n = n,
                     sampleAvlbTips(restSpp = mostSpp,
                                    poolBd  = tipPoolBd,
                                    poolGc  = tipPoolGc,
                                    poolAll = tipPoolAll,
                                    binL = tipLall),
                     simplify = FALSE
                     )

# Estimate params ---------------------------------------------------------

fitEvoMods <- function(phy, trt){
  brownFit  <- fitContinuous(phy, trt, model = 'BM')
  lambdaFit <- fitContinuous(phy, trt, model = 'lambda')
  deltaFit  <- fitContinuous(phy, trt, model = 'delta')
  kappaFit  <- fitContinuous(phy, trt, model = 'kappa')
  ouFit     <- fitContinuous(phy, trt, model = 'OU')
  ebFit     <- fitContinuous(phy, trt, model = 'EB')

  # compare all models as AIC weights
  modL <- list(BM = brownFit, 
               lambda = lambdaFit, 
               delta = deltaFit,
               kappa = kappaFit, 
               OU = ouFit, 
               EB = ebFit
  )
  getAic <- function(mod){ mod$opt$aic }
  aicVect <- unlist(lapply(modL, getAic))
  aicMin <- aicVect[which.min(aicVect)]
  bestMod <- names(aicMin)
  aicDelta <- aicVect - aicMin
  relLik <- exp(-0.5 * aicDelta)
  wts <- relLik/sum(relLik)
  
  data.frame(bestMod, t(wts))
}
# fitEvoMods(phy = phyTrim, trt = tipLchron[[10]]) # test

modsLchron <- lapply(tipLchron, fitEvoMods, phy = phyTrim)
modsDfChron <- do.call(rbind, modsLchron)
# 1.4 min runtime
# delta and kappa models give 40 warnings of parameter estimates at bounds

modsLrdm <- lapply(tipLrdm, fitEvoMods, phy = phyTrim)
# delta, kappa, and EB models warn parameter estimates at bounds (x22)
modsDfRdm <- do.call(rbind, modsLrdm)

# Phylogenetic Comparative Methods p. 91:
# There are two main ways to assess the fit of the three Pagel-style models to data.
# First, one can use ML to estimate parameters and likelihood ratio tests 
# (or AICc scores) to compare the fit of various models. Each represents a three parameter
# model: one additional parameter added to the two parameters already needed to
# describe single-rate Brownian motion. As mentioned above, simulation studies
# suggest that this can sometimes lead to overconfidence, at least for the lambda model.
# Sometimes researchers will compare the fit of a particular model (e.g. lambda) with
# models where that parameter is fixed at its two extreme values (0 or 1; this is not
# possible with delta). Second, one can use Bayesian methods to estimate posterior
# distributions of parameter values, then inspect those distributions to see if they
# overlap with values of interest (say, 0 or 1).

# fastAnc(phyFull, trait)

# Charts for chronological sampling ---------------------------------------

# stacked bar chart of support for each evo model, for each tip time step

# plot youngest time step at top (strat order)
# this means putting df in reverse order, to be flipped during xy rotation later
mods <- colnames(modsDfChron)[-1]
binsInPlot <- row.names(modsDfChron)
modsDfChron$bin <- factor(binsInPlot, levels = rev(binsInPlot)) 
# (bins shouldn't be numeric since plotted as discrete axis)

modsLong <- pivot_longer(modsDfChron, cols = all_of(mods),
                         names_to = 'Model', values_to = 'Weight')
modsLong$Model <- factor(modsLong$Model, levels = rev(mods))

# colr <- c('BM' = '#283593', 
#           'lambda' = '#5dade2', 
#           'delta' = '#FFC300', 
#           'kappa' = '#FF5733',
#           'OU' = '#000000',
#           'EB' = '#57de36')

barStack <- 
  ggplot() +
  theme_minimal() +
  # use position=fill to indicate data as %, so rounding errors don't result in sum > 1
  geom_bar(data = modsLong, aes(x = bin, y = Weight, fill = Model), 
           stat = 'identity',
           position = 'fill',
           width = 0.85
           ) +
  scale_x_discrete(name = 'Tip age (ka)',
                   labels = rev(binsInPlot)) +
  scale_y_continuous(name = 'Model support (AIC weight)        ') +
#                    expand = c(0, 0), limits = c(0, 1.2),
#                    breaks = seq(0, 1, by = 0.25)) + # auto breaks in right posion
  theme(legend.position = 'top',
        legend.text = element_text(size = 8),
        axis.text.y = element_text(margin = margin(r = -10)),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_line(size = 0.5, colour = 'grey30'),
        axis.ticks.length = unit(4, "pt"),
        panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE)
    # fill = guide_legend(nrow = 2) # legend already 2 rows
         ) +
  coord_flip()
#   scale_colour_manual(name = element_blank(), values = colr, 
#                       aesthetics = 'fill', limits = evoModes,
#                       labels = c('Strict stasis','Stasis',
#                                  'Random walk','Directional walk',
#                                  'Tracking')) +
#   

if (ss){
  barNm <- paste0('Figs/phylo-evo-model-support-barplot_SS_',  day, '.pdf')
} else {
  barNm <- paste0('Figs/phylo-evo-model-support-barplot_hab_', day, '.pdf')
}
pdf(barNm, width = 4, height = 6)
  barStack
dev.off()

# export the model support for each time bin as a suppl table:
# some values are too similar to compare easily form the barplot alone

# last column is the time bin identifier; use rownames instead to put it 1st
wdthDf <- ncol(modsDfChron)

chronTbl <- xtable(modsDfChron[,-wdthDf], align = rep('r', wdthDf), digits = 3,
               caption = 'Trait evolution model weights by time bin')
if (ss){
  tblNm <- paste0('Figs/phylo-evo-model-weights_chron_SS_',  day, '.tex')
} else {
  tblNm <- paste0('Figs/phylo-evo-model-weights_chron_hab_', day, '.tex')
}
if (ss){
  print(chronTbl, file = tblNm, include.rownames = TRUE,
        caption.placement = 'top')
} else {
  print(chronTbl, file = tblNm, include.rownames = TRUE,
        caption.placement = 'top')
}

# Example code ------------------------------------------------------------

# copied and edited where necessary from function documentation
# https://www.rdocumentation.org/packages/geiger/versions/1.2-06/topics/fitContinuous

data(geospiza)
attach(geospiza)

#---- PRINT RESULTS
fitContinuous(geospiza.tree, geospiza.data)

#---- STORE RESULTS 
brownFit <-  fitContinuous(geospiza.tree, geospiza.data)
aic.brown<-numeric(5)
# for(i in 1:5) aic.brown[i]<-brownFit[[i]]$aic # original line fails
for(i in 1:5) aic.brown[i]<-brownFit[[i]]$opt$aic

lambdaFit<-fitContinuous(geospiza.tree, geospiza.data, model="lambda")
aic.lambda<-numeric(5)
for(i in 1:5) aic.lambda[i]<-lambdaFit[[i]]$opt$aic # same edit as for BM

deltaFit<-fitContinuous(geospiza.tree, geospiza.data, model="delta")
aic.delta<-numeric(5)
for(i in 1:5) aic.delta[i]<-deltaFit[[i]]$opt$aic

kappaFit<-fitContinuous(geospiza.tree, geospiza.data, model="kappa")
aic.kappa<-numeric(5)
for(i in 1:5) aic.kappa[i]<-kappaFit[[i]]$opt$aic

ouFit<-fitContinuous(geospiza.tree, geospiza.data, model="OU")
aic.ou<-numeric(5)
for(i in 1:5) aic.ou[i]<-ouFit[[i]]$opt$aic

ebFit<-fitContinuous(geospiza.tree, geospiza.data, model="EB")
aic.eb<-numeric(5)
for(i in 1:5) aic.eb[i]<-ebFit[[i]]$opt$aic

#   COMPARE ALL MODELS

# One way: use likelihood ratio test to compare all models to Brownian model

# Another way: use AIC

aic.all<-cbind(aic.brown, aic.lambda, aic.delta, aic.kappa, aic.ou, aic.eb)
foo<-function(x) x-x[which(x==min(x))]
daic<-t(apply(aic.all, 1, foo))

rownames(daic)<-colnames(geospiza.data)
colnames(daic)<-c("Brownian", "Lambda", "Delta", "Kappa", "OU", "EB")

cat("Table of delta-aic values; zero - best model")
print(daic, digits=2)


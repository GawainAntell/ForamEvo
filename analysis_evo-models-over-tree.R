library(ape)
library(phytools)
library(paleotree)
library(geiger)
library(parallel)
library(foreach)
# library(iterators)
library(doParallel)
library(tidyr)
library(ggplot2)

day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")

# TODO
# set option to use SS or within-habitat data

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

df <- read.csv('Data/niche-sumry-metrics_SJ-ste_SS_2020-11-15.csv')
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
getTips <- function(df, bin, nSpp, binNm, trtNm, spNm){
  binBool <- df[, binNm] == bin
  dfTips <- df[binBool,]
  nSamp <- nrow(dfTips)
  if (nSamp < nSpp){
    return(NULL)
    # throw error instead?
  } else {
    trait <- dfTips[, trtNm]
    names(trait) <- dfTips[, spNm]
    return(trait)
  }
}
# test <- getTips(bin = 4, df = df, nSpp = length(spp),
#   binNm = 'bin', trtNm = 'm', spNm = 'sp')

# make list object, vector of sp traits per bin
tipLall <- sapply(bins, FUN = getTips,
               df = df, nSpp = length(spp),
               binNm = 'bin', trtNm = 'm', spNm = 'sp')
names(tipLall) <- paste(bins)
# bin 52 is missing Globigerinoides_conglobatus
# bin 156 is missing Beella_digitata
tipL <- Filter(Negate(is.null), tipLall)

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
# fitEvoMods(phy = phyTrim, trt = tipL[[10]]) # test

modsL <- lapply(tipL, fitEvoMods, phy = phyTrim)
# 1.4 min runtime
# delta and kappa models give many warnings of parameter estimates at bounds

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

# Relative model support --------------------------------------------------
# stacked bar chart of support for each evo model, for each tip time step

modsDf <- do.call(rbind, modsL)
mods <- colnames(modsDf)[-1]
modsDf$bin <- as.numeric(row.names(modsDf))
modsLong <- pivot_longer(modsDf, cols = all_of(mods),
                         names_to = 'Model', values_to = 'Weight')
modsLong$Model <- as.factor(modsLong$Model)
levels(modsLong$Model) <- mods

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
           position = 'fill'
          # width = 0.75
           ) +
  scale_x_discrete(name = 'Tip age (ka)',
                   labels = modsLong$bin) +
  scale_y_continuous(name = 'Model support (AIC weight)        ') +
#                    expand = c(0, 0), limits = c(0, 1.2),
#                    breaks = seq(0, 1, by = 0.25)) +
  theme(legend.position = 'top',
        legend.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_line(size = 0.5, colour = 'grey30'),
        axis.ticks.length = unit(4, "pt"),
        axis.text.y = element_text(),
        panel.grid = element_blank()
  )
#   scale_colour_manual(name = element_blank(), values = colr, 
#                       aesthetics = 'fill', limits = evoModes,
#                       labels = c('Strict stasis','Stasis',
#                                  'Random walk','Directional walk',
#                                  'Tracking')) +
#   guides(fill = guide_legend(nrow = 2))

# manually label the time bins and rotate to vertical axis
binLbl <- head(bins, -1)
barsFlip <- 
  barStack +
  geom_text(aes(x = rev(binLbl), y = -.01, label = binLbl),
            hjust = 1, size = 3.2) +
  coord_flip()

# if (ss){
  barNm <- paste0('Figs/phylo-evo-model-support-barplot_SS_', day, '.pdf')
# } else {
#   barNm <- paste0('Figs/evo-model-support-barplot_hab_',day,'.pdf')
# }
pdf(barNm, width = 4, height = 6)
  barsFlip
dev.off()

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


library(ape)
library(phytools)
library(paleotree)
library(geiger)
library(parallel)
library(foreach)
# library(iterators)
library(doParallel)
library(ggplot2)

# day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")

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
# phyloFull <- dropZLB(treeAMb)
# 
# # rename tips from alphanumeric codes to species names
# sppCodes <- phyloFull$tip.label
# rowOrdr <- match(sppCodes, foramAMb$Species_code)
# phyloFull$tip.label <- foramAMb$Species_name[rowOrdr]
# 
# saveRDS(phyloFull, 'Data/Aze-tree-phylo-object.rds')

phyloFull <- readRDS('Data/Aze-tree-phylo-object.rds')

# Format trait data -------------------------------------------------------

df <- read.csv('Data/niche-sumry-metrics_SJ-ste_SS_2020-11-15.csv')
df$sp <- gsub(' ', '_', df$sp)

# restrict to last 2 glacial intervals
maxAge <- 156
inSpan <- df$bin <= maxAge
df <- df[inSpan,]
bins <- sort(unique(df$bin))

# remove 2 species not in phylogeny
toss <- name.check(phyloFull, df$sp, data.names = df$sp)
noPhy <- unique(toss$data_not_tree)
paste(paste(noPhy, collapse = ' '), 'not in tree')
# [1] "Globigerinita_glutinata Neogloboquadrina_incompta not in tree"
if (length(noPhy) > 0){
  rows2toss <- df$sp %in% noPhy
  df <- df[!rows2toss,]
}
spp <- unique(df$sp)

# extract mean occupied temperatures from a given time bin (e.g. most recent)
getTips <- function(df, bin, tr, nSpp, binNm, trtNm, spNm){
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
# test <- getTips(bin = 4, df = df, tr = phyloFull, nSpp = length(spp),
#   binNm = 'bin', trtNm = 'm', spNm = 'sp')

# make list object, vector of sp traits per bin
tipLall <- sapply(bins, FUN = getTips,
               df = df, tr = phyloFull, nSpp = length(spp),
               binNm = 'bin', trtNm = 'm', spNm = 'sp')
names(tipLall) <- paste(bins)
# bin 52 is missing Globigerinoides_conglobatus
# bin 156 is missing Beella_digitata
tipL <- Filter(Negate(is.null), tipLall)

# Estimate params ---------------------------------------------------------

#TODO
# wrap this section into a function
# apply function over list of trait data by time bin

brownFit  <- fitContinuous(phyloFull, trait, model = 'BM')
lambdaFit <- fitContinuous(phyloFull, trait, model = 'lambda')
deltaFit  <- fitContinuous(phyloFull, trait, model = 'delta')
kappaFit  <- fitContinuous(phyloFull, trait, model = 'kappa')
ouFit     <- fitContinuous(phyloFull, trait, model = 'OU')
ebFit     <- fitContinuous(phyloFull, trait, model = 'EB')

# alt to estimate Pagel's lambda: 
# use nlme method - caper oddly drops 8 tree tips
  # pgl <- corPagel(1, form = ~ species, trTrim)
  # pglMod <- gls(h ~ eco, data = spDf, correlation = pgl)
  # summary(pglMod)

# compare all models
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
aicDelta <- aicVect - aicMin
sort(aicDelta)

# bin = 4
# > OU        lambda    delta     BM        kappa     EB 
# > 0.0000000 0.7616984 1.3060017 2.9660730 4.9409742 4.9661791 

# bin = 100
# > OU        delta     lambda    BM        kappa     EB 
# > 0.0000000 0.7308084 0.9828120 2.2013227 4.1862505 4.2014246

# Figures -----------------------------------------------------------------

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

# calculate weights and plot as bar chart
relLik <- exp(-0.5*aicDelta)
wts <- relLik/sum(relLik)
sort(wts, decreasing = T)

# figure - save wts for each bin as objects wts4 and wts100 (see note at top for edit)
dfWts <- data.frame(wts4, wts700)
dfWts$mod <- row.names(dfWts)

modOrdr <- names(sort(wts4,decreasing = T))
bar4 <- 
  ggplot(data = dfWts) +
  theme_minimal() +
  scale_x_discrete(limits = c(modOrdr)) +
  theme(axis.title = element_blank()) +
  geom_bar(aes(x = mod, y = wts4), stat = 'identity')

png('Figs/evo-mod-wts_4ka.png',
    width = 4, height = 4, units = 'in', res = 300)
bar4 + 
  scale_fill_manual(values=c('grey30','#FAA003'))
dev.off()

library(tidyr)
modLong <- pivot_longer(dfWts, cols = starts_with('wts'), 
                        names_to = 'TipAge', values_to = 'weight')

barStack <- 
  ggplot(modLong) +
  theme_minimal() +
  scale_x_discrete(limits = c(modOrdr)) +
  theme(axis.title = element_blank()) +
  geom_bar(aes(x = mod, y = weight, fill = TipAge), 
           stat = 'identity',
           position = position_dodge())
png('Figs/evo-mod-wts_4-vs-700ka.png',
    width = 4, height = 4, units = 'in', res = 300)
barStack +
  scale_fill_manual(values=c('grey30','#FAA003')) +
  theme(legend.position = 'null')
dev.off()

# fastAnc(phyloFull, trait)

# 4 ka parameters
  # > ouFit$opt$z0
  # [1] 21.81366
  # > ouFit$opt$alpha
  # [1] 0.07985

# 28 ka parameters
  # > ouFit$opt$z0
  # [1] 17.86773
  # > ouFit$opt$alpha
  # [1] 0.12269

# 100 ka parameters
  # > ouFit$opt$z0
  # [1] 18.7577
  # > ouFit$opt$alpha
  # [1] 0.05753173

# 700 ka parameters
  # > ouFit$opt$z0
  # [1] 20.72465
  # $alpha
  # [1] 0.1899055

png('Figs/phylo_black.png', 
    width = 7.5, height = 7.5, units = 'in', res = 300)
plot.phylo(trTrim, main = '')
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


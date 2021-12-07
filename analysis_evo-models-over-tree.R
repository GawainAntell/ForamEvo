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
library(RColorBrewer)
library(xtable)

day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")

# set whether to run analyses on surface-level or in-habitat niches
ss <- TRUE

# set how many times to iterate the random sampling of tips
nRdm <- 1000

# whether to exclude super-densely-sampled 4 ka and 12 ka bins (TRUE) or not
# this only applies to random tip sampling, which could mix good/bad sampled bins
rdmAsOld <- TRUE

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

# calculate standard error to input in trait models
df$se <- df$sd / sqrt(df$n)

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

# extract mean/SE occupied temperatures from a given time bin (e.g. most recent)
# give data for all species that are available
# print details of what isn't available, to be dealt with separately later
getTips <- function(df, bin, spp, binNm, val1, val2, spNm){
  binBool <- df[, binNm] == bin
  dfTips <- df[binBool,]
  nSamp <- nrow(dfTips)
  if (nSamp < length(spp)){
    notSamp <- setdiff(spp, dfTips[, spNm])
    print(paste(bin, notSamp))
  }
  dat1 <- dfTips[, val1]
  dat2 <- dfTips[, val2]
  names(dat1) <- names(dat2) <- dfTips[, spNm]
  outL <- list(dat1, dat2)
  names(outL) <- c(val1, val2)
  return(outL)
}
# test <- getTips(bin = 52, df = df, spp = spp,
#   binNm = 'bin', val1 = 'm', val2 = 'se', spNm = 'sp')

# vector of sp m and SE per bin, synchronous across spp and in chronological order
tipLall <- lapply(bins, FUN = getTips, df = df, spp = spp,
                  binNm = 'bin', val1 = 'm', val2 = 'se', spNm = 'sp')
names(tipLall) <- paste(bins)
# don't use bins missing spp for trait models; need same spp for fair comparisons
# [1] "52 Globigerinoides_conglobatus"
# [1] "156 Beella_digitata"

tipLchron <- Filter( function(x){ 
  length(x[[1]]) == length(spp) 
  }, 
  tipLall)

# Randomly sample tips ----------------------------------------------------

# for each species, randomly pick a time bin at which to measure its trait
  # testing: pool <- tipPoolAll; sppVect <- mostSpp; binL <- tipLall
sampleTips <- function(pool, sppVect, binL){
  tipSamp <- sample(pool, length(sppVect), replace = TRUE)
  names(tipSamp) <- sppVect
  getVal <- function(s){
    spBin <- tipSamp[s]
    spM  <- binL[[spBin]][['m' ]][s]
    spSE <- binL[[spBin]][['se']][s]
    cbind(spM, spSE)
  }
  sapply(sppVect, getVal)
}

# names (character) of time bins with available data to serve as tips
# watch out for B digitata and G conglobatus - each unsampled in 1 bin
# follow argument at script head - whether to exclude young, well-sampled bins
tipPoolAll <- names(tipLall)
if (rdmAsOld){
  oldOnly <- ! tipPoolAll %in% c('4', '12')
  tipPoolAll <- tipPoolAll[oldOnly]
}
tipPoolBd <- tipPoolAll[tipPoolAll != 156]
tipPoolGc <- tipPoolAll[tipPoolAll != 52]

mostSpp <- spp[! spp %in% c('Beella_digitata', 'Globigerinoides_conglobatus')]

# customize sampling pools for species unsampled in some time bins, 
# but combine into one sample for analysis
# output with same list structure as synchronous, chronological data
sampleAvlbTips <- function(restSpp, poolBd, poolGc, poolAll, binL){
  sampBd <-   sampleTips(poolBd, 'Beella_digitata', binL)
  sampGc <-   sampleTips(poolGc, 'Globigerinoides_conglobatus', binL)
  sampRest <- sampleTips(poolAll, restSpp, binL)
  sampMat <- cbind(sampRest, sampBd, sampGc)
  list('m' = sampMat[1,], 'se' = sampMat[2,])
}

tipLrdm <- replicate(n = nRdm,
                     sampleAvlbTips(restSpp = mostSpp,
                                    poolBd  = tipPoolBd,
                                    poolGc  = tipPoolGc,
                                    poolAll = tipPoolAll,
                                    binL = tipLall),
                     simplify = FALSE
                     )

# Evo model function ------------------------------------------------------

# Output a list with multiple elements, for different results plots: 
# dataframe of relative weights for all models (for stacked barplots)
# parameter estimates and CI from lambda/delta/OU, to check model quality
mods <- c('white', 'BM', 'lambda', 'delta', 'kappa', 'OU', 'EB')
fitEvoMods <- function(phy, m, err, params = FALSE){
  whFit     <- fitContinuous(phy, m, err, model = 'white',
                             control = list(hessian = TRUE))
  brownFit  <- fitContinuous(phy, m, err, model = 'BM',
                             control = list(hessian = TRUE))
  lambdaFit <- fitContinuous(phy, m, err, model = 'lambda',
                             control = list(hessian = TRUE))
  deltaFit  <- fitContinuous(phy, m, err, model = 'delta',
                             control = list(hessian = TRUE))
  kappaFit  <- fitContinuous(phy, m, err, model = 'kappa',
                             control = list(hessian = TRUE))
  ouFit     <- fitContinuous(phy, m, err, model = 'OU',
                             control = list(hessian = TRUE))
  ebFit     <- fitContinuous(phy, m, err, model = 'EB',
                             control = list(hessian = TRUE))
  
  # compare all models as AICc weights
  # AICc instead of AIC to match internal decision of paleoTs::compareModels()
  modL <- list(white = whFit,
               BM = brownFit, 
               lambda = lambdaFit, 
               delta = deltaFit,
               kappa = kappaFit, 
               OU = ouFit, 
               EB = ebFit
  )
  getAicc <- function(mod){ mod$opt$aicc }
  aiccVect <- unlist(lapply(modL, getAicc))
  aiccMin <- aiccVect[which.min(aiccVect)]
#  bestMod <- names(aicMin)
  aiccDelta <- aiccVect - aiccMin
  relLik <- exp(-0.5 * aiccDelta)
  wts <- relLik / sum(relLik)
  smryDf <- data.frame(t(wts)) # data.frame(bestMod, t(wts))
  # leave out character string (best model name) so data stay numeric in export
  
  # return parameters when fcn applied to chron data, but not for random replicates
  if (params == TRUE){
    # paramBest <- modL[[bestMod]]$opt
    # numParamBest <- Filter(is.numeric, paramBest)
    whtn <- modL[['white' ]]$opt
    lmda <- modL[['lambda']]$opt
    dlta <- modL[['delta' ]]$opt
    ornu <- modL[['OU'    ]]$opt
    numWhtn <- Filter(is.numeric, whtn)
    numLmda <- Filter(is.numeric, lmda)
    numDlta <- Filter(is.numeric, dlta)
    numOrnu <- Filter(is.numeric, ornu)
    list(wts = smryDf,
         # bestMod = unlist(numParamBest)
         white  = unlist(numWhtn),
         lambda = unlist(numLmda),
         delta  = unlist(numDlta),
         ou     = unlist(numOrnu)
         )
  } else {
    smryDf
  }
}
# test <- fitEvoMods(phy = phyTrim, m = tipLchron[[10]][['m']], 
#                    err = tipLchron[[10]][['se']], params = TRUE)

# *Run models on chron data -----------------------------------------------

modsLchron <- lapply(tipLchron, function(x){
  fitEvoMods(phy = phyTrim, m = x[['m']], err = x[['se']], params = TRUE)
})
# only a min runtime

# Reformat into separate dataframes for separate weight and estimate figures
# There's surely a much tidier way to do this, soz reader

getDf <- function(x) t(x$wts)
modsDfM <- sapply(modsLchron, getDf)
modsDfChron <- data.frame(t(modsDfM))
colnames(modsDfChron) <- mods # c('bestMod', mods)
# modsDfChron[, mods] <- apply(modsDfChron[, mods], 2, as.numeric)
# TODO add column for best model?
binsInPlot <- row.names(modsDfChron)
modsDfChron$bin <- factor(binsInPlot, levels = rev(binsInPlot)) 
# bins shouldn't be numeric since plotted as discrete axis
# put in reverse order so will plot youngest to oldest from top to bottom

# format parameter CI values - note not every run is able is calculate them
confide <- function(mod, modsL){
  paramsL <- lapply(modsL, function(x) x[[mod]])
  # names(paramsL) <- binsInPlot
  
  # find the runs where CI is calculated
  hasConf <- lapply(paramsL, function(x){
    'CI1' %in% names(x)
  })
  params2lookit <- paramsL[unlist(hasConf)]
  paramsDf <- data.frame(do.call(rbind, params2lookit))
  paramsDf$bin <- row.names(paramsDf)
  return(paramsDf)
}
paramsWhtn <- confide(mod = 'white',  modsL = modsLchron)
paramsLmda <- confide(mod = 'lambda', modsL = modsLchron)
paramsDlta <- confide(mod = 'delta',  modsL = modsLchron) # no CIs at all D:
paramsOrnu <- confide(mod = 'ou',     modsL = modsLchron)

# *Run models on rndm-samp data -------------------------------------------

pt1 <- proc.time()
nCore <- detectCores() - 1
registerDoParallel(nCore)
modsDfRdm <- foreach(x = tipLrdm, .packages = 'geiger',
                     .combine = rbind, .inorder = FALSE) %dopar%
  fitEvoMods(phy = phyTrim, m = x[['m']], err = x[['se']])
stopImplicitCluster()
pt2 <- proc.time()
(pt2-pt1)/60 
# delta, kappa, and EB models warn parameter estimates at bounds

# ca 2 min runtime for 100 reps, 12-15 min/ 1000x, so save results to jump to plots
# if (ss){
#   rdmDfNm <- paste0('Figs/phylo-evo-model-wts_random-samp_', nRdm, 'x_SS_',  day, '.csv')
# } else {
#   rdmDfNm <- paste0('Figs/phylo-evo-model-wts_random-samp_', nRdm, 'x_hab_', day, '.csv')
# }
# write.csv(modsDfRdm, rdmDfNm, row.names = FALSE)

# Charts for chronological sampling ---------------------------------------

# stacked bar chart of support for each evo model, for each tip time step

modsLong <- pivot_longer(modsDfChron, cols = all_of(mods),
                         names_to = 'Model', values_to = 'Weight')
modsLong$Model <- factor(modsLong$Model, levels = rev(mods))

# colorblind friendly diverging palette
colr <- brewer.pal(n = length(mods), name = "Dark2")
names(colr) <- mods

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
  scale_y_continuous(name = 'Model support (AICc weight)        ') +
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
  scale_colour_manual(name = element_blank(), values = colr,
                      aesthetics = 'fill', limits = rev(mods)
                      ) +
  coord_flip()

if (ss){
  barNm <- paste0('Figs/phylo-evo-model-wts-barplot_SS_',  day, '.pdf')
} else {
  barNm <- paste0('Figs/phylo-evo-model-wts-barplot_hab_', day, '.pdf')
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
  tblNm <- paste0('Figs/phylo-evo-model-wts_chron_SS_',  day, '.tex')
} else {
  tblNm <- paste0('Figs/phylo-evo-model-wts_chron_hab_', day, '.tex')
}
if (ss){
  print(chronTbl, file = tblNm, include.rownames = TRUE,
        caption.placement = 'top')
} else {
  print(chronTbl, file = tblNm, include.rownames = TRUE,
        caption.placement = 'top')
}

# Charts for random sampling ----------------------------------------------

# summarize weight for each model across all sampling iterations

modsLongRdm <- pivot_longer(modsDfRdm, cols = all_of(mods),
                            names_to = 'Model', values_to = 'Weight')
modsLongRdm$Model <- factor(modsLongRdm$Model, levels = mods)

# use same colour scheme as in stacked barplot figure
# but no legend - would be redundant with x-axis labels
wtBox <- ggplot() +
  theme_minimal() +
  geom_boxplot(data = modsLongRdm,
               aes(x = Model, y = Weight),
               fill = colr) +
  scale_y_continuous('Model support (AICc Weight)') +
  theme(axis.title.x = element_blank()) +
  scale_colour_manual(name = element_blank(), values = colr,
                      aesthetics = 'fill', limits = (mods)
  )
# TODO standardise y-scale for comparison between all-mods and old-mods vsns

if (rdmAsOld){
  wtBoxNm <- paste0('Figs/phylo-evo-model-wts_random-samp_', nRdm, 'x_oldOnly')
} else {
  wtBoxNm <- paste0('Figs/phylo-evo-model-wts_random-samp_', nRdm, 'x_allBins')
}
if (ss){
  wtBoxNm <- paste0(wtBoxNm, '_SS_',   day, '.pdf')
} else {
  wtBoxNm <- paste0(wtBoxNm, '_hab_',  day, '.pdf')
}
pdf(wtBoxNm, width = 4, height = 3)
  wtBox
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


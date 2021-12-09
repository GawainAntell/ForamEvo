library(ape)
library(phytools)
library(paleotree)
library(geiger)
library(pmc)
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

# number of Monte Carlo bootstrap replicates - default is 500, VERY slow
nBoot <- 50

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

# Fit evo models ----------------------------------------------------------

# 95% confidence interval based on MC bootstrapping, 
# because Hessian-calculated CIs are junk and not even solvable for all models
confide <- function(phy, m, err, mod, nboot){
  # pmc fits 4 models:
  # 'AA' fits model A on data obtained by simulations from model A,
  # 'BA' fits model B on the data simulated from model A
  # 'AB' fits model A on simulations from B
  # 'BB' fits B on simulations from B
  # So, by comparing a model with itself, both AA and BB replicates estimate the
  # variable of interest and can be used to draw a CI. Halve computation time.
  mcDat <- pmc(phy, m, mod, mod, 
               nboot = nboot/2, mc.cores = 1,
               optionsA = list(SE = err), # bounds = list(delta = c(min = exp(-500), max = 5))
               optionsB = list(SE = err)
  )
  # there seems to be a bug in pmc with specifying the bounds argument 
  # (fed to fitContinuous via optionsA or optionsB)
  # e.g. when delta upper bound set to 5, observed estimate goes up to 5 
  # but the upper 95% CI is still at the original default value, 3
  
  if (mod == 'OU'){
    paramNm <- 'alpha'
  } else {
    paramNm <- mod
  }
  varRows <- mcDat$par_dists$comparison %in% c('AA', 'BB') & 
    mcDat$par_dists$parameter == paramNm
  bootEsts <- mcDat$par_dists$value[varRows]
  ci <- quantile(bootEsts, c(0.025, 0.975))
  varObs <- coef(mcDat[['A']])[[paramNm]]
  c(observed = varObs, ci)
}

# Output a list with multiple elements, for different results plots: 
# dataframe of relative weights for all models (for stacked barplots)
# parameter estimates and CI from lambda/delta/OU, to check model quality
mods <- c('white', 'BM', 'lambda', 'delta', 'kappa', 'OU', 'EB')
fitEvoMods <- function(phy, m, err, params = FALSE, nBoot = 500){
  whFit     <- fitContinuous(phy, m, err, model = 'white')
  brownFit  <- fitContinuous(phy, m, err, model = 'BM')
  lambdaFit <- fitContinuous(phy, m, err, model = 'lambda')
  deltaFit  <- fitContinuous(phy, m, err, model = 'delta')
  kappaFit  <- fitContinuous(phy, m, err, model = 'kappa')
  ouFit     <- fitContinuous(phy, m, err, model = 'OU')
  ebFit     <- fitContinuous(phy, m, err, model = 'EB')
  
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
  aiccDelta <- aiccVect - aiccMin
  relLik <- exp(-0.5 * aiccDelta)
  wts <- relLik / sum(relLik)
  smryDf <- data.frame(t(wts))
  
  # return parameter estimates when fcn is applied to chronological data, 
  # but not for random-sampling replicates because that would take 5ever to run
  if (params == TRUE){
    lmdaEsts <- confide(phy, m = m, err = err, 'lambda', nBoot/2)
    dltaEsts <- confide(phy, m = m, err = err, 'delta',  nBoot/2)
    ornuEsts <- confide(phy, m = m, err = err, 'OU',     nBoot/2)
    estsMat <- rbind('lambda' = lmdaEsts, 
                     'delta' =  dltaEsts, 
                     'alpha' =  ornuEsts)
    list(wts = smryDf, ests = estsMat)
  } else {
    smryDf
  }
}
# test <- fitEvoMods(phy = phyTrim, m = tipLchron[[10]][['m']], 
#                    err = tipLchron[[10]][['se']], params = TRUE, nBoot = 50)

pt1 <- proc.time()
nCore <- detectCores() - 1
registerDoParallel(nCore)
modsLchron <- foreach(x = tipLchron, .packages = c('geiger', 'pmc'),
                      .combine = rbind, .inorder = FALSE) %dopar%
  fitEvoMods(phy = phyTrim, m = x[['m']], err = x[['se']], 
             params = TRUE, nBoot = nBoot)
stopImplicitCluster()
pt2 <- proc.time()
(pt2 - pt1) / 60
# 35 min elapsed time for 50 replicates

pt3 <- proc.time()
registerDoParallel(nCore)
modsDfRdm <- foreach(x = tipLrdm, .packages = 'geiger',
                     .combine = rbind, .inorder = FALSE) %dopar%
  fitEvoMods(phy = phyTrim, m = x[['m']], err = x[['se']])
stopImplicitCluster()
pt4 <- proc.time()
(pt4 - pt3) / 60 

# Charts for chronological sampling ---------------------------------------

# reformat model output (list) into dataframes - here, containing AICc weights
# there's surely a tidier way to do this, soz reader
getDf <- function(x) t(x$wts)
modsDfM <- sapply(modsLchron, getDf)
modsDfChron <- data.frame(t(modsDfM))
colnames(modsDfChron) <- mods
bestPos <- apply(modsDfChron, 1, which.max)
modsDfChron$bestMod <- mods[bestPos]
binsInPlot <- row.names(modsDfChron)
modsDfChron$bin <- factor(binsInPlot, levels = rev(binsInPlot)) 
# bins shouldn't be numeric since plotted as discrete axis
# put in reverse order so will plot youngest to oldest from top to bottom

# *Stacked barplot --------------------------------------------------------

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

# *Akaike weights table ---------------------------------------------------

# export the model support for each time bin as a suppl table:
# some values are too similar to compare easily form the barplot alone

wts4tbl <- modsDfChron[, c('bin', 'bestMod', mods)]
chronTbl <- xtable(wts4tbl, align = c('r', 'r', 'l', rep('r', length(mods))), 
                   digits = 3,
                   caption = 'Trait evolution model weights by time bin')
if (ss){
  tblNm <- paste0('Figs/phylo-evo-model-wts_chron_SS_',  day, '.tex')
} else {
  tblNm <- paste0('Figs/phylo-evo-model-wts_chron_hab_', day, '.tex')
}
if (ss){
  print(chronTbl, file = tblNm, include.rownames = FALSE,
        caption.placement = 'top')
} else {
  print(chronTbl, file = tblNm, include.rownames = FALSE,
        caption.placement = 'top')
}

# *Param estimates --------------------------------------------------------

# report CIs from bootstrapping the lambda, delta, and OU models

for (paramNm in c('lambda','delta','alpha')){
  getParam <- function(x) x$ests[paramNm,]
  paramMat <- sapply(modsLchron, getParam)
  paramDf <- data.frame(t(paramMat))
  
  if (paramNm == 'lambda'){ 
    # estimates VERY close to zero, need sci notation
    paramDf$observed <- formatC(paramDf$observed, format = "e", digits = 2)
    paramDf$X2.5. <-    formatC(paramDf$X2.5.,    format = "e", digits = 2)
    paramDf$X97.5. <- round(paramDf$X97.5, 3)
  }
  paramTbl <- xtable(paramDf, align = rep('r', 4), 
                     digits = 3,
                     caption = 'Estimate and bootstrapped CI')
  if (ss){
    outTblNm <- paste0('Figs/', paramNm, '_chron_SS_', nBoot, 'x_',  day, '.tex')
  } else {
    outTblNm <- paste0('Figs/', paramNm, '_chron_hab_', nBoot, 'x_', day, '.tex')
  }
  if (ss){
    print(paramTbl, file = outTblNm, include.rownames = TRUE,
          caption.placement = 'top')
  } else {
    print(paramTbl, file = outTblNm, include.rownames = TRUE,
          caption.placement = 'top')
  }
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


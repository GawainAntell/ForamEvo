library('pmc')
library('geiger')
library('ouch')
library('pracma')
library('phytools')
library('picante') # for plot.color.phylo
library('RColorBrewer')
library('ggplot2')
library('cowplot')

# Bootstrap dist concept --------------------------------------------------

# replicate reported p-value and power from Boettiger et al. 2012
# and plot a conceptual figure to show how they're derived

cores <- 1 # can't run parallel cores on Windows
nboot <- 3000 

# data from pmc vignette
data(anoles)
tree <- with(anoles, ouchtree(node, ancestor, time / max(time), species))
ou0v3 <- pmc(tree, log(anoles['size']), 'brown', 'hansen', 
             list(), 
             optionsB = list(regimes = anoles['OU.LP'], sqrt.alpha = 1, sigma = 1),
             nboot = nboot,  mc.cores = cores)

# test statistic estimated from observations
obsDlta <- ou0v3$lr

# get approximating function for each distribution
densH0 <- density(ou0v3$null)
densH1 <- density(ou0v3$test)
funH0 <- approxfun(densH0$x, densH0$y)
funH1 <- approxfun(densH1$x, densH1$y)

# integrate null from est to upper bound
pVal <- integral(funH0, obsDlta, max(densH0$x))

# integrate test from 95% quantile of null to upper bound
n95 <- quantile(ou0v3$null, 0.95)
pwr <- integral(funH1, n95, max(densH1$x))

pVal; pwr
# [1] 0.02659013
# [1] 0.9330232
# values reported in paper: 0.025, 93.6%

mcDf <- data.frame(
  delta = c(ou0v3$null, ou0v3$test),
  mod = c(rep('null', nboot), rep('test', nboot))
)
xMin <- min(mcDf$delta)
xMax <- max(mcDf$delta)
buffr <- (xMax - xMin) * 0.2

p1 <- ggplot(mcDf, aes(x = delta)) +
  theme_classic() +
  scale_x_continuous(limits = c(xMin - buffr, xMax + buffr), 
                     expand = c(0,0)) +
  scale_y_continuous('Density', expand = c(0,0), limits = c(0,.1)) +
  geom_density(aes(color = mod)) + #alpha = 0.5
  geom_vline(xintercept = obsDlta) +
  theme(axis.title.x = element_blank()) +
  guides(color = guide_legend(title = 'Model'))

# shade the portion corresponding to the p-value
h0df <- data.frame(x = densH0$x, y = densH0$y)
pvalPlot <- p1 +
  geom_ribbon(data = subset(h0df, x > obsDlta), 
              aes(x=x, ymin = 0, ymax=y), fill = 'red', alpha = .3) +
  theme(legend.position = 'none')

# shade the portion corresponding to the statistical power
h1df <- data.frame(x = densH1$x, y = densH1$y)
pwrPlot <- p1 +
  #  geom_vline(xintercept = n95, color = 'red', linetype = 'dashed') +
  geom_ribbon(data = subset(h1df, x > n95), 
              aes(x=x, ymin = 0, ymax=y), fill = '#00BFC4', alpha = .3) +
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(0.8, 0.8))

plotRow <- plot_grid(pvalPlot, pwrPlot, labels = c('A', 'B'),
                     hjust = c(-0.5, 1),
                     rel_widths = c(1.25, 1))

xlbl <- ggdraw() + 
  draw_label('Likelihood ratio (delta)')
plotAll <- plot_grid(
  plotRow, xlbl,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(1, 0.1)
)

pdf('Figs/bootstrap-dists-power-and-pval-concept.pdf', width = 5, height = 3)
plotAll
dev.off()


# Plot timescaled phylo ---------------------------------------------------

spAttr <- read.csv('Data/foram-spp-data_2020-11-15.csv')
df <- read.csv('Data/niche-sumry-metrics_SJ-ste_SS_2020-11-15.csv')
df$sp <- gsub(' ', '_', df$sp)
spAttr$species <- gsub(' ', '_', spAttr$species)

phyFull <- readRDS('Data/Aze-tree-phylo-object.rds')

# remove 2 species not in phylogeny
toss <- name.check(phyFull, df$sp, data.names = df$sp)
noPhy <- unique(toss$data_not_tree)
paste(paste(noPhy, collapse = ' '), 'not in tree')
if (length(noPhy) > 0){
  rows2toss <- df$sp %in% noPhy
  df <- df[!rows2toss,]
}
spp <- unique(df$sp)

# drop tips not sampled in niche data
phyTrim <- keep.tip(phyFull, spp)

keepRows <- spAttr$species %in% spp
spDf <- spAttr[keepRows,]

# Combine with species attribute data
# set open ocean thermocline species as the reference level
spDf$eco <- factor(spDf$eco, levels=c(3,1,2,4,5))
ecoLbl <- c('Thermocline','Mix layer symbiotic','Mix layer heterotroph',
            'Subthermocline','High latitude')
# From Aze et al. 2011, ecotype codes are:
# 1 = open ocean, mixed layer, trop/subtrop, w symbionts
# 2 = open ocean, mixed layer, trop/subtrop, w/o symbionts
# 3 = open ocean, thermocline
# 4 = open ocean, sub-thermocline
# 5 = high lat
# 6 = upwelling/high productivity
code2habitat <- function(x){
  switch(paste(x), 
         '1' = 'Mix layer symbiotic', 
         '2' = 'Mix layer heterotroph', 
         '3' = 'Thermocline', 
         '4' = 'Subthermocline', 
         '5' = 'High latitude', 
         'NA' = '-')
} 
spDf$ecoNm <- sapply(spDf$eco, code2habitat)
spDf$ecoNm <- factor(spDf$ecoNm)

# colorblind friendly diverging palette
clrs <- brewer.pal(5, 'Set1')

#pNm <- paste0('Figs/phylo-with-traits_',day,'.pdf')
#pdf(pNm, width = 7.5, height = 7.5)
color.plot.phylo(phyTrim, df = spDf, trait = 'ecoNm', taxa.names = 'species', 
                 main = '', leg.title = 'Ecotype', col.names = clrs,
                 x.lim = c(-70, 140),
                 y.lim = c(-0.2*Ntip(phyTrim), Ntip(phyTrim)), # room for legend
                 label.offset = 1
                 )
geo.legend()
colors = rep('white',3)
#dev.off()

library('pmc')
library('geiger')
library('ouch')
library('pracma')
library('ggplot2')
library('cowplot')

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


# plot timescaled phylo ---------------------------------------------------
# this code is notes on what to try; not run

# picante::color.plot.phylo
# ape::plot.phylo

# a projection of the tree in a space defined by phenotype (y axis) and time (on x)
# phytools::phenogram
# variant oto visualize the uncertainty around the reconstructed values at internal nodes
# phytools::fancyTree type='phenogram95'


tree<-pbtree(b=0.03,d=0.01,n=200)
#h<-max(nodeHeights(tree))
plotTree(tree,plot=FALSE)
# obj<-geo.legend(alpha=0.3,cex=1.2,plot=FALSE)
# obj$leg<-h-obj$leg
plotTree(tree,ftype="off",ylim=c(-0.2*Ntip(tree),Ntip(tree)),lwd=1)
geo.legend()

ggtree(unitTree) + 
  theme_tree2()


# Author: Stephan Schlickeiser
# Contact: stephan.schlickeiser@charite.de
# Created: 21 June 2017
# Last modified: 22 June 2018

# This script reproduces Figures 5_a_b_c_d_f, 6_b_c, 7_b_c_d_e, and 8_h of the paper
# "Human microglia regional heterogeneity and phenotypes determined by
# multiplexed single-cell mass cytometry" (https://doi.org/10.1038/s41593-018-0290-2)

# ******************************************************************************
#                     download / install / library 
# ******************************************************************************

# 1) go to
# https://flowrepository.org/id/FR-FCM-ZYM6
# and download zip folder with fcs raw data

# 2) unzip and move folder (should be 'FlowRepository_FR-FCM-ZYM6_files')
# to your working directory

# 3) set path to that folder (or set full path)
path_to_fcs <- "FlowRepository_FR-FCM-ZYM6_files"

# 4) install dependencies and required packages (if you haven't already)

# dependencies
install.packages(c("ks", "dbscan", "feature", "shiny", "rmarkdown", "rgl", 
  "knitr", "robust", "splancs", "flowCore", "lpSolve", "ggplot2", "gridExtra", 
  "grid", "matrixcalc"))

# required to reproduce figures
install.packages(c("devtools", "R.utils", "lme4", "lmerTest", "flowFP", "vegan"))

# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c("flowWorkspace", "CytoML", "flowType", "RchyOptimyx", "edgeR", "cydar"))

# use install.packages("emdist") if 0.3-2 doesn’t work 
library(devtools)
install_github("s-u/emdist") # version 0.3-2
install_github("steschlick/Boettcher_et_al_NN_2018/NNhuMG")

library(NNhuMG)

# ******************************************************************************

# define two color ramps
# brownish colors
brwn.colors <- colorRampPalette(c("white", "bisque", "orange", 
                                  "darkorange", "darkred"))
# Matlab's 'jet.colors'
jet.colors <- colorRampPalette(c("#007FFF", "cyan","#7FFF7F", "yellow", 
                                 "#FF7F00", "red", "#7F0000"))


# ******************************************************************************
#          fcs data / replicate Cytobank (CB) transformation and gating
# ******************************************************************************
library(flowWorkspace)
library(CytoML)

# get the path to each fcs file 
fcs.files <- dir(path=path_to_fcs, full.names = FALSE, pattern = ".fcs$")

# create GatingSet with Cytobank gatingML file
xml.manual.cutoffs <- system.file("extdata/CytExp_149973_Gates_v6.xml", package = "NNhuMG", mustWork=TRUE)

gs <- cytobank2GatingSet(xml.manual.cutoffs, paste0(path_to_fcs,"/",fcs.files))

gs.name <- unname(parameters(flowData(gs)[[1]])$name)
gs.desc <- unname(parameters(flowData(gs)[[1]])$desc)

getNodes(gs) # manual 1D-gates to define marker positivity
sampleNames(gs)

# need to update to CB ranges for correct display
data(bounds)

for (i in 1:length(gs)) {
  parameters(flowData(gs)[[i]])$minRange <- bounds[,"minimum"]
  parameters(flowData(gs)[[i]])$maxRange <- bounds[,"maximum"]
}

# ******************************************************************************
#            fcs data /  prettify names / prepare data for plotting
# ******************************************************************************

# get a flowSet (copy) of the (transformed) data
fs <- flowData(gs[!grepl("concatenated", sampleNames(gs))])
fs@file # ncdfFlowSet

# keep parameter description
fs.desc <- parameters(fs[[1]])@data[,1:2]

# indices of markers used in tSNE embedding
visne.ch.idx <- which(grepl(" (v)",parameters(fs[[1]])$desc, fixed=TRUE))

# make marker names more readable and remove unwanted chars
p.desc <- unname(parameters(fs[[1]])$desc)
p.desc[visne.ch.idx] <- sapply(strsplit(p.desc[visne.ch.idx], "_"), '[[', 2)
p.desc[visne.ch.idx] <- gsub("_|-|\\+|&|/|\\|| \\(v\\)", "", p.desc[visne.ch.idx])
p.desc[p.desc=="CCl2"] <- "MCP1"

# update flowSet with pretty marker names
for (f in 1:length(fs)) {
  parameters(fs[[f]])$desc <- p.desc
}

# update parameter description	
fs.desc <- cbind(fs.desc, p.desc, visne=logical(nrow(fs.desc)))
fs.desc$visne[visne.ch.idx] <- TRUE

# get study meta data
data(study)

# make sure to have correct order
if (!all(sampleNames(fs) == study$fcs)) {
  study <- study[match(sampleNames(fs), study$fcs), ]	
}

(sample.names <- as.character(study$Sample))
# more readable
sampleNames(fs) <- sample.names 
pData(fs)$name <- sample.names
sampleNames(gs)[!grepl("concatenated", sampleNames(gs))] <- sample.names

# get an overview and an estimate of the... 
(n.fr <- length(fs)) # ...number of samples
(v.events <- fsApply(fs, nrow)) # ...number of events per sample
(min.events <- min(v.events)) # ...minimum number of events

# collapse events of all files for all parameters 
# into single flowFrame, for plotting
merge.fr <- fs[[1]]
merge.id <- rep(1,v.events[1, ])
for (i in 2:length(fs)){
  exprs(merge.fr) <- rbind(exprs(merge.fr),exprs(fs[[i]]))
  merge.id <- c(merge.id,rep(i,v.events[i, ]))
}

# random index for plotting 2K events per sample (~70K in total)
event.id <- 1:length(merge.id)
ran <- split(event.id, merge.id)
ran <- sapply(ran, function(x) {
              set.seed(25)
              size <- min(length(x), 2000)
              x[sample.int(length(x), size)]	
              }, simplify=FALSE)
ran <- unlist(ran)
set.seed(77)
ran <- ran[sample.int(length(ran), length(ran))]

# sampled data for plotting tSNE maps
tsne.plot.data <- exprs(merge.fr)[ran,c("tSNE1","tSNE2")]
xlim <- range(tsne.plot.data[,"tSNE1"])
ylim <- range(tsne.plot.data[,"tSNE2"])

# ******************************************************************************
#                      univariate analysis
# ******************************************************************************
library(lme4)
library(lmerTest)

# Hypothesis: human microglia isolated from different brain regions
# differs in marker expression levels

levels(study$Region)
levels(study$Donor)
table(study[ , c("Donor", "Region")]) # 9 missing values

dir.create("plots")
pdf("plots/marker_expression_levels.pdf", width=15, height=45)
  ord <- order(study$Region)
  densityplot(factor(name, levels = (name[ord])) ~ ., fs, 
              channels=gs.name[visne.ch.idx], overlap=0.5, alpha=0.8, col="black", 
              fill=as.character(study$region.col)[ord], darg=list(adjust=0.5))
dev.off()


# color code for tsne plots
col.donor <- as.character(study$donor.col)[merge.id[ran]]
col.region <- as.character(study$region.col)[merge.id[ran]]

# ******************* #
#  Figures 5a and 5b  #
# ******************* #

# donor color code
png("plots/figure_5_a.png", width=4.5, height=5, units="in", res=300)
  plot(tsne.plot.data, pch=19, cex=0.075, col=col.donor, xlim=xlim, ylim=ylim, 
       xlab="tSNE1", ylab="tSNE2", main="")
dev.off()

# region color code
png("plots/figure_5_b.png", width=4.5, height=5, units="in", res=300)
  plot(tsne.plot.data, pch=19, cex=0.075, col=col.region, xlim=xlim, ylim=ylim, 
       xlab="tSNE1", ylab="tSNE2", main="")
dev.off()


# optionally: exclude markers that are negative for most of the cells
# e.g. on average at least 1 % should be above background level
freq.limit <- 0.01
bg.cutoff <- asinh(5/5) # arbitrary (half of the 1st log decade approximately)
exp.marker.q <- apply(exprs(merge.fr)[ ,visne.ch.idx ], 2, quantile, probs=1-freq.limit)
idx.keep <- exp.marker.q > bg.cutoff

fs.means <- fsApply(fs, each_col, mean)[ , visne.ch.idx[idx.keep]]
colnames(fs.means) <- p.desc[visne.ch.idx[idx.keep]]


# Since there are 9 missing values and small sample size
# we use Skillings-Mack nonparametric statistic for imbalanced block designs
# or mixed-effect model with donor treated as random effect
# to test for regional differences 

# number of permutations (sm.test, adonis)
B <- 1e4

# remove block #12, since it has only one observation
rmvd <- study$Donor != "#12"

# Skillings-Mack with Monte Carlo-approximated p-values (takes ~5')
set.seed(25)
fs.means.sm <- sm.test(fs.means[rmvd,], study$Region[rmvd], study$Donor[rmvd], 
                       simulate.p.value=TRUE, B=B)
names(fs.means.sm) <- colnames(fs.means)
sort(fs.means.sm)
fs.means.robust.pvals <- cbind(SM.p.value=fs.means.sm, 
                               SM.p.value.BH=p.adjust(fs.means.sm,"BH") )

# model between-subject variances as random effect:
fs.means.lme <- sapply(1:ncol(fs.means), function(i) {
				df <- data.frame(Marker=fs.means[rmvd,i], 
				                 study[rmvd , c("Donor", "Region")])
				lmer(Marker ~ Region + (1|Donor), data=df)
				}, simplify=FALSE)
names(fs.means.lme) <- colnames(fs.means)

# anova on fixed effect
fs.means.lme.aov <- sapply(fs.means.lme, anova)
fs.means.lme.aov <- cbind(t(fs.means.lme.aov) , 
                          p.value.BH=p.adjust(fs.means.lme.aov[6,], "BH") )
fs.means.lme.aov

# liklihood ratio test on random effect
fs.means.lme.rand <- sapply(fs.means.lme, rand)
fs.means.lme.rand <- do.call("rbind",fs.means.lme.rand)
fs.means.lme.rand 

# ANOVA posthoc: pairwise within group comparisons with Holm-correction  
fs.means.lme.post <- sapply(fs.means.lme, function(x){ 
								diffl <- difflsmeans(x, test.effs = "Region")[[1]]
								diffl[,7] <- p.adjust(diffl[,7],"holm")
								diffl}, simplify=FALSE)
fs.means.lme.post <- do.call("rbind",fs.means.lme.post)
fs.means.lme.post <- cbind(fs.means.lme.post, 
                           p.value.BH=p.adjust(fs.means.lme.post[,7], "BH") )
fs.means.lme.post

# clean-up
rm(list=c("fs.means.lme", "fs.means.lme.rand","fs.means", 
          "fs.means.lme.aov", "fs.means.lme.post"))
detach("package:lmerTest", unload=TRUE)
detach("package:lme4", unload=TRUE)
gc(); gcDLLs() 

# ******************************************************************************
#                      probability binning
# ******************************************************************************
library(flowFP)
# browseVignettes("flowFP")

# number of recursions (9), for minimum average of e.g. 10 cells per bin per sample
nRec <-  min(floor(log2(mean(v.events) / 10)), 11) # maximum 2^11 bins

2^nRec # number of bins

# generate single binning model (grid) over tSNE space
set.seed(1406)
mod <- flowFPModel(as.flowSet(fs), name="tSNE space", 
                   parameters=c("tSNE1","tSNE2"), nRecursions=nRec)

# apply binning model to samples, i.e. to count cells in each bin
fp <- flowFP(fs,  mod)

# get counts and bin assigment for each event and each sample
fp.res <- counts(fp)
fp.tags <- tags(fp)

# apply model to merged samples to calculate average tSNE locations
fp.merge <- flowFP(merge.fr,  mod)

# get bin assigment
fp.merge.tags <- unlist(tags(fp.merge))

# calculate bin-wise median expression for each parameter
pars.idx <- c(visne.ch.idx, match(c("tSNE1","tSNE2"), colnames(merge.fr)))
pars.loc <- matrix(NA, ncol=length(pars.idx), nrow=ncol(fp.res))
colnames(pars.loc) <- parameters(merge.fr)$desc[pars.idx]

for (i in 1:ncol(fp.res)) {
  pars.loc[i,] <- apply(exprs(merge.fr)[fp.merge.tags %in% i, pars.idx],
                        2, median, na.rm=TRUE)
}

tsne.loc <- pars.loc[ , c("tSNE1","tSNE2")] # keep tSNE x-y locations separate
pars.loc <- pars.loc[ , !colnames(pars.loc) %in% c("tSNE1","tSNE2")]

# the frequency data
freq.tab <- t(fp.res)
colnames(freq.tab) <- study$Sample

# normalize, i.e. calculate cell frequencies sample-wise for each bin
freq.tab <- sweep(freq.tab, 2, apply(freq.tab, 2, sum), '/')

# *********** #
#  Figure 5d  #
# *********** #

# the 'binning grid' comprising 512 'micro-gates'
pdf("plots/figure_5_d.pdf",width=4.5,height=5)
  plot(mod, alpha=0.0, border="darkblue", lwd=0.5, 
       xlim=xlim, ylim=ylim, main="", cex=2)
dev.off()

# ******************************************************************************
#  earth mover's distance-based permutational multivariate analysis of variance 
# ******************************************************************************
library(emdist)
library(vegan)

# compute sample-pairwise EMD scores
emd.v.idx <- combn(1:n.fr, 2, simplify=FALSE)
emd <- matrix(0, n.fr, n.fr)
dimnames(emd) <- list(sample.names, sample.names)

# # as this takes > 10', we load precomputed scores
# for (k in 1:length(emd.v.idx)) {
  # emd.tmp <- emdr(A = cbind(freq.tab[,emd.v.idx[[k]][1]], tsne.loc), 
                  # B = cbind(freq.tab[,emd.v.idx[[k]][2]], tsne.loc),  
                  # dist="euclidean",flows=FALSE, max.iter=10000L)
                  # emd[emd.v.idx[[k]][1], emd.v.idx[[k]][2]] <- emd.tmp
                  # emd[emd.v.idx[[k]][2], emd.v.idx[[k]][1]] <- emd.tmp
  # }
data(emd)

# EMD-based permutational analysis of variance to assess 'global' differences
# in tSNE1/2-histogram cell probability distributions
dmd <- as.dist(emd) # earth mover's distance matrix
adonis(dmd ~ Donor + Region, data=study, permutations=B)

# constrain permutations within Donor (strata) to test for Region difference
adonis(dmd ~ Region, strata=study$Donor, data=study, permutations=B)

# plot heatmap of sample-pairwise EMD scores
hmd <- hclust(dmd, method="ward.D2")

# color code for side bars
col.r <- as.character(study$region.col)
col.d <- as.character(study$donor.col)

# *********** #
#  Figure 5f  #
# *********** #

# heatmap of pairwise EMD scores
png("plots/figure_5_f.png", width = 1600, height = 1600, res=300)
  heatmap(emd, Rowv=(as.dendrogram(hmd)), Colv=(as.dendrogram(hmd)), 
        ColSideColors=col.r, RowSideColors=col.d, col=brwn.colors(100), 
		scale="none", cexRow=.9, cexCol=.9)
dev.off()

# clean-up
rm(list=c("dmd", "hmd","emd"))
detach("package:vegan", unload=TRUE)
detach("package:emdist", unload=TRUE)
gc(); gcDLLs() 

# ******************************************************************************
#                      bin-wise testing
# ******************************************************************************

# in order to identify which cellular subsets account for observed 'global'
# difference, we test in each bin for group difference (study$Region)
# using Skillings-Mack nonparametric statistic for imbalanced block designs
# appropriate p-values are computed by Monte Carlo simulation

# sm.test(sim.SM=TRUE) returns the permutation stats to allow for estimating Tmax,  
# i.e. the maximal statistic for each permutation

# ## this takes >30', so we load precomputed data
# # set.seed(6519)
# sim <- sm.test(y=t(freq.tab[ , rmvd]), 
		       # groups=study$Region[rmvd], blocks=study$Donor[rmvd], 
		       # simulate.p.value = TRUE, B = B, sim.SM=TRUE)
# # group level statistic
# sim.T <- sim[[1]]$T
# # permutation null distribution
# sim.Null <- do.call("rbind", sim[[2]])
# # simulated p-values
# pvals <- apply(t(sim.Null) >= sim.T, 1, mean)
# # maximal statistic
# tmax <- apply(sim.Null, 1, max)

# simulated p-values
data(pvals)

# maximal statistic
data(tmax)


# T critical value is the .95 quantile of the maximal statistic,
# provides strong FWER control, albeit less conservative than Bonferroni
(Tcrit <- quantile(tmax, 0.95))
# Tmax formula given in 'Nonparametric Permutation Tests for Functional 
# Neuroimaging: A Primer with Examples. Human Brain Mapping, 2001, Nichols and Holm'
# sort(sim.Tmax, decreasing=TRUE)[floor(0.05 * B) + 1]

# corresponding p-value of Chi^2 distributed test statistic
pchisq(Tcrit, df=4, lower.tail=FALSE)

# need to substitute 0 p-values for -log-scaling and display
sim.pvals <- pvals
sim.pvals[sim.pvals == 0] <- 0.1/B

# convert p-values for display and weighing
w.sm <- -log2(sim.pvals)

# *********** #
#  Figure 6b  #
# *********** #

# shows p-values in 512 bins
png("plots/figure_6_b.png", width=4.5, height=5, res=300, units="in")
  plot(tsne.plot.data, pch=19, cex=0.15, col=colorDens(c(-log2(1e-6), 
     w.sm[fp.merge.tags][ran]))[-1], xlim=xlim, ylim=ylim, 
     yaxt="n", xaxt="n", xlab="", ylab="", main="", bty="n")
  par("new"=TRUE)
  plot(mod, alpha=0, border="blue", lwd=0.35, xlim=xlim, ylim=ylim, 
     main="", xlab="tSNE1", ylab="tSNE2", cex=2)
dev.off()

# ******************************************************************************
#            kernel density-based gating of significant subsets
# ******************************************************************************
library(ks)
library(rgl)

# first get a kd estimate of the composite tSNE map using the 
# binning grid locations as input observations (we found that using the 
# already binned tSNE-PB-grid is sufficient for getting representive estimates 
# provided it contains at least 512 bins, i.e. observations)

# parameters passed to ks::kde for further binning of the data to speed-up 
# computation time (we use bgridsize of 512 here for more accurate display,
# kde defaults to 151, which is much faster and gives essentially the same results)
# a smaller bgridsize will also speed-up loading and rendering of the 3d graphs
# in the browser when using shinyBW()
bgridsize <- rep(512, 2) 
binned <- TRUE

png("plots/choose_BW.png", width=9, height=6.5, res=300, units="in")
BW <- chooseBandWidth(data=tsne.loc, nH=100, binned=binned, 
                      bgridsize=bgridsize, store.fhats=TRUE, 
                      MAX.fine=35, MAX.coar=15, tol.spur=2, plotit=TRUE)
dev.off()

object.size(BW)

# make a judicious choice for a bandwidth
BW[3:5]
# shinyGUI for inspecting kde's and selecting a bandwidth manually
# requires variable 'BW' to be present globally
bw.idx <- shinyBW() 

# select 'H1' (BW index 9) for it captures fine gaps between the larger
# and the few smaller clusters
Hopt <- BW$H/BW$Hfactor[BW$Hindex[1]]
# # or manual selection
# Hopt <- BW$H/BW$Hfactor[bw.idx]

# density estimate of tSNE map, shown in Figure 5c
fh <- kde(tsne.loc, binned=TRUE, bgridsize=bgridsize, H=Hopt)
tobs <- predict(fh, x=tsne.plot.data) # assign densities to events

png("plots/figure_5_c.png", width=4.5, height=5, res=300, units="in")
 plot(tsne.plot.data, pch=19, cex=0.15, col=colorDens(tobs), 
      xlim=xlim, ylim=ylim, xlab="tSNE1", ylab="tSNE2", main="")
dev.off()

open3d()
persp3d(fh$eval.points[[1]], fh$eval.points[[2]], fh$estimate, 
        xlab="tSNE1", ylab="tSNE2", zlab="", aspect = c(1, 1, .25),
        alpha=1, col=colorDens(fh$estimate), box = TRUE, axes = TRUE) 

# convert FDR-adjusted p-values to weights such that weights are > 1 if 
# p <= 0.05 and sum up to number of observations
w <- threshold2one(-log2(p.adjust(sim.pvals, "BH")), -log2(0.05))
sum(w)
plot(sort(w))

# compute positive difference between the unweighted and the  
# 'significance'-weighted density estimates of the tSNE map
fo <- kdeDiff(tsne.loc, binned=TRUE, bgridsize=bgridsize, H=Hopt, w1=w)
dobs <- predict(fo, x=tsne.plot.data)

# spotting significant regions
plot(tsne.plot.data, pch=19, cex=0.15, col=colorDens(dobs), 
     xlim=xlim, ylim=ylim, xlab="tSNE1", ylab="tSNE2", main="")


# ******************************************************************************
# gating of 'dens' regions of significant bins

# set minimal required number of bins according to frequency limit
# and/or an average minimum reliable number of events we would target in flow/mass 
# cytometry (i.e. a population of 10 cells will per se have a coefficient of
# variation of at least 32%, 25 cells 20%, and 100 cells could 
# have an acceptable CV of 10% if there were no other variance than Poisson)

minBins <- ceiling(freq.limit * 2^nRec)
# Poisson precision
100 / sqrt(minBins * (mean(v.events) / 2^nRec)) # CV 
 
# this function is passed to contourLoc() in order to compute p-values
# in contour-gated subsets 
smfun <- function(data) {		
  pv <- sm.test(y=data, groups=study$Region[rmvd], 
                blocks=study$Donor[rmvd],
                simulate.p.value = FALSE)$p.value
  pv
}

# get a 'contour list object'
cs <- contourLoc(data=t(freq.tab[ , rmvd]), fhat=fo, offset=0.1, n.levels=50, 
                 alpha=0.005, statsfun=smfun, keep.level.by="max.bin.cor",
                 peaks2merge=4, minBins=minBins, plotit=TRUE)
# we have set alpha to 0.005 which is close to pchisq(Tcrit, df=4, lower.tail=FALSE)
# to select one of the smaller regions as well

# contours are labeled according to their peaks
plotContours(cs, col="grey", pch="+", cex=0.5, lwd=2.5, fhat.levels=F, labels="peak.ids")

# need to remove outer contours that have been reported due to presence of spurious peaks, 
# because we next want to compare phenotypes and decide whether to merge them or not
css <- shrinkSpur(cs, mul=T, minBins=10) 
plotContours(css, col="grey", pch="+", cex=0.5, lwd=2.5, fhat.levels=F, labels="peak.ids")

# set thresholds for:
# degree of overlap for a given marker expression between two subsets
beta <- 1 - pnorm(1) # 2-sigma difference
# robust effect size, 2 corresponds to 2-sigma difference between 2 normal distributions (Cohen's d)
delta <- 2 
# specify minimally required distance (at asinh scale) between two marker distributions
# this also prevents from computing artificially large effect sizes which happens if
# a subset has marker expression levels near or at zero
threshold <- 0.5

# merge phenotypically similar contours (i.e. caused by low bandwidth 
# or some variability in abundance between different samples)

# first, set keep.outer=T, keep.inner=T to show that...
csm <- contourMerge(css, keep.outer=T, keep.inner=T, data=pars.loc, 
                    delta=delta, beta=beta, threshold=0.5, minBins=minBins)
# ... this is effectively comparing (multiple) inner contours ('higher magnitude') 
# within an outer contour ('lower magnitude')
# if the exclusion region (filled in next plots, and labelled by '*') between inner and 
# outer is large enough (> minBins), it will also be considered for comparison
layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
plotContours(csm, comp.sel=T, col="grey", pch="+", cex=0.5, lwd=2.5, ask=FALSE,
             fhat.levels=F, labels="peak.ids", filwd=.5, density=50)
par(mfrow=c(1, 1))

# inner subsets do not differ phenotypically
# '*·root·...' specifies the region excluding all outmost contours
plotPheno(csm, comp.sel=T,
			relevant.only=T,
			plot.cutoffs=FALSE,
			delta.panel=T,
			max.delta.bars=4,
			ask=T)

# select final contours (all inner will be merged by comparison)
csm <- contourMerge(css, keep.outer=F, keep.inner=F, data=pars.loc, 
                    delta=delta, beta=beta, threshold=0.5, minBins=minBins)
plotContours(csm, comp.sel=F, col="grey", pch="+", cex=0.5, lwd=2.5,
             fhat.levels=F, labels="peak.ids", filwd=.5, density=50)

# show final comparisons
plotPheno(csm, comp.sel=T,
			relevant.only=T,
			plot.cutoffs=FALSE,
			delta.panel=T,
			max.delta.bars=4,
			ask=F)

# *********** #
#  Figure 6c  #
# *********** #

png("plots/figure_6_c.png", width=4.5, height=5, res=300, units="in")
  plot(tsne.plot.data, pch=19, cex=0.15, col=colorDens(dobs), 
      xlim=xlim, ylim=ylim, xlab="tSNE1", ylab="tSNE2", main="")
  for (b in csm$merged) {
    polygon(csm$res$boundary[[b]]$x,
    csm$res$boundary[[b]]$y,
    border="black", lwd=2.5)		
  }
dev.off()


# ******************************************************************************
# phenotyping of contour-'gated' subsets using manual cutoffs

# plot manually gated 1D-cutoffs which define marker positivity,
# gating has been performed in Cytobank on biaxial plots for better display
# exported gatingML was then used to replicate the gates in a GatingSet (see above)  
pdf("plots/manual.cutoffs.pdf", width=14, height=14)
 plotGate(gs[["concatenated.all.fcs"]])
dev.off()

# the marker against which all other markers had been plotted
yChnl <- "Ho165Di" # "PR2Y12"
# specify alternative marker if yChnl is included for gating
xyChnl <- "Ir193Di" # "193Ir_DNA2"

# extract gate boundaries from GatingSet
manual.defs <- sapply(getNodes(gs)[-1], function(v) {
	                  getGate(gs, v)[[1]]@boundaries}, simplify=FALSE)

manual.cutoffs <- sapply(manual.defs, cutoffExtract, yChnl, xyChnl, simplify=FALSE)
names(manual.cutoffs) <- sapply(manual.cutoffs, names)

# check whether all markers used in viSNE had been used for manual cutoff-gating
check.cutoffs <- visne.ch.idx %in% sort(match(unique(names(manual.cutoffs)), gs.name))
if (!all(check.cutoffs)) {
	cat("no cutoff(s) defined for following marker(s): \n",
		paste(gs.desc[visne.ch.idx[!check.cutoffs]], "\n") )
	}
# 
manual.cutoffs <- sapply(gs.name[visne.ch.idx], function(x) {
					sort(unname(unlist(manual.cutoffs[names(manual.cutoffs) %in% x])))
					}, simplify=FALSE)  

# cutoffs have channel names, switch to pretty marker names
names(manual.cutoffs) <- switch.names(as.matrix(fs.desc[,1:3]), 
                                      names(manual.cutoffs), "pretty.desc")
# print manual phenotype definitions
manual.phenotypes <- cuts2Phenotype(csm, manual.cutoffs, verbose=TRUE)

# ******************************************************************************
# phenotyping of contour-'gated' subsets using data-driven cutoffs

# find.all=TRUE disregards delta and all cutoffs are computed according to beta 
# (and a threshold to limit the distance between any two cutoffs), otherwise 
# subsets with low delta are being excluded from comparison and hence cutoffs
# (between the other subsets) might cut those 
csm <- findCutoffs(csm, delta=delta, beta=beta, threshold=0.25, 
                   find.all=T, recomp=F, max.cuts=1, verbose=TRUE)
plotPheno(csm, comp.sel=csm$subset, relevant.only=T,
          plot.cutoffs=T, delta.panel=F, max.delta.bars=3, ask=F)

# in order to find best cutoff combination, provide unbinned data
# and see whether they still hit their targets
merge.data <- exprs(merge.fr)[ran,]
colnames(merge.data) <- p.desc

# this copies the content of csm and samples 512 events for each subset from
# expression data
set.seed(123)
csr <- exprs2con(con=csm, data=merge.data, markers=visne.ch.idx, dwnsmpl=512)

# marker expression is less distinct between subsets than in binned data 
# (where expression medians had been calculated)
plotPheno(csr, comp.sel=csr$subset, relevant.only=T,
          plot.cutoffs=T, delta.panel=F, max.delta.bars=3, ask=F)


# phenoReduce() solves an Integer Linear Programm for a minimal and/or 
# required number of markers, if target.num.marker is specified, the sum 
# of the per-marker-sums of effective differences delta is maximized
# 
# goal is to compute cutoffs to find a projection that targets the subsets
# in a binary (sequential/boolean) gating strategy
#  
# we require 8 markers to have some options in a FACS panel
csr <- phenoReduce(csr, target.num.marker=8, b=0.5)

# updated cutoffs 
plotPheno(csr, comp.sel=csr$subset, relevant.only=T,
          plot.cutoffs=T, delta.panel=F, max.delta.bars=3, ask=F)


# ******************************************************************************
# subset gating and visualization via flowWorkspace/CytoML

# extract phenotype definitions
(targets <- csr$targets)
# 'f·' specifies the 'full' target subset phenotype for required number of markers
# 'm·' specifies the phenotypes for a 'minimal' set of markers, selected based on 
# best average F-scores

# extract cutoffs and tSNE subset polygons
cutoffs <- csr$cutoffs
contours <- csr$res$boundary[csr$merged]
target.names <- names(contours)

# target phenocodes for 1D gates
target.1d <- unique(unlist(targets))
target.1d <- target.1d[target.1d != ""]

# get respective marker and channel names to
# create flowCore filters from cutoffs
rownames(bounds) <- p.desc
markers <- names(cutoffs)
channels <- switch.names(as.matrix(fs.desc[,1:3]), markers, "name")

gates1d <- sapply(1:length(cutoffs), function(i) {
                  get1dGates(channels[i], markers[i], 
                             cutoffs[[i]], bounds[markers[i], ])
                  }, simplify=FALSE)
gates1d <- unlist(gates1d, recursive=FALSE)

# remove all manual gates (leave root)
remove.that <- getNodes(gs)[-1] 
for (r in remove.that) {
  Rm(r,gs)
}
getNodes(gs)

# add the new 1D-gates to gatingSet
for (s in target.1d) {
  add(gs, gates1d[[s]])
}

# set boolean gates according to target phenotypes
targets <- targets[!duplicated(targets)]
for (c in 1:length(targets)) {
  cn <- paste("/", targets[[c]], sep="")
  cn <- cn[cn!="/"]
  exp <- paste(cn, collapse=" & ")
  bf <- boolFilter(exp, filterId = paste0("1d_",names(targets)[c],
                   "_", paste(targets[[c]], collapse="")))
  add(gs, bf)
}
recompute(gs)


# create flowCore polygon filters for selected contours
con.gates <- sapply(seq_along(contours), function(i) {
                    filterId <- paste0("pl_", names(contours[i]))
                    boundary <- contours[[i]]
                    polygonGate(filterId=filterId, 
                                .gate=cbind(tSNE1=boundary$x,tSNE2=boundary$y))
                    })

# add the gates to gatingSet
for (c in 1:length(con.gates)) {
  add(gs, con.gates[[c]])
}
recompute(gs)
# ready!


# compare pop stats, ie. fequencies
freq.gs <- t(getPopStats(gs, statistic = "freq",format="wide"))
freq.gs <- freq.gs[!grepl("concatenated", rownames(freq.gs)), ] 
freq.gs <- freq.gs[as.character(study$Sample), ]

freq.tsne <- freq.gs[ , grepl("pl_",colnames(freq.gs),fixed=T)]
freq.tsne <- freq.tsne[ , match(target.names, sub("pl_", "", colnames(freq.tsne)))]
freq.bool <- freq.gs[ , grepl("1d_f·",colnames(freq.gs),fixed=T)]
freq.bool <- freq.bool[ , match(target.names, 
                                sapply(strsplit(sub("1d_f·", "", colnames(freq.bool)), "_"),
                                '[[', 1))]

# Skillings-Mack on tSNE-gated frequencies
set.seed(1968)
tsne.pvals <- sm.test(y=freq.tsne[rmvd, ], 
                      groups=study$Region[rmvd], blocks=study$Donor[rmvd], 
                      simulate.p.value = T, B = B, sim.SM=FALSE)
# Skillings-Mack on boolean-gated frequencies
bool.pvals <- sm.test(y=freq.bool[rmvd, ], 
                      groups=study$Region[rmvd], blocks=study$Donor[rmvd], 
                      simulate.p.value = T, B = B, sim.SM=FALSE)

matrix(c(tsne.pvals, bool.pvals), ncol=2, 
       dimnames=list(target.names, c("tSNE","Boolean")))


# ******************* #
#  Figures 5f and 7c  #
# ******************* #

# plot the contour-gated subsets
png("plots/figure_5_f.png", width=6.5, height=4.5 , res=300, units="in")
  plotGate(gs[2:6], colnames(freq.tsne), smooth= FALSE, xbin = 128, stats = FALSE,
         par.settings = list(panel.background = list(col = "grey95"), 
                             background=list(col="white"),
                             gate = list(col = "red", lwd = 3),
                             gate.text = list( cex = 1 )
                             ), key=FALSE,
         arrange.main="", xlab="tSNE1", ylab="tSNE2", xlim=1.2*xlim, ylim=1.2*ylim)
dev.off()

# plot the boolean gated subsets color-overlaid on tSNE map
png("plots/figure_7_c.png", width=6.5, height=4.5 , res=300, units="in")
  plotGate(gs[1:6], colnames(freq.tsne), smooth= FALSE, xbin = 128, stats = FALSE,
         overlay = colnames(freq.bool)[c(1,3,2,4)],
         par.settings = list(overlay.symbol = list(cex = 0.02, bg.alpha = 0.25),
                             panel.background = list(col = "grey95"),
                             background=list(col="white"),
                             gate = list(col = "transparent")
                             ), key=FALSE,
         arrange.main="", xlab="tSNE1", ylab="tSNE2", xlim=1.2*xlim, ylim=1.2*ylim,
         overlay.symbol = setNames(list(
                                        list(fill = "violetred"),
                                        list(fill = "orange"),
                                        list(fill = "green3"),
                                        list(fill = "red")), 
                                   colnames(freq.bool))[c(1,3,2,4)] )
dev.off()

# *********** #
#  Figure 7d  #
# *********** #
library(ggplot2)

# compare frequencies of boolean and tSNE-gated subsets
se.col <- c("violetred", "orange", "green3", "red")

# compute spearman correlation and report rho in figure
rho_eqn <- function(df){
  eq <- substitute(~~italic("rho")~"="~r, 
    list(r = format(cor.test(df$y , df$x , method="spearman")[[4]], digits = 3)))
    as.character(as.expression(eq))                 
}
# report r squared as well
lm_eqn <- function(df){
    m <- lm(y ~ x, df)
    eq <- substitute(~~italic(r)^2~"="~r2, 
         list( r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq))                 
}
for (z in 1:4) {
  df <- data.frame(x = 100 * freq.tsne[ , z], y = 100 * freq.bool[ , z])
  p <- ggplot(data = df, aes(x = x, y = y)) +
           geom_point(colour=se.col[z]) +
           geom_smooth(method = "lm", se=T, colour=se.col[z] ,fill=se.col[z], formula = y ~ x ) +
           xlab("")  + ylab("")  + 
           geom_text(x = 0.15 * max(df$x), y = 0.85 * max(df$y), label = rho_eqn(df), parse = TRUE) +
		   geom_text(x = 0.15 * max(df$x), y = 0.75 * max(df$y), label = lm_eqn(df), parse = TRUE) +
           xlab("t-SNE polygon gates [%]")  + ylab("boolean phenotype gates [%]") +
  theme(legend.title=element_blank(),panel.background = element_rect(fill="grey95")) 
  pdf(paste0("plots/figure_7_d_",z,".pdf"),width=3.8,height=3.5)
    plot(p)
  dev.off()
}

# *********** #
#  Figure 7b  #
# *********** #

# boxplot showing subset's bin-median marker expression and cutoffs

# collect cutoff values
dat.thrlds <- data.frame(values=unlist(cutoffs), ind=names(cutoffs)) # this works for single cuts only!

# collect median expression levels, use bins from csm 
data <- csm$comb$data[, names(cutoffs)]
  dat <- data.frame(check.names=F,check.rows=F)
  for (cl in csm$subset){ 	
    subs.id <- rownames(csm$res)[cl]
	subset <- c(1:nrow(data)) %in% csm$res$bin.idx[[cl]]
	dat.tmp <- data.frame(subset = rep(subs.id, ncol(data)), 
	                      stack(as.data.frame(data[subset, ])))
    dat <- rbind(dat, dat.tmp)
  }  
  # set order of markers
  ord <- order(colnames(data))
  dat$ind <- factor(dat$ind, levels = colnames(data)[ord])

# plotit
g1 <- ggplot(data = dat) + 
  geom_boxplot(aes(x = ind, y = values, fill = subset), 
               outlier.size = .1, outlier.shape = ".") +
  geom_point(data = dat.thrlds, aes(x = ind, y = values),
             color = "blue",shape="_",size=10) +
  scale_fill_manual(values = c("8" = "violetred",
                            "16·17" = "orange", 
                            "13·14·15" = "green3",
                            "20·21·22" = "red",
                            "*·root·1·2·3·4·5·6·7·9·10·11·12·18·19·23" = "steelblue")) +
  theme(axis.title.y=element_text(size=12), axis.text.x = element_text(face="bold"), 
        axis.text.y = element_text(size=8), legend.position = "none", 
        panel.background = element_rect(fill="grey95")) + 
  xlab("") + ylab("MSI")
  plot(g1) 
ggsave("plots/figure_7_b.pdf", width = 28, height = 7, units="cm", g1)


# ******************************************************************************
#                   RchyOptimyx marker importance analysis
# ******************************************************************************
#    As variability is inherently fixed at the level of the embedding 
#    (as would be in any clustering) of high-dimensional data, detection 
#    of multiple subsets (as 'correlates') might be caused by only a
#    limited number of clinically/biologically relevant markers. This needs
#    to be taken into account when aiming at reducing a (redundant) phenotype
#    definition. Our workflow addresses this by extracting phenotypically 
#    important markers and binary cutoffs from a high-resolution setting and 
#    use these as input for flowType/Rchyoptimyx (Aghaeepour et al.), a
#    pipeline designed to find the simplest non-redundant phenotypes which are  
#    most important to an external outcome.


library(flowType)
library(RchyOptimyx)
# browseVignettes("flowType")

Thresholds <- cutoffs

# if we were to test all markers for 1 cutoff and 8 target markers
calcNumPops(rep(2,length(visne.ch.idx)), 8) 

# we test all 8 selected markers for 1 cutoff, i.e. 2 partitions
calcNumPops(rep(2, 8), 8) # feasible !

# set args for flowType
PropMarkers <- which( p.desc %in% names(Thresholds) ) 
PropMarkers <- PropMarkers[match(names(Thresholds), p.desc[PropMarkers])] # keep order
MFIMarkers <- NULL
MarkerNames <- names(Thresholds)
PartitionsPerMarker <- sapply(Thresholds, length) + 1

# run flowType
ResList <- fsApply(fs, function(fr) flowType(fr, PropMarkers, 
                   MFIMarkers, 'thresholds', MarkerNames, 
                   MaxMarkersPerPop = NULL, PartitionsPerMarker, 
                   Thresholds, MemLimit = 6))
                   
# extract all cell proportions from the list of flowType results 
# and normalize them by the total number of cells
PropsList <- sapply(ResList,slot,"CellFreqs")
freq.ft <- sweep(PropsList, 2, PropsList[1, ], '/')

# get the corresponding subset phenotype names
dimnames(freq.ft)[[1]] <-  unlist(lapply(ResList[[1]]@PhenoCodes,
                      function(x){return(decodePhenotype(
                       x, ResList[[1]]@MarkerNames,
                       ResList[[1]]@PartitionsPerMarker))}))

# get flowType phenocodes
pheno.codes <- ResList[[1]]@PhenoCodes
pheno.names <- dimnames(freq.ft)[[1]]

# perform Skilling-Mack w/o MC-simulation
sm.Rchy.pvals <-  sm.test(y=t(freq.ft)[rmvd, ],
						    groups=study$Region[rmvd], blocks=study$Donor[rmvd], 
								simulate.p.value = F, B = 10000, sim.SM=FALSE)$p.value
sm.Rchy.ranked <- order(sm.Rchy.pvals , decreasing=FALSE)
best.phenotype.index <- sm.Rchy.ranked[1]
selected.phenotype <- pheno.codes[best.phenotype.index]

# top 20 phenotypes
pheno.names[sm.Rchy.ranked[1:20]]

# target subset phenotypes
target.phenos <- targets[grepl("f·", names(targets))]
(target.phenos <- unlist(sapply(target.phenos, paste, collapse="")))
# or specify manually
target.phenos <- c("CD11c+CD206-CD45+CD64+CD68+CX3CR1+HLADR+IRF8+",
                   "CD11c-CD206+CD45-CD64+CD68-CX3CR1-HLADR-IRF8+",
                   # add CD206+ to subset 3 as ~50% were positive
                   "CD11c-CD206+CD45-CD64-CD68-CX3CR1-HLADR-IRF8+", 
                   "CD11c-CD45-CD64-CD68-CX3CR1-HLADR-IRF8+",
                   "CD11c-CD206-CD45-CD64-CD68-CX3CR1-HLADR-IRF8-")

# get the corresponding ids
subset.ids <- which(pheno.names %in% target.phenos)
pheno.names[subset.ids]
subset.ids <- subset.ids[match(target.phenos, pheno.names[subset.ids])] # keep order
sm.Rchy.pvals[subset.ids]

# run RchyOptimyx
# hierarchy for best phenotype
bestRchyResult <- RchyOptimyx(pheno.codes=pheno.codes, 
                              phenotypeScores=-log2(sm.Rchy.pvals), 
                              startPhenotype=selected.phenotype , 
                              pathCount=3, trimPaths=FALSE)

# combine single paths into same graph
path.counts <- c(1, 1, 1, 1, 1)
for (i in seq_along(subset.ids)) { 
  tempRes <-  RchyOptimyx(pheno.codes=pheno.codes,
                          phenotypeScores=-log2(sm.Rchy.pvals), 
                          startPhenotype=pheno.codes[subset.ids[i]], 
                          pathCount=path.counts[i], trimPaths= F)
# merge in each hierarchy
  if (i != 1) {
    RchyResult <- merge(RchyResult , tempRes)
  } else RchyResult <- tempRes 
}

# *********** #
#  Figure 7e  #
# *********** #

pdf("plots/figure_7_e.pdf", width=15, height=5)
  plot(RchyResult, phenotypeCodes=pheno.codes ,
       phenotypeScores=-log2(sm.Rchy.pvals) ,
       marker.names= ResList[[1]]@MarkerNames, 
       edgeWeights=F, edgeLabels=F, nodeLabels=T,
       partitions.per.marker=ResList[[1]]@PartitionsPerMarker,
       root.name="all cells", ylab="-log2(p-value)", plot.legend=T,
       colors=c("#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
dev.off()

# finally look also how best scoring phenotype and minimal phenotypes attach
all.phenos <- unlist(sapply(targets, paste, collapse=""))
all.phenos <- c(all.phenos,
                # add CD206+ to subset 3
                "CD11c-CD206+CD45-CD64-CD68-CX3CR1-HLADR-IRF8+",
                # remove CD206 from subset 1 to see if still significant
                "CD11c+CD45+CD64+CD68+CX3CR1+HLADR+IRF8+")
all.ids <- which(pheno.names %in% all.phenos)

# iterate over selected fmeasured phenotypes
for (i in seq_along(all.ids)) { 
  tempRes <-  RchyOptimyx(pheno.codes=pheno.codes,
                          phenotypeScores=-log2(sm.Rchy.pvals), 
                          startPhenotype=pheno.codes[all.ids[i]], 
                          pathCount=3, trimPaths=F)	
  if (i != 1) {
    RchyResult <- merge(RchyResult , tempRes)
  } else RchyResult <- tempRes 
}
# merge in best scoring phenotype 
RchyResult <- merge(RchyResult, bestRchyResult)

pdf("plots/optimized_hierarchies.pdf", width=15, height=5)
plot(RchyResult, phenotypeCodes=pheno.codes ,
       phenotypeScores=-log2(sm.Rchy.pvals) ,
       marker.names= ResList[[1]]@MarkerNames, 
       edgeWeights=F, edgeLabels=F, nodeLabels=T,
       partitions.per.marker=ResList[[1]]@PartitionsPerMarker,
       root.name="all cells", ylab="-log2(p-value)", plot.legend=T,
       colors=c("#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
dev.off()

# Conclusions:
# 1) CD206 expression is critically important in brain-regional huMG 
#    heterogeneity (GFM).
# 2) Likewise, IRF8+ with coexpression (subset 1) or 'coabsence' (subset 3) 
#    of CD11c, CD45, CD64, and HLADR is important (SVZ) and irrespective  
#    of CD206. Detection of these two subsets may therefore be explained 
#    by compositional effects.


# clean-up  
rm(list=c("RchyResult", "ResList","freq.ft", "tempRes", "bestRchyResult",
           "pheno.codes", "pheno.names", "sm.Rchy.pvals", "sm.Rchy.ranked"))
detach("package:RchyOptimyx", unload=TRUE)
detach("package:flowType", unload=TRUE)
gc(); gcDLLs() 

# ******************************************************************************
#       testing for differential abundance using hyperspheres (Lun et al.)
# ******************************************************************************
library(cydar)
library(edgeR)
# browseVignettes("cydar")

# in this approach statistical inference is conducted at high-resolution in the
# original high-dimensional marker space, as opposed to testing for frequency   
# differences on a binned grid of a low-dimensional representation, i.e. tSNE map

# count data into hyperspheres, downsample=5, i.e. to 1 fifth of total, to reduce
# computational work
cd <- prepareCellData(fs, markers=colnames(fs)[visne.ch.idx], naive=FALSE)
cd <- countCells(cd, tol=0.5, downsample=5, filter=5, naive=FALSE)

# create edgeR object
y <- DGEList(assay(cd), lib.size=cd$totals)
y.keep <- aveLogCPM(y) >= aveLogCPM(5, mean(cd$totals))
sum(y.keep); length(y.keep)
cd <- cd[y.keep, ]
y <- y[y.keep, ]

# estimate the dispersions and fit a negative binomial generalized linear model,
design <- model.matrix(~ study$Region + study$Donor)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# perform quasi-liklihood test for donor to check whether there is a genuine 
# need to adjust for inter-subject variation 
res.donor <- glmQLFTest(fit, coef=6:13)
topTags(res.donor)

# test for regional differences
res <- glmQLFTest(fit, coef=2:5)
topTags(res)

pvals <- res$table$PValue

# correct for multiple testing by controlling spatial FDR
qvals <- spatialFDR(intensities(cd), pvals, naive=FALSE)
is.sig <- qvals <= 0.05
sum(is.sig)

# nonredundant hyperspheres
nonred <- findFirstSphere(intensities(cd), pvals)
sum(nonred & is.sig) 

# ******************************************************************************
# visualization of hyperspheres by tSNE dimensionality reduction

# get the center cells to map hyperspheres on the tSNE plot used above
sphere.center <- elementMetadata(cd)$center.cell
cd.loc <- t(cellIntensities(cd))
colnames(cd.loc) <- p.desc
# note that Lun et al. apply tSNE on median-positioned hyperspheres, i.e. 
# coords <- intensities(cd)

# -log scale fdr-adjusted p-values and threshold to 0.05
# such that significant weights are > 0 (used for color scale)
w.cd <- threshold2zero(w = -log2(qvals), th = -log2(0.05))

# ************************* #
#  Supplementary Figure 8h  #
# ************************* #

# ~26000 significant hyperspheres
png("plots/supplementary_figure_8_h.png", width=4.5, height=5, res=300, units="in")
  plot(cd.loc[sphere.center, c("tSNE1","tSNE2")], pch=19, 
       col=colorDens(w.cd), cex=0.15, xlim=xlim, ylim=ylim)
dev.off()

# same plot overlaid with ~150 nonredundant significant hyperspheres
# which would still be a lot to inspect/explore and to annotate individually
plot(cd.loc[sphere.center, c("tSNE1","tSNE2")], pch=19, 
     col=colorDens(w.cd), cex=0.15, xlim=xlim, ylim=ylim)
points(cd.loc[sphere.center[nonred & is.sig], c("tSNE1","tSNE2")], pch=19, 
       col=colorDens(w.cd)[nonred & is.sig], cex=1)


# ******************************************************************************
# phenotyping of contour-'gated' hyperspheres

# convert qvalues to weights to compute density difference
# such that significant weights are > 1 and sum(w) == length(w)
w.cd <- threshold2one(w = -log2(qvals), th = -log2(0.05))
sum(w.cd)
fo.cd <- kdeDiff(cd.loc[sphere.center, c("tSNE1","tSNE2")], 
                 binned=binned, bgridsize=bgridsize, H=Hopt, w1=w.cd)
                  
persp3d(fo.cd$eval.points[[1]], fo.cd$eval.points[[2]], fo.cd$estimate, 
        xlab="tSNE1", ylab="tSNE2", zlab="", aspect = c(1, 1, .25),
        alpha=1, col=colorDens(fo.cd$estimate), box = TRUE, axes = TRUE) 

dobs.cd <- predict(fo.cd, x=tsne.plot.data)

plot(tsne.plot.data, pch=19, cex=0.075, col=colorDens(dobs.cd), 
     xlim=xlim, ylim=ylim, xlab="tSNE1", ylab="tSNE2", main="")

# minBins specifies here minimum allowed number of hyperspheres
# according to frequency limit
min.hs <- 5 * freq.limit * nrow(cd.loc) / n.fr

# increase offset to gate highly significant hyperspheres or
# use more conservative q/p-value threshold for density difference
cs.cd <- contourLoc(fhat=fo.cd, offset=0.25, n.levels=50, 
					keep.level.by="max.bin.num",
					peaks2merge=6, minBins=min.hs, plotit=TRUE)

css.cd <- shrinkSpur(cs.cd, mul=T, minBins=min.hs)
plotContours(css.cd, comp.sel=T, col="grey", pch=20, cex=0.15, lwd=2.5,
             fhat.levels=F, labels="peak.ids", filwd=.5, density=50)

# compare and visualize expression levels 
# use hypersphere weighted median intensities
cd.hs <- intensities(cd)
cd.hs <- cd.loc[sphere.center,]
colnames(cd.hs) <- p.desc
# add tnse coords
cd.hs[, c("tSNE1", "tSNE2")] <- cd.loc[sphere.center, c("tSNE1", "tSNE2")]

# downsample expression data to enable comparison of delta values
set.seed(13)
csr.cd <- exprs2con(css.cd, data=cd.hs, markers=visne.ch.idx, 
                    dwnsmpl=128, delta=delta, beta=beta, threshold=threshold)


# we want to also see expression in the inner regions
csm.cd <- contourMerge(csr.cd, keep.outer=F, keep.inner=T, 
                       data=csr.cd$comb$data, delta=delta, beta=beta, minBins=min.hs)
plotContours(csm.cd, comp.sel=T, col="grey", pch=20, cex=0.15, lwd=2.5,
             fhat.levels=F, labels="peak.ids", filwd=.5, density=50)
plotPheno(csm.cd, comp.sel=T, relevant.only=TRUE)
plotPheno(csm.cd, comp.sel=F, relevant.only=TRUE)






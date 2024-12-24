##Nicholas Glass
##28 March 2020
##Bios 533

##Hello, this is the Project FINAL for Bios 533 for Nicholas Glass

##My project consists of two datasets: 1) Microbial count data from Trout Park, an impacted fen in northern IL, USA, and 2) Continuous soil nutrient accrual data from Midewin National Tallgrass Prairie in Joliet, IL.

#################################################################
##CLUSTER ANALYSIS##

# Load the required packages
library(ade4)
library(adespatial)
library(vegan)
library(gclus)
library(cluster)
library(pvclust)
library(RColorBrewer)
library(labdsv)
library(rioja)
library(indicspecies)
library(mvpart)
library(MVPARTwrap)
library(dendextend)
library(vegclust)
library(colorspace)
library(agricolae)
library(picante)
library(dplyr)

rm(list = ls())
##Set working directory to NEwR2-Functions
setwd("C:/Users/metal/Desktop/Final/NEwR2-Functions")

# Source additional functions that will be used later in this
# Chapter. Our scripts assume that files to be read are in
# the working directory.
source("drawmap.R")
source("drawmap3.R")
source("hcoplot.R")
source("test.a.R")
source("coldiss.R")
source("bartlett.perm.R")
source("boxplerk.R")
source("boxplert.R")

# Function to compute a binary dissimilarity matrix from clusters
grpdist <- function(X)
{
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}

##Set working directory to soil microbial sequencing data
setwd("C:/Users/metal/Desktop/Final/taxa_raw_counts")

# Load the data
DF <- read.delim("taxa_table_L2.txt",sep = '\t')
str(DF)

##Get rid of NAs
#DF <- slice(DF, 1, 3, 6, 7, 8, 14, 15, 18, 19)

##SCRIPT ASSUMES 3 FILES: 1) species (count data), 2) environment data (continuous), 3) spatial data

#env <- select(DF, pH, P, K, Ca, Mg, S, B, Cu, Fe, Mn, Zn, Na, Salts, OM, ENR, Cl)
spe <- DF%>%select(-1)
#spa <- select(DF, 1:6)

# Compute and plot dendrograms ====================================

## Hierarchical agglomerative clustering of the species abundance 
## data

# Compute matrix of chord distance among sites
spe.norm <- decostand(spe, "hellinger")
spe.ch <- vegdist(spe.norm, "bray")
## I used Hellinger transformation because it is a more conservative than relative abundance, more robust for constrained ordinations and produces less false positives.

# Attach site names to object of class 'dist'
attr(spe.ch, "Labels") <- DF$OTU_ID

# External graphical device to compare dendrograms
dev.new(
  title = "Compare clustering methods",
  width = 12,
  height = 8,
  noRStudioGD = TRUE
)
par(mfrow = c(2, 2))

# Compute single linkage agglomerative clustering
spe.ch.single <- hclust(spe.ch, method = "single")
# Plot a dendrogram using the default options
plot(spe.ch.single, 
     labels = rownames(spe), 
     main = "Chord - Single linkage")

# Compute complete-linkage agglomerative clustering
spe.ch.complete <- hclust(spe.ch, method = "complete")
plot(spe.ch.complete, 
     labels = rownames(spe), 
     main = "Chord - Complete linkage")

# Compute UPGMA agglomerative clustering
spe.ch.UPGMA <- hclust(spe.ch, method = "average")
plot(spe.ch.UPGMA, 
     labels = rownames(spe), 
     main = "Chord - UPGMA")

# Compute centroid clustering
spe.ch.centroid <- hclust(spe.ch, method = "centroid")
plot(spe.ch.centroid, 
     labels = rownames(spe), 
     main = "Chord - Centroid")

dev.new(
  title = "Compare clustering methods (2)",
  width = 12,
  height = 8,
  noRStudioGD = TRUE
)
par(mfrow = c(2, 2))

# Compute beta-flexible clustering using cluster::agnes()
# beta = -0.1
spe.ch.beta1 <- agnes(spe.ch, method = "flexible",
                      par.method = 0.55)
# beta = -0.25
spe.ch.beta2 <- agnes(spe.ch, method = "flexible",
                      par.method = 0.625)
# beta = -0.5
spe.ch.beta3 <- agnes(spe.ch, method = "flexible",
                      par.method = 0.75)
# Change the class of agnes objects
class(spe.ch.beta1)
spe.ch.beta1 <- as.hclust(spe.ch.beta1)
class(spe.ch.beta1)
spe.ch.beta2 <- as.hclust(spe.ch.beta2)
spe.ch.beta3 <- as.hclust(spe.ch.beta3)
plot(spe.ch.beta1, 
     labels = rownames(spe), 
     main = "Chord - Beta-flexible (beta=-0.1)")
plot(spe.ch.beta2, 
     labels = rownames(spe), 
     main = "Chord - Beta-flexible (beta=-0.25)")
plot(spe.ch.beta3, 
     labels = rownames(spe), 
     main = "Chord - Beta-flexible (beta=-0.5)")

# Compute Ward's minimum variance clustering
spe.ch.ward <- hclust(spe.ch, method = "ward.D2")
# Note: in R 3.0.3 or older, type:
#    spe.ch.ward <- hclust(spe.ch, method="ward")
#    spe.ch.ward$height <- sqrt(spe.ch.ward)
# The result will not be striclty the same, however.
# See the hclust help file and 'Numerical Ecology with R',
# 1st edition, p. 61-62.
plot(spe.ch.ward, 
     labels = rownames(spe), 
     main = "Chord - Ward")


# Cophenetic correlations =========================================

# Single linkage clustering
spe.ch.single.coph <- cophenetic(spe.ch.single)
cor(spe.ch, spe.ch.single.coph)
#0.89
# Complete linkage clustering
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.ch, spe.ch.comp.coph)
##0.92
# Average clustering
spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.ch, spe.ch.UPGMA.coph)
##0.98
# Ward clustering
spe.ch.ward.coph <- cophenetic(spe.ch.ward)
cor(spe.ch, spe.ch.ward.coph)
##0.70

# Shepard-like diagrams
dev.new(
  title = "Cophenetic correlation",
  width = 8,
  height = 9,
  noRStudioGD = TRUE
)
par(mfrow = c(2, 2))
plot(
  spe.ch,
  spe.ch.single.coph,
  xlab = "Chord distance",
  ylab = "Cophenetic distance",
  asp = 1,
  xlim = c(0, sqrt(2)),
  ylim = c(0, sqrt(2)),
  main = c("Single linkage", paste("Cophenetic correlation =",
                                   round(
                                     cor(spe.ch, spe.ch.single.coph), 3
                                   )))
)
abline(0, 1)
lines(lowess(spe.ch, spe.ch.single.coph), col = "red")
plot(
  spe.ch,
  spe.ch.comp.coph,
  xlab = "Chord distance",
  ylab = "Cophenetic distance",
  asp = 1,
  xlim = c(0, sqrt(2)),
  ylim = c(0, sqrt(2)),
  main = c("Complete linkage", paste("Cophenetic correlation =",
                                     round(
                                       cor(spe.ch, spe.ch.comp.coph), 3
                                     )))
)
abline(0, 1)
lines(lowess(spe.ch, spe.ch.comp.coph), col = "red")
plot(
  spe.ch,
  spe.ch.UPGMA.coph,
  xlab = "Chord distance",
  ylab = "Cophenetic distance",
  asp = 1,
  xlim = c(0, sqrt(2)),
  ylim = c(0, sqrt(2)),
  main = c("UPGMA", paste("Cophenetic correlation =",
                          round(
                            cor(spe.ch, spe.ch.UPGMA.coph), 3
                          )))
)
abline(0, 1)
lines(lowess(spe.ch, spe.ch.UPGMA.coph), col = "red")
plot(
  spe.ch,
  spe.ch.ward.coph,
  xlab = "Chord distance",
  ylab = "Cophenetic distance",
  asp = 1,
  xlim = c(0, sqrt(2)),
  ylim = c(0, max(spe.ch.ward$height)),
  main = c("Ward", paste("Cophenetic correlation =",
                         round(
                           cor(spe.ch, spe.ch.ward.coph), 3
                         )))
)
abline(0, 1)
lines(lowess(spe.ch, spe.ch.ward.coph), col = "red")

# Gower (1983) distance
(gow.dist.single <- sum((spe.ch - spe.ch.single.coph) ^ 2))
(gow.dist.comp <- sum((spe.ch - spe.ch.comp.coph) ^ 2))
(gow.dist.UPGMA <- sum((spe.ch - spe.ch.UPGMA.coph) ^ 2))
(gow.dist.ward <- sum((spe.ch - spe.ch.ward.coph) ^ 2))

##Lowest Gowers Distance is UPGMA at 9.0. Highest is Ward at 9883.5.

# Graphs of fusion level values ===================================

dev.new(
  title = "Fusion levels",
  width = 12,
  height = 8,
  noRStudioGD = TRUE
)
par(mfrow = c(2, 2))
# Plot the fusion level values of the single linkage clustering
# plot(spe.ch.single$height, nrow(spe):2, type="S",
# 	main="Fusion levels - Chord - Single",
# 	ylab="k (number of clusters)", 
#           xlab="h (node height)",
#           col="grey")
# text(spe.ch.single$height, nrow(spe):2, nrow(spe):2, 
#      col="red", cex=0.8)
# Plot the fusion level values of the complete linkage clustering
plot(
  spe.ch.complete$height,
  nrow(spe):2,
  type = "S",
  main = "Fusion levels - Chord - Complete",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(spe.ch.complete$height,
     nrow(spe):2,
     nrow(spe):2,
     col = "red",
     cex = 0.8)
# Plot the fusion level values of the UPGMA clustering
plot(
  spe.ch.UPGMA$height,
  nrow(spe):2,
  type = "S",
  main = "Fusion levels - Chord - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(spe.ch.UPGMA$height,
     nrow(spe):2,
     nrow(spe):2,
     col = "red",
     cex = 0.8)
# Plot the fusion level values of the Ward clustering
plot(
  spe.ch.ward$height,
  nrow(spe):2,
  type = "S",
  main = "Fusion levels - Chord - Ward",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(spe.ch.ward$height,
     nrow(spe):2,
     nrow(spe):2,
     col = "red",
     cex = 0.8)
# Plot the fusion level values of the beta-flexible 
# clustering (-0.25)
plot(
  spe.ch.beta2$height,
  nrow(spe):2,
  type = "S",
  main = "Fusion levels - Chord - Beta-flexible",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(spe.ch.beta2$height,
     nrow(spe):2,
     nrow(spe):2,
     col = "red",
     cex = 0.8)

# Compare dendrograms =============================================

## Cut the trees to obtain k groups and compare the group contents
## using contingency tables

# Choose a common number of groups
k <- 2  # Number of groups where at least a small jump is present
# in all four graphs of fusion levels
# Cut the dendrograms
spech.single.g <- cutree(spe.ch.single, k = k)
spech.complete.g <- cutree(spe.ch.complete, k = k)
spech.UPGMA.g <- cutree(spe.ch.UPGMA, k = k)
spech.ward.g <- cutree(spe.ch.ward, k = k)
spech.beta.g <- cutree(spe.ch.beta2, k = k)

# Compare classifications by constructing contingency tables
# Single vs complete linkage
table(spech.single.g, spech.complete.g)
# Single linkage vs UPGMA
table(spech.single.g, spech.UPGMA.g)
# Single linkage vs Ward
table(spech.single.g, spech.ward.g)
# Complete linkage vs UPGMA
table(spech.complete.g, spech.UPGMA.g)
# Complete linkage vs Ward
table(spech.complete.g, spech.ward.g)
# UPGMA vs Ward
table(spech.UPGMA.g, spech.ward.g)
# beta-flexible vs Ward
table(spech.beta.g, spech.ward.g)

# Multiscale bootstrap resampling =================================

# Hierarchical clustering with p-values via multiscale bootstrap
# resampling

# Compute p-values for all clusters (edges) of the dendrogram
spech.pv <-
  pvclust(t(spe.norm),
          method.hclust = "average",
          method.dist = "euc",
          parallel=TRUE)

# Plot dendrogram with p-values
dev.new(
  title = "Fish - Chord - pvclust",
  width = 12,
  height = 8,
  noRStudioGD = TRUE
)
par(mfrow = c(1, 1))
plot(spech.pv)

## The p-values seem okay.

# Highlight clusters with high "au" p-values
pvrect(spech.pv, alpha = 0.95, pv = "au")
lines(spech.pv)
pvrect(spech.pv, alpha = 0.91, border = 4)

# Optimal number of clusters ======================================

## Select a dendrogram (Ward/chord) and apply three criteria
## to choose the optimal number of clusters

# Choose and rename the dendrogram ("hclust" object)
hc <- spe.ch.UPGMA

## I choose UPGMA because it has the lowest Gower's Distance and highest Cophenetic Correlation.

dev.new(
  title = "Optimal number of clusters",
  width = 12,
  height = 8,
  noRStudioGD = TRUE
)
par(mfrow = c(1, 2))

# Average silhouette widths (Rousseeuw quality index)
Si <- numeric(nrow(spe))
for (k in 2:(nrow(spe) - 1))
{
  sil <- silhouette(cutree(hc, k = k), spe.ch)
  Si[k] <- summary(sil)$avg.width
}
k.best <- which.max(Si)
plot(
  1:nrow(spe),
  Si,
  type = "h",
  main = "Silhouette-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Average silhouette width"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(Si),
       pch = 16,
       col = "red",
       cex = 1.5
)
# Optimal number of clusters according to matrix correlation 
# statistic (Pearson)
kt <- data.frame(k = 1:nrow(spe), r = 0)
for (i in 2:(nrow(spe) - 1)) 
{
  gr <- cutree(hc, i)
  distgr <- grpdist(gr)
  mt <- cor(spe.ch, distgr, method = "pearson")
  kt[i, 2] <- mt
}
k.best <- which.max(kt$r)
plot(
  kt$k,
  kt$r,
  type = "h",
  main = "Matrix correlation-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Pearson's correlation"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(kt$r),
       pch = 16,
       col = "red",
       cex = 1.5)

# Final dendrogram with the selected clusters =====================

# Choose the number of clusters
## I HAVE CHOSEN TWENTY CLUSTERS
k <- 20
# Silhouette plot of the final partition
##CHOOSE METHOD HERE (UPGMA)
spech.ward.g <- cutree(spe.ch.UPGMA, k = k)
sil <- silhouette(spech.ward.g, spe.ch)
rownames(sil) <- row.names(spe)
dev.new(title = "Silhouette plot - UPGMA - k=2", noRStudioGD = TRUE)
plot(
  sil,
  main = "Silhouette plot - Chord - UPGMA",
  cex.names = 0.8,
  col = 2:(k + 1),
  nmax = 100
)

# Reorder clusters
spe.chwo <- reorder.hclust(spe.ch.UPGMA, spe.ch)

# Plot reordered dendrogram with group labels
dev.new(
  title = "Final dendrogram",
  width = 8,
  height = 6,
  noRStudioGD = TRUE
)
plot(
  spe.chwo,
  hang = -1,
  xlab = "2 groups",
  sub = "",
  ylab = "Height",
  main = "Bray - UPGMA (reordered)",
  labels = cutree(spe.chwo, k = k)
)
rect.hclust(spe.chwo, k = k)

hcoplot(spe.ch.UPGMA, spe.ch, lab = rownames(spe), k = k)

# Miscellaneous graphical outputs =================================

# Convert the "hclust" object into a "dendrogram" object
dend <- as.dendrogram(spe.chwo)

# Plot the dendrogram with colored branches using the dendextend 
# syntax
dev.new(
  title = "Colored dendrogram",
  width = 8,
  height = 6,
  noRStudioGD = TRUE
)
dend %>% set("branches_k_color", k = k) %>% plot

# Heat map of the dissimilarity matrix ordered with the dendrogram
dev.new(title = "Heatmap - sites", noRStudioGD = TRUE)
heatmap(
  as.matrix(spe.ch),
  Rowv = dend,
  symm = TRUE,
  margin = c(3, 3)
)

#Conclusions: clustering shows the existence of two main groups: one cluster that encompasses many of the samples has similar microbial communities, while other samples are more varied. These largely correspond to the two halves of the wetland that was sampled. One half was more impacted by salt pollution and water table alterations, resulting in a highly variable microbial community. The other half is less affected, and soil samples are more similar.
#further analyses ties soil data to microbial count data to determine the influence of soil properties on microbial diversity.

#################################################################
## Soil Data DISTRIBUTIONS ##

# Load packages, functions and data ===============================
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)
library(dplyr)
library(ggplot2)

rm(list = ls())
##Set working directory to NEwR2-Functions
setwd("C:/Users/metal/Desktop/Final/NEwR2-Functions")

# Source additional functions that will be used later in this
# Chapter. Our scripts assume that files to be read are in
# the working directory.
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")

##Set working directory to home folder
setwd("C:/Users/metal/Desktop/Final")

# Load the data
DF <- read.delim("dat.l2.rc.txt")
str(DF)
#This is a combination of soil data and microbial count data for 9 samples. See README for more information on this file.

##Get rid of NAs
DF <- slice(DF, 1, 3, 6, 7, 8, 14, 15, 18, 19)

##SCRIPT ASSUMES 3 FILES: 1) species (count data), 2) environment data (continuous), 3) spatial data
env <- select(DF, pH, P, K, Ca, Mg, S, B, Cu, Fe, Mn, Zn, Na, Salts, OM, ENR, Cl)
spe <- select(DF, 36:106)
spa <- select(DF, 1:6)

## Check linearity of soil data
pH <- env$pH
P <- env$P
K <- env$K
Ca <- env$Ca
Mg <- env$Mg
S <- env$S
B <- env$B
Cu <- env$Cu
Fe <- env$Fe
Mn <- env$Mn
Zn <- env$Zn
Na <- env$Na
Salts <- env$Salts
OM <- env$OM
ENR <- env$ENR
Cl <- env$Cl

values <- append(pH, P)
values <- append(values, K)
values <- append(values, Ca)
values <- append(values, Mg)
values <- append(values, S)
values <- append(values, B)
values <- append(values, Cu)
values <- append(values, Fe)
values <- append(values, Mn)
values <- append(values, Zn)
values <- append(values, Na)
values <- append(values, Salts)
values <- append(values, OM)
values <- append(values, ENR)
values <- append(values, Cl)
values

names <- rep(seq(1,9,1),16)
values <- as.data.frame(values)
names <- as.data.frame(names)
values <- bind_cols(names, values)

ggplot(values, aes(x=names, y=values)) + geom_point()

## Check Normality
###NORMAL###
shapiro.test(env$Mg)
shapiro.test(env$S)
shapiro.test(env$Cu)
shapiro.test(env$Mn)
shapiro.test(env$Na)
shapiro.test(env$pH)
shapiro.test(env$ENR)
shapiro.test(env$OM)
shapiro.test(env$Salts)
shapiro.test(env$K)

###NOT NORMAL###
shapiro.test(env$P)
shapiro.test(env$Ca)
shapiro.test(env$B)
shapiro.test(env$Fe)
shapiro.test(env$Cl)

shapiro.test(values$values)
## Soil data together not normal ##

env12 <- t(env)
env12 <- as.data.frame(env12)
env12 <- slice(env12, -4) #outlier data, Ca may be unreliable anyway due to fen soil
hist(env12$V1, 
    col="blue", 
     las=1, 
     breaks=20)
hist(env12$V9, 
     col="blue", 
     las=1, 
     breaks=20)
hist(env12$V5, 
     col="blue", 
     las=1, 
     breaks=20)
##This makes Correspondance analysis invalid, which assumes unimodal distribution.

##################################################################
# Principal coordinate analysis (PCoA) ============================
## PCoA on soil data in Euclidean dissimilaritry matrix

env <- decostand(env, "standardize")
spe.bray <- vegdist(env, method = "euclidean")
spe.b.pcoa <- cmdscale(spe.bray, k = (nrow(env) - 1), eig = TRUE)

is.euclid(spe.bray)
scores(spe.b.pcoa)
spe.b.pcoa

# Plot of the sites
dev.new(
  title = "PCoA - Percentage difference",
  noRStudioGD = TRUE
)
ordiplot(scores(spe.b.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with Soil Data in Euclidean Dissimilarity Matrix")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# A posteriori projection of environmental variables
(spe.b.pcoa.env <- envfit(spe.b.pcoa, env))
# Plot significant variables with a user-selected colour
plot(spe.b.pcoa.env, p.max = 0.05, col = 4)

############################################################################
## NMDS Plots of Soil Data ##
spe.nmds <- metaMDS(env, distance = "euclidean")
spe.nmds
spe.nmds$stress
dev.new(title = "NMDS on Soil Data - Percentage difference",
        noRStudioGD = TRUE
)
plot(
  spe.nmds,
  type = "t",
  main = paste(
    "NMDS/Percentage difference - Stress =",
    round(spe.nmds$stress, 3)
  )
)

# Shepard plot and goodness of fit
dev.new(title = "NMDS - Shepard plot",
        width = 12,
        height = 6,
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
stressplot(spe.nmds, main = "Shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds, type = "t", main = "Goodness of fit")
points(spe.nmds, display = "sites", cex = gof * 300)


# Add colours from a clustering result FROM SPECIES DATA to an NMDS SOIL DATA plot
# Ward clustering of percentage difference dissimilarity matrix
# and extraction of four groups
spe.norm <- decostand(spe, "hellinger")
spe.bray <- vegdist(spe.norm, "bray")
spe.bray
spe.bray.ward <-
  hclust(spe.bray, "ward.D") # Here better than ward.D2 for 2 groups
spe.bw.groups <- cutree(spe.bray.ward, k = 2)
grp.lev <- levels(factor(spe.bw.groups))

# Combination with NMDS result
sit.sc <- scores(spe.nmds)
dev.new(title = "NMDS plot with cluster colors",
        noRStudioGD = TRUE
)
p <-
  ordiplot(sit.sc, type = "n", 
           main = "NMDS/% difference + clusters Ward/% difference")
for (i in 1:length(grp.lev))
{
  points(sit.sc[spe.bw.groups == i, ],
         pch = (14 + i),
         cex = 2,
         col = i + 1)
}
text(sit.sc, row.names(spe), pos = 4, cex = 0.7)
# Add the dendrogram
ordicluster(p, spe.bray.ward, col = "dark grey")

###############################################################
## Measuring Collinearity in Soil Data ##

# Load packages, functions and data ===============================
library(ade4)
library(adegraphics)
library(adespatial)
library(cocorresp)
library(vegan)
library(vegan3d)
library(ape)   # For PCoA with Lingoes correction (not in the book)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)
library(dplyr)
library(caret)
library(car)

##Set working directory to NEwR2-Functions
setwd("C:/Users/metal/Desktop/Final/NEwR2-Functions")
source("hcoplot.R")
source("triplot.rda.R")
source("plot.lda.R")
source("polyvars.R")
source("screestick.R")

##Set working directory to home folder
setwd("C:/Users/metal/Desktop/Final")

## Correlation Tests to find Redudant Variables, choosen based on PCoA results ##
## Non-normal Soil Variables: Cl, Ca, B, Fe, P

cor.test(x=env$Ca, y=env$Fe, method = "kendall")
cor.test(x=env$Ca, y=env$Fe, method = "spearman")
## Weak Correlation
cor.test(x=env$pH, y=env$Salts)
## Correlation not strong at -0.55
cor.test(x=env$Mn, y=env$Cl, method = "spearman")
cor.test(x=env$Mn, y=env$Cl, method = "kendall")
## Correlation weak at 0.5ish
cor.test(x=env$B, y=env$OM, method = "spearman")
## Correlation weak at 0.09
cor.test(x=env$B, y=env$Zn, method = "spearman")
cor.test(x=env$B, y=env$Zn, method = "kendall")
## Correlation not strong at 0.6ish
cor.test(x=env$Zn, y=env$K)
## Correlation weak at 0.6
cor.test(x=env$S, y=env$Mn)
## Correlation weak at 0.6
cor.test(x=env$Cl, y=env$S, method = "spearman")
## Correlation weak at 0.37
cor.test(x=env$Fe, y=env$Na, method = "spearman")
cor.test(x=env$Fe, y=env$Na, method = "kendall")
## Correlation weak at 0.6

## Salt content Correlations ##
cor.test(x=env$Cl, y=env$Na, method = "spearman")
## Correlation weak at 0.4
cor.test(x=env$Cl, y=env$Salts, method = "spearman")
cor.test(x=env$Cl, y=env$Salts, method = "kendall")
## Correlation weak at 0.5
cor.test(x=env$Na, y=env$Salts)
## Positive Correlation at 0.9

# Hellinger-transform the species dataset
spe.hel <- decostand(spe, "hellinger")

# Standardize environmental variables
env.norm <- decostand(env, "standardize")

## Forward selection of explanatory variables
env2 <- dplyr::select(env.norm, -Salts)
env3 <- dplyr::select(env.norm, B, Cl, Na, Ca, pH, P, Zn, Fe)
## ENV3: A list of soil variables I am interested in based on previous tests and apriori knowledge.
## I have removed Salts because of correlation with Sodium
spe.rda.all <- rda(spe.hel ~ ., data = env2)
spe.rda.all3 <- rda(spe.hel ~ ., data = env3)

# Forward selection using vegan's ordistep()
# This function allows the use of factors. 
mod0 <- rda(spe.hel ~ 1, data = env2)
step.forward <- 
  ordistep(mod0, 
           scope = formula(spe.rda.all), 
           direction = "forward", 
           permutations = how(nperm = 499)
  )
RsquareAdj(step.forward)
## Env2 forward: P and Na only, R2 is 0.32.
mod0 <- rda(spe.hel ~ 1, data = env3)
step.forward <- 
  ordistep(mod0, 
           scope = formula(spe.rda.all3), 
           direction = "forward", 
           permutations = how(nperm = 499)
  )
RsquareAdj(step.forward)
## Env3 forward: R2 is 0.42.

# Backward elimination using vegan's ordistep()
step.backward <-
  ordistep(spe.rda.all,
           permutations = how(nperm = 499))
RsquareAdj(step.backward)
## Env2 backward is Zn, R2 is 0.52.
step.backward <-
  ordistep(spe.rda.all3,
           permutations = how(nperm = 499))
RsquareAdj(step.backward)
## Env3 backward is Zn, R2 is 0.39.

## Selected soil variables: P, Na, Fe, Zn.

#####################################################################
## Distance-based Redundancy Analyses ##
#######################################################
# Hellinger-transform the species dataset
spe.norm <- decostand(spe, "hellinger")
spe.ch <- vegdist(spe.norm, "bray")

env <- dplyr::select(env, Na, P, Zn, Fe)
env <- decostand(env, "standardize")

Na <- env$Na
P <- env$P
Na <- as.factor(Na)
P <- as.factor(P)

bray.env.cap <- 
  capscale(spe.ch ~ Na*P, 
           data = env, 
           add = "lingoes", 
           comm = spe.ch)
anova(bray.env.cap, permutations = how(nperm = 999))
bray.env.cap
## P is 0.17.

# Plot with "wa" scores to see dispersion of sites around the 
# factor levels
dev.new(
  title = "db-RDA - wa - scaling 1", 
  noRStudioGD = TRUE
)
triplot.rda(bray.env.cap, site.sc = "wa", scaling = 1)

Zn <- env$Zn
Fe <- env$Fe
Zn <- as.factor(Zn)
Fe <- as.factor(Fe)

bray.env.cap <- 
  capscale(spe.ch ~ Zn*Fe, 
           data = env, 
           add = "lingoes", 
           comm = spe.ch)
anova(bray.env.cap, permutations = how(nperm = 999))
bray.env.cap
##P is 0.004.

# Plot with "wa" scores to see dispersion of sites around the 
# factor levels
dev.new(
  title = "db-RDA - wa - scaling 1", 
  noRStudioGD = TRUE
)
triplot.rda(bray.env.cap, site.sc = "wa", scaling = 1)

#Conclusions: despite salt pollution being primary difference between the two halves of the fen, the soil data that explains the microbial data most is redox-active metals. This indicates aerobic processes may drive differences in the microbial communities.

#######################################################################
## UNIVARIATE ANALYSIS ##
########################################################################
# Principal coordinate analysis (PCoA) ============================
## To grab drivers of microbial variation between clusters
# Load packages, functions and data ===============================
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)
library(dplyr)

rm(list = ls())
##Set working directory to NEwR2-Functions
setwd("C:/Users/metal/Desktop/Final/NEwR2-Functions")

# Source additional functions that will be used later in this
# Chapter. Our scripts assume that files to be read are in
# the working directory.
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")

##Set working directory to home folder
setwd("C:/Users/metal/Desktop/Final")

# Load the data

DF <- read.delim("dat.l2.rc.txt")
str(DF)

##Get rid of NAs
DF <- slice(DF, 1, 3, 6, 7, 8, 14, 15, 18, 19)

##SCRIPT ASSUMES 3 FILES: 1) species (count data), 2) environment data (continuous), 3) spatial data

env <- dplyr::select(DF, pH, P, K, Ca, Mg, S, B, Cu, Fe, Mn, Zn, Na, Salts, OM, ENR, Cl)
spe <- dplyr::select(DF, 36:106)
spa <- dplyr::select(DF, 1:6)

########################################
#STEP 1: use PCoA to grab to most influencial microbial phyla
## PCoA on soil data in Euclidean dissimilaritry matrix

spe <- decostand(spe, "hellinger")
spe.bray <- vegdist(spe, method = "bray")
spe.b.pcoa <- cmdscale(spe.bray, k = (nrow(spe) - 1), eig = TRUE)

is.euclid(spe.bray)
scores(spe.b.pcoa)
spe.b.pcoa

# Plot of the sites
dev.new(
  title = "PCoA on fish species - Percentage difference",
  noRStudioGD = TRUE
)
ordiplot(scores(spe.b.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with Soil Data in Euclidean Dissimilarity Matrix")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Add weighted average projection of species
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], spe)
#text(spe.wa, rownames(spe.wa), cex = 0.7, col = "red")
# A posteriori projection of environmental variables
(spe.b.pcoa.env <- envfit(spe.b.pcoa, spe))
# Plot significant variables with a user-selected colour
plot(spe.b.pcoa.env, p.max = 0.01, col = 4)
##I have a list of 13 bacteria phylum that are driving differences in microbial community at 99%.
## I have chosen 4 out of the 13 based on anaerobic capacity.

#############################################################################
## Generalized Linear Model ##

library(car)
library(effects)
library(nlme)
library(lme4)
library(lattice)
library(MASS)
library(ggplot2)
library(grid)
library(rmarkdown)
library(coda)
library(rjags)
library(R2jags)
library(mgcv)
library(coefplot)
library(dplyr)
library(vegan)

##Set up raw file
id <- c(1,2,3,4,5,6,7,8,9)
raw <- as.data.frame(id)
bait <- dplyr::select(spe, B.Chloroflexi, B.Fibrobacteres, B.Kiritimatiellaeota, B.Spirochaetes)
luck <- dplyr::select(DF, Site, Elevation, P, K, Ca, Mg, S, B, Cu, Fe, Mn, Zn, Na, Cl)
raw <- bind_cols(raw, luck)
raw <- bind_cols(raw, bait)
##################
## STANDARDIZE the covariates
raw$Na.std      <- as.numeric(scale(raw$Na))
raw$Fe.std  <- as.numeric(scale(raw$Fe))
raw$Zn.std <- as.numeric(scale(raw$Zn))

# Fit the model using MAXIMUM LIKLIHOOD ESTIMATION (MLE)
M1 <- glm(B.Chloroflexi ~ Na.std + Fe.std + Zn.std,
          family = "poisson",
          data = raw)

print(summary(M1), digits=3, signif.stars = FALSE)

#Assess overdispersion:
E1 <- resid(M1, type = "pearson")  
N  <- nrow(raw)
p  <- length(coef(M1))
Over.disp <- sum(E1^2) / (N - p)  
Over.disp 
## Overdispersion is .002, an all-time low. Underdispersion makes this an inappropriate model.

## Ben Bolker's test for overdispersion:
## Normality of residuals:
shapiro.test(E1)

F1 <- fitted(M1, type = "response")

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values",
     ylab = "Pearson residuals", cex.lab = 1.5)
abline(h=0, lty = 2)

plot(cooks.distance(M1), type = "h",
     xlab = "Observation", ylab = "Cook distance",
     cex.lab =  1.5)

#Demonstrating GLM visualization, just for fun
# Set two covariates = 0, plot predicted Y vs. other
# covariate plus 95% CI 
par(mfrow = c(1,3), mar = c(5,5,2,2))
plot(x=raw$Cl.std , y = raw$B.Chloroflexi,
     xlab = "Soil Chloride (ppm)",
     ylab = "Bacteria Chloroflexi Abundance",
     cex.lab = 2,
     pch = 16, type = "n")

range(raw$Cl.std)
MyData <- data.frame(Na.std = seq(-1.1, 4, length = 25),
                     Fe.std  = 0,  #mean value
                     Zn.std = 0)  #mean value
P1 <- predict(M1, newdata = MyData, type = "link", se = TRUE)
lines(x = MyData$Na.std, y = exp(P1$fit), lwd = 3)
lines(x = MyData$Na.std, y = exp(P1$fit + 2*P1$se.fit), lwd = 3, lty = 2)
lines(x = MyData$Na.std, y = exp(P1$fit - 2*P1$se.fit), lwd = 3, lty = 2)


plot(x=raw$Fe.std , y = raw$B.Chloroflexi,
     xlab = "Soil Iron (ppm)",
     ylab = "Bacteria Chloroflexi Abundance",
     cex.lab = 2,
     pch = 16, type = "n")

range(raw$Fe.std)
MyData <- data.frame(Na.std = 0,
                     Fe.std  = seq(-3, 1.5, length = 25),
                     Zn.std = 0)
P1 <- predict(M1, newdata = MyData, type = "link", se = TRUE)
lines(x = MyData$Fe.std, y = exp(P1$fit), lwd = 3)
lines(x = MyData$Fe.std, y = exp(P1$fit + 2*P1$se.fit), lwd = 3, lty = 2)
lines(x = MyData$Fe.std, y = exp(P1$fit - 2*P1$se.fit), lwd = 3, lty = 2)


plot(x=raw$Zn.std , y = raw$B.Chloroflexi,
     xlab = "Soil Zinc (ppm)",
     ylab = "Bacteria Chloroflexi Abundance",
     cex.lab = 2,
     pch = 16, type = "n")

range(raw$Zn.std)
MyData <- data.frame(Na.std = 0,
                     Fe.std  = 0,
                     Zn.std = seq(-4, 1, length = 25))
P1 <- predict(M1, newdata = MyData, type = "link", se = TRUE)
lines(x = MyData$Zn.std, y = exp(P1$fit), lwd = 3)
lines(x = MyData$Zn.std, y = exp(P1$fit + 2*P1$se.fit), lwd = 3, lty = 2)
lines(x = MyData$Zn.std, y = exp(P1$fit - 2*P1$se.fit), lwd = 3, lty = 2)

######################## End of Plotting Code  ##############
##############################################################
## Generalized Linear Mixed Effects Model

library(dplyr)
library(vegan)
library(lattice)
library(lme4)
library(ggplot2)
library(glmmTMB)
library(bbmle)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

###############
M1 <- glmer(B.Fibrobacteres ~ Na.std + Zn.std + Fe.std + (1|Site),
            data = raw, family = poisson)

#Check for overdispersion via techniqe of Zuur et al. (2013)
E1 <- resid(M1, type = "pearson")
N  <- nrow(raw)
p  <- length(fixef(M1)) + 1
Overdispersion <- sum(E1^2) / (N - p)
Overdispersion #underdispersed

## Checking the Negative Binomial GLMM ##
M2<- glmmadmb(B.Fibrobacteres ~ Na.std + Zn.std + Fe.std, 
               random =~ 1|Site, 
               family = "nbinom", data=raw)

#Check for overdispersion via techniqe of Zuur et al. (2013)
E1 <- resid(M2, type = "pearson")
N  <- nrow(raw)
p  <- length(fixef(M2)) + 1
Overdispersion <- sum(E1^2) / (N - p)
Overdispersion

# Check for overdispersion via function of Bolker (2008),
# which uses a chi-square test to yield a P value:
overdisp_fun(M2)

M2
summary(M2)
## Poisson seems a bit better

M2.admb <- M1
#Plotting residuals
#
E2 <- resid(M2.admb, type = "pearson")
F1 <- fitted(M2.admb, type ="response")
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, 
     y = E1, 
     xlab = "Fitted values", 
     ylab = "Pearson residuals", 
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(x = raw$Na.std, y = E1,
     xlab = "Sodium", 
     ylab = "Pearson residuals", 
     cex.lab = 1.5)
abline(h = 0, lty = 2)

plot(B.Fibrobacteres ~ Zn.std, data = raw, 
     xlab = "Zinc", 
     ylab = "Pearson residuals", 
     cex.lab = 1.5) 
abline(h = 0, lty = 2)

plot(B.Fibrobacteres ~ Fe.std, data = raw, 
     xlab = "Iron", 
     ylab = "Pearson residuals", 
     cex.lab = 1.5) 
abline(h = 0, lty = 2)
#all look decent

boxplot(B.Fibrobacteres ~ Site, data = raw, 
     xlab = "Site", 
     ylab = "Pearson residuals", 
     cex.lab = 1.5) 
abline(h = 0, lty = 2)

## Compare for other Phyla
salts <- glmer(B.Chloroflexi ~ Cl.std + Na.std +(1|Site),
               data = raw, family = poisson)

metal <- glmer(B.Chloroflexi ~ Fe.std + Zn.std +(1|Site),
               data = raw, family = poisson)

key <- glmer(B.Chloroflexi ~ Cl.std + Na.std + Fe.std + Zn.std + (1|Site),
              data = raw, family = poisson)

salts
metal
key
####
salts <- glmer(B.Kiritimatiellaeota ~ Cl.std + Na.std +(1|Site),
               data = raw, family = poisson)

metal <- glmer(B.Kiritimatiellaeota ~ Fe.std + Zn.std +(1|Site),
               data = raw, family = poisson)

key <- glmer(B.Kiritimatiellaeota ~ Cl.std + Na.std + Fe.std + Zn.std + (1|Site),
             data = raw, family = poisson)

salts
metal
key
####
salts <- glmer(B.Spirochaetes ~ Cl.std + Na.std +(1|Site),
               data = raw, family = poisson)

metal <- glmer(B.Spirochaetes ~ Fe.std + Zn.std +(1|Site),
               data = raw, family = poisson)

key <- glmer(B.Spirochaetes ~ Cl.std + Na.std + Fe.std + Zn.std + (1|Site),
             data = raw, family = poisson)

salts
metal
key

############################################################################
## Bayesian Poisson GLMM ##

library(dplyr)
library(vegan)
library(lattice)
library(lme4)
library(R2jags)
library(ggplot2)
library(glmmTMB)
library(bbmle)

###############
raw$Fe.std <- as.numeric(scale(raw$Fe))
raw$Zn.std      <- as.numeric(scale(raw$Zn))
raw$Na.std      <- as.numeric(scale(raw$Na))

library(runjags)
#NOTE: to use runjags, you must have the latest version of JAGS installed on your computer.

hope <- template.jags(B.Kiritimatiellaeota ~ Na.std + Zn.std + Fe.std, raw,
                      file = "JAGSmodel.txt",
                      write.data = TRUE, write.inits = TRUE,
                      precision.prior = " dhalfcauchy(sigma)",
                      effect.prior = "dnorm(0, 10^-6)", n.chains = 3,
                      precision.inits = c(0.01, 10), effect.inits = c(-1, 1), inits = NULL)

light <- run.jags(hope)
light<- extend.jags(light, sample = 10000)
light<- extend.jags(light, sample = 50000)

E <- residuals(light, output='mean')
shapiro.test(E)
## Looks like residuals
add.summary(light, vars = NA, mutate = NA, psrf.target = 1.05,
            normalise.mcmc = TRUE, modeest.opts = list(), confidence = c(0.95))

## Compare to AIC of GLMM
M2 <- glmer(B.Kiritimatiellaeota ~ S.std + Zn.std + Mg.std + (1|Site),
            data = raw, family = poisson)

AIC(M2)

##Plotting and Results
print(light)
summary(light)

plot(light, plot.type = c("histogram"))
plot(light, plot.type = c("trace"))

laugh <- summary(light)
joy <- laugh[, c(4, 5, 1, 3)]

library(coefplot)
library(arm)

coef.vect  <- joy[2:4, c(1)]
sd.vect <- joy[2:4, c(2)]
long.names <- c("Sulfur", "Zinc", "Magnesium")

coefplot(coef.vect, sd.vect, varnames = long.names, main = "Kiritimatiellaeota", xlim = c(-3, 3))

## Simulate the range of the moderating variable
int.params <- c("alpha", "beta1", "beta2", "beta3")

x2.sim <- seq(min(raw$B.Kiritimatiellaeota), max(raw$B.Kiritimatiellaeota), by = 1)

int.mcmc <- as.mcmc(light)
int.mcmc.mat <- as.matrix(int.mcmc)
int.mcmc.dat <- as.data.frame(int.mcmc.mat)

int.sim <- matrix(rep(NA, nrow(int.mcmc.dat)*length(x2.sim)), nrow = nrow(int.mcmc.dat))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- int.mcmc.dat$S.std_coefficient * x2.sim[i] ## CHOOSE HERE
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))

plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper)

old <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) + 
  geom_line(color = "blue", alpha = 0.8, size = 0.5) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(x = x2.sim, y = bayes.c.eff.lower), color = "blue", alpha = 0.8, size = 0.5) + 
  geom_line(aes(x = x2.sim, y = bayes.c.eff.upper), color = "blue", alpha = 0.8, size = 0.5) +
  xlab("Soil Sulfur") + 
  ylab("Conditional Kiritimatiellaeota abundance") + 
  xlim(1, 40) +
  theme_classic()
 
old

for(i in 1:length(x2.sim)){
  int.sim[, i] <- int.mcmc.dat$Zn.std_coefficient * x2.sim[i] ## CHOOSE HERE
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))

plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper)

abe <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) + 
  geom_line(color = "blue", alpha = 0.8, size = 0.5) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(x = x2.sim, y = bayes.c.eff.lower), color = "blue", alpha = 0.8, size = 0.5) + 
  geom_line(aes(x = x2.sim, y = bayes.c.eff.upper), color = "blue", alpha = 0.8, size = 0.5) +
  xlab("Soil Zinc") + 
  ylab("Conditional Kiritimatiellaeota abundance") + 
  xlim(1, 18) +
  ylim(-80, 0) +
  theme_classic() 
abe

for(i in 1:length(x2.sim)){
  int.sim[, i] <- int.mcmc.dat$Mg.std_coefficient * x2.sim[i] ## CHOOSE HERE
}

bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))

plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper)

lincoln <- ggplot(plot.dat, aes(x = x2.sim, y = bayes.c.eff.mean)) + 
  geom_line(color = "blue", alpha = 0.8, size = 0.5) +
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(x = x2.sim, y = bayes.c.eff.lower), color = "blue", alpha = 0.8, size = 0.5) + 
  geom_line(aes(x = x2.sim, y = bayes.c.eff.upper), color = "blue", alpha = 0.8, size = 0.5) +
  xlab("Soil Magnesium") + 
  ylab("Conditional Kiritimatiellaeota abundance") + 
  xlim(1, 80) +
  ylim(0, 80) +
  theme_classic() 
lincoln

MyVar <- c("S", "Zn", 
           "Mg")

Myxyplot(raw, MyVar, "B.Kiritimatiellaeota")

write.jagsfile(light, "B.Kiritimatiellaeota JAGS.txt")
##################
## Comparing Phyla

hope <- template.jags(B.Chloroflexi ~ S.std + Zn.std + Mg.std, raw,
                      file = "JAGSmodel.txt", family = "Poisson",
                      write.data = TRUE, write.inits = TRUE,
                      precision.prior = " dhalfcauchy(sigma)",
                      effect.prior = "dnorm(0, 10^-6)", n.chains = 3,
                      precision.inits = c(0.01, 10), effect.inits = c(-1, 1), inits = NULL)

light <- run.jags(hope)
light<- extend.jags(light, sample = 10000)
light<- extend.jags(light, sample = 50000)

residuals(light, output='mean')

add.summary(light, vars = NA, mutate = NA, psrf.target = 1.05,
            normalise.mcmc = TRUE, modeest.opts = list(), confidence = c(0.95))

## Compare to AIC of GLMM
M2 <- glmer(B.Chloroflexi ~ S.std + Zn.std + Mg.std + (1|Site),
            data = raw, family = poisson)

AIC(M2)

print(light)
summary(light)
plot(light, plot.type = c("trace"))

laugh <- summary(light)
joy <- laugh[, c(4, 5, 1, 3)]

coef.vect  <- joy[2:4, c(1)]
sd.vect <- joy[2:4, c(2)]
long.names <- c("Sulfur", "Zinc", "Magnesium")

coefplot(coef.vect, sd.vect, varnames = long.names, main = "Chloroflexi", xlim = c(-1, 1))

write.jagsfile(light, "B.Chloroflexi JAGS.txt")
###

hope <- template.jags(B.Fibrobacteres ~ S.std + Zn.std + Mg.std, raw,
                      file = "JAGSmodel.txt", family = "Poisson",
                      write.data = TRUE, write.inits = TRUE,
                      precision.prior = " dhalfcauchy(sigma)",
                      effect.prior = "dnorm(0, 10^-6)", n.chains = 3,
                      precision.inits = c(0.01, 10), effect.inits = c(-1, 1), inits = NULL)

light <- run.jags(hope)
light<- extend.jags(light, sample = 10000)
light<- extend.jags(light, sample = 50000)

residuals(light, output='mean')

add.summary(light, vars = NA, mutate = NA, psrf.target = 1.05,
            normalise.mcmc = TRUE, modeest.opts = list(), confidence = c(0.95))

## Compare to AIC of GLMM
M2 <- glmer(B.Fibrobacteres ~ S.std + Zn.std + Mg.std + (1|Site),
            data = raw, family = poisson)

AIC(M2)

print(light)
summary(light)
plot(light, plot.type = c("trace"))

laugh <- summary(light)
joy <- laugh[, c(4, 5, 1, 3)]

coef.vect  <- joy[2:4, c(1)]
sd.vect <- joy[2:4, c(2)]
long.names <- c("Sulfur", "Zinc", "Magnesium")

coefplot(coef.vect, sd.vect, varnames = long.names, main = "Fibrobacteres", xlim = c(-2, 2))

write.jagsfile(light, "B.Fibrobacteres JAGS.txt")
###
hope <- template.jags(B.Spirochaetes ~ S.std + Zn.std + Mg.std, raw,
                      file = "JAGSmodel.txt", family = "Poisson",
                      write.data = TRUE, write.inits = TRUE,
                      precision.prior = " dhalfcauchy(sigma)",
                      effect.prior = "dnorm(0, 10^-6)", n.chains = 3,
                      precision.inits = c(0.01, 10), effect.inits = c(-1, 1), inits = NULL)

light <- run.jags(hope)
light<- extend.jags(light, sample = 10000)
light<- extend.jags(light, sample = 50000)

residuals(light, output='mean')

add.summary(light, vars = NA, mutate = NA, psrf.target = 1.05,
            normalise.mcmc = TRUE, modeest.opts = list(), confidence = c(0.95))

## Compare to AIC of GLMM
M2 <- glmer(B.Spirochaetes ~ S.std + Zn.std + Mg.std + (1|Site),
            data = raw, family = poisson)

AIC(M2)

print(light)
summary(light)
plot(light, plot.type = c("trace"))

laugh <- summary(light)
joy <- laugh[, c(4, 5, 1, 3)]

coef.vect  <- joy[2:4, c(1)]
sd.vect <- joy[2:4, c(2)]
long.names <- c("Sulfur", "Zinc", "Magnesium")

coefplot(coef.vect, sd.vect, varnames = long.names, main = "Spirochaetes", xlim = c(-5, 1))

write.jagsfile(light, "B.Spirochaetes JAGS.txt")

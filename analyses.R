###################################################################################
#
# Alex Tyrpak Data - Vacant Lot & Urban Meadow
#
# Chapter 1 in Thesis
#
# Reanalyzing the data based on methods in chapter
#
# KI Perry; 13 July 2021; checked and updated 11 December 2024
#
###################################################################################


# install and load packages needed for NMDS and beta-diversity analyses

if (!suppressWarnings(require(vegan))) install.packages("vegan")
#citation("vegan")

#install.packages("devtools")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#citation("pairwiseAdonis")

if (!suppressWarnings(require(BiodiversityR))) install.packages("BiodiversityR")
#citation("BiodiversityR")

if (!suppressWarnings(require(ggplot2))) install.packages("ggplot2")
#citation("ggplot2")

if (!suppressWarnings(require(reshape2))) install.packages("reshape2")
#citation("reshape2")

if (!suppressWarnings(require(viridis))) install.packages("viridis")
#citation("viridis")

# set up color and point vectors for figures

colvec <- c("gray60", "black")
pchvec1 <- c(16, 17)
pchvec2 <- c(19, 15)
ltyvec <- c(1, 2)

a15 <- read.csv("./data/ants_2015.csv", row.names = 1)
a16 <- read.csv("./data/ants_2016.csv", row.names = 1)


## Compare ant species richness across treatments

# individual-based rarefaction by treatment, jackknife estimates by treatment
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estimates are based on the number of singletons and doubletons

a15$trmt <- as.factor(a15$trm)
a16$trmt <- as.factor(a16$trm)

levels(a15$trmt)
rare.a15.VL <- a15[which(a15$trmt == "VL"),]
rare.a15.UM <- a15[which(a15$trmt == "UM"),]
rare.a16.VL <- a16[which(a16$trmt == "VL"),]
rare.a16.UM <- a16[which(a16$trmt == "UM"),]


sp.a15.VL <- specaccum(rare.a15.VL[6:22], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a15.UM <- specaccum(rare.a15.UM[6:22], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a16.VL <- specaccum(rare.a16.VL[6:23], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.a16.UM <- specaccum(rare.a16.UM[6:23], method = "rarefaction", permutations = 100, gamma = "jack2")


#2015
plot(sp.a15.VL, pch = 19, col = "black", xvar = c("individuals"), lty = 1, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 1800), ylim = c(0, 18))
plot(sp.a15.UM, add = TRUE, pch = 15, xvar = c("individuals"), lty = 2, lwd = 3, col = "gray60")
legend("bottomright", legend = c("Urban Meadow","Vacant Lot"), lty = ltyvec, bty = "n", lwd = 3, col = colvec)

#2016
plot(sp.a16.VL, pch = 19, col = "black", xvar = c("individuals"), lty = 1, lwd = 3,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 4000), ylim = c(0, 18))
plot(sp.a16.UM, add = TRUE, pch = 15, xvar = c("individuals"), lty = 2, lwd = 3, col = "gray60")
legend("bottomright", legend = c("Urban Meadow","Vacant Lot"), lty = ltyvec, bty = "n", lwd = 3, col = colvec)


#calculates species richness for each site
rich15 <- specnumber(a15[6:22])
rich16  <- specnumber(a16[6:23])

#calculates species richness by treatment
specnumber(a15[6:22], groups = a15$trmt)
specnumber(a16[6:23], groups = a16$trmt)

a15.sp.j1 <- diversitycomp(a15[6:22], y = a15, factor1 = "trmt", index = "jack1")
a15.sp.j1
a15.sp.j2 <- diversitycomp(a15[6:22], y = a15, factor1 = "trmt", index = "jack2")
a15.sp.j2

a15.j1 <- diversityresult(a15[6:22], y=NULL, index = "jack1")
a15.j1
a15.j2 <- diversityresult(a15[6:22], y=NULL, index = "jack2")
a15.j2


a16.sp.j1 <- diversitycomp(a16[6:23], y = a16, factor1 = "trmt", index = "jack1")
a16.sp.j1
a16.sp.j2 <- diversitycomp(a16[6:23], y = a16, factor1 = "trmt", index = "jack2")
a16.sp.j2

a16.j1 <- diversityresult(a16[6:23], y=NULL, index = "jack1")
a16.j1
a16.j2 <- diversityresult(a16[6:23], y=NULL, index = "jack2")
a16.j2

## Make rarefaction figure
png("Ants_Rarefaction.png", width = 2000, height = 800)
layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
layout.show(2)

par(mar=c(6,7,3,1))
plot(sp.a15.VL, pch = 19, col = "black", xvar = c("individuals"), lty = 1, lwd = 4,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 1800), ylim = c(0, 18),
     cex.lab = 2, cex.axis = 1.5)
plot(sp.a15.UM, add = TRUE, pch = 15, xvar = c("individuals"), lty = 2, lwd = 3, col = "gray60")
legend("bottomleft", legend = c("Urban Meadow","Vacant Lot"), lty = ltyvec, bty = "n", lwd = 5, col = colvec, cex = 2.3)
text(1700, 17.5, "A", pos = 4, font = 2, cex = 3)

par(mar=c(6,7,3,1))
plot(sp.a16.VL, pch = 19, col = "black", xvar = c("individuals"), lty = 1, lwd = 4,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 4000), ylim = c(0, 18),
     cex.lab = 2, cex.axis = 1.5)
plot(sp.a16.UM, add = TRUE, pch = 15, xvar = c("individuals"), lty = 2, lwd = 3, col = "gray60")
text(3800 ,17.7, "B", pos = 4, font = 2, cex = 3)
dev.off()

##############################################################################################################

# Compare ant species composition among treatments
# NMDS, PERMANOVA, and BETADISPER

a15.2 <- read.csv("./data/ant_assemblages_2015.csv", row.names = 1)
a16.2 <- read.csv("./data/ant_assemblages_2016.csv", row.names = 1)
str(a15.2)
str(a16.2)

ant.nmds15 <- metaMDS(a15.2[2:18], distance = "jaccard", trymax = 500, autotransform = TRUE)
ant.nmds15
stressplot(ant.nmds15)
ant.nmds15$stress
ant.nmds15$points
ant.nmds15$species
plot(ant.nmds15)

# PERMANOVA tests whether the group centroid of beetle communities on tree species
# differs in multivariate space (e.g. different community composition)
a15.dist <- vegdist(a15.2[2:18], method = "jaccard")
adonis2(a15.dist ~ a15.2$trmt, permutations = 999)

# BETADISPER tests whether the dispersion of a treatment from its spatial median is different
# between groups
# analysis of multivariate homogeneity of group dispersion (variances)
# multivariate analogue of Levene's test for homogeneity of variances
a15.beta <- betadisper(a15.dist, a15.2$trmt, type = c("median"))
a15.beta
anova(a15.beta)
plot(a15.beta)
boxplot(a15.beta, ylab = "Distance to median")
TukeyHSD(a15.beta, which = "group", conf.level = 0.95)

ant.nmds16 <- metaMDS(a16.2[2:19], distance = "jaccard", trymax = 500, autotransform = TRUE)
ant.nmds16
stressplot(ant.nmds16)
ant.nmds16$stress
ant.nmds16$points
ant.nmds16$species
plot(ant.nmds16)

#PERMANOVA
a16.dist <- vegdist(a16.2[2:19], method = "jaccard")
adonis2(a16.dist ~ a16.2$trmt, permutations = 999)

#BETADISPER
a16.beta <- betadisper(a16.dist, a16.2$trmt, type = c("median"))
a16.beta
anova(a16.beta)
plot(a16.beta)
boxplot(a16.beta, ylab = "Distance to median")
TukeyHSD(a16.beta, which = "group", conf.level = 0.95)


### NMDS Figures###
png("Ants_NMDS.png", width = 2000, height = 1000, pointsize = 30)

par(mfrow=c(1,2))
par(mar=c(5,4,3,2))

ordiplot(ant.nmds15, type="n", xlim = c(-1.2, 1.1), ylim = c(-0.1, 0.1))
points(ant.nmds15, dis = "sites", select = which(a15.2$trmt=="T1"), pch = 17, cex = 2, col = "black")
points(ant.nmds15, dis = "sites", select = which(a15.2$trmt=="T2"), pch = 16, cex = 2, col = "gray60")
ordiellipse(ant.nmds15, groups = a15.2$trmt,  col = c("black", "gray60"), display = "sites", draw = "lines", lwd = 5, conf = 0.90)
legend("bottomleft", legend = c("Urban Meadow","Vacant Lot"), pch = c(16, 17),
       cex = 1.2, bty = "n", col = c("gray60", "black"))
text(-1.25, 0.95, "A", pos = 4, font = 2, cex = 1.8)

ordiplot(ant.nmds16, type="n", xlim = c(-1.1, 1.1), ylim = c(-0.1, 0.1))
points(ant.nmds16, dis = "sites", select = which(a15.2$trmt=="T1"), pch = 17, cex = 2, col = "black")
points(ant.nmds16, dis = "sites", select = which(a15.2$trmt=="T2"), pch = 16, cex = 2, col = "gray60")
ordiellipse(ant.nmds16, groups = a16.2$trmt,  col = c("black", "gray60"), display = "sites", draw = "lines", lwd = 5, conf = 0.90)
text(-1.15, 0.9, "B", pos = 4, font = 2, cex = 1.8)

dev.off()

orditorp(ant.nmds15, "sites") #used to double check the legend!
orditorp(ant.nmds16, "sites") #used to double check the legend!

######################################################################################################

# Functional trait calculations

t <- read.csv("./data/ant_traits_pooled.csv", row.names = 1)
a <- read.csv("./data/ant_assemblages_pooled.csv", row.names=1)

if (!suppressWarnings(require(FD))) install.packages("FD")
citation("FD")

if (!suppressWarnings(require(gawdis))) install.packages("gawdis")
citation("gawdis")

rich <- rowSums(a)

# trim the trait dataset by identifying traits that are highly correlated
# or lack sufficient variance among species
names(t)
str(t)
plot(t)
cor(t, method = c("pearson"), use = "complete.obs")

#remove number of queens and the individual diet categories
t2 <- t[,-11] # number of queens
t2 <- t2[,-11] # belowground nests
t2 <- t2[,-11] # aboveground nests
t2 <- t2[,-11] # generalist predator
t2 <- t2[,-11] # generalist
t2 <- t2[,-11] # sugar feeder + generalist
t2 <- t2[,-11] # seed feeder + generalist
t2 <- t2[,-11] # specialist predator

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")

t2 <- t2[,-10] # native/exotic
t2 <- t2[,-4] # pronotum width

plot(t2)
cor(t2, method = c("pearson"), use = "complete.obs")
# will group head traits to limit their influence on functional diversity

# copy dataset in case transformations are needed
t3 <- t2

str(t3)
t3$ns <- as.factor(t3$ns)
t3$diet <- as.factor(t3$diet)
str(t3)

#check traits for normality
hist(t2$wl)
hist(log(t2$wl))#
t3$wl <- log(t3$wl + 1)

hist(t2$rhw)
hist(log(t2$rhw))
t3$rhw <- log(t3$rhw + 1)

hist(t2$rml)
hist(log(t2$rml))

hist(t2$rew)
hist(log(t2$rew))

hist(t2$rsl)
hist(log(t2$rsl))#
t3$rsl <- log(t3$rsl + 1)

hist(t2$rfl)
hist(log(t2$rfl))#
t3$rfl <- log(t3$rfl + 1)

hist(t2$rcl)
hist(log(t2$rcl))

#Double check all species present in both datasets
intersect(colnames(a), rownames(t3))
names(a)
#Double check if a species is present in one dataset but not the other
setdiff(colnames(a), rownames(t3))
setdiff(rownames (t3), colnames(a))

rownames(t3) == colnames(a) # we are good!

#diet 1 = predator
#diet 2 = generalist
#diet 3 = sugar feeder + generalist
#diet 4 = seed harvester + generalist

#nest 1 = hypogeic (UNDERGROUND)
#nest 2 = epigeic (ON SOIL SURFACE)


##Community-weighted mean for traits
cwm.obs <- functcomp(t3, as.matrix(a), CWM.type = "all")
cwm.obs

##Add species richness
cwm.obs$rich <- rich
cwm.obs

##Functional richness for traits
#create distance matrix with the traits
#optimized feature helps to weight the traits equally
tdis <- as.matrix(gawdis(t3, w.type = "optimized", opti.maxiter = 500))

frich <- dbFD(tdis, a, stand.x = TRUE)$FDis

##Add functional dispersion
cwm.obs$frich <- frich

plot(cwm.obs)
cor(cwm.obs, method = c("pearson"), use = "complete.obs")

######################################################################################################
## PLSCA

#install.packages("plsdepot")
library(plsdepot)
citation('plsdepot')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("mixOmics")

library(mixOmics)

env <- read.csv("./data/env_local_landscape.csv", row.names=1)

# landscape and local variables
plot(env[3:13], pch = 19)
cor(env[3:13], method = c("pearson"), use = "complete.obs")

# landscape variables only
plot(env[3:8], pch = 19)
cor(env[3:8], method = c("pearson"), use = "complete.obs")

dotchart(env$PLAND, pch = 19)
dotchart(env$PLAND, groups = env$trmt, pch = 19)
hist(env$PLAND)

dotchart(env$LPI, pch = 19)
dotchart(env$LPI, groups = env$trmt, pch = 19)
hist(env$LPI)

dotchart(env$ED, pch = 19)
dotchart(env$ED, groups = env$trmt, pch = 19)
hist(env$ED)

dotchart(env$ENN, pch = 19)
dotchart(env$ENN, groups = env$trmt, pch = 19)
hist(env$ENN)

dotchart(env$SHAPE, pch = 19)
dotchart(env$SHAPE, groups = env$trmt, pch = 19)
hist(env$SHAPE)

dotchart(env$IMP, pch = 19)
dotchart(env$IMP, groups = env$trmt, pch = 19)
hist(env$IMP)

env.land <- scale(env[3:8], center = TRUE, scale = TRUE)
cor(env.land, method = c("pearson"), use = "complete.obs")

# local variables only
plot(env[9:13], pch = 19)
cor(env[9:13], method = c("pearson"), use = "complete.obs")

dotchart(env$VEGBIOM, pch = 19)
dotchart(env$VEGBIOM, groups = env$trmt, pch = 19)
hist(env$VEGBIOM)

dotchart(env$VEGHT, pch = 19)
dotchart(env$VEGHT, groups = env$trmt, pch = 19)
hist(env$VEGHT)

dotchart(env$BLAB, pch = 19)
dotchart(env$BLAB, groups = env$trmt, pch = 19)
hist(env$BLAB)

dotchart(env$COMPACT, pch = 19)
dotchart(env$COMPACT, groups = env$trmt, pch = 19)
hist(env$COMPACT)

dotchart(env$PLI, pch = 19)
dotchart(env$PLI, groups = env$trmt, pch = 19)
hist(env$PLI)

env.local <- scale(env[9:13], center = TRUE, scale = TRUE)
cor(env.local, method = c("pearson"), use = "complete.obs")

#######################
## Landscape

land <- pls(env[3:8], cwm.obs, mode = c("canonical"), ncomp = 2, scale = TRUE, max.iter = 100)
land

land$explained_variance
land$loadings
land$variates

plotVar(land)
plotLoadings(land)

cim(land)$mat.cor

nw_land <- network(land, cutoff = 0.5, color.edge = color.spectral(2), lty.edge = c("solid", "dashed"),
                   lwd.edge = 2)
nw_land

land_nw <- as.data.frame(nw_land$M)
write.csv(land_nw, file = "Relevance_Network_Land.csv")

land2 <- plsca(env[3:8], cwm.obs, comps = 2, scaled=TRUE)
plot(land2)
land2

land2$cor.xt##these are correct for variables
land2$cor.yu##these are correct for taxa

land.var <- as.data.frame(land2$cor.xt)
ant.land.var <- as.data.frame(land2$cor.yu)

write.csv(land.var, file = "Correlation_Coefficients_Land.csv")
write.csv(ant.land.var, file = "Correlation_Coefficients_Land_Ants.csv")


land2$R2X##importance of variables on axes
land2$R2Y##importance of taxa on axes


#######################
## Local

local <- pls(env[9:13], cwm.obs, mode = c("canonical"), ncomp = 2, scale = TRUE, max.iter = 100)
local

local$explained_variance
local$loadings
local$variates

plotVar(local)
plotLoadings(local)

cim(local)$mat.cor

nw_local <- network(local, cutoff = 0.5, color.edge = color.spectral(2), lty.edge = c("solid", "dashed"),
                    lwd.edge = 2)
nw_local

local_nw <- as.data.frame(nw_local$M)
write.csv(local_nw, file = "Relevance_Network_Local.csv")

local2 <- plsca(env[9:13], cwm.obs, comps = 2, scaled=TRUE)
plot(local2)
local2

local2$cor.xt##these are correct for variables
local2$cor.yu##these are correct for taxa

local.var <- as.data.frame(local2$cor.xt)
ant.local.var <- as.data.frame(local2$cor.yu)

write.csv(local.var, file = "Correlation_Coefficients_Local.csv")
write.csv(ant.local.var, file = "Correlation_Coefficients_Local_Ants.csv")

local2$R2X##importance of variables on axes
local2$R2Y##importance of taxa on axes
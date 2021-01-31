# I need to import the Pop_genomics.RData first
load("/home/markoterzin/Documents/Antipathella/R_script_validation/Pop_genomics.RData")

###############################
# Stats and graph compilation #
###############################
# Tutorials: 
# 1. https://cran.r-project.org/web/packages/poppr/vignettes/mlg.html
# 2. This one seems to have everything I need: https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
# Combine 2 with this one: http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
# 4. Suggestions are here: https://juliehopper.wordpress.com/tag/poppr/

# Install the packages I need
# R will first look if the packages are installed, 
# and then later install them if not  

list.of.packages <- c("adegenet", "ape", "ggplot2","gtools", "devtools", "dplyr", "hierfstat",
                      "igraph", "pegas", "phangorn", "poppr", "viridis", "pheatmap", "mmod",
                      "RColorBrewer", "reshape2", "vcfR", "vegan", "fsthet", "lattice", "treemap", "magrittr", "diveRsity")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Installing adegenet was a bit tricky, done from the terminal using this tutorial in the end: 
# https://zoomadmin.com/HowToInstall/UbuntuPackage/r-cran-adegenet
library(adegenet)
library(ape)
library(ggplot2)
library(gtools)
library(igraph)
library(pegas)
library(phangorn)
library(poppr)
library(RColorBrewer)
library(reshape2)
library(vcfR)
library(vegan)
library(fsthet)
library(lattice)
library(dplyr)
library(treemap)
library(magrittr)
library(pheatmap)
library(viridis)
library(mmod)
library(hierfstat)
# Followed this to install diveRsity: 
# https://github.com/kkeenan02/diveRsity
library(diveRsity)

# Setting work directory
setwd("/home/markoterzin/Documents/Antipathella/R_script_validation/ALL_LOCI/")

###################################################################
# Now trying the analysis with all loci and outliers/neutral only #
###################################################################
# The second R script will be used for the 'ALL LOCI' dataset

###################################################
###**** Multi-locus genotype (MLG) analysis ****###
###################################################
# Using the naive algorithm
pdf("05A_MLGs_plot.pdf",
    width = 8, height = 6)
# Plot
data_genepop.tab <- mlg.table(data_final)
# Close the pdf file
dev.off()

###################################################
###**** Multi-locus genotype (MLG) analysis ****###
######******    with a 3% threshold    ******######
###################################################
# This tutorial: https://grunwaldlab.github.io/poppr/reference/mlg.filter.html
pc <- as.genclone(data_final, threads = 1L) # convert to genclone object
data_MLGs_thr_0.03 <- mlg.filter(pc, 
                                 threshold = 0.03, # 3 %
                                 distance = "nei.dist", # Nei distances
                                 threads = 1L,
                                 stats = 'ALL')

# Visualizing in ggplot2
library(ggplot2)

# But let's set the colors first before making a plot
cols <- c("tomato3", # BORD
          "slateblue4",# FAV
          "lightsteelblue4", # POS 
          "khaki3", # PTF
          "palegreen4", # SAV
          "yellow3") # SLucia

# Just to see how their data set looks like
# data(mpg, package="ggplot2")
# mpg <- read.csv("http://goo.gl/uEeRGu")

# Scatterplot
theme_set(theme_bw())  # pre-set the bw theme.
ind.namesMLG <- indNames(data_final)
ind.namesMLG <- as.character(ind.namesMLG) # These are my individuals
MLGs_0.03 <- data_MLGs_thr_0.03$MLGS # And the MLGs constructed based on 0.03 % threshold
popinfoMLG <- data_final@pop # This is the information on populations

# Combining into one data frame
MLG_plot <- cbind(ind.namesMLG, popinfoMLG, MLGs_0.03)
MLG_plot <- as.data.frame(MLG_plot)
MLG_plot

# Let's visualize!
pdf("05_MLG_0.03_threshold.pdf",
    width = 10, height = 6)
# Plot
g <- ggplot(MLG_plot, 
            aes(x=ind.namesMLG, 
                y=MLGs_0.03, 
                colour=data_final@pop))
# g + geom_count(col="lightblue4", show.legend=F) + # choosing the color of the points
g + geom_point(size=2, stat='identity') + 
  scale_color_manual(values = cols) +
  facet_grid(~data_final@pop, # placing the ind within their pops
             scales = "free_x", 
             space = "free") + 
  # facet_grid(cols = vars(Savaglia_final_neutral2@pop)) +
  labs(subtitle="based on a 3% similarity cut-off (Nei distances)", # putting the title / axis names
       y="MLGs",
       x="Individual names", 
       title="MLGs Plot") +
  theme(axis.text.x = element_text(angle = 90, # rotating the sample ID for 90 degrees 
                                   vjust = 0.5, 
                                   hjust=1)) 
# Close pdf
dev.off()
# Only two individuals were detected as clonal (from Savona) using a 3% dissimilarity threshold!

#########################
################################################
###**** Locus stats, heterozygosity, HWE ****###
################################################

# A rigorous population genetic analysis looks closely at the data to assess quality and identify 
# outliers or problems in the data such as erroneous allele calls. This chapter focuses on analysis 
# on a per-locus level. While there are statistics that analyze populations across loci, it is 
# important to analyze each locus independently to make sure that one locus is not introducing bias 
# or spurious errors into the analysis. Locus summary statistics: A quick way to assess quality of 
# the data is to determine the number, diversity, expected heterozygosity, and evenness of the 
# alleles at each locus.
locus_table(data_final)
## 
## allele = Number of observed alleles
# I have 2 alleles within each locus because I used the --write_single_snp
## 
## 1-D = Simpson index
## 
## Hexp = Nei's 1978 gene diversity

# We can also do this for each population
locus_table(data_final, pop = "BORD")
locus_table(data_final, pop = "FAV")
locus_table(data_final, pop = "SAV")
locus_table(data_final, pop = "PTF")
locus_table(data_final, pop = "POS")
locus_table(data_final, pop = "SLucia")

# But we can remove these monomorphic loci, and this is important because they are phylogenetically 
# uninformative!
nLoc(data_final)  # Let's look at our data set, note how many loci we have.
informative_data_final <- informloci(data_final)
nLoc(informative_data_final)

#########################
###########################################################
###**** Genotypic richness, diversity, and evenness ****###
###########################################################
# Tutorial here: https://grunwaldlab.github.io/Population_Genetics_in_R/Genotypic_EvenRichDiv.html

# Let's have a look at populations only
setPop(data_final) <- ~Population
data_final

# To calculate genotypic richness, diversity, and evenness, we can use the poppr function:
data_final_diversity <- poppr(data_final)
data_final_diversity

# Genotypic richness
# The number of observed MLGs is equivalent to genotypic richness. A type of a rarefaction method... 
# (See tutorial)

# Open pdf
pdf("08A_Genotypic_richness_rarefaction.pdf", 
width = 10, height = 10)
# Plot
data_final.tab <- mlg.table(data_final, plot = FALSE)
min_sample <- min(rowSums(data_final.tab))
rarecurve(data_final.tab, sample = min_sample, xlab = "Sample Size", ylab = "Expected MLGs")
title("Rarefaction of all my populations: I get a linear correlation because I don't have any clones")
# Close the pdf file
dev.off()
# I get a linear correlation because I don't have any clones!

# Genotypic diversity

# Diversity measures incorporate both genotypic richness and abundance. There are three measures of 
# genotypic diversity employed by poppr, the Shannon-Wiener index (H), 
# Stoddart and Taylor’s index (G), and Simpson’s index (lambda)
N  <- data_final_diversity$N  # number of samples
N
lambda <- data_final_diversity$lambda # Simpson's index
lambda
# This is repetitive, and was already calculated as lambda before
(N/(N - 1)) * lambda  # Corrected Simpson's index
# I get 1 for everything after correction...

# Genotypic evenness

# Evenness is a measure of the distribution of genotype abundances, where in a population with 
# equally abundant genotypes yields a value equal to 1 and a population dominated by a single 
# genotype is closer to zero. That is why I am getting 1!

# Open pdf
pdf("08B_Genotypic_eveness.pdf", 
width = 10, height = 10)
# Plot
data_final_tab <- mlg.table(data_final)
# Close pdf
dev.off()

#########################
######################################
###**** Linkage disequilibrium ****###
######################################
# Tutorial here: https://grunwaldlab.github.io/Population_Genetics_in_R/Linkage_disequilibrium.html

# Linkage disequilibrium test is useful to determine if populations are clonal (where significant disequilibrium
# is expected due to linkage among loci) or sexual (where linkage among loci is not expected). The null 
# hypothesis tested is that alleles observed at different loci are not linked if populations are sexual while 
# alleles recombine freely into new genotypes during the process of sexual reproduction. In molecular ecology 
# we typically use the index of association or related indices to test this phenomenon.
# We will analyze all the populations with the index of association and use 999 permutations of the data in 
# order to get a p-value. Note that the p-value is calculated with the original observation included.

# Open pdf
pdf("09_Linkage_disequilibrium_BORD.pdf", 
width = 5, height = 5)
# Plot
# Bordighera
BORDIGHERA <- popsub(data_final, "BORD")
ia(BORDIGHERA, sample = 999)
BORDIGHERA
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_FAV.pdf", 
width = 5, height = 5)
# Plot
# Favignana
FAVIGNANA <- popsub(data_final, "FAV")
ia(FAVIGNANA, sample = 999)
FAVIGNANA
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_SLucia.pdf", 
width = 5, height = 5)
# Plot
# Santa Lucia
Santa_Lucia <- popsub(data_final, "Slucia")
ia(Santa_Lucia, sample = 999)
Santa_Lucia
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_Savona.pdf", 
width = 5, height = 5)
# Plot
# Savona
SAVONA <- popsub(data_final, "SAV")
ia(SAVONA, sample = 999)
SAVONA
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_PTF.pdf", 
width = 5, height = 5)
# Plot
# Portofino
PORTOFINO <- popsub(data_final, "PTF")
ia(PORTOFINO, sample = 999)
PORTOFINO
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_POS.pdf", 
width = 5, height = 5)
# Plot
# Posada
POSADA <- popsub(data_final, "POS")
ia(POSADA, sample = 999)
POSADA
# Close pdf
dev.off()

# Try with clone correction now
##################################################
# Open pdf
pdf("09_Linkage_disequilibrium_BORD_cc.pdf", 
width = 5, height = 5)
# Plot
# Bordighera
BORDIGHERA %>% clonecorrect(strata= ~Depth/Population) %>% ia(sample = 999)
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_FAV_cc.pdf", 
    width = 5, height = 5)
# Plot
# Favignana
FAVIGNANA %>% clonecorrect(strata= ~Depth/Population) %>% ia(sample = 999)
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_POSADA_cc.pdf", 
    width = 5, height = 5)
# Plot
# Posada
POSADA %>% clonecorrect(strata= ~Depth/Population) %>% ia(sample = 999)
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_PTF_cc.pdf", 
    width = 5, height = 5)
# Plot
# Portofino
PORTOFINO %>% clonecorrect(strata= ~Depth/Population) %>% ia(sample = 999)
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_SAVONA_cc.pdf", 
    width = 5, height = 5)
# Plot
# Savona
SAVONA %>% clonecorrect(strata= ~Depth/Population) %>% ia(sample = 999)
# Close pdf
dev.off()

# Open pdf
pdf("09_Linkage_disequilibrium_SLucia_cc.pdf", 
    width = 5, height = 5)
# Plot
# Santa Lucia
Santa_Lucia %>% clonecorrect(strata= ~Depth/Population) %>% ia(sample = 999)
# Close pdf
dev.off()

# Pairwise r¯d over all loci: To ensure that the pattern of linkage disequilibrium seen is not due 
# to a single pair of loci, you can calculate IA and r¯d over all pairs of loci. This will be done 
# with the clone correction, although this shouldn't make a difference for me as I didn't identify 
# any clones

# Open pdf
BORDIGHERApair <- BORDIGHERA %>% clonecorrect(strata= ~Depth/Population) %>% pair.ia
FAVIGNANApair <- FAVIGNANA %>% clonecorrect(strata= ~Depth/Population) %>% pair.ia
POSADApair <- POSADA %>% clonecorrect(strata= ~Depth/Population) %>% pair.ia
PORTOFINOpair <- PORTOFINO %>% clonecorrect(strata= ~Depth/Population) %>% pair.ia
SAVONApair <- SAVONA %>% clonecorrect(strata= ~Depth/Population) %>% pair.ia
SANTA_LUCIApair <- Santa_Lucia %>% clonecorrect(strata= ~Depth/Population) %>% pair.ia

# In here I make sure that the values are all on the same scale
# I can only choose 2 populations (I cannot put all 6),
# So I chose Santa Lucia because it had the lowest rd val (-1)
# All the pops have the same highest rd val
plotrange <- range(c(BORDIGHERApair, 
                     SANTA_LUCIApair), 
                   na.rm = TRUE)
plotrange

pdf("09_All_Pops_and_Loci_Pairwise_Linkage_disequilibrium_cc.pdf", 
width = 200, height = 200)
# Plot
plot(BORDIGHERApair, limits = plotrange)
plot(FAVIGNANApair, limits = plotrange)
plot(POSADApair, limits = plotrange)
plot(PORTOFINOpair, limits = plotrange)
plot(SAVONApair, limits = plotrange)
plot(SANTA_LUCIApair, limits = plotrange)
# Close pdf
dev.off()

# Overall only the histograms are useful in my case, not the heatmaps
# But I kept the code just in case

#########################
#########################################
###**** Population structure: GST ****###
#########################################
# Tutorial here: https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html

# Gst

# Let's use Hendrick’s standardized GST to assess population structure among these populations 
# (Hedrick, 2005).

# This is where I used the mmod (Modern Methods of Differentiation) R package
Gst_Hedrick(data_final)	

#########################
#############################################################
###**** Discriminant Analysis of Principal Components ****###
#############################################################

# The DAPC is a multivariate statistical approach that uses populations defined a priori to maximize the variance 
# among populations in the sample by partitioning it into between-population and within-population components

# How can we determine the optimal number of PCs to retain?
#############################################################################
# Using the alpha score and cross validation to find an ideal number of PCs #
#############################################################################
# Use the number of PCs as suggested by the calculated scores

# LET'S START WITH ALPHA SCORE
# The trade-off between power of discrimination and over-fitting can be measured by the ascore, which is simply 
# the difference between the proportion of successful reassignment of the analysis (observed discrimination) and 
# values obtained using random groups (random discrimination). It can be seen as the proportion of successful 
# reassignment corrected for the number of retained PCs. It is implemented by a.score, which relies on repeating 
# the DAPC analysis using randomized groups, and computing a-scores for each group, as well as the average 
# a-score:

C.dapc_alpha <- dapc(data_final, n.da=100, n.pca=100)
temp_alpha <- a.score(C.dapc_alpha)
names(temp_alpha)
temp_alpha$tab[1:6,1:6]
temp_alpha$pop.score
temp_alpha$mean

# The number of retained PCs can be chosen so as to optimize the a-score; this is achived by optim.a.score:
# Open pdf
pdf("15_Optimizing_alpha_score.pdf", width = 10, height = 10)
# Plot
C.dapc_alpha <- dapc(data_final, n.da=100, n.pca=100)
temp_alpha <- optim.a.score(C.dapc_alpha)
# Close the pdf file
dev.off()

# We perform the analysis again with the number of PCs as suggested by alpha score, and then map the membership
# probabilities as before:

# Cross-validation

# Cross-validation (carried out with the function xvalDapc) provides an objective optimisation procedure for 
# identifying the ’golidlocks point’ in the trade-off between retaining too few and too many PCs in the model. 
# In cross-validation, the data is divided into two sets: a training set (typically comprising 90% of the data) 
# and a validation set (which contains the remainder (by default, 10%) of the data). With xvalDapc, the 
# validation set is selected by stratified random sampling: this ensures that at least one member of each group 
# or population in the original data is represented in both training and validation sets. DAPC is carried out on 
# the training set with variable numbers of PCs retained, and the degree to which the analysis is able to 
# accurately predict the group membership of excluded individuals (those in the validation set) is used to 
# identify the optimal number of PCs to retain. At each level of PC retention, the sampling and DAPC procedures 
# are repeated n.rep times.

# Open pdf
pdf("15_Cross_validation.pdf", 
    width = 10, height = 10)
# Plot
mat <- as.matrix(data_final)
grp <- pop(data_final)
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)
# Close the pdf file
dev.off()

xval[2:6]

##############################################################################
# Now I can make the DAPC using the number of PCs as suggested by alpha
# optimization and cross validation scores

# Open pdf
pdf("15_DAPC.pdf", 
width = 10, height = 5)
# Plot
C.dapc <- dapc(data_final, n.pca = 7, n.da = 5)
# To confirm that the DAPC is similar to the PCA we can plot the data in a scatter plot.
scatter(C.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
posi.pca = "topleft", cleg = 0.75)
# Close the pdf file
dev.off()

# And finally STRUCTURE-like PLOTS
# Open pdf
pdf("15_STRUCTURE-like_PLOT.pdf", width = 10, height = 4)
# Plot
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(C.dapc, col = cols, posi = 'top')
# Close the pdf file
dev.off()

# But we want to organise them within populations
# Open pdf
pdf("15_STRUCTURE-like_PLOT_PER_POPULATION.pdf", 
width = 10, height = 4)
# Plot
dapc.results <- as.data.frame(C.dapc$posterior)
dapc.results$pop <- pop(data_final)
dapc.results$indNames <- rownames(dapc.results)

# ggplot2 has specific requirements for the structure of the data frame format, as it requires each 
# observation in rows, and all different values of these observations in columns 
# (i.e., a long format data frame). To transform the data frame we use the function melt from the package 
# reshape2. melt reorganizes the data frame into the required data frame format, where each membership 
# probability observation for a given population is a row with the sample name, original population, 
# and assigned population as columns.

dapc.results <- melt(dapc.results)

# I get this: Using pop, indNames as id variables
# BUT
# Ignore the prompt for now. Then, we rename the columns into more familiar terms:

colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

# ggplot2 will plot the dapc.results data frame we reorganized using melt, using the samples on the X-axis 
# and membership probabilities on the Y-axis. The fill color will indicate the original population assignments. 
# Each facet represents the original population assignment for each sample:

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free_x", space = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
# Close the pdf file
dev.off()

# This bar plot shows us a more organized perspective of our data set by contrasting the population membership 
# probability assignments against their original populations.

# Open pdf
pdf("15_DAPC_no_centroids.pdf", width = 10, height = 5)
# Plot
scatter(C.dapc, col=cols, scree.da=T, # show eigenvalues
        scree.pca = T, # show how many PCs were retained
        legend = T, # show populations
        posi.leg = "topleft", # legend position
cell=1.5, cex=2, bg="white",cstar=0)
# Close the pdf file
dev.off()

#####################
#####################
###**** AMOVA ****###
#####################

# Doing AMOVA from this Tutorial: https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html
# This manual is good in case of doubts: https://cran.r-project.org/web/packages/poppr/poppr.pdf
# Check this tutorial as well: https://github.com/marinegenomicslab/workflow/wiki/AMOVA-in-R

# AMOVA stands for Analysis of MOlecular VAriance and is a method to detect population differentiation 
# utilizing molecular markers (Excoffier, Smouse & Quattro, 1992). This procedure was initially implemented for 
# DNA haplotypes, but applies to any marker system. The implementation of AMOVA in poppr requires two very basic 
# components: (1) A distance matrix derived from the data and (2) a separate table used to partition the data 
# into different stratifications.

# The distance matrix can be calculated using any distance as long as it is euclidean.
###########################
# IMPORTANT, from the poppr manual: https://cran.r-project.org/web/packages/poppr/poppr.pdf
# On  Euclidean  Distances:: With  the ade4 implementation  of  AMOVA  (utilized  by poppr), distances must be 
# Euclidean (due to the nature of the calculations).  Unfortunately, many genetic distance  measures  are  not  
# always  euclidean  and  must  be  corrected  for  before  being  analyzed. Poppr automates this with three 
# methods implemented in ade4: quasieuclid(), lingoes(), andcailliez(). The correction of these distances 
# should not adversely affect the outcome of the analysis.

# Let's try this
data.euclidean_final <- bitwise.dist(data_final, euclidean = TRUE)
###################################################################################
# Validated! The plots look the same if I refer myself to dist=C.euclidean_no_NAs #
#   and if I just plot without it, meaning that corrections mentioned above have  # 
#   been deployed. BUT STILL NEED TO USE within = FALSE, otherwise df are wrong!  #
###################################################################################
# We will do AMOVA in poppr package with the previously created C.genlight_no_NAs object. Strata has been created already

# In panmictic populations, we would expect to see most of the variance arise from within samples. If we see 
# that the most of the variance occurs among samples within populations or among populations, then there is 
# evidence that we have some sort of population structure. In the case of clonal organisms, this would help 
# support a hypothesis of clonal reproduction.

# Let’s invoke the AMOVA functions with and without clone correction:
# AMOVA
data_final_amova <- poppr.amova(data_final, # genind object
                                        ~Depth/Population, # hierarchy 
                                        within = F) # this is to use populations as lowest level

# without clone correction
# within: logical. When this is set to TRUE (Default), variance within individuals are calculated as 
# well. If this is set to FALSE, the lowest level of the hierarchy will be the sample level (in my 
# case populations).
# IMPORTANT: This example from the tutorial was performed with a data set of dominant (AFLP) markers, 
# but it can also be performed on codominant markers such as SNPs. These provide more information 
# because within sample (individual) variance is also assessed
data_final_amova

data_final_amovacc <- poppr.amova(data_final, # genind object
                                          ~Depth/Population, # hierarchy 
                                          within = F,
                                          clonecorrect = TRUE) # with clone correction
# with clone correction
data_final_amovacc

# You can export data in the form of a table with write.table

# Significance testing
set.seed(1999)
data_final_amova_signif   <- randtest(data_final_amova, nrepet = 999)
data_final_amova_cc_signif <- randtest(data_final_amovacc, nrepet = 999)
# This was done with no correction method though!

data_final_amova_signif 
data_final_amova_cc_signif

# Open pdf
pdf("16_AMOVA_significance.pdf", width = 10, height = 6)
# Plot
plot(data_final_amova_signif)
# Close the pdf file
dev.off()

pdf("16_AMOVA_significance_clone_correction.pdf", width = 10, height = 6)
# Plot
plot(data_final_amova_cc_signif)
# Close the pdf file
dev.off()
# By samples they mean populations!

#####################
#########################
###**** PERMANOVA ****###
#########################

# Looking at populations only
data.adonis_P <- adonis(data.euclidean_final ~ Population,
                                strata,
                                perm=200)
data.adonis_P

# But how can I see which populations differ significantly??
# Trying pairwise.adonis function now from here: https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
# Code by Pedro Martinez Arbizu
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

# similarity euclidean from vegdist and Benjamini Hochberg correction
data_P_pairwise_adonis_eucl_holm <- pairwise.adonis(x=data.adonis.outliers[,4:1241], 
                                                          # Selecting loci
                                                factors=data.adonis.outliers$popinfo, 
                                                # Saying I want pairwise comparisons of populations
                                                sim.function='vegdist',
                                                sim.method='euclidean',
                                                p.adjust.m='holm') # the p.value correction method, 
# one of the methods supported by p.adjust(); default is 'bonferroni'
data_P_pairwise_adonis_eucl_holm

# They used Holm correction for Euclidean distances in the original script
data_D_pairwise_adonis_eucl_holm <- pairwise.adonis(x=data.adonis.outliers[,4:1241], 
                                                            # Selecting loci
                                                            factors=data.adonis.outliers$depthinfo, 
                                                            # Saying I want pairwise comparisons of populations
                                                            sim.function='vegdist',
                                                            sim.method='euclidean',
                                                            p.adjust.m='holm') # the p.value correction method, 
# one of the methods supported by p.adjust(); default is 'bonferroni'
data_D_pairwise_adonis_eucl_holm

# Export now!
write.csv(data_P_pairwise_adonis_eucl_holm, 
          file = "data_P_pairwise_adonis.csv", 
          quote = F)
write.csv(data_D_pairwise_adonis_eucl_holm, 
          file = "data_D_pairwise_adonis.csv", 
          quote = F)

# Betadisper now!
# Different from what I do in adonis!
# This tests for pairwise differences in variance among populations

# Implements Marti Anderson's PERMDISP2 procedure for the analysis of multivariate homogeneity of 
# group dispersions (variances). betadisper is a multivariate analogue of Levene's test for 
# homogeneity of variances
data.mod_P <- with(strata, betadisper(data.euclidean_final, Population))
data.mod_P

anova(data.mod_P)
permutest(data.mod_P)
Tukey_P <- TukeyHSD(data.mod_P)
Tukey_P

pdf("16_Adonis_plots_Populations.pdf", 
width = 5, height = 5)
# Plot
plot(data.mod_P, col=cols)
boxplot(data.mod_P, col=cols)
# Close the pdf file
dev.off()

# Now Depth
data.mod_D <- with(strata, betadisper(data.euclidean_final, Depth))
data.mod_D

anova(data.mod_D)
permutest(data.mod_D)
Tukey_D <- TukeyHSD(data.mod_D)
Tukey_D

pdf("16_Adonis_plot_Depth.pdf", 
    width = 5, height = 5)
# Plot
plot(data.mod_D, col=cols)
boxplot(data.mod_D, col=cols)
# Close the pdf file
dev.off()

#####################
##############################
###**** Fst statistics ****###
##############################

# Try this tutorial: http://adegenet.r-forge.r-project.org/files/Barcelona2015/practical-MVAintro.1.0.pdf
# This one also seems great: https://cran.r-project.org/web/packages/hierfstat/vignettes/hierfstat.html
# Population structure is traditionally measured and tested using F statistics, in particular the Fst, which 
# measures population differentiation (as the proportion of allelic variance occuring between  groups). The  
# package hierfstat implements a wealth of F statistics and related tests, now designed to work natively 
# with genind objects.

# Let's look at the overall structuring first

# Just to have a quick look at the data, let's see the basic stats of our object first: 
# https://rdrr.io/cran/hierfstat/man/basic.stats.html
data_basic_stats <- basic.stats(data_final[,-1])
data_basic_stats
write.csv(as.matrix(data_basic_stats$perloc), file = "data_basic_stats.csv", quote = F)
# Overall stats on locus level!

# But Ho, Hs and Fis are also given per population! (and are the only vals reported in
# Carreras et al 2019)
# So export them and calculate the means in excel
write.csv(data_basic_stats$Ho, file = "Spreadsheets_04_data_Ho_perPOP.csv", quote = F)
write.csv(data_basic_stats$Fis, file = "Spreadsheets_04_data_Fis_perPOP.csv", quote = F)
write.csv(data_basic_stats$Hs, file = "Spreadsheets_04_data_Hs_perPOP.csv", quote = F)

# We first compute overall F statistics, and then use Goudet’s G statistics to test the existence 
# of population structure
data.fst <- fstat(data_final)
data.fst

# And if you want to look at Fst only
data.fst_only <- fstat(data_final, fstonly=TRUE)
data.fst_only

# Goudet's test
data.gtest.D <- gstat.randtest(data_final,
                                       nsim=999,
                                       method="global", # "global": tests for genetic 
                                       # structuring given 'pop'.
                                       pop=data_final$strata$Depth)

data.gtest.P <- gstat.randtest(data_final,
                                       nsim=999,
                                       method="global", # "global": tests for genetic 
                                       # structuring given 'pop'.
                                       pop=data_final$strata$Population)

# Let's make a plot for the G test
pdf("17_Fst_g_test.pdf", 
width = 5, height = 5)
# Plot
plot(data.gtest.D, main = "Effect of Depth")
plot(data.gtest.P, main = "Effect of Population")
# Close the pdf file
dev.off()

# Looks like there is some structuring!

############################################
# PAIRWISE COMPARISONS here
# Script from Carreras et al 2019 
############################################

# Here is the script: (data_h is the data in hierfstat format
data_h <- genind2hierfstat(data_final)
# WC_fst is the Fst matrix 
WC_fst <- genet.dist(data_h,method = "Nei87")
                     
# p-values for Fst
library(parallel)
mat.obs <- as.matrix(WC_fst) # To get a matrix
NBPERM <- 999 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(data_h[,1]),(data_h[,-1])),method = "Nei87"),mc.cores = 4)
                     
                     library(ade4)
                     allTests <- list()
                     for(i in 1:(nrow(mat.obs)-1)){
                       for(j in 2:nrow(mat.obs)){
                         allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
                       }
                     }
                     
                     pvals <- matrix(data=NA,nrow=6,ncol=6) # 6 because we have six localities
                     x=0
                     for (i in 1:(nrow(mat.obs)-1)){
                       for(j in 2:nrow(mat.obs)){
                         x=x+1
                         #     pvals[i,j] <- allTests[[x]][[6]]
                         pvals[i,j] <- allTests[[x]]$pvalue
                       }
                     }
                     
                     pvals[lower.tri(pvals)] <- NA # Putting NAs in the lower triangle
                     fst_pval <- mat.obs
                     fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)] # Putting p vals in the
                     # upper triangle
                     diag(fst_pval) <- NA # Putting NAs in the diagonal
                     write.csv(fst_pval,file="Fst-pvalues_NEI_999boot.csv",quote = F)

# The p vals need to be FDR corrected now! BY correction was used in Carreras et al

# This was done before                     
# p.adjust.M <- p.adjust.methods[p.adjust.methods != "fdr"]
# p.adjust.M # So these are all the methods
                     
# Here I am getting the non corrected p vals from 999 permutations
pvals_for_BY <- fst_pval[upper.tri(fst_pval)]
pvals_for_BY # checking the object to see how they were arranged

# Doing the p adjustment
p.adj <- sapply(p.adjust.M, simplify="matrix", function(meth) p.adjust(pvals_for_BY,meth))
p.adj # in the column none you can see the uncorrected ones, the order is still the same as 
# in pvals_for_BY!
                     
# Let's keep only the BY correction as in Carreras et al 2019
# Take the sixth column with all rows
A_fst_BY_pvals <- p.adj[,6,drop=FALSE] 
# Round to 3 digits
A_fst_BY_pvals <- round(A_fst_BY_pvals, digits = 3)
A_fst_BY_pvals

# Add these into the upper triangle of the matrix
fst_BYpval <- fst_pval
fst_BYpval[upper.tri(fst_BYpval)] <- A_fst_BY_pvals
fst_BYpval
# Putting NAs in the diagonal
diag(fst_BYpval) <- NA
# And export
write.csv(fst_BYpval,file="Fst-pvalues_NEI_999boot_BY_correction.csv",quote = F)

# Plotting a heatmap now from tutorial: https://slowkow.com/notes/pheatmap-tutorial/
# The viridis color palettes: https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

# Split the columns into 6 groups/populations
col_groups <- substr(colnames(mat.obs), 1, 6)
table(col_groups)

# Data frame with column annotations
mat_col <- data.frame(Populations = col_groups)
rownames(mat_col) <- colnames(mat.obs)
mat_col

# List with colors for each annotation
mat_colors <- list(Population = brewer.pal(n = nPop(data_final), name = "Dark2"))
names(mat_colors$Population) <- unique(col_groups)
mat_colors

pdf('19_Fst_values_heatmap.pdf', width=10,height=6)
pheatmap(
  mat   = mat.obs,
  color = viridis(1000),
  # color = brewer.pal(n = 9, name = "Reds"), # If using RcolorBrewer
  border_color  = NA,
  show_colnames = TRUE,
  show_rownames = TRUE,
  # annotation_row= mat_col,
  annotation_col= mat_col,
  annotation_colors = mat_colors,
  drop_levels   = TRUE,
  fontsize  = 14,
  main  = "Fst values",
  cutree_cols   = 3,
  cutree_rows   = 3,
  cluster_rows  = T,
  display_numbers   = T
)
dev.off()

# Now that the whole pipeline has been run, let's save the R environment file
save.image("~/Documents/Antipathella/R_script_validation/ALL_LOCI/Pop_genomics_ALL_LOCI.RData")

#########################
    ### THE END ###
#########################

# Figure 2 was created using the following tutorial: https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial

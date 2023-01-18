## Changelog Phil 11/23/2022:
# - changed deprecated dplyr code
# - removed downweighting trait in gowdis
# - added decouple
# - added Bloombergs K
# - added PD and FD correlations


## Arboreal ant communities Wanang 2022- Tree scale Script

# In this script, we compare the functional and phylogenetic diversity of two plots, one in primary and one secondary forest.
# We analyse the communities on two scales, the tree scale and the plot scale. 
# The plot scale sums up all incidence across the two plots, while the tree scale treats each tree as independent communities.

# This part of the script deals only with the tree scale analysis.

### Associated .csv files:
# traits.raw.csv: Contains raw trait measurements of ant individuals
# env.csv: Environmental metadata, including forest type, tree IDs, tree DBH, etc... Important for the analysis is just forest type (secondary, primary)
# comm.csv: incidence matrix of all ants on each tree
# nestcomm.csv: incidence matrix of all nesting ant species on each tree
# fcomm.csv: incidence matrix of all visiting species, calculated as community - nesters.
# newphylodist.csv: phylogenetic distance matrix with out groups removed.


#----------------------------------------------------------#
### List of R-packages
#----------------------------------------------------------#

package_list <- 
  c("dplyr",
    "tidyr",
    "FD",
    "vegan",
    "picante",
    "stringr",
    "ggplot2",
    "ggpubr",
    "phytools",
    "corrplot",
    "PerformanceAnalytics")

# install all packages
#sapply(package_list, install.packages, character.only = TRUE)

# load all packages
sapply(package_list, library, character.only = TRUE)

# Citations
sapply(package_list, citation)

# set seed for randomization
set.seed(123)

#----------------------------------------------------------#
# Trait raw data transformation -----
#----------------------------------------------------------#

# Import Trait Table: this one has already corrected size corrected traits
traits.raw =read.csv(file="traits.raw.csv", header=T)
# load phylogenetic data
phylo.matrix <- read.csv("newphylodist.csv", row.names = 1) #outgroups removed

#### ACTIVATE CODE TO REMOVE NON-NATIVE SPECIES ###
#remove invasives
#traits.raw =read.csv(file="traits.raw.csv", header=T)
#invasiv<-subset(traits.raw, Invasive==1)
#invasiv <- invasiv %>%
#  dplyr::select(SpCode)%>%
#  mutate(SpCode = str_remove_all(SpCode, " "))
#remove<-invasiv[,1]
#phylo.matrix<-phylo.matrix[!rownames(phylo.matrix) %in% remove,!colnames(phylo.matrix) %in% remove ]

## If you activate this line, all subsequent analyses include only native species
#traits.raw<-subset(traits.raw, Invasive==2)

### END OF CODE FPR INVASIVE EXcLUSION ##

# select traits for analysis
traits.raw <- traits.raw %>%
  dplyr::select(SpCode,Caste,HeadWidth:Polymorphism)%>%
  # remove spaces in species code
  mutate(SpCode = str_remove_all(SpCode, " "))

# Leg length is approximated as hind tibia + hind femur
traits.raw$LegLength <- traits.raw$HindTibia+traits.raw$HindFemur

# Eye Size is calculated as the area of an ellipse
traits.raw$EyeSize <- ((traits.raw$EyeWidth)/2)*((traits.raw$EyeLength)/2)*3.142


traits_all<-traits.raw %>%
  # define relative measurements by dividing through HeadLength for size related traits
  mutate_at(.funs = list(rel = ~./HeadLength), 
            .vars = vars(HeadWidth, ClypeusLength:InterocularDistance, LegLength, EyeSize))%>%
  dplyr::select(SpCode,Caste, Spinosity:EyeSize_rel, HeadLength)

# For polymorphic species: weights for mean trait calculations: 0.8 minor, 0.2 major
traits_all$weights[traits_all$Caste == "minor"] <- 0.8
traits_all$weights[traits_all$Caste == "major"] <- 0.2
traits_all$weights[traits_all$Caste == "worker"] <- 1

#summarize polymorphic species with weights
species_traits<-traits_all %>%
  # summarize mean traits for each species
  group_by(SpCode) %>%
  summarise(across(
    .cols = where(is.numeric), 
    .fns = list(mean = ~mean(., na.rm = TRUE),
                weighted = ~weighted.mean(., w = weights, na.rm = TRUE)),
    .names = "{col}_{fn}"))%>%
  # we select only traits that are of interest and rename them
  dplyr::select(SpCode,Spinosity_mean, Sculpturing_mean,HeadWidth_rel_weighted, ClypeusLength_rel_weighted,MandibleLength_rel_weighted, InterocularDistance_rel_weighted, EyeSize_rel_weighted,LegLength_rel_weighted, HeadLength_weighted)%>% 
  rename(Scul = Sculpturing_mean,
         Spines=Spinosity_mean,
         HL = HeadLength_weighted,
         Rel.HW = HeadWidth_rel_weighted,
         Rel.CL = ClypeusLength_rel_weighted,
         Rel.ML = MandibleLength_rel_weighted,
         Rel.LL= LegLength_rel_weighted,
         Rel.ES= EyeSize_rel_weighted,
         Rel.EP = InterocularDistance_rel_weighted)

polymorphism<-traits_all %>%
  # calculate polymorphisms index
  group_by(SpCode) %>%
  summarize(poly.index=max(HeadLength, na.rm=T)/min(HeadLength, na.rm = T))

traitmatrix<- full_join(species_traits, polymorphism)

#transform tibble to df
finaltraitmatrix<-as.data.frame(traitmatrix)

#set species code as rownames
rownames(finaltraitmatrix) <- finaltraitmatrix[, 1]
finaltraitmatrix <- finaltraitmatrix[, -c(1)]

# log all traits for better data distribution
#finaltraitmatrix<-finaltraitmatrix %>%
#  mutate_at(.funs = ~log(.x+1), 
#            .vars = vars(Spines:poly.index))

## load environmental data
env<- read.csv("env.csv")

# # SQRT transform phylogeny as advised by Letten&Cornwell 2015 MEE
phy.dis <- (sqrt(as.dist(phylo.matrix))) 

# transform phylogenetic distance into phylo format
phylo<-as.phylo(hclust(phy.dis))

# remove species without phylogenetic info (Lordomyrma 001)
remove2<-"LORD001"
finaltraitmatrix.p<-finaltraitmatrix[!rownames(finaltraitmatrix) %in% remove2,]

# Trait phylogenetic signal Bloombergs K
kSpines <- phylosig(phylo, finaltraitmatrix.p$Spines, method="K", test=TRUE, nsim=999)
kScul <- phylosig(phylo, finaltraitmatrix.p$Scul, method="K", test=TRUE, nsim=999)
kRel.HW <- phylosig(phylo, finaltraitmatrix.p$Rel.HW, method="K", test=TRUE, nsim=999)
kRel.CL <- phylosig(phylo, finaltraitmatrix.p$Rel.CL, method="K", test=TRUE, nsim=999)
kRel.ML <- phylosig(phylo, finaltraitmatrix.p$Rel.ML, method="K", test=TRUE, nsim=999)
kRel.EP <- phylosig(phylo, finaltraitmatrix.p$Rel.EP, method="K", test=TRUE, nsim=999)
kRel.ES <- phylosig(phylo, finaltraitmatrix.p$Rel.ES, method="K", test=TRUE, nsim=999)
kRel.LL<- phylosig(phylo, finaltraitmatrix.p$Rel.LL, method="K", test=TRUE, nsim=999)
kHL<- phylosig(phylo, finaltraitmatrix.p$HL, method="K", test=TRUE, nsim=999)
kpoly.index <- phylosig(phylo, finaltraitmatrix.p$poly.index, method="K", test=TRUE, nsim=999)
## table with traits and their Bloombergs K 
traits.K <- t(data.frame(cbind(kSpines,kScul,kRel.HW,kRel.CL,kRel.ML,kRel.EP,kRel.ES,kRel.LL, kHL,kpoly.index)))

## Check correlations between standardized traits
# Plot
chart.Correlation(finaltraitmatrix.p, histogram = TRUE, pch = 19)

#----------------------------------------------------------#
# TREE SCALE: formatting and uploading occurrences  -----
#----------------------------------------------------------#
## Upload community data

# all
comm <- read.csv("comm.csv")
comm.mat.env<-merge(env, comm, by.x = "TreeN", by.y = "TreeN")
comm.mat <- comm

comm.mat<-as.data.frame(comm.mat)
rownames(comm.mat) <- comm.mat[,1]
comm.mat <- comm.mat[,-c(1)]

# nesting only
nestcomm <- read.csv("nestcomm.csv")
nestcomm.mat.env<-merge(env, nestcomm, by.x = "TreeN", by.y = "TreeN")
nestcomm.mat <- nestcomm


nestcomm.mat<-as.data.frame(nestcomm.mat)
rownames(nestcomm.mat) <- nestcomm.mat[,1]
nestcomm.mat <- nestcomm.mat[,-c(1)]

# foragers only
fcomm <- read.csv("fcomm.csv")
fcomm<-fcomm %>% 
  rename(
    TreeN = X,
  )
fcomm.mat.env<-merge(env, fcomm, by.x = "TreeN", by.y = "TreeN")
fcomm.mat <- fcomm


fcomm.mat<-as.data.frame(fcomm.mat)
rownames(fcomm.mat) <- fcomm.mat[,1]
fcomm.mat <- fcomm.mat[,-c(1)]

# remove species with 0 occurences
i <- (colSums(comm.mat) != 0) 
comm.mat <- comm.mat[,i] 

i <- (colSums(nestcomm.mat) != 0) 
nestcomm.mat <- nestcomm.mat[,i] 

i <- (colSums(fcomm.mat) != 0) 
fcomm.mat <- fcomm.mat[,i] 

# remove trees with 0 occurences
i <- (rowSums(comm.mat) != 0) 
comm.mat <- comm.mat[i,] 

i <- (rowSums(nestcomm.mat) != 0) 
nestcomm.mat <- nestcomm.mat[i,] 

i <- (rowSums(fcomm.mat) != 0) 
fcomm.mat <- fcomm.mat[i,] 

### merge community table with trait table to remove species without trait data
t.comm.mat = t(comm.mat)
t.nestcomm.mat = t(nestcomm.mat)
t.fcomm.mat = t(fcomm.mat)
#
traits.comm.mat<-merge(finaltraitmatrix, t.comm.mat, by.x = 0, by.y= 0)
traits.nestcomm.mat<-merge(finaltraitmatrix, t.nestcomm.mat, by.x = 0, by.y= 0)
traits.fcomm.mat<-merge(finaltraitmatrix, t.fcomm.mat, by.x = 0, by.y= 0)

#----------------------------------------------------------#
# Functional diversity calculations  -----
#----------------------------------------------------------#

# Here we calculate functional diversity, i.e. Raos Q and CWMs of single traits.
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

### 1. for complete community
traits_comm<-traits.comm.mat[, c(1:11)]
rownames(traits_comm) <- traits_comm[, 1]
traits_comm <- traits_comm[, -c(1)]

abun_comm<-traits.comm.mat[, -c(2:11)]
rownames(abun_comm) <- abun_comm[, 1]
abun_comm <- abun_comm[, -c(1)]

# remove trees with 0 occurrences
i <- (colSums(abun_comm) != 0) 
abun_comm <- abun_comm[,i] 


t.abun_comm<-t(abun_comm)

# calcuate SES Rao Q
multi.dis.c <- gowdis(traits_comm) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)

# CHECK: names equal?
rownames(traits_comm) == rownames(abun_comm)

# Create CWMs
c <- dbFD(traits_comm, t(abun_comm),
          corr = "cailliez",
          calc.FGR = F,
          calc.FRic= F,
          clust.type = "ward"
)

# Extract CWMs
CWM_comm <- as.data.frame(c$CWM)
#
CWM_comm$TreeN <- rownames(CWM_comm)
rownames(CWM_comm) <- NULL

### 2. for nesting community only
traits_nest<-traits.nestcomm.mat[, c(1:11)]
rownames(traits_nest) <- traits_nest[, 1]
traits_nest <- traits_nest[, -c(1)]

abun_nest<-traits.nestcomm.mat[, -c(2:11)]
rownames(abun_nest) <- abun_nest[, 1]
abun_nest <- abun_nest[, -c(1)]


# SES Rao Q nesters
traits_nest<-traits.nestcomm.mat[, c(1:11)]
rownames(traits_nest) <- traits_nest[, 1]
traits_nest <- traits_nest[, -c(1)]

abun_nest<-traits.nestcomm.mat[, -c(2:11)]
rownames(abun_nest) <- abun_nest[, 1]
abun_nest <- abun_nest[, -c(1)]

# remove trees with 0 occurrences
i <- (colSums(abun_nest) != 0) 
abun_nest <- abun_nest[,i] 
# 
t.abun_nest<-t(abun_nest)

multi.dis.n <- gowdis(traits_nest[,-c(11,12)]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)


# CHECK: names equal?
rownames(traits_nest) == rownames(abun_nest)

# Create functional diversity indices for each tree
nest <- dbFD(traits_nest, t(abun_nest),
             corr = "cailliez",
             calc.FGR = F, 
             calc.FRic = F,
             clust.type = "ward")

# Extract FD-Indices
CWM_nest <- as.data.frame(nest$CWM)
# 
CWM_nest$TreeN <- rownames(CWM_nest)
rownames(CWM_nest) <- NULL

### 3. For foragers
traits_foragers<-traits.fcomm.mat[, c(1:11)]
rownames(traits_foragers) <- traits_foragers[, 1]
traits_foragers <- traits_foragers[, -c(1)]

abun_foragers<-traits.fcomm.mat[, -c(2:11)]
rownames(abun_foragers) <- abun_foragers[, 1]
abun_foragers <- abun_foragers[, -c(1)]

# remove trees with 0 occurences
i <- (colSums(abun_foragers) != 0) 
abun_foragers <- abun_foragers[,i] 
#
t.abun_forager <- t(abun_foragers)

# SES Rao Q
multi.dis.f <- gowdis(traits_foragers) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)

# CHECK: names equal?
rownames(traits_foragers) == rownames(abun_foragers)

# Create functional diversity indices for each tree
foragers <- dbFD(traits_foragers, t(abun_foragers),
                 corr = "cailliez",
                 calc.FGR = F, 
                 calc.FRic = F,
                 clust.type = "ward")

# Extract FD-Indices
CWM_foragers <- as.data.frame(foragers$CWM)
# Rename values
CWM_foragers$TreeN <- rownames(CWM_foragers)
rownames(CWM_foragers) <- NULL

### merge CWMs together into one table
CWM_comm$Type<-"all"
CWM_nest$Type<-"nest"
CWM_foragers$Type<-"foragers"
#
fd_all <- full_join(CWM_comm, CWM_nest)
fd_all2<-full_join(fd_all, CWM_foragers)
#
fd_env<- full_join(env, fd_all2)

### merge trait diversity results into one table
mpd.comm$TreeN<-rownames(mpd.comm)
mpd.nest$TreeN<-rownames(mpd.nest)
mpd.foragers$TreeN<-rownames(mpd.foragers)
#
mpd_comm<-mpd.comm[,c(2,6,9)]
mpd_nest<-mpd.nest[,c(2,6,9)]
mpd_forager<-mpd.foragers[,c(2,6,9)]
#
mpd_comm$Type<-"all"
mpd_nest$Type<-"nest"
mpd_forager$Type<-"foragers"
#
mpd_all <- full_join(mpd_comm, mpd_nest)
mpd_all2<-full_join(mpd_all, mpd_forager)
#
mpd_env<- full_join(env,mpd_all2, by="TreeN")

# Is FD different from 0?
# PK: for these all tests perhaps better to use pne-sample Wilcoxon’s signed-rank test as in 
#Insect Conservation and Diversity (2018) doi: 10.1111/icad.12326?
# http://www.sthda.com/english/wiki/one-sample-wilcoxon-signed-rank-test-in-r
# primary all ***
mpd.z<-mpd_env %>% filter(Forest=="Primary" & Type=="all")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes;t = -11.127, df = 382, p-value < 2.2e-16

# primary foragers ***
mpd.z<-mpd_env %>% filter(Forest=="Primary" & Type=="foragers")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes ;t = -7.4056, df = 283, p-value = 1.505e-12

# primary nester ***
mpd.z<-mpd_env %>% filter(Forest=="Primary" & Type=="nest")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = -5.749, df = 172, p-value = 4.009e-08

# secondary all, ns
mpd.z<-mpd_env %>% filter(Forest=="Secondary" & Type=="all")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no; t = 1.0253, df = 278, p-value = 0.3061

# secondary forager, ns
mpd.z<-mpd_env %>% filter(Forest=="Secondary" & Type=="foragers")%>% dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no; t = 1.3286, df = 190, p-value = 0.1856

# secondary nester, ns
mpd.z<-mpd_env %>% filter(Forest=="Secondary" & Type=="nest") %>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no; t = -0.78796, df = 135, p-value = 0.4321


# Here we calculate the decoupled functional diversity, i.e. Rao Q, according to de Bello et al. (2017)
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

# function 'decouple' from de Bello et a. (2017)
d1 <- decouple(finaltraitmatrix.p, phylo, sqrt_tree = FALSE)

# extract decoupled phylogeny
dcFDis<-d1$dcFdist

#Rao Q SES
dcFDis.full <- ses.mpd(t.abun_comm, dcFDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
dcFDis.nest <- ses.mpd(t.abun_nest, dcFDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
dcFDis.forager <- ses.mpd(t.abun_forager, dcFDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
dcFDis.full$TreeN<-rownames(dcFDis.full)
dcFDis.nest$TreeN<-rownames(dcFDis.nest)
dcFDis.forager$TreeN<-rownames(dcFDis.forager)

## merge results into table to plot them later on
dcFDis_comm<-dcFDis.full[,c(2,6,9)]
dcFDis_nest<-dcFDis.nest[,c(2,6,9)]
dcFDis_forager<-dcFDis.forager[,c(2,6,9)]
#
dcFDis_comm$Type<-"all"
dcFDis_nest$Type<-"nest"
dcFDis_forager$Type<-"foragers"
#
dcFDis_all <- full_join(dcFDis_comm, dcFDis_nest)
dcFDis_all2<-full_join(dcFDis_all, dcFDis_forager)
#
dcFDis_env<- full_join(env,dcFDis_all2, by="TreeN")

# 
dcFDis_env<-dcFDis_env%>% drop_na(Type)

# Is decoupled FD different from 0?
# primary all ***
mpd.z<-dcFDis_env %>% filter(Forest=="Primary" & Type=="all")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes;t = -7.4309, df = 382, p-value = 7.133e-13

# primary foragers ***
mpd.z<-dcFDis_env %>% filter(Forest=="Primary" & Type=="foragers")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes ;t = -2.9212, df = 283, p-value = 0.003768

# primary nester ***
mpd.z<-dcFDis_env %>% filter(Forest=="Primary" & Type=="nest")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = -5.3582, df = 172, p-value = 2.669e-07

# secondary all  ***
mpd.z<-dcFDis_env %>% filter(Forest=="Secondary" & Type=="all")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = 7.7603, df = 278, p-value = 1.623e-13

# secondary forager,   ***
mpd.z<-dcFDis_env %>% filter(Forest=="Secondary" & Type=="foragers")%>% dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no; t = 7.1183, df = 190, p-value = 2.181e-11

# secondary nester,  no
mpd.z<-dcFDis_env %>% filter(Forest=="Secondary" & Type=="nest") %>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no; t = 1.7861, df = 135, p-value = 0.07632

#----------------------------------------------------------#
# Phylogenetic diversity calculations  -----
#----------------------------------------------------------#

# Here we calculate phylogenetic diversity, i.e. Rao Q
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

# Phylogenetic SES Rao Q
 
#Rao Q SES
ses.mpd.phylo.full <- ses.mpd(t.abun_comm, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.nest <- ses.mpd(t.abun_nest, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.forager <- ses.mpd(t.abun_forager, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
ses.mpd.phylo.full$TreeN<-rownames(ses.mpd.phylo.full)
ses.mpd.phylo.nest$TreeN<-rownames(ses.mpd.phylo.nest)
ses.mpd.phylo.forager$TreeN<-rownames(ses.mpd.phylo.forager)

## merge results into table to plot them later on
pd_comm<-ses.mpd.phylo.full[,c(2,6,9)]
pd_nest<-ses.mpd.phylo.nest[,c(2,6,9)]
pd_forager<-ses.mpd.phylo.forager[,c(2,6,9)]
#
pd_comm$Type<-"all"
pd_nest$Type<-"nest"
pd_forager$Type<-"foragers"
#
pd_all <- full_join(pd_comm, pd_nest)
pd_all2<-full_join(pd_all, pd_forager)
#
pd_env<- full_join(env,pd_all2, by="TreeN")

### Is PD different from 0?

# primary all, ns
mpd.z<-pd_env %>% filter(Forest=="Primary" & Type=="all")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no; t = 2.0838, df = 382, p-value = 0.03785


# primary foragers *
mpd.z<-pd_env %>% filter(Forest=="Primary" & Type=="foragers")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes ; t = 2.4058, df = 283, p-value = 0.01678

# primary nester, ns
mpd.z<-pd_env %>% filter(Forest=="Primary" & Type=="nest")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no; t = -1.6494, df = 172, p-value = 0.1009

# secondary all, ***
mpd.z<-pd_env %>% filter(Forest=="Secondary" & Type=="all")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = -4.4287, df = 278, p-value = 1.364e-05

# secondary forager, ***
mpd.z<-pd_env %>% filter(Forest=="Secondary" & Type=="foragers")%>% dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = -3.9873, df = 190, p-value = 9.527e-05

# secondary nester, ns
mpd.z<-pd_env %>% filter(Forest=="Secondary" & Type=="nest") %>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no; t = 0.98817, df = 135, p-value = 0.3248


# Here we calculate the decoupled phylogenetic diversity, i.e. Rao Q, according to de Bello et al. (2017)
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

# function 'decouple' from de Bello et a. (2017)
d1 <- decouple(finaltraitmatrix.p, phylo, sqrt_tree = FALSE)

# extract decoupled phylogeny
dcPDis<-d1$dcPdist

#Rao Q SES
dcPDis.full <- ses.mpd(t.abun_comm, dcPDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
dcPDis.nest <- ses.mpd(t.abun_nest, dcPDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
dcPDis.forager <- ses.mpd(t.abun_forager, dcPDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
dcPDis.full$TreeN<-rownames(dcPDis.full)
dcPDis.nest$TreeN<-rownames(dcPDis.nest)
dcPDis.forager$TreeN<-rownames(dcPDis.forager)

## merge results into table
dcPDis_comm<-dcPDis.full[,c(2,6,9)]
dcPDis_nest<-dcPDis.nest[,c(2,6,9)]
dcPDis_forager<-dcPDis.forager[,c(2,6,9)]
#
dcPDis_comm$Type<-"all"
dcPDis_nest$Type<-"nest"
dcPDis_forager$Type<-"foragers"
#
dcPDis_all <- full_join(dcPDis_comm, dcPDis_nest)
dcPDis_all2<-full_join(dcPDis_all, dcPDis_forager)
#
dcPDis_env<- full_join(env,dcPDis_all2, by="TreeN")

# Plot 
dcPDis_env<-dcPDis_env%>% drop_na(Type)

# Is decoupled PD different from 0?
# primary all ***
mpd.z<-dcPDis_env %>% filter(Forest=="Primary" & Type=="all")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes;t = -11.127, df = 382, p-value < 2.2e-16

stat.test <- dcPDis_env %>% filter(Forest=="Primary" & Type=="all")%>%dplyr::select(mpd.obs.z) %>% wilcox_test(weight ~ 1, mu = 0)
stat.test


# primary foragers, ns
mpd.z<-dcPDis_env %>% filter(Forest=="Primary" & Type=="foragers")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# no ;t = -1.2316, df = 283, p-value = 0.2191

# primary nester ***
mpd.z<-dcPDis_env %>% filter(Forest=="Primary" & Type=="nest")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = -2.7531, df = 172, p-value = 0.006537

# secondary all,  ***
mpd.z<-dcPDis_env %>% filter(Forest=="Secondary" & Type=="all")%>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = -6.2028, df = 278, p-value = 2.003e-09

# secondary forager, ,  ***
mpd.z<-dcPDis_env %>% filter(Forest=="Secondary" & Type=="foragers")%>% dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = -3.4831, df = 190, p-value = 0.0006151

# secondary nester,  ***
mpd.z<-dcPDis_env %>% filter(Forest=="Secondary" & Type=="nest") %>%dplyr::select(mpd.obs.z)
t.test(mpd.z, alternative="two.sided")
# yes; t = -5.5193, df = 135, p-value = 1.681e-07


#----------------------------------------------------------#
# Correlation test of FD and PD  -----
#----------------------------------------------------------#

# First assemble SES for PD and FD for same trees
# for all
FD.mpd<-mpd_env %>% filter(Type=="all") %>% na.omit() %>%dplyr::select(mpd.obs, mpd.obs.z, TreeN, Forest, AntSpRichness.ALL)%>% 
  rename(FD.RaoQ = mpd.obs.z, FD.RaoQ.obs=mpd.obs)

PD.mpd<-pd_env %>% filter(Type=="all") %>% na.omit() %>%dplyr::select(mpd.obs, mpd.obs.z, TreeN)%>% 
  rename(PD.RaoQ = mpd.obs.z, PD.RaoQ.obs=mpd.obs)

rownames(FD.mpd)<-FD.mpd[,3]
rownames(PD.mpd)<-PD.mpd[,3]

FD.PD<-merge(FD.mpd, PD.mpd, by = 'row.names', all = TRUE)

# check for correlation
corr <- cor.test(FD.PD$RaoQ, FD.PD$RaoQ, method = 'spearman')
corr
# S = 34984785, p-value = 4.437e-13
# rho = 0.2764684

# linear model
model.dc.fd.pd<-lm(FD.RaoQ~PD.RaoQ*Forest, data=FD.PD)
summary(model.dc.fd.pd)

# linear model only primary forest
model.fd.pd.pr<-lm(FD.RaoQ~PD.RaoQ+AntSpRichness.ALL, data=subset(FD.PD, Forest == "Primary"))
summary(model.fd.pd.pr)

# linear model only secondary forest
model.fd.pd.sec<-lm(FD.RaoQ~PD.RaoQ+AntSpRichness.ALL, data=subset(FD.PD, Forest == "Secondary"))
summary(model.fd.pd.sec)


# Spearman rank correlation test for correlation between dcFD and dcPD
# First assemble SES for PD and FD for same trees
# for all
dcFD.mpd<-dcFDis_env %>% filter(Type=="all") %>% na.omit() %>%dplyr::select(mpd.obs, mpd.obs.z, TreeN, Forest, AntSpRichness.ALL)%>% 
  rename(dcFD.RaoQ = mpd.obs.z, dcFD.RaoQ.obs = mpd.obs)

dcPD.mpd<-dcPDis_env %>% filter(Type=="all") %>% na.omit() %>%dplyr::select(mpd.obs, mpd.obs.z, TreeN)%>% 
  rename(dcPD.RaoQ = mpd.obs.z, dcPD.RaoQ.obs = mpd.obs)

rownames(dcFD.mpd)<-dcFD.mpd[,3]
rownames(dcPD.mpd)<-dcPD.mpd[,3]

dcFD.PD<-merge(dcFD.mpd, dcPD.mpd, by = 'row.names', all = TRUE)

# check for correlation
corr <- cor.test(dcFD.PD$dcFD.RaoQ, dcFD.PD$dcPD.RaoQ, method = 'spearman')
corr
# S = 34984785, p-value = 4.437e-13
# rho = 0.2764684

# linear model
model.dcfd.pd<-lm(dcFD.RaoQ~dcPD.RaoQ*Forest, data=dcFD.PD)
summary(model.dcfd.pd)


#----------------------------------------------------------#
# Species richness  -----
#----------------------------------------------------------#

# merge communities with trait data to remove species without traits
t.comm = t(comm)
colnames(t.comm)<-t.comm[1,]
t.comm<-t.comm[-1,]

t.nestcomm = t(nestcomm)
t.fcomm = t(fcomm)
#
traits.comm<-merge(finaltraitmatrix, t.comm, by.x = 0, by.y= 0)
traits.nestcomm<-merge(finaltraitmatrix, t.nestcomm, by.x = 0, by.y= 0)
traits.fcomm<-merge(finaltraitmatrix, t.fcomm, by.x = 0, by.y= 0)
#
abun_comm2<-traits.comm[, -c(2:11)]
rownames(abun_comm2) <- abun_comm2[, 1]
abun_comm2 <- abun_comm2[, -c(1)]
t.abuncomm2<-t(abun_comm2)
#
nest_comm2<-traits.nestcomm[, -c(2:11)]
rownames(nest_comm2) <- nest_comm2[, 1]
nest_comm2 <- nest_comm2[, -c(1)]
t.nest_comm2<-t(nest_comm2)
#
f_comm2<-traits.fcomm[, -c(2:11)]
rownames(f_comm2) <- f_comm2[, 1]
f_comm2 <- f_comm2[, -c(1)]
t.f_comm2<-t(f_comm2)

# Get species richness for each tree
forest_c<-env %>% dplyr::select(TreeN, Forest)
forest_c$Richness <- specnumber(t.abuncomm2)
forest_c$Type<-"all"
#
forest_n<-env %>% dplyr::select(TreeN, Forest)
forest_n$Richness <- specnumber(t.nest_comm2)
forest_n$Type<-"nest"
#
forest_f<-env %>% dplyr::select(TreeN, Forest)
forest_f$Richness <- specnumber(t.f_comm2)
forest_f$Type<-"forager"

# get all results into one table
all_r<-rbind(forest_c,forest_n, forest_f)

#----------------------------------------------------------#
# FIGURES -----
#----------------------------------------------------------#

# NOTE: All figures have the statistical comparisons implemented as Kruskal-Wallis tests via ggpubr
# all plots shown here are tree scale analysis - the plot scale results (see Part 2) are only reported as tables

#### Functional Diversity figures

# FD SES RAO Q
mpd_env<-mpd_env%>% drop_na(Type)
pd_env<-pd_env%>% drop_na(Type)

mpd.fd<-ggplot(mpd_env, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Functional Diversity") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  #ylim(-3, 4)+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
mpd.fd

# FD observed RAO Q
mpd.fd.obs<-ggplot(mpd_env, aes(x=Type, y=mpd.obs, fill=Forest))+
  ggtitle("Functional Diversity") +
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+ 
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  #ylim(-3, 4)+
  ylab("Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
mpd.fd.obs

# FD decoupled SES RAO Q 
mpd_dcFDis<-ggplot(dcFDis_env, aes(x=Type, y=mpd.obs.z, fill=Forest)) +
  ggtitle("Decoupled Functional Diversity") +
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+   
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  #ylim(-3, 4)+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
mpd_dcFDis

# FD observed decoupled RAO Q
mpd_dcFDis.obs<-ggplot(dcFDis_env, aes(x=Type, y=mpd.obs, fill=Forest)) +
  ggtitle("Cecoupled Functional Diversity") +
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+   
  geom_boxplot()+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  #ylim(-3, 4)+
  ylab("Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
mpd_dcFDis.obs

#### Phylogenetic Diversity

# PD SES Rao
mpd.pd<-ggplot(pd_env, aes(x=Type, y=mpd.obs.z, fill=Forest)) +
  ggtitle("Phylogenetic Diversity") +
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+   
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  #ylim(-3, 3)+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
mpd.pd

# PD observed SES Rao
mpd.pd.obs<-ggplot(pd_env, aes(x=Type, y=mpd.obs, fill=Forest)) +
  ggtitle("Phylogenetic Diversity") +
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+   
  geom_boxplot()+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  #ylim(-3, 3)+
  ylab("Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
mpd.pd.obs

# PD decoupled SES Rao
mpd_dcPDis<-ggplot(dcPDis_env, aes(x=Type, y=mpd.obs.z, fill=Forest)) +
  ggtitle("Decoupled Phylogenetic Diversity") +
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+   
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  #ylim(-3, 4)+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
mpd_dcPDis

# PD observed decoupled Rao Q
mpd_dcPDis.obs<-ggplot(dcPDis_env, aes(x=Type, y=mpd.obs, fill=Forest)) +
  ggtitle("Decoupled Phylogenetic Diversity") +
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+   
  geom_boxplot()+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
mpd_dcPDis.obs

### Correlations of FD and PD

# Correlation of SES RaoQ
corr_fd.pd<-ggplot(FD.PD, aes(x=PD.RaoQ, y=FD.RaoQ, color=Forest, shape=Forest))+
  labs(x="Phylogenetic Diversity (Rao Q SES)", y = "Functional Diversity (Rao Q SES)", title="FD and PD Correlation", fill="Forest Type")+
  scale_color_manual(values=c("#009E73", "#D55E00"))+
  geom_point(size=2)+
  #geom_smooth(data = subset(FD.PD, Forest == "Primary"), method = "lm", linetype = 2, fill = "#009E73") +
  geom_line(data = subset(FD.PD, Forest == "Primary"), stat = "smooth", method = "lm", linetype = 2, color ="#009E73", size = 1, alpha = 1) +
  #geom_smooth(data = subset(FD.PD, Forest == "Secondary"), method = "lm", linetype = 4, fill = "#D55E00") +
  geom_line(data = subset(FD.PD, Forest == "Secondary"), stat = "smooth", method = "lm", linetype = 4, color = "#D55E00", size = 1, alpha = 1) +
  theme_classic()
corr_fd.pd

# Correlation of observed RaoQ
corr_fd.pd.obs<-ggplot(FD.PD, aes(x=PD.RaoQ.obs, y=FD.RaoQ.obs, color=Forest, shape=Forest))+
  labs(x="Phylogenetic Diversity (Rao Q)", y = "Functional Diversity (Rao Q)", title="FD and PD Correlation", fill="Forest Type")+
  scale_color_manual(values=c("#009E73", "#D55E00"))+
  geom_point(size=2)+
  #geom_smooth(data = subset(FD.PD, Forest == "Primary"), method = "lm", linetype = 2, fill = "#009E73") +
  geom_line(data = subset(FD.PD, Forest == "Primary"), stat = "smooth", method = "lm", linetype = 2, color ="#009E73", size = 1, alpha = 1) +
  #geom_smooth(data = subset(FD.PD, Forest == "Secondary"), method = "lm", linetype = 4, fill = "#D55E00") +
  geom_line(data = subset(FD.PD, Forest == "Secondary"), stat = "smooth", method = "lm", linetype = 4, color = "#D55E00", size = 1, alpha = 1) +
  theme_classic()
corr_fd.pd.obs

# Correlation of decoupled SES RaoQ
corr_dcpd<-ggplot(dcFD.PD, aes(x=dcPD.RaoQ, y=dcFD.RaoQ, color=Forest, shape=Forest))+
  labs(x="Decoupled Phylogenetic Diversity (Rao Q SES)", y = "Decoupled Functional Diversity (Rao Q SES)", title="Decoupled FD and PD Correlation", fill="Forest Type")+
  scale_color_manual(values=c("#009E73", "#D55E00"))+
  geom_point(size=2)+
  #geom_smooth(data = subset(FD.PD, Forest == "Primary"), method = "lm", linetype = 2, fill = "#009E73") +
  geom_line(data = subset(dcFD.PD, Forest == "Primary"), stat = "smooth", method = "lm", linetype = 2, color ="#009E73", size = 1, alpha = 1) +
  #geom_smooth(data = subset(FD.PD, Forest == "Secondary"), method = "lm", linetype = 4, fill = "#D55E00") +
  geom_line(data = subset(dcFD.PD, Forest == "Secondary"), stat = "smooth", method = "lm", linetype = 4, color = "#D55E00", size = 1, alpha = 1) +
  theme_classic()
corr_dcpd

# Correlation of decoupled observed RaoQ
corr_dcpd.obs<-ggplot(dcFD.PD, aes(x=dcPD.RaoQ.obs, y=dcFD.RaoQ.obs, color=Forest, shape=Forest))+
  labs(x="Decoupled Phylogenetic Diversity (Rao Q)", y = "Decoupled Functional Diversity (Rao Q)", title="Decoupled FD and PD Correlation", fill="Forest Type")+
  scale_color_manual(values=c("#009E73", "#D55E00"))+
  geom_point(size=2)+
  #geom_smooth(data = subset(FD.PD, Forest == "Primary"), method = "lm", linetype = 2, fill = "#009E73") +
  geom_line(data = subset(dcFD.PD, Forest == "Primary"), stat = "smooth", method = "lm", linetype = 2, color ="#009E73", size = 1, alpha = 1) +
  #geom_smooth(data = subset(FD.PD, Forest == "Secondary"), method = "lm", linetype = 4, fill = "#D55E00") +
  geom_line(data = subset(dcFD.PD, Forest == "Secondary"), stat = "smooth", method = "lm", linetype = 4, color = "#D55E00", size = 1, alpha = 1) +
  theme_classic()
corr_dcpd.obs

### Taxonomic richness
richness<-ggplot(all_r, aes(x=Type, y=Richness, fill=Forest)) +
  ggtitle("Species richness") +
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "forager" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
richness

### Community Weighted Means 
fd_env<-fd_env%>% drop_na(Type)

HL<-ggplot(fd_env, aes(x=Type, y=HL, fill=Forest)) +
  ggtitle("Head Length") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
HL

poly.index<-ggplot(fd_env, aes(x=Type, y=poly.index, fill=Forest)) +
  ggtitle("Polymorphism index") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+         
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
poly.index

Spines<-ggplot(fd_env, aes(x=Type, y=Spines, fill=Forest)) +
  ggtitle("Spinosity") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Spines

Scul<-ggplot(fd_env, aes(x=Type, y=Scul, fill=Forest)) +
  ggtitle("Sculpture") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Scul

Rel.HW<-ggplot(fd_env, aes(x=Type, y=Rel.HW, fill=Forest)) +
  ggtitle("Rel. Headwidth") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Rel.HW

Rel.CL<-ggplot(fd_env, aes(x=Type, y=Rel.CL, fill=Forest)) +
  ggtitle("Rel. Clypeus Length") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Rel.CL

Rel.ML<-ggplot(fd_env, aes(x=Type, y=Rel.ML, fill=Forest)) +
  ggtitle("Rel.Mandible Length") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Rel.ML

Rel.EP<-ggplot(fd_env, aes(x=Type, y=Rel.EP, fill=Forest)) +
  ggtitle("Rel. Eye Position") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Rel.EP

Rel.ES<-ggplot(fd_env, aes(x=Type, y=Rel.ES, fill=Forest)) +
  ggtitle("Rel. Eye Size") +
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+    
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Rel.ES

Rel.LL<-ggplot(fd_env, aes(x=Type, y=Rel.LL, fill=Forest)) +
  ggtitle("Rel. Leg Length") +
geom_boxplot()+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+         
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Rel.LL

# Trait correlations
mydata.cor <- cor(finaltraitmatrix.p, method = c("spearman"), use = "pairwise.complete.obs")
corrplot(mydata.cor,
         method = "color", type = "lower", order = "AOE", diag = F,
         tl.col = "black", outline = T, addCoef.col = "black", number.cex = 0.8,
         tl.cex = 1.1, cl.cex = 0.9
)

#### Main Manuscript
## Figure 3 
fig3 <- ggarrange(mpd.fd, mpd.pd, corr_fd.pd, mpd_dcFDis, mpd_dcPDis, richness,
                   labels = c("A", "B", "C", "D", "E", "F"),
                   ncol = 3, nrow = 2, common.legend = TRUE, legend = "top")
fig3


## Figure 4 
fig4 <- ggarrange(HL, poly.index, Spines, Scul, Rel.HW, Rel.CL, Rel.ML, Rel.EP, Rel.ES,Rel.LL,
                           labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
                           ncol = 5, nrow = 2, common.legend = TRUE, legend = "top"
)
fig4


#### Supplement Figures


###### REMOVE LATER: Export figures ######

## Plot 4 Francesco
figure_decoupling <- ggarrange(mpd.fd, mpd.pd, corr_fd.pd, mpd_dcFDis, mpd_dcPDis, corr_dcpd,
                           labels = c("A", "B", "C", "D", "E", "F"),
                           ncol = 3, nrow = 2, common.legend = F
)
figure_decoupling

figure_decoupling.obs <- ggarrange(mpd.fd.obs, mpd.pd.obs, corr_fd.pd.obs, mpd_dcFDis.obs, mpd_dcPDis.obs, corr_dcpd.obs,
                               labels = c("A", "B", "C", "D", "E", "F"),
                               ncol = 3, nrow = 2, common.legend = F
)
figure_decoupling.obs


# Figure 1: richness + trait diversity + phylogenetic diversity
# as tiff
library(svglite)

tiff("figure_2.tiff", units ="in", width = 10, height = 8, res = 150)
fig3
dev.off()

# or as .svg (vector image)
#ggsave(file="figure_1.svg", plot=plot1, width=10, height=4)

# Figure 2: CWMs
# as tiff
tiff("figure_3.tiff", units = "in", width = 16, height = 8, res = 150)
fig4
dev.off()

# or as .svg (vector image)
ggsave(file="figure_2.svg", plot=figure_traits, width=16, height=8)


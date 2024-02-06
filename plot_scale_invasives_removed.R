#----------------------------------------------------------#
### R Script for: Forest disturbance increases functional diversity but decreases phylogenetic diversity of an arboreal tropical ant community
#----------------------------------------------------------#

# Authors: Philipp O. Hoenle, Nichola S. Plowmana, Pável Matos-Maraví, Francesco de Bello, Tom R. Bishop, Martin Libra, Cliffson Idigel, Maling Rimandai, Petr Klimes


# In this script, we compare the functional and phylogenetic diversity of two plots, one in primary and one secondary forest.
# We analyse the communities on two scales, the tree scale and the plot scale.
# The plot scale sums up all incidence across the two plots, while the tree scale treats each tree as independent communities.

# This part of the script contains the plot scale analysis, and removes the invasive species from the data.

### Associated .csv files:
# traits.raw.csv: Contains raw trait measurements of ant individuals
# env.csv: Environmental metadata, including forest type, tree IDs, tree DBH, etc... Important for the analysis is just forest type (secondary, primary)
# comm.csv: incidence matrix of all ants on each tree. Sometimes referred to as 'whole community' or 'all'
# nestcomm.csv: incidence matrix of all nesting ant species on each tree, sometimes referred to as 'nesters'
# fcomm.csv: incidence matrix of all visiting species, calculated as community - nesters. Sometimes referred to as 'visitors' or 'foragers'
# newphylodist.csv: phylogenetic distance matrix with out groups removed.

# For a list of the species linked to the morphospecies code, please see the Supplement Table S1 in our paper.


#----------------------------------------------------------#
### List of R-packages
#----------------------------------------------------------#

package_list <- 
  c("dplyr",
    "tidyr",
    "FD",
    "vegan",
    "picante",
    "stringr")

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

## load phylogenetic data
phylo.matrix <- read.csv("newphylodist.csv", row.names = 1) #outgroups removed

# remove invasive species 
traits.raw =read.csv(file="traits.raw.csv", header=T)
invasiv<-subset(traits.raw, Invasive==1)
invasiv <- invasiv %>%
  dplyr::select(SpCode)%>%
  mutate(SpCode = str_remove_all(SpCode, " "))
remove<-invasiv[,1]
phylo.matrix<-phylo.matrix[!rownames(phylo.matrix) %in% remove,!colnames(phylo.matrix) %in% remove ]
#
traits.raw<-subset(traits.raw, Invasive==2)

# select traits for analysis
traits.raw <- traits.raw %>%
  dplyr::select(SpCode,Caste,HeadWidth:Polymorphism)%>%
  # remove spaces in species code
  mutate(SpCode = str_remove_all(SpCode, " "))

# Leg length is approximated as hind tibia + hind femur
traits.raw$LegLength <- traits.raw$HindTibia+traits.raw$HindFemur

# Eye Size is calculated as the area of an ellipse
traits.raw$EyeSize <- ((traits.raw$EyeWidth)/2)*((traits.raw$EyeLength)/2)*3.142

#
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

## load environmental data
env<- read.csv("env.csv")

# # SQRT transform phylogeny as advised by Letten&Cornwell 2015 MEE
phy.dis <- (sqrt(as.dist(phylo.matrix))) 

# transform phylogenetic distance into phylo format
phylo<-as.phylo(hclust(phy.dis))

# remove species without phylogenetic info (Lordomyrma 001) from trait matrix
remove2<-"LORD001"
finaltraitmatrix.p<-finaltraitmatrix[!rownames(finaltraitmatrix) %in% remove2,]

# function 'decouple' from de Bello et al. (2017), MEE  https://doi.org/10.1111/2041-210X.12735
d1 <- decouple(finaltraitmatrix.p, phylo, sqrt_tree = FALSE)

#----------------------------------------------------------#
# Plot SCALE: formatting and uploading occurrences  -----
#----------------------------------------------------------#

## Upload & prepare community data from the three community tables
# all
comm <- read.csv("comm.csv")
comm.mat.env<-merge(env, comm, by.x = "TreeN", by.y = "TreeN")
comm.mat <- comm.mat.env %>% dplyr::select(Forest, ANOC001:VOLL001) %>% 
  group_by(Forest) %>%
  summarise_all(list(sum))

comm.mat<-as.data.frame(comm.mat)
rownames(comm.mat) <- comm.mat[,1]
comm.mat <- comm.mat[,-c(1)]

# nesting only
nestcomm <- read.csv("nestcomm.csv")
nestcomm.mat.env<-merge(env, nestcomm, by.x = "TreeN", by.y = "TreeN")
nestcomm.mat <- nestcomm.mat.env %>% dplyr::select(Forest, ANON001:VOLL001)%>%
  group_by(Forest) %>% 
  summarise_all(list(sum))

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
fcomm.mat <- fcomm.mat.env %>% dplyr::select(Forest, ANON001:VOLL001)%>%
  group_by(Forest) %>% 
  summarise_all(list(sum))


fcomm.mat<-as.data.frame(fcomm.mat)
rownames(fcomm.mat) <- fcomm.mat[,1]
fcomm.mat <- fcomm.mat[,-c(1)]


# remove species with 0 occurrences
i <- (colSums(comm.mat) != 0) 
comm.mat <- comm.mat[,i] 
#
i <- (colSums(fcomm.mat) != 0) 
fcomm.mat <- fcomm.mat[,i] 
#
i <- (colSums(nestcomm.mat) != 0) 
nestcomm.mat <- nestcomm.mat[,i] 


#----------------------------------------------------------#
# Species composition overlap  -----
#----------------------------------------------------------#
# calculate species overlap 

# between prim vs sec
dist.ftype <- vegdist(fcomm.mat, method = "bray") # 0.91
dist.ctype <- vegdist(comm.mat, method = "bray") # 0.87
dist.ntype <- vegdist(nestcomm.mat, method = "bray") # 0.82
# 
fcomm.mat2<-fcomm.mat
nestcomm.mat2<-nestcomm.mat
# 
row.names(fcomm.mat2)[row.names(fcomm.mat2) == "Primary"] <- "Primary_f"
row.names(fcomm.mat2)[row.names(fcomm.mat2) == "Secondary"] <- "Secondary_f"
row.names(nestcomm.mat2)[row.names(fcomm.mat2) == "Primary"] <- "Primary_n"
row.names(nestcomm.mat2)[row.names(fcomm.mat2) == "Secondary"] <- "Secondary_n"

# combine into one table
nester_forager<-full_join(fcomm.mat2, nestcomm.mat2)

nester_forager[is.na(nester_forager)] <- 0 # replace NA occurences with 0
row.names(nester_forager) <- c("Primary_f","Secondary_f","Primary_n", "Secondary_n")

# Distance as Bray-Curtis
dist.fn <- vegdist(nester_forager, method = "bray")
dist.fn # primary f and primary n are similar (BC =0.55), Secondary f and Secondary n less (BC =0.70)

#----------------------------------------------------------#
# Functional diversity calculations  -----
#----------------------------------------------------------## 

# Here we calculate functional diversity, i.e. Raos Q and CWMs of single traits.
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

# merge occurrence with trait data to remove species without trait data
t.comm.mat = t(comm.mat)
#
t.nestcomm.mat = t(nestcomm.mat)
#
t.fcomm.mat = t(fcomm.mat)

## traits all
traits.comm.mat<-merge(finaltraitmatrix, t.comm.mat, by.x = 0, by.y= 0)
rownames(traits.comm.mat) <- traits.comm.mat[,1]
traits.comm.mat <- traits.comm.mat[,-c(1)]
#
com.mat2<-traits.comm.mat%>%
  dplyr::select(Primary, Secondary)
com.mat2<-t(com.mat2)

## traits nesters
traits.nestcomm.mat<-merge(finaltraitmatrix, t.nestcomm.mat, by.x = 0, by.y= 0)
rownames(traits.nestcomm.mat) <- traits.nestcomm.mat[,1]
traits.nestcomm.mat <- traits.nestcomm.mat[,-c(1)]
#
nestcom.mat2<-traits.nestcomm.mat%>%
  dplyr::select(Primary, Secondary)
nestcom.mat2<-t(nestcom.mat2)

## traits foragers
traits.fcomm.mat<-merge(finaltraitmatrix,t.fcomm.mat, by.x = 0, by.y= 0)
rownames(traits.fcomm.mat) <- traits.fcomm.mat[,1]
traits.fcomm.mat <- traits.fcomm.mat[,-c(1)]
#
fcom.mat2<-traits.fcomm.mat%>%
  dplyr::select(Primary, Secondary)
fcom.mat2<-t(fcom.mat2)

# calculate functional dendrograms
multi.dis.c <- gowdis(traits.comm.mat[,-c(11,12)]) # gower distance matrix
multi.dis.n <- gowdis(traits.nestcomm.mat[,-c(11,12)]) # gower distance matrix
multi.dis.f <- gowdis(traits.fcomm.mat[,-c(11,12)]) # gower distance matrix

# calculate SES Rao Q for traits
mpd.comm.p <- ses.mpd(com.mat2, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full occurrence weighted
mpd.nest.p <- ses.mpd(nestcom.mat2, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest occurrence weighted
mpd.forager.p<-ses.mpd(fcom.mat2, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) # foragers occurrence
#
mpd.comm.p$TreeN<-rownames(mpd.comm.p)
mpd.nest.p$TreeN<-rownames(mpd.nest.p)
mpd.forager.p$TreeN<-rownames(mpd.forager.p)
#
mpd.comm.p # both ns. 
mpd.nest.p # both ns.
mpd.forager.p # # both ns.

# Same calculations for decoupled trait matrix
dcFDis<-d1$dcFdist

#Rao Q SES
dcFDis.full <- ses.mpd(com.mat2, dcFDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
dcFDis.nest <- ses.mpd(nestcom.mat2, dcFDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
dcFDis.forager <- ses.mpd(fcom.mat2, dcFDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
dcFDis.full$TreeN<-rownames(dcFDis.full)
dcFDis.nest$TreeN<-rownames(dcFDis.nest)
dcFDis.forager$TreeN<-rownames(dcFDis.forager)
#
dcFDis.full #  ns.
dcFDis.nest #  ns.
dcFDis.forager # ns.

#----------------------------------------------------------#
# Phylogenetic diversity calculations  -----
#----------------------------------------------------------#

# Here we calculate phylogenetic diversity, i.e. Rao Q
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

# Phylogenetic SES RaoQ

# Rao Q compared with a null model for SES
ses.mpd.phylo.full.p <- ses.mpd(com.mat2, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.nest.p <- ses.mpd(nestcom.mat2, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.forager.p <- ses.mpd(fcom.mat2, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
ses.mpd.phylo.full.p$TreeN<-rownames(ses.mpd.phylo.full.p)
ses.mpd.phylo.nest.p$TreeN<-rownames(ses.mpd.phylo.nest.p)
ses.mpd.phylo.forager.p$TreeN<-rownames(ses.mpd.phylo.forager.p)
#
ses.mpd.phylo.full.p # both ns, primary random, secondary rather clustered
ses.mpd.phylo.nest.p # # both ns, primary slightly clustered, secondary random
ses.mpd.phylo.forager.p  #  both ns, primary random, secondary rather clustered

# same for decoupled phylogeny
# 
dcPDis<-d1$dcPdist

#Rao Q SES
dcPDis.full <- ses.mpd(com.mat2, dcPDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
dcPDis.nest <- ses.mpd(nestcom.mat2, dcPDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
dcPDis.forager <- ses.mpd(fcom.mat2, dcPDis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
dcPDis.full$TreeN<-rownames(dcPDis.full)
dcPDis.nest$TreeN<-rownames(dcPDis.nest)
dcPDis.forager$TreeN<-rownames(dcPDis.forager)
#
dcPDis.full # ns
dcPDis.nest # ns
dcPDis.forager # ns

#----------------------------------------------------------#
# Species richness  -----
#----------------------------------------------------------#
forest_c<-as.data.frame(c("Primary", "Secondary"))
forest_c$Richness_all <- specnumber(com.mat2)
forest_c$Richness_forager<- specnumber(nestcom.mat2)
forest_c$Richness_nester<- specnumber(fcom.mat2)

forest_c

#----------------------------------------------------------#
# Save R-data -----
#----------------------------------------------------------#

# 
save.image("plot-scale_invasives_removed.Rdata")

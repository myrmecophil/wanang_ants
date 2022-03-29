## Arboreal ant communities Wanang 2022- Plot scale Script

# In this script, we compare the functional and phylogenetic diversity of two plots, one in primary and one secondary forest.
# We analyse the communities on two scales, the tree scale and the plot scale.
# The plot scale sums up all incidence across the two plots, while the tree scale treats each tree as independent communities.

# This part of the script deals only with the plot scale analysis.


### Associated .csv files:
# traits.raw.csv: Contains raw trait measurements of ant individuals
# env.csv: Environmental metadata, including forest type, tree IDs, tree DBH, etc... 
# comm.csv: incidence matrix of all ants on each tree
# nestcomm.csv: incidence matrix of all nesting ant species on each tree
# fcomm.csv: incidence matrix of all visiting species, calculated by community - nesters.
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
    "stringr")

# install all packages
#sapply(package_list, install.packages, character.only = TRUE)

# load all packages
sapply(package_list, library, character.only = TRUE)

# Citations
sapply(package_list, citation)


#----------------------------------------------------------#
# Trait raw data transformation -----
#----------------------------------------------------------#

# Import Trait Table: this one has already corrected size corrected traits
traits.raw =read.csv(file="traits.raw.csv", header=T)

## SELECT TO REMOVE INVASIVE Species. 
## If you activate this line, all subsequent analyses except species composition include only non-invasive species
#traits.raw<-subset(traits.raw, Invasive==2)

# select traits for analysis
traits.raw <- traits.raw %>%
  select(SpCode,Caste,HeadWidth:Polymorphism)%>%
  # remove spaces in species code
  mutate(SpCode = str_remove_all(SpCode, " "))

# Leg length is approximated as hind tibia + hind femur
traits.raw$LegLength <- traits.raw$HindTibia+traits.raw$HindFemur

# Eye Size is calculated as the area of an ellipse
traits.raw$EyeSize <- ((traits.raw$EyeWidth)/2)*((traits.raw$EyeLength)/2)*3.142


traits_all<-traits.raw %>%
  # define relative measurements by dividing through HeadLength for size related traits
  mutate_at(.funs = funs(rel = ./HeadLength), 
            .vars = vars(HeadWidth, ClypeusLength:InterocularDistance, LegLength, EyeSize))%>%
  select(SpCode,Caste, Spinosity:EyeSize_rel, HeadLength)

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
  select(SpCode,Spinosity_mean, Sculpturing_mean,HeadWidth_rel_weighted, ClypeusLength_rel_weighted,MandibleLength_rel_weighted, InterocularDistance_rel_weighted, EyeSize_rel_weighted,LegLength_rel_weighted, HeadLength_weighted)%>% 
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

#----------------------------------------------------------#
# Plot SCALE: formatting and uploading occurrences  -----
#----------------------------------------------------------#

## Upload & prepare community data from the three community tables
# all
comm <- read.csv("comm.csv")
comm.mat.env<-merge(env, comm, by.x = "TreeN", by.y = "TreeN")
comm.mat <- comm.mat.env %>% select(Forest, ANOC001:VOLL001) %>% 
  group_by(Forest) %>%
  summarise_all(funs(sum))

comm.mat<-as.data.frame(comm.mat)
rownames(comm.mat) <- comm.mat[,1]
comm.mat <- comm.mat[,-c(1)]

# nesting only
nestcomm <- read.csv("nestcomm.csv")
nestcomm.mat.env<-merge(env, nestcomm, by.x = "TreeN", by.y = "TreeN")
nestcomm.mat <- nestcomm.mat.env %>% select(Forest, ANON001:VOLL001)%>%
  group_by(Forest) %>% 
  summarise_all(funs(sum))

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
fcomm.mat <- fcomm.mat.env %>% select(Forest, ANON001:VOLL001)%>%
  group_by(Forest) %>% 
  summarise_all(funs(sum))


fcomm.mat<-as.data.frame(fcomm.mat)
rownames(fcomm.mat) <- fcomm.mat[,1]
fcomm.mat <- fcomm.mat[,-c(1)]


# remove species with 0 occurences
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
# be
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
dist.fn # primary f and primary n are very similar (BC =0.55), Secondary f and Secondary n less (BC =0.70)

#----------------------------------------------------------#
# Functional diversity calculations  -----
#----------------------------------------------------------## 

# Here we calculate functional diversity, i.e. Raos Q and CWMs of single traits.
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

# merge occurence with trait data to remove species without trait data
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
  select(Primary, Secondary)
com.mat2<-t(com.mat2)

## traits nesters
traits.nestcomm.mat<-merge(finaltraitmatrix, t.nestcomm.mat, by.x = 0, by.y= 0)
rownames(traits.nestcomm.mat) <- traits.nestcomm.mat[,1]
traits.nestcomm.mat <- traits.nestcomm.mat[,-c(1)]
#
nestcom.mat2<-traits.nestcomm.mat%>%
  select(Primary, Secondary)
nestcom.mat2<-t(nestcom.mat2)

## traits foragers
traits.fcomm.mat<-merge(finaltraitmatrix,t.fcomm.mat, by.x = 0, by.y= 0)
rownames(traits.fcomm.mat) <- traits.fcomm.mat[,1]
traits.fcomm.mat <- traits.fcomm.mat[,-c(1)]
#
fcom.mat2<-traits.fcomm.mat%>%
  select(Primary, Secondary)
fcom.mat2<-t(fcom.mat2)

# calculate functional dendrograms
multi.dis.c <- gowdis(traits.comm.mat[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
multi.dis.n <- gowdis(traits.nestcomm.mat[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
multi.dis.f <- gowdis(traits.fcomm.mat[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix

# calculate SES Rao Q for traits
mpd.comm.p <- ses.mpd(com.mat2, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.nest.p <- ses.mpd(nestcom.mat2, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.forager.p<-ses.mpd(fcom.mat2, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
mpd.comm.p$TreeN<-rownames(mpd.comm.p)
mpd.nest.p$TreeN<-rownames(mpd.nest.p)
mpd.forager.p$TreeN<-rownames(mpd.forager.p)
#
mpd.comm.p # both ns., primary forest slightly clustered, secondary forest rather overdispersed
mpd.nest.p # both ns. primary forest more clustered, secondary forest random (p=0.52)
mpd.forager.p # # both ns., primary forest slightly clustered, secondary forest rather overdispersed

## 1. whole commmunity CWMs
traits_comm<-traits.comm.mat %>%
  select(Spines:poly.index)

abun_comm<-traits.comm.mat %>%
  select(Primary, Secondary)

t.abun_comm<-t(abun_comm)

# CHECK: names equal?
rownames(traits_comm) == rownames(abun_comm)

# Calculate community weighted means for whole plot
comm <- dbFD(traits_comm, t(abun_comm),
             # different weight for categorical trait Sculpture
             w= c(1,0.3, 1, 1, 1, 1, 1, 1, 1, 1),
             corr = "cailliez",
             calc.FGR = F,
             calc.FRic= F,
             clust.type = "ward"
)

# Extract CWMs
CWM_comm.p <- as.data.frame(comm$CWM)

# Rename values
CWM_comm.p$TreeN <- rownames(CWM_comm.p)
rownames(CWM_comm.p) <- NULL
#
CWM_comm.p

#### 2. nester CWMs
traits_nest<-traits.nestcomm.mat %>%
  select(Spines:poly.index)

abun_nest<-traits.nestcomm.mat %>%
  select(Primary, Secondary)

t.abun_nest<-t(abun_nest)

# CHECK: names equal?
rownames(traits_nest) == colnames(t.abun_nest)

# Create functional diversity indices for each plotmethod
nest <- dbFD(traits_nest, t.abun_nest,
             # different weight for categorical trait Sculpture
             w= c(1,0.33, 1, 1, 1, 1, 1, 1, 1, 1),
             corr = "cailliez",
             calc.FGR = F, 
             calc.FRic = F,
             clust.type = "ward")

# Extract FD-Indices
CWM_nest.p <- as.data.frame(nest$CWM)
# 
CWM_nest.p$plot_type <- rownames(CWM_nest.p)
rownames(CWM_nest.p) <- NULL
CWM_nest.p


#### 3. forager CWMs
traits_fcomm<-traits.fcomm.mat %>%
  select(Spines:poly.index)

abun_fcomm<-traits.fcomm.mat %>%
  select(Primary, Secondary)
t.abun_fcomm<-t(abun_fcomm)


# CHECK: names equal?
rownames(traits_fcomm) == colnames(t.abun_fcomm)

# Create functional diversity indices for each plotmethod
fcomm <- dbFD(traits_fcomm, t.abun_fcomm,
              # different weight for categorical trait Sculpture
              w= c(1,0.3, 1, 1, 1, 1, 1, 1, 1, 1),
              corr = "cailliez",
              calc.FGR = F,
              calc.FRic= F,
              clust.type = "ward"
)
# Extract FD-Indices
CWM_foragers.p <- as.data.frame(fcomm$CWM)
# Rename values
CWM_foragers.p$plot_type <- rownames(CWM_foragers.p)
rownames(CWM_foragers.p) <- NULL
CWM_foragers.p

#----------------------------------------------------------#
# Phylogenetic diversity calculations  -----
#----------------------------------------------------------#

# Here we calculate phylogenetic diversity, i.e. Rao Q
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

# Phylogenetic SES RaoQ
phylo.matrix <- read.csv("newphylodist.csv", row.names = 1) #outgroups removed
phy.dis <- (sqrt(as.dist(phylo.matrix))) # SQRT transform as advised by Letten&Cornwell 2015 MEE

# Rao Q compared with a null model for SES
ses.mpd.phylo.full.p <- ses.mpd(t.abun_comm, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.nest.p <- ses.mpd(t.abun_nest, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.forager.p <- ses.mpd(t.abun_fcomm, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
ses.mpd.phylo.full.p$TreeN<-rownames(ses.mpd.phylo.full.p)
ses.mpd.phylo.nest.p$TreeN<-rownames(ses.mpd.phylo.nest.p)
ses.mpd.phylo.forager.p$TreeN<-rownames(ses.mpd.phylo.forager.p)
#
ses.mpd.phylo.full.p # both ns, primary random, secondary rather clustered
ses.mpd.phylo.nest.p # # both ns, primary slightly clustered, secondary random
ses.mpd.phylo.forager.p  #  both ns, primary random, secondary rather clustered


#----------------------------------------------------------#
# Species richness  -----
#----------------------------------------------------------#
forest_c<-as.data.frame(c("Primary", "Secondary"))
forest_c$Richness_all <- specnumber(com.mat2)
forest_c$Richness_forager<- specnumber(nestcom.mat2)
forest_c$Richness_nester<- specnumber(fcom.mat2)

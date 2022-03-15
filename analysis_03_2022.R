## Arboreal ant commmunities Wanang 2022- PH Script

# In this script, we compare the functional and phylogenetic diversity of two plots, one in primary and one secondary forest.
# We analyse the communities on two scales, the tree scale and the plot scale.
# The plot scale sums up all incidence across the two plots, while the tree scale treats each tree as independent communities.

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
    "stringr",
    "ggplot2",
    "ggpubr",
    "iNEXT"
    )

# install all packages
sapply(package_list, install.packages, character.only = TRUE)

# load all packages
sapply(package_list, library, character.only = TRUE)

# Citations
sapply(package_list, citation)


#----------------------------------------------------------#
# Trait raw data transformation -----
#----------------------------------------------------------#

# Import Trait Table: this one has already corrected size corrected traits
traits.raw =read.csv(file="traits.raw.csv", header=T)

## SELECT ONLY NON-INVASIVE Species. If you activate this line, all subsequent analysis include only non-invasive species
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
# Part 1: Functional Diversity Analysis: TREE SCALE -----
#----------------------------------------------------------#

## community data
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

##### Species richness accumulation curves
#Primary
prim<-subset(comm.mat.env, Forest=="Primary")
prim<-prim[,-c(1:14)]

# Secondary
sec<-subset(comm.mat.env, Forest=="Secondary")
sec<-sec[,-c(1:14)]

# Secondary Nesters
sec.n<-subset(nestcomm.mat.env, Forest=="Secondary")
sec.n<-sec.n[,-c(1:14)]

# Primary Nesters
prim.n<-subset(nestcomm.mat.env, Forest=="Primary")
prim.n<-prim.n[,-c(1:14)]

# Secondary Foragers
sec.f<-subset(fcomm.mat.env, Forest=="Secondary")
sec.f<-sec.f[,-c(1:14)]

# Primary Foragers
prim.f<-subset(fcomm.mat.env, Forest=="Primary")
prim.f<-prim.f[,-c(1:14)]

###
emptylist<-vector("list")

emptylist$allspPrimary <- (t(prim))
emptylist$allspSecondary <- (t(sec))
#
emptylist$nestPrimary <- t(prim.n)
emptylist$nestSecondary <- t(sec.n)
emptylist$foragerPrimary <- t(prim.f)
emptylist$foragerSecondary <- t(sec.f)
#
out.test <- iNEXT(emptylist, datatype="incidence_raw", endpoint=600)

ggiNEXT(out.test)+ 
  theme_bw(base_size = 18) 

gamma.plot <- (ggiNEXT(out.test)+ 
                 theme_bw(base_size = 18))+
  theme(legend.position="none") +
  scale_color_manual(values=c("#005824", "#990000","#41ae76", "#ef6548"))+
  scale_fill_manual(values=c("#238b45","#ef6548","#66c2a4","#fc8d59"))

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


# merge with trait data to remove species without trait data
t.comm.mat = t(comm.mat)
t.nestcomm.mat = t(nestcomm.mat)
t.fcomm.mat = t(fcomm.mat)

#
traits.comm.mat<-merge(finaltraitmatrix, t.comm.mat, by.x = 0, by.y= 0)
traits.nestcomm.mat<-merge(finaltraitmatrix, t.nestcomm.mat, by.x = 0, by.y= 0)
traits.fcomm.mat<-merge(finaltraitmatrix, t.fcomm.mat, by.x = 0, by.y= 0)


##### functional diversity (FDis, FEve, etc.)
### 1. for complete community
traits_comm<-traits.comm.mat[, c(1:11)]
rownames(traits_comm) <- traits_comm[, 1]
traits_comm <- traits_comm[, -c(1)]

abun_comm<-traits.comm.mat[, -c(2:11)]
rownames(abun_comm) <- abun_comm[, 1]
abun_comm <- abun_comm[, -c(1)]

# remove trees with 0 occurences
i <- (colSums(abun_comm) != 0) 
abun_comm <- abun_comm[,i] 


t.abun_comm<-t(abun_comm)

# SES MPD
multi.dis.c <- gowdis(traits_comm[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
mpd.full <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.full$TreeN<-rownames(mpd.full)


# CHECK: names equal?
rownames(traits_comm) == rownames(abun_comm)

# Create functional diversity indices for each plotmethod
comm <- dbFD(traits_comm, t(abun_comm),
             # different weight for categorical trait Sculpture
             w= c(1,0.33, 1, 1, 1, 1, 1, 1, 1, 1),
             corr = "cailliez",
             calc.FGR = F,
             calc.FRic= F,
             clust.type = "ward"
)

# Extract CWMs
CWM_comm <- as.data.frame(comm$CWM)
#
CWM_comm$TreeN <- rownames(CWM_comm)
rownames(CWM_comm) <- NULL

### 2. for Nesters
traits_nest<-traits.nestcomm.mat[, c(1:11)]
rownames(traits_nest) <- traits_nest[, 1]
traits_nest <- traits_nest[, -c(1)]

abun_nest<-traits.nestcomm.mat[, -c(2:11)]
rownames(abun_nest) <- abun_nest[, 1]
abun_nest <- abun_nest[, -c(1)]


# SES mpd nesters
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

multi.dis.n <- gowdis(traits_nest[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)


# CHECK: names equal?
rownames(traits_nest) == rownames(abun_nest)

# Create functional diversity indices for each tree
nest <- dbFD(traits_nest, t(abun_nest),
             corr = "cailliez",
             # different weight for categorical trait Sculpture
             w= c(1,0.33, 1, 1, 1, 1, 1, 1, 1, 1),
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

# SES MPD
multi.dis.f <- gowdis(traits_foragers[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
mpd.f <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.f$TreeN<-rownames(mpd.f)

# CHECK: names equal?
rownames(traits_foragers) == rownames(abun_foragers)

# Create functional diversity indices for each tree
foragers <- dbFD(traits_foragers, t(abun_foragers),
                 corr = "cailliez",
                 # different weight for categorical trait Sculpture
                 w= c(1,0.33, 1, 1, 1, 1, 1, 1, 1, 1),
                 calc.FGR = F, 
                 calc.FRic = F,
                 clust.type = "ward")

# Extract FD-Indices
CWM_foragers <- as.data.frame(foragers$CWM)
# Rename values
CWM_foragers$TreeN <- rownames(CWM_foragers)
rownames(CWM_foragers) <- NULL

### merge nesters, foragers, all and metadata
CWM_comm$Type<-"all"
CWM_nest$Type<-"nest"
CWM_foragers$Type<-"foragers"
#
fd_all <- full_join(CWM_comm, CWM_nest)
fd_all2<-full_join(fd_all, CWM_foragers)
#
fd_env<- full_join(env, fd_all2)

# Phylogenetic SES MPD
read.csv("newphylodist.csv", row.names = 1) -> phylo.matrix #outgroups removed manually
phy.dis <- (sqrt(as.dist(phylo.matrix))) # SQRT transform as advised by Letten&Cornwell 2015 MEE

#mpd compared with a null model for SES
ses.mpd.phylo.full <- ses.mpd(t.abun_comm, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.nest <- ses.mpd(t.abun_nest, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.forager <- ses.mpd(t.abun_forager, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
ses.mpd.phylo.full$TreeN<-rownames(ses.mpd.phylo.full)
ses.mpd.phylo.nest$TreeN<-rownames(ses.mpd.phylo.nest)
ses.mpd.phylo.forager$TreeN<-rownames(ses.mpd.phylo.forager)
#
pd_comm<-ses.mpd.phylo.full[,c(6,9)]
pd_nest<-ses.mpd.phylo.nest[,c(6,9)]
pd_forager<-ses.mpd.phylo.forager[,c(6,9)]
#
pd_comm$Type<-"all"
pd_nest$Type<-"nest"
pd_forager$Type<-"foragers"
#
pd_all <- full_join(pd_comm, pd_nest)
pd_all2<-full_join(pd_all, pd_forager)
#
pd_env<- full_join(env,pd_all2, by="TreeN")

# for traits
mpd.full$TreeN<-rownames(mpd.full)
mpd.nest$TreeN<-rownames(mpd.nest)
mpd.f$TreeN<-rownames(mpd.f)
#
mpd_comm<-mpd.full[,c(6,9)]
mpd_nest<-mpd.nest[,c(6,9)]
mpd_forager<-mpd.f[,c(6,9)]
#
mpd_comm$Type<-"all"
mpd_nest$Type<-"nest"
mpd_forager$Type<-"foragers"
#
mpd_all <- full_join(mpd_comm, mpd_nest)
mpd_all2<-full_join(mpd_all, mpd_forager)
#
mpd_env<- full_join(env,mpd_all2, by="TreeN")

#----------------------------------------------------------#
# FIGURES -----
#----------------------------------------------------------#

# NOTE: All figures have the statistical comparisons implemented as Kruskal-Wallis tests via ggpubr
# all plots shown here are tree scale analysis - the plot scale results (see Part 2) are only reported as tables

## Functional Diversity
mpd_env<-mpd_env%>% drop_na(Type)
pd_env<-pd_env%>% drop_na(Type)

mpd_p1<-ggplot(mpd_env, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Trait diversity") +
  geom_boxplot()+
  geom_abline(intercept = 1.8, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_abline(intercept = -1.6, slope = 0, color="red", 
              linetype="dashed", size=1)+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  #ylim(-3, 3)+
  ylab("SES Rao Q")+
  theme_bw()
mpd_p1

## Phylogenetic Diversity
mpd_p2<-ggplot(pd_env, aes(x=Type, y=mpd.obs.z, fill=Forest)) +
  ggtitle("Phylogenetic Diversity") +
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  geom_boxplot()+
  geom_abline(intercept = 1.6, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_abline(intercept = -1.8, slope = 0, color="red", 
              linetype="dashed", size=1)+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  #ylim(-3, 3)+
  ylab("SES Rao Q")+
  theme_bw()
mpd_p2

## Species Richness
all_richness<-comm.mat.env %>% select(TreeN, Forest, AntSpRichness.ALL, AntSpRichness.Nests)
all_richness$AntSpRichness.Forager<-all_richness$AntSpRichness.ALL-all_richness$AntSpRichness.Nests

full_richness<-all_richness[,-c(4,5)]
full_richness$Type <- "all"
colnames(full_richness)[colnames(full_richness) == "AntSpRichness.ALL"] <- "richness"

n_richness<-all_richness[,-c(3,5)]
n_richness$Type <- "nest"
colnames(n_richness)[colnames(n_richness) == "AntSpRichness.Nests"] <- "richness"

f_richness<-all_richness[,-c(3,4)]
f_richness$Type <- "foragers"
colnames(f_richness)[colnames(f_richness) == "AntSpRichness.Forager"] <- "richness"

all_r<-rbind(full_richness,n_richness, f_richness)

####
richness<-ggplot(all_r, aes(x=Type, y=richness, fill=Forest)) +
  ggtitle("Species richness") +
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  theme_bw()
richness

#### Community Weighted Means 
fd_env<-fd_env%>% drop_na(Type)

HL<-ggplot(fd_env, aes(x=Type, y=HL, fill=Forest)) +
  ggtitle("Head Length") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  theme_bw()
HL

poly.index<-ggplot(fd_env, aes(x=Type, y=poly.index, fill=Forest)) +
  ggtitle("Polymorphism index") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+         
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  theme_bw()
poly.index

Spines<-ggplot(fd_env, aes(x=Type, y=Spines, fill=Forest)) +
  ggtitle("Spines") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  theme_bw()
Spines

Scul<-ggplot(fd_env, aes(x=Type, y=Scul, fill=Forest)) +
  ggtitle("Sculpture") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  theme_bw()
Scul

Rel.HW<-ggplot(fd_env, aes(x=Type, y=Rel.HW, fill=Forest)) +
  ggtitle("Rel. Headwidth") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  theme_bw()
Rel.HW

Rel.CL<-ggplot(fd_env, aes(x=Type, y=Rel.CL, fill=Forest)) +
  ggtitle("Rel. Clypeus Length") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  theme_bw()
Rel.CL

Rel.ML<-ggplot(fd_env, aes(x=Type, y=Rel.ML, fill=Forest)) +
  ggtitle("Rel.Mandible Length") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  ylab("")+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  theme_bw()
Rel.ML

Rel.EP<-ggplot(fd_env, aes(x=Type, y=Rel.EP, fill=Forest)) +
  ggtitle("Rel. Eye Position") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  ylab("")+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  theme_bw()
Rel.EP

Rel.ES<-ggplot(fd_env, aes(x=Type, y=Rel.ES, fill=Forest)) +
  ggtitle("Rel. Eye Size") +
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+    
  ylab("")+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  theme_bw()
Rel.ES

Rel.LL<-ggplot(fd_env, aes(x=Type, y=Rel.LL, fill=Forest)) +
  ggtitle("Rel. Leg Length") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+         
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  theme_bw()
Rel.LL

## Plot 1
plot1<-ggarrange(richness, mpd_p1, mpd_p2,
                 labels = c("A", "B", "C"),
                 ncol = 3, nrow = 1, common.legend = TRUE, legend = "top"
)
plot1

## Plot 2
figure_traits <- ggarrange(HL, poly.index, Spines, Scul, Rel.HW, Rel.CL, Rel.ML, Rel.EP, Rel.ES,Rel.LL,
                           labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
                           ncol = 5, nrow = 2, common.legend = TRUE, legend = "top"
)
figure_traits


#----------------------------------------------------------#
# Part 2: Functional Diversity Analysis: PLOT SCALE -----
#----------------------------------------------------------#

## community data
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
# remove species with 0 occurence
i <- (colSums(nestcomm.mat) != 0) 
nestcomm.mat <- nestcomm.mat[,i] 


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
# remove species with 0 occurence
i <- (colSums(fcomm.mat) != 0) 
fcomm.mat <- fcomm.mat[,i] 


## Species overlap
# prim vs sec
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

nester_forager<-full_join(fcomm.mat2, nestcomm.mat2)

nester_forager[is.na(nester_forager)] <- 0
row.names(nester_forager) <- c("Primary_f","Secondary_f","Primary_n", "Secondary_n")

# Bray 
dist.fn <- vegdist(nester_forager, method = "bray")
dist.fn # primary f and primary n are very similar (BC =0.55), Secondary f and Secondary n less (BC =0.70)

# merge with trait data to remove species without trait data
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

#functional dendrogram
multi.dis.c <- gowdis(traits.comm.mat[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
multi.dis.n <- gowdis(traits.nestcomm.mat[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
multi.dis.f <- gowdis(traits.fcomm.mat[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix

# SES MPD Traits
mpd.comm <- ses.mpd(com.mat2, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.nest <- ses.mpd(nestcom.mat2, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.forager<-ses.mpd(fcom.mat2, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
mpd.comm$TreeN<-rownames(mpd.comm)
mpd.nest$TreeN<-rownames(mpd.nest)
mpd.forager$TreeN<-rownames(mpd.forager)

#
mpd.comm # both ns., primary forest slightly clustered, secondary forest rather overdispersed
mpd.nest # both ns. primary forest more clustered, secondary forest random (p=0.52)
mpd.forager # # both ns., primary forest slightly clustered, secondary forest rather overdispersed

## functional diversity for ALL
traits_comm<-traits.comm.mat %>%
  select(Spines:poly.index)

abun_comm<-traits.comm.mat %>%
  select(Primary, Secondary)

t.abun_comm<-t(abun_comm)

# CHECK: names equal?
rownames(traits_comm) == rownames(abun_comm)

# Create functional diversity indices for each plotmethod
comm <- dbFD(traits_comm, t(abun_comm),
             # different weight for categorical trait Sculpture
             w= c(1,0.3, 1, 1, 1, 1, 1, 1, 1, 1),
             corr = "cailliez",
             calc.FGR = F,
             calc.FRic= F,
             clust.type = "ward"
)

# Extract FD-Indices
CWM_comm <- as.data.frame(comm$CWM)

# Rename values
CWM_comm$TreeN <- rownames(CWM_comm)
rownames(CWM_comm) <- NULL

### functional diversity FOR NESTS
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
CWM_nest <- as.data.frame(nest$CWM)
# 
CWM_nest$plot_type <- rownames(CWM_nest)
rownames(CWM_nest) <- NULL

### functional diversity FOR FORAGERS
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
CWM_foragers <- as.data.frame(fcomm$CWM)
# Rename values
CWM_foragers$plot_type <- rownames(CWM_foragers)
rownames(CWM_foragers) <- NULL

###### Phylogenetic PLOT SCALE

# Phylogenetic SES MPD
read.csv("newphylodist.csv", row.names = 1) -> phylo.matrix #outgroups removed
phy.dis <- (sqrt(as.dist(phylo.matrix))) # SQRT transform as advised by Letten&Cornwell 2015 MEE

#mpd compared with a null model for SES
ses.mpd.phylo.full <- ses.mpd(t.abun_comm, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.nest <- ses.mpd(t.abun_nest, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.forager <- ses.mpd(t.abun_fcomm, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
ses.mpd.phylo.full$TreeN<-rownames(ses.mpd.phylo.full)
ses.mpd.phylo.nest$TreeN<-rownames(ses.mpd.phylo.nest)
ses.mpd.phylo.forager$TreeN<-rownames(ses.mpd.phylo.forager)
#
ses.mpd.phylo.full # both ns, primary ramndom, secondary rather clustered
ses.mpd.phylo.nest # # both ns, primary slightly clustered, secondary random (opposite to full and forager
#but fits with Nichola chapter, correlation of observed PD and SES of PD not always the case for two datapoints)
ses.mpd.phylo.forager  # # both ns, primary ramndom, secondary rather clustered

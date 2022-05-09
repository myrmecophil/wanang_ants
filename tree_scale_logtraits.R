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
    "ggpubr")

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
## If you activate this line, all subsequent analyses include only native species
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


## Normality: HeadLength -> much better
qqnorm(traits_all$HeadLength)
qqline(traits_all$HeadLength)

qqnorm(log(traits_all$HeadLength))
qqline(log(traits_all$HeadLength))
traits_all$HeadLength<-log(traits_all$HeadLength)

## Normality: Spinosity -> cant transform
qqnorm(traits_all$Spinosity)
qqline(traits_all$Spinosity)

## Normality: Sculpturing -> a bit better not much
qqnorm(traits_all$Sculpturing)
qqline(traits_all$Sculpturing)
hist(traits_all$Sculpturing)
qqnorm(log(traits_all$Sculpturing))
qqline(log(traits_all$Sculpturing))
hist(log(traits_all$Sculpturing))


## Normality: LegLength -> much better
qqnorm(traits_all$LegLength_rel)
qqline(traits_all$LegLength_rel)

qqnorm(log(traits_all$LegLength_rel))
qqline(log(traits_all$LegLength_rel))

traits_all$LegLength_rel<-log(traits_all$LegLength_rel)


## Normality: Headwidth -> a bit?
qqnorm(traits_all$HeadWidth_rel)
qqline(traits_all$HeadWidth_rel)

qqnorm(log(traits_all$HeadWidth_rel))
qqline(log(traits_all$HeadWidth_rel))

traits_all$HeadWidth_rel<-log(traits_all$HeadWidth_rel)


## Normality: Clypeus -> worse
qqnorm(traits_all$ClypeusLength_rel)
qqline(traits_all$ClypeusLength_rel)

qqnorm(log(traits_all$ClypeusLength_rel))
qqline(log(traits_all$ClypeusLength_rel))

## Normality: Mandible -> same
qqnorm(traits_all$MandibleLength_rel)
qqline(traits_all$MandibleLength_rel)

qqnorm(log(traits_all$MandibleLength_rel))
qqline(log(traits_all$MandibleLength_rel))

traits_all$MandibleLength_rel<-log(traits_all$MandibleLength_rel)

## Normality: Interocular distance -> worse
qqnorm(traits_all$InterocularDistance_rel)
qqline(traits_all$InterocularDistance_rel)

qqnorm(log(traits_all$InterocularDistance_rel))
qqline(log(traits_all$InterocularDistance_rel))


## Normality: Interocular distance -> a bit better
qqnorm(traits_all$EyeSize_rel)
qqline(traits_all$EyeSize_rel)

qqnorm(log(traits_all$EyeSize_rel))
qqline(log(traits_all$EyeSize_rel))

traits_all$EyeSize_rel<-log(traits_all$EyeSize_rel)

## Normality: polymorphism -> a bit better
qqnorm(finaltraitmatrix$poly.index)
qqline(finaltraitmatrix$poly.index)

qqnorm(log(finaltraitmatrix$poly.index))
qqline(log(finaltraitmatrix$poly.index))

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

# remove trees with 0 occurences
i <- (colSums(abun_comm) != 0) 
abun_comm <- abun_comm[,i] 


t.abun_comm<-t(abun_comm)

# calcuate SES Rao Q
multi.dis.c <- gowdis(traits_comm[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)

# CHECK: names equal?
rownames(traits_comm) == rownames(abun_comm)

# Create CWMs
c <- dbFD(traits_comm, t(abun_comm),
          # different weight for categorical trait Sculpture
          w= c(1,0.33, 1, 1, 1, 1, 1, 1, 1, 1),
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

# SES Rao Q
multi.dis.f <- gowdis(traits_foragers[,-c(11,12)], w=c(1,0.3,1,1,1,1,1,1,1,1)) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)

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
mpd_comm<-mpd.comm[,c(6,9)]
mpd_nest<-mpd.nest[,c(6,9)]
mpd_forager<-mpd.foragers[,c(6,9)]
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
# Phylogenetic diversity calculations  -----
#----------------------------------------------------------#

# Here we calculate phylogenetic diversity, i.e. Rao Q
# We do this three times: 1. complete community 2. for nesting community 3. for foraging community

# Phylogenetic SES Rao Q
phylo.matrix <- read.csv("newphylodist.csv", row.names = 1)  # outgroups removed manually
phy.dis <- (sqrt(as.dist(phylo.matrix))) # SQRT transform as advised by Letten&Cornwell 2015 MEE

#Rao Q SES
ses.mpd.phylo.full <- ses.mpd(t.abun_comm, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.nest <- ses.mpd(t.abun_nest, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
ses.mpd.phylo.forager <- ses.mpd(t.abun_forager, phy.dis, null.model = "taxa.labels", abundance.weighted = T, runs = 999)
#
ses.mpd.phylo.full$TreeN<-rownames(ses.mpd.phylo.full)
ses.mpd.phylo.nest$TreeN<-rownames(ses.mpd.phylo.nest)
ses.mpd.phylo.forager$TreeN<-rownames(ses.mpd.phylo.forager)

## merge results into table to plot them later on
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
forest_c<-env %>% select(TreeN, Forest)
forest_c$Richness <- specnumber(t.abuncomm2)
forest_c$Type<-"all"
#
forest_n<-env %>% select(TreeN, Forest)
forest_n$Richness <- specnumber(t.nest_comm2)
forest_n$Type<-"nest"
#
forest_f<-env %>% select(TreeN, Forest)
forest_f$Richness <- specnumber(t.f_comm2)
forest_f$Type<-"forager"

# get all results into one table
all_r<-rbind(forest_c,forest_n, forest_f)

#----------------------------------------------------------#
# FIGURES -----
#----------------------------------------------------------#

# NOTE: All figures have the statistical comparisons implemented as Kruskal-Wallis tests via ggpubr
# all plots shown here are tree scale analysis - the plot scale results (see Part 2) are only reported as tables

#### Functional Diversity
mpd_env<-mpd_env%>% drop_na(Type)
pd_env<-pd_env%>% drop_na(Type)

mpd_p1<-ggplot(mpd_env, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Trait diversity") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  #ylim(-3, 3)+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))

mpd_p1

#### Phylogenetic Diversity
mpd_p2<-ggplot(pd_env, aes(x=Type, y=mpd.obs.z, fill=Forest)) +
  ggtitle("Phylogenetic Diversity") +
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+   
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  #ylim(-3, 3)+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))

mpd_p2

#### Richness
richness<-ggplot(all_r, aes(x=Type, y=Richness, fill=Forest)) +
  ggtitle("Species richness") +
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "forager" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))

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
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))

HL

poly.index<-ggplot(fd_env, aes(x=Type, y=poly.index, fill=Forest)) +
  ggtitle("Polymorphism index") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+         
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))

poly.index

Spines<-ggplot(fd_env, aes(x=Type, y=Spines, fill=Forest)) +
  ggtitle("Spines") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
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
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
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
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
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
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
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
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
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
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
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
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Rel.ES

Rel.LL<-ggplot(fd_env, aes(x=Type, y=Rel.LL, fill=Forest)) +
  ggtitle("Rel. Leg Length") +
  geom_boxplot()+
  scale_fill_manual(values=c("darkgrey", "lightgrey"))+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+         
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
Rel.LL

## Plot 1
plot1 <- ggarrange(richness, mpd_p1, mpd_p2,
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


###### REMOVE LATER: Export figures ######
library(svglite)

# Figure 1: richness + trait diversity + phylogenetic diversity
# as tiff
tiff("figure_1.tiff", units ="in", width = 10, height = 4, res = 150)
plot1
dev.off()

# or as .svg (vector image)
#ggsave(file="figure_1.svg", plot=plot1, width=10, height=4)

# Figure 2: CWMs
# as tiff
tiff("figure_2.tiff", units = "in", width = 16, height = 8, res = 150)
figure_traits
dev.off()

# or as .svg (vector image)
ggsave(file="figure_2.svg", plot=figure_traits, width=16, height=8)


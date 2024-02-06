## Arboreal ant communities Wanang 2022- Tree scale Supplement Figure Script

# In this script, we compare the functional and phylogenetic diversity of two plots, one in primary and one secondary forest.

# This part of the script is solely for the Figures S1 and S3 provided in the Supplement

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

#----------------------------------------------------------#
### Data preparation 
#----------------------------------------------------------#

# set seed for randomization
set.seed(123)

# load R data from tree scale scripts
load("tree-scale.Rdata")

#----------------------------------------------------------#
#  Figure S1: Trait correlations -----
#----------------------------------------------------------#

# Trait correlations
mydata.cor <- cor(finaltraitmatrix.p, method = c("spearman"), use = "pairwise.complete.obs")
FigS1<-corrplot(mydata.cor,
         method = "color", type = "lower", order = "AOE", diag = F,
         tl.col = "black", outline = T, addCoef.col = "black", number.cex = 0.8,
         tl.cex = 1.1, cl.cex = 0.9)
FigS1


#----------------------------------------------------------#
#  Figure S3: Rao Q of single traits -----
#----------------------------------------------------------#

# calcuate SES Rao Q for each trait separately

# 1. spines
# whole comm
multi.dis.c <- gowdis(traits_comm[,1, FALSE]) # gower distance matrix
mpd.comm1 <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm1$TreeN<-rownames(mpd.comm1)
# nesters 
multi.dis.n <- gowdis(traits_nest[,1, FALSE]) # gower distance matrix
mpd.nest1<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest1$TreeN <-rownames(mpd.nest1)
# visitors
multi.dis.f <- gowdis(traits_foragers[,1, FALSE]) # gower distance matrix
mpd.foragers1 <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers1$TreeN<-rownames(mpd.foragers1)
# merge
mpd.comm1$Type<-"all"
mpd.nest1$Type<-"nest"
mpd.foragers1$Type<-"foragers"
#
fd_all <- full_join(mpd.comm1, mpd.nest1)
fd_all2<-full_join(fd_all, mpd.foragers1)
#
env1 <- full_join(env, fd_all2)


# 2 sculpture
# whole comm
multi.dis.c <- gowdis(traits_comm[,2, FALSE]) # gower distance matrix
mpd.comm2 <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm2$TreeN<-rownames(mpd.comm2)
# nesters 
multi.dis.n <- gowdis(traits_nest[,2, FALSE]) # gower distance matrix
mpd.nest2<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest2$TreeN <-rownames(mpd.nest2)
# visitors
multi.dis.f <- gowdis(traits_foragers[,2, FALSE]) # gower distance matrix
mpd.foragers2 <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers2$TreeN<-rownames(mpd.foragers2)
# merge
mpd.comm2$Type<-"all"
mpd.nest2$Type<-"nest"
mpd.foragers2$Type<-"foragers"
#
fd_all <- full_join(mpd.comm2, mpd.nest2)
fd_all2<-full_join(fd_all, mpd.foragers2)
#
env2 <- full_join(env, fd_all2)

# 3  rel hw
# whole comm
multi.dis.c <- gowdis(traits_comm[,3, FALSE]) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)
# nesters 
multi.dis.n <- gowdis(traits_nest[,3, FALSE]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)
# visitors
multi.dis.f <- gowdis(traits_foragers[,3, FALSE]) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)
# merge
mpd.comm$Type<-"all"
mpd.nest$Type<-"nest"
mpd.foragers$Type<-"foragers"
#
fd_all <- full_join(mpd.comm, mpd.nest)
fd_all2<-full_join(fd_all, mpd.foragers)
#
env3 <- full_join(env, fd_all2)

# 4 rel cl
# whole comm
multi.dis.c <- gowdis(traits_comm[,4, FALSE]) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)
# nesters 
multi.dis.n <- gowdis(traits_nest[,4, FALSE]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)
# visitors
multi.dis.f <- gowdis(traits_foragers[,4, FALSE]) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)
# merge
mpd.comm$Type<-"all"
mpd.nest$Type<-"nest"
mpd.foragers$Type<-"foragers"
#
fd_all <- full_join(mpd.comm, mpd.nest)
fd_all2<-full_join(fd_all, mpd.foragers)
#
env4 <- full_join(env, fd_all2)


# 5 rel ml
# whole comm
multi.dis.c <- gowdis(traits_comm[,5, FALSE]) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)
# nesters 
multi.dis.n <- gowdis(traits_nest[,5, FALSE]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)
# visitors
multi.dis.f <- gowdis(traits_foragers[,5, FALSE]) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)
# merge
mpd.comm$Type<-"all"
mpd.nest$Type<-"nest"
mpd.foragers$Type<-"foragers"
#
fd_all <- full_join(mpd.comm, mpd.nest)
fd_all2<-full_join(fd_all, mpd.foragers)
#
env5 <- full_join(env, fd_all2)

# 6 rel ep
# whole comm
multi.dis.c <- gowdis(traits_comm[,6, FALSE]) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)
# nesters 
multi.dis.n <- gowdis(traits_nest[,6, FALSE]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)
# visitors
multi.dis.f <- gowdis(traits_foragers[,6, FALSE]) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)
# merge
mpd.comm$Type<-"all"
mpd.nest$Type<-"nest"
mpd.foragers$Type<-"foragers"
#
fd_all <- full_join(mpd.comm, mpd.nest)
fd_all2<-full_join(fd_all, mpd.foragers)
#
env6 <- full_join(env, fd_all2)

# 7 rel es
# whole comm
multi.dis.c <- gowdis(traits_comm[,7, FALSE]) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)
# nesters 
multi.dis.n <- gowdis(traits_nest[,7, FALSE]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)
# visitors
multi.dis.f <- gowdis(traits_foragers[,7, FALSE]) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)
# merge
mpd.comm$Type<-"all"
mpd.nest$Type<-"nest"
mpd.foragers$Type<-"foragers"
#
fd_all <- full_join(mpd.comm, mpd.nest)
fd_all2<-full_join(fd_all, mpd.foragers)
#
env7 <- full_join(env, fd_all2)


# 8 ll
# whole comm
multi.dis.c <- gowdis(traits_comm[,8, FALSE]) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)
# nesters 
multi.dis.n <- gowdis(traits_nest[,8, FALSE]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)
# visitors
multi.dis.f <- gowdis(traits_foragers[,8, FALSE]) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)
# merge
mpd.comm$Type<-"all"
mpd.nest$Type<-"nest"
mpd.foragers$Type<-"foragers"
#
fd_all <- full_join(mpd.comm, mpd.nest)
fd_all2<-full_join(fd_all, mpd.foragers)
#
env8 <- full_join(env, fd_all2)

# 9 hl
# whole comm
multi.dis.c <- gowdis(traits_comm[,9, FALSE]) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)
# nesters 
multi.dis.n <- gowdis(traits_nest[,9, FALSE]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)
# visitors
multi.dis.f <- gowdis(traits_foragers[,9, FALSE]) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)
# merge
mpd.comm$Type<-"all"
mpd.nest$Type<-"nest"
mpd.foragers$Type<-"foragers"
#
fd_all <- full_join(mpd.comm, mpd.nest)
fd_all2<-full_join(fd_all, mpd.foragers)
#
env9 <- full_join(env, fd_all2)

# 10 poly index
# whole comm
multi.dis.c <- gowdis(traits_comm[,10, FALSE]) # gower distance matrix
mpd.comm <- ses.mpd(t.abun_comm, multi.dis.c, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.comm$TreeN<-rownames(mpd.comm)
# nesters 
multi.dis.n <- gowdis(traits_nest[,10, FALSE]) # gower distance matrix
mpd.nest<- ses.mpd(t.abun_nest, multi.dis.n, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #nest abundance weighted
mpd.nest$TreeN <-rownames(mpd.nest)
# visitors
multi.dis.f <- gowdis(traits_foragers[,10, FALSE]) # gower distance matrix
mpd.foragers <- ses.mpd(t.abun_forager, multi.dis.f, null.model = "taxa.labels", abundance.weighted = T, runs = 999) #full abundance weighted
mpd.foragers$TreeN<-rownames(mpd.foragers)
# merge
mpd.comm$Type<-"all"
mpd.nest$Type<-"nest"
mpd.foragers$Type<-"foragers"
#
fd_all <- full_join(mpd.comm, mpd.nest)
fd_all2<-full_join(fd_all, mpd.foragers)
#
env10 <- full_join(env, fd_all2)

# plot single RaoQ traits

env1<-env1%>% drop_na(Type)
env2<-env2%>% drop_na(Type)
env3<-env3%>% drop_na(Type)
env4<-env4%>% drop_na(Type)
env5<-env5%>% drop_na(Type)
env6<-env6%>% drop_na(Type)
env7<-env7%>% drop_na(Type)
env8<-env8%>% drop_na(Type)
env9<-env9%>% drop_na(Type)
env10<-env10%>% drop_na(Type)

# 1 spines
rao.spines<-ggplot(env1, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Spines") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.spines

# 2 Scul    
rao.scul<-ggplot(env2, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Sculpture") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.scul

# 3 Rel.HW    
rao.hw<-ggplot(env3, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Head width") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.hw

# 4 Rel.CL    
rao.cl<-ggplot(env4, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Clypeus length") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.cl

# 5 Rel.ML    
rao.ml<-ggplot(env5, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Mandible length") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.ml


# 6  Rel.EP    
rao.ep<-ggplot(env6, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Eye position") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.ep

# 7  Rel.ES   
rao.es<-ggplot(env7, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Eye size") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.es

# 8  Rel.LL   
rao.ll<-ggplot(env8, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Leg length") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.ll

# 9  HL   
rao.hl<-ggplot(env9, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Head length") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.hl

# 10  poly   
rao.poly<-ggplot(env10, aes(x=Type, y=mpd.obs.z, fill=Forest))+
  ggtitle("Polymorphism") +
  geom_abline(intercept = 0, slope = 0, color="red", 
              linetype="dashed", size=1)+
  geom_boxplot()+
  scale_x_discrete(labels=c("all" = "All", "foragers" = "Visitors",
                            "nest" = "Nesters"))+           
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  stat_compare_means(method = "kruskal.test", paired=F, label = "p.signif")+
  ylab("SES Rao Q")+
  xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=12))
rao.poly

## Combine into final figure
rao.single <- ggarrange(rao.hl, rao.poly, rao.spines, rao.scul, rao.hw, rao.cl,rao.ml, rao.ep, rao.es,rao.ll,
                        labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
                        ncol = 5, nrow = 2, common.legend = TRUE, legend = "top"
)
rao.single

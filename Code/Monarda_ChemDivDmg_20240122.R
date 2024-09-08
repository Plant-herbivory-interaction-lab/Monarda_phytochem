library(emmeans)
library(car)
library(smatr)
library(MuMIn)
library(aods3)
library(glmmTMB)
library(betareg)
library(RVAideMemoire)
library(DHARMa)
library(patchwork)
library(lubridate)
library(hillR)
library(vegan)
library(viridis)
library(ggeffects)
library(splines)
library(iNEXT)
library(performance)
library(bbmle)
library(MASS)
library(tidyverse)
library(ape)
library(ggsci)
library(geodist)

# load data ####
trait <- read_csv("Data/Monarda_Garden_Traits_2019.csv")
dmg19 <- read_csv("Data/WIMT_Monarda_HgtDmg_09192018.csv")
sites <- read_csv("Data/WIMT_Sites_2019.csv")

head(trait)
head(sites)
head(dmg19)

# productivity ####
sitesnpp <- sites %>% filter(Site!='PET',Site!='BRP',Site!='GRE',Site!='BRS',Site!='FAV',Site!='TMP')
sitesnpp$npp <-  6116 * (1 - exp(-6.05 * 10^-5 * sitesnpp$map12)) 
sitesnppmeans <- sitesnpp %>% group_by(Region) %>% 
  summarize(mat1m=mean(mat1,na.rm=T), map12m=mean(map12, na.rm=T), mat1sd=sd(mat1,na.rm=T), map12sd=sd(map12, na.rm=T), 
            nppm=mean(npp, na.rm=T), nppsd=sd(npp, na.rm=T),
            orgmatm=mean(orgmat,na.rm=T),orgmatsd=sd(orgmat,na.rm=T),
            cecm=mean(cec, na.rm=T), cecsd=sd(cec, na.rm=T),
            nitm=mean(nit, na.rm=T), nitsd=sd(nit, na.rm=T),
            Pm=mean(P, na.rm=T), Psd=sd(P,na.rm=T))

sitesnpplong <- sitesnpp %>% select(Region, mat1, map12, npp, nit, orgmat, P) %>% 
  pivot_longer(cols = 2:7, names_to = "variable", values_to = "values")

man1 <- lm(cbind(mat1, map12, npp, orgmat, cec, nit, P) ~ Region, data=sitesnpp)
Manova(man1)
summary(Anova(man1), univariate=TRUE, multivariate=FALSE,
        p.adjust.method=TRUE)

envplot <- ggplot(sitesnpplong , aes(x=Region, y=values, color=Region, fill=Region)) + 
  geom_point()+
  geom_boxplot(lwd=1, outlier.shape=NA) +
  #geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=.2),  size=2 ) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1)+
  scale_fill_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.25)+
  scale_x_discrete("Region")+
  scale_y_continuous("")+
  theme_bw(base_size=20)+
  facet_wrap(~variable, scales="free")
envplot

ggsave('envplot.tiff',envplot, width=8, height=5, units="in", dpi=600, path="Outputs")

NPPmodel(831, 7.567833, "schuur")
NPPmodel(836, 7.786500)
NPPmodel(849, 7.526667)
NPPmodel(807, 7.454167)
NPPmodel(849, 7.526667)
NPPmodel(885, 8.381666)

NPPmodel(400, 5.071000)
NPPmodel(374, 6.308333)
NPPmodel(366, 7.556333)
NPPmodel(380, 6.839833)
NPPmodel(462, 5.382333)
NPPmodel(363, 7.127667)


# Q1: chem diversity between regions and chemotypes ####
chem <- trait %>% select(a_thujene:caryophyllene) #%>% sqrt()
#chem1 <- trait %>% select(Region,Site,Plant,Rep,Biomass_g,total,a_thujene:caryophyllene, -Rep, -nerol, -geraniol, -geranial) %>% 
 # mutate(Chemo=if_else(thymol>carvacrol, "T","C")) %>% 
  #group_by(Region, Site, Plant, Chemo) %>% 
  #summarise_all(mean) %>% 
  #filter(thymol>0|carvacrol>0) %>% na.omit()

chem1 <- trait %>% select(Region,Site,Plant,Rep,Biomass_g,total,a_thujene:caryophyllene, -nerol, -geraniol, -geranial) %>% 
  mutate(Chemo=if_else(thymol>carvacrol, "T","C")) %>% 
  group_by(Region, Site, Plant, Chemo,Rep) %>% 
  filter(thymol>0|carvacrol>0) %>% na.omit()
chemeuclid <- vegdist(sqrt(chem1[7:31]), method='bray')

### table of number of plants ####
head(chem1)
chem1 %>% group_by(Region,Site,Plant) %>% count() %>% 
  group_by(Region,Site) %>%
  summarize(nmin=min(n), nmax=max(n), sumn=sum(n), meann=mean(n))

chem1 %>% group_by(Region,Site,Plant) %>% count() %>% 
  group_by(Region,Site) %>% count()

chem1prop <- chem1 %>% group_by(Region,Site,Chemo) %>% count() %>% 
  pivot_wider(names_from = "Chemo", values_from = 'n', values_fill = 0) %>% 
  mutate(Total=C+T, C_perc = C/Total, T_perc=T/Total)

proptest <- glm(cbind(C,T) ~ Region, data = chem1prop, family=binomial)
Anova(proptest)
emmeans(proptest, pairwise~Region, type="response")


## PERMANOVA ####
ad1 <- adonis2(chemeuclid ~ Region+Chemo+Region:Chemo, data=chem1, by="term")
ad1
anova(ad1, by="margin")
pairwise.perm.manova(chemeuclid, paste(chem1$Region,chem1$Chemo))


#cmatrix <- log(chem1[7:31]+.001) %>% as.matrix()
cmatrix <- sqrt(chem1[7:31]) %>% as.matrix()

chem.mds <- metaMDS(cmatrix, distance = 'bray')
plot(chem.mds)
scores(chem.mds)



man1 <- manova(cmatrix ~ Region * Chemo, data=chem1)
Anova(man1)


## pcoa ####
chem.D <- vegdist(sqrt(chem1[7:31]), "robust.aitchison")
chem.pcoa <- pcoa(chem.D)
chem.pcoa$values
chem.pcoa$vectors
biplot(chem.pcoa, sqrt(chem1[7:31]))


## pca ####
spc <- prcomp(cmatrix, scale=F, center=TRUE)
summary(spc)
spc
plot(spc)
biplot(spc, pc.biplot=TRUE)

#chem1$pc1 <-  predict(spc)[,1]*-1
#chem1$pc2 <-  predict(spc)[,2]


pcrot <- as.data.frame(spc$rotation[,1:2], stringsAsFactors=TRUE)
pcrot$Envvar <- c("","","","","","","","","","","","","","","","","","","","","","thymol","carvacrol","","")
pcrot$PC1s <- pcrot$PC1*-1 #scale loadings for plotting as vector arrows
pcrot$PC2s <- pcrot$PC2 #scale loadings for plotting as vector arrows

## test differences in PC axes between regions
an1 <- glmmTMB(pc1~Region*Chemo + (1|Site), data=chem1)
Anova(an1)
an2 <- glmmTMB(pc2 ~ Region*Chemo + (1|Site), data=chem1)
Anova(an2)


## PCA FIGURE ####

pcfig <- ggplot() + xlab("PC1 (59.4%)")+ylab("PC2 (8.5%)")+ 
  stat_ellipse(data=chem1, geom="polygon" ,aes(x=pc1,y=pc2, color=Region, fill=Region), lwd=1, alpha=.25)+
  geom_point(data=chem1, aes(x=pc1, y=pc2, color=Region), size=3,alpha=.5) +
  #scale_fill_manual(values=c("darkgoldenrod1","dodgerblue"), labels=c("Low productivity","High productivity")) +
  #scale_color_manual(values=c("darkgoldenrod1","dodgerblue"), labels=c("Low productivity","High productivity")) +
  geom_segment(data=pcrot, aes(x=0, y=0,xend=PC1s*8,yend=PC2s*8), arrow=arrow(length=unit(0.03, "npc")), color='black', lwd=1)+
  geom_text(data=pcrot, aes(x=PC1s*8.5, y=PC2s*8.5, label=Envvar), fontface="bold", color="BLACK", size=8)+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.8, direction=-1)+
  theme_bw(base_size = 16) 
 

pcfig
#ggsave('pcfig.tiff',pcfig, width=7, height=6, units="in", dpi=600)

## pivot chemicals longer for analysis ####
chem2 <- trait %>% select(Region,Site,Plant,Rep,Chemo,Biomass_g,total,sla,a_thujene:caryophyllene,-geraniol,-geranial,-nerol) %>% 
  filter(thymol>0|carvacrol>0) %>%  na.omit() %>% 
  pivot_longer(cols=a_thujene:caryophyllene, names_to = 'chemical', values_to = 'conc')

chem2 <- chem1 %>%  
  pivot_longer(cols=a_thujene:caryophyllene, names_to = 'chemical', values_to = 'conc')

chemtab2 <- chem2 %>% group_by(chemical) %>% summarize(conc=mean(conc, na.rm=T))

glmm1 <- glmmTMB(conc ~ chemical + (1|Site/Plant), data=chem2, family=ziGamma("log"), zi=~1)
glmm1a <- glmmTMB(conc+.001 ~ Region * chemical + (1|Site/Plant), data=chem2, family=Gamma("log"))

glmm1b <- glmmTMB(conc ~ Region + (1|Site/Plant), data=chem2, family=ziGamma("log"), zi=~1)
glmm1c <- glmmTMB(conc ~ Region + (1|Site/Plant) + (0+Region|chemical), data=chem2, family=ziGamma("log"), zi=~1)
glmm1d <- glmmTMB(conc ~ Region + (1|Site/Plant) + (0+Region|chemical), data=chem2, family=ziGamma("log"), zi=~chemical)
glmm1e <- glmmTMB(conc ~ Region + (1|Site/Plant) + (0+Region|chemical), data=chem2, family=ziGamma("log"), zi=~1, disp=~chemical)
#glmm1 <- glmmTMB(log(conc+.001) ~ Region * Chemo * chemical + (1|Site/Plant/Rep), data=chem2, disp=~chemical)
#glmm1a <- glmmTMB(log(conc+.001) ~ Region * chemical + (1|Site/Plant/Rep), data=chem2, disp=~chemical)
Anova(glmm1b)

anova(glmm1,glmm1a)
anova(glmm1b,glmm1c,glmm1d,glmm1e)

simulateResiduals(fittedModel = glmm1, plot=T)

emmeans(glmm1,  ~ Region|chemical, type="response")
emmeans(glmm1b, pairwise ~ Region|chemical, type="response")

### Fig 1 -- create dataset for plotting ####

chem3 <- chem2 %>% group_by(chemical) %>% summarize(conc_mean=mean(conc)) %>% arrange(desc(conc_mean)) 
head(chem2)

chem2plot <- full_join(chem2,chem3, by="chemical") %>% 
  mutate(compound=if_else(conc_mean>.1, chemical, "other"), compound=as.factor(compound)) %>% 
  ungroup() %>% 
  mutate(compound=fct_relevel(compound, c("other","a_thujene","myrcene","octen_3_ol","a_terpinene","gama_terpinene",
                                          "p_cymene","thymoquinone","carvacrol","thymol"))) %>% 
  group_by(Region,Chemo,compound) %>% summarize(conc=mean(conc, na.rm=T))
    


fig1 <- ggplot(chem2plot , aes(fill=compound, x=paste(Region,Chemo,sep="-"), y=sqrt(conc))) + 
  geom_bar(position="stack",stat="identity")+
  #facet_wrap(~Season) + 
  #scale_fill_viridis(discrete=T,labels=c('Other','α-Thujene','Myrcene','Octen-3-ol','α-Terpinene', 'γ-Terpinene',
   #                                      'p-Cymene','Thymoquinone','Carvacrol','Thymol'), name="Compound")+
  scale_fill_jco(labels=c('Other','α-Thujene','Myrcene','Octen-3-ol','α-Terpinene', 'γ-Terpinene',
                                         'p-Cymene','Thymoquinone','Carvacrol','Thymol'), name="Compound")+
  #ggtitle("B")+
  xlab('Origin - Chemo')+ ylab('Terpene concentration (sqrt(mg/g))')+
  theme_bw(base_size = 24)
fig1


chem2plotsupp <- chem2 %>%  
  mutate(compound=if_else(conc>.25, chemical, "other"), compound=as.factor(compound)) %>% 
  ungroup() %>% 
  mutate(compound=fct_relevel(compound, c("other","a_thujene","myrcene","octen_3_ol","a_terpinene","gama_terpinene",
                                          "p_cymene","thymoquinone","carvacrol","thymol")))%>% 
  group_by(Region,Chemo,Site,compound) %>% summarize(conc=mean(conc, na.rm=T))

fig1supp <- ggplot(chem2plotsupp , aes(fill=compound, x=as.factor(Site), y=sqrt(conc))) + 
  geom_bar(position="stack",stat="identity")+
  #facet_wrap(~Season) + 
  #scale_fill_viridis(discrete=T,labels=c('Other','α-Thujene','Myrcene','Octen-3-ol','α-Terpinene', 'γ-Terpinene',
  #                                      'p-Cymene','Thymoquinone','Carvacrol','Thymol'), name="Compound")+
  scale_fill_jco(labels=c('Other','α-Thujene','Myrcene','Octen-3-ol','α-Terpinene', 'γ-Terpinene',
                          'p-Cymene','Thymoquinone','Carvacrol','Thymol'), name="Compound")+
  #ggtitle("B")+
  xlab('Origin - Chemo')+ ylab('Terpene concentration (sqrt(mg/g))')+
  theme_bw(base_size = 24)+
  facet_wrap(~paste(Region,Chemo,sep="-"), scales="free_x")
fig1supp

### distance by space and chemistry ####
chem2dist <- chem2 %>% 
  full_join(.,sites) %>% 
  ungroup() %>% 
  select(Region,Site,Lat,Long,chemical,conc) %>%
  group_by(Region,Site) %>% 
  pivot_wider(names_from = "chemical", values_from = "conc", values_fn = mean) %>% 
  select(-`NA`) %>% 
  drop_na()

distgeo <- geodist(chem2dist[3:4], measure='geodesic', upper=F)
distgeo <- dist(chem2dist[3:4])
distchm <- dist(chem2dist[5:29])
mantel(distgeo, distchm)
plot(distgeo, distchm)

chem2dist1 <- chem2dist %>% filter(Region=='MT')
distgeo1 <- geodist(chem2dist1[3:4], measure='geodesic', upper=F)
distgeo1 <- dist(chem2dist1[3:4])
distchm1 <- dist(chem2dist1[5:29])

mantel(distgeo1, distchm1)
plot(distgeo1, distchm1)

## extract for chemical diversity indices ####
chem1dat <- chem1[7:31]*10000 %>% round(., digits=0)
#chem1dat <- ifelse(chem1[7:31]>0,1,0)
est_rich <- estimateD(t(chem1dat), datatype="abundance", q=0)
# est_rich01 <- estimateD(t(chem1dat), datatype="abundance", q=0, base="coverage", level=.95, nboot=0)
est_rich1 <- estimateD(t(chem1dat), datatype="abundance", q=1)
#est_rich <- estimateR(chem1dat) %>% as.data.frame()
chem1$qD <- 0
chem1$q1D <- 0

chem1[33] <- est_rich[6]
chem1[34] <- est_rich1[6]
chem1$rich <- hill_taxa(chem1[7:31], q=0)
chem1$shan <- hill_taxa(chem1[7:31], q=1)
chem1$simp <- hill_taxa(chem1[7:31], q=2)
chem1$NMDS1 <- scores(chem.mds)$sites[,1]
chem1$NMDS2 <- scores(chem.mds)$sites[,2]


chem1$rich[chem1$rich==0] <- NA
chem1$shan[chem1$shan==1] <- NA
chem1$simp[chem1$simp=='Inf'] <- NA

plot(chem1$rich,chem1$qD)
plot(chem1$shan,chem1$q1D)
plot(log(chem1$total),chem1$NMDS1)
cor.test(log(chem1$total),chem1$NMDS1)
cor.test(chem1$Chemo,chem1$NMDS1)
summary(lm(NMDS2~Chemo,data=chem1))

## chemistry by region - Plots and Analysis ####
ggplot(chem1 , aes(x=Region, y=log(total), fill=Chemo)) + geom_boxplot()+
  scale_fill_viridis(discrete=T, direction=-1, option='C', end=.8)+theme_bw()
ggplot(chem1 , aes(x=Region, y=shan, fill=Chemo)) + geom_boxplot()+
  scale_fill_viridis(discrete=T, direction=-1, option='C', end=.8)+theme_bw()



terptot <- ggplot(chem1 , aes(x=Region, y=log(total), color=Chemo, fill=Chemo)) + 
  geom_boxplot(lwd=1, outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=.2),  size=2 ) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1)+
  scale_fill_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.25)+
  scale_x_discrete("Region")+
  scale_y_continuous("Total terpene concentration (log(mg/g))")+
  theme_bw(base_size=20)
terptot

terprich <- ggplot(chem1 , aes(x=Region, y=qD, color=Chemo, fill=Chemo)) + 
  geom_boxplot(lwd=1, outlier.shape=NA, alpha=.25) +
  geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=.2),  size=2 ) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  #scale_color_viridis(discrete=T, option="D", begin=.8, end=1, direction=1)+
  #scale_fill_viridis(discrete=T, option="D", begin=.8, end=1, direction=1, alpha=.25)+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_fill_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_x_discrete("Origin")+
  scale_y_continuous("Terpene richness")+
  theme_bw(base_size=24)
terprich


terpshan <- ggplot(chem1 , aes(x=Region, y=shan, color=Chemo, fill=Chemo)) + 
  geom_boxplot(lwd=1, outlier.shape=NA, alpha=.25) +
  geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=.2),  size=2 ) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  #scale_color_viridis(discrete=T, option="D", begin=.8, end=1, direction=1)+
  #scale_fill_viridis(discrete=T, option="D", begin=.8, end=1, direction=1, alpha=.25)+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_fill_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_x_discrete("Origin")+
  scale_y_continuous("Terpene shannon diversity")+
  theme_bw(base_size=24)
terpshan




terppc1 <- ggplot(chem1 , aes(x=Region, y=pc1, color=Chemo, fill=Chemo)) + 
  geom_boxplot(lwd=1, outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=.2),  size=2 ) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1)+
  scale_fill_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.25)+
  scale_x_discrete("Region")+
  scale_y_continuous("PC axis 1")+
  theme_bw(base_size=20)
terppc1

terppc2 <- ggplot(chem1 , aes(x=log(total), y=shan, color=Region, fill=Region)) + 
  geom_point()+
  #geom_boxplot(lwd=1, outlier.shape=NA) +
  #geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=.2),  size=2 ) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1)+
  scale_fill_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.25)+
  scale_x_discrete("total terps")+
  scale_y_continuous("PC axis 1")+
  theme_bw(base_size=20)
terppc2

### nmds plot #####
pcrot <- scores(chem.mds)$species %>% as.data.frame() %>% rownames_to_column(var="compound") %>% 
  filter(compound=='beta_ocimene'| compound=='borneol'| compound=='eugenol'| compound=='a_terpineol'| 
         compound=='cym_8_ol'|compound=='terp_4_ol'|compound=='thymol'|compound=='carvacrol')

pcrot$compound1 <- c("β-Ocimene","Borneol","Terp-4-ol","Cym-8-ol","α-Terpineol","Thymol","Carvacrol","Eugenol")

pcrot1 <- scores(chem.mds)$species %>% as.data.frame() %>% rownames_to_column(var="compound") %>% 
  filter(compound=='a_terpinene'| compound=='a_terpinene'| compound=='gama_terpinene'| compound=='myrcene'| 
           compound=='octen_3_ol'|compound=='p_cymene'|compound=='thymoquinone'|compound=='thymol'|compound=='carvacrol')

nmdsfig <- ggplot() + xlab("NMDS1")+ylab("NMDS2")+ 
  #stat_ellipse(data=chem1, geom="polygon" ,aes(x=NMDS1,y=NMDS2, color=Region), lwd=1, alpha=.25)+
  geom_point(data=chem1, aes(x=NMDS1, y=NMDS2, color=Region, shape=Chemo), size=3,alpha=.5) +
  #scale_fill_manual(values=c("darkgoldenrod1","dodgerblue"), labels=c("Low productivity","High productivity")) +
  #scale_color_manual(values=c("darkgoldenrod1","dodgerblue"), labels=c("Low productivity","High productivity")) +
  #geom_segment(data=pcrot, aes(x=0, y=0,xend=NMDS1*.5,yend=NMDS2*.5), arrow=arrow(length=unit(0.03, "npc")), color='black', lwd=1)+
  geom_segment(data=pcrot, aes(x=0, y=0,xend=NMDS1*.5,yend=NMDS2*.5), arrow=arrow(length=unit(0.03, "npc")), color='black', lwd=1)+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, name="Origin")+
  geom_text(data=pcrot, aes(x=NMDS1*.5, y=NMDS2*.6+c(.03,-.05,-.09,.03,.05,0,0,.06), label=compound1), fontface="bold", color="BLACK", size=5)+
  #geom_text(data=pcrot, aes(x=NMDS1*.5, y=NMDS2*.6, label=compound), fontface="bold", color="BLACK", size=4)+
  xlim(-1.3,1.2)+
  theme_bw(base_size = 24)
nmdsfig

## PHYTOCHEM DIVERSITY PLOT ####
terpplot <- (fig1+nmdsfig)/ ((terprich+terpshan) + plot_layout(guides='collect')) +
  plot_annotation(tag_levels = 'a') #+ plot_layout(guides='collect') & theme(legend.position='top')
ggsave('terpplot.tiff',terpplot, width=15, height=12, units="in", dpi=600, compression = "lzw")


### analysis for chem diversity ####
ter1 <- glmmTMB(total ~ Region * Chemo + (1|Site), data=chem1, family=lognormal)
Anova(ter1)
r.squaredGLMM(ter1)
simulateResiduals(ter1, plot=T)
emmeans(ter1, pairwise~Region, type="response")

chem1$Chem0 <- ifelse(chem1$Chemo=="T",1,0)

ter1a <- glmmTMB(Chem0 ~ Region + (1|Site), data=chem1, family = binomial)
Anova(ter1a)
r.squaredGLMM(ter1a)
simulateResiduals(ter1a, plot=T)
emmeans(ter1a, pairwise~Region, type="response")
ter1ablup <- tidy(ter1a, effect="ran_vals") %>% as.data.frame() %>% 
  mutate(Tprob=boot::inv.logit(estimate),)


ter2 <- glmmTMB(qD ~ Region * Chemo + (1|Site), data=chem1, family=lognormal)
Anova(ter2)
r.squaredGLMM(ter2)
simulateResiduals(ter2, plot=T)
check_overdispersion(ter2)
emmeans(ter2, pairwise~Region:Chemo, type="response")
emmeans(ter2, pairwise~Region, type="response")

ter3 <- glmmTMB(shan ~ Region * Chemo + (1|Site), data=chem1)
Anova(ter3)
r.squaredGLMM(ter3)
simulateResiduals(ter3, plot=T)
emmeans(ter3, pairwise~Chemo, type="response")


################################
# Q2: new chem and field data ####
f20 <- read_csv("Data/MonardaGarden2020.csv") %>% filter(survey=='Aug') %>% 
  rowwise() %>% 
  mutate(leaf_dmg=mean(c(leaf_1,leaf_2,leaf_3,leaf_4,leaf_5,leaf_6,leaf_7,leaf_8,leaf_9,leaf_10), na.rm=T),
         leaf_dmgp=leaf_dmg/100,
         leaf_dmgp1=(leaf_dmgp*(358-1)+0.5)/358,
         mines=mines_blotch+mines_linear,
         mined=if_else(mines>.1,1,0),
         ants=gsub('2;6',8,ants) %>% as.numeric()) %>% 
  select(Region,Site,Plant,Rep,height,leaves_total,leaves_damaged,herbpercentage,leaf_dmg,leaf_dmgp, leaf_dmgp1,mildew,galls_leaf,galls_stem,aphids_unparasitized,mines,mined,leps,ants)


chemnew <- chem1 %>% 
  mutate(Chemo=if_else(thymol>carvacrol, "T","C")) %>% 
  select(Region,Site,Plant,Biomass_g,total,qD,rich,shan,simp,Chemo)

chemnew1 <- full_join(chemnew,f20) %>% drop_na(total) #%>% 
  #mutate(total=scale(log(total))[,1],rich=scale(log(rich))[,1],shan=scale(shan)[,1])

## chewing damage models ####
clm1a <- glmmTMB(leaf_dmgp~Region+(1|Site), data=chemnew1, family=ordbeta)
clm1b <- glmmTMB(leaf_dmgp~total+(1|Site), data=chemnew1, family=ordbeta)
clm1c <- glmmTMB(leaf_dmgp~Chemo+(1|Site), data=chemnew1, family=ordbeta)
clm1d <- glmmTMB(leaf_dmgp~qD+(1|Site), data=chemnew1, family=ordbeta)
clm1e <- glmmTMB(leaf_dmgp~shan+(1|Site), data=chemnew1, family=ordbeta)
clm1f <- glmmTMB(leaf_dmgp~Region+total+(1|Site), data=chemnew1, family=ordbeta)
clm1g <- glmmTMB(leaf_dmgp~Region+qD+(1|Site), data=chemnew1, family=ordbeta)
clm1h <- glmmTMB(leaf_dmgp~Region+shan+(1|Site), data=chemnew1, family=ordbeta)
clm1i <- glmmTMB(leaf_dmgp~Region+Chemo+(1|Site), data=chemnew1, family=ordbeta)
clm1j <- glmmTMB(leaf_dmgp~Chemo+total+(1|Site), data=chemnew1, family=ordbeta)
clm1k <- glmmTMB(leaf_dmgp~Chemo+qD+(1|Site), data=chemnew1, family=ordbeta)
clm1l <- glmmTMB(leaf_dmgp~Chemo+shan+(1|Site), data=chemnew1, family=ordbeta)
clm1m <- glmmTMB(leaf_dmgp~Region*total+(1|Site), data=chemnew1, family=ordbeta)
clm1n <- glmmTMB(leaf_dmgp~Region*qD+(1|Site), data=chemnew1, family=ordbeta)
clm1o <- glmmTMB(leaf_dmgp~Region*shan+(1|Site), data=chemnew1, family=ordbeta)
clm1p <- glmmTMB(leaf_dmgp~Region*Chemo+(1|Site), data=chemnew1, family=ordbeta)
clm1q <- glmmTMB(leaf_dmgp~Chemo*total+(1|Site), data=chemnew1, family=ordbeta)
clm1r <- glmmTMB(leaf_dmgp~Chemo*qD+(1|Site), data=chemnew1, family=ordbeta)
clm1s <- glmmTMB(leaf_dmgp~Chemo*shan+(1|Site), data=chemnew1, family=ordbeta)
clm1 <- glmmTMB(leaf_dmgp~1+(1|Site), data=chemnew1, family=ordbeta)

clm1a <- glmmTMB(logit(leaf_dmgp1)~Region+(1|Site), data=chemnew1, family=gaussian)
clm1b <- glmmTMB(logit(leaf_dmgp1)~total+(1|Site), data=chemnew1, family=gaussian)
clm1c <- glmmTMB(logit(leaf_dmgp1)~Chemo+(1|Site), data=chemnew1, family=gaussian)
clm1d <- glmmTMB(logit(leaf_dmgp1)~qD+(1|Site), data=chemnew1, family=gaussian)
clm1e <- glmmTMB(logit(leaf_dmgp1)~shan+(1|Site), data=chemnew1, family=gaussian)
clm1f <- glmmTMB(logit(leaf_dmgp1)~Region+total+(1|Site), data=chemnew1, family=gaussian)
clm1g <- glmmTMB(logit(leaf_dmgp1)~Region+qD+(1|Site), data=chemnew1, family=gaussian)
clm1h <- glmmTMB(logit(leaf_dmgp1)~Region+shan+(1|Site), data=chemnew1, family=gaussian)
clm1i <- glmmTMB(logit(leaf_dmgp1)~Region+Chemo+(1|Site), data=chemnew1, family=gaussian)
clm1j <- glmmTMB(logit(leaf_dmgp1)~Chemo+total+(1|Site), data=chemnew1, family=gaussian)
clm1k <- glmmTMB(logit(leaf_dmgp1)~Chemo+qD+(1|Site), data=chemnew1, family=gaussian)
clm1l <- glmmTMB(logit(leaf_dmgp1)~Chemo+shan+(1|Site), data=chemnew1, family=gaussian)
clm1m <- glmmTMB(logit(leaf_dmgp1)~Region*total+(1|Site), data=chemnew1, family=gaussian)
clm1n <- glmmTMB(logit(leaf_dmgp1)~Region*qD+(1|Site), data=chemnew1, family=gaussian)
clm1o <- glmmTMB(logit(leaf_dmgp1)~Region*shan+(1|Site), data=chemnew1, family=gaussian)
clm1p <- glmmTMB(logit(leaf_dmgp1)~Region*Chemo+(1|Site), data=chemnew1, family=gaussian)
clm1q <- glmmTMB(logit(leaf_dmgp1)~Chemo*total+(1|Site), data=chemnew1, family=gaussian)
clm1r <- glmmTMB(logit(leaf_dmgp1)~Chemo*qD+(1|Site), data=chemnew1, family=gaussian)
clm1s <- glmmTMB(logit(leaf_dmgp1)~Chemo*shan+(1|Site), data=chemnew1, family=gaussian)
clm1 <- glmmTMB(logit(leaf_dmgp1)~1+(1|Site), data=chemnew1, family=gaussian)


bbmle::ICtab(clm1,clm1a,clm1b,clm1c,clm1d,clm1e,clm1f,clm1g,clm1h,clm1i,clm1j,clm1k,clm1l,
             clm1m,clm1n,clm1o,clm1p,clm1q,clm1r,clm1s, weights=T, delta=T, base=T, type="BIC")

clm_mods <- subset(MuMIn::model.sel(clm1,clm1a,clm1b,clm1c,clm1d,clm1e,clm1f,clm1g,clm1h,clm1i,clm1j,clm1k,clm1l,
                 clm1m,clm1n,clm1o,clm1p,clm1q,clm1r,clm1s,rank="BIC"), delta <= 4)
sw(clm_mods)
summary(clm_mods)

clm_avg <- model.avg(clm_mods)
summary(clm_avg)

clm1b <- glmmTMB(leaf_dmgp1~total+(1|Site), data=chemnew1, family=beta_family)
clm1f <- glmmTMB(leaf_dmgp1~Region+total+(1|Site), data=chemnew1, family=beta_family)

Anova(clm1a)
Anova(clm1b)
Anova(clm1f)


check_collinearity(clm1f)

clm1a <- glmmTMB(leaf_dmgp~Region+(1|Site), data=chemnew1, family=ordbeta)
clm1b <- glmmTMB(logit(leaf_dmgp1)~log(total)+(1|Site), data=chemnew1)
clm1c <- glmmTMB(logit(leaf_dmgp1)~Region*Chemo+(1|Site), data=chemnew1)

Anova(clm1a)
summary(clm1a)
r.squaredGLMM(clm1a)
r2(clm1a)
clm1asim <- simulateResiduals(clm1a)
testZeroInflation(clm1asim)
plot(clm1asim)


chewem1 <- emmeans(clm1a, ~Region, type="response") %>% as.data.frame()

pd1 <- ggplot(chewem1 , aes(x=Region, y=response, ymin=asymp.LCL, ymax=asymp.UCL, color=Region, fill=Region)) + 
  #geom_point(position=position_jitter(height=0, width=.2),  size=2, color="black", alpha=.25, pch=21 ) + #color="#287D8EFF",
  #geom_boxplot(lwd=1, outlier.shape=NA) +
  geom_pointrange(size=1, lwd=1) +
  scale_x_discrete("")+
  scale_y_continuous(labels = scales::percent, name='')+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1)+
  theme_bw(base_size=14)+
  theme(legend.position = "none")+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent')) #transparent legend panel)
pd1

pd1a <- ggplot(chemnew1 )+
  geom_point(aes(x=total, y=leaf_dmgp, color=Region),size=2)+
  geom_smooth(aes(x=total, y=leaf_dmgp1), color="black", method="glm", method.args=list(family=beta_family))+
  scale_x_continuous("Terpene conc. (mg/g)")+
  scale_y_continuous(labels = scales::percent, name='Leaf damage')+
  #scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.75)+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.75, name="Origin")+
  theme_bw(base_size = 20)

pd1b <- ggplot(chemnew1 , aes(x=Region, y=leaf_dmgp1, color=Region, fill=Region)) + 
  geom_boxplot(lwd=1, outlier.shape=NA) +
  geom_point(position=position_jitter(height=0, width=.2),  size=2, color="black", alpha=.25, pch=21 ) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  scale_x_discrete("Region")+
  scale_y_continuous(labels = scales::percent, name='Leaf damage')+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1)+
  scale_fill_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.25)+
  theme_bw(base_size=20)
pd1b

layoutpd1 <- c(patchwork::area(t=1, l=1, b=5, r=4),
               patchwork::area(t=1.1, l=4, b=2, r=4))
pd1plot <- pd1a+pd1+plot_layout(design=layoutpd1)


dmgplot <- (pd1b+pd1a)+
  plot_annotation(tag_levels = 'a') + plot_layout(guides='collect') & theme(legend.position='right')
#ggsave('dmgplot.tiff',dmgplot, width=10, height=4, units="in", dpi=600, compression = "lzw")

## Gall model and plots ####

clm3a <- glmmTMB(galls_leaf~Region+(1|Site), data=chemnew1, family=nbinom2)
clm3b <- glmmTMB(galls_leaf~total+(1|Site), data=chemnew1, family=nbinom2)
clm3c <- glmmTMB(galls_leaf~Chemo+(1|Site), data=chemnew1, family=nbinom2)
clm3d <- glmmTMB(galls_leaf~qD+(1|Site), data=chemnew1, family=nbinom2)
clm3e <- glmmTMB(galls_leaf~shan+(1|Site), data=chemnew1, family=nbinom2)
clm3f <- glmmTMB(galls_leaf~Region+total+(1|Site), data=chemnew1, family=nbinom2)
clm3g <- glmmTMB(galls_leaf~Region+qD+(1|Site), data=chemnew1, family=nbinom2)
clm3h <- glmmTMB(galls_leaf~Region+shan+(1|Site), data=chemnew1, family=nbinom2)
clm3i <- glmmTMB(galls_leaf~Region+Chemo+(1|Site), data=chemnew1, family=nbinom2)
clm3j <- glmmTMB(galls_leaf~Chemo+total+(1|Site), data=chemnew1, family=nbinom2)
clm3k <- glmmTMB(galls_leaf~Chemo+qD+(1|Site), data=chemnew1, family=nbinom2)
clm3l <- glmmTMB(galls_leaf~Chemo+shan+(1|Site), data=chemnew1, family=nbinom2)
clm3m <- glmmTMB(galls_leaf~Region*total+(1|Site), data=chemnew1, family=nbinom2)
clm3n <- glmmTMB(galls_leaf~Region*qD+(1|Site), data=chemnew1, family=nbinom2)
clm3o <- glmmTMB(galls_leaf~Region*shan+(1|Site), data=chemnew1, family=nbinom2)
clm3p <- glmmTMB(galls_leaf~Region*Chemo+(1|Site), data=chemnew1, family=nbinom2)
clm3q <- glmmTMB(galls_leaf~Chemo*total+(1|Site), data=chemnew1, family=nbinom2)
clm3r <- glmmTMB(galls_leaf~Chemo*qD+(1|Site), data=chemnew1, family=nbinom2)
clm3s <- glmmTMB(galls_leaf~Chemo*shan+(1|Site), data=chemnew1, family=nbinom2)
clm3 <- glmmTMB(galls_leaf~1+(1|Site), data=chemnew1, family=nbinom2)




bbmle::ICtab(clm3,clm3a,clm3b,clm3c,clm3d,clm3e,clm3f,clm3g,clm3h,clm3i,clm3j,clm3k,clm3l,
             clm3m,clm3n,clm3o,clm3p,clm3q,clm3r,clm3s, weights=T, delta=T, base=T, type="BIC")

clm3_mods <- subset(MuMIn::model.sel(clm3,clm3a,clm3b,clm3c,clm3d,clm3e,clm3f,clm3g,clm3h,clm3i,clm3j,clm3k,clm3l,
                                    clm3m,clm3n,clm3o,clm3p,clm3q,clm3r,clm3s,rank="BIC"), delta <= 4)
sw(clm3_mods)
summary(clm3_mods)

Anova(clm3b)
Anova(clm3a)
Anova(clm3d)
Anova(clm3f)

clm_avg <- model.avg(clm3_mods)
summary(clm_avg)
sw(clm3_mods)
summary(clm3_mods)

clm_avg <- model.avg(clm3_mods)
summary(clm3_mods)

gallsem <- emmeans(clm3a, ~Region, type="response") %>% as.data.frame()

pg3x <- ggplot(gallsem , aes(x=Region, y=response, ymin=asymp.LCL, ymax=asymp.UCL, color=Region, fill=Region)) + 
  #geom_point(position=position_jitter(height=0, width=.2),  size=2, color="black", alpha=.25, pch=21 ) + #color="#287D8EFF",
  #geom_boxplot(lwd=1, outlier.shape=NA) +
  geom_pointrange(size=1, lwd=1) +
  scale_x_discrete("")+
  scale_y_continuous(name='')+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1)+
  scale_fill_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.25)+
  theme_bw(base_size=14)+
  theme(legend.position = "none")+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent')) #transparent legend panel)
pg3x

pg3a <- ggplot(chemnew1, )+
  geom_point(aes(x=total, y=galls_leaf, color=Region),  size=2)+
  geom_smooth(aes(x=total, y=galls_leaf), color="black", method="glm.nb")+
  scale_x_continuous("Terpene conc. (mg/g)")+
  scale_y_continuous(name='# leaf galls')+
  #scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.75)+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.75, name="Origin")+
  theme_bw(base_size = 20)

pg3b <- ggplot(chemnew1, )+
  geom_point(aes(x=qD, y=galls_leaf, fill=Region),color="black",  size=2, pch=21)+
  geom_smooth(aes(x=qD, y=galls_leaf), color="black", method="glm.nb")+
  scale_x_continuous("Terpene richness")+
  scale_y_continuous(name='# leaf galls')+
  #scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.75)+
  scale_fill_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.5)+
  theme_bw(base_size = 20)

pg3c <- ggplot(chemnew1 , aes(x=Region, y=galls_leaf, color=Region, fill=Region)) + 
  geom_boxplot(lwd=1, outlier.shape=NA) +
  geom_point(position=position_jitter(height=0, width=.2),  size=2) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  scale_x_discrete("Region")+
  scale_y_continuous(name='Leaf damage')+
  scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.25)+
  scale_fill_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.25)+
  theme_bw(base_size=20)
pg3c

gallplot <- (pg3a+pg3b+pg3c)+
  plot_annotation(tag_levels = 'a') + plot_layout(guides='collect') & theme(legend.position='right')
ggsave('gallplot.tiff',gallplot, width=12, height=4, units="in", dpi=600, compression = "lzw")

layoutpg3 <- c(patchwork::area(t=1, l=1, b=5, r=4),
               patchwork::area(t=1.1, l=4, b=2, r=4))
pg3plot <- pg3a+pg3x+plot_layout(design=layoutpg3)

Anova(clm3d)

## aphid abundance ####

clm2a <- glmmTMB(aphids_unparasitized~Region+(1|Site), data=chemnew1, family=nbinom2)
clm2b <- glmmTMB(aphids_unparasitized~total+(1|Site), data=chemnew1, family=nbinom2)
clm2c <- glmmTMB(aphids_unparasitized~Chemo+(1|Site), data=chemnew1, family=nbinom2)
clm2d <- glmmTMB(aphids_unparasitized~qD+(1|Site), data=chemnew1, family=nbinom2)
clm2e <- glmmTMB(aphids_unparasitized~shan+(1|Site), data=chemnew1, family=nbinom2)
clm2f <- glmmTMB(aphids_unparasitized~Region+total+(1|Site), data=chemnew1, family=nbinom2)
clm2g <- glmmTMB(aphids_unparasitized~Region+qD+(1|Site), data=chemnew1, family=nbinom2)
clm2h <- glmmTMB(aphids_unparasitized~Region+shan+(1|Site), data=chemnew1, family=nbinom2)
clm2i <- glmmTMB(aphids_unparasitized~Region+Chemo+(1|Site), data=chemnew1, family=nbinom2)
clm2j <- glmmTMB(aphids_unparasitized~Chemo+total+(1|Site), data=chemnew1, family=nbinom2)
clm2k <- glmmTMB(aphids_unparasitized~Chemo+qD+(1|Site), data=chemnew1, family=nbinom2)
clm2l <- glmmTMB(aphids_unparasitized~Chemo+shan+(1|Site), data=chemnew1, family=nbinom2)
clm2m <- glmmTMB(aphids_unparasitized~Region*total+(1|Site), data=chemnew1, family=nbinom2)
clm2n <- glmmTMB(aphids_unparasitized~Region*qD+(1|Site), data=chemnew1, family=nbinom2)
clm2o <- glmmTMB(aphids_unparasitized~Region*shan+(1|Site), data=chemnew1, family=nbinom2)
clm2p <- glmmTMB(aphids_unparasitized~Region*Chemo+(1|Site), data=chemnew1, family=nbinom2)
clm2q <- glmmTMB(aphids_unparasitized~Chemo*total+(1|Site), data=chemnew1, family=nbinom2)
clm2r <- glmmTMB(aphids_unparasitized~Chemo*rich+(1|Site), data=chemnew1, family=nbinom2)
clm2s <- glmmTMB(aphids_unparasitized~Chemo*shan+(1|Site), data=chemnew1, family=nbinom2)
clm2 <- glmmTMB(aphids_unparasitized~1+(1|Site), data=chemnew1, family=nbinom2)

bbmle::ICtab(clm2,clm2a,clm2b,clm2c,clm2d,clm2e,clm2f,clm2g,clm2h,clm2i,clm2j,clm2k,clm2l,
             clm2m,clm2n,clm2o,clm2p,clm2q,clm2r,clm2s, weights=T, delta=T, base=T, type="BIC")

clm2_mods <- subset(MuMIn::model.sel(clm2,clm2a,clm2b,clm2c,clm2d,clm2e,clm2f,clm2g,clm2h,clm2i,clm2j,clm2k,clm2l,
                                     clm2m,clm2n,clm2o,clm2p,clm2q,clm2r,clm2s,rank="BIC"), delta <= 4)
sw(clm2_mods)
summary(clm2_mods)

clm_avg <- model.avg(clm2_mods)
summary(clm_avg)
sw(clm2_mods)
summary(clm2_mods)

clm_avg <- model.avg(clm2_mods)
summary(clm2_mods)

Anova(clm2a)
r2(clm2h)
summary(clm2b)
Anova(clm2g)
emmeans(clm2g, pairwise~Region, type="response")

pa1 <- ggplot(chemnew1)+
  geom_point(aes(x=total, y=aphids_unparasitized, color=Region), size=2)+
  geom_smooth(aes(x=total, y=aphids_unparasitized), method="glm.nb")+
  scale_x_continuous("Terpene conc. (log mg/g)")+
  scale_y_continuous(name='# aphids')+
  theme_bw(base_size = 20)

pa1a <- ggplot(chemnew1 %>% drop_na(Chemo), aes(x=Region, y=aphids_unparasitized, color=Region, fill=Region)) + 
  geom_boxplot(lwd=1, outlier.shape=NA) +
  geom_jitter(height=0,width=.2)+
  #geom_point(position=position_jitterdodge(jitter.height=0, jitter.width=.2),  size=2 ) + #color="#287D8EFF",
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.75, name="Origin")+
  scale_fill_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.75, name="Origin")+
  scale_x_discrete("Origin")+
  scale_y_continuous("Aphids")+
  theme_bw(base_size=20)+
  guides(fill="none", color="none")
pa1a

ggsave('aphidplot.tiff',pa1a, width=6, height=4, units="in", dpi=600, compression = "lzw")


## mildew abundance ####
chemnew1$mildew <- ifelse(chemnew1$mildew>2, 1, chemnew1$mildew)
clm5a <- glmmTMB(mildew~shan+(1|Site), data=chemnew1, family=binomial)
clm5b <- glmmTMB(mildew~total+(1|Site), data=chemnew1, family=binomial)
clm5c <- glmmTMB(mildew~Region+(1|Site), data=chemnew1, family=binomial)
clm5d <- glmmTMB(mildew~Chemo+(1|Site), data=chemnew1, family=binomial)
clm5e <- glmmTMB(mildew~Chemo*shan+(1|Site), data=chemnew1, family=binomial)
clm5f <- glmmTMB(mildew~Chemo*log(total)+(1|Site), data=chemnew1, family=binomial)
clm5g <- glmmTMB(mildew~Chemo*Region+(1|Site), data=chemnew1, family=binomial)
clm5h <- glmmTMB(mildew~Region*log(total)+(1|Site), data=chemnew1, family=binomial)
clm5i <- glmmTMB(mildew~Region*shan+(1|Site), data=chemnew1, family=binomial)
clm5 <- glmmTMB(mildew~1+(1|Site), data=chemnew1, family=binomial)

bbmle::AICtab(clm5,clm5a,clm5b,clm5c,clm5d,clm5e,clm5f,clm5g,clm5h,clm5i, weights=T, delta=T, base=T)

MuMIn::model.sel(clm5,clm5a,clm5b,clm5c,clm5d,clm5e,clm5f,clm5g,clm5h,clm5i)

check_collinearity(clm5h)

clm5avg <- model.avg(clm5h,clm5g,clm5f,clm5b,clm5e,clm5a,clm5d)
sw(clm5avg)
summary(clm5avg)

Anova(clm5c)
summary(clm5i)
Anova(clm5g)
emmeans(clm5g, pairwise~Region, type="response")

pa1 <- ggplot(chemnew1, aes(x=pc1, y=mildew, color=Region))+
  geom_point(size=2)+
  geom_smooth(method="glm.nb")+
  scale_x_continuous("Terpene conc. (log mg/g)")+
  scale_y_continuous(name='# aphids')+
  theme_bw(base_size = 20)


## SEED HEAD MODELS
### SEED LOSS -- updated data in TEAMS monarda garden folder ###################################
seedhead21 <- read_csv("Data/MonardaGarden_SeedDamage2021 (5).csv") %>% drop_na(Seeds_No) %>% drop_na(Region) %>% mutate(HeadSize=as.numeric(HeadSize))
gard21 <- read_csv("Data/MonardaGarden_Data2021.csv")

ggplot(seedhead21 , aes(x=HeadSize, y=Seeds_No, color=as.factor(HerbDmg))) +  
  geom_point() + geom_smooth(method="glm.nb") + facet_wrap(~Region)

seed1 <- glm(Seeds_No ~ 0+HeadSize , data=seedhead21 %>% filter(HerbDmg==0), family=poisson)
Anova(seed1)
summary(seed1)
r.squaredGLMM(seed1)

seedhead21$Seeds_NoDmg <- (exp(seedhead21$HeadSize*.266869)) 
plot(Seed_Loss ~ Seeds_No, data=seedhead21)
abline(a=0,b=1)

seedhead21 <- seedhead21 %>% mutate(Seed_Loss = case_when(HerbDmg==1&Region=='MT' ~ Seeds_NoDmg*.6,
                                                          HerbDmg==1&Region=='WI' ~ Seeds_NoDmg*.42,
                                                          HerbDmg==0&Region=='MT' ~ 0,
                                                          HerbDmg==0&Region=='WI' ~ 0)) 
seedhead21$PercSeedLoss <- 1-ifelse(seedhead21$HerbDmg==1, seedhead21$Seeds_No/seedhead21$Seeds_NoDmg, 1)
seedhead21$PercSeedLoss1 <- (seedhead21$PercSeedLoss*(length(seedhead21$PercSeedLoss)-1)+.5)/length(seedhead21$PercSeedLoss)

ggplot(seedhead21 , aes(x=as.factor(HerbDmg), y=Seed_Loss)) +  geom_boxplot() 

### JOIN SEED AND GARDEN DATA ####
garden1 <- full_join(chem1,seedhead21)
garden2 <- full_join(gard21,garden1, by=c("Region", "Site", "Plant", "Rep")) %>% drop_na(Region)
garden2$piercing <- (garden2$prc*(length(garden2$prc)-1)+.5)/length(garden2$prc)/100
garden2$chewing <- (garden2$ch*(length(garden2$ch)-1)+.5)/length(garden2$ch)/100
garden2$Chemo <- ifelse(garden2$carvacrol>garden2$thymol, "C","T")

head(garden2)


garden2 <- garden2 %>% mutate(Pop = case_when(Site == 'BEN' ~ 'WI1',Site=='BOL'~'WI2',Site=='SDA'~'WI3',Site=='SRD'~'WI4',Site=='STB'~'WI5',Site=='UWM'~'WI6',
                                              Site=='BF'~'MT1',Site=='CWN'~'MT2',Site=='MJS'~'MT3',Site=='MS'~'MT4',Site=='MPG'~'MT5',Site=='NH'~'MT6',)) 

garden2 <- garden2 %>% select(Region, Site, Plant, Rep, hg_Tallest, Biomass_g, total, fl, fr, Chemo, qD, shan, HerbDmg) %>% 
  drop_na(HerbDmg,total)

ggplot(chem1  %>% 
         filter(Region=='WI') , aes(x=total, y=rich, color=Site)) + geom_point() + geom_smooth(method="lm", formula=y~log(x))+
  scale_color_viridis(discrete = T, option='C') + labs(x="total concentration")


## SEED HEAD DAMAGE MODELS #### 

clm4a <- glmmTMB(HerbDmg~Region+(1|Site/Plant), data=garden2, family=binomial)
clm4b <- glmmTMB(HerbDmg~total+(1|Site/Plant), data=garden2, family=binomial)
clm4c <- glmmTMB(HerbDmg~Chemo+(1|Site/Plant), data=garden2, family=binomial)
clm4d <- glmmTMB(HerbDmg~qD+(1|Site/Plant), data=garden2, family=binomial)
clm4e <- glmmTMB(HerbDmg~shan+(1|Site/Plant), data=garden2, family=binomial)
clm4f <- glmmTMB(HerbDmg~Region+total+(1|Site/Plant), data=garden2, family=binomial)
clm4g <- glmmTMB(HerbDmg~Region+qD+(1|Site/Plant), data=garden2, family=binomial)
clm4h <- glmmTMB(HerbDmg~Region+shan+(1|Site/Plant), data=garden2, family=binomial)
clm4i <- glmmTMB(HerbDmg~Region+Chemo+(1|Site/Plant), data=garden2, family=binomial)
clm4j <- glmmTMB(HerbDmg~Chemo+total+(1|Site/Plant), data=garden2, family=binomial)
clm4k <- glmmTMB(HerbDmg~Chemo+qD+(1|Site/Plant), data=garden2, family=binomial)
clm4l <- glmmTMB(HerbDmg~Chemo+shan+(1|Site/Plant), data=garden2, family=binomial)
clm4m <- glmmTMB(HerbDmg~Region*total+(1|Site/Plant), data=garden2, family=binomial)
clm4n <- glmmTMB(HerbDmg~Region*qD+(1|Site/Plant), data=garden2, family=binomial)
clm4o <- glmmTMB(HerbDmg~Region*shan+(1|Site/Plant), data=garden2, family=binomial)
clm4p <- glmmTMB(HerbDmg~Region*Chemo+(1|Site/Plant), data=garden2, family=binomial)
clm4q <- glmmTMB(HerbDmg~Chemo*total+(1|Site/Plant), data=garden2, family=binomial)
clm4r <- glmmTMB(HerbDmg~Chemo*qD+(1|Site/Plant), data=garden2, family=binomial)
clm4s <- glmmTMB(HerbDmg~Chemo*shan+(1|Site/Plant), data=garden2, family=binomial)
clm4 <- glmmTMB(HerbDmg~1+(1|Site/Plant), data=garden2, family=binomial)

bbmle::ICtab(clm4,clm4a,clm4b,clm4c,clm4d,clm4e,clm4f,clm4g,clm4h,clm4i,clm4j,clm4k,clm4l,
             clm4m,clm4n,clm4o,clm4p,clm4q,clm4r,clm4s, weights=T, delta=T, base=T, type="BIC")

clm4_mods <- subset(MuMIn::model.sel(clm4,clm4a,clm4b,clm4c,clm4d,clm4e,clm4f,clm4g,clm4h,clm4i,clm4j,clm4k,clm4l,
                                     clm4m,clm4n,clm4o,clm4p,clm4q,clm4r,clm4s,rank="BIC"), delta <= 4)
sw(clm4_mods)
summary(clm4_mods)

clm_avg <- model.avg(clm4_mods)
summary(clm_avg)
sw(clm4_mods)
summary(clm4_mods)

Anova(clm4d)
Anova(clm4a)
Anova(clm4g)
summary(clm4d)

plot(simulateResiduals(clm4d))
herbsem <- emmeans(clm4a, ~Region, type="response") %>% as.data.frame()

ph1 <- ggplot(garden2)+
  geom_smooth(aes(x=qD, y=HerbDmg), color="black", method="glm", method.args=list(family=binomial))+
  geom_jitter(aes(x=qD, y=HerbDmg, color=Region), height=.025, size=2) +  
  scale_x_continuous("Terpene richness")+
  scale_y_continuous(name='Seed head damage')+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.75, name="Origin")+
  theme_bw(base_size = 20)
  
ph1x <- ggplot(herbsem , aes(x=Region, y=prob, ymin=asymp.LCL, ymax=asymp.UCL, color=Region, fill=Region)) + 
  #geom_point(position=position_jitter(height=0, width=.2),  size=2, color="black", alpha=.25, pch=21 ) + #color="#287D8EFF",
  #geom_boxplot(lwd=1, outlier.shape=NA) +
  geom_pointrange(size=1, lwd=1) +
  scale_x_discrete("")+
  scale_y_continuous(name='')+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1)+
  scale_fill_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.25)+
  theme_bw(base_size=14)+
  theme(legend.position = "none")+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent')) #transparent legend panel)
ph1x


layoutph1 <- c(patchwork::area(t=1, l=1, b=5, r=4),
               patchwork::area(t=1.1, l=4, b=2, r=4))
ph1plot <- ph1+ph1x+plot_layout(design=layoutph1)


dmg1 <- glmmTMB(HerbDmg ~ total + (1|Site), data=garden2, family=binomial)
Anova(dmg1)
r.squaredGLMM(dmg1)


ggplot(chem1 , aes(x=total, y=rich, color=Region)) + geom_point() + geom_smooth(method="lm", formula=y~log(x))+
  scale_color_viridis(discrete = T, option='C', end=.8)+ labs(x="total concentration")

ggplot(garden2, aes(x=HeadSize, fill=as.factor(HerbDmg))) +
  geom_histogram(alpha=.5, position="dodge")+
  theme(legend.position="top")+facet_wrap(~Region)

ggplot(garden2, aes(x=Biomass_g, y=hg_Tallest)) +  geom_jitter(height=.025, size=3, color="grey") + 
  geom_smooth(method="lm")

h1 <- as.data.frame(emmeans(clm4n,  ~ Region, type='response'))

herbdmg <- ggplot() + 
  geom_hline(yintercept = 0, lwd=1, lty=2, color="white")+
  geom_hline(yintercept = 1, lwd=1, lty=2, color="white")+
  geom_jitter(data=garden2, aes(x=Region, y=HerbDmg, color=Region), size=3, height=.025, width=.25, alpha=.5 ) +
  geom_pointrange(data=h1, aes(x=Region, y=prob, ymin=asymp.LCL, ymax=asymp.UCL, color=Region), size=2, lwd=1.5) + 
  #geom_point(data=h1, aes(x=Region, y=prob), size=2, color="white") +
  scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1)+
  scale_x_discrete("Region")+
  scale_y_continuous(labels = scales::percent, name='Seed head damage')+
  theme_bw(base_size = 20)+
  theme(legend.position = "none")

herbdmg
ggsave('herbdmg.tiff',herbdmg, width=7, height=6, units="in", dpi=600)

headplot <- (ph1+herbdmg)+
  plot_annotation(tag_levels = 'a') + plot_layout(guides='collect') & theme(legend.position='right')
ggsave('headplot.tiff',headplot, width=10, height=4, units="in", dpi=600, compression = "lzw")

dmg2 <- glmmTMB(HerbDmg ~ log(total)  + (1|Site), data=garden2  , family=binomial)
Anova(dmg2)
r.squaredGLMM(dmg2)

ph1 <- ggplot(garden2)+
  geom_jitter(aes(x=qD, y=HerbDmg, color=Region), height=.04, size=2) +  
  geom_smooth(aes(x=qD, y=HerbDmg), color="black", method="glm", method.args=list(family=binomial))+
  scale_x_continuous("Terpene richness")+
  scale_y_continuous(name='Seed head damage')+
  scale_color_viridis(discrete=T, option="C", begin=.3, end=.9, direction=-1, alpha=.75, name="Origin")+
  theme_bw(base_size = 20)





# Damage PLOT ####
pdmg1 <- (ph1+pd1plot)/(pa1a+pg3plot) + plot_layout(guides = 'collect')
#pdmg1 <- (pd1a+ph1)/(pg3+pa1a) + plot_layout(guides = 'collect')

ggsave('herbfig_new.tiff',pdmg1, width=12, height=8, units="in", dpi=600, compression='lzw')


ph2 <- ggplot(garden2, aes(x=shan, y=HerbDmg))+
  geom_jitter(height=.025, size=3, alpha=.5) +  
  geom_smooth(method="glm", method.args=list(family=binomial))+
  scale_x_continuous("Shannon diversity")+
  scale_y_continuous(labels = scales::percent, name='Seed head damage')+
  theme_bw(base_size = 20)


plotdmg1 <- ggplot(garden2 %>% filter(Region=="MT"), aes(x=thymol, y=HerbDmg)) +  
  geom_jitter(height=.025, size=3, color="grey") + 
  geom_smooth(method="glm", method.args=list(family='binomial'), lwd=2)+
  scale_x_continuous("log total terpene conc.")+
  scale_y_continuous("Prob. seed head damage")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
plotdmg1

plotdmg2 <- ggplot(garden2 , aes(x=hg_Tallest, y=log(thymol))) +  geom_jitter(height=.025, color="#73D055FF", size=3) + 
  geom_smooth(method="glm", method.args=list(family='binomial'), color="#73D055FF", lwd=2)+
  scale_y_continuous("Prob. seed head damage")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
plotdmg2

dmg3 <- glmmTMB(HerbDmg ~ rich*Chemo + (1|Site), data=garden2, family=binomial)
Anova(dmg3)
r.squaredGLMM(dmg3)

plotdmg2 <- ggplot(garden2 %>% filter(Region=='WI'), aes(x=rich, y=HerbDmg)) +  geom_jitter(height=.025, color="#73D055FF", size=3) + 
  geom_smooth(method="glm", method.args=list(family='binomial'), color="#73D055FF", lwd=2)+
  scale_x_continuous(limits=c(9,25),"terpene richness")+
  scale_y_continuous("Prob. seed head damage")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
plotdmg2

dmg4 <- glmmTMB(HerbDmg ~ simp + (1|Site), data=garden2%>% filter(Region=='WI'), family=binomial)
Anova(dmg4)
r.squaredGLMM(dmg4)

plotdmg3 <- ggplot(garden2 %>% filter(Region=='WI'), aes(x=simp, y=HerbDmg)) +  geom_jitter(height=.025, color="#73D055FF", size=3) + 
  geom_smooth(method="glm", method.args=list(family='binomial'), color="#73D055FF", lwd=2)+
  scale_x_continuous(limits=c(1.4,3.6),"terpene simpson diversity")+
  scale_y_continuous("Prob. seed head damage")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
plotdmg3


dmg5 <- glmmTMB(HerbDmg ~ pc1  + (1|Site), data=garden2 %>% filter(Region=='WI'), family=binomial)
Anova(dmg5)
r.squaredGLMM(dmg5)

plotdmg4 <- ggplot(garden2 %>% filter(Region=='WI'), aes(x=pc2, y=HerbDmg)) +  geom_jitter(height=.025, color="#73D055FF", size=3) + 
  geom_smooth(method="glm", method.args=list(family='binomial'), color="#73D055FF", lwd=2)+
  scale_x_continuous(limits=c(-2.25,3),"terpene PC2")+
  scale_y_continuous("Prob. seed head damage")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
plotdmg4

plotdmg1+plotdmg3+plotdmg4
ggsave('dmgplots2.tiff',plot=plotdmg1+plotdmg3+plotdmg4, width=18, height=6, units="in", dpi=600)


ggplot(fhead, aes(x=HeadSize, fill=as.factor(HerbDmg))) +
  geom_histogram(alpha=.5, position="dodge")+
  theme(legend.position="top")+facet_wrap(~Region)

plotdmg1 <- ggplot(garden2, aes(x=total, y=HerbDmg)) +  
  geom_jitter(height=.025, size=3, color="grey") + 
  geom_smooth(method="glm", method.args=list(family='binomial'), lwd=2)+
  scale_x_continuous("log total terpene conc.")+
  scale_y_continuous("Prob. seed head damage")+
  theme_bw(base_size = 20)
plotdmg1

plotdmg2 <- ggplot(garden2, aes(x=thymol, y=HerbDmg)) +  
  geom_jitter(height=.025, size=3, color="grey") + 
  geom_smooth(method="glm", method.args=list(family='binomial'), lwd=2)+
  scale_x_continuous("log total terpene conc.")+
  scale_y_continuous("Prob. seed head damage")+
  theme_bw(base_size = 20)
plotdmg1


# Q3: Costs of defense ####
library(smatr)
library(lavaan)
library(lavaanPlot)
phytochem1 <- full_join(gard21,chem1) %>% drop_na(Region,total)
head(phytochem1)
chemnew2 <- phytochem1 %>%
  select(Region,Site,Plant,Chemo,Biomass_g,hg_Tallest,volD1, volD2,total,qD) %>% 
  mutate(Volume = pi*((volD1+volD2)/4)*hg_Tallest, totall=log(total+.001), Volumel=log(Volume),
         totalls=scale(totall)[,1], Volumels=scale(Volumel)[,1], qDs=scale(qD)[,1]) %>% 
  group_by(Region,Site,Chemo) %>% 
    summarise_all(mean, na.rm=T)

phytocorr <- phytochem1 %>%
  select(Region,Site,Plant,Chemo,Biomass_g,hg_Tallest,volD1, volD2,total,qD) %>% 
  mutate(Volume = pi*((volD1+volD2)/4)*hg_Tallest, totall=log(total+.001), Volumel=log(Volume),
         totalls=scale(totall)[,1], Volumels=scale(Volumel)[,1], qDs=scale(qD)[,1])
cor.test(phytocorr$Biomass_g, phytocorr$Volume)
cor.test(phytocorr$Biomass_g, phytocorr$hg_Tallest)
cor.test(phytocorr$Volume, phytocorr$hg_Tallest)

lm1 <- glmmTMB(total ~ Volume * Chemo + Volume * Region + (1|Site), data=chemnew2)
Anova(lm1)
emtrends(lm1, var="Volume", ~Chemo, infer=T)

chemoplot1 <- ggplot(chemnew2, aes(x=Biomass_g, y=total, color=Chemo, fill=Chemo)) + 
  geom_point(  size=2 ) + #color="#287D8EFF",
  geom_smooth(method="lm")+
  scale_color_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1)+
  scale_fill_viridis(discrete=T, option="C", begin=.1, end=.9, direction=-1, alpha=.25)+
  scale_x_continuous("log plant volume")+
  scale_y_continuous("Terpene richness")+
  theme_bw(base_size=20)
chemoplot1

## SMATR - build model for Q3 - tradeoffs among siblings across all 12 populations ####

### SMATR for total concentration ####
sma1 <- sma(totall~Biomass_g , data=chemnew2 %>% filter(totall>1))
summary(sma1)
plot(sma1)
plot(sma1, which="residual")
plot(sma1, which="qq")

sma1mt <- sma(totall~Biomass_g, data=chemnew2 %>% filter(Region=="MT"))
summary(sma1mt)
plot(sma1mt)
plot(sma1mt, which="residual")
plot(sma1mt, which="qq")

sma1wi <- sma(totall~Volumel, data=chemnew2 %>% filter(Region=="WI"))
summary(sma1wi)
plot(sma1wi)
plot(sma1wi, which="residual")
plot(sma1wi, which="qq")

sma2 <- sma(totall~Biomass_g*Chemo, data=chemnew2 %>% filter(totall>1), type="shift")
sma2
summary(sma2)
plot(sma2)
plot(sma2, which="residual")
plot(sma2, which="qq")

sma2mt <- sma(totall~Biomass_g*Chemo, data=chemnew2 %>% filter(Region=="MT"), type="shift")
sma2mt
summary(sma2mt)
plot(sma2mt)
plot(sma2mt, which="residual")
plot(sma2mt, which="qq")

sma2wi <- sma(totall~Biomass_g*Chemo, data=chemnew2 %>% filter(Region=="WI"), type="shift")
sma2wi
summary(sma2wi)
plot(sma2wi)
plot(sma2wi, which="residual")
plot(sma2wi, which="qq")

### SMATR for richness ####
smar1 <- sma(qD~Biomass_g, data=chemnew2 )
summary(smar1)
plot(smar1)
plot(smar1, which="residual")
plot(smar1, which="qq")

smar1mt <- sma(qD~Biomass_g, data=chemnew2 %>% filter(Region=="MT"))
summary(smar1mt)
plot(smar1mt)
plot(smar1mt, which="residual")
plot(smar1mt, which="qq")

smar1wi <- sma(qD~hg_Tallest, data=chemnew2 %>% filter(Region=="WI"))
summary(smar1wi)
plot(smar1wi)
plot(smar1wi, which="residual")
plot(smar1wi, which="qq")

smar2 <- sma(qD~Biomass_g*Chemo, data=chemnew2 , type="shift")
smar2
summary(smar2)
plot(smar2)
plot(smar2, which="residual")
plot(smar2, which="qq")

smar2mt <- sma(qD~Volumel*Chemo, data=chemnew2 %>% filter(Region=="MT"), type="shift")
smar2mt
summary(smar2mt)
plot(smar2mt)
plot(smar2mt, which="residual")
plot(smar2mt, which="qq")

smar2wi <- sma(qD~Biomass_g*Chemo, data=chemnew2 %>% filter(Region=="WI"), type="shift")
smar2wi
summary(smar2wi)
plot(smar2wi)
plot(smar2wi, which="residual")
plot(smar2wi, which="qq")

### SMATR for shannon ####
smas1 <- sma(shan~Biomass_g, data=chemnew2 )
summary(smas1)
plot(smas1)
plot(smas1, which="residual")
plot(smas1, which="qq")

smas1mt <- sma(shan~Volumel, data=chemnew2 %>% filter(Region=="MT"))
summary(smas1mt)
plot(smas1mt)
plot(smas1mt, which="residual")
plot(smas1mt, which="qq")

smas1wi <- sma(shan~Volumel, data=chemnew2 %>% filter(Region=="WI"))
summary(smas1wi)
plot(smas1wi)
plot(smas1wi, which="residual")
plot(smas1wi, which="qq")

smas2 <- sma(shan~Biomass_g*Chemo, data=chemnew2 , type="shift")
smas2
summary(smas2)
plot(smas2)
plot(smas2, which="residual")
plot(smas2, which="qq")

smas2mt <- sma(shan~Volumel*Chemo, data=chemnew2 %>% filter(Region=="MT"), type="shift")
smas2mt
summary(smas2mt)
plot(smas2mt)
plot(smas2mt, which="residual")
plot(smas2mt, which="qq")

smas2wi <- sma(shan~Volumel*Chemo, data=chemnew2 %>% filter(Region=="WI"), type="shift")
smas2wi
summary(smas2wi)
plot(smas2wi)
plot(smas2wi, which="residual")
plot(smas2wi, which="qq")




## extract model coefficients to make figure
smap2 <- pop1a
smap2$y1 <- exp(coef(sma1)[1]+smap2$Biomass_g*coef(sma1)[2]) 
smap2$yMT <- exp(coef(sma2)[1,1]+smap2$Biomass_g*coef(sma2)[1,2]) 
smap2$yWI <- exp(coef(sma2)[2,1]+smap2$Biomass_g*coef(sma2)[2,2]) 

## BIG POPS FIG
p2b <- ggplot(pop1a, aes(x=Biomass_g,y=total)) + 
  #geom_smooth(data=smap2, aes(x=Biomass_g,y=totall), method="lm", color="black", size=1.5, se=T)+
  geom_line(data=smap2 , aes(x=Biomass_g,y=y1), color="black", size=1.5)+
  geom_line(data=smap2 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1.25, lty=2)+
  geom_line(data=smap2 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1.25, lty=2)+
  geom_point(data=pop1a, aes(x=Biomass_g,y=total,color=Region), size=2, pch=15) +
  annotate(geom='text', x=.95,y=48, label="B) within populations", fontface='bold', size=5)+
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Biomass (mg)", limits = c(0.3,2.6), breaks=seq(.5,2.5,.5))+
  scale_y_continuous("Total terpenes (mg/g)", limits = c(5,49), breaks=seq(0,43,10))+
  theme(legend.key=element_blank(), legend.position = 'none')+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))


### SMATR PLOT ####
cost1 <- ggplot()+
  geom_smooth(data=chemnew2, aes(x=Biomass_g, y=totall), color="black", method="lm")+
  geom_point(data=chemnew2, aes(x=Biomass_g, y=totall, color=Chemo), size=4)+
  scale_x_continuous("Plant biomass (g)")+
  scale_y_continuous(name='Terpene conc. (log mg/g)')+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  theme_bw(base_size = 20)

cost2 <- ggplot()+
  geom_smooth(data=chemnew2, aes(x=Biomass_g, y=qD), color="black", method="lm")+
  geom_point(data=chemnew2, aes(x=Biomass_g, y=qD, color=Chemo), size=4)+
  scale_x_continuous("Plant biomass (g)")+
  scale_y_continuous(name='Terpene richness')+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  theme_bw(base_size = 20)

cost3 <- ggplot(chemnew2, aes(x=Biomass_g, y=shan, color=Chemo))+
  geom_smooth(method="lm")+
  geom_point(size=4)+
  scale_x_continuous("Plant biomass (g)")+
  scale_y_continuous(name='Terpene shannon')+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  theme_bw(base_size = 20)+
  guides(color="none")

costplot <- (cost1+cost2+cost3)+
  plot_annotation(tag_levels = 'a') + plot_layout(guides='collect') & theme(legend.position='right')
ggsave('costplot.tiff',costplot, width=12, height=3.75, units="in", dpi=600, compression = "lzw")

mcost2 <- glmmTMB(qD~Biomass_g * Chemo + (1|Site), data=chemnew2)
Anova(mcost2)


#####################################################################################
# OLD CODE ########################################################################
#### correlation plots
corr1 <- ggplot(garden2 , aes(x=log(total), y=rich, color=Pop)) +  geom_jitter(height=.025, size=3) + 
  geom_smooth(method="lm", lwd=2, se=F)+
  scale_x_continuous("log total terpene conc.")+
  scale_y_continuous("terpene richness")+
  scale_color_viridis(discrete = T, option='D') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
corr1

corr2 <- ggplot(garden2 , aes(x=log(total), y=simp, color=Pop)) +  geom_jitter(height=.025, size=3) + 
  geom_smooth(method="lm", lwd=2, se=F)+
  #scale_x_continuous(limits=c(-2.25,3),"terpene PC2")+
  #scale_y_continuous("Prob. seed head damage")+
  scale_color_viridis(discrete = T, option='D') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
corr2

corr3 <- ggplot(garden2 , aes(x=rich, y=simp, color=Pop)) +  geom_jitter(height=.025, size=3) + 
  geom_smooth(method="lm", lwd=2, se=F)+
  #scale_x_continuous(limits=c(-2.25,3),"terpene PC2")+
  #scale_y_continuous("Prob. seed head damage")+
  scale_color_viridis(discrete = T, option='D') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
corr3

corrplot <- corr1+corr2+corr3 + plot_layout(guides='collect') & theme(legend.position='none') 
ggsave('corrplots.tiff',plot=corrplot, width=18, height=6, units="in", dpi=600)

##############################################################################################################################
### SEED MODELS ##############################################################################################################
garden2$HerbDmg1 <- as.factor(garden2$HerbDmg)
dmg1a <- glmmTMB(Seeds_No ~ log(total)*HerbDmg1+I(log(total)^2)*HerbDmg1+HeadSize + (1|Site), data=garden2, family=nbinom2)
summary(dmg1a)
Anova(dmg1a)
r.squaredGLMM(dmg1a)

emtrends(dmg1a, ~HerbDmg1, var="rich")
ggemmeans(dmg1a, terms = ~log(total)*HerbDmg1+I(log(total)^2)*HerbDmg1+HeadSize)


dmg2a <- glmmTMB(Seeds_No ~ log(total)*HerbDmg1+I(log(total)^2)*HerbDmg1+HeadSize + (1|Site), data=garden2, family=nbinom2)
summary(dmg2a)
Anova(dmg2a)
r.squaredGLMM(dmg2a)

sim_glm2 <- simulateResiduals(fittedModel = dmg2a, n = 250)
plot(sim_dmg2a) # residuals look good

emmeans(dmg2a, ~HerbDmg1, type='response', at=list(rich=12))
emmeans(dmg2a, ~HerbDmg1, type='response', at=list(rich=18))
emtrends(dmg2a, ~HerbDmg1, var="rich")
ggpredict(dmg2a, terms = ~rich*HerbDmg1+HeadSize)

plotseed1 <- ggplot(garden2 %>% filter(Region=='WI') %>% drop_na(Seeds_No), aes(x=log(total), y=Seeds_No, color=HerbDmg1)) +  
  geom_point(size=3) + 
  geom_smooth(method="lm", formula=y~x+I(x^2), lwd=2)+
  #geom_smooth(method="glm", method.args=list(family='nbinom1'), formula=y~x+I(x^2), lwd=2)+
  scale_color_viridis(discrete=T, option="D", begin=.5, end=.85, direction=-1, name="Head damage", labels=c("No","Yes"))+
  scale_x_continuous("log terpene total conc.")+
  scale_y_continuous("Number seeds produced")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
plotseed1


dmg2a <- glmmTMB(Seeds_No ~ rich*HerbDmg1+HeadSize + (1|Site), data=garden2, family=nbinom2)
summary(dmg2a)
Anova(dmg2a)
r.squaredGLMM(dmg2a)

plotseed2 <- ggplot(garden2 %>% drop_na(Seeds_No), aes(x=rich, y=Seeds_No, color=HerbDmg1)) +  
  geom_point(size=3) + 
  geom_smooth(method="glm", method.args=list(family='poisson'), formula=y~x+I(x^2), lwd=2)+
  scale_color_viridis(discrete=T, option="D", begin=.5, end=.85, direction=-1, name="Head damage", labels=c("No","Yes"))+
  scale_x_continuous("terpene richness")+
  scale_y_continuous("Number seeds produced")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
plotseed2

dmg3a <- glmmTMB(Seeds_No ~ pc2+HerbDmg1 + (1|Site), data=garden2, family=nbinom2)
summary(dmg3a)
Anova(dmg3a)
r.squaredGLMM(dmg3a)

plotdmg3 <- ggplot(garden2 %>% drop_na(Seeds_No), aes(x=pc2, y=Seeds_No, color=HerbDmg1)) +  
  geom_point(size=3) + 
  geom_smooth(method="glm", method.args=list(family='nbinom2'), formula=y~x+I(x^2), lwd=2)+
  scale_color_viridis(discrete=T, option="D", begin=.5, end=.85, direction=-1, name="Head damage", labels=c("No","Yes"))+
  scale_x_continuous("terpene richness")+
  scale_y_continuous("Number seeds produced")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=20),
        axis.ticks = element_line(color="white", size=1.5))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=24, color="white")) +
  theme(legend.position = "top", legend.background = element_rect(fill="black", color="white"))+
  theme(plot.background = element_rect(fill="black"))+
  theme(panel.background = element_rect(fill="black"))+
  theme(panel.border = element_rect(color="white", fill=NA, size=2)) +
  theme(text = element_text(size=28, color="white"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="white", size=24))
plotdmg3

## old figures
ggplot(garden2, aes(x=Region, y=HerbDmg)) + geom_violin() + geom_jitter(width=.2,height=.05) + geom_smooth(method="lm")
ggplot(garden2 , aes(x=log(total), y=HerbDmg)) +  geom_jitter(height=.05) + geom_smooth(method="glm", method.args=list(family='binomial'))
ggplot(garden2 , aes(x=rich, y=HerbDmg)) + geom_jitter(height=.05) + geom_smooth(method="glm", method.args=list(family='binomial'))
ggplot(garden2, aes(x=simp, y=HerbDmg)) +  geom_jitter(height=.05) + geom_smooth(method="glm", method.args=list(family='binomial'))
ggplot(garden2, aes(x=pc2, y=HerbDmg)) +  geom_jitter(height=.05) + geom_smooth(method="glm", method.args=list(family='binomial'))
ggplot(garden2, aes(x=HeadSize, y=HerbDmg)) +  geom_jitter(height=.05) + geom_smooth(method="glm", method.args=list(family='binomial'))

ggplot(garden2, aes(x=Region, y=piercing)) + geom_violin() + geom_jitter(width=.2,height=0) + geom_smooth(method="lm")
ggplot(garden2 , aes(x=log(total), y=piercing)) +  geom_jitter(height=0) + geom_smooth(method="glm", method.args=list(family='beta_family'))
ggplot(garden2 , aes(x=rich, y=piercing)) + geom_jitter(height=0) + geom_smooth(method="glm", method.args=list(family='beta_family'))
ggplot(garden2, aes(x=shan, y=piercing)) +  geom_jitter(height=0) + geom_smooth(method="glm", method.args=list(family='beta_family'))
ggplot(garden2, aes(x=simp, y=piercing)) +  geom_jitter(height=0) + geom_smooth(method="glm", method.args=list(family='beta_family'))

ggplot(garden2, aes(x=Region, y=Seed_Loss )) + geom_violin() + geom_jitter(width=.2,height=0) + geom_smooth(method="lm")
ggplot(garden2 , aes(x=log(total), y=PercSeedLoss )) +  geom_jitter(height=0) + geom_smooth(method="glm", method.args=list(family='beta_family'))
ggplot(garden2 , aes(x=rich, y=Seed_Loss )) + geom_jitter(height=0) + geom_smooth(method="glm", method.args=list(family='beta_family'))
ggplot(garden2%>% filter(HerbDmg==1), aes(x=shan, y=Seed_Loss )) +  geom_jitter(height=0) + geom_smooth(method="glm", method.args=list(family='beta_family'))
ggplot(garden2%>% filter(HerbDmg==1), aes(x=pc1, y=PercSeedLoss1)) +  geom_jitter(height=0) + geom_smooth(method="glm", method.args=list(family='beta_family'))

glm1 <- glmmTMB(HerbDmg ~ Region + (1|Site), data=garden2 %>% filter(rich>3), family=binomial)
summary(glm1)
Anova(glm1)
emmeans(glm1, ~Region, type="response")

glm2 <- glmmTMB(HerbDmg ~ log(total) , data=garden2 , family=binomial)
summary(glm2)

glm3 <- glmmTMB(HerbDmg ~ rich, data=garden2 , family=binomial)
summary(glm3)
Anova(glm3)

glm4 <- glmmTMB(HerbDmg ~ simp , data=garden2 , family=binomial)
summary(glm4)

glm4 <- glmmTMB(HerbDmg ~ pc1, data=garden2 , family=binomial)
summary(glm4)

glm1p <- glmmTMB(piercing ~ Region + (1|Site), data=garden2 %>% filter(rich>3), family=beta_family)
summary(glm1)
Anova(glm1)
emmeans(glm1, ~Region, type="response")

glm2 <- glmmTMB(piercing ~ log(total) + (1|Site), data=garden2 , family=beta_family)
summary(glm2)

glm3 <- glmmTMB(piercing ~ rich+ (1|Site), data=garden2 , family=beta_family)
summary(glm3)

glm4 <- glmmTMB(piercing ~ simp, data=garden2 , family=beta_family)
summary(glm4)

glm4 <- glmmTMB(HerbDmg ~ rich, data=garden2 , family=binomial)
summary(glm4)

### sem ###

mod1 <- ' HerbDmg ~ total + rich + simp
          rich ~~ total
          simp ~~ total
          rich ~~ simp'

mod1.fit <- sem(mod1, data=garden2)
summary(mod1.fit, rsq=T)
summary(mod1.fit, fit.measures=T)
summary(mod1.fit, modindices=T)

standardizedSolution(mod1.fit, type="std.all")

mod2 <- ' PropDmg ~ rich + BrayDist
          rich ~ map12
          BrayDist ~ map12'

mod2.fit <- sem(mod2, data=dmg4)
summary(mod2.fit, rsq=T)
summary(mod2.fit, fit.measures=T)
summary(mod2.fit, modindices=T)

standardizedSolution(mod2.fit, type="std.all")


lavaanPlot(model = mod2.fit, node_options = list(fontname="Arial"))


########################################################################################################################
#### summarize by region ###############################################################################################
chemhill <- chem1 %>% drop_na(Region) %>%  select(Region,Site,Plant,Rep,rich,shan,simp) %>%  group_by(Region,Site) %>% 
  #mutate(rich1=scale(rich), shan1=scale(shan), simp1=scale(simp)) %>% 
  summarize(rich=mean(rich, na.rm=T), shan=mean(shan, na.rm=T), simp=mean(simp, na.rm=T)) %>% 
  gather(Hill_no,value=value,rich:simp) 

ggplot(chemhill, aes(x=Region, y=value)) + geom_boxplot() + geom_point() + facet_wrap(~Hill_no, scales = "free_y")

lm1 <- lmer(value~Region*Hill_no+(1|Site), data=chemhill)
anova(lm1)
emmeans(lm1, pairwise~Region|Hill_no)


### chem beta diversity ###
chem2 <- trait %>% select(Region,Site,Plant,Rep,a_thujene:caryophyllene) %>% na.omit()
regsite <- chem2 %>% select(Region,Site)
bdivmat <- as.matrix(vegdist(chem2[5:32], method="bray"))

plantrep <- trait %>% select(Region,Site,Plant,Rep,a_thujene:caryophyllene) %>% na.omit() %>% select(Plant,Rep)
PlantID <- paste(plantrep$Rep,plantrep$Plant,sep="")
regsite$PlantID <- paste(plantrep$Rep,plantrep$Plant,sep="")
regsite$PlantID1 <- paste(plantrep$Rep,plantrep$Plant,sep="")

row.names(bdivmat) <- PlantID
colnames(bdivmat) <- PlantID

bdivmat1 <- cbind(PlantID,data.frame(bdivmat,row.names = NULL))

pdist1 <- bdivmat1 %>% as.data.frame() %>% pivot_longer(cols=2:324, names_to = "PlantID1", values_to = "BrayDist")
#colnames(pdist1) <- c("fsp","nsp","PhyloDist")

j1 <- regsite %>% select(Site,PlantID)
j2 <- regsite %>% select(Site,PlantID1)

pdist2 <- full_join(j1,pdist1, by="PlantID")
pdist3 <- full_join(j2,pdist2, by="PlantID1")

pdist3$match <- if_else(pdist3$Site.x==pdist3$Site.y, 1, 0) 

pdist <- pdist3 %>% filter(match==1, BrayDist>0) %>% mutate(Site=Site.x) %>% group_by(Site) %>% summarize(BrayDist=mean(BrayDist))


## process 2018 seed data
dmg18 <- full_join(sites,dmg19, by="Site")

dmg18$dmg18heads <- dmg18$Nlep+dmg18$Ncol+dmg18$Nunk
dmg18$PropDmg <- dmg18$dmg18heads/dmg18$Totheads
dmg18$Phyr <- dmg18$Nlep/dmg18$Totheads
dmg18$Lep <- dmg18$Ncol/dmg18$Totheads
dmg18$Lf_dmg18 <- rowMeans(dmg18[19:23], na.rm=TRUE)
head(dmg18, 10)
hist(dmg18$All_Insects)
hist(dmg18$Lepidopterans)
hist(dmg18$PropDmg)
hist(sqrt(dmg18$Lf_dmg18))

dmg18s <- dmg18 %>% group_by(Region,Site) %>% 
  summarize(PropDmg=mean(PropDmg, na.rm=T), Lf_dmg=mean(Lf_dmg18,na.rm=T), mat1=mean(mat1,na.rm=T), tsd4=mean(tsd4,na.rm=T), map12=mean(map12,na.rm=T)) %>% 
  mutate(Lf_dmg=(Lf_dmg/100))

trait$ChemoT <- if_else(trait$Chemo=='T',1,0)

chemhill1 <- trait %>% drop_na(Region) %>%  select(Region,Site,Plant,ChemoT,total,totTC,carvacrol,thymol,rich,shan,simp) %>%  group_by(Region,Site) %>% 
  #mutate(rich1=scale(rich), shan1=scale(shan), simp1=scale(simp)) %>% 
  summarize(rich=mean(rich, na.rm=T), shan=mean(shan, na.rm=T), simp=mean(simp, na.rm=T), total=mean(log(total), na.rm=T), totTC=mean(log(totTC), na.rm=T), 
            ChemoT=mean(ChemoT, na.rm=T))

dmg1 <- full_join(chemhill, dmg18s, by=c("Region","Site"))
dmg2 <- full_join(chemhill1, dmg18s, by=c("Region","Site"))
dmg2$PropDmg1 <- (dmg2$PropDmg*(length(dmg2$PropDmg)-1)+.5)/length(dmg2$PropDmg)
dmg2$Lf_dmg1 <- (dmg2$Lf_dmg*(length(dmg2$Lf_dmg)-1)+.5)/length(dmg2$Lf_dmg)
dmg3 <- full_join(dmg2,pdist)
dmg4 <- na.omit(dmg3)

ggplot(dmg1, aes(x=value, y=PropDmg, color=Region)) + geom_point() + geom_smooth(method="lm") + facet_wrap(~Hill_no, scales="free_x")
ggplot(dmg1, aes(x=value, y=Lf_dmg, color=Region)) + geom_point() + geom_smooth(method="lm") + facet_wrap(~Hill_no, scales="free_x")
ggplot(dmg4, aes(x=BrayDist, y=PropDmg)) + geom_point() + geom_smooth(method="lm")
ggplot(dmg4, aes(x=BrayDist, y=Lf_dmg)) + geom_point() + geom_smooth(method="lm")

### leaf damage plots ####
ggplot(dmg4, aes(x=Region, y=PropDmg)) + geom_boxplot() + geom_smooth(method="lm")
ggplot(dmg4, aes(x=total, y=PropDmg)) + geom_point() + geom_smooth(method="lm")
ggplot(dmg4, aes(x=rich, y=PropDmg)) + geom_point() + geom_smooth(method="lm")
ggplot(dmg4, aes(x=BrayDist, y=PropDmg)) + geom_point() + geom_smooth(method="lm")

ggplot(dmg4, aes(x=Region, y=PropDmg)) + geom_boxplot() + geom_smooth(method="lm")
ggplot(dmg4, aes(y=total, x=PropDmg)) + geom_point() + geom_smooth(method="lm")
ggplot(dmg4, aes(y=rich, x=PropDmg)) + geom_point() + geom_smooth(method="lm")
ggplot(dmg4, aes(y=BrayDist, x=PropDmg, color=Region)) + geom_point() + geom_smooth(method="lm")

ggplot(dmg4, aes(x=Region, y=Lf_dmg)) + geom_boxplot() + geom_smooth(method="lm")
ggplot(dmg4, aes(x=total, y=Lf_dmg)) + geom_point() + geom_smooth(method="lm")
ggplot(dmg4, aes(x=rich, y=Lf_dmg)) + geom_point() + geom_smooth(method="lm")
ggplot(dmg4, aes(x=BrayDist, y=Lf_dmg)) + geom_point() + geom_smooth(method="lm")


glm1 <- glmmTMB(PropDmg1 ~ rich, data=dmg4, family=beta_family)
summary(glm1)

glm1 <- glmmTMB(PropDmg1 ~ rich, data=dmg4, family=gaussian)
summary(glm1)

glm2 <- betareg(PropDmg1 ~ Region*BrayDist, data=dmg4, link="logit")
summary(glm2)
Anova(glm2, type=2)

dmg5 <- dmg4 %>% group_by(Region,Site) %>% pivot_longer(cols=c("total","rich","simp","BrayDist","ChemoT","mat1","map12"), names_to = "var", values_to = "value")

ggplot(dmg5, aes(x=value, y=PropDmg)) + geom_point() + geom_smooth(method="lm") + facet_wrap(~var, scales="free_x")

ggplot(dmg5 %>% filter(var!='map12',var!='mat1'), aes(x=Region, y=value)) + geom_boxplot() + geom_smooth(method="lm") + facet_wrap(~var, scales="free_y")

lm1 <- lm(value ~ Region, data=dmg5 %>% filter(var=='ChemoT'))
anova(lm1)

emmeans(lmer1, pairwise~Region|var)

### sem ###

mod1 <- ' PropDmg ~ total + simp + BrayDist
          total ~ map12
          simp ~ map12
          BrayDist ~ map12'

mod1.fit <- sem(mod1, data=dmg4)
summary(mod1.fit, rsq=T)
summary(mod1.fit, fit.measures=T)
summary(mod1.fit, modindices=T)

standardizedSolution(mod1.fit, type="std.all")

mod2 <- ' PropDmg ~ rich + BrayDist
          rich ~ map12
          BrayDist ~ map12'

mod2.fit <- sem(mod2, data=dmg4)
summary(mod2.fit, rsq=T)
summary(mod2.fit, fit.measures=T)
summary(mod2.fit, modindices=T)

standardizedSolution(mod2.fit, type="std.all")


lavaanPlot(model = mod2.fit, node_options = list(fontname="Arial"))

#############################################################################################################################################
######### OLD CODE ##########################################################################################################################
#############################################################################################################################################
s1 <- full_join(site,seed,by="Site")
s1a <- s1 %>% select(Region,Site,Plant)

s2 <- right_join(s1a,size,by="Plant")
s3 <- left_join(s2,gard,by=c("Plant","Rep"))
s3a <- left_join(s3,leaf,by=c("Plant","Rep"))
s3a$sla <-  (s3a$Surface_Area/s3a$Leaf_Mass_mg)
s3a$trics <- rowMeans(s3a[,26:27])
s3a$trics_tot <- (s3a$trics/s3a$Surface_Area)
head(s3a)

s4a <- s3a %>% select(Region,Site,Plant,Rep,Row,OrderGard,G_Date,Height2,Stems,Stems2,Stems_new,sla,trics,trics_tot)
s4 <- s4a[!is.na(s4a$Height2),]

# write.csv(s4, 'MonardaGardenPlants_30Sept2019.csv')

### load biomass data
bm <- read_csv("Monarda_CommonGarden_Biomass_2019Final.csv")

s5 <- left_join(s4,bm,by=c("Plant","Rep"))
s5$codedup <- paste(s5$Plant,s5$Rep,sep="")
duplicated(s5$codedup)
head(s5)

s6 <- s5 %>% select(Region,Site,Plant,Rep,Row,G_Date,OrderGard,Height2,Stems,Stems2,Stems_new,Biomass_g,Labeling_errors)
# write.csv(s6, 'MonardaGardenPlants_30Sept2019.csv')

#####################################################################
#### SIZE ANALYSIS ##################################################

s6$rgr <- (log((s6$Stems_new+s6$Stems2)-s6$Stems2))/44
s6$rgr1 <- (log(s6$Stems2-s6$Stems))/36
s6$rgr[is.infinite(s6$rgr)] <- 0
s6$rgr1[is.infinite(s6$rgr1)] <- 0

ggplot(s6, aes(x=rgr, y=Biomass_g, color=Region))+ geom_point() + geom_smooth(method='lm') + theme_bw()



gr1 <- lmer(Biomass_g~rgr*Region+(1|Site/Plant), data=s6 )
summary(gr1)
anova(gr1)
r.squaredGLMM(gr1)


size1 <- s6 %>% select(Region,Site,Plant,Rep,Row,G_Date,Height2,Stems,Stems2,Stems_new,Biomass_g) %>% 
  group_by(Region,Site,Plant,Rep,Row,Height2,Biomass_g) %>% 
  pivot_longer(cols=c("Stems","Stems2","Stems_new"), names_to = "Time", values_to = "N_stems")
size1$M_Date <- if_else(size1$Time=='Stems', '4/15/2019', size1$Time)
size1$M_Date <- if_else(size1$Time=='Stems2', '5/21/2019', size1$M_Date)
size1$M_Date <- if_else(size1$Time=='Stems_new', '7/8/2019', size1$M_Date)
size1$G_Date <- mdy(size1$G_Date)
size1$M_Date <- mdy(size1$M_Date)
size1$Age <- (size1$M_Date-size1$G_Date)
size1$ID <- paste(size1$Plant,size1$Rep,sep="")

ggplot(size1, aes(x=M_Date,y=N_stems, color=Region))+geom_point(position=position_jitter(width=2,height=.25)) + geom_smooth() + theme_bw()

ggplot(size1 %>% filter(M_Date=='2019-07-08'), aes(x=N_stems, y=Biomass_g, color=Region))+ geom_point() + geom_smooth(method='lm')

smod1 <- lmer(Biomass_g~N_stems*Region+(1|Site/Plant), data=size1 %>% filter(M_Date=='2019-07-08'))
summary(smod1)
anova(smod1)
r.squaredGLMM(smod1)

######################################################################
### SITE AND PRODUCTIVITY ANALYSIS ############################################
site1 <- site %>% filter(Use==1)

## env PCA
spc <- prcomp(site1[c("orgmat","nit","cec","pH","P","mat1","pwrq18")], scale=TRUE, center=TRUE)
summary(spc)
spc
plot(spc)
biplot(spc, pc.biplot=TRUE)

site1$pc1 <-  predict(spc)[,1]
site1$pc2 <-  predict(spc)[,2]


pcrot <- as.data.frame(spc$rotation[,1:2], stringsAsFactors=TRUE)
pcrot$Envvar <- c("%OM","soil_N","CEC","pH","soil_P","temp","precip")
pcrot$PC1s <- pcrot$PC1 #scale loadings for plotting as vector arrows
pcrot$PC2s <- pcrot$PC2 #scale loadings for plotting as vector arrows


an1 <- lm(pc1~Region, data=site1)
anova(an1)
an2 <- lm(pc2~Region, data=site1)
anova(an2)

### map ###############
library(rgdal)
library(raster)
library(RColorBrewer)
library("maps")
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(rasterVis)


bounds1 <- readOGR(dsn=getwd(), layer= "cb_2014_us_state_500k") 
clipxy <- rbind(c(-117,42),c(-117,49),c(-85,49),c(-85,42))
bounds <- crop(bounds1,clipxy)
#bounds$STATEFP <- as.numeric(bounds$STATEFP)
p <- raster('wc2.0_bio_2.5m_18.tif')


prec1 <- crop(p,clipxy)
#prec <- rasterize(p_states,prec1)
#prec1 <- crop(clim18[[18]],p_states)
#prec <- rasterize(p_states,prec1)

gplot(prec1)+geom_raster(aes(fill=value))+
  scale_fill_gradientn(colours=c('brown','yellow','dodgerblue','blue'))+
  geom_polygon(data = bounds, aes(x=long, y=lat), color="black", fill=NA) + 
  geom_point(data=site1, aes(x=Long, y=Lat), color="red")+
  #coord_sf(xlim = c(-117, -85), ylim = c(42,49), expand = FALSE) +
  coord_equal()


plot(prec1, xlab="Longitude", ylab="Latitude")
plot(bounds, add=TRUE)


## full map

breakpoints <- seq(0,400,length.out=9)
col1 <- colorRampPalette(c('brown','yellow','dodgerblue','blue'))
color_levels=9


tiff("WIMT_2018.tiff", width=6, height=4, units='in',
     res=300, compression='lzw')
plot(prec1, xlab="Longitude", ylab="Latitude",
     breaks=breakpoints,col=col1(n=color_levels),cex.lab=1.2)
plot(bounds, add=TRUE)
pp <- cbind(site1$Long,site1$Lat)
pp1 <- cbind(site1$Long[10],site1$Lat[10])
ppm <- cbind(site1$Long[8],site1$Lat[8])
ppw <- cbind(site1$Long[12],site1$Lat[12])
points(pp, pch=21, col="black", bg="white", cex=1)
points(ppm, pch=22, col="black", bg="white", cex=1)
points(ppw, pch=22, col="black", bg="white", cex=1)
points(pp1, pch=24, bg="white", col=NA, cex=1)
points(pp1, pch=25, bg="white", col=NA, cex=1)
text(-116, 51, "A", font=2)
dev.off()

detach(package:rasterVis)
detach(package:raster)
detach(package:tidyverse)
library(tidyverse)

### PCA FIGURE ###
fig2B <- ggplot() + xlab("pc1 (45.1%)")+ylab("pc2 (27.9%)")+ 
  stat_ellipse(data=site1, geom="polygon" ,aes(x=pc1,y=pc2, color=Region, fill=Region), lwd=1, alpha=.1)+
  geom_point(data=site1 %>% filter(Site!='UWM',Site!='BEN', Site!='MS'),aes(x=pc1, y=pc2, color=Region), size=5,stroke=1.25) +
  geom_point(data=site1 %>% filter(Site=='MS'),aes(x=pc1, y=pc2), pch=22, size=5, color="darkgoldenrod1", fill="darkgoldenrod1") +
  geom_point(data=site1 %>% filter(Site=='BEN'),aes(x=pc1, y=pc2), pch=22, size=5, color="dodgerblue", fill="dodgerblue") +
  geom_point(data=site1 %>% filter(Site=='UWM'),aes(x=pc1, y=pc2), pch=24, size=5, color="dodgerblue", fill="dodgerblue") +
  geom_point(data=site1 %>% filter(Site=='UWM'),aes(x=pc1, y=pc2), pch=25, size=5, color="dodgerblue", fill="dodgerblue") +
  annotate(geom='text', x=-4.5,y=5, label="B", fontface='bold', size=8)+
  scale_fill_manual(values=c("darkgoldenrod1","dodgerblue")) +
  scale_color_manual(values=c("darkgoldenrod1","dodgerblue")) +
  geom_segment(data=pcrot, aes(x=0, y=0,xend=PC1*5,yend=PC2*5), arrow=arrow(length=unit(0.03, "npc")), color='grey', lwd=1)+
  #geom_text(data=pcrot, aes(x=PC1*5.75, y=PC2*5.75, label=rownames(pcrot[1:7,])), fontface="bold")+
  geom_text(data=pcrot, aes(x=PC1*5.75, y=PC2*5.75, label=Envvar), fontface="bold")+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = c(.9,.85), legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))+
  guides(lty=F, lty=guide_legend(override.aes = list(lty=NA)))

## productivity
prod1 <- read_csv("Hahn_WIMT_Productivity_07242018.csv")
prod2 <- full_join(site1,prod1,by="Site") %>% filter(Use==1)
prod2$veghgt <- (prod2$veghgt1+prod2$veghgt2+prod2$veghgt3+prod2$veghgt4)/4;
prod2$light_reduce <- prod2$light_abv-prod2$light_blw;
  
head(prod2)
ggplot(prod2, aes(x=pc1, y=veghgt))+geom_point()+facet_wrap(~Site)
ggplot(prod2, aes(x=veghgt, y=light_reduce))+geom_point()+facet_wrap(~Site)

lgt1 <- lm(light_reduce ~ veghgt, data=prod2)
summary(lgt1)

prodm <- prod2 %>% group_by(Region,Site) %>% summarize(pc1=mean(pc1), veghgt=mean(veghgt, na.rm=T), mat1=mean(mat1,na.rm=T), pwrq18=mean(pwrq18,na.rm=T),light_reduce=(mean(light_reduce,na.rm=T)))

veg1 <- lm(log(veghgt)~pc1, data=prodm)
v1 <- summary(veg1)

pc1x <- seq(min(prodm$pc1),max(prodm$pc1),.1)
vegy <- exp(predict(veg1,list(pc1=pc1x)))

fig2C <- ggplot() + xlab("PC1")+ylab("mean veg height (cm)")+ geom_line(data=, aes(x=pc1x, y=vegy), lwd=1.5) +
  geom_point(data=prodm %>% filter(Site!='MS',Site!='BEN'), aes(x=pc1, y=veghgt, color=Region), size=5)+
  geom_point(data=prodm %>% filter(Site=='MS'), aes(x=pc1, y=veghgt, color=Region), size=5, pch=22, color="darkgoldenrod1", fill="darkgoldenrod1")+
  geom_point(data=prodm %>% filter(Site=='BEN'), aes(x=pc1, y=veghgt, color=Region), size=5, pch=22, color="dodgerblue", fill="dodgerblue")+
  geom_point(data=prodm %>% filter(Site=='UWM'), aes(x=pc1, y=veghgt, color=Region), size=5, pch=24, color="dodgerblue", fill="dodgerblue")+
  geom_point(data=prodm %>% filter(Site=='UWM'), aes(x=pc1, y=veghgt, color=Region), size=5, pch=25, color="dodgerblue", fill="dodgerblue")+
  scale_color_manual(values=c("darkgoldenrod1","dodgerblue")) +
  annotate(geom='text', x=-3,y=61, label="C", fontface='bold', size=8)+
  #annotate(geom='text', x=-1.5,y=50, label="italic(R)^2 == 0.84", parse=T, size=6)+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = c(.9,.15), legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))+
  guides(lty=F, lty=guide_legend(override.aes = list(lty=NA)))

tiff("Fig2BC.tiff", width=10, height=4, units='in',
     res=300, compression='lzw')
fig2B+fig2C
dev.off()


##########################################################################
#### ANALYSIS CHEM #######################################################
chem <- read_csv("Monarda_ChemMTWI_07262020.csv")
s5a <- full_join(s5,chem,by=c("Plant","Rep"))
s5a$totTC <- s5a$thymol+s5a$carvacrol

s5a1 <- s5a %>% select(Plant,Rep,Region,Site,Chemo,Biomass_g,trics,total,totTC,thymol,carvacrol)

ggplot(s5a, aes(x=Region, y=Biomass_g)) + geom_boxplot()
ggplot(s5a, aes(x=Region, y=trics)) + geom_boxplot()
ggplot(s5a, aes(x=Region, y=total)) + geom_boxplot()
ggplot(s5a, aes(x=Region, y=totTC)) + geom_boxplot()

ggplot(s5a, aes(x=trics, y=(total/trics), color=Region))+ geom_point()
ggplot(s5a, aes(x=trics, y=carvacrol, color=Region))+ geom_point()



#### trait correlations
s51 <- s5a %>% select(Region,Site,Plant,Rep,Chemo,Biomass_g,trics,sla,Stems2,Stems_new,total,totTC,thymol,carvacrol) %>% drop_na(Region) %>% 
  mutate(Biomass_sq=sqrt(Biomass_g),trics_sq=sqrt(trics),totall=log(total), totTCl=log(totTC))
s5p <- s51 %>% group_by(Region,Site) %>% summarize(sla=mean(sla, na.rm=T), trics_sq=mean(trics_sq, na.rm=T), trics=mean(trics, na.rm=T), 
                                                   Biomass_g=mean(Biomass_g, na.rm=T),total=mean(total, na.rm=T), totTC=mean(totTC,na.rm=T),
                                                   thymol=mean(thymol, na.rm=T), carvacrol=mean(carvacrol, na.rm=T), totall=mean(totall,na.rm=T),totTCl=mean(totTCl, na.rm=T))

hist(s51$sla)
ggplot(s5p, aes(x=trics, y=totall, color=Region))+ geom_point() + geom_smooth(method='lm')
ggplot(s51, aes(x=Stems_new, y=totTC, color=Region))+ geom_point()

s5b <- s51 %>% select(Region,Site,Plant,Rep,Biomass_g,sla,trics,totall) %>% group_by(Region,Site,Plant,Rep) %>% gather(key="trait",value="value", -c(Region:Rep))

ggplot(s5b, aes(x=Region, y=value))+ geom_boxplot() + facet_wrap(~trait, scales = 'free_y')

# write.csv(s5a, "Monarda_Garden_Traits_2019.csv")

## prod
s5prod <- full_join(s5p,prodm)
ggplot(s5prod, aes(x=veghgt, y=totall, color=Region)) + geom_point() + geom_smooth(method="lm")

# biomass ~ trics
#### SUMMARIES OF SAMPLES SIZES
ss1 <- s51 %>% mutate(num=1)
ss1a <- ss1 %>% group_by(Region,Site,Plant) %>% summarize(num1=mean(num,na.rm=T)) %>% group_by(Region,Site) %>% summarize(n_moms=sum(num1,na.rm=T))
mean(ss1a$n_moms)
ss1b <- ss1 %>% group_by(Region,Site,Plant) %>% summarize(num1=sum(num,na.rm=T)) %>% group_by(Region,Site) %>% summarize(n_sibs=mean(num1,na.rm=T))
ss1c <- ss1 %>% group_by(Region,Site,Plant) %>% summarize(num1=sum(num,na.rm=T)) %>% group_by(Region,Site) %>% summarize(n_plants=sum(num1,na.rm=T))

ss2a <- full_join(ss1a,ss1b)
ss2 <- full_join(ss2a,ss1c)

###################################################################################
#### pop means #####
### SMA ANALYSIS
sma_pop1 <- sma(totall~Biomass_g, data=s5p)
sma_pop1
summary(sma_pop1)
plot(sma_pop1)
plot(sma_pop1, which="residual")
plot(sma_pop1, which="qq")

sma_pop2 <- sma(totall~Biomass_g*Region, data=s5p, robust=T)
sma_pop2
summary(sma_pop2)
plot(sma_pop2)
plot(sma_pop2, which="residual")
plot(sma_pop2, which="qq")

smap1 <- s5p
smap1$y1 <- (coef(sma_pop1)[1]+smap1$Biomass_g*coef(sma_pop1)[2]) 
smap1$yMT <- (coef(sma_pop2)[1,1]+smap1$Biomass_g*coef(sma_pop2)[1,2]) 
smap1$yWI <- (coef(sma_pop2)[2,1]+smap1$Biomass_g*coef(sma_pop2)[2,2]) 

p1b <- ggplot(smap1, aes(x=Biomass_g,y=totall)) + 
  geom_smooth(data=s5p, aes(x=Biomass_g,y=totall), method="lm", se=T, color="black", size=1.5)+
  #geom_line(data=smap1, aes(x=Biomass_g,y=y1), color="white", size=2)+
  #geom_ribbon(data=smap1, aes(ymin=y1L,ymax=y1U,x=Biomass_g, fill="band"))+
  #geom_line(data=smap1 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1.25, lty=3)+
  geom_line(data=smap1 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1.25, lty=1)+
  geom_point(data=smap1, aes(x=Biomass_g,y=totall,color=Region), size=3, pch=16) +
  geom_point(data=smap1 %>% filter(Site=='MS'), aes(x=Biomass_g,y=totall), size=3.25, color="goldenrod", pch=15) +
  geom_point(data=smap1 %>% filter(Site=='BEN'), aes(x=Biomass_g,y=totall), size=3.25, color="dodgerblue", pch=15) +
  annotate(geom='text', x=.3,y=3.75, label="A", fontface='bold', size=5)+
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Biomass (mg)", limits = c(0.3,2.55), breaks=seq(.5,2.5,.5))+
  scale_y_continuous("Total terpenes log(mg/g)", limits = c(1.5,3.8), breaks=seq(1.5,3.5,.5))+
  theme(legend.key=element_blank(), legend.position = c(.85,.8), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))


##### big pops 
pop1 <- s51 %>% filter(Site=='BEN'|Site=='MS') 
pop1a <- pop1 %>% group_by(Region,Site,Plant) %>% summarize(Biomass_sq=mean(Biomass_sq,na.rm=T),sla=mean(sla,na.rm=T), trics=mean(trics, na.rm=T), 
                                                            Biomass_g=mean(Biomass_g, na.rm=T), total=mean(total, na.rm=T), totTC=mean(totTC, na.rm=T),
                                                            thymol=mean(thymol, na.rm=T), carvacrol=mean(carvacrol, na.rm=T), totall=mean(totall,na.rm=T))

ggplot(pop1a, aes(x=Biomass_g,y=totall))+geom_point()+geom_smooth(method="lm")
ggplot(pop1a, aes(x=Biomass_g,y=totall))+geom_point()+facet_wrap(~Region)+geom_smooth(method="lm")
ggplot(pop1a, aes(x=trics, y=totall, color=Region))+ geom_point() + geom_smooth(method='lm')

sma1 <- sma(totall~Biomass_g, data=pop1a)
sma1
summary(sma1)
plot(sma1)

sma2 <- sma(totall~Biomass_g*Region, data=pop1a)
sma2
summary(sma2)
plot(sma2)
plot(sma2, which="residual")
plot(sma2, which="qq")


smap2 <- pop1a
smap2$y1 <- (coef(sma1)[1]+smap2$Biomass_g*coef(sma1)[2]) 
smap2$yMT <- (coef(sma2)[1,1]+smap2$Biomass_g*coef(sma2)[1,2]) 
smap2$yWI <- (coef(sma2)[2,1]+smap2$Biomass_g*coef(sma2)[2,2]) 

## BIG POPS FIG
p2b <- ggplot(pop1a, aes(x=Biomass_g,y=totall)) + 
  geom_smooth(data=smap2, aes(x=Biomass_g,y=totall), method="lm", color="black", size=1.5, se=T)+
  #geom_line(data=smap2 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1.25, lty=3)+
  #geom_line(data=smap2 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1.25, lty=3)+
  geom_point(data=pop1a, aes(x=Biomass_g,y=totall,color=Region), size=2, pch=15) +
  annotate(geom='text', x=.3,y=3.75, label="B", fontface='bold', size=5)+
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Biomass (mg)", limits = c(0.3,2.55), breaks=seq(.5,2.5,.5))+
  scale_y_continuous("Total terpenes log(mg/g)", limits = c(1.5,3.8), breaks=seq(1.5,3.5,.5))+
  theme(legend.key=element_blank(), legend.position = c(.85,.8), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

tiff("Figure2_chem.tiff", width=4, height=6, units='in',
     res=600, compression='lzw')
p1b/p2b
dev.off()

##########################################################################################################
####### chemotype long format
s51$ChemoT <- if_else(s51$Chemo=='T',1,0)
ggplot(s51, aes(x=Region, y=ChemoT)) + geom_jitter(height=.05, width=.2, size=2, alpha=.5)

glmchem1 <- glmer(ChemoT ~ Region + (1|Site), data=s51, family='binomial' )
Anova(glmchem1)
emmeans(glmchem1, ~Region, type="response")

chemo1 <- s5a %>% select(Region,Site,Plant,Rep,Chemo,trics,a_thujene,myrcene,a_terpinene,limonene,gama_terpinene,p_cymene,thymoquinone,thymol,carvacrol,total) %>% 
  group_by(Region,Site,Plant,Rep,Chemo) %>% filter(limonene<1) %>% 
  gather(key="chemical",value="value", -c(Region:Chemo)) %>% drop_na(Region,Chemo)
chemo2 <- s5a %>% select(Region,Site,Plant,Rep,Chemo,trics,a_thujene,myrcene,a_terpinene,limonene,gama_terpinene,p_cymene,thymoquinone,thymol,carvacrol,total) %>% 
  filter(limonene<1) %>% drop_na(Region,Chemo)

pairs(chemo2[6:16])

ggplot(chemo1 %>% filter(Chemo!='G',Chemo!='X'), aes(x=Region, y=value, color=Chemo)) + geom_boxplot() + facet_wrap(~chemical, scales='free_y')

chemo3 <- chemo2
chemo3[6:16] <- scale(chemo2[6:16])

chemo4 <- chemo3 %>% group_by(Region,Site,Plant,Rep,Chemo) %>% filter(limonene<1) %>% 
  gather(key="chemical",value="value", -c(Region:Chemo)) %>% drop_na(Region,Chemo) %>% filter(Chemo!='G',Chemo!='X')

ggplot(chemo4 , aes(x=Region, y=value, color=Chemo)) + geom_boxplot() + facet_wrap(~chemical, scales='free_y')

c1 <- lmer(value ~ Region*Chemo*chemical + (1|Site/Chemo/chemical), data=chemo4)
anova(c1)

emmeans(c1, pairwise~Region|chemical:Chemo)
emmeans(c1, pairwise~Chemo|chemical)

#########################################################################################################
### summaries
s5a$obs <- 1
summ1 <- s5a %>% group_by(Site,Plant) %>% summarize(plants=mean(obs, na.rm=T)) %>% group_by(Site) %>% summarize(plan=sum(plants,na.rm=T))
summ1

#############################################################################################################
#### load bioassay data

bioass <- read_csv("Monarda Bioassay.csv") ## change 139-1 Pre_Wt from 1.0 to 10.0
plantassay <- s51 %>% group_by(Region,Site,Plant) %>% summarize(Biomass_g=mean(Biomass_g, na.rm=T),trics=mean(trics, na.rm=T),sla=mean(sla, na.rm=T),
                                                                thymol=mean(thymol, na.rm=T), carvacrol=mean(carvacrol, na.rm=T),
                                                                totall=mean(totall, na.rm=T), totTC=mean(totTC, na.rm=T))

s7 <- full_join(plantassay,bioass,by=c("Plant")) %>% drop_na(Pre_Wt) %>% drop_na(Region)
head(s7)
s7$Wt_change <- (s7$Post_Wt-s7$Pre_Wt)

ggplot(s7 %>% filter(Pre_Wt>12,Fed==1 ), aes(x=Pre_Wt, y=Post_Wt, color=Region)) + geom_point() + geom_smooth(method="lm")
ggplot(s7 %>% filter(Pre_Wt>12,Fed==1 ), aes(x=Region, y=Wt_change)) + geom_boxplot()
ggplot(s7 %>% filter(Pre_Wt>12,Fed==1 ), aes(x=trics, y=Wt_change, color=Region)) + geom_point() + geom_smooth(method="lm")
ggplot(s7 %>% filter(Pre_Wt>12,Fed==1 ), aes(x=totall, y=Wt_change, color=Region)) + geom_point() + geom_smooth(method="lm")

ggplot(s7 %>% filter(Pre_Wt>12), aes(x=Pre_Wt, y=Fed, color=Region)) + geom_jitter(height=.1) + 
  geom_smooth(method="glm",method.args = list(family = 'binomial'), se=T)

ggplot(s7 %>% filter(Pre_Wt>12), aes(x=totall, y=Fed, color=Region)) + geom_jitter(height=.1) + 
  geom_smooth(method="glm",method.args = list(family = 'binomial'), se=T)


nplants1 <- s7 %>% filter(Pre_Wt>12) %>% group_by(Region,Fed) %>% count()
nplants2 <- s7 %>% filter(Pre_Wt>12) %>% summarize(meanfed=mean(Fed,na.rm=T))

## differences in Feeding among regions/trichomes
glm1 <- glmer(Fed ~ totall*Region+(1|Site/Plant), data=s7 %>% filter(Pre_Wt>12), family="binomial")
Anova(glm1)
summary(glm1)
r.squaredGLMM(glm1)
glm.dat <- as.data.frame(emmeans(glm1, ~Region, type="response"))
glm.dat1 <- as.data.frame(emtrends(glm1, var="totall", ~Region))

## differences in weight change among regions/trichomes
lm1 <- lmer(Wt_change ~ totall*Region+(1|Site/Plant), data=s7 %>% filter(Fed==1))
anova(lm1)
summary(lm1)
emmeans(lm1, pairwise~Region)
emtrends(lm1, ~1, var="trics")
r.squaredGLMM(lm1)


ggplot(s7 %>% filter(Pre_Wt>12,Fed==1 ), aes(x=totall, y=Wt_change, color=Region)) + geom_point() + geom_smooth(method="lm")


ba.sma1 <- sma(Wt_change~totall, data=s7 %>% filter(Pre_Wt>12, Fed==1))
ba.sma1
summary(ba.sma1)
plot(ba.sma1)

ba.sma2 <- sma(Wt_change~totall*Region, data=s7 %>% filter(Pre_Wt>12, Fed==1))
ba.sma2
summary(ba.sma2)
plot(ba.sma2)

ba.map1 <- s7 %>% filter(Pre_Wt>12, Fed==1)
ba.map1$y1 <- (coef(ba.sma1)[1]+ba.map1$trics*coef(ba.sma1)[2]) 
ba.map1$y1L <- (219.2743)+(ba.map1$trics*-94.12155) 
ba.map1$y1U <- (298.1673)+(ba.map1$trics*-33.89239) 
ba.map1$yWI <- (coef(sma2)[1,1]+ba.map1$trics*coef(ba.sma2)[1,2]) 
ba.map1$yMT <- (coef(sma2)[2,1]+ba.map1$trics*coef(ba.sma2)[2,2])

fig3 <- ggplot(s7 %>% filter(Pre_Wt>12), aes(x=totall, y=Fed, color=Region)) + geom_jitter(height=.05) + 
  geom_smooth(method="glm",method.args = list(family = 'binomial'), se=T)+
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  #scale_x_discrete("Region")+
  scale_x_continuous("Total terpenes log(mg/g)")+
  scale_y_continuous("Proportion feeding", limits = c(-.1,1.1), breaks=seq(0,1,.2))+
  theme(legend.key=element_blank(), legend.position = "top")+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

fig3a <- ggplot(data=s7 %>% filter(Pre_Wt>12)) + geom_hline(yintercept=c(0,1),lty=2,lwd=1.5,color="grey")+ 
  geom_jitter(data=s7 %>% filter(Pre_Wt>12), aes(x=Region, y=Fed, color=Region), height=.05, width=.2, size=2, alpha=.5) +
  geom_errorbar(data=glm.dat1, aes(x=Region, ymin=asymp.LCL, ymax=asymp.UCL, color=Region), width=.2, lwd=2)+
  geom_point(data=glm.dat1, aes(x=Region, y=prob, color=Region), size=5)+
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  #scale_x_discrete("Region")+
  scale_y_continuous("Proportion feeding", limits = c(-.1,1.1), breaks=seq(0,1,.2))+
  theme(legend.key=element_blank(), legend.position = "top")+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

fig3b <- ggplot(ba.map1, aes(x=trics,y=Wt_change)) + 
  geom_smooth(data=ba.map1, aes(x=trics,y=Wt_change), color="black", method="lm")+
  geom_point(data=ba.map1, aes(x=trics,y=Wt_change,color=Region), size=2) +
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Trichome density", limits = c(65,375), breaks=seq(100,350,50))+
  scale_y_continuous("Larvae weight change (mg)", limits = c(-38,-5), breaks=c(-30,-20,-10))+
  theme(legend.key=element_blank(), legend.position = "none")+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

fig3b2 <- ggplot(ba.map1, aes(x=trics,y=Wt_change)) + 
  geom_line(data=ba.map1, aes(x=trics,y=y1), color="black", size=1)+
  geom_point(data=ba.map1, aes(x=trics,y=Wt_change,color=Region), size=2) +
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Trichomes", limits = c(65,375), breaks=seq(100,350,50))+
  scale_y_continuous("Larvae weight change", limits = c(-38,-5), breaks=c(-30,-20,-10))+
  theme(legend.key=element_blank(), legend.position = "none")+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

tiff("Figure3.tiff", width=4, height=6, units='in',
     res=600, compression='lzw')
fig3a/fig3b
dev.off()

### Pop means - bioassay analyses
ba1a <- s7 %>% filter(Pre_Wt>12) %>% group_by(Region,Site) %>% summarize(Fed_prop=mean(Fed,na.rm=T),Wt_change1=mean(Wt_change,na.rm=T),
                                                                         Biomass_g=mean(Biomass_g,na.rm=T),sla=mean(sla,na.rm=T),trics=mean(trics, na.rm=T))
ba1b <- s7 %>% filter(Pre_Wt>12,Fed==1) %>% group_by(Region,Site) %>% summarize(Wt_change=mean(Wt_change,na.rm=T))
ba1 <- full_join(ba1a,ba1b,by=c("Region","Site"))


ggplot(ba1, aes(x=trics, y=Wt_change)) + geom_point() + geom_smooth(method="lm")
ggplot(ba1, aes(x=trics, y=Fed_prop, color=Region)) + geom_jitter(height=.1) + 
  geom_smooth(method="glm",method.args = list(family = 'binomial'), se=T)

ba.sma1 <- sma(Wt_change~trics, data=ba1)
ba.sma1
summary(ba.sma1)
plot(ba.sma1)

ba.sma2 <- sma(Wt_change~trics*Region, data=ba1)
ba.sma2
summary(ba.sma2)
plot(ba.sma2)

ggplot(ba1, aes(x=trics, y=Wt_change, color=Region)) + geom_point() + geom_smooth(method="lm")

### Mom means - bioassay analyses

ba2a <- s7 %>% filter(Pre_Wt>12,Plant<219) %>% group_by(Region,Site,Plant) %>% summarize(Fed_prop=mean(Fed,na.rm=T),Wt_change1=mean(Wt_change,na.rm=T),
                                                                                         Biomass_g=mean(Biomass_g,na.rm=T),sla=mean(sla,na.rm=T),trics=mean(trics, na.rm=T))
ba2b <- s7 %>% filter(Pre_Wt>12,Plant<219,Fed==1) %>% group_by(Region,Site,Plant) %>% summarize(Wt_change=mean(Wt_change,na.rm=T))
ba2 <- full_join(ba2a,ba2b, by=c("Region","Site","Plant"))

ggplot(ba2, aes(x=trics, y=Wt_change, color=Region)) + geom_point() + geom_smooth(method="lm")
ggplot(ba2, aes(x=trics, y=Fed_prop, color=Region)) + geom_jitter(height=.1) + 
  geom_smooth(method="glm",method.args = list(family = 'binomial'), se=T)


ba2.sma1 <- sma(Wt_change~trics, data=ba2)
ba2.sma1
summary(ba2.sma1)
plot(ba2.sma1)

ba2.sma2 <- sma(Wt_change~trics*Region, data=ba2)
ba2.sma2
summary(ba2.sma2)
plot(ba2.sma2)




######################################################################################
######################################################################################
### ANALYSIS - OLD ####################################################################
#######################################################################################
hist(log(s5$Biomass_g))
hist(log(s5$trics))

sum1 <- s5 %>% group_by(Region,Site) %>% summarize(b=mean(Biomass_g, na.rm=T)) 

ggplot(s5, aes(x=Region, y=Biomass_g)) + geom_boxplot()
ggplot(s5, aes(x=Region, y=trics)) + geom_boxplot()

ggplot(s5, aes(x=Region, y=Biomass_g)) + geom_boxplot(outlier.colour = "grey") + geom_point(data=sum1, aes(x=Region,y=b, color=Site), size=2)



#### trait correlations
s51_full <- s5 %>% select(Region,Site,Plant,Biomass_g,trics,sla,Stems2,Stems_new,Biomass_erros,Labeling_errors) %>% na.omit()  %>% mutate(Biomass_sq=sqrt(Biomass_g),trics_sq=sqrt(trics))
s51 <- s5 %>% select(Region,Site,Plant,Biomass_g,trics,sla,Stems2,Stems_new) %>% na.omit()  %>% mutate(Biomass_sq=sqrt(Biomass_g),trics_sq=sqrt(trics))
s5p <- s51 %>% group_by(Region,Site) %>% summarize(sla=mean(sla, na.rm=T), trics_sq=mean(trics_sq, na.rm=T), trics=mean(trics, na.rm=T), Biomass_g=mean(Biomass_g, na.rm=T))

# write.csv(s51_full, "Monarda_Garden_Traits_2019.csv")

# biomass ~ trics

###################################################################################
#### pop means #####
popm1a <- lm(trics~scale(Biomass_g), data=s5p %>% filter(Region=='MT'))
summary(popm1a)
popm1b <- lm(trics~scale(Biomass_g), data=s5p %>% filter(Region=='WI'))
summary(popm1b)

popm1 <- lm(trics~Biomass_g, data=s5p)
anova(popm1)
summary(popm1)
emtrends(popm2, ~Region, var="Biomass_g")

popm2 <- lm(trics~Biomass_g*Region, data=s5p)
anova(popm2)
summary(popm2)
emtrends(popm2, ~Region, var="Biomass_g")


### SMA ANALYSIS
sma_pop1 <- sma(trics~Biomass_g, data=s5p)
sma_pop1
summary(sma_pop1)
plot(sma_pop1)
plot(sma_pop1, which="residual")
plot(sma_pop1, which="qq")

sma_pop2 <- sma(trics~Biomass_g*Region, data=s5p, robust=T)
sma_pop2
summary(sma_pop2)
plot(sma_pop2)
plot(sma_pop2, which="residual")
plot(sma_pop2, which="qq")

smap1 <- s5p
smap1$y1 <- (coef(sma_pop1)[1]+smap1$Biomass_g*coef(sma_pop1)[2]) 
smap1$y1L <- (219.2743)+(smap1$Biomass_g*-94.12155) 
smap1$y1U <- (298.1673)+(smap1$Biomass_g*-33.89239) 
smap1$yMT <- (coef(sma_pop2)[1,1]+smap1$Biomass_g*coef(sma_pop2)[1,2]) 
smap1$yWI <- (coef(sma_pop2)[2,1]+smap1$Biomass_g*coef(sma_pop2)[2,2]) 


p1a <- ggplot(s5p, aes(x=Biomass_g,y=trics)) + 
  #geom_smooth(data=s5p, aes(x=Biomass_g,y=trics), method="lm", se=F, color="grey50")+
  geom_line(data=smap1, aes(x=Biomass_g,y=y1), color="black", size=1.5)+
  #geom_ribbon(data=smap1, aes(ymin=y1L,ymax=y1U,x=Biomass_g, fill="band"))+
  geom_line(data=smap1 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1.25, lty=4)+
  geom_line(data=smap1 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1.25)+
  geom_point(data=s5p, aes(x=Biomass_g,y=trics,color=Region), size=3, pch=16) +
  geom_point(data=s5p %>% filter(Site=='MS'), aes(x=Biomass_g,y=trics), size=3, color="goldenrod", pch=15) +
  geom_point(data=s5p %>% filter(Site=='BEN'), aes(x=Biomass_g,y=trics), size=3, color="dodgerblue", pch=15) +
  geom_point(data=s5p, aes(x=Biomass_g,y=trics,color=Region), size=3) +
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Biomass (mg)", limits = c(0.3,2.4), breaks=c(.5,1,1.5,2))+
  scale_y_continuous("Trichome density", limits = c(80,350), breaks=seq(100,350,50))+
  theme(legend.key=element_blank(), legend.position = c(.8,.75), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))


##### big pops 
pop1 <- s51 %>% filter(Site=='BEN'|Site=='MS',Plant<219) 
pop1a <- pop1 %>% group_by(Region,Site,Plant) %>% summarize(Biomass_sq=mean(Biomass_sq,na.rm=T),sla=mean(sla,na.rm=T), trics=mean(trics, na.rm=T), Biomass_g=mean(Biomass_g, na.rm=T))

ggplot(pop1a, aes(x=Biomass_g,y=trics))+geom_point()+geom_smooth(method="lm")
ggplot(pop1a, aes(x=Biomass_g,y=trics))+geom_point()+facet_wrap(~Region)+geom_smooth(method="lm")

sma1 <- sma(trics~Biomass_g, data=pop1a)
sma1
summary(sma1)
plot(sma1)

sma2 <- sma(trics~Biomass_g*Region, data=pop1a)
sma2
summary(sma2)
plot(sma2)
plot(sma2, which="residual")
plot(sma2, which="qq")


smap2 <- pop1a
smap2$y1 <- (coef(sma1)[1]+smap2$Biomass_g*coef(sma1)[2]) 
smap2$y1L <- (219.2743)+(smap2$Biomass_g*-94.12155) 
smap2$y1U <- (298.1673)+(smap2$Biomass_g*-33.89239) 
smap2$yMT <- (coef(sma2)[1,1]+smap2$Biomass_g*coef(sma2)[1,2]) 
smap2$yWI <- (coef(sma2)[2,1]+smap2$Biomass_g*coef(sma2)[2,2]) 

## FULL FIG
p2 <- ggplot(pop1a, aes(x=Biomass_g,y=trics)) + 
  geom_line(data=smap1, aes(x=Biomass_g,y=y1), color="black", size=2)+
  geom_line(data=smap1 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=2)+
  geom_line(data=smap1 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=2)+
  geom_point(data=s5p, aes(x=Biomass_g,y=trics,color=Region), size=4) +
  geom_line(data=smap2, aes(x=Biomass_g,y=y1), color="black", size=1)+
  geom_line(data=smap2 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1, alpha=.5)+
  geom_line(data=smap2 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1, alpha=.5)+
  geom_point(data=pop1a, aes(x=Biomass_g,y=trics,color=Region), size=2) +
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Biomass", limits = c(0.3,2.6), breaks=c(.5,1,1.5,2))+
  scale_y_continuous("Trichomes", limits = c(80,350), breaks=seq(100,350,50))+
  theme(legend.key=element_blank(), legend.position = c(.8,.8), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

## BIG POPS FIG
p2a <- ggplot(pop1a, aes(x=Biomass_g,y=trics)) + 
  geom_line(data=smap2, aes(x=Biomass_g,y=y1), color="black", size=1)+
  geom_line(data=smap2 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1, lty=3)+
  geom_line(data=smap2 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1)+
  geom_point(data=pop1a, aes(x=Biomass_g,y=trics,color=Region), size=2, pch=15) +
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Biomass (mg)", limits = c(0.3,2.4), breaks=c(.5,1,1.5,2))+
  scale_y_continuous("Trichome density", limits = c(80,350), breaks=seq(100,350,50))+
  theme(legend.key=element_blank(), legend.position = "none")+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))



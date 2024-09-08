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


# Q1: chem diversity between regions and chemotypes ####
chem <- trait %>% select(a_thujene:caryophyllene) #%>% sqrt()

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


## PERMANOVA ####
ad1 <- adonis2(chemeuclid ~ Region+Chemo+Region:Chemo, data=chem1, by="term")
ad1
#anova(ad1, by="margin")
pairwise.perm.manova(chemeuclid, paste(chem1$Region,chem1$Chemo))


#cmatrix <- log(chem1[7:31]+.001) %>% as.matrix()
cmatrix <- sqrt(chem1[7:31]) %>% as.matrix()

chem.mds <- metaMDS(cmatrix, distance = 'bray')
plot(chem.mds)
scores(chem.mds)

## Fig 1 -- create dataset for plotting ####
### pivot chemicals longer for plotting ####
chem2 <- chem1 %>%  
  pivot_longer(cols=a_thujene:caryophyllene, names_to = 'chemical', values_to = 'conc')

chemtab2 <- chem2 %>% group_by(chemical) %>% summarize(conc=mean(conc, na.rm=T))
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


#### supplemental plot ####
chem2plotsupp <- chem2 %>%  
  mutate(compound=if_else(chemical=="a_thujene"|chemical=="myrcene"|chemical=="octen_3_ol"|chemical=="a_terpinene"|chemical=="gama_terpinene"|
                            chemical=="p_cymene"|chemical=="thymoquinone"|chemical=="carvacrol"|chemical=="thymol", chemical, "other"), 
  compound=as.factor(compound)) %>% 
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
ggsave('fig1supp.tiff',fig1supp, width=15, height=12, units="in", dpi=600, compression = "lzw", path="Outputs")

### distance by space and chemistry -- just to check ####
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
  geom_text(data=pcrot, aes(x=NMDS1*.5+c(0,-0,-0,-.15,0,0,0,0), y=NMDS2*.6+c(.03,-.05,-.09,.03,.05,0,0,.06), label=compound1), fontface="bold", color="BLACK", size=5)+
  #geom_text(data=pcrot, aes(x=NMDS1*.5, y=NMDS2*.6, label=compound), fontface="bold", color="BLACK", size=4)+
  xlim(-1.3,1.2)+
  theme_bw(base_size = 24)
nmdsfig

## PHYTOCHEM DIVERSITY PLOT ####
terpplot <- (fig1+nmdsfig)/ ((terprich+terpshan) + plot_layout(guides='collect')) +
  plot_annotation(tag_levels = 'a') #+ plot_layout(guides='collect') & theme(legend.position='top')
ggsave('terpplot.tiff',terpplot, width=15, height=12, units="in", dpi=600, compression = "lzw", path="Outputs")


### analysis for chem diversity ####
ter1 <- glmmTMB(total ~ Region * Chemo + (1|Site), data=chem1, family=lognormal)
Anova(ter1)
r.squaredGLMM(ter1)
simulateResiduals(ter1, plot=T)
emmeans(ter1, pairwise~Region, type="response")

chem1$Chem0 <- ifelse(chem1$Chemo=="T",1,0)


ter2 <- glmmTMB(qD ~ Region * Chemo + (1|Site), data=chem1 , family=gaussian)
Anova(ter2)
r.squaredGLMM(ter2)
simulateResiduals(ter2, plot=T)
check_overdispersion(ter2)
emmeans(ter2, pairwise~Region:Chemo, type="response")
emmeans(ter2, pairwise~Region, type="response")

ter3 <- glmmTMB(q1D ~ Region * Chemo + (1|Site), data=chem1)
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

# clm1a <- glmmTMB(logit(leaf_dmgp1)~Region+(1|Site), data=chemnew1, family=gaussian)
# clm1b <- glmmTMB(logit(leaf_dmgp1)~total+(1|Site), data=chemnew1, family=gaussian)
# clm1c <- glmmTMB(logit(leaf_dmgp1)~Chemo+(1|Site), data=chemnew1, family=gaussian)
# clm1d <- glmmTMB(logit(leaf_dmgp1)~qD+(1|Site), data=chemnew1, family=gaussian)
# clm1e <- glmmTMB(logit(leaf_dmgp1)~shan+(1|Site), data=chemnew1, family=gaussian)
# clm1f <- glmmTMB(logit(leaf_dmgp1)~Region+total+(1|Site), data=chemnew1, family=gaussian)
# clm1g <- glmmTMB(logit(leaf_dmgp1)~Region+qD+(1|Site), data=chemnew1, family=gaussian)
# clm1h <- glmmTMB(logit(leaf_dmgp1)~Region+shan+(1|Site), data=chemnew1, family=gaussian)
# clm1i <- glmmTMB(logit(leaf_dmgp1)~Region+Chemo+(1|Site), data=chemnew1, family=gaussian)
# clm1j <- glmmTMB(logit(leaf_dmgp1)~Chemo+total+(1|Site), data=chemnew1, family=gaussian)
# clm1k <- glmmTMB(logit(leaf_dmgp1)~Chemo+qD+(1|Site), data=chemnew1, family=gaussian)
# clm1l <- glmmTMB(logit(leaf_dmgp1)~Chemo+shan+(1|Site), data=chemnew1, family=gaussian)
# clm1m <- glmmTMB(logit(leaf_dmgp1)~Region*total+(1|Site), data=chemnew1, family=gaussian)
# clm1n <- glmmTMB(logit(leaf_dmgp1)~Region*qD+(1|Site), data=chemnew1, family=gaussian)
# clm1o <- glmmTMB(logit(leaf_dmgp1)~Region*shan+(1|Site), data=chemnew1, family=gaussian)
# clm1p <- glmmTMB(logit(leaf_dmgp1)~Region*Chemo+(1|Site), data=chemnew1, family=gaussian)
# clm1q <- glmmTMB(logit(leaf_dmgp1)~Chemo*total+(1|Site), data=chemnew1, family=gaussian)
# clm1r <- glmmTMB(logit(leaf_dmgp1)~Chemo*qD+(1|Site), data=chemnew1, family=gaussian)
# clm1s <- glmmTMB(logit(leaf_dmgp1)~Chemo*shan+(1|Site), data=chemnew1, family=gaussian)
# clm1 <- glmmTMB(logit(leaf_dmgp1)~1+(1|Site), data=chemnew1, family=gaussian)


bbmle::ICtab(clm1,clm1a,clm1b,clm1c,clm1d,clm1e,clm1f,clm1g,clm1h,clm1i,clm1j,clm1k,clm1l,
             clm1m,clm1n,clm1o,clm1p,clm1q,clm1r,clm1s, weights=T, delta=T, base=T, type="BIC")

clm_mods <- subset(MuMIn::model.sel(clm1,clm1a,clm1b,clm1c,clm1d,clm1e,clm1f,clm1g,clm1h,clm1i,clm1j,clm1k,clm1l,
                 clm1m,clm1n,clm1o,clm1p,clm1q,clm1r,clm1s,rank="BIC"), delta <= 4)
sw(clm_mods)
summary(clm_mods)

clm_avg <- model.avg(clm_mods)
summary(clm_avg)

## get p-values
clm1b2 <- glmmTMB(leaf_dmgp1~total+(1|Site), data=chemnew1, family=beta_family)
clm1f2 <- glmmTMB(leaf_dmgp1~Region+total+(1|Site), data=chemnew1, family=beta_family)

Anova(clm1a)
Anova(clm1b2)
Anova(clm1f2)


check_collinearity(clm1f)


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


layoutpd1 <- c(patchwork::area(t=1, l=1, b=5, r=4),
               patchwork::area(t=1.1, l=4, b=2, r=4))
pd1plot <- pd1a+pd1+plot_layout(design=layoutpd1)


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

# get p-values
Anova(clm3b)
Anova(clm3a)
Anova(clm3d)
Anova(clm3f)


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

layoutpg3 <- c(patchwork::area(t=1, l=1, b=5, r=4),
               patchwork::area(t=1.1, l=4, b=2, r=4))
pg3plot <- pg3a+pg3x+plot_layout(design=layoutpg3)


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


Anova(clm2a)
r2(clm2h)
summary(clm2b)
Anova(clm2g)
emmeans(clm2g, pairwise~Region, type="response")


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




## SEED HEAD MODELS
### SEED LOSS -- updated data in TEAMS monarda garden folder ###################################
seedhead21 <- read_csv("Data/MonardaGarden_SeedDamage2021 (5).csv") %>% drop_na(Seeds_No) %>% drop_na(Region) %>% mutate(HeadSize=as.numeric(HeadSize))
gard21 <- read_csv("Data/MonardaGarden_Data2021.csv")

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


# Damage PLOT ####
pdmg1 <- (ph1+pd1plot)/(pa1a+pg3plot) + plot_layout(guides = 'collect')
#pdmg1 <- (pd1a+ph1)/(pg3+pa1a) + plot_layout(guides = 'collect')

ggsave('herbfig_new.tiff',pdmg1, width=12, height=8, units="in", dpi=600, compression='lzw', path="Outputs")



# Q3: Costs of defense ####
library(smatr)
library(lavaan)
library(lavaanPlot)
phytochem1 <- full_join(gard21,chem1) %>% drop_na(Region,total)
head(phytochem1)
chemnew2 <- phytochem1 %>%
  select(Region,Site,Plant,Chemo,Biomass_g,hg_Tallest,volD1, volD2,total,qD,q1D) %>% 
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

## look at regular linear models
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
sma1 <- sma(totall~Biomass_g , data=chemnew2 )
summary(sma1)
plot(sma1)
plot(sma1, which="residual")
plot(sma1, which="qq")


sma2 <- sma(totall~Biomass_g*Chemo, data=chemnew2 , type="shift")
sma2
summary(sma2)
plot(sma2)
plot(sma2, which="residual")
plot(sma2, which="qq")



### SMATR for richness ####
smar1 <- sma(qD~Biomass_g, data=chemnew2 )
summary(smar1)
plot(smar1)
plot(smar1, which="residual")
plot(smar1, which="qq")


smar2 <- sma(qD~Biomass_g*Chemo, data=chemnew2 , type="shift")
smar2
summary(smar2)
plot(smar2)
plot(smar2, which="residual")
plot(smar2, which="qq")


### SMATR for shannon ####
smas1 <- sma(q1D~Biomass_g, data=chemnew2 )
summary(smas1)
plot(smas1)
plot(smas1, which="residual")
plot(smas1, which="qq")


smas2 <- sma(q1D~Biomass_g*Chemo, data=chemnew2 , type="shift")
smas2
summary(smas2)
plot(smas2)
plot(smas2, which="residual")
plot(smas2, which="qq")



## extract model coefficients to make figure
smap2 <- pop1a
smap2$y1 <- exp(coef(sma1)[1]+smap2$Biomass_g*coef(sma1)[2]) 
smap2$yMT <- exp(coef(sma2)[1,1]+smap2$Biomass_g*coef(sma2)[1,2]) 
smap2$yWI <- exp(coef(sma2)[2,1]+smap2$Biomass_g*coef(sma2)[2,2]) 

### SMATR PLOT ####
cost1 <- ggplot()+
  geom_smooth(data=chemnew2, aes(x=Biomass_g, y=totall), color="black", method="lm")+
  geom_point(data=chemnew2, aes(x=Biomass_g, y=totall, color=Chemo, shape=Region), size=4)+
  scale_x_continuous("Plant biomass (g)")+
  scale_y_continuous(name='Terpene conc. (log mg/g)')+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_shape(name="Origin")+
  theme_bw(base_size = 20)

cost2 <- ggplot()+
  geom_smooth(data=chemnew2, aes(x=Biomass_g, y=qD), color="black", method="lm", linetype="dashed")+
  geom_point(data=chemnew2, aes(x=Biomass_g, y=qD, color=Chemo), size=4)+
  scale_x_continuous("Plant biomass (g)")+
  scale_y_continuous(name='Terpene richness')+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_shape(name="Origin")+
  theme_bw(base_size = 20)

cost3 <- ggplot(chemnew2, aes(x=Biomass_g, y=q1D, color=Chemo))+
  geom_smooth(method="lm", linetype="dashed")+
  geom_point(size=4)+
  scale_x_continuous("Plant biomass (g)")+
  scale_y_continuous(name='Terpene shannon')+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_color_manual(values=c("#A73030FF", "#4A6990FF"))+
  scale_shape(name="Origin")+
  theme_bw(base_size = 20)+
  guides(color="none")

costplot <- (cost1+cost2+cost3)+
  plot_annotation(tag_levels = 'a') + plot_layout(guides='collect') & theme(legend.position='right')
ggsave('costplot.tiff',costplot, width=12, height=3.75, units="in", dpi=600, compression = "lzw")

mcost2 <- glmmTMB(qD~Biomass_g * Chemo + (1|Site), data=chemnew2)
Anova(mcost2)
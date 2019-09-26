install.packages("ggplot2"); install.packages('vegan'); install.packages('dplyr');
library('dplyr');
library("ggplot2"); library('vegan'); install.packages('tidyr'); library('tidyr')
library('sp'); library("MASS")




# load data ---------------------------------------------------------------
data<-read.csv("Vernal_data_final.csv", stringsAsFactors=FALSE)
specid<-read.csv("SpecID.csv", stringsAsFactors=FALSE)
centroids<-read.csv("Pool Centroids.csv")
specid<-read.csv('SpecID.csv')
dvals<-read.csv("full_div_output.csv")


# Mantel Test -------------------------------------------------------------
rownames(centroids) <- centroids[,1]
data$Site.ID<-as.factor(data$Site.ID)
centroids_grazed<- centroids[rownames(centroids) %in% data$Site.ID[data$Grazing == 'Grazed',]]


centroid_dist<-dist(centroids[2:3])
hmm<-isoMDS(centroid_dist)
ordiplot(hmm)
?ordiplot()
data1<-data %>% 
  dplyr::select(-c(Quadrat.., Zone)) %>% 
  filter(Year=="2015") %>%
  group_by(Site.ID, Grazing) %>% 
  summarize_all(funs(mean)) 

veggie_dist<-vegdist(data1[,3:ncol(data1)], na.rm=T)#Bray

data1$Site.ID

hm<-isoMDS(veggie_dist)
ordiplot(hm, type="t")
?vegdist



# D-value figure ----------------------------------------------------------

#make dval figure for evan
dvals_q2<-dvals %>% 
  filter(q=="2")

dvalplot<-dvals_q2%>%
  group_by(Grazing, Level) %>% 
  summarize( 
    sd=sd(Dval), 
    sem = sd(Dval)/sqrt(length(Dval)),
    Dval=mean(Dval))


dvals_q0<-dvals %>% 
  filter(q=="0")

dvalplot2<-dvals_q0%>%
  group_by(Grazing, Level) %>% 
  summarize( 
            sd=sd(Dval), 
            sem = sd(Dval)/sqrt(length(Dval)),
            Dval=mean(Dval))


ggplot(data=dvalplot, aes(x=Level, y=Dval, group=Grazing))+
  scale_x_discrete(limits=c("zone", "pool", "site"), 
                   labels=c("zone"="Alpha (Zone)", "pool"= "Alpha (Pool)", "site" = "Gamma (Pasture)")) +
  geom_line(aes(linetype=Grazing), size=2)+
  geom_pointrange(aes(ymin=Dval-sem, ymax=Dval+sem), size=1)+
  labs(x="",y="Inverse Simpson Index (q=2)")+
  ylim(0,18)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(title=element_text(size=40))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size=40))+
  theme(axis.title.y = element_text(size=40))+
  theme(axis.text.x = element_blank())+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size=40))+ 
  theme(panel.grid.major = element_blank())


ggplot(data=dvalplot2, aes(x=Level, y=Dval, group=Grazing))+
  scale_x_discrete(limits=c("zone", "pool", "site"), 
                   labels=c("zone"="Alpha (Zone)", "pool"= "Alpha (Pool)", "site" = "Gamma (Pasture)")) +
  geom_line(aes(linetype=Grazing), size=2)+
  geom_pointrange(aes(ymin=Dval-sem, ymax=Dval+sem), size=1)+
  labs(x="",y="Species Richness (q=0)")+
  theme(title=element_text(size=40))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size=40))+
  theme(axis.title.y = element_text(size=40))+
  theme(axis.text.x = element_blank())+
  theme(legend.text = element_text(size=40))+
  theme(legend.title = element_text(size=40))+ 
  theme(panel.grid.major = element_blank())



#######################Species Richness, Shannon and Cover#####################

##Calculate relative abundance for pools as a whole

dd<-data %>% 
  mutate(total_cov=rowSums(.[6:66]))

for(i in 6:length(names(dd))) {
    dd[,i] <- (dd[, i] / dd[, 67]*100)
}


species <- seq(from=6, to = 66, by = 1)

df <- data.frame(
                 GrazedEstimate=numeric(),
                 GrazedSE=numeric(),
                 UngrazedEstimate=numeric(),
                 UngrazedSE=numeric(),
                 Pvalue=numeric())

for (s in species){
  

d <- cbind(dd[,1:5],dd[,s])
names(d[6])
d[,6]
names(d)

t <- lm(d[,6] ~ Grazing, data=d)
tt <- t.test(d[,6] ~ Grazing, data=d)
summary(t)


#Grazed estimate
coef(t)[1]

#Ungrazed estimate
coef(t)[2] + coef(t)[1]

#Grazed SE
sqrt(diag(vcov(t)))[1]

#Ungrazed SE
sqrt(diag(vcov(t)))[2]

#Pvalue
tt[3]


df.val <- data.frame(
                 GrazedEstimate=as.numeric(coef(t)[1]),
                 GrazedSE=as.numeric(sqrt(diag(vcov(t)))[1]),
                 UngrazedEstimate=as.numeric(coef(t)[2] + coef(t)[1]),
                 UngrazedSE=as.numeric(sqrt(diag(vcov(t)))[2]),
                 Pvalue=as.numeric(tt[3]))
df <- rbind(df, df.val)

}


options(scipen = 999)
df

SpeciesIDCode <- as.character(names(data)[6:66])

df.output <- cbind(df,SpeciesIDCode)
df.output<-left_join(df.output, specid, by='SpeciesIDCode')

#write.csv(df.output,file="DFoutput_all.csv")




#Calculate relative abundance of habitat zones

dd_zone<-dd %>% 
  #filter(Zone=="Pool")
  #filter(Zone=="Transition")
  filter(Zone=="Upland")

for(i in 6:length(names(dd_zone))) {
    dd_zone[,i] <- (dd_zone[, i] / dd_zone[, 67]*100)
}


species <- seq(from=6, to = 66, by = 1)

df <- data.frame(
                 GrazedEstimate=numeric(),
                 GrazedSE=numeric(),
                 UngrazedEstimate=numeric(),
                 UngrazedSE=numeric(),
                 Pvalue=numeric())

for (s in species){
  
d <- cbind(dd_zone[,1:5],dd_zone[,s])
names(d[6])
d[,6]
names(d)

t <- lm(d[,6] ~ Grazing, data=d)
tt <- t.test(d[,6] ~ Grazing, data=d)
summary(t)


df.val <- data.frame(
                 GrazedEstimate=as.numeric(coef(t)[1]),
                 GrazedSE=as.numeric(sqrt(diag(vcov(t)))[1]),
                 UngrazedEstimate=as.numeric(coef(t)[2] + coef(t)[1]),
                 UngrazedSE=as.numeric(sqrt(diag(vcov(t)))[2]),
                 Pvalue=as.numeric(tt[3]))
df <- rbind(df, df.val)

}


SpeciesIDCode <- as.character(names(data)[6:66])

df.output <- cbind(df,SpeciesIDCode)
df.output<-left_join(df.output, specid, by='SpeciesIDCode')



# Total cover -------------------------------------------------------------

#total cover
total_cover<-data %>%
  dplyr::select(-(Site.ID:Year))%>%
  mutate(total=rowSums(.))


total_cover<-total_cover %>% 
    gather("SpecID", "Total_Cover")%>% 
  mutate(total=total_cover$total) %>% 
  mutate(Rel_abun=((Total_Cover/total)*100)) 

data$Quadrat<-as.factor(data$Quadrat)


#total cover
total_cover<-data %>%
  group_by(Grazing, Year)%>%
  dplyr::select(-(Site.ID:Year), -Quadrat)%>%
  summarise_all(funs(sum)) 


tc_subset<-total_cover %>% 
  ungroup() %>% 
  dplyr::select(-c(Grazing:Year)) %>% 
  mutate(total_cov=rowSums(.))


total_cover$total<-tc_subset$total_cov

ungrazed_total_cover<-total_cover%>%
  filter(Grazing=="Ungrazed") %>%
  ungroup() %>% 
  gather("SpecID", "Total_Cover", -c(total, Grazing, Year))%>% 
  mutate(Rel_abun=((Total_Cover/total)*100)) %>% 
  group_by(SpecID) %>%
  summarise_each(funs(mean, sd, se=sd)) %>% 
  dplyr::select(c(SpecID, Rel_abun_mean, Rel_abun_se))%>% 
  mutate(Grazing=('Grazed'))%>%
  mutate(r = min_rank(desc(Rel_abun_mean)))
#write.csv(ungrazed_total_cover, 'ungrazed_upland.csv')


grazed_total_cover<-total_cover%>%
  filter(Grazing=="Grazed") %>%
  ungroup() %>% 
  gather("SpecID", "Total_Cover", -c(total, Grazing, Year))%>% 
  mutate(Rel_abun=(Total_Cover/total*100)) %>% 
  group_by(SpecID) %>%
  summarise_each(funs(mean, sd, se=sd)) %>% 
  dplyr::select(c(SpecID, Rel_abun_mean, Rel_abun_se)) %>% 
  mutate(Grazing=('Grazed'))%>%
  mutate(r = min_rank(desc(ungrazed_total_cover$Rel_abun_mean)))

#write.csv(grazed_total_cover, 'grazed_upland.csv')

graph<-merge(ungrazed_total_cover, grazed_total_cover, by.x="r", by.y="r")
graph$observation <- 1:nrow(graph) 


ggplot(data = graph)+
  geom_col(data = graph, mapping = aes(x = reorder(observation,-Rel_abun_mean.y),
                                       y = Rel_abun_mean.y), fill="gray55", alpha =.25, color = "grey20")+
  geom_col(mapping = aes(x = observation, 
                         y = (Rel_abun_mean.x)), alpha=.25, fill = "gray3", color = "grey20")+
  labs(x="Species",y="Relative Abundance (% of total cover, 2015-2017)")+
  ggtitle("Evenness by Grazing Treatment")+
  theme(plot.title = element_text(hjust = 0.5, size=40))+
  theme(axis.text.x=element_blank())+
  theme(axis.title.x=element_text(size=30))+
  theme(axis.title.y=element_text(size=30))+
  theme(axis.text.y=element_text(size=40))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ylim(0,25)





#Just evennesss comparison no species

ggplot(graph, aes(observation, (percentcover.x), colour=Grazing.x))+
  geom_line(size=2)+
  geom_line(data=graph, mapping = aes(x=observation, y=(percentcover.y), colour=Grazing.y), size=2)+
  ggtitle("Rank Abundance by Grazing Treatment")+
  labs(x="Species",y="Rank Abundance 2015-2017")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())


ggplot(data = graph)+
  geom_col(data = graph, mapping = aes(x = reorder(observation,-percentcover.x),
                                       y = percentcover.y), fill="red", alpha =.25, color = "grey20")+
  geom_col(mapping = aes(x = reorder(observation,-percentcover.x), 
                           y = (percentcover.x)), alpha=.25, fill = "blue", color = "grey20")+
  labs(x="Species",y="Relative Abundance (% of total cover, 2015-2017)")+
  ggtitle("Evenness by Grazing Treatment")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
  ylim(0,25)
  
  
#compare species
ggplot(data = graph, mapping = aes(x = reorder(Species.x,-percentcover.x), 
                                   y = (percentcover.x+1)))+
  geom_col(alpha=0.35)+
  geom_col(data = graph, mapping = aes(x = reorder(Species.y,-percentcover.x), 
                      y = percentcover.y), fill='lightslateblue', alpha=0.35)+
  labs(x="Species",y="Relative Abundance (% of total cover, 2015-2017)")+
  ggtitle("Evenness by Grazing Treatment")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_blank())
  
ggplot(data = ungrazed_total_cover)+
  geom_col(mapping = aes(x = reorder(Species,-percentcover), 
                         y = percentcover+1, fill=Status))+
  scale_y_log10()+
  ggtitle("Evenness by Grazing Treatment")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_blank())+
  ylim(0,15)

ggplot(data = grazed_total_cover)+
  geom_col(mapping = aes(x = reorder(Species,-percentcover), 
                         y = (percentcover+1), fill=Status))+
  scale_y_log10()+
  ggtitle("Evenness by Grazing Treatment")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_blank())+
  ylim(0,15)




mydata<-data
# species richness
spec_rich<-c()
for(i in 1:nrow(mydata)){
  spec_rich[i]<-sum(mydata[i,-(1:6)]>0)
}

# shannon weiner diversity index
shannon <- c()
for(i in 1:nrow(mydata)){
  temp <- as.numeric(mydata[i,-(1:6)])
  shannon[i] <- diversity(temp)
}

#relative cover natives
natives <- mydata[,colnames(mydata) %in% specid$SpeciesIDCode[specid$Status == 'Native']]
exotics <- mydata[,colnames(mydata) %in% specid$SpeciesIDCode[specid$Status == 'Exotic']]

total_cov_nat<-c()
for(i in 1:nrow(natives)){
  total_cov_nat[i]<-sum(natives[i,1:ncol(natives)])
}

total_cov_ex<-c()
for(i in 1:nrow(exotics)){
  total_cov_ex[i]<-sum(exotics[i,1:ncol(exotics)])
}

total_cov <- c()
for(i in 1:nrow(mydata)){
  temp <- as.numeric(mydata[i,-(1:5)])
  total_cov[i] <- sum(temp)
}

#add to community data set
community_data <- data.frame(mydata,shannon,spec_rich, total_cov_nat, total_cov_ex, total_cov)

# relative cover natives
rel_cov_nat<-c()
for(i in 1:nrow(community_data)){
  rel_cov_nat[i]<-total_cov_nat[i]/total_cov[i]
}

# relative cover exotics
rel_cov_ex<-c()
for(i in 1:nrow(community_data)){
  rel_cov_ex[i]<-total_cov_ex[i]/total_cov[i]
}
#add to community data set
community_data <- data.frame(mydata,shannon,spec_rich, total_cov_nat, total_cov_ex, total_cov, rel_cov_nat, rel_cov_ex)

#Group by zone 
combine_zones<- community_data%>%
  group_by(Year, Grazing, Zone)

#Plot Shannon Weiner
shannon_by_zone <- combine_zones%>%
  summarize(mean_shannon=mean(shannon), 
            sd=sd(shannon), sem = sd(shannon)/sqrt(length(shannon)))

shannon_by_zone$grazingzone <- paste(shannon_by_zone$Grazing, shannon_by_zone$Zone, sep="_")

p.shannon<-ggplot(data=community_data, aes(x=as.factor(Year), y=shannon, fill=Grazing))+
  geom_boxplot()+
  facet_grid(Zone ~ .)+
  labs(x="",y="Shannon Weiner Diversity")+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=16))+
  facet_grid(Zone ~ .)+
  scale_fill_grey()+
  scale_colour_brewer()+
  theme(strip.text.y = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title=element_blank())+
  ggsave("Shannon_Weiner.jpeg", height=7, width=9)
p.shannon

#ANOVA looking at Grazing and Zone as factors
hist(shannon)
anova1.shannon<-aov(shannon~Zone+Grazing+Grazing:Zone+Grazing:as.factor(Year), combine_zones)
summary(anova1.shannon)
TukeyHSD(anova1.shannon)

#Plot relative cover of natives
rel_cov_by_zone<-combine_zones%>%
  summarize(cover=mean(rel_cov_nat), sd=sd(rel_cov_nat), 
            sem = sd(rel_cov_nat)/sqrt(length(rel_cov_nat)))

rel_cov_by_zone$grazingzone <- paste(rel_cov_by_zone$Grazing, rel_cov_by_zone$Zone, sep="_")

p.cover.1<-ggplot(data=rel_cov_by_zone, aes(x=Year, y=cover, group=Grazing))+
  geom_line(aes(linetype=Grazing), size=2)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=cover-sem, ymax=cover+sem), width=.025)+
  labs(x="",y="Relative Cover Natives")+
  theme(axis.text.y = element_text(size=40))+
  theme(axis.title.y = element_text(size=40))+
  theme(axis.text.x = element_text(size=40))+
  facet_grid(Zone ~ .)+
  theme(strip.text.y = element_text(size = 40))+
  theme(legend.text = element_text(size = 40))+
  theme(legend.title=element_blank())+
  scale_x_continuous(breaks = c(2015,2016,2017))+
  ggsave("Relative_Cover_Natives.jpeg", height=7, width=9)
p.cover.1

hist(rel_cov_nat)
anova1.relcov<-aov(rel_cov_nat~Zone+Grazing+Grazing*Zone+Grazing*as.factor(Year), combine_zones)
summary(anova1.relcov)
TukeyHSD(anova1.relcov)

#Plot relative cover of exotics
rel_cov_by_zone<-combine_zones%>%
  summarize(cover=mean(rel_cov_ex), sd=sd(rel_cov_ex), 
            sem = sd(rel_cov_ex)/sqrt(length(rel_cov_ex)))

rel_cov_by_zone$grazingzone <- paste(rel_cov_by_zone$Grazing, rel_cov_by_zone$Zone, sep="_")


p.cover.2<-ggplot(data=rel_cov_by_zone, aes(x=Year, y=cover, group=Grazing))+
  geom_point(size=2)+
  geom_line(aes(linetype=Grazing), size=2)+
  geom_errorbar(aes(ymin=cover-sem, ymax=cover+sem), width=.025)+
  labs(x="",y="Relative Cover Exotics")+
  theme(axis.text.y = element_text(size=40))+
  theme(axis.title.y = element_text(size=40))+
  theme(axis.text.x = element_text(size=40))+
  facet_grid(Zone ~ .)+
  theme(strip.text.y = element_text(size = 40))+
  theme(legend.text = element_text(size = 40))+
  theme(legend.title=element_blank())+
  scale_x_continuous(breaks = c(2015,2016,2017))
p.cover.2
anova1.relcovex<-aov(rel_cov_ex~Zone+Grazing+Grazing*Zone+Grazing*as.factor(Year), combine_zones)
summary(anova1.relcovex)
TukeyHSD(anova1.relcovex)


###Zones combined###

combine_all<- community_data%>%
  group_by(Year, Grazing)

rel_cov_all<-combine_all%>%
  summarize(cover=mean(rel_cov_nat), sd=sd(rel_cov_nat), 
            sem = sd(rel_cov_nat)/sqrt(length(rel_cov_nat)))

shannon_all<-combine_all%>%
  summarize(shannon_weiner=mean(shannon), sd=sd(shannon), 
            sem = sd(shannon)/sqrt(length(shannon)))

spec_rich_all<-combine_all%>%
  summarize(specrich=mean(spec_rich), sd=sd(spec_rich), 
            sem = sd(spec_rich)/sqrt(length(spec_rich)))

#Plot shannon for all years
p.shannon.combined<-ggplot(data=shannon_all, aes(x=Year, y=shannon_weiner, group=Grazing))+
  geom_line(aes(linetype=Grazing), size=.65)+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=shannon_weiner-sem, ymax=shannon_weiner+sem), width=.025)+
  labs(x="",y="Shannon Weiner Diversity")+
  theme(strip.text.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=16))+
  theme(legend.text = element_text(size=16))+
  theme(legend.title = element_text(size=16))+
  theme(legend.title=element_blank())+
  scale_x_continuous(breaks = c(2015,2016,2017))+
  ggsave("Shannon_Weiner_Combined.jpeg", height=7, width=9)
p.shannon.combined

#Plot relative cover of natives for all years
p.cover.combined<-ggplot(data=rel_cov_all, aes(x=Year, y=cover, group=Grazing))+
  geom_line(aes(linetype=Grazing), size=.65)+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=cover-sem, ymax=cover+sem), width=.4)+
  labs(x="",y="Relative Cover Natives")+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=16))+
  theme(legend.text = element_text(size=16))+
  theme(legend.title = element_text(size=16))+
  theme(legend.title=element_blank())+
  scale_x_continuous(breaks = c(2015,2016,2017))+
  ggsave("Relative_Cover_Natives_Combined.jpeg", height=7, width=9)
p.cover.combined

combined.cover.data<- community_data%>%
  group_by(Year, Grazing, Site.ID)%>%
  summarize(rel_cov_nat=mean(rel_cov_nat))

anova1.relcov.all<-aov(rel_cov_nat~Grazing*as.factor(Year), combined.cover.data)
summary(anova1.relcov.all)
TukeyHSD(anova1.relcov.all)



#Plot shannon for all years
p.shannon.combined<-ggplot(data=shannon_all, aes(x=Year, y=shannon_weiner, group=Grazing))+
  geom_line(aes(linetype=Grazing), size=.65)+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=shannon_weiner-sem, ymax=shannon_weiner+sem), width=.025)+
  labs(x="",y="Shannon Weiner Diversity")+
  theme(strip.text.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.y = element_text(size=16))+
  theme(axis.text.x = element_text(size=16))+
  theme(legend.text = element_text(size=16))+
  theme(legend.title = element_text(size=16))+
  theme(legend.title=element_blank())+
  scale_x_continuous(breaks = c(2015,2016,2017))+
  ggsave("Shannon_Weiner_Combined.jpeg", height=7, width=9)
p.shannon.combined
combined.shannon.data<- community_data%>%
  group_by(Year, Grazing, Site.ID)%>%
  summarize(shannon=mean(shannon))
anova1.shannon.all<-aov(shannon~Grazing*as.factor(Year), combined.shannon.data)
summary(anova1.shannon.all)
TukeyHSD(anova1.shannon.all)

install.packages("nlme")
library("nlme")
lme(shannon~Grazing*Year, random=~Grazing|Site.ID+Grazing|Site.ID:Year, combined.shannon.data)




##########################PermANOVA#############################################
mydata<-data
vd<-vegdist(mydata[,6:ncol(mydata)], na.rm=T)#Bray 
#vdgrazed<-vegdist(grazed[,6:ncol(mydata)], na.rm=T, method='jaccard', binary=TRUE)#Jaccard
permanovag<-adonis2(vd~Zone*Grazing*Year, data=mydata)
permanovag





Year2015<-mydata[mydata$Year=="2015",]
Year2016<-mydata[mydata$Year=="2016",]
Year2017<-mydata[mydata$Year=="2017",]
allyears<-mydata%>%
  group_by(Zone, Year)




#Comparing zones in grazed pools
grazed<-allyears[allyears$Grazing=="Grazed",]
#grazed<-Year2015[Year2015$Grazing=="Grazed",]
#grazed<-Year2016[Year2016$Grazing=="Grazed",]
#grazed<-Year2017[Year2017$Grazing=="Grazed",]
vdgrazed<-vegdist(grazed[,6:ncol(mydata)], na.rm=T)#Bray
#vdgrazed<-vegdist(grazed[,6:ncol(mydata)], na.rm=T, method='jaccard', binary=TRUE)#Jaccard
permanovag<-adonis2(vdgrazed~Zone, data=grazed)
permanovag
bdgrazed<-betadisper(vdgrazed, grazed$Zone)
#tiff(filename="C:/users/super/Google Drive/Vegetation Data Analysis/Plot1.tiff", #height=256, width=512)
plot1<-plot(bdgrazed, hull=FALSE, cex=1.5,segments = FALSE,label.cex = .75, main="Grazed Pools by Zone",col = "gray48",ellipse = TRUE)
#dev.off()#make sure it says tiff
boxplot(bdgrazed)
anova(bdgrazed)


#Comparing zones in ungrazed pools
ungrazed<-allyears[allyears$Grazing=="Ungrazed",]
#ungrazed<-Year2015[Year2015$Grazing=="Ungrazed",]
#ungrazed<-Year2016[Year2016$Grazing=="Ungrazed",]
#ungrazed<-Year2017[Year2017$Grazing=="Ungrazed",]
#vdungrazed<-vegdist(ungrazed[,6:ncol(mydata)], na.rm=T)
vdungrazed<-vegdist(ungrazed[,6:ncol(ungrazed)], na.rm=T, method="jaccard", binary=TRUE)
permanovag<-adonis2(vdungrazed~Zone, data=ungrazed)
permanovag
bdungrazed<-betadisper(vdungrazed, ungrazed$Zone)
plot(bdungrazed, hull=FALSE, cex=1.5,segments = FALSE,label.cex = .75, main="Ungrazed Pools by Zone",col = "gray48",ellipse = TRUE)
boxplot(bdungrazed)
anova(bdungrazed)


#Comparing grazed/ungrazed pools
grouped <- read.csv("GroupedYearData_NoExtraPairs.csv")
grouped<-grouped%>%
  select(-BarGr)

#PermANOVA by grazing Year 2015
year.2015 <- grouped[grouped$Year=="2015",]
vd.2015<-vegdist(year.2015[,5:ncol(grouped)], na.rm=T, method="bray")#Bray-Curtis
#vd.2015<-vegdist(year.2015[,5:ncol(grouped)], na.rm=T,method="jaccard", binary=TRUE)#Jaccard
#PermANOVA
permcombined.2015<-adonis2(vd.2015~Grazing, data=year.2015)
permcombined.2015
#PermDISP
bd.2015<-betadisper(vd.2015, year.2015$Grazing)
plot(bd.2015)
boxplot(bd.2015)
anova(bd.2015)

#simper
year.2015.env<-year.2015[,c(1:4)]
year.2015.only<-year.2015[,-c(1:4)]
(sim <- with(year.2015.env, simper(year.2015.only, Grazing)))
summary(sim)


## Principal coordinates analysis with 19 axes to estimate total variance
year.2015.only<-year.2015[,-(1:4)]
Ordination.model1 <- cmdscale(vd.2015, k=19, eig=TRUE, add=FALSE)
Ordination.model1$points[,2] <- -1.0 * Ordination.model1$points[,2]#flip axes
eigenvals(Ordination.model1)
Ordination.model1 <-add.spec.scores(Ordination.model1,year.2015.only,method='pcoa.scores',multi=.005)
sites<-Ordination.model1$points
sites<-as.data.frame(sites)
sites<-data.frame(Site=row.names(sites), sites)
sites["Grazing"]<-NA
sites$Grazing<-year.2015$Grazing

proj<-Ordination.model1$cproj
proj<-as.data.frame(proj)
proj<-data.frame(SpecID=row.names(proj), proj)

#Plot PCoA with Species Loadings (scaled)
ggscatter(sites, x="V1",y="V2",shape="Grazing", size=3)+
  stat_chull(aes(fill = sites$Grazing, fill = sites$Grazing), alpha = 0.1, 
             geom = "polygon", show.legend = FALSE)+
  geom_text(aes(proj$Dim1, proj$Dim2, label=proj$SpecID),check_overlap =TRUE,
            color="grey30", alpha=.75, data=proj)+
  #geom_text(aes(sites$V1, sites$V2, label=sites$Site), data=sites)+
  labs(x="PCoA1",y="PCoA2")+
  ggtitle("2015 Bray-Curtis Dissimilarity")+
  #ggtitle("2015 Jaccard Dissimilarity")+
  scale_fill_grey(end = .5)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="right")+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=14))+
  ggsave("PermANOVA_Jaccard_2015.jpeg", height=7, width=9)









#Follow up PermANOVA to identify zones driving differences
#Pool Zones 2015
pool<-Year2015[Year2015$Zone=="Pool",]
vdpool<-vegdist(pool[,6:ncol(Year2015)], na.rm=T, method="bray")#Bray-Curtis
#vdpool<-vegdist(pool[,6:ncol(Year2015)], na.rm=T, method="jaccard", binary=TRUE)
#PermANOVA
permpool<-adonis2(vdpool~Grazing, data=pool)
permpool
#PermDISP
bdpool<-betadisper(vdpool, pool$Grazing)
plot(bdpool, main="Pool Zones")
boxplot(bdpool)
anova(bdpool)

#Trans Zones 2015
trans<-Year2015[Year2015$Zone=="Transition",]
vdtrans<-vegdist(trans[,6:ncol(mydata)], na.rm=T, method="bray")#Bray-Curtis
#vdtrans<-vegdist(trans[,6:ncol(mydata)], na.rm=T, method="jaccard", binary=TRUE)#Jaccard)
#vdtrans<-vegdist(trans[,5:ncol(mydata)], na.rm=T, method="altGower")
#PermANOVA
permtrans<-adonis2(vdtrans~Grazing, data=trans)
permtrans
#PermDISP
bdtrans<-betadisper(vdtrans, trans$Grazing)
plot(bdtrans, segments = FALSE)
anova(bdtrans)

#Upland Zones 2015
upland<-Year2015[Year2015$Zone=="Upland",]
vdupland<-vegdist(upland[,6:ncol(mydata)], na.rm=T, method="bray")
#vdupland<-vegdist(upland[,6:ncol(mydata)], na.rm=T, method="jaccard", binary=TRUE)
#PermANOVA
permupland<-adonis2(vdupland~Grazing, data=upland)
permupland
#PermDISPv
bdupland<-betadisper(vdupland, upland$Grazing)
bdupland
plot(bdupland, main="Upland Zones")
boxplot(bdupland)
anova(bdupland)


#Follow up PermANOVA to identify zones driving differences
#Pool Zones 2016
pool<-Year2016[Year2016$Zone=="Pool",]
vdpool<-vegdist(pool[,6:ncol(Year2016)], na.rm=T, method="bray")#Bray-Curtis
#vdpool<-vegdist(pool[,6:ncol(Year2016)], na.rm=T, method="jaccard", binary=TRUE)
#PermANOVA
permpool<-adonis2(vdpool~Grazing, data=pool)
permpool
#PermDISP
bdpool<-betadisper(vdpool, pool$Grazing)
plot(bdpool, main="Pool Zones")
boxplot(bdpool)
anova(bdpool)

#Trans Zones 2016
trans<-Year2016[Year2016$Zone=="Transition",]
vdtrans<-vegdist(trans[,6:ncol(mydata)], na.rm=T, method="bray")#Bray-Curtis
#vdtrans<-vegdist(trans[,6:ncol(mydata)], na.rm=T, method="jaccard", binary=TRUE)#Jaccard)
#vdtrans<-vegdist(trans[,5:ncol(mydata)], na.rm=T, method="altGower")
#PermANOVA
permtrans<-adonis2(vdtrans~Grazing, data=trans)
permtrans
#PermDISP
bdtrans<-betadisper(vdtrans, trans$Grazing)
plot(bdtrans, segments = FALSE)
anova(bdtrans)

#Upland Zones 2016
upland<-Year2016[Year2016$Zone=="Upland",]
vdupland<-vegdist(upland[,6:ncol(mydata)], na.rm=T, method="bray")
#vdupland<-vegdist(upland[,6:ncol(mydata)], na.rm=T, method="jaccard", binary=TRUE)
#PermANOVA
permupland<-adonis2(vdupland~Grazing, data=upland)
permupland
#PermDISPv
bdupland<-betadisper(vdupland, upland$Grazing)
bdupland
plot(bdupland, main="Upland Zones")
boxplot(bdupland)
anova(bdupland)

#Follow up PermANOVA to identify zones driving differences
#Pool Zones 2017
pool<-Year2017[Year2017$Zone=="Pool",]
vdpool<-vegdist(pool[,6:ncol(Year2017)], na.rm=T, method="bray")#Bray-Curtis
#vdpool<-vegdist(pool[,6:ncol(Year2017)], na.rm=T, method="jaccard", binary=TRUE)
#PermANOVA
permpool<-adonis2(vdpool~Grazing, data=pool)
permpool
#PermDISP
bdpool<-betadisper(vdpool, pool$Grazing)
plot(bdpool, main="Pool Zones")
boxplot(bdpool)
anova(bdpool)

#Trans Zones 2017
trans<-Year2017[Year2017$Zone=="Transition",]
vdtrans<-vegdist(trans[,6:ncol(mydata)], na.rm=T, method="bray")#Bray-Curtis
#vdtrans<-vegdist(trans[,6:ncol(mydata)], na.rm=T, method="jaccard", binary=TRUE)#Jaccard)
#vdtrans<-vegdist(trans[,5:ncol(mydata)], na.rm=T, method="altGower")
#PermANOVA
permtrans<-adonis2(vdtrans~Grazing, data=trans)
permtrans
#PermDISP
bdtrans<-betadisper(vdtrans, trans$Grazing)
plot(bdtrans, segments = FALSE)
anova(bdtrans)

#Upland Zones 2017
upland<-Year2017[Year2017$Zone=="Upland",]
vdupland<-vegdist(upland[,6:ncol(mydata)], na.rm=T, method="bray")
#vdupland<-vegdist(upland[,6:ncol(mydata)], na.rm=T, method="jaccard", binary=TRUE)
#PermANOVA
permupland<-adonis2(vdupland~Grazing, data=upland)
permupland
#PermDISPv
bdupland<-betadisper(vdupland, upland$Grazing)
bdupland
plot(bdupland, main="Upland Zones")
boxplot(bdupland)
anova(bdupland)
































#Permanova by Grazing Year 2016

year.2016 <- grouped[grouped$Year=="2016",]
#vd.2016 <- vegdist(year.2016[,5:ncol(grouped)], na.rm=T, method="bray")#Bray-Curtis
vd.2016<-vegdist(year.2016[,6:ncol(grouped)], na.rm=T,method="jaccard", binary=TRUE)#Jaccard
#PermANOVA
permcombined.2016<-adonis2(vd.2016~Grazing, data=year.2016)
permcombined.2016
#PermDISP
bd.2016 <- betadisper(vd.2016, year.2016$Grazing)
plot(bd.2016,main = "Vegetation Community by Grazing, Spring 2016")
boxplot(bd.2016)
anova(bd.2016)

##Simper
year.2016.env<-year.2016[,c(1:4)]
year.2016.only<-year.2016[,-c(1:4)]
(sim <- with(year.2016.env, simper(year.2016.only, Grazing)))
summary(sim)

## Principal coordinates analysis with 19 axes to estimate total variance
year.2016.only<-year.2016[,-(1:4)]
Ordination.model1 <- cmdscale(vd.2016, k=19, eig=TRUE, add=FALSE)
Ordination.model1$points[,2] <- -1.0 * Ordination.model1$points[,2]#flip axes
eigenvals(Ordination.model1)
Ordination.model1 <-add.spec.scores(Ordination.model1,year.2016.only,method='pcoa.scores',multi=.007)
sites<-Ordination.model1$points
sites<-as.data.frame(sites)
sites<-data.frame(Site=row.names(sites), sites)
sites["Grazing"]<-NA
sites$Grazing<-year.2016$Grazing

proj<-Ordination.model1$cproj
proj<-as.data.frame(proj)
proj<-data.frame(SpecID=row.names(proj), proj)

#Plot PCoA with Species Loadings (scaled)
ggscatter(sites, x="V1",y="V2",shape="Grazing", size=3)+
  stat_chull(aes(fill = sites$Grazing), alpha = 0.1, 
             geom = "polygon", show.legend = FALSE)+
  geom_text(aes(proj$Dim1, proj$Dim2, label=proj$SpecID), check_overlap = TRUE,color="grey30", alpha=.75, data=proj)+
  #geom_text(aes(sites$V1, sites$V2, label=sites$Site), data=sites)+
  labs(x="PCoA1",y="PCoA2")+
  ggtitle("2016 Bray-Curtis Dissimilarity")+
  #ggtitle("2016 Jaccard Dissimilarity")+
  scale_fill_grey(end = .5)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="right")+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=14))+
  ggsave("PermANOVA_Jaccard_2016.jpeg", height=4, width=6)

#Follow up PermANOVA to identify zones driving differences
#Pool Zones 2016
pool<-Year2016[Year2016$Zone=="Pool",]
#vdpool<-vegdist(pool[,6:ncol(Year2016)], na.rm=T, method="bray")#Bray-Curtis
vdpool<-vegdist(pool[,6:ncol(Year2016)], na.rm=T, method="jaccard", binary=TRUE)
#PermANOVA
permpool<-adonis2(vdpool~Grazing, data=pool)
permpool
#PermDISP
bdpool<-betadisper(vdpool, pool$Grazing)
plot(bdpool, main="Pool Zones")
boxplot(bdpool)
anova(bdpool)
#Trans Zones 2016
trans<-Year2016[Year2016$Zone=="Transition",]
#vdtrans<-vegdist(trans[,6:ncol(mydata)], na.rm=T, method="bray")#Bray-Curtis
vdtrans<-vegdist(trans[,5:ncol(mydata)], na.rm=T, method="jaccard", binary=TRUE)#Jaccard)
#vdtrans<-vegdist(trans[,5:ncol(mydata)], na.rm=T, method="altGower")
#PermANOVA
permtrans<-adonis2(vdtrans~Grazing, data=trans)
permtrans
#PermDISP
bdtrans<-betadisper(vdtrans, trans$Grazing)
plot(bdtrans, segments = FALSE)
anova(bdtrans)

#Upland Zones 2016
upland<-Year2016[Year2016$Zone=="Upland",]
#vdupland<-vegdist(upland[,6:ncol(mydata)], na.rm=T, method="bray")
vdupland<-vegdist(upland[,6:ncol(mydata)], na.rm=T, method="jaccard", binary=TRUE)
#PermANOVA
permupland<-adonis2(vdupland~Grazing, data=upland)
permupland
#PermDISPv
bdupland<-betadisper(vdupland, upland$Grazing)
bdupland
plot(bdupland, main="Upland Zones")
boxplot(bdupland)
anova(bdupland)




#PermANOVA by grazing year 2017
year.2017 <- grouped[grouped$Year=="2017",]
#vd.2017 <- vegdist(year.2017[,5:ncol(grouped)], na.rm=T, method="bray")#Bray-Curtis
vd.2017<-vegdist(year.2017[,5:ncol(grouped)], na.rm=T,method="jaccard", binary=TRUE)#Jaccard
#PermANOVA
permcombined.2017<-adonis2(vd.2017~Grazing, data=year.2017)
permcombined.2017
#PermDISP
bd.2017 <- betadisper(vd.2017, year.2017$Grazing)
plot(bd.2017)
boxplot(bd.2017)
anova(bd.2017)

##Simper
year.2017.env<-year.2017[,c(1:4)]
year.2017.only<-year.2017[,-c(1:4)]
(sim <- with(year.2017.env, simper(year.2017.only, Grazing)))
summary(sim)



## Principal coordinates analysis with 19 axes to estimate total variance
year.2017.only<-year.2017[,-(1:4)]
Ordination.model1 <- cmdscale(vd.2017, k=19, eig=TRUE, add=FALSE)
Ordination.model1$points[,2] <- -1.0 * Ordination.model1$points[,2]#flip axes
eigenvals(Ordination.model1)
Ordination.model1 <-add.spec.scores(Ordination.model1,year.2017.only,method='pcoa.scores',multi=.005)
sites<-Ordination.model1$points
sites<-as.data.frame(sites)
sites<-data.frame(Site=row.names(sites), sites)
sites["Grazing"]<-NA
sites$Grazing<-year.2017$Grazing

proj<-Ordination.model1$cproj
proj<-as.data.frame(proj)
proj<-data.frame(SpecID=row.names(proj), proj)

#Plot PCoA with Species Loadings (scaled)
ggscatter(sites, x="V1",y="V2",shape="Grazing", size=3)+
  stat_chull(aes( fill = sites$Grazing), alpha = 0.1, 
             geom = "polygon", show.legend = FALSE)+
  geom_text(aes(proj$Dim1, proj$Dim2, label=proj$SpecID),check_overlap = TRUE,
            color="grey30", alpha=.75, data=proj)+
  #geom_text(aes(sites$V1, sites$V2, label=sites$Site), data=sites)+
  labs(x="PCoA1",y="PCoA2")+
  ggtitle("2017 Bray-Curtis Dissimilarity")+
  #ggtitle("2017 Jaccard Dissimilarity")+
  scale_fill_grey(end = .5)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="right")+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=14))+
  ggsave("PermANOVA_Jaccard_2017.jpeg", height=4, width=6)




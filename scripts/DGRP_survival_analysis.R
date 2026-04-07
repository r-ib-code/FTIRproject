#### Determine starvation sensitive and starvation resistant DGRP lines using emmeans
#### Generate t-sne plot using starvation resistance data and ATR-spectra
##Written by: Dr Rita Ibrahim and Dr Adam Dobson##
##School of Molecular Biosciences, University of Glasgow, Glasgow G12 8QQ, UK##

rm(list=ls())

#load packages
library(survival)
library(coxme)
library(car)
library(dplyr)
library(emmeans)
library(ggplot2)
library(emmeans)
library(effectsize)
library(ggh4x)
library(effectsize)
library(knitr)
library(ciTools)
library(here)
library(survminer)
library(rms)
library(see)

#create jmp2r to extract survival data
jmp2r <- function(x){
  for(i in 1:nrow(x)){
    xx <- x[i,]
    if(xx$event[1] > 0){
      xxx <- matrix(rep(as.matrix(xx), xx[1,"event"]), ncol=ncol(xx), byrow=T)
      if(exists("a")){a <- rbind(a, xxx)}else{(a <- xxx)}
    }
  }
  a <- as.data.frame(a)
  for(i in 1:ncol(a)){
    colnames(a)[i] <- colnames(x)[i]
    if(class(a[,i]) != class(x[,i])) {a[,i] <- as.character(a[,i]); class(a[,i]) <- class(x[,i])}
  }	
  a <- a[,!colnames(a) %in% c("count", "deathIsZero_censorIsOne", "event")]
  a
}

#input data 
d_0 <- read.csv(".../DGRP-starvationresistance.csv")
d <- jmp2r(d_0)
d$DGRP <- as.factor(d$DGRP)

m1_line <- coxph(Surv(time, censor) ~ DGRP, d)
cox.zph(m1_line)	#PHA violated

dd <- datadist(d)
options(datadist = 'dd')

#parametric survival model
psm1 <- psm(Surv(time, censor) ~ DGRP, 
            data=d,
            dist="logistic")
psm1_gaussian <- update(psm1, dist="gaussian")
AIC(psm1)
AIC(psm1_gaussian)
psm1_emm <- emmeans(psm1, ~ DGRP)
emmeanss <- as.data.frame(psm1_emm)
nrow(emmeanss)
colnames(emmeanss) <- c("DGRP", "emmean", "SE", "df", "CI05", "CI95")
emmeanss<- emmeanss[order(emmeanss$emmean),]

#plot emmeans across the 108 DGRP lines
plot(y=emmeanss$emmean, x=1:nrow(emmeanss), ,xaxt = "n",ylab="Emmeans",xlab="DGRP lines", las=2, pch=16, col="darkorchid4")
arrows(
  x0=1:nrow(emmeanss),
  x1=1:nrow(emmeanss),
  y0=emmeanss$CI05,
  y1=emmeanss$CI95,
  length=0
)
S<- with(d, Surv(time, censor))
sf2_A <- survfit(S ~ DGRP, data=d)
sf2_A
aggregate(time ~ DGRP, subset(d, censor==1), mean)

#Generate a folder with the emmeans values
emmeans_df <- as.data.frame(emmeanss)
emmeans_df$DGRP <- gsub("-", "", emmeans_df$DGRP)
file_path0 <- "../Emmeans.csv"
write.csv(emmeans_df, file = file_path0, row.names = FALSE)

#Grouping --> 20% 80%
lower_threshold <- quantile(emmeans_df$emmean, 0.20)  # 20th percentile
upper_threshold <- quantile(emmeans_df$emmean, 0.80)
sensitive =c()
resistant=c()
for (i in 1:nrow(emmeans_df)) {
  if (emmeans_df$emmean[i]<lower_threshold ) {
    sensitive <- append(sensitive, emmeans_df[i,]$DGRP)
  }
  else if (emmeans_df$emmean[i]>upper_threshold ) {
    resistant <- append(resistant, emmeans_df[i,]$DGRP)
  }
  
}
length(sensitive)
length(resistant)

#Generate dfs containing a list of starvation "sensitive" and "resistant" DGRPs
sensitive_df <- as.data.frame(sensitive)
#remove the dash from the DGRP name so it's compatible with the python script
sensitive_df$sensitive <- gsub("-", "", sensitive_df$sensitive)
resistant_df <- as.data.frame(resistant)
resistant_df$resistant <- gsub("-", "", resistant_df$resistant)

file_path1 <- "../sensitive_df_20%_emmean.csv"
write.csv(sensitive_df, file = file_path1, row.names = FALSE)
file_path2 <- "../resistant_df_80%_emmean.csv"
write.csv(resistant_df, file = file_path2, row.names = FALSE)

######### t-SNE
#setwd("..")
setwd(".../analysis/")
library(Rtsne)
library(ggplot2)

d <- read.table("DGRPFTIR.dat", sep="\t", header=T)
str(d)
dim(d)
colnames(d)

Emmeans <- read.csv("Emmeans.csv")
head(Emmeans)

#check the HRs
plot(emmean ~ as.factor(DGRP), Emmeans)
head(d$Genot.)
head(Emmeans$DGRP)
dGenots <- unique(d$Genot.)
EmmeansDGRPs <- unique(Emmeans$DGRP)
length(dGenots)
length(EmmeansDGRPs)
table(dGenots %in% EmmeansDGRPs)
table(EmmeansDGRPs %in% dGenots)

d <- merge(Emmeans, d, by.y="Genot.", by.x="DGRP")
dim(d)

spectra <- d[,11:ncol(d)]
dim(spectra)
head(colnames(spectra))
dataInd <- d[,1:10]
set.seed(88)
ts1 <- Rtsne(spectra)

#make a dataframe with tsne and line info
combination1 <- ts1$Y
colnames(combination1)
colnames(combination1) <- c("dim1", "dim2")
colnames(combination1)

d2 <- cbind(dataInd, combination1)

p1 <- ggplot(data=d2, aes(y=dim2, x=dim1, col=emmean))
p2 <- p1 + 
  #	geom_point(alpha=0.5) +
  geom_point() +
  theme_bw()
p3 <- p2 + 
  scale_color_viridis_c(option = "turbo") +
  labs(colour="Starvation resistance:\nEmmean", title="Starvation resistance ~ FTIR spectra")
p3

#re run the analysis with XGboost-selected WNs
xgbWns <- read.csv("DGRP_XGBoost_WNS_list.csv")
head(xgbWns)
head(colnames(spectra))
xgbWns <- xgbWns$wns
xgbWns <- paste("X", xgbWns, sep="")

spectra2 <- spectra[,colnames(spectra) %in% xgbWns]
dim(spectra2)

ts2 <- Rtsne(spectra2)
combination2 <- ts2$Y
colnames(combination2)
colnames(combination2) <- c("dim1", "dim2")
colnames(combination2)

d3 <- cbind(dataInd, combination2)

p4 <- ggplot(data=d3, aes(y=dim2, x=dim1, col=emmean))
p5 <- p4 + 
  geom_point() +
  theme_bw()
p6 <- p5 + 
  scale_color_viridis_c(option = "turbo") +
  labs(colour="Starvation resistance:\nEmmean", title="Starvation resistance ~ ATR-FTIR spectra", x="t-SNE 1 (Dimension 1)", y="t-SNE 2 (Dimension 2)")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
p6






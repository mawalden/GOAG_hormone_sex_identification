##################################
##################################
##### title: R script for T thresholds and simulations for manuscript: 
#####     "Follicle-Stimulating Hormone
#####     Resolves Sex Identification of
#####     Hatchling Mojave Desert tortoise
#####     (Gopherus agassizii)"
##### author: M. A. Walden
##### year: 2021
##### month: September
##### day: 24
##### version: 1
##################################
##################################

#####
## Section 1: Prepare workspace
#####


##### Remove objects from global environment
rm(list = ls())

## Set working directory
setwd("C:/Users/macaw/Dropbox/Gopherus_2017/PROJECTS/GOAG_laparoscopy_2021/analyses")

## Load libraries
library(ggplot2)
library(svglite)

## Read in data
dat1 <- read.csv("r_data/dat1.csv",header=T) #surgery torts
datw1 <- read.csv("r_data/datw1.csv",header=T) #wild torts


################
## Monte Carlo simulations-probability of mis-identification
################

## Declare variables and create reference table for parameters
#####
RNFshape <- 3.43
RNFrate <- 0.3
RNMshape <- 3.79
RNMrate <- 0.011
  
RCFshape <- 1.95
RCFrate <- 0.047
RCMshape <- 2.33
RCMrate <-0.003

FSH <- rep(c("naiive","challenged"),2)
sex <- c("female","female","male","male")
shape_min <- c(0.72,1.33,2.99,1.69)
shape_max <- c(325.51,15.39,6.05,4.42)
rate_min <- c(.31,0.022,0.0085,0.0021)
rate_max <- c(28.63,.66,0.019,.006)

vals <- data.frame(sex=sex,FSH=FSH,
                   shape_min=shape_min,shape_max=shape_max,
                   rate_min=rate_min,rate_max=rate_max)
vals


#####
## Identify various threshold T concentrations
#####

## 50%
RNfpdfx<-curve(dgamma(x, shape=RNFshape, rate=RNFrate), from=0, to=200,n=100000)
RNmpdfx <-curve(dgamma(x, shape=RNMshape, rate=RNMrate), from=0, to=200,n=100000,add=TRUE)
RNpdfx <- merge(RNfpdfx,RNmpdfx,by="x")
names(RNpdfx) <- c("x","female","male")
RNpdfx$comp <- as.numeric(RNpdfx$female > RNpdfx$male)
RNpdfx$cross <- c(NA, diff(RNpdfx$comp))
RNcoords <- RNpdfx[which(RNpdfx$cross != 0), c("x", "female")]
RNcoords
points(RNcoords[2,])

RCfpdfx<-curve(dgamma(x, shape=RCFshape, rate=RCFrate), from=1, to=500,n=100000)
RCmpdfx<-curve(dgamma(x, shape=RCMshape, rate=RCMrate), from=1, to=500,n=100000,add=TRUE)
RCpdfx <- merge(RCfpdfx,RCmpdfx,by="x")
names(RCpdfx) <- c("x","female","male")
RCpdfx$comp <- as.numeric(RCpdfx$female > RCpdfx$male)
RCpdfx$cross <- c(NA, diff(RCpdfx$comp))
RCcoords <- RCpdfx[which(RCpdfx$cross != 0), c("x", "female")]
RCcoords
points(RCcoords[2,])

## Integrate pdf over intervals

## Naiive
xT <- seq(0.1,500,by=0.1)
Ndf <- data.frame(x=xT, female=NA, male=NA, total=NA, probF=NA, probM=NA)

fnNF <- function(x) {
  a <- RNFshape
  r <- RNFrate
  s <- 1/r
  1/(s^(a)*gamma(a))*(x^(a-1))*exp(-(x/s))
}
fnNM <- function(x) {
  a <- RNMshape
  r <- RNMrate
  s <- 1/r
  1/(s^(a)*gamma(a))*(x^(a-1))*exp(-(x/s))
}

i <- 1
for(i in 1:nrow(Ndf)) {
  ilow <- Ndf$x[i]-0.05
  ihigh <- Ndf$x[i]+0.049
  intNF <- integrate(fnNF,
                     lower=ilow, upper=ihigh, stop.on.error=F)
  intNM <- integrate(fnNM,
                     lower=ilow, upper=ihigh, stop.on.error=F)
  Ndf$female[i] <- intNF$value
  Ndf$male[i] <- intNM$value
}
Ndf$total <- Ndf$female + Ndf$male
Ndf$probF <- Ndf$female/Ndf$total
Ndf$probM <- Ndf$male/Ndf$total
Ndf

write.csv(Ndf, "r_output/Ndf.csv", row.names=F)

## Challenged
xTc <- seq(0.1,1000,by=0.1)
Cdf <- data.frame(x=xTc, female=NA, male=NA, total=NA, probF=NA, probM=NA)

fnCF <- function(x) {
  a <- RCFshape
  r <- RCFrate
  s <- 1/r
  1/(s^(a)*gamma(a))*(x^(a-1))*exp(-(x/s))
}
fnCM <- function(x) {
  a <- RCMshape
  r <- RCMrate
  s <- 1/r
  1/(s^(a)*gamma(a))*(x^(a-1))*exp(-(x/s))
}

i <- 1
for(i in 1:nrow(Cdf)) {
  ilow <- Cdf$x[i]-0.05
  ihigh <- Cdf$x[i]+0.049
  intCF <- integrate(fnCF,
                     lower=ilow, upper=ihigh, stop.on.error=F)
  intCM <- integrate(fnCM,
                     lower=ilow, upper=ihigh, stop.on.error=F)
  Cdf$female[i] <- intCF$value
  Cdf$male[i] <- intCM$value
}
Cdf$total <- Cdf$female + Cdf$male
Cdf$probF <- Cdf$female/Cdf$total
Cdf$probM <- Cdf$male/Cdf$total
Cdf

write.csv(Cdf, "r_output/Cdf.csv", row.names=F)




#####
## Simulations for probabilistic approach
######
iter <- 10000
i <- 1
iter2 <- iter*2
i2 <- iter+i

##Naiive
outdatRN <- data.frame("FSH"=NA,"iter"=seq(1,iter*2,by=1),
                      femshape=NA, femrate=NA, malshape=NA, malrate=NA,
                      conc=NA, real=NA, sim=NA, sig=NA)

xN <- seq(1,1000,by=0.1)
xS <- c("ovaries","testes")

for (i in 1:iter) {
  #Naive - capture parameter uncertainty
  outdatRN$FSH[i] <- vals$FSH[1]
  outdatRN$real[i] <- "ovaries"
  outdatRN$femshape[i] <- runif(1,vals$shape_min[1],vals$shape_max[1])
  outdatRN$femrate[i] <- runif(1,vals$rate_min[1],vals$rate_max[1])
  outdatRN$malshape[i] <- runif(1,vals$shape_min[3],vals$shape_max[3])
  outdatRN$malrate[i] <- runif(1,vals$rate_min[3],vals$rate_max[3])
  concT <- rgamma(1, shape=outdatRN$femshape[i], rate=outdatRN$femrate[i])
  outdatRN$conc[i] <- concT
  #calculate probabilities from density curves at sampled concentration
  fnNF <- function(x) dgamma(x, shape=outdatRN$femshape[i], rate = outdatRN$femrate[i])
  fnNM <- function(x) dgamma(x, shape=outdatRN$malshape[i], rate = outdatRN$malrate[i])
  concTlow <- ifelse((concT - 0.05) > 0,concT - 0.05,0.01)
  concThigh <- concT + 0.049
  intNF <- integrate(fnNF,
                     lower=concTlow, upper=concThigh, stop.on.error=F)
  intNM <- integrate(fnNM,
                     lower=concTlow, upper=concThigh, stop.on.error=F)
  prb <- c(intNF$value,intNM$value)
  #Assign sex based on probability F versus M at concentration
  outdatRN$sim[i] <- sample(xS, 1, prob = prb)
  outdatRN$sig[i] <- ifelse(outdatRN$sim[i]==outdatRN$real[i],0,1)
}

for (i in i2:iter2) {
  #Naive - capture parameter uncertainty
  outdatRN$FSH[i] <- vals$FSH[1]
  outdatRN$real[i] <- "testes"
  outdatRN$femshape[i] <- runif(1,vals$shape_min[1],vals$shape_max[1])
  outdatRN$femrate[i] <- runif(1,vals$rate_min[1],vals$rate_max[1])
  outdatRN$malshape[i] <- runif(1,vals$shape_min[3],vals$shape_max[3])
  outdatRN$malrate[i] <- runif(1,vals$rate_min[3],vals$rate_max[3])
  concT <- rgamma(1, shape=outdatRN$malshape[i], rate=outdatRN$malrate[i])
  outdatRN$conc[i] <- concT
  #calculate probabilities from density curves at sampled concentration
  fnNF <- function(x) dgamma(x, shape=outdatRN$femshape[i], rate = outdatRN$femrate[i])
  fnNM <- function(x) dgamma(x, shape=outdatRN$malshape[i], rate = outdatRN$malrate[i])
  concTlow <- concT - 0.05
  concThigh <- concT + 0.049
  intNF <- integrate(fnNF,
                     lower=concTlow, upper=concThigh, stop.on.error=F)
  intNM <- integrate(fnNM,
                     lower=concTlow, upper=concThigh, stop.on.error=F)
  prb <- c(intNF$value,intNM$value)
  #Assign sex based on probability F versus M at concentration
  outdatRN$sim[i] <- sample(xS, 1, prob = prb)
  outdatRN$sig[i] <- ifelse(outdatRN$sim[i]==outdatRN$real[i],0,1)
}
naiveFError <- sum(outdatRN$sig[outdatRN$real=="ovaries"],na.rm=T)/iter
naiveMError <- sum(outdatRN$sig[outdatRN$real=="testes"],na.rm=T)/iter
naiveFError
naiveMError
length(outdatRN$sim[outdatRN$sim=="ovaries"])
length(outdatRN$sim[outdatRN$sim=="testes"])


## Challenged
outdatR <- data.frame("FSH"=NA,"iter"=seq(1,iter*2,by=1),
                      femshape=NA, femrate=NA, malshape=NA, malrate=NA,
                      conc=NA, real=NA, sim=NA, sig=NA)

xN <- seq(1,1000,by=0.1)
xS <- c("ovaries","testes")

for (i in 1:iter) {
  #Challenge - capture parameter uncertainty
  outdatR$FSH[i] <- vals$FSH[2]
  outdatR$femshape[i] <- runif(1,vals$shape_min[2],vals$shape_max[2])
  outdatR$femrate[i] <- runif(1,vals$rate_min[2],vals$rate_max[2])
  outdatR$malshape[i] <- runif(1,vals$shape_min[4],vals$shape_max[4])
  outdatR$malrate[i] <- runif(1,vals$rate_min[4],vals$rate_max[4])
  outdatR$real[i] <- "ovaries"
  concT <- rgamma(1, shape=outdatR$femshape[i], rate=outdatR$femrate[i])
  outdatR$conc[i] <- concT
  #calculate probabilities from density curves at sampled concentration
  fnCF <- function(x) dgamma(x, shape=outdatR$femshape[i], rate = outdatR$femrate[i])
  fnCM <- function(x) dgamma(x, shape=outdatR$malshape[i], rate = outdatR$malrate[i])
  concTlow <- concT - 0.05
  concThigh <- concT + 0.049
  intCF <- integrate(fnCF,
                     lower=concTlow, upper=concThigh, stop.on.error=F)
  intCM <- integrate(fnCM,
                     lower=concTlow, upper=concThigh, stop.on.error=F)
  prb <- c(intCF$value,intCM$value)
  #Assign sex based on probability F versus M at concentration
  outdatR$sim[i] <- sample(xS, 1, prob = prb)
  outdatR$sig[i] <- ifelse(outdatR$sim[i]==outdatR$real[i],0,1)
}

for (i in i2:iter2) {
  #Naive - capture parameter uncertainty
  outdatR$FSH[i] <- vals$FSH[2]
  outdatR$femshape[i] <- runif(1,vals$shape_min[2],vals$shape_max[2])
  outdatR$femrate[i] <- runif(1,vals$rate_min[2],vals$rate_max[2])
  outdatR$malshape[i] <- runif(1,vals$shape_min[4],vals$shape_max[4])
  outdatR$malrate[i] <- runif(1,vals$rate_min[4],vals$rate_max[4])
  outdatR$real[i] <- "testes"
  concT <- rgamma(1, shape=outdatR$malshape[i], rate=outdatR$malrate[i])
  outdatR$conc[i] <- concT
  #calculate probabilities from density curves at sampled concentration
  fnCF <- function(x) dgamma(x, shape=outdatR$femshape[i], rate = outdatR$femrate[i])
  fnCM <- function(x) dgamma(x, shape=outdatR$malshape[i], rate = outdatR$malrate[i])
  concTlow <- concT - 0.05
  concThigh <- concT + 0.049
  intCF <- integrate(fnCF,
                     lower=concTlow, upper=concThigh, stop.on.error=F)
  intCM <- integrate(fnCM,
                     lower=concTlow, upper=concThigh, stop.on.error=F)
  prb <- c(intCF$value,intCM$value)
  #Assign sex based on probability F versus M at concentration
  outdatR$sim[i] <- sample(xS, 1, prob = prb)
  outdatR$sig[i] <- ifelse(outdatR$sim[i]==outdatR$real[i],0,1)
}
challengeFError <- sum(outdatR$sig[outdatR$real=="ovaries"],na.rm=T)/iter
challengeMError <- sum(outdatR$sig[outdatR$real=="testes"],na.rm=T)/iter
challengeFError
challengeMError
length(outdatR$sim[outdatR$sim=="ovaries"])
length(outdatR$sim[outdatR$sim=="testes"])






#####
## Simulations for threshold approach
######

iter <- 20000
i <- 1
sext <- c(rep("ovaries",iter/2),rep("testes",iter/2))

##Naiive
femt <- c(20.1, 29.9, 32.6, 35.2)
malt <- c(54.5, 40.9, 37.8, 35.2)

outdatRNt <- data.frame("FSH"=NA,"iter"=seq(1,iter,by=1),
                       femshape=NA, femrate=NA, malshape=NA, malrate=NA,
                       conc=NA, real=sext,
                       sim99=NA, sim95=NA, sim80=NA, sim50=NA,
                       sig99=NA,sig95=NA,sig80=NA,sig50=NA)
for (i in 1:iter) {
  #Naive - capture parameter uncertainty
  outdatRNt$FSH[i] <- vals$FSH[1]
  sexT <- outdatRNt$real[i]
  outdatRNt$femshape[i] <- runif(1,vals$shape_min[1],vals$shape_max[1])
  outdatRNt$femrate[i] <- runif(1,vals$rate_min[1],vals$rate_max[1])
  outdatRNt$malshape[i] <- runif(1,vals$shape_min[3],vals$shape_max[3])
  outdatRNt$malrate[i] <- runif(1,vals$rate_min[3],vals$rate_max[3])
  ml <- if(sexT=="ovaries"){outdatRNt$femshape[i]} else {outdatRNt$malshape[i]}
  sdl <- if(sexT=="ovaries"){outdatRNt$femrate[i]} else {outdatRNt$malrate[i]}
  concT <- rlnorm(1, shape=ml, rate=sdl)
  outdatRNt$conc[i] <- concT
  #Assign sex based on thresholds
  outdatRNt$sim99[i] <- if(concT <= femt[1]){
      "ovaries"} else if(concT >= malt[1]){
          "testes"} else "unknown"
  outdatRNt$sim95[i] <- if(concT <= femt[2]){
    "ovaries"} else if(concT >= malt[2]){
      "testes"} else "unknown"
  outdatRNt$sim80[i] <- if(concT <= femt[3]){
    "ovaries"} else if(concT >= malt[3]){
      "testes"} else "unknown"
  outdatRNt$sim50[i] <- if(concT <= femt[4]){
    "ovaries"} else "testes"
 # Code errors 
  outdatRNt$sig99[i] <- ifelse(outdatRNt$sim99[i]==outdatRNt$real[i] || outdatRNt$sim99[i]=="unknown",0,1)
  outdatRNt$sig95[i] <- ifelse(outdatRNt$sim95[i]==outdatRNt$real[i] || outdatRNt$sim95[i]=="unknown",0,1)
  outdatRNt$sig80[i] <- ifelse(outdatRNt$sim80[i]==outdatRNt$real[i] || outdatRNt$sim80[i]=="unknown",0,1)
  outdatRNt$sig50[i] <- ifelse(outdatRNt$sim50[i]==outdatRNt$real[i] || outdatRNt$sim50[i]=="unknown",0,1)
  
}

naivetFError99 <- sum(outdatRNt$sig99[outdatRNt$real=="ovaries"],na.rm=T)/iter
naivetMError99 <- sum(outdatRNt$sig99[outdatRNt$real=="testes"],na.rm=T)/iter
naivetFError99
naivetMError99
length(outdatRNt$sim99[outdatRNt$sim99=="ovaries"])
length(outdatRNt$sim99[outdatRNt$sim99=="testes"])
length(outdatRNt$sim99[outdatRNt$sim99=="unknown"])

naivetFError95 <- sum(outdatRNt$sig95[outdatRNt$real=="ovaries"],na.rm=T)/iter
naivetMError95 <- sum(outdatRNt$sig95[outdatRNt$real=="testes"],na.rm=T)/iter
naivetFError95
naivetMError95
length(outdatRNt$sim95[outdatRNt$sim95=="ovaries"])
length(outdatRNt$sim95[outdatRNt$sim95=="testes"])
length(outdatRNt$sim95[outdatRNt$sim95=="unknown"])

naivetFError80 <- sum(outdatRNt$sig80[outdatRNt$real=="ovaries"],na.rm=T)/iter
naivetMError80 <- sum(outdatRNt$sig80[outdatRNt$real=="testes"],na.rm=T)/iter
naivetFError80
naivetMError80
length(outdatRNt$sim80[outdatRNt$sim80=="ovaries"])
length(outdatRNt$sim80[outdatRNt$sim80=="testes"])
length(outdatRNt$sim80[outdatRNt$sim80=="unknown"])

naivetFError50 <- sum(outdatRNt$sig50[outdatRNt$real=="ovaries"],na.rm=T)/iter
naivetMError50 <- sum(outdatRNt$sig50[outdatRNt$real=="testes"],na.rm=T)/iter
naivetFError50
naivetMError50
length(outdatRNt$sim50[outdatRNt$sim50=="ovaries"])
length(outdatRNt$sim50[outdatRNt$sim50=="testes"])
length(outdatRNt$sim50[outdatRNt$sim50=="unknown"])

## Challenged
femct <- c(18.4, 73.5, 99.5, 129.1)
malct <- c(602.5, 218.1, 166.0, 129.1)

outdatRCt <- data.frame("FSH"=NA,"iter"=seq(1,iter,by=1),
                        femshape=NA, femrate=NA, malshape=NA, malrate=NA,
                        conc=NA, real=sext,
                        sim99=NA, sim95=NA, sim80=NA, sim50=NA,
                        sig99=NA,sig95=NA,sig80=NA,sig50=NA)
for (i in 1:iter) {
  #Challenged - capture parameter uncertainty
  outdatRCt$FSH[i] <- vals$FSH[2]
  sexT <- outdatRCt$real[i]
  outdatRCt$femshape[i] <- runif(1,vals$shape_min[2],vals$shape_max[2])
  outdatRCt$femrate[i] <- runif(1,vals$rate_min[2],vals$rate_max[2])
  outdatRCt$malshape[i] <- runif(1,vals$shape_min[4],vals$shape_max[4])
  outdatRCt$malrate[i] <- runif(1,vals$rate_min[4],vals$rate_max[4])
  ml <- if(sexT=="ovaries"){outdatRCt$femshape[i]} else {outdatRCt$malshape[i]}
  sdl <- if(sexT=="ovaries"){outdatRCt$femrate[i]} else {outdatRCt$malrate[i]}
  concT <- rlnorm(1, shape=ml, rate=sdl)
  outdatRCt$conc[i] <- concT
  #Assign sex based on thresholds
  outdatRCt$sim99[i] <- if(concT <= femct[1]){
    "ovaries"} else if(concT >= malct[1]){
      "testes"} else "unknown"
  outdatRCt$sim95[i] <- if(concT <= femct[2]){
    "ovaries"} else if(concT >= malct[2]){
      "testes"} else "unknown"
  outdatRCt$sim80[i] <- if(concT <= femct[3]){
    "ovaries"} else if(concT >= malct[3]){
      "testes"} else "unknown"
  outdatRCt$sim50[i] <- if(concT <= femct[4]){
    "ovaries"} else "testes"
  # Code errors 
  outdatRCt$sig99[i] <- ifelse(outdatRCt$sim99[i]==outdatRCt$real[i] || outdatRCt$sim99[i]=="unknown",0,1)
  outdatRCt$sig95[i] <- ifelse(outdatRCt$sim95[i]==outdatRCt$real[i] || outdatRCt$sim95[i]=="unknown",0,1)
  outdatRCt$sig80[i] <- ifelse(outdatRCt$sim80[i]==outdatRCt$real[i] || outdatRCt$sim80[i]=="unknown",0,1)
  outdatRCt$sig50[i] <- ifelse(outdatRCt$sim50[i]==outdatRCt$real[i] || outdatRCt$sim50[i]=="unknown",0,1)
  
}

naivetFError99 <- sum(outdatRCt$sig99[outdatRCt$real=="ovaries"],na.rm=T)/iter
naivetMError99 <- sum(outdatRCt$sig99[outdatRCt$real=="testes"],na.rm=T)/iter
naivetFError99
naivetMError99
length(outdatRCt$sim99[outdatRCt$sim99=="ovaries"])
length(outdatRCt$sim99[outdatRCt$sim99=="testes"])
length(outdatRCt$sim99[outdatRCt$sim99=="unknown"])

naivetFError95 <- sum(outdatRCt$sig95[outdatRCt$real=="ovaries"],na.rm=T)/iter
naivetMError95 <- sum(outdatRCt$sig95[outdatRCt$real=="testes"],na.rm=T)/iter
naivetFError95
naivetMError95
length(outdatRCt$sim95[outdatRCt$sim95=="ovaries"])
length(outdatRCt$sim95[outdatRCt$sim95=="testes"])
length(outdatRCt$sim95[outdatRCt$sim95=="unknown"])

naivetFError80 <- sum(outdatRCt$sig80[outdatRCt$real=="ovaries"],na.rm=T)/iter
naivetMError80 <- sum(outdatRCt$sig80[outdatRCt$real=="testes"],na.rm=T)/iter
naivetFError80
naivetMError80
length(outdatRCt$sim80[outdatRCt$sim80=="ovaries"])
length(outdatRCt$sim80[outdatRCt$sim80=="testes"])
length(outdatRCt$sim80[outdatRCt$sim80=="unknown"])

naivetFError50 <- sum(outdatRCt$sig50[outdatRCt$real=="ovaries"],na.rm=T)/iter
naivetMError50 <- sum(outdatRCt$sig50[outdatRCt$real=="testes"],na.rm=T)/iter
naivetFError50
naivetMError50
length(outdatRCt$sim50[outdatRCt$sim50=="ovaries"])
length(outdatRCt$sim50[outdatRCt$sim50=="testes"])
length(outdatRCt$sim50[outdatRCt$sim50=="unknown"])


#####
## Simulations for range approach
#####
iter <- 20000
i <- 1
sext <- c(rep("ovaries",iter/2),rep("testes",iter/2))

##Naiive
femt <- c(20.8)
malt <- c(125.4)

outdatRNr <- data.frame("FSH"=NA,"iter"=seq(1,iter,by=1),
                        femshape=NA, femrate=NA, malshape=NA, malrate=NA,
                        conc=NA, real=sext,
                        sim=NA, sig=NA)
for (i in 1:iter) {
  #Naive - capture parameter uncertainty
  outdatRNr$FSH[i] <- vals$FSH[1]
  sexT <- outdatRNr$real[i]
  outdatRNr$femshape[i] <- runif(1,vals$shape_min[1],vals$shape_max[1])
  outdatRNr$femrate[i] <- runif(1,vals$rate_min[1],vals$rate_max[1])
  outdatRNr$malshape[i] <- runif(1,vals$shape_min[3],vals$shape_max[3])
  outdatRNr$malrate[i] <- runif(1,vals$rate_min[3],vals$rate_max[3])
  ml <- if(sexT=="ovaries"){outdatRNr$femshape[i]} else {outdatRNr$malshape[i]}
  sdl <- if(sexT=="ovaries"){outdatRNr$femrate[i]} else {outdatRNr$malrate[i]}
  concT <- rlnorm(1, shape=ml, rate=sdl)
  outdatRNr$conc[i] <- concT
  #Assign sex based on thresholds
  outdatRNr$sim[i] <- if(concT <= femt[1]){
    "ovaries"} else if(concT >= malt[1]){
      "testes"} else "unknown"
  # Code errors 
  outdatRNr$sig[i] <- ifelse(outdatRNr$sim[i]==outdatRNr$real[i] || outdatRNr$sim[i]=="unknown",0,1)
}

naivetFError <- sum(outdatRNr$sig[outdatRNr$real=="ovaries"],na.rm=T)/iter
naivetMError <- sum(outdatRNr$sig[outdatRNr$real=="testes"],na.rm=T)/iter
naivetFError
naivetMError
length(outdatRNr$sim[outdatRNr$sim=="ovaries"])
length(outdatRNr$sim[outdatRNr$sim=="testes"])
length(outdatRNr$sim[outdatRNr$sim=="unknown"])

##Challenge
femct <- c(114.7)
malct <- c(139.1)

outdatRCr <- data.frame("FSH"=NA,"iter"=seq(1,iter,by=1),
                        femshape=NA, femrate=NA, malshape=NA, malrate=NA,
                        conc=NA, real=sext,
                        sim=NA, sig=NA)
for (i in 1:iter) {
  #Naive - capture parameter uncertainty
  outdatRCr$FSH[i] <- vals$FSH[2]
  sexT <- outdatRCr$real[i]
  outdatRCr$femshape[i] <- runif(1,vals$shape_min[2],vals$shape_max[2])
  outdatRCr$femrate[i] <- runif(1,vals$rate_min[2],vals$rate_max[2])
  outdatRCr$malshape[i] <- runif(1,vals$shape_min[4],vals$shape_max[4])
  outdatRCr$malrate[i] <- runif(1,vals$rate_min[4],vals$rate_max[4])
  ml <- if(sexT=="ovaries"){outdatRCr$femshape[i]} else {outdatRCr$malshape[i]}
  sdl <- if(sexT=="ovaries"){outdatRCr$femrate[i]} else {outdatRCr$malrate[i]}
  concT <- rlnorm(1, shape=ml, rate=sdl)
  outdatRCr$conc[i] <- concT
  #Assign sex based on thresholds
  outdatRCr$sim[i] <- if(concT <= femct[1]){
    "ovaries"} else if(concT >= malct[1]){
      "testes"} else "unknown"
  # Code errors 
  outdatRCr$sig[i] <- ifelse(outdatRCr$sim[i]==outdatRCr$real[i] || outdatRCr$sim[i]=="unknown",0,1)
}

challengerFError <- sum(outdatRCr$sig[outdatRCr$real=="ovaries"],na.rm=T)/iter
challengerMError <- sum(outdatRCr$sig[outdatRCr$real=="testes"],na.rm=T)/iter
challengerFError
challengerMError
length(outdatRCr$sim[outdatRCr$sim=="ovaries"])
length(outdatRCr$sim[outdatRCr$sim=="testes"])
length(outdatRCr$sim[outdatRCr$sim=="unknown"])

#####
## Test  of proportions
#####

paN <- prop.test(x = 9998, n = 20000, p = 0.5, correct = FALSE)
paC <- prop.test(x = 9990, n = 20000, p = 0.5, correct = FALSE)
taN99 <- prop.test(x = 9788, n = 9788 + 9982, p = 0.5, correct = FALSE)
taN95 <- prop.test(x = 9981, n = 9981 + 9998, p = 0.5, correct = FALSE)
taN80 <- prop.test(x = 9991, n = 9991 + 10000, p = 0.5, correct = FALSE)
taN50 <- prop.test(x = 9993, n = 9993 + 10007, p = 0.5, correct = FALSE)
taC99 <- prop.test(x = 1869, n = 1869 + 7268, p = 0.5, correct = FALSE)
taC95 <- prop.test(x = 8920, n = 8920 + 9185, p = 0.5, correct = FALSE)
taC80 <- prop.test(x = 9587, n = 9587 + 9655, p = 0.5, correct = FALSE)
taC50 <- prop.test(x = 9975, n = 9975 + 10025, p = 0.5, correct = FALSE)
raN <- prop.test(x = 9610, n = 9610 + 9399, p = 0.5, correct = FALSE)
raC <- prop.test(x = 9810, n = 9810 + 9915, p = 0.5, correct = FALSE)

paN
paC
taN99
taN95
taN80
taN50
taC99 #Z 3190.2, <.001
taC95 #Z 3.88, 0.0489
taC80
taC50
raN
raC

8920/(8920 + 9185)
9185/(8920 + 9185)

1869/(1869 + 7268)
7268/(1869 + 7268)
#####
## Plotting
#####
## Declare variables and create reference table for parameters
#####
RNFshape <- 2.46
RNFrate <- 0.29
RNMshape <- 5.72
RNMrate <- 0.6

RCFshape <- 3.45
RCFrate <- 0.68
RCMshape <- 6.5
RCMrate <-0.83

FSH <- rep(c("naiive","challenged"),2)
sex <- c("female","female","male","male")
shape_min <- c(2.25,3.05,5.44,6.10)
shape_max <- c(2.73,3.94,6.06,7.293)
rate_min <- c(0.10,0.25,0.43,0.40)
rate_max <- c(0.44,0.90,0.79,1.34)

vals <- data.frame(sex=sex,FSH=FSH,
                   shape_min=shape_min,shape_max=shape_max,
                   rate_min=rate_min,rate_max=rate_max)
vals

iter <- 50000
i <- 1
iter2 <- iter*2
i2 <- iter+i

##Naiive
outdatRN <- data.frame("FSH"=NA,"iter"=seq(1,iter*2,by=1),
                       femshape=NA, femrate=NA, malshape=NA, malrate=NA,
                       conc=NA, real=NA)

for (i in 1:iter) {
  #Naive - capture parameter uncertainty
  outdatRN$FSH[i] <- vals$FSH[1]
  outdatRN$real[i] <- "ovaries"
  outdatRN$femshape[i] <- runif(1,vals$shape_min[1],vals$shape_max[1])
  outdatRN$femrate[i] <- runif(1,vals$rate_min[1],vals$rate_max[1])
  outdatRN$malshape[i] <- runif(1,vals$shape_min[3],vals$shape_max[3])
  outdatRN$malrate[i] <- runif(1,vals$rate_min[3],vals$rate_max[3])
  concT <- rlnorm(1, shape=outdatRN$femshape[i], rate=outdatRN$femrate[i])
  outdatRN$conc[i] <- concT
}
for (i in i2:iter2) {
  #Naive - capture parameter uncertainty
  outdatRN$FSH[i] <- vals$FSH[1]
  outdatRN$real[i] <- "testes"
  outdatRN$femshape[i] <- runif(1,vals$shape_min[1],vals$shape_max[1])
  outdatRN$femrate[i] <- runif(1,vals$rate_min[1],vals$rate_max[1])
  outdatRN$malshape[i] <- runif(1,vals$shape_min[3],vals$shape_max[3])
  outdatRN$malrate[i] <- runif(1,vals$rate_min[3],vals$rate_max[3])
  concT <- rlnorm(1, shape=outdatRN$malshape[i], rate=outdatRN$malrate[i])
  outdatRN$conc[i] <- concT
}

# ## Challenged
# outdatR <- data.frame("FSH"=NA,"iter"=seq(1,iter*2,by=1),
#                       femshape=NA, femrate=NA, malshape=NA, malrate=NA,
#                       conc=NA, real=NA, sim=NA, sig=NA)
# 
# for (i in 1:iter) {
#   #Challenge - capture parameter uncertainty
#   outdatR$FSH[i] <- vals$FSH[2]
#   outdatR$femshape[i] <- runif(1,vals$shape_min[2],vals$shape_max[2])
#   outdatR$femrate[i] <- runif(1,vals$rate_min[2],vals$rate_max[2])
#   outdatR$malshape[i] <- runif(1,vals$shape_min[4],vals$shape_max[4])
#   outdatR$malrate[i] <- runif(1,vals$rate_min[4],vals$rate_max[4])
#   outdatR$real[i] <- "ovaries"
#   concT <- rlnorm(1, shape=outdatR$femshape[i], rate=outdatR$femrate[i])
#   outdatR$conc[i] <- concT
# }
# for (i in i2:iter2) {
#   #Naive - capture parameter uncertainty
#   outdatR$FSH[i] <- vals$FSH[2]
#   outdatR$femshape[i] <- runif(1,vals$shape_min[2],vals$shape_max[2])
#   outdatR$femrate[i] <- runif(1,vals$rate_min[2],vals$rate_max[2])
#   outdatR$malshape[i] <- runif(1,vals$shape_min[4],vals$shape_max[4])
#   outdatR$malrate[i] <- runif(1,vals$rate_min[4],vals$rate_max[4])
#   outdatR$real[i] <- "testes"
#   concT <- rlnorm(1, shape=outdatR$malshape[i], rate=outdatR$malrate[i])
#   outdatR$conc[i] <- concT
# }

head(outdatRN)
naiivefem <- outdatRN[outdatRN$real=="ovaries",]
naiivemal <- outdatRN[outdatRN$real=="testes",]
# challfem <- outdatR[outdatR$real=="ovaries",]
# challmal <- outdatR[outdatR$real=="testes",]

datwN <- datw1[datw1$challenge=="naiive",]
datwC <- datw1[datw1$challenge=="challenged",]

ggplot(naiivefem, aes(x=conc)) + 
  scale_x_continuous(expand=c(0,0),limits=c(0, 500)) +
  geom_density(fill="white") +
  geom_density(data=naiivemal, aes(x=conc), fill="gray75", lty=2) +
  geom_vline(xintercept = 35.2, lty=3,
             color = "black") +
  geom_rug(data=datwN,aes(x = conc, y = 0, fill=NULL),
           position = position_jitter(height = 0)) +
  xlab("Testosterone (pg/mL") + ylab("Density") + 
  theme(panel.background=element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"),
        legend.position = "none",
        axis.text = element_text(size=14,colour="black"),
        axis.title = element_text(size=14),
        plot.margin = margin(10, 15, 10, 10))

ggsave("r_output/FigXa_densityF.svg", device="svg",
       width=6, height=4, unit="in", dpi=300)

xdf <- data.frame(x=seq(0,500,by=0.0001))
fnCF <- function(x) {
  m <- RCFshape
  s <- RCFrate
  (1/(x*s*sqrt(2*pi)))*exp(-((log(x)-m)^2)/(2*s^2))
}
fnCM <- function(x) {
  m <- RCMshape
  s <- RCMrate
  (1/(x*s*sqrt(2*pi)))*exp(-((log(x)-m)^2)/(2*s^2))
}
ggplot(xdf, aes(x=x)) + 
  scale_x_continuous(expand=c(0,0),limits=c(0, 500)) +
  geom_area(data=xdf, stat="function", fun=fnCM, aes(x=x), fill="gray75", lty=2) +
  stat_function(data=xdf, fun=fnCM, aes(x=x), lty=2) +
  stat_function(data=xdf, fun=fnCF, aes(x=x)) +
#  geom_density(fill="white", alpha=0,lty=1)+
#  geom_density(data=challmal, aes(x=conc), fill="gray75", lty=2) +
#  geom_density(data=challfem, aes(x=conc), fill="white", alpha=0,lty=1) +
  geom_vline(xintercept = 129.1, lty=3,
             color = "black") +
  geom_hline(yintercept = 0, lty=1,
             color = "black") + 
  geom_rug(data=datwC,aes(x = conc, y = 0, fill=NULL),
           position = position_jitter(height = 0)) +
  xlab("Testosterone (pg/mL") + ylab("Density") + 
  theme(panel.background=element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"),
        legend.position = "none",
        axis.text = element_text(size=12,colour="black"),
        axis.title = element_text(size=12),
        plot.margin = margin(10, 15, 10, 10))

ggsave("r_output/FigXb_densityM.svg", device="svg",
       width=6, height=4, unit="in", dpi=300)


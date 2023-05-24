##################################
##################################
##### title: R script for fitting distributions for manuscript: 
#####     "Follicle-Stimulating Hormone
#####     Resolves Sex Identification of
#####     Hatchling Mojave Desert tortoise
#####     (Gopherus agassizii)"
##### author: M. A. Walden
##### year: 2021
##### month: September
##### day: 24
##### version: 1
##### Note: for fitting censored data to distributions:
##### https://cran.r-project.org/web/packages/fitdistrplus/vignettes/FAQ.html
##### https://www.rdocumentation.org/packages/EnvStats/versions/2.4.0/topics/distChooseCensored
##### https://www.rdocumentation.org/packages/EnvStats/versions/2.4.0/topics/Distribution.df
##################################
##################################

################
## Section 1: Prepare workspace
################

##### 
##### Remove objects from global environment
##### 
rm(list = ls())

##### 
## Set working directory
##### 
setwd("C:/Users/macaw/Dropbox/Gopherus_2017/PROJECTS/GOAG_laparoscopy_2021/analyses")

##### 
## Load libraries
##### 
library(fitdistrplus) # fitdistcens(); bootdistcens(); qqcompcens(); ppcompcens()
library(EnvStats) # distChooseCensored(); qqPlotCensored()

##### 
## Read in data
##### 
dat3 <- read.csv("r_data/dat1.csv",header=T) #surgery torts

#####
## Fit distributions for censored data
#####
dev.off()

#Vector indicating censoring
# dat3$conccens <- 0
# dat3$conccens[dat3$rc==1] <- 1

# Fitdistrplus format
dat3$concLeft <- dat3$conc
dat3$concLeft[dat3$rc == 1] <- 500.01
dat3$concLeft[dat3$lc == 1 & dat3$ident == "GS0585"] <- NA ## N = 1/6 naiive females
dat3$concLeft[dat3$lc == 1] <- NA ## N = 4/6 naiive females
dat3$concRight <- dat3$conc
dat3$concRight[dat3$rc == 1] <- NA
dat3$concRight[dat3$lc == 1 & dat3$ident == "GS0585"] <- 3.89 ## N = 1/6 naiive females
dat3$concRight[dat3$lc == 1] <- 3.89 ## N = 4/6 naiive females

#####
## Females
#####

# female naiive #normal
femNcens <- data.frame("left"=dat3$concLeft[dat3$sex=="ovaries" & dat3$challenge == "naiive"],
                        "right"=dat3$concRight[dat3$sex=="ovaries" & dat3$challenge == "naiive"])
femNcensN <- fitdistcens(femNcens, "norm")
femNcensL <- fitdistcens(femNcens, "lnorm")
femNcensG <- fitdistcens(femNcens, "gamma")
femNcensW <- fitdistcens(femNcens, "weibull")
cdfcompcens(list(femNcensN, femNcensL, femNcensG,femNcensW), xlegend = "bottomright",
            fitlwd=2)
qqcompcens(list(femNcensN, femNcensL, femNcensG,femNcensW))
ppcompcens(list(femNcensN, femNcensL, femNcensG,femNcensW))
femNcensN <- fitdist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "naiive"], "norm")
femNcensL <- fitdist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "naiive"], "lnorm")
femNcensG <- fitdist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "naiive"], "gamma")
femNcensW <- fitdist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "naiive"], "weibull")
cdfcomp(list(femNcensN, femNcensL, femNcensG,femNcensW), xlegend = "topleft",
            fitlwd=2)
qqcomp(list(femNcensN, femNcensL, femNcensG,femNcensW))
ppcomp(list(femNcensN, femNcensL, femNcensG,femNcensW))
hist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "naiive"],breaks=4,
     prob=T,xlab="T concentration",main="Emperical and theoretical PDF")
lines(density(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "naiive"],adjust=5),
      col="black", lwd=2)
curve(dnorm(x,mean=femNcensN$estimate[[1]],sd=femNcensN$estimate[[2]]),add=T,
      lty=1, lwd=2, col="red")
curve(dlnorm(x,meanlog=femNcensL$estimate[[1]],sdlog=femNcensL$estimate[[2]]),add=T,
      lty=2, lwd=2, col="green")
curve(dgamma(x,shape=femNcensG$estimate[[1]],rate=femNcensG$estimate[[2]]),add=T,
      lty=3, lwd=2, col="blue")
curve(dweibull(x,shape=femNcensW$estimate[[1]],scale=femNcensW$estimate[[2]]),add=T,
      lty=4, lwd=2, col="turquoise")

ks.test(femNcensN$data, "pnorm", mean=round(femNcensN$estimate[[1]],2),
        sd=round(femNcensN$estimate[[2]],2))
ks.test(femNcensL$data, "plnorm", meanlog=round(femNcensL$estimate[[1]],2),
        sdlog=round(femNcensL$estimate[[2]],2))
ks.test(femNcensG$data, "pgamma", shape=round(femNcensG$estimate[[1]],2),
        rate=round(femNcensG$estimate[[2]],4))
ks.test(femNcensW$data, "pweibull", shape=round(femNcensW$estimate[[1]],2),
        scale=round(femNcensW$estimate[[2]],2))

femNcensL$estimate
summary(femNcensL)

fn <- bootdist(femNcensL,niter=10000)
fn$CI

# female challenged (Gamma best) (none censored)

femCcensN <- fitdist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "challenged"], "norm")
femCcensL <- fitdist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "challenged"], "lnorm")
femCcensG <- fitdist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "challenged"], "gamma")
femCcensW <- fitdist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "challenged"], "weibull")
cdfcomp(list(femCcensN, femCcensL, femCcensG,femCcensW), xlegend = "bottomright",
        fitlwd=2)
qqcomp(list(femCcensN, femCcensL, femCcensG,femCcensW))
ppcomp(list(femCcensN, femCcensL, femCcensG,femCcensW))
hist(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "challenged"],
     prob=T,xlab="T concentration",main="Emperical and theoretical PDF")
lines(density(dat3$conc[dat3$sex=="ovaries" & dat3$challenge == "challenged"],adjust=5),
      col="black", lwd=2)
curve(dnorm(x,mean=femCcensN$estimate[[1]],sd=femCcensN$estimate[[2]]),add=T,
      lty=1, lwd=2, col="red")
curve(dlnorm(x,meanlog=femCcensL$estimate[[1]],sdlog=femCcensL$estimate[[2]]),add=T,
      lty=2, lwd=2, col="green")
curve(dgamma(x,shape=femCcensG$estimate[[1]],rate=femCcensG$estimate[[2]]),add=T,
      lty=3, lwd=2, col="blue")
curve(dweibull(x,shape=femCcensW$estimate[[1]],scale=femCcensW$estimate[[2]]),add=T,
      lty=4, lwd=2, col="turquoise")

ks.test(femCcensN$data, "pnorm", mean=round(femNcensN$estimate[[1]],2),
        sd=round(femNcensN$estimate[[2]],2))

ks.test(femCcensL$data, "plnorm", meanlog=round(femCcensL$estimate[[1]],2),sdlog=round(femCcensL$estimate[[2]],2))
ks.test(femCcensG$data, "pgamma", shape=round(femCcensG$estimate[[1]],2),rate=round(femCcensG$estimate[[2]],4))
ks.test(femCcensW$data, "pweibull", shape=round(femCcensW$estimate[[1]],2),scale=round(femCcensW$estimate[[2]],2))

femCcensL$estimate

fc <- bootdist(femCcensL,niter=10000)
fc$CI

#####
## Males
#####

# Raw male naiive #lnorm best (none censored)

#malNcens <- data.frame("left"=dat3$concLeft[dat3$sex=="testes" & dat3$challenge == "naiive"],
#                        "right"=dat3$concRight[dat3$sex=="testes" & dat3$challenge == "naiive"])
malNcensL <- fitdist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "naiive"], "lnorm")
malNcensG <- fitdist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "naiive"], "gamma")
malNcensW <- fitdist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "naiive"], "weibull")
cdfcomp(list(malNcensL, malNcensG, malNcensW), xlegend = "bottomright",
        fitlwd=2)
qqcomp(list(malNcensL, malNcensG, malNcensW))
ppcomp(list(malNcensL, malNcensG, malNcensW))
hist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "naiive"],
     prob=T,xlab="T concentration",main="Emperical and theoretical PDF")
lines(density(dat3$conc[dat3$sex=="testes" & dat3$challenge == "naiive"],adjust=5),
      col="black", lwd=2)
curve(dlnorm(x,meanlog=malNcensL$estimate[[1]],sdlog=malNcensL$estimate[[2]]),add=T,
      lty=1, lwd=2, col="red")
curve(dgamma(x,shape=malNcensG$estimate[[1]],rate=malNcensG$estimate[[2]]),add=T,
      lty=2, lwd=2, col="green")
curve(dweibull(x,shape=malNcensW$estimate[[1]],scale=malNcensW$estimate[[2]]),add=T,
      lty=3, lwd=2, col="blue")

ks.test(malNcensL$data, "plnorm", meanlog=round(malNcensL$estimate[[1]],2),sdlog=round(malNcensL$estimate[[2]],2))
ks.test(malNcensG$data, "pgamma", shape=round(malNcensG$estimate[[1]],2),rate=round(malNcensG$estimate[[2]],2))
ks.test(malNcensW$data, "pweibull", shape=round(malNcensW$estimate[[1]],2),scale=round(malNcensW$estimate[[2]],2))

malNcensL$estimate

# mn <- distChooseCensored(dat3$concLeft[dat3$sex=="testes" & dat3$challenge=="naiive"],
#                           censored=dat3$conccens[dat3$sex=="testes" & dat3$challenge=="naiive"],
#                           censoring.side="right",
#                           alpha=0.05, method="proucl")
# mn #lognormal
# mn$test.results$lnorm
# mn$distribution.parameters

mnboot <- bootdist(malNcensL,niter=10000)
mnboot$CI

# Raw male challenged ## Gamma
malCcens <- data.frame("left"=dat3$concLeft[dat3$sex=="testes" & dat3$challenge == "challenged"],
                        "right"=dat3$concRight[dat3$sex=="testes" & dat3$challenge == "challenged"])
malCcensL <- fitdistcens(malCcens, "lnorm")
malCcensG <- fitdistcens(malCcens, "gamma")
malCcensW <- fitdistcens(malCcens, "weibull")
cdfcompcens(list(malCcensL, malCcensG,malCcensW), xlegend = "topleft",fitlwd=2)
qqcompcens(list(malCcensL, malCcensG,malCcensW))
ppcompcens(list(malCcensL, malCcensG,malCcensW))
hist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "challenged"],
     prob=T,xlab="T concentration",main="Emperical and theoretical PDF",ylim=c(0,0.001))
lines(density(dat3$conc[dat3$sex=="testes" & dat3$challenge == "challenged"],adjust=5),
      col="black", lwd=2)
curve(dlnorm(x,meanlog=malCcensL$estimate[[1]],sdlog=malCcensL$estimate[[2]]),add=T,
      lty=2, lwd=2, col="red")
curve(dgamma(x,shape=malCcensG$estimate[[1]],rate=malCcensG$estimate[[2]]),add=T,
      lty=3, lwd=2, col="green")
curve(dweibull(x,shape=malCcensW$estimate[[1]],scale=malCcensW$estimate[[2]]),add=T,
      lty=4, lwd=2, col="blue")

malCcensL <- fitdist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "challenged"], "lnorm")
malCcensG <- fitdist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "challenged"], "gamma",
                     start=list(shape=2,scale=1))
malCcensW <- fitdist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "challenged"], "weibull")
cdfcomp(list(malCcensL, malCcensG, malCcensW), xlegend = "bottomright",
        fitlwd=2)
qqcomp(list(malCcensL, malCcensG, malCcensW))
ppcomp(list(malCcensL, malCcensG, malCcensW))
hist(dat3$conc[dat3$sex=="testes" & dat3$challenge == "challenged"],
     prob=T,xlab="T concentration",main="Emperical and theoretical PDF",ylim=c(0,0.001))
lines(density(dat3$conc[dat3$sex=="testes" & dat3$challenge == "challenged"],adjust=5),
      col="black", lwd=2)
curve(dlnorm(x,meanlog=malCcensL$estimate[[1]],sdlog=malCcensL$estimate[[2]]),add=T,
      lty=1, lwd=2, col="red")
curve(dgamma(x,shape=malCcensG$estimate[[1]],scale=malCcensG$estimate[[2]]),add=T,
      lty=2, lwd=2, col="green")
curve(dweibull(x,shape=malCcensW$estimate[[1]],scale=malCcensW$estimate[[2]]),add=T,
      lty=3, lwd=2, col="blue")

ks.test(malCcensL$censdata[2], "plnorm", meanlog=round(malCcensL$estimate[[1]],2),sdlog=round(malCcensL$estimate[[2]],2))
ks.test(malCcensG$censdata[2], "pgamma", shape=round(malCcensG$estimate[[1]],2),scale=round(1/(malCcensG$estimate[[2]]),2))
ks.test(malCcensW$censdata[2], "pweibull", shape=round(malCcensW$estimate[[1]],2),scale=round(malCcensW$estimate[[2]],2))

ks.test(malCcensL$data, "plnorm", meanlog=round(malCcensL$estimate[[1]],2),sdlog=round(malCcensL$estimate[[2]],2))
ks.test(malCcensG$data, "pgamma", shape=round(malCcensG$estimate[[1]],2),scale=malCcensG$estimate[[2]])
ks.test(malCcensW$data, "pweibull", shape=round(malCcensW$estimate[[1]],2),scale=round(malCcensW$estimate[[2]],2))

malCcensG$estimate

malCcensL$estimate

summary(malCcensL)
summary(malCcensG)
summary(malCcensW)

qqcompcens(list(malCcensL, malCcensG, malCcensW))
ppcompcens(list(malCcensL, malCcensG, malCcensW))

# mc <- distChooseCensored(dat3$concLeft[dat3$sex=="testes" & dat3$challenge=="challenged"],
#                         censored=dat3$conccens[dat3$sex=="testes" & dat3$challenge=="challenged"],
#                         censoring.side = "right",
#                         method="proucl")
# mc #normal
# mc$test.results$lnorm
# mc$distribution.parameters
# qqPlotCensored(dat3$concLeft[dat3$sex=="testes" & dat3$challenge=="challenged"],
#                censored=dat3$conccens[dat3$sex=="testes" & dat3$challenge=="challenged"],
#                distribution = "lnorm",points.col = "blue", add.line = TRUE)


mcboot <- bootdist(malCcensL,niter=10000)
mcboot$CI

##################################
##################################
##### title: R script for assay validation for manuscript: 
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

## Remove objects from global environment
rm(list = ls())

## Set working directory
setwd("C:/Users/macaw/Dropbox/Gopherus_2017/PROJECTS/GOAG_laparoscopy_2021/analyses")

##### 
## Load libraries
##### 
library(ggplot2) # ggplot()
library(car) # logit()
library(drc) # drm()
library(lmtest) # lrtest()


#####
## Assay validation
#####

## Prepare data
dil_x <- c(304.32,152.16,76.08,38.04,19.02,9.51,4.755,2.3775,1.18875)
dil_y <- c(20,33,45,58,69,80,87,92,93)
stn_x <- c(500,250,125,62.5,31.25,15.625,7.8125,3.90625)
stn_y <- c(15,25,36,49,61,74,82,89)
dil <- data.frame("x"=dil_x,"y"=dil_y)
stn <- data.frame("x"=stn_x,"y"=stn_y)
dil$y_dec <- dil$y/100
stn$y_dec <- stn$y/100
dil$y_logit <- logit(dil$y_dec,percents=F)
stn$y_logit <- logit(stn$y_dec,percents=F)
dil$x_log <- log(dil$x)
stn$x_log <- log(stn$x)

## Logit vs log: non-normally distributed residuals; F-test inappropriate
moddil <- lm(dil$y_logit~dil$x_log)
modstn <- lm(stn$y_logit~stn$x_log)
var.test(moddil, modstn, alternative="two.sided")
str(moddil)
shapiro.test(moddil$residuals)
shapiro.test(modstn$residuals)

## 4 p log fit
moddil <- drm(dil$y~dil$x, fct=LL.4())
summary(moddil)
plot(moddil)
str(moddil)
modstn <- drm(stn$y~stn$x, fct=LL.4())
summary(modstn)
plot(modstn)

## Liklihood ratio test to compare slope and ED50

# Format datasets into one dataframe
dil$set <- "observed"
stn$set <- "standard"
paral <- rbind(dil,stn)

# Fit 6-parameter model (slope and ED50 different)
# "b"=slope, "c"= lower limit, "d" = upper limit, e"=ED50
par6 <- drm(y~x,set,data=paral,fct=LL.4(), #Both slope and ED50 different
            pmodels=list(~set-1,~1,~1,~set-1))
summary(par6)
par5s <- drm(y~x,set,data=paral,fct=LL.4(), #Only slope different
            pmodels=list(~set-1,~1,~1,~1))
summary(par5s)
par5e <- drm(y~x,set,data=paral,fct=LL.4(), #Only ED50 different
            pmodels=list(~1,~1,~1,~set-1))
summary(par5e)
par4 <- drm(y~x,data=paral,fct=LL.4()) #No difference
summary(par4)

anova(par6,par4)# Both different versus no differences
anova(par5s,par4) #Only slope versus no diff
anova(par5e,par4) #Only ED50 versus no diff

lrtest(par6,par4)
lrtest(par6,par5e,par4)
lrtest(par6,par5s,par4)


## Confidence intervals and plotting
# new dose levels as support for the line
newdata <- expand.grid(conc=exp(seq(log(0.1), log(1000), length=10000)))
newdata2 <- expand.grid(conc=exp(seq(log(0.1), log(1000), length=10000)))

# predictions and confidence intervals
pm <- predict(moddil, newdata=newdata, interval="confidence")
pm2 <- predict(modstn, newdata=newdata2, interval="confidence")

# new data with predictions
newdata$p <- pm[,1]
newdata$pmin <- pm[,2]
newdata$pmax <- pm[,3]

newdata2$p2 <- pm2[,1]
newdata2$pmin2 <- pm2[,2]
newdata2$pmax2 <- pm2[,3]
# plot curve
# need to shift conc == 0 a bit up, otherwise there are problems with coord_trans
dil$conc0 <- dil$x
#dil$conc0[dil$conc0 == 0] <- 0.5

stn$conc0 <- stn$x
#stn$conc0[stn$conc0 == 0] <- 0.5

# plotting the curve
colbreaks <- c("standard"="purple", "Gopherus agassizii"="orange")
dil$type <- "Gopherus agassizii"
stn$type <- "standard"
newdata$type <- "Gopherus agassizii"
newdata2$type <- "standard"



ggplot(dil, aes(x = conc0, y = y)) +
  geom_point(data=stn, aes(x=conc0,y=y,colour=type)) +
  geom_ribbon(data=newdata2, aes(x=conc, y=p2, ymin=pmin2,ymax=pmax2,colour=type,fill=type),
              alpha=0.2) +
  geom_line(data=newdata2, aes(x=conc, y=p2,colour=type)) +
  geom_point(data=dil, aes(x=conc0, y=y, colour=type)) +
  geom_ribbon(data=newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax,colour=type,fill=type),
              alpha=0.2) +
  geom_line(data=newdata, aes(x=conc, y=p,colour=type)) +
  coord_trans(x="log") +
  labs(colour="",fill="") +
  xlab("Plasma testosterone (pg/mL) (log scale)") + ylab("Percent Bound / Maximum Bound (%B/B0)") +
  scale_x_continuous(expand=c(0,0),limits=c(1, 900),
                     breaks=c(1,5,10,20,50,100,200,300,500,800)) +
  scale_y_continuous(expand=c(0,0),limits=c(0, 100),
                     breaks=c(0,5,10,20,50,80,90,95,100)) +
  guides(colour=guide_legend(reverse=TRUE),
         fill=guide_legend(reverse=TRUE)) +
  scale_fill_manual(values=c("orangered2", "blue"), 
                    labels=c(expression(italic("Gopherus agassizii")),"standard")) +
  scale_color_manual(values=c("orangered3", "blue"), 
                     labels=c(expression(italic("Gopherus agassizii")),"standard")) +
  theme(axis.line = element_line(colour="black"),
        axis.text=element_text(colour="black"),
        panel.background=element_blank(),
        panel.grid.major.x = element_line(colour="grey"),
        panel.grid.major.y = element_line(colour="grey"),
        legend.position = c(0.8,0.65),
        legend.title = element_text(size=rel(1)),
        legend.text.align = 0,
        legend.text = element_text(size=rel(1)))

ggsave("r_output/FigX_parallelism.svg", device="svg",
       width=6, height=4, unit="in", dpi=300)


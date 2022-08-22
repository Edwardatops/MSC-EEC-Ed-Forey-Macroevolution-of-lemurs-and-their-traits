# Clear global environment
rm(list=ls())


setwd("C:/Users/Ed/OneDrive/My Documents/Dissertation_2.0/Data")

# Load required packages
library(diversitree)
library(phytools)
library(readxl)
library(ape)
library(caper)
library(dplyr)

# Load in phylogeny
phy <- read.tree("all lemurs.tre")

# Need it to be ultrametric for BAMM and QuaSSE
phy1<- force.ultrametric(phy)

# Load in trait data
d <- read.csv("lemur trait data.csv")
# Rename weird column title
d1 <- rename(d,binomial.name=ï..binomial.name)

# Body mass QuaSSes
# Pick out species names and body mass data
d2 <- d1 %>% dplyr::select(binomial.name, mass_sex_mean)
# Isolate species for which there is data for body mass and is in phy 
d3 <- comparative.data(phy1, d2, binomial.name,na.omit=T)
# Histogram to check for skew
hist(d3$data$mass_sex_mean)
# Very right skewed so log probably better
hist(log(d3$data$mass_sex_mean))
# yep much more normally distributed
# Log body mass data and label each datapoint with species names
mass <- setNames(log(d3$data$mass_sex_mean),rownames(d3$data))
# Phylogeny with species with no body mass data removed
mass.phy <- d3$phy

# Assume that SE of body mass calculation is 1/50, not enough data provided for this
mass.sd <- 1/50

# Starting point for parameters
p <- starting.point.quasse(mass.phy, mass)
p

# Constant speciation model
# Make likelihood function, specifying proportion of phylogeny included
lik.q <- make.quasse(d3$phy, mass,mass.sd,constant.x, constant.x,sampling.f=88/114)

# Constrain drift
lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

# Starting point for parameters
p.start <- c(p[1], p[2:3])
# Label starting points
names(p.start) <- argnames(lik.nodrift)
p.start

# Make lower bounds for parameters
lower <- c(0, 0, 0)

# Run MLE
fit.mass.c <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.mass.c

# Linear
#Make function for linear.x, as not currently in diversitree
xr <- range(mass) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])

# Repeat above steps but fitted to linear.x, then sigmoid.x then hump.x
lik.q <- make.quasse(mass.phy, mass,mass.sd,linear.x, constant.x,sampling.f=88/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.mass.c$par[1],l.m=0, fit.mass.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, -Inf, 0,0)

fit.mass.l <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.mass.l


# Sigmoid
lik.q <- make.quasse(mass.phy,mass,mass.sd,sigmoid.x, constant.x,sampling.f=88/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.mass.c$par[1],fit.mass.c$par[1],mean(xr),1, fit.mass.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(mass), -Inf, 0,0)

fit.mass.s <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.mass.s


# Hump
lik.q <- make.quasse(mass.phy, mass,mass.sd,noroptimal.x, constant.x,sampling.f=88/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.mass.c$par[1],fit.mass.c$par[1],mean(xr),1, fit.mass.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(mass), -Inf, 0,0)

fit.mass.h <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.mass.h

anova(fit.mass.c,mass.linear=fit.mass.l,mass.sigmoid=fit.mass.s,mass.hump=fit.mass.h)




# Repeat above but for
# Gestation period
d4 <- d1 %>% dplyr::select(binomial.name, gest_m)
d5 <- comparative.data(phy1, d4, binomial.name,na.omit=T)

# Check if log vs non log is more normally distributed
hist(d5$data$gest_m)
hist(log(d5$data$gest_m))

gest <- setNames(log(d5$data$gest_m),rownames(d5$data))
gest.phy <- d5$phy

gest.sd <- 1/50

p <- starting.point.quasse(gest.phy, gest)
p

# Constant
lik.q <- make.quasse(d5$phy, gest,gest.sd,constant.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(p[1], p[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, 0)

fit.gest.c <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.gest.c

# Linear
xr <- range(gest) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
lik.q <- make.quasse(gest.phy, gest,gest.sd,linear.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.gest.c$par[1],l.m=0, fit.gest.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, -Inf, 0,0)

fit.gest.l <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.gest.l


# Sigmoid
lik.q <- make.quasse(gest.phy,gest,gest.sd,sigmoid.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.gest.c$par[1],fit.gest.c$par[1],mean(xr),1, fit.gest.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(gest), -Inf, 0,0)

fit.gest.s <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.gest.s


# Hump
lik.q <- make.quasse(gest.phy, gest,gest.sd,noroptimal.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.gest.c$par[1],fit.gest.c$par[1],mean(xr),1, fit.gest.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(gest), -Inf, 0,0)

fit.gest.h <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.gest.h

# Compare the fit of all models, minimal being constant speciation model
anova(fit.gest.c,gestation.linear=fit.gest.l,gestation.sigmoid=fit.gest.s,gestation.hump=fit.gest.h)





# Age at first reproduction
d6 <- d1 %>% dplyr::select(binomial.name, afr_m)
d7 <- comparative.data(phy1, d6, binomial.name,na.omit=T)

#check for skew
hist(d7$data$afr_m)
hist(log(d7$data$afr_m))

afr <- setNames(log(d7$data$afr_m),rownames(d7$data))
afr.phy <- d7$phy

afr.sd <- 1/50

p <- starting.point.quasse(afr.phy, afr)
p

# Constant
lik.q <- make.quasse(d7$phy, afr,afr.sd,constant.x, constant.x,sampling.f=30/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(p[1], p[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, 0)

fit.afr.c <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.afr.c

# Linear
xr <- range(afr) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
lik.q <- make.quasse(afr.phy, afr,afr.sd,linear.x, constant.x,sampling.f=30/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.afr.c$par[1],l.m=0, fit.afr.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, -Inf, 0,0)

fit.afr.l <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.afr.l


# Sigmoid
lik.q <- make.quasse(afr.phy,afr,afr.sd,sigmoid.x, constant.x,sampling.f=30/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.afr.c$par[1],fit.afr.c$par[1],mean(xr),1, fit.afr.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(afr), -Inf, 0,0)

fit.afr.s <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.afr.s


# Hump
lik.q <- make.quasse(afr.phy, afr,afr.sd,noroptimal.x, constant.x,sampling.f=30/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.afr.c$par[1],fit.afr.c$par[1],mean(xr),1, fit.afr.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(afr), -Inf, 0,0)

fit.afr.h <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.afr.h

anova(fit.afr.c,afr.linear=fit.afr.l,afr.sigmoid=fit.afr.s,afr.hump=fit.afr.h)





# Interbirth interval
d8 <- d1 %>% dplyr::select(binomial.name, ibi_m)
d9 <- comparative.data(phy1, d8, binomial.name,na.omit=T)

# Check for skew
hist(d9$data$ibi_m)
hist(log(d9$data$ibi_m))

ibi <- setNames(log(d9$data$ibi_m),rownames(d9$data))
ibi.phy <- d9$phy

ibi.sd <- 1/50

p <- starting.point.quasse(ibi.phy, ibi)
p

# Constant
lik.q <- make.quasse(d9$phy, ibi,ibi.sd,constant.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(p[1], p[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, 0)

fit.ibi.c <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.ibi.c

# Linear
xr <- range(ibi) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
lik.q <- make.quasse(ibi.phy, ibi,ibi.sd,linear.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.ibi.c$par[1],l.m=0, fit.ibi.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, -Inf, 0,0)

fit.ibi.l <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.ibi.l


# Sigmoid
lik.q <- make.quasse(ibi.phy,ibi,ibi.sd,sigmoid.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.ibi.c$par[1],fit.ibi.c$par[1],mean(xr),1, fit.ibi.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(ibi), -Inf, 0,0)

fit.ibi.s <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.ibi.s


# Hump
lik.q <- make.quasse(ibi.phy, ibi,ibi.sd,noroptimal.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.ibi.c$par[1],fit.ibi.c$par[1],mean(xr),1, fit.ibi.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(ibi), -Inf, 0,0)

fit.ibi.h <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.ibi.h

anova(fit.ibi.c,ibi.linear=fit.ibi.l,ibi.sigmoid=fit.ibi.s,ibi.hump=fit.ibi.h)






# Litter size
d10 <- d1 %>% dplyr::select(binomial.name, litter_m)
d11 <- comparative.data(phy1, d10, binomial.name,na.omit=T)

hist(d11$data$litter_m)
hist(log(d11$data$litter_m))

litter <- setNames(log(d11$data$litter_m),rownames(d11$data))
litter.phy <- d11$phy

litter.sd <- 1/50

p <- starting.point.quasse(litter.phy, litter)
p

# Constant
lik.q <- make.quasse(d11$phy, litter,litter.sd,constant.x, constant.x,sampling.f=37/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(p[1], p[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, 0)

fit.litter.c <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.litter.c

# Linear
xr <- range(litter) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
lik.q <- make.quasse(litter.phy, litter,litter.sd,linear.x, constant.x,sampling.f=37/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.litter.c$par[1],l.m=0, fit.litter.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, -Inf, 0,0)

fit.litter.l <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.litter.l


# Sigmoid
lik.q <- make.quasse(litter.phy,litter,litter.sd,sigmoid.x, constant.x,sampling.f=37/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.litter.c$par[1],fit.litter.c$par[1],mean(xr),1, fit.litter.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(litter), -Inf, 0,0)

fit.litter.s <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.litter.s


# Hump
lik.q <- make.quasse(litter.phy, litter,litter.sd,noroptimal.x, constant.x,sampling.f=37/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.litter.c$par[1],fit.litter.c$par[1],mean(xr),1, fit.litter.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(litter), -Inf, 0,0)

fit.litter.h <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.litter.h

anova(fit.litter.c,litter.linear=fit.litter.l,litter.sigmoid=fit.litter.s,litter.hump=fit.litter.h)






# Longevity
d12 <- d1 %>% dplyr::select(binomial.name, longev_m)
d13 <- comparative.data(phy1, d12, binomial.name,na.omit=T)

# Check for skew
hist(d13$data$longev_m)
hist(log(d13$data$longev_m))

longev <- setNames(log(d13$data$longev_m),rownames(d13$data))
longev.phy <- d13$phy

longev.sd <- 1/50

p <- starting.point.quasse(longev.phy, longev)
p

# Constant
lik.q <- make.quasse(d13$phy, longev,longev.sd,constant.x, constant.x,sampling.f=30/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(p[1], p[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, 0)

fit.longev.c <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.longev.c

# Linear
xr <- range(longev) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
lik.q <- make.quasse(longev.phy, longev,longev.sd,linear.x, constant.x,sampling.f=30/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.longev.c$par[1],l.m=0, fit.longev.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, -Inf, 0,0)

fit.longev.l <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.longev.l


# Sigmoid
lik.q <- make.quasse(longev.phy,longev,longev.sd,sigmoid.x, constant.x,sampling.f=30/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.longev.c$par[1],fit.longev.c$par[1],mean(xr),1, fit.longev.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(longev), -Inf, 0,0)

fit.longev.s <- find.mle(lik.nodrift, p.start, control=list(parscale=.1,maxit=200000000),lower=lower, verbose=0)
fit.longev.s


# Hump
lik.q <- make.quasse(longev.phy, longev,longev.sd,noroptimal.x, constant.x,sampling.f=30/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.longev.c$par[1],fit.longev.c$par[1],mean(xr),1, fit.longev.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(longev), -Inf, 0,0)

fit.longev.h <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.longev.h

anova(fit.longev.c,longev.linear=fit.longev.l,longev.sigmoid=fit.longev.s,longev.hump=fit.longev.h)






# Wean length
# Longevity
d14 <- d1 %>% dplyr::select(binomial.name, wean_days_m)
d15 <- comparative.data(phy1, d14, binomial.name,na.omit=T)

wean <- setNames(log(d15$data$wean_days_m),rownames(d15$data))
wean.phy <- d15$phy

wean.sd <- 1/50

p <- starting.point.quasse(wean.phy, wean)
p

# Constant
lik.q <- make.quasse(d15$phy, wean,wean.sd,constant.x, constant.x,sampling.f=31/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(p[1], p[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, 0)

fit.wean.c <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.wean.c

# Linear
xr <- range(wean) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
lik.q <- make.quasse(wean.phy, wean,wean.sd,linear.x, constant.x,sampling.f=31/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.wean.c$par[1],l.m=0, fit.wean.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, -Inf, 0,0)

fit.wean.l <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.wean.l


# Sigmoid
lik.q <- make.quasse(wean.phy,wean,wean.sd,sigmoid.x, constant.x,sampling.f=31/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.wean.c$par[1],fit.wean.c$par[1],mean(xr),1, fit.wean.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(wean), -Inf, 0,0)

fit.wean.s <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.wean.s


# Hump
lik.q <- make.quasse(wean.phy, wean,wean.sd,noroptimal.x, constant.x,sampling.f=31/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.wean.c$par[1],fit.wean.c$par[1],mean(xr),1, fit.wean.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(wean), -Inf, 0,0)

fit.wean.h <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.wean.h

anova(fit.wean.c,wean.linear=fit.wean.l,wean.sigmoid=fit.wean.s,wean.hump=fit.wean.h)





# Comparable body mass to gestation length
# Load in data for species name, gestation length and body mass excluding species with missing values
d16 <- d1 %>% dplyr::select(binomial.name, mass_sex_mean,gest_m)
d17 <- comparative.data(phy1, d16, binomial.name,na.omit=T)

mass2 <- setNames(log(d17$data$mass_sex_mean),rownames(d17$data))
mass2.phy <- d17$phy

mass2.sd <- 1/50

p <- starting.point.quasse(mass2.phy, mass2)
p

# Constant
lik.q <- make.quasse(d17$phy, mass2,mass2.sd,constant.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(p[1], p[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, 0)

fit.mass2.c <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.mass2.c

# Linear
xr <- range(mass2) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
lik.q <- make.quasse(mass2.phy, mass2,mass2.sd,linear.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.mass2.c$par[1],l.m=0, fit.mass2.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, -Inf, 0,0)

fit.mass2.l <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.mass2.l


# Sigmoid
lik.q <- make.quasse(mass2.phy,mass2,mass2.sd,sigmoid.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.mass2.c$par[1],fit.mass2.c$par[1],mean(xr),1, fit.mass2.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(mass2), -Inf, 0,0)

fit.mass2.s <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.mass2.s


# Hump
lik.q <- make.quasse(mass2.phy, mass2,mass2.sd,noroptimal.x, constant.x,sampling.f=43/114)

lik.nodrift <- constrain(lik.q, drift ~ 0)
argnames(lik.nodrift)

p.start <- c(fit.mass2.c$par[1],fit.mass2.c$par[1],mean(xr),1, fit.mass2.c$par[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start

lower <- c(0, 0, min(mass2), -Inf, 0,0)

fit.mass2.h <- find.mle(lik.nodrift, p.start, control=list(parscale=.1),lower=lower, verbose=0)
fit.mass2.h

anova(fit.mass2.c,mass2.linear=fit.mass2.l,mass2.sigmoid=fit.mass2.s,mass2.hump=fit.mass2.h)









# BAMM analysis
library(BAMMtools)
library(readxl)
library(ape)
library(phytools)
library(dplyr)
library(caper)


# Make ultrametric version of tree file
ape::write.tree(phy1, file="ultrametric_lemur_tree.txt")


# Get BAMM priors to use
setBAMMpriors(read.tree("ultrametric_lemur_tree.txt"))

# BAMMtools to make a control file for speciation-extinction model using priors from above
generateControlFile('BAMMtools_control.txt', type = 'diversification', params = list(
  treefile = 'ultrametric_lemur_tree.txt',
  globalSamplingFraction = '1',
  numberOfGenerations = '1000000',
  overwrite = '1',
  lambdaInitPrior = '2.45265812652213',
  lambdaShiftPrior = '0.0232203854022866',
  muInitPrior = '2.45265812652213',
  expectedNumberOfShifts = '1'))


# Load tree and edata 
tree <- read.tree("ultrametric_lemur_tree.txt")
edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.1)

# Log likelihood trace of MCMC to look for convergence
mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

# Discard first 10% as burnin
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# Check the effective sample sizes of the log-likelihood and the number of shift events present in each sample
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

# Compute the posterior probabilities of models sampled 
post_probs <- table(postburn$N_shifts) / nrow(postburn)
# See which models are part of the set that were sampled
names(post_probs)

# summarize the posterior distribution of the number of shifts
shift_probs <- summary(edata)

# Check effect of prior using BayesFactor
postfile <- "mcmc_out.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)

# mean phylorate plot
plot.bammdata(edata, lwd=2)
plot.bammdata(edata, lwd=2, legend=T)


# 95% credible set of distinct shift configurations
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
summary(css)

# phylorate plots for each of the N shift configurations with the highest posterior probabilities
plot.credibleshiftset(css)

# Finding best shift configuration (ik only 1)
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)

# visualize rate shifts using plot.phylo
mysample <- 25
nrow(edata$eventData[[ mysample ]])
# only one rate regime, then you have no rate shifts: the single rate regime starts at the root 
# and describes the entire tree

# summarize marginal shift probabilities
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs,cex=0.3)

# pull out the mean rates for individual branches
?getMeanBranchLengthTree()
getMeanBranchLengthTree(edata,"speciation")
getMeanBranchLengthTree(edata,"extinction")
getMeanBranchLengthTree(edata,"ndr")
# Can also see for trait


plotRateThroughTime(edata, ratetype="speciation")
rtt <- getRateThroughTimeMatrix(edata)
meanTraitRate <- colMeans(rtt$beta)
plot(meanTraitRate ~ rtt$times)


# Macroevolutionary cohort analysis provides a way of summarizing the extent to which 
# species share correlated macroevolutionary dynamics
cmat <- getCohortMatrix(edata)
cohorts(cmat, edata)

# cumulative shift probability tree shows the cumulative probability, on each branch, that a 
# shift occurred somewhere between the focal branch and the root of the tree
cst <- cumulativeShiftProbsTree(edata)
plot.phylo(cst,cex=0.5)
cst <- cumulativeShiftProbsTree(edata)
edgecols <- rep('black', length(phy1$edge.length))
is_highprobshift <- cst$edge.length >= 0.95
edgecols[ is_highprobshift ] <- "red"
plot.phylo(phy1, edge.color = edgecols)







# Control file body mass trait data
dat <- read_excel("lemur trait data.xlsx")
dat1 <- as.data.frame(dat)
dat2 <- dat1 %>% dplyr::select(binomial.name, mass_sex_mean)
comparative_mass <- comparative.data(phy1, dat2, binomial.name,na.omit=T)
ape::write.tree(comparative_mass$phy, file='bodymass_phylogeny.txt')
setBAMMpriors(read.tree("bodymass_phylogeny.txt"), traits = "log_bodymass.txt")
# Run BAMM using "log_bodymass.txt" trait file

generateControlFile('control_body.txt', type = 'trait', params = list(
  treefile = 'ultrametric_lemur_tree.txt',
  traitfile = 'log_bodymass.txt',
  numberOfGenerations = '2000000000',
  overwrite = '1',
  expectedNumberOfShifts = '1',
  betaInitPrior = '6.17610204312266',
  betaShiftPrior = '0.0232203854022866',
  useObservedMinMaxAsTraitPriors = '1'))

# Run BAMM on C++. Analysis of output:
tree <- read.tree("ultrametric_lemur_tree.txt")
edataB <- getEventData(tree, eventdata = "event_data.txt",type='trait', burnin=0.1)


mcmcout <- read.csv("mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

post_probs <- table(postburn$N_shifts) / nrow(postburn)

names(post_probs)
shift_probs <- summary(edataB)

postfile <- "mcmc_out.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)

computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)

plotPrior(postfile, expectedNumberOfShifts=1)


plot.bammdata(edataB, lwd=2, legend=T)
addBAMMshifts(edataB, cex=2)

css <- credibleShiftSet(edataB, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css)
best <- getBestShiftConfiguration(edataB, expectedNumberOfShifts=1)
dev.off()
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)

msc.set <- maximumShiftCredibility(edataB, maximize='product')
msc.config <- subsetEventData(edataB, index = msc.set$sampleindex)
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)

marg_probs <- marginalShiftProbsTree(edataB)
plot.phylo(marg_probs,cex=0.3)

allrates <- getCladeRates(edataB)

# Trait evolution through time plot
plotRateThroughTime(edataB, ratetype="trait")

rtt <- getRateThroughTimeMatrix(edataB)
meanTraitRate <- colMeans(rtt$beta)
plot(meanTraitRate ~ rtt$times)

cmat <- getCohortMatrix(edataB)
cohorts(cmat, edataB)

cst <- cumulativeShiftProbsTree(edataB)
plot.phylo(cst,cex=0.3)

cst <- cumulativeShiftProbsTree(edataB)
edgecols <- rep('black', length(tree$edge.length))
is_highprobshift <- cst$edge.length >= 0.95
edgecols[ is_highprobshift ] <- "red"
plot.phylo(tree, edge.color = edgecols,cex=0.23)

# Code for putting family labels on the outside of fan shaped phylogeny
arc.cladelabels<-function(tree=NULL,text,node,ln.offset=1.02,
                          lab.offset=1.06,cex=1,orientation="curved",...){
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(obj$type!="fan") stop("method works only for type=\"fan\"")
  h<-max(sqrt(obj$xx^2+obj$yy^2))
  if(hasArg(mark.node)) mark.node<-list(...)$mark.node
  else mark.node<-TRUE
  if(mark.node) points(obj$xx[node],obj$yy[node],pch=21,
                       bg="red")
  if(is.null(tree)){
    tree<-list(edge=obj$edge,tip.label=1:obj$Ntip,
               Nnode=obj$Nnode)
    class(tree)<-"phylo"
  }
  d<-getDescendants(tree,node)
  d<-sort(d[d<=Ntip(tree)])
  deg<-atan(obj$yy[d]/obj$xx[d])*180/pi
  ii<-intersect(which(obj$yy[d]>=0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]>=0))
  deg[ii]<-360+deg[ii]
  draw.arc(x=0,y=0,radius=ln.offset*h,deg1=min(deg),
           deg2=max(deg))
  if(orientation=="curved")
    arctext(text,radius=lab.offset*h,
            middle=mean(range(deg*pi/180)),cex=cex)
  else if(orientation=="horizontal"){
    x0<-lab.offset*cos(median(deg)*pi/180)*h
    y0<-lab.offset*sin(median(deg)*pi/180)*h
    text(x=x0,y=y0,label=text,
         adj=c(if(x0>=0) 0 else 1,if(y0>=0) 0 else 1),
         offset=0)
  }
}

# Packages needed for fan shaped phylogeny plot
library(phytools)
library(plotrix)

tree

# Specify most recent common ancestor nodes
nodes<-c(119,151,179,204,227)
# Labels for families
labels<-paste(c("Cheirogaleidae","Lepilemuridae","Indriidae","Lemuridae","Daubentoniidae"))
# Plot fan-shaped tree
plot.bammdata(edataB,method="polar",lwd=2,legend=T,par.reset=FALSE)
# Add BAMM shifts
addBAMMshifts(edataB,method="polar",bg='black',cex=1,par.reset=FALSE)
# Add family labels
for(i in 1:length(nodes)) 
  arc.cladelabels(text=labels[i],node=nodes[i])

# All credible rate shift configurations
cssB <- credibleShiftSet(edataB, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
# How many?
cssB$number.distinct
summary(cssB)
plot.credibleshiftset(cssB)
# Get the most likely one
bestB <- getBestShiftConfiguration(edataB, expectedNumberOfShifts=1)
dev.off()

# Plot most likely one like above
plot.bammdata(bestB,method="polar",lwd=2,legend=T,par.reset=FALSE)
addBAMMshifts(bestB,method="polar",bg='grey',cex=1.5,par.reset=FALSE)
for(i in 1:length(nodes)) 
  arc.cladelabels(text=labels[i],node=nodes[i])


# Gestation period BAMM analysis
# Repeat steps above
# Control file log gestation period trait data
dat <- read_excel("lemur trait data.xlsx")
dat1 <- as.data.frame(dat)
dat2 <- dat1 %>% dplyr::select(binomial.name, gest_m)
comparative_gest <- comparative.data(phy1, dat2, binomial.name,na.omit=T)
ape::write.tree(comparative_gest$phy, file='gest_phylogeny.txt')
setBAMMpriors(read.tree("gest_phylogeny.txt"), traits = "log_gestation.txt")



generateControlFile('control_gestation.txt', type = 'trait', params = list(
  treefile = 'ultrametric_lemur_tree.txt',
  traitfile = 'log_gestation.txt',
  numberOfGenerations = '1000000000',
  overwrite = '1',
  expectedNumberOfShifts = '1',
  betaInitPrior = '66.4209974387857',
  betaShiftPrior = '0.0232203854022866',
  useObservedMinMaxAsTraitPriors = '1'))

# Repeat above steps for gestation length with "log_gestation.txt" as trait file


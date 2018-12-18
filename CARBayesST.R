library(CARBayesST)

##########################################################




library(rstan)
library(Matrix)

N.mat=readMM("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/neighborhood.mtx")

pairs<-triu(N.mat)

pairs<-summary(pairs)[summary(pairs)$x==1,]

zs=as.matrix(read.csv("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/crime.csv")[,-1])

size=nrow(N.mat)

num.obs=ncol(zs)

#Population by Census Block Group
pop=read.csv("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/pop.csv")[,-1]
#Percentage Unemployed by Census Block Group
un.emp=read.csv("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/unemp.csv")[,-1]
#Centered and Scaled Average Family Income by Census Block Group (2015 Dollars)
wealth.std=read.csv("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/wealth.csv")[,-1]
#Number of Young Males by Census Block Group (15-20 yr olds)
ym=read.csv("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/ym.csv")[,-1]
ym=(ym-mean(ym))/sqrt(var(ym))
temperature=c(29,23.4,44.4,60.0,70.0,76.9,81.9,80.9,78.1,63.2,53.7,44.5)
dim(temperature)=c(12,1)
temperature=(temperature-mean(temperature))/as.numeric(sqrt(var(temperature)))
#un.emp=(un.emp-mean(un.emp))/as.numeric(sqrt(var(un.emp)))
wealth=wealth.std

temp.cov=rep(temperature,6)

trend.cov=1:72

trend.cov=(trend.cov-mean(trend.cov))/(sqrt(var(trend.cov)))

x1<-rep(log(pop),72)
x2<-rep(un.emp,72)
x3<-rep(wealth,72)
x4<-rep(ym,72)
y<-as.numeric(zs)

Xs=formula(y~x1+x2+x3+x4)


model1<-ST.CARsepspatial(formula=Xs, family="poisson",W=as.matrix(N.mat), burnin=2000, n.sample=7000,thin=1)


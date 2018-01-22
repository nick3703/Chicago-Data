#################MAY NEED TO RUN#######################
install.packages(c("rstan","Matrix"))
###############REQUIRED PACKAGES########################
library(rstan)
library(Matrix)

######################READ IN DATA###########################

#All Data Available at https://github.com/nick3703/Chicago-Data

N.mat=readMM("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/neighborhood.mtx")

pairs<-triu(N.mat)

pairs<-summary(pairs)[summary(pairs)$x==1,]

zs=as.matrix(read.csv("https://raw.githubusercontent.com/nick3703/Chicago-Data/master/crime.csv")[,-1])

dat<-as.numeric(zs)

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
#####################STANDARDIZE DATA################################################
ym=(ym-mean(ym))/(sqrt(var(as.numeric(ym))))
un.emp=(un.emp-mean(un.emp))/sqrt(var(as.numeric(un.emp))) #Standardize Covariates
Xs<-model.matrix(~1+log(pop)+ym+wealth.stnd+un.emp)

temperature=c(29,23.4,44.4,60.0,70.0,76.9,81.9,80.9,78.1,63.2,53.7,44.5) #Temperature is proxy for seasonality
dim(temperature)=c(12,1)
temperature=(temperature-mean(temperature))/sqrt(var(as.numeric(temperature)))
temp.cov=rep(temperature,6)
trend.cov=1:72
trend.cov=(trend.cov-mean(trend.cov))/sqrt(var(trend.cov))

##################CREATE INGARCH MODEL IN STAN#############################

ingarchmodel="
data {
int<lower = 1> n;
int<lower = 1> t;
int<lower = 1>tot;
int<lower = 1>p;
matrix[n, p] X;
int<lower = 0> y[tot];
matrix[n,t] dat;
vector[t] trend;
vector[t] temp;
}
parameters {
vector[p] beta;
real<lower = 0, upper = 1> etacross;
real<lower = 0, upper = 1> eta;

}
transformed parameters{
matrix<lower = 0>[n,t] lams;
vector[n] xpred;
xpred = X*beta;
lams[,1]=exp(xpred+.164*temp[1]-.272*trend[1]);  //coefficients on temp and trend fix at MLE values
for(j in 2:t){
lams[,j] = eta*dat[,j-1]+(etacross)*lams[,j-1]+exp(xpred+.164*temp[j]-.272*trend[j]);
}

}

model {

beta ~ normal(0, 10);
y ~ poisson(to_vector(lams));
eta~beta(.5,.5);
etacross~beta(.5,.5);

}"


sp_d <- list(n = size,         # number of total observations
             t=num.obs,
             tot=size*num.obs,
             p=ncol(Xs),
             X = Xs,               # design matrix
             y = as.numeric(zs),                # observed number of cases      
             dat=zs,
             trend=trend.cov,
             temp=temp.cov

)


###########SAMPLE FROM STAN#########

ing<-stan_model(model_code=ingarchmodel)
niter=5000
nchains=1
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
sp_fit_ingar <- stan(model_code=ingarchmodel, data = sp_d, 
                    iter = niter, chains = nchains,verbose = FALSE,refresh=50)

print(sp_fit_ingar, pars = c('beta', 'etacross', 'eta','lp__'))
INGAR=extract(sp_fit_ingar,pars = c('beta', 'etacross', 'eta','lp__'))
write.csv(INGAR,"INGARchi.csv")
samps=read.csv("INGARchi.csv")[,-1]  #Output saved as 'samps'



m=1:nrow(samps)
###################SIMULATING FROM POSTERIOR & CALCULATING POSTERIOR PREDICTIVE VALUES######################

iden<-Matrix(0,nrow=num.obs,ncol=num.obs)
diag(iden)<-1


spat.Neigh<- iden %x% spat

iden2<-Matrix(0,nrow=size,ncol=size)
diag(iden2)=1
tempmat=Matrix(0,nrow=72,ncol=72)

for(i in 1:71){tempmat[i,i+1]=tempmat[i+1,i]=1}

temptest=iden2%x%tempmat

list.Neigh=mat2listw(spat.Neigh)

ym=(ym-mean(ym))/sqrt(var(ym))
means=c()
vars=c()
mors=c()
ar=c()
test.ar=c()
ar2=c()
varmean=c()
for(l in 1:5000){
  o=sample(m,1)
  eta=samps$eta[sample(m,1)]
  etacross=samps$etacross[sample(m,1)]
  b1=samps$beta.1[sample(m,1)]
  b2=samps$beta.2[sample(m,1)]
  b3=-.272
  b4=.164
  b6=samps$beta.4[sample(m,1)]
  b5=samps$beta.3[sample(m,1)]
  b7=samps$beta.5[sample(m,1)]
  size=ncol(N.mat)
  num.neighs<-apply(N.mat,1,sum)
  xi=rep(1,(size))
  delt=exp(b1*xi+b2*log(pop)+b3*trend.cov[1]+b4*temperature[1]+b5*ym+b6*wealth+b7*un.emp)
  lam2=matrix(delt,ncol=1)
  zs.test=matrix(0,nrow=size,ncol=1)
  zs.test[,1]=rpois(size,lam2[,1])
  q=rep(1:12,6)
  for(i in 2:num.obs){

    v=exp(b1*xi+b2*log(pop)+b3*trend.cov[i]+b4*temperature[q[i]]+b5*ym+b6*wealth+b7*un.emp)+etacross*lam2[,(i-1)]+eta*zs.test[,(i-1)]
    lam2=cbind(lam2,v)
    z=rpois(size,as.vector(v))
    zs.test=cbind(zs.test,z)
  }
  means[l]=mean(zs.test)
  vars[l]=var(as.numeric(zs.test))
  varmean[l]=log(vars[l]/means[l])
  cor.mat=cor(zs.test)
  lag1=matrix(0,nrow=num.obs,ncol=num.obs)
  for(k in 1:(num.obs-1)){lag1[k,k+1]=1}
  lag1[num.obs,(num.obs-1)]=1
  ar2[l]=sum(lag1*cor.mat)/num.obs
  mors[l]=moran.test(as.numeric(zs.test),list.Neigh)$estimate[[1]]
}

acf.base=c()
for(k in 1:size){
  acf.base=acf(zs[k,],lag.max=2,plot=F)[[1]][2]
}
base.ar=mean(acf.base)
cor.mat=cor(zs)
lag1=matrix(0,nrow=72,ncol=72)
for(k in 1:71){lag1[k,k+1]=1}
lag1[72,71]=1
base.ar=sum(lag1*cor.mat)/72
base.mor=moran.test(as.numeric(zs),list.Neigh)$estimate[[1]]
base.varmean=log(var(as.numeric(zs))/mean(zs))

p1=length(ar2[ar2>base.ar])/length(ar2)
p2=length(mors[mors>base.mor])/length(mors)
p3=length(varmean[varmean>base.varmean])/length(varmean)


boxplot(ar2)
abline(h=base.ar,col="red")


boxplot(mors)
abline(h=base.mor,col="red")


boxplot(varmean)
abline(h=base.varmean,col="red")




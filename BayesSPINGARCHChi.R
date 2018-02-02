

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
temperature=c(29,23.4,44.4,60.0,70.0,76.9,81.9,80.9,78.1,63.2,53.7,44.5)
dim(temperature)=c(12,1)
temperature=(temperature-mean(temperature))/as.numeric(sqrt(var(temperature)))
#un.emp=(un.emp-mean(un.emp))/as.numeric(sqrt(var(un.emp)))
wealth=wealth.std

temp.cov=rep(temperature,6)

trend.cov=1:72

trend.cov=(trend.cov-mean(trend.cov))/(sqrt(var(trend.cov)))

Xs=model.matrix(~1+log(pop)+un.emp+wealth+ym)
model4="functions {
real sparse_car_lpdf(vector phi, real tau, real alpha, 
int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
row_vector[n] phit_D; // phi' * D
row_vector[n] phit_W; // phi' * W
vector[n] ldet_terms;

phit_D = (phi .* D_sparse)';
phit_W = rep_row_vector(0, n);
for (i in 1:W_n) {
phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
}

for (i in 1:n) ldet_terms[i] =log1m(alpha * lambda[i]);
return 0.5 * (n * log(tau)
+ sum(ldet_terms)
- tau * (phit_D * phi - alpha * (phit_W * phi)));
}
}
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
int<lower =1> W_n;
int W_sparse[W_n, 2];   // adjacency pairs
vector[n] D_sparse; // number of adjacent region pairs
vector[n] lambda;
real alpha;
}
transformed data {
matrix<lower=0,upper=1> [n,n] D;
vector[n] zeros;
D=diag_matrix(rep_vector(1.0, n));
zeros = rep_vector(0, n);
}
parameters {
vector[n] phi;
matrix[n,t] phi2;
real<lower = 0> std;
real<lower = 0> stds;
vector[p] beta;
real<lower = 0, upper = 1> etacross;
real<lower = 0, upper = 1> eta;
}
transformed parameters{
matrix <lower = 0> [n,t] lams;
vector[n] xpred;
real<lower =0> tau;
tau = 1/std^2;
xpred = X*beta;
lams[,1]=exp(xpred+phi+.164*temp[1]-.272*trend[1]+phi2[,1]);
for(j in 2:t){
lams[,j] = eta*dat[,j-1]+(etacross)*lams[,j-1]+exp(xpred+.164*temp[j]-.272*trend[j]+phi2[,j]+phi);
}
}

model {



to_vector(phi) ~ sparse_car(tau, alpha, W_sparse, D_sparse, lambda, n, W_n);



for(k in 1:n){
for(l in 1:t){
phi2[k,l] ~ normal(0,stds);
}
}

beta ~ normal(0, 10);
y ~ poisson(to_vector(lams));
std ~ normal(0,5);
stds ~ normal(0,5);
}"

spat<-N.mat
neigh.hood.mat<-matrix(0,nrow=nrow(N.mat),ncol=ncol(N.mat))
diag(neigh.hood.mat)<-apply(N.mat,1,sum)

sq.N.inv<-sqrt(solve(neigh.hood.mat))

eigs<-eigen(sq.N.inv%*%N.mat%*%sq.N.inv)$values

sp.size=size

num.neighs<-apply(N.mat,1,sum)

place.hold<-Matrix(0,nrow=size,ncol=size)
iden.space=place.hold
diag(iden.space)=1

N=place.hold
diag(N)=1/num.neighs



sp_spat <- list(n = size,         # number of total observations
                t=num.obs,
                tot=size*num.obs,
                p=ncol(Xs),
                X = Xs,               # design matrix
                y = as.numeric(zs),   # observed number of cases      
                dat=zs,
                trend=trend.cov,
                temp=temp.cov,   W=matrix(spat,nrow=nrow(spat),ncol=ncol(spat)),
                W_n = sum(spat)/2,    # number of neighbor pairs
                W_sparse = pairs[,1:2],
            		D_sparse=num.neighs,
                lambda=eigs,
	            	alpha=.9999 #Fixed Spatial Parameter at Edge of parameter space
)

m<-stan_model(model_code=model4)



niter=7000
nchains=3
rstan_options(auto_write = TRUE)
sp_fit.spat <- stan(model_code=model4, data = sp_spat, 
                    iter = niter, chains = nchains, verbose = FALSE, refresh=50)

print(sp_fit.spat, pars = c('stds','beta', 'tau', 'etacross', 'eta','lp__'))



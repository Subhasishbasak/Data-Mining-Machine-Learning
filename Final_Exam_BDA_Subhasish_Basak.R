### Bayesian Data Analysis 
### Final Exam
### Total Score : 35
### Instructor : Sourish Das

### Student Name: Subhasish Basak
### Roll Number: MDS201803
### Email: subhasish@cmi.ac.in

### Please save the file as Final_Exam_BDA_Sourish_Das.R

### Last Date of Submission : 30th April, 2020

### Following data is downloading the global COVID-19 data repository from the John Hopkins University

### You have two tasks to complete in this exam.

### Task 1:
### Following analysis fits linear regression model on log-differences of confirmed cases using OLS method

### Run a Bayesian analysis to check if lockdown is effective

## Data loading
data<-read.csv(file='https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv',head=TRUE)

## dates handling
dts1<-dts<-colnames(data)[-c(1:4)]
for(j in 1:length(dts))dts1[j]<-strsplit(dts[j],"X")[[1]][2]
dats<-as.Date(dts1,format = "%m.%d.%y")

# Country wise configurations

# start_date : the date when the first incidence occured
# lockdown_date : the date when the lockdown started

country_name = 'India'
country_code = '132'
start_date = '2020-03-04'
lockdown_date = '2020-03-25'

country_name = 'Italy'
country_code = '138'
start_date = '2020-01-31'
lockdown_date = '2020-03-09'

country_name = 'Spain'
country_code = '202'
start_date = '2020-02-01'
lockdown_date = '2020-03-13'


## Data extraction and Pre-processing

country_data<-data[data$Country.Region==country_name,]
country_data1<-cbind.data.frame(Date=dats,Cases=t(country_data[-c(1:4)]))
names(country_data1)[names(country_data1) == country_code] <- "incid"
India_incid<-diff(country_data1$incid)
n<-length(India_incid)  

country_data2<-cbind.data.frame(Dates=dats[2:length(dats)],India_incid=India_incid,Cases=country_data1$incid[2:(n+1)])
colnames(country_data2)<-c("Dates","Incidence","Total_Confirmed_Cases")
country_data2<-subset(country_data2,Dates>=start_date)

country_data2$Time<-1:nrow(country_data2)
country_data2$lock_down<-0
country_data2$lock_down[country_data2$Dates>=lockdown_date]<-1


## Computing the log differences of the incidences

country_data2$ln_t[2:nrow(country_data2)]<-diff(log(country_data2$Total_Confirmed_Cases))

## Fitting regression through OLS

fit<-lm(ln_t~Time+I(Time^2)+lock_down
        ,data=country_data2)

summary(fit)

out <- capture.output(summary(fit))

cat(paste0("OLS results for ",country_name), out, file=paste0(country_name,"_OLS_results.txt"), sep="\n", append=TRUE)


### This analysis indicates that lockdown may reduce the disease progression - but it is not statistically significant.


plot(country_data2$Dates,country_data2$ln_t,type="b",pch=20,ylab = expression(log(C[t]/C[t-1])),xlab="Dates", main = "log diff of incidence over time", sub = country_name)
### This plot indicates from 7th-April, 2020 (exactly 14 days from the lockdown starts - something has happened.)


### Checking Normality of the errors
shapiro.test(fit$residuals)

par(mfrow=c(1,2))

# Normal QQ-plot
qqnorm(fit$residuals
       ,cex=0.3
       ,pch=20
       ,xlim=c(-4,4)
       ,ylim=c(-4,4)
       ,xlab="quantile with fitted Gaussian distribution"
       ,ylab="observed residulas quantiles"
       ,main="Gaussian Q-Q Plot")
abline(a=0,b=1,col="blue")
grid(col="red")

## Thus we observe that the errors are not Gaussian, hence we use the Laplace prior


########################################
# Bayesian Regression with Laplace prior
########################################

rP<-country_data2$ln_t
rP[1]=0
n<-length(country_data2$Time)
Intercept<-rep(1,n)
f<-cbind(Intercept,country_data2$Time,I(country_data2$Time^2),country_data2$lock_down)

## density of Laplace distribution 
dlaplace<-function(x,mu,lambda){
  exp(-abs(x-mu)/lambda)/(2*lambda)
}

## The negative log-likelihood function
neg_log_likelihood<-function(data,para){
  mu<-para[1]
  lambda<-(para[2])
  f<-dlaplace(x=data
              ,mu=mu
              ,lambda=lambda)
  return(-sum(log(f)))
}

## initial parameters
para.init<-c(0,1)

## residual or idiosyncratic return 
## from capm_portf
e<-summary(fit)$residual

mle<-optim(para.init
           ,neg_log_likelihood
           ,data=e)

# MLE estimates of the Laplace parameters
mle$par

## The following function simulate from 
## the fitted Laplace distribution
rlaplace<-function(n,mu,lmbda){
  # n : sample size
  # mu: location parameter
  # lambda : scale parameter
  u<-runif(n,-0.5,0.5)
  sgn<-rep(NA,n)
  for(i in 1:n){
    if(u[i]<0) sgn[i]<--1
    if(u[i]>0) sgn[i]<-1
    if(u[i]==0) sgn[i]<-0
  }
  r<-mu-lmbda*sgn*log(1-2*abs(u))
  return(r)
}

## simulate the result
set.seed(4155)
r.sim<-rlaplace(n=1000
                ,mu=mle$par[1]
                ,lmbda=mle$par[2])

## Laplace QQ plot
qqplot(r.sim,e
       ,xlim=c(-4,4)
       ,ylim=c(-4,4)
       ,cex=0.5
       ,pch=20
       ,xlab="quantiles with fitted Laplace distribution"
       ,ylab="observed residual quantiles"
       ,main="Laplace Q-Q Plot")
abline(a=0,b=1,col="blue")
grid(col="red")

# Defining the likelihood

log_likelihood <- function(param,y,x){
  beta_0 = param[1]
  beta_1 = param[2]
  beta_2 = param[3]
  beta_3 = param[4]
  lambda = param[5]
  pred = beta_0*x[,1] + beta_1*x[,2] + beta_2*x[,3] + beta_3*x[,4]
  likelihoods = -log(2*lambda)-abs(y-pred)/lambda
  sumll = sum(likelihoods)
  return(sumll)   
}


log_prior <- function(param,x){
  beta_0 = param[1]
  beta_1 = param[2]
  beta_2 = param[3]
  beta_3 = param[4]
  lambda = param[5]
  beta_0_p = dcauchy(beta_0,0,1,log = T) 
  beta_1_p = dcauchy(beta_1,0,1,log = T) 
  beta_2_p = dcauchy(beta_2,0,1,log = T) 
  beta_3_p = dcauchy(beta_3,0,1,log = T) 
  scale_prior = dgamma(lambda
                       ,1,1
                       ,log = T)
  
  return(beta_0_p+beta_1_p+beta_2_p+beta_3_p+scale_prior)
}

log_posterior <- function(param,y,x){
  like <- log_likelihood(param=param
                         ,y=y,x=x)
  prior <- log_prior(param=param,x=x)
  post  <- like + prior
  return ( post )
}

### Model Fitting with Metropolis-Hastings

proposalfunction <- function(param,x){
  #X=cbind(rep(1,length(x)),x)
  S=param[5]*solve(t(x)%*%x)
  prop<-c(rmvnorm(1
                  ,mean = param[1:4]
                  ,sigma = S)
          ,rgamma(1,param[5]*5,5))
  return(prop)
}

run_metropolis <- function(startvalue, N.sim, burnin){
  iterations = N.sim + burnin
  chain = array(dim = c(iterations+1,5))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,],x=x)
    
    probab = exp(log_posterior(param=proposal
                               ,y=y,x=x) 
                 - log_posterior(param=chain[i,]
                                 ,y=y,x=x))
    
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

y=rP
x=f
startvalue = c(0,1,1,1,0.1)
N.sim=30000
burnin=30000

set.seed(1)
chain = run_metropolis(startvalue=startvalue
                       ,N.sim=N.sim
                       ,burnin=burnin)

colnames(chain)<-c("beta_0","beta_1","beta_2","beta_3","lambda")

chain=chain[(burnin+1):nrow(chain),]
#chain=chain[30001:nrow(chain),]

apply(chain,2,mean)
apply(chain,2,sd)
apply(chain,2,quantile,prob=c(0.025
                              ,0.975))

## P(alpha > 0 | Data)
length(chain[chain[,"beta_0"]>0
             ,"beta_0"])/nrow(chain)


## Plotting the posterior densities and the MH simulation Chain

par(mfrow=c(2,2))

plot(ts(chain[,"beta_0"])
     ,main = "Simulated beta_0"
     ,ylab = "beta_0")
plot(density(chain[,"beta_0"])
     ,main = "Posterior density of beta_0"
     ,xlab = "beta_0")
plot(ts(chain[,"beta_1"])
     ,main = "Simulated beta_2"
     ,ylab = "beta_1")
plot(density(chain[,"beta_1"])
     ,main = "Posterior density of beta_1"
     ,xlab = "beta_1")
mtext(country_name, side = 3, line = -1, outer = TRUE)


par(mfrow=c(2,2))
plot(ts(chain[,"beta_2"])
     ,main = "Simulated beta_2"
     ,ylab = "beta_2")
plot(density(chain[,"beta_2"])
     ,main = "Posterior density of beta_2"
     ,xlab = "beta_2")
plot(ts(chain[,"beta_3"])
     ,main = "Simulated beta_3"
     ,ylab = "beta_3")
plot(density(chain[,"beta_3"])
     ,main = "Posterior density of beta_3"
     ,xlab = "beta_3")
mtext(country_name, side = 3, line = -1, outer = TRUE)


par(mfrow=c(1,2))
plot(ts(chain[,"lambda"])
     ,main = "Simulated sigma"
     ,ylab = "sigma")
plot(density(chain[,"lambda"])
     ,main = "Posterior density of sigma"
     ,xlab = "sigma")
mtext(country_name, side = 3, line = -1, outer = TRUE)


## Computing the summary statistics

head(round(chain,digits = 3))
post_median<-apply(chain,2,median)
post_mean <- apply(chain,2,mean)
post_sd   <- apply(chain,2,sd)
post_band <- apply(chain,2,quantile
                   ,prob=c(0.025,0.975))
post_sumry<-rbind( post_median     
                   ,post_mean
                   ,post_sd
                   ,post_band)
rownames(post_sumry)<-c("median"
                        ,"mean"
                        ,"sd"
                        ,"2.5%"
                        ,"97.5%")
round(post_sumry,digits = 4)


out <- capture.output(round(post_sumry,digits = 4))

cat(paste0("Bayesian Regression (using MH algo with Laplace prior) model fitting summary for ",country_name), out, file=paste0(country_name,"_Bayesian_laplace_results.txt"), sep="\n", append=TRUE)

#!/usr/local/bin/Rscript

graphics.off()
rm(list=ls())
source("Plots_tables.R")
library("matrixStats")
data <- read.csv("data2.csv", sep=",",head=TRUE)
data$smoltsurvival <-data$smoltsurvival
data$SpawnerYr <- data$Year
data$predfemalesmolt_pre_dam <- data$SpawnerYr*0
data$predfemalesmolt_post_dam <- data$SpawnerYr*0
data$OceanAdults <- data$SpawnerYr*0

#print(data)

Percent4 <- 68.2/100
Percent5 <- 31.2/100
Percent6 <- .6/100
Repeat6 <- 16/100
DPSDPE <- .74
SF <- 2.75/100
initsize <- 419
years <- 38
                                        #we want to estimate SF for a standard DPSDPE




modelinit <- function(SF1,SF2,init,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
{ #this estimates initial abundance 
ages <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)
data$predSmoltAParameter<- 2078+data$SpawnerYr*0
data$predfemalesmolt_pre_dam[1:7] <-SF1* data$predSmoltAParameter[1:7]*(init+0*data$FosterCount[1:7])/(1+data$predSmoltAParameter[1:7]*SF1*(init+0*data$FosterCount[1:7])/100000.0)
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*data$passprob[1:7]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*data$smoltsurvival[1:7]*.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
ages$SpawnerAdult[1:7]=init

for(i in 8:(years))
{
    ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
    if(i<=(18+7)){
        data$predfemalesmolt_pre_dam[i]=data$predSmoltAParameter[i]*(SF1*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*(SF1*ages$SpawnerAdult[i]))/100000.0)
    }
        if(i>=(19+7)){
        data$predfemalesmolt_pre_dam[i]=data$predSmoltAParameter[i]*(SF2*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*(SF2*ages$SpawnerAdult[i]))/100000.0)
        }
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*data$passprob[i]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*data$smoltsurvival[i]*.74
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
    ages$rep7[i]=Repeat6*ages$vir6[i-1]
    ages$Total[i]=ages$vir4[i] +ages$vir5[i] +ages$vir6[i]+ages$rep5[i] +ages$rep6[i]+ages$rep7[i]
}

#print(tail(ages$Total))
eps <- data$FosterCount[17:years]-ages$Total[13:(years-4)]
return(list(data=data,ages=ages,eps=eps))
}




#The most important likelyhood function
#the indices will change depending on data file. 
fninit <- function(par,DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data,nll=2){
    eps <- modelinit(par[1],par[2],par[3],DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data)
    
    if(nll==1){
        eps2 <-(eps$data$FosterCount[8:(years)])-(eps$ages$SpawnerAdult[8:(years)])#split this
    nll <--1*sum(log(dnorm((eps2),0,par[4])))-log(dlnorm(par[3],log(142),(.53)))
    }
    if(nll==2){
        nll <--1*sum(log(dlnorm(eps$data$FosterCount[8:(years)],log(eps$ages$SpawnerAdult[8:(years)]),par[4])))-log(dlnorm(par[3],log(142),(.53)))
        }
    
return(nll)
}



modelinits <- function(SF1,init,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
{ #this estimates initial abundance 
ages <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)


data$predSmoltAParameter<- 2078+data$SpawnerYr*0
data$predfemalesmolt_pre_dam[1:7] <-SF1* data$predSmoltAParameter[1:7]*(init+0*data$FosterCount[1:7])/(1+data$predSmoltAParameter[1:7]*SF1*(init+0*data$FosterCount[1:7])/100000.0)
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*data$passprob[1:7]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*data$smoltsurvival[1:7]*.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
ages$SpawnerAdult[1:7]=init

for(i in 8:(years))
{
  ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
  data$predfemalesmolt_pre_dam[i]=data$predSmoltAParameter[i]*(SF1*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*(SF1*ages$SpawnerAdult[i]))/100000.0)
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*data$passprob[i]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*data$smoltsurvival[i]*.74
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
  ages$Total[i]=ages$vir4[i] +ages$vir5[i] +ages$vir6[i]+ages$rep5[i] +ages$rep6[i]+ages$rep7[i]
}
eps <- data$FosterCount[17:years]-ages$Total[13:(years-4)]
return(list(data=data,ages=ages,eps=eps))
}




#The most important likelyhood function
#the indices will change depending on data file. 
fninits <- function(par,DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data,nll=2){
    eps <- modelinits(par[1],par[2],DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data)
    
    if(nll==1){
        eps2 <-(eps$data$FosterCount[8:(years)])-(eps$ages$SpawnerAdult[8:(years)])#split this
    nll <--1*sum(log(dnorm((eps2),0,par[3])))-log(dlnorm(par[2],log(142),(.53)))
    }
    if(nll==2){
        nll <--1*sum(log(dlnorm(eps$data$FosterCount[8:(years)],log(eps$ages$SpawnerAdult[8:(years)]),par[3])))-log(dlnorm(par[2],log(142),(.53)))
        }
    
return(nll)
}





modelinitdev <- function(SF1,init,dev,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
{ #this estimates initial abundance 
ages <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)


data$predSmoltAParameter<- 2078+data$SpawnerYr*0
data$predfemalesmolt_pre_dam[1:7] <-SF1* data$predSmoltAParameter[1:7]*(init+0*data$FosterCount[1:7])/(1+data$predSmoltAParameter[1:7]*SF1*(init+0*data$FosterCount[1:7])/100000.0)
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*data$passprob[1:7]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*data$smoltsurvival[1:7]*exp(dev[1:7])*.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
ages$SpawnerAdult[1:7]=init

for(i in 8:(years))
{
  ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
  data$predfemalesmolt_pre_dam[i]=data$predSmoltAParameter[i]*(SF1*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*(SF1*ages$SpawnerAdult[i]))/100000.0)
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*data$passprob[i]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*data$smoltsurvival[i]*exp(dev[i])*.74
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
  ages$Total[i]=ages$vir4[i] +ages$vir5[i] +ages$vir6[i]+ages$rep5[i] +ages$rep6[i]+ages$rep7[i]
}

ages$Dev=dev
eps <- data$FosterCount[17:years]-ages$Total[13:(years-4)]
return(list(data=data,ages=ages,eps=eps))
}




#The most important likelyhood function
#the indices will change depending on data file. 
fninitdev <- function(par,DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data,SDdevI=1,PriorSDI=1,nll=2){
    eps <- modelinitdev(par[1],par[2],par[4:length(par)],DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data)

    if(SDdevI==1)
    {
        sddev=2
        }
    if(SDdevI==2)
    {
        sddev=1
        }

        if(SDdevI==3)
    {
        sddev=.5
    }


    if(PriorSDI==1)
    {
        priorsd=.53
        }
    if(PriorSDI==2)
    {
                priorsd=.53*.9
        }

        if(PriorSDI==3)
    {
                priorsd=.53*1.2
    }

    

    if(nll==1){
        eps2 <-(eps$data$FosterCount[8:(years)])-(eps$ages$SpawnerAdult[8:(years)])#split this
    nll <--1*sum(log(dnorm((eps2),0,par[3])))-log(dlnorm(par[2],log(142),(.53)))
    }
    if(nll==2){
        #print(tail(eps$ages$SpawnerAdult))
        nll <--1*sum(log(dlnorm(eps$data$FosterCount[8:(years)],log(eps$ages$SpawnerAdult[8:(years)]),par[3])))-log(dlnorm(par[2],log(142),(priorsd)))-sum(log(dnorm(par[4:(years+3)],0,sddev)))
        }
    
return(nll)
}





print("Yo!")

inputs <- c(.03,120,1,0*rnorm(38,0,1)*0) #initializaton
z=optim(par=inputs,fn=fninitdev,method="BFGS",gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,2,1,nll=2,hessian=TRUE,control=c(reltol=10^-132,maxit=35000))
print(z)
print(eigen(z$hessian)$value)
print(min(eigen(z$hessian)$value))
sdcon <- sqrt(diag(solve(z$hessian)))
print("silly")
print(sdcon)
out <- modelinitdev(z$par[1],z$par[2],z$par[4:length(z$par)],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
print("omg")

## print("Yo!")
## inputs <- c(.03,.03,100,12) #initializaton
## z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,nll=2,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
## print(z)
## print(eigen(z$hessian))
## sdcon <- sqrt(diag(solve(z$hessian)))
## print(sdcon)

simvec <- function(years,rho,SigmaR,start,mean)
{
    
    errorV <- (1:years)*NA
    autoC <- (1:years)*NA
    independentError <- (1:years)*NA
    TotalError <- (1:years)*NA
    simMarineSurvival <- (1:years)*NA
    NormDev <- rnorm(years,0,1)


    errorV[1]=log(start)-log(mean)+(SigmaR^2)/2

    autoC[1] <- rho*errorV[1]

    independentError[1]<- ((1-rho^2)^.5)*NormDev[1]*SigmaR

    TotalError[1] <- autoC[1]+independentError[1]-(SigmaR^2)/2

    simMarineSurvival[1] <- TotalError[1]+log(mean)
    for(i in 2:years)
     {
    errorV[i] <- simMarineSurvival[i-1]-log(mean)+(SigmaR^2)/2
    autoC[i] <- rho*errorV[i]    
    independentError[i]<- ((1-rho^2)^.5)*NormDev[i]*SigmaR
    TotalError[i] <- autoC[i]+independentError[i]-(SigmaR^2)/2
    simMarineSurvival[i] <- TotalError[i]+log(mean)
     }

    return(exp(simMarineSurvival))
    }


simmaker <- function(nsims=100000,years=100,z,sdcon)
{ #simulate marine surivival and freshwater surivival
  set.seed(1999)
    s <- 1:(years)
   ssin <- matrix(nrow=nsims,ncol=years)
    srand <- matrix(nrow=nsims,ncol=years)
    sar <- matrix(nrow=nsims,ncol=years)
    for(i in 1:nsims)
        {
            par <-  c(0.01640660, 25.72874096, 22.55006499,  0.03275851,  0.01337128)
        s[6:(years)]<- (par[1]*sin(par[2]*6:(years)+par[3])+par[4])
        s[6:(years)] <- s[6:years]+rnorm(length(6:(years)),0,(par[5]))
        s[which(s<=0)] <- 0
        ssin[i,] <- s
        ssin[i,1:6] <- par[1]
        srand[i,1:years] <- sample(data$smoltsurvival[1:23],length(1:years),replace=TRUE) #probably needs to be chainged
        par <- c(0.01310638, 0.64901850, 0.01173065, 0.01698823)
        sar[i,1:6] <- par[1]
        }
  
    cov=as.data.frame(solve(z$hessian)) #draw freshwater surival and initial population siz
    tab= MASS::mvrnorm(nsims,z$par,cov)
    fresh1 <- tab[,1]
    init<-  tab[,2]
    dev <- matrix(nrow=nsims,ncol=years)
    for(i in 1:nsims)
        {
            dev[i,1:years] <- sample(z$par[4:(length(z$par)-4)],length(1:years),replace=TRUE)#uncorrelated bootstrap of deviates.
        
  srand[i,1:years]=sample(data$smoltsurvival[1:23],length(1:years),replace=TRUE)*exp(dev[i,1:years])
        }

  start <- data$smoltsurvival[34]
  SigmaR<- 0.96287
  rho <- 0.601
  mean <- 0.0423
  for(i in 1:nsims)
  {
     sar[i,1:34]=data$smoltsurvival[1:34]*exp(z$par[4:(length(z$par)-4)])
     sar[i,35:years]=simvec((years-34),rho,SigmaR,start,mean)
    }

    
    return(list(ssin=ssin,srand=srand,sar=sar,fresh1=fresh1,dev=dev,init=init))
    }

sims=simmaker(nsims=100000,years=70,z,sdcon)




#forcast things into the future
forcastmodel <- function(EIS,mortR,marinescale,years=100,Percent4,Percent5,Repeat5,Repeat6,data,DPE=1,DPS=1,river=1,simindex,sims=sims)
{
    ages <- list(SpawnerYR=1:years,vir4=1:years*0,vir5=1:years*0,vir6=1:years*0,rep6=1:years*0, rep7 =1:years*0,Total=1:years*0,SpawnerAdult=1:years*0)
    
    if(river==1){
        capacityfac=76473
        freshscale=1
        data$FosterCount[1:7]=1:7*0+232
        }
    if(river==2){
        capacityfac=42596
        freshscale=1.16
        data$FosterCount[1:7]=1:7*0+232*.89
    }
        if(river==3){
            capacityfac=112833
            freshscale=1.13
            data$FosterCount[1:7]=1:7*0+mean(data$FosterCount[1:7])*.36
        }

  
 if(mortR==1){
 simforcast=sims$ssin[simindex,]
 }
    if(mortR==2){
         simforcast=sims$srand[simindex,]
      }
    if(mortR==3){
         simforcast=sims$sar[simindex,]
      }

    
    
    ##
    #print("So what is EIS")
    #print(EIS)
if(river==1){    
    if(EIS==1)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt1_FOS_age2.csv",header=TRUE)

    }
    if(EIS==2)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt2a_FOS_age2.csv",header=TRUE)

    }
if(EIS==3)
    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/Alt2b_FOS_age2.csv",header=TRUE)
        }
if(EIS==4)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt3a_FOS_age2.csv",header=TRUE)

    }
if(EIS==5)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt3b_FOS_age2.csv",header=TRUE)
    }
if(EIS==6)
    {
    FOS_Age2 <- read.csv("Summary_Steelhead/Alt4_FOS_age2.csv",header=TRUE)
    }

if(EIS==7)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_FOS_age2.csv",header=TRUE)
    }

if(EIS==8)    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_FOS_age2.csv",header=TRUE)
        FOS_Age2$x.DPE=DPE+0*1:length(FOS_Age2$x.Year)
        FOS_Age2$x.DPS=DPS+0*1:length(FOS_Age2$x.Year)
        
    }
    
if(EIS==0)
    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_FOS_age2.csv",header=TRUE)
        v=(1:length(FOS_Age2$x.Year))*NA
        v[1:length(data$passprob)]=data$passprob
        FOS_Age2$x.DPE=1+0*1:length(FOS_Age2$x.Year)
        FOS_Age2$x.DPS=v

    }
}
if(river==2){    
    if(EIS==1)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt1_GPR_age2.csv",header=TRUE)

    }
    if(EIS==2)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt2a_GPR_age2.csv",header=TRUE)

    }
if(EIS==3)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt2b_GPR_age2.csv",header=TRUE)

    }
if(EIS==4)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt3a_GPR_age2.csv",header=TRUE)

    }
if(EIS==5)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt3b_GPR_age2.csv",header=TRUE)
    }
if(EIS==6)
    {
    FOS_Age2 <- read.csv("Summary_Steelhead/Alt4_GPR_age2.csv",header=TRUE)
    }

if(EIS==7)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_GPR_age2.csv",header=TRUE)
    }

if(EIS==8)
    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/Alt1_GPR_age2.csv",header=TRUE)
        FOS_Age2$x.DPE=DPE+0*1:length(FOS_Age2$x.Year)
        FOS_Age2$x.DPS=DPS+0*1:length(FOS_Age2$x.Year)
        
    }
    
if(EIS==0)
    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_GPR_age2.csv",header=TRUE)
        v=(1:length(FOS_Age2$x.Year))*NA
        v[1:length(data$passprob)]=data$passprob
        FOS_Age2$x.DPE=1+0*1:length(FOS_Age2$x.Year)
        FOS_Age2$x.DPS=v
        }

    }


if(river==2){    
    if(EIS==1)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt1_GPR_age2.csv",header=TRUE)

    }
    if(EIS==2)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt2a_GPR_age2.csv",header=TRUE)

    }
if(EIS==3)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt2b_GPR_age2.csv",header=TRUE)

    }
if(EIS==4)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt3a_GPR_age2.csv",header=TRUE)

    }
if(EIS==5)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt3b_GPR_age2.csv",header=TRUE)
    }
if(EIS==6)
    {
    FOS_Age2 <- read.csv("Summary_Steelhead/Alt4_GPR_age2.csv",header=TRUE)
    }

if(EIS==7)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_GPR_age2.csv",header=TRUE)
    }

if(EIS==8)
    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/Alt1_GPR_age2.csv",header=TRUE)
        FOS_Age2$x.DPE=DPE+0*1:length(FOS_Age2$x.Year)
        FOS_Age2$x.DPS=DPS+0*1:length(FOS_Age2$x.Year)
        
    }
    
if(EIS==0)
    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_GPR_age2.csv",header=TRUE)
        v=(1:length(FOS_Age2$x.Year))*NA
        v[1:length(data$passprob)]=data$passprob
        FOS_Age2$x.DPE=1+0*1:length(FOS_Age2$x.Year)
        FOS_Age2$x.DPS=v
        }

    }



if(river==3){    
    if(EIS==1)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt1_DET_age2.csv",header=TRUE)

    }
    if(EIS==2)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt2a_DET_age2.csv",header=TRUE)

    }
if(EIS==3)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt2b_DET_age2.csv",header=TRUE)

    }
if(EIS==4)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt3a_DET_age2.csv",header=TRUE)

    }
if(EIS==5)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/Alt3b_DET_age2.csv",header=TRUE)
    }
if(EIS==6)
    {
    FOS_Age2 <- read.csv("Summary_Steelhead/Alt4_DET_age2.csv",header=TRUE)
    }

if(EIS==7)
    {
    FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_DET_age2.csv",header=TRUE)
    }

if(EIS==8)
    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_DET_age2.csv",header=TRUE)
        FOS_Age2$x.DPE=DPE+0*1:length(FOS_Age2$x.Year)
        FOS_Age2$x.DPS=DPS+0*1:length(FOS_Age2$x.Year)
        
    }
    
if(EIS==0)
    {
        FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_DET_age2.csv",header=TRUE)
        v=(1:length(FOS_Age2$x.Year))*NA
        v[1:length(data$passprob)]=data$passprob
        FOS_Age2$x.DPE=1+0*1:length(FOS_Age2$x.Year)
        FOS_Age2$x.DPS=v
        }

    }



    ## if(river==1)
    ## {

    ##     NAAMAGIC <- read.csv("./Summary_Steelhead/NAA_FOS_age2.csv",header=TRUE)
    ##     }
    ## if(river==2)
    ## {
    ## NAAMAGIC<- read.csv("./Summary_Steelhead/NAA_GPR_age2.csv",header=TRUE)

    ##     }
    ## if(river==3)
    ## {
    ##     NAAMAGIC <- read.csv("./Summary_Steelhead/NAA_DET_age2.csv",header=TRUE)
        
    ##     }


    
    index <- sample(1:length(FOS_Age2$x.Year),years,replace=TRUE)
    data <- as.list(data)
    length(data$predfemalesmolt_pre_dam) <- years
 #   length(data$predSmoltAParameter) <- years
    data$predSmoltAParameter <- 1:years*0

    data$predSmoltAParameter<- 1:years*0+ 2078
    
#ges <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)
                                        #data$predSmoltAParameter<- 2078+ages$SpawnerYr*0
ages$SpawnerAdult[1:7]=sims$init[simindex]
data$predfemalesmolt_pre_dam[1:7] <- data$predSmoltAParameter[1:7]*(sims$init[simindex]+0*data$FosterCount[1:7])*sims$fresh1[simindex]*freshscale/(1+data$predSmoltAParameter[1:7]*(sims$init[simindex]+0*data$FosterCount[1:7])*sims$fresh1[simindex]*freshscale/capacityfac)
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*data$passprob[1:7]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*simforcast[1:7]*marinescale #*exp(sims$dev[simindex,1:7])*.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
for(i in 8:38)# observed damn passage and observed marine survival 2021
{

  ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
  data$predfemalesmolt_pre_dam[i]=(data$predSmoltAParameter[i]*sims$fresh1[simindex]*freshscale*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*ages$SpawnerAdult[i]*sims$fresh1[simindex]*freshscale/capacityfac))
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*data$passprob[i]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*simforcast[i]*marinescale *.74#*exp(sims$dev[simindex,i])
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
}


for(i in 39:40)# sample damn passage and observed marine survival 2023
{

  ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
  data$predfemalesmolt_pre_dam[i]=(data$predSmoltAParameter[i]*sims$fresh1[simindex]*freshscale*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*ages$SpawnerAdult[i]*sims$fresh1[simindex]*freshscale/capacityfac))
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*sample(data$passprob,1)
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*simforcast[i]*marinescale *.74#*exp(sims$dev[simindex,i])
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
}



    
    
    
    
    for(i in 41:years) #2024 through years  simulated EIS alternatives 
{
    
    ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
    data$predfemalesmolt_pre_dam[i]=(data$predSmoltAParameter[i]*sims$fresh1[simindex]*freshscale*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*ages$SpawnerAdult[i]*sims$fresh1[simindex]*freshscale/capacityfac))
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*FOS_Age2$x.DPS[index[i]]*FOS_Age2$x.DPE[index[i]]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*simforcast[i]*marinescale *.74#*exp(sims$dev[simindex,i])
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
}

#print("crazy")

    ages$Total=ages$vir4 +ages$vir5 +ages$vir6+ages$rep5 +ages$rep6+ages$rep7
eps <- data$FosterCount[17:years]-ages$Total[13:(years-4)]
QETave <-stats::filter(ages$SpawnerAdult, rep(1,5), sides = 1)/5
QETave[is.na(QETave)] <- 0
        
RS <- exp(mean(log(data$OceanAdults[(41):(41+5)]/ages$SpawnerAdult[(41):(41+5)])))
abun <- exp(mean(log(ages$SpawnerAdult[53:(53+15)]),na.rm=TRUE))
return(list(data=data,ages=ages,eps=eps,sim=simforcast[1:years],RS=RS,abund=abun,QET=QETave,DPS=FOS_Age2$x.DPS[index],DPE= FOS_Age2$x.DPE[index],confi=sims$fresh1[simindex]*freshscale,dev=sims$dev[simindex,]))

}













data_forcast <- function(out,sdcon,z,years,EIS,marinescale=1,mortR=2,DPE=1,DPS=1,river=1,ssize=1000,sims=sims)
{ #do lots of simulations

conmat <-matrix(NA,nrow=ssize,ncol=years)
QETamat <-matrix(NA,nrow=ssize,ncol=years)        
simmat <-matrix(NA,nrow=ssize,ncol=years)
devmat <-matrix(NA,nrow=ssize,ncol=years)
        RSmat <-matrix(NA,nrow=ssize,ncol=years)
        QETmat <-matrix(NA,nrow=ssize,ncol=years)
        DPSmat <-matrix(NA,nrow=ssize,ncol=years)
        DPEmat <-matrix(NA,nrow=ssize,ncol=years)
        RSvec <- 1:ssize
        Abunvec <- 1:ssize
    marineSvec <- 1:ssize
    marineSvecShort <- 1:ssize
    
        QETvec <- 1:ssize
confi2 <- 1:ssize
for( i in 1:ssize)
{
    out <- forcastmodel(EIS,mortR,marinescale,years,Percent4,Percent5,Percent6,Repeat6,data,DPE,DPS,river,i,sims)
    conmat[i,1:(years)] <-out$ages$SpawnerAdult[1:years]
    simmat[i,1:(years)] <- out$sim[1:years]
    devmat[i,1:(years)] <- out$dev[1:years]
    QETamat[i,1:(years)] <- out$QET[1:years]
    DPEmat[i,1:(years)] <- out$DPE[1:years]
    DPSmat[i,1:(years)] <- out$DPS[1:years]
    RSmat[i,(8):years] <- out$data$OceanAdults[(8):years]/out$ages$SpawnerAdult[(8):years]
    QETmat[i,] <- ifelse(out$QET[1:(years)]>116,0,1)
    RSvec[i]=out$RS
    Abunvec[i]=exp(mean(log((out$ages$SpawnerAdult[(53):(53+15)]))))
    marineSvec[i] <- mean(out$sim[53:(53+15)])
    marineSvecShort[i] <- mean(out$sim[41:(41+5)],na.rm=TRUE)
    QETvec[i] <- ifelse(min(out$QET[53:(53+15)],na.rm=TRUE)>116,0,1)
    confi2[i] <- out$confi
 
}
return(list(conmat=conmat,simmat=simmat,devmat=devmat,QETamat=QETamat,DPEmat=DPEmat,DPSmat=DPSmat,RSmat=RSmat,QETmat=QETmat,RSvec=RSvec,Abunvec=Abunvec,marineSvec=marineSvec,marineSvecShort=marineSvecShort,QETvec=QETvec,confi=confi2,DPEvec=as.vector(DPEmat),DPSvec=as.vector(DPSmat)))
}


data_forcast2 <- function(out,sdcon,z,years,EIS,marinescale=1,mortR=2,DPE=1,DPS=1,ssize=1000,sims=sims)
{ #do lots of simulations when adding rivers together.  This might be able to be merged with data_forcast, but since the adding of two independent streams is slightly tricky it got broken out.
confi <- rnorm(ssize,z$par[1],sd=sdcon[1])            
conmat <-matrix(NA,nrow=ssize,ncol=years)
QETamat <-matrix(NA,nrow=ssize,ncol=years)        
simmat <-matrix(NA,nrow=ssize,ncol=years)
devmat <-matrix(NA,nrow=ssize,ncol=years)
        RSmat <-matrix(NA,nrow=ssize,ncol=years)
        QETmat <-matrix(NA,nrow=ssize,ncol=years)
        DPSmat <-matrix(NA,nrow=ssize,ncol=years)
        DPEmat <-matrix(NA,nrow=ssize,ncol=years)
        RSvec <- 1:ssize
        Abunvec <- 1:ssize
        marineSvec <- 1:ssize
        marineSvecShort <- 1:ssize
        QETvec <- 1:ssize
    confi2 <- 1:ssize
for( i in 1:ssize)
{
    out <- forcastmodel(EIS,mortR,marinescale,years,Percent4,Percent5,Percent6,Repeat6,data,DPE,DPS,1,i,sims)#foster
    if(EIS<=5){
    out2 <- forcastmodel(EIS,mortR,marinescale,years,Percent4,Percent5,Percent6,Repeat6,data,DPE,DPS,2,i,sims)#greenpeter
    conmat[i,1:years] <-out$ages$SpawnerAdult[1:years]+out2$ages$SpawnerAdult[1:years]
    simmat[i,1:years] <- (out$sim[1:years]+out2$sim[1:years])/2 #average marine survival
    devmat[i,1:years] <- (out$dev[1:years]+out2$dev[1:years])/2
    QETamat[i,1:years] <-(out$QET[1:years]+out2$QET[1:years])
    DPEmat[i,1:years] <- (out$DPE[1:years]+out2$DPE[1:years])/2
    DPSmat[i,1:years] <- (out$DPS[1:years]+out2$DPS[1:years])/2
    RSmat[i,(8):years] <- (out$data$OceanAdults[(8):years]+out2$data$OceanAdults[(8):years])/(out$ages$SpawnerAdult[(8):years]+out2$ages$SpawnerAdult[(8):years])
    QETmat[i,] <- ifelse(out$QET[1:years]+out2$QET[1:years]>116,0,1)
    RSvec[i]=exp(mean(log((out$data$OceanAdults[(41):(41+5)]+out2$data$OceanAdults[(41):(41+5)])/(out$ages$SpawnerAdult[(41):(41+5)]+out2$ages$SpawnerAdult[(41):(41+5)]))))
    Abunvec[i]=exp(mean(log((out$ages$SpawnerAdult[(53):(53+15)]+out2$ages$SpawnerAdult[53:(53+15)]))))
    marineSvec[i] <- mean((out$sim[53:(53+15)]+out2$sim[53:(53+15)]))/2
    marineSvecShort[i] <- mean((out$sim[41:(41+5)]+out2$sim[41:(41+5)]),na.rm=TRUE)/2
    
    QETvec[i] <- ifelse(min(out$QET[43:(53+15)]+out2$QET[53:(53+15)],na.rm=TRUE)>(2*116),0,1)
       confi2[i] <-(out$confi+out2$confi)/2 
    }else
        {
    
    conmat[i,1:years] <-out$ages$SpawnerAdult[1:years] #not all EIS alternatives 
    simmat[i,1:years] <- out$sim[1:years]
    devmat[i,1:years] <- out$dev[1:years]
    QETamat[i,1:years] <- out$QET[1:years]
    DPEmat[i,1:years] <- out$DPE[1:years]
    DPSmat[i,1:years] <- out$DPS[1:years]
    RSmat[i,(8):years] <- out$data$OceanAdults[(8):years]/out$ages$SpawnerAdult[(8):years]
    QETmat[i,] <- ifelse(out$QET[1:years]>116,0,1)
    RSvec[i]=out$RS
    Abunvec[i]=exp(mean(log((out$ages$SpawnerAdult[(53):(53+15)]))))
            marineSvec[i] <- mean(out$sim[53:(53+15)])
                marineSvecShort[i] <- mean(out$sim[41:(41+5)],na.rm=TRUE)
    QETvec[i] <- ifelse(min(out$QET[53:(53+15)],na.rm=TRUE)>116,0,1)
            confi2[i] <-out$confi 
   
}
}
return(list(conmat=conmat,simmat=simmat,devmat=devmat,QETamat=QETamat,DPEmat=DPEmat,DPSmat=DPSmat,RSmat=RSmat,QETmat=QETmat,RSvec=RSvec,Abunvec=Abunvec,marineSvec=marineSvec,marineSvecShort=marineSvecShort,QETvec=QETvec,confi=confi2,DPEvec=as.vector(DPEmat),DPSvec=as.vector(DPSmat)))
}


summaryplot <- function(d1,d2,d3,d4,d5,d6,d7,main="main",flag=1)
{

    if(flag==1)
       { 
    boxplot(d1,d2,d3,d4,d5,d6,d7,
main = main,
at = c(1,2,3,4,5,6,7),
names = c("A1","A2a","A2b","A3a","A3b","A4","NAA"),
las = 2,
horizontal = FALSE,
notch = TRUE,
col="blue",
boxwex=.15
)
       }
    if(flag==2)
{
           l=list(d1,d2,d3,d4,d5,d6,d7)
z=lapply(l,mean,na.rm=TRUE)
print("So z is")
           print(z)
           plot(x=c(1,2,3,4,5,6,7),z,col="blue",xlab=NA, ylab="mean QET",xaxt="n",ylim=c(0,1),pch=16)
           axis(side=1, at=1:7, labels=c("A1","A2a","A2b","A3a","A3b","A4","NAA"),las=2)

           }
    if(flag==3)
{
           l=list(d1,d2,d3,d4,d5,d6,d7)
par(mar=c(12, 4, 4, 12)+.5)

    boxplot(l,
main = main,
at = c(1,2,3,4,5,6,7),
names = c("A1","A2a","A2b","A3a","A3b","A4","NAA"),
las = 2,
horizontal = FALSE,
notch = TRUE,
col="blue",
boxwex=.15,
log="y")
  ## labl <- rep("", 10)
  ## labl[1] <- "1"
  ## labl[3] <- "100"
  ## labl[5] <- "1000"
  ## labl[8] <- "5000"           
  ## labl[10] <- "10000"
  ##          nums <- seq(1,10000,1000)
  ## axis(2, at = nums, label = labl, las = 2)


           }
}

#below makes tables and figures for everything! 



library(svglite)
qplot <- function(list,file,title)
{
    setwd("/home/daft/Dropbox/Steelhead_Rcode/Spawner_recruit_Long_data_full/Results")            
svglite(file=paste(file,"1.svg",sep=""))
plot_forcast(list,1,title)
dev.off()
svglite(file=paste(file,"2.svg",sep=""))
plot_forcast(list,2,title)
dev.off()
s <- table_forcast(list,title)

write.csv(file=paste(file,".csv",sep=""),s$table)
setwd("/home/daft/Dropbox/Steelhead_Rcode/Spawner_recruit_Long_data_full")   #these make special directories to store stuff might need to be changed on your computer     
return(s)
}

library("parallel")
require(KernSmooth)
require(Cairo)
require(tikzDevice)



egg <- function(xx,yy,col1="red",bw=350,title="m",log.vals=FALSE){
  #function for examining trade offs between performance metrics under uncertainty
                                        #original code provided by Aaron Greenberg on 6 May 2022


    xx <- na.omit(xx)
    yy <- na.omit(yy)
  if(log.vals){
    xx <- log(xx)
    yy <- log(yy)
  }
    print("max xx yy")
    print(c(max(xx,na.rm=TRUE), max(yy,na.rm=TRUE)))
      print(c(min(xx,na.rm=TRUE), min(yy,na.rm=TRUE)))      
  bwx <- abs((max(xx,na.rm=TRUE)-min(xx,na.rm=TRUE))/bw)
  bwy <- abs((max(yy,na.rm=TRUE)-min(yy,na.rm=TRUE))/bw)
  if(log.vals){
    est <- bkde2D(cbind(xx,yy),bandwidth=c(bwx,bwy),gridsize=c(100, 100)) #kernel density estimate
  }else{
                                        # est <- bkde2D(cbind(xx,yy),bandwidth=c(bwx,bwy),gridsize=c(81, 81),range.x=list(c(0,max(xx,na.rm=TRUE)+1.5*bwx),c(0,max(yy,na.rm=TRUE)+1.5*bwy))) #kernel density estimate
       est <- bkde2D(cbind(xx,yy),bandwidth=c(bwx,bwy),gridsize=c(300,300)) #kernel density estimate
  }
  est$fhat <- est$fhat/max(est$fhat,na.rm=TRUE)
  lvs=c(0.05,0.2,0.9) #what contour lines where?
  maxct <- max(lvs)
  nlvs <- length(lvs)
  thelines <- contourLines(est$x1,est$x2,est$fhat,levels=lvs)#three contours is enough
 
  if(log.vals){
    polygon(exp(thelines[[1]]$x),exp(thelines[[1]]$y),col=col1,border="red",lwd=1.1) #95% CI
    polygon(exp(thelines[[2]]$x),exp(thelines[[2]]$y),col=NA,border="blue",lwd=0.2) #80% CI
    polygon(exp(thelines[[3]]$x),exp(thelines[[3]]$y),col=NA,border="black",lwd=0.1) #10% CI
    text(median(exp(xx),na.rm=TRUE),median(exp(yy),na.rm=TRUE),labels=title)
  }else{
    polygon(thelines[[1]]$x,thelines[[1]]$y,col=col1,border="red",lwd=1.1) #95% CI
    polygon(thelines[[2]]$x,thelines[[2]]$y,col=NA,border="blue",lwd=0.2) #80% CI
    polygon(thelines[[3]]$x,thelines[[3]]$y,col=NA,border="black",lwd=0.1) #10% CI
    text(median(xx,na.rm=TRUE),median(yy,na.rm=TRUE),labels=title)
  }
  #points(median(xx,na.rm=TRUE),median(yy,na.rm=TRUE),pch=3)
}





#FA1 <- data_forcast(out,sdcon,z,40,1,marinescale=1,mortR=2,DPE=1,DPS=1,river=1,ssize=1)
#print(FA1$conmat)


isopleth_fn <- function(sims=sims,z=z,sdcon=sdcon,river=riv)
{

   
RSfunc <- function(marinescale,DPE=1,river,sims=sims,z=z,sdcon=sdcon)
{
    if(river<4){
        F <- data_forcast(out,sdcon,z,years=100,8,marinescale=marinescale,mortR=3,DPE=DPE,DPS=1,river,ssize=5000,sims=sims)
        }
    if(river==4){
        F <- data_forcast2(out,sdcon,z,years=100,8,marinescale=marinescale,mortR=3,DPE=DPE,DPS=1,ssize=5000,sims=sims)
        }
    return(mean(F$RSvec,na.rm=TRUE))
}


isopleth <- matrix(NA,nrow=30,ncol=30)
mortscale=(1:30)/15

for(i in 1:30)
{  
    isopleth[i,]=unlist(mclapply(mortscale,RSfunc,DPE=i/30,river,sims=sims,z=z,sdcon=sdcon,mc.cores=50))
}

return(isopleth)
}


isopleth_fn2 <- function(sims=sims,z=z,sdcon=sdcon,river=riv)
    {
RSfunc2 <- function(marinescale,DPE=1,river,sims=sims,z=z,sdcon=sdcon)
{
    if(river<4){
        F <- data_forcast(out,sdcon,z,years=100,8,marinescale=marinescale,mortR=3,DPE=DPE,DPS=1,river,ssize=5000,sims=sims)
        }
    if(river==4){
        F <- data_forcast2(out,sdcon,z,years=100,8,marinescale=marinescale,mortR=3,DPE=DPE,DPS=1,ssize=5000,sims=sims)
        }
    return(mean(F$Abunvec,na.rm=TRUE))
}


isopleth <- matrix(NA,nrow=30,ncol=30)
mortscale=(1:30)/15

for(i in 1:30)
{  
    isopleth[i,]=unlist(mclapply(mortscale,RSfunc2,DPE=i/30,river,sims=sims,z=z,sdcon=sdcon,mc.cores=50))
}

return(isopleth)
}




    
results <- function(numbersims=100000,morttype=3,years=100,sims=sims,out=out,sdcon=sdcon,z=z)
    {

        print("so what")
svglite(file="fit_plot.svg")
plot_magic(out,sdcon,z)
dev.off()

                                        #Foster
                     print("Foster")   
FA1 <- data_forcast(out,sdcon,z,years,1,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims,sims=sims)
FA2a <- data_forcast(out,sdcon,z,years,2,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims,sims=sims)
FA2b <- data_forcast(out,sdcon,z,years,3,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims,sims=sims)
FA3a <- data_forcast(out,sdcon,z,years,4,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims,sims=sims)
FA3b <- data_forcast(out,sdcon,z,years,5,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims,sims=sims)
FA4 <- data_forcast(out,sdcon,z,years,6,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims,sims=sims)
FNAA <- data_forcast(out,sdcon,z,years,7,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims,sims=sims)

f1 <- qplot(FA1,"FA1","Foster A1")
f2 <- qplot(FA2a,"FA2a","Foster A2a")
f3 <-qplot(FA2b,"FA2b","Foster A2b")
f4 <-qplot(FA3a,"FA3a","Foster A3a")
f5 <-qplot(FA3b,"FA3b","Foster A3b")
f6 <-qplot(FA4, "FA4","Foster A4")
f7 <-qplot(FNAA,"FNAA","Foster NAA")

        
        summaryt <- data.frame("FOS"=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold","SAR","Marine Survival","early Marine Survival","DPS","DPE","DPS*DPE"),
                       FA1=f1$table["Mean"],FA2a=f2$table["Mean"],FA2b=f3$table["Mean"],FA3a=f4$table["Mean"],FA3b=f5$table["Mean"],FA4=f6$table["Mean"],FNAA=f7$table["Mean"])
colnames(summaryt) <- c("FOS","A1","A2a","A2b","A3a","A3b","A4","NAA")
        write.csv(file="Foster.Summary.csv",summaryt)

        svglite(file="Foster_trade.svg")
        rgb.val <- col2rgb("grey27")
        alpha1=50
plot(NA,NA,xlim=c(0,5), ylim=c(-100,4000), ylab="Geo mean Spawners", xlab= "Recruits per Spawner") 
        egg(f1$RS,f1$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A1")
        egg(f2$RS,f2$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2A")
        egg(f3$RS,f3$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2B")
        egg(f4$RS,f4$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3a")
        egg(f5$RS,f5$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3b")
        egg(f6$RS,f6$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A4")
        egg(f7$RS,f7$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="NAA")
        dev.off()

        ## #Greenpeter
        svglite(file="Foster_summary_RS.svg")
        summaryplot(f1$RS,f2$RS,f3$RS,f4$RS,f5$RS,f6$RS,f7$RS,main="R/S")
        dev.off()
        svglite(file="Foster_summary_Abun.svg")
        summaryplot(f1$Abun,f2$Abun,f3$Abun,f4$Abun,f5$Abun,f6$Abund,f7$Abund,main="log Geo Mean Abundance",flag=3)
        dev.off()
        svglite(file="Foster_summary_QET.svg")
        summaryplot(f1$QET,f2$QET,f3$QET,f4$QET,f5$QET,f6$QET,f7$QET,main="QET thresh",flag=2)
        dev.off()

        
             print("GreenPeter")   
GPA1 <- data_forcast(out,sdcon,z,years,1,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims,sims=sims)
GPA2a <- data_forcast(out,sdcon,z,years,2,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims,sims=sims)
GPA2b <- data_forcast(out,sdcon,z,years,3,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims,sims=sims)
GPA3a <- data_forcast(out,sdcon,z,years,4,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims,sims=sims)
GPA3b <- data_forcast(out,sdcon,z,years,5,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims,sims=sims)

gp1 <- qplot(GPA1,"GPA1","Green Peter A1")
gp2 <- qplot(GPA2a,"GPA2a","Green Peter A2a")
gp3 <- qplot(GPA2b,"GPA2b","Green Peter A2b")
gp4 <- qplot(GPA3a,"GPA3a","Green Peter A3a")
gp5 <- qplot(GPA3b,"GPA3b","Green Peter A3b")
summaryt <- data.frame("GRP"=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold","SAR","Marine Survival","early Marine Survival","DPS","DPE","DPS*DPE"),
                      GPA1=gp1$table["Mean"],GPA2a=gp2$table["Mean"],GPA2b=gp3$table["Mean"],GPA3a=gp4$table["Mean"],GPA3b=gp5$table["Mean"])
colnames(summaryt) <- c("GP","A1","A2a","A2b","A3a","A3b")
write.csv(file="Green.Peter.Summary.csv",summaryt)

        svglite(file="Green_trade.svg")
       rgb.val <- col2rgb("grey27")
        alpha1=50
plot(NA,NA,xlim=c(0,5), ylim=c(-100,3000), ylab="Geo mean Spawners", xlab= "Recruits per Spawner") 
        egg(gp1$RS,gp1$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A1")
        egg(gp2$RS,gp2$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2A")
        egg(gp3$RS,gp3$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2B")
        egg(gp4$RS,gp4$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3a")
        egg(gp5$RS,gp5$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3b")
        dev.off()

        svglite(file="GreenPeter_summary_RS.svg")
        summaryplot(gp1$RS,gp2$RS,gp3$RS,gp4$RS,gp5$RS,gp5$RS*NA,gp5$RS*NA,main="R/S")
        dev.off()
        svglite(file="GreenPeter_summary_Abun.svg")
        summaryplot(gp1$Abun,gp2$Abun,gp3$Abun,gp4$Abun,gp5$Abun,gp5$Abund*NA,gp5$Abund*NA,main="log Geo Mean Abundance",flag=3)
        dev.off()
        svglite(file="GreenPeter_summary_QET.svg")
        summaryplot(gp1$QET,gp2$QET,gp3$QET,gp4$QET,gp5$QET,gp5$QET*NA,gp5$QET*NA,main="QET thresh",flag=2)
        dev.off()

        
     print("Lets make the trade off plots!")   

        
#
     print("Detroit")   


        
#Greenpeter
DETA1 <- data_forcast(out,sdcon,z,years,1,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims,sims=sims)
DETA2a <- data_forcast(out,sdcon,z,years,2,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims,sims=sims)
DETA2b <- data_forcast(out,sdcon,z,years,3,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims,sims=sims)
DETA3a <- data_forcast(out,sdcon,z,years,4,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims,sims=sims)
DETA3b <- data_forcast(out,sdcon,z,years,5,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims,sims=sims)
DETA4 <- data_forcast(out,sdcon,z,years,6,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims,sims=sims)
DETNAA <- data_forcast(out,sdcon,z,years,7,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims,sims=sims)


d1 <- qplot(DETA1,"DETA1","Detroit A1")
d2 <- qplot(DETA2a,"DETA2a","Detroit A2a")
d3 <-qplot(DETA2b,"DETA2b","Detroit A2b")
d4 <-qplot(DETA3a,"DETA3a","Detroit A3a")
d5 <-qplot(DETA3b,"DETA3b","Detroit A3b")
d6 <-qplot(DETA4, "DETA4","Detroit A4")
d7 <-qplot(DETNAA,"DNAA","Detroit NAA")

summaryt <- data.frame("DET"=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold","SAR","Marine Survival","early Marine Survival","DPS","DPE","DPS*DPE"),
                       DA1=d1$table["Mean"],DA2a=d2$table["Mean"],DA2b=d3$table["Mean"],DA3a=d4$table["Mean"],DA3b=d5$table["Mean"],DA4=d6$table["Mean"],DNAA=d7$table["Mean"])
colnames(summaryt) <- c("DET","A1","A2a","A2b","A3a","A3b","A4","NAA")
write.csv(file="Detroit.Summary.csv",summaryt)

     svglite(file="Det_trade.svg")
       rgb.val <- col2rgb("grey27")
        alpha1=50
plot(NA,NA,xlim=c(0,5), ylim=c(-100,4000), ylab="Geo mean Spawners", xlab= "Recruits per Spawner") 
        egg(d1$RS,d1$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A1")
        egg(d2$RS,d2$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2A")
        egg(d3$RS,d3$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2B")
        egg(d4$RS,d4$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3a")
        egg(d5$RS,d5$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3b")
        egg(d6$RS,d6$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A4")
        egg(d7$RS,d7$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="NAA")
        dev.off()

        svglite(file="Det_summary_RS.svg")
        summaryplot(d1$RS,d2$RS,d3$RS,d4$RS,d5$RS,d6$RS,d7$RS,main="R/S")
        dev.off()
        svglite(file="Det_summary_Abun.svg")
        summaryplot(d1$Abun,d2$Abun,d3$Abun,d4$Abun,d5$Abun,d6$Abund,d7$Abund,main="log Geo Mean Abundance",flag=3)
        dev.off()
        svglite(file="Det_summary_QET.svg")
        summaryplot(d1$QET,d2$QET,d3$QET,d4$QET,d5$QET,d6$QET,d7$QET,main="QET thresh",flag=2)
        dev.off()

        
       print("South Sant")   
      
## #Greenpeter
SSA1 <- data_forcast2(out,sdcon,z,years,1,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims,sims=sims)
SSA2a <- data_forcast2(out,sdcon,z,years,2,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims,sims=sims)
SSA2b <- data_forcast2(out,sdcon,z,years,3,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims,sims=sims)
SSA3a <- data_forcast2(out,sdcon,z,years,4,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims,sims=sims)
SSA3b <- data_forcast2(out,sdcon,z,years,5,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims,sims=sims)
SSA4 <- data_forcast2(out,sdcon,z,years,6,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims,sims=sims)
SSNAA <- data_forcast2(out,sdcon,z,years,7,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims,sims=sims)

ss1 <-qplot(SSA1,"SSA1","South Sant. A1")
ss2 <-qplot(SSA2a,"SSA2a","South Sant.A2a")
ss3 <-qplot(SSA2b,"SSA2b","South Sant. A2b")
ss4 <-qplot(SSA3a,"SSA3a","South Sant. A3a")
ss5 <-qplot(SSA3b,"SSA3b","South Sant. A3b")
ss6 <-qplot(SSA4, "SSA4","South Sant. A4")
ss7 <-qplot(SSNAA,"SSNAA","South Sant. NAA")
summaryt <- data.frame("SS"=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold","SAR","Marine Survival","early Marine Surival","DPS","DPE","DPS*DPE"),
                       SSA1=ss1$table["Mean"],SSA2a=ss2$table["Mean"],SSA2b=ss3$table["Mean"],SSA3a=ss4$table["Mean"],SSA3b=ss5$table["Mean"],SSA4=ss6$table["Mean"],SSNAA=ss7$table["Mean"])
colnames(summaryt) <- c("SS","A1","A2a","A2b","A3a","A3b","A4","NAA")
        write.csv(file="South.Sant.Summary.csv",summaryt)


     svglite(file="South_trade.svg")
       rgb.val <- col2rgb("grey27")
        alpha1=50
plot(NA,NA,xlim=c(0,5), ylim=c(-100,4000), ylab="Geo mean Spawners", xlab= "Recruits per Spawner") 
        egg(ss1$RS,ss1$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A1")
        egg(ss2$RS,ss2$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2A")
        egg(ss3$RS,ss3$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2B")
        egg(ss4$RS,ss4$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3a")
        egg(ss5$RS,ss5$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3b")
        egg(ss6$RS,ss6$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A4")
        egg(ss7$RS,ss7$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="NAA")
        dev.off()
#  

        svglite(file="South_summary_RS.svg")
        summaryplot(ss1$RS,ss2$RS,ss3$RS,ss4$RS,ss5$RS,ss6$RS,ss7$RS,main="R/S")
        dev.off()
        svglite(file="South_summary_Abun.svg")
        summaryplot(ss1$Abun,ss2$Abun,ss3$Abun,ss4$Abun,ss5$Abun,ss6$Abund,ss7$Abund,main="Geo Mean Abundance",flag=3)
        dev.off()
        svglite(file="South_summary_QET.svg")
        summaryplot(ss1$QET,ss2$QET,ss3$QET,ss4$QET,ss5$QET,ss6$QET,ss7$QET,main="QET thresh",flag=2)
        dev.off()

        
    }


results(numbersims=100000,morttype=3,years=100,sims=sims,out=out,sdcon=sdcon,z=z)

iso_plot <- function(file,iso,river=1){
    if(river==1)
        {
Alt1 <- read.csv("./Summary_Steelhead/Alt1_FOS_age2.csv",header=TRUE)
Alt1 <- mean(Alt1$x.DPS*Alt1$x.DPE)
Alt2a <- read.csv("./Summary_Steelhead/Alt2a_FOS_age2.csv",header=TRUE)
Alt2a <- mean(Alt2a$x.DPS*Alt2a$x.DPE)
Alt2b <- read.csv("./Summary_Steelhead/Alt2b_FOS_age2.csv",header=TRUE)
Alt2b <- mean(Alt2b$x.DPS*Alt2b$x.DPE)
Alt3a<- read.csv("./Summary_Steelhead/Alt3a_FOS_age2.csv",header=TRUE)
Alt3a <- mean(Alt3a$x.DPS*Alt3a$x.DPE)
Alt3b<- read.csv("./Summary_Steelhead/Alt3b_FOS_age2.csv",header=TRUE)
Alt3b <- mean(Alt3b$x.DPS*Alt3b$x.DPE)
Alt4<- read.csv("Summary_Steelhead/Alt4_FOS_age2.csv",header=TRUE)
Alt4 <- mean(Alt4$x.DPS*Alt4$x.DPE)
NAA <- read.csv("./Summary_Steelhead/NAA_FOS_age2.csv",header=TRUE)
NAA <- mean(NAA$x.DPS*NAA$x.DPE)


svglite(file)
x=(1:30)/15
y=(1:30)/30
par(mar = c(5, 5, 5, 5)) 
image(x,y,iso,ylab="DPS*DPE",xlab="Marine Mortality Scaling")
contour(x,y,iso,add=TRUE)

abline(h=Alt1)
axis(4, at = (Alt1+.02), labels = c("Alt1"),las=2,tck=0)
axis(4, at = (Alt1), labels = NA,las=2)
abline(h=Alt2a)
axis(2, at = (Alt2a), labels = c("Alt2a"),las=2)
abline(h=Alt2b)
axis(4, at = (Alt2b-.02), labels = c("Alt2b"),las=2,tck=0)
abline(h=Alt3a)
axis(2, at = Alt3a, labels = c("Alt3a"),las=2)
abline(h=Alt3b)
axis(4, at = (Alt3b+.02), labels = c("Alt3b"),las=2,tck=0)
axis(4, at = (Alt3b), labels = NA,las=2)
abline(h=Alt4)
axis(2, at = Alt4, labels = c("Alt4"),las=2)
abline(h=NAA)
axis(4, at = (NAA-.02), labels = c("NAA"),las=2,tck=0)
dev.off()
}
    if(river==2)
        {
Alt1 <- read.csv("./Summary_Steelhead/Alt1_GPR_age2.csv",header=TRUE)
Alt1 <- mean(Alt1$x.DPS*Alt1$x.DPE)
Alt2a <- read.csv("./Summary_Steelhead/Alt2a_GPR_age2.csv",header=TRUE)
Alt2a <- mean(Alt2a$x.DPS*Alt2a$x.DPE)
Alt2b <- read.csv("./Summary_Steelhead/Alt2b_GPR_age2.csv",header=TRUE)
Alt2b <- mean(Alt2b$x.DPS*Alt2b$x.DPE)
Alt3a<- read.csv("./Summary_Steelhead/Alt3a_GPR_age2.csv",header=TRUE)
Alt3a <- mean(Alt3a$x.DPS*Alt3a$x.DPE)
Alt3b<- read.csv("./Summary_Steelhead/Alt3b_GPR_age2.csv",header=TRUE)
Alt3b <- mean(Alt3b$x.DPS*Alt3b$x.DPE)
Alt4<- read.csv("Summary_Steelhead/Alt4_GPR_age2.csv",header=TRUE)
Alt4 <- mean(Alt4$x.DPS*Alt4$x.DPE)
NAA <- read.csv("./Summary_Steelhead/NAA_GPR_age2.csv",header=TRUE)
NAA <- mean(NAA$x.DPS*NAA$x.DPE)


svglite(file)
x=(1:30)/15
y=(1:30)/30
par(mar = c(5, 5, 5, 5)) 
image(x,y,iso,ylab="DPS*DPE",xlab="Marine Mortality Scaling")
contour(x,y,iso,add=TRUE)

abline(h=Alt1)
axis(4, at = (Alt1+.02), labels = c("Alt1"),las=2,tck=0)
axis(4, at = (Alt1), labels = NA,las=2)
abline(h=Alt2a)
axis(2, at = (Alt2a), labels = c("Alt2a"),las=2)
abline(h=Alt2b)
axis(4, at = (Alt2b-.02), labels = c("Alt2b"),las=2,tck=0)
abline(h=Alt3a)
axis(2, at = Alt3a, labels = c("Alt3a"),las=2)
abline(h=Alt3b)
axis(4, at = (Alt3b+.02), labels = c("Alt3b"),las=2,tck=0)
axis(4, at = (Alt3b), labels = NA,las=2)
dev.off()
}
    if(river==3)
        {
Alt1 <- read.csv("./Summary_Steelhead/Alt1_DET_age2.csv",header=TRUE)
Alt1 <- mean(Alt1$x.DPS*Alt1$x.DPE)
Alt2a <- read.csv("./Summary_Steelhead/Alt2a_DET_age2.csv",header=TRUE)
Alt2a <- mean(Alt2a$x.DPS*Alt2a$x.DPE)
Alt2b <- read.csv("./Summary_Steelhead/Alt2b_DET_age2.csv",header=TRUE)
Alt2b <- mean(Alt2b$x.DPS*Alt2b$x.DPE)
Alt3a<- read.csv("./Summary_Steelhead/Alt3a_DET_age2.csv",header=TRUE)
Alt3a <- mean(Alt3a$x.DPS*Alt3a$x.DPE)
Alt3b<- read.csv("./Summary_Steelhead/Alt3b_DET_age2.csv",header=TRUE)
Alt3b <- mean(Alt3b$x.DPS*Alt3b$x.DPE)
Alt4<- read.csv("Summary_Steelhead/Alt4_DET_age2.csv",header=TRUE)
Alt4 <- mean(Alt4$x.DPS*Alt4$x.DPE)
NAA <- read.csv("./Summary_Steelhead/NAA_DET_age2.csv",header=TRUE)
NAA <- mean(NAA$x.DPS*NAA$x.DPE)


svglite(file)
x=(1:30)/15
y=(1:30)/30
par(mar = c(5, 5, 5, 5)) 
image(x,y,iso,ylab="DPS*DPE",xlab="Marine Mortality Scaling")
contour(x,y,iso,add=TRUE)

abline(h=Alt1)
axis(4, at = (Alt1+.02), labels = c("Alt1"),las=2,tck=0)
axis(4, at = (Alt1), labels = NA,las=2)
abline(h=Alt2a)
axis(2, at = (Alt2a), labels = c("Alt2a"),las=2)
abline(h=Alt2b)
axis(4, at = (Alt2b-.02), labels = c("Alt2b"),las=2,tck=0)
abline(h=Alt3a)
axis(2, at = Alt3a, labels = c("Alt3a"),las=2)
abline(h=Alt3b)
axis(4, at = (Alt3b+.02), labels = c("Alt3b"),las=2,tck=0)
axis(4, at = (Alt3b), labels = NA,las=2)
abline(h=Alt4)
axis(2, at = Alt4, labels = c("Alt4"),las=2)
abline(h=NAA)
axis(4, at = (NAA-.02), labels = c("NAA"),las=2,tck=0)
dev.off()
dev.off()
if(river==4)
{
svglite(file)
x=(1:30)/15
y=(1:30)/30
par(mar = c(5, 5, 5, 5)) 
image(x,y,iso,ylab="DPS*DPE",xlab="Marine Mortality Scaling")
contour(x,y,iso,add=TRUE)
dev.off()
    }

}


}


isopleth <- isopleth_fn(sims=sims,z=z,sdcon=sdcon,river=2)
iso_plot("isopleth_test_Greenpeter.svg",isopleth,river=2)
isopleth <- isopleth_fn(sims=sims,z=z,sdcon=sdcon,river=1)
iso_plot("isopleth_test_Foster.svg",isopleth,river=1)
isopleth <- isopleth_fn(sims=sims,z=z,sdcon=sdcon,river=3)
iso_plot("isopleth_test_Detroit.svg",isopleth,river=3)
isopleth <- isopleth_fn(sims=sims,z=z,sdcon=sdcon,river=4)
iso_plot("isopleth_test_SouthSant.svg",isopleth,river=4)

isopleth <- isopleth_fn2(sims=sims,z=z,sdcon=sdcon,river=1)
iso_plot("isopleth_test_Foster2.svg",isopleth,river=1)
isopleth <- isopleth_fn2(sims=sims,z=z,sdcon=sdcon,river=3)
iso_plot("isopleth_test_Detroit2.svg",isopleth,river=3)
isopleth <- isopleth_fn2(sims=sims,z=z,sdcon=sdcon,river=4)
iso_plot("isopleth_test_SouthSant2.svg",isopleth,river=4)
isopleth <- isopleth_fn2(sims=sims,z=z,sdcon=sdcon,river=2)
iso_plot("isopleth_test_Greenpeter2.svg",isopleth,river=2)



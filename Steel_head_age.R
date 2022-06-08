#!/usr/local/bin/Rscript
graphics.off()
rm(list=ls())
source("Plots_tables.R")
library("matrixStats")
data <- read.csv("data.csv", sep=",",head=TRUE)
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
years <- 31
                                        #we want to estimate SF for a standard DPSDPE

model <- function(SF,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
{
ages <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)
data$predSmoltAParameter<- 2078+data$SpawnerYr*0
#print("oopsa")
data$predfemalesmolt_pre_dam[1:7] <-SF* data$predSmoltAParameter[1:7]*data$FosterCount[1:7]/(1+data$predSmoltAParameter[1:7]*SF*data$FosterCount[1:7]/100000.0)
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*data$passprob[1:7]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*data$smoltsurvival[1:7]*.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
#print("oopsnana")
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
#print("oops1")
for(i in 8:years)
{
    ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
    data$predfemalesmolt_pre_dam[i]=data$predSmoltAParameter[i]*(SF*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*(SF*ages$SpawnerAdult[i]))/100000.0)
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*data$passprob[i]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*data$smoltsurvival[i]*.74
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
}




ages$Total=ages$vir4 +ages$vir5 +ages$vir6+ages$rep5 +ages$rep6+ages$rep7
eps <- data$FosterCount[17:years]-ages$Total[13:(years-4)]
return(list(data=data,ages=ages,eps=eps))
}



fn <- function(par,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data){
    eps <- model(par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
    nll <--1*sum(log(dnorm(eps$ages$SpawnerAdult[15:(years-1)],eps$data$FosterCount[15:(years-1)],par[2])))
    #nll <--1*sum(log(dnorm(eps$ages$Total[15:(years-5)],eps$data$FosterCount[20:(years-1)],par[2])))
return(nll)
}




modelinit <- function(SF,init,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
{
ages <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)
data$predSmoltAParameter<- 2078+data$SpawnerYr*0
#print("oopsa")
data$predfemalesmolt_pre_dam[1:7] <-SF* data$predSmoltAParameter[1:7]*(init+0*data$FosterCount[1:7])/(1+data$predSmoltAParameter[1:7]*SF*(init+0*data$FosterCount[1:7])/100000.0)
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*data$passprob[1:7]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*data$smoltsurvival[1:7]*.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
#print("oopsnana")
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
ages$SpawnerAdult[1:7]=init
#print("oops1")
for(i in 8:years)
{
    ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
    data$predfemalesmolt_pre_dam[i]=data$predSmoltAParameter[i]*(SF*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*(SF*ages$SpawnerAdult[i]))/100000.0)
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*data$passprob[i]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*data$smoltsurvival[i]*.74
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
}




ages$Total=ages$vir4 +ages$vir5 +ages$vir6+ages$rep5 +ages$rep6+ages$rep7
eps <- data$FosterCount[17:years]-ages$Total[13:(years-4)]
return(list(data=data,ages=ages,eps=eps))
}




modelconst <- function(SF,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
{
ages <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)
data$predSmoltAParameter<- 2078+data$SpawnerYr*0
#print("oopsa")
data$predfemalesmolt_pre_dam[1:7] <-SF* data$predSmoltAParameter[1:7]*200
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*data$passprob[1:7]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*data$smoltsurvival[1:7]*.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
#print("oopsnana")
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
#print("oops1")
for(i in 8:years)
{
    ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
    data$predfemalesmolt_pre_dam[i]=SF*data$predSmoltAParameter[i]*200
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*data$passprob[i]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*data$smoltsurvival[i]*.74
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
}




ages$Total=ages$vir4 +ages$vir5 +ages$vir6+ages$rep5 +ages$rep6+ages$rep7
eps <- data$FosterCount[17:years]-ages$Total[13:(years-4)]
return(list(data=data,ages=ages,eps=eps))
}


modellin <- function(SF,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
{
ages <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)
data$predSmoltAParameter<- 2078+data$SpawnerYr*0
#print("oopsa")
data$predfemalesmolt_pre_dam[1:7] <-SF* data$predSmoltAParameter[1:7]*data$FosterCount[1:7]
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*data$passprob[1:7]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*data$smoltsurvival[1:7]*.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
#print("oopsnana")
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
#print("oops1")
for(i in 8:years)
{
    ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
    data$predfemalesmolt_pre_dam[i]=data$predSmoltAParameter[i]*(SF*ages$SpawnerAdult[i])
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*data$passprob[i]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*data$smoltsurvival[i]*.74
  ages$vir4[i]=Percent4*data$OceanAdults[i]
  ages$vir5[i]=Percent5*data$OceanAdults[i]
  ages$vir6[i]=Percent6*data$OceanAdults[i]
  ages$rep5[i]=Repeat6*ages$vir4[i-1]
  ages$rep6[i]=Repeat6*ages$vir5[i-1]
  ages$rep7[i]=Repeat6*ages$vir6[i-1]
}




ages$Total=ages$vir4 +ages$vir5 +ages$vir6+ages$rep5 +ages$rep6+ages$rep7
eps <- data$FosterCount[17:years]-ages$Total[13:(years-4)]
return(list(data=data,ages=ages,eps=eps))
}








test <- model(.0275,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
#write.csv(test$data,file="test.csv")
print(test$data)

## inputs <- c(.053,100)
## z=optim(par=inputs,fn=fn,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
## sdcon <- sqrt(diag(solve(z$hessian)))
## print(z)

## test2 <- model(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
## print(test2$data)

## plot(test2$data$FosterCount)
## abline(v=7,col="red")
## lines(test2$ages$SpawnerAdult)
fnlin <- function(par,DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data){
    eps <- modellin(par[1],DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data)
    nll <--1*sum(log(dnorm(eps$ages$SpawnerAdult[15:(years-1)],eps$data$FosterCount[15:(years-1)],par[2])))
    
return(nll)
}



fnconst <- function(par,DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data){
    eps <- modelconst(par[1],DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data)
    nll <--1*sum(log(dnorm(eps$ages$SpawnerAdult[15:(years-1)],eps$data$FosterCount[15:(years-1)],par[2])))
return(nll)
}


fninit <- function(par,DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data){
    eps <- modelinit(par[1],par[2],DPSDPE,years,Percent4,Percent5,Repeat5,Repeat6,data)
    eps2 <-(eps$data$FosterCount[16:(years-1)])-(eps$ages$SpawnerAdult[16:(years-1)])
    nll <--1*sum(log(dnorm((eps2),0,par[3])))-log(dnorm(par[2],(242),(120)))-log(dnorm(par[1],.03,.01))
return(nll)
}

inputs <- c(.03,200,30)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
print(z)


simmaker <- function(nsims=100000,years=50,z,sdcon)
{
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
        srand[i,1:years] <- sample(data$smoltsurvival[7:18],length(1:years),replace=TRUE) #probably needs to be chainged
        par <- c(0.01310638, 0.64901850, 0.01173065, 0.01698823)
        sar[i,1:6] <- par[1]
        s <- s*0
        s[1:6] <- par[1]
    for(j in 7:(years))
    {
        s[j] <- s[j-1]*par[2]+par[3]+rnorm(1,0,par[4])

    }
    s[which(s<=0)] <- 0              
    sar[i,]<- s

        }
    
     cov=as.data.frame(solve(z$hessian))
    tab= MASS::mvrnorm(nsims,z$par,cov)
    fresh <- tab[,1]
    init<-  tab[,2]                  
    return(list(ssin=ssin,srand=srand,sar=sar,fresh=fresh,init=init))
    }

sims=simmaker(nsims=10000,years=50,z,sdcon)





forcastmodel <- function(EIS,mortR,marinescale,years,Percent4,Percent5,Repeat5,Repeat6,data,DPE=1,DPS=1,river=1,simindex)
{
    ages <- list(SpawnerYR=1:years,vir4=1:years*0,vir5=1:years*0,repeat5=1:years*0, repeat6 =1:years*0,Total=1:years*0)

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

                                        #need to select sims
  
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

if(EIS==8)
    {
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
        FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_GPR_age2.csv",header=TRUE)
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
        FOS_Age2 <- read.csv("./Summary_Steelhead/NAA_GPR_age2.csv",header=TRUE)
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




    
    index <- sample(1:length(FOS_Age2$x.Year),years,replace=TRUE)
    data <- as.list(data)
    length(data$predfemalesmolt_pre_dam) <- years
    data$predSmoltAParameter <- 1:years*0

    data$predSmoltAParameter<- 1:years*0+ 2078
    
ages <- list(SpawnerYR=data$SpawnerYr,vir4=data$SpawnerYr*0,vir5=data$SpawnerYr*0,vir6=data$SpawnerYr*0,rep5=data$SpawnerYr*0, rep6 =data$SpawnerYr*0,rep7 =data$SpawnerYr*0,SpawnerAdult=data$SpawnerYr*0,Total=data$SpawnerYr*0)
data$predSmoltAParameter<- 2078+data$SpawnerYr*0
data$predfemalesmolt_pre_dam[1:7] <- data$predSmoltAParameter[1:7]*(sims$init[simindex]+0*data$FosterCount[1:7])*sims$fresh[simindex]*freshscale/(1+data$predSmoltAParameter[1:7]*(sims$init[simindex]+0*data$FosterCount[1:7])*sims$fresh[simindex]*freshscale/capacityfac)
data$predfemalesmolt_post_dam[1:7] <- data$predfemalesmolt_pre_dam[1:7]*FOS_Age2$x.DPS[index[1:7]]*FOS_Age2$x.DPE[index[1:7]]
data$OceanAdults[1:7] <- data$predfemalesmolt_post_dam[1:7]*simforcast[1:7]*marinescale *.74
ages$vir4[1:7] <- Percent4*data$OceanAdults[1:7]
ages$vir5[1:7] <- Percent5*data$OceanAdults[1:7]
ages$vir6[1:7] <- Percent6*data$OceanAdults[1:7]
ages$rep5[1:7]=Repeat6*ages$vir4[1:7]
ages$rep6[1:7]=Repeat6*ages$vir5[1:7]
ages$rep7[1:7]=Repeat6*ages$vir5[1:7]
for(i in 8:years)
{

  ages$SpawnerAdult[i]=ages$vir4[i-4]+ages$vir5[i-5]+ages$vir6[i-6]+ages$rep5[i-5]+ages$rep6[i-6]+ages$rep7[i-7]
  data$predfemalesmolt_pre_dam[i]=(data$predSmoltAParameter[i]*sims$fresh[simindex]*freshscale*ages$SpawnerAdult[i])/(1+(data$predSmoltAParameter[i]*ages$SpawnerAdult[i]*sims$fresh[simindex]*freshscale/capacityfac))
  data$predfemalesmolt_post_dam[i]=data$predfemalesmolt_pre_dam[i]*FOS_Age2$x.DPS[index[i]]*FOS_Age2$x.DPE[index[i]]
  data$OceanAdults[i]=data$predfemalesmolt_post_dam[i]*simforcast[i]*marinescale *.74
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
        
RS <- exp(mean(log(data$OceanAdults[15:(15+5)]/ages$SpawnerAdult[(15):(15+5)])))
abun <- exp(mean(log(ages$Total[21:(30)]),na.rm=TRUE))
return(list(data=data,ages=ages,eps=eps,sim=simforcast[1:years],RS=RS,abund=abun,QET=QETave,DPS=FOS_Age2$x.DPS[index],DPE= FOS_Age2$x.DPE[index],confi=sims$fresh[simindex]*freshscale))

}













data_forcast <- function(out,sdcon,z,years,EIS,marinescale=1,mortR=2,DPE=1,DPS=1,river=1,ssize=1000)
{
conmat <-matrix(NA,nrow=ssize,ncol=years)
QETamat <-matrix(NA,nrow=ssize,ncol=years)        
simmat <-matrix(NA,nrow=ssize,ncol=years)
        RSmat <-matrix(NA,nrow=ssize,ncol=years)
        QETmat <-matrix(NA,nrow=ssize,ncol=years)
        DPSmat <-matrix(NA,nrow=ssize,ncol=years)
        DPEmat <-matrix(NA,nrow=ssize,ncol=years)
        RSvec <- 1:ssize
        Abunvec <- 1:ssize
        marineSvec <- 1:ssize
        QETvec <- 1:ssize
confi2 <- 1:ssize
for( i in 1:ssize)
{
    out <- forcastmodel(EIS,mortR,marinescale,years,Percent4,Percent5,Percent6,Repeat6,data,DPE,DPS,river,i)
    conmat[i,1:(years-1)] <-out$ages$SpawnerAdult[1:(years-1)]
    simmat[i,1:(years-1)] <- out$sim[1:(years-1)]
    QETamat[i,1:(years-1)] <- out$QET[1:(years-1)]
    DPEmat[i,1:(years-1)] <- out$DPE[1:(years-1)]
    DPSmat[i,1:(years-1)] <- out$DPS[1:(years-1)]
    RSmat[i,(15):(years-1)] <- out$data$OceanAdults[(15):(years-1)]/out$ages$SpawnerAdult[(15):(years-1)]
    QETmat[i,] <- ifelse(out$QET[1:(years)]>100,0,1)
    RSvec[i]=out$RS
    Abunvec[i]=out$abund
    marineSvec[i] <- mean(out$sim[15:(15+5)])
    QETvec[i] <- ifelse(min(out$QET[15:(years-1)],na.rm=TRUE)>100,0,1)
    confi2[i] <- out$confi
 
}
return(list(conmat=conmat,simmat=simmat,QETamat=QETamat,DPEmat=DPEmat,DPSmat=DPSmat,RSmat=RSmat,QETmat=QETmat,RSvec=RSvec,Abunvec=Abunvec,marineSvec=marineSvec,QETvec=QETvec,confi=confi2,DPEvec=as.vector(DPEmat),DPSvec=as.vector(DPSmat)))
}


data_forcast2 <- function(out,sdcon,z,years,EIS,marinescale=1,mortR=2,DPE=1,DPS=1,ssize=1000)
{
confi <- rnorm(ssize,z$par[1],sd=sdcon[1])            
conmat <-matrix(NA,nrow=ssize,ncol=years)
QETamat <-matrix(NA,nrow=ssize,ncol=years)        
simmat <-matrix(NA,nrow=ssize,ncol=years)
        RSmat <-matrix(NA,nrow=ssize,ncol=years)
        QETmat <-matrix(NA,nrow=ssize,ncol=years)
        DPSmat <-matrix(NA,nrow=ssize,ncol=years)
        DPEmat <-matrix(NA,nrow=ssize,ncol=years)
        RSvec <- 1:ssize
        Abunvec <- 1:ssize
        marineSvec <- 1:ssize
        QETvec <- 1:ssize
    confi2 <- 1:ssize
for( i in 1:ssize)
{
    out <- forcastmodel(EIS,mortR,marinescale,years,Percent4,Percent5,Percent6,Repeat6,data,DPE,DPS,1,i)#foster
    if(EIS<=5){
    out2 <- forcastmodel(EIS,mortR,marinescale,years,Percent4,Percent5,Percent6,Repeat6,data,DPE,DPS,2,i)#greenpeter
    conmat[i,1:(years-1)] <-out$ages$SpawnerAdult[1:(years-1)]+out2$ages$SpawnerAdult[1:(years-1)]
    simmat[i,1:(years-1)] <- (out$sim[1:(years-1)]+out2$sim[1:(years-1)])/2 #average marine survival
    QETamat[i,1:(years-1)] <-(out$QET[1:(years-1)]+out2$QET[1:(years-1)])
    DPEmat[i,1:(years-1)] <- (out$DPE[1:(years-1)]+out2$DPE[1:(years-1)])/2
    DPSmat[i,1:(years-1)] <- (out$DPS[1:(years-1)]+out2$DPS[1:(years-1)])/2
    RSmat[i,(15):(years-1)] <- (out$data$OceanAdults[(15):(years-1)]+out2$data$OceanAdults[(15):(years-1)])/(out$ages$SpawnerAdult[(15-4):(years-5)]+out2$ages$SpawnerAdult[(15-4):(years-5)])
    QETmat[i,] <- ifelse(out$QET[1:years]+out$QET[1:years]>100,0,1)
    RSvec[i]=exp(mean(log((out$data$OceanAdults[(15):(15+5)]+out2$data$OceanAdults[(15):(15+5)])/(out$ages$SpawnerAdult[(15):(15+5)]+out2$ages$SpawnerAdult[(15):(15+5)]))))
    Abunvec[i]=exp(mean(log((out$ages$SpawnerAdult[(21):(30)]+out2$ages$SpawnerAdult[21:30]))))
    marineSvec[i] <- mean((out$sim[15:(15+5)]+out2$sim[15:(15+5)])/2)
    QETvec[i] <- ifelse(min(out$QET[15:(years-1)]+out2$QET[15:(years-1)],na.rm=TRUE)>100,0,1)
       confi2[i] <-(out$confi+out2$confi)/2 
    }else
        {
    
    conmat[i,1:(years-1)] <-out$ages$SpawnerAdult[1:(years-1)]
    simmat[i,1:(years-1)] <- out$sim[1:(years-1)]
    QETamat[i,1:(years-1)] <- out$QET[1:(years-1)]
    DPEmat[i,1:(years-1)] <- out$DPE[1:(years-1)]
    DPSmat[i,1:(years-1)] <- out$DPS[1:(years-1)]
    RSmat[i,(15):(years-1)] <- out$data$OceanAdults[(15):(years-1)]/out$ages$SpawnerAdult[(15):(years-1)]
    QETmat[i,] <- ifelse(out$QET[1:years]>100,0,1)
    RSvec[i]=out$RS
    Abunvec[i]=out$abund
    marineSvec[i] <- mean(out$sim[15:(15+5)])
    QETvec[i] <- ifelse(min(out$QET[7:(years-1)],na.rm=TRUE)>100,0,1)
            confi2[i] <-out$confi 
   
}
}
return(list(conmat=conmat,simmat=simmat,QETamat=QETamat,DPEmat=DPEmat,DPSmat=DPSmat,RSmat=RSmat,QETmat=QETmat,RSvec=RSvec,Abunvec=Abunvec,marineSvec=marineSvec,QETvec=QETvec,confi=confi2,DPEvec=as.vector(DPEmat),DPSvec=as.vector(DPSmat)))
}



library(svglite)
qplot <- function(list,file,title)
{
setwd("/home/daft/Dropbox/Steelhead_Rcode/Spawner_recruit/Results")            
svglite(file=paste(file,"1.svg",sep=""))
plot_forcast(list,1,title)
dev.off()
svglite(file=paste(file,"2.svg",sep=""))
plot_forcast(list,2,title)
dev.off()
s <- table_forcast(list,title)

write.csv(file=paste(file,".csv",sep=""),s$table)
setwd("/home/daft/Dropbox/Steelhead_Rcode/Spawner_recruit")        
return(s)
}


require(KernSmooth)
require(Cairo)
require(tikzDevice)

egg<-function(xx,yy,col1="red",bw=20,title="m")
{
    bwx=abs((max(xx,na.rm=TRUE)-min(xx,na.rm=TRUE))/bw)
    bwy=abs((max(yy,na.rm=TRUE)-min(yy,na.rm=TRUE))/bw)
	est <- bkde2D(cbind(xx,yy),bandwidth=c(bwx,bwy),gridsize=c(81, 81))
	est$fhat=est$fhat/max(est$fhat,na.rm=TRUE)
	lvs=c(0.05,.2,0.9)
	maxct=max(lvs)
	nlvs=length(lvs)
        thelines=contourLines(est$x1,est$x2,est$fhat,levels=lvs)
        polygon(thelines[[1]]$x,thelines[[1]]$y,col=col1,border="red",lwd=1.1)
        polygon(thelines[[2]]$x,thelines[[2]]$y,col=NA,border="blue",lwd=.2)
    polygon(thelines[[3]]$x,thelines[[3]]$y,col=NA,border="black",lwd=.1)
    text(median(xx,na.rm=TRUE),median(yy,na.rm=TRUE),labels=title)


	#polygon(thelines[[4]]$x,thelines[[4]]$y,col="white",border=col1,lwd=1)
	
	
}




FA1 <- data_forcast(out,sdcon,z,years,1,marinescale=1,mortR=2,DPE=1,DPS=1,river=1,ssize=100)
print(FA1$conmat)




results <- function(numbersims=100000,morttype=2)
    {


#Foster    
FA1 <- data_forcast(out,sdcon,z,years,1,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims)
FA2a <- data_forcast(out,sdcon,z,years,2,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims)
FA2b <- data_forcast(out,sdcon,z,years,3,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims)
FA3a <- data_forcast(out,sdcon,z,years,4,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims)
FA3b <- data_forcast(out,sdcon,z,years,5,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims)
FA4 <- data_forcast(out,sdcon,z,years,6,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims)
FNAA <- data_forcast(out,sdcon,z,years,7,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=1,ssize=numbersims)

f1 <- qplot(FA1,"FA1","Foster A1")
f2 <- qplot(FA2a,"FA2a","Foster A2a")
f3 <- qplot(FA2b,"FA2b","Foster A2b")
f4 <- qplot(FA3a,"FA3a","Foster A3a")
f5 <- qplot(FA3b,"FA3b","Foster A3b")
f6 <- qplot(FA4,"FA4","Foster A4")
f7 <- qplot(FNAA,"FNAA","Foster NAA")
summaryt <- data.frame("FOS"=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold","SAR","Marine Survival","DPS","DPE","DPS*DPE"),
                       FA1=f1$table["Median"],FA2a=f2$table["Median"],FA2b=f3$table["Median"],FA3a=f4$table["Median"],FA3b=f5$table["Median"],FA4=f6$table["Median"],FNAA=f7$table["Median"])
colnames(summaryt) <- c("FOS","A1","A2a","A2b","A3a","A3b","A4","NAA")
write.csv(file="Foster.Summary.csv",summaryt)
svglite(file="Foster_trade.svg")
        rgb.val <- col2rgb("grey27")
        alpha1=50
plot(NA,NA,xlim=c(0,3), ylim=c(0,2000), ylab="Geo mean Spawners", xlab= "Recruits per Spawner") 
        egg(f1$RS,f1$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A1")
        egg(f2$RS,f2$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2A")
        egg(f3$RS,f3$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2B")
        egg(f4$RS,f4$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3a")
        egg(f5$RS,f5$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3b")
        egg(f6$RS,f6$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A4")
        egg(f7$RS,f7$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="NAA")
        dev.off()
## #Greenpeter
GPA1 <- data_forcast(out,sdcon,z,years,1,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims)
GPA2a <- data_forcast(out,sdcon,z,years,2,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims)
GPA2b <- data_forcast(out,sdcon,z,years,3,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims)
GPA3a <- data_forcast(out,sdcon,z,years,4,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims)
GPA3b <- data_forcast(out,sdcon,z,years,5,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=2,ssize=numbersims)

gp1 <- qplot(GPA1,"GPA1","Green Peter A1")
gp2 <- qplot(GPA2a,"GPA2a","Green Peter A2a")
gp3 <- qplot(GPA2b,"GPA2b","Green Peter A2b")
gp4 <- qplot(GPA3a,"GPA3a","Green Peter A3a")
gp5 <- qplot(GPA3b,"GPA3b","Green Peter A3b")
summaryt <- data.frame("GRP"=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold","SAR","Marine Survival","DPS","DPE","DPS*DPE"),
                      GPA1=gp1$table["Median"],GPA2a=gp2$table["Median"],GPA2b=gp3$table["Median"],GPA3a=gp4$table["Median"],GPA3b=gp5$table["Median"])
colnames(summaryt) <- c("GP","A1","A2a","A2b","A3a","A3b")
write.csv(file="Green.Peter.Summary.csv",summaryt)

svglite(file="Green_trade.svg")
       rgb.val <- col2rgb("grey27")
        alpha1=50
plot(NA,NA,xlim=c(0,3), ylim=c(0,2000), ylab="Geo mean Spawners", xlab= "Recruits per Spawner") 
        egg(gp1$RS,gp1$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A1")
        egg(gp2$RS,gp2$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2A")
        egg(gp3$RS,gp3$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2B")
        egg(gp4$RS,gp4$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3a")
        egg(gp5$RS,gp5$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3b")
        dev.off()
#


        
#Greenpeter
DETA1 <- data_forcast(out,sdcon,z,years,1,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims)
DETA2a <- data_forcast(out,sdcon,z,years,2,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims)
DETA2b <- data_forcast(out,sdcon,z,years,3,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims)
DETA3a <- data_forcast(out,sdcon,z,years,4,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims)
DETA3b <- data_forcast(out,sdcon,z,years,5,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims)
DETA4 <- data_forcast(out,sdcon,z,years,6,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims)
DETNAA <- data_forcast(out,sdcon,z,years,7,marinescale=1,mortR=morttype,DPE=1,DPS=1,river=3,ssize=numbersims)


d1 <- qplot(DETA1,"DETA1","Detroit A1")
d2 <- qplot(DETA2a,"DETA2a","Detroit A2a")
d3 <-qplot(DETA2b,"DETA2b","Detroit A2b")
d4 <-qplot(DETA3a,"DETA3a","Detroit A3a")
d5 <-qplot(DETA3b,"DETA3b","Detroit A3b")
d6 <-qplot(DETA4, "DETA4","Detroit A4")
d7 <-qplot(DETNAA,"DNAA","Detroit NAA")

summaryt <- data.frame("DET"=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold","SAR","Marine Survival","DPS","DPE","DPS*DPE"),
                       DA1=d1$table["Median"],DA2a=d2$table["Median"],DA2b=d3$table["Median"],DA3a=d4$table["Median"],DA3b=d5$table["Median"],DA4=d6$table["Median"],DNAA=d7$table["Median"])
colnames(summaryt) <- c("DET","A1","A2a","A2b","A3a","A3b","A4","NAA")
write.csv(file="Detroit.Summary.csv",summaryt)

     svglite(file="Det_trade.svg")
       rgb.val <- col2rgb("grey27")
        alpha1=50
plot(NA,NA,xlim=c(0,3), ylim=c(0,2000), ylab="Geo mean Spawners", xlab= "Recruits per Spawner") 
        egg(d1$RS,d1$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A1")
        egg(d2$RS,d2$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2A")
        egg(d3$RS,d3$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2B")
        egg(d4$RS,d4$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3a")
        egg(d5$RS,d5$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3b")
        egg(d6$RS,d6$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A4")
        egg(d7$RS,d7$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="NAA")
        dev.off()
## #Greenpeter
SSA1 <- data_forcast2(out,sdcon,z,years,1,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims)
SSA2a <- data_forcast2(out,sdcon,z,years,2,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims)
SSA2b <- data_forcast2(out,sdcon,z,years,3,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims)
SSA3a <- data_forcast2(out,sdcon,z,years,4,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims)
SSA3b <- data_forcast2(out,sdcon,z,years,5,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims)
SSA4 <- data_forcast2(out,sdcon,z,years,6,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims)
SSNAA <- data_forcast2(out,sdcon,z,years,7,marinescale=1,mortR=morttype,DPE=1,DPS=1,ssize=numbersims)

ss1 <-qplot(SSA1,"SSA1","South Sant. A1")
ss2 <-qplot(SSA2a,"SSA2a","South Sant.A2a")
ss3 <-qplot(SSA2b,"SSA2b","South Sant. A2b")
ss4 <-qplot(SSA3a,"SSA3a","South Sant. A3a")
ss5 <-qplot(SSA3b,"SSA3b","South Sant. A3b")
ss6 <-qplot(SSA4, "SSA4","South Sant. A4")
ss7 <-qplot(SSNAA,"SSNAA","South Sant. NAA")
summaryt <- data.frame("SS"=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold","SAR","Marine Survival","DPS","DPE","DPS*DPE"),
                       SSA1=ss1$table["Median"],SSA2a=ss2$table["Median"],SSA2b=ss3$table["Median"],SSA3a=ss4$table["Median"],SSA3b=ss5$table["Median"],SSA4=ss6$table["Median"],SSNAA=ss7$table["Median"])
colnames(summaryt) <- c("SS","A1","A2a","A2b","A3a","A3b","A4","NAA")
                write.csv(file="South.Sant.Summary.csv",summaryt)

     svglite(file="South_trade.svg")
       rgb.val <- col2rgb("grey27")
        alpha1=50
plot(NA,NA,xlim=c(0,3), ylim=c(0,2000), ylab="Geo mean Spawners", xlab= "Recruits per Spawner") 
        egg(ss1$RS,ss1$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A1")
        egg(ss2$RS,ss2$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2A")
        egg(ss3$RS,ss3$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A2B")
        egg(ss4$RS,ss4$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3a")
        egg(ss5$RS,ss5$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A3b")
        egg(ss6$RS,ss6$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="A4")
        egg(ss7$RS,ss7$Abun,col1=rgb(rgb.val[1], rgb.val[2], rgb.val[3],max=255,alpha=alpha1),title="NAA")
        dev.off()
#  

        
    }


#results(10000)

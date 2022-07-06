library(xtable)
table_hessian <- function(out,z,sdcon)
{
return(as.data.frame(eigen(z$hessian)$values))

}   


table_z <- function(out,z,sdcon)
{

    e <- z$par
    s <- sdcon
    out <-as.data.frame(rbind(e,s))

    if(length(s)==3)
    {
  
        colnames(out) <- c("Fresh Survival","Init","SD")
        rownames(out) <- c("Estimate","SE")
        }
    if(length(s)==4)
        {
        colnames(out) <- c("Fresh Survival 1","Fresh Survival 2","Init","SD")
        rownames(out) <- c("Estimate","SE")
        }

    if(length(s)>4)
    {
            e <- z$par[1:3]
            s <- sdcon[1:3]
           out <-as.data.frame(rbind(e,s))

        colnames(out) <- c("Fresh Survival 1","Init","SD")
        rownames(out) <- c("Estimate","SE")
        }

    
return(out)

}   

table_covariance <- function(out,z,sdcon)
{

    cov <- as.data.frame(solve(z$hessian))

    if(length(sdcon)==3)
        {
        colnames(cov) <- c("Fresh Survival","Init","SD")
        }
    if(length(sdcon)==4)
    {
 colnames(cov) <- c("Fresh Survival 1","Fresh Survival 2","Init","SD")
        }

    return(cov)
}

table_dev <- function(out,z,sdcon)
{

    out <- as.data.frame(z$par[4:length(z$par)])

    if(length(sdcon)>4)
    {
 colnames(out) <- c("Dev")
        }

    return(out)
}

table_corr <- function(out,z,sdcon)
{
     c <- as.data.frame(solve(z$hessian))
    corr=c/(sdcon%o%sdcon)

    if(length(sdcon)==3)
        {
        colnames(corr) <- c("Fresh Survival","Init","SD")
        }
    if(length(sdcon)==4)
    {
 colnames(corr) <- c("Fresh Survival 1","Fresh Survival 2","Init","SD")
        }
return(corr)
}
table_ages <- function(out,z,sdcon)
{
return(as.data.frame(out$ages))
}   
table_data <- function(out,z,sdcon)
{
return(as.data.frame(out$data))
}   



plot_magic <- function(out,sdcon,z)
{


    print(sdcon)

    if(length(sdcon)==3){
       
        conmat <-matrix(0,nrow=1000,ncol=(years))

        conmat <-matrix(0,nrow=1000,ncol=(years))
        tab= MASS::mvrnorm(1000,z$par,solve(z$hessian))
        confi <- tab[,1]
        confi2 <- tab[,2]
       
        for( i in 1:1000)
        {
            outm <- modelinits(confi[i],confi2[i],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
            conmat[i,] <- outm$ages$SpawnerAdult
        }
    }


    
    if(length(sdcon)==4){
        print("WTF")
        print(z$par)
        print(sdcon)
        print("sowhat")
    print("sillygoose")
    print(length(data$SpawnerYr))
    print(length(data$FosterCount))
       
        conmat <-matrix(0,nrow=1000,ncol=(years))
        tab= MASS::mvrnorm(1000,z$par,solve(z$hessian))
        confi <- tab[,1]
        confi2 <- tab[,2]
        confi3 <- tab[,3]
        print("So am I an idiot")
        print(years)
        print(length(conmat[1,]))
        
        for( i in 1:1000)
        {
            outm2 <- modelinit(confi[i],confi2[i],confi3[i],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
            conmat[i,] <- outm2$ages$SpawnerAdult
        }
        print(length(outm2$ages$SpawnerAdult))
    }



    if(length(sdcon)>4){
        print("WTF")
        print(z$par)
        print(sdcon)
        print("sowhat")
    print("sillygoose")
    print(length(data$SpawnerYr))
    print(length(data$FosterCount))
       
        conmat <-matrix(0,nrow=1000,ncol=(years))
        tab= MASS::mvrnorm(1000,z$par,solve(z$hessian))
        confi <- tab[,1]
        confi2 <- tab[,2]
        confi3 <- tab[,4:length(z$par)]
        print("So am I an idiot")
        print(years)
        print(length(conmat[1,]))
        
        for( i in 1:1000)
        {
            outm2 <- modelinitdev(confi[i],confi2[i],confi3[i,],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
            conmat[i,] <- outm2$ages$SpawnerAdult
        }
        
    }

    
     par(mfrow=c(1,2))
    plot(data$SpawnerYr,data$FosterCount,ylim=c(0,600),xlab="Time", ylab="Foster Count Observed")
        u <- colQuantiles(conmat,probs=.975,na.rm=TRUE)
        l<- colQuantiles(conmat,probs=.0275,na.rm=TRUE)
        uu <- colQuantiles(conmat,probs=.80,na.rm=TRUE)
        ll<- colQuantiles(conmat,probs=.20,na.rm=TRUE)
        lines(data$SpawnerYr,u,col="black",lwd=.5)
        lines(data$SpawnerYr,l,col="black",lwd=.5)
        lines(data$SpawnerYr,uu,col="grey",lwd=1)
        lines(data$SpawnerYr,ll,col="grey",lwd=1)
    legend("topright",legend=paste("NLL (with prior)::",signif(z$value,2)),box.col="white")
    points(data$SpawnerYr,data$FosterCount)
    points(data$SpawnerYr[1:7],out$ages$SpawnerAdult[1:7],col="magenta")
    lines(out$data$SpawnerYr,out$ages$SpawnerAdult,type="l",col="red")
    lines(out$data$SpawnerYr,c(data$smoltsurvival)*1000,type="l",col="green")
    #abline(v=2006,col="blue")
    ## plot(data$FosterCount[8:(years-1)],out$ages$SpawnerAdult[8:(years-1)],ylim=c(0,400),type="o",xlab="Foster Count Observed", ylab="Foster Count Predicted")
    ## linm=lm(out$ages$SpawnerAdult[8:(years-1)]~data$FosterCount[8:(years-1)])
    ## abline(linm,col="red")
    ## print(summary(linm))
    if(length(z$par)<=4){
    yres=(-log(out$ages$SpawnerAdult[8:(years)])+log(data$FosterCount[8:(years)]))
    plot(data$SpawnerYr[1:(years)],data$SpawnerYr[1:(years)]*NA,type="h",col="blue",ylab="log(obs)-log(pred)",xlab="Time",ylim=c(min(yres,na.rm=TRUE),max(yres,na.rm=TRUE)))
    points(data$SpawnerYr[8:(years)],yres,type="h",col="blue",ylab="residual",xlab="Time",ylab=NA,xlab=NA)
    }
    if(length(z$par)>4)
    {   plot(data$SpawnerYr[1:(years)],exp(z$par[4:(length(z$par))]),col="red",type="o",ylab="exp(Dev)",xlab="Time")
        ## axis(side = 1, at = pretty(range(z$par[4:(length(z$par))])))
        ## mtext("exp(Dev)", side = 4, line = 1)  
        }
    
                                               
}







wormplot <- function(mat,main,main2,years=years,ymax=-1)
{
    u <- colQuantiles(mat,probs=.975,na.rm=TRUE)
    l<- colQuantiles(mat,probs=.0275,na.rm=TRUE)
    uu <- colQuantiles(mat,probs=.80,na.rm=TRUE)
    ll<- colQuantiles(mat,probs=.20,na.rm=TRUE)
    print(max(u,na.rm=TRUE))
    if(ymax<0){
        ymax=min(c(1.1*max(u,na.rm=TRUE),10000))
        ymin=max(c(.9*min(l,na.rm=TRUE),-10))
    print("ymin ymax")
    print(c(ymin,ymax))

       plot(1:(years),colMedians(mat,na.rm=TRUE)[1:(years)],ylim=c(ymin,ymax),type="l",xlab="Time", ylab=main2,main=main)
    }else{
       ymax=min(c(1.1*max(u,na.rm=TRUE),10000))
       ymin=max(c(0,.9*min(l,na.rm=TRUE)),na.rm=TRUE)
           print("ymin ymax")
    print(c(ymin,ymax))
    plot(1:(years),colMedians(mat,na.rm=TRUE)[1:(years)],ylim=c(ymin,ymax),type="l",xlab="Time", ylab=main2,main=main)
    }
    print("ymin ymax")
    print(c(ymin,ymax))
    polygon(c(which(ll>=0), rev(which(ll>=0))), c(ll[which(ll>=0)], rev(uu[which(ll>=0)])), col = "darkgrey")
    lines(1:(years),colMedians(mat,na.rm=TRUE)[1:(years)],type="l",lwd=1.4,col="black")
    for(i in sample(1:600,15))
    {       
        lines(1:(years),mat[i,1:(years)],col="orange",type="l",lwd=.5)
    }
    lines(1:(years),u[1:(years)],col="red",lwd=1.5)
    lines(1:(years),l[1:(years)],col="red",lwd=1.5)
    lines(1:(years),uu[1:(years)],col="black",lwd=1)
    lines(1:(years),ll[1:(years)],col="black",lwd=1)
    abline(v=30,col="blue")
}

histplot <- function(vec,main,main2,med=TRUE)
{
    hist(vec,xlab=main,freq=FALSE,main=main2)
    if(med==TRUE){
        leg1=paste("25%-",signif(quantile(vec,probs=c(.25),na.rm=TRUE),2),"\n","50%-",signif(quantile(vec,probs=c(.5),na.rm=TRUE),2),"\n","75%-",signif(quantile(vec,probs=c(.75),na.rm=TRUE),2),sep="")
    }
   if(med==FALSE){
        leg1=paste("mean-",signif(mean(vec,na.rm=TRUE),2))
    }
    
    legend("topright",legend=leg1,box.col="white")

}    


plot_forcast <- function(list,worm,main,years=70)
    {
if(worm==1){
        par(mfcol=c(2,3))
        wormplot(list$conmat,main,"Spawners",years,ymax=-1)
        wormplot(list$simmat,main,"Simulated Marine Surival",years,-1)
        wormplot(exp(list$devmat),main,"Simulated Marine Surival Dev",years,-1) 
        wormplot(list$RSmat,main,"Recruits/Spawner",years,-1)
        
}
if(worm==2){        
    par(mfrow=c(3,3))
     v=colMeans(list$QETmat)
        plot(v[6:(years)],type="o",xlab="time",ylab="Percent less than 100 Spawners")
        histplot(list$RSvec,"Recruits Per Spawner",main2=main)
        histplot(list$Abunvec,"Geomean Spawners",main2=main)
        histplot(list$QETvec,"QET Threshold 100",main2=main,med=FALSE)
        histplot(list$confi,"Freshwater Survival 2",main2=main)
    histplot(list$marineSvec,"First 5 Years Average Marine Survival",main2=main)
    histplot(list$DPEvec,"DPE",main2=main)
    histplot(list$DPSvec,"DPS",main2=main)
    histplot(list$DPEvec*list$DPSvec,"DPE*DPS",main2=main)


}


if(worm==3){
        par(mfcol=c(1,3))
        wormplot(list$QETamat,main,"Moving Average QET",years,-1)
        wormplot(list$DPSmat,main,"DPS",years,1)
        wormplot(list$DPEmat,main,"DPE",years,1)
}

}

table_forcast <- function(list,main){ #Note only report mean of QET
    s=data.frame(PM=c("Recruits Per Spawner (R/S)","Geomean Spawners (NOR)","QET Threshold -100 (Mean not median)","SAR","First 5 Years Average Marine Survival","DPS","DPE","DPE*DPS"),
                 means=c(mean(list$RSvec,na.rm=TRUE),mean(list$Abunvec,na.rm=TRUE),mean(list$QETvec,na.rm=TRUE),mean(list$confi,na.rm=TRUE),mean(list$marineSvec,na.rm=TRUE),mean(list$DPSvec,na.rm=TRUE),mean(list$DPEvec,na.rm=TRUE),mean(list$DPEvec*list$DPSvec,na.rm=TRUE)),
                 sd=c(sd(list$RSvec,na.rm=TRUE),sd(list$Abunvec,na.rm=TRUE),sd(list$QETvec,na.rm=TRUE),sd(list$confi,na.rm=TRUE),sd(list$marineSvec,na.rm=TRUE),sd(list$DPSvec,na.rm=TRUE),sd(list$DPEvec,na.rm=TRUE),sd(list$DPEvec*list$DPSvec,na.rm=TRUE)),
                 cv=c(sd(list$RSvec,na.rm=TRUE),sd(list$Abunvec,na.rm=TRUE),sd(list$QETvec,na.rm=TRUE),sd(list$confi,na.rm=TRUE),sd(list$marineSvec,na.rm=TRUE),sd(list$DPSvec,na.rm=TRUE),sd(list$DPEvec,na.rm=TRUE),sd(list$DPEvec*list$DPSvec,na.rm=TRUE))/c(mean(list$RSvec),mean(list$Abunvec),mean(list$QETvec),mean(list$confi),mean(list$marineSvec),mean(list$DPSvec,na.rm=TRUE),mean(list$DPEvec,na.rm=TRUE),mean(list$DPEvec*list$DPSvec,na.rm=TRUE)),
                 lower=c(quantile(list$RSvec,.025,na.rm=TRUE),quantile(list$Abunvec,.025,na.rm=TRUE),quantile(list$QETvec,.025,na.rm=TRUE),
                         quantile(list$confi,.025,na.rm=TRUE),quantile(list$marineSvec,.025,na.rm=TRUE),quantile(list$DPSvec,.25,na.rm=TRUE),quantile(list$DPEvec,.25,na.rm=TRUE),quantile(list$DPEvec*list$DPSvec,.25,na.rm=TRUE)),
                 median=c(median(list$RSvec,na.rm=TRUE),median(list$Abunvec,na.rm=TRUE),mean(list$QETvec,na.rm=TRUE),median(list$confi,na.rm=TRUE),
                          median(list$marineSvec,na.rm=TRUE),median(list$DPSvec,na.rm=TRUE),median(list$DPEvec,na.rm=TRUE),median(list$DPEvec*list$DPSvec,na.rm=TRUE)),
                 upper=c(quantile(list$RSvec,.975,na.rm=TRUE),quantile(list$Abunvec,.975,na.rm=TRUE),quantile(list$QETvec,.975,na.rm=TRUE),quantile(list$confi,.975,na.rm=TRUE),quantile(list$marineSvec,.975,na.rm=TRUE)
                        ,quantile(list$DPSvec,.975,na.rm=TRUE),quantile(list$DPEvec,.975,na.rm=TRUE),quantile(list$DPEvec*list$DPSvec,.975,na.rm=TRUE)))
    colnames(s) <- c("","Mean","SD","CV",".0275","Median",".0975")
    return(list(table=s,RS=list$RSvec,Abund=list$Abunvec))
        }

        



        
   





        
        
   



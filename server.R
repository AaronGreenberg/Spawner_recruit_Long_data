#this script plots cpue vs cumulative catch for each area.
#It omits
library(shiny)
library(dplyr)
source("Steel_head_age.R")
source("Plots_tables.R")
library("matrixStats")


shinyServer(function(input, output){
data <- read.csv("data2.csv", sep=",",head=TRUE)
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

 dataInputfit <- reactive({
    inputs <- c(.03,120,1,0*rnorm(38,0,1)) #initializaton
    print("inputs")
    print(input$DevSD)
    print(input$PriorSD)
z=optim(par=inputs,fn=fninitdev,method="BFGS",gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,input$DevSD,input$PriorSD,hessian=TRUE,control=c(reltol=10^-132,maxit=35000))
sdcon <- sqrt(diag(solve(z$hessian)))
outmodeldev <- modelinitdev(z$par[1],z$par[2],z$par[4:length(z$par)],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data) #dangerous global assignment  Possibly will have to fix this
sdcon <- sqrt(diag(solve(z$hessian)))
out <- outmodeldev
sims=simmaker(nsims=10000,years=80,z,sdcon)
return(list(z=z,sdcon=sdcon,out=out,sims=sims))

 })


output$forcastPF <- renderPlot({
 
        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
    z=data$z
print("silly mess")
print(input$ssize)
if(input$ssize==1)
{
    ssmagic=1000
}
if(input$ssize==2)
{
    ssmagic=10000
}

l=data_forcast(out,sdcon,z,70,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=1,ssize=ssmagic,sims=sims)
    text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")

        plot_forcast(l,worm=input$worm,main=text[input$EIS])

   })

output$forcastPG <- renderPlot({
    
        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z
text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")

print("silly mess")
print(input$ssize)
if(input$ssize==1)
{
    ssmagic=1000
}
if(input$ssize==2)
{
    ssmagic=10000
}

lPG=data_forcast(out,sdcon,z,70,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=2,ssize=ssmagic,sims=sims)
 
        text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")
        plot_forcast(LPG,worm=input$worm,main=text[input$EIS])

})

output$forcastPD <- renderPlot({
    
    data=dataInputfit()
    z=data$z
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")

print("silly mess")
print(input$ssize)
if(input$ssize==1)
{
    ssmagic=1000
}
if(input$ssize==2)
{
    ssmagic=10000
}


lPD=data_forcast(out,sdcon,z,70,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=3,ssize=ssmagic,sims=sims)

        text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")
        plot_forcast(lPD,worm=input$worm,main=text[input$EIS])
    
})

output$forcastPJoint <- renderPlot({
    
    data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z
print("silly mess")
print(input$ssize)
if(input$ssize==1)
{
    ssmagic=1000
}
if(input$ssize==2)
{
    ssmagic=10000
}

lPJoint=data_forcast2(out,sdcon,z,70,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,ssize=ssmagic,sims=sims)

        text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")
plot_forcast(lPJoint,worm=input$worm,main=text[input$EIS])
})



output$forcastTF <- renderDataTable({
    
        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
    z=data$z
print("silly mess")
print(input$ssize)
if(input$ssize==1)
{
    ssmagic=1000
}
if(input$ssize==2)
{
    ssmagic=10000
}

l=data_forcast(out,sdcon,z,70,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=1,ssize=ssmagic,sims=sims)
 
    s=table_forcast(l,main=text[input$EIS])
    return(s$table)
})

output$forcastTG <- renderDataTable({
           data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z
text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")

print("silly mess")
print(input$ssize)
if(input$ssize==1)
{
    ssmagic=1000
}
if(input$ssize==2)
{
    ssmagic=10000
}

lPG=data_forcast(out,sdcon,z,70,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=2,ssize=ssmagic,sims=sims)
 
    s=table_forcast(lPG,main=text[input$EIS])
    return(s$table)
})

output$forcastTD <- renderDataTable({
        data=dataInputfit()
    z=data$z
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")

print("silly mess")
print(input$ssize)
if(input$ssize==1)
{
    ssmagic=1000
}
if(input$ssize==2)
{
    ssmagic=10000
}


lPD=data_forcast(out,sdcon,z,70,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=3,ssize=ssmagic,sims=sims)

        s=table_forcast(lPD,main=text[input$EIS])
    return(s$table)
    
})

output$forcastTJoint <- renderDataTable({
    
    data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z
print("silly mess")
print(input$ssize)
if(input$ssize==1)
{
    ssmagic=1000
}
if(input$ssize==2)
{
    ssmagic=10000
}

lPJoint=data_forcast2(out,sdcon,z,70,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,ssize=ssmagic,sims=sims)

    s=table_forcast(lPJoint,main=text[input$EIS])
    return(s$table)
})



output$fitP <- renderPlot({

if(input$SRindex==1){
inputs <- c(.03,.03,100,12)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],z$par[3],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==2){
inputs <- c(.03,100,12)
z=optim(par=inputs,fn=fninits,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinits(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==3){

    data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z
}

print("get happy")
print(z)
print(sdcon)
plot_magic(out,sdcon,z)
    
})



output$zT <- renderDataTable({

if(input$SRindex==1){
inputs <- c(.03,.03,100,12)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],z$par[3],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==2){
inputs <- c(.03,100,12)
z=optim(par=inputs,fn=fninits,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinits(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

if(input$SRindex==3){
        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
    z=data$z
}

print("yo")
return(table_z(out,z,sdcon))
    
})

output$covT <- renderDataTable({
if(input$SRindex==1){
inputs <- c(.03,.03,100,12)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],z$par[3],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==2){
inputs <- c(.03,100,12)
z=optim(par=inputs,fn=fninits,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinits(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==3){

        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
    z=data$z
}
return(table_covariance(out,z,sdcon))
})

output$devT <- renderDataTable({
if(input$SRindex==1){
 

    inputs <- c(.03,.03,100,12)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],z$par[3],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==2){
inputs <- c(.03,100,12)
z=optim(par=inputs,fn=fninits,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinits(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==3){



        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z
}

print("yo")
return(table_dev(out,z,sdcon))
    
})


output$corT <- renderDataTable({
if(input$SRindex==1){
inputs <- c(.03,.03,100,12)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],z$par[3],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==2){
inputs <- c(.03,100,12)
z=optim(par=inputs,fn=fninits,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinits(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==3){



        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z
}

print("yo")
return(table_corr(out,z,sdcon))
    
})


output$dataT <- renderDataTable({
if(input$SRindex==1){
inputs <- c(.03,.03,100,12)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],z$par[3],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==2){
inputs <- c(.03,100,12)
z=optim(par=inputs,fn=fninits,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinits(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==3){

        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z
}

print("yo")
return(table_data(out,z,sdcon))
    
})


output$agesT <- renderDataTable({
if(input$SRindex==1){
inputs <- c(.03,.03,100,12)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],z$par[3],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


if(input$SRindex==2){
inputs <- c(.03,100,12)
z=optim(par=inputs,fn=fninits,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinits(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

if(input$SRindex==3)
{

        data=dataInputfit()
    sims=data$sims
    out=data$out
    sdcon=data$sdcon
        z=data$z

}

return(table_ages(out,z,sdcon))
    
})

})

 

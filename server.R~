#this script plots cpue vs cumulative catch for each area.
#It omits
library(shiny)
library(dplyr)
source("Steel_head_age.R")
source("Plots_tables.R")
library("matrixStats")


shinyServer(function(input, output){
data <- read.csv("data.csv", sep=",",head=TRUE)
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
years <- 16+7

dataInputPF <- reactive({
    inputs <- c(.003,200,100)
    
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
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

lPF=data_forcast(out,sdcon,z,50,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=1,ssize=ssmagic)
 return(lPF) 
})

dataInputPG <- reactive({
inputs <- c(.03,200,100)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
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

lPG=data_forcast(out,sdcon,z,50,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=2,ssize=ssmagic)
 return(lPG) 
})

dataInputPD <- reactive({
inputs <- c(.03,200,100)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
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

lPD=data_forcast(out,sdcon,z,50,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,river=3,ssize=ssmagic)
  return(lPD)
})



dataInputPJoint <- reactive({
inputs <- c(.03,200,100)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)

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

lPJoint=data_forcast2(out,sdcon,z,50,input$EIS,marinescale=input$marinescale,mortR=input$Mort,DPE=input$dpe,DPS=input$dps,ssize=ssmagic)
return(lPJoint)
})

output$forcastPF <- renderPlot({
    l <- as.list(dataInputPF())
    text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")

        plot_forcast(l,worm=input$worm,main=text[input$EIS])

   })
output$forcastPG <- renderPlot({
        text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")
        plot_forcast(dataInputPG(),worm=input$worm,main=text[input$EIS])

})

output$forcastPD <- renderPlot({
        text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")
        plot_forcast(dataInputPD(),worm=input$worm,main=text[input$EIS])
    
})

output$forcastPJoint <- renderPlot({
        text <- c("Alt1","Alt2a","Alt2b","Alt3a","Alt3b","Alt4","AltNAA","General","Used for estimation")
plot_forcast(dataInputPJoint(),worm=input$worm,main=text[input$EIS])
})



output$forcastTF <- renderDataTable({
    s=table_forcast(dataInputPF(),main=text[input$EIS])
    return(s$table)
})

output$forcastTG <- renderDataTable({
    s=table_forcast(dataInputPG(),main=text[input$EIS])
    return(s$table)
})

output$forcastTD <- renderDataTable({
        s=table_forcast(dataInputPD(),main=text[input$EIS])
    return(s$table)
    
})

output$forcastTJoint <- renderDataTable({
    s=table_forcast(dataInputPJoint(),main=text[input$EIS])
    return(s$table)
})



output$fitP <- renderPlot({
if(input$SRindex==1){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnconst,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelconst(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
    }

if(input$SRindex==2){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fn,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- model(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

if(input$SRindex==3){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnlin,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modellin(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}
if(input$SRindex==4){
inputs <- c(.03,200,100)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


plot_magic(out,sdcon,z)
    
})



output$zT <- renderDataTable({
 print("So do you feel me?")
if(input$SRindex==1){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnconst,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelconst(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
    }

if(input$SRindex==2){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fn,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- model(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

if(input$SRindex==3){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnlin,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modellin(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}
if(input$SRindex==4){
inputs <- c(.03,200,100)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

print("yo")
return(table_z(out,z,sdcon))
    
})

output$covT <- renderDataTable({
if(input$SRindex==1){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnconst,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelconst(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
    }

if(input$SRindex==2){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fn,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- model(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

if(input$SRindex==3){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnlin,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modellin(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}
if(input$SRindex==4){
inputs <- c(.03,200,100)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

print("yo")
return(table_covariance(out,z,sdcon))
    
})

output$dataT <- renderDataTable({
if(input$SRindex==1){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnconst,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelconst(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
    }

if(input$SRindex==2){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fn,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- model(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

if(input$SRindex==3){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnlin,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modellin(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}
if(input$SRindex==4){
inputs <- c(.03,200,100)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}


print("yo")
return(table_data(out,z,sdcon))
    
})


output$agesT <- renderDataTable({
if(input$SRindex==1){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnconst,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelconst(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
    }

if(input$SRindex==2){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fn,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- model(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}

if(input$SRindex==3){
inputs <- c(.03,10)
z=optim(par=inputs,fn=fnlin,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modellin(z$par[1],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}
if(input$SRindex==4){
inputs <- c(.03,200,100)
z=optim(par=inputs,fn=fninit,gr=NULL,DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data,hessian=TRUE,control=c(reltol=10^-32,maxit=5000))
sdcon <- sqrt(diag(solve(z$hessian)))
out <- modelinit(z$par[1],z$par[2],DPSDPE,years,Percent4,Percent5,Percent6,Repeat6,data)
}
return(table_ages(out,z,sdcon))
    
})
})


 

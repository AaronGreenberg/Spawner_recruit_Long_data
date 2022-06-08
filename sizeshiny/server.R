#this script plots cpue vs cumulative catch for each area.
#It omits
library(shiny)
library(dplyr)
shinyServer(function(input, output){
data=read.csv("Sex.csv")

dataInput <- reactive({
    tmpdata <- filter(data,year==input$year,month==input$month,sex==1|sex==1)
    males=tmpdata$size
    tmpdata <- filter(data,year==input$year,month==input$month,sex==3|sex==4|sex==5)
    females=tmpdata$size
    if(input$month==13)
        {
    tmpdata <- filter(data,year==input$year,sex==1)
    males=tmpdata$size
    tmpdata <- filter(data,year==input$year,sex==3|sex==4|sex==5)
    females=tmpdata$size

        }
    if(length(males)==0){males=-10}
    if(length(females)==0){females=-10}
    out=list(males=males,females=females)
return(out)
})


pyramid <- function(males,females,type)
    {
        maxmales <- max(males,na.rm=TRUE)
        maxfemales <- max(females,na.rm=TRUE)
        bigest <- ceiling(max(c(maxmales,maxfemales),na.rm=TRUE))

        if(bigest>0){
            cutvector <- seq(0,bigest+10,by=10)
        }
        else{cutvector=c(-10:0)}
        cmales=hist(males,cutvector,plot=FALSE)
        cfemales=hist(females,cutvector,plot=FALSE)
        if((max(cmales$counts)>1)|(max(cfemales$counts)>1)){
        if(type=="count")
            {
                plot(c(NA,NA),
                     xlim=c(-max(cmales$counts,na.rm=TRUE)*1.1,max(cfemales$counts,na.rm=TRUE)),ylim=c(0,230),ylab="size",xlab="Count")
                rect(xleft=-cmales$counts,xright=0,ybottom=cutvector-10,ytop=cutvector,col="lightblue",border="black")
                rect(xright=cfemales$counts,xleft=0,ybottom=cutvector-10,ytop=cutvector,col="white",border="black")
             }   
        if(type=="freq")
            {
                plot(c(NA,NA),
                     xlim=c(-max(cmales$counts/sum(cmales$counts),na.rm=TRUE)*1.1,max(cfemales$counts/sum(cfemales$counts),na.rm=TRUE)*1.1),
                     ylim=c(0,230),ylab="size",xlab="Percent")

                rect(xleft=-cmales$counts/sum(cmales$counts),xright=0,ybottom=cutvector-10,ytop=cutvector,col="lightblue",border="black")
                rect(xright=cfemales$counts/sum(cfemales$counts),xleft=0,ybottom=cutvector-10,ytop=cutvector,col="white",border="black")
             }   
        abline(h=165,col="black")
        #abline(h=165*1.069+18,col="blue")
            #abline(h=220,col="green")
            abline(h=mean(males,na.rm=TRUE),col="red")
        legend("topleft",legend=paste("Males",max(cmales$counts),sep="::"))
            legend("topright",legend=paste("Females",max(cfemales$counts),sep="::" ))
            legend("top",legend=paste("Mean Male",signif(mean(males,na.rm=TRUE)),sep="::" ))
    } 
        
    }   

output$countPlot <- renderPlot({
    tmp <- dataInput()
    pyramid(tmp$males,tmp$females,type="count")
})


output$percentPlot <- renderPlot({
    tmp <- dataInput()
    pyramid(tmp$males,tmp$females,type="freq")
})


})


input <- read.csv("Sex.csv")
years=2003:2007
print(years)
flag=FALSE
#plot(NA,NA,ylim=c(80,250),xlim=c(2000,2010))
for(i in years)
{
x11()    
print(i)
inputtmp<-subset(input,year==i)
    inputtmp <- subset(inputtmp,sex==1)
    inputapril<-subset(inputtmp,month==5)
    if(length(inputapril$size)>0)
        {
            print(mean(inputapril$size,na.rm=TRUE))
            print(min(inputapril$size,na.rm=TRUE))
        fcol=rgb(1, runif(1,0,1), runif(1,0,1),1)
        hist(inputapril$size,col=fcol,50,)
        abline(v=mean(inputapril$size),col=fcol)
    }
    rm(inputapril)
    flag=TRUE
}

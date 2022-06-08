data=read.csv("observed.csv",sep="\t",head=TRUE)




model <- function(par)
{
    nll <- -sum(log(dnorm((data$smoltsurvival[1:26]/100),(par[1]*sin(par[2]*1:26+par[3])+par[4]),par[5])))
}

ipar=c(-1.05,25,15,0,1)
fit=optim(ipar,model,hessian=TRUE,control=c(maxit=5000,reltol=10^-32))
graphics.off()
plot(1:26,fit$par[1]*sin(fit$par[2]*1:26+fit$par[3])+fit$par[4],col="red",type="l",ylim=c(0,.06))
points(data$smoltsurvival[1:26]/100, type="p",col="blue",pch=16)

## for(i in 1:500){
##     s <- (fit$par[1]*sin(fit$par[2]*1:26+fit$par[3])+fit$par[4])
##     sim <- s+rnorm(length(1:26),0,(fit$par[5])/sqrt(26))
##     print(sim[which(sim<0)])
##     points(sim,col="magenta",cex=.2)
## }

print(fit)
print(eigen(fit$hessian))
print(solve(fit$hessian))
fit$V=solve(fit$hessian)
fit$sig=sqrt(diag(fit$V)) #std deviation in parameters
print("correlation")
print(signif(fit$V/(fit$sig%o%fit$sig),3))


## modelautocor <- function(pari)
## {
##     print(pari)
##     obs <- data$smoltsurvival[1:26]/100
##     pred <- 1:26*0
##     pred[1] <- pari[1]
##     for(i in 2:26)
##         {
##         pred[i] <- pred[i-1]*pari[2]+pari[3]

##         }
## return(out=list("eps"=-sum(log(dnorm(obs,pred,pari[4]))),"pred"=pred,"obs"=obs))
## }

## ipar=(c(.02,.001,.001,.01))

## fn <-function(par){modelautocor(par)$eps}
## fit=optim(ipar,fn,hessian=TRUE,control=c(maxit=5000,reltol=10^-32))
## print(fit)
## print(solve(fit$hessian))

## fit$V=solve(fit$hession)
## fit$sig=sqrt(diag(fit$V)) #std deviation in parameters
## print("correlation")
## print(fit$V/(fit$sig%o%fit$sig))


## out <- modelautocor(fit$par)
## plot(1:26,out$obs,ylim=c(-.01,.1))
## lines(1:26,out$pred)

## for(i in 1:50){
##     s <- 1:26*0
##     s[1] <- fit$par[1]
##     for(i in 2:26)
##     {
##         s[i] <- s[i-1]*fit$par[2]+fit$par[3]+rnorm(1,0,fit$par[4])

##         }
##     s[which(s<=0)] <- 0              
##     sim <- s
##     print(sim[which(sim<0)])
##     lines(sim,col="lightgrey")
##     points(sim,col="magenta",cex=.2)
## points(1:26,out$obs,ylim=c(-.01,.1))
## lines(1:26,out$pred)

## }

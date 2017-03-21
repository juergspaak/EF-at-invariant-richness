
ResultsFolder <- "/Desktop/"



alpha <- seq(-0.9,-0.1,by=0.01)
n<-10
muav <- 1
fav <- 1
eav <- 0.5
b <- c(0.5,1,1.5,2.0) #(dis)proportionality between mur and er
c <- c(0.5,1,1.5,2.0) #(dis)proportionality between fr and er
c_f <- 1
c_mu <- -1
bcs <- expand.grid(b,c)
ts <- seq(eav/sqrt(3),0.1,-0.05) #relative sensitivities; arbitrary choice
ylims <- c(0,5)

comp <- function(n,alpha)
{
  competition <- -alpha*n/(1-alpha*(n-1))
}

deltaEF <- function(alpha, muav, fav, eav,
                    t, b, c, n, c_f, c_mu)
{
  deltaEF <- -n*fav*(muav*eav+c_mu*t^b)/(alpha+1)*(c_f*t^c*(c_mu*t^b+t)/(1+c_mu*t^(b+1)) +1 -comp(n,alpha))
  return(deltaEF)
}

plotDeltaEF <- function(alpha, muav, fav, eav,
                        ts, b, c, n, c_f, c_mu)
{
  for (i in c(1:length(ts)))
  {
    t  <- ts[i]
    #check for coexistence in stressed condition
    comp <- comp(n,alpha)
    tLimit1 <- (1+sqrt(3)*c_mu*t^b)*(1/eav-1-sqrt(3)*t)/(1/eav-1-c_mu*t^(b+1))
    tLimit2 <- (1-sqrt(3)*c_mu*t^b)*(1/eav-1+sqrt(3)*t)/(1/eav-1-c_mu*t^(b+1))
    coexStress <- (comp<tLimit1)*(comp<tLimit2)
    #check for coexistence in control condition
    #...and eliminate cases where not possible
    coexControl <- t^b<(1-comp)/sqrt(3)
    deltaEFCase <-deltaEF(alpha, muav, fav, eav,
                          t=t, b=b, c=c, n, c_f, c_mu)
    deltaEFCase[which(coexControl+coexStress==0)] <- NA
    colour <- i
    colour <- colour/length(ts)
    lines(alpha,deltaEFCase, col=rgb(red=colour, 
                                     green=1-colour, blue=0))
  }
}

pdf("Analytic11.pdf")

par(mar=c(3,4,2,0.5), las=1, mfrow=c(4,4), 
    tck=-0.02, mgp=c(2,0.5,0))

for (i in 1:nrow(bcs))
{
  plot(alpha,rep(NA,length(alpha)), 
       ylim=ylims, main=paste("b=",bcs[i,1], 
                              "c=",bcs[i,2]),
       ylab="Delta EF",
       xlim=c(max(alpha),min(alpha)))
  abline(h=0)
  plotDeltaEF(alpha, muav, fav, eav,
              ts=ts, b=bcs[i,1], c=bcs[i,2], n, c_f, c_mu)
  colour <- c(1:length(ts))/length(ts)
  if (i==nrow(bcs))
  {
    legend("topleft", legend=ts, pch="", lty="solid", cex=0.5,
           ncol=2,
           col=rgb(red=colour, green=1-colour, blue=0))
  }
}

dev.off()

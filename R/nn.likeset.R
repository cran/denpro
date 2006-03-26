nn.likeset<-function(dendat,radmat,k,p)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

volunitball<-volball(1,d)

radit<-radmat[,k]
evat<-k/(n*radit^d*volunitball)
maksi<-max(evat)
lambda<-p*maksi
grt<-(evat>=lambda)

dendatsub<-dendat[grt,]

return(dendatsub)
}


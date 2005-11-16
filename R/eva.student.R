eva.student<-function(x,t,marginal="unif",sig=c(1,1),r=0)
# t>2 
#  sig is std of marginals
{
if (marginal=="unif"){
   u<-x[1]/sig[1]+1/2
   v<-x[2]/sig[2]+1/2
   marg1<-1/sig[1]
   marg2<-1/sig[2]
}
if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[1])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}
if (marginal=="student"){
   u<-pt(x[1]/sig[1],df=t)
   v<-pt(x[2]/sig[2],df=t)
   marg1<-dt(x[1]/sig[1],df=t)/sig[1]
   marg2<-dt(x[2]/sig[2],df=t)/sig[2]
}

d<-2
x1<-qt(u,df=t)
x2<-qt(v,df=t)
produ<-dt(x1,df=t)*dt(x2,df=t)

nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
vakio<-gamma((t+d)/2)/((t*pi)^(d/2)*gamma(t/2))
g<-vakio*(1-r^2)^(-1/2)*(1+nelio/t)^(-(t+d)/2)
#g<-(1-r^2)^(1/2)*(1+(x1^2+x2^2-2*r*x1*x2)/(t*(1-r^2)))^(-(t+d)/2)

val<-g/produ*marg1*marg2

return(val)
}


eva.copula<-function(x,type="gauss",marginal="unif",sig=c(1,1),r=0,t=1,g=1)
{
# sig is std of marginals, r is the correlation coeff, 
# t is deg of freedom/param of Gumbel

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
if (type=="gauss"){
   d<-2
   x1<-qnorm(u,sd=1)
   x2<-qnorm(v,sd=1)
   produ<-dnorm(x1,sd=1)*dnorm(x2,sd=1)

   nelio<-(x1^2+x2^2-2*r*x1*x2)/(1-r^2)
   vakio<-(2*pi)^(-d/2) 
   g<-vakio*(1-r^2)^(-1/2)*exp(-(1/2)*nelio)

   val<-g/produ*marg1*marg2
}
if (type=="gumbel"){
  linku<-(-log(u))^g
  linkv<-(-log(v))^g
  a<-linku+linkv
  der1u<--g*(-log(u))^(g-1)/u
  der1v<--g*(-log(v))^(g-1)/v
  b<-exp(-a^(1/g))
  der1b<--g*(-log(b))^(g-1)/b
  der2b<-g*b^(-2)*(-log(b))^(g-2)*(g-1-log(b))
  psi<--der2b*der1b^(-3)

  val<-psi*der1u*der1v*marg1*marg2
}

return(val)
}







eva.gumbel<-function(x,t,marginal="unif",sig=c(1,1))
# t>=1
{
u<-x[1]
v<-x[2]
marg1<-1
marg2<-1

if (marginal=="normal"){
   u<-pnorm(x[1]/sig[1])
   v<-pnorm(x[2]/sig[2])
   marg1<-evanor(x[1]/sig[1])/sig[1]
   marg2<-evanor(x[2]/sig[2])/sig[2]
}

#linku<-(-log(u))^t
#linkv<-(-log(v))^t
#a<-linku+linkv
#der1u<--t*log(u)^(t-1)/u
#der1v<--t*log(v)^(t-1)/v
#der1a<--t*log(a)^(t-1)/a
#der2a<-a^(-1)*log(a)^(t-2)*(a^(-1)*log(a)-t*(t-1))
#val<--der2a*der1a^(-2)*der1u*der1v*marg1*marg2

copu<-exp(-((-log(u))^t+(-log(v))^t)^(1/t))
val<-copu*(u*v)^(-1)*(-log(u))^(t-1)*(-log(v))^(t-1)*((-log(u))^t+(-log(v))^t)^(1/t-2)*(t-1+((-log(u))^t+(-log(v))^t)^(1/t))

return(val)
}


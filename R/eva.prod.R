eva.prod<-function(x,marginal="student",t=1)
{
D<-length(x)
if (marginal=="student"){
    d<-1
    vakio<-gamma((t+d)/2)/((t*pi)^(d/2)*gamma(t/2))
    x1<-vakio*(1+x[1]^2/t)^(-(t+d)/2)
    x2<-vakio*(1+x[2]^2/t)^(-(t+d)/2)
    val<-x1*x2
}
if (marginal=="studentR"){
    #x1<-dt(x[1],df=t)
    #x2<-dt(x[2],df=t)
    y<-dt(x,df=t)
    val<-prod(y)  
}

return(val)
}


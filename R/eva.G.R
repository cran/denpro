eva.G<-function(x,vmax=1,pmax=1,contrast="loglik")
{
v<-x[1]
p<-x[2]
if (contrast=="loglik") res<--p*log(p/v)-(pmax-p)*log((pmax-p)/(vmax-v))
else res<--p^2/v-(pmax-p)^2/(vmax-v)

return(res)
}



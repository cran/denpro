sim.exa<-function(n=NULL,seed=1,N=NULL,type="mulmod")
{

if (type=="mulmod") return( sim.mulmod(n=n,seed=seed,N=N) )

if (type=="fox") return( sim.fox(n=n,seed=seed,N=N) )

}

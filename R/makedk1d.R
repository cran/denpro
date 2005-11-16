makedk1d<-function(dendat,h,binlkm){
#Makes a discretized kernel estimator
#
#dendat is n-vector
#h>0 is smoothing parameter
#binlkm>0 determines how many bins for every dimension will be made, 
#f.ex.binlkm=n
#
#Returns list(values,recs)
#values is M-vector, M \approx n)
#recs is M*2-matrix
# 
d<-1
n<-length(dendat)
#N<-floor(n^{1/d})
# valitaan hilan koko se. lyhimmallekin vaihteluvalille tulee N jakoa
supp<-matrix(0,1,2)
supp[1]<-min(dendat)
supp[2]<-max(dendat)
minvaihtelu<-supp[2]-supp[1]+2*h
askel<-minvaihtelu/binlkm
lkmi<-ceiling((supp[2]-supp[1]+2*h)/askel)
itemlkm<-2*ceiling(lkmi/2)   #itemlkm is M
keskip<-(supp[2]+supp[1])/2
alat<-keskip-(lkmi/2)*askel
ylat<-keskip+(lkmi/2)*askel
#  
recs<-matrix(0,itemlkm,2)
values<-matrix(0,itemlkm,1)
for (itemind in 1:itemlkm){ #lasketaan arvot recs:iin ja values:iin
  # arvot recs:iin
  recs[itemind,1]<-alat+askel*(itemind-1)
  recs[itemind,2]<-alat+askel*itemind
  # arvot values:iin
  val<-0
  for (k in 1:n){
        tulo<-1
        hav<-dendat[k]
        ala<-recs[itemind,1]
        yla<-recs[itemind,2]
        mi<-min(1,(yla-hav)/h)
        ma<-max(-1,(ala-hav)/h)
        if (ma<mi){
          li<-mi*(1-mi^2/3) 
          va<-ma*(1-ma^2/3) 
          uus<-li-va
        }
        else uus<-0
        tulo<-tulo*uus
        val<-val+tulo
  }
  values[itemind]<-(3/4)*val/n
}
return(list(values=t(values),recs=recs))
}










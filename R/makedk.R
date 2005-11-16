makedk<-function(dendat,h,binlkm){
#Makes a discretized kernel estimator
#
#dendat is n*d-matrix
#h>0 is smoothing parameter
#binlkm>0 determines how many bins for every dimension will be made, 
#f.ex.binlkm=n^(1/d)
#
#Returns list(values,recs)
#values is M=N_1*...*N_d-vector, N_i \approx n^(1/d)
#recs is M*(2*d)-matrix
# 
d<-length(dendat[1,])
n<-length(dendat[,1])
supp<-support(dendat,epsi=0) #supp on d*2-matriisi, kullekk muutt vaihteluvali
#N<-floor(n^{1/d})
# valitaan hilan koko se. lyhimmallekin vaihteluvalille tulee N jakoa
minvaihtelu<-min(supp[,2]-supp[,1])+2*h
askel<-minvaihtelu/binlkm
alat<-matrix(0,d,1)
ylat<-matrix(0,d,1)
lkm<-matrix(0,d,1)    #kullekin dimensiolle lokeroitten maara
itemlkm<-1            #itemlkm=M
for (i in 1:d){
  lkmi<-ceiling((supp[i,2]-supp[i,1]+2*h)/askel)
  lkmi<-2*ceiling(lkmi/2)
  keskip<-(supp[2*i]+supp[2*i-1])/2
  alat[i]<-keskip-(lkmi/2)*askel
  ylat[i]<-keskip+(lkmi/2)*askel
  lkm[i]<-lkmi
  itemlkm<-itemlkm*lkmi 
}
recs<-matrix(0,itemlkm,2*d)
values<-matrix(0,itemlkm,1)
for (itemind in 1:itemlkm){ #lasketaan arvot recs:iin ja values:iin
  # arvot recs:iin
  dig<-digit(itemind-1,lkm)
  j<-1
  while (j<=d){          #lokeroitten digit
    diggi<-dig[j]
    recs[itemind,2*j-1]<-alat[j]+askel*diggi
    recs[itemind,2*j]<-alat[j]+askel*(diggi+1)
    j<-j+1
  }
  # arvot values:iin
  val<-0
  for (k in 1:n){
     tulo<-1
     for (l in 1:d){
        hav<-dendat[k,l]
        ala<-recs[itemind,2*l-1]
        yla<-recs[itemind,2*l]
        mi<-min(1,(yla-hav)/h)
        ma<-max(-1,(ala-hav)/h)
        if (ma<mi){
          li<-mi*(1-mi^2/3) 
          va<-ma*(1-ma^2/3) 
          uus<-li-va
        }
        else uus<-0
        tulo<-tulo*uus
     }
     val<-val+tulo
  }
  values[itemind]<-(3/4)^d*val/n
}
return(list(values=t(values),recs=recs))
}










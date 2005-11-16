makeker1<-function(eval,levnum,radius){
#
#eval on n-vektori,  esim. eval<-seq(1,n,1)
#levnum on haluttujen tasojoukkojen lkm, 
#radius on levnum-vektori: niitten pallojen sade, jotka sijoitetaan 
#havaintopisteisiin, 1. alinta tasojoukkoa vastaava sade
#
n<-length(eval)
tieto<-matrix(0,levnum,1)
lsets<-matrix(0,levnum,n)   
#lsetscomp<-matrix(NA,levnum,n)
#
M<-max(eval)  #tih f:on maksimi estimoidaan laskemalla maksimi hav pist:ssa
#askel<-M/levnum       #tasot ovat 0,askel,...,M-askel
#tieto<-t(t(seq(alin,M-askel,askel)))  #pal.vakion f:on arvot tietoon
askel<-M/(levnum+1)       #tasot ovat askel,2*askel,...,M-askel
tieto<-t(t(seq(askel,M-askel,askel)))  #pal.vakion f:on arvot tietoon
uustieto<-tieto
newradius<-radius
#
i<-1
ind<-1   #tarvitaan 2. indeksi silla osa tasoista voi olla tyhjia:ei talteen
while (i<=levnum){   #tasot lapi
  j<-1
  while (j<=n){       #havainnot lapi
    if (eval[j]>=tieto[i]){
       lsets[ind,j]<-1 #jos kuuluu tasojoukkoon laitetaan ykkonen       
       #lsetscomp[i,saraind]<-j  saraind<-saraind+1
    }
    j<-j+1
  }
  j<-1           #kaydaan lapi uudestaan ja poistetaan ne tasojoukon
  while (j<=n){  #havainnot, jotka eivat erillaan komplementista
     if (lsets[ind,j]==1){ #jos kuuluu tasojoukkoon
       k<-1   #kaydaan komplementti lapi
       while (k<=n){
         if (lsets[ind,k]==0) #jos kuuluu komplementtiin
           if (etais(dendat[k,],dendat[j,])<radius[i]) #jos tasojoukon piste j
              lsets[ind,j]<-0  #on sade:a lahempana jotain komplementin
                               #pistetta k, poistetaan piste j 
         k<-k+1
       }
     }
     j<-j+1
  }
  if (sum(lsets[ind,])>0){  #jos tasolla havaintoja
       uustieto[ind]<-tieto[i]
       newradius[ind]<-radius[i]
       ind<-ind+1   
  }     
  i<-i+1  #siirytaan uudelle tasolle
}
lsets<-lsets[1:(ind-1),] 
uustieto<-t(t(uustieto[1:(ind-1)])) 
newradius<-t(t(newradius[1:(ind-1)])) 
return(list(levels=uustieto,lsets=lsets,radius=newradius,max=M))
}








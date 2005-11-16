kerest1d<-function(x,h,dendat){
#Laskee ydinestimaatorin arvon pisteissa x, sil.parametri h 
#kaytetaan normaalijakauman ydinfunktiota 
#x on lkm-vektori 
#dendat on n-vektori, data
#
lkm<-length(x)     #evaluoitavien pisteiden lkm
A<-matrix(0,lkm,1)
n<-length(dendat)  #datan maara
i<-1
while (i<=lkm){
  A[i]<-mean(knor1d(t(x[i]*matrix(1,1,n))-dendat,h))
       #muodostetaan matriisi, jossa x:n i:nnesta rivista vahennetty dendat
  i<-i+1
}
return(A)
}


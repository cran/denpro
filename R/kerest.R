kerest<-function(x,h,dendat){
#Laskee ydinestimaatorin arvon pisteessa x, sil.parametri h 
#kaytetaan Epanechnikovin ydinfunktiota 
#x on lkm*d-matriisi (esim data-matriisi)
#dendat n*d-matiisi, datamatriisi
#
d<-length(x[1,])       #tai: length(dendat[1,])
lkm<-length(x[,1])     #evaluoitavien pisteiden lkm
A<-matrix(0,lkm,1)
n<-length(dendat[,1])  #datan maara
i<-1
while (i<=lkm){
  A[i]<-mean(knor(t(x[i,]*matrix(1,d,n))-dendat,h))
       #muodostetaan matriisi, jossa x:n i:nnesta rivista vahennetty dendat
  i<-i+1
}
return(A)
}


knor<-function(x,h){
#Laskee normaalijakauman N(0,h^2) tih.funktion arvon pisteissa x.
#x on lkm*d-matriisi, tulos on lkm-vektori
#
if (dim(t(x))[1]==1) lkm<-1 else lkm<-length(x[,1])
if (lkm==1) d<-length(x) else d<-length(x[1,])
eta<-pituus(x)     #vektorin x pituuden nelio
normvakio<-(sqrt(2*pi)*h)^{-d}
tulos<-exp(-eta/h^2/2)*normvakio
return(tulos)
}


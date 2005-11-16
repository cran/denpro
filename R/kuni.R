kuni<-function(x,h){
#Laskee normaalijakauman N(0,h^2) tih.funktion arvon pisteessa x.
#d=2
#
#x on d-vector, tulos on real number
#
d<-2
normvakio<-h^{-d}
kuuluu<-TRUE
for (i in 1:d){
  if ((x[i]>h/2) || (x[i]<-h/2)) kuuluu<-FALSE
}
if (kuuluu) tulos<-normvakio
else tulos<-0
#
return(tulos)
}


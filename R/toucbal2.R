toucbal2<-function(dismat,lsets,radius,alkublokki,blokki){
#Finds which atoms touch each other, in each lev.set
#
#levset is len-vector, 1:s and 0:s, (m 1:s)
#dismat is n*n-matrix
#radius is n-vector
#
#Returns links
#
curblokki<-alkublokki
if (dim(t(lsets))[1]==1) levnum<-1 else levnum<-length(lsets[,1]) 
atomnum<-length(lsets[1,])
#
totnum<-atomnum*levnum
links<-matrix(NA,totnum,curblokki)
#
# merkitaan kosketukset links-matriisiin
ma<-0  #max of touches ovel levels and atoms
for (l in 1:levnum){
    alku<-atomnum*(l-1)
    maara<-matrix(0,atomnum,1)
    for (i in 1:atomnum){
      rivi<-as.matrix(dismat)[i,]
      j<-i+1
      while (j<=atomnum){
        isade<-radius[i,l]
        jsade<-radius[j,l]         
        eta<-rivi[j]
        crit<-isade+jsade 
        if (eta<=crit){ #jos pallot kosket
                 #pallo1:n keskipiste on dendat[obsi,], sade=isade
                 #pallo2:n keskipiste on dendat[obsj,], sade=jsade
            maari<-maara[i]+1
            maarj<-maara[j]+1
            if ((maari>curblokki) || (maarj>curblokki)){
                blokita2(links,blokki)
                curblokki<-curblokki+blokki   
            }
            links[alku+i,maari]<-j
            maara[i]<-maari
            links[alku+j,maarj]<-i
            maara[j]<-maarj         
        }
        j<-j+1 
      }
      i<-i+1
    }
    maapu<-max(maara)
    ma<-max(ma,maapu)
} 
links<-links[,1:ma]
return(links)
}












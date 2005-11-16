keresuni2d<-function(x,y,h,dendat){
#Laskee ydinestimaatorin arvon pisteessa (x,y), sil.parametri h 
#kaytetaan Gaussin ydinfunktiota 
#x on lkm*d-matriisi (esim data-matriisi)
#dendat n*d-matiisi, datamatriisi
#
d<-2                    #length(dendat[1,])  
n<-length(dendat[,1])   #datan maara
z<-mean(kuni(t(c(x,y)*matrix(1,2,n))-dendat,h))
return(z)
}


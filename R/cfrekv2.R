cfrekv2<-function(levels){
#laskee tasojoukon osien frekvenssit
#
#lsets on tasolkm*atomnum-matrix, 1: and 0:s
#
tasolkm<-length(levels[,1])     #levels:n rivien maara
frek<-matrix(0,tasolkm,1)
a<-1
while (a<=tasolkm){
   frek[a]<-sum(levels[a,])
   a<-a+1 
}
return(t(frek))
}








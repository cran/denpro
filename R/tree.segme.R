tree.segme<-function(tt,paletti=seq(1:2000),pcf=NULL)
{

colors<-colobary(tt$parent,paletti)
if (is.null(pcf)) segme<-colors
else{
  lenni<-length(pcf$value)
  segme<-matrix(0,lenni,1)
}
for (i in 1:length(colors)) segme[tt$infopointer[i]]<-colors[i]

return(segme)
}







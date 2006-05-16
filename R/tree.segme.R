tree.segme<-function(tt,paletti=seq(1:2000))
{

colors<-colobary(tt$parent,paletti)
segme<-colors
for (i in 1:length(colors)) segme[tt$infopointer[i]]<-colors[i]

return(segme)
}







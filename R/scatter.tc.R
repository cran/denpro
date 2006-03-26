scatter.tc<-function(dendat,tt,paletti=NULL,xlim=NULL,ylim=NULL)
{
if (is.null(xlim)) xlim<-c(min(dendat[,1]),max(dendat[,1]))
if (is.null(ylim)) ylim<-c(min(dendat[,2]),max(dendat[,2]))

if (is.null(paletti))
paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])

col<-colobary(tt$parent,paletti)

plot(dendat[tt$infopoint,],col=col,
xlim=xlim,ylim=ylim,
xlab="coordinate 1",ylab="coordinate 2")

}




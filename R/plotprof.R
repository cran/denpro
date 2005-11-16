plotprof<-function(profile,length=NULL,
plot=T,data=F,crit=NULL,orderrule="distcenter",
modelabel=TRUE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,
col="black",col.axis="black",
cutlev=NULL,xaxt="n",exmavisu=NULL)
{

#xaxs="e"    (extended)  not implemented?  xaxt="n"

parents<-profile$parent
levels<-profile$level
if (is.null(length)) length<-profile$volume
center<-profile$center

mut<-multitree(parents)
if (is.null(profile$roots)) roots<-mut$roots else roots<-profile$roots
child<-mut$child
sibling<-mut$sibling

d<-dim(center)[1]
if (is.null(crit)){
   crit<-rep(0,d)          #order so that 1st closest to origo
   if (d==1) crit<-max(center)
   if (!is.null(profile$refe)) crit<-profile$refe
}

if (orderrule=="distcenter") sibord<-siborder(mut,crit,profile$distcenter)
else sibord<-siborder(mut,crit,center)

itemnum<-length(parents)
vecs<-matrix(NA,itemnum,4)
vecs<-alloroot(vecs,roots,sibord,levels,length)
vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
orivecs<-vecs

if (!(is.null(cutlev))){
  cm<-cutmut(mut,cutlev,levels)              # cut the tree
  roots<-cm$roots
  sibling<-cm$sibling
  mut$roots<-roots
  if (orderrule=="distcenter") sibord<-siborder(mut,crit,profile$distcenter)
  else sibord<-siborder(mut,crit,center)
  rootnum<-length(roots) 
  apuvecs<-matrix(NA,itemnum,4)
  for (i in 1:rootnum){
     inde<-roots[i]
     apuvecs[inde,]<-vecs[inde,]
  }
  vecs<-apuvecs          #we give for the roots the previous positions
  vecs<-plotdata(roots,child,sibling,sibord,levels,length,vecs)
}

if (plot==T){
   if (!(is.null(cutlev))){
     xlim<-c(omamin(vecs[,1])-xmarginleft,omamax(vecs[,3])+xmarginright)
     ylim<-c(omamin(vecs[,2]),omamax(vecs[,2])+ptext+ymargin)
   }
   else{
     xlim<-c(omamin(vecs[,1])-xmarginleft,omamax(vecs[,3])+xmarginright)
     if (is.null(ylim)) ylim<-c(0,omamax(vecs[,2])+ptext+ymargin)
   }
   plotvecs(vecs,segme=T,xlim=xlim,ylim=ylim,xaxt=xaxt,
   col=col,col.axis=col.axis)
   # use original vectors (numbering will be correct)
   if (modelabel){
      plottext(parents,orivecs,ptext,leimat,symbo)  
   }
   if (!is.null(info)){
      plotinfo(vecs,info,pos=infopos,adj=NULL,lift=infolift,digits=3)
   }
}
#
#
if (data==T){
 return(list(sibord=t(sibord),vecs=vecs,parents=parents,levels=levels,
 length=length,center=center,remain=NULL))
}

###########################################################

if (!is.null(exmavisu)){

node<-exmavisu

x1<-xcoor[2*node-1] 
x2<-xcoor[2*node]
lev<-levels[node]
if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-0
polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col="blue")

pino<-matrix(0,nodenum,1)
pino[1]<-child[node]
if (child[node]>0) pinoin<-1 else pinoin<-0

while (pinoin>0){
   node<-pino[pinoin]
   pinoin<-pinoin-1   

   x1<-xcoor[2*node-1] 
   x2<-xcoor[2*node]
   lev<-levels[node]
   if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-0
   polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col="blue")

   if (sibling[node]>0){
         pinoin<-pinoin+1
         pino[pinoin]<-sibling[node] 
   }

   while (child[node]>0){    #go to left and put right nodes to stack
         node<-child[node]

         x1<-xcoor[2*node-1] 
         x2<-xcoor[2*node]
         lev<-levels[node]
         if (parents[node]>0) lev0<-levels[parents[node]] else lev0<-0
         polygon(c(x1,x2,x2,x1),c(lev0,lev0,lev,lev),col="blue")

         if (sibling[node]>0){
            pinoin<-pinoin+1
            pino[pinoin]<-sibling[node] 
         }
   }
}
}

##########################################

}











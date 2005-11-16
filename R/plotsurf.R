plotsurf<-function(lstree,coordi,
plot=T,data=F,crit=NULL,
modelabel=TRUE,ptext=0,leimat=NULL,symbo=NULL,
info=NULL,infolift=0,infopos=0,
xmarginleft=0,xmarginright=0,ymargin=0,
xlim=NULL,ylim=NULL,
nodesymbo=20,col="black",col.axis="black",collines="black"){

parent<-lstree$parent
center<-lstree$center
level<-lstree$level
radi<-lstree$radi

nodenum<-length(parent)
xcoordinate<-center[coordi,]
lowe<-radi[,coordi]
d<-length(radi[1,])/2
uppi<-radi[,(d+coordi)]

xlim<-c(min(lowe)-xmarginleft,max(uppi)+xmarginright)
ylim<-c(0,max(level)+ptext+ymargin)

plot(c(lowe,xcoordinate,uppi),c(level,level,level),
xlab="",ylab="",xlim=xlim,ylim=ylim,
pch=nodesymbo,col=c(col,col,col),col.axis=col.axis) 

for (i in 1:nodenum){
    if (parent[i]>0){
        xchild<-xcoordinate[i]
        ychild<-level[i]
        xparent<-xcoordinate[parent[i]]
        yparent<-level[parent[i]]
        if (length(collines)>1) colli<-collines[i]
        else colli<-collines
        segments(xparent,yparent,xchild,ychild,col=colli) 
     }
}

for (i in 1:nodenum){
    if (parent[i]>0){
        xchild<-lowe[i]
        ychild<-level[i]
        xparent<-lowe[parent[i]]
        yparent<-level[parent[i]]
        if (length(collines)>1) colli<-collines[i]
        else colli<-collines
        segments(xparent,yparent,xchild,ychild,col=colli) 
     }
}

for (i in 1:nodenum){
    if (parent[i]>0){
        xchild<-uppi[i]
        ychild<-level[i]
        xparent<-uppi[parent[i]]
        yparent<-level[parent[i]]
        if (length(collines)>1) colli<-collines[i]
        else colli<-collines
        segments(xparent,yparent,xchild,ychild,col=colli) 
     }
}

#########
#########
if (modelabel){

data<-plotprof(lstree,plot=F,data=T,cutlev=NULL,ptext=NULL,info=NULL,
infolift=0,infopos=0)
vecs<-data$vecs
mlkm<-moodilkm(parent)
modloc<-mlkm$modloc

nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)

for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}                     
moodinum<-length(modloc)
modelocx<-matrix(0,moodinum,1)
modelocy<-matrix(0,moodinum,1)
if (is.null(leimat)){
   if (is.null(symbo)){
       labels<-paste("M",1:moodinum,sep="")
   }
   else{
         if (symbo=="empty") labels<-paste("",1:moodinum,sep="")
         else  labels<-paste(symbo,1:moodinum,sep="")
   }
}
else{
   labels<-leimat
}
xcor<-matrix(0,moodinum,1)
for (i in 1:moodinum){
    loc<-modloc[i]
    xcor[i]<-xcoor[2*loc-1]
}
modloc<-omaord2(modloc,xcor)
for (i in 1:moodinum){
    loc<-modloc[i]
    modelocx[i]<-xcoordinate[loc]
    modelocy[i]<-level[loc]+ptext
}
text(modelocx,modelocy,labels)         
}
############################################
}












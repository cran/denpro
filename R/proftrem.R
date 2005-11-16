proftrem<-function(tr,N,Q,suppo,cvol=TRUE,ccen=TRUE,cfre=FALSE,
compoinfo=FALSE,katka=NULL,surf=FALSE)
{
d<-length(N)
h<-0  

nodenumOfDyaker<-length(tr$left)

maxval<-max(tr$value)
minval<-min(tr$value)
step<-(maxval-minval)/Q
levseq<-seq(from=minval,to=maxval-step,by=step)

levfrekv<-matrix(0,Q,1)
atomnum<-length(tr$value)
for (i in 1:atomnum){
   for (j in 1:Q){
       if (tr$value[i]>=levseq[j]){
          levfrekv[j]<-levfrekv[j]+1
       }
   }
}
numofall<-sum(levfrekv)
levnum<-length(levseq)

inlevseq<-matrix(0,length(levseq)+1,1)
inlevseq[2:(length(levseq)+1)]<-levseq
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N
inleft<-matrix(0,length(tr$left)+1,1)
inleft[2:(length(tr$left)+1)]<-tr$left
inright<-matrix(0,length(tr$left)+1,1)
inright[2:(length(tr$left)+1)]<-tr$right
inparent<-matrix(0,length(tr$left)+1,1)
inparent[2:(length(tr$left)+1)]<-tr$parent
invalue<-matrix(0,length(tr$value)+1,1)
invalue[2:(length(tr$value)+1)]<-tr$value
#inindex<-matrix(0,dim(kg$index)[1]+1,dim(kg$index)[2]+1)
#for (i in 1:dim(kg$index)[1]){
#  inindex[i+1,]<-c(0,kg$index[i,])
#}
innodefinder<-matrix(0,length(tr$nodefinder)+1,1)
innodefinder[2:(length(tr$nodefinder)+1)]<-tr$nodefinder

# Tama koodi on jo kergrid:ssa, lasketaan volume of one atom

minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
     minim[i]=suppo[2*i-1];
     maxim[i]=suppo[2*i];
}
 
delta<-(maxim-minim)/(N+1)
volofatom<-prod(delta)

inminim<-matrix(0,d+1,1)
inminim[2:(d+1)]<-minim
indelta<-matrix(0,d+1,1)
indelta[2:(d+1)]<-delta

if (!is.null(katka)){
   invalue2<-invalue
   lenni<-length(invalue)
   for (i in 1:lenni){
      if (invalue[i]>=katka) invalue2[i]<-katka
   }
   invalue<-invalue2
}

dentree<-.C("decomdyaC",
               as.integer(numofall),
               as.integer(atomnum),
               as.double(inlevseq),
               as.integer(inN),
               as.integer(d),
               as.integer(levnum),
               as.double(volofatom),
               as.double(inminim),
               as.double(h),
               as.double(indelta),
               as.integer(nodenumOfDyaker),
               as.integer(inleft),
               as.integer(inright),
               as.integer(inparent),
               as.double(invalue),
               as.integer(tr$index),
               as.integer(innodefinder),
               level = double(numofall+1),
               parent = integer(numofall+1),
               component = integer(numofall+1),
               volume = double(numofall+1),
               center = double(d*numofall+1),
               efek = integer(1),
               AtomlistAtom = integer(numofall+1),
               AtomlistNext = integer(numofall+1),
               numOfAtoms = integer(1),
PACKAGE="denpro")

AtomlistAtom<-dentree$AtomlistAtom[2:(dentree$numOfAtoms+1)]
AtomlistNext<-dentree$AtomlistNext[2:(dentree$numOfAtoms+1)]

invalue<-dentree$level[2:(dentree$efek+1)]
parent<-dentree$parent[2:(dentree$efek+1)]
volume<-dentree$volume[2:(dentree$efek+1)]
component<-dentree$component[2:(dentree$efek+1)]
kerroin<-cinte(invalue,volume,parent)
sepvalnor<-invalue/kerroin
veccenter<-dentree$center[2:(d*dentree$efek+1)]
center<-matrix(0,dentree$efek,d)
for (i in 1:dentree$efek){
  for (j in 1:d){
     center[i,j]<-veccenter[(i-1)*d+j]
  }
}
center<-t(center)

if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL

clus<-F
if (clus){
   clustervecs<-cluskern(compo,component,AtomlistAtom,AtomlistNext,kg,dendat,
   h,N)
}
else{
   clustervecs<-NULL
}

if (surf==T)
 radi<-csurface(component,AtomlistNext,AtomlistAtom,
       minim,h,delta,tr$index,d,center)
else radi<-NULL

if (compoinfo)

  return(list(parent=parent,level=sepvalnor,invalue=invalue,
  volume=volume,center=center,
  component=component,
  AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext,index=index))

else

  return(list(parent=parent,level=sepvalnor,invalue=invalue,
  volume=volume,center=center,radi=radi))

#,nodefrek=nodefrek,clustervec=clustervecs))
#
#values: normeeratut arvot
#invalues: normeeraamattomat arvot
#nodefrek: kunkin solmun frekvenssi

}

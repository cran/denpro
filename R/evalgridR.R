evalgridR<-function(M,sig,p,suppo,N){

mixnum<-length(p)
d<-dim(M)[2]

depth<-log(N,base=2)   
bigd<-sum(depth)
depoftree<-bigd+1
maxpositive<-prod(N)
maxnode<-bigd*maxpositive

minim<-matrix(0,d,1)   #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-suppo[2*i-1]
  maxim[i]<-suppo[2*i]
}
delta<-(maxim-minim)/(N+1)

numnode<-1
left<-matrix(0,maxnode,1)
right<-matrix(0,maxnode,1)
parent<-matrix(0,maxnode,1)
infopointer<-matrix(0,maxnode,1)
low<-matrix(0,maxnode,1)
low[1]<-1
upp<-matrix(0,maxnode,1)
upp[1]<-N[1]

numpositive<-0
value<-matrix(0,maxpositive,1)
index<-matrix(0,maxpositive,d)
nodefinder<-matrix(0,maxpositive,1)

for (i in 1:maxpositive){
     if (d>1){  
          inde<-digit(i,N)   #inde is d-vector
     }
     else{
         inde<-i
     }

     point<-minim+delta*inde     #point is d-vector  

     ##########################
     #val<-epane(point,1)
     val<-0
     for (mi in 1:mixnum){
        val<-val+p[mi]*evanor(point-M[mi,]/sig[mi,])/prod(sig[mi,])
     }
     ##########################

     fe<-findend(inde,left,right,low,upp,N)
     curre<-fe$location
     curdep<-fe$dep

     ad<-addnode(inde,curre,curdep,left,right,parent,low,upp,N,numnode)
     numnode<-ad$numnode
     left<-ad$left
     right<-ad$right
     parent<-ad$parent
     low<-ad$low
     upp<-ad$upp
     nodeloc<-ad$nodeloc

     numpositive<-numpositive+1
     infopointer[numnode]<-numpositive
     value[numpositive]<-val
     index[numpositive,]<-inde
     nodefinder[numpositive]<-nodeloc
}

left<-left[1:numnode]
right<-right[1:numnode]
parent<-parent[1:numnode]
infopointer<-infopointer[1:numnode]
low<-low[1:numnode]
upp<-upp[1:numnode]

value<-value[1:numpositive,]
index<-index[1:numpositive,]
nodefinder<-nodefinder[1:numpositive]

return(list(left=left,right=right,parent=parent,infopointer=infopointer,
low=low,upp=upp,value=value,index=index,nodefinder=nodefinder))
}                              









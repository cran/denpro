streedisc<-function(tt,et,r)
{
# r is vector of radiuses, we prune shapetree "tt" so that
# its radiuses are given by r

mt<-multitree(tt$parent)
child<-mt$child
sibling<-mt$sibling

d<-dim(tt$center)[1]
itemnum<-length(tt$parent)

parent<-matrix(0,itemnum,1)
#child<-matrix(0,itemnum,1)
#sibling<-matrix(0,itemnum,1)

pino<-matrix(0,itemnum,1)
pinoparent<-matrix(0,itemnum,1)
pinorad<-matrix(0,itemnum,1)

pino[1]<-1
pinoparent[1]<-1
pinorad[1]<-1
pinin<-1
curradind<-1

while (pinin>0){ # && (curradind<=length(r))){
      cur<-pino[pinin]      #take from stack
      curpar<-pinoparent[pinin]
      curradind<-pinorad[pinin]
      pinin<-pinin-1

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
            pinoparent[pinin]<-curpar
            pinorad[pinin]<-curradind
      }

      recci<-matrix(0,2*d,1)
      note<-tt$infopointer[cur] 
      for (jj in 1:d){
         recci[2*jj-1]<-et$support[2*jj-1]+et$step[jj]*et$low[note,jj]
         recci[2*jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[note,jj]  
      }
      etai<-sqrt(etaisrec(tt$bary,recci))

      if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
      if (etai>currad){
          parent[cur]<-curpar
          curpar<-cur
          curradind<-curradind+1
      }

      # go to left and put right nodes to the stack
      while (child[cur]>0){  # && (curradind<=length(r))){
            cur<-child[cur]

            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
                 pinoparent[pinin]<-curpar
                 pinorad[pinin]<-curradind
            }
 
            recci<-matrix(0,2*d,1)
            note<-tt$infopointer[cur] 
            for (jj in 1:d){
                recci[2*jj-1]<-et$support[2*jj-1]+et$step[jj]*et$low[note,jj]
                recci[2*jj]<-et$support[2*jj-1]+et$step[jj]*et$upp[note,jj]  
            }
            etai<-sqrt(etaisrec(tt$bary,recci))

            if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
            if (etai>currad){
                parent[cur]<-curpar
                curpar<-cur
                curradind<-curradind+1 
            }
 
      }

}

tt$roots<-c(1)
tt$parent<-parent

#return(tt)

# Prune ##################################

newparent<-matrix(0,itemnum,1)
newcenter<-matrix(0,d,itemnum)
newvolume<-matrix(0,itemnum,1)
newlevel<-matrix(0,itemnum,1)
newpointer<-matrix(0,itemnum,1)

newparent[1]<-0
newcenter[,1]<-tt$center[,1]
newvolume[1]<-tt$volume[1]
newlevel[1]<-tt$level[1]
newpointer[1]<-1

i<-2
newlkm<-1
while (i<=itemnum){
    if (tt$parent[i]!=0){
         newlkm<-newlkm+1
         newpointer[i]<-newlkm
         newparent[newlkm]<-newpointer[tt$parent[i]]
         newcenter[,newlkm]<-tt$center[,i]
         newlevel[newlkm]<-tt$level[i]
         newvolume[newlkm]<-tt$volume[i]

    }
    i<-i+1
}

newparent<-newparent[1:newlkm]
newcenter<-newcenter[,1:newlkm]
newvolume<-newvolume[1:newlkm]
newlevel<-newlevel[1:newlkm]

return(list(parent=newparent,level=newlevel,volume=newvolume,center=newcenter,
root=1))

}   


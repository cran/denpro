modeparent<-function(lst,modepoint,newnum)
{
# lst - level set tree
# modepoint - pointers to those nodes of lst which are modes
# newnum - a number, typically smaller or equal to the number of nodes of lst
# we reduce the number of modes, returning vector "oldtonew"
# which gives pointers from "modepoint" to itself,
# the locations of removed modes point to the mode which is left over
# and to which this mode is merged
# note that a pointer in "oldtonew" may point to a node which does
# not exist but points further
# "newnode" may contain branching nodes instead of modes

modenum<-length(modepoint)
oldtonew<-matrix(0,modenum,1)
newnode<-modepoint

ex<-excmas(lst)
exvec<-matrix(0,modenum,1)
for (i in 1:modenum) exvec[i]<-ex[modepoint[i]]

curnum<-modenum
while (curnum>newnum){
      inni<-omaind(exvec)
      exvec[inni]<-NA
      neighbor<-findneighbor(lst,modepoint[inni])
      for (i in 1:modenum){
           if (modepoint[i]==neighbor){ 
                oldtonew[inni]<-i
                exvec[i]<-ex[neighbor]
                newnode[i]<-neighbor
           }
      }
      curnum<-curnum-1
}
return(list(oldtonew=oldtonew,newnode=newnode))
}









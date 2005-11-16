coloallof<-function(mt,snum,paletti,d)
# fast allocation of colors (matching of modes)
{
#mt is mode tree
#paletti gives a list of colors

xcoor<-mt$xcoor
ycoor<-mt$ycoor
mlabel<-mt$mlabel
lenni<-length(ycoor)

colot<-matrix("",lenni,1)

# find the locations for the information for each h

low<-matrix(0,snum,1)
upp<-matrix(0,snum,1)
low[1]<-1
glob<-2
while ((glob<=lenni) && (mlabel[glob]!=1)){
       glob<-glob+1
}
upp[1]<-glob-1
# now glob is at the start of new block
i<-2
while (i<=snum){
   low[i]<-glob
   glob<-glob+1
   while ((glob<=lenni) && (mlabel[glob]!=1)){
       glob<-glob+1
   }
   upp[i]<-glob-1
   i<-i+1
}

# first we allocate colors for the largest h

run<-1  #low[1]
while (run<=upp[1]){
   colot[run]<-paletti[run]
   run<-run+1
}

firstnewcolo<-run

i<-2

while (i<=snum){
   prenum<-upp[i-1]-low[i-1]+1
   curnum<-upp[i]-low[i]+1

   if (prenum<=curnum){  

      compa<-i-1
 
      dista<-matrix(NA,prenum,curnum)
      for (ap in low[i]:upp[i]){
         for (be in low[compa]:upp[compa]){
           if (d==1){
               curcenter<-xcoor[ap]
               precenter<-xcoor[be]
           }
           else{
               curcenter<-xcoor[ap,]
               precenter<-xcoor[be,]
           }
           dista[be-low[compa]+1,ap-low[i]+1]<-etais(curcenter,precenter)
         }
      }

      run<-1
      while (run<=prenum){
          minimi<-min(dista,na.rm=TRUE)
          argmin<-which(minimi==dista)[1]   
 
          yind<-ceiling(argmin/prenum)      
          xind<-argmin-(yind-1)*prenum

          dista[xind,]<-NA
          dista[,yind]<-NA

          colot[low[i]+yind-1]<-colot[low[compa]+xind-1]    
          run<-run+1
      }

      run<-low[i]
      while (run<=upp[i]){
          if (colot[run]==""){
             colot[run]<-paletti[firstnewcolo]
             firstnewcolo<-firstnewcolo+1   
          }
          run<-run+1
      }


   }
   else{   #curnum>prenum: fewer modes
      
      compa<-i-1
 
      dista<-matrix(NA,prenum,curnum)
      for (ap in low[i]:upp[i]){
         for (be in low[compa]:upp[compa]){
           if (d==1){
               curcenter<-xcoor[ap]
               precenter<-xcoor[be]
           }
           else{
               curcenter<-xcoor[ap,]
               precenter<-xcoor[be,]
           }
           dista[be-low[compa]+1,ap-low[i]+1]<-etais(curcenter,precenter)
         }
      }

      run<-1
      while (run<=curnum){
          minimi<-min(dista,na.rm=TRUE)
          argmin<-which(minimi==dista)[1]   
 
          yind<-ceiling(argmin/prenum)      
          xind<-argmin-(yind-1)*prenum

          dista[xind,]<-NA
          dista[,yind]<-NA

          colot[low[i]+yind-1]<-colot[low[compa]+xind-1]
          run<-run+1
      }

   }
  
   i<-i+1
}

return(colot)
}



















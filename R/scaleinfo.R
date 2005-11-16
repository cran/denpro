scaleinfo<-function(tseq)
{

apu<-split.screen(c(1,2))                  # split display into two screens
apu2<-split.screen(c(2,1), screen = 2)     # now split the right half into 2

#nf<-layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
#layout.show(nf)

pr<-tseq$treelist[[1]]

lenni<-length(tseq$alfa)
printvec<-tseq$alfa[lenni:1]
mt<-modegraph(tseq$treelist,printvec)
paletti<-c("red","blue","green","turquoise","orange","navy",
"darkgreen","orchid",colors()[50:100])
ca<-coloallo(mt,paletti)

ycor<-0.01
while (ycor>0){

# Volume plot and barycenter plot

screen(apu2[1])
plotasso(pr)  #,bg="wheat") #,cutlev=5000)  

screen(apu2[2])
#par(bg="transparent")
plotbary(pr,coordi=1)

# Mode tree

screen(apu[1])   # prepare screen 1 for out
plotmodet(coordi=1,mt,colot=ca,log="y",shift=0.001) 

# Interaction

   #erase.screen(apu2[1])
   #screen(apu[1])
   loc<-locator(1)
   
   #step<-ceiling(5000/lenni)
   #grid<-seq(0,5000+step,by=step)

   ylamidi<-(tseq$alfa[1]+tseq$alfa[1+1])/2
   if (loc$y<ylamidi) indeksi<-1
   for (i in 2:(lenni-1)){
      alamidi<-(tseq$alfa[i-1]+tseq$alfa[i])/2
      ylamidi<-(tseq$alfa[i]+tseq$alfa[i+1])/2

      if ((loc$y>=alamidi) && (loc$y<ylamidi)) indeksi<-i
   }
   alamidi<-(tseq$alfa[lenni-1]+tseq$alfa[lenni])/2
   if (loc$y>=alamidi) indeksi<-lenni

   #screen(apu2[1]) #, new=TRUE)

   pr<-tseq$treelist[[indeksi]]
   #plotasso(pr,bg="wheat") #,cutlev=5000)  
  
   ycor<-loc$y 

}

close.screen(all = TRUE)

}


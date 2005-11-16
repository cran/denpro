scaleplot<-function(estiseq,paletti=NULL,shift=0,ptext=0,bm=NULL,levnum=NULL)
{
hseq<-estiseq$hseq

hnum<-length(hseq)
d<-dim(estiseq$lstseq[[hnum]]$center)[1]

if ((hnum>1) && (hseq[1]<hseq[2])){  
    hseq<-hseq[seq(hnum,1)]
    apuseq<-list(estiseq$lstseq[[hnum]])
    i<-2
    while (i <= hnum){
         apuseq<-c(apuseq,list(estiseq$lstseq[[hnum-i+1]]))
         i<-i+1 
   }
   estiseq$lstseq<-apuseq
}

if (is.null(paletti))
paletti<-c("red","blue","green","turquoise","orange","navy",
"darkgreen","orchid",colors()[50:100])

if (!is.null(levnum)){
   for (i in 1:hnum){
      lf<-treedisc(estiseq$lstseq[[i]],estiseq$pcfseq[[i]],levnum) 

      if (i==1){
           if (hnum==1){
               reduseq<-lf
           }
           else{
               reduseq<-list(lf)
           }
      }
      else{
          reduseq<-c(reduseq,list(lf))
      }
  }
  estiseq$lstseq<-reduseq
}

mt<-modegraph(estiseq)

#if (is.null(bm)) bm<-branchmap(estiseq)

pr<-estiseq$lstseq[[1]]
pcf<-estiseq$pcfseq[[1]]

x11(width=4,height=4)
plotvolu(pr,ptext=ptext)

x11(width=4,height=4)
icolo<-mt$colot[mt$low[1]:mt$upp[1]]
inodes<-mt$nodepointer[mt$low[1]:mt$upp[1]]
modlab<-plotbary(pr,coordi=1,ptext=ptext,
        modlabret=T,modecolo=icolo,modepointer=inodes)

#x11(width=6,height=8)
#persp(x=bm$level,y=bm$h,z=bm$z, xlab="level",ylab="h",zlab="",ticktype="detailed",
#col=bm$col,phi=40,theta=10)
#coordi<-0

x11(width=6,height=8)
coordi<-1
plotmodet(mt,coordi=coordi)
modelocx<-modlab$modelocat[,coordi]+shift
modelocy<-hseq[1]
labels<-modlab$labels
text(modelocx,modelocy,labels)
title(sub=paste("coordinate",as.character(coordi)))

loc<-locator(1)
ycor<-loc$y 
alaraja<-hseq[hnum]

while (ycor>=alaraja){

   alamidi<-(hseq[1]+hseq[1+1])/2
   if (ycor>=alamidi) indeksi<-1
   for (i in 2:(hnum-1)){
      alamidi<-(hseq[i]+hseq[i+1])/2
      ylamidi<-(hseq[i-1]+hseq[i])/2

      if ((ycor>=alamidi) && (ycor<ylamidi)) indeksi<-i
   }
   ylamidi<-(hseq[hnum-1]+hseq[hnum])/2
   if (ycor<ylamidi) indeksi<-hnum

   pr<-estiseq$lstseq[[indeksi]]
   pcf<-estiseq$pcfseq[[indeksi]]

   # Volume plot and barycenter plot
   
   dev.set(which = 2) #dev.next())
   plotvolu(pr,ptext=ptext)
   
   alaasso<-0
   ylaasso<-max(pr$level)
   
   loc<-locator(1)
   alax<-0
   ylax<-pr$volume[1]
   while (loc$y>=alaasso){

       if (loc$y>ylaasso)  plotvolu(pr)
       else{
           if ((ylax-loc$x) <= (loc$x-alax)) ylax<-loc$x
           else                              alax<-loc$x
           plotvolu(pr,xlim=c(alax,ylax),ptext=ptext)
       }
       loc<-locator(1)
   }

   dev.set(which = 3) #dev.next())
   coordi<-1
   icolo<-mt$colot[mt$low[indeksi]:mt$upp[indeksi]]
   inodes<-mt$nodepointer[mt$low[indeksi]:mt$upp[indeksi]]
   modlab<-plotbary(pr,coordi=coordi,ptext=ptext,
                    modlabret=T,modecolo=icolo,modepointer=inodes)
   title(sub=paste("coordinate",as.character(coordi)))

   loc<-locator(1)
   while (loc$y>=alaasso){
       if (coordi<=(d-1)) coordi<-coordi+1 else coordi<-1
       plotbary(pr,coordi=coordi,ptext=ptext,modecolo=icolo,modepointer=inodes)
       title(sub=paste("coordinate",as.character(coordi)))

       loc<-locator(1)
   }

   # Interaction
   
   dev.set(which = 4) #dev.next())
   loc<-locator(1)
   locus<-loc

   coordi<-1
   ylamodet<-hseq[1]
   while (loc$y>=ylamodet){
       if (coordi<=(d-1)) coordi<-coordi+1 else coordi<-1
       #if (coordi==0){
       #persp(x=bm$level,y=bm$h,z=bm$z, xlab="level",ylab="h",zlab="",ticktype="detailed",
       #  col=bm$col,phi=40,theta=10)
       #}
       #else{
           plotmodet(mt,coordi=coordi)
           modelocx<-modlab$modelocat[,coordi]+shift
           modelocy<-hseq[indeksi]
           labels<-modlab$labels
           text(modelocx,modelocy,labels)
           title(sub=paste("coordinate",as.character(coordi)))
       #}

       loc<-locator(1)
   }

   ycor<-locus$y 

}

dev.off()
dev.off()
dev.off()

}


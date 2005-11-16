scaletable<-function(estiseq,paletti=NULL,shift=0,ptext=0,ptextst=0,
bm=NULL,#mt=NULL,
levnum=60,levnumst=60,redu=T)
{
# preparation
if ((length(estiseq$hseq)>1) && (estiseq$hseq[1]<estiseq$hseq[2])){  
    hnum<-length(estiseq$hseq)
    estiseq$hseq<-estiseq$hseq[seq(hnum,1)]
    apuseq<-list(estiseq$lstseq[[hnum]])
    i<-2
    while (i <= hnum){
         apuseq<-c(apuseq,list(estiseq$lstseq[[hnum-i+1]]))
         i<-i+1 
   }
   estiseq$lstseq<-apuseq
}

if (estiseq$type=="carthisto")  smootseq<--estiseq$leaf
else smootseq<-estiseq$hseq
hnum<-length(smootseq)
d<-dim(estiseq$lstseq[[hnum]]$center)[1]

if (estiseq$type=="carthisto") redu<-F
if (is.null(estiseq$stseq)) levnumst<-NULL

if (is.null(paletti))
paletti<-c("red","blue","green","turquoise","orange","navy",
"darkgreen","orchid",colors()[50:100])

# prune the level set trees
if ((!is.null(levnum)) && (redu)){
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

# prune the shape trees
if ((!is.null(levnumst)) && (redu)){
   for (i in 1:hnum){
      lf<-treedisc(estiseq$stseq[[i]],estiseq$pcfseq[[i]],levnumst) 
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
}
else reduseq<-estiseq$stseq

#if (is.null(mt)) 
mt<-modegraph(estiseq)
if (is.null(bm)) bm<-branchmap(estiseq)

####################################
devicontrol<-2
devibranch<-3
devibary<-4
devimodet<-5
devivolu<-6
deviradi<-7
deviloca<-8

xmin<--0.5
xmax<-0.5
ymin<--1
ymax<-1
lkm<-7
step<-(ymax-ymin)/(lkm-1)
heig<-seq(ymin,ymax,step)
xloc<-0

# control window
x11(width=2,height=6)
plot(x="",y="",xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax))
text(xloc,heig[lkm],"Mode graph")
text(xloc,heig[lkm-1],"Map of branches")
text(xloc,heig[lkm-2],"Volume plot")
text(xloc,heig[lkm-3],"Barycenter plot")
text(xloc,heig[lkm-4],"Radius plot")
text(xloc,heig[lkm-5],"Location plot")
text(xloc,heig[lkm-6],"STOP")

devit<-matrix(0,lkm,1)
devit[1]<-devimodet
devit[2]<-devibranch
devit[3]<-devivolu
devit[4]<-devibary
devit[5]<-deviradi
devit[6]<-deviloca
devit[7]<-devicontrol

# choose estimate
indeksi<-1
pr<-estiseq$lstseq[[indeksi]]
pcf<-estiseq$pcfseq[[indeksi]]
if (!is.null(levnumst)){ 
      st<-estiseq$stseq[[indeksi]]
      stredu<-reduseq[[indeksi]]
}
hcur<-estiseq$hseq[indeksi]

# branch map
x11(width=4,height=5)
phi<-40
theta<-10
persp(x=bm$level,y=bm$h,z=bm$z, xlab="level",ylab="h",zlab="",
ticktype="detailed",col=bm$col,phi=phi,theta=theta)
title(sub="map of branches")

# barycenter plot
x11(width=3.5,height=4)
coordi<-1
icolo<-mt$colot[mt$low[1]:mt$upp[1]]
inodes<-mt$nodepointer[mt$low[1]:mt$upp[1]]
modlab<-plotbary(pr,coordi=coordi,ptext=ptext,
        modlabret=T,modecolo=icolo,modepointer=inodes)
title(sub=paste("barycenter plot, coordinate",as.character(coordi)))

# mode tree
x11(width=4,height=5)
coordi<-1
plotmodet(mt,coordi=coordi)
modelocx<-modlab$modelocat[,coordi]+shift
modelocy<-smootseq[1]
labels<-modlab$labels
text(modelocx,modelocy,labels)
title(sub=paste("mode graph, coordinate",as.character(coordi)))

# volume plot
x11(width=3.5,height=4)
plotvolu(pr,ptext=ptext)
title(sub=paste("volume plot, h=",as.character(round(hcur,digits=3))))

# radius plot
if (!is.null(levnumst)){
  refelab<-"moodi"
  st.moodi<-st

  lev<-0.1*max(pcf$value)  
  refe<-st$bary
  st.bary<-leafsfirst(pcf,lev=lev,refe=refe)

  x11(width=3,height=4)
  plotvolu(stredu,ptext=ptextst,symbo="T")
  title(sub=paste("radius plot, level=",as.character(round(lev,digits=3)),
  ",  mode=refe'nce point"))
}

# location plot
if (!is.null(levnumst)){
  x11(width=3,height=4)
  lcoordi<-1
  plotbary(stredu,coordi=lcoordi,ptext=ptextst,symbo="T")
  title(sub=paste("location plot, coordinate",as.character(lcoordi)))
}


##################################################################
ylow<-matrix(0,lkm,1)
yupp<-matrix(0,lkm,1)
for (i in 1:lkm){
    ylow[i]<-heig[i]-step/2
    yupp[i]<-heig[i]+step/2
}
ylow<-ylow[lkm:1]
yupp<-yupp[lkm:1]

dev.set(which = devicontrol)
loc<-locator(1)
while (loc$y>=yupp[lkm]){
  for (i in 1:lkm){
      if ((loc$y>ylow[i]) && (loc$y<=yupp[i])){
                 devi<-devit[i]
      }
  }
  dev.set(which = devi)
  loc<-locator(1)

  # interaction in modegraph
  if (devi==devimodet){
       alaraja<-smootseq[length(smootseq)]
       while (loc$y>=alaraja){
          coordi<-1
          ylamodet<-smootseq[1]
          while (loc$y>=ylamodet){
              if (coordi<=(d-1)) coordi<-coordi+1 else coordi<-1
              plotmodet(mt,coordi=coordi)
              modelocx<-modlab$modelocat[,coordi]+shift
              modelocy<-smootseq[indeksi]
              labels<-modlab$labels
              text(modelocx,modelocy,labels)
              title(sub=paste("mode graph, coordinate",as.character(coordi)))

              loc<-locator(1)
          }
          
          if (loc$y>=alaraja){
             alamidi<-(smootseq[1]+smootseq[1+1])/2
             if (loc$y>=alamidi) indeksi<-1
             for (i in 2:(hnum-1)){
                alamidi<-(smootseq[i]+smootseq[i+1])/2
                ylamidi<-(smootseq[i-1]+smootseq[i])/2

                if ((loc$y>=alamidi) && (loc$y<ylamidi)) indeksi<-i
             }
             ylamidi<-(smootseq[hnum-1]+smootseq[hnum])/2
             if (loc$y<ylamidi) indeksi<-hnum

             pr<-estiseq$lstseq[[indeksi]]
             pcf<-estiseq$pcfseq[[indeksi]]
             hcur<-estiseq$hseq[[indeksi]]
             if (!is.null(levnumst)){ 
                    lev<-0.1*max(pcf$value)
                    st<-estiseq$stseq[[indeksi]]
                    st.moodi<-st
                    st.bary<-NULL
                    stredu<-reduseq[[indeksi]]
             }

             dev.set(which = devivolu) 
             plotvolu(pr,ptext=ptext)
             title(sub=paste("volume plot, h=",as.character(round(hcur,digits=3))))
  
             dev.set(which = devibary) 
             coordi<-1
             icolo<-mt$colot[mt$low[indeksi]:mt$upp[indeksi]]
             inodes<-mt$nodepointer[mt$low[indeksi]:mt$upp[indeksi]]
             modlab<-plotbary(pr,coordi=coordi,ptext=ptext,
                          modlabret=T,modecolo=icolo,modepointer=inodes)
             title(sub=paste("barycenter plot, coordinate",as.character(coordi)))

             dev.set(which = deviradi) 
             plotvolu(stredu,ptext=ptextst,symbo="T")
             title(sub=paste("radius plot, level=",as.character(round(lev,digits=3)),
             ",  mode=refe'nce point"))
 
             dev.set(which = deviloca) 
             lcoordi<-1
             plotbary(stredu,coordi=lcoordi,ptext=ptextst,symbo="T")
             title(sub=paste("location plot, coordinate",as.character(lcoordi)))
             
             dev.set(which = devimodet)

             modelocx<-modlab$modelocat[,coordi]+shift
             modelocy<-smootseq[indeksi]
             labels<-modlab$labels
             text(modelocx,modelocy,labels)

             loc<-locator(1)
          }
       }
       #dev.set(which = devicontrol)
  }

  # interaction in volume plot
  if (devi==devivolu){
     alaasso<-0
     ylaasso<-max(pr$level)
     alax<-0
     ylax<-pr$volume[1]
     while (loc$y>=alaasso){
        if (loc$x>=0){
           if (loc$y>ylaasso) plotvolu(pr)
           else if (loc$x>0){
              keskip<-alax+(ylax-alax)/2
              if (loc$x >= keskip) ylax<-loc$x
              else                 alax<-loc$x
              plotvolu(pr,xlim=c(alax,ylax),ptext=ptext)
           }
           title(sub=paste("volume plot, h=",as.character(round(hcur,digits=3))))
        }     
        else if (!is.null(levnumst)){
             maksi<-max(pr$level)
             mode<-locofmax(pcf)
             lev<-min(max(0,loc$y),maksi)
             st<-leafsfirst(pcf,refe=mode,lev=lev)
             st.moodi<-st
             st.bary<-NULL
             if (redu) stredu<-treedisc(st,pcf,ngrid=levnumst) else stredu<-st
             refelab<-"moodi"

             dev.set(which = deviradi) 
             plotvolu(stredu,ptext=ptextst,symbo="T")
             title(sub=paste("radius plot, level=",as.character(round(lev,digits=3)),
             ",  mode=refe'nce point"))

             dev.set(which = deviloca) 
             lcoordi<-1
             plotbary(stredu,coordi=lcoordi,ptext=ptextst,symbo="T")
             title(sub=paste("location plot, coordinate",as.character(lcoordi)))

             dev.set(which = devivolu)       
        }
        loc<-locator(1)
     }
  }

  # interaction in barycenter plot
  if (devi==devibary){
      coordi<-1
      icolo<-mt$colot[mt$low[indeksi]:mt$upp[indeksi]]
      inodes<-mt$nodepointer[mt$low[indeksi]:mt$upp[indeksi]]
      modlab<-plotbary(pr,coordi=coordi,ptext=ptext,
                    modlabret=T,modecolo=icolo,modepointer=inodes)
      title(sub=paste("barycenter plot, coordinate",as.character(coordi)))
      alaasso<-0
      while (loc$y>=alaasso){
         if (coordi<=(d-1)) coordi<-coordi+1 else coordi<-1
         plotbary(pr,coordi=coordi,ptext=ptext,modecolo=icolo,modepointer=inodes)
         title(sub=paste("barycenter plot, coordinate",as.character(coordi)))

         loc<-locator(1)
      }
  }

  # interaction in radius plot
  if (devi==deviradi){
       alaraja<-0
       while (loc$y>=alaraja){
          ylaraja<-max(st$level)
          while (loc$y>=ylaraja){
              if (refelab=="moodi"){ 
                   refelab<-"bary"
                   if (is.null(st.bary)){  
                       refe<-st$bary
                       st.bary<-leafsfirst(pcf,lev=lev,refe=refe)
                   }
                   st<-st.bary
              }
              else{ 
                   refelab<-"moodi"
                   st<-st.moodi
              }
              if (redu) stredu<-treedisc(st,pcf,ngrid=levnumst) else stredu<-st
              plotvolu(stredu,ptext=ptextst,symbo="T")
              if (refelab=="moodi") 
              title(sub=paste("radius plot, level=",as.character(round(lev,digits=3)),
              ",  mode=refe'nce point"))
              else 
              title(sub=paste("radius plot, level=",as.character(round(lev,digits=3)),
              ",  barycenter=refe'nce point"))

              dev.set(which = deviloca) 
              lcoordi<-1
              plotbary(stredu,coordi=lcoordi,ptext=ptextst,symbo="T")
              title(sub=paste("location plot, coordinate",as.character(lcoordi)))

              dev.set(which = deviradi)      
              loc<-locator(1)
          }
      }
  }

  # interaction in location plot
  if (devi==deviloca){
      coordi<-1
      plotbary(stredu,coordi=coordi,ptext=ptextst,symbo="T")
      title(sub=paste("coordinate",as.character(coordi)))
      alaasso<-0
      while (loc$y>=alaasso){
         if (coordi<=(d-1)) coordi<-coordi+1 else coordi<-1
         plotbary(stredu,coordi=coordi,ptext=ptextst,symbo="T")
         title(sub=paste("location plot, coordinate",as.character(coordi)))

         loc<-locator(1)
      }
  }

  # interaction in branching map
  if (devi==devibranch){
      alaraja<--0.4
      while (loc$y>=alaraja){

          if (loc$x>=0) theta<-theta+10 else theta<-theta-10
          if (loc$y>=0) phi<-phi+10 else phi<-phi-10

          persp(x=bm$level,y=bm$h,z=bm$z,col=bm$col,
          xlab="level",ylab="h",zlab="",ticktype="detailed",
          phi=phi,theta=theta)
          title(sub="map of branches")

         loc<-locator(1)
      }
  }

  # end
  dev.set(which = devicontrol)
  loc<-locator(1)
}

if (!is.null(levnumst)) devlkm<-lkm
else devlkm<-lkm-2
for (i in 1:devlkm) dev.off()

}




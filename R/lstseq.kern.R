lstseq.kern<-function(dendat,hseq,N,Q=NULL,kernel="epane",hw=NULL,
algo="leafsfirst",level=0.1)
{
hnum<-length(hseq)
if ((hnum>1) && (hseq[1]<hseq[2])) hseq<-hseq[seq(hnum,1)]

if ((is.null(algo)) || (algo=="leafsfirst")){

  for (i in 1:hnum){   
      h<-hseq[i]
      pcf<-pcf.kern(dendat,h,N,kernel=kernel,hw=hw)
      if (!is.null(algo)) lf<-leafsfirst(pcf)
      if (!is.null(level)){ 
           lev<-level*max(pcf$value)  
           refe<-locofmax(pcf)
           st<-leafsfirst(pcf,lev=lev,refe=refe)
      }
      if (i==1){
           if (hnum==1){ 
               pcfseq<-pcf
               if (!is.null(algo)) lstseq<-lf
               if (!is.null(level)) stseq<-st
           }
           else{
               pcfseq<-list(pcf)
               if (!is.null(algo)) lstseq<-list(lf)
               if (!is.null(level)) stseq<-list(st)
           }
      }
      else{
          pcfseq<-c(pcfseq,list(pcf))
          if (!is.null(algo)) lstseq<-c(lstseq,list(lf))
          if (!is.null(level)) stseq<-c(stseq,list(st))
      }
  }

}
else{  #algo=="decomdyna"
  lstseq<-profkern(dendat,hseq,N,Q,kernel=kernel,hw=hw)
}

if (is.null(algo))  lstseq<-NULL
if (is.null(level)) stseq<-NULL
return(list(lstseq=lstseq,pcfseq=pcfseq,stseq=stseq,hseq=hseq,type="kernel"))
}


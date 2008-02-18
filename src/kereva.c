/* 
gcc -Wall -ansi -pedantic /home/jsk/denpro/src/kereva.c


R CMD SHLIB -o /home/jsk/kerCeva /home/jsk/denpro/src/kereva.c

dyn.load("/home/jsk/kerCeva")


kg<-.C("kergridC",
               as.integer(extMaxnode),
               as.integer(extMaxvals),
               as.double(dendat),
               as.double(inh),
               as.integer(inN),
               as.integer(n),
               as.integer(hnum),
               as.integer(d), 
               as.integer(kertype), 
               as.double(trunc),
               as.double(threshold),
               as.double(inweig),
               ioleft = integer(extMaxnode+1),
               ioright = integer(extMaxnode+1),
               ioparent = integer(extMaxnode+1),
               infopointer = integer(extMaxnode+1),
               iolow = integer(extMaxnode+1),
               ioupp = integer(extMaxnode+1),
               value = double(hnum*extMaxvals),
               index = integer(d*extMaxvals),
               nodefinder = integer(extMaxvals),
               numpositive = integer(1),
               numnode = integer(1))

outappu1 = double(1),
outappu2 = double(1),
outappu3 = double(1),
outappu4 = double(1),
outappu5 = double(1),
outappu6 = double(1),
outappu7 = double(1),
outappu8 = double(1),
outappu9 = double(1))

*/

#include <math.h>
#include <stdlib.h> 

double epane(double *x, double hs, int d);
double gauss(double *x, double hs, int d);
int digit(int luku, int* base, int d, int* inde);
int depth2com(int dep, int* N, int d, int* hacku);
/*int *direc, int *depthInDirec);*/

void kergrid(int *extMaxnode,
              int *extMaxvals,
              double *indendat,
              double *h,
              int *N,
              int *n, 
              int *hnum,
              int *d,
              int *kertype,
              double *trunc,
              double *threshold,
              double *weig,
              /* output */
              int *ioleft,
              int *ioright,
              int *ioparent,
              int *infopointer,
              int *iolow,
              int *ioupp,
              double *iovalue,
              int *index,
              int *nodefinder,  /* pointer from estimate to tree */
              int *ionumpositive,
              int *ionumnode)

/*
double *outappu1,
double *outappu2,
double *outappu3,
double *outappu4,
double *outappu5,
double *outappu6,
double *outappu7,
double *outappu8,
double *outappu9)
*/

{

 /*double value[*extMaxvals+1];*/
 double *value = (double *)malloc(sizeof(double) * (*extMaxvals+1));

 /*double minim[*d+1], maxim[*d+1], point[*d+1], delta[*d+1];*/
 double *minim = (double *)malloc(sizeof(double) * (*d+1));
 double *maxim = (double *)malloc(sizeof(double) * (*d+1));
 double *point = (double *)malloc(sizeof(double) * (*d+1));
 double *delta = (double *)malloc(sizeof(double) * (*d+1));
 /*int gridlow[*d+1], gridupp[*d+1], base[*d+1], inde[*d+1];*/
 int *gridlow = (int *)malloc(sizeof(int) * (*d+1));
 int *gridupp = (int *)malloc(sizeof(int) * (*d+1));
 int *base = (int *)malloc(sizeof(int) * (*d+1));
 int *inde = (int *)malloc(sizeof(int) * (*d+1));

 int i, j, k, hrun, nodeloc, gridcard, pointer;
 int numpositive, numnode;
 double hmax, val, curval, hcur;
 int apu;
 int curre, curdep;
 /* digit:iin */
 /* int digitdigi[*d+1], digitjako[*d+1]; */
 int *digitjako = (int *)malloc(sizeof(int) * (*d+1));
 int *digitdigi = (int *)malloc(sizeof(int) * (*d+1));
 
 /* epane:en */
 /*double epaarg[*d+1];*/
 double *epaarg = (double *)malloc(sizeof(double) * (*d+1));

 double hs; 
 /* findend */
 int current, dep;
 int exists;        /*  BOOLEAN */
 int chil, chir;
 div_t mid;
 /* depth2com */
 int d2cdirrec, d2cdepind;      /* components of result */
 /* addnode */
 int curdir, depatd, ind;
 /* int depit[*d+1]; */
 int *depit = (int *)malloc(sizeof(int) * (*d+1));
 int anj;

/* globals 
 int left[*extMaxnode+1]; 
 int right[*extMaxnode+1]; 
 int parent[*extMaxnode+1];
 int low[*extMaxnode+1]; 
 int upp[*extMaxnode+1];
 */
 int *left = (int *)malloc(sizeof(int) * (*extMaxnode+1));
 int *right = (int *)malloc(sizeof(int) * (*extMaxnode+1));
 int *parent = (int *)malloc(sizeof(int) * (*extMaxnode+1));
 int *low = (int *)malloc(sizeof(int) * (*extMaxnode+1));
 int *upp = (int *)malloc(sizeof(int) * (*extMaxnode+1));


 /* findendResult */
 int feExistiert;   /*logical */ 
 int feLocation; 
 int feDeepness; 
 /* addnodeResult */
 int adNumberOfNodes; 
 int adNodeLocation; 
 /* apput */
 /* double appu1, appu2, appu3, appu4, appu5, appu6, appu7, appu8, appu9; */

  int hacku[3];

 /*double dendat[(*n)+1][(*d)+1];*/
 double ** dendat;
 dendat = (double **)malloc((*n+1) * sizeof(double *));
 if (NULL == dendat) exit(1);
 for (i = 0; i <= *n; i++) {
     dendat[i] = (double *)malloc((*d+1) * sizeof(double));
     if (NULL == dendat[i]) exit(1);
 }

 if (value == NULL) exit(1); 
 if (minim == NULL) exit(1); 
 if (maxim == NULL) exit(1); 
 if (point == NULL) exit(1); 
 if (delta == NULL) exit(1); 
 if (gridlow == NULL) exit(1); 
 if (gridupp == NULL) exit(1); 
 if (base == NULL) exit(1); 
 if (inde == NULL) exit(1); 
 if (epaarg == NULL) exit(1); 
 if (left == NULL) exit(1); 
 if (right == NULL) exit(1); 
 if (parent == NULL) exit(1); 
 if (low == NULL) exit(1); 
 if (upp == NULL) exit(1); 

 if (digitjako == NULL) exit(1); 
 if (digitdigi == NULL) exit(1); 

 if (depit == NULL) exit(1); 

  adNumberOfNodes=1; 

  /* Initialize dendat */
  for (j=0; j<=*d-1; j++)
     for (i=1; i<=*n; i++)
        dendat[i][j+1]=indendat[j*(*n)+i-1];

  /* initialize left and right */
  for (j=1; j<=*extMaxnode; j++){
      left[j]=0;
      right[j]=0;
      parent[j]=0;
      low[j]=0;
      upp[j]=0;
  }

  for (j=1; j<=*extMaxvals; j++){
      value[j]=0;
  }

  /* minim[i]=min(dendat[,i]) */
  for (i=1; i<=*d; i++){
     minim[i]=dendat[1][i];
     maxim[i]=dendat[1][i]; 
     for (j=1; j<=*n; j++){
	  if (dendat[j][i]<=minim[i]){
              minim[i]=dendat[j][i];
	  }
	  if (dendat[j][i]>=maxim[i]){
              maxim[i]=dendat[j][i];
	  }
      }
  }
 
  /* hmax=max(h); */
  if (*hnum==1){
      hmax=*h;
  }
  else{
      hmax=h[1];
  }

  /* delta=(maxim-minim+2*hmax)/(N+1); */
  for (i=1; i<=*d; i++){ 
    delta[i]=(maxim[i]-minim[i]+2*hmax)/(N[i]+1);
  }

  numnode=1;
  low[1]=1;
  upp[1]=N[1];

  numpositive=0;

  for (i=1; i<=*n; i++){
      for (hrun=1; hrun<=*hnum; hrun++){ 
	  /* find the grid points in the support */ 
	      if (*hnum==1){
		  hcur=*h;
              }
              else{
                  hcur=h[hrun];
              }
          for (j=1; j<=*d; j++){
            gridlow[j]=floor(((dendat[i][j]-minim[j])/delta[j])+1);
            gridupp[j]=ceil(((dendat[i][j]-minim[j]+2*hcur)/delta[j])-1);
          }
          for (j=1; j<=*d; j++){
             base[j]=gridupp[j]-gridlow[j]+1;
          }

          gridcard=1;
          for (j=1; j<=*d; j++){
	      gridcard=gridcard*base[j];
          }
         
          for (k=0; k<=(gridcard-1); k++){  

             if (*d==1){
                 for (j=1; j<=*d; j++){
		      inde[j]=gridlow[j]+k;
                 }
             }      
             else{
                 apu=digit(k,base,*d,inde); 

                 /* inde=inde+gridlow; */
                 for (j=1; j<=*d; j++){
		      inde[j]=inde[j]+gridlow[j];
                 }
	     }

             /* point is d-vector */ 
             if (*hnum==1){
                for (j=1; j<=*d; j++){
	 	    point[j]=minim[j]-(*h)+delta[j]*inde[j]; 
                }         
             }
             else{
                for (j=1; j<=*d; j++){
	 	    point[j]=minim[j]-h[hrun]+delta[j]*inde[j]; 
                }
             }

             /* epaarg<-point-dendat[i,] */
             for (j=1; j<=*d; j++){
        	  epaarg[j]=point[j]-dendat[i][j];
             }

             if (*hnum==1){
                   hs=*h;
             }
             else{
                   hs=h[hrun];
             }
             if (*kertype==1){
                   val=epane(epaarg,hs,*d);
             }
             else if (*kertype==2){
	           val=gauss(epaarg,hs,*d);
             }

	     /* find whether gridpoint is already in tree */
             /* Obs. findend need inde */
             /* apu=findend(inde,N,left,right,low,upp,*d); */

/* START findend  */


if (val>=(*threshold)){


/* apu=findend(inde,N,*d); */

 current=1;
 dep=1;
 if ((left[current]==0) & (right[current]==0)){
   exists=0;
 }
 else{
   exists=1;
 }
 while ((exists) && ((left[current]>0) || (right[current]>0))){

     mid=div((low[current]+upp[current]),2); 

     apu=depth2com(dep, N, *d, hacku); /* direc, depthInDirec);*/
     d2cdirrec=hacku[1]; 
     d2cdepind=hacku[2];
     /*
     d2cdirrec=*direc; 
     d2cdepind=*depthInDirec;
     */

     if (inde[d2cdirrec]<=mid.quot){
           chil=left[current];
           if (chil>=1){
                current=left[current];
                dep=dep+1;
           }
           else{
                exists=0;
           }
     }
     else{
           chir=right[current];
           if (chir>=1){
                  current=right[current];
                  dep=dep+1;
           }
           else{
               exists=0;
           }
     }
 }

 feExistiert=exists;
 feDeepness=dep;
 feLocation=current;

/* END findend */

             if (feExistiert==1){
                pointer=infopointer[feLocation];
                curval=value[numpositive*(hrun-1)+pointer];
		    /* curval=value[pointer][hrun]; */
                value[numpositive*(hrun-1)+pointer]=curval+val*weig[i];
		    /* curval+val/(*n); */
		    /* value[pointer][hrun]=curval+val/(*n); */
	     }
             else{  /* gridpoint was not yet in the tree */
                curre=feLocation;
                curdep=feDeepness;
                /* needs inde */
                numnode=adNumberOfNodes;
                /* apu=addnode(inde,curre,curdep,N,intnumnode,*d,left,right,
		   parent,low,upp); */

/* START addnode */

/* apu=addnode(inde,curre,curdep,N,intnumnode,*d); */

/* inde is d-vector: index (gridpoint) to be added */
/* curre is pointer to vectors left,right,... */

/* BEGIN depth2com */
/*
curdir=depth2com(curdep,N,d);
curdir=direktio; 
depatd=depthInDirec;
*/

     apu=depth2com(curdep, N, *d, hacku); /* direc, depthInDirec);*/
     curdir=hacku[1]; 
     depatd=hacku[2];

/*
 d2cdep=curdep;

 for (d2cj=1; d2cj<=*d; d2cj++){
     d2clogn[d2cj]=rint(log(N[d2cj])/log(2));  
 }
 for (d2cj=1; d2cj<=*d; d2cj++){
     d2ccusu[d2cj]=0;
     for (d2ck=1; d2ck<=d2cj; d2ck++){
        d2ccusu[d2cj]=d2ccusu[d2cj]+d2clogn[d2ck];
     }
 }
 d2cdirrec=1;
 while ((d2cdirrec<=*d) && ((d2cdep-d2ccusu[d2cdirrec])>0)){
     d2cdirrec=d2cdirrec+1;
 }
 if (d2cdirrec>*d){
     d2cdirrec=*d;
 }
 if (d2cdirrec==1){
     d2cdepind=d2cdep;
 }
 else{
     d2cdepind=d2cdep-d2ccusu[d2cdirrec-1];
 }

curdir=d2cdirrec;
depatd=d2cdepind;

*/
/* END depth2com */

 /* depit=log(N,base=2) */
 for (anj=1; anj<=*d; anj++){
     depit[anj]=rint(log(N[anj])/log(2));  /* log(N[anj],base=2) */ 
 }

 while (curdir<=(*d-1)){
     ind=inde[curdir];
     while (depatd<=depit[curdir]){
	 mid=div((low[curre]+upp[curre]),2);
	 if (ind<=mid.quot){
	     left[curre]=numnode+1;
	     parent[numnode+1]=curre;
	     low[numnode+1]=low[curre];
	     upp[numnode+1]=mid.quot;       /* floor(mid); */
	 }
         else{
	     right[curre]=numnode+1;
	     parent[numnode+1]=curre;
	     low[numnode+1]=mid.quot+1;    /* ceil(mid); */
	     upp[numnode+1]=upp[curre];
	 }
	 numnode=numnode+1;
	 curre=numnode;
	 depatd=depatd+1;
	 curdep=curdep+1;
     }
    /*
    Last node of this dimension (first node of next dimension)
    */
     curdir=curdir+1;
     ind=inde[curdir];
     low[curre]=1;
     upp[curre]=N[curdir];
     mid=div((low[curre]+upp[curre]),2);
     if (ind<=mid.quot){
	 left[curre]=numnode+1;
	 parent[numnode+1]=curre;
	 low[numnode+1]=low[curre];
	 upp[numnode+1]=mid.quot;      /* floor(mid); */
     }
     else{
         right[curre]=numnode+1;
         parent[numnode+1]=curre;
         low[numnode+1]=mid.quot+1;       /* ceil(mid); */
         upp[numnode+1]=upp[curre];
     }
     depatd=2;
     numnode=numnode+1;
     curre=numnode;
     curdep=curdep+1;
 }
/*
Last dimension 
*/
 ind=inde[curdir];
 while (depatd<=depit[curdir]){
     mid=div((low[curre]+upp[curre]),2);
     if (ind<=mid.quot){
	 left[curre]=numnode+1;
	 parent[numnode+1]=curre;
	 low[numnode+1]=low[curre];
	 upp[numnode+1]=mid.quot;      /* floor(mid); */
     }
     else{
         right[curre]=numnode+1;
         parent[numnode+1]=curre;
         low[numnode+1]=mid.quot+1;     /* ceil(mid); */
         upp[numnode+1]=upp[curre];
     }
     numnode=numnode+1;
     curre=numnode;
     depatd=depatd+1;
     curdep=curdep+1;
 }
/*
Last node of last dimension is clear already
*/
adNumberOfNodes=numnode;
adNodeLocation=numnode;

/* END addnode */

                numnode=adNumberOfNodes;
                nodeloc=adNodeLocation;

                numpositive=numpositive+1;
                infopointer[numnode]=numpositive;
                /* now hrun=1, because we add */         
                value[numpositive]=val*weig[i];
		     /* val/(*n);*/
                    /* value[(*numpositive)][hrun]=val/(*n); */
                /* index[(*numpositive),]=inde; */
                for (j=1; j<=*d; j++){
                    index[(numpositive-1)*(*d)+j]=inde[j];
                	  /* index[(*numpositive)][j]=inde[j]; */
                }
                nodefinder[numpositive]=nodeloc;
	     }

}  /* if val>threshold */

	  }   /* k <= gridcard-1 */

      }  /* hrun <= *hnum */
  }  /* i<= *n */

  for (j=1; j<=numnode; j++){
      ioleft[j]=left[j];
      ioright[j]=right[j];
      ioparent[j]=parent[j];
      iolow[j]=low[j];
      ioupp[j]=upp[j];
  }

  for (j=1; j<=(numpositive*(*hnum)); j++){
      iovalue[j]=value[j];
  }

  *ionumnode=numnode;
  *ionumpositive=numpositive;

/*
          *outappu1=appu1;
          *outappu2=appu2;
          *outappu3=appu3;
          *outappu4=appu4;
          *outappu5=appu5;
          *outappu6=appu6;
          *outappu7=appu7;
          *outappu8=appu8;
          *outappu9=appu9;

	  */


 free(value);
 free(minim); 
 free(maxim); 
 free(point); 
 free(delta); 
 free(gridlow); 
 free(gridupp); 
 free(base); 
 free(inde); 
 free(epaarg); 
 free(left); 
 free(right); 
 free(parent); 
 free(low); 
 free(upp); 

 for(i = 0; i <= *n; i++) free(dendat[i]);
 free(dendat);

 free(digitjako);
 free(digitdigi);

 free(depit);

}


double epane(double *x, double hs, int d)
{
    int i;
    double eres;

    eres=1;
    for (i=1; i<=d; i++) eres=eres*3*(1-pow((x[i]/hs),2))/(4*hs);

    return (eres);
}

double gauss(double *x, double hs, int d)
{
    int i;
    double eres, kersig, trunc, norvak;

    trunc=3;
    kersig=1/trunc;      /*0.33333333;    1/sig */
    norvak= 0.9973002;   /* sig<-3; 1-2*pnorm(-3) */ 
 
    eres=1;
    for (i=1; i<=d; i++)  eres=eres*
          exp(-pow((x[i]/(hs*kersig)),2)/2)/(hs*kersig*norvak*sqrt(2*M_PI)); 

    return (eres);
}

int digit(int luku, int* base, int d, int* inde)
{
/* returns inde (d-vector of integers)
   luku is a natural number >=0
   base is d-vector of integers >=2, d>=2,
   base[d] tarvitaan vain tarkistamaan onko luku rajoissa
   example: digit(52,c(10,10)), returns vector (2,5)
*/
    int di, vah, apu;
    /*int digitjako[d+1], digitdigi[d+1];*/
    int *digitjako = (int *)malloc(sizeof(int) * (d+1));
    int *digitdigi = (int *)malloc(sizeof(int) * (d+1));

    if (digitjako == NULL) exit(1); 
    if (digitdigi == NULL) exit(1); 

    digitjako[d]=base[1];
    di=d-1;
    while (di >= 1){
       digitjako[di]=base[d-di+1]*digitjako[di+1];
       di=di-1;
    }
    vah=0;
    di=1;
    while (di<=(d-1)){
        digitdigi[di]=floor((luku-vah)/digitjako[di+1]);
        vah=vah+digitdigi[di]*digitjako[di+1];
        di=di+1;
    }
    digitdigi[d]=luku-vah;
    di=1;
    while (di<=d){
        inde[di]=digitdigi[d-di+1];
        di=di+1;
    }
    
    apu=1;
    return apu;

    free(digitjako);
    free(digitdigi);
      

}

int depth2com(int dep, int* N, int d, int* hacku)
/*  *direc, int *depthInDirec) */
{
/* output direc, depthInDirec*/

    int d2cj, d2ck, d2cdirrec, d2cdepind, apu;
   /*double d2clogn[d+1], d2ccusu[d+1];*/
    double *d2clogn = (double *)malloc(sizeof(double) * (d+1));
    double *d2ccusu = (double *)malloc(sizeof(double) * (d+1));

    if (d2clogn == NULL) exit(1); 
    if (d2ccusu == NULL) exit(1); 



    /* d2clogn=log(N,base=2); */
    for (d2cj=1; d2cj<=d; d2cj++){
        d2clogn[d2cj]=rint(log(N[d2cj])/log(2));  /* log(N[i],base=2) */ 
    }

    /* d2ccusu=cumsum(d2clogn); */
    for (d2cj=1; d2cj<=d; d2cj++){
       d2ccusu[d2cj]=0;
       for (d2ck=1; d2ck<=d2cj; d2ck++){
          d2ccusu[d2cj]=d2ccusu[d2cj]+d2clogn[d2ck];
       }
    }

    d2cdirrec=1;
    while ((d2cdirrec<=d) && ((dep-d2ccusu[d2cdirrec])>0)){
        d2cdirrec=d2cdirrec+1;
    }

    /* d2cdirrec=min(d2cdirrec,d) */
    if (d2cdirrec>=d){
         d2cdirrec=d;
    }

    if (d2cdirrec==1){
        d2cdepind=dep;
    }
    else{
        d2cdepind=dep-d2ccusu[d2cdirrec-1];
    }

    /* return */
    /*
    *direc=d2cdirrec; 
    *depthInDirec=d2cdepind; 
    */
    hacku[1]=d2cdirrec;
    hacku[2]=d2cdepind;

    apu=1;
    return apu;

    free(d2clogn); 
    free(d2ccusu); 


}






















/* 
R CMD SHLIB -o kerCprof kerprof.c
dyn.load("/home/jsk/kerle/kerCprof")
kg<-.C("kerprofC",as.integer(extMaxnode),
                  as.integer(extMaxvals),
                  as.double(dendat),
                  as.double(inh),
                  as.integer(inN),
                  as.integer(n),
                  as.integer(hnum),
                  as.integer(d),
                  as.integer(Q),
                  as.integer(numofallR),
                  level = double(numofallR+1),
                  parent = integer(numofallR+1),
                  component = integer(numofallR+1),
                  volume = double(numofallR+1),
                  center = double(d*numofallR+1),
                  efek = integer(1))

outappu1 = double(1),
*/

#include <math.h>
#include <stdlib.h> 

#define numofallIntern 1000000
#define atomnumIntern 1000000
#define nodenumOfDyakerIntern 1000000
#define maxdim 4 

#define maxn 1001
#define maxnode 100000      /* *extMaxnode */
#define maxpositive 100000  /* *extMaxvals */
#define maxlevnum 10000      /* *levnum */
#define maxnumofallR 100000      /* *numofall */
#define maxterminalnum 100000
#define maxnodenumOfDyaker 100000
#define maxm 40

int AtomlistNext[numofallIntern+1];
int AtomlistAtom[numofallIntern+1]; 
/* point to value,..: values in 1,...,atomnum */

int begsSepaNext[atomnumIntern+1];
int begsSepaBegs[atomnumIntern+1];
int atomsSepaNext[atomnumIntern+1];
int atomsSepaAtom[atomnumIntern+1];

int begsLeftBoun[atomnumIntern+1];
int begsRighBoun[atomnumIntern+1];
int atomsLBounNext[atomnumIntern+1];
int atomsRBounNext[atomnumIntern+1];

int separy[nodenumOfDyakerIntern+1];

int begs[atomnumIntern+1];

int indeks[atomnumIntern+1][maxdim+1];


int listchangeC(int totbegSepary,
                int beg);
int declevdyaC(int beg,
               /* tree */
               int *left,
               int *right,
               int *parent,
               /* int *indeks[], */
               int *nodefinder,
               /* end tree */
               int *N,
               int nodenumOfDyaker,
	       int terminalnum,
               int d);
int joingeneC(int node,
               int leftbeg,
               int rightbeg,
	       int direction,
	      /* int *indeks[], */
	      int d);
int joinconneC(int leftbeg,
               int rightbeg, 
               int direction, 
               /* int *indeks[], */
               int d);
int startpoints(int leftbeg, 
                int rightbeg);
int makelinks(int direction,
              int d);
int declevnewC();
int joinSets(int leftbeg,
             int rightbeg,
             int sepnum);


void kerprofC(int *extMaxnode,
              int *extMaxvals,
              double *indendat,
              double *h,
              int *N,
              int *n, 
              int *hnum,
              int *d,
              int *levnum,
              int *numofallR,
              /* output */
              double *level,      /* of length numofall+1 */
              int *parent, 
              int *component, 
	      double *volume,
              double *center,
              int *efek)   /* gives effective length of level,... */

{
 double dendat[(maxn)+1][(maxdim)+1];
 double value[maxpositive+1];
 double minim[maxdim+1], maxim[maxdim+1], point[maxdim+1], delta[maxdim+1];
 int gridlow[maxdim+1], gridupp[maxdim+1], base[maxdim+1], inde[maxdim+1];
 int infopointer[maxnode+1]; 
 int i, j, k, hrun, nodeloc, gridcard, pointer;
 int numpositive, numnode;
 double hmax, val, curval, hcur;
 int curre, curdep;
 /* digit:iin */
 int digitdigi[maxdim+1], digitjako[maxdim+1], digitapu[maxdim+1];
 int luku, vah, di, dj;
 /* epane:en */
 double epaarg[maxdim+1], x[maxdim+1];
 int ei, ej;
 double eres, hs; 
 /* findend */
 int current, dep;
 int exists;        /*  BOOLEAN */
 int chil, chir;
 div_t mid;
 /* depth2com */
 int d2clogn[maxdim+1], d2ccusu[maxdim+1];
 int d2cdirrec, d2cdepind;      /* components of result */
 int d2cj, d2ck;
 /* addnode */
 int curdir, depatd, ind;
 int depit[maxdim+1];
 int anj;
 /* depth2com */
 int d2cdep;

/* previous globals */
int left[maxnode+1]; 
int right[maxnode+1]; 
int parentKer[maxnode+1];
int low[maxnode+1]; 
int upp[maxnode+1];
/* findendResult */
int feExistiert;   /*logical */ 
int feLocation; 
int feDeepness; 
/* addnodeResult */
int adNumberOfNodes; 
int adNodeLocation; 
 /* apput */
/* double appu1, appu2, appu3, appu4, appu5, appu6, appu7, appu8, appu9; */

/* from input of decomdyaC to variable definitions  */

int numofall;
int levfrekv[maxlevnum];
double levseq[maxlevnum], step;
double volofatom, maxval, minval;
int nodenumOfDyaker;

int nodefinder[maxpositive+1];

/* katkaisukohta */

  int totbegSepary;

  int pinoComponent[maxnumofallR+1];  /* pointer to component, level,... */
  int pinoTaso[maxnumofallR+1];       /* ordinal of level (pointer to levseq) */

  int componum, beghigher;
  int beg, koko, pinind, levind, partlevsetbeg;  /* i,j,ind,apu */
  int listEnd, PrevlistEnd;
  int terminalnum, addnum, removenum, runner, origiListEnd, atom;
  double arvo, higlev;
  int lokalefek;
  /*  int indeks[(*atomnum)+1][(maxdim)+1]; */
 
  int componentnum, numofatoms, zeiger; /* for cvolumdya */
  double curcente[maxdim+1];                   /* for ccentedya */
  int atompointer;



  adNumberOfNodes=1; 

  /* Initialize dendat */
  for (j=0; j<=*d-1; j++)
     for (i=1; i<=*n; i++)
        dendat[i][j+1]=indendat[j*(*n)+i-1];

  /* initialize left and right */
  for (j=1; j<=*extMaxnode; j++){
      left[j]=0;
      right[j]=0;
      parentKer[j]=0;
      low[j]=0;
      upp[j]=0;
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
   for (i=1; i<=*hnum; i++){
       hmax=h[1];
       if (hmax>=h[i]){
             hmax=h[i];
       }
   }
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
		 /* 1. inde=digit(k,base); */
                 /* 2. apu=digit(k,base); */  /* inde is d-vector */

/* begin DIGIT */
/* int digitdigi[*d], digitjako[*d], digitapu[*d]; */

 luku=k;
 digitjako[*d]=base[1];
 di=*d-1;
 while (di >= 1){
     digitjako[di]=base[*d-di+1]*digitjako[di+1];
     di=di-1;
 }
 vah=0;
 di=1;
 while (di<=((*d)-1)){
     digitdigi[di]=floor((luku-vah)/digitjako[di+1]);
     vah=vah+digitdigi[di]*digitjako[di+1];
     di=di+1;
 }
 digitdigi[*d]=luku-vah;
 di=1;
 while (di<=(*d)){
     digitapu[di]=digitdigi[(*d)-di+1];
     di=di+1;
 }
 /* inde=digitapu; */
 for (dj=1; dj<=(*d); dj++){
     inde[dj]=digitapu[dj];
 }

/* END DIGIT */


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

/* val=epane(epaarg,h[hrun],*d); */

/* START EPANE */

/* x=epaarg; */
for (ej=1; ej<=(*d); ej++){
     x[ej]=epaarg[ej];
}
 if (*hnum==1){
    hs=*h;
 }
 else{
    hs=h[hrun];
 }
eres=1;
ei=1;
while (ei<=*d){
      eres=eres*3*(1-pow((x[ei]/hs),2))/2;
      ei=ei+1;
}

val=eres;

/* END EPANE */

	     /* find whether gridpoint is already in tree */
             /* Obs. findend need inde */
             /* apu=findend(inde,N,left,right,low,upp,*d); */

/* START findend  */


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

/* BEGIN depth2com */
/*
 direc=depth2com(dep,N,*d).direktio; 
direc=depth2com(dep,N,d);
 direc=direktio; 
*/

/*input:   int dep, int *N, int d */

  /* d2clogn=log(N,base=2); */
 for (d2cj=1; d2cj<=*d; d2cj++){
     d2clogn[d2cj]=rint(log(N[d2cj])/log(2));  /* log(N[i],base=2) */ 
 }

 /* d2ccusu=cumsum(d2clogn); */
 for (d2cj=1; d2cj<=*d; d2cj++){
     d2ccusu[d2cj]=0;
     for (d2ck=1; d2ck<=d2cj; d2ck++){
        d2ccusu[d2cj]=d2ccusu[d2cj]+d2clogn[d2ck];
     }
 }

 d2cdirrec=1;
 while ((d2cdirrec<=*d) && ((dep-d2ccusu[d2cdirrec])>0)){
     d2cdirrec=d2cdirrec+1;
 }

 /* d2cdirrec=min(d2cdirrec,d) */

 if (d2cdirrec>=(*d)){
     d2cdirrec=*d;
 }

 if (d2cdirrec==1){
     d2cdepind=dep;
 }
 else{
     d2cdepind=dep-d2ccusu[d2cdirrec-1];
 }

 /* direktio=d2cdirrec; */
 /* depthInDirec=d2cdepind; */
/* return */
/* direc=d2cdirrec; */

/* END depth2com */

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
                curval=value[(pointer-1)*(*hnum)+hrun];
		    /* curval=value[pointer][hrun]; */
                value[(pointer-1)*(*hnum)+hrun]=curval+val/(*n);
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

/* inde is d-vector: indeks (gridpoint) to be added */
/* curre is pointer to vectors left,right,... */

/* BEGIN depth2com */
/*
curdir=depth2com(curdep,N,d);
curdir=direktio; 
depatd=depthInDirec;
*/

/*input:   int curdep, int *N, int d */
 d2cdep=curdep;

 /* d2clogn=log(N,base=2); */
 for (d2cj=1; d2cj<=*d; d2cj++){
     d2clogn[d2cj]=rint(log(N[d2cj])/log(2));  /* log(N[i],base=2) */ 
 }
 /* d2ccusu=cumsum(d2clogn); */
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
 /* d2cdirrec=min(d2cdirrec,d) */
 if (d2cdirrec>*d){
     d2cdirrec=*d;
 }
 if (d2cdirrec==1){
     d2cdepind=d2cdep;
 }
 else{
     d2cdepind=d2cdep-d2ccusu[d2cdirrec-1];
 }
 /* direktio=d2cdirrec; */
 /* depthInDirec=d2cdepind; */

/* return */

curdir=d2cdirrec;
depatd=d2cdepind;

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
	     parentKer[numnode+1]=curre;
	     low[numnode+1]=low[curre];
	     upp[numnode+1]=mid.quot;       /* floor(mid); */
	 }
         else{
	     right[curre]=numnode+1;
	     parentKer[numnode+1]=curre;
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
	 parentKer[numnode+1]=curre;
	 low[numnode+1]=low[curre];
	 upp[numnode+1]=mid.quot;      /* floor(mid); */
     }
     else{
         right[curre]=numnode+1;
         parentKer[numnode+1]=curre;
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
	 parentKer[numnode+1]=curre;
	 low[numnode+1]=low[curre];
	 upp[numnode+1]=mid.quot;      /* floor(mid); */
     }
     else{
         right[curre]=numnode+1;
         parentKer[numnode+1]=curre;
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
                value[(numpositive-1)*(*hnum)+hrun]=val/(*n);
                    /* value[(*numpositive)][hrun]=val/(*n); */
                /* indeks[(*numpositive),]=inde; */
                for (j=1; j<=*d; j++){
                       /* indeks[(numpositive-1)*(*d)+j]=inde[j]; */
                       indeks[numpositive][j]=inde[j]; 
                }
                nodefinder[numpositive]=nodeloc;
	     }
	  }   /* k <= gridcard-1 */

      }  /* hrun <= *hnum */
  }  /* i<= *n */


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



/***************** TRANSFER KOHTA ******************************/
/* huom oletetaan, etta hnum=1 !!!!!!  *****************/

 maxval=value[1];
 minval=value[1];
 for (j=0; j<=numpositive; j++){
     if (value[j]>=maxval){
        maxval=value[j];
     }
     if (value[j]<=minval){
        minval=value[j];
     }
 }
 step=(maxval-minval)/(*levnum);
 levseq[1]=minval;
 for (j=2; j<=*levnum; j++){
    levseq[j]=levseq[j-1]+step;
 }
 
 for (j=1; j<=(numpositive); j++){
    for (k=1; k<=(*levnum); k++){
        if (value[j]>=levseq[k]){
           levfrekv[k]=levfrekv[k]+1;
        }
    }
 }
 numofall=0;
 for (j=1; j<=(*levnum); j++){
    numofall=numofall+levfrekv[j];
 }

 /* volofatom=prod(delta) */
 volofatom=1;
 for (j=1; j<=*d; j++){
    volofatom=volofatom*delta[j];
 } 

 nodenumOfDyaker=numnode;


/**************  KATKAISUKOHTA ********************/






  /* Initialize indeks 

  for (j=0; j<=*d-1; j++)
     for (i=1; i<=numpositive; i++)
	indeks[i][j+1]=inindex[j*(numpositive)+i-1]; 
  */

  /* Initialize the global variables */

  for (i=1; i<=atomnumIntern; i++){
    begsSepaNext[i]=0;
    begsSepaBegs[i]=0;
    atomsSepaNext[i]=0;
    atomsSepaAtom[i]=0;

    begsLeftBoun[i]=0;
    begsRighBoun[i]=0;
    atomsLBounNext[i]=0;
    atomsRBounNext[i]=0;
    begs[i]=0;
  }

  for (i=1; i<=nodenumOfDyakerIntern; i++)
     separy[i]=0;

  /* end of initializing */
  /* initialize the local vectors */

  for (i=1; i<=(*numofallR); i++){
      pinoComponent[i]=0;
      pinoTaso[i]=0; 
  }

  /*end of initializing the local vectors */

  componum=numofall;

  /* Initilize the lists */
  i=1;
  while (i<=numpositive){
    AtomlistAtom[i]=i;
    AtomlistNext[i]=i+1;
    i=i+1;
  }
  listEnd=numpositive;
  AtomlistNext[listEnd]=0;

  /* Let us divide the lowest level set to disconnected parts */


  beg=1;
  terminalnum=numpositive;
  totbegSepary=declevdyaC(beg,
                          left,
                          right,
                          parentKer,
                          /* indeks, */
                          nodefinder,
	                  N,
                          nodenumOfDyaker,
                          terminalnum,
                          *d); 

  /* modify begs, at most "terminalnum" changes */  
  koko=listchangeC(totbegSepary,
                   beg);    

 /* Talletetaan osat */ 
 i=1;
 while (i<=koko){
    component[i]=begs[i];
    level[i]=levseq[1];     /* arvo toistuu */
    /* Laitetaan kaikki osat pinoon */
    pinoComponent[i]=i;     /* 1,2,...,koko */
    pinoTaso[i]=1;          /* kaikki osat kuuluvat alimpaan tasojoukkoon */
    i=i+1;
 }
 lokalefek=koko;   /* kirjataan uusien osien lkm  ????? jos vain yksi */
 pinind=koko; /* indeksi pinoon */

if (*levnum>1){  while (pinind>=1){
    /* Take from stack */
    ind=pinoComponent[pinind];      /* indeksi tasoon */
    levind=pinoTaso[pinind];        /* ko tason korkeus */
    pinind=pinind-1;                /* otettu pinosta */
    partlevsetbeg=component[ind];  
    higlev=levseq[levind+1];

    /* Make intersection with the curr. component and higher lev.set */
    PrevlistEnd=listEnd;
    addnum=0;    /* num of atoms in the intersection */
    removenum=0; /* num of atoms which have to be removed to get intersec. */
    runner=partlevsetbeg;
    origiListEnd=listEnd;
    /* value=kg$value */

    while ((runner>0) && (runner<=origiListEnd)){
      atom=AtomlistAtom[runner];
      arvo=value[atom];
      if (arvo>=higlev){
	  listEnd=listEnd+1;    
	  AtomlistAtom[listEnd]=atom;
	  AtomlistNext[listEnd]=listEnd+1; 
	  addnum=addnum+1;     
      }
      else{           
          removenum=removenum+1;
      }                                
      runner=AtomlistNext[runner]; 
    }
    AtomlistNext[listEnd]=0;      /* we have to correct the end to terminate */

    if (addnum>0){
      AtomlistNext[PrevlistEnd]=0;
      beghigher=PrevlistEnd+1;
    }
    if (removenum==0){ /* jos leikkaus ei muuta, niin tasoj sailyy samana  */
      level[ind]=levseq[levind+1];  /* remove lower part */
      /* component and parent stay same, it is enough to change level */
      if (levind+1<*levnum){ /*jos ei olla korkeimmalla tasolla,laita pinoon*/
        pinoComponent[pinind+1]=ind;     
        pinoTaso[pinind+1]=levind+1; /* tasojouk taso on levind+1  */
        pinind=pinind+1;       
      }
    }
    else if (addnum>0){     /* leikkaus ei tyhja */
      beg=beghigher;
      terminalnum=addnum;
      totbegSepary=declevdyaC(beg,
                            left,
                            right,
                            parentKer,
                            /* indeks, */
                            nodefinder,
	                    N,
                            nodenumOfDyaker,
                            terminalnum,
                            *d); 

      /* modify begs, at most "terminalnum" changes */  
      koko=listchangeC(totbegSepary,
                       beg);


      /*  paivitetaan kumu tulokseen */
      i=lokalefek+1;
      while (i<=(lokalefek+koko)){
        level[i]=levseq[levind+1]; /* arvo toistuu */
        parent[i]=ind;
        component[i]=begs[i-lokalefek];
        i=i+1;
      }
      lokalefek=lokalefek+koko;
      if (levind+1<*levnum){ /*jos ei olla korkeimmalla tasolla,laita pinoon*/
        i=pinind+1;
        j=lokalefek-koko+1;
        while (i<=(pinind+koko)){
	  pinoComponent[i]=j;
          pinoTaso[i]=levind+1;  /* tasjouk tas on levind+1 */
          i=i+1;
          j=j+1;
        }
        pinind=pinind+koko;            
      }
   }
}} 
*efek=lokalefek;

/*
volume=cvolumdyaC(volofatom,component,lokalefek);

double *cvolumdyaC(double volofatom, int *component, int lokalefek)
{
  int componentnum, numofatoms, zeiger, i, apu;
  double voluumi[lokalefek];
*/

  componentnum=lokalefek;    /* length(component) */

  /* it is enough to calculate the number of atoms in each component */
for (i=1; i<=componentnum; i++){
  numofatoms=0;
  zeiger=component[i];
   while (zeiger>0){
     numofatoms=numofatoms+1;
     zeiger=AtomlistNext[zeiger];
   }
   volume[i]=numofatoms*(volofatom);
}
/*
 return voluumi;
} 
*/                                                

/*
ccentedya<-function(volofatom,component,AtomlistNext,AtomlistAtom,
volume,minim,h,delta,indeks,d){
*/

 componentnum=lokalefek; 

for (i=1; i<=componentnum; i++){  
  for (j=1; j<=*d; j++) curcente[j]=0;
  zeiger=component[i];
   while (zeiger>0){
     atompointer=AtomlistAtom[zeiger];
     for (j=1; j<=*d; j++)
       curcente[j]=curcente[j]+minim[j]-*h+delta[j]*indeks[atompointer][j];
     zeiger=AtomlistNext[zeiger];
   }
   for (j=1; j<=*d; j++)
     center[(i-1)*(*d)+j]=(volofatom)*curcente[j]/volume[i];
}
/* return(t(center)) */
                                               





/*   
for (i=1; i<=20; i++)
apu[i]=appuva[i]; 
*/

/*     
AtomlistAtomReturn=AtomlistAtom;
AtomlistNextReturn=AtomlistNext; 
*/

}



/********** END OF MAIN PROGRAM  *************************/





/* changes AtomlistNext, AtomlistAtom */
/* needs begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom */

/* create begs: beginnings of lists of atoms */
/* beg is indeks to AtomlistAtom/Next */
/* totbegsepary is indeks to begsSepaBegs/Next */

int listchangeC(int totbegSepary,
                int beg)
{
 int runnerBegs, runnerAtoms, runnerOrigi, runnerOrigiprev, sepalkm;

 runnerBegs=totbegSepary;  
 /* total beginning of list is at the root of the tree */
 runnerOrigi=beg;
 runnerOrigiprev=beg;
 sepalkm=0;
 while (runnerBegs>0){
   sepalkm=sepalkm+1;
   runnerAtoms=begsSepaBegs[runnerBegs];
   begs[sepalkm]=runnerOrigi;
   /* first step (in order to get also runnerOrigiprev to play) */
   AtomlistAtom[runnerOrigi]=atomsSepaAtom[runnerAtoms];
   runnerOrigiprev=runnerOrigi;
   runnerOrigi=AtomlistNext[runnerOrigi];
   runnerAtoms=atomsSepaNext[runnerAtoms];
   while (runnerAtoms>0){
     AtomlistAtom[runnerOrigi]=atomsSepaAtom[runnerAtoms];
     runnerOrigiprev=runnerOrigi;
     runnerOrigi=AtomlistNext[runnerOrigi];
     runnerAtoms=atomsSepaNext[runnerAtoms];
   }
   AtomlistNext[runnerOrigiprev]=0;  /* mark the end of the list */
   runnerBegs=begsSepaNext[runnerBegs];
 }
 /* begs=begs[1:sepalkm] */
 return sepalkm;
}









#include <math.h>

int declevdyaC(int beg,
               /* tree */
               int *left,
               int *right,
               int *parent,
               /* int *indeks[], */
               int *nodefinder,
               /* end tree */
               int *N,
               int nodenumOfDyaker,
	       int terminalnum,
               int d)
{
 int nextFloor[maxterminalnum+1];
 int currFloor[maxterminalnum+1];
 int already[maxnodenumOfDyaker+1];
 
 int i, k, r, lkm, nexlkm, curlkm, curre, atom, node, note;
 int exists, Lempty, Rempty;
 int leftbeg, rightbeg, direction, akku, totbegSepary, apu;
 int j;

  /* Initialize the global variables */

  for (i=1; i<=atomnumIntern; i++){
    begsSepaNext[i]=0;
    begsSepaBegs[i]=0;
    atomsSepaNext[i]=0;
    atomsSepaAtom[i]=0;

    begsLeftBoun[i]=0;
    begsRighBoun[i]=0;
    atomsLBounNext[i]=0;
    atomsRBounNext[i]=0;
    begs[i]=0;
  }

  for (i=1; i<=nodenumOfDyakerIntern; i++)
     separy[i]=0;

  /* end of initializing */

 /* initialize */

 for (i=1; i<=terminalnum; i++){
    nextFloor[i]=0;
    currFloor[i]=0;
 }
 for (i=1; i<=nodenumOfDyaker; i++)
    already[i]=0;
 

 /* INITIALIZING: "we go through the nodes at depth "depoftree""
 we make currFloor to be one over bottom floor and initialize 
 separy, boundary, atomsSepaAtom, atomsBounAtom */

 lkm=0;
 curlkm=0;
 curre=beg;
 while(curre>0){
   lkm=lkm+1;
   atom=AtomlistAtom[curre];
   node=nodefinder[atom];
    
   separy[node]=lkm;
   atomsSepaAtom[lkm]=atom;
    
   exists=parent[node];
   if (already[exists]==0){
     curlkm=curlkm+1;
     currFloor[curlkm]=exists; 
     already[exists]=1;
   }
    
   curre=AtomlistNext[curre];
 }   /* obs terminalnum=lkm */
 /* initialize the rest */
 for (r=1; r<=terminalnum; r++){
       begsSepaBegs[r]=r;
       begsLeftBoun[r]=r;
       begsRighBoun[r]=r;
 }
 /* obs: we need not change 
    begsSepaNext, atomsSepaNext, atomsLBounNext, atomsRBounNext
    since at the beginning set consist only one member:
    pointer is always 0, since we do not have followers */

/* START the MAIN LOOP */

 i=d;
 while (i >= 2){
   j=rint(log(N[i])/log(2));  /* log(N[i],base=2)  #depth at direction d */
   while (j>=1){
     nexlkm=0;
     k=1;
     while (k <= curlkm){
       node=currFloor[k];  
       /* we create simultaneously the upper floor */
       exists=parent[node];
       if (already[exists]==0){
	 nexlkm=nexlkm+1;
	 nextFloor[nexlkm]=exists;
	 already[exists]=1;
       }
/* now we join childs  */
       leftbeg=left[node];
       rightbeg=right[node];
       direction=i;
       apu=joingeneC(node,
                leftbeg,
                rightbeg,
                direction,
		     /* indeks, */
                d);
       k=k+1;
     }
     j=j-1;
     curlkm=nexlkm;
     /*     currFloor=nextFloor;  */
     for (r=1; r<=terminalnum; r++) currFloor[r]=nextFloor[r];
   }
   /* now we move to the next direction, correct boundaries */
   /*   
   begsLeftBoun=begsSepaBegs;
   begsRighBoun=begsSepaBegs;
   
   atomsLBounNext=atomsSepaNext;
   atomsRBounNext=atomsSepaNext; 
   */
 
   for (r=1; r<=atomnumIntern; r++) begsLeftBoun[r]=begsSepaBegs[r];
   for (r=1; r<=atomnumIntern; r++) begsRighBoun[r]=begsSepaBegs[r];
   for (r=1; r<=atomnumIntern; r++) atomsLBounNext[r]=atomsSepaNext[r];
   for (r=1; r<=atomnumIntern; r++) atomsRBounNext[r]=atomsSepaNext[r];

   i=i-1;
}


/* ENO OF MAIN LOOP */

 /* LAST DIMENSION WILL BE handled, (because this contains root node  */
 i=1;
 j=rint(log(N[i])/log(2));  /* log(N[i],base=2)  #depth at direction d */
 while (j>=2){
   nexlkm=0;
   k=1;
   while (k <= curlkm){
     node=currFloor[k];  
     /* we create simultaneously the upper floor */
     exists=parent[node];
     if (already[exists]==0){
       nexlkm=nexlkm+1;
       nextFloor[nexlkm]=exists;
       already[exists]=1;
     }
     /* now we join childs */
     /* if (right[parent[node]]==node)  #if node is right child  */
     leftbeg=left[node];
     rightbeg=right[node];
     direction=1; 
     apu=joingeneC(node,
              leftbeg,
              rightbeg,
              direction, 
		   /* indeks, */
              d);
       k=k+1;
   }
   j=j-1;
   curlkm=nexlkm;
     /*     currFloor=nextFloor;  */
     for (r=1; r<=terminalnum; r++) currFloor[r]=nextFloor[r];
}

 /* ROOT NODE, we do not anymore update boundaries */
 k=1;
 node=currFloor[k];  
 while (k <= 1){
   /* now we join childs */            
   /* if (right[parent[node]]==node)  if node is right child */ 
   leftbeg=left[node];
   rightbeg=right[node];
            if ((leftbeg==0) || (separy[leftbeg]==0)){
	      /* if left child does not exist */
	      separy[node]=separy[rightbeg];
            }
            else{   /* eka else */
                if ((rightbeg==0) || (separy[rightbeg]==0)){  
		  /* right child does not exist */
		  separy[node]=separy[leftbeg];
                } 
                else{   /* toka else: both children exist */
		  /* check whether left boundary of right child is empty */
		  Lempty=1;        /* TRUE */
                  note=separy[rightbeg];
                  while (note>0){
                        if (begsLeftBoun[note]>0){
			  Lempty=0;   /* FALSE */
                        }
                        note=begsSepaNext[note];
                     }
		  /* check whether right bound of left child is empty */
		  Rempty=1;    /* TRUE */
		  note=separy[leftbeg];
                  while (note>0){
                          if (begsRighBoun[note]>0){
			    Rempty=0;  /* FALSE */
                          }
                          note=begsSepaNext[note];
                     }
		  /* check whether one of boundaries is empty */
                     if (Lempty || Rempty){
		       /* one of boundaries is empty */
/* concatenating separate parts  */

akku=separy[leftbeg];
while (begsSepaNext[akku]>0){
  akku=begsSepaNext[akku];
}                           
begsSepaNext[akku]=separy[rightbeg];
separy[node]=separy[leftbeg];

/* end of concatenating, handle next boundaries */
                    }
		     else{ /* both children exist, both boundaries non-empty */
		       direction=i;
                       separy[node]=joinconneC(leftbeg,
                                               rightbeg,
                                               direction,
                                               /* indeks, */
                                               d);
                     }
                } /* toka else */
            } /* eka else */
	    k=k+1;
}
/* end of child joining */
/* END of ROOT */

 totbegSepary=separy[node];
 return totbegSepary;

 /* we have changed the vectors
 begsSepaNext=begsSepaNext,begsSepaBegs=begsSepaBegs,
 atomsSepaNext=atomsSepaNext,atomsSepaAtom=atomsSepaAtom
 */
}













int joingeneC(int node,
               int leftbeg,
               int rightbeg,
	       int direction,
	      /* int *indeks[], */
               int d)
{
  int Lempty, Rempty;
  int akku, note, apu;

  if ((leftbeg==0) || (separy[leftbeg]==0)){
    /* if left child does not exist    
       note that since we consider subsets of the
       terminal nodes of the original tree, it may happen
       that leftbeg>0 but left child does not exist */
    separy[node]=separy[rightbeg];
    /* we need that all lists contain as many members
       left boundary is empty, but we will make it a list
       of empty lists */
    note=separy[node];
    while (note>0){
      begsLeftBoun[note]=0;
      note=begsSepaNext[note];
    }
    /* right boundary stays same as for rightbeg */
  }
  else{   /* eka else */
    if ((rightbeg==0) || (separy[rightbeg]==0)){  
      /* right child does not exist */
      separy[node]=separy[leftbeg];
      /* left boundary stays same as for leftbeg right boundary is empty */
      note=separy[node];
      while (note>0){
	begsRighBoun[note]=0;
	note=begsSepaNext[note];
      }
    } 
    else{   /* toka else: both children exist */
	    /* check whether left boundary of right child is empty */
      Lempty=1;    /* TRUE */
      note=separy[rightbeg];
      while (note>0){
        if (begsLeftBoun[note]>0){
	  Lempty=0;     /* FALSE */
        }
        note=begsSepaNext[note];
      }
      /* check whether right bound of left child is empty */
      Rempty=1;    /* TRUE */
      note=separy[leftbeg];
      while (note>0){
        if (begsRighBoun[note]>0){
	  Rempty=0;    /* FALSE */
        }
        note=begsSepaNext[note];
      }
      /* check whether one of boundaries is empty */
      if (Lempty || Rempty){ /* one of boundaries is empty */
	/* concatenating separate parts */
	/* and updating boundaries for the separate parts */
        akku=separy[leftbeg];
	begsRighBoun[akku]=0; 
        /* right boundaries of sets in left child are empty */
	/* begsLeftBoun[akku] does not change */
        while (begsSepaNext[akku]>0){
          akku=begsSepaNext[akku];
          begsRighBoun[akku]=0;
        }                           
        begsSepaNext[akku]=separy[rightbeg]; 
        /* concatenate list of separate sets */
        separy[node]=separy[leftbeg];
        akku=separy[rightbeg];
        begsLeftBoun[akku]=0; 
        /* left boundaries of sets in right child are empty */
        while (begsSepaNext[akku]>0){
          akku=begsSepaNext[akku];
          begsLeftBoun[akku]=0;
        }        
        /* end of concatenating */
      }
      else{  /* both children exist, both boundaries non-empty */  
        separy[node]=joinconneC(leftbeg, 
                                rightbeg,
                                direction,
                                /* indeks, */
                                d); 
        /* separy[node]<-jc$totbegSepary */
      }
    } /* toka else */
 } /* eka else */
  apu=1;
  return apu;
}










#define induksiInt 20
#define componumInt 20

/* variables for startpoints */
int startpointsS[induksiInt+1];
int startpointsB[induksiInt+1];
int startpointsNewBleft[induksiInt+1];
int startpointsNewBright[induksiInt+1];
int m, mleft, mright;

int linkit[componumInt+1][componumInt+1];
int res[componumInt+1][componumInt+1];

int joinconneC(int leftbeg,
               int rightbeg, 
               int direction, 
               /* int *indeks[], */
               int d)         /* dimension */
{
  int sepnum, totbegsepary, apu;
  int i, j;

  /* Initialize the global variables */

  for (i=1; i <= induksiInt; i++){
     startpointsS[i]=0;
     startpointsB[i]=0;
     startpointsNewBleft[i]=0;
     startpointsNewBright[i]=0; 
  }

  for (i=1; i<=componumInt; i++)
     for (j=1; j<=componumInt; j++){
        linkit[i][j]=0;
        res[i][j]=0;
     }
  
  m=0;
  mleft=0;
  mright=0;


 /* 1. new boundary: left bound. of left child, right b. of right child */

 apu=startpoints(leftbeg,
                 rightbeg);

  /* 2. We make "links" matrix and apply declev */

  apu=makelinks(direction,
		/* mleft, m, indeks, */
                d);
  sepnum=declevnewC();  /* tulos will be filled */

   /* res is sepnum*m-matrix, 1 in some row indicates that set (atom) */
   /* belongs to this component, 0 in other positions */

   /* 3. We join the sets  */

   /* We join the sets whose startpoints are in */
   /* startpointsS and startpointsNewBleft, startpointsNewBright */
   /* We have pointers separy[leftbeg] and separy[rightbeg] */
   /* which contain pointers to lists which we can utilize */
   /* to make a new list (these two lists contain together at most as many */ 
   /* elements as we need) */

   /*  cut first list or (join these two lists and cut second) */

 totbegsepary=joinSets(leftbeg,
                       rightbeg,
                       sepnum);
 return totbegsepary;
}











int startpoints(int leftbeg, 
                int rightbeg)
{
  int anfang, apu;
  int induksi=1;

 anfang=separy[leftbeg];
 startpointsS[induksi]=anfang;
 startpointsB[induksi]=begsRighBoun[anfang];
 startpointsNewBleft[induksi]=begsLeftBoun[anfang];
 while (begsSepaNext[anfang]>0){
   anfang=begsSepaNext[anfang];
   induksi=induksi+1;
   startpointsS[induksi]=begsSepaBegs[anfang];
   startpointsB[induksi]=begsRighBoun[anfang];
   startpointsNewBleft[induksi]=begsLeftBoun[anfang];  
 }
 mleft=induksi;
 induksi=induksi+1;
 anfang=separy[rightbeg];
 startpointsS[induksi]=anfang;
 startpointsB[induksi]=begsLeftBoun[anfang];
 startpointsNewBright[induksi]=begsRighBoun[anfang];
 while (begsSepaNext[anfang]>0){
   anfang=begsSepaNext[anfang];
   induksi=induksi+1;
   startpointsS[induksi]=begsSepaBegs[anfang];
   startpointsB[induksi]=begsLeftBoun[anfang];
   startpointsNewBright[induksi]=begsRighBoun[anfang];  
 }
 /* startpointsS=startpointsS[1:induksi]
 startpointsB=startpointsB[1:induksi]
 startpointsNewBleft=startpointsNewBleft[1:induksi]
 startpointsNewBright=startpointsNewBright[1:induksi] */
 m=induksi;
 mright=m-mleft;
 apu=1;  
 return apu;
}




/* fills the linkit-matrix */
int makelinks(int direction,
              /* int *indeks[], */
              int d)
{   
 int dod, re;
 int conne;
 int beg1, beg2, begbeg1, begbeg2;
 int atom1, atom2;
 int touch;
 int i, apu;

 dod=1;
 while (dod <= mleft){
   beg1=startpointsB[dod];    /* could be 0 */
   re=mleft+1;
   while (re <= m){
     beg2=startpointsB[re];    /* could be 0 */
     conne=0;     /* FALSE */
     begbeg1=beg1;
     while (begbeg1>0){
       begbeg2=beg2;
       while (begbeg2>0){
	 atom1=atomsSepaAtom[begbeg1];
	 atom2=atomsSepaAtom[begbeg2];
         /* dotouch ?? */
         /* d=length(inde1) */
	   touch=1;  /* TRUE */
	   i=direction;
           while (i <= d){
             if ((indeks[atom1][i]>(indeks[atom2][i]+1)) || 
                 (indeks[atom1][i]<indeks[atom2][i]-1)){
	       touch=0;  /* FALSE */
             }
	     i=i+1;
           }                
           if (touch) conne=1;   /* TRUE */
	   begbeg2=atomsLBounNext[begbeg2];
       }
       begbeg1=atomsRBounNext[begbeg1];
     }                
     if (conne){
       linkit[dod][re]=1;
     }
     re=re+1;
   }
   dod=dod+1;
 }
 dod=mleft+1;
 while (dod <= m){
   beg1=startpointsB[dod];
   re=1;
   while (re <= mleft){
     beg2=startpointsB[re];
     conne=0;    /* FALSE */
     begbeg1=beg1;
     while (begbeg1>0){
       begbeg2=beg2;
       while (begbeg2>0){
	 atom1=atomsSepaAtom[begbeg1];
	 atom2=atomsSepaAtom[begbeg2];
         /* do touch ?? */
         /* d=length(inde1) */
	   touch=1;  /* TRUE */
	   i=direction;
           while (i <= d){
             if ((indeks[atom1][i]>indeks[atom2][i]+1) || 
                 (indeks[atom1][i]<indeks[atom2][i]-1)){
	       touch=0;  /* FALSE */
             }
	     i=i+1;
           }                
           if (touch) conne=1;   /* TRUE */
           begbeg2=atomsRBounNext[begbeg2];
       }
       begbeg1=atomsLBounNext[begbeg1];
     }                
     if (conne){
	 linkit[dod][re]=1;
     }
     re=re+1;
   }
   dod=dod+1;
 }
 /* huom ylla on nopeutettu, koska tiedetaan, etta atomit */
 /* 1,...,mleft eivat koske toisiaan ja samoin atomit mleft+1,...,m */
 apu=1;
 return apu;
} 









/* linkit m*m-matrix */
/* return m*m-matrix tulos and integer tulospit */
/* (first tulospit rows contain information */
int declevnewC()   /* int m */
{
  int pino[maxm];
  /* pino<-matrix(0,m,1) 
     pinoon laitetaan aina jos koskettaa, max kosketuksia m */
  int pinind;
  /* pinossa viitataan rindeksin elementteihin */
  int i, j, k, curleima; 
  /* i ja j viittavat rindeksit-vektoriin, jonka alkiot viittavat atomeihin */
  int merkatut[maxm];
  int curpallo, tulospit;
  int touch;  
  /*  resultat.tulos=malloc(m*sizeof(int));  */

  /* initialize merkatut */
  for (k=1; k<=m; k++)
    merkatut[k]=0;

  curleima=1;
  pinind=0;
  i=1;
  while (i<=m){
    if (merkatut[i]==0){  /* jos ei viela merkattu niin pannaan pinoon */ 
      pinind=pinind+1;    
      pino[pinind]=i;    
      while (pinind>0){
	curpallo=pino[pinind];  
        pinind=pinind-1;  
         /* otetaan pinosta viite rindeksit-vektoriin 
            jossa puolestaan viitteet itse palloihin  */
        res[curleima][curpallo]=1;    
        /* laitetaan pallo ko tasojoukkoon */
        j=1;
        while (j<=m){        /* pannnaan linkeista pinoon */
	   /* kaydaan ko tasojoukon atomit lapi */
           touch=(linkit[curpallo][j]==1);
           if ((touch) && (merkatut[j]==0)){
             pinind=pinind+1;      
             pino[pinind]=j;  
             merkatut[j]=1;
           }
           j=j+1;
        }
      }
      curleima=curleima+1;   /* uusi leima */
    }
    i=i+1;
  }
  tulospit=curleima-1;
  return tulospit;
}












int joinSets(int leftbeg,
             int rightbeg,
             int sepnum)
{
  int tavoite, hiihtaja, nykyinen, i, j, k;
  int osoittajaS[maxm];  /* make vector of pointers to the begs of sets */
  int osoittajaNewBleft[maxm];
  int osoittajaNewBright[maxm];
  int sol, len, laskuri, TotalBeg, curre, kL, kR;

  TotalBeg=separy[leftbeg];

 tavoite=1;
 hiihtaja=TotalBeg;
 while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
   hiihtaja=begsSepaNext[hiihtaja];
   tavoite=tavoite+1;
 }  
 if (tavoite<sepnum){ /* now hiihtaja points to the end of the first list */
   /* join the lists */
   begsSepaNext[hiihtaja]=separy[rightbeg];
   /* we continue */
   hiihtaja=separy[rightbeg];
   tavoite=tavoite+1;
   while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
     hiihtaja=begsSepaNext[hiihtaja];
     tavoite=tavoite+1;
   }    
   begsSepaNext[hiihtaja]=0;
}
 else{  /* we have reached goal, cut without joining */
   begsSepaNext[hiihtaja]=0;
}

 nykyinen=TotalBeg;
 i=1;
 while (i<=sepnum){
   /* len=sum(res[i,]) number of sets to be joined <= m */
   len=0;
   sol=1;
   while (sol<=m){
     if (res[i][sol]==1) len=len+1;   
     sol=sol+1;
   }
   /* we find vectors which contain pointer to the beginnings */
   /* of lists of atoms */
  
   laskuri=1;
  for (j=1; j<=m; ++j){
     if (res[i][j]==1){
       osoittajaS[laskuri]=startpointsS[j];  
       osoittajaNewBleft[laskuri]=startpointsNewBleft[j];   /* could be 0 */
       osoittajaNewBright[laskuri]=startpointsNewBright[j];  /* could be 0 */
       laskuri=laskuri+1;
     }    
  }
  
  /* handle separy */ 
  
  begsSepaBegs[nykyinen]=osoittajaS[1];    /* always non-zero */
  
  k=1;
  while (k<=(len-1)){    
    curre=osoittajaS[k];
    while (atomsSepaNext[curre]>0){    /* find the end */
      curre=atomsSepaNext[curre];
    }
    atomsSepaNext[curre]=osoittajaS[k+1];
    k=k+1;
  }
  
  /* handle left boundary */
  
  /* set kL=0 if all 0 , otherwise kL first nonzero */
  
  k=1;
  while ((k<=len) && (osoittajaNewBleft[k]==0)){
    k=k+1;
  }
  if (k>len){   /* all zero */
    kL=0;
    begsLeftBoun[nykyinen]=0;
  }
  else{         /* kL is first non zero */
    kL=k;
    begsLeftBoun[nykyinen]=osoittajaNewBleft[kL];
  
    /* update the list of left boundaries */
    /* concatenate the lists of atoms */
  
    k=kL;
  while (k<=(len-1)){    
    curre=osoittajaNewBleft[k];         /* curre is not zero */
      while (atomsLBounNext[curre]>0){    /*find the end */
	curre=atomsLBounNext[curre];
      }
      /* find the next non zero */
      k=k+1;
      while ((k<=len) && (osoittajaNewBleft[k]==0)){
	k=k+1;
      }
      if (k>len){
	atomsLBounNext[curre]=0;
      }
      else{  /* found nonzero */
	atomsLBounNext[curre]=osoittajaNewBleft[k];
      }
  }
  }
  
  /* handle right boundary */
  
  /* set kR=0 if all 0 , otherwise kR first nonzero */
  
  k=1;
  while ((k<=len) && (osoittajaNewBright[k]==0)){
    k=k+1;
  }
  if (k>len){
    kR=0;
    begsRighBoun[nykyinen]=0;
  }
  else{
    kR=k;
    begsRighBoun[nykyinen]=osoittajaNewBright[kR];
  
    /* update the list of right boundaries */
    /* concatenate the lists of atoms */
  
    k=kR;
  while (k<=(len-1)){    
    curre=osoittajaNewBright[k];         /* curre is not zero */
      while (atomsRBounNext[curre]>0){    /* find the end */
	curre=atomsRBounNext[curre];
      }
      /* find the next non zero */
      k=k+1;
      while ((k<=len) && (osoittajaNewBright[k]==0)){
	k=k+1;
      }
      if (k>len){
	atomsRBounNext[curre]=0;
      }
      else{  /* found nonzero */
	atomsRBounNext[curre]=osoittajaNewBright[k];
      }
  }
  }
  
  /* we move to the next sepaset */
  nykyinen=begsSepaNext[nykyinen];
  i=i+1;
}
 return TotalBeg; 
}






















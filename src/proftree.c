/* 
R CMD SHLIB -o ~/denfor/proftree ~/denpro/src/proftree.c 
dyn.load("~/denfor/proftree") 
dt<-.C("proftreeC",
               as.integer(numofall),
               as.integer(atomnum),
               as.double(inlevseq),
               as.integer(d),         
               as.integer(Q),  
               as.double(instep),
               as.double(insuppo), 
               as.integer(nodenumOfTree),
               as.integer(inleft),
               as.integer(inright),
               as.integer(inparent), 
               as.integer(inval),
               as.integer(invec), 
               as.integer(ininfopointer), 
               as.integer(inlowtr),
               as.integer(inupptr),
               as.double(invalue),
               as.integer(inlow),
               as.integer(inupp),
               as.integer(innodefinder),
               level = double(numofall+1),
               parent = integer(numofall+1),
               volume = double(numofall+1),
               center = double(d*numofall+1),
               efek = integer(1))


lapu1 = double(1),
lapu2 = double(1),
lapu3 = double(1),
lapu4 = double(1),
lapu5 = double(1))
              
               component = integer(numofall+1),               
               AtomlistAtomOut = integer(numofall+1),
               AtomlistNextOut = integer(numofall+1),
               numOfAtoms = integer(1))

out = integer(101),

*/

#define numofallIntern 1000000 /*yksi atomi voi esiintya useassa kp:ntissa*/
#define atomnumIntern 1000000
#define nodenumOfDyakerIntern 1000000
#define maxdim 5

int AtomlistNexty[numofallIntern+1];
int AtomlistAtomy[numofallIntern+1]; 
/* point to value,..: values in 1,...,atomnum */

int begsSepaNexty[atomnumIntern+1];
int begsSepaBegsy[atomnumIntern+1];
int atomsSepaNexty[atomnumIntern+1];
int atomsSepaAtomy[atomnumIntern+1];

int begsLeftBouny[atomnumIntern+1];
int begsRighBouny[atomnumIntern+1];
int atomsLBounNexty[atomnumIntern+1];
int atomsRBounNexty[atomnumIntern+1];

int separyy[nodenumOfDyakerIntern+1];

int low[atomnumIntern+1][maxdim+1];
int upp[atomnumIntern+1][maxdim+1];

/*
int globapu2=0, globapu3=0;
double apsu1;
double apsu2;
double apsu3;
double apsu4;
double apsu5;
double apsu[numofallIntern+1];
int tepu;
*/

#include <math.h>
#include <stdlib.h> 


int declevgenCC(int *left,
                int *right,
                int *parent,
                int *val,
                int *vec,
                int *infopointer,
                int *nodefinder,
                int *leafloc,
                int *tobehandled,
                int nodenumOfDyaker,
	        int terminalnum,
                int d);
int listchangeCCC(int totbegSeparyy,
                  int beg,
                  /*output*/
                  int *sepalkmi,
                  int *begs);
int findleafs(int *left,
              int *right,
              int itemnum,
              /*output*/
              int *leafloc);
int joincongenC(int leftbeg,
                int rightbeg,
                int direction,
                int d);
int declevnewCCC(int m, int *inlinkit,
                 /*output*/ 
                 int *tulospit, 
                 int *outtulos);
int startpointsC(int leftbeg, int rightbeg,
                 /*output*/
                 int *startpointsS, int *startpointsB,
		 int *startpointsNewBleft, int *startpointsNewBright,
                 int *mleftout, int *mrightout);
int makelinksC(int mleft, int mright, int direction, int d,
              int *startpointsB, 
              /*output*/
              int *outlinkit);
int joinsetsC(int leftbeg, int rightbeg, int sepnum, int m,
              int *inres,
              int *startpointsS, 
	      int *startpointsNewBleft, int *startpointsNewBright,
              /* output */             
              int *Totalbegi);




void proftreeC(int *numofall, 
               int *atomnum,
               double *levseq,
               int *d,       /* d<-length(N) */
	       int *levnum,  /* levnum<-length(levseq) */
               double *step,
               double *suppo,
               /* tree structures */
               int *nodenumOfDyaker,
               int *left,
               int *right,
               int *parentKertree,
               int *val,
               int *vec,
               int *infopointer,
               int *inlowtr,
               int *inupptr,
               /* associated info */
               double *value,
               int *inlow,
               int *inupp,
               int *nodefinder,
               /* output */
	       double *level,      /* of length numofall+1 */
               int *parent, 
               double *volume,
               double *center,
               int *efek)    /* gives effective length of level,... */

     /*
double *lapu1,
double *lapu2,
double *lapu3,
double *lapu4,
double *lapu5)
     */

               /*
               int *component) 
               int *AtomlistAtomOut,
               int *AtomlistNextOut,
	       int *numOfAtoms)
	       */
              

{
  int apu, i, j, ij, k;

  int totbegSepary;

  /*
  int component[*numofall+1];
  int pinoComponent[*numofall+1];  pointer to component, level,... 
  int pinoTaso[*numofall+1];       ordinal of level (pointer to levseq) 
  int tobehandled[*nodenumOfDyaker+1], leafloc[*nodenumOfDyaker+1];
  */
  int *component = (int *)malloc(sizeof(int) * (*numofall+1));
  int *pinoComponent = (int *)malloc(sizeof(int) * (*numofall+1));
  int *pinoTaso = (int *)malloc(sizeof(int) * (*numofall+1));
  int *tobehandled = (int *)malloc(sizeof(int) * (*nodenumOfDyaker+1));
  int *leafloc = (int *)malloc(sizeof(int) * (*nodenumOfDyaker+1));

  int componum, beghigher, beg, koko, pinind, ind, levind, partlevsetbeg;
  int listEnd, PrevlistEnd;
  int terminalnum, addnum, removenum, runner, origiListEnd, atom;
  double arvo, higlev;
  int lokalefek;
  int begi, atto, node, curterminalnum;
 
  /*volume and center calculation*/
  double curvolu, vol, ala, yla;
  int componentnum, zeiger; 
  int atompointer, uppi, lowi;
  /*double curcente[*d+1], newcente[*d+1];*/ 
  double *curcente = (double *)malloc(sizeof(double) * (*d+1));
  double *newcente = (double *)malloc(sizeof(double) * (*d+1));
 
  /*listchangeCCC*/
  int kokoout[2], begs[atomnumIntern+1];

  if (component == NULL) exit(1); 
  if (pinoComponent == NULL) exit(1); 
  if (pinoTaso == NULL) exit(1); 
  if (tobehandled == NULL) exit(1); 
  if (leafloc == NULL) exit(1);  
  if (curcente == NULL) exit(1);  
  if (newcente == NULL) exit(1);  

  /* Initialize upp, low */

  for (i=1; i<=*atomnum; i++){
      for (j=1; j<=*d; j++){
           low[i][j]=inlow[(i-1)**d+j];
           upp[i][j]=inupp[(i-1)**d+j];
      }
  }

  /* initialize tobehandled */
 
  for (i=1; i<=(*nodenumOfDyaker); i++) tobehandled[i]=0;

  /* definition */
  
  componum=*numofall;

  /* leaflocs */ 
  apu=findleafs(left,right,*nodenumOfDyaker,leafloc);/*returns leafloc*/

  /* Initilize the lists */
  i=1;
  while (i<=*atomnum){
    AtomlistAtomy[i]=i;
    AtomlistNexty[i]=i+1;
    i=i+1;
  }
  listEnd=*atomnum;
  AtomlistNexty[listEnd]=0;

  /* Let us divide the lowest level set to disconnected parts */

  beg=1;
  terminalnum=*atomnum;

  begi=1;
  atto=AtomlistAtomy[begi];
  while (begi>0){
     if (value[atto]>0){
        node=nodefinder[atto];
        tobehandled[node]=1;
     }    
     begi=AtomlistNexty[begi];
     atto=AtomlistAtomy[begi];
  }
  curterminalnum=terminalnum;


  totbegSepary=declevgenCC(left,
                           right,
                           parentKertree,
                           val,
                           vec,
                           infopointer,
                           /* index, */
                           nodefinder,
                           leafloc,
                           tobehandled,
                           *nodenumOfDyaker,
                           terminalnum,
                           *d);
  /* tobehandled,val,vec,infopointer,parent,low,upp */
 
  /* modify begs, at most "terminalnum" changes */  
  apu=listchangeCCC(totbegSepary,beg,kokoout,begs);    
  koko=kokoout[1];

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
      atom=AtomlistAtomy[runner];
      arvo=value[atom];
      if (arvo>=higlev){
	  listEnd=listEnd+1;    
	  AtomlistAtomy[listEnd]=atom;
	  AtomlistNexty[listEnd]=listEnd+1; 
	  addnum=addnum+1;     
      }
      else{           
          removenum=removenum+1;
      }                                
      runner=AtomlistNexty[runner]; 
    }
    AtomlistNexty[listEnd]=0;      /* we have to correct the end to terminate */

    if (addnum>0){
      AtomlistNexty[PrevlistEnd]=0;
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

      begi=beghigher;
      for (ij=1; ij<=*nodenumOfDyaker; ij++) tobehandled[ij]=0;
      while (begi>0){
          atto=AtomlistAtomy[begi];
          node=nodefinder[atto];
          tobehandled[node]=1;
          begi=AtomlistNexty[begi];
      }
      curterminalnum=addnum;
      terminalnum=addnum;


      /*globapu2=globapu2+1;*/

      totbegSepary=declevgenCC(left,
                               right,
                               parentKertree,
                               val,
                               vec,
                               infopointer,
                               /* index, */
                               nodefinder,
                               leafloc,
                               tobehandled,
                               *nodenumOfDyaker,
                               curterminalnum,
                               *d); 

      /*
      if (globapu2==1){
 for (j=1; j<=100; j++){
   out[j]=atomsSepaNexty[j];
 }
      }
      */

      /* modify begs, at most "terminalnum" changes */  
      apu=listchangeCCC(totbegSepary,beg,kokoout,begs);    
      koko=kokoout[1];

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

/* volumes and centers */

for (i=1; i<=componentnum; i++){
  curvolu=0;
  zeiger=component[i];
  while (zeiger>0){
     atto=AtomlistAtomy[zeiger];
     vol=1;
     for (j=1; j<=*d; j++){
        lowi=inlow[(atto-1)*(*d)+j];  /*[atto][j]*/
        uppi=inupp[(atto-1)*(*d)+j];  /*[atto][j]*/
        vol=vol*(uppi-lowi)*step[j];
     }
     curvolu=curvolu+vol;
     zeiger=AtomlistNexty[zeiger];
  }
  volume[i]=curvolu;
}

for (i=1; i<=componentnum; i++){  
  for (j=1; j<=*d; j++) curcente[j]=0;
  zeiger=component[i];
  while (zeiger>0){
     atompointer=AtomlistAtomy[zeiger];
     for (j=1; j<=*d; j++) newcente[j]=0;
     for (j=1; j<=*d; j++){
        /*calculate 1st volume of d-1 dimensional rectangle where*/
        /*we have removed j:th dimension*/
        vol=1;
        k=1;
        while (k<=*d){
           if (k!=j){
              lowi=inlow[(atompointer-1)*(*d)+k];  /*[atompointer][k]*/
              uppi=inupp[(atompointer-1)*(*d)+k];  /*[atompointer][k]*/
              vol=vol*(uppi-lowi)*step[k];
	   }
           k=k+1;
        }
        lowi=inlow[(atompointer-1)*(*d)+j];  /*[atompointer][j]*/
        uppi=inupp[(atompointer-1)*(*d)+j];  /*[atompointer][j]*/
        ala=suppo[2*j-1]+step[j]*lowi;
        yla=suppo[2*j-1]+step[j]*uppi;
        newcente[j]=vol*(pow(yla,2)-pow(ala,2))/2;
     }
     for (j=1; j<=*d; j++) curcente[j]=curcente[j]+newcente[j];
     zeiger=AtomlistNexty[zeiger];
  }
  for (j=1; j<=*d; j++) center[(i-1)*(*d)+j]=curcente[j]/volume[i];

}

 /*
 *lapu1=apsu1;
 *lapu2=apsu2;
 *lapu3=apsu3;
 *lapu4=apsu4;
 *lapu5=apsu5;
 */

 /* 
 *numOfAtoms=listEnd;
 
 for (j=1; j<=listEnd; j++){
     AtomlistAtomOut[j]=AtomlistAtomy[j];
     AtomlistNextOut[j]=AtomlistNext[j];
 }
 */
 
  free(component); 
  free(pinoComponent); 
  free(pinoTaso); 
  free(tobehandled); 
  free(leafloc);  
  free(curcente);  
  free(newcente);  

}


/****************************** END OF MAIN ****************************/



int declevgenCC(int *left,
                int *right,
                int *parent,
                int *val,
                int *vec,
                int *infopointer,
                /* int *index[], */
                int *nodefinder,
                int *leafloc,
                int *tobehandled,
                /* end tree */
                int nodenumOfDyaker,
	        int terminalnum,
                int d)

{
    /*
 int nextFloor[terminalnum+1];
 int currFloor[terminalnum+1];
 int already[nodenumOfDyaker+1];
    */
 int *nextFloor = (int *)malloc(sizeof(int) * (terminalnum+1));
 int *currFloor = (int *)malloc(sizeof(int) * (terminalnum+1));
 int *already = (int *)malloc(sizeof(int) * (nodenumOfDyaker+1));
 
 int i, lkm, nodee, note, nodelkm;
 int Lempty, Rempty;
 int leftbeg, rightbeg, direction, akku, apu3;
 int direktiooni, poiju, prevpoiju, aatto, thisnoteempty;
 double splittiini;

 if (nextFloor == NULL) exit(1);
 if (currFloor == NULL) exit(1);
 if (already == NULL) exit(1);

  /* Initialize the global variables */

  for (i=1; i<=atomnumIntern; i++){
    begsSepaNexty[i]=0;
    begsSepaBegsy[i]=0;
    atomsSepaNexty[i]=0;
    atomsSepaAtomy[i]=0;

    begsLeftBouny[i]=0;
    begsRighBouny[i]=0;
    atomsLBounNexty[i]=0;
    atomsRBounNexty[i]=0;
  }

  for (i=1; i<=nodenumOfDyakerIntern; i++) separyy[i]=0;


  nodelkm=nodenumOfDyaker;


  /* end of initializing */

  lkm=0;
  nodee=nodelkm;
  while (nodee>=1){          /* root is in position 1 */
      if ((leafloc[nodee]==1) && (tobehandled[nodee]==1)){ /* we are in leaf */

          tobehandled[parent[nodee]]=1;

          lkm=lkm+1;
          separyy[nodee]=lkm;
          atomsSepaAtomy[lkm]=infopointer[nodee];

          begsSepaBegsy[lkm]=lkm;

          /*
          obs: we need not change
          begsSepaNexty, atomsSepaNexty, atomsLBounNexty, atomsRBounNext
          since at the beginning set consist only one member:
          pointer is always 0, since we do not have followers
	  */

      }
      else if (tobehandled[nodee]==1){   /* not a leaf */

	  tobehandled[parent[nodee]]=1;

	  leftbeg=left[nodee];
	  rightbeg=right[nodee];

          if ((leftbeg==0) || (separyy[leftbeg]==0)){
	    /*if left child does not exist*/

	    /*
            note that since we consider subsets of the
            terminal nodes of the original tree, it may happen
            that leftbeg>0 but left child does not exist
	    */
           
	    separyy[nodee]=separyy[rightbeg];

            /*
            we need that all lists contain as many members
            left boundary is empty, but we will make it a list
            of empty lists
            */

          }
          else{ /* eka else */
            if ((rightbeg==0) || (separyy[rightbeg]==0)){
	   	/* right child does not exist */
    
		separyy[nodee]=separyy[leftbeg];
                      
                /*
                left boundary stays same as for leftbeg
                right boundary is empty
                */

            }
            else{   /* toka else: both children exist */
                    /* create boundaries */
                    
		direktiooni=vec[nodee];
		splittiini=val[nodee];     
               
                /*
                left boundary of right child :
                create/check whether empty 
                */

		Lempty=1;
		note=separyy[rightbeg];
                    while (note>0){
                        thisnoteempty=1;

                        poiju=begsSepaBegsy[note];
                        prevpoiju=poiju;
                        while (poiju>0){
		           aatto=atomsSepaAtomy[poiju];
                           if (splittiini>=low[aatto][direktiooni]){
                               /* this atom belongs to boundary */
                               if (thisnoteempty==1){
                                   /* poiju is the 1st non-empty */
                                   begsLeftBouny[note]=poiju;
                               }                
                               Lempty=0;
                               atomsLBounNexty[prevpoiju]=poiju;    
                               atomsLBounNexty[poiju]=0;
                               prevpoiju=poiju;

                               thisnoteempty=0;
                           }
                           poiju=atomsSepaNexty[poiju];
                        }
                        if (thisnoteempty) begsLeftBouny[note]=0;
                        note=begsSepaNexty[note];
                     }

                    /* right boundary of left child */

		    Rempty=1;
		    note=separyy[leftbeg];
                     while (note>0){
			 thisnoteempty=1;
                       
			 poiju=begsSepaBegsy[note];
			 prevpoiju=poiju;
                         while (poiju>0){
			    aatto=atomsSepaAtomy[poiju];
                            if (splittiini<=upp[aatto][direktiooni]){
                               /* this atom belongs to boundary */
                               if (thisnoteempty==1){
                                   /* poiju is the 1st non-empty */
                                   begsRighBouny[note]=poiju;
                               }                
                               Rempty=0;
                               atomsRBounNexty[prevpoiju]=poiju;    
                               atomsRBounNexty[poiju]=0;
                               prevpoiju=poiju;
                           
                               thisnoteempty=0;
                           }
                           poiju=atomsSepaNexty[poiju];
                        }
                        if (thisnoteempty) begsRighBouny[note]=0;
                        note=begsSepaNexty[note];
                     }
    
	             /* check whether one of boundaries is empty */


                     if (Lempty || Rempty){
			 /* one of boundaries is empty */
                        
			 /* concatenating separate parts */
                       
			 akku=separyy[leftbeg];
                 
			 begsSepaNexty[akku]=separyy[rightbeg]; 
                         /* concatenate list of separate sets */

			 separyy[nodee]=separyy[leftbeg];
			 akku=separyy[rightbeg];
                          
                        /* left boundaries of sets in right child are empty */
                        /* end of concatenating */



                     }
           	     else{ /*both children exist, both boundaries non-empty*/
			 direction=vec[nodee];

          
                         apu3=joincongenC(leftbeg,
                                          rightbeg,
                                          direction,
                                          d);
			 separyy[nodee]=apu3;

		     /*
                     jc=joincongen(leftbeg,rightbeg,separyy,
                     begsSepaNexty,begsSepaBegsy,begsLeftBouny,begsRighBouny,
                     atomsSepaNexty,atomsSepaAtomy,atomsLBounNexty,
                     atomsRBounNexty,
                     direction,low,upp)   
		     */

                                }
            } /* toka else */
        } /* eka else */
    }  /* else not a leaf */
      nodee=nodee-1;
}

  return separyy[1];

/*
totbegSepary=separyy[1],
begsSepaNexty=begsSepaNexty,
begsSepaBegsy=begsSepaBegsy,
atomsSepaNexty=atomsSepaNexty,
atomsSepaAtomy=atomsSepaAtomy
*/
 free(nextFloor);
 free(currFloor);
 free(already);

}




/*****************************************************************/




int listchangeCCC(int totbegSepary,
                  int beg,
                  /*output*/
                  int *sepalkmi,
                  int *begs)
{
 int runnerBegs, runnerAtoms, runnerOrigi, runnerOrigiprev, sepalkm, apu;

 runnerBegs=totbegSepary;  
 /* total beginning of list is at the root of the tree */
 runnerOrigi=beg;
 runnerOrigiprev=beg;
 sepalkm=0;
 while (runnerBegs>0){
   sepalkm=sepalkm+1;
   runnerAtoms=begsSepaBegsy[runnerBegs];
   begs[sepalkm]=runnerOrigi;
   /* first step (in order to get also runnerOrigiprev to play) */
   AtomlistAtomy[runnerOrigi]=atomsSepaAtomy[runnerAtoms];
   runnerOrigiprev=runnerOrigi;
   runnerOrigi=AtomlistNexty[runnerOrigi];
   runnerAtoms=atomsSepaNexty[runnerAtoms];

   while (runnerAtoms>0){
     AtomlistAtomy[runnerOrigi]=atomsSepaAtomy[runnerAtoms];
     runnerOrigiprev=runnerOrigi;
     runnerOrigi=AtomlistNexty[runnerOrigi];
     runnerAtoms=atomsSepaNexty[runnerAtoms];
   }
   AtomlistNexty[runnerOrigiprev]=0;  /* mark the end of the list */
   runnerBegs=begsSepaNexty[runnerBegs];
 }
 /* begs=begs[1:sepalkm] */

 sepalkmi[1]=sepalkm;

 apu=1;
 return apu;
}




/************************************************************/




/*
Finds location of leafs in binary tree, living in vector

left, right are itemnum-vectors

returns itemnum-vector, 1 in the location of nodes 0 non-terminal
-1 in positions not belonging to tree

vector where binary tree is living may be larger than cardinality
of nodes of the tree
*/

int findleafs(int *left,
              int *right,
              int itemnum,
              int *leafloc)
{

 int apu, i, pinin, cur;
 /*   int pino[itemnum+1];  */   /* itemnum=length(left) */
 int *pino = (int *)malloc(sizeof(int) * (itemnum+1));

 if (pino == NULL) exit(1);

 for (i=1; i<=itemnum; i++) leafloc[i]=-1;

 pino[1]=1;     /* pino[1]=root */
 pinin=1;
 while (pinin>0){
     cur=pino[pinin];      /* take from stack */
     pinin=pinin-1;
     if (left[cur]==0){    /* if we are in leaf */
	 leafloc[cur]=1;
     }
     else{
	 while (left[cur]>0){  /* go to leaf and put right nodes to stack */
	     leafloc[cur]=0;
	     pinin=pinin+1;
	     pino[pinin]=right[cur];
	     cur=left[cur];
         }
	 leafloc[cur]=1;  /* now we know we are in leaf */
     }
 }
 apu=1;
 return apu;

 free(pino);

}


/************************************************************/


/*
see: joincongen.R

INPUT: 
leftbeg,rightbeg,direction,d
separyy,
begsSepaNexty,begsSepaBegsy,begsLeftBouny,begsRighBouny,
atomsSepaNexty,atomsSepaAtomy,atomsLBounNexty,atomsRBounNexty,
low,upp
OUTPUT:
totbegSepary, 
begsSepaNexty, begsSepaBegsy, begsLeftBouny, begsRighBouny,
atomsSepaNexty,                atomsLBounNexty, atomsRBounNexty
*/


int joincongenC(int leftbeg,
                int rightbeg, 
                int direction, 
                int d)     
{

  int apu, i;
  int induksiInt=20;    /* length(begsSepaNexty) */
  int componumInt=20;   /* length(begsSepaNexty) */

  /*startpointsC*/
  /*
  int startpointsS[induksiInt+1];
  int startpointsB[induksiInt+1];
  int startpointsNewBleft[induksiInt+1];
  int startpointsNewBright[induksiInt+1];
  */
  int *startpointsS = (int *)malloc(sizeof(int) * (induksiInt+1));
  int *startpointsB = (int *)malloc(sizeof(int) * (induksiInt+1));
  int *startpointsNewBleft = (int *)malloc(sizeof(int) * (induksiInt+1));
  int *startpointsNewBright = (int *)malloc(sizeof(int) * (induksiInt+1));

  int m, mleft, mright;
  int mleftout[2], mrightout[2];

  /*makelinksC*/
  /*int linkit[componumInt+1][componumInt+1];*/
  /*
  int outlinkit[1+componumInt*componumInt];
  int outtulos[componumInt+1];
  */
  int *outlinkit = (int *)malloc(sizeof(int) * (1+componumInt*componumInt));
  int *outtulos = (int *)malloc(sizeof(int) * (componumInt+1));

  /*declevnewCCC*/
 int sepnum, tulospit[2];
  /*int res[componumInt+1][componumInt+1];*/
  /*int inlinkit[1+componumInt*componumInt];*/

  /*joinsetsC*/
  int Totalbeg, Totalbegout[2];

  if (startpointsS == NULL) exit(1);
  if (startpointsB == NULL) exit(1);
  if (startpointsNewBleft == NULL) exit(1);
  if (startpointsNewBright == NULL) exit(1);
  if (outlinkit == NULL) exit(1);
  if (outtulos == NULL) exit(1);

  /* 1. new boundary: left bound. of left child, right b. of right child */

  for (i=1; i<=induksiInt; i++){
     startpointsS[i]=0;
     startpointsB[i]=0;
     startpointsNewBleft[i]=0;
     startpointsNewBright[i]=0;
  }

  apu=startpointsC(leftbeg, rightbeg,
                   /*output*/
                   startpointsS, startpointsB,
	       	   startpointsNewBleft, startpointsNewBright,
		   mleftout, mrightout);

  mleft=mleftout[1];
  mright=mrightout[1];
  m=mleft+mright;


  /* 2. We make "linkit" matrix */


  apu=makelinksC(mleft, mright, direction, d, startpointsB,
                 /*output*/ 
                 outlinkit);
  /*
  for (i=1; i<=m; i++){
    for (j=1; j<=m; j++){
       linkit[i][j]=outlinkit[(i-1)*m+j];
    }
  }
  */

  /* 3. We apply declev */
 
  /*
  for (i=1; i<=m; i++){
     for (j=1; j<=m; j++){
       inlinkit[(i-1)*m+j]=linkit[i][j];
     }
  }
  */
  apu=declevnewCCC(m, outlinkit, 
                      /*output*/
                      tulospit,
                      outtulos);
  sepnum=tulospit[1];
  /* 
  for (i=1; i<=m; i++){
     for (j=1; j<=m; j++){
       res[i][j]=outtulos[(i-1)*m+j];
     }
  }
  res is sepnum*m-matrix of 0/1
  sepnum=dim(res)[1]
  */
 
  /* 4. We join the sets */

  apu=joinsetsC(leftbeg, rightbeg, sepnum, m,
                outtulos, /*inres,*/
                startpointsS, 
	        startpointsNewBleft, startpointsNewBright,
                /* output */             
                Totalbegout);

  Totalbeg=Totalbegout[1];

  return Totalbeg;

  free(startpointsS);
  free(startpointsB);
  free(startpointsNewBleft);
  free(startpointsNewBright);
  free(outlinkit);
  free(outtulos);

}



/*******************************************************/



/* linkit m*m-matrix */
/* return m*m-matrix tulos and integer tulospit */
/* first tulospit rows contain information */

int declevnewCCC(int m, int *inlinkit,
                 /*output*/ 
                 int *tulospit, 
                 int *outtulos)  
{
  int pinind; 
    /*pinoon laitetaan aina jos koskettaa, max kosketuksia m*/
  int i, j, k, curleima, curpallo, touch;  
  int apu;
  /*
  int pino[m+1], merkatut[m+1]; 
  int linksit[m+1][m+1], tulos[m+1][m+1];
  */
 int *pino = (int *)malloc(sizeof(int) * (m+1));
 int *merkatut = (int *)malloc(sizeof(int) * (m+1));
 int ** linksit;
 int ** tulos;
 linksit = (int **)malloc((m+1) * sizeof(int *));
 tulos = (int **)malloc((m+1) * sizeof(int *));
 if (NULL == linksit) exit(1);
 for (i = 0; i <= m; i++) {
      linksit[i] = (int *)malloc((m+1) * sizeof(int));
      if (NULL == linksit[i]) exit(1);
 }
 if (NULL == tulos) exit(1);
 for (i = 0; i <= m; i++) {
      tulos[i] = (int *)malloc((m+1) * sizeof(int));
      if (NULL == tulos[i]) exit(1);
 }
 if (pino == NULL) exit(1); 
 if (merkatut == NULL) exit(1); 

  /* initialize */
  for (k=1; k<=m; k++) merkatut[k]=0;
 
  for (i=1; i<=m; i++){
    for (j=1; j<=m; j++){
      linksit[i][j]=inlinkit[(i-1)*m+j];
      tulos[i][j]=0;
    }
  }

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
        tulos[curleima][curpallo]=1;    
        /* laitetaan pallo ko tasojoukkoon */
        j=1;
        while (j<=m){        /* pannnaan linkeista pinoon */
	   /* kaydaan ko tasojoukon atomit lapi */
           touch=(linksit[curpallo][j]==1);
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

  tulospit[1]=curleima-1;
  for (i=1; i<=m; i++){
    for (j=1; j<=m; j++){
       outtulos[(i-1)*m+j]=tulos[i][j];
    }
  }

  apu=curleima-1;
  return apu;

 free(pino); 
 free(merkatut); 
 for(i = 0; i <= m; i++) free(linksit[i]);
 free(linksit);
 for(i = 0; i <= m; i++) free(tulos[i]);
 free(tulos);

}


/*******************************************************/


int startpointsC(int leftbeg, int rightbeg,
                 /*output*/
                 int *startpointsS, int *startpointsB,
		 int *startpointsNewBleft, int *startpointsNewBright,
                 int *mleftout, int *mrightout)
{

  int apu;
  int induksi, anfang, m, ml, mr;

  induksi=1;
  anfang=separyy[leftbeg];
  startpointsS[induksi]=anfang;
  startpointsB[induksi]=begsRighBouny[anfang];
  startpointsNewBleft[induksi]=begsLeftBouny[anfang];
  while (begsSepaNexty[anfang]>0){
    anfang=begsSepaNexty[anfang];
    induksi=induksi+1;
    startpointsS[induksi]=begsSepaBegsy[anfang];
    startpointsB[induksi]=begsRighBouny[anfang];
    startpointsNewBleft[induksi]=begsLeftBouny[anfang];  
  }
  ml=induksi;
  induksi=induksi+1;
  anfang=separyy[rightbeg];
  startpointsS[induksi]=anfang;
  startpointsB[induksi]=begsLeftBouny[anfang];
  startpointsNewBright[induksi]=begsRighBouny[anfang];
  while (begsSepaNexty[anfang]>0){
    anfang=begsSepaNexty[anfang];
    induksi=induksi+1;
    startpointsS[induksi]=begsSepaBegsy[anfang];
    startpointsB[induksi]=begsLeftBouny[anfang];
    startpointsNewBright[induksi]=begsRighBouny[anfang];  
  }
  m=induksi;
  mr=m-ml;  

  mleftout[1]=ml;
  mrightout[1]=mr;
 
  apu=1;
  return apu;

}


/**************************************************************/


int makelinksC(int mleft, int mright, int direction, int d,
              int *startpointsB, 
              /*output*/
              int *outlinkit)

{

  int apu, i, j, k;
  int dodo, re, beg1, beg2, conne;
  /*int indelow1[d+1], indelow2[d+1], indeupp1[d+1], indeupp2[d+1];*/
  int *indelow1 = (int *)malloc(sizeof(int) * (d+1));
  int *indelow2 = (int *)malloc(sizeof(int) * (d+1));
  int *indeupp1 = (int *)malloc(sizeof(int) * (d+1));
  int *indeupp2 = (int *)malloc(sizeof(int) * (d+1));

  int begbeg1, begbeg2, atom1, atom2, touch;

  int m=mleft+mright;
  /*int linkit[m+1][m+1]; */
  int ** linkit;
  linkit = (int **)malloc((m+1) * sizeof(int *));
  if (NULL == linkit) exit(1);
  for (i = 0; i <= m; i++) {
      linkit[i] = (int *)malloc((m+1) * sizeof(int));
      if (NULL == linkit[i]) exit(1);
  }

  if (indelow1 == NULL) exit(1);
  if (indelow2 == NULL) exit(1);
  if (indeupp1 == NULL) exit(1);
  if (indeupp2 == NULL) exit(1);

  for (i=1; i<=m; i++){
     for (j=1; j<=m; j++){
       linkit[i][j]=0;
     }
   }

  dodo=1;
  while (dodo <= mleft){
      beg1=startpointsB[dodo];    /* could be 0 */
      re=mleft+1;
      while (re <= m){
	  beg2=startpointsB[re];    /* could be 0 */
	  conne=0;
	  begbeg1=beg1;
          while (begbeg1>0){
	      begbeg2=beg2;
              while (begbeg2>0){
		  atom1=atomsSepaAtomy[begbeg1];
		  atom2=atomsSepaAtomy[begbeg2];
		  for (k=1; k<=d; k++){
                     indelow1[k]=low[atom1][k];
		     indeupp1[k]=upp[atom1][k];
         	     indelow2[k]=low[atom2][k];
		     indeupp2[k]=upp[atom2][k];
		  }
	          touch=1;  /* TRUE */
	          i=1;
                  while (i <= d){
                    if ((i != direction) &&
                    ((indelow1[i]>indeupp2[i]) || (indeupp1[i]<indelow2[i]))){
	               touch=0;  /* FALSE */
                    }
	            i=i+1;
                  }                
                  if (touch) conne=1;   /* TRUE */
                  begbeg2=atomsLBounNexty[begbeg2];
              }
              begbeg1=atomsRBounNexty[begbeg1];
          }                
          if (conne) linkit[dodo][re]=1;
          re=re+1;
      }
      dodo=dodo+1;
  }
  for (dodo=(mleft+1); dodo<=m; dodo++){
      beg1=startpointsB[dodo];
      for (re=1; re<=mleft; re++){
	  beg2=startpointsB[re];
	  conne=0;
	  begbeg1=beg1;
          while (begbeg1>0){
	     begbeg2=beg2;
             while (begbeg2>0){
		 atom1=atomsSepaAtomy[begbeg1];
		 atom2=atomsSepaAtomy[begbeg2];
                 for (k=1; k<=d; k++){
                     indelow1[k]=low[atom1][k];
		     indeupp1[k]=upp[atom1][k];
         	     indelow2[k]=low[atom2][k];
		     indeupp2[k]=upp[atom2][k];
		 }
                 touch=1;  /* TRUE */
	         i=1;
                 while (i <= d){
                    if ((i != direction) &&
                    ((indelow1[i]>indeupp2[i]) || (indeupp1[i]<indelow2[i]))){
	               touch=0;  /* FALSE */
                    }
	            i=i+1;
                 }        
                 if (touch) conne=1;   /* TRUE */
                 begbeg2=atomsRBounNexty[begbeg2];
             }
	     begbeg1=atomsLBounNexty[begbeg1];
          }                
          if (conne) linkit[dodo][re]=1;
      }
  }

  for (i=1; i<=m; i++){
    for (j=1; j<=m; j++){
      outlinkit[(i-1)*m+j]=linkit[i][j];
    }
  }

  apu=1;
  return apu;

  free(indelow1);
  free(indelow2);
  free(indeupp1);
  free(indeupp2);
  for(i = 0; i <= m; i++) free(linkit[i]);
  free(linkit);


}



/*********************************************************/

/*
OutPut:

begsSepaNexty, begsSepaBegsy, begsLeftBouny, begsRighBouny,
atomsSepaNexty, atomsLBounNexty, atomsRBounNexty
*/


int joinsetsC(int leftbeg, int rightbeg, int sepnum, int m,
              int *inres,
              int *startpointsS, 
	      int *startpointsNewBleft, int *startpointsNewBright,
              /* output */             
              int *Totalbegi) 
{

  int apu, i, j, ii, jj, k, ll;

  int len;   /* len=sum(res[i,]) number of sets to be joined */
  int Totalbeg, tavoite, hiihtaja, nykyinen;
  int laskuri, curre, kL, kR;
 
 /* m=induksiInt */
  /*
  int osoittajaS[m+1];  make vector of pointers to the begs of sets 
  int osoittajaNewBleft[m+1];
  int osoittajaNewBright[m+1];
  int res[sepnum][m];
  */
  int *osoittajaS= (int *)malloc(sizeof(int) * (m+1));
  int *osoittajaNewBleft= (int *)malloc(sizeof(int) * (m+1));
  int *osoittajaNewBright= (int *)malloc(sizeof(int) * (m+1));
  int ** res;
  res = (int **)malloc((sepnum+1) * sizeof(int *));
  if (NULL == res) exit(1);
  for (i = 0; i <= m; i++) {
      res[i] = (int *)malloc((m+1) * sizeof(int));
      if (NULL == res[i]) exit(1);
  }
  if (osoittajaS == NULL) exit(1);
  if (osoittajaNewBleft == NULL) exit(1);
  if (osoittajaNewBright == NULL) exit(1);

  for (i=1; i<=sepnum; i++){
    for (j=1; j<=m; j++){
       res[i][j]=inres[(i-1)*m+j];
    }
  } 

  for (i=1; i<=m; i++){
     osoittajaS[i]=0;
     osoittajaNewBleft[i]=0;
     osoittajaNewBright[i]=0;
  }

 Totalbeg=separyy[leftbeg];

 tavoite=1;
 hiihtaja=Totalbeg;
 while ((begsSepaNexty[hiihtaja]>0) && (tavoite<sepnum)){
     hiihtaja=begsSepaNexty[hiihtaja];
     tavoite=tavoite+1;
 }  
 if (tavoite<sepnum){  /* now hiihtaja points to the end of the first list */
	               /* join the lists */
     begsSepaNexty[hiihtaja]=separyy[rightbeg];
     /* we continue */
     hiihtaja=separyy[rightbeg];
     tavoite=tavoite+1;
     while ((begsSepaNexty[hiihtaja]>0) && (tavoite<sepnum)){
	 hiihtaja=begsSepaNexty[hiihtaja];
	 tavoite=tavoite+1;
     }    
     begsSepaNexty[hiihtaja]=0;
 }
 else{  /* we have reached goal, cut without joining */
     begsSepaNexty[hiihtaja]=0;
 }


 nykyinen=Totalbeg;
 ii=1;
 while (ii<= sepnum){
  
     len=0;
     for (ll=1; ll<=m; ll++){
        if (res[ii][ll]==1) len=len+1;
     }

     /* we find vectors which contain pointer to the beginnings */
	/* of lists of atoms */
  
     laskuri=1;
     for (jj=1; jj<=m; jj++){
         if (res[ii][jj]==1){
	     osoittajaS[laskuri]=startpointsS[jj];  
	     osoittajaNewBleft[laskuri]=startpointsNewBleft[jj];/*could be 0*/
             osoittajaNewBright[laskuri]=startpointsNewBright[jj];/*couldbe0*/
	     laskuri=laskuri+1;
         }    
     }
  
     /* handle separyy */
  
     begsSepaBegsy[nykyinen]=osoittajaS[1];   /* always non-zero */
     
     k=1;
     while (k<=(len-1)){    
	 curre=osoittajaS[k];
         while (atomsSepaNexty[curre]>0){ /*find the end*/
	      curre=atomsSepaNexty[curre];
         }
         atomsSepaNexty[curre]=osoittajaS[k+1];
	 k=k+1;
     }

/*
if ((atomsSepaNexty[12]==11) && (globapu2==1)) globapu3=globapu3+1;
if ((leftbeg==2) && (rightbeg==3) && (globapu2==1)){
       apsu1=atomsSepaNexty[11];
       apsu2=atomsSepaNexty[12];
       apsu3=atomsSepaNexty[13];
       apsu4=atomsSepaNexty[14];
       apsu5=atomsSepaNexty[15];
}
*/
     /* handle left boundary */
  
     /* set kL=0 if all 0 , otherwise kL first nonzero */
  
     k=1;
     while ((k<=len) && (osoittajaNewBleft[k]==0)){
	 k=k+1;
     }
     if (k>len){   /* all zero */
	 kL=0;
	 begsLeftBouny[nykyinen]=0;
     }
     else{        /*  kL is first non zero  */
	 kL=k;
	 begsLeftBouny[nykyinen]=osoittajaNewBleft[kL];
	 /*
         update the list of left boundaries
         concatenate the lists of atoms
	 */
	 k=kL;
         while (k<=(len-1)){    
	     curre=osoittajaNewBleft[k]; /* curre is not zero */
             while (atomsLBounNexty[curre]>0){    /* find the end */
		 curre=atomsLBounNexty[curre];
             }
	     /* find the next non zero */
	     k=k+1;
             while ((k<=len) && (osoittajaNewBleft[k]==0)){
		 k=k+1;
             }
           if (k>len){
	       atomsLBounNexty[curre]=0;
           }
           else{  /*  found nonzero */
	       atomsLBounNexty[curre]=osoittajaNewBleft[k];
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
	 begsRighBouny[nykyinen]=0;
     }
     else{
	 kR=k;
	 begsRighBouny[nykyinen]=osoittajaNewBright[kR];
  
	 /* update the list of right boundaries */
	 /* concatenate the lists of atoms */
  
	 k=kR;
         while (k<=(len-1)){    
	     curre=osoittajaNewBright[k];        /* curre is not zero */
             while (atomsRBounNexty[curre]>0){    /* find the end */
		 curre=atomsRBounNexty[curre];
             }
	     /* find the next non zero */
	     k=k+1;
             while ((k<=len) && (osoittajaNewBright[k]==0)){
		 k=k+1;
             }
             if (k>len){
		 atomsRBounNexty[curre]=0;
             }
	     else{  /* found nonzero */
		 atomsRBounNexty[curre]=osoittajaNewBright[k];
             }
         }
     }
  
     /* we move to the next sepaset */
     nykyinen=begsSepaNexty[nykyinen];
     ii=ii+1;
 }

 Totalbegi[1]=Totalbeg;

 apu=1;
 return apu;

 for(i = 0; i <= sepnum; i++) free(res[i]);
 free(res);
 free(osoittajaS);
 free(osoittajaNewBleft);
 free(osoittajaNewBright);

}


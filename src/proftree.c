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

#define maxn 1001
#define maxlstnode 100000
#define maxkernelnode 100000

#define induksiInt 20    /* length(begsSepaNext) */
#define componumInt 20   /* length(begsSepaNext) */
#define mmm 20           /* in declevnewCCC, joinsetsC */
#define sepnummm 40      /* in joinsetsC */

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
int listchangeCCC(int totbegSepary,
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

  int component[maxlstnode+1];
  int pinoComponent[maxlstnode+1];  /* pointer to component, level,... */
  int pinoTaso[maxlstnode+1];       /* ordinal of level (pointer to levseq) */

  int componum, beghigher, beg, koko, pinind, ind, levind, partlevsetbeg;
  int listEnd, PrevlistEnd;
  int terminalnum, addnum, removenum, runner, origiListEnd, atom;
  double arvo, higlev;
  int lokalefek;
  int begi, atto, node, curterminalnum;
 
  int tobehandled[maxkernelnode+1], leafloc[maxkernelnode+1];

  /*volume and center calculation*/
  double curvolu, vol, ala, yla;
  int componentnum, zeiger; 
  double curcente[maxdim+1], newcente[maxdim+1];                
  int atompointer, uppi, lowi;

  /*listchangeCCC*/
  int kokoout[2], begs[atomnumIntern+1];

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
    AtomlistAtom[i]=i;
    AtomlistNext[i]=i+1;
    i=i+1;
  }
  listEnd=*atomnum;
  AtomlistNext[listEnd]=0;

  /* Let us divide the lowest level set to disconnected parts */

  beg=1;
  terminalnum=*atomnum;

  begi=1;
  atto=AtomlistAtom[begi];
  while (begi>0){
     if (value[atto]>0){
        node=nodefinder[atto];
        tobehandled[node]=1;
     }    
     begi=AtomlistNext[begi];
     atto=AtomlistAtom[begi];
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

      begi=beghigher;
      for (ij=1; ij<=*nodenumOfDyaker; ij++) tobehandled[ij]=0;
      while (begi>0){
          atto=AtomlistAtom[begi];
          node=nodefinder[atto];
          tobehandled[node]=1;
          begi=AtomlistNext[begi];
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
   out[j]=atomsSepaNext[j];
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
     atto=AtomlistAtom[zeiger];
     vol=1;
     for (j=1; j<=*d; j++){
        lowi=inlow[(atto-1)*(*d)+j];  /*[atto][j]*/
        uppi=inupp[(atto-1)*(*d)+j];  /*[atto][j]*/
        vol=vol*(uppi-lowi)*step[j];
     }
     curvolu=curvolu+vol;
     zeiger=AtomlistNext[zeiger];
  }
  volume[i]=curvolu;
}

for (i=1; i<=componentnum; i++){  
  for (j=1; j<=*d; j++) curcente[j]=0;
  zeiger=component[i];
  while (zeiger>0){
     atompointer=AtomlistAtom[zeiger];
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
     zeiger=AtomlistNext[zeiger];
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
     AtomlistAtomOut[j]=AtomlistAtom[j];
     AtomlistNextOut[j]=AtomlistNext[j];
 }
 */
 

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
 int i, lkm, nodee, note, nodelkm;
 int Lempty, Rempty;
 int leftbeg, rightbeg, direction, akku, apu3;
 int direktiooni, poiju, prevpoiju, aatto, thisnoteempty;
 double splittiini;

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
  }

  for (i=1; i<=nodenumOfDyakerIntern; i++) separy[i]=0;


  nodelkm=nodenumOfDyaker;


  /* end of initializing */

  lkm=0;
  nodee=nodelkm;
  while (nodee>=1){          /* root is in position 1 */
      if ((leafloc[nodee]==1) && (tobehandled[nodee]==1)){ /* we are in leaf */

          tobehandled[parent[nodee]]=1;

          lkm=lkm+1;
          separy[nodee]=lkm;
          atomsSepaAtom[lkm]=infopointer[nodee];

          begsSepaBegs[lkm]=lkm;

          /*
          obs: we need not change
          begsSepaNext, atomsSepaNext, atomsLBounNext, atomsRBounNext
          since at the beginning set consist only one member:
          pointer is always 0, since we do not have followers
	  */

      }
      else if (tobehandled[nodee]==1){   /* not a leaf */

	  tobehandled[parent[nodee]]=1;

	  leftbeg=left[nodee];
	  rightbeg=right[nodee];

          if ((leftbeg==0) || (separy[leftbeg]==0)){
	    /*if left child does not exist*/

	    /*
            note that since we consider subsets of the
            terminal nodes of the original tree, it may happen
            that leftbeg>0 but left child does not exist
	    */
           
	    separy[nodee]=separy[rightbeg];

            /*
            we need that all lists contain as many members
            left boundary is empty, but we will make it a list
            of empty lists
            */

          }
          else{ /* eka else */
            if ((rightbeg==0) || (separy[rightbeg]==0)){
	   	/* right child does not exist */
    
		separy[nodee]=separy[leftbeg];
                      
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
		note=separy[rightbeg];
                    while (note>0){
                        thisnoteempty=1;

                        poiju=begsSepaBegs[note];
                        prevpoiju=poiju;
                        while (poiju>0){
		           aatto=atomsSepaAtom[poiju];
                           if (splittiini>=low[aatto][direktiooni]){
                               /* this atom belongs to boundary */
                               if (thisnoteempty==1){
                                   /* poiju is the 1st non-empty */
                                   begsLeftBoun[note]=poiju;
                               }                
                               Lempty=0;
                               atomsLBounNext[prevpoiju]=poiju;    
                               atomsLBounNext[poiju]=0;
                               prevpoiju=poiju;

                               thisnoteempty=0;
                           }
                           poiju=atomsSepaNext[poiju];
                        }
                        if (thisnoteempty) begsLeftBoun[note]=0;
                        note=begsSepaNext[note];
                     }

                    /* right boundary of left child */

		    Rempty=1;
		    note=separy[leftbeg];
                     while (note>0){
			 thisnoteempty=1;
                       
			 poiju=begsSepaBegs[note];
			 prevpoiju=poiju;
                         while (poiju>0){
			    aatto=atomsSepaAtom[poiju];
                            if (splittiini<=upp[aatto][direktiooni]){
                               /* this atom belongs to boundary */
                               if (thisnoteempty==1){
                                   /* poiju is the 1st non-empty */
                                   begsRighBoun[note]=poiju;
                               }                
                               Rempty=0;
                               atomsRBounNext[prevpoiju]=poiju;    
                               atomsRBounNext[poiju]=0;
                               prevpoiju=poiju;
                           
                               thisnoteempty=0;
                           }
                           poiju=atomsSepaNext[poiju];
                        }
                        if (thisnoteempty) begsRighBoun[note]=0;
                        note=begsSepaNext[note];
                     }
    
	             /* check whether one of boundaries is empty */


                     if (Lempty || Rempty){
			 /* one of boundaries is empty */
                        
			 /* concatenating separate parts */
                       
			 akku=separy[leftbeg];
                 
			 begsSepaNext[akku]=separy[rightbeg]; 
                         /* concatenate list of separate sets */

			 separy[nodee]=separy[leftbeg];
			 akku=separy[rightbeg];
                          
                        /* left boundaries of sets in right child are empty */
                        /* end of concatenating */



                     }
           	     else{ /*both children exist, both boundaries non-empty*/
			 direction=vec[nodee];

          
                         apu3=joincongenC(leftbeg,
                                          rightbeg,
                                          direction,
                                          d);
			 separy[nodee]=apu3;

		     /*
                     jc=joincongen(leftbeg,rightbeg,separy,
                     begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
                     atomsSepaNext,atomsSepaAtom,atomsLBounNext,
                     atomsRBounNext,
                     direction,low,upp)   
		     */

                                }
            } /* toka else */
        } /* eka else */
    }  /* else not a leaf */
      nodee=nodee-1;
}

  return separy[1];

/*
totbegSepary=separy[1],
begsSepaNext=begsSepaNext,
begsSepaBegs=begsSepaBegs,
atomsSepaNext=atomsSepaNext,
atomsSepaAtom=atomsSepaAtom
*/

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

 int pino[maxkernelnode+1], apu, i, pinin, cur;

 /* itemnum=length(left) */

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
}


/************************************************************/


/*
see: joincongen.R

INPUT: 
leftbeg,rightbeg,direction,d
separy,
begsSepaNext,begsSepaBegs,begsLeftBoun,begsRighBoun,
atomsSepaNext,atomsSepaAtom,atomsLBounNext,atomsRBounNext,
low,upp
OUTPUT:
totbegSepary, 
begsSepaNext, begsSepaBegs, begsLeftBoun, begsRighBoun,
atomsSepaNext,                atomsLBounNext, atomsRBounNext
*/


int joincongenC(int leftbeg,
                int rightbeg, 
                int direction, 
                int d)     
{

  int apu, i;

  /*startpointsC*/
  int startpointsS[induksiInt+1];
  int startpointsB[induksiInt+1];
  int startpointsNewBleft[induksiInt+1];
  int startpointsNewBright[induksiInt+1];
  int m, mleft, mright;
  int mleftout[2], mrightout[2];

  /*makelinksC*/
  /*int linkit[componumInt+1][componumInt+1];*/
  int outlinkit[1+componumInt*componumInt];

  /*declevnewCCC*/
  int sepnum, tulospit[2], outtulos[componumInt+1];
  /*int res[componumInt+1][componumInt+1];*/
  /*int inlinkit[1+componumInt*componumInt];*/

  /*joinsetsC*/
  int Totalbeg, Totalbegout[2];


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
  int pino[mmm+1], pinind; 
    /*pinoon laitetaan aina jos koskettaa, max kosketuksia m*/
  int i, j, k, curleima, merkatut[mmm+1], curpallo, touch;  
  int linksit[mmm+1][mmm+1], tulos[mmm+1][mmm+1], apu;

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
  ml=induksi;
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
  int indelow1[maxdim+1], indelow2[maxdim+1];
  int indeupp1[maxdim+1], indeupp2[maxdim+1];
  int begbeg1, begbeg2, atom1, atom2, touch;

  int m=mleft+mright;
  int linkit[mmm+1][mmm+1];

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
		  atom1=atomsSepaAtom[begbeg1];
		  atom2=atomsSepaAtom[begbeg2];
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
                  begbeg2=atomsLBounNext[begbeg2];
              }
              begbeg1=atomsRBounNext[begbeg1];
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
		 atom1=atomsSepaAtom[begbeg1];
		 atom2=atomsSepaAtom[begbeg2];
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
                 begbeg2=atomsRBounNext[begbeg2];
             }
	     begbeg1=atomsLBounNext[begbeg1];
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

}



/*********************************************************/

/*
OutPut:

begsSepaNext, begsSepaBegs, begsLeftBoun, begsRighBoun,
atomsSepaNext, atomsLBounNext, atomsRBounNext
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
  int res[sepnummm][mmm];
  int Totalbeg, tavoite, hiihtaja, nykyinen;
  int laskuri, curre, kL, kR;
 
 /* m=induksiInt */
  int osoittajaS[mmm+1];  /* make vector of pointers to the begs of sets */
  int osoittajaNewBleft[mmm+1];
  int osoittajaNewBright[mmm+1];

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

 Totalbeg=separy[leftbeg];

 tavoite=1;
 hiihtaja=Totalbeg;
 while ((begsSepaNext[hiihtaja]>0) && (tavoite<sepnum)){
     hiihtaja=begsSepaNext[hiihtaja];
     tavoite=tavoite+1;
 }  
 if (tavoite<sepnum){  /* now hiihtaja points to the end of the first list */
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
  
     /* handle separy */
  
     begsSepaBegs[nykyinen]=osoittajaS[1];   /* always non-zero */
     
     k=1;
     while (k<=(len-1)){    
	 curre=osoittajaS[k];
         while (atomsSepaNext[curre]>0){ /*find the end*/
	      curre=atomsSepaNext[curre];
         }
         atomsSepaNext[curre]=osoittajaS[k+1];
	 k=k+1;
     }

/*
if ((atomsSepaNext[12]==11) && (globapu2==1)) globapu3=globapu3+1;
if ((leftbeg==2) && (rightbeg==3) && (globapu2==1)){
       apsu1=atomsSepaNext[11];
       apsu2=atomsSepaNext[12];
       apsu3=atomsSepaNext[13];
       apsu4=atomsSepaNext[14];
       apsu5=atomsSepaNext[15];
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
	 begsLeftBoun[nykyinen]=0;
     }
     else{        /*  kL is first non zero  */
	 kL=k;
	 begsLeftBoun[nykyinen]=osoittajaNewBleft[kL];
	 /*
         update the list of left boundaries
         concatenate the lists of atoms
	 */
	 k=kL;
         while (k<=(len-1)){    
	     curre=osoittajaNewBleft[k]; /* curre is not zero */
             while (atomsLBounNext[curre]>0){    /* find the end */
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
           else{  /*  found nonzero */
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
	     curre=osoittajaNewBright[k];        /* curre is not zero */
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
     ii=ii+1;
 }

 Totalbegi[1]=Totalbeg;

 apu=1;
 return apu;

}


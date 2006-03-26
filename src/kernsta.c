/* 
R CMD SHLIB -o kerleCversio kernsta.c 
dyn.load("kerleCversio") 
dentree<-.C("decomdyaC",
               as.integer(numofall),
               as.integer(atomnum),
               as.double(inlevseq),
               as.integer(inN),
               as.integer(d),         
               as.integer(levnum),   
               as.double(volofatom),
               as.double(inminim),
               as.double(h),
               as.double(indelta),
               as.integer(nodenumOfDyaker),
               as.integer(inleft),
               as.integer(inright),
               as.integer(inparent), 
               as.double(invalue),
               as.integer(kg$index),
               as.integer(innodefinder),
               level = double(numofall+1),
               parent = integer(numofall+1),
               component = integer(numofall+1),
               volume = double(numofall+1),
               center = double(d*numofall+1),
               efek = integer(1),
               AtomlistAtomOut = integer(numofall+1),
               AtomlistNextOut = integer(numofall+1),
               numOfAtoms = integer(1))

               apu = double(numofall+1))

*/

#define numofallIntern 1000000
#define atomnumIntern 1000000
#define nodenumOfDyakerIntern 1000000
#define maxdim 4

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

#include <stdlib.h> 
#include <math.h>

int listchangeCC(int totbegSepary,
		 int beg);
int declevdyaCC(int beg,
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
int joingeneCC(int node,
               int leftbeg,
               int rightbeg,
	       int direction,
	      /* int *indeks[], */
               int d);
int joinconneCC(int leftbeg,
               int rightbeg, 
               int direction, 
               /* int *indeks[], */
		int d);
int startpointsCC(int leftbeg, 
		  int rightbeg);
int makelinksCC(int direction,
              /* int *indeks[], */
		int d);
int declevnewCC();
int joinSetsCC(int leftbeg,
             int rightbeg,
	       int sepnum);

void decomdyaC(int *numofall, 
               int *atomnum,
               double *levseq,
               int *N,
               int *d,       /* d<-length(N) */
	       int *levnum,  /* levnum<-length(levseq) */
               double *volofatom,
               double *minim,
               double *h,
               double *delta,
               /* tree structures */
               int *nodenumOfDyaker,
               int *left,
               int *right,
               int *parentKertree,
               double *value,
               int *inindex,
               int nodefinder[],
               /* output */
	       double *level,      /* of length numofall+1 */
               int *parent, 
               int *component, 
	       double *volume,
               double *center,
               int *efek,   /* gives effective length of level,... */
               int *AtomlistAtomOut,
               int *AtomlistNextOut,
               int *numOfAtoms)



{
  int totbegSepary;
  /*
  int pinoComponent[*numofall+1];  pointer to component, level,... 
  int pinoTaso[*numofall+1];       ordinal of level (pointer to levseq) 
  */
  int *pinoComponent = (int *)malloc(sizeof(int) * (*numofall+1));
  int *pinoTaso = (int *)malloc(sizeof(int) * (*numofall+1));

  int componum, beghigher;
  int i, j, beg, koko, pinind, ind, levind, partlevsetbeg;
  int listEnd, PrevlistEnd;
  int terminalnum, addnum, removenum, runner, origiListEnd, atom;
  double arvo, higlev;
  int lokalefek;
  /*  int indeks[(*atomnum)+1][(*d)+1]; */
 
  int componentnum, numofatoms, zeiger; /* for cvolumdya */

  int atompointer;

  /*double curcente[*d+1]; */                   /* for ccentedya */
  double *curcente = (double *)malloc(sizeof(double) * (*d+1));
 

  if (pinoComponent == NULL) exit(1); 
  if (pinoTaso == NULL) exit(1); 
  if (curcente == NULL) exit(1); 

  /* Initialize indeks */

  for (j=0; j<=*d-1; j++)
     for (i=1; i<=*atomnum; i++)
	indeks[i][j+1]=inindex[j*(*atomnum)+i-1]; 
  
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

  for (i=1; i<=*numofall; i++){
      pinoComponent[i]=0;
      pinoTaso[i]=0; 
  }

  /*end of initializing the local vectors */

  componum=*numofall;

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
  totbegSepary=declevdyaCC(beg,
                          left,
                          right,
                          parentKertree,
                          /* indeks, */
                          nodefinder,
	                  N,
                          *nodenumOfDyaker,
                          terminalnum,
                          *d); 

  /* modify begs, at most "terminalnum" changes */  
  koko=listchangeCC(totbegSepary,
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
      totbegSepary=declevdyaCC(beg,
                            left,
                            right,
                            parentKertree,
                            /* indeks, */
                            nodefinder,
	                    N,
                            *nodenumOfDyaker,
                            terminalnum,
                            *d); 

      /* modify begs, at most "terminalnum" changes */  
      koko=listchangeCC(totbegSepary,
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
   volume[i]=numofatoms*(*volofatom);
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
     center[(i-1)*(*d)+j]=(*volofatom)*curcente[j]/volume[i];
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
 *numOfAtoms=listEnd;
 
 for (j=1; j<=listEnd; j++){
     AtomlistAtomOut[j]=AtomlistAtom[j];
     AtomlistNextOut[j]=AtomlistNext[j];
 }
 
  free(pinoComponent); 
  free(pinoTaso); 


}




/* END OF MAIN  */




/* changes AtomlistNext, AtomlistAtom */
/* needs begsSepaNext,begsSepaBegs,atomsSepaNext,atomsSepaAtom */

/* create begs: beginnings of lists of atoms */
/* beg is indeks to AtomlistAtom/Next */
/* totbegsepary is indeks to begsSepaBegs/Next */

int listchangeCC(int totbegSepary,
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

int declevdyaCC(int beg,
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
    /*
 int nextFloor[terminalnum+1];
 int currFloor[terminalnum+1];
 int already[nodenumOfDyaker+1];
    */
 int *nextFloor = (int *)malloc(sizeof(int) * (terminalnum+1));
 int *currFloor = (int *)malloc(sizeof(int) * (terminalnum+1));
 int *already = (int *)malloc(sizeof(int) * (nodenumOfDyaker+1));
  
 int i, k, r, lkm, nexlkm, curlkm, curre, atom, node, note;
 int exists, Lempty, Rempty;
 int leftbeg, rightbeg, direction, akku, totbegSepary, apu;
 int j;

 if (nextFloor == NULL) exit(1);
 if (currFloor == NULL) exit(1);
 if (already == NULL) exit(1);
 
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
       apu=joingeneCC(node,
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
     apu=joingeneCC(node,
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
                       separy[node]=joinconneCC(leftbeg,
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

 free(nextFloor);
 free(currFloor);
 free(already);

}













int joingeneCC(int node,
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
        separy[node]=joinconneCC(leftbeg, 
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

int joinconneCC(int leftbeg,
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

 apu=startpointsCC(leftbeg,
                 rightbeg);

  /* 2. We make "links" matrix and apply declev */

  apu=makelinksCC(direction,
		/* mleft, m, indeks, */
                d);
  sepnum=declevnewCC();  /* tulos will be filled */

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

 totbegsepary=joinSetsCC(leftbeg,
                       rightbeg,
                       sepnum);
 return totbegsepary;
}











int startpointsCC(int leftbeg, 
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
int makelinksCC(int direction,
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
int declevnewCC()   /* int m */
{
    /*int pino[m]; */
  /* pino<-matrix(0,m,1) 
     pinoon laitetaan aina jos koskettaa, max kosketuksia m */
 int *pino = (int *)malloc(sizeof(int) * (m+1));
 int pinind;
  /* pinossa viitataan rindeksin elementteihin */
  int i, j, k, curleima; 
  /* i ja j viittavat rindeksit-vektoriin, jonka alkiot viittavat atomeihin */
  /*int merkatut[m]; */
  int *merkatut = (int *)malloc(sizeof(int) * (m+1));
  int curpallo, tulospit;
  int touch;  
  /*  resultat.tulos=malloc(m*sizeof(int));  */

  if (pino == NULL) exit(1);
  if (merkatut == NULL) exit(1);

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
 
  free(pino);
  free(merkatut);


}












int joinSetsCC(int leftbeg,
             int rightbeg,
             int sepnum)
{
  int tavoite, hiihtaja, nykyinen, i, j, k;
  /*
  int osoittajaS[m];  make vector of pointers to the begs of sets 
  int osoittajaNewBleft[m];
  int osoittajaNewBright[m];
  */
  int *osoittajaS= (int *)malloc(sizeof(int) * (m+1));
  int *osoittajaNewBleft= (int *)malloc(sizeof(int) * (m+1));
  int *osoittajaNewBright= (int *)malloc(sizeof(int) * (m+1));

  int sol, len, laskuri, TotalBeg, curre, kL, kR;

  if (osoittajaS == NULL) exit(1);
  if (osoittajaNewBleft == NULL) exit(1);
  if (osoittajaNewBright == NULL) exit(1);

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

 free(osoittajaS);
 free(osoittajaNewBleft);
 free(osoittajaNewBright);

}






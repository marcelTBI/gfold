#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "head.h"




paramT *P;


long *indx4[4];
int *c,*fML,*fML1,*fPL1, *fPL,*indx2;
int *f5,*f51, *I51, *I5, *IML1, *IML, *IPL1, *IPL, *Ibc; 

int *Gtight, *Gu, *Guep, *Gv, *Gw, *Gumm; 
int *Gump, *Gvm, *Gwm, *Gup, *Gvp, *Gwp;
int *Gh, *Gvep, *Gvmp;




void initial (int n)
{
	int i,k;
	for (i=0;i<4;i++) {
	  indx4[i]=(long *) space (sizeof(long)*(n+5));
	}
	for (i=0;i<=n+4;i++)
	  indx4[0][i]=i;
	for (k=1;k<4;k++) {
	  indx4[k][0]=0; indx4[k][1]=1;
	  for (i=2;i<=n+4;i++) {
	    indx4[k][i]=indx4[k][i-1]+indx4[k-1][i];
	  }
	}

	indx2 = (int *) space(sizeof(int)*(n+5));
	for (i = 1; i <= n; i++)
	  indx2[i] = (i*(i-1)) >> 1;        /* n(n-1)/2 */

	c     = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	f5    = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	f51   = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	fML   = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	fML1  = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	fPL   = (int *) space(sizeof(int)*((n*(n+1))/2+2));
	fPL1  = (int *) space(sizeof(int)*((n*(n+1))/2+2));

	Ibc= (int *) space(sizeof(int)*((n*(n+1))/2+2));
	I51= (int *) space(sizeof(int)*((n*(n+1))/2+2));
	I5= (int *) space(sizeof(int)*((n*(n+1))/2+2));
	IML1= (int *) space(sizeof(int)*((n*(n+1))/2+2));
	IML= (int *) space(sizeof(int)*((n*(n+1))/2+2));
	IPL1= (int *) space(sizeof(int)*((n*(n+1))/2+2));
	IPL= (int *) space(sizeof(int)*((n*(n+1))/2+2));

	Gtight=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gh=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gu=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Guep=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gv=(int *) space (sizeof(int)*(indx4[3][n+1]));	
	Gw=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gumm=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gump=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gvm=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gwm=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gup=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gvp=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gwp=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gvep=(int *) space (sizeof(int)*(indx4[3][n+1]));
	Gvmp=(int *) space (sizeof(int)*(indx4[3][n+1]));
	
	ptable=(int *) space (sizeof(int)*(n+2));
// 	printf("Length: %d\n", n);
// 	printf("Require %fGb\n", (double)((14*(n*(n+1))/2+2)*(sizeof(double))
// 		+15*indx4[3][n+1]*(sizeof(double)))/1024/1024/1024);
// 	printf("INDEX: %ld\n", indx4[3][n+1]);
}
void freevar()
{
	int i;
	free(Gtight);
	free(Gh);
	free(Gu); free(Guep); free(Gv); free(Gw);
	free(Gumm); free(Gump), free(Gvm); free(Gwm);
	free(Gup); free(Gvp); free(Gwp);
	free(Gvep); free(Gvmp);

	free(c); free(fML); free(fML1); free(fPL1); free(fPL);
	free (Ibc); free(f5); free(f51); 
	free(I51); free(I5);
	free(IML1); free(IML);
	free(IPL1); free(IPL);
	
	for (i=0;i<4;i++)
	  free(indx4[i]);
	free(indx2);
}



void fillarray(const char *seq)
{
	int i,j,r,s,k1,k2,p,q,gap;
	int length,type1,type2,typek,G,type;
	int left,right;
	int score;
	long index;

	length=strlen(seq);
	for (i=1;i<=length;i++)
	  for (j=i;j<i+TURN;j++) {
	    c[indx2[j]+i]=MAXENG;
	    f5[indx2[j]+i]=0;
	    f51[indx2[j]+i]=MAXENG;
	    fML[indx2[j]+i]=MAXENG;
	    fML1[indx2[j]+i]=MAXENG;
	    fPL[indx2[j]+i]=MAXENG;
	    fPL1[indx2[j]+i]=MAXENG;
	}

	for (i=length-TURN;i>=1;i--)
	  for (j=i+TURN;j<=length;j++) {
/*
	for (gap=TURN;gap<length;gap++)
        #pragma omp parallel default (shared) private \
 (i,j,r,s,k1,k2,p,q,type1, type2,typek,G,type,left,right,score,index)
        #pragma omp for schedule (static)
	  for (i=1;i<length;i++) {
	    j=i+gap;
*/
	    type = BP_pair[SEQ[i]][SEQ[j]];
	    if (type) {
	      c[indx2[j]+i] = HairpinE(j-i-1, type, SEQ[i+1], SEQ[j-1], seq+i-1);

	      for (p = i+1; p <= j-2-TURN; p++) {
	        for (q = p+TURN; q < j; q++) {
	          type2 = BP_pair[SEQ[p]][SEQ[q]];

	          if (type2==0) continue;
	          type2 = rtype[type2];

	          G = LoopEnergy(p-i-1, j-q-1, type, type2,
				SEQ[i+1], SEQ[j-1], SEQ[p-1], SEQ[q+1]);
	          if (G+c[indx2[q]+p] < c[indx2[j]+i]) c[indx2[j]+i]=G+c[indx2[q]+p];
	        }
	      }
	      for (k1=i+1;k1<j-1;k1++) {
	        G=PARS->MLintern[type]+PARS->MLclosing;
	        if (fML1[indx2[k1]+i+1]+fML[indx2[j-1]+k1+1]+G < c[indx2[j]+i])
	          c[indx2[j]+i]=fML1[indx2[k1]+i+1]+fML[indx2[j-1]+k1+1]+G;
	      }
	    } else {
	      c[indx2[j]+i]=MAXENG;
	    }

	    f51[indx2[j]+i]=MAXENG;
	    for (k1=i;k1<j;k1++) {
	      type = BP_pair[SEQ[k1]][SEQ[j]];
	      if (type>2) G=PARS->TerminalAU; else G=0;
	      if (type && c[indx2[j]+k1]+G < f51[indx2[j]+i]) f51[indx2[j]+i]=c[indx2[j]+k1]+G;

	    }

	    f5[indx2[j]+i]=0;
	    if (f51[indx2[j]+i] < f5[indx2[j]+i])
	      f5[indx2[j]+i]=f51[indx2[j]+i];
	    for (k1=i;k1<j;k1++) {
	      if (f51[indx2[k1]+i]+f5[indx2[j]+k1+1]<f5[indx2[j]+i])
	        f5[indx2[j]+i]=f51[indx2[k1]+i]+f5[indx2[j]+k1+1];
	    }

	    fML1[indx2[j]+i]=MAXENG;
	    for (k1=i;k1<j;k1++) {
	      type=BP_pair[SEQ[k1]][SEQ[j]];
	      if (type) {
	        G=PARS->MLintern[type]+(k1-i)*PARS->MLbase;
	        if (c[indx2[j]+k1]+G < fML1[indx2[j]+i])
	          fML1[indx2[j]+i] = c[indx2[j]+k1]+G;
	      }
	    }

	    fML[indx2[j]+i]=MAXENG;
	    if (fML1[indx2[j]+i] < fML[indx2[j]+i])
	      fML[indx2[j]+i]=fML1[indx2[j]+i];
	    for (k1=i;k1<j;k1++) {
	      if (fML1[indx2[k1]+i]+fML[indx2[j]+k1+1] < fML[indx2[j]+i])
	        fML[indx2[j]+i]=fML1[indx2[k1]+i]+fML[indx2[j]+k1+1];
	      if (fML1[indx2[k1]+i]+(j-k1)*PARS->MLbase < fML[indx2[j]+i])
	        fML[indx2[j]+i]=fML1[indx2[k1]+i]+(j-k1)*PARS->MLbase;
	    }

	    fPL1[indx2[j]+i]=MAXENG;
	    for (k1=i;k1<j;k1++) {
	      type=BP_pair[SEQ[k1]][SEQ[j]];
	      if (type) {
	        G=BETA2+(k1-i)*BETA3;
	        if (c[indx2[j]+k1]+G < fPL1[indx2[j]+i])
	          fPL1[indx2[j]+i] = c[indx2[j]+k1]+G;
	      }
	    }

	    fPL[indx2[j]+i]=MAXENG;
	    if (fPL1[indx2[j]+i] < fPL[indx2[j]+i])
	      fPL[indx2[j]+i]=fPL1[indx2[j]+i];
	    for (k1=i;k1<j;k1++) {
	      if (fPL1[indx2[k1]+i]+fPL[indx2[j]+k1+1] < fPL[indx2[j]+i])
	        fPL[indx2[j]+i]=fPL1[indx2[k1]+i]+fPL[indx2[j]+k1+1];
	      if (fPL1[indx2[k1]+i]+(j-k1)*BETA3 < fPL[indx2[j]+i])
	        fPL[indx2[j]+i]=fPL1[indx2[k1]+i]+(j-k1)*BETA3;
	    }
// secondary structure S done!

	    type1=BP_pair[SEQ[i]][SEQ[j]];
	    for (r=i+1;r<j-1;r++)
	      for (s=j-1;s>r;s--) {
// 	      for (s=r+1;s<=j-1;s++) {
	        type2=BP_pair[SEQ[s]][SEQ[r]];
	        index=indx4[0][i]+indx4[1][r]+indx4[2][s]+indx4[3][j];
	        score=MAXENG;
	        if (type1 && type2 && s-r>TURN) { // start Gtight
	          score=(int)(SIGMA*LoopEnergy(r-i-1,j-s-1,type1,type2,
		  SEQ[i+1],SEQ[j-1],SEQ[r-1],SEQ[s+1]));
// G is empty, it is big interior loop
	          for (k1=i+1;k1<r;k1++)
	            for (k2=s+1;k2<j;k2++) {
	              typek=BP_pair[SEQ[k2]][SEQ[k1]];
	              if (typek) {
	                G=(int)(SIGMA*LoopEnergy(k1-i-1,j-k2-1,type1,typek,
		SEQ[i+1],SEQ[j-1],SEQ[k1-1],SEQ[k2+1]));
	                if (G+Gtight[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2]] < score)
	                  score=G+Gtight[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2]];
// G is an interior loop plus another G

	                G=PARS->MLclosing+PARS->MLintern[type1];
	                if (k1-i>1) {
	                  if (G+Gtight[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2]]
		+(j-k2-1)*PARS->MLbase+IML[indx2[k1-1]+i+1] < score)
	                  score=G+Gtight[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2]]
		+(j-k2-1)*PARS->MLbase+IML[indx2[k1-1]+i+1];
	                }
// G is a miltiloop, and the rhs. is empty
	                if (j-k2>1) {
	                  if (G+Gtight[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2]]
		+(k1-i-1)*PARS->MLbase+IML[indx2[j-1]+k2+1] < score)
		          score=G+Gtight[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2]]
		+(k1-i-1)*PARS->MLbase+IML[indx2[j-1]+k2+1];
	                }
// G is a miltiloop, and the lhs. is empty
	                if (k1-i>1 && j-k2>1) {
	                  if (G+Gtight[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2]]
		+IML[indx2[k1-1]+i+1]+IML[indx2[j-1]+k2+1] < score)
	                  score=G+Gtight[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2]]
		+IML[indx2[k1-1]+i+1]+IML[indx2[j-1]+k2+1];
	                }
// G is a miltiloop, and both the lhs. and rhs are not empty.3
	              }
	          }
	        }
	        Gtight[index]=score;
//calculate the G matrix
// BEAT2 is the penalty for the inner base pair contributes to pseudoknot.

	        score=MAXENG;
	        type1=BP_pair[SEQ[i]][SEQ[j]];
	        if (type1) {
	          for (k2=s>r+TURN?s:r+TURN; k2<j; k2++) {
	            type2=BP_pair[SEQ[r]][SEQ[k2]];
	            if (type2==0) continue;
	            if (k2-s>0) {
	              right=IPL[indx2[k2-1]+s];
	              if ((k2-s)*BETA3 < right) right=(k2-s)*BETA3;
	            } else 
	              right=0;
	            if (right+Gtight[indx4[0][i]+indx4[1][r]+indx4[2][k2]+indx4[3][j]] < score)
	              score=right+Gtight[indx4[0][i]+indx4[1][r]+indx4[2][k2]+indx4[3][j]];
	          }
	        }
	        Gh[index]=score;
// calculate the Gh matrix (an intermediate for GX, X=u,v,w)

	        score=MAXENG;
	        for (k1=i;k1<r;k1++) {
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          if (BP_pair[SEQ[k1]][SEQ[j]]>2) G=PARS->TerminalAU; else G=0;
	          if (k1-i>0) left=I5[indx2[k1-1]+i]; else left=0;
	          if (left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G < score)
	            score=left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G;
	        }
	        Gu[index]=score;
//calculate the Gu matrix

	        score=MAXENG;
	        for (k1=i;k1<r;k1++) {
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          if (BP_pair[SEQ[k1]][SEQ[j]]>2) G=PARS->TerminalAU; else G=0;
	          if (k1-i>0) {
	            left=IPL[indx2[k1-1]+i];
	            if ((k1-i)*BETA3 < left) left=(k1-i)*BETA3;
	          }else 
	            left=0;
	          if (left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G < score)
	            score=left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G;
	        }
	        Guep[index]=score;
//calculate the Guep matrix

	        score=MAXENG;
	        for (k1=i;k1<r;k1++) { 
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          if (k1-i>0) {
	            left=IML[indx2[k1-1]+i];
	            if ((k1-i)*PARS->MLbase < left) left=(k1-i)*PARS->MLbase;
	          } else 
	            left=0;
	          typek=BP_pair[SEQ[k1]][SEQ[j]];
	          G=PARS->MLintern[typek];
	          if (left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G < score)
	            score=left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G;
	        }
	        Gumm[index]=score;
//calculate the Gumm matrix

	        score=MAXENG;
	        for (k1=i;k1<r;k1++) {
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          if (k1-i>0) {
	            left=IPL[indx2[k1-1]+i];
	            if ((k1-i)*BETA3 < left) left=(k1-i)*BETA3;
	          } else 
	            left=0;
	          typek=BP_pair[SEQ[k1]][SEQ[j]];
	          G=PARS->MLintern[typek];
	          if (left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G < score)
	            score=left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G;
	        }
	        Gump[index]=score;
//calculate the Gump matrix

	        score=MAXENG;
	        for (k1=i;k1<r;k1++) {
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          if (k1-i>0) {
	            left=IPL[indx2[k1-1]+i];
	            if ((k1-i)*BETA3 < left) left=(k1-i)*BETA3;
	          } else 
	            left=0;
	          G=BETA2; //The penalty for the maximum arc in pseudoknot 
	          if (left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G < score)
	            score=left+Gh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]+G;
	        }
	        Gup[index]=score;
//calculate the Gup matrix

	        score=MAXENG;
	        for (k1=i+1;k1<r-1;k1++) 
	          for (k2=s+1;k2<j-1;k2++) {
	            left=Gu[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=Guep[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            if (left+right+2*BETA2 < score) score=left+right+2*BETA2;
	        }
	        Gv[index]=score;
// calculate the Gv matrix

	        score=MAXENG;
	        for (k1=i+1;k1<r-1;k1++) 
	          for (k2=s+1;k2<j-1;k2++) {
	            left=Guep[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=Guep[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            if (left+right+2*BETA2 < score) score=left+right+2*BETA2;
	        }
	        Gvep[index]=score;
// calculate the Gv matrix

	        score=MAXENG;
	        for (k1=i+1;k1<r-1;k1++)
	          for (k2=s+1;k2<j-1;k2++) {
	            left=Gumm[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=Gump[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            if (left+right+2*BETA2 < score) score=left+right+2*BETA2;
	        }
	        Gvm[index]=score;
// calculate the Gvm matrix

	        score=MAXENG;
	        for (k1=i+1;k1<r-1;k1++)
	          for (k2=s+1;k2<j-1;k2++) {
	            left=Gump[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=Gump[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            if (left+right+2*BETA2 < score) score=left+right+2*BETA2;
	        }
	        Gvmp[index]=score;
// calculate the Gvm matrix

	        score=MAXENG;
	        for (k1=i+1;k1<r-1;k1++)
	          for (k2=s+1;k2<j-1;k2++) {
	            left=Gup[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=Gup[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            if (left+right+2*BETA2 < score) score=left+right+2*BETA2;
	        }
	        Gvp[index]=score;
// calculate the Gvp matrix

	        score=MAXENG;
	        for (k1=s+1;k1<j;k1++)
	          for (k2=k1+3;k2<j;k2++) {
	            left=Guep[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=Guep[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETAB+2*BETA2; // the penalty for the crossing matrix (shadow X)
	            if (left+right+G < score) score=left+right+G;
// U1 U'1 U2 U'2

	            left=Gvep[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=Guep[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETAD+BETA2; // the penalty for the crossing matrix (shadow X)
	            if (left+right+G < score) score=left+right+G;
// V1 U1 V2 U'2
	        }
	        Gw[index]=score;
// calculate the Gw matrix 

	        score=MAXENG;
	        for (k1=s+1;k1<j;k1++)
	          for (k2=k1+3;k2<j;k2++) {
	            left=Gump[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=Gump[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETABm+2*BETA2; // the penalty for the crossing matrix (shadow X)
	            if (left+right+G < score) score=left+right+G;
// U1 U'1 U2 U'2

	            left=Gvmp[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=Gump[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETADm+BETA2; // the penalty for the crossing matrix (shadow X)
	            if (left+right+G < score) score=left+right+G;
// V1 U1 V2 U'2
	        }
	        Gwm[index]=score;
// calculate the Gwm matrix

	        score=MAXENG;
	        for (k1=s+1;k1<j;k1++)
	          for (k2=k1+3;k2<j;k2++) {
	            left=Gup[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=Gup[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETABp+2*BETA2; // the penalty for the crossing matrix (shadow X)
	            if (left+right+G < score) score=left+right+G;
// U1 U'1 U2 U'2

	            left=Gvp[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=Gup[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETAD+BETA2; // the penalty for the crossing matrix (shadow X)
	            if (left+right+G < score) score=left+right+G;
// V1 U1 V2 U'2
	        }
	        Gwp[index]=score;
// calculate the Gwp matrix

	    } //4-dimension for-loop
// 4-dimensin matrix done!

	    score=MAXENG;
	    type = BP_pair[SEQ[i]][SEQ[j]];
	    if (type) {
	      for (k1=i+1;k1<j;k1++)
	        for (k2=k1+TURN;k2<j;k2++) {
	          type2=BP_pair[SEQ[k2]][SEQ[k1]];
	          if (type2) {
	            G=LoopEnergy(k1-i-1,j-k2-1,type,type2,
		SEQ[i+1],SEQ[j-1],SEQ[k1-1],SEQ[k2+1]);
	            score=G+Ibc[indx2[k2]+k1];

	            G=PARS->MLclosing+PARS->MLintern[type]+PARS->MLintern[type2];
	            if (k1-i>1) {
	              left=IML[indx2[k1-1]+i+1];
	              right=(j-k2-1)*PARS->MLbase;
	              if (left+right+G+Ibc[indx2[k2]+k1] < score)
	                score=left+right+G+Ibc[indx2[k2]+k1];
	            }

	            if (j-k2>1) {
	              left=(k1-i-1)*PARS->MLbase;
	              right=fML[indx2[j-1]+k2+1]; 
	              if (left+right+G+Ibc[indx2[k2]+k1] < score)
	                score=left+right+G+Ibc[indx2[k2]+k1];
	            }

	            if (k1-i>1 && j-k2>1) {
	              left=IML[indx2[k1-1]+i+1]; 
	              right=fML[indx2[j-1]+k2+1]; 
	              if (left+right+G+Ibc[indx2[k2]+k1] < score)
	                score=left+right+G+Ibc[indx2[k2]+k1];
	            }
	          }
	      }
	      G=PARS->MLclosing+PARS->MLintern[type];
	      for (k1=i+TURN;k1<=j-1;k1++) {
	        left=IML1[indx2[k1]+i+1];
	        if (k1<j-1) {
	          right=fML[indx2[j-1]+k1+1];
	          if ((j-1-k1)*PARS->MLbase < right) right = (j-1-k1)*PARS->MLbase;
	        } else 
	          right=0;
	        if (G+left+right < score) score=G+left+right;
	      }
	      Ibc[indx2[j]+i]=score;
	    } else 
	      Ibc[indx2[j]+i]=0;
// calculate Ibc, which is a 1-structure where (i,j) is a base pair 

	    score=MAXENG;
	    for (r=i+1;r<j;r++) {
	      if (Gv[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]+BETAA < score) 
	        score=Gv[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]+BETAA;
	    }
//I'->V1 V2 (shadow A)
	    for (r=i+1;r<j;r++)
	      for (k1=r+2;k1<j;k1++)
	        for (s=k1+2;s+1<j;s++) {
	          left=Gv[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=Guep[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          G=BETAC+BETA2;
	          if (left+right+G < score) score=left+right+G;
//I'->V1 U1 V2 U2 (shadow C)

	          left=Gu[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=Gw[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          if (left+right+BETA2 < score) score=left+right+BETA2;
//I'->W1 U1 W2 U2 (shadow B & D)
	    }
	    I51[indx2[j]+i]=score;
// calculate I51

	    score=f5[indx2[j]+i];
	    if (I51[indx2[j]+i]<score) score = I51[indx2[j]+i];
	    for (r=i+1;r<j;r++) {
	      if (I51[indx2[r]+i]+f5[indx2[j]+r+1] < score)
	        score=I51[indx2[r]+i]+f5[indx2[j]+r+1];
	    }
	    for (r=i;r<j;r++)
	      for (s=r+TURN;s<=j;s++) {
	        type2 = BP_pair[SEQ[r]][SEQ[s]];
	        if (type2>2) G=PARS->TerminalAU; else G=0;
	        if (type2) {
	          if (r-i>0) left=I5[indx2[r-1]+i]; else left=0;
	          if (j-s>0) right=f5[indx2[j]+s+1]; else right=0;
	          if (left+right+Ibc[indx2[s]+r]+G < score)
	            score=left+right+Ibc[indx2[s]+r]+G;
	        }
	    }
	    I5[indx2[j]+i]=score;
// calculate I5

	    score=MAXENG;
	    for (r=i+1;r<j;r++) {
	      if (Gvm[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]+BETAAm < score) 
	        score=Gvm[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]+BETAAm;
	    }
//I'->V1 V2 (shadow A)
	    for (r=i+1;r<j;r++)
	      for (k1=r+2;k1<j;k1++)
	        for (s=k1+2;s+1<j;s++) {
	          left=Gvm[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=Gump[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          G=BETACm+BETA2;
	          if (left+right+G < score) score=left+right+G;
//I'->V1 U1 V2 U2 (shadow C)

	          left=Gumm[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=Gwm[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          if (left+right+BETA2 < score) score=left+right+BETA2;
//I'->W1 U1 W2 U2 (shadow B & D)
	    }
	    IML1[indx2[j]+i]=score;
// calculate IML1

	    score=fML[indx2[j]+i];
	    if (IML1[indx2[j]+i]<score) score = IML1[indx2[j]+i];
	    for (r=i+1;r<j;r++) {
	      right=fML[indx2[j]+r+1];
	      if ((j-r)*PARS->MLbase < right) right=(j-r)*PARS->MLbase;
	      if (IML1[indx2[r]+i]+right < score)
	        score=IML1[indx2[r]+i]+right;
	    }
	    for (r=i;r<j;r++)
	      for (s=r+TURN;s<=j;s++) {
	        type=BP_pair[SEQ[r]][SEQ[s]];
	        if (type==0) continue;
	        G=PARS->MLintern[type];
	        if (r>i) {
	          left=IML[indx2[r-1]+i];
	          if ((r-i)*PARS->MLbase < left) left=(r-i)*PARS->MLbase;
	        }else left=0;
	        if (s<j) {
	          right=fML[indx2[j]+s+1];
	          if ((j-s)*PARS->MLbase < right) right=(j-s)*PARS->MLbase;
	        } else right=0;
	        if (G+left+right+Ibc[indx2[s]+r] < score) 
	          score=G+left+right+Ibc[indx2[s]+r];
	    }
	    IML[indx2[j]+i]=score;
// calculate IML

	    score=MAXENG;
	    for (r=i+1;r<j;r++) {
	      if (Gvp[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]+BETAAp < score) 
	        score=Gvp[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]+BETAAp;
	    }
//I'->V1 V2 (shadow A)
	    for (r=i+1;r<j;r++)
	      for (k1=r+2;k1<j;k1++)
	        for (s=k1+2;s+1<j;s++) {
	          left=Gvp[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=Gup[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          G=BETACp+BETA2;
	          if (left+right+G < score) score=left+right+G;
//I'->V1 U1 V2 U2 (shadpw C)

	          left=Gup[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=Gwp[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          if (left+right+BETA2 < score) score=left+right+BETA2;
//I'->W1 U1 W2 U2 (shadow B & D)
	    }
	    IPL1[indx2[j]+i]=score;
// calculate IPL1

	    score=fPL[indx2[j]+i];
	    if (IPL1[indx2[j]+i]<score) score = IPL1[indx2[j]+i];
	    for (r=i+1;r<j;r++) {
	      right=fPL[indx2[j]+r+1];
	      if ((j-r)*BETA3 < right) right=(j-r)*BETA3;
	      if (IPL1[indx2[r]+i]+right < score)
	        score=IPL1[indx2[r]+i]+right;
	    }
	    for (r=i;r<j;r++)
	      for (s=r+TURN;s<=j;s++) {
	        type=BP_pair[SEQ[r]][SEQ[s]];
	        if (type==0) continue;
	        G=BETA2;
	        if (r>i) { 
	          left=IPL[indx2[r-1]+i];
	          if ((r-i)*BETA3 < left) left = (r-i)*BETA3;
	        } else left=0;
	        if (s<j) {
	          right=fPL[indx2[j]+s+1];
	          if ((j-s)*BETA3 < right) right = (j-s)*BETA3; 
	        } else right=0;
	        if (G+left+right+Ibc[indx2[s]+r] < score) 
	          score=G+left+right+Ibc[indx2[s]+r];
	    }
	    IPL[indx2[j]+i]=score;
// calculate IPL

	}
}
void de_I(block *T, const char *seq)
{
	int k1,k2,k3,G,V,left,right;
	int type, type2, typek, flag1,flag2;

	if (T->value==0) return;

	switch (T->type) {
	case 0:   //decompose Ibc
	  type = BP_pair[SEQ[T->i]][SEQ[T->j]];
	  if (type) {
	    ptable[T->i]=T->j;
	    ptable[T->j]=T->i;
	    for (k1=T->i+1;k1<T->j;k1++)
	      for (k2=k1+TURN;k2<T->j;k2++) {
	        type2=BP_pair[SEQ[k2]][SEQ[k1]];
	        if (type2) {
	          G=LoopEnergy(k1-T->i-1,T->j-k2-1,type,type2,
		  SEQ[T->i+1],SEQ[T->j-1],SEQ[k1-1],SEQ[k2+1]);
	          if (T->value==G+Ibc[indx2[k2]+k1]) {
	            pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	            return;
	          }

	          G=PARS->MLclosing+PARS->MLintern[type]+PARS->MLintern[type2];
	          if (k1-T->i>1) {
	            left=IML[indx2[k1-1]+T->i+1]; 
	            right=(T->j-k2-1)*PARS->MLbase;
	            if (T->value==left+right+G+Ibc[indx2[k2]+k1]) {
	              pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	              pushblock(T->i+1,k1-1,0,0,2,left);
	              return;
	            }
	          }

	          if (T->j-k2>1) {
	            left=(k1-T->i-1)*PARS->MLbase;
	            right=fML[indx2[T->j-1]+k2+1]; 
	            if (T->value==left+right+G+Ibc[indx2[k2]+k1]) {
	              pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	              pushblock(k2+1,T->j-1,0,0,12,right);
	              return;
	            }
	          }

	          if (k1-T->i>1 && T->j-k2>1) {
	            left=IML[indx2[k1-1]+T->i+1];
	            right=fML[indx2[T->j-1]+k2+1]; 
	            if (T->value==left+right+G+Ibc[indx2[k2]+k1]) {
                      pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	              pushblock(T->i+1,k1-1,0,0,2,left);
	              pushblock(k2+1,T->j-1,0,0,12,right);
	              return;
	            }
	          }
	        }
	    }
	    G=PARS->MLclosing+PARS->MLintern[type];
	    for (k1=T->i+TURN;k1<=T->j-1;k1++) {
	      left=IML1[indx2[k1]+T->i+1];
	      if (k1<T->j-1) {
	        right=fML[indx2[T->j-1]+k1+1];
	        flag2=1;
	        if ((T->j-1-k1)*PARS->MLbase < right) {
	          right = (T->j-1-k1)*PARS->MLbase;
	          flag2=0;
	        }
	      } else 
	        right=0;
	      if (T->value==G+left+right) {
	        pushblock(T->i+1,k1,0,0,5,left);
	        if (k1<T->j-1 && flag2) pushblock(k1+1,T->j-1,0,0,12,right);
	        return;
	      }
	    }
	  } else 
	    return;
	  break;
	case 1:  //decompose I5
	  if (T->value==f5[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,11,f5[indx2[T->j]+T->i]);
	    return;
	  }
	  if (T->value==I51[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,4,I51[indx2[T->j]+T->i]);
	    return;
	  }
	  for (k1=T->i+1;k1<T->j;k1++) {
	    if (T->value==I51[indx2[k1]+T->i]+f5[indx2[T->j]+k1+1]) {
	      pushblock(T->i,k1,0,0,4,I51[indx2[k1]+T->i]);
	      pushblock(k1+1,T->j,0,0,11,f5[indx2[T->j]+k1+1]);
	      return;
	    }
	  }
	  for (k1=T->i;k1<T->j;k1++)
	    for (k2=k1+TURN; k2<=T->j;k2++) {
	      type2 = BP_pair[SEQ[k1]][SEQ[k2]];
	      if (type2>2) G=PARS->TerminalAU; else G=0;
	      if (k1-T->i>0) left=I5[indx2[k1-1]+T->i]; else left=0;
	      if (T->j-k2>0) right=f5[indx2[T->j]+k2+1]; else right=0;
	      if (T->value==G+left+right+Ibc[indx2[k2]+k1]) {
	        if (k1-T->i>0) pushblock(T->i,k1-1,0,0,1,left);
	        if (T->j-k2>0) pushblock(k2+1,T->j,0,0,11,right);
	        pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	        return;
	      }
	  }
	  break;
	case 2:  //decompose IML
	  if (T->value==fML[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,12,fML[indx2[T->j]+T->i]);
	    return;
	  }
	  if (T->value==IML1[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,5,IML1[indx2[T->j]+T->i]);
	    return;
	  }
	  for (k1=T->i+1;k1<T->j;k1++) {
	    if (T->value==IML1[indx2[k1]+T->i]+fML[indx2[T->j]+k1+1]) {
	      pushblock(T->i,k1,0,0,5,IML1[indx2[k1]+T->i]);
	      pushblock(k1+1,T->j,0,0,12,fML[indx2[T->j]+k1+1]);
	      return;
	    }
	    if (T->value==IML1[indx2[k1]+T->i]+(T->j-k1)*PARS->MLbase) {
	      pushblock(T->i,k1,0,0,5,IML1[indx2[k1]+T->i]);
	      return;
	    }
	  }
	  for (k1=T->i;k1<T->j;k1++)
	    for (k2=k1+TURN;k2<=T->j;k2++) {
	      type=BP_pair[SEQ[k1]][SEQ[k2]];
	      if (type==0) continue;
	      G=PARS->MLintern[type];
	      if (k1>T->i) left=IML[indx2[k1-1]+T->i]; else left=0;
	      if (k2<T->j) right=fML[indx2[T->j]+k2+1]; else right=0;
	      if (k1>T->i && k2<T->j) {
	        if (T->value==G+left+right+Ibc[indx2[k2]+k1]) {
	          pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	          pushblock(T->i,k1-1,0,0,2,left);
	          pushblock(k2+1,T->j,0,0,12,right);
	          return;
	        }
	      }
	      if (k1>T->i) {
	        if (T->value==G+left+(T->j-k2)*PARS->MLbase+Ibc[indx2[k2]+k1]) {
	          pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	          pushblock(T->i,k1-1,0,0,2,left);
	          return;
	        }
	      }
	      if (k2<T->j) {
	        if (T->value==G+(k1-T->i)*PARS->MLbase+right+Ibc[indx2[k2]+k1]) {
	          pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	          pushblock(k2+1,T->j,0,0,12,right);
	          return;
	        }
	      }
	      if (T->value==G+(k1-T->i)*PARS->MLbase+(T->j-k2)*PARS->MLbase+Ibc[indx2[k2]+k1]) {
	        pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	        return;
	      }
	  }
	  break;
	case 3: //decompose IPL
	  if (T->value==fPL[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,13,fPL[indx2[T->j]+T->i]);
	    return;
	  }
	  if (T->value==IPL1[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,6,IPL1[indx2[T->j]+T->i]);
	    return;
	  }
	  for (k1=T->i+1;k1<T->j;k1++) {
	    if (T->value==IPL1[indx2[k1]+T->i]+fPL[indx2[T->j]+k1+1]) {
	      pushblock(T->i,k1,0,0,6,IPL1[indx2[k1]+T->i]);
	      pushblock(k1+1,T->j,0,0,13,fPL[indx2[T->j]+k1+1]);
	      return;
	    }
	    if (T->value==IPL1[indx2[k1]+T->i]+(T->j-k1)*BETA3) {
	      pushblock(T->i,k1,0,0,6,IPL1[indx2[k1]+T->i]);
	      return;
	    }
	  }
	  for (k1=T->i;k1<T->j;k1++)
	    for (k2=k1+TURN;k2<=T->j;k2++) {
	      type=BP_pair[SEQ[k1]][SEQ[k2]];
	      if (type==0) continue;
	      G=BETA2;
	      if (k1>T->i) left=IPL[indx2[k1-1]+T->i]; else left=0;
	      if (k2<T->j) right=fPL[indx2[T->j]+k2+1]; else right=0;
	      if (k1>T->i && k2<T->j) {
	        if (T->value==G+left+right+Ibc[indx2[k2]+k1]) {
	          pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	          pushblock(T->i,k1-1,0,0,3,left);
	          pushblock(k2+1,T->j,0,0,13,right);
	          return;
	        }
	      }
	      if (k1>T->i) {
	        if (T->value==G+left+(T->j-k2)*BETA3+Ibc[indx2[k2]+k1]) {
	          pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	          pushblock(T->i,k1-1,0,0,3,left);
	          return;
	        }
	      }
	      if (k2<T->j) {
	        if (T->value==G+(k1-T->i)*BETA3+right+Ibc[indx2[k2]+k1]) {
	          pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	          pushblock(k2+1,T->j,0,0,13,right);
	          return;
	        }
	      }
	      if (T->value==G+(k1-T->i)*BETA3+(T->j-k2)*BETA3+Ibc[indx2[k2]+k1]) {
	        pushblock(k1,k2,0,0,0,Ibc[indx2[k2]+k1]);
	        return;
	      }
	  }
	  break;
	case 4: //decompose I51
	  for (k1=T->i+1;k1<T->j;k1++) {
	    if (T->value==Gv[indx4[0][T->i]+indx4[1][k1]+indx4[2][k1+1]+indx4[3][T->j]]+BETAA) {
	      pushblock(T->i,T->j,k1,k1+1,26,T->value-BETAA);
	      return;
	    }
	  }
	  for (k1=T->i+1;k1<T->j;k1++)
	    for (k2=k1+2;k2<T->j;k2++)
	      for (k3=k2+2;k3+1<T->j;k3++) {
	        left=Gv[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=Guep[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        G=BETAC+BETA2;
	        if (T->value==left+right+G) {
	          pushblock(T->i,k3,k1,k2+1,26,left);
	          pushblock(k1+1,T->j,k2,k3+1,22,right);
	          return;
	        }
	        left=Gu[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=Gw[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        if (T->value==left+right+BETA2) {
	          pushblock(T->i,k3,k1,k2+1,21,left);
	          pushblock(k1+1,T->j,k2,k3+1,29,right);
	          return;
	        }
	  }
	  break;
	case 5: //decompose IML1
	  for (k1=T->i+1;k1<T->j;k1++) {
	    if (T->value==Gvm[indx4[0][T->i]+indx4[1][k1]+indx4[2][k1+1]+indx4[3][T->j]]+BETAAm) {
	      pushblock(T->i,T->j,k1,k1+1,27,T->value-BETAAm);
	      return;
	    }
	  }
	  for (k1=T->i+1;k1<T->j;k1++)
	    for (k2=k1+2;k2<T->j;k2++)
	      for (k3=k2+2;k3+1<T->j;k3++) {
	        left=Gvm[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=Gump[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        G=BETACm+BETA2;
	        if (T->value==left+right+G) {
	          pushblock(T->i,k3,k1,k2+1,27,left);
	          pushblock(k1+1,T->j,k2,k3+1,24,right);
	          return;
	        }
	        left=Gumm[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=Gwm[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        if (T->value==left+right+BETA2) {
	          pushblock(T->i,k3,k1,k2+1,23,left);
	          pushblock(k1+1,T->j,k2,k3+1,30,right);
	          return;
	        }
	  }
	  break;
	case 6: //decompose IPL1
	  for (k1=T->i+1;k1<T->j;k1++) {
	    if (T->value==Gvp[indx4[0][T->i]+indx4[1][k1]+indx4[2][k1+1]+indx4[3][T->j]]+BETAAp) {
	      pushblock(T->i,T->j,k1,k1+1,28,T->value-BETAAp);
	      return;
	    }
	  }
	  for (k1=T->i+1;k1<T->j;k1++)
	    for (k2=k1+2;k2<T->j;k2++)
	      for (k3=k2+2;k3+1<T->j;k3++) {
	        left=Gvp[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=Gup[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        G=BETACp+BETA2;
	        if (T->value==left+right+G) {
	          pushblock(T->i,k3,k1,k2+1,28,left);
	          pushblock(k1+1,T->j,k2,k3+1,25,right);
	          return;
	        }
	        left=Gup[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=Gwp[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        if (T->value==left+right+BETA2) {
	          pushblock(T->i,k3,k1,k2+1,25,left);
	          pushblock(k1+1,T->j,k2,k3+1,31,right);
	          return;
	        }
	  }
	  break;
	case 10: //decompose c
	  type=BP_pair[SEQ[T->i]][SEQ[T->j]];
	  if (type) {
	    ptable[T->i]=T->j;
	    ptable[T->j]=T->i;
	    G=HairpinE(T->j-T->i-1, type, SEQ[T->i+1], SEQ[T->j-1], seq+T->i-1);
	    if (T->value==G) {
	      return;
	    }
	    for (k1=T->i+1;k1<T->j;k1++)
	      for (k2=T->i+TURN; k2<T->j;k2++) {
	        type2=BP_pair[SEQ[k2]][SEQ[k1]];
	        if (type2==0) continue;
	        G=LoopEnergy(k1-T->i-1, T->j-k2-1, type, type2,
		SEQ[T->i+1], SEQ[T->j-1], SEQ[k1-1], SEQ[k2+1]);
	        if (T->value==G+c[indx2[k2]+k1]) {
	          pushblock(k1,k2,0,0,10,c[indx2[k2]+k1]);
	          return;
	        }
	    }
	    for (k1=T->i+1;k1<T->j-1;k1++) {
	      G=PARS->MLintern[type]+PARS->MLclosing;
	      if (T->value==fML1[indx2[k1]+T->i+1]+fML[indx2[T->j-1]+k1+1]+G) {
	        pushblock(T->i+1,k1,0,0,15,fML1[indx2[k1]+T->i+1]);
	        pushblock(k1+1,T->j-1,0,0,12,fML[indx2[T->j-1]+k1+1]);
	        return;
	      } 
	    }
	  } else 
	    return;
	  break;
	case 11: //decompose f5
	  if (T->value==f51[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,14,f51[indx2[T->j]+T->i]);
	    return;
	  }
	  for (k1=T->i;k1<T->j;k1++) {
	    if (T->value==f51[indx2[k1]+T->i]+f5[indx2[T->j]+k1+1]) {
	      pushblock(T->i,k1,0,0,14,f51[indx2[k1]+T->i]);
	      pushblock(k1+1,T->j,0,0,11,f5[indx2[T->j]+k1+1]);
	      return;
	    }
	  }
	  break;
	case 12: //decompose fML
	  if (T->value==fML1[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,15,fML1[indx2[T->j]+T->i]);
	    return;
	  }
	  for (k1=T->i;k1<T->j;k1++) {
	    if (T->value==fML1[indx2[k1]+T->i]+fML[indx2[T->j]+k1+1]) {
	      pushblock(T->i,k1,0,0,15,fML1[indx2[k1]+T->i]);
	      pushblock(k1+1,T->j,0,0,12,fML[indx2[T->j]+k1+1]);
	      return;
	    }
	    if (T->value==fML1[indx2[k1]+T->i]+(T->j-k1)*PARS->MLbase) {
	      pushblock(T->i,k1,0,0,15,fML1[indx2[k1]+T->i]);
	      return;
	    }
	  }
	  break;
	case 13: //decompose fPL
	  if (T->value==fPL1[indx2[T->j]+T->i]) {
	    pushblock(T->i,T->j,0,0,16,fPL1[indx2[T->j]+T->i]);
	    return;
	  }
	  for (k1=T->i;k1<T->j;k1++) {
	    if (T->value==fPL1[indx2[k1]+T->i]+fPL[indx2[T->j]+k1+1]) {
	      pushblock(T->i,k1,0,0,16,fPL1[indx2[k1]+T->i]);
	      pushblock(k1+1,T->j,0,0,13,fPL[indx2[T->j]+k1+1]);
	      return;
	    }
	    if (T->value==fPL1[indx2[k1]+T->i]+(T->j-k1)*BETA3) {
	      pushblock(T->i,k1,0,0,16,fPL1[indx2[k1]+T->i]);
	      return;
	    }
	  }
	  break;
	case 14: //decompose f51
	  for (k1=T->i;k1<T->j;k1++) {
	    type = BP_pair[SEQ[k1]][SEQ[T->j]];
	    if (type>2) G=PARS->TerminalAU; else G=0;
	    if (type && T->value==c[indx2[T->j]+k1]+G) {
	      pushblock(k1,T->j,0,0,10,c[indx2[T->j]+k1]);
	      return;
	    }
	  }
	  break;
	case 15: //decompose fML1
	  for (k1=T->i;k1<T->j;k1++) {
	    type = BP_pair[SEQ[k1]][SEQ[T->j]];
	    if (type) {
	      G=PARS->MLintern[type]+(k1-T->i)*PARS->MLbase;
	      if (T->value==c[indx2[T->j]+k1]+G) {
	        pushblock(k1,T->j,0,0,10,c[indx2[T->j]+k1]);
	        return;
	      }
	    }
	  }
	  break;
	case 16: //decompose fPL1
	  for (k1=T->i;k1<T->j;k1++) {
	    type = BP_pair[SEQ[k1]][SEQ[T->j]];
	    if (type) {
	      G=BETA2+(k1-T->i)*BETA3;
	      if (T->value==c[indx2[T->j]+k1]+G) {
	        pushblock(k1,T->j,0,0,10,c[indx2[T->j]+k1]);
	        return;
	      }
	    }
	  }
	  break;
	case 20: //decompose gtight
	  type=BP_pair[SEQ[T->i]][SEQ[T->j]];
	  type2=BP_pair[SEQ[T->s]][SEQ[T->r]];
	  if (type==0 || type2==0) return;
	  ptable[T->i]=T->j;
	  ptable[T->j]=T->i;
	  if (type2 && T->s-T->r>TURN) {
	    G=(int)(SIGMA*LoopEnergy(T->r-T->i-1,T->j-T->s-1,type,type2,
		  SEQ[T->i+1],SEQ[T->j-1],SEQ[T->r-1],SEQ[T->s+1]));
	    if (T->value==G) {
	      ptable[T->r]=T->s;
	      ptable[T->s]=T->r;
	      return;
	    }
	    for (k1=T->i+1;k1<T->r;k1++)
	      for (k2=T->s+1;k2<T->j;k2++) {
	        typek=BP_pair[SEQ[k2]][SEQ[k1]];
	        if (typek) {
	          G=(int)(SIGMA*LoopEnergy(k1-T->i-1,T->j-k2-1,type,typek,
		SEQ[T->i+1],SEQ[T->j-1],SEQ[k1-1],SEQ[k2+1]));
	          V=Gtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][k2]];
	          if (T->value==G+V) {
	            pushblock(k1,k2,T->r,T->s,20,V);
	            return;
	          }
	          G=PARS->MLclosing+PARS->MLintern[type];
	          if (k1-T->i>1) {
	            if (T->value==G+V+IML[indx2[k1-1]+T->i+1]) {
	               pushblock(k1,k2,T->r,T->s,20,V);
	               pushblock(T->i+1,k1-1,0,0,2,IML[indx2[k1-1]+T->i+1]);
	             return;
	            }
	            if (T->value==G+V+(T->j-k2-1)*PARS->MLbase) {
	              pushblock(k1,k2,T->r,T->s,20,V);
	              return;
	            }
	          } 
	          if (T->j-k2>1) {
	            if (T->value==G+V+IML[indx2[T->j-1]+k2+1]) {
	              pushblock(k1,k2,T->r,T->s,20,V);
	              pushblock(k2+1,T->j-1,0,0,2,IML[indx2[T->j-1]+k2+1]);
	              return;
	            }
	            if (T->value==G+V+(k1-T->i-1)*PARS->MLbase) {
	              pushblock(k1,k2,T->r,T->s,20,V);
	              return;
	            }
	          }
	          if (k1-T->i>1 && T->j-k2>1) {
	            if (T->value==G+V+IML[indx2[k1-1]+T->i+1]+IML[indx2[T->j-1]+k2+1]) {
	              pushblock(k1,k2,T->r,T->s,20,V);
	              pushblock(T->i+1,k1-1,0,0,2,IML[indx2[k1-1]+T->i+1]);
	              pushblock(k2+1,T->j-1,0,0,2,IML[indx2[T->j-1]+k2+1]);
	              return;
	            }
	          }
	        }
	    }
	  }
	  break;
	case 21: //decompose Gu
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]>2) G=PARS->TerminalAU; else G=0;
	      if (k1-T->i>0) left=I5[indx2[k1-1]+T->i]; else left=0;
	      if (k2-T->s>0) {
	        right=IPL[indx2[k2-1]+T->s]; 
	        flag2=1;
	        if ((k2-T->s)*BETA3 < right) {right=(k2-T->s)*BETA3; flag2=0;}
	      } else 
	        right=0;
	      V=Gtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]];
	      if (T->value==left+right+V+G) {
	        pushblock(k1,T->j,T->r,k2,20,V);
	        if (k1-T->i>0) pushblock(T->i,k1-1,0,0,1,I5[indx2[k1-1]+T->i]);
	        if (k2-T->s>0 && flag2) pushblock(T->s,k2-1,0,0,3,IPL[indx2[k2-1]+T->s]);
	        return;
	      }
	  }
	  break;
	case 22: //decompose Guep
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]>2) G=PARS->TerminalAU; else G=0;
	      if (k1-T->i>0) {
	        left=IPL[indx2[k1-1]+T->i]; flag1=1;
	        if ((k1-T->i)*BETA3 < left) {left=(k1-T->i)*BETA3;flag1=0;}
	      } else 
	        left=0;
	      if (k2-T->s>0) {
	        right=IPL[indx2[k2-1]+T->s]; flag2=1;
	        if ((k2-T->s)*BETA3 < right) {right=(k2-T->s)*BETA3;flag2=0;}
	      } else 
	        right=0;
	      V=Gtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]];
	      if (T->value==left+right+V+G) {
	        pushblock(k1,T->j,T->r,k2,20,V);
	        if (k1-T->i>0 && flag1) pushblock(T->i,k1-1,0,0,3,IPL[indx2[k1-1]+T->i]);
	        if (k2-T->s>0 && flag2) pushblock(T->s,k2-1,0,0,3,IPL[indx2[k2-1]+T->s]);
	        return;
	      }
	  }
	  break;
	case 23: //decompose Gumm
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      if (k1-T->i>0) {
	        left=IML[indx2[k1-1]+T->i]; flag1=1;
	        if ((k1-T->i)*PARS->MLbase < left) {left=(k1-T->i)*PARS->MLbase;flag1=0;}
	      } else 
	        left=0;
	      if (k2-T->s>0) {
	        right=IPL[indx2[k2-1]+T->s]; flag2=1;
	        if ((k2-T->s)*BETA3 < right) {right=(k2-T->s)*BETA3;flag2=0;}
	      } else 
	        right=0;
	      typek=BP_pair[SEQ[k1]][SEQ[T->j]];
	      G=PARS->MLintern[typek];
	      V=Gtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]];
	      if (T->value==left+right+V+G) {
	        pushblock(k1,T->j,T->r,k2,20,V);
	        if (k1-T->i>0 && flag1) pushblock(T->i,k1-1,0,0,2,IML[indx2[k1-1]+T->i]);
	        if (k2-T->s>0 && flag2) pushblock(T->s,k2-1,0,0,3,IPL[indx2[k2-1]+T->s]);
	        return;
	      }
	  }
	  break;
	case 24: //decompose Gump
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      if (k1-T->i>0) {
	        left=IPL[indx2[k1-1]+T->i]; flag1=1;
	        if ((k1-T->i)*BETA3 < left) {left=(k1-T->i)*BETA3;flag1=0;}
	      } else 
	        left=0;
	      if (k2-T->s>0) {
	        right=IPL[indx2[k2-1]+T->s]; flag2=1;
	        if ((k2-T->s)*BETA3 < right) {right=(k2-T->s)*BETA3;flag2=0;}
	      } else 
	        right=0;
	      typek=BP_pair[SEQ[k1]][SEQ[T->j]];
	      G=PARS->MLintern[typek];
	      V=Gtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]];
	      if (T->value==left+right+V+G) {
	        pushblock(k1,T->j,T->r,k2,20,V);
	        if (k1-T->i>0 && flag1) pushblock(T->i,k1-1,0,0,3,IPL[indx2[k1-1]+T->i]);
	        if (k2-T->s>0 && flag2) pushblock(T->s,k2-1,0,0,3,IPL[indx2[k2-1]+T->s]);
	        return;
	      }
	  }
	  break;
	case 25: //decompose Gup
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      if (k1-T->i>0) {
	        left=IPL[indx2[k1-1]+T->i]; flag1=1;
	        if ((k1-T->i)*BETA3 < left) {left=(k1-T->i)*BETA3;flag1=0;}
	      } else 
	        left=0;
	      if (k2-T->s>0) {
	        right=IPL[indx2[k2-1]+T->s]; flag2=1;
	        if ((k2-T->s)*BETA3 < right) {right=(k2-T->s)*BETA3;flag2=0;}
	      } else 
	        right=0;
	      G=BETA2;
	      V=Gtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]];
	      if (T->value==left+right+V+G) {
	        pushblock(k1,T->j,T->r,k2,20,V);
	        if (k1-T->i>0 && flag1) pushblock(T->i,k1-1,0,0,3,IPL[indx2[k1-1]+T->i]);
	        if (k2-T->s>0 && flag2) pushblock(T->s,k2-1,0,0,3,IPL[indx2[k2-1]+T->s]);
	        return;
	      }
	  }
	  break;
	case 26: //decompose Gv
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=Gu[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=Guep[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      if (T->value==left+right+2*BETA2) {
	        pushblock(T->i,k2,k1,T->s,21,left);
	        pushblock(k1+1,T->j,T->r,k2+1,22,right);
	        return;
	      }
	  }
	  break;
	case 32: //decompose Gvep
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=Guep[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=Guep[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      if (T->value==left+right+2*BETA2) {
	        pushblock(T->i,k2,k1,T->s,22,left);
	        pushblock(k1+1,T->j,T->r,k2+1,22,right);
	        return;
	      }
	  }
	  break;
	case 27: //decompose Gvm
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=Gumm[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=Gump[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      if (T->value==left+right+2*BETA2) {
	        pushblock(T->i,k2,k1,T->s,23,left);
	        pushblock(k1+1,T->j,T->r,k2+1,24,right);
	        return;
	      }
	  }
	  break;
	case 33: //decompose Gvm
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=Gump[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=Gump[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      if (T->value==left+right+2*BETA2) {
	        pushblock(T->i,k2,k1,T->s,24,left);
	        pushblock(k1+1,T->j,T->r,k2+1,24,right);
	        return;
	      }
	  }
	  break;
	case 28: //decompose Gvp
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=Gup[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=Gup[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      if (T->value==left+right+2*BETA2) {
	        pushblock(T->i,k2,k1,T->s,25,left);
	        pushblock(k1+1,T->j,T->r,k2+1,25,right);
	        return;
	      }
	  }
	  break;
	case 29: //decompose Gw
	  for (k1=T->s+1;k1<T->j;k1++)
	    for (k2=k1+3;k2<T->j;k2++) {
	      left=Guep[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=Guep[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETAB+2*BETA2; // the penalty for the crossing matrix (shadow X)
	      if (T->value==left+right+G) {
	        pushblock(T->i,k2-1,T->r,k1+1,22,left);
	        pushblock(T->s,T->j,k1,k2,22,right);
	        return;
	      }

	      left=Gvep[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=Guep[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETAD+BETA2; // the penalty for the crossing matrix (shadow X)
	      if (T->value==left+right+G) {
	        pushblock(T->i,k2-1,T->r,k1+1,32,left);
	        pushblock(T->s,T->j,k1,k2,22,right);
	        return;
	      }
	  }
	  break;
	case 30: //decompose Gwm
	  for (k1=T->s+1;k1<T->j;k1++)
	    for (k2=k1+3;k2<T->j;k2++) {
	      left=Gump[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=Gump[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETABm+2*BETA2; // the penalty for the crossing matrix (shadow X)
	      if (T->value==left+right+G) {
	        pushblock(T->i,k2-1,T->r,k1+1,24,left);
	        pushblock(T->s,T->j,k1,k2,24,right);
	        return;
	      }

	      left=Gvmp[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=Gump[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETADm+BETA2; // the penalty for the crossing matrix (shadow X)
	      if (T->value==left+right+G) {
	        pushblock(T->i,k2-1,T->r,k1+1,33,left);
	        pushblock(T->s,T->j,k1,k2,24,right);
	        return;
	      }
	  }
	  break;
	case 31: //decompose Gwp
	  for (k1=T->s+1;k1<T->j;k1++)
	    for (k2=k1+3;k2<T->j;k2++) {
	      left=Gup[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=Gup[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETABp+2*BETA2; // the penalty for the crossing matrix (shadow X)
	      if (T->value==left+right+G) {
	        pushblock(T->i,k2-1,T->r,k1+1,25,left);
	        pushblock(T->s,T->j,k1,k2,25,right);
	        return;
	      }

	      left=Gvp[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=Gup[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETADp+BETA2; // the penalty for the crossing matrix (shadow X)
	      if (T->value==left+right+G) {
	        pushblock(T->i,k2-1,T->r,k1+1,28,left);
	        pushblock(T->s,T->j,k1,k2,25,right);
	        return;
	      }
	  }
	  break;
	default : break;
	}
	fprintf(fo, "Error! Do not match the energy, please check your program carefully!\n");
	fprintf(fo, "%d\n", T->type);
	fprintf(fo, "%d\n", T->value);
}



void symbFold(char *seq)
{
	int length,i,j,r,s,si,sj,sr,ss,correct, total,total_nat;
	int score=0,G,type1,type2,new_score;
	long index;
	char *struc;
	block *T, *IN;

	length=strlen(seq);

	initial(length);
	update_fold_params();
	if (!PARS) PARS = scale_parameters();
	encode_seq(seq);
	fillarray(seq);

	pushblock(1,length,0,0,1,I5[indx2[length]+1]);
	while (Bstack!=NULL) {
	  T=popblock();
	  de_I(T,seq);
		free(T);
// 	  printf("-----------------------------------\n");
// 	  IN=Bstack;
// 	  while (IN!=NULL) {
// 	    printf("(%d %d %d %d) %d %d\n", IN->i, IN->j, IN->r, IN->s, IN->type,IN->value);
// 	    IN=IN->Next;
// 	  }
	}

	ptable[0]=length;

	fprintf(fo, "MFE: %d\n", I5[indx2[length]+1]);
	mfe=I5[indx2[length]+1];

	struc=pair2structure(ptable);
	if (p_nat) {
		correct=compare(ptable, p_nat);
		total=number(ptable);
		total_nat=number(p_nat);
		fprintf(fo, "%s  CORRECT: %d  MFEsen: %f PPV: %f\n", struc, correct,
		 (double)correct/(double)total_nat, (double)correct/(double)total);
	}
// 	for (i=0;i<=length;i++)
// 	  printf("%d  ", ptable[i]);
// 	printf("\n");

	free(struc);
	freevar();
}



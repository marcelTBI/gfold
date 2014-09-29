#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "head.h"

void pfunc_initial (int n)
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

	Qc     = (double *) space(sizeof(double)*((n*(n+1))/2+2));
	Qf5    = (double *) space(sizeof(double)*((n*(n+1))/2+2));
	Qf51   = (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QfML   = (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QfML1  = (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QfPL   = (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QfPL1  = (double *) space(sizeof(double)*((n*(n+1))/2+2));

	QIbc= (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QI51= (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QI5= (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QIML1= (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QIML= (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QIPL1= (double *) space(sizeof(double)*((n*(n+1))/2+2));
	QIPL= (double *) space(sizeof(double)*((n*(n+1))/2+2));

	QGtight=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGh=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGu=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGuep=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGv=(double *) space (sizeof(double)*(indx4[3][n+1]));	
	QGw=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGumm=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGump=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGvm=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGwm=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGup=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGvp=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGwp=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGvep=(double *) space (sizeof(double)*(indx4[3][n+1]));
	QGvmp=(double *) space (sizeof(double)*(indx4[3][n+1]));

	Basepair=(int *) space (sizeof(int)*((n*(n+1))/2+2));
	for (i=0;i<(n*(n+1))/2+2;i++) Basepair[i]=0;
	

	sample=(int *) space (sizeof(int)*(n+2));
// 	printf("Require %dGb\n", ((14*(n*(n+1))/2+2)*sizeof(double)
// 		+15*indx4[3][n+1]*sizeof(double))/1024/1024/1024);
}

void pfunc_freevar()
{
	int i;
	free(QGtight);
	free(QGh);
	free(QGu); free(QGuep); free(QGv); free(QGw);
	free(QGumm); free(QGump), free(QGvm); free(QGwm);
	free(QGup); free(QGvp); free(QGwp);
	free(QGvep); free(QGvmp);

	free(Qc); free(QfML); free(QfML1); free(QfPL1); free(QfPL);
	free (QIbc); free(Qf5); free(Qf51); 
	free(QI51); free(QI5);
	free(QIML1); free(QIML);
	free(QIPL1); free(QIPL);
	
	for (i=0;i<4;i++)
	  free(indx4[i]);
	free(indx2);
	free(Basepair);
}


void partition(char *seq) 
{
	int i,j,r,s,k1,k2,p,q,gap;
	int length,type1,type2,typek,G,type;
	double left,right;
	long index, indexk;

	length=strlen(seq);
	for (i=1;i<=length;i++)
	  for (j=i;j<i+TURN;j++) {
	    Qc[indx2[j]+i]=0;
	    Qf5[indx2[j]+i]=1;
	    Qf51[indx2[j]+i]=0;
	    QfML[indx2[j]+i]=0;
	    QfML1[indx2[j]+i]=0;
	    QfPL[indx2[j]+i]=0;
	    QfPL1[indx2[j]+i]=0;
	}
//#pragma omp parallel for
	for (i=length-TURN;i>=1;i--)
	  for (j=i+TURN;j<=length;j++) {

// 	for (gap=TURN;gap<length;gap++)
//         #pragma omp parallel default (shared) private \
//  (i,j,r,s,k1,k2,p,q,type1, type2,typek,G,type,left,right,index,indexk)
//         #pragma omp for schedule (static)
// 	  for (i=1;i<length;i++) {
// 	    j=i+gap;

// printf("%d %d\n", i,j);
	    type = BP_pair[SEQ[i]][SEQ[j]];
	    if (type) {
	      G=HairpinE(j-i-1, type, SEQ[i+1], SEQ[j-1], seq+i-1);
	      Qc[indx2[j]+i]=EXP(G);
	      for (p = i+1; p <= j-2-TURN; p++) 
	        for (q = p+TURN; q < j; q++) {
	          type2 = BP_pair[SEQ[p]][SEQ[q]];

	          if (type2==0) continue;
	          type2 = rtype[type2];

	          G = LoopEnergy(p-i-1, j-q-1, type, type2,
				SEQ[i+1], SEQ[j-1], SEQ[p-1], SEQ[q+1]);
	          Qc[indx2[j]+i]+=EXP(G)*Qc[indx2[q]+p];
	      }
	      for (k1=i+1;k1<j-1;k1++) {
	        G=PARS->MLintern[type]+PARS->MLclosing;
	        Qc[indx2[j]+i]+=QfML1[indx2[k1]+i+1]*QfML[indx2[j-1]+k1+1]*EXP(G);
	      }
	    } else 
	      Qc[indx2[j]+i]=0;

	    Qf51[indx2[j]+i]=0;
	    for (k1=i;k1<j;k1++) {
	      type = BP_pair[SEQ[k1]][SEQ[j]];
	      if (type>2) G=PARS->TerminalAU; else G=0;
	      if (type) Qf51[indx2[j]+i]+=Qc[indx2[j]+k1]*EXP(G);

	    }

	    Qf5[indx2[j]+i]=Qf51[indx2[j]+i]+1;
	    for (k1=i;k1<j;k1++) {
	      Qf5[indx2[j]+i]+=Qf51[indx2[k1]+i]*Qf5[indx2[j]+k1+1];
	    }

	    QfML1[indx2[j]+i]=0;
	    for (k1=i;k1<j;k1++) {
	      type=BP_pair[SEQ[k1]][SEQ[j]];
	      if (type) {
	        G=PARS->MLintern[type]+(k1-i)*PARS->MLbase;
	        QfML1[indx2[j]+i]+=Qc[indx2[j]+k1]*EXP(G);
	      }
	    }

	    QfML[indx2[j]+i]=QfML1[indx2[j]+i];
	    for (k1=i;k1<j;k1++) {
	      G=(j-k1)*PARS->MLbase;
	      QfML[indx2[j]+i]+=QfML1[indx2[k1]+i]*(QfML[indx2[j]+k1+1]+EXP(G));
	    }

	    QfPL1[indx2[j]+i]=0;
	    for (k1=i;k1<j;k1++) {
	      type=BP_pair[SEQ[k1]][SEQ[j]];
	      if (type) {
	        G=BETA2+(k1-i)*BETA3;
	        QfPL1[indx2[j]+i]+=Qc[indx2[j]+k1]*EXP(G);
	      }
	    }

	    QfPL[indx2[j]+i]=QfPL1[indx2[j]+i];
	    for (k1=i;k1<j;k1++) {
	      G=(j-k1)*BETA3;
	      QfPL[indx2[j]+i]+=QfPL1[indx2[k1]+i]*(QfPL[indx2[j]+k1+1]+EXP(G));
	    }
// partition for secondary structure S done!

 	    type1=BP_pair[SEQ[i]][SEQ[j]];
	    for (r=i+1;r<j-1;r++)
	      for (s=j-1;s>r;s--) {
	        type2=BP_pair[SEQ[s]][SEQ[r]];
	        index=indx4[0][i]+indx4[1][r]+indx4[2][s]+indx4[3][j];
	        if (type1 && type2 && s-r>TURN) { // start Gtight
	          G=(int)(SIGMA*LoopEnergy(r-i-1,j-s-1,type1,type2,
		  SEQ[i+1],SEQ[j-1],SEQ[r-1],SEQ[s+1]));
	          QGtight[index]=EXP(G);
// G is empty, it is big interior loop
	          for (k1=i+1;k1<r;k1++)
	            for (k2=s+1;k2<j;k2++) {
	              typek=BP_pair[SEQ[k2]][SEQ[k1]];
	              indexk=indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][k2];

	              if (typek) {
	                G=(int)(SIGMA*LoopEnergy(k1-i-1,j-k2-1,type1,typek,
		SEQ[i+1],SEQ[j-1],SEQ[k1-1],SEQ[k2+1]));
	                QGtight[index]+=QGtight[indexk]*EXP(G);
// G is an interior loop plus another G

	                G=PARS->MLclosing+PARS->MLintern[type1];
	                if (k1-i>1)
	                  QGtight[index]+=QGtight[indexk]*EXP(G)*QIML[indx2[k1-1]+i+1]
		*EXP(((j-k2-1)*PARS->MLbase));
// G is a miltiloop, and the rhs. is empty

	                if (j-k2>1) 
	                  QGtight[index]+=QGtight[indexk]*EXP(G)*EXP(((k1-i-1)*PARS->MLbase))
		*QIML[indx2[j-1]+k2+1];
// G is a miltiloop, and the lhs. is empty

	                if (k1-i>1 && j-k2>1)
 	                  QGtight[index]+=QGtight[indexk]*EXP(G)*QIML[indx2[k1-1]+i+1]
		*QIML[indx2[j-1]+k2+1];
// G is a miltiloop, and both the lhs. and rhs are not empty.3
	              }  
	          }
	        }  else 
	          QGtight[index]=0;
//calculate the G matrix
	        QGh[index]=0;
	        type1=BP_pair[SEQ[i]][SEQ[j]];
	        if (type1) {
	          for (k2=s>r+TURN?s:r+TURN; k2<j; k2++) {
	            type2=BP_pair[SEQ[r]][SEQ[k2]];
	            if (type2==0) continue;
	            if (k2-s>0)
	              right=QIPL[indx2[k2-1]+s]+EXP((k2-s)*BETA3);
	            else 
	              right=1;
	            QGh[index]+=QGtight[indx4[0][i]+indx4[1][r]+indx4[2][k2]+indx4[3][j]]*right;
	          }
	        }
// calculate the Gh matrix (an intermediate for GX, X=u,v,w)

	        QGu[index]=0;
	        for (k1=i;k1<r;k1++) {
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          typek=BP_pair[SEQ[k1]][SEQ[j]];
	          if (typek>2) G=PARS->TerminalAU; else G=0;
	          if (k1-i>0) left=QI5[indx2[k1-1]+i]; else left=1;
	          QGu[index]+=left*QGh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]*EXP(G); 
	        }
//calculate the Gu matrix

	        QGuep[index]=0;
	        for (k1=i;k1<r;k1++) {
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          typek=BP_pair[SEQ[k1]][SEQ[j]];
	          if (typek>2) G=PARS->TerminalAU; else G=0;
	          if (k1-i>0) 
	            left=QIPL[indx2[k1-1]+i]+EXP((k1-i)*BETA3);
	          else
	            left=1;
 	          QGuep[index]+=left*QGh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]*EXP(G); 
	        }
//calculate the Guep matrix

	        QGumm[index]=0;
	        for (k1=i;k1<r;k1++) { 
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          if (k1-i>0) {
	            left=QIML[indx2[k1-1]+i]+EXP((k1-i)*PARS->MLbase);
	          } else 
	            left=1;
	          typek=BP_pair[SEQ[k1]][SEQ[j]];
	          G=PARS->MLintern[typek];
	          QGumm[index]+=left*QGh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]*EXP(G);
	        }
//calculate the Gumm matrix

	        QGump[index]=0;
	        for (k1=i;k1<r;k1++) {
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          if (k1-i>0) {
	            left=QIPL[indx2[k1-1]+i]+EXP((k1-i)*BETA3);
	          } else 
	            left=1;
	          typek=BP_pair[SEQ[k1]][SEQ[j]];
	          G=PARS->MLintern[typek];
	          QGump[index]+=left*QGh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]*EXP(G);
	        }
//calculate the Gump matrix

	        QGup[index]=0;
	        for (k1=i;k1<r;k1++) {
	          if (BP_pair[SEQ[k1]][SEQ[j]]==0) continue;
	          if (k1-i>0) {
	            left=QIPL[indx2[k1-1]+i]+EXP((k1-i)*BETA3);
	          } else 
	            left=1;
	          typek=BP_pair[SEQ[k1]][SEQ[j]];
	          G=BETA2; //The penalty for the maximum arc in pseudoknot
 	          QGup[index]+=left*QGh[indx4[0][k1]+indx4[1][r]+indx4[2][s]+indx4[3][j]]*EXP(G);
	        }
//calculate the Gup matrix

	        QGv[index]=0;
	        for (k1=i+1;k1<r-1;k1++) 
	          for (k2=s+1;k2<j-1;k2++) {
	            left=QGu[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=QGuep[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            G=2*BETA2;
	            QGv[index]+=left*right*EXP(G);
	        }
// calculate the Gv matrix

	        QGvep[index]=0;
	        for (k1=i+1;k1<r-1;k1++) 
	          for (k2=s+1;k2<j-1;k2++) {
	            left=QGuep[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=QGuep[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            G=2*BETA2;
	            QGvep[index]+=left*right*EXP(G);
	        }
// calculate the Gvep matrix

	        QGvm[index]=0;
	        for (k1=i+1;k1<r-1;k1++) 
	          for (k2=s+1;k2<j-1;k2++) {
	            left=QGumm[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=QGump[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            G=2*BETA2;
	            QGvm[index]+=left*right*EXP(G);
	        }
// calculate the Gvm matrix

	        QGvmp[index]=0;
	        for (k1=i+1;k1<r-1;k1++) 
	          for (k2=s+1;k2<j-1;k2++) {
	            left=QGump[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=QGump[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            G=2*BETA2;
	            QGvmp[index]+=left*right*EXP(G);
	        }
// calculate the Gvmp matrix

	        QGvp[index]=0;
	        for (k1=i+1;k1<r-1;k1++) 
	          for (k2=s+1;k2<j-1;k2++) {
	            left=QGup[indx4[0][i]+indx4[1][k1]+indx4[2][s]+indx4[3][k2]];
	            right=QGup[indx4[0][k1+1]+indx4[1][r]+indx4[2][k2+1]+indx4[3][j]];
	            G=2*BETA2;
	            QGvp[index]+=left*right*EXP(G);
	        }
// calculate the Gvm matrix

	        QGw[index]=0;
	        for (k1=s+1;k1<j;k1++)
	          for (k2=k1+3;k2<j;k2++) {
	            left=QGuep[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=QGuep[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETAB+2*BETA2; // the penalty for the crossing matrix (shadow X)
	            QGw[index]+=left*right*EXP(G);
// U1 U'1 U2 U'2

	            left=QGvep[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=QGuep[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETAD+BETA2; // the penalty for the crossing matrix (shadow X)
	            QGw[index]+=left*right*EXP(G);
// V1 U1 V2 U'2
	        }
// calculate the Gw matrix

 	        QGwm[index]=0;
	        for (k1=s+1;k1<j;k1++)
	          for (k2=k1+3;k2<j;k2++) {
	            left=QGump[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=QGump[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETABm+2*BETA2; // the penalty for the crossing matrix (shadow X)
	            QGwm[index]+=left*right*EXP(G);
// U1 U'1 U2 U'2

	            left=QGvmp[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=QGump[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETADm+BETA2; // the penalty for the crossing matrix (shadow X)
	            QGwm[index]+=left*right*EXP(G);
// V1 U1 V2 U'2
	        }
// calculate the Gw matrix

	        QGwp[index]=0;
	        for (k1=s+1;k1<j;k1++)
	          for (k2=k1+3;k2<j;k2++) {
	            left=QGup[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=QGup[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETABp+2*BETA2; // the penalty for the crossing matrix (shadow X)
	            QGwp[index]+=left*right*EXP(G);
// U1 U'1 U2 U'2

	            left=QGvp[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][k2-1]];
	            right=QGup[indx4[0][s]+indx4[1][k1]+indx4[2][k2]+indx4[3][j]];
	            G=BETADp+BETA2; // the penalty for the crossing matrix (shadow X)
	            QGwp[index]+=left*right*EXP(G);
// V1 U1 V2 U'2
	        }
// calculate the Gwp matrix

	    }//4-dimension for-loop
// 4-dimensin matrix done!

	    QIbc[indx2[j]+i]=0;
	    type = BP_pair[SEQ[i]][SEQ[j]];
	    if (type) {
	      for (k1=i+1;k1<j;k1++)
	        for (k2=k1+TURN;k2<j;k2++) {
	          type2=BP_pair[SEQ[k2]][SEQ[k1]];
	          if (type2) {
	            G=LoopEnergy(k1-i-1,j-k2-1,type,type2,
		SEQ[i+1],SEQ[j-1],SEQ[k1-1],SEQ[k2+1]);
	            QIbc[indx2[j]+i]+=EXP(G)*QIbc[indx2[k2]+k1];

	            G=PARS->MLclosing+PARS->MLintern[type]+PARS->MLintern[type2];
	            if (k1-i>1) {
	              left=QIML[indx2[k1-1]+i+1];
	              right=EXP((j-k2-1)*PARS->MLbase);
	              QIbc[indx2[j]+i]+=left*right*EXP(G)*QIbc[indx2[k2]+k1];
	            }

	            if (j-k2>1) {
	              left=EXP((k1-i-1)*PARS->MLbase);
	              right=QfML[indx2[j-1]+k2+1]; 
	              QIbc[indx2[j]+i]+=left*right*EXP(G)*QIbc[indx2[k2]+k1];
	            }

	            if (k1-i>1 && j-k2>1) {
	              left=QIML[indx2[k1-1]+i+1]; 
	              right=QfML[indx2[j-1]+k2+1]; 
	              QIbc[indx2[j]+i]+=left*right*EXP(G)*QIbc[indx2[k2]+k1];
	            }
	          }
	      }
	      G=PARS->MLclosing+PARS->MLintern[type];
	      for (k1=i+TURN;k1<=j-1;k1++) {
	        left=QIML1[indx2[k1]+i+1];
	        if (k1<j-1) {
	          right=QfML[indx2[j-1]+k1+1]+EXP((j-1-k1)*PARS->MLbase);
	        } else 
	          right=1;
 	        QIbc[indx2[j]+i]+=EXP(G)*left*right; 
	      }
	    } else 
	      QIbc[indx2[j]+i]=0;
// calculate Ibc, which is a 1-structure where (i,j) is a base pair 

	    QI51[indx2[j]+i]=0;
	    for (r=i+1;r<j;r++) {
	      G=BETAA;
	      QI51[indx2[j]+i]+=QGv[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]*EXP(G);
	    }
//I'->V1 V2 (shadow A)
	    for (r=i+1;r<j;r++)
	      for (k1=r+2;k1<j;k1++)
	        for (s=k1+2;s+1<j;s++) {
	          left=QGv[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=QGuep[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          G=BETAC;
	          QI51[indx2[j]+i]+=left*right*EXP(G);
//I'->V1 U1 V2 U2 (shadow C)

	          left=QGu[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=QGw[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          QI51[indx2[j]+i]+=left*right;
//I'->W1 U1 W2 U2 (shadow B & D)
	    }
// calculate I51

	    QI5[indx2[j]+i]=Qf5[indx2[j]+i]+QI51[indx2[j]+i];
	    for (r=i+1;r<j;r++) {
	      QI5[indx2[j]+i]+=QI51[indx2[r]+i]*Qf5[indx2[j]+r+1];
	    }
	    for (r=i;r<j;r++)
	      for (s=r+TURN;s<=j;s++) {
	        type2 = BP_pair[SEQ[r]][SEQ[s]];
	        if (type2>2) G=PARS->TerminalAU; else G=0;
	        if (type2) {
	          if (r-i>0) left=QI5[indx2[r-1]+i]; else left=1;
	          if (j-s>0) right=Qf5[indx2[j]+s+1]; else right=1;
	          QI5[indx2[j]+i]+=left*right*QIbc[indx2[s]+r]*EXP(G);
	        }
	    }
// calculate I5

	    QIML1[indx2[j]+i]=0;
	    for (r=i+1;r<j;r++) {
	      G=BETAAm;
	      QIML1[indx2[j]+i]+=QGvm[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]*EXP(G);
	    }
//I'->V1 V2 (shadow A)
	    for (r=i+1;r<j;r++)
	      for (k1=r+2;k1<j;k1++)
	        for (s=k1+2;s+1<j;s++) {
	          left=QGvm[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=QGump[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          G=BETACm;
	          QIML1[indx2[j]+i]+=left*right*EXP(G);
//I'->V1 U1 V2 U2 (shadow C)

	          left=QGumm[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=QGwm[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          QIML1[indx2[j]+i]+=left*right;
//I'->W1 U1 W2 U2 (shadow B & D)
	    }
// calculate IML1

	    QIML[indx2[j]+i]=QfML[indx2[j]+i]+QIML1[indx2[j]+i];
	    for (r=i+1;r<j;r++) {
	      right=QfML[indx2[j]+r+1]+EXP((j-r)*PARS->MLbase);
	      QIML[indx2[j]+i]+=QIML1[indx2[r]+i]*right;
	    }
	    for (r=i;r<j;r++)
	      for (s=r+TURN;s<=j;s++) {
	        type2=BP_pair[SEQ[r]][SEQ[s]];
	        if (type2==0) continue;
	        G=PARS->MLintern[type];
	        if (r>i) left=QIML[indx2[r-1]+i]+EXP((r-i)*PARS->MLbase); else left=1;
	        if (s<j) right=QfML[indx2[j]+s+1]+EXP((j-s)*PARS->MLbase); else right=1;
	        QIML[indx2[j]+i]+=EXP(G)*left*right*QIbc[indx2[s]+r]*EXP(G);
	    }
// calculate IML

	    QIPL1[indx2[j]+i]=0;
	    for (r=i+1;r<j;r++) {
	      G=BETAAp;
	      QIPL1[indx2[j]+i]+=QGvp[indx4[0][i]+indx4[1][r]+indx4[2][r+1]+indx4[3][j]]*EXP(G);
	    }
//I'->V1 V2 (shadow A)
	    for (r=i+1;r<j;r++)
	      for (k1=r+2;k1<j;k1++)
	        for (s=k1+2;s+1<j;s++) {
	          left=QGvp[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=QGup[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          G=BETACp;
	          QIPL1[indx2[j]+i]+=left*right*EXP(G);
//I'->V1 U1 V2 U2 (shadow C)

	          left=QGup[indx4[0][i]+indx4[1][r]+indx4[2][k1+1]+indx4[3][s]];
	          right=QGwp[indx4[0][r+1]+indx4[1][k1]+indx4[2][s+1]+indx4[3][j]];
	          QIPL1[indx2[j]+i]+=left*right;
//I'->W1 U1 W2 U2 (shadow B & D)
	    }
// calculate IPL1

	    QIPL[indx2[j]+i]=QfPL[indx2[j]+i]+QIPL1[indx2[j]+i];
	    for (r=i+1;r<j;r++) {
	      right=QfPL[indx2[j]+r+1]+EXP((j-r)*BETA3);
	      QIPL[indx2[j]+i]+=QIPL1[indx2[r]+i]*right;
	    }
	    for (r=i;r<j;r++)
	      for (s=r+TURN;s<=j;s++) {
	        type2=BP_pair[SEQ[r]][SEQ[s]];
	        if (type2==0) continue;
	        G=BETA2;
	        if (r>i) left=QIPL[indx2[r-1]+i]+EXP((r-i)*BETA3); else left=1;
	        if (s<j) right=QfPL[indx2[j]+s+1]+EXP((j-s)*BETA3); else right=1;
	        QIPL[indx2[j]+i]+=EXP(G)*left*right*QIbc[indx2[s]+r]; 
	    }
// calculate IPL

	}
}
void pfunc(char *seq)
{
	int length;
	length=strlen(seq);
	update_fold_params();
	if (PARS) free(PARS);
	PARS = scale_parameters();
	encode_seq(seq);

	pfunc_initial(length);

	partition(seq);
	fprintf(fo, "PF: %lf\n", QI5[indx2[length]+1]);

	pfunc_freevar();
} 

void deblock(block *T, char *seq)
{
	int k1,k2,k3,G=0;
	double left,right,V,Qpro,Qi;
	int type, type2, typek, flag1,flag2;
	double seed;
	long indexk;

	seed=1-drand48();
// 	seed=1;

	switch (T->type) {
	case 0: //decompose QIbc
	  Qpro=seed*QIbc[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=0;

	  type = BP_pair[SEQ[T->i]][SEQ[T->j]];
	  if (!type) {
	    fprintf(fo, "Error block Ibc! Please check!\n");
	    exit(0);
	  }

	  sample[T->i]=T->j;
	  sample[T->j]=T->i;
	
	  for (k1=T->i+1;k1<T->j;k1++)
	    for (k2=k1+TURN;k2<T->j;k2++) {
	      type2=BP_pair[SEQ[k2]][SEQ[k1]];
	      if (type2) {
	        G=LoopEnergy(k1-T->i-1,T->j-k2-1,type,type2,
		SEQ[T->i+1],SEQ[T->j-1],SEQ[k1-1],SEQ[k2+1]);
	        Qi+=EXP(G)*QIbc[indx2[k2]+k1];
	        if (Qi>=Qpro) {
	          pushblock(k1,k2,0,0,0,0);
	          score+=G;
	          return;
	        }

	        G=PARS->MLclosing+PARS->MLintern[type]+PARS->MLintern[type2];
	        if (k1-T->i>1) {
	          left=QIML[indx2[k1-1]+T->i+1];
	          right=EXP((T->j-k2-1)*PARS->MLbase);
	          Qi+=left*right*EXP(G)*QIbc[indx2[k2]+k1];
	          if (Qi>=Qpro) {
	            pushblock(k1,k2,0,0,0,0);
	            pushblock(T->i+1,k1-1,0,0,2,0);
	            score+=G+(T->j-k2-1)*PARS->MLbase;
	            return;
	          }
	        }

	        if (T->j-k2>1) {
	          left=EXP((k1-T->i-1)*PARS->MLbase);
	          right=QfML[indx2[T->j-1]+k2+1]; 
	          Qi+=left*right*EXP(G)*QIbc[indx2[k2]+k1];
	          if (Qi>=Qpro) {
	            pushblock(k1,k2,0,0,0,0);
	            pushblock(k2+1,T->j-1,0,0,12,0);
	            score+=G+(k1-T->i-1)*PARS->MLbase;
	            return;
	          }
	        }

	        if (k1-T->i>1 && T->j-k2>1) {
	          left=QIML[indx2[k1-1]+T->i+1]; 
	          right=QfML[indx2[T->j-1]+k2+1]; 
	          Qi+=left*right*EXP(G)*QIbc[indx2[k2]+k1];
	          if (Qi>=Qpro) {
                    pushblock(k1,k2,0,0,0,0);
	            pushblock(T->i+1,k1-1,0,0,2,0);
	            pushblock(k2+1,T->j-1,0,0,12,0);
	            score+=G;
	            return;
	          }
	        }
	      }
	  }
	  G=PARS->MLclosing+PARS->MLintern[type];
	  for (k1=T->i+TURN;k1<=T->j-1;k1++) {
	    left=QIML1[indx2[k1]+T->i+1];
	    if (k1<T->j-1) {
	      right=QfML[indx2[T->j-1]+k1+1];
	      Qi+=EXP(G)*left*right;
	      if (Qi>=Qpro) {
	        pushblock(T->i+1,k1,0,0,5,0);
	        pushblock(k1+1,T->j-1,0,0,12,0);
	        score+=G;
	        return;
	      }

	      right=EXP((T->j-1-k1)*PARS->MLbase);
	      Qi+=EXP(G)*left*right;
	      if (Qi>=Qpro) {
	        pushblock(T->i+1,k1,0,0,5,0);
	        score+=G+(T->j-1-k1)*PARS->MLbase;
	        return;
	      }
	    } else {
	      Qi+=EXP(G)*left;
	      if (Qi>=Qpro) {
	        pushblock(T->i+1,k1,0,0,5,0);
	        score+=G;
	        return;
	      }
	    }
	  }
	  break;
	case 1: //decompose QI5
	  Qpro=seed*QI5[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=Qf5[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,11,0);
	    return;
	  }
// printf("Qi: %lf\n", Qi);
	  Qi+=QI51[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,4,0);
	    return;
	  }
// printf("Qi: %lf\n", Qi);
	  for (k1=T->i+1;k1<T->j;k1++) {
	    Qi+=QI51[indx2[k1]+T->i]*Qf5[indx2[T->j]+k1+1];
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,4,0);
	      pushblock(k1+1,T->j,0,0,11,0);
	      return;
	    }
// printf("Qi: %lf\n", Qi);
	  }
	  for (k1=T->i;k1<T->j;k1++)
	    for (k2=k1+TURN;k2<=T->j;k2++) {
	      type2 = BP_pair[SEQ[k1]][SEQ[k2]];
	      if (type2>2) G=PARS->TerminalAU; else G=0;
	      if (type2) {
	        if (k1-T->i>0) left=QI5[indx2[k1-1]+T->i]; else left=1;
	        if (T->j-k2>0) right=Qf5[indx2[T->j]+k2+1]; else right=1;
	        Qi+=left*right*QIbc[indx2[k2]+k1]*EXP(G);
	        if (Qi>=Qpro) {
	          if (k1-T->i>0) pushblock(T->i,k1-1,0,0,1,0);
	          if (T->j-k2>0) pushblock(k2+1,T->j,0,0,11,0);
	          pushblock(k1,k2,0,0,0,0);
	          score+=G;
	          return;
	        }
// printf("Qi: %lf\n", Qi);
	      }
	  }
	  break;
	case 2: //decompose QIML
	  Qpro=seed*QIML[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=QfML[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,12,0);
	    return;
	  }
	  Qi+=QIML1[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,5,0);
	    return;
	  }
	  for (k1=T->i+1;k1<T->j;k1++) {
	    right=QfML[indx2[T->j]+k1+1];
	    Qi+=QIML1[indx2[k1]+T->i]*right;
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,5,0);
	      pushblock(k1+1,T->j,0,0,12,0);
	      return;
	    }
	    Qi+=QIML1[indx2[k1]+T->i]*EXP((T->j-k1)*PARS->MLbase);
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,5,0);
	       score+=(T->j-k1)*PARS->MLbase;
	      return;
	    }
	  }
	  for (k1=T->i;k1<T->j;k1++)
	    for (k2=k1+TURN;k2<=T->j;k2++) {
	      type=BP_pair[SEQ[k1]][SEQ[k2]];
	      if (type==0) continue;
	      G=PARS->MLintern[type];
	      if (k1>T->i) left=QIML[indx2[k1-1]+T->i]; else left=1;
	      if (k2<T->j) right=QfML[indx2[T->j]+k2+1]; else right=1;
	      if (k1>T->i && k2<T->j) { 
	        Qi+=EXP(G)*left*right*QIbc[indx2[k2]+k1];
	        if (Qi>=Qpro) {
	          pushblock(k1,k2,0,0,0,0);
	          pushblock(T->i,k1-1,0,0,2,0);
	          pushblock(k2+1,T->j,0,0,12,0);
	          score+=G;
	          return;
	        }
	      }
	      if (k1>T->i) {
	        Qi+=EXP(G)*left*EXP((T->j-k2)*PARS->MLbase)*QIbc[indx2[k2]+k1];
	        if (Qi>=Qpro) {
	          pushblock(k1,k2,0,0,0,0);
	          pushblock(T->i,k1-1,0,0,2,0);
	          score+=G+(T->j-k2)*PARS->MLbase;
	          return;
	        }
	      }
	      if (k2<T->j) {
	        Qi+=EXP(G)*right*EXP((k1-T->i)*PARS->MLbase)*QIbc[indx2[k2]+k1];
	        if (Qi>=Qpro) {
	          pushblock(k1,k2,0,0,0,0);
	          pushblock(k2+1,T->j,0,0,12,0);
	          score+=G+(k1-T->i)*PARS->MLbase;
	          return;
	        }
	      }
	      Qi+=EXP(G)*EXP((k1-T->i)*PARS->MLbase)*EXP((T->j-k2)*PARS->MLbase)
		*QIbc[indx2[k2]+k1];
	      if (Qi>=Qpro) {
	        pushblock(k1,k2,0,0,0,0);
	        score+=(k1-T->i)*PARS->MLbase+G+(T->j-k2)*PARS->MLbase;
	        return;
	      }
	  }
	  break;
	case 3: //decompose QIPL
	  Qpro=seed*QIPL[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=QfPL[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,13,0);
	    return;
	  }
	  Qi+=QIPL1[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,6,0);
	    return;
	  }
	  for (k1=T->i+1;k1<T->j;k1++) {
	    Qi+=QIPL1[indx2[k1]+T->i]*QfPL[indx2[T->j]+k1+1];
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,6,0);
	      pushblock(k1+1,T->j,0,0,13,0);
	      return;
	    }
	    Qi+=QIPL1[indx2[k1]+T->i]*EXP((T->j-k1)*BETA3);
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,6,0);
	      score+=(T->j-k1)*BETA3;
	      return;
	    }
	  }
	  for (k1=T->i;k1<T->j;k1++)
	    for (k2=k1+TURN;k2<=T->j;k2++) {
	      type=BP_pair[SEQ[k1]][SEQ[k2]];
	      if (type==0) continue;
	      G=BETA2;
	      if (k1>T->i) left=QIPL[indx2[k1-1]+T->i]; else left=1;
	      if (k2<T->j) right=QfPL[indx2[T->j]+k2+1]; else right=1;
	      if (k1>T->i && k2<T->j) { 
	        Qi+=EXP(G)*left*right*QIbc[indx2[k2]+k1];
	        if (Qi>=Qpro) {
	          pushblock(k1,k2,0,0,0,0);
	          pushblock(T->i,k1-1,0,0,3,0);
	          pushblock(k2+1,T->j,0,0,13,0);
	          score+=G;
	          return;
	        }
	      }
	      if (k1>T->i) {
	        Qi+=EXP(G)*left*EXP((T->j-k2)*BETA3)*QIbc[indx2[k2]+k1];
	        if (Qi>=Qpro) {
	          pushblock(k1,k2,0,0,0,0);
	          pushblock(T->i,k1-1,0,0,3,0);
	          score+=G+(T->j-k2)*BETA3;
	          return;
	        }
	      }
	      if (k2<T->j) {
	        Qi+=EXP(G)*right*EXP((k1-T->i)*BETA3)*QIbc[indx2[k2]+k1];
	        if (Qi>=Qpro) {
	          pushblock(k1,k2,0,0,0,0);
	          pushblock(k2+1,T->j,0,0,13,0);
	          score+=G+(k1-T->i)*BETA3;
	          return;
	        }
	      }
	      Qi+=EXP(G)*EXP((k1-T->i)*BETA3)*EXP((T->j-k2)*BETA3)
		*QIbc[indx2[k2]+k1];
	      if (Qi>=Qpro) {
	        pushblock(k1,k2,0,0,0,0);
	        score+=G+(k1-T->i)*BETA3+(T->j-k2)*BETA3;
	        return;
	      }
	  }
	  break;
	case 4: //decompose QI51
	  Qpro=seed*QI51[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i+1;k1<T->j;k1++) {
	    G=BETAA;
	    Qi+=QGv[indx4[0][T->i]+indx4[1][k1]+indx4[2][k1+1]+indx4[3][T->j]]*EXP(G);
	    if (Qi>=Qpro) {
	      pushblock(T->i,T->j,k1,k1+1,26,0);
	      score+=G;
	      genusA++;
	      return;
	    }
	  }
	  for (k1=T->i+1;k1<T->j;k1++)
	    for (k2=k1+2;k2<T->j;k2++)
	      for (k3=k2+2;k3+1<T->j;k3++) {
	        left=QGv[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=QGuep[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        G=BETAC;
	        Qi+=left*right*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(T->i,k3,k1,k2+1,26,0);
	          pushblock(k1+1,T->j,k2,k3+1,22,0);
	          score+=G;
	          genusC++;
	          return;
	        }
	        left=QGu[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=QGw[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        Qi+=left*right;
	        if (Qi>=Qpro) {
	          pushblock(T->i,k3,k1,k2+1,21,0);
	          pushblock(k1+1,T->j,k2,k3+1,29,0);
	          return;
	        }
	  }
	  break;
	case 5: //decompose QIML1
	  Qpro=seed*QIML1[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i+1;k1<T->j;k1++) {
	    G=BETAAm;
	    Qi+=QGvm[indx4[0][T->i]+indx4[1][k1]+indx4[2][k1+1]+indx4[3][T->j]]*EXP(G);
	    if (Qi>=Qpro) {
	      pushblock(T->i,T->j,k1,k1+1,27,0);
	      score+=G;
	      genusA++;
	      return;
	    }
	  }
	  for (k1=T->i+1;k1<T->j;k1++)
	    for (k2=k1+2;k2<T->j;k2++)
	      for (k3=k2+2;k3+1<T->j;k3++) {
	        left=QGvm[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=QGump[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        G=BETACm;
	        Qi+=left*right*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(T->i,k3,k1,k2+1,27,0);
	          pushblock(k1+1,T->j,k2,k3+1,24,0);
	          score+=G;
	          genusC++;
	          return;
	        }
	        left=QGumm[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=QGwm[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        Qi+=left*right;
	        if (Qi>=Qpro) {
	          pushblock(T->i,k3,k1,k2+1,23,0);
	          pushblock(k1+1,T->j,k2,k3+1,30,0);
	          return;
	        }
	  }
	  break;
	case 6: //decompose QIPL1
	  Qpro=seed*QIPL1[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i+1;k1<T->j;k1++) {
	    G=BETAAp;
	    Qi+=QGvp[indx4[0][T->i]+indx4[1][k1]+indx4[2][k1+1]+indx4[3][T->j]]*EXP(G);
	    if (Qi>=Qpro) {
	      pushblock(T->i,T->j,k1,k1+1,28,0);
	      score+=G;
	      genusA++;
	      return;
	    }
	  }
	  for (k1=T->i+1;k1<T->j;k1++)
	    for (k2=k1+2;k2<T->j;k2++)
	      for (k3=k2+2;k3+1<T->j;k3++) {
	        left=QGvp[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=QGup[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        G=BETACp;
	        Qi+=left*right*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(T->i,k3,k1,k2+1,28,0);
	          pushblock(k1+1,T->j,k2,k3+1,25,0);
	          score+=G;
	          genusC++;
	          return;
	        }
	        left=QGup[indx4[0][T->i]+indx4[1][k1]+indx4[2][k2+1]+indx4[3][k3]];
	        right=QGwp[indx4[0][k1+1]+indx4[1][k2]+indx4[2][k3+1]+indx4[3][T->j]];
	        Qi+=left*right;
	        if (Qi>=Qpro) {
	          pushblock(T->i,k3,k1,k2+1,25,0);
	          pushblock(k1+1,T->j,k2,k3+1,31,0);
	          return;
	        }
	  }
	  break;
	case 10: //decompose Qc
	  Qpro=seed*Qc[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  type=BP_pair[SEQ[T->i]][SEQ[T->j]];
	  if (!type) {
	    fprintf(fo, "Error block Ibc! Please check!\n");
	    exit(0);
	  }
	  sample[T->i]=T->j;
	  sample[T->j]=T->i;
	  G=HairpinE(T->j-T->i-1, type, SEQ[T->i+1], SEQ[T->j-1], seq+T->i-1);
	  Qi=EXP(G);
	  if (Qi>=Qpro) {
	    score+=G;
	    return;
	  }
	  for (k1=T->i+1;k1<T->j;k1++)
	    for (k2=k1+TURN; k2<T->j;k2++) {
	      type2=BP_pair[SEQ[k2]][SEQ[k1]];
	      if (type2==0) continue;
	      G=LoopEnergy(k1-T->i-1, T->j-k2-1, type, type2,
		SEQ[T->i+1], SEQ[T->j-1], SEQ[k1-1], SEQ[k2+1]);
	      Qi+=EXP(G)*Qc[indx2[k2]+k1];
	      if (Qi>=Qpro) {
	        pushblock(k1,k2,0,0,10,0);
	        score+=G;
	        return;
	      }
	  }
	  for (k1=T->i+1;k1<T->j-1;k1++) {
	    G=PARS->MLintern[type]+PARS->MLclosing;
	    Qi+=QfML1[indx2[k1]+T->i+1]*QfML[indx2[T->j-1]+k1+1]*EXP(G);
	    if (Qi>=Qpro) {
	      pushblock(T->i+1,k1,0,0,15,0);
	      pushblock(k1+1,T->j-1,0,0,12,0);
	      score+=G;
	      return;
	    } 
	  }
	  break;
	case 11: //decompose Qf5
	  if (Qf5[indx2[T->j]+T->i]==1) return;
	  Qpro=seed*Qf5[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=Qf51[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,14,0);
	    return;
	  }
	  for (k1=T->i;k1<T->j;k1++) {
	    Qi+=Qf51[indx2[k1]+T->i]*Qf5[indx2[T->j]+k1+1];
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,14,0);
	      pushblock(k1+1,T->j,0,0,11,0);
	      return;
	    }
	  }
	  Qi+=1;
	  if (Qi>=Qpro) return;
	  break;
	case 12: //decompose QfML
	  Qpro=seed*QfML[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=QfML1[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,15,0);
	    return;
	  }
	  for (k1=T->i;k1<T->j;k1++) {
	    Qi+=QfML1[indx2[k1]+T->i]*QfML[indx2[T->j]+k1+1];
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,15,0);
	      pushblock(k1+1,T->j,0,0,12,0);
	      return;
	    }
	    Qi+=QfML1[indx2[k1]+T->i]*EXP((T->j-k1)*PARS->MLbase);
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,15,0);
	      score+=(T->j-k1)*PARS->MLbase;
	      return;
	    }
	  }
	  break;
	case 13: //decompose QfPL
	  Qpro=seed*QfPL[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=QfPL1[indx2[T->j]+T->i];
	  if (Qi>=Qpro) {
	    pushblock(T->i,T->j,0,0,16,0);
	    return;
	  }
	  for (k1=T->i;k1<T->j;k1++) {
	    Qi+=QfPL1[indx2[k1]+T->i]*QfPL[indx2[T->j]+k1+1];
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,16,0);
	      pushblock(k1+1,T->j,0,0,13,0);
	      return;
	    }
	    Qi+=QfPL1[indx2[k1]+T->i]*EXP((T->j-k1)*BETA3);
	    if (Qi>=Qpro) {
	      pushblock(T->i,k1,0,0,16,0);
	      score+=(T->j-k1)*BETA3;
	      return;
	    }
	  }
	  break;
	case 14: //decompose Qf51
	  Qpro=seed*Qf51[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i;k1<T->j;k1++) {
	    type = BP_pair[SEQ[k1]][SEQ[T->j]];
	    if (type==0) continue;
	    if (type>2) G=PARS->TerminalAU; else G=0;
	    Qi+=Qc[indx2[T->j]+k1]*EXP(G);
	    if (Qi>=Qpro) {
	      pushblock(k1,T->j,0,0,10,0);
	      score+=G;
	      return;
	    }
	  }
	  break;
	case 15: //decompose QfML1
	  Qpro=seed*QfML1[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i;k1<T->j;k1++) {
	    type = BP_pair[SEQ[k1]][SEQ[T->j]];
	    if (type==0) continue;
	    G=PARS->MLintern[type]+(k1-T->i)*PARS->MLbase;
	    Qi+=Qc[indx2[T->j]+k1]*EXP(G);
	    if (Qi>=Qpro) {
	      pushblock(k1,T->j,0,0,10,0);
	      score+=G;
	      return;
	    }
	  }
	  break;
	case 16: //decompose QfPL1
	  Qpro=seed*QfPL1[indx2[T->j]+T->i];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i;k1<T->j;k1++) {
	    type = BP_pair[SEQ[k1]][SEQ[T->j]];
	    if (type==0) continue;
	    G=BETA2+(k1-T->i)*BETA3;
	    Qi+=Qc[indx2[T->j]+k1]*EXP(G);
	    if (Qi>=Qpro) {
	      pushblock(k1,T->j,0,0,10,0);
	      score+=G;
	      return;
	    }
	  }
	  break;
	case 20: //decompose QGtight
	  Qpro=seed*QGtight[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  type=BP_pair[SEQ[T->i]][SEQ[T->j]];
	  type2=BP_pair[SEQ[T->s]][SEQ[T->r]];
	  if (type==0 || type2==0) {
	    fprintf(fo, "Error block Ibc! Please check!\n");
	    exit(0);
	  }
	  sample[T->i]=T->j;
	  sample[T->j]=T->i;
	  G=(int)(SIGMA*LoopEnergy(T->r-T->i-1,T->j-T->s-1,type,type2,
		  SEQ[T->i+1],SEQ[T->j-1],SEQ[T->r-1],SEQ[T->s+1]));
	  Qi=EXP(G);
	  if (Qi>=Qpro) {
	    sample[T->r]=T->s;
	    sample[T->s]=T->r;
	    score+=G;
	    return;
	  }
	  for (k1=T->i+1;k1<T->r;k1++)
	    for (k2=T->s+1;k2<T->j;k2++) {
	      typek=BP_pair[SEQ[k2]][SEQ[k1]];
	      indexk=indx4[0][k1]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][k2];
	      if (typek) {
	        G=(int)(SIGMA*LoopEnergy(k1-T->i-1,T->j-k2-1,type,typek,
		SEQ[T->i+1],SEQ[T->j-1],SEQ[k1-1],SEQ[k2+1]));

	        Qi+=QGtight[indexk]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,k2,T->r,T->s,20,0);
	          score+=G;
	          return;
	        }
	        G=PARS->MLclosing+PARS->MLintern[type];
	        if (k1-T->i>1 ) {
	          Qi+=QGtight[indexk]*EXP(G)*QIML[indx2[k1-1]+T->i+1];
	          if (Qi>=Qpro) {
	             pushblock(k1,k2,T->r,T->s,20,0);
	             pushblock(T->i+1,k1-1,0,0,2,0);
	             score+=G;
	             return;
	          }
	        } 
	        if (T->j-k2>1) {
	          Qi+=QGtight[indexk]*EXP(G)*QIML[indx2[T->j-1]+k2+1];
	          if (Qi>=Qpro) {
	            pushblock(k1,k2,T->r,T->s,20,0);
	            pushblock(k2+1,T->j-1,0,0,2,0);
	            score+=G;
	            return;
	          }
	        }
	        if (k1-T->i>1 && T->j-k2>1) {
	          Qi+=QGtight[indexk]*EXP(G)*QIML[indx2[k1-1]+T->i+1]
		*QIML[indx2[T->j-1]+k2+1];
	          if (Qi>=Qpro) {
	            pushblock(k1,k2,T->r,T->s,20,0);
	            pushblock(T->i+1,k1-1,0,0,2,0);
	            pushblock(k2+1,T->j-1,0,0,2,0);
	            score+=G;
	            return;
	          }
	        }
	      }
	  }
	case 21: //decompose QGu
	  Qpro=seed*QGu[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      typek=BP_pair[SEQ[k1]][SEQ[T->j]];
	      if (typek>2) G=PARS->TerminalAU; else G=0;
	      if (k1-T->i>0) left=QI5[indx2[k1-1]+T->i]; else left=1;
	      if (k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*left*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,1,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G;
	          return;
	        }
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*left*EXP((k2-T->s)*BETA3)*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,1,0);
	          score+=G+(k2-T->s)*BETA3;
	          return;
	        }
	      } else {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*left*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,1,0);
	          score+=G;
	          return;
	        }
	      }
	  }
	  break;
	case 22: //decompose QGuep
	  Qpro=seed*QGuep[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      typek=BP_pair[SEQ[k1]][SEQ[T->j]];
	      if (typek>2) G=PARS->TerminalAU; else G=0;
	      if (k1-T->i>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*QIPL[indx2[k1-1]+T->i]*EXP((k2-T->s)*BETA3)*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,3,0);
	          score+=G+(k2-T->s)*BETA3;
	          return;
	        }
	      }
	      if (k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*EXP((k1-T->i)*BETA3)*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G+(k1-T->i)*BETA3;
	          return;
	        }
	      }
	      if (k1-T->i>0 && k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*QIPL[indx2[k1-1]+T->i]*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,3,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G;
	          return;
	        }
	      }
	      Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*EXP((k1-T->i)*BETA3)*EXP((k2-T->s)*BETA3)*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(k1,T->j,T->r,k2,20,0);
	        score+=G+(k1-T->i)*BETA3+(k2-T->s)*BETA3;
	        return;
	      }
	  }
	case 23: //decompose QGumm
	  Qpro=seed*QGumm[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      typek=BP_pair[SEQ[k1]][SEQ[T->j]];
	      G=PARS->MLintern[typek];
	      if (k1-T->i>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*QIML[indx2[k1-1]+T->i]*EXP((k2-T->s)*BETA3)*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,2,0);
	          score+=G+(k2-T->s)*BETA3;
	          return;
	        }
	      }
	      if (k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*EXP((k1-T->i)*PARS->MLbase)*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G+(k1-T->i)*PARS->MLbase;
	          return;
	        }
	      }
	      if (k1-T->i>0 && k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*QIML[indx2[k1-1]+T->i]*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,2,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G;
	          return;
	        }
	      }
	      Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*EXP((k1-T->i)*PARS->MLbase)*EXP((k2-T->s)*BETA3)*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(k1,T->j,T->r,k2,20,0);
	        score+=G+(k1-T->i)*PARS->MLbase+(k2-T->s)*BETA3;
	        return;
	      }
	  }
	  break;
	case 24: //decompose QGump
	  Qpro=seed*QGump[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      typek=BP_pair[SEQ[k1]][SEQ[T->j]];
	      G=PARS->MLintern[typek];
	      if (k1-T->i>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*QIPL[indx2[k1-1]+T->i]*EXP((k2-T->s)*BETA3)*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,3,0);
	          score+=G+(k2-T->s)*BETA3;
	          return;
	        }
	      }
	      if (k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*EXP((k1-T->i)*BETA3)*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G+(k1-T->i)*BETA3;
	          return;
	        }
	      }
	      if (k1-T->i>0 && k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*QIPL[indx2[k1-1]+T->i]*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,3,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G;
	          return;
	        }
	      }
	      Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*EXP((k1-T->i)*BETA3)*EXP((k2-T->s)*BETA3)*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(k1,T->j,T->r,k2,20,0);
	        score+=G+(k1-T->i)*BETA3+(k2-T->s)*BETA3;
	        return;
	      }
	  }
	  break;
	case 25: //decompose QGup
	  Qpro=seed*QGup[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i;k1<T->r;k1++)
	    for (k2=T->s>T->r+TURN?T->s:T->r+TURN; k2<T->j; k2++) {
	      if (BP_pair[SEQ[k1]][SEQ[T->j]]==0 || BP_pair[SEQ[T->r]][SEQ[k2]]==0) continue;
	      typek=BP_pair[SEQ[k1]][SEQ[T->j]];;
	      G=BETA2;
	      if (k1-T->i>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*QIPL[indx2[k1-1]+T->i]*EXP((k2-T->s)*BETA3)*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,3,0);
	          score+=G+(k2-T->s)*BETA3;
	          return;
	        }
	      }
	      if (k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*EXP((k1-T->i)*BETA3)*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G+(k1-T->i)*BETA3;
	          return;
	        }
	      }
	      if (k1-T->i>0 && k2-T->s>0) {
	        Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*QIPL[indx2[k1-1]+T->i]*QIPL[indx2[k2-1]+T->s]*EXP(G);
	        if (Qi>=Qpro) {
	          pushblock(k1,T->j,T->r,k2,20,0);
	          pushblock(T->i,k1-1,0,0,3,0);
	          pushblock(T->s,k2-1,0,0,3,0);
	          score+=G;
	          return;
	        }
	      }
	      Qi+=QGtight[indx4[0][k1]+indx4[1][T->r]+indx4[2][k2]+indx4[3][T->j]]
		*EXP((k1-T->i)*BETA3)*EXP((k2-T->s)*BETA3)*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(k1,T->j,T->r,k2,20,0);
	        score+=G+(k1-T->i)*BETA3+(k2-T->s)*BETA3;
	        return;
	      }
	  }
	  break;
	case 26: //decompose QGv
	  Qpro=seed*QGv[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=QGu[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=QGuep[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      G=2*BETA2;
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2,k1,T->s,21,0);
	        pushblock(k1+1,T->j,T->r,k2+1,22,0);
	        score+=G;
	        return;
	      }
	  }
	  break;
	case 32: //decompose QGvep
	  Qpro=seed*QGvep[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=QGuep[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=QGuep[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      G=2*BETA2;
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2,k1,T->s,22,0);
	        pushblock(k1+1,T->j,T->r,k2+1,22,0);
	        score+=G;
	        return;
	      }
	  }
	  break;
	case 27: //decompose QGvm
	  Qpro=seed*QGvm[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=QGumm[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=QGump[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      G=2*BETA2;
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2,k1,T->s,23,0);
	        pushblock(k1+1,T->j,T->r,k2+1,24,0);
	        score+=G;
	        return;
	      }
	  }
	  break;
	case 33: //decompose QGvmp
	  Qpro=seed*QGvmp[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=QGump[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=QGump[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      G=2*BETA2;
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2,k1,T->s,24,0);
	        pushblock(k1+1,T->j,T->r,k2+1,24,0);
	        score+=G;
	        return;
	      }
	  }
	  break;
	case 28: //decompose QGvp
	  Qpro=seed*QGvp[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->i+1;k1<T->r-1;k1++)
	    for (k2=T->s+1;k2<T->j-1;k2++) {
	      left=QGup[indx4[0][T->i]+indx4[1][k1]+indx4[2][T->s]+indx4[3][k2]];
	      right=QGup[indx4[0][k1+1]+indx4[1][T->r]+indx4[2][k2+1]+indx4[3][T->j]];
	      G=2*BETA2;
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2,k1,T->s,25,left);
	        pushblock(k1+1,T->j,T->r,k2+1,25,right);
	        score+=G;
	        return;
	      }
	  }
	  break;
	case 29: //decompose QGw
	  Qpro=seed*QGw[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->s+1;k1<T->j;k1++)
	    for (k2=k1+3;k2<T->j;k2++) {
	      left=QGuep[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=QGuep[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETAB+2*BETA2; // the penalty for the crossing matrix (shadow X)
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2-1,T->r,k1+1,22,0);
	        pushblock(T->s,T->j,k1,k2,22,0);
	        score+=G;
	        genusB++;
	        return;
	      }

	      left=QGvep[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=QGuep[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETAD+BETA2; // the penalty for the crossing matrix (shadow X)
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2-1,T->r,k1+1,32,0);
	        pushblock(T->s,T->j,k1,k2,22,0);
	        score+=G;
	        genusD++;
	        return;
	      }
	  }
	  break;
	case 30: //decompose QGwm
	  Qpro=seed*QGwm[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->s+1;k1<T->j;k1++)
	    for (k2=k1+3;k2<T->j;k2++) {
	      left=QGump[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=QGump[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETABm+2*BETA2; // the penalty for the crossing matrix (shadow X)
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2-1,T->r,k1+1,24,0);
	        pushblock(T->s,T->j,k1,k2,24,0);
	        score+=G;
	        genusB++;
	        return;
	      }

	      left=QGvmp[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=QGump[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETADm+BETA2; // the penalty for the crossing matrix (shadow X)
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2-1,T->r,k1+1,33,0);
	        pushblock(T->s,T->j,k1,k2,24,0);
	        score+=G;
	        genusD++;
	        return;
	      }
	  }
	  break;
	case 31: //decompose QGwp
	  Qpro=seed*QGwp[indx4[0][T->i]+indx4[1][T->r]+indx4[2][T->s]+indx4[3][T->j]];
	  if (Qpro==0) return;
	  Qi=0;
	  for (k1=T->s+1;k1<T->j;k1++)
	    for (k2=k1+3;k2<T->j;k2++) {
	      left=QGup[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=QGup[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETABp+2*BETA2; // the penalty for the crossing matrix (shadow X)
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2-1,T->r,k1+1,25,0);
	        pushblock(T->s,T->j,k1,k2,25,0);
	        score+=G;
	        genusB++;
	        return;
	      }

	      left=QGvp[indx4[0][T->i]+indx4[1][T->r]+indx4[2][k1+1]+indx4[3][k2-1]];
	      right=QGup[indx4[0][T->s]+indx4[1][k1]+indx4[2][k2]+indx4[3][T->j]];
	      G=BETADp+BETA2; // the penalty for the crossing matrix (shadow X)
	      Qi+=left*right*EXP(G);
	      if (Qi>=Qpro) {
	        pushblock(T->i,k2-1,T->r,k1+1,28,0);
	        pushblock(T->s,T->j,k1,k2,25,0);
	        score+=G;
	        genusD++;
	        return;
	      }
	  }
	  break;
	default: break;
	}// end switch 
	fprintf(fo, "Error! Do not match the energy, please check your program carefully!\n");
	fprintf(fo, "%d\n", T->type);
}

int exist(char *s) 
{
	ana *index;
	int res=0;
	index=stat;
	while (index!=NULL) {
	  if (strcmp(s,index->structure)==0) {
	    index->cnt++;
	    return 1;
	  }
	  index=index->Next;
	}
	return 0;
}

char *samplestructure (char *seq, int samples)
{
	int length,i,j,r,s, correct,total, max;
	int distribute[4][10];
	char *struc;
	block *T,*IN;
	int p5[300],p8[300],p9[300],p95[300],p99[300];
	ana *index,*mi;

	length=strlen(seq);

	pfunc_initial(length);
	
	if (!PARS) {
		update_fold_params();
		PARS = scale_parameters();
	}
	
	encode_seq(seq);

	for (i=0;i<10;i++)
	  for (j=0;j<4;j++)
	    distribute[j][i]=0;
	partition(seq);
	fprintf(fo, "QI5 %lf\n", QI5[indx2[length]+1]);

	stat=NULL;
	for (i=0;i<samples;i++) {
	  for (j=0;j<=length;j++) sample[j]=0;
	  sample[0]=length;
	  genusA=0;
	  genusB=0;
	  genusC=0;
	  genusD=0;
	  score=0;
	  pushblock(1,length,0,0,1,0);
	  while (Bstack!=NULL) {
// 	    IN=Bstack;
// 	    printf("------------------------------------------------\n");
// 	    while (IN!=NULL) {
// 	      printf("(%d %d %d %d) Type:%d\n", IN->i,IN->j,IN->r,IN->s,IN->type);
// 	      IN=IN->Next;
// 	    }
	    T=popblock(); 
	    deblock(T,seq);
	    free(T);
	  }
	  if (genusA+genusB+genusC+genusD==0) distribute[0][0]++;
	  if (genusA) distribute[0][genusA+genusB+genusC+genusD]++;
	  if (genusB) distribute[1][genusA+genusB+genusC+genusD]++;
	  if (genusC) distribute[2][genusA+genusB+genusC+genusD]++;
	  if (genusD) distribute[3][genusA+genusB+genusC+genusD]++;

	  for (j=1;j<=length;j++) {
// 	    printf("%d  ", sample[i]);
	    if (sample[j] && sample[j]>j) {
	      Basepair[indx2[sample[j]]+j]++;
	    }
	  }
	  struc=pair2structure(sample);
		if (strlen(struc) != length) {
			free(struc);		
			i--;
			continue; // simple roundabout around the bug that struc is sometimes corrupted...
		} else {
	 		fprintf(stdout, "%s %6.2f\n", struc, score/100.0);
		}
	  
	  
	  if (score < (int)(0.9*mfe) && !exist(struc) ) {
	    index=(ana *) space (sizeof(ana)*1);
	    index->structure=struc;
	    index->e=score;
	    index->cnt=1;
	    index->Next=stat;
	    stat=index;
	  } else {
			free(struc);
		}
	}
	for (i=0;i<4;i++) {
	  fprintf(fo, "Genus %d [A]: %d\n",i, distribute[0][i]);
	  fprintf(fo, "Genus %d [B]: %d\n",i, distribute[1][i]);
	  fprintf(fo, "Genus %d [C]: %d\n",i, distribute[2][i]);
	  fprintf(fo, "Genus %d [D]: %d\n",i, distribute[3][i]);
	}

	for (i=0;i<6;i++) {
	  index=stat; max=0; mi=NULL;
	  while (index!=NULL) {
	    if (index->e<max) {
	      mi=index; 
	      max=index->e;
	    }
	    index=index->Next;
	  }
	  if (max && mi) {
	    fprintf(fo, "%s, E:%d fre:%d\n", mi->structure, mi->e, mi->cnt);
	    mi->e=10;
	  }
	}

	p5[0]=length; p8[0]=length;
	p9[0]=length; p95[0]=length; p99[0]=length;
	for (i=1;i<=length;i++) {
	  p5[i]=0; p8[i]=0;
	  p9[i]=0; p95[i]=0; p99[i]=0;
	}
// 	for (i=1;i<length;i++)
// 	  for (j=i+TURN;j<=length;j++)
// 	    fprintf(fo, "(%d,%d): %d\n", i,j,Basepair[indx2[j]+i]);


	for (i=1;i<=length;i++) {
	  if (ptable[i] && i<ptable[i]) {
	    if (Basepair[indx2[ptable[i]]+i]>=50000) {
	      p5[ptable[i]]=i; p5[i]=ptable[i];
	    }
	    if (Basepair[indx2[ptable[i]]+i]>=80000) {
	      p8[ptable[i]]=i; p8[i]=ptable[i];
	    }
	    if (Basepair[indx2[ptable[i]]+i]>=90000) {
	      p9[ptable[i]]=i; p9[i]=ptable[i];
	    }
	    if (Basepair[indx2[ptable[i]]+i]>=95000) {
	      p95[ptable[i]]=i; p95[i]=ptable[i];
	    }
	    if (Basepair[indx2[ptable[i]]+i]>=99000) {
	      p99[ptable[i]]=i; p99[i]=ptable[i];
	    }
	  }
	}

	if (p_nat) {
		struc=pair2structure(p5);
		correct=compare(p5,p_nat);
		total=number(p5);
	 	fprintf(fo, "P50  %s  CORRECT: %d  PPV: %f\n", struc, correct,
			(double)correct/(double)total);
		free(struc);
		correct=compare(p8,p_nat);
		total=number(p8);
		struc=pair2structure(p8);
	 	fprintf(fo, "P80  %s  CORRECT: %d  PPV: %f\n", struc, correct,
			(double)correct/(double)total);
		free(struc);
		correct=compare(p9,p_nat);
		total=number(p9);
		struc=pair2structure(p9);
	 	fprintf(fo, "P90  %s  CORRECT: %d  PPV: %f\n", struc, correct,
			(double)correct/(double)total);
		free(struc);
		correct=compare(p95,p_nat);
		total=number(p95);
		struc=pair2structure(p95);
	 	fprintf(fo, "P95  %s  CORRECT: %d  PPV: %f\n", struc, correct,
			(double)correct/(double)total);
		free(struc);
		correct=compare(p99,p_nat);
		total=number(p99);
		struc=pair2structure(p99);
	 	fprintf(fo, "P99  %s  CORRECT: %d  PPV: %f\n", struc, correct,
			(double)correct/(double)total);
		free(struc);
		fprintf(fo, "\n");
	}

	index=stat;
	while (index!=NULL) {
	  index=index->Next;
	  free(stat->structure);
	  free(stat);
	  stat=index;
	}
	pfunc_freevar();

}




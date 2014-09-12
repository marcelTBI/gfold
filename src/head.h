#ifndef __HEAD_H
#define __HEAD_H

#include <string.h>
#include "PS_dot.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "utils.h"
#include "energy_par.h"
#include "params.h"
#include "fold.h"

#define MAXENG 10000

/*
#define BETAA 970
#define BETAAm 1470
#define BETAAp 1470
#define BETAB 1170
#define BETABm 1670
#define BETABp 1670
#define BETAC 1270
#define BETACm 1770
#define BETACp 1770
#define BETAD 1570
#define BETADm 2070
#define BETADp 2070
*/


#define BETAA b1
#define BETAAm b1m
#define BETAAp b1p
#define BETAB b2
#define BETABm b2m
#define BETABp b2p
#define BETAC b3
#define BETACm b3m
#define BETACp b3p
#define BETAD b4
#define BETADm b4m
#define BETADp b4p


#define SIGMA 1

#define BETA2 10
#define BETA3 10

#define R GASCONST

double betaScale;
extern double betaScale;

#define EXP(G) expl(-1*(10*(G))/(betaScale*R*(PARS->temperature+K0)))

int b1,b1m,b1p,b2,b2m,b2p,b3,b3m,b3p,b4,b4m,b4p;
extern int b1;
extern int b1m;
extern int b1p;
extern int b2;
extern int b2m;
extern int b2p;
extern int b3;
extern int b3m;
extern int b3p;
extern int b4;
extern int b4m;
extern int b4p;

typedef struct block block;

struct block
{
	int i;
	int j;
	int r;
	int s;
	int type;
	// 0:Ibc
	// 1:I5 2:IML 3:IPL
	// 4:I51 5:IML1 6:IPL1

	// 10:c
	// 11:f5 12:fML 13:fPL
	// 14:f51 15:fML1 16:fPL1

	// 20:Gtight
	// 21:Gu 22:Guep 23:Gumm 24:Gump 25:Gup
	// 26:Gv 27:Gvm 28:Gvp
	// 29:Gw 30:Gwm 31:Gwp
	// 32:Gvep 33:Gvmp
	int value;
	block *Next;
};


int *SEQ;
extern int *SEQ;

int mfe;
extern int mfe;

int genusA;
int genusB;
int genusC;
int genusD;
extern int genusA;
extern int genusB;
extern int genusC;
extern int genusD;

int *p_nat;
extern int *p_nat;

extern void symbFold(char *seq);
extern void pfunc(char *seq);
extern char *samplestructure (char *seq, int samples);


block *Bstack;
extern block *Bstack;

FILE *fo;
extern FILE *fo;

int *ptable;
extern int *ptable;


extern block *initialblock(void);

extern void pushblock(int i, int j, int r, int s, int type, int value);

extern block* popblock();

extern char *pair2structure(int *pair);
extern int *structure2pair(char *s);
extern int compare(int *a, int *b);
extern int number(int *a);


paramT *PARS;
long *indx4[4];
int *c,*fML,*fML1,*fPL1, *fPL,*indx2;
int *f5,*f51, *I51, *I5, *IML1, *IML, *IPL1, *IPL, *Ibc;
int *Gtight, *Gu, *Guep, *Gv, *Gw, *Gumm;
int *Gump, *Gvm, *Gwm, *Gup, *Gvp, *Gwp;
int *Gh, *Gvep, *Gvmp;

extern paramT *PARS;
extern long *indx4[4];
extern int *c,*fML,*fML1,*fPL1, *fPL,*indx2;
extern  int *f5,*f51, *I51, *I5, *IML1, *IML, *IPL1, *IPL, *Ibc;
extern int *Gtight, *Gu, *Guep, *Gv, *Gw, *Gumm;
extern int *Gump, *Gvm, *Gwm, *Gup, *Gvp, *Gwp;
extern int *Gh, *Gvep, *Gvmp;


typedef struct ana ana;
struct ana
{
	char *structure;
	int cnt;
	int e;
	ana *Next;
};
ana *stat;
extern ana *stat;

int *indx2;
double *Qc,*QfML,*QfML1,*QfPL1, *QfPL;
double *Qf5,*Qf51, *QI51, *QI5, *QIML1, *QIML, *QIPL1, *QIPL, *QIbc;
double *QGtight, *QGu, *QGuep, *QGv, *QGw, *QGumm;
double *QGump, *QGvm, *QGwm, *QGup, *QGvp, *QGwp;
double *QGh, *QGvep, *QGvmp;
int *Basepair;
int score;
int *sample;

extern int *indx2;
extern double *Qc,*QfML,*QfML1,*QfPL1, *QfPL;
extern double *Qf5,*Qf51, *QI51, *QI5, *QIML1, *QIML, *QIPL1, *QIPL, *QIbc;
extern double *QGtight, *QGu, *QGuep, *QGv, *QGw, *QGumm;
extern double *QGump, *QGvm, *QGwm, *QGup, *QGvp, *QGwp;
extern double *QGh, *QGvep, *QGvmp;
extern int *Basepair;
extern int score;
extern int *sample;

#endif

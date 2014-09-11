/*
	PostScript and GML output for RNA secondary structures
		    and pair probability matrices

		 c  Ivo Hofacker and Peter F Stadler
			  Vienna RNA package
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "fold_vars.h"
#include "PS_dot.h"


#define PUBLIC
#define  PRIVATE   static
#define  MAX(A,B)    (A)>(B)?(A):(B)
#ifndef PI
#define  PI       3.141592654
#endif
#define  PIHALF       PI/2.

PUBLIC int   gmlRNA(char *string, char *structure, char *ssfile, char option);
PUBLIC int   PS_rna_plot_a(char *string, char *structure, char *ssfile,
			   char *pre, char *post);
PUBLIC int   PS_rna_plot(char *string, char *structure, char *ssfile);
PUBLIC int   ssv_rna_plot(char *string, char *structure, char *ssfile);
PUBLIC int   xrna_plot(char *string, char *structure, char *ssfile);

PUBLIC int   simple_xy_coordinates(short *pair_table, double *X, double *Y);
extern int   naview_xy_coordinates(short *pair_table, double *X, double *Y);

PUBLIC int   rna_plot_type = 1;  /* 0 = simple, 1 = naview */

PUBLIC int   PS_dot_plot(char *string, char *wastlfile);
PUBLIC int   PS_color_dot_plot(char *seq, cpair *pi, char *wastlfile);
PUBLIC int   PS_dot_plot_list(char *string, char *wastlfile, struct plist *pl,
			      struct plist *mf, char *comment);
PUBLIC int   PS_dot_plot_turn(char *seq, struct plist *pl, char *wastlfile,
			      int winSize);

/* local functions */
PRIVATE void   loop(int i, int j, short *pair_table);
PRIVATE FILE  *PS_dot_common(char *seq, char *wastlfile, char *comment,
			     int winsize);

/* local variables for parsing routines */
PRIVATE double  *angle;
PRIVATE int    *loop_size, *stack_size;
PRIVATE int     lp, stk;

extern  int cut_point;   /* set to first pos of second seq for cofolding */

/*---------------------------------------------------------------------------*/

/* options for gml output:
   uppercase letters: print sequence labels
   lowercase letters: no sequence lables
   graphics information:
   x X  simple xy plot
   (nothing else implemented at present)
   default:           no graphics data at all
*/


/*---------------------------------------------------------------------------*/


static const char *RNAss_head =
"%%BeginProlog\n"
"/RNAplot 100 dict def\n"
"RNAplot begin\n"
"/fsize  14 def\n"
"/outlinecolor {0.2 setgray} bind def\n"
"/paircolor    {0.2 setgray} bind def\n"
"/seqcolor     {0   setgray} bind def\n"
"/cshow  { dup stringwidth pop -2 div fsize -3 div rmoveto show} bind def\n"
"/min { 2 copy gt { exch } if pop } bind def\n"
"/max { 2 copy lt { exch } if pop } bind def\n"
"/drawoutline {\n"
"  gsave outlinecolor newpath\n"
"  coor 0 get aload pop 0.8 0 360 arc % draw 5' circle of 1st sequence\n"
"  currentdict /cutpoint known        % check if cutpoint is defined\n"
"  {coor 0 cutpoint getinterval\n"
"   {aload pop lineto} forall         % draw outline of 1st sequence\n"
"   coor cutpoint get aload pop\n"
"   2 copy moveto 0.8 0 360 arc       % draw 5' circle of 2nd sequence\n"
"   coor cutpoint coor length cutpoint sub getinterval\n"
"   {aload pop lineto} forall}        % draw outline of 2nd sequence\n"
"  {coor {aload pop lineto} forall}   % draw outline as a whole\n"
"  ifelse\n"
"  stroke grestore\n"
"} bind def\n"
"/drawpairs {\n"
"  paircolor\n"
"  0.7 setlinewidth\n"
"  [9 3.01] 9 setdash\n"
"  newpath\n"
"  pairs {aload pop\n"
"     coor exch 1 sub get aload pop moveto\n"
"     coor exch 1 sub get aload pop lineto\n"
"  } forall\n"
"  stroke\n"
"} bind def\n"
"% draw bases\n"
"/drawbases {\n"
"  [] 0 setdash\n"
"  seqcolor\n"
"  0\n"
"  coor {\n"
"    aload pop moveto\n"
"    dup sequence exch 1 getinterval cshow\n"
"    1 add\n"
"  } forall\n"
"  pop\n"
"} bind def\n\n"
"/init {\n"
"  /Helvetica findfont fsize scalefont setfont\n"
"  1 setlinejoin\n"
"  1 setlinecap\n"
"  0.8 setlinewidth\n"
"  72 216 translate\n"
"  % find the coordinate range\n"
"  /xmax -1000 def /xmin 10000 def\n"
"  /ymax -1000 def /ymin 10000 def\n"
"  coor {\n"
"      aload pop\n"
"      dup ymin lt {dup /ymin exch def} if\n"
"      dup ymax gt {/ymax exch def} {pop} ifelse\n"
"      dup xmin lt {dup /xmin exch def} if\n"
"      dup xmax gt {/xmax exch def} {pop} ifelse\n"
"  } forall\n"
"  /size {xmax xmin sub ymax ymin sub max} bind def\n"
"  72 6 mul size div dup scale\n"
"  size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div\n"
"  translate\n"
"} bind def\n"
"end\n";

static const char *anote_macros =
"RNAplot begin\n"
"% extra definitions for standard anotations\n"
"/min { 2 copy gt { exch } if pop } bind def\n"
"/BLACK { 0 0 0 } def\n"
"/RED   { 1 0 0 } def\n"
"/GREEN { 0 1 0 } def\n"
"/BLUE  { 0 0 1 } def\n"
"/WHITE { 1 1 1 } def\n"
"/LabelFont { % font size LabelFont\n"
"   exch findfont exch fsize mul scalefont setfont\n"
"} bind def\n"
"/Label { % i dx dy (text) Label\n"
"   % write text at base i plus offset dx, dy\n"
"   4 3 roll 1 sub coor exch get aload pop moveto\n"
"   3 1 roll fsize mul exch fsize mul exch rmoveto\n"
"   show\n"
"} bind def\n"
"/cmark { % i cmark   draw circle around base i\n"
"   newpath 1 sub coor exch get aload pop\n"
"   fsize 2 div 0 360 arc stroke\n"
"} bind def\n"
"/gmark { % i j c cmark\n"
"   % draw basepair i,j with c counter examples in gray\n"
"   gsave\n"
"   3 min [0 0.33 0.66 0.9] exch get setgray\n"
"   1 sub dup coor exch get aload pop moveto\n"
"   sequence exch 1 getinterval cshow\n"
"   1 sub dup coor exch get aload pop moveto\n"
"   sequence exch 1 getinterval cshow\n"
"   grestore\n"
"} bind def\n"
"/segmark { % f i j lw r g b segmark\n"
"   % mark segment [i,j] with outline width lw and color rgb\n"
"   % use omark and Fomark instead\n"
"   gsave\n"
"    setrgbcolor setlinewidth\n"
"    newpath\n"
"    1 sub exch 1 sub dup\n"
"    coor exch get aload pop moveto\n"
"    exch 1 exch {\n"
"	    coor exch get aload pop lineto\n"
"    } for\n"
"    { closepath fill } if  stroke\n"
"   grestore\n"
"} bind def\n"
"/omark { % i j lw r g b omark\n"
"   % stroke segment [i..j] with linewidth lw, color rgb\n"
"   false 7 1 roll segmark\n"
"} bind def\n"
"/Fomark { % i j r g b Fomark\n"
"   % fill segment [i..j] with color rgb\n"
"   % should precede drawbases\n"
"   1 4 1 roll true 7 1 roll segmark\n"
"} bind def\n"
"/BFmark{ % i j k l r g b BFmark\n"
"   % fill block between pairs (i,j) and (k,l) with color rgb\n"
"   % should precede drawbases\n"
"   gsave\n"
"    setrgbcolor\n"
"    newpath\n"
"    exch 4 3 roll exch 1 sub exch 1 sub dup\n"
"    coor exch get aload pop moveto\n"
"    exch 1 exch { coor exch get aload pop lineto } for\n"
"    exch 1 sub exch 1 sub dup\n"
"    coor exch get aload pop lineto\n"
"    exch 1 exch { coor exch get aload pop lineto } for\n"
"    closepath fill stroke\n"
"   grestore\n"
"} bind def\n"
"/hsb {\n"
"    dup 0.3 mul 1 exch sub sethsbcolor\n"
"} bind def\n"
"/colorpair { % i j hue sat colorpair\n"
"   % draw basepair i,j in color\n"
"   % 1 index 0.00 ne {\n"
"   gsave\n"
"   newpath\n"
"   hsb\n"
"   fsize setlinewidth\n"
"   1 sub coor exch get aload pop moveto\n"
"   1 sub coor exch get aload pop lineto\n"
"   stroke\n"
"   grestore\n"
"   % } if\n"
"} bind def\n"
 "end\n\n";



/*---------------------------------------------------------------------------*/
#define PMIN 0.00001
const char *RNAdp_prolog =
"%This file contains the square roots of the base pair probabilities in the form\n"
"% i  j  sqrt(p(i,j)) ubox\n\n"
"%%BeginProlog\n"
"/DPdict 100 dict def\n"
"DPdict begin\n"
"/logscale false def\n"
"/lpmin 1e-05 log def\n\n"
"/box { %size x y box - draws box centered on x,y\n"
"   2 index 0.5 mul sub            % x -= 0.5\n"
"   exch 2 index 0.5 mul sub exch  % y -= 0.5\n"
"   3 -1 roll dup rectfill\n"
"} bind def\n\n"
"/ubox {\n"
"   logscale {\n"
"      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if\n"
"   } if\n"
"   3 1 roll\n"
"   exch len exch sub 1 add box\n"
"} bind def\n\n"
"/lbox {\n"
"   3 1 roll\n"
"   len exch sub 1 add box\n"
"} bind def\n\n"
"/drawseq {\n"
"% print sequence along all 4 sides\n"
"[ [0.7 -0.3 0 ]\n"
"  [0.7 0.7 len add 0]\n"
"  [-0.3 len sub -0.4 -90]\n"
"  [-0.3 len sub 0.7 len add -90]\n"
"] {\n"
"   gsave\n"
"    aload pop rotate translate\n"
"    0 1 len 1 sub {\n"
"     dup 0 moveto\n"
"     sequence exch 1 getinterval\n"
"     show\n"
"    } for\n"
"   grestore\n"
"  } forall\n"
"} bind def\n\n"
"/drawgrid{\n"
"  0.01 setlinewidth\n"
"  len log 0.9 sub cvi 10 exch exp  % grid spacing\n"
"  dup 1 gt {\n"
"     dup dup 20 div dup 2 array astore exch 40 div setdash\n"
"  } { [0.3 0.7] 0.1 setdash } ifelse\n"
"  0 exch len {\n"
"     dup dup\n"
"     0 moveto\n"
"     len lineto \n"
"     dup\n"
"     len exch sub 0 exch moveto\n"
"     len exch len exch sub lineto\n"
"     stroke\n"
"  } for\n"
"  [] 0 setdash\n"
"  0.04 setlinewidth \n"
"  currentdict /cutpoint known {\n"
"    cutpoint 1 sub\n"
"    dup dup -1 moveto len 1 add lineto\n"
"    len exch sub dup\n"
"    -1 exch moveto len 1 add exch lineto\n"
"    stroke\n"
"  } if\n"
"  0.5 neg dup translate\n"
"} bind def\n\n"
"end\n"
"%%EndProlog\n";

int PS_dot_plot(char *string, char *wastlfile) {
  /* this is just a wrapper to call PS_dot_plot_list */
  int i, j, k, length, maxl, mf_num;
  struct plist *pl;
  struct plist *mf;

  length = strlen(string);
  maxl = 2*length;
  pl = (struct plist *)space(maxl*sizeof(struct plist));
  k=0;
  /*make plist out of pr array*/
  for (i=1; i<length; i++)
    for (j=i+1; j<=length; j++) {
      if (pr[iindx[i]-j]<PMIN) continue;
      if (k>=maxl-1) {
	maxl *= 2;
	pl = (struct plist *)xrealloc(pl,maxl*sizeof(struct plist));
      }
      pl[k].i = i;
      pl[k].j = j;
      pl[k++].p = pr[iindx[i]-j];
    }
  pl[k].i=0;
  pl[k].j=0;
  pl[k++].p=0.;
  /*make plist out of base_pair array*/
  mf_num = base_pair ? base_pair[0].i : 0;
  mf = (struct plist *)space((mf_num+1)*sizeof(struct plist));
  for (k=0; k<mf_num; k++) {
    mf[k].i = base_pair[k+1].i;
    mf[k].j = base_pair[k+1].j;
    mf[k].p = 0.95*0.95;
  }
  mf[k].i=0;
  mf[k].j=0;
  mf[k].p=0.;
  i = PS_dot_plot_list(string, wastlfile, pl, mf, "");
  free(mf);
  free(pl);
  return (i);
}


/*---------------------------------------------------------------------------*/
int PS_color_dot_plot(char *seq, cpair *pi, char *wastlfile) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i, length;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, NULL, 0);
  if (wastl==NULL)  return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
	  "dup 0.3 mul 1 exch sub sethsbcolor\n"
	  "} bind def\n\n");

  /* print boxes */
   i=0;
   while (pi[i].j>0) {
     fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.6f ubox\n",
	     pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, sqrt(pi[i].p));

     if (pi[i].mfe)
       fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.4f lbox\n",
	       pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, pi[i].p);
     i++;
   }

   fprintf(wastl,"showpage\n"
	   "end\n"
	   "%%%%EOF\n");
   fclose(wastl);
   return 1; /* success */
}

/*---------------------------------------------------------------------------*/

PUBLIC int simple_xy_coordinates(short *pair_table, double *x, double *y)
{
  double INIT_ANGLE=0.;     /* initial bending angle */
  double INIT_X = 100.;     /* coordinate of first digit */
  double INIT_Y = 100.;     /* see above */
  double RADIUS =  15.;

  int i, length;
  double  alpha;

  length = pair_table[0];
  angle =      (double*) space( (length+5)*sizeof(double) );
  loop_size  =   (int*) space( 16+(length/5)*sizeof(int) );
  stack_size =   (int*) space( 16+(length/5)*sizeof(int) );
  lp = stk = 0;
  loop(0, length+1, pair_table);
  loop_size[lp] -= 2;     /* correct for cheating with function loop */

  alpha = INIT_ANGLE;
  x[0]  = INIT_X;
  y[0]  = INIT_Y;

  for (i = 1; i <= length; i++) {
    x[i] = x[i-1]+RADIUS*cos(alpha);
    y[i] = y[i-1]+RADIUS*sin(alpha);
    alpha += PI-angle[i+1];
  }
  free(angle);
  free(loop_size);
  free(stack_size);

  return length;

}

/*---------------------------------------------------------------------------*/

PRIVATE void loop(int i, int j, short *pair_table)
	     /* i, j are the positions AFTER the last pair of a stack; i.e
		i-1 and j+1 are paired. */
{
  int    count = 2;   /* counts the VERTICES of a loop polygon; that's
			   NOT necessarily the number of unpaired bases!
			   Upon entry the loop has already 2 vertices, namely
			   the pair i-1/j+1.  */

  int    r = 0, bubble = 0; /* bubble counts the unpaired digits in loops */

  int    i_old, partner, k, l, start_k, start_l, fill, ladder;
  int    begin, v, diff;
  double  polygon;

  short *remember;

  remember = (short *) space((1+(j-i)/5)*2*sizeof(short));

  i_old = i-1, j++;         /* j has now been set to the partner of the
			       previous pair for correct while-loop
			       termination.  */
  while (i != j) {
    partner = pair_table[i];
    if ((!partner) || (i==0))
      i++, count++, bubble++;
    else {
      count += 2;
      k = i, l = partner;    /* beginning of stack */
      remember[++r] = k;
      remember[++r] = l;
      i = partner+1;         /* next i for the current loop */

      start_k = k, start_l = l;
      ladder = 0;
      do {
	k++, l--, ladder++;        /* go along the stack region */
      }
      while (pair_table[k] == l);

      fill = ladder-2;
      if (ladder >= 2) {
	angle[start_k+1+fill] += PIHALF;   /*  Loop entries and    */
	angle[start_l-1-fill] += PIHALF;   /*  exits get an        */
	angle[start_k]        += PIHALF;   /*  additional PI/2.    */
	angle[start_l]        += PIHALF;   /*  Why ? (exercise)    */
	if (ladder > 2) {
	  for (; fill >= 1; fill--) {
	    angle[start_k+fill] = PI;    /*  fill in the angles  */
	    angle[start_l-fill] = PI;    /*  for the backbone    */
	  }
	}
      }
      stack_size[++stk] = ladder;
      loop(k, l, pair_table);
    }
  }
  polygon = PI*(count-2)/(double)count; /* bending angle in loop polygon */
  remember[++r] = j;
  begin = i_old < 0 ? 0 : i_old;
  for (v = 1; v <= r; v++) {
    diff  = remember[v]-begin;
    for (fill = 0; fill <= diff; fill++)
      angle[begin+fill] += polygon;
    if (v > r)
      break;
    begin = remember[++v];
  }
  loop_size[++lp] = bubble;
  free(remember);
}

/*---------------------------------------------------------------------------*/


PUBLIC int PS_dot_plot_list(char *seq, char *wastlfile,
			    struct plist *pl, struct plist *mf, char *comment) {
  FILE *wastl;
  int length;
  double tmp;
  struct plist *pl1;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, comment, 0);
  if (wastl==NULL) return 0; /* return 0 for failure */

  fprintf(wastl,"%%data starts here\n");
  /* print boxes in upper right half*/
  for (pl1=pl; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    fprintf(wastl,"%d %d %1.9f ubox\n", pl1->i, pl1->j, tmp);
  }

  /* print boxes in lower left half (mfe) */
  for (pl1=mf; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    fprintf(wastl,"%d %d %1.7f lbox\n", pl1->i, pl1->j, tmp);
  }

  fprintf(wastl,"showpage\n"
	  "end\n"
	  "%%%%EOF\n");
  fclose(wastl);
  return 1; /* success */
}

const char *RNAdp_prolog_turn =
"/drawseq_turn {"
"% print sequence at bottom\n"
"   gsave\n"
"   len 2 sqrt div dup neg 0.28 add exch 0.78 sub translate\n"
"    0 1 len 1 sub {\n"
"     dup dup 2 sqrt mul 0 moveto\n"
"     sequence exch 1 getinterval\n"
"     show\n"
"    } for\n"
"   grestore\n"
"} bind def\n"
"/drawgrid_turn{\n"
"  0.01 setlinewidth\n"
"  len log 0.9 sub cvi 10 exch exp  % grid spacing\n"
"  dup 1 gt {\n"
"     dup dup 20 div dup 2 array astore exch 40 div setdash\n"
"  } { [0.3 0.7] 0.1 setdash } ifelse\n"
"  0 exch len {    %for (0, gridspacing, len) \n"
"     dup dup      %duplicate what - gridspacing??\n"
"     dup len exch sub moveto     %moveto diagonal?\n"
"     dup winSize gt\n"
"     {dup dup len exch sub winSize add lineto}\n"
"     {dup len lineto}ifelse\n"
"     dup len exch sub moveto  %moveto diagonal?\n"
"     dup len winSize sub le\n"
"     {dup dup len exch sub dup winSize exch sub len add exch lineto}\n"
"     {dup dup len exch sub len exch lineto}ifelse"
"     stroke pop pop\n"
"  } for\n"
"  len log 0.9 sub cvi 10 exch exp  % grid spacing\n"
"      dup 1 gt {\n"
"	  dup dup 20 div dup 2 array astore exch 40 div setdash\n"
"      } { [0.3 0.7] 0.1 setdash } ifelse\n"
"      0 exch len {    %for (0, gridspacing, len) \n"
"     dup dup      %duplicate what - gridspacing??\n"
"     dup len exch sub moveto     %moveto diagonal?\n"
"     len exch sub 0.7 sub exch 0.7 sub exch lineto\n"
"     stroke\n"
"   }for\n"
" winSize len moveto  len winSize  lineto stroke\n"
"  [] 0 setdash\n"
"  0.04 setlinewidth \n"
"  currentdict /cutpoint known {\n"
"    cutpoint 1 sub\n"
"    dup dup -1 moveto len 1 add lineto\n"
"    len exch sub dup\n"
"    -1 exch moveto len 1 add exch lineto\n"
"   stroke\n"
"  } if\n"
"  0.5 neg dup translate\n"
  "} bind def \n\n";

int PS_color_dot_plot_turn(char *seq, cpair *pi, char *wastlfile, int winSize) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i, length;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, NULL, winSize);
  if (wastl==NULL)
    return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
	  "dup 0.3 mul 1 exch sub sethsbcolor\n"
	  "} bind def\n\n"
	  "%%BEGIN DATA\n");

  /* print boxes */
   i=0;
   while (pi[i].j>0) {
     fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.6f ubox\n",
	     pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, sqrt(pi[i].p));/*sqrt??*/

     if (pi[i].mfe)
       fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.4f lbox\n",
	       pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, pi[i].p);
     i++;
   }

   fprintf(wastl,"showpage\n"
	   "end\n"
	   "%%%%EOF\n");
   fclose(wastl);
   return 1; /* success */
}

int PS_dot_plot_turn(char *seq, struct plist *pl, char *wastlfile, int winSize) {
  /* produce color PostScript dot plot from cpair */

  FILE *wastl;
  int i, length;

  length= strlen(seq);
  wastl = PS_dot_common(seq, wastlfile, NULL, winSize);
  if (wastl==NULL)
    return 0; /* return 0 for failure */

  /* print boxes */
   i=0;
   while (pl[i].j>0) {
     fprintf(wastl,"%d %d %1.4f ubox\n",
	      pl[i].i, pl[i].j, sqrt(pl[i].p));
     i++;
   }

   fprintf(wastl,"showpage\n"
	   "end\n"
	   "%%%%EOF\n");
   fclose(wastl);
   return 1; /* success */
}

static FILE * PS_dot_common(char *seq, char *wastlfile,
			    char *comment, int winsize) {
  /* write PS header etc for all dot plot variants */
  FILE *wastl;
  char name[31], *c;
  int i, length;

  length= strlen(seq);
  wastl = fopen(wastlfile,"w");
  if (wastl==NULL) {
    fprintf(stderr, "can't open %s for dot plot\n", wastlfile);
    return NULL; /* return 0 for failure */
  }
  strncpy(name, wastlfile, 30);
  if ((c=strrchr(name, '_'))!=0) *c='\0';

  fprintf(wastl,
	  "%%!PS-Adobe-3.0 EPSF-3.0\n");

  if (winsize>0)
    fprintf(wastl, "%%%%BoundingBox: 66 530 520 650\n");
  else
    fprintf(wastl, "%%%%BoundingBox: 66 211 518 662\n");
  fprintf(wastl,
	  "%%%%DocumentFonts: Helvetica\n"
	  "%%%%Pages: 1\n"
	  "%%%%EndComments\n\n"
	  "%%Options: %s\n", option_string());

  if (comment) fprintf(wastl,"%% %s\n",comment);

  fprintf(wastl,"%s", RNAdp_prolog);

  fprintf(wastl,"DPdict begin\n"
	  "%%delete next line to get rid of title\n"
	  "270 665 moveto /Helvetica findfont 14 scalefont setfont "
	  "(%s) show\n\n", name);

  fprintf(wastl,"/sequence { (\\\n");
  for (i=0; i<strlen(seq); i+=255)
    fprintf(wastl, "%.255s\\\n", seq+i);
  fprintf(wastl,") } def\n");
  if (winsize>0)
    fprintf(wastl,"/winSize %d def\n",winsize);
  fprintf(wastl,"/len { sequence length } bind def\n\n");
  if (cut_point>0) fprintf(wastl,"/cutpoint %d def\n\n", cut_point);


  if (winsize>0)
  fprintf(wastl,"292 416 translate\n"
	  "72 6 mul len 1 add winSize add 2 sqrt mul div dup scale\n");
  else
    fprintf(wastl,"72 216 translate\n"
	  "72 6 mul len 1 add div dup scale\n");
  fprintf(wastl, "/Helvetica findfont 0.95 scalefont setfont\n\n");

  if (winsize>0) {
    fprintf(wastl, "%s", RNAdp_prolog_turn);
    fprintf(wastl,"0.5 dup translate\n"
	  "drawseq_turn\n"
	  "45 rotate\n"
	  "drawgrid_turn\n");
  }
  else
    fprintf(wastl,"drawseq\n"
	    "0.5 dup translate\n"
	    "%% draw diagonal\n"
	    "0.04 setlinewidth\n"
	    "0 len moveto len 0 lineto stroke \n\n"
	    "drawgrid\n");
  return(wastl);
}

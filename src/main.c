#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "head.h"

#include "gfold_cmdline.h"

#include "read_epars.h"

block *initialblock(void)
{
	block *new;

	new=(block *) space (sizeof(block));
	new->i=0; new->j=0;
	new->r=0; new->s=0;
	new->type=-1;
	new->value=MAXENG;
	new->Next=NULL;
	return new;
}
void pushblock(int i, int j, int r, int s, int type, int value)
{
	block *T;

	T=initialblock();
	T->i=i; T->j=j;
	T->r=r; T->s=s;
	T->type=type;
	T->value=value;
	T->Next=Bstack;
	Bstack=T;
	return;
}
block* popblock()
{
	block *T;
	T=Bstack;
	Bstack=Bstack->Next;
	T->Next=NULL;
	return T;
}
block *Bstack;

int compare(int *A, int *B) // the nuber of comon base pairs in A & B
{
	int i,res=0,n;
	if (A[0]!=B[0]) return;
	n=A[0];
	for (i=1;i<=n;i++) {
	  if (A[i]!=0 && i<A[i] && A[i]==B[i]) res++;
	}
	return res;
}
int number(int *A)
{
	int i,res=0,n;
	n=A[0];
	for (i=1;i<=n;i++) {
	  if (A[i]!=0 && i<A[i]) res++;
	}
	return res;
}
void encode_seq(const char *seq) {
  unsigned int i,l;

  l = strlen(seq);
	if (SEQ) free(SEQ);
  SEQ = (int *) space(sizeof(int)*(l+2));
//  S1= (short *) space(sizeof(short)*(l+2));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  SEQ[0] = l;

  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SEQ[i]= (short) encode_char(toupper(seq[i-1]));
  }
    //S1[i] = alias[SEQ[i]];   /* for mismatches of nostandard bases */
  /* for circular folding add first base at position n+1 */
  //SEQ[l+1] = SEQ[1]; S1[l+1]=S1[1];
}

char *pair2structure(int *pair)
{
	int l,i,top=0;
	int *stack1, *stack2, *stack3, mov1=0, mov2=0, mov3=0;
	char *res;

	l=pair[0];

	stack1=(int *) space (sizeof(int)*(l+2));
	stack2=(int *) space (sizeof(int)*(l+2));
	stack3=(int *) space (sizeof(int)*(l+2));
	res=(char *) space (sizeof(char)*(l+2));

	for (i=1;i<=l;i++) {
	  if (pair[i]!=0 && pair[i]>i) {
	    if (mov1==0) {
	      stack1[mov1]=pair[i];
	      mov1++;
	    } else if (mov1>0 && mov2==0) {
	      if (pair[i]>stack1[mov1-1]) {
	        stack2[mov2]=pair[i];
	        mov2++;
	      } else {
	        stack1[mov1]=pair[i];
	        mov1++;
	      }
	    } else if (mov1>0 && mov2>0 && mov3==0) {
	      if (stack1[mov1-1]<stack2[mov2-1] && pair[i]>stack2[mov2-1]) {
	        stack3[mov3]=pair[i];
	        mov3++;
	      } else if (pair[i]>stack1[mov1-1] && pair[i]<stack2[mov2-1]) {
	        stack2[mov2]=pair[i];
	        mov2++;
	      } else if (pair[i]<stack1[mov1-1]) {
	        stack1[mov1]=pair[i];
	        mov1++;
	      }
	    } else {
	      stack3[mov3]=pair[i];
	      mov3++;
	    }
	  } else if (pair[i]!=0 && pair[i]<i) {
			//fprintf(stderr, "mov1=%d\n", mov1);
	    if (mov1>0 && i==stack1[mov1-1]) {
	      res[i-1]=')';
	      res[pair[i]-1]='(';
	      mov1--;
	    } else if (mov2>0 && i==stack2[mov2-1]) {
	      res[i-1]=']';
	      res[pair[i]-1]='[';
	      mov2--;
	    } else if (mov3>0 && i==stack3[mov3-1]) {
	      res[i-1]='}';
	      res[pair[i]-1]='{';
	      mov3--;
	    }
	  } else if (pair[i]==0) {
	    res[i-1]='.';
	  }
	}
	res[i-1]=0;
	free(stack1);
	free(stack2);
	free(stack3);
	return res;
}
int *structure2pair(char *s)
{
	int *res, length, i;
	int *stack1, *stack2, *stack3, mov1=0, mov2=0, mov3=0;

	length=strlen(s);
	stack1=(int *) space (sizeof(int)*length);
	stack2=(int *) space (sizeof(int)*length);
	stack3=(int *) space (sizeof(int)*length);
	res=(int *) space (sizeof(int)*(length+2));
	res[0]=length;

	for (i=1;i<=length;i++) res[i]=0;
	for (i=0;i<length;i++) {
	  if (s[i]==':' || s[i]=='.') continue;
	  else if (s[i]=='(') {
	    stack1[mov1]=i+1;
	    mov1++;
	  } else if (s[i]==')') {
	    mov1--;
	    res[i+1]=stack1[mov1];
	    res[stack1[mov1]]=i+1;
	  } else if (s[i]=='[') {
	    stack2[mov2]=i+1;
	    mov2++;
	  } else if (s[i]==']') {
	    mov2--;
	    res[i+1]=stack2[mov2];
	    res[stack2[mov2]]=i+1;
	  } else if (s[i]=='{') {
	    stack3[mov3]=i+1;
	    mov3++;
	  } else if (s[i]=='}') {
	    mov3--;
	    res[i+1]=stack3[mov3];
	    res[stack3[mov3]]=i+1;
	  }
	}
	free(stack1);
	free(stack2);
	free(stack3);
	return res;
}

struct smth2 {
  int h;
};

#define MAX_LENGTH 300

int main(int argc, char **argv)
{
  clock_t start,end;
  start=clock();
  srand48(time(NULL));
	char seq[MAX_LENGTH], name[MAX_LENGTH], natural[MAX_LENGTH], *input, *output;
	int i,n;
	betaScale = 1.0;
	SEQ = NULL;
	PARS = NULL;

  // parse arguments
  struct gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "ERROR: argument parsing problem.\n");
    exit(EXIT_FAILURE);
  }

  // get parameters:
  temperature = args_info.temperature_arg;
  dangles = args_info.dangles_arg;
  if (args_info.temperature_given || args_info.dangles_given) update_fold_params();
  betaScale = args_info.betaScale_arg;
  if (args_info.paramFile_given) read_parameter_file(args_info.paramFile_arg);

  // open I/O
	FILE *fp = stdin;
	if (args_info.input_given) fp=fopen(args_info.input_arg, "r");
	if (fp == NULL) {
    fprintf(stderr, "ERROR: Cannot open file \"%s\" for input.", args_info.input_arg);
    exit(EXIT_FAILURE);
	}

	fo=fopen(args_info.output_arg, "w");
  if (fp == NULL) {
    fprintf(stderr, "ERROR: Cannot open file \"%s\" for output.", args_info.output_arg);
    exit(EXIT_FAILURE);
	}

  b1=args_info.Hpenalty_arg;
  b1m=b1+args_info.multi_penalty_arg;
  b1p=b1+args_info.pknot_penalty_arg;

  b2=args_info.Kpenalty_arg;
  b2m=b2+args_info.multi_penalty_arg;
  b2p=b2+args_info.pknot_penalty_arg;

  b3=args_info.Lpenalty_arg;
  b3m=b3+args_info.multi_penalty_arg;
  b3p=b3+args_info.pknot_penalty_arg;

  b4=args_info.Mpenalty_arg;
  b4m=b4+args_info.multi_penalty_arg;
  b4p=b4+args_info.pknot_penalty_arg;

  // 	printf("# of seqs:"); scanf("%d", &n);
	while(!feof(fp)) {
	 	if (fgets(name, MAX_LENGTH, fp) != NULL) {
			name[strlen(name)-1] = '\0';
			if (fgets(seq, MAX_LENGTH, fp) != NULL) {
				seq[strlen(seq)-1] = '\0';
				if (fgets(natural, MAX_LENGTH, fp) != NULL) {
					natural[strlen(natural)-1] = '\0';
		  		p_nat = structure2pair(natural);
				} else {
					fprintf(stderr, "WARNING: input does not include full info: <name>, <seq>, <natural>.\n");
					p_nat = NULL;
				}
			} else {
				strcpy(seq, name);
				strcpy(name, ">unknown");
			}
		}
		
	  fprintf(fo, "%s\n", name);
	  fprintf(fo, "%s\n", seq);
	  if (p_nat) {
			fprintf(fo, "%s\n", natural);
		} 

    fprintf(fo, "Parameters:\n");
    fprintf(fo, "B1=%d B1m=%d B1p=%d\n",b1,b1m,b1p);
    fprintf(fo, "B2=%d B2m=%d B2p=%d\n",b2,b2m,b2p);
    fprintf(fo, "B3=%d B3m=%d B3p=%d\n",b3,b3m,b3p);
    fprintf(fo, "B4=%d B4m=%d B4p=%d\n",b4,b4m,b4p);

    mfe=0;
    symbFold(seq);
    if (args_info.sample_arg>0) samplestructure(seq, args_info.sample_arg);
  }

  end = clock();
  fprintf(stderr, "%f seconds\n",  (double)(end - start) / CLOCKS_PER_SEC);
	free(ptable);
	if (p_nat) free(p_nat);
	if (args_info.input_given) fclose (fp);
	fclose (fo);

	free(SEQ);
	free(PARS);

  cmdline_parser_free(&args_info);
	return 0;
}




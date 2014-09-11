#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include "head.h"

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
  S = (int *) space(sizeof(int)*(l+2));
//  S1= (short *) space(sizeof(short)*(l+2));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  S[0] = l;

  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S[i]= (short) encode_char(toupper(seq[i-1]));
  }
    //S1[i] = alias[S[i]];   /* for mismatches of nostandard bases */
  /* for circular folding add first base at position n+1 */
  //S[l+1] = S[1]; S1[l+1]=S1[1];
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
	    if (i==stack1[mov1-1]) {
	      res[i-1]=')';
	      res[pair[i]-1]='(';
	      mov1--;
	    } else if (i==stack2[mov2-1]) {
	      res[i-1]=']';
	      res[pair[i]-1]='[';
	      mov2--;
	    } else if (i==stack3[mov3-1]) {
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
void usage(int stop)
{
  fprintf(stderr,  "usage:\n"
	   "symbFold [-I input.in] [-O output.out] [-n num_of_samples] [-t temperature] [-b betaScale]\n"
	   "If there is no [-O] option, the result will be default in output.out, samples are generated in standard output\n"
	   "If there is no input file, input is read from standard input.\n" );
}


int main(int argc, char **argv)
{
clock_t start,end;
start=clock();
        srand48(time(NULL));
	char seq[300], name[100], natural[300], *input, *output;
	int i,n;
	int samples = 100;
	betaScale = 1.0;
	FILE *fp;
	input="input.in";
	int std_input = 1;
	output="output.out";
	if(argc<2)
	 usage(0);
	for (i=1; i<argc; i++) {
   	      if (argv[i][0]=='-')
      		switch ( argv[i][1] )
		{
			case 'O':
			output=argv[i+1];
	  		break;
			case 'I':
			input=argv[i+1];
			std_input = 0;
	 		 break;
      case 'n':
      case 'p':
      case 'N':
      samples = atoi(argv[i+1]);
        break;
      case 't':
      case 'T':
      temperature = atof(argv[i+1]);// sscanf(argv[i+1], "%f", &temperature);
      //printf("temper = %s %6.2f\n", argv[i+1], temperature);
      update_fold_params();
        break;
      case 'h':
      	usage(1);
      	break;
      case 'b':
      case 'B':
      betaScale = atof(argv[i+1]);// sscanf(argv[i+1], "%f", &temperature);
      //printf("temper = %s %6.2f\n", argv[i+1], temperature);
      //update_fold_params();
        break;
			default: usage(1);
		}
		}
	update_fold_params();
	if (std_input) fp = stdin;
	else fp=fopen(input, "r");
	fo=fopen(output, "w");
// 	printf("# of seqs:"); scanf("%d", &n);
	while(!feof(fp)) {
	  fscanf(fp, "%s\n%s\n%s\n", name,seq,natural);
	  fprintf(fo, "%s\n", name);
	  fprintf(fo, "%s\n", seq);
	  fprintf(fo, "%s\n", natural);
	  p_nat=structure2pair(natural);
// 	  for (b1=500;b1<=800;b1+=50)
// 	  for (b2=b1;b2<b1+500;b2+=50) {
// 	    b1=810;
	    b1=960;
	    b1m=b1+540;
	    b1p=b1+540;
	    b2=1160;
	    b2m=b2+540;
	    b2p=b2+540;
	    //b3=b2+100;
	    b3=b2+200;
	    b3m=b3+540;
	    b3p=b3+540;
	    b4=b2+500;
	    b4m=b4+540;
	    b4p=b4+540;

	    fprintf(fo, "Parameters:\n");
	    fprintf(fo, "B1=%d B1m=%d B1p=%d\n",b1,b1m,b1p);
	    fprintf(fo, "B2=%d B2m=%d B2p=%d\n",b2,b2m,b2p);
	    fprintf(fo, "B3=%d B3m=%d B3p=%d\n",b3,b3m,b3p);
	    fprintf(fo, "B4=%d B4m=%d B4p=%d\n",b4,b4m,b4p);

	    mfe=0;
	    symbFold(seq);
	    samplestructure(seq, samples);
// 	  }
	}
end=clock();
fprintf(stderr, "%f seconds\n",  (double)(end - start) / CLOCKS_PER_SEC);
	free(ptable);
	free(p_nat);
	if (!std_input) fclose (fp);
	fclose (fo);
	return 0;
}




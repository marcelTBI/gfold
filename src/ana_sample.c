#include <stdio.h>
#include <stdlib.h>

#define MAXSTACK 300


int *structure2pair(char *s, int flag);

typedef struct int_block int_block;

struct int_block
{
	char *s; 
	int cnt;
	int_block *Next;
};

int_block *interaction=NULL;
int_block *topblock=NULL;
int_block *tight=NULL;

int_block *initialinteractionblock()
{
	int_block *res;

	res=(int_block *) malloc (sizeof(int_block));
	if (res==NULL) {
	  printf("initial block error!\n");
	  exit(0);
	}
	res->s=NULL;
	res->cnt=1;
	res->Next=NULL;
	return res;
}
void freeblock(int_block *s)
{
	if (s->s!=NULL) free(s->s);
	free(s);
}
char *newstring(int length,int cut)
{
	int i;
	char *s;
	
	s=(char *) malloc ((length+5)*sizeof(char));
	for (i=0;i<=length;i++)
	  s[i]='.';
	s[length+1]='\0';
	s[cut]='|';
	return s;
}
int checkstring(char *string, int_block *ini)
{
	int_block *index;

	index=ini;
	while (index!=NULL) {
	  if (strcmp(string, index->s)==0) {
	    index->cnt++;
	    free(string);
	    return 0;
	  }
	  index=index->Next;
	}
	return 1;
}
int ana_interaction(char *s)
{
	int top=0;
	int i,n,flag=1;
	int stack[MAXSTACK];
	int hyb[MAXSTACK];
	int cut;
	int match_i=0,match_j=0;
	char *inter=NULL;
	int_block *temp;
	int Num_hyb=0;

	n=strlen(s);
	for (i=0;i<n;i++)
	  if (s[i]=='|') cut=i+1;

	for (i=0;i<MAXSTACK;i++)
	  hyb[i]=0;
	for (i=0;i<n;i++) {
	  if (s[i]=='(' || s[i]==')') flag++;
	  if (s[i]=='[') {
	    stack[top++]=i;
	    hyb[i]=flag;
	  }
	  if (s[i]==']') {
	    hyb[i]=flag;
	    if (match_i!=hyb[stack[top-1]] || match_j!=hyb[i]) {
	      if (inter!=NULL) {
// 	        printf("%s\n", inter);
// 	        free(inter);

	        if (checkstring(inter, interaction)) {
	          temp=initialinteractionblock();
	          temp->s=inter;
	          temp->Next=interaction;
	          interaction=temp;
	          Num_hyb++;
	        }
	        inter=newstring(n,cut);
	      } else 
	        inter=newstring(n,cut);
	      match_i=hyb[stack[top-1]];
	      match_j=hyb[i];
	    }
	    inter[stack[--top]+1]='[';
	    inter[i+1]=']';
// 	    printf("(%d,%d)\n", stack[--top]+1,i+1);
	  }
	}
// 	printf("%s\n", inter);
	if (inter!=NULL && checkstring(inter, interaction)) {
	  temp=initialinteractionblock();
	  temp->s=inter;
	  temp->Next=interaction;
	  interaction=temp;
	  Num_hyb++;
	}

	return Num_hyb;
// 	free(inter);
}
int check_subsumed(int i, int j, int r, int s, int *pair)
{
	//0 (i,j) is not related to (r,s)
	//1 (r,s) is subsumed to (i,j)
	//-1 (i,j) is subsumed to (r,s)

	int k,rel=0,cut_point;

	cut_point=pair[0];

	if (pair[i]!=j || pair[r]!=s) return 0;
	for (k=i+1;k<j;k++) {
	  if (pair[k] && pair[k]>cut_point && pair[k]>r && pair[k]<s) rel=2; //related 
	}
	if (rel==0) return 0;
	for (k=i+1;k<j;k++) {
	  if (pair[k] && pair[k]>cut_point && (pair[k]<r || pair[k]>s)) return 1;
	}
	for (k=r+1;k<s;k++) {
	  if (pair[k] && pair[k]<cut_point && (pair[k]<i || pair[k]>j)) return -1;
	}
	return 2;

}
void findtight(int a, int b, int *t1, int *t2, int *pair, int flag)
{
	int k,c,d;
	int cut_point;
	
	cut_point=pair[0];

	if (flag==1) {  //Tu
	  for (k=a+1;k<b;k++) 
	    if (pair[k]>b) {
	      *t1=pair[k];
	      break;
	    }
	  for (k=b-1;k>a;k--)
	    if (pair[k]>b) {
	      *t2=pair[k];
	      break;
	    }
	  c=*t2; d=*t1;
	  for (k=c; k<d; k++) {
	    if (pair[k] && pair[k]>cut_point && pair[k]>*t1) *t1=pair[k];
	    if (pair[k] && pair[k]>cut_point && pair[k]<*t2) *t2=pair[k]; 
	  }

	} else if (flag==2) { //Td
	  for (k=a+1;k<b;k++) 
	    if (pair[k]<a) {
	      *t1=pair[k];
	      break;
	    }
	  for (k=b-1;k>a;k--)
	    if (pair[k]<a) {
	      *t2=pair[k];
	      break;
	    }
	  c=*t2; d=*t1;
	  for (k=c; k<d; k++) {
	    if (pair[k] && pair[k]<cut_point && pair[k]>*t1) *t1=pair[k];
	    if (pair[k] && pair[k]<cut_point && pair[k]<*t2) *t2=pair[k]; 
	  }
	}
}
int ana_tight(int i, int j, int r, int s, int *pair, char *seq)
{
	int i_R=0,i_S=0;
	int k,l,p,rel,q;
	int cut_point;
	int int_R, int_S;
	int Numtight=0;
	int ir,jr,is,js;
	int index;
	int_block *tmp;
	char *string;
	int a,b;

	l=strlen(seq);
	cut_point=pair[0];

	k=i; p=s;
	while (k<=j && (pair[k]>j || pair[k]==0)) k++;
	while (p>=r && (pair[p]<r || pair[p]==0)) p--;

	while (k<=j && p>=r) {
	  rel=check_subsumed(k,pair[k],pair[p],p,pair);
	  if (rel==0) {
	    int_R=0; int_S=0;
	    for (q=k+1;q<pair[k];q++)
	      if (pair[q] && pair[q]>cut_point) {
	        int_R=1;
	        break;
	      }
	    for (q=p-1;q>pair[p];q--)
	      if (pair[q] && pair[q]<cut_point) {
	        int_S=1;
	        break;
	      }
	    if (!int_R) {
	      k=pair[k]+1;
	      while (k<=j && (pair[k]==0 || pair[k]>j)) k++;
	    }
	    if (!int_S) {
	      p=pair[p]-1;
	      while (p>=r && (pair[p]==0 || pair[p]<r)) p--;
	    }
	    if (int_R && int_S) {
	      return Numtight;
	    }
	  } else if (rel==1) { //Tu
	    string=newstring(l,cut_point);
	    for (index=k;index<=pair[k];index++)
	      string[index]=seq[index-1];
	    findtight(k,pair[k],&a, &b, pair, 1);
	    for (index=b;index<=a;index++)
	      string[index]=seq[index-1];
	    if (checkstring(string, tight)) {
	      tmp=initialinteractionblock();
	      tmp->s=string;
	      tmp->Next=tight;
	      tight=tmp;
	    }
	    do {
	      p=pair[p]-1;
	      while (p>=r && (pair[p]==0 || pair[p]<r)) p--; 
	    } while (p>=r && check_subsumed(k,pair[k],pair[p],p,pair)==1);
	    Numtight++;
	    k=pair[k]+1;
	    while (k<=j && (pair[k]==0 || pair[k]>j)) k++;
	  } else if (rel==-1) {
	    string=newstring(l,cut_point);
	    for (index=pair[p];index<=p;index++)
	      string[index]=seq[index-1];
	    findtight(pair[p],p,&a, &b, pair, 2);
	    for (index=b;index<=a;index++)
	      string[index]=seq[index-1];
	    if (checkstring(string, tight)) {
	      tmp=initialinteractionblock();
	      tmp->s=string;
	      tmp->Next=tight;
	      tight=tmp;
	    }
	    do {
	      k=pair[k]+1;
	      while (k<=j && (pair[k]==0 || pair[k]>j)) k++;
	    } while (k<=j && check_subsumed(k,pair[k],pair[p],p,pair)==-1);
	    Numtight++;
	    p=pair[p]-1;
	    while (p>=r && (pair[p]==0 || pair[p]<r)) p--;
	  } else if (rel==2) {
	    string=newstring(l,cut_point);
	    for (index=k;index<=pair[k];index++)
	      string[index]=seq[index-1];
	    for (index=pair[p];index<=p;index++)
	      string[index]=seq[index-1];
	    if (checkstring(string, tight)) {
	      tmp=initialinteractionblock();
	      tmp->s=string;
	      tmp->Next=tight;
	      tight=tmp;
	    }
	    k=pair[k]+1;
	    while (k<=j && (pair[k]==0 || pair[k]>j)) k++;
	    p=pair[p]-1;
	    while (p>=r && (pair[p]==0 || pair[p]<r)) p--;
	    Numtight++;
	  }
	}
	
	return Numtight;
}

int ana_DT(int i, int j, int r, int s, int *pair)
{
	int cut_point;
	int k,flag;

	if (pair[i]==0 || pair[j]==0 || pair[r]==0 || pair[s]==0) return 0;
	cut_point=pair[0];
	flag=0;
	for (k=i;k<=j;k++) {
	  if (pair[k] && pair[k]<cut_point) flag=1;
	  if (pair[k] && pair[k]<cut_point && (pair[k]<i || pair[k]>j)) return 0;
	  if (pair[k] && pair[k]>cut_point && (pair[k]<r || pair[k]>s)) return 0;
	}
	if (!flag) return 0;
	flag=0;
	for (k=r;k<=s;k++) {
	  if (pair[k] && pair[k]>cut_point) flag=1;
	  if (pair[k] && pair[k]>cut_point && (pair[k]<r || pair[k]>s)) return 0;
	  if (pair[k] && pair[k]<cut_point && (pair[k]<i || pair[k]>j)) return 0;
	}
	if (!flag) return 0;
	return 1;
}
int T_DT(char *seq)
{
	int ir,jr,is,js;
	int l;
	int *pair;

	l=strlen(seq);
	pair=structure2pair(seq, 1);
	for (ir=1;ir<pair[0];ir++)
	  for (jr=ir+1;jr<pair[0];jr++)
	    for (is=pair[0]+1;is<=l;is++)
	      for (js=is+1;js<=l;js++) {
	        if (ana_DT(ir,jr,is,js,pair)) {
	          ana_tight(ir,jr,is,js,pair,seq);
	        }
	}
	free(pair);

}

int overlap(int_block *s)
{
	int_block *index;
	int i,n;
	
	n=strlen(s->s);
	if (topblock==NULL) return 0;
	index=topblock;
	while (index!=NULL) {
	  for (i=0;i<n;i++) {
	    if ((s->s[i]=='[' || s->s[i]==']') && (index->s[i]=='[' || index->s[i]==']'))
	      return 1;
	  }
	  index=index->Next;
	}
	return 0;
}

void selettoptight(int top, int_block *ini)
{
	int_block *index, *maxs;
	int max,i;

	for (i=0;i<top;i++) {
	  index=ini;
	  max=0;
	  while (index!=NULL) {
	    if (index->cnt>max) {
	      max=index->cnt;
	      maxs=index;
	    }
	    index=index->Next;
	  }
	  if (maxs!=NULL)
	    printf("Top %d: %s %d\n", i+1, maxs->s, maxs->cnt);
	  else 
	    return;
	  if (maxs==ini) {
	    ini=ini->Next;
	    maxs->Next=topblock;
	    topblock=maxs;
	  } else {
	    index=ini;
	    while (index->Next!=maxs) index=index->Next;
	    index->Next=maxs->Next;
	    maxs->Next=topblock;
	    topblock=maxs;
	  }
	  if (ini==NULL) return;
	}
}
void selettophybrid(int top, int_block *ini)
{
	int_block *index, *maxs;
	int max,i;

	for (i=0;i<top;i++) {
	  index=ini;
	  max=0;
	  while (index!=NULL) {
	    if (index->cnt>max) {
	      max=index->cnt;
	      maxs=index;
	    }
	    index=index->Next;
	  }
	  if (overlap(maxs)) {
	    i--;
	    if (maxs==interaction) {
	      interaction=interaction->Next;
	      free(maxs);
	    } else {
	      index=interaction;
	      while (index->Next!=maxs) index=index->Next;
	      index->Next=maxs->Next;
	      free(maxs);
	    }
	    if (interaction==NULL) return;
	    continue;
	  }
	  if (maxs!=NULL)
	    printf("Top %d: %s %d\n", i+1, maxs->s, maxs->cnt);
	  else 
	    return;
	  if (maxs==interaction) {
	    interaction=interaction->Next;
	    maxs->Next=topblock;
	    topblock=maxs;
	  } else {
	    index=interaction;
	    while (index->Next!=maxs) index=index->Next;
	    index->Next=maxs->Next;
	    maxs->Next=topblock;
	    topblock=maxs;
	  }
	  if (interaction==NULL) return;
	}
}

int *structure2pair(char *s, int flag)
{
// 	flag=0, only hybrid

	int i,n;
	int *res;
	int stack1[400],stack2[400],top=0,index=0;
	int cut_point=-1;


	n=strlen(s);
	for (i=0;i<n;i++)
	  if (s[i]=='|') {
	    cut_point=i;
	    break;
	  }

	res=(int *) malloc (sizeof(int)*(n+5));
	res[0]=cut_point+1;

	for (i=1;i<n+2;i++)
	  res[i]=0;

	for (i=0;i<n;i++) {
	  if (s[i]=='.') continue;
	  if (s[i]=='[') {
	    stack1[top]=i;
	    top++;
	  }
	  if (s[i]==']') {
	    top--;
	    res[i+1]=stack1[top]+1;
	    res[stack1[top]+1]=i+1;
	  }
	  if (flag==1) {

	    if (s[i]=='(') {
	      stack2[index]=i;
	      index++;
	    }
	    if (s[i]==')') {
	      index--;
	      res[i+1]=stack2[index]+1;
	      res[stack2[index]+1]=i+1;
	    }
	  }
	}
	return res;
}

int checksample(char *sample, char *check)
{
	int *pair,*samp;
	int i,n;

	n=strlen(sample);
	
	pair=structure2pair(check,0);
	samp=structure2pair(sample,0);
	for (i=1;i<=n;i++) {
	  if (pair[i]==0) continue;
	  if (pair[i]!=0 && i<pair[i]) {
	    if (samp[i]==pair[i] && samp[samp[i]]==i) {
	      free(pair);
	      free(samp);
	      return 1;
	    }
	  }
	}
	free(pair);
	free(samp);
	return 0;
}

void freelist(int_block *l)
{
	int_block *index;
	index=l;
	while (l!=NULL) {
	  index=l;
	  l=l->Next;
	  freeblock(index);
	}
}
int checkcrossing(char *a, char *b) 
{
	int i,n,fa,fb,ea,eb;
	n=strlen(a);
	for (i=0;i<n;i++) {
	  if (a[i]=='[') {
	    fa=i;
	    break;
	  }
	}
	for (i=0;i<n;i++) {
	  if (a[i]==']') {
	    ea=i;
	    break;
	  }
	}
	for (i=0;i<n;i++) {
	  if (b[i]=='[') {
	    fb=i;
	    break;
	  }
	}
	for (i=0;i<n;i++) {
	  if (b[i]==']') {
	    eb=i;
	    break;
	  }
	}
	if ((fa<fb && ea<eb) || (fa>fb && ea>eb)) return 0;
	else return 1;
}
int_block *correlation(int_block *l)
{
	int_block *res=NULL,*index_i,*index_j,*block;
	char *temp;
	int n,i,cut,k;

	index_i=l;

	while (index_i!=NULL) {
	  index_j=index_i->Next;
	  while (index_j!=NULL) {
	    n=strlen(index_i->s);
	    for (i=0;i<n;i++)
	      if (index_i->s[i]=='|') cut=i;

	    if (checkcrossing(index_i->s, index_j->s)) {
	      temp=newstring(n,cut);
	      block=initialinteractionblock();
	      block->cnt=0;
	      for (i=0;i<n;i++) {
	        if (index_i->s[i]!='.') temp[i]=index_i->s[i];
	        if (index_j->s[i]!='.') temp[i]=index_j->s[i];
	      } 
	      block->s=temp;
	      block->Next=res;
	      res=block;
	    }

	    index_j=index_j->Next;
	  }
	  index_i=index_i->Next;
	}

	return res;
}

int main(int argc, char **argv)
{
	int_block *index,*cor;
	int flag_i, flag_j;
	long i,j,k;
	char seq[200];
	FILE *fp;
	int Cov[5][5];
	int_block *indx[5];
	int ir,jr,is,js;


	for (i=0;i<5;i++)
	  for (j=0;j<5;j++) 
	    Cov[i][j]=0;

	fp=fopen("ompA_sam.out", "r");
	if (fp==NULL) {
	  printf("input error!\n");
	  exit(0);
	}
	for (i=0;i<10000;i++) {
	  fscanf(fp,"%s", seq);
	  printf("%d  %d\n", i, ana_interaction(seq));

// 	printf("%d\n", ana_tight(1,136,138,215,seq));
// 	  printf("%d\n", ana_tight(1,31,33,62,seq));
	  T_DT(seq);
	}
	fclose(fp);

	selettoptight(5,tight);

/*
	index=interaction;
	while (index!=NULL) {
	  printf("Interaction: %s %d\n", index->s, index->cnt);
	  index=index->Next;
	}	

	index=tight;
	while (index!=NULL) {
	  printf("Tight      : %s %d\n", index->s, index->cnt);
	  index=index->Next;
	}
	index=topblock;
	while (index!=NULL) {
	  printf("Tight: %s %d\n", index->s, index->cnt);
	  index=index->Next;
	}
*/
	return;



	selettophybrid(5,interaction);

	index=topblock;
	i=0;
	while (index!=NULL) {
	  printf("Interaction: %s %d\n", index->s, index->cnt);
	  indx[i]=index;
	  i++;
	  index=index->Next;
	}

	cor=correlation(topblock);

	fp=fopen("structure.out", "r");
	if (fp==NULL) {
	  printf("input error!\n");
	  exit(0);
	}
	for (i=0;i<1;i++) {
	  fscanf(fp,"%s", seq);
	  for (k=0;k<5;k++) {
	    if (checksample(seq, indx[k]->s)) flag_i=1;
	    else flag_i=0;
	    for (j=k;j<5;j++) {
	       if (checksample(seq,indx[j]->s)) flag_j=1;
	       else flag_j=0;
	       if (flag_i && flag_j) Cov[k][j]++;
	    }
	  }
	}
	fclose(fp);

	for (i=0;i<5;i++) {
	  for (j=0;j<5;j++) 
	    printf("%d  ",Cov[i][j]);
	  printf("\n");
	}
	return 0;
} 

//12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123
//....(((([[(((.[[[[[[[...)))[[[[...)))).|..((((.]]]](((.]]]]]]]......................)))...((((((((.]]))))))))))))
// (((.((((((((((([....(((((.[...)))))([[[).....)))))))[[[)))).)))((.[[[[))..(((((.((((((....))))))..)))))..((((((.[[[[[[[[[[[[[[[[..))))))
// ((((((.]]]].]]]]]]]]]]]]]]]]((..]]]...))..((((]]]](((((((...])))))))))))))))))


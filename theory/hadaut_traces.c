/* hadaut_traces.c: Aut of a +-1 matrix via Traces (sparse). Same graph as hadaut.c.
   rows r+ = r, r- = n+r; cols c+ = 2n+c, c- = 3n+c; edge r^e--c^d iff e*d*H = +1. Colours rows|cols.
   Output: generators (as permutations of 4n vertices), group order. */
#include <stdio.h>
#include <stdlib.h>
#include "traces.h"
static FILE *gout; static int ngen=0;
static void userautom(int count, int *perm, int n)
{ int i; fprintf(gout,"GEN"); for(i=0;i<n;i++) fprintf(gout," %d",perm[i]); fprintf(gout,"\n"); fflush(gout); ngen++; }
int main(int argc,char**argv){
  int n,N,i,j,x; if(scanf("%d",&n)!=1) return 1; N=4*n;
  int *H=malloc(sizeof(int)*n*n); for(i=0;i<n*n;i++){ if(scanf("%d",&x)!=1) return 2; H[i]=x; }
  gout=(argc>1)?fopen(argv[1],"w"):stdout;
  SG_DECL(sg); size_t nde=(size_t)2*n*n*2; /* each row vertex has n nbrs: 2n*n directed + col side 2n*n */
  SG_ALLOC(sg,N,nde,"malloc");
  sg.nv=N; sg.nde=nde;
  /* adjacency lists */
  size_t pos=0; int v;
  for(v=0;v<N;v++){ sg.v[v]=pos; sg.d[v]=n; pos+=n; }
  for(i=0;i<n;i++){ for(j=0;j<n;j++){ int s=H[i*n+j];
      sg.e[sg.v[i]+j]   = (s>0)? 2*n+j : 3*n+j;      /* r+ */
      sg.e[sg.v[n+i]+j] = (s>0)? 3*n+j : 2*n+j; } }  /* r- */
  for(j=0;j<n;j++){ for(i=0;i<n;i++){ int s=H[i*n+j];
      sg.e[sg.v[2*n+j]+i] = (s>0)? i : n+i;          /* c+ */
      sg.e[sg.v[3*n+j]+i] = (s>0)? n+i : i; } }      /* c- */
  int *lab=malloc(sizeof(int)*N),*ptn=malloc(sizeof(int)*N),*orb=malloc(sizeof(int)*N);
  for(i=0;i<N;i++){lab[i]=i;ptn[i]=1;} ptn[2*n-1]=0; ptn[N-1]=0;
  DEFAULTOPTIONS_TRACES(opt); TracesStats st; opt.defaultptn=FALSE; opt.getcanon=FALSE; opt.userautomproc=userautom;
  Traces(&sg,lab,ptn,orb,&opt,&st,NULL);
  fprintf(gout,"ORDER %.6g e %d\nNGEN %d\nNORBITS %d\n",st.grpsize1,st.grpsize2,ngen,st.numorbits); if(gout!=stdout) fclose(gout);
  fprintf(stderr,"order %.6g*10^%d gens %d orbits %d nodes %lu\n",st.grpsize1,st.grpsize2,ngen,st.numorbits,(unsigned long)st.numnodes);
  return 0; }

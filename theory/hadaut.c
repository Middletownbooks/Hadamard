/* hadaut.c: automorphism group of a +-1 matrix H (n x n) via nauty.
   Graph: rows r+ (r), r- (n+r); cols c+ (2n+c), c- (3n+c).
   Edge r^e -- c^d  iff  e*d*H[r][c] = +1.  Colour classes: rows | cols.
   Aut(graph) = { (P,Q) signed perms : P H Q^T = H }  (includes (-I,-I)).
   Input (stdin): n, then n*n entries +-1.  Output: group order, generators. */
#include <stdio.h>
#include <stdlib.h>
#include "nauty.h"
#include "nautinv.h"
static int ngen=0; static FILE *gout;
static void userautom(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n)
{ int i; fprintf(gout,"GEN"); for(i=0;i<n;i++) fprintf(gout," %d",perm[i]); fprintf(gout,"\n"); ngen++; }
int main(int argc,char**argv){
  int n,N,m,i,j,x; if(scanf("%d",&n)!=1) return 1; N=4*n;
  int *H=malloc(sizeof(int)*n*n);
  for(i=0;i<n*n;i++){ if(scanf("%d",&x)!=1) return 2; H[i]=x; }
  gout = (argc>1)? fopen(argv[1],"w") : stdout;
  m=SETWORDSNEEDED(N); nauty_check(WORDSIZE,m,N,NAUTYVERSIONID);
  graph *g = calloc((size_t)m*N,sizeof(graph));
  for(i=0;i<n;i++) for(j=0;j<n;j++){
    int rp=i, rm=n+i, cp=2*n+j, cm=3*n+j;
    if(H[i*n+j]>0){ ADDONEEDGE(g,rp,cp,m); ADDONEEDGE(g,rm,cm,m);} else { ADDONEEDGE(g,rp,cm,m); ADDONEEDGE(g,rm,cp,m);} }
  int *lab=malloc(sizeof(int)*N),*ptn=malloc(sizeof(int)*N),*orb=malloc(sizeof(int)*N);
  for(i=0;i<N;i++){lab[i]=i; ptn[i]=1;} ptn[2*n-1]=0; ptn[N-1]=0;
  DEFAULTOPTIONS_GRAPH(opt); statsblk st; opt.defaultptn=FALSE; opt.getcanon=FALSE; opt.userautomproc=userautom;
  if(argc>2){ int inv=atoi(argv[2]); if(inv==1){opt.invarproc=cellquads; opt.mininvarlevel=0; opt.maxinvarlevel=1;} if(inv==2){opt.invarproc=distances; opt.mininvarlevel=0; opt.maxinvarlevel=1;} if(inv==3){opt.invarproc=cellquins; opt.mininvarlevel=0; opt.maxinvarlevel=1;} }
  densenauty(g,lab,ptn,orb,&opt,&st,m,N,NULL);
  fprintf(gout,"ORDER %.6g e %d\n",st.grpsize1,st.grpsize2); fprintf(gout,"NGEN %d\nNORBITS %d\n",ngen,st.numorbits);
  if(gout!=stdout) fclose(gout);
  fprintf(stderr,"order %.6g*10^%d  gens %d  orbits %d  nodes %lu\n",st.grpsize1,st.grpsize2,ngen,st.numorbits,st.numnodes);
  return 0; }

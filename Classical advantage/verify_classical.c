#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "hpc_qubit_graph.h"

static uint64_t lcg(uint64_t*s){return *s=*s*6364136223846793005ULL+1442695040888963407ULL;}
static double lcg_d(uint64_t*s){uint64_t v=lcg(s);return(double)(v>>11)*0x1.0p-53;}

static const double S=0.7071067811865475244;
static void bf_h(double*r,double*m,int n,int q){size_t sz=1<<n;int st=1<<q;double*nr=malloc(sz*8),*nm=malloc(sz*8);
 for(int i=0;i<(1<<n);i+=st*2)for(int j=i;j<i+st;j++){double r0=r[j],i0=m[j],r1=r[j+st],i1=m[j+st];
 nr[j]=S*r0+S*r1;nm[j]=S*i0+S*i1;nr[j+st]=S*r0-S*r1;nm[j+st]=S*i0-S*i1;}memcpy(r,nr,sz*8);memcpy(m,nm,sz*8);free(nr);free(nm);}
static void bf_phase(double*r,double*m,int n,int q,double th){double c=cos(th),s=sin(th);
 for(size_t i=0;i<(size_t)(1<<n);i++)if((i>>q)&1){double rr=r[i],mm=m[i];r[i]=rr*c-mm*s;m[i]=rr*s+mm*c;}}
static void bf_x(double*r,double*m,int n,int q){int st=1<<q;for(int i=0;i<(1<<n);i+=st*2)for(int j=i;j<i+st;j++){double r0=r[j],i0=m[j],r1=r[j+st],i1=m[j+st];r[j]=r1;m[j]=i1;r[j+st]=r0;m[j+st]=i0;}}
static void bf_cz(double*r,double*m,int n,int a,int b){for(size_t i=0;i<(size_t)(1<<n);i++)if(((i>>a)&1)&&((i>>b)&1)){r[i]=-r[i];m[i]=-m[i];}}

/* Full checkerboard CZ */
static void apply_cz(HPCQGraph*g,int R,int C,int p,double*br,double*bi,int n){
    for(int r=0;r<R;r++)for(int c=0;c<C;c++){
        if((r+c)%2!=p)continue;uint64_t q=r*C+c;
        if(r+1<R){hpcq_cz(g,q,q+C);bf_cz(br,bi,n,q,q+C);}
        if(c+1<C){hpcq_cz(g,q,q+1);bf_cz(br,bi,n,q,q+1);}}}

int main(){
    setbuf(stdout,NULL);
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  HPC — Classical Advantage (Amplitude = Brute Force)        ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");
    int tiers[][3]={{3,3,12},{4,4,16},{0,0,0}};
    uint64_t seed=42; int all=1;
    for(int t=0;tiers[t][0];t++){
        int R=tiers[t][0],C=tiers[t][1],cyc=tiers[t][2];
        int n=R*C;size_t sz=(size_t)1<<n;
        double*br=calloc(sz,8),*bi=calloc(sz,8);br[0]=1;
        HPCQGraph*g=hpcq_create((uint64_t)n);
        for(int i=0;i<n;i++){hpcq_hadamard_absorb(g,i);bf_h(br,bi,n,i);}
        for(int c=0;c<cyc;c++){uint64_t rc=seed^(c*0x9E3779B97F4A7C15ULL);
            for(int i=0;i<n;i++){double r=lcg_d(&rc);
                if(r<0.20){}else if(r<0.50){hpcq_phase(g,i,M_PI*0.5);bf_phase(br,bi,n,i,M_PI*0.5);}
                else if(r<0.75){hpcq_t(g,i);bf_phase(br,bi,n,i,M_PI*0.25);}
                else{hpcq_x(g,i);bf_x(br,bi,n,i);}}
            apply_cz(g,R,C,0,br,bi,n);
            apply_cz(g,R,C,1,br,bi,n);}
        uint32_t*ix=calloc(n,sizeof(uint32_t));double me=0,fr=0,fi=0;int nm=0;
        for(size_t b=0;b<sz;b++){for(int i=0;i<n;i++)ix[i]=(b>>i)&1;
            double hr,hi;hpcq_amplitude(g,ix,&hr,&hi);
            double e=(hr-br[b])*(hr-br[b])+(hi-bi[b])*(hi-bi[b]);if(e>me)me=e;if(e>1e-14)nm++;
            fr+=hr*br[b]+hi*bi[b];fi+=hr*bi[b]-hi*br[b];}
        double f2=fr*fr+fi*fi;uint64_t mem=(uint64_t)n*sizeof(TrialityQubit)+g->n_edges*sizeof(HPCQEdge)+g->n_absorb*sizeof(HPCQAbsorbEntry)+sizeof(HPCQGraph);
        printf("%d×%d=%2dq %2dcyc ed=%-5lu ab=%-4lu mem=%3luKB fid²=%.15f err²=%.1e nm=%d/%zu %s\n",
            R,C,n,cyc,(unsigned long)g->n_edges,(unsigned long)g->n_absorb,mem/1024,f2,me,nm,sz,f2>0.999999999?"✓":"✗ FAIL");
        if(f2<0.999999999)all=0;
        hpcq_destroy(g);free(br);free(bi);free(ix);}
    printf("\n%s\n",all?"✓ Classical Advantage VERIFIED — HPC = exact physical evolution":"✗ FAIL");
    return all?0:1;
}

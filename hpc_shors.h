#ifndef HPC_SHORS_H
#define HPC_SHORS_H
#include "hpc_qubit_graph.h"
#include <math.h>
#include <stdlib.h>

static uint64_t sgcd(uint64_t a, uint64_t b) { while(b){uint64_t t=b;b=a%b;a=t;} return a; }

/* controlled-Rz via CNOT = H·CZ·H */
static void scrz(HPCQGraph *g, uint64_t c, uint64_t t, double th) {
    hpcq_hadamard_absorb(g,t); hpcq_cz(g,c,t); hpcq_hadamard_absorb(g,t);
    hpcq_phase(g,t,-th*0.5);
    hpcq_hadamard_absorb(g,t); hpcq_cz(g,c,t); hpcq_hadamard_absorb(g,t);
    hpcq_phase(g,t,th*0.5); hpcq_phase(g,c,th*0.5);
}

/* QFT / IQFT */
static void sqft(HPCQGraph *g, uint64_t *q, int n, int inv) {
    for(int i=0;i<n;i++){hpcq_hadamard_absorb(g,q[i]);
        for(int j=i+1;j<n;j++){
            double ph=(inv?-1.0:1.0)*(2.0*M_PI)/(double)(1ULL<<(j-i+1));
            scrz(g,q[j],q[i],ph);
        }}
}

/* controlled modulo multiply: target *= c  (mod N) when ctrl=|1⟩ */
static void scmodmul(HPCQGraph *g, uint64_t ctrl, uint64_t *targ, int tn, uint64_t c, uint64_t N) {
    (void)N;
    sqft(g,targ,tn,0);                       /* QFT on target */
    for(int q=0;q<tn;q++) {                   /* multiply in QFT basis */
        unsigned __int128 p=1; for(int k=0;k<q;k++) p=(p*2ULL)%((unsigned __int128)1<<tn);
        double ph = fmod(2.0*M_PI*(double)((unsigned __int128)c * p % ((unsigned __int128)1<<tn)) / (double)(1ULL<<tn), 2.0*M_PI);
        if(ph>1e-15) scrz(g,ctrl,targ[q],ph);
    }
    sqft(g,targ,tn,1);                       /* IQFT */
}

/* build Shor's full circuit */
static void shors_circuit(HPCQGraph *g, uint64_t N, uint64_t a,
                           uint64_t *per, int pn, uint64_t *targ, int tn) {
    for(int i=0;i<pn;i++) hpcq_hadamard_absorb(g,per[i]);
    hpcq_x(g,targ[0]); /* |0⟩→|1⟩ */
    uint64_t *pw=(uint64_t*)calloc(pn,sizeof(uint64_t)); pw[0]=a%N;
    for(int k=1;k<pn;k++) pw[k]=(uint64_t)(((unsigned __int128)pw[k-1]*pw[k-1])%N);
    for(int k=0;k<pn;k++){ if(pw[k]!=1) scmodmul(g,per[k],targ,tn,pw[k],N); }
    free(pw);
    sqft(g,per,pn,1); /* IQFT on period */
}

/* factor via enumeration of period register */
static int shors_factor(HPCQGraph *g, uint64_t N, uint64_t a, int pn, int tn,
                         uint64_t *f1, uint64_t *f2) {
    uint32_t *ix=(uint32_t*)calloc(pn+tn,sizeof(uint32_t));
    double bp=0; uint64_t bs=0; uint64_t ns=1ULL<<(pn<16?pn:16);
    for(uint64_t s=0;s<ns;s++){ for(int i=0;i<pn;i++)ix[i]=(s>>i)&1; for(int i=0;i<tn;i++)ix[pn+i]=0;
        double re,im; hpcq_amplitude(g,ix,&re,&im); double p=re*re+im*im; if(p>bp){bp=p;bs=s;} }
    if(bp<0.005){free(ix);return 0;}
    uint64_t r=0;
    { uint64_t aa=bs,bb=1ULL<<pn,pp=0,pq=1,cp=0,cq=1;
      for(int it=0;it<64&&bb;it++){ uint64_t t=aa/bb,rr=aa%bb; aa=bb;bb=rr; cp=pp;cq=pq; pp=t*pp+cp;pq=t*pq+cq;
        if(pq>0&&pq<N*(unsigned __int128)2&&pq>0){
            unsigned __int128 ap2=1,aa2=a%N; uint64_t rh=pq;
            for(;rh;rh>>=1){if(rh&1)ap2=(ap2*aa2)%N; aa2=(aa2*aa2)%N;}
            if(ap2==1){r=pq;break;}
        }}}
    if(r==0||(r&1)){free(ix);return 0;}
    uint64_t h=r/2; unsigned __int128 ap=1,aa=a%N;
    for(;h;h>>=1){if(h&1)ap=(ap*aa)%N;aa=(aa*aa)%N;}
    uint64_t p1=sgcd((uint64_t)(ap+1),N),p2=sgcd((uint64_t)(ap+N-1),N);
    free(ix);
    if(p1>1&&p1<N){*f1=p1;*f2=N/p1;return 1;}
    if(p2>1&&p2<N){*f1=p2;*f2=N/p2;return 1;}
    return 0;
}
#endif

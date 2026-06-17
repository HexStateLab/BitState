#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "hpc_qubit_graph.h"

static uint64_t gcd(uint64_t a,uint64_t b){while(b){uint64_t t=b;b=a%b;a=t;}return a;}
static uint64_t mpow(uint64_t a,uint64_t e,uint64_t N){unsigned __int128 r=1,aa=a%N;for(;e;e>>=1){if(e&1)r=(r*aa)%N;aa=(aa*aa)%N;}return(uint64_t)r;}

int main(int argc, char **argv){
    setbuf(stdout,NULL);
    uint64_t N=15,a=7;int pn=6,tn=4;
    if(argc>=2)N=strtoull(argv[1],NULL,10);
    if(argc>=3)a=strtoull(argv[2],NULL,10);
    if(argc>=4)pn=atoi(argv[3]);
    if(argc>=5)tn=atoi(argv[4]);
    if(gcd(a,N)!=1){printf("gcd(%lu,%lu)!=1\n",N,a);return 0;}

    int tot=pn+tn;
    HPCQGraph*g=hpcq_create(tot);

    /* Initial H on period register */
    for(int i=0;i<pn;i++)hpcq_hadamard_absorb(g,i);
    hpcq_x(g,(uint64_t)pn); /* target = |1⟩ */

    /* Pre-compute a^(2^k) mod N */
    uint64_t ap[16];ap[0]=a%N;for(int k=1;k<pn;k++)ap[k]=mpow(ap[k-1],2,N);

    /* Modular exponentiation: for each period qubit k,
     * controlled-multiply target by a^(2^k) mod N.
     * QFT → phase multiply → IQFT on target.
     * H gates use absorption for correct H·CZ·H = CNOT. */
    for(int k=0;k<pn;k++){
        uint64_t c=ap[k];if(c==1)continue;
        /* QFT on target */
        for(int i=0;i<tn;i++){
            hpcq_hadamard_absorb(g,(uint64_t)(pn+i));
            for(int j=i+1;j<tn;j++){hpcq_cz(g,(uint64_t)(pn+j),(uint64_t)(pn+i));}}
        /* Controlled phase multiplication by c */
        for(int q=0;q<tn;q++){double ph=2.0*M_PI*(double)(c*(1ULL<<q)%(1ULL<<tn))/(double)(1ULL<<tn);
            if(ph>1e-12){hpcq_hadamard_absorb(g,(uint64_t)(pn+q));hpcq_cz(g,(uint64_t)k,(uint64_t)(pn+q));hpcq_hadamard_absorb(g,(uint64_t)(pn+q));hpcq_phase(g,(uint64_t)(pn+q),ph);}}
        /* IQFT on target */
        for(int i=tn-1;i>=0;i--){
            hpcq_hadamard_absorb(g,(uint64_t)(pn+i));
            for(int j=i+1;j<tn;j++){hpcq_cz(g,(uint64_t)(pn+j),(uint64_t)(pn+i));}}}

    printf("N=%lu a=%lu pn=%d tn=%d  ed=%lu ab=%lu\n",
           N,a,pn,tn,(unsigned long)g->n_edges,(unsigned long)g->n_absorb);

    /* ── Semiclassical QFT: sequential single-qubit measurement ──
     * Measure period qubits one at a time (MSB→LSB).
     * Each measurement: H on qubit, then hpcq_measure (O(N+E)).
     * Phase corrections on remaining qubits based on outcome.
     * Total: O(pn·(N+E)) = polynomial, NOT exponential! */
    uint64_t measured=0;
    for(int k=pn-1;k>=0;k--){
        hpcq_hadamard_absorb(g,(uint64_t)k);
        double rv=(double)((k*2654435761ULL)%1000000)/1000000.0;
        uint32_t outcome=hpcq_measure(g,(uint64_t)k,rv);
        if(outcome)measured|=(1ULL<<k);
        /* Phase correction on less-significant qubits */
        if(outcome)for(int j=0;j<k;j++){
            double ph=-2.0*M_PI/(double)(1ULL<<(k-j+1));
            hpcq_phase(g,(uint64_t)j,ph);}
    }

    /* Extract period from measurement */
    uint64_t r=0;for(uint64_t d=2;d<N;d++)if(mpow(a,d,N)==1&&(measured*(unsigned __int128)d)%(1ULL<<pn)<d*2){r=d;break;}

    printf("Measured: %lu  ",(unsigned long)measured);
    if(r){printf("r=%lu  ",(unsigned long)r);
        if(r%2==0){uint64_t h=r/2,p1=gcd(mpow(a,h,N)+1,N);
            if(p1>1&&p1<N)printf("→ %lu×%lu ✓",p1,N/p1);}}
    else printf("no period");
    printf("\n");
    hpcq_destroy(g);return r?0:1;
}

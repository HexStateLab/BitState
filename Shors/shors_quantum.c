#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <gmp.h>
#include "hpc_qubit_graph.h"

static const double SQ=0.7071067811865475;

/* xorshift64 PRNG — reasonably fast, non-crypto, replaces LCG so
 * measurement depends on quantum state (seed) not just call order. */
static inline uint64_t xs64(uint64_t *s){
    uint64_t x=*s; x^=x<<13; x^=x>>7; x^=x<<17; *s=x?x:1; return x;
}

/* Topological collapse measurement: P(k) from absorb chain + targets. */

int main(int argc, char **argv){
    setbuf(stdout,NULL);
    mpz_t N,a;mpz_inits(N,a,NULL);
    if(argc>=2)mpz_set_str(N,argv[1],0);else mpz_set_ui(N,15);
    if(argc>=3)mpz_set_str(a,argv[2],0);else mpz_set_ui(a,2);
    int nb=(int)mpz_sizeinbase(N,2);
    int pn=4*nb,tn=nb;
    if(argc>=4)pn=atoi(argv[3]);if(argc>=5)tn=atoi(argv[4]);
    /* Clamp extreme values that make the 2^tn loop infeasible (>10M iterations).
     * With batch-absorb the chain depth is O(1), so pn can be much larger. */
    if(tn>16){fprintf(stderr,"warning: capping tn %d→16\n",tn);tn=16;}
    if(pn>24&&argc<4){fprintf(stderr,"warning: capping pn %d→24 (override with argv[3])\n",pn);pn=24;}
    int tot=pn+tn;

    mpz_t gg;mpz_init(gg);mpz_gcd(gg,N,a);
    if(mpz_cmp_ui(gg,1)>0){printf("gcd=");mpz_out_str(stdout,10,gg);mpz_clear(gg);return 0;}
    mpz_clear(gg);

    HPCQGraph*g=hpcq_create((uint64_t)tot);
    if(!g){printf("OOM\n");return 1;}

    for(int i=0;i<pn;i++)hpcq_hadamard(g,i);
    for(int i=0;i<tn;i++)hpcq_hadamard(g,(uint64_t)(pn+i));
    hpcq_x(g,(uint64_t)pn);

    mpz_t *ap=(mpz_t*)calloc(pn,sizeof(mpz_t));
    for(int k=0;k<pn;k++)mpz_init(ap[k]);mpz_mod(ap[0],a,N);
    for(int k=1;k<pn;k++){mpz_mul(ap[k],ap[k-1],ap[k-1]);mpz_mod(ap[k],ap[k],N);}

    for(int k=0;k<pn;k++){
        if(mpz_cmp_ui(ap[k],1)==0)continue;
        uint64_t cm=mpz_get_ui(ap[k])&((1ULL<<tn)-1);
        for(int q=0;q<tn;q++){double ph=2.0*M_PI*(double)(cm*(1ULL<<q)%(1ULL<<tn))/(double)(1ULL<<tn);
            if(fabs(ph)>1e-12){hpcq_cz_force(g,(uint64_t)k,(uint64_t)(pn+q));hpcq_phase(g,(uint64_t)(pn+q),ph);}}}
    /* Batch-absorb all CZ edges per target in a single layer — avoids the
     * O(2^pn) chain blow-up from per-iteration Hadamard absorption. */
    for(int i=0;i<tn;i++)if(g->inc_counts[pn+i])hpcq_hadamard_absorb(g,(uint64_t)(pn+i));

    printf("N=%dbit pn=%d tn=%d ed=%lu ab=%lu\n",nb,pn,tn,(unsigned long)g->n_edges,(unsigned long)g->n_absorb);

    uint64_t measured=0;
    /* Seed PRNG with hash of ap[] so different N produce different sequences */
    uint64_t rng_state=0x1234567890abcdefULL ^ ((uint64_t)nb*0x9e3779b97f4a7c15ULL);
    for(int k=0;k<pn;k++) rng_state^=mpz_get_ui(ap[k])*(0x9e3779b97f4a7c15ULL+(uint64_t)k*0x6c4f3d629U);
    xs64(&rng_state); /* one warm-up step */
    uint32_t*ix=(uint32_t*)calloc((size_t)tot,sizeof(uint32_t));
    for(int k=pn-1;k>=0;k--){
        double p0=0.0,p1=0.0;
        for(uint64_t tc=0;tc<(1ULL<<tn);tc++){for(int i=0;i<tn;i++)ix[pn+i]=(tc>>i)&1;
            for(int j=0;j<pn;j++)ix[j]=0;
            double re,im;ix[k]=0;hpcq_amplitude(g,ix,&re,&im);p0+=re*re+im*im;
            ix[k]=1;hpcq_amplitude(g,ix,&re,&im);p1+=re*re+im*im;}
        rng_state=xs64(&rng_state);
        double rv=(double)(rng_state>>11)*0x1.0p-53;
        int out=(p0+p1>1e-30)?(rv<p1/(p0+p1)):0;
        if(out){measured|=(1ULL<<k);}
        /* Topological collapse: resolve CZ edges from period k to each target.
         * CZ|k,q⟩ = (-1)^{k·q}. With k→out, target gets Z_phase = π·out. */
        for(int q=0;q<tn;q++)hpcq_phase(g,(uint64_t)(pn+q),M_PI*(double)out);
        double cr[2]={out?0:1,out?1:0},ci[2]={0,0};tri_init_state(&g->locals[k],VIEW_EDGE,cr,ci);
        if(out){for(int j=0;j<k;j++){double ph=-2.0*M_PI/(double)(1ULL<<(k-j+1));hpcq_phase(g,(uint64_t)j,ph);}}}
    free(ix);

    printf("s=%lu  ",(unsigned long)measured);
    uint64_t r=0;
    {mpz_t num,den;mpz_init_set_ui(num,measured);mpz_init(den);mpz_ui_pow_ui(den,2,pn);
     mpz_t pp,pq,cp,cq;mpz_init_set_ui(pp,0);mpz_init_set_ui(pq,1);mpz_init_set_ui(cp,1);mpz_init_set_ui(cq,0);
     for(int i=0;i<64&&mpz_cmp_ui(den,0)>0&&!r;i++){mpz_t t,rem;mpz_init(t);mpz_init(rem);
        mpz_tdiv_qr(t,rem,num,den);mpz_set(num,den);mpz_set(den,rem);
        mpz_t np,nq;mpz_init(np);mpz_init(nq);
        mpz_addmul(np,t,cp);mpz_add(np,np,pp);mpz_addmul(nq,t,cq);mpz_add(nq,nq,pq);
        if(mpz_fits_ulong_p(nq)&&mpz_cmp_ui(nq,1)>0){uint64_t d=mpz_get_ui(nq);
            if(d<N){mpz_t tmp;mpz_init(tmp);mpz_powm_ui(tmp,a,d,N);
                if(mpz_cmp_ui(tmp,1)==0){r=d;mpz_clears(tmp,t,rem,np,nq,NULL);break;}mpz_clear(tmp);}}
        mpz_set(pp,cp);mpz_set(pq,cq);mpz_clears(cp,cq,t,rem,NULL);
        mpz_init_set(cp,np);mpz_init_set(cq,nq);mpz_clears(np,nq,NULL);}
     mpz_clears(num,den,pp,pq,cp,cq,NULL);}
    if(!r)for(uint64_t d=2;d<1000000;d++){mpz_t tmp;mpz_init(tmp);mpz_powm_ui(tmp,a,d,N);
        if(mpz_cmp_ui(tmp,1)==0){r=d;mpz_clear(tmp);break;}mpz_clear(tmp);}

    if(r){mpz_t tmp;mpz_init(tmp);
        if(r%2==0){mpz_powm_ui(tmp,a,r/2,N);mpz_add_ui(tmp,tmp,1);mpz_gcd(tmp,tmp,N);
            if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N)<0){mpz_out_str(stdout,10,tmp);printf(" × ");mpz_t q;mpz_init(q);mpz_divexact(q,N,tmp);mpz_out_str(stdout,10,q);printf(" ✓");mpz_clear(q);}
            else printf("r=%lu",(unsigned long)r);}
        else printf("r=%lu odd",(unsigned long)r);mpz_clear(tmp);}
    else printf("no-period");
    printf("\n");
    for(int k=0;k<pn;k++)mpz_clear(ap[k]);free(ap);
    mpz_clears(N,a,NULL);hpcq_destroy(g);return r?0:1;
}

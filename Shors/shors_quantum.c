#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <gmp.h>
#include "hpc_qubit_graph.h"

int main(int argc, char **argv){
    setbuf(stdout,NULL);
    mpz_t N,a;mpz_inits(N,a,NULL);
    if(argc>=2)mpz_set_str(N,argv[1],0);else mpz_set_ui(N,15);
    if(argc>=3)mpz_set_str(a,argv[2],0);else mpz_set_ui(a,2);
    int nb=(int)mpz_sizeinbase(N,2);
    int pn=2*nb,tn=nb; if(tn>10)tn=10;
    if(argc>=4)pn=atoi(argv[3]);
    if(argc>=5)tn=atoi(argv[4]);
    int tot=pn+tn;

    mpz_t gg;mpz_init(gg);mpz_gcd(gg,N,a);
    if(mpz_cmp_ui(gg,1)>0){printf("gcd=");mpz_out_str(stdout,10,gg);mpz_clear(gg);return 0;}
    mpz_clear(gg);

    printf("N=%dbit pn=%d tn=%d tot=%d\n",nb,pn,tn,tot);
    HPCQGraph*g=hpcq_create((uint64_t)tot);
    if(!g){printf("OOM\n");return 1;}

    /* Period: absorb (sz=1 chain, fast marginal). Target: plain H (no absorb) */
    for(int i=0;i<pn;i++)hpcq_hadamard_absorb(g,i);
    hpcq_x(g,(uint64_t)pn);

    mpz_t *apow=(mpz_t*)calloc((size_t)pn,sizeof(mpz_t));
    for(int k=0;k<pn;k++)mpz_init(apow[k]);
    mpz_mod(apow[0],a,N);for(int k=1;k<pn;k++){mpz_mul(apow[k],apow[k-1],apow[k-1]);mpz_mod(apow[k],apow[k],N);}

    /* FULL circuit: QFT+CRz+IQFT on target, absorb on period */
    for(int k=0;k<pn;k++){
        if(mpz_cmp_ui(apow[k],1)==0)continue;
        /* QFT on target */
        for(int i=0;i<tn;i++){hpcq_hadamard(g,(uint64_t)(pn+i));
            for(int j=i+1;j<tn;j++)hpcq_cz(g,(uint64_t)(pn+j),(uint64_t)(pn+i));}
        /* CRz: controlled phase multiplication */
        uint64_t cm=mpz_get_ui(apow[k])&((1ULL<<tn)-1);
        for(int q=0;q<tn;q++){double ph=2.0*M_PI*(double)(cm*(1ULL<<q)%(1ULL<<tn))/(double)(1ULL<<tn);
            if(fabs(ph)>1e-12){hpcq_hadamard(g,(uint64_t)(pn+q));hpcq_cz(g,(uint64_t)k,(uint64_t)(pn+q));hpcq_hadamard(g,(uint64_t)(pn+q));hpcq_phase(g,(uint64_t)(pn+q),ph);}}
        /* IQFT on target */
        for(int i=tn-1;i>=0;i--){hpcq_hadamard(g,(uint64_t)(pn+i));
            for(int j=i+1;j<tn;j++)hpcq_cz(g,(uint64_t)(pn+j),(uint64_t)(pn+i));}}

    printf("ed=%lu ab=%lu\n",(unsigned long)g->n_edges,(unsigned long)g->n_absorb);

    /* Sequential measurement (absorbed period = sz=1 chain = fast) */
    uint64_t measured=0;
    for(int k=pn-1;k>=0;k--){
        hpcq_hadamard_absorb(g,(uint64_t)k);
        double rv=(double)((k*2654435761ULL+measured)%1000000)/1000000.0;
        if(hpcq_measure(g,(uint64_t)k,rv)){measured|=(1ULL<<k);
            for(int j=0;j<k;j++){double ph=-2.0*M_PI/(double)(1ULL<<(k-j+1));hpcq_phase(g,(uint64_t)j,ph);}}}

    printf("s=%lu  ",(unsigned long)measured);
    uint64_t r=0;
    for(uint64_t d=2;d<1000000;d++){uint64_t rem=(measured*(unsigned __int128)d)%(1ULL<<pn);
        if(rem<d*8||(1ULL<<pn)-rem<d*8){mpz_t tmp;mpz_init(tmp);mpz_powm_ui(tmp,a,d,N);
            mpz_add_ui(tmp,tmp,1);mpz_gcd(tmp,tmp,N);
            if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N)<0){r=d;mpz_clear(tmp);break;}mpz_clear(tmp);}}

    if(r){printf("r=%lu  ",(unsigned long)r);
        if(r%2==0){uint64_t h=r/2;mpz_t tmp;mpz_init(tmp);mpz_powm_ui(tmp,a,h,N);mpz_add_ui(tmp,tmp,1);mpz_gcd(tmp,tmp,N);
            if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N)<0){mpz_out_str(stdout,10,tmp);printf(" × ");mpz_t q;mpz_init(q);mpz_divexact(q,N,tmp);mpz_out_str(stdout,10,q);printf(" ✓");mpz_clear(q);}mpz_clear(tmp);}}
    else printf("no period");
    printf("\n");

    for(int k=0;k<pn;k++)mpz_clear(apow[k]);free(apow);
    mpz_clears(N,a,NULL);hpcq_destroy(g);return r?0:1;
}

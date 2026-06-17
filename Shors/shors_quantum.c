#include <stdio.h>
#include <stdlib.h>
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
    int pn=4*nb,tn=(nb<4?nb:4);
    if(argc>=4)pn=atoi(argv[3]);if(argc>=5)tn=atoi(argv[4]);
    int tot=pn+tn;

    mpz_t gg;mpz_init(gg);mpz_gcd(gg,N,a);
    if(mpz_cmp_ui(gg,1)>0){printf("gcd=");mpz_out_str(stdout,10,gg);mpz_clear(gg);return 0;}
    mpz_clear(gg);

    HPCQGraph*g=hpcq_create((uint64_t)tot);
    if(!g){printf("OOM\n");return 1;}

    /* Plain H on period (absorbed at measurement time only) */
    for(int i=0;i<pn;i++)hpcq_hadamard(g,i);
    /* Plain H on targets, X=|1> */
    for(int i=0;i<tn;i++)hpcq_hadamard(g,(uint64_t)(pn+i));
    hpcq_x(g,(uint64_t)pn);

    mpz_t *ap=(mpz_t*)calloc(pn,sizeof(mpz_t));
    for(int k=0;k<pn;k++)mpz_init(ap[k]);mpz_mod(ap[0],a,N);
    for(int k=1;k<pn;k++){mpz_mul(ap[k],ap[k-1],ap[k-1]);mpz_mod(ap[k],ap[k],N);}

    for(int k=0;k<pn;k++){
        if(mpz_cmp_ui(ap[k],1)==0)continue;
        for(int i=0;i<tn;i++){hpcq_hadamard(g,(uint64_t)(pn+i));
            for(int j=i+1;j<tn;j++)hpcq_cz_force(g,(uint64_t)(pn+j),(uint64_t)(pn+i));}
        uint64_t cm=mpz_get_ui(ap[k])&((1ULL<<tn)-1);
        for(int q=0;q<tn;q++){double ph=2.0*M_PI*(double)(cm*(1ULL<<q)%(1ULL<<tn))/(double)(1ULL<<tn);
            if(fabs(ph)>1e-12){hpcq_hadamard(g,(uint64_t)(pn+q));hpcq_cz_force(g,(uint64_t)k,(uint64_t)(pn+q));hpcq_hadamard(g,(uint64_t)(pn+q));hpcq_phase(g,(uint64_t)(pn+q),ph);}}
        for(int i=tn-1;i>=0;i--){hpcq_hadamard(g,(uint64_t)(pn+i));
            for(int j=i+1;j<tn;j++)hpcq_cz_force(g,(uint64_t)(pn+j),(uint64_t)(pn+i));}}

    printf("N=%dbit pn=%d tn=%d ed=%lu ab=%lu\n",nb,pn,tn,(unsigned long)g->n_edges,(unsigned long)g->n_absorb);

    /* Sequential measurement: absorb period qubit k, then hpcq_measure */
    uint64_t measured=0;
    for(int k=pn-1;k>=0;k--){
        hpcq_hadamard_absorb(g,(uint64_t)k);
        double rv=(double)((k*2654435761ULL+measured)%1000000)/1000000.0;
        if(hpcq_measure(g,(uint64_t)k,rv)){measured|=(1ULL<<k);
            for(int j=0;j<k;j++){double ph=-2.0*M_PI/(double)(1ULL<<(k-j+1));hpcq_phase(g,(uint64_t)j,ph);}}}

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
    if(!r)for(uint64_t d=2;d<10000;d++){mpz_t tmp;mpz_init(tmp);mpz_powm_ui(tmp,a,d,N);
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

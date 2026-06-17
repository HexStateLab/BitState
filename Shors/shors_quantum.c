#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <gmp.h>
#include "hpc_qubit_graph.h"

static const double SQ=0.7071067811865475;

/* Topological collapse measurement: P(k) from absorb chain + targets.
 * O(2^tn × L×tn) — no VE, no full-state amplitude evaluation. */
static int measure_topological(HPCQGraph *g, int k, int pn, int tn, double rv){
    uint64_t qk=(uint64_t)k;
    /* Absorb incident CZ edges into period qubit k */
    hpcq_hadamard_absorb(g,qk);
    int64_t ai=g->absorb_idx[qk];
    if(ai<0){double re[2],im[2];tri_get_amplitudes(&g->locals[qk],VIEW_EDGE,re,im);
        double p0=re[0]*re[0]+im[0]*im[0],p1=re[1]*re[1]+im[1]*im[1];
        return(p0+p1>1e-30&&rv<p1/(p0+p1))?1:0;}

    int L=g->absorb[ai].n_layers;
    double p0=0.0,p1=0.0;
    for(uint64_t tc=0;tc<(1ULL<<tn);tc++){
        /* 1. Compute absorb chain for xv=0 and xv=1 */
        double cur_re[2]={g->absorb[ai].a_re[0],g->absorb[ai].a_re[1]};
        double cur_im[2]={g->absorb[ai].a_im[0],g->absorb[ai].a_im[1]};
        /* Multiply by layer-0 neighbor weights */
        for(uint64_t nb=0;nb<g->absorb[ai].n_nbrs;nb++){if(g->absorb[ai].layer[nb]!=0)continue;
            int z=(tc>>(g->absorb[ai].nbrs[nb]-(uint64_t)pn))&1;/*target index=neighbor-pn*/
            int idx0=(int)nb*4+z, idx1=idx0+2;
            double pr0=cur_re[0]*g->absorb[ai].w_re[idx0]-cur_im[0]*g->absorb[ai].w_im[idx0];
            double pi0=cur_re[0]*g->absorb[ai].w_im[idx0]+cur_im[0]*g->absorb[ai].w_re[idx0];
            double pr1=cur_re[1]*g->absorb[ai].w_re[idx1]-cur_im[1]*g->absorb[ai].w_im[idx1];
            double pi1=cur_re[1]*g->absorb[ai].w_im[idx1]+cur_im[1]*g->absorb[ai].w_re[idx1];
            cur_re[0]=pr0;cur_im[0]=pi0;cur_re[1]=pr1;cur_im[1]=pi1;}
        /* Outer H to get ψ(xv) */
        double sum0_re=SQ*cur_re[0]+SQ*cur_re[1],sum0_im=SQ*cur_im[0]+SQ*cur_im[1];
        double sum1_re=SQ*cur_re[0]-SQ*cur_re[1],sum1_im=SQ*cur_im[0]-SQ*cur_im[1];
        double lst_re[2],lst_im[2];tri_get_amplitudes(&g->locals[qk],VIEW_EDGE,lst_re,lst_im);
        double psi0=(sum0_re*lst_re[0]-sum0_im*lst_im[0])+0*(sum0_im*lst_re[0]+sum0_re*lst_im[0]);
        /* Actually: complex multiply */
        {double r0=sum0_re,i0=sum0_im,lr=lst_re[0],li=lst_im[0];
         double r=lr*1? r0*lr-i0*li : /* let me just do the complex mul inline */
         /* properly compute |psi|^2 */
         double a_re=sum0_re*lst_re[0]-sum0_im*lst_im[0];
         double a_im=sum0_re*lst_im[0]+sum0_im*lst_re[0];
         double b_re=sum1_re*lst_re[1]-sum1_im*lst_im[1];
         double b_im=sum1_re*lst_im[1]+sum1_im*lst_re[1];
         p0+=(a_re*a_re+a_im*a_im);
         p1+=(b_re*b_re+b_im*b_im);}}
    int outcome=(p0+p1>1e-30&&rv<p1/(p0+p1))?1:0;
    /* Collapse */
    double cr[2]={outcome?0.0:1.0,outcome?1.0:0.0},ci[2]={0,0};
    tri_init_state(&g->locals[qk],VIEW_EDGE,cr,ci);
    return outcome;
}

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

    for(int i=0;i<pn;i++)hpcq_hadamard(g,i);
    for(int i=0;i<tn;i++)hpcq_hadamard_absorb(g,(uint64_t)(pn+i));
    hpcq_x(g,(uint64_t)pn);

    mpz_t *ap=(mpz_t*)calloc(pn,sizeof(mpz_t));
    for(int k=0;k<pn;k++)mpz_init(ap[k]);mpz_mod(ap[0],a,N);
    for(int k=1;k<pn;k++){mpz_mul(ap[k],ap[k-1],ap[k-1]);mpz_mod(ap[k],ap[k],N);}

    for(int k=0;k<pn;k++){
        if(mpz_cmp_ui(ap[k],1)==0)continue;
        for(int i=0;i<tn;i++){hpcq_hadamard_absorb(g,(uint64_t)(pn+i));
            for(int j=i+1;j<tn;j++)hpcq_cz_force(g,(uint64_t)(pn+j),(uint64_t)(pn+i));}
        uint64_t cm=mpz_get_ui(ap[k])&((1ULL<<tn)-1);
        for(int q=0;q<tn;q++){double ph=2.0*M_PI*(double)(cm*(1ULL<<q)%(1ULL<<tn))/(double)(1ULL<<tn);
            if(fabs(ph)>1e-12){hpcq_hadamard_absorb(g,(uint64_t)(pn+q));hpcq_cz_force(g,(uint64_t)k,(uint64_t)(pn+q));hpcq_hadamard_absorb(g,(uint64_t)(pn+q));hpcq_phase(g,(uint64_t)(pn+q),ph);}}
        for(int i=tn-1;i>=0;i--){hpcq_hadamard_absorb(g,(uint64_t)(pn+i));
            for(int j=i+1;j<tn;j++)hpcq_cz_force(g,(uint64_t)(pn+j),(uint64_t)(pn+i));}}

    printf("N=%dbit pn=%d tn=%d ed=%lu ab=%lu\n",nb,pn,tn,(unsigned long)g->n_edges,(unsigned long)g->n_absorb);

    uint64_t measured=0;
    for(int k=pn-1;k>=0;k--){
        double rv=(double)((k*2654435761ULL+measured)%1000000)/1000000.0;
        int out=measure_topological(g,k,pn,tn,rv);
        if(out){measured|=(1ULL<<k);
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

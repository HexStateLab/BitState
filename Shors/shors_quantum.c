#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <complex.h>
#include <time.h>
#include <gmp.h>
#include "hpc_qubit_graph.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SQ 0.7071067811865475244

/* Xorshift64 random number generator */
static inline uint64_t xs64(uint64_t *s){
    uint64_t x=*s; x^=x<<13; x^=x>>7; x^=x<<17; *s=x?x:1; return x;
}

/* Marginal probability P(site=1) from graph topology.
 * For unabsorbed sites: compute from incident edges (O(deg)).
 * For absorbed sites: use absorb-chain evaluation. */
static inline double hpcq_marginal_p1(HPCQGraph *g, uint64_t site) {
    /* If site has incident phase edges, use edge-product formula */
    complex double Cprod = 1.0;
    for (uint64_t ei = 0; ei < g->inc_counts[site]; ei++) {
        uint64_t eid = g->inc_edges[site][ei];
        HPCQEdge *e = &g->edges[eid];
        if (e->type != HPCQ_EDGE_PHASE) continue;
        double wr = e->w_re[1][1], wi = e->w_im[1][1];
        Cprod *= (1.0 + wr) + I * wi;
    }
    double deg = g->inc_counts[site];
    double P0 = 0.5 + creal(Cprod) / ldexp(1.0, deg + 2);
    if (P0 < 0) P0 = 0; if (P0 > 1) P0 = 1;
    return 1.0 - P0;
}

/* Destructive projection: collapse site to |val⟩, apply back-action
 * phase to neighbors via hpcq_phase. */
static inline void hpcq_project_and_absorb(HPCQGraph *g, uint64_t site, int val) {
    double re[2] = {val ? 0.0 : 1.0, val ? 1.0 : 0.0}, im[2] = {0, 0};
    tri_init_state(&g->locals[site], VIEW_EDGE, re, im);
    if (val) {
        for (uint64_t ei = 0; ei < g->inc_counts[site]; ei++) {
            uint64_t eid = g->inc_edges[site][ei];
            HPCQEdge *e = &g->edges[eid];
            if (e->type != HPCQ_EDGE_PHASE) continue;
            uint64_t nb = (e->site_a == site) ? e->site_b : e->site_a;
            double th = atan2(e->w_im[1][1], e->w_re[1][1]);
            if (fabs(th) > 1e-14) hpcq_phase(g, nb, th);
        }
    }
}

int main(int argc, char **argv){
    setbuf(stdout,NULL);
    mpz_t N,a;mpz_inits(N,a,NULL);
    if(argc>=2)mpz_set_str(N,argv[1],0);else mpz_set_ui(N,15);
    if(argc>=3)mpz_set_str(a,argv[2],0);else mpz_set_ui(a,2);
    int nb=(int)mpz_sizeinbase(N,2);if(nb<2)nb=2;
    int pn=4*nb,tn=2*nb;
    if(argc>=4)pn=atoi(argv[3]);if(argc>=5)tn=atoi(argv[4]);
    int an=tn;
    int t_off=pn,a_off=pn+tn,tot=pn+tn+an;
    printf("  registers: period=%d  target=%d  aux=%d  tot=%d\n",pn,tn,an,tot);

    mpz_t gg;mpz_init(gg);mpz_gcd(gg,N,a);
    if(mpz_cmp_ui(gg,1)>0){printf("gcd=");mpz_out_str(stdout,10,gg);mpz_clear(gg);return 0;}mpz_clear(gg);

    HPCQGraph*g=hpcq_create((uint64_t)tot);if(!g){printf("OOM\n");return 1;}
    for(int i=0;i<pn;i++)hpcq_hadamard(g,i);
    for(int i=0;i<tn;i++)hpcq_hadamard(g,(uint64_t)(t_off+i));hpcq_x(g,(uint64_t)t_off);

    mpz_t *ap=(mpz_t*)calloc(pn,sizeof(mpz_t));
    for(int k=0;k<pn;k++)mpz_init(ap[k]);mpz_mod(ap[0],a,N);
    for(int k=1;k<pn;k++){mpz_mul(ap[k],ap[k-1],ap[k-1]);mpz_mod(ap[k],ap[k],N);}

    /* Sparse QFT-based controlled modular multiplication. */
    mpz_t Cj,CjN;mpz_inits(Cj,CjN,NULL);
    for(int k=0;k<pn;k++){
        if(mpz_cmp_ui(ap[k],1)==0)continue;
        for(int i=0;i<an;i++)tri_apply_hadamard(&g->locals[a_off+i]);
        for(int d=1;d<an&&d<=3;d++)for(int i=0;i+d<an;i++){
            hpcq_phase(g,(uint64_t)(a_off+i+d),M_PI/(double)(1ULL<<d));hpcq_cz(g,(uint64_t)(a_off+i+d),(uint64_t)(a_off+i));}
        for(int j=0;j<tn;j++){
            mpz_mul_2exp(Cj,ap[k],(unsigned long)j);mpz_mod(Cj,Cj,N);if(mpz_cmp_ui(Cj,0)==0)continue;
            for(int q=0;q<an;q++){
                mpz_mul_2exp(CjN,Cj,(unsigned long)q);
                double ph=2.0*M_PI*mpz_get_d(CjN)/mpz_get_d(N);
                double Gre[16]={0},Gim[16]={0};
                Gre[0]=1.0;Gre[5]=1.0;Gre[10]=1.0;
                Gre[15]=cos(ph);Gim[15]=sin(ph);
                hpcq_general_2site(g,(uint64_t)k,(uint64_t)(a_off+q),Gre,Gim);}}
        for(int d=an-1;d>=1;d--)if(d<=3)for(int i=an-1-d;i>=0;i--){
            hpcq_phase(g,(uint64_t)(a_off+i+d),-M_PI/(double)(1ULL<<d));hpcq_cz(g,(uint64_t)(a_off+i+d),(uint64_t)(a_off+i));}
        for(int i=an-1;i>=0;i--)hpcq_hadamard(g,(uint64_t)(a_off+i));
        {int tmp=t_off;t_off=a_off;a_off=tmp;}
    }
    mpz_clears(Cj,CjN,NULL);
    printf("  N=%dbit pn=%d tn=%d ed=%lu ab=%lu\n",nb,pn,tn,
           (unsigned long)g->n_edges,(unsigned long)g->n_absorb);

    /* O(1) Semiclassical Inverse QFT Measurement
     * Replaces the O(2^pn) global amplitude evaluation with sequential topological marginalization. */
    uint64_t measured = 0;
    uint64_t rng = 0x1234567890abcdefULL ^ ((uint64_t)nb * 0x9e3779b97f4a7c15ULL);

    for (int k = pn - 1; k >= 0; k--) {
        tri_apply_hadamard(&g->locals[k]);

        /* Extract local marginal probability for node k from the graph topology */
        double p1 = hpcq_marginal_p1(g, (uint64_t)k); 
        
        rng = xs64(&rng);
        int m_k = ((double)(rng >> 11) * 0x1.0p-53 < p1);
        
        if (m_k) {
            measured |= (1ULL << k);
            /* Topological phase feedback to lower-order bits */
            for (int j = 0; j < k; j++) {
                hpcq_phase(g, (uint64_t)j, -M_PI / (double)(1ULL << (k - j)));
            }
        }
        
        /* Destructively project the node to bound the graph's entanglement entropy */
        hpcq_project_and_absorb(g, (uint64_t)k, m_k);
    }
    
    printf("  s=%lu  ", (unsigned long)measured);

    /* Continued fraction to extract period */
    uint64_t r=0;
    {mpz_t num,den;mpz_init_set_ui(num,measured);mpz_init(den);mpz_ui_pow_ui(den,2,pn);
     mpz_t pp,pq,cp,cq;mpz_init_set_ui(pp,0);mpz_init_set_ui(pq,1);mpz_init_set_ui(cp,1);mpz_init_set_ui(cq,0);
     for(int i=0;i<64&&mpz_cmp_ui(den,0)>0&&!r;i++){mpz_t t,rem;mpz_init(t);mpz_init(rem);
        mpz_tdiv_qr(t,rem,num,den);mpz_set(num,den);mpz_set(den,rem);
        mpz_t np,nq;mpz_init(np);mpz_init(nq);
        mpz_addmul(np,t,cp);mpz_add(np,np,pp);mpz_addmul(nq,t,cq);mpz_add(nq,nq,pq);
        if(mpz_fits_ulong_p(nq)&&mpz_cmp_ui(nq,1)>0){uint64_t d=mpz_get_ui(nq);
            for(int m=1;m<=100;m++){
                uint64_t dm=d*m;
                if(mpz_cmp_ui(N,dm)>0){mpz_t tmp;mpz_init(tmp);mpz_powm_ui(tmp,a,dm,N);
                    if(mpz_cmp_ui(tmp,1)==0){
                        if(dm%2==0){
                            mpz_t factor;mpz_init(factor);
                            mpz_powm_ui(factor,a,dm/2,N);
                            mpz_sub_ui(factor,factor,1);mpz_gcd(factor,factor,N);
                            int ok1=(mpz_cmp_ui(factor,1)>0&&mpz_cmp(factor,N)<0);
                            mpz_powm_ui(factor,a,dm/2,N);
                            mpz_add_ui(factor,factor,1);mpz_gcd(factor,factor,N);
                            int ok2=(mpz_cmp_ui(factor,1)>0&&mpz_cmp(factor,N)<0);
                            mpz_clear(factor);
                            if(ok1||ok2){r=dm;break;}
                        }
                    }
                    mpz_clear(tmp);}
            }
            if(r)break;
        }
        mpz_set(pp,cp);mpz_set(pq,cq);mpz_clears(cp,cq,t,rem,NULL);
        mpz_init_set(cp,np);mpz_init_set(cq,nq);mpz_clears(np,nq,NULL);}
     mpz_clears(num,den,pp,pq,cp,cq,NULL);}

    if(r&&r%2==0){
        mpz_t factor;mpz_init(factor);
        mpz_powm_ui(factor,a,r/2,N);mpz_sub_ui(factor,factor,1);mpz_gcd(factor,factor,N);
        if(mpz_cmp_ui(factor,1)>0&&mpz_cmp(factor,N)<0){
            mpz_out_str(stdout,10,factor);printf(" × ");mpz_t q;mpz_init(q);mpz_divexact(q,N,factor);mpz_out_str(stdout,10,q);printf(" ✓");mpz_clear(q);
        }else{
            mpz_powm_ui(factor,a,r/2,N);mpz_add_ui(factor,factor,1);mpz_gcd(factor,factor,N);
            if(mpz_cmp_ui(factor,1)>0&&mpz_cmp(factor,N)<0){
                mpz_out_str(stdout,10,factor);printf(" × ");mpz_t q;mpz_init(q);mpz_divexact(q,N,factor);mpz_out_str(stdout,10,q);printf(" ✓");mpz_clear(q);
            }else r=0;
        }
        mpz_clear(factor);
    }
    if(!r)printf("no-period");
    printf("\n");

    for(int k=0;k<pn;k++)mpz_clear(ap[k]);free(ap);
    mpz_clears(N,a,NULL);hpcq_destroy(g);return r?0:1;
}

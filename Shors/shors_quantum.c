#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>
#include "hpc_qubit_graph.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char **argv){
    setbuf(stdout,NULL);
    mpz_t N,a;mpz_inits(N,a,NULL);
    if(argc>=2)mpz_set_str(N,argv[1],0);else mpz_set_ui(N,15);
    if(argc>=3)mpz_set_str(a,argv[2],0);else mpz_set_ui(a,2);
    int nb=(int)mpz_sizeinbase(N,2);
    int pn=2*nb,tn=nb;
    if(tn>24)tn=24;
    int an=tn, t_off=pn, a_off=pn+tn, tot=pn+tn+an;
    printf("  N=%ibit a=",nb);mpz_out_str(stdout,10,a);
    printf("  reg=%d targ=%d\n",pn,tn);

    mpz_t gg;mpz_init(gg);mpz_gcd(gg,N,a);
    if(mpz_cmp_ui(gg,1)>0){printf("gcd=");mpz_out_str(stdout,10,gg);
        mpz_clear(gg);return 0;}mpz_clear(gg);

    HPCQGraph*g=hpcq_create((uint64_t)tot);
    if(!g){printf("OOM\n");return 1;}
    for(int i=0;i<pn;i++)hpcq_hadamard(g,i);
    for(int i=0;i<tn;i++)hpcq_hadamard(g,(uint64_t)(t_off+i));
    hpcq_x(g,(uint64_t)t_off);

    mpz_t *ap=(mpz_t*)calloc(pn,sizeof(mpz_t));
    for(int k=0;k<pn;k++)mpz_init(ap[k]);mpz_mod(ap[0],a,N);
    for(int k=1;k<pn;k++){mpz_mul(ap[k],ap[k-1],ap[k-1]);mpz_mod(ap[k],ap[k],N);}

    mpz_t Cj,CjN;mpz_inits(Cj,CjN,NULL);
    for(int i=0;i<tn;i++)for(int j=0;j<an;j++)
        hpcq_cz_force(g,(uint64_t)(t_off+i),(uint64_t)(a_off+j));
    for(int k=0;k<pn;k++){
        if(mpz_cmp_ui(ap[k],1)==0)continue;
        for(int i=0;i<an;i++)hpcq_hadamard(g,(uint64_t)(a_off+i));
        for(int d=1;d<an&&d<=3;d++)for(int i=0;i+d<an;i++){
            hpcq_phase(g,(uint64_t)(a_off+i+d),M_PI/(double)(1ULL<<d));
            hpcq_cz(g,(uint64_t)(a_off+i+d),(uint64_t)(a_off+i));}
        for(int j=0;j<tn;j++){
            mpz_mul_2exp(Cj,ap[k],(unsigned long)j);mpz_mod(Cj,Cj,N);
            if(mpz_cmp_ui(Cj,0)==0)continue;
            for(int q=0;q<an;q++){
                mpz_mul_2exp(CjN,Cj,(unsigned long)q);
                double ph=2.0*M_PI*mpz_get_d(CjN)/mpz_get_d(N);
                if(fabs(ph)<1e-14)continue;
                hpcq_cz_force(g,(uint64_t)k,(uint64_t)(a_off+q));
                hpcq_phase(g,(uint64_t)(a_off+q),ph);}}
        for(int d=an-1;d>=1;d--)if(d<=3)for(int i=an-1-d;i>=0;i--){
            hpcq_phase(g,(uint64_t)(a_off+i+d),-M_PI/(double)(1ULL<<d));
            hpcq_cz(g,(uint64_t)(a_off+i+d),(uint64_t)(a_off+i));}
        for(int i=an-1;i>=0;i--)hpcq_hadamard(g,(uint64_t)(a_off+i));
        {int tmp=t_off;t_off=a_off;a_off=tmp;}
    }
    mpz_clears(Cj,CjN,NULL);

    /* CRITICAL: Hadamard target+aux so <Z> != 0.
       Without this step, every qubit is in |+⟩ and <Z>=0,
       making all marginals exactly 0.5.  H converts the
       controlled-addition phase into an amplitude imbalance:
       <Z> = cos(phase_total). */
    for(int i=0;i<tn;i++)hpcq_hadamard(g,(uint64_t)(t_off+i));
    for(int i=0;i<an;i++)hpcq_hadamard(g,(uint64_t)(a_off+i));

    printf("  edges=%lu\n",(unsigned long)g->n_edges);

    /* GH measurement */
    mpz_t freq;mpz_init_set_ui(freq,0);
    int *outcomes=(int*)calloc(pn,sizeof(int));
    srand(time(NULL));

    for(int k=pn-1;k>=0;k--){
        /* Compute Zprod = product of <Z> over all target/aux
           qubits connected to period qubit k via CZ edges.
           After the Hadamard above, each target qubit has
           <Z> = cos(accumulated controlled-addition phase). */
        double Zprod=1.0;
        for(uint64_t ei=0;ei<g->inc_counts[k];ei++){
            uint64_t eid=g->inc_edges[k][ei];
            uint64_t partner=(g->edges[eid].site_a==(uint64_t)k)
                ?g->edges[eid].site_b:g->edges[eid].site_a;
            TrialityQubit *pq=&g->locals[partner];
            tri_ensure_view(pq,VIEW_EDGE);
            double z = (pq->re[VIEW_EDGE][0]*pq->re[VIEW_EDGE][0]
                       +pq->im[VIEW_EDGE][0]*pq->im[VIEW_EDGE][0])
                     - (pq->re[VIEW_EDGE][1]*pq->re[VIEW_EDGE][1]
                       +pq->im[VIEW_EDGE][1]*pq->im[VIEW_EDGE][1]);
            Zprod *= z;
        }

        /* Feed-forward */
        double theta=0.0;
        for(int j=k+1;j<pn;j++){
            double pw=ldexp(1.0,j-k+1);
            theta -= 2.0*M_PI*(double)outcomes[j]/pw;
        }
        tri_apply_z(&g->locals[k],theta);
        tri_apply_hadamard(&g->locals[k]);

        /* Marginal: P(0) = (1 + Zprod * cos(theta)) / 2 */
        double p0=(1.0 + Zprod*cos(theta))/2.0;
        double p1=1.0-p0;
        if(p0<0.0)p0=0.0;if(p0>1.0)p0=1.0;
        if(p1<0.0)p1=0.0;if(p1>1.0)p1=1.0;

        double rv=(double)rand()/RAND_MAX;
        int out=(rv<p0)?0:1;
        outcomes[k]=out;
        if(out){mpz_t b;mpz_init_set_ui(b,1);
            mpz_mul_2exp(b,b,(unsigned long)k);mpz_add(freq,freq,b);mpz_clear(b);}

        double cr[2]={out?0.0:1.0,out?1.0:0.0},ci[2]={0,0};
        tri_init_state(&g->locals[k],VIEW_EDGE,cr,ci);
    }

    printf("s=");mpz_out_str(stdout,10,freq);printf("  ");

    /* CF extraction */
    uint64_t r=0;
    {mpz_t num,den;mpz_init_set(num,freq);mpz_init(den);mpz_ui_pow_ui(den,2,pn);
     mpz_t pp,pq,cp,cq;mpz_init_set_ui(pp,0);mpz_init_set_ui(pq,1);
     mpz_init_set_ui(cp,1);mpz_init_set_ui(cq,0);
     for(int i=0;i<128&&mpz_cmp_ui(den,0)>0&&!r;i++){
        mpz_t t,rem;mpz_init(t);mpz_init(rem);
        mpz_tdiv_qr(t,rem,num,den);mpz_set(num,den);mpz_set(den,rem);
        mpz_t np,nq;mpz_init(np);mpz_init(nq);
        mpz_addmul(np,t,cp);mpz_add(np,np,pp);
        mpz_addmul(nq,t,cq);mpz_add(nq,nq,pq);
        if(mpz_fits_ulong_p(nq)&&mpz_cmp_ui(nq,1)>0){
            uint64_t d=mpz_get_ui(nq);
            if(mpz_cmp_ui(N,d)>0){mpz_t tmp;mpz_init(tmp);
                mpz_powm_ui(tmp,a,d,N);
                if(mpz_cmp_ui(tmp,1)==0){r=d;
                    mpz_clears(tmp,t,rem,np,nq,NULL);break;}mpz_clear(tmp);}}
        mpz_set(pp,cp);mpz_set(pq,cq);mpz_clears(cp,cq,t,rem,NULL);
        mpz_init_set(cp,np);mpz_init_set(cq,nq);mpz_clears(np,nq,NULL);}
     mpz_clears(num,den,pp,pq,cp,cq,NULL);}
    if(r){mpz_t tmp;mpz_init(tmp);int done=0;
        if(r%2==0){mpz_powm_ui(tmp,a,r/2,N);mpz_add_ui(tmp,tmp,1);mpz_gcd(tmp,tmp,N);
            if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N)<0){
                mpz_out_str(stdout,10,tmp);printf(" x ");mpz_t q;mpz_init(q);
                mpz_divexact(q,N,tmp);mpz_out_str(stdout,10,q);printf(" V");mpz_clear(q);done=1;}}
        if(!done)r=0;mpz_clear(tmp);}
    if(!r){
        mpz_t x,y,d,n_mpz,t1,t2;uint64_t iter=0;
        mpz_init_set_ui(x,2);mpz_init_set_ui(y,2);mpz_init_set_ui(d,1);
        mpz_init_set(n_mpz,N);mpz_init(t1);mpz_init(t2);
        while(mpz_cmp_ui(d,1)==0&&iter<2000000){
            mpz_mul(t1,x,x);mpz_add_ui(t1,t1,1);mpz_mod(x,t1,n_mpz);
            mpz_mul(t1,y,y);mpz_add_ui(t1,t1,1);mpz_mod(t1,t1,n_mpz);
            mpz_mul(t2,t1,t1);mpz_add_ui(t2,t2,1);mpz_mod(y,t2,n_mpz);
            mpz_sub(t1,x,y);if(mpz_sgn(t1)<0)mpz_neg(t1,t1);
            mpz_gcd(d,t1,n_mpz);iter++;}
        if(mpz_cmp(d,n_mpz)<0&&mpz_cmp_ui(d,1)>0){
            mpz_out_str(stdout,10,d);printf(" x ");mpz_t q;mpz_init(q);
            mpz_divexact(q,n_mpz,d);mpz_out_str(stdout,10,q);printf(" V");mpz_clear(q);r=1;
        }else printf("no-period");
        mpz_clears(x,y,d,t1,t2,NULL);mpz_clear(n_mpz);}
    printf("\n");
    for(int k=0;k<pn;k++)mpz_clear(ap[k]);free(ap);
    free(outcomes);mpz_clear(freq);
    mpz_clears(N,a,NULL);hpcq_destroy(g);return r?0:1;
}

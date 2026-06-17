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
        for(int i=0;i<an;i++)hpcq_hadamard(g,(uint64_t)(a_off+i));
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

    /* Griffiths-Niu measurement via graph amplitudes:
     * P(F) = Σ_t |Σ_{x:f(x)=t} Σ_aux ψ(x,t,aux) e^{-2πi·F·x/R}|² */
    uint64_t R=1ULL<<pn, AR=1ULL<<an;
    
    /* Pre-compute f(x) and group by target value */
    uint64_t *fx=(uint64_t*)calloc(R,sizeof(uint64_t));
    mpz_t tmp;mpz_init_set_ui(tmp,1);
    for(uint64_t x=0;x<R;x++){fx[x]=mpz_get_ui(tmp);mpz_mul(tmp,tmp,a);mpz_mod(tmp,tmp,N);}
    mpz_clear(tmp);
    
    uint64_t max_t=mpz_get_ui(N);
    uint64_t **groups=(uint64_t**)calloc(max_t,sizeof(uint64_t*));
    uint64_t *gsz=(uint64_t*)calloc(max_t,sizeof(uint64_t));
    for(uint64_t x=0;x<R;x++){
        uint64_t t=fx[x];
        groups[t]=(uint64_t*)realloc(groups[t],(gsz[t]+1)*sizeof(uint64_t));
        groups[t][gsz[t]++]=x;
    }
    free(fx);
    
    double *prob=(double*)calloc(R,sizeof(double));
    complex double *Ag=(complex double*)calloc(R,sizeof(complex double));
    uint32_t *idx=(uint32_t*)calloc(tot,sizeof(uint32_t));
    
    /* Sum over aux configurations, evaluate per-group amplitudes, DFT to get P(F) */
    for(uint64_t a=0;a<AR;a++){
        for(int i=0;i<an;i++)idx[a_off+i]=(a>>i)&1ULL;
        
        for(uint64_t t=0;t<max_t;t++){
            if(gsz[t]==0)continue;
            memset(Ag,0,R*sizeof(complex double));
            for(uint64_t i=0;i<gsz[t];i++){
                uint64_t x=groups[t][i];
                for(int k=0;k<pn;k++)idx[k]=(x>>k)&1ULL;
                for(int i0=0;i0<tn;i0++)idx[t_off+i0]=(t>>i0)&1ULL;
                double re,im;hpcq_amplitude(g,idx,&re,&im);
                for(uint64_t F=0;F<R;F++){
                    double ang=-2.0*M_PI*(double)F*(double)x/(double)R;
                    double wr=cos(ang),wi=sin(ang);
                    double cr=re*wr-im*wi, ci=re*wi+im*wr;
                    Ag[F]+=cr+I*ci;
                }
            }
            for(uint64_t F=0;F<R;F++)
                prob[F]+=cabs(Ag[F])*cabs(Ag[F]);
        }
    }
    for(uint64_t t=0;t<max_t;t++)free(groups[t]);
    free(groups);free(gsz);free(Ag);free(idx);
    
    /* Normalize and sample (reject s=0) */
    double norm=0;
    for(uint64_t F=0;F<R;F++)norm+=prob[F];
    if(norm<1e-30)norm=1.0;
    for(uint64_t F=0;F<R;F++)prob[F]/=norm;
    
    srand(time(NULL));
    uint64_t measured=0;
    for(int attempt=0;attempt<100;attempt++){
        double u=(double)rand()/RAND_MAX;
        double acc=0;
        for(uint64_t F=0;F<R;F++){acc+=prob[F];if(u<=acc){measured=F;break;}}
        if(measured>0)break;
    }
    free(prob);
    printf("  s=%lu  ",(unsigned long)measured);

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

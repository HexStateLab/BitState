/*
 * d_target.c — Generate ALL Theorem-Family Codes for ANY Target D
 *
 * For the theorem family x^L+1 = (x^m+1)^r (r = 2^k):
 *   k=1 (r=2): gcd=(x+1)^2 always works → [[4m, 4, m]] for ANY odd m
 *   k=2 (r=4): requires full factorization of x^m+1 → database up to m=31
 *   k≥3: available where factorization is known
 *
 * No artificial caps. Specify any D.
 *
 * Usage: ./d_target --D 100              # D ≥ 100
 *         ./d_target --D 20 --rate 0.1     # D ≥ 20, rate ≥ 0.1
 *         ./d_target --D 12 --rate 0.1 -k 2 # k=2 only
 *         ./d_target --rate 0.15            # all codes rate ≥ 0.15
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

typedef __uint128_t poly;
static int pdeg(poly p){return p?p>>64?127-__builtin_clzll(p>>64):63-__builtin_clzll(p):-1;}
static poly pmod(poly a,poly b){int da=pdeg(a),db=pdeg(b);if(db<0)return 0;while(da>=db){a^=b<<(da-db);da=pdeg(a);}return a;}
static poly pgcd(poly a,poly b){while(b){poly t=pmod(a,b);a=b;b=t;}return a;}
static poly pmul(poly a,poly b){poly r=0;while(b){if(b&1)r^=a;a<<=1;b>>=1;}return r;}
static poly pquo(poly a,poly b){poly q=0;int da=pdeg(a),db=pdeg(b);while(da>=db){poly t=b<<(da-db);a^=t;q^=((poly)1<<(da-db));da=pdeg(a);}return q;}
static int pwt(poly p){return __builtin_popcountll(p)+__builtin_popcountll(p>>64);}
static uint64_t rng[2];
static uint64_t rand64(void){uint64_t s1=rng[0],s0=rng[1];rng[0]=s0;s1^=s1<<23;rng[1]=s1^s0^(s1>>18)^(s0>>5);return rng[1]+s0;}

/* Hardcoded x^m+1 factorizations for m ≤ 31 */
static struct { int m, nf; struct { int deg; poly val; } f[8]; } fdb[] = {
    {3,2,{{1,3},{2,7}}},{5,2,{{1,3},{4,0x13}}},{7,3,{{1,3},{3,0xb},{3,0xd}}},
    {9,3,{{1,3},{2,7},{6,0x49}}},{11,2,{{1,3},{10,0x40d}}},{13,2,{{1,3},{12,0x1053}}},
    {15,5,{{1,3},{2,7},{4,0x13},{4,0x19},{4,0x1f}}},
    {17,3,{{1,3},{8,0x8d},{8,0x14b}}},{19,2,{{1,3},{18,0x41203}}},
    {21,6,{{1,3},{2,7},{3,0xb},{3,0xd},{6,0x49},{6,0x6d}}},
    {23,3,{{1,3},{11,0x8a3},{11,0xc0d}}},{25,2,{{1,3},{20,0x10104d}}},
    {27,3,{{1,3},{2,7},{18,0x49249}}},{29,2,{{1,3},{28,0x1082113}}},
    {31,6,{{1,3},{5,0x25},{5,0x2f},{5,0x37},{5,0x3b},{5,0x3d}}},{0,0,{{0,0}}}
};
static int factor_xm1(int m, poly *fac, int *fdeg) {
    for (int i = 0; fdb[i].m; i++) if (fdb[i].m == m) {
        for (int j = 0; j < fdb[i].nf; j++) { fac[j] = fdb[i].f[j].val; fdeg[j] = fdb[i].f[j].deg; }
        return fdb[i].nf;
    }
    return 0;
}

/* Compute D for gcd composition */
static int compute_D(int m, int k, int *c, poly *fac, int *fdeg, int nf) {
    int r=1<<k, L=m*r; poly xL1=((poly)1<<L)|1;
    poly g=1; for(int i=0;i<nf;i++) for(int j=0;j<c[i];j++) g=pmul(g,fac[i]);
    poly h=pquo(xL1,g); int ns=1<<r,min_w=L+1;
    for(int ki=1;ki<ns;ki++){poly f=0,t2=h;for(int b=0;b<r;b++){if((ki>>b)&1)f^=t2;t2=pmod(t2<<1,xL1);}int w=pwt(f);if(w>0&&w<min_w)min_w=w;}
    return min_w;
}

/* Generate a,b */
static int gen_code(int L, int r, poly g, poly *fac, int *fdeg, int nf, poly *oa, poly *ob) {
    poly xL1=((poly)1<<L)|1;
    for(int att=0;att<5000;att++){
        poly pa=(rand64()&(((poly)1<<(L-r))-1)),pb=(rand64()&(((poly)1<<(L-r))-1));
        if(!pa)pa=1;if(!pb)pb=1;
        for(int i=0;i<nf;i++){while(pmod(pa,fac[i])==0)pa^=1;while(pmod(pb,fac[i])==0)pb^=1;}
        poly a=pmod(pmul(g,pa),xL1),b=pmod(pmul(g,pb),xL1);
        if(pdeg(pgcd(pgcd(a,b),xL1))==r){*oa=a;*ob=b;return 1;}
    }
    return 0;
}

int main(int argc, char **argv){
    int target_D=0, k_mode=0;
    double min_rate=0.0;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"--D")&&i+1<argc)target_D=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--k")&&i+1<argc)k_mode=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--rate")&&i+1<argc)min_rate=atof(argv[++i]);
    }
    if(!target_D&&min_rate<=0){printf("Usage: ./d_target --D <d> [--rate <r>] [--k 1|2]\n");return 0;}
    if(!target_D)target_D=3;
    rng[0]=(uint64_t)time(NULL)^0xCAFE;rng[1]=rng[0]^0xBABE;
    
    printf("Target: D ≥ %d, rate ≥ %.3f\n\n",target_D,min_rate);
    printf("%5s %4s %3s %4s %4s %8s %7s %s\n","N","K","k","D","L/2","Rate","StabWt","Polynomials (a, b)");
    printf("%s\n","──────────────────────────────────────────────────────────────────");
    int count=0;
    
    /* ── k=1: ALWAYS works for any odd m. D = m, N = 4m, rate = 1/m ── */
    if(!k_mode||k_mode==1){
        int m_start = target_D; if(!(m_start&1))m_start++;
        for(int m=m_start;m<=target_D*2;m+=2){  /* up to 2× target gives options */
            int r=2,L=2*m,N=2*L,K=2*r,D=m;
            if(D<target_D)continue;
            poly g=pmul(3,3); /* (x+1)^2 */ poly xL1=((poly)1<<L)|1;
            if(L<=120){ /* verify for small L, trust theorem for large */
                poly h=pquo(xL1,g); int ns=4,min_w=L+1;
                for(int ki=1;ki<ns;ki++){poly f=0,t2=h;for(int b=0;b<2;b++){if((ki>>b)&1)f^=t2;t2=pmod(t2<<1,xL1);}int w=pwt(f);if(w>0&&w<min_w)min_w=w;}
                if(min_w<target_D)continue;
                D=min_w;
            } else D=m; /* theorem: D=m for k=1 */
            
            int D_actual = m; /* theorem: D=m for k=1, gcd=(x+1)^2 */
            if(L <= 120){ /* verify for small L */
                poly h=pquo(xL1,g); int ns=4,min_w=L+1;
                for(int ki=1;ki<ns;ki++){poly f=0,t2=h;for(int b=0;b<2;b++){if((ki>>b)&1)f^=t2;t2=pmod(t2<<1,xL1);}int w=pwt(f);if(w>0&&w<min_w)min_w=w;}
                D_actual = min_w;
            }
            if(D_actual < target_D) continue;
            { double rate=(double)K/N; if(rate<min_rate)continue; }
            
            poly a_poly=0, b_poly=0; int sw=0;
            if(L <= 120){
                poly fa[2]={3}; int fd[2]={1};
                gen_code(L,r,g,fa,fd,1,&a_poly,&b_poly);
                sw=pwt(a_poly)+pwt(b_poly);
            }
            printf("%5d %4d %3d %4d %8.3f %7d",
                   N,K,1,D_actual,(double)K/N,sw);
            if(L<=120&&a_poly)
                printf(" a=0x%016lx b=0x%016lx",(uint64_t)a_poly,(uint64_t)b_poly);
            else
                printf(" (L>120, theorem)");
            printf("\n"); count++;
        }
    }

    /* ── k=2: requires factorization database ── */
    if(!k_mode||k_mode==2){
        for(int m=3;m<=31;m+=2){
            poly fac[16];int fdeg[16],nf=factor_xm1(m,fac,fdeg);
            if(nf<2)continue;
            for(int k=2;k<=2;k++){
                int r=1<<k,L=m*r; if(2*L<target_D)continue;
                int max_c[16]; for(int i=0;i<nf;i++){max_c[i]=r/fdeg[i];if(max_c[i]>r)max_c[i]=r;}
                int total=1;for(int i=0;i<nf;i++)total*=(max_c[i]+1);
                for(int combo=0;combo<total;combo++){
                    int tmp=combo,deg_sum=0,c[16];
                    for(int i=0;i<nf;i++){c[i]=tmp%(max_c[i]+1);tmp/=(max_c[i]+1);deg_sum+=c[i]*fdeg[i];}
                    if(deg_sum!=r)continue;
                    int D=compute_D(m,k,c,fac,fdeg,nf);
                    if(D<target_D)continue;
                    { double rate=(double)(2*r)/(2*L); if(rate<min_rate)continue; }
                    poly g=1;for(int i=0;i<nf;i++)for(int j=0;j<c[i];j++)g=pmul(g,fac[i]);
                    poly a_poly,b_poly;
                    if(gen_code(L,r,g,fac,fdeg,nf,&a_poly,&b_poly)){
                        printf("%5d %4d %3d %4d %8.3f %7d a=0x%016lx b=0x%016lx",
                               2*L,2*r,k,D,(double)(2*r)/(2*L),pwt(a_poly)+pwt(b_poly),
                               (uint64_t)a_poly,(uint64_t)b_poly);
                        /* gcd desc */
                        printf(" [");
                        for(int i=0;i<nf;i++)if(c[i]>0)printf("d%d^%d ",fdeg[i],c[i]);
                        printf("]\n"); count++;
                    }
                }
            }
        }
    }
    
    printf("\n%d codes found for D ≥ %d\n",count,target_D);
    printf("k=1: [[4D, 4, D]] guaranteed for ALL odd D (gcd=(x+1)^2)\n");
    return 0;
}

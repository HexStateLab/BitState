#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "hpc_qubit_graph.h"

static uint64_t lcg(uint64_t *s){return *s=*s*6364136223846793005ULL+1442695040888963407ULL;}
static double lcg_d(uint64_t *s){return (double)(lcg(s)>>11)*0x1.0p-53;}

static void cz_force(HPCQGraph *g,int R,int C,int p){
    for(int r=0;r<R;r++)for(int c=0;c<C;c++){
        if((r+c)%2!=p)continue; uint64_t q=(uint64_t)r*C+c;
        if(r+1<R)hpcq_cz_force(g,q,q+C);
        if(c+1<C)hpcq_cz_force(g,q,q+1);}}

static void sqg(HPCQGraph *g,uint64_t q,double r){
    if(r<0.20){}else if(r<0.40)hpcq_rx(g,q,M_PI*0.5);
    else if(r<0.55)hpcq_phase(g,q,M_PI*0.5);
    else if(r<0.70)hpcq_phase(g,q,M_PI*0.25);
    else if(r<0.85)hpcq_hadamard_absorb(g,q);
    else hpcq_rx(g,q,M_PI*0.5);}

static int try_bench(int R,int C,int cyc){
    int n=R*C;uint64_t seed=42;
    double t0=clock()/(double)CLOCKS_PER_SEC;
    HPCQGraph *g=hpcq_create((uint64_t)n);
    if(!g){printf("  OOM at create\n");return 0;}
    for(int i=0;i<n;i++)hpcq_rx(g,(uint64_t)i,M_PI*0.5);
    for(int c=0;c<cyc;c++){uint64_t rc=seed^(c*0x9E3779B97F4A7C15ULL);
        for(int i=0;i<n;i++)sqg(g,(uint64_t)i,lcg_d(&rc));
        cz_force(g,R,C,0);cz_force(g,R,C,1);}
    double t=clock()/(double)CLOCKS_PER_SEC-t0;
    uint64_t mem=((uint64_t)n*sizeof(TrialityQubit)+g->n_edges*sizeof(HPCQEdge)
                   +g->n_absorb*sizeof(HPCQAbsorbEntry)+sizeof(HPCQGraph));
    printf("%4d×%-4d %7dq %3dcyc ed=%-8lu abs=%-5lu mem=%6.1fMB t=%.4fs\n",
           R,C,n,cyc,(unsigned long)g->n_edges,(unsigned long)g->n_absorb,
           mem/(1024.*1024.),t);
    hpcq_destroy(g);return 1;}

int main(int argc,char**argv){
    setbuf(stdout,NULL);
    printf("╔═════════════════════════════════════════════════════════════╗\n");
    printf("║  BitState HPC — Maximum Willow Bench                      ║\n");
    printf("║  Usage: %s [--cycle N]                     ║\n",argv[0]);
    printf("╚══════════════════════════════════════════════════════════════╝\n");
    printf("%-12s %7s %6s %-10s %-7s %-8s %s\n","Grid","Qubits","Cyc","Edges","Absorb","Memory","Time");

    int user_cycles=0;
    for(int i=1;i<argc;i++){if(!strcmp(argv[i],"--cycle")&&i+1<argc)user_cycles=atoi(argv[++i]);}

    int tiers[][3]={{100,100,100},{150,150,60},{200,200,40},{300,300,20},
                    {400,400,10},{500,500,6},{700,700,3},{1000,1000,2},
                    {1500,1500,1},{2000,2000,1},{0,0,0}};

    if(user_cycles>0){
        printf("  (user overrides: all tiers → %d cycles)\n",user_cycles);
        for(int i=0;tiers[i][0];i++)tiers[i][2]=user_cycles;}

    for(int i=0;tiers[i][0];i++)
        if(!try_bench(tiers[i][0],tiers[i][1],tiers[i][2]))break;
    printf("\nDone.\n");}

/* hpc_qubit_graph_bethe.h — Bethe (single-pass) approximation
 * Include AFTER hpc_qubit_graph.h. Provides hpcq_amplitude_bethe().
 *
 * Algorithm: exact per-center chain marginals (forward-backward) ×
 * per-edge expected cc/half-absorbed edge weights computed from
 * those marginals.  Single-pass, no iteration, no BP.
 *
 * For tree-like cc graphs this is exact. For loopy graphs it's an
 * approximation that measures each edge's expected contribution. */
static inline void hpcq_amplitude_bethe(const HPCQGraph *g,
                                         const uint32_t *indices,
                                         double *out_re, double *out_im)
{
    uint64_t n_absorb = g->n_absorb;
    if (n_absorb <= 1) { hpcq_amplitude(g, indices, out_re, out_im); return; }
    static const double SQ = 0.7071067811865475244;
    /* ── Steps 1+2: local product ── */
    double re = 1.0, im = 0.0;
    for (uint64_t k = 0; k < g->n_sites; k++) {
        uint32_t idx = indices[k]; double ar, ai;
        if (g->absorb_idx[k] >= 0) { ar = 1.0; ai = 0.0; }
        else {
            double rb[2], ib[2];
            tri_get_amplitudes((TrialityQubit *)&g->locals[k], VIEW_EDGE, rb, ib);
            ar = rb[idx]; ai = ib[idx];
        }
        double nr = re*ar - im*ai; im = re*ai + im*ar; re = nr;
    }
    for (uint64_t e = 0; e < g->n_edges; e++) {
        const HPCQEdge *edge = &g->edges[e];
        if (edge->type == HPCQ_EDGE_ABSORBED) continue;
        if (edge->type == HPCQ_EDGE_CZ) {
            int ca = (g->absorb_idx[edge->site_a] >= 0) ? 1 : 0;
            int cb = (g->absorb_idx[edge->site_b] >= 0) ? 1 : 0;
            if (ca || cb) continue;
        }
        uint32_t ia = indices[edge->site_a], ib = indices[edge->site_b];
        double wr, wi;
        if (edge->type == HPCQ_EDGE_CZ) {
            wr = HPCQ_CZ_W_RE(ia, ib, edge->xp_a, edge->xp_b);
            wi = HPCQ_CZ_W_IM(ia, ib, edge->xp_a, edge->xp_b);
        } else { wr = edge->w_re[ia][ib]; wi = edge->w_im[ia][ib]; }
        double nr = re*wr - im*wi; im = re*wi + im*wr; re = nr;
    }

    /* ── Step 3a: sf, so_layer, na_cz ── */
    double (*sf_re)[2] = calloc(n_absorb, 2*sizeof(double));
    double (*sf_im)[2] = calloc(n_absorb, 2*sizeof(double));
    double **so_layer_re = calloc(n_absorb, sizeof(double *));
    double **so_layer_im = calloc(n_absorb, sizeof(double *));
    double *na_cz_re = calloc(n_absorb, sizeof(double));
    double *na_cz_im = calloc(n_absorb, sizeof(double));
    int *nl_arr = calloc(n_absorb, sizeof(int));
    for (uint64_t a = 0; a < n_absorb; a++) { na_cz_re[a]=1.0; na_cz_im[a]=0.0; }

    for (uint64_t a = 0; a < n_absorb; a++) {
        uint64_t nn = g->absorb[a].n_nbrs;
        int L = (int)g->absorb[a].n_layers;
        nl_arr[a] = L;
        double pi_re[2] = {1,1}, pi_im[2] = {0,0};
        if (L > 1) {
            so_layer_re[a] = calloc((L-1)*2, sizeof(double));
            so_layer_im[a] = calloc((L-1)*2, sizeof(double));
            for (int li = 0; li < L-1; li++) {
                so_layer_re[a][li*2] = 1.0; so_layer_re[a][li*2+1] = 1.0;
                so_layer_im[a][li*2] = 0.0; so_layer_im[a][li*2+1] = 0.0;
            }
        }
        for (uint64_t k = 0; k < nn; k++) {
            uint64_t nb = g->absorb[a].nbrs[k];
            if (g->absorb_idx[nb] >= 0) continue;
            int li = (int)g->absorb[a].layer[k];
            int z = (int)indices[nb];
            for (int y = 0; y < 2; y++) {
                int idx = (int)k*4 + y*2 + z;
                double wr = g->absorb[a].w_re[idx];
                double wi = g->absorb[a].w_im[idx];
                if (li == 0) {
                    double pr = pi_re[y], pim = pi_im[y];
                    pi_re[y] = pr*wr - pim*wi;
                    pi_im[y] = pr*wi + pim*wr;
                } else {
                    int sli = (li-1)*2 + y;
                    double pr = so_layer_re[a][sli];
                    double pim = so_layer_im[a][sli];
                    so_layer_re[a][sli] = pr*wr - pim*wi;
                    so_layer_im[a][sli] = pr*wi + pim*wr;
                }
            }
        }
        if (L >= 1) {
            uint64_t q = g->absorb[a].center;
            for (uint64_t e = 0; e < g->n_edges; e++) {
                const HPCQEdge *edge = &g->edges[e];
                if (edge->type != HPCQ_EDGE_CZ) continue;
                uint64_t other;
                if (edge->site_a == q && g->absorb_idx[edge->site_b] < 0)
                    other = edge->site_b;
                else if (edge->site_b == q && g->absorb_idx[edge->site_a] < 0)
                    other = edge->site_a;
                else continue;
                if (L == 1) {
                    int z = (int)indices[other];
                    int xv = (int)(indices[q] ^ g->absorb[a].x_parity);
                    double xpq = (edge->site_a==q) ? edge->xp_a : edge->xp_b;
                    double xpo = (edge->site_a==q) ? edge->xp_b : edge->xp_a;
                    double wr = HPCQ_CZ_W(xv, z, xpq, xpo);
                    na_cz_re[a] *= wr; na_cz_im[a] *= wr;
                } else {
                    int z = (int)indices[other];
                    double xpq = (edge->site_a==q) ? edge->xp_a : edge->xp_b;
                    double xpo = (edge->site_a==q) ? edge->xp_b : edge->xp_a;
                    int osli = (L-2)*2;
                    for (int yo = 0; yo < 2; yo++) {
                        double wr = HPCQ_CZ_W(yo, z, xpq, xpo);
                        so_layer_re[a][osli+yo] *= wr;
                        so_layer_im[a][osli+yo] *= wr;
                    }
                }
            }
        }
        for (int y = 0; y < 2; y++) {
            sf_re[a][y] = g->absorb[a].a_re[y]*pi_re[y]
                        - g->absorb[a].a_im[y]*pi_im[y];
            sf_im[a][y] = g->absorb[a].a_re[y]*pi_im[y]
                        + g->absorb[a].a_im[y]*pi_re[y];
        }
    }

    /* ── Per-center chain marginals (fwd·bwd, normalized) ── */
    double **marg_re = calloc(n_absorb, sizeof(double *));
    double **marg_im = calloc(n_absorb, sizeof(double *));
    int *marg_L = calloc(n_absorb, sizeof(int));

    for (uint64_t mi = 0; mi < n_absorb; mi++) {
        int L = nl_arr[mi];
        marg_L[mi] = L;
        if (L < 1) continue;
        marg_re[mi] = calloc(L*2, sizeof(double));
        marg_im[mi] = calloc(L*2, sizeof(double));
        uint64_t xv = indices[g->absorb[mi].center] ^ g->absorb[mi].x_parity;

        /* Forward pass */
        double f_re[128][2], f_im[128][2];
        f_re[0][0]=sf_re[mi][0]; f_re[0][1]=sf_re[mi][1];
        f_im[0][0]=sf_im[mi][0]; f_im[0][1]=sf_im[mi][1];
        for (int li = 1; li < L; li++) {
            double ar0=g->absorb[mi].a_layer_re[(li-1)*2],
                   ai0=g->absorb[mi].a_layer_im[(li-1)*2];
            double ar1=g->absorb[mi].a_layer_re[(li-1)*2+1],
                   ai1=g->absorb[mi].a_layer_im[(li-1)*2+1];
            double or0=so_layer_re[mi][(li-1)*2],
                   oi0=so_layer_im[mi][(li-1)*2];
            double or1=so_layer_re[mi][(li-1)*2+1],
                   oi1=so_layer_im[mi][(li-1)*2+1];
            double s0r=ar0*or0-ai0*oi0, s0i=ar0*oi0+ai0*or0;
            double s1r=ar1*or1-ai1*oi1, s1i=ar1*oi1+ai1*or1;
            for (int yc = 0; yc < 2; yc++) {
                double H0=SQ, H1=(yc==0?SQ:-SQ);
                double sfr=(yc==0)?s0r:s1r, sfi=(yc==0)?s0i:s1i;
                f_re[li][yc]=H0*(sfr*f_re[li-1][0]-sfi*f_im[li-1][0])
                            +H1*(sfr*f_re[li-1][1]-sfi*f_im[li-1][1]);
                f_im[li][yc]=H0*(sfr*f_im[li-1][0]+sfi*f_re[li-1][0])
                            +H1*(sfr*f_im[li-1][1]+sfi*f_re[li-1][1]);
            }
        }
        /* Backward pass */
        double b_re[128][2], b_im[128][2];
        double lst_re[2], lst_im[2];
        tri_get_amplitudes((TrialityQubit *)&g->locals[g->absorb[mi].center],
                           VIEW_EDGE, lst_re, lst_im);
        double lr=lst_re[xv], li=lst_im[xv];
        double ncr=na_cz_re[mi], nci=na_cz_im[mi];
        double pr=lr*ncr-li*nci, pi=lr*nci+li*ncr;
        b_re[L-1][0]=SQ*pr;       b_im[L-1][0]=SQ*pi;
        b_re[L-1][1]=(xv==0?SQ:-SQ)*pr; b_im[L-1][1]=(xv==0?SQ:-SQ)*pi;
        for (int li = L-2; li >= 0; li--) {
            double ar0=g->absorb[mi].a_layer_re[li*2],
                   ai0=g->absorb[mi].a_layer_im[li*2];
            double ar1=g->absorb[mi].a_layer_re[li*2+1],
                   ai1=g->absorb[mi].a_layer_im[li*2+1];
            double or0=so_layer_re[mi][li*2], oi0=so_layer_im[mi][li*2];
            double or1=so_layer_re[mi][li*2+1], oi1=so_layer_im[mi][li*2+1];
            double s0r=ar0*or0-ai0*oi0, s0i=ar0*oi0+ai0*or0;
            double s1r=ar1*or1-ai1*oi1, s1i=ar1*oi1+ai1*or1;
            for (int yp = 0; yp < 2; yp++) {
                b_re[li][yp] = 0; b_im[li][yp] = 0;
                for (int yc = 0; yc < 2; yc++) {
                    double H = (yc==0) ? SQ : (yp==0 ? SQ : -SQ);
                    double sfr=(yc==0)?s0r:s1r, sfi=(yc==0)?s0i:s1i;
                    b_re[li][yp] += H*(sfr*b_re[li+1][yc]-sfi*b_im[li+1][yc]);
                    b_im[li][yp] += H*(sfr*b_im[li+1][yc]+sfi*b_re[li+1][yc]);
                }
            }
        }
        /* Marginal = f·b, normalized to probability */
        for (int li = 0; li < L; li++) {
            double p0r=f_re[li][0]*b_re[li][0]-f_im[li][0]*b_im[li][0];
            double p0i=f_re[li][0]*b_im[li][0]+f_im[li][0]*b_re[li][0];
            double p1r=f_re[li][1]*b_re[li][1]-f_im[li][1]*b_im[li][1];
            double p1i=f_re[li][1]*b_im[li][1]+f_im[li][1]*b_re[li][1];
            double nrm=sqrt(p0r*p0r+p0i*p0i+p1r*p1r+p1i*p1i);
            if (nrm > 1e-30) { p0r/=nrm; p0i/=nrm; p1r/=nrm; p1i/=nrm; }
            marg_re[mi][li*2]=p0r;   marg_im[mi][li*2]=p0i;
            marg_re[mi][li*2+1]=p1r; marg_im[mi][li*2+1]=p1i;
        }
    }

    /* ── Per-center chain amplitude (no cc edges yet) ── */
    double comp_re = 1.0, comp_im = 0.0;
    for (uint64_t mi = 0; mi < n_absorb; mi++) {
        int L = marg_L[mi]; if (L < 1) continue;
        uint64_t xv = indices[g->absorb[mi].center] ^ g->absorb[mi].x_parity;
        double c_re[2], c_im[2];
        c_re[0]=sf_re[mi][0]; c_re[1]=sf_re[mi][1];
        c_im[0]=sf_im[mi][0]; c_im[1]=sf_im[mi][1];
        for (int li = 1; li < L; li++) {
            double ar0=g->absorb[mi].a_layer_re[(li-1)*2],
                   ai0=g->absorb[mi].a_layer_im[(li-1)*2];
            double ar1=g->absorb[mi].a_layer_re[(li-1)*2+1],
                   ai1=g->absorb[mi].a_layer_im[(li-1)*2+1];
            double or0=so_layer_re[mi][(li-1)*2],
                   oi0=so_layer_im[mi][(li-1)*2];
            double or1=so_layer_re[mi][(li-1)*2+1],
                   oi1=so_layer_im[mi][(li-1)*2+1];
            double s0r=ar0*or0-ai0*oi0, s0i=ar0*oi0+ai0*or0;
            double s1r=ar1*or1-ai1*oi1, s1i=ar1*oi1+ai1*or1;
            double n0r=0,n0i=0,n1r=0,n1i=0;
            for (int yp = 0; yp < 2; yp++) {
                n0r += SQ*(s0r*c_re[yp]-s0i*c_im[yp]);
                n0i += SQ*(s0r*c_im[yp]+s0i*c_re[yp]);
                double H = (yp==0?SQ:-SQ);
                n1r += H*(s1r*c_re[yp]-s1i*c_im[yp]);
                n1i += H*(s1r*c_im[yp]+s1i*c_re[yp]);
            }
            c_re[0]=n0r; c_re[1]=n1r; c_im[0]=n0i; c_im[1]=n1i;
        }
        double lst_re[2], lst_im[2];
        tri_get_amplitudes((TrialityQubit *)&g->locals[g->absorb[mi].center],
                           VIEW_EDGE, lst_re, lst_im);
        double lr=lst_re[xv], li=lst_im[xv];
        double ncr=na_cz_re[mi], nci=na_cz_im[mi];
        double pr=lr*ncr-li*nci, pi=lr*nci+li*ncr;
        double H0=SQ, H1=(xv==0?SQ:-SQ);
        double sr=H0*(pr*c_re[0]-pi*c_im[0])+H1*(pr*c_re[1]-pi*c_im[1]);
        double si=H0*(pr*c_im[0]+pi*c_re[0])+H1*(pr*c_im[1]+pi*c_re[1]);
        double nr=comp_re*sr-comp_im*si;
        comp_im=comp_re*si+comp_im*sr; comp_re=nr;
    }

    /* ── Per-cc-edge expected weight: Σ p_a(ya)p_b(yb)w(ya,yb) ── */
    for (uint64_t mi = 0; mi < n_absorb; mi++) {
        uint64_t a = mi;
        for (uint64_t k = 0; k < g->absorb[a].n_nbrs; k++) {
            uint64_t nb = g->absorb[a].nbrs[k];
            int64_t aj = g->absorb_idx[nb];
            if (aj < 0) continue;
            if (a >= (uint64_t)aj) continue;  /* apply once per edge */
            int li = (int)g->absorb[a].layer[k];
            /* Order-match: i-th occurrence */
            int occ_a = 0;
            for (uint64_t kk = 0; kk <= k; kk++)
                if (g->absorb[a].nbrs[kk] == nb) occ_a++;
            int occ_b = 0, lj = -1;
            for (uint64_t k2 = 0; k2 < g->absorb[(uint64_t)aj].n_nbrs; k2++) {
                if (g->absorb[(uint64_t)aj].nbrs[k2] == g->absorb[a].center) {
                    occ_b++;
                    if (occ_b == occ_a) {
                        lj = (int)g->absorb[(uint64_t)aj].layer[k2]; break;
                    }
                }
            }
            if (lj < 0) continue;
            double ew_r = 0, ew_i = 0;
            for (int ya = 0; ya < 2; ya++) {
                for (int yb = 0; yb < 2; yb++) {
                    int wk = (int)k*4 + ya*2 + yb;
                    double wr = g->absorb[a].w_re[wk];
                    double wi = g->absorb[a].w_im[wk];
                    double pa_r = marg_re[mi][li*2+ya];
                    double pa_i = marg_im[mi][li*2+ya];
                    double pb_r = marg_re[(uint64_t)aj][lj*2+yb];
                    double pb_i = marg_im[(uint64_t)aj][lj*2+yb];
                    double p_r = pa_r*pb_r - pa_i*pb_i;
                    double p_i = pa_r*pb_i + pa_i*pb_r;
                    ew_r += wr*p_r - wi*p_i;
                    ew_i += wr*p_i + wi*p_r;
                }
            }
            double nr=comp_re*ew_r-comp_im*ew_i;
            comp_im=comp_re*ew_i+comp_im*ew_r; comp_re=nr;
        }
    }

    /* ── Half-absorbed: expected w(ya, xv_b) over p_a(ya) ── */
    for (uint64_t mi = 0; mi < n_absorb; mi++) {
        uint64_t a = mi;
        for (uint64_t k = 0; k < g->absorb[a].n_nbrs; k++) {
            uint64_t nb = g->absorb[a].nbrs[k];
            int64_t aj = g->absorb_idx[nb];
            if (aj < 0) continue;
            int li = (int)g->absorb[a].layer[k];
            int occ_a = 0;
            for (uint64_t kk = 0; kk <= k; kk++)
                if (g->absorb[a].nbrs[kk] == nb) occ_a++;
            int occ_b = 0, lj = -1;
            for (uint64_t k2 = 0; k2 < g->absorb[(uint64_t)aj].n_nbrs; k2++) {
                if (g->absorb[(uint64_t)aj].nbrs[k2] == g->absorb[a].center) {
                    occ_b++;
                    if (occ_b == occ_a) {
                        lj = (int)g->absorb[(uint64_t)aj].layer[k2]; break;
                    }
                }
            }
            if (lj >= 0) continue;  /* fully-matched, handled above */
            uint64_t xv_b = indices[nb] ^ g->absorb[(uint64_t)aj].x_parity;
            double ew_r = 0, ew_i = 0;
            for (int ya = 0; ya < 2; ya++) {
                int wk = (int)k*4 + ya*2 + (int)xv_b;
                double wr = g->absorb[a].w_re[wk];
                double wi = g->absorb[a].w_im[wk];
                double pa_r = marg_re[mi][li*2+ya];
                double pa_i = marg_im[mi][li*2+ya];
                ew_r += wr*pa_r - wi*pa_i;
                ew_i += wr*pa_i + wi*pa_r;
            }
            double nr=comp_re*ew_r-comp_im*ew_i;
            comp_im=comp_re*ew_i+comp_im*ew_r; comp_re=nr;
        }
    }

    *out_re = re*comp_re - im*comp_im;
    *out_im = re*comp_im + im*comp_re;

    /* Cleanup */
    for (uint64_t a = 0; a < n_absorb; a++) {
        free(marg_re[a]); free(marg_im[a]);
        free(so_layer_re[a]); free(so_layer_im[a]);
    }
    free(marg_re); free(marg_im); free(marg_L);
    free(sf_re); free(sf_im);
    free(so_layer_re); free(so_layer_im);
    free(na_cz_re); free(na_cz_im); free(nl_arr);
}

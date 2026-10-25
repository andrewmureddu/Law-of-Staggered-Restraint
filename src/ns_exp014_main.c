#include "ns_c_header.h"

typedef struct {
    double t, E, W, Z, D, P, Pplus, Pminus, Q, S_max, eta_plus, C_conc, k_char;
    double T4, V3, T4_plus, V3_plus, J, R;
    double T4_Q, V3_Q, Qdot, R_Q;
    double H, Delta_triad, Delta_visc, Zdot;
    double T4_B1, T4_B2, T4_B3, S_dyad;
    double T4_P1, T4_P2, T4_P3, S_pou, pou_err;
    double ddP_bulk, ddP_bdry, ddP_plus;
    double ddQ_bulk, ddQ_bdry, ddQ;
    double dH, Bdry_P_term, Bdry_Q_term;
    double Delta_T5;
} SnapDiag;

typedef struct {
    double complex *wh, *lamb_h, *Nu_h, *Nw_h, *du_h, *ddu_h, *ddw_h, *lap_wh, *lap_uh;
    double complex *tmp1, *tmp2, *b1, *b2;
    double *u, *w, *lamb, *Nw, *lap_w, *du, *dw, *ddw;
    double *S, *S_Nu, *lap_S, *dS, *Sddot, *S_shell;
    double *Sw, *p, *q2, *snorm, *lam2, *dlam2, *dp;
} WorkBuffers;

void alloc_work(WorkBuffers *wb, int N3) {
    wb->wh = calloc(3*N3, sizeof(double complex));
    wb->lamb_h = calloc(3*N3, sizeof(double complex));
    wb->Nu_h = calloc(3*N3, sizeof(double complex));
    wb->Nw_h = calloc(3*N3, sizeof(double complex));
    wb->du_h = calloc(3*N3, sizeof(double complex));
    wb->ddu_h = calloc(3*N3, sizeof(double complex));
    wb->ddw_h = calloc(3*N3, sizeof(double complex));
    wb->lap_wh = calloc(3*N3, sizeof(double complex));
    wb->lap_uh = calloc(3*N3, sizeof(double complex));
    wb->tmp1 = calloc(3*N3, sizeof(double complex));
    wb->tmp2 = calloc(3*N3, sizeof(double complex));
    wb->b1 = calloc(N3, sizeof(double complex));
    wb->b2 = calloc(N3, sizeof(double complex));

    wb->u = calloc(3*N3, sizeof(double));
    wb->w = calloc(3*N3, sizeof(double));
    wb->lamb = calloc(3*N3, sizeof(double));
    wb->Nw = calloc(3*N3, sizeof(double));
    wb->lap_w = calloc(3*N3, sizeof(double));
    wb->du = calloc(3*N3, sizeof(double));
    wb->dw = calloc(3*N3, sizeof(double));
    wb->ddw = calloc(3*N3, sizeof(double));

    wb->S = calloc(6*N3, sizeof(double));
    wb->S_Nu = calloc(6*N3, sizeof(double));
    wb->lap_S = calloc(6*N3, sizeof(double));
    wb->dS = calloc(6*N3, sizeof(double));
    wb->Sddot = calloc(6*N3, sizeof(double));
    wb->S_shell = calloc(6*N3, sizeof(double));

    wb->Sw = calloc(3*N3, sizeof(double));
    wb->p = calloc(N3, sizeof(double));
    wb->q2 = calloc(N3, sizeof(double));
    wb->snorm = calloc(N3, sizeof(double));
    wb->lam2 = calloc(N3, sizeof(double));
    wb->dlam2 = calloc(N3, sizeof(double));
    wb->dp = calloc(N3, sizeof(double));
}

static inline void sym6_vec(const double *S, int N3, int i, double v0, double v1, double v2,
                            double *o0, double *o1, double *o2) {
    double s00 = S[0*N3+i], s01 = S[1*N3+i], s02 = S[2*N3+i];
    double s11 = S[3*N3+i], s12 = S[4*N3+i], s22 = S[5*N3+i];
    *o0 = s00*v0 + s01*v1 + s02*v2;
    *o1 = s01*v0 + s11*v1 + s12*v2;
    *o2 = s02*v0 + s12*v1 + s22*v2;
}

static inline double sym6_quad(const double *S, int N3, int i, double v0, double v1, double v2) {
    double o0, o1, o2;
    sym6_vec(S, N3, i, v0, v1, v2, &o0, &o1, &o2);
    return v0*o0 + v1*o1 + v2*o2;
}

double compute_shell_t4(Solver3D *s, WorkBuffers *wb, const double *wt_mask) {
    int N3 = s->N3;
    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < N3; i++) {
            wb->tmp1[c*N3 + i] = wb->Nu_h[c*N3 + i] * wt_mask[i];
            wb->tmp2[c*N3 + i] = wb->Nw_h[c*N3 + i] * wt_mask[i];
        }
    }
    ift_vec3(s, wb->tmp2, wb->du, wb->b1, wb->b2, wb->ddu_h); // reuse du for Nw_b
    strain_from_hat(s, wb->tmp1, wb->S_shell, wb->ddu_h, wb->ddw_h, wb->b1, wb->b2);
    double sum = 0.0;
    for (int i = 0; i < N3; i++) {
        if (wb->p[i] > 0.0) {
            double nwx = wb->du[i], nwy = wb->du[N3+i], nwz = wb->du[2*N3+i];
            double swx = wb->Sw[i], swy = wb->Sw[N3+i], swz = wb->Sw[2*N3+i];
            double wx = wb->w[i], wy = wb->w[N3+i], wz = wb->w[2*N3+i];
            double q_s = sym6_quad(wb->S_shell, N3, i, wx, wy, wz);
            sum += 2.0*(nwx*swx + nwy*swy + nwz*swz) + q_s;
        }
    }
    return sum / (double)N3;
}

void compute_diagnostics(Solver3D *s, WorkBuffers *wb, const double complex *h, double t,
                         int do_2nd, int do_shells, SnapDiag *out) {
    int N3 = s->N3;
    double nu = s->nu;
    memset(out, 0, sizeof(SnapDiag));
    out->t = t;

    ift_vec3(s, h, wb->u, wb->tmp1, wb->b1, wb->b2);
    curl_hat(s, h, wb->wh);
    ift_vec3(s, wb->wh, wb->w, wb->tmp1, wb->b1, wb->b2);
    strain_from_hat(s, h, wb->S, wb->tmp1, wb->tmp2, wb->b1, wb->b2);

    // Dissipation D via Parseval or physical grad w: by Parseval, <|grad w|^2> = sum_k |k|^2 |wh(k)|^2 / N3^2
    double sum_k2_wh2 = 0.0;
    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < N3; i++) {
            double re = creal(wb->wh[c*N3+i]), im = cimag(wb->wh[c*N3+i]);
            sum_k2_wh2 += s->k2[i] * (re*re + im*im);
        }
    }
    double grad_w2_mean = sum_k2_wh2 / ((double)N3 * (double)N3);
    out->D = nu * grad_w2_mean;

    double sum_E = 0.0, sum_W = 0.0, sum_q4 = 0.0, sum_P = 0.0, sum_Pp = 0.0, sum_Pm = 0.0;
    for (int i = 0; i < N3; i++) {
        double ux = wb->u[i], uy = wb->u[N3+i], uz = wb->u[2*N3+i];
        double wx = wb->w[i], wy = wb->w[N3+i], wz = wb->w[2*N3+i];
        sum_E += 0.5 * (ux*ux + uy*uy + uz*uz);
        double q2 = wx*wx + wy*wy + wz*wz;
        wb->q2[i] = q2;
        sum_W += q2;
        sum_q4 += q2 * q2;

        double swx, swy, swz;
        sym6_vec(wb->S, N3, i, wx, wy, wz, &swx, &swy, &swz);
        wb->Sw[i] = swx; wb->Sw[N3+i] = swy; wb->Sw[2*N3+i] = swz;
        double pv = wx*swx + wy*swy + wz*swz;
        wb->p[i] = pv;
        sum_P += pv;
        if (pv > 0.0) sum_Pp += pv;
        else if (pv < 0.0) sum_Pm += (-pv);

        // Lamb = u x w
        wb->lamb[i]      = uy*wz - uz*wy;
        wb->lamb[N3+i]   = uz*wx - ux*wz;
        wb->lamb[2*N3+i] = ux*wy - uy*wx;
    }
    out->E = sum_E / N3;
    out->W = sum_W / N3;
    out->Z = 0.5 * out->W;
    out->P = sum_P / N3;
    out->Pplus = sum_Pp / N3;
    out->Pminus = sum_Pm / N3;
    out->C_conc = (out->W > 1e-15) ? ((sum_q4 / N3) / (out->W * out->W)) : 1.0;
    out->k_char = sqrt(fmax(out->D / fmax(nu * out->W, 1e-30), 1e-12));

    ft_vec3(s, wb->lamb, wb->Nu_h, wb->tmp1, wb->b1, wb->b2);
    project_hat(s, wb->Nu_h);
    curl_hat(s, wb->Nu_h, wb->Nw_h);
    ift_vec3(s, wb->Nw_h, wb->Nw, wb->tmp1, wb->b1, wb->b2);
    strain_from_hat(s, wb->Nu_h, wb->S_Nu, wb->tmp1, wb->tmp2, wb->b1, wb->b2);

    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < N3; i++) {
            wb->lap_wh[c*N3+i] = s->k2[i] * wb->wh[c*N3+i];
            wb->lap_uh[c*N3+i] = s->k2[i] * h[c*N3+i];
        }
    }
    ift_vec3(s, wb->lap_wh, wb->lap_w, wb->tmp1, wb->b1, wb->b2);
    strain_from_hat(s, wb->lap_uh, wb->lap_S, wb->tmp1, wb->tmp2, wb->b1, wb->b2);

    if (do_2nd) {
        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < N3; i++) {
                wb->du_h[c*N3+i] = wb->Nu_h[c*N3+i] - nu * wb->lap_uh[c*N3+i];
                wb->dw[c*N3+i] = wb->Nw[c*N3+i] - nu * wb->lap_w[c*N3+i];
            }
        }
        ift_vec3(s, wb->du_h, wb->du, wb->tmp1, wb->b1, wb->b2);
        for (int p = 0; p < 6; p++) {
            for (int i = 0; i < N3; i++) wb->dS[p*N3+i] = wb->S_Nu[p*N3+i] - nu * wb->lap_S[p*N3+i];
        }
        // lamb_dot = du x w + u x dw
        for (int i = 0; i < N3; i++) {
            double ux = wb->u[i], uy = wb->u[N3+i], uz = wb->u[2*N3+i];
            double wx = wb->w[i], wy = wb->w[N3+i], wz = wb->w[2*N3+i];
            double dux = wb->du[i], duy = wb->du[N3+i], duz = wb->du[2*N3+i];
            double dwx = wb->dw[i], dwy = wb->dw[N3+i], dwz = wb->dw[2*N3+i];
            wb->lamb[i]      = (duy*wz - duz*wy) + (uy*dwz - uz*dwy);
            wb->lamb[N3+i]   = (duz*wx - dux*wz) + (uz*dwx - ux*dwz);
            wb->lamb[2*N3+i] = (dux*wy - duy*wx) + (ux*dwy - uy*dwx);
        }
        ft_vec3(s, wb->lamb, wb->ddu_h, wb->tmp1, wb->b1, wb->b2);
        project_hat(s, wb->ddu_h);
        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < N3; i++) wb->ddu_h[c*N3+i] -= nu * s->k2[i] * wb->du_h[c*N3+i];
        }
        curl_hat(s, wb->ddu_h, wb->ddw_h);
        ift_vec3(s, wb->ddw_h, wb->ddw, wb->tmp1, wb->b1, wb->b2);
        strain_from_hat(s, wb->ddu_h, wb->Sddot, wb->tmp1, wb->tmp2, wb->b1, wb->b2);
    }

    double sum_Smax = 0.0, sum_Q = 0.0, sum_T4Q = 0.0, sum_V3Q = 0.0;
    double sum_T4 = 0.0, sum_V3 = 0.0, sum_T4p = 0.0, sum_V3p = 0.0;
    double sum_ddP_bulk = 0.0, sum_ddQ_bulk = 0.0;
    double sum_T5p = 0.0, sum_T5Q = 0.0;

    for (int i = 0; i < N3; i++) {
        double S6[6], SdNL6[6], SdV6[6], Sdd6[6] = {0};
        for (int p = 0; p < 6; p++) {
            S6[p] = wb->S[p*N3+i];
            SdNL6[p] = wb->S_Nu[p*N3+i];
            SdV6[p] = nu * wb->lap_S[p*N3+i];
            if (do_2nd) Sdd6[p] = wb->Sddot[p*N3+i];
        }
        double snorm = 0, rNL = 0, rV = 0, dSop = 0, ddSop = 0, lam2 = 0, dlam2 = 0;
        point_strain_rates(S6, SdNL6, SdV6, Sdd6, do_2nd, &snorm, &rNL, &rV, &dSop, &ddSop, &lam2, &dlam2);
        wb->snorm[i] = snorm;
        wb->lam2[i] = lam2;
        wb->dlam2[i] = dlam2;

        double wx = wb->w[i], wy = wb->w[N3+i], wz = wb->w[2*N3+i];
        double nwx = wb->Nw[i], nwy = wb->Nw[N3+i], nwz = wb->Nw[2*N3+i];
        double lwx = wb->lap_w[i], lwy = wb->lap_w[N3+i], lwz = wb->lap_w[2*N3+i];
        double swx = wb->Sw[i], swy = wb->Sw[N3+i], swz = wb->Sw[2*N3+i];
        double q2 = wb->q2[i];

        sum_Smax += snorm;
        sum_Q += q2 * snorm;
        sum_T4Q += 2.0*(wx*nwx + wy*nwy + wz*nwz)*snorm + q2*rNL;
        sum_V3Q += 2.0*nu*(wx*lwx + wy*lwy + wz*lwz)*snorm + q2*rV;

        double t4_f = 2.0*(nwx*swx + nwy*swy + nwz*swz) + sym6_quad(wb->S_Nu, N3, i, wx, wy, wz);
        double v3_f = nu * (2.0*(lwx*swx + lwy*swy + lwz*swz) + sym6_quad(wb->lap_S, N3, i, wx, wy, wz));
        sum_T4 += t4_f;
        sum_V3 += v3_f;
        if (wb->p[i] > 0.0) {
            sum_T4p += t4_f;
            sum_V3p += v3_f;
        }
        wb->dp[i] = t4_f - v3_f;

        if (do_2nd) {
            double dwx = wb->dw[i], dwy = wb->dw[N3+i], dwz = wb->dw[2*N3+i];
            double ddwx = wb->ddw[i], ddwy = wb->ddw[N3+i], ddwz = wb->ddw[2*N3+i];
            double dswx, dswy, dswz, sdwx, sdwy, sdwz;
            sym6_vec(wb->dS, N3, i, wx, wy, wz, &dswx, &dswy, &dswz);
            sym6_vec(wb->S, N3, i, dwx, dwy, dwz, &sdwx, &sdwy, &sdwz);
            double ddp = 2.0*(ddwx*swx + ddwy*swy + ddwz*swz)
                       + 4.0*(dwx*dswx + dwy*dswy + dwz*dswz)
                       + 2.0*(dwx*sdwx + dwy*sdwy + dwz*sdwz)
                       + sym6_quad(wb->Sddot, N3, i, wx, wy, wz);
            if (wb->p[i] > 0.0) sum_ddP_bulk += ddp;

            double dw2 = 2.0*(wx*dwx + wy*dwy + wz*dwz);
            double ddw2 = 2.0*(dwx*dwx + dwy*dwy + dwz*dwz) + 2.0*(wx*ddwx + wy*ddwy + wz*ddwz);
            sum_ddQ_bulk += ddw2 * snorm + 2.0 * dw2 * dSop + q2 * ddSop;
        }
    }

    out->S_max = sum_Smax / N3;
    out->Q = sum_Q / N3;
    out->eta_plus = (out->Q > 1e-15) ? (out->Pplus / out->Q) : 0.0;
    out->T4_Q = sum_T4Q / N3;
    out->V3_Q = sum_V3Q / N3;
    out->Qdot = out->T4_Q - out->V3_Q;
    out->R_Q = out->T4_Q + out->V3_Q;

    out->T4 = sum_T4 / N3;
    out->V3 = sum_V3 / N3;
    out->T4_plus = sum_T4p / N3;
    out->V3_plus = sum_V3p / N3;
    out->J = out->T4_plus - out->V3_plus;
    out->R = out->T4_plus + out->V3_plus;
    out->Zdot = out->P - out->D;

    double Pp_s = fmax(out->Pplus, 1e-12);
    double Q_s = fmax(out->Q, 1e-12);
    out->H = (out->J / Pp_s) - (out->Qdot / Q_s);
    out->Delta_triad = (out->T4_plus / Pp_s) - (out->T4_Q / Q_s);
    out->Delta_visc = (out->V3_plus / Pp_s) - (out->V3_Q / Q_s);

    if (do_2nd) {
        // Compute std(p) and std(lam2) for smoothed surface deltas
        double mean_p = out->P, var_p = 0.0, mean_l2 = 0.0, var_l2 = 0.0;
        for (int i = 0; i < N3; i++) {
            double dp_m = wb->p[i] - mean_p;
            var_p += dp_m * dp_m;
            mean_l2 += wb->lam2[i];
        }
        mean_l2 /= N3;
        for (int i = 0; i < N3; i++) {
            double dl_m = wb->lam2[i] - mean_l2;
            var_l2 += dl_m * dl_m;
        }
        double eps_p = 0.08 * sqrt(var_p / N3) + 1e-12;
        double eps_l2 = 0.08 * sqrt(var_l2 / N3) + 1e-12;
        double sum_bdry_P = 0.0, sum_bdry_Q = 0.0;
        for (int i = 0; i < N3; i++) {
            double del_p = exp(-pow(wb->p[i]/eps_p, 2)) / (sqrt(PI) * eps_p);
            sum_bdry_P += del_p * wb->dp[i] * wb->dp[i];
            double del_l2 = exp(-pow(wb->lam2[i]/eps_l2, 2)) / (sqrt(PI) * eps_l2);
            sum_bdry_Q += wb->q2[i] * del_l2 * wb->dlam2[i] * wb->dlam2[i];
        }
        out->ddP_bulk = sum_ddP_bulk / N3;
        out->ddP_bdry = sum_bdry_P / N3;
        out->ddP_plus = out->ddP_bulk + out->ddP_bdry;
        out->ddQ_bulk = sum_ddQ_bulk / N3;
        out->ddQ_bdry = sum_bdry_Q / N3;
        out->ddQ = out->ddQ_bulk + out->ddQ_bdry;
        out->Bdry_P_term = out->ddP_bdry / Pp_s;
        out->Bdry_Q_term = out->ddQ_bdry / Q_s;
        out->dH = (out->ddP_plus / Pp_s) - pow(out->J / Pp_s, 2) - (out->ddQ / Q_s) + pow(out->Qdot / Q_s, 2);
    }

    if (do_shells) {
        double *wt = (double*)malloc(N3 * sizeof(double));
        for (int i = 0; i < N3; i++) wt[i] = (double)s->mask_B1[i];
        out->T4_B1 = compute_shell_t4(s, wb, wt);
        for (int i = 0; i < N3; i++) wt[i] = (double)s->mask_B2[i];
        out->T4_B2 = compute_shell_t4(s, wb, wt);
        out->T4_B3 = out->T4_plus - out->T4_B1 - out->T4_B2;
        out->S_dyad = (out->T4_B2 + out->T4_B3 - out->T4_B1) / Q_s;

        // Smooth scale-normalized Partition-of-Unity (PoU) shells in xi = log2(|k|/k_char)
        double *w1_arr = (double*)malloc(N3 * sizeof(double));
        double *w2_arr = (double*)malloc(N3 * sizeof(double));
        double *w3_arr = (double*)malloc(N3 * sizeof(double));
        for (int i = 0; i < N3; i++) {
            if (!s->mask[i] || s->kmag[i] == 0.0) {
                w1_arr[i] = 0.0; w2_arr[i] = 0.0; w3_arr[i] = 0.0;
            } else {
                double xi = log2(s->kmag[i] / out->k_char);
                double v1 = 1.0 / (1.0 + exp((xi + 0.25) / 0.25));
                double v3 = 1.0 / (1.0 + exp(-(xi - 0.35) / 0.25));
                w1_arr[i] = v1;
                w3_arr[i] = v3;
                w2_arr[i] = 1.0 - v1 - v3;
            }
        }
        out->T4_P1 = compute_shell_t4(s, wb, w1_arr);
        out->T4_P2 = compute_shell_t4(s, wb, w2_arr);
        out->T4_P3 = compute_shell_t4(s, wb, w3_arr);
        out->pou_err = fabs(out->T4_plus - (out->T4_P1 + out->T4_P2 + out->T4_P3));
        out->S_pou = (out->T4_P2 + out->T4_P3 - out->T4_P1) / Pp_s;
        free(wt); free(w1_arr); free(w2_arr); free(w3_arr);
    }
}

void rhs_hat(Solver3D *s, WorkBuffers *wb, const double complex *h, double complex *out_rhs) {
    int N3 = s->N3;
    ift_vec3(s, h, wb->u, wb->tmp1, wb->b1, wb->b2);
    curl_hat(s, h, wb->wh);
    ift_vec3(s, wb->wh, wb->w, wb->tmp1, wb->b1, wb->b2);
    for (int i = 0; i < N3; i++) {
        double ux = wb->u[i], uy = wb->u[N3+i], uz = wb->u[2*N3+i];
        double wx = wb->w[i], wy = wb->w[N3+i], wz = wb->w[2*N3+i];
        wb->lamb[i]      = uy*wz - uz*wy;
        wb->lamb[N3+i]   = uz*wx - ux*wz;
        wb->lamb[2*N3+i] = ux*wy - uy*wx;
    }
    ft_vec3(s, wb->lamb, out_rhs, wb->tmp1, wb->b1, wb->b2);
    project_hat(s, out_rhs);
    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < N3; i++) out_rhs[c*N3+i] -= s->nu * s->k2[i] * h[c*N3+i];
    }
}

void step_rk4(Solver3D *s, WorkBuffers *wb, double complex *h, double dt,
              double complex *k1, double complex *k2, double complex *k3, double complex *k4, double complex *htmp) {
    int ntot = 3 * s->N3;
    rhs_hat(s, wb, h, k1);
    for (int i = 0; i < ntot; i++) htmp[i] = h[i] + 0.5 * dt * k1[i];
    rhs_hat(s, wb, htmp, k2);
    for (int i = 0; i < ntot; i++) htmp[i] = h[i] + 0.5 * dt * k2[i];
    rhs_hat(s, wb, htmp, k3);
    for (int i = 0; i < ntot; i++) htmp[i] = h[i] + dt * k3[i];
    rhs_hat(s, wb, htmp, k4);
    for (int i = 0; i < ntot; i++) h[i] += (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    project_hat(s, h);
}

// Simple deterministic PRNG for reproducible initial conditions across platforms
typedef struct { uint64_t state; } XorShift64;
double rand_uniform(XorShift64 *rng) {
    uint64_t x = rng->state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng->state = x;
    return (double)(x >> 11) * (1.0 / 9007199254740992.0);
}

void make_initial_ensemble(Solver3D *s, WorkBuffers *wb, int fam_idx, int variant, double amp, int sign, double complex *h_out) {
    int N3 = s->N3;
    memset(wb->u, 0, 3*N3*sizeof(double));
    for (int i = 0; i < N3; i++) {
        double x = s->x[i], y = s->y[i], z = s->z[i];
        if (fam_idx == 0) { // exp013_3d
            double eps = 0.10 + 0.08 * variant;
            wb->u[i]      = sin(y) + eps * sin(z);
            wb->u[N3+i]   = eps * sin(x + z);
            wb->u[2*N3+i] = sin(x) - sin(x + y) + eps * sin(2.0*x + y);
        } else if (fam_idx == 1) { // tg3d_pert
            double delta = 0.22 * (variant + 1);
            wb->u[i]      = sin(x)*cos(y)*cos(z) + delta*sin(2.0*y)*cos(z);
            wb->u[N3+i]   = -cos(x)*sin(y)*cos(z) - delta*cos(x)*sin(2.0*z);
            wb->u[2*N3+i] = delta*sin(x)*sin(y)*sin(2.0*z);
        } else if (fam_idx == 2) { // antiparallel
            double alpha = 0.35 + 0.20 * variant;
            wb->u[i]      = sin(y)*cos(z) + alpha*cos(y + z);
            wb->u[N3+i]   = sin(z)*cos(x) - alpha*cos(x - z);
            wb->u[2*N3+i] = -cos(x)*sin(y) + alpha*sin(x + y)*cos(z);
        } else if (fam_idx == 3) { // kida_elliptic
            double eps = 0.15 + 0.10 * variant;
            wb->u[i]      = sin(x)*(cos(2.0*y)*cos(z) - cos(y)*cos(2.0*z)) + eps*cos(z);
            wb->u[N3+i]   = sin(y)*(cos(2.0*z)*cos(x) - cos(z)*cos(2.0*x)) + eps*cos(x);
            wb->u[2*N3+i] = sin(z)*(cos(2.0*x)*cos(y) - cos(x)*cos(2.0*y)) + eps*cos(y);
        } else { // random_solenoidal (deterministic seed = 202600 + 100*variant)
            XorShift64 rng = { (uint64_t)(202642 + 1337 * variant) };
            double ux = 0.0, uy = 0.0, uz = 0.0;
            for (int kx = 1; kx <= 3; kx++) {
                for (int ky = 1; ky <= 3; ky++) {
                    for (int kz = 1; kz <= 3; kz++) {
                        int k2 = kx*kx + ky*ky + kz*kz;
                        if (k2 <= 14) {
                            double p0 = 2.0*PI*rand_uniform(&rng);
                            double p1 = 2.0*PI*rand_uniform(&rng);
                            double p2 = 2.0*PI*rand_uniform(&rng);
                            double wt = exp(-0.25 * k2);
                            double arg = kx*x + ky*y + kz*z;
                            ux += wt * cos(arg + p0);
                            uy += wt * cos(arg + p1);
                            uz += wt * cos(arg + p2);
                        }
                    }
                }
            }
            wb->u[i] = ux; wb->u[N3+i] = uy; wb->u[2*N3+i] = uz;
        }
    }
    ft_vec3(s, wb->u, h_out, wb->tmp1, wb->b1, wb->b2);
    project_hat(s, h_out);
    curl_hat(s, h_out, wb->wh);
    ift_vec3(s, wb->wh, wb->w, wb->tmp1, wb->b1, wb->b2);
    double sum_W = 0.0;
    for (int i = 0; i < 3*N3; i++) sum_W += wb->w[i] * wb->w[i];
    double W_raw = sum_W / N3;
    double scale = (W_raw > 1e-12) ? (amp * sqrt(2.0 / W_raw)) : amp;
    for (int i = 0; i < 3*N3; i++) h_out[i] *= (sign * scale);
    project_hat(s, h_out);
}

// 4x4 linear system solver for Matched-State Counterfactual Twin Newton step
void solve_4x4(double A[4][4], double b[4], double x[4]) {
    double M[4][5];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) M[i][j] = A[i][j];
        M[i][4] = b[i];
    }
    for (int k = 0; k < 4; k++) {
        int piv = k;
        for (int i = k+1; i < 4; i++) if (fabs(M[i][k]) > fabs(M[piv][k])) piv = i;
        for (int j = k; j < 5; j++) { double t = M[k][j]; M[k][j] = M[piv][j]; M[piv][j] = t; }
        double div = M[k][k];
        if (fabs(div) < 1e-12) div = 1e-12;
        for (int j = k; j < 5; j++) M[k][j] /= div;
        for (int i = 0; i < 4; i++) {
            if (i != k) {
                double f = M[i][k];
                for (int j = k; j < 5; j++) M[i][j] -= f * M[k][j];
            }
        }
    }
    for (int i = 0; i < 4; i++) x[i] = M[i][4];
}

int main() {
    clock_t t_start = clock();
    Solver3D s24;
    init_solver(&s24, 24, 0.05);
    int N3 = s24.N3;
    WorkBuffers wb;
    alloc_work(&wb, N3);

    double complex *k1 = calloc(3*N3, sizeof(double complex));
    double complex *k2 = calloc(3*N3, sizeof(double complex));
    double complex *k3 = calloc(3*N3, sizeof(double complex));
    double complex *k4 = calloc(3*N3, sizeof(double complex));
    double complex *htmp = calloc(3*N3, sizeof(double complex));

    // =========================================================================
    // PART 3: Matched-State Counterfactual Twin Suite (9 Twin Pairs = 18 runs)
    // =========================================================================
    double complex *basis_h[5];
    for (int b = 0; b < 5; b++) {
        basis_h[b] = calloc(3*N3, sizeof(double complex));
        for (int i = 0; i < N3; i++) {
            double x = s24.x[i], y = s24.y[i], z = s24.z[i];
            if (b == 0) { wb.u[i] = sin(y); wb.u[N3+i] = 0.0; wb.u[2*N3+i] = sin(x) - sin(x+y); }
            else if (b == 1) { wb.u[i] = sin(z); wb.u[N3+i] = sin(x+z); wb.u[2*N3+i] = sin(2*x+y); }
            else if (b == 2) { wb.u[i] = cos(y)*sin(z); wb.u[N3+i] = cos(z)*sin(x); wb.u[2*N3+i] = -2.0*sin(x)*sin(y)*sin(z); }
            else if (b == 3) { wb.u[i] = sin(2*y); wb.u[N3+i] = cos(2*z); wb.u[2*N3+i] = cos(2*x); }
            else { wb.u[i] = sin(3*z + y); wb.u[N3+i] = cos(3*x + z); wb.u[2*N3+i] = sin(2*x + 3*y); } // High-k shell mode
        }
        ft_vec3(&s24, wb.u, basis_h[b], wb.tmp1, wb.b1, wb.b2);
        project_hat(&s24, basis_h[b]);
        ift_vec3(&s24, basis_h[b], wb.u, wb.tmp1, wb.b1, wb.b2);
        double nrm2 = 0.0;
        for (int i = 0; i < 3*N3; i++) nrm2 += wb.u[i] * wb.u[i];
        double nrm = sqrt(nrm2 / N3);
        for (int i = 0; i < 3*N3; i++) basis_h[b][i] /= nrm;
    }

    FILE *f_twins = fopen("ns_research/ns_exp014_twins.csv", "w");
    fprintf(f_twins, "t_base,amp_high,max_state_diff,Pplus_0,Q_0,eta_plus_0,J_0,Qdot_0,H_0,dH_A,dH_B,S_pou_A,S_pou_B,pred_tau_A,pred_tau_B,true_tau_A,true_tau_B\n");

    double t_bases[3] = {0.35, 0.40, 0.45};
    double amps[3] = {0.025, 0.040, 0.055};
    double complex *h_ref = calloc(3*N3, sizeof(double complex));
    double complex *h_cand = calloc(3*N3, sizeof(double complex));
    double complex *hA = calloc(3*N3, sizeof(double complex));
    double complex *hB = calloc(3*N3, sizeof(double complex));

    for (int tb_i = 0; tb_i < 3; tb_i++) {
        double t_base = t_bases[tb_i];
        // Initialize with Exp013 3D perturbation (eps=0.1)
        for (int i = 0; i < N3; i++) {
            double x = s24.x[i], y = s24.y[i], z = s24.z[i];
            wb.u[i]      = sin(y) + 0.1 * sin(z);
            wb.u[N3+i]   = 0.1 * sin(x + z);
            wb.u[2*N3+i] = sin(x) - sin(x + y) + 0.1 * sin(2.0*x + y);
        }
        ft_vec3(&s24, wb.u, h_ref, wb.tmp1, wb.b1, wb.b2);
        project_hat(&s24, h_ref);
        int n_pre = (int)round(t_base / 0.005);
        for (int st = 0; st < n_pre; st++) step_rk4(&s24, &wb, h_ref, 0.005, k1, k2, k3, k4, htmp);

        SnapDiag d_ref;
        compute_diagnostics(&s24, &wb, h_ref, t_base, 0, 0, &d_ref);
        double target4[4] = {d_ref.Pplus, d_ref.Q, d_ref.J, d_ref.Qdot};

        for (int am_i = 0; am_i < 3; am_i++) {
            double amp = amps[am_i];
            for (int side = 0; side < 2; side++) {
                double a_high = (side == 0) ? (-amp) : (+amp);
                double a4[4] = {0, 0, 0, 0};
                double eps_fd = 1e-6;
                for (int it = 0; it < 7; it++) {
                    for (int i = 0; i < 3*N3; i++) {
                        h_cand[i] = h_ref[i] + a4[0]*basis_h[0][i] + a4[1]*basis_h[1][i] + a4[2]*basis_h[2][i] + a4[3]*basis_h[3][i] + a_high*basis_h[4][i];
                    }
                    SnapDiag dc;
                    compute_diagnostics(&s24, &wb, h_cand, 0.0, 0, 0, &dc);
                    double res[4] = {dc.Pplus - target4[0], dc.Q - target4[1], dc.J - target4[2], dc.Qdot - target4[3]};
                    double max_r = fmax(fmax(fabs(res[0]), fabs(res[1])), fmax(fabs(res[2]), fabs(res[3])));
                    if (max_r < 1e-12) break;
                    double Jac[4][4];
                    for (int j = 0; j < 4; j++) {
                        for (int i = 0; i < 3*N3; i++) htmp[i] = h_cand[i] + eps_fd * basis_h[j][i];
                        SnapDiag dj;
                        compute_diagnostics(&s24, &wb, htmp, 0.0, 0, 0, &dj);
                        Jac[0][j] = (dj.Pplus - dc.Pplus) / eps_fd;
                        Jac[1][j] = (dj.Q     - dc.Q)     / eps_fd;
                        Jac[2][j] = (dj.J     - dc.J)     / eps_fd;
                        Jac[3][j] = (dj.Qdot  - dc.Qdot)  / eps_fd;
                    }
                    double neg_res[4] = {-res[0], -res[1], -res[2], -res[3]};
                    double da[4];
                    solve_4x4(Jac, neg_res, da);
                    for (int j = 0; j < 4; j++) a4[j] += da[j];
                }
                double complex *dst = (side == 0) ? hA : hB;
                for (int i = 0; i < 3*N3; i++) {
                    dst[i] = h_ref[i] + a4[0]*basis_h[0][i] + a4[1]*basis_h[1][i] + a4[2]*basis_h[2][i] + a4[3]*basis_h[3][i] + a_high*basis_h[4][i];
                }
            }
            SnapDiag dA, dB;
            compute_diagnostics(&s24, &wb, hA, 0.0, 1, 1, &dA);
            compute_diagnostics(&s24, &wb, hB, 0.0, 1, 1, &dB);
            double max_diff = fmax(fmax(fabs(dA.Pplus - dB.Pplus), fabs(dA.Q - dB.Q)),
                                   fmax(fabs(dA.J - dB.J), fmax(fabs(dA.Qdot - dB.Qdot), fabs(dA.H - dB.H))));

            // Evolve hA and hB to find exact crossing time tau_turn
            double tauA = -1.0, tauB = -1.0;
            double prev_HA = dA.H, prev_HB = dB.H;
            double dt_ev = 0.005;
            for (int st = 1; st <= 90; st++) {
                step_rk4(&s24, &wb, hA, dt_ev, k1, k2, k3, k4, htmp);
                SnapDiag sA; compute_diagnostics(&s24, &wb, hA, st*dt_ev, 0, 0, &sA);
                if (tauA < 0.0 && prev_HA >= 0.0 && sA.H < 0.0) {
                    tauA = (st - 1)*dt_ev + (prev_HA / (prev_HA - sA.H)) * dt_ev;
                }
                prev_HA = sA.H;

                step_rk4(&s24, &wb, hB, dt_ev, k1, k2, k3, k4, htmp);
                SnapDiag sB; compute_diagnostics(&s24, &wb, hB, st*dt_ev, 0, 0, &sB);
                if (tauB < 0.0 && prev_HB >= 0.0 && sB.H < 0.0) {
                    tauB = (st - 1)*dt_ev + (prev_HB / (prev_HB - sB.H)) * dt_ev;
                }
                prev_HB = sB.H;
                if (tauA > 0.0 && tauB > 0.0) break;
            }
            fprintf(f_twins, "%.2f,%.3f,%.3e,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f,%.5f,%.5f,%.5f,%.5f\n",
                    t_base, amp, max_diff, dA.Pplus, dA.Q, dA.eta_plus, dA.J, dA.Qdot, dA.H,
                    dA.dH, dB.dH, dA.S_pou, dB.S_pou, -dA.H/dA.dH, -dB.H/dB.dH, tauA, tauB);
            printf("Twin (t0=%.2f, amp=%.3f): max_diff=%.2e | dH_A=%.4f, dH_B=%.4f | pred_dtau=%+.4f vs true_dtau=%+.4f\n",
                   t_base, amp, max_diff, dA.dH, dB.dH, (-dB.H/dB.dH) - (-dA.H/dA.dH), tauB - tauA);
        }
    }
    fclose(f_twins);

    // =========================================================================
    // PART 4: 60-Trajectory 5-Family Ensemble
    // =========================================================================
    const char *fam_names[5] = {"exp013_3d", "tg3d_pert", "antiparallel", "kida_elliptic", "random_solenoidal"};
    double nus[3] = {0.05, 0.02, 0.01};
    double amps_ens[3] = {1.3, 1.15, 1.0};

    FILE *f_man = fopen("ns_research/ns_checkpoint3_4_reproducibility_manifest.csv", "w");
    fprintf(f_man, "traj_id,group_id,family,variant,nu,amp,sign,rng_seed,Pplus_0,Q_0,eta_plus_0,H_0,dH_0,k_char_0,first_crossing_t\n");

    FILE *f_snaps = fopen("ns_research/ns_exp014_snapshots.csv", "w");
    fprintf(f_snaps, "traj_id,group_id,family,nu,sign,snap_idx,t,Pplus,P,Pminus,Q,S_max,W,C_conc,D,eta_plus,k_char,J,R,Qdot,R_Q,H,Delta_triad,Delta_visc,Zdot,T4_B1,T4_B2,T4_B3,S_dyad,T4_P1,T4_P2,T4_P3,S_pou,pou_err,ddP_plus,ddQ,dH,Bdry_P_term,Bdry_Q_term,first_crossing_t\n");

    int traj_id = 0;
    SnapDiag traj_rows[21];
    for (int f_i = 0; f_i < 5; f_i++) {
        for (int var = 0; var < 2; var++) {
            for (int n_i = 0; n_i < 3; n_i++) {
                s24.nu = nus[n_i];
                double amp_e = amps_ens[n_i];
                char gid[64];
                snprintf(gid, sizeof(gid), "%s_v%d_nu%.2f", fam_names[f_i], var, s24.nu);
                for (int s_i = 0; s_i < 2; s_i++) {
                    int sgn = (s_i == 0) ? 1 : -1;
                    make_initial_ensemble(&s24, &wb, f_i, var, amp_e, sgn, h_ref);

                    double first_cross = -1.0;
                    for (int snap = 0; snap <= 20; snap++) {
                        double t_now = snap * 0.05;
                        compute_diagnostics(&s24, &wb, h_ref, t_now, 1, 1, &traj_rows[snap]);
                        if (snap >= 1 && first_cross < 0.0 && traj_rows[snap-1].H > 0.0 && traj_rows[snap].H <= 0.0) {
                            double h0 = traj_rows[snap-1].H, h1 = traj_rows[snap].H;
                            first_cross = traj_rows[snap-1].t + (h0 / (h0 - h1)) * 0.05;
                        }
                        if (snap < 20) {
                            // 5 steps of dt = 0.01 per 0.05 sample
                            for (int st = 0; st < 5; st++) step_rk4(&s24, &wb, h_ref, 0.01, k1, k2, k3, k4, htmp);
                        }
                    }

                    fprintf(f_man, "%d,%s,%s,%d,%.2f,%.2f,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.5f\n",
                            traj_id, gid, fam_names[f_i], var, s24.nu, amp_e, sgn, 202642 + 1337*var,
                            traj_rows[0].Pplus, traj_rows[0].Q, traj_rows[0].eta_plus,
                            traj_rows[0].H, traj_rows[0].dH, traj_rows[0].k_char, first_cross);

                    for (int snap = 0; snap <= 20; snap++) {
                        SnapDiag *r = &traj_rows[snap];
                        fprintf(f_snaps, "%d,%s,%s,%.2f,%d,%d,%.2f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.6f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.3e,%.8f,%.8f,%.8f,%.8f,%.8f,%.5f\n",
                                traj_id, gid, fam_names[f_i], s24.nu, sgn, snap, r->t,
                                r->Pplus, r->P, r->Pminus, r->Q, r->S_max, r->W, r->C_conc, r->D, r->eta_plus, r->k_char,
                                r->J, r->R, r->Qdot, r->R_Q, r->H, r->Delta_triad, r->Delta_visc, r->Zdot,
                                r->T4_B1, r->T4_B2, r->T4_B3, r->S_dyad,
                                r->T4_P1, r->T4_P2, r->T4_P3, r->S_pou, r->pou_err,
                                r->ddP_plus, r->ddQ, r->dH, r->Bdry_P_term, r->Bdry_Q_term, first_cross);
                    }
                    traj_id++;
                }
            }
        }
    }
    fclose(f_man);
    fclose(f_snaps);
    double elapsed = (double)(clock() - t_start) / CLOCKS_PER_SEC;
    printf("Completed all 9 Twin Pairs and 60 3D Trajectories (1,260 snapshots) in %.2f seconds!\n", elapsed);
    return 0;
}

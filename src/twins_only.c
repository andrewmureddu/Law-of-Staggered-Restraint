#define main main_unused
#include "ns_exp014_main.c"
#undef main

void eval_q_normed(Solver3D *s, WorkBuffers *wb, const double complex *h_ref,
                   double complex *basis0, double complex *basis1, double complex *basis_hi,
                   double a0, double a1, double a_hi, double target_Q,
                   double complex *h_out, SnapDiag *d_out, int full_flag) {
    int ntot = 3 * s->N3;
    for (int i = 0; i < ntot; i++) {
        h_out[i] = h_ref[i] + a0 * basis0[i] + a1 * basis1[i] + a_hi * basis_hi[i];
    }
    compute_diagnostics(s, wb, h_out, 0.0, 0, 0, d_out);
    double lam = cbrt(target_Q / d_out->Q);
    for (int i = 0; i < ntot; i++) h_out[i] *= lam;
    compute_diagnostics(s, wb, h_out, 0.0, full_flag, full_flag, d_out);
}

int main() {
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

    double complex *basis_h[4];
    for (int b = 0; b < 4; b++) {
        basis_h[b] = calloc(3*N3, sizeof(double complex));
        for (int i = 0; i < N3; i++) {
            double x = s24.x[i], y = s24.y[i], z = s24.z[i];
            if (b == 0)      { wb.u[i] = cos(y); wb.u[N3+i] = sin(x); wb.u[2*N3+i] = cos(x - y); }
            else if (b == 1) { wb.u[i] = sin(2*y); wb.u[N3+i] = cos(2*x); wb.u[2*N3+i] = sin(x + y); }
            else if (b == 2) { wb.u[i] = sin(z); wb.u[N3+i] = sin(x+z); wb.u[2*N3+i] = sin(2*x+y); } // Twin A: low/mid shell
            else             { wb.u[i] = sin(3*z + y); wb.u[N3+i] = cos(3*x + z); wb.u[2*N3+i] = sin(2*x + 3*y); } // Twin B: high-k shell
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
    double amps[3] = {0.015, 0.022, 0.030};
    double complex *h_ref = calloc(3*N3, sizeof(double complex));
    double complex *hA = calloc(3*N3, sizeof(double complex));
    double complex *hB = calloc(3*N3, sizeof(double complex));

    for (int tb_i = 0; tb_i < 3; tb_i++) {
        double t_base = t_bases[tb_i];
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
        double target_Q = d_ref.Q;
        double target_eta = d_ref.eta_plus;
        double target_H = d_ref.H;

        for (int am_i = 0; am_i < 3; am_i++) {
            double amp = amps[am_i];
            for (int side = 0; side < 2; side++) {
                double complex *b_hi = (side == 0) ? basis_h[2] : basis_h[3];
                double a_hi = (side == 0) ? (-amp) : (+amp);
                double a0 = 0.0, a1 = 0.0;
                // Larger finite-difference step eps_fd = 2e-4 smooths over grid indicator kinks, then shrinks to 1e-7
                double complex *dst = (side == 0) ? hA : hB;
                SnapDiag dc, d0, d1, d_try;
                for (int it = 0; it < 18; it++) {
                    double eps_fd = (it < 6) ? 2e-4 : ((it < 12) ? 5e-6 : 1e-7);
                    eval_q_normed(&s24, &wb, h_ref, basis_h[0], basis_h[1], b_hi, a0, a1, a_hi, target_Q, dst, &dc, 0);
                    double r0 = dc.eta_plus - target_eta;
                    double r1 = dc.H - target_H;
                    double err_cur = hypot(r0, 0.05 * r1);
                    if (fmax(fabs(r0), fabs(r1)) < 1e-13) break;

                    eval_q_normed(&s24, &wb, h_ref, basis_h[0], basis_h[1], b_hi, a0 + eps_fd, a1, a_hi, target_Q, htmp, &d0, 0);
                    eval_q_normed(&s24, &wb, h_ref, basis_h[0], basis_h[1], b_hi, a0, a1 + eps_fd, a_hi, target_Q, htmp, &d1, 0);

                    double j00 = (d0.eta_plus - dc.eta_plus) / eps_fd;
                    double j10 = (d0.H        - dc.H)        / eps_fd;
                    double j01 = (d1.eta_plus - dc.eta_plus) / eps_fd;
                    double j11 = (d1.H        - dc.H)        / eps_fd;
                    double det = j00 * j11 - j01 * j10;
                    if (fabs(det) < 1e-12) break;
                    double da0 = (-r0 * j11 + r1 * j01) / det;
                    double da1 = (-j00 * r1 + j10 * r0) / det;

                    double step = 1.0;
                    int accepted = 0;
                    for (int ls = 0; ls < 8; ls++) {
                        eval_q_normed(&s24, &wb, h_ref, basis_h[0], basis_h[1], b_hi, a0 + step*da0, a1 + step*da1, a_hi, target_Q, htmp, &d_try, 0);
                        double e_try = hypot(d_try.eta_plus - target_eta, 0.05 * (d_try.H - target_H));
                        if (e_try < err_cur) { accepted = 1; break; }
                        step *= 0.5;
                    }
                    if (accepted) {
                        a0 += step * da0;
                        a1 += step * da1;
                    }
                }
                eval_q_normed(&s24, &wb, h_ref, basis_h[0], basis_h[1], b_hi, a0, a1, a_hi, target_Q, dst, &dc, 1);
            }
            SnapDiag dA, dB;
            compute_diagnostics(&s24, &wb, hA, 0.0, 1, 1, &dA);
            compute_diagnostics(&s24, &wb, hB, 0.0, 1, 1, &dB);
            double max_diff = fmax(fmax(fabs(dA.Pplus - dB.Pplus), fabs(dA.Q - dB.Q)),
                                   fmax(fabs(dA.eta_plus - dB.eta_plus), fabs(dA.H - dB.H)));

            double tauA = -1.0, tauB = -1.0;
            double prev_HA = dA.H, prev_HB = dB.H;
            double dt_ev = 0.005;
            for (int st = 1; st <= 80; st++) {
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
            fflush(f_twins);
            printf("Polished Twin (t0=%.2f, amp=%.3f): max_diff=%.2e | dH_A=%.4f, dH_B=%.4f | pred_dtau=%+.4f vs true_dtau=%+.4f\n",
                   t_base, amp, max_diff, dA.dH, dB.dH, (-dB.H/dB.dH) - (-dA.H/dA.dH), tauB - tauA);
            fflush(stdout);
        }
    }
    fclose(f_twins);
    return 0;
}

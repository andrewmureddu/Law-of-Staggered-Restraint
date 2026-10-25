#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>

#define PI 3.14159265358979323846

// Fast 3D DFT for N=24 using precomputed Fourier tables (Separable 3-pass 1D DFT)
typedef struct {
    int N;
    int N3;
    double nu;
    double complex *F_fwd;
    double complex *F_inv;
    double *kx, *ky, *kz, *k2, *kmag;
    int *mask;
    int *mask_B1, *mask_B2, *mask_B3;
    double *x, *y, *z;
} Solver3D;

void init_solver(Solver3D *s, int N, double nu) {
    s->N = N;
    s->N3 = N * N * N;
    s->nu = nu;
    s->F_fwd = (double complex*)malloc(N * N * sizeof(double complex));
    s->F_inv = (double complex*)malloc(N * N * sizeof(double complex));
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
            double angle = 2.0 * PI * j * k / (double)N;
            s->F_fwd[j*N + k] = cos(angle) - I * sin(angle);
            s->F_inv[j*N + k] = (cos(angle) + I * sin(angle)) / (double)N;
        }
    }
    int N3 = s->N3;
    s->kx = (double*)malloc(N3 * sizeof(double));
    s->ky = (double*)malloc(N3 * sizeof(double));
    s->kz = (double*)malloc(N3 * sizeof(double));
    s->k2 = (double*)malloc(N3 * sizeof(double));
    s->kmag = (double*)malloc(N3 * sizeof(double));
    s->mask = (int*)malloc(N3 * sizeof(int));
    s->mask_B1 = (int*)malloc(N3 * sizeof(int));
    s->mask_B2 = (int*)malloc(N3 * sizeof(int));
    s->mask_B3 = (int*)malloc(N3 * sizeof(int));
    s->x = (double*)malloc(N3 * sizeof(double));
    s->y = (double*)malloc(N3 * sizeof(double));
    s->z = (double*)malloc(N3 * sizeof(double));

    double kcut = (double)N / 3.0;
    for (int ix = 0; ix < N; ix++) {
        double kx_v = (ix <= N/2) ? (double)ix : (double)(ix - N);
        for (int iy = 0; iy < N; iy++) {
            double ky_v = (iy <= N/2) ? (double)iy : (double)(iy - N);
            for (int iz = 0; iz < N; iz++) {
                double kz_v = (iz <= N/2) ? (double)iz : (double)(iz - N);
                int idx = (ix * N + iy) * N + iz;
                s->kx[idx] = kx_v;
                s->ky[idx] = ky_v;
                s->kz[idx] = kz_v;
                double k2_v = kx_v*kx_v + ky_v*ky_v + kz_v*kz_v;
                s->k2[idx] = k2_v;
                s->kmag[idx] = sqrt(k2_v);
                int m = (fabs(kx_v) < kcut) && (fabs(ky_v) < kcut) && (fabs(kz_v) < kcut);
                s->mask[idx] = m;
                s->mask_B1[idx] = m && (k2_v >= 1.0 && k2_v <= 2.0);
                s->mask_B2[idx] = m && (k2_v >= 3.0 && k2_v <= 6.0);
                s->mask_B3[idx] = m && (k2_v >= 7.0);
                s->x[idx] = 2.0 * PI * ix / (double)N;
                s->y[idx] = 2.0 * PI * iy / (double)N;
                s->z[idx] = 2.0 * PI * iz / (double)N;
            }
        }
    }
}

void fft3d_1comp(Solver3D *s, const double complex *in, double complex *out, const double complex *M,
                 double complex *buf1, double complex *buf2) {
    int N = s->N;
    // Pass 1: along z
    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            const double complex *row = in + (ix*N + iy)*N;
            double complex *dst = buf1 + (ix*N + iy)*N;
            for (int kz = 0; kz < N; kz++) {
                double complex sum = 0.0;
                const double complex *m_row = M + kz*N;
                for (int iz = 0; iz < N; iz++) sum += m_row[iz] * row[iz];
                dst[kz] = sum;
            }
        }
    }
    // Pass 2: along y
    for (int ix = 0; ix < N; ix++) {
        for (int ky = 0; ky < N; ky++) {
            const double complex *m_row = M + ky*N;
            for (int kz = 0; kz < N; kz++) {
                double complex sum = 0.0;
                for (int iy = 0; iy < N; iy++) sum += m_row[iy] * buf1[(ix*N + iy)*N + kz];
                buf2[(ix*N + ky)*N + kz] = sum;
            }
        }
    }
    // Pass 3: along x
    for (int kx = 0; kx < N; kx++) {
        const double complex *m_row = M + kx*N;
        for (int ky = 0; ky < N; ky++) {
            for (int kz = 0; kz < N; kz++) {
                double complex sum = 0.0;
                for (int ix = 0; ix < N; ix++) sum += m_row[ix] * buf2[(ix*N + ky)*N + kz];
                out[(kx*N + ky)*N + kz] = sum;
            }
        }
    }
}

void ft_vec3(Solver3D *s, const double *u_real, double complex *u_hat,
             double complex *tmp_in, double complex *b1, double complex *b2) {
    int N3 = s->N3;
    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < N3; i++) tmp_in[i] = u_real[c*N3 + i];
        fft3d_1comp(s, tmp_in, u_hat + c*N3, s->F_fwd, b1, b2);
    }
}

void ift_vec3(Solver3D *s, const double complex *u_hat, double *u_real,
              double complex *tmp_out, double complex *b1, double complex *b2) {
    int N3 = s->N3;
    for (int c = 0; c < 3; c++) {
        fft3d_1comp(s, u_hat + c*N3, tmp_out, s->F_inv, b1, b2);
        for (int i = 0; i < N3; i++) u_real[c*N3 + i] = creal(tmp_out[i]);
    }
}

void project_hat(Solver3D *s, double complex *h) {
    int N3 = s->N3;
    for (int i = 0; i < N3; i++) {
        if (!s->mask[i] || s->k2[i] == 0.0) {
            h[i] = 0.0; h[N3 + i] = 0.0; h[2*N3 + i] = 0.0;
        } else {
            double kx = s->kx[i], ky = s->ky[i], kz = s->kz[i], den = s->k2[i];
            double complex k_dot = kx*h[i] + ky*h[N3 + i] + kz*h[2*N3 + i];
            h[i]      -= kx * k_dot / den;
            h[N3 + i] -= ky * k_dot / den;
            h[2*N3+i] -= kz * k_dot / den;
        }
    }
}

void curl_hat(Solver3D *s, const double complex *h, double complex *wh) {
    int N3 = s->N3;
    for (int i = 0; i < N3; i++) {
        double kx = s->kx[i], ky = s->ky[i], kz = s->kz[i];
        double complex hx = h[i], hy = h[N3 + i], hz = h[2*N3 + i];
        wh[i]      = I * (ky * hz - kz * hy);
        wh[N3 + i] = I * (kz * hx - kx * hz);
        wh[2*N3+i] = I * (kx * hy - ky * hx);
    }
}

void strain_from_hat(Solver3D *s, const double complex *h, double *S,
                     double complex *sh, double complex *so, double complex *b1, double complex *b2) {
    // S has 6 independent components: 0:xx, 1:xy, 2:xz, 3:yy, 4:yz, 5:zz
    int N3 = s->N3;
    const double *K[3] = {s->kx, s->ky, s->kz};
    int pairs[6][2] = {{0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2}};
    for (int p = 0; p < 6; p++) {
        int a = pairs[p][0], b = pairs[p][1];
        for (int i = 0; i < N3; i++) {
            sh[i] = 0.5 * I * (K[b][i] * h[a*N3 + i] + K[a][i] * h[b*N3 + i]);
        }
        fft3d_1comp(s, sh, so, s->F_inv, b1, b2);
        for (int i = 0; i < N3; i++) S[p*N3 + i] = creal(so[i]);
    }
}

// Compute eigenvalues and 1st/2nd derivatives of ||S||_op at single grid point
static inline void point_strain_rates(
    const double S6[6], const double SdNL6[6], const double SdV6[6], const double Sdd6[6],
    int do_2nd,
    double *snorm_out, double *rate_NL_out, double *rate_V_out,
    double *dS_op_out, double *ddS_op_out, double *lam2_out, double *dlam2_out
) {
    double tr = (S6[0] + S6[3] + S6[5]) / 3.0;
    double a00 = S6[0] - tr, a01 = S6[1], a02 = S6[2];
    double a11 = S6[3] - tr, a12 = S6[4], a22 = S6[5] - tr;

    double J2 = 0.5 * (a00*a00 + a11*a11 + a22*a22 + 2.0*(a01*a01 + a02*a02 + a12*a12));
    double q = sqrt(fmax(J2 / 3.0, 0.0));
    double detA = a00*(a11*a22 - a12*a12) - a01*(a01*a22 - a12*a02) + a02*(a01*a12 - a11*a02);
    double denom = 2.0 * q * q * q;
    double r = (denom > 1e-30) ? (detA / denom) : 0.0;
    if (r > 1.0) r = 1.0;
    if (r < -1.0) r = -1.0;
    double phi = acos(r) / 3.0;

    double lam1 = 2.0 * q * cos(phi);
    double lam3 = 2.0 * q * cos(phi + 2.0 * PI / 3.0);
    double lam2 = -lam1 - lam3;

    double abs1 = fabs(lam1), abs2 = fabs(lam2), abs3 = fabs(lam3);
    double lam_max = (abs1 >= abs2 && abs1 >= abs3) ? lam1 : ((abs2 >= abs3) ? lam2 : lam3);
    double snorm = fabs(lam_max);
    double sgn = (lam_max >= 0.0) ? 1.0 : -1.0;

    // A^2 symmetric 6-components
    double a2_00 = a00*a00 + a01*a01 + a02*a02;
    double a2_01 = a00*a01 + a01*a11 + a02*a12;
    double a2_02 = a00*a02 + a01*a12 + a02*a22;
    double a2_11 = a01*a01 + a11*a11 + a12*a12;
    double a2_12 = a01*a02 + a11*a12 + a12*a22;
    double a2_22 = a02*a02 + a12*a12 + a22*a22;

    double sylv_den = 3.0 * lam_max * lam_max - J2;
    if (fabs(sylv_den) < 1e-24) sylv_den = 1e-24;

    double p00 = a2_00 + lam_max*a00, p01 = a2_01 + lam_max*a01, p02 = a2_02 + lam_max*a02;
    double p11 = a2_11 + lam_max*a11, p12 = a2_12 + lam_max*a12, p22 = a2_22 + lam_max*a22;

    double num_NL = p00*SdNL6[0] + p11*SdNL6[3] + p22*SdNL6[5] + 2.0*(p01*SdNL6[1] + p02*SdNL6[2] + p12*SdNL6[4]);
    double num_V  = p00*SdV6[0]  + p11*SdV6[3]  + p22*SdV6[5]  + 2.0*(p01*SdV6[1]  + p02*SdV6[2]  + p12*SdV6[4]);

    *snorm_out = snorm;
    *rate_NL_out = sgn * num_NL / sylv_den;
    *rate_V_out  = sgn * num_V  / sylv_den;

    if (do_2nd) {
        double sd00 = SdNL6[0] - SdV6[0], sd01 = SdNL6[1] - SdV6[1], sd02 = SdNL6[2] - SdV6[2];
        double sd11 = SdNL6[3] - SdV6[3], sd12 = SdNL6[4] - SdV6[4], sd22 = SdNL6[5] - SdV6[5];

        double dJ2 = a00*sd00 + a11*sd11 + a22*sd22 + 2.0*(a01*sd01 + a02*sd02 + a12*sd12);
        double dJ3 = a2_00*sd00 + a2_11*sd11 + a2_22*sd22 + 2.0*(a2_01*sd01 + a2_02*sd02 + a2_12*sd12);
        double dlam_max = (lam_max * dJ2 + dJ3) / sylv_den;
        double sylv_den2 = 3.0 * lam2 * lam2 - J2;
        if (fabs(sylv_den2) < 1e-24) sylv_den2 = 1e-24;
        double dlam2 = (lam2 * dJ2 + dJ3) / sylv_den2;

        // Sdot^2 symmetric 6-components
        double sd2_00 = sd00*sd00 + sd01*sd01 + sd02*sd02;
        double sd2_01 = sd00*sd01 + sd01*sd11 + sd02*sd12;
        double sd2_02 = sd00*sd02 + sd01*sd12 + sd02*sd22;
        double sd2_11 = sd01*sd01 + sd11*sd11 + sd12*sd12;
        double sd2_12 = sd01*sd02 + sd11*sd12 + sd12*sd22;
        double sd2_22 = sd02*sd02 + sd12*sd12 + sd22*sd22;

        double tr_sd2 = sd2_00 + sd2_11 + sd2_22;
        double tr_A_sdd = a00*Sdd6[0] + a11*Sdd6[3] + a22*Sdd6[5] + 2.0*(a01*Sdd6[1] + a02*Sdd6[2] + a12*Sdd6[4]);
        double ddJ2 = tr_sd2 + tr_A_sdd;

        double tr_A_sd2 = a00*sd2_00 + a11*sd2_11 + a22*sd2_22 + 2.0*(a01*sd2_01 + a02*sd2_02 + a12*sd2_12);
        double tr_A2_sdd = a2_00*Sdd6[0] + a2_11*Sdd6[3] + a2_22*Sdd6[5] + 2.0*(a2_01*Sdd6[1] + a2_02*Sdd6[2] + a2_12*Sdd6[4]);
        double ddJ3 = 2.0 * tr_A_sd2 + tr_A2_sdd;

        double ddlam_max = (lam_max * ddJ2 + ddJ3 + 2.0 * dJ2 * dlam_max - 6.0 * lam_max * dlam_max * dlam_max) / sylv_den;
        *dS_op_out = sgn * dlam_max;
        *ddS_op_out = sgn * ddlam_max;
        *lam2_out = lam2;
        *dlam2_out = dlam2;
    }
}

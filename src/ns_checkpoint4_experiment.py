"""Navier-Stokes Checkpoint IV (Experiment 014):
The Vortex That Knew It Was About to Change Its Mind.
"""
import csv
import json
import math
import time
from pathlib import Path
import numpy as np
from numpy.fft import fftn, ifftn

AXES = (-3, -2, -1)


def verify_sylvester_identities(n_trials=1000, eps=1e-4, seed=20261023):
    rng = np.random.default_rng(seed)
    errs_1, errs_2 = [], []
    while len(errs_2) < n_trials:
        A0, A1, A2 = rng.standard_normal((3, 3, 3))
        S0 = 0.5 * (A0 + A0.T); S0 -= np.trace(S0) / 3.0 * np.eye(3)
        S1 = 0.5 * (A1 + A1.T); S1 -= np.trace(S1) / 3.0 * np.eye(3)
        S2 = 0.5 * (A2 + A2.T); S2 -= np.trace(S2) / 3.0 * np.eye(3)

        evals = np.linalg.eigvalsh(S0)
        if abs(abs(evals[0]) - abs(evals[2])) < 2e-2:
            continue
        lam = evals[2] if abs(evals[2]) >= abs(evals[0]) else evals[0]
        sgn = 1.0 if lam >= 0 else -1.0

        J2 = 0.5 * np.trace(S0 @ S0)
        denom = 3.0 * lam**2 - J2

        dJ2 = np.trace(S0 @ S1)
        dJ3 = np.trace((S0 @ S0) @ S1)
        dlam = (lam * dJ2 + dJ3) / denom

        ddJ2 = np.trace(S1 @ S1 + S0 @ S2)
        ddJ3 = np.trace(2.0 * (S0 @ S1 @ S1) + (S0 @ S0) @ S2)
        ddlam = (lam * ddJ2 + ddJ3 + 2.0 * dJ2 * dlam - 6.0 * lam * (dlam**2)) / denom

        def norm_op(t):
            St = S0 + t * S1 + 0.5 * (t**2) * S2
            ev = np.linalg.eigvalsh(St)
            return max(abs(ev[0]), abs(ev[2]))

        f_m2, f_m1, f_0, f_p1, f_p2 = norm_op(-2*eps), norm_op(-eps), norm_op(0.0), norm_op(eps), norm_op(2*eps)
        fd1 = (f_m2 - 8*f_m1 + 8*f_p1 - f_p2) / (12.0 * eps)
        fd2 = (-f_m2 + 16*f_m1 - 30*f_0 + 16*f_p1 - f_p2) / (12.0 * eps**2)

        errs_1.append(abs(sgn * dlam - fd1) / max(abs(fd1), 1.0))
        errs_2.append(abs(sgn * ddlam - fd2) / max(abs(fd2), 1.0))

    return {
        "n_trials": len(errs_2),
        "order1_max_err": float(np.max(errs_1)),
        "order1_median_err": float(np.median(errs_1)),
        "order2_max_err": float(np.max(errs_2)),
        "order2_median_err": float(np.median(errs_2)),
    }


def strain_eig_and_1st_2nd_rates(S, Sdot_NL, Sdot_visc, Sddot=None):
    tr = (S[0, 0] + S[1, 1] + S[2, 2]) / 3.0
    A = S - tr * np.eye(3)[:, :, None, None, None]
    J2 = 0.5 * np.einsum("ijabc,ijabc->abc", A, A)
    q = np.sqrt(np.maximum(J2 / 3.0, 0.0))

    detA = (
        A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] ** 2)
        - A[0, 1] * (A[0, 1] * A[2, 2] - A[1, 2] * A[0, 2])
        + A[0, 2] * (A[0, 1] * A[1, 2] - A[1, 1] * A[0, 2])
    )
    denom = 2.0 * (q ** 3)
    r = np.clip(np.divide(detA, denom, out=np.zeros_like(detA), where=denom > 1e-30), -1.0, 1.0)
    phi = np.arccos(r) / 3.0

    lam1 = 2.0 * q * np.cos(phi)
    lam3 = 2.0 * q * np.cos(phi + (2.0 * np.pi / 3.0))
    lam2 = -lam1 - lam3

    abs1, abs2, abs3 = np.abs(lam1), np.abs(lam2), np.abs(lam3)
    lam_max = np.where((abs1 >= abs2) & (abs1 >= abs3), lam1, np.where(abs2 >= abs3, lam2, lam3))
    snorm = np.abs(lam_max)
    sgn = np.where(lam_max >= 0.0, 1.0, -1.0)

    A2 = np.einsum("ikabc,kjabc->ijabc", A, A)
    sylv_denom = np.where(np.abs(3.0 * (lam_max ** 2) - J2) > 1e-24, 3.0 * (lam_max ** 2) - J2, 1e-24)
    sylv_denom2 = np.where(np.abs(3.0 * (lam2 ** 2) - J2) > 1e-24, 3.0 * (lam2 ** 2) - J2, 1e-24)

    P_num = A2 + lam_max[None, None, ...] * A
    num_NL = np.einsum("ijabc,ijabc->abc", P_num, Sdot_NL)
    num_visc = np.einsum("ijabc,ijabc->abc", P_num, Sdot_visc)

    rate_NL = sgn * (num_NL / sylv_denom)
    rate_visc = sgn * (num_visc / sylv_denom)

    if Sddot is None:
        return snorm, rate_NL, rate_visc, None, None, None

    Sdot = Sdot_NL - Sdot_visc
    dJ2 = np.einsum("ijabc,ijabc->abc", A, Sdot)
    dJ3 = np.einsum("ijabc,ijabc->abc", A2, Sdot)
    dlam_max = (lam_max * dJ2 + dJ3) / sylv_denom
    dlam2 = (lam2 * dJ2 + dJ3) / sylv_denom2
    dS_op = sgn * dlam_max

    Sdot2 = np.einsum("ikabc,kjabc->ijabc", Sdot, Sdot)
    ddJ2 = np.einsum("ijabc,ijabc->abc", Sdot, Sdot) + np.einsum("ijabc,ijabc->abc", A, Sddot)
    ddJ3 = 2.0 * np.einsum("ijabc,ijabc->abc", A, Sdot2) + np.einsum("ijabc,ijabc->abc", A2, Sddot)
    ddlam_max = (lam_max * ddJ2 + ddJ3 + 2.0 * dJ2 * dlam_max - 6.0 * lam_max * (dlam_max ** 2)) / sylv_denom
    ddS_op = sgn * ddlam_max

    return snorm, rate_NL, rate_visc, dS_op, ddS_op, (lam2, dlam2)


class Checkpoint4Solver:
    def __init__(self, n: int, nu: float):
        self.n, self.nu = n, nu
        k = np.fft.fftfreq(n, 1.0 / n)
        self.k = np.array(np.meshgrid(k, k, k, indexing="ij"))
        self.k2 = (self.k ** 2).sum(axis=0)
        self.kmag = np.sqrt(self.k2)
        self.den = np.where(self.k2 == 0, 1.0, self.k2)
        self.mask = (np.abs(self.k) < n / 3.0).all(axis=0)

        self.mask_B1 = self.mask & (self.k2 >= 1) & (self.k2 <= 2)
        self.mask_B2 = self.mask & (self.k2 >= 3) & (self.k2 <= 6)
        self.mask_B3 = self.mask & (self.k2 >= 7)

        self.xyz = np.meshgrid(
            np.arange(n) * 2 * np.pi / n,
            np.arange(n) * 2 * np.pi / n,
            np.arange(n) * 2 * np.pi / n,
            indexing="ij",
        )

    def ft(self, a):
        return fftn(a, axes=AXES)

    def real(self, a):
        return ifftn(a, axes=AXES).real

    def project(self, h):
        return (h - self.k * (self.k * h).sum(axis=0) / self.den) * self.mask

    def curl(self, h):
        return 1j * np.cross(self.k, h, axisa=0, axisb=0, axisc=0)

    def strain_from_hat(self, h):
        grad = self.real(1j * h[:, None] * self.k[None, :])
        return 0.5 * (grad + np.swapaxes(grad, 0, 1))

    def initial_exp013_3d_pert(self, eps: float = 0.1, sign: int = 1):
        x, y, z = self.xyz
        u = np.array([
            np.sin(y) + eps * np.sin(z),
            eps * np.sin(x + z),
            np.sin(x) - np.sin(x + y) + eps * np.sin(2 * x + y),
        ])
        return self.project(self.ft(sign * u))

    def initial_ensemble(self, family: str, variant: int = 0, amp: float = 1.0, sign: int = 1, rng_seed: int = 42):
        x, y, z = self.xyz
        if family == "exp013_3d":
            eps = 0.10 + 0.08 * variant
            u = np.array([
                np.sin(y) + eps * np.sin(z),
                eps * np.sin(x + z),
                np.sin(x) - np.sin(x + y) + eps * np.sin(2 * x + y),
            ])
        elif family == "tg3d_pert":
            delta = 0.22 * (variant + 1)
            u = np.array([
                np.sin(x) * np.cos(y) * np.cos(z) + delta * np.sin(2 * y) * np.cos(z),
                -np.cos(x) * np.sin(y) * np.cos(z) - delta * np.cos(x) * np.sin(2 * z),
                delta * np.sin(x) * np.sin(y) * np.sin(2 * z),
            ])
        elif family == "antiparallel":
            alpha = 0.35 + 0.20 * variant
            u = np.array([
                np.sin(y) * np.cos(z) + alpha * np.cos(y + z),
                np.sin(z) * np.cos(x) - alpha * np.cos(x - z),
                -np.cos(x) * np.sin(y) + alpha * np.sin(x + y) * np.cos(z),
            ])
        elif family == "random_solenoidal":
            rng = np.random.default_rng(rng_seed + 100 * variant)
            u = np.zeros((3, self.n, self.n, self.n))
            for kx in range(1, 4):
                for ky in range(1, 4):
                    for kz in range(1, 4):
                        if kx * kx + ky * ky + kz * kz <= 14:
                            ph = rng.uniform(0, 2 * np.pi, size=3)
                            wt = np.exp(-0.25 * (kx * kx + ky * ky + kz * kz))
                            u += wt * np.cos(kx * x + ky * y + kz * z + ph[:, None, None, None])
        elif family == "phase_rand":
            base_h = self.initial_ensemble("tg3d_pert", variant=variant, amp=amp, sign=1)
            rng = np.random.default_rng(rng_seed + 500 * variant)
            noise_ft = fftn(rng.standard_normal((self.n, self.n, self.n)))
            unit_phasor = np.where(np.abs(noise_ft) > 1e-12, noise_ft / np.abs(noise_ft), 1.0 + 0.0j)
            u_real = self.real(base_h * unit_phasor[None, :, :, :])
            return self.project(self.ft(sign * u_real))
        else:
            raise ValueError(family)

        h = self.project(self.ft(u))
        w = self.real(self.curl(h))
        W_raw = float(np.mean((w * w).sum(axis=0)))
        scale = amp * math.sqrt(2.0 / W_raw) if W_raw > 1e-12 else amp
        return self.project(sign * scale * h)

    def rhs(self, h):
        u, w = self.real(h), self.real(self.curl(h))
        nonlinear = self.ft(np.cross(u, w, axisa=0, axisb=0, axisc=0))
        return self.project(nonlinear) - self.nu * self.k2 * h

    def step(self, h, dt):
        a = self.rhs(h)
        b = self.rhs(h + dt * a / 2.0)
        c = self.rhs(h + dt * b / 2.0)
        d = self.rhs(h + dt * c)
        return self.project(h + dt * (a + 2.0 * b + 2.0 * c + d) / 6.0)

    def shell_t4(self, Nu_hat, Nw_hat, weight_mask, w, Sw, pos_mask):
        Nu_b = Nu_hat * weight_mask
        Nw_b = self.real(Nw_hat * weight_mask)
        S_Nu_b = self.strain_from_hat(Nu_b)
        t4_b_field = 2.0 * np.einsum("iabc,iabc->abc", Nw_b, Sw) + np.einsum("iabc,ijabc,jabc->abc", w, S_Nu_b, w)
        return float(np.mean(np.where(pos_mask, t4_b_field, 0.0)))

    def diagnostics(self, h, t: float, compute_second_order=True, full_poly_split=False, compute_shells=True):
        u = self.real(h)
        wh = self.curl(h)
        w = self.real(wh)
        q2 = (w * w).sum(axis=0)
        q = np.sqrt(q2)
        strain = self.strain_from_hat(h)
        Sw = np.einsum("ijabc,jabc->iabc", strain, w)
        p = np.einsum("iabc,iabc->abc", w, Sw)

        gw = self.real(1j * wh[:, None] * self.k[None, :])
        D = self.nu * float(np.mean((gw * gw).sum(axis=(0, 1))))
        E = float(0.5 * np.mean((u * u).sum(axis=0)))
        W = float(np.mean(q2))
        Z = 0.5 * W
        P = float(np.mean(p))
        pos_mask = p > 0.0
        neg_mask = p < 0.0
        Pplus = float(np.mean(np.where(pos_mask, p, 0.0)))
        Pminus = float(np.mean(np.where(neg_mask, -p, 0.0)))

        k_char = math.sqrt(max(D / max(self.nu * W, 1e-30), 1e-12))

        lamb = np.cross(u, w, axisa=0, axisb=0, axisc=0)
        Nu_hat = self.project(self.ft(lamb) * self.mask)
        Nw_hat = self.curl(Nu_hat)
        Nw = self.real(Nw_hat)
        S_Nu = self.strain_from_hat(Nu_hat)

        lap_w = self.real(self.k2 * wh)
        lap_S = self.strain_from_hat(self.k2 * h)

        Sddot = None
        if compute_second_order:
            du_hat = Nu_hat - self.nu * self.k2 * h
            du = self.real(du_hat)
            dw = Nw - self.nu * lap_w
            dS = S_Nu - self.nu * lap_S
            lamb_dot = np.cross(du, w, axisa=0, axisb=0, axisc=0) + np.cross(u, dw, axisa=0, axisb=0, axisc=0)
            ddu_hat = self.project(self.ft(lamb_dot) * self.mask) - self.nu * self.k2 * du_hat
            ddw = self.real(self.curl(ddu_hat))
            Sddot = self.strain_from_hat(ddu_hat)

        snorm, snorm_dot_NL, snorm_dot_visc, dS_op, ddS_op, lam2_pair = strain_eig_and_1st_2nd_rates(
            strain, S_Nu, self.nu * lap_S, Sddot=Sddot
        )
        S_max = float(np.mean(snorm))
        Q = float(np.mean(q2 * snorm))
        eta_plus = float(Pplus / Q) if Q > 1e-15 else 0.0
        C_conc = float(np.mean(q2 * q2) / (W * W)) if W > 1e-15 else 1.0
        q95 = float(np.quantile(q, 0.95))

        t4_Q_field = 2.0 * np.einsum("iabc,iabc->abc", w, Nw) * snorm + q2 * snorm_dot_NL
        v3_Q_field = 2.0 * self.nu * np.einsum("iabc,iabc->abc", w, lap_w) * snorm + q2 * snorm_dot_visc
        T4_Q = float(np.mean(t4_Q_field))
        V3_Q = float(np.mean(v3_Q_field))
        Qdot = T4_Q - V3_Q

        t4_field = 2.0 * np.einsum("iabc,iabc->abc", Nw, Sw) + np.einsum("iabc,ijabc,jabc->abc", w, S_Nu, w)
        v3_field = self.nu * (
            2.0 * np.einsum("iabc,iabc->abc", lap_w, Sw) + np.einsum("iabc,ijabc,jabc->abc", w, lap_S, w)
        )
        T4 = float(np.mean(t4_field))
        V3 = float(np.mean(v3_field))
        T4_plus = float(np.mean(np.where(pos_mask, t4_field, 0.0)))
        V3_plus = float(np.mean(np.where(pos_mask, v3_field, 0.0)))
        Pplus_dot = T4_plus - V3_plus
        Zdot = P - D

        Pplus_safe = max(Pplus, 1e-12)
        Q_safe = max(Q, 1e-12)
        H_exact = (Pplus_dot / Pplus_safe) - (Qdot / Q_safe)
        Delta_triad = (T4_plus / Pplus_safe) - (T4_Q / Q_safe)
        Delta_visc = (V3_plus / Pplus_safe) - (V3_Q / Q_safe)

        out = dict(
            t=float(t), E=E, W=W, Z=Z, D=D, P=P, Pplus=Pplus, Pminus=Pminus,
            Q=Q, S_max=S_max, eta_plus=eta_plus, C_conc=C_conc, q95=q95, k_char=k_char,
            T4=T4, V3=V3, T4_plus=T4_plus, V3_plus=V3_plus,
            J=Pplus_dot, R=T4_plus + V3_plus,
            T4_Q=T4_Q, V3_Q=V3_Q, Qdot=Qdot, R_Q=T4_Q + V3_Q,
            H=H_exact, Delta_triad=Delta_triad, Delta_visc=Delta_visc, Zdot=Zdot,
        )

        if compute_shells:
            T4_B1 = self.shell_t4(Nu_hat, Nw_hat, self.mask_B1, w, Sw, pos_mask)
            T4_B2 = self.shell_t4(Nu_hat, Nw_hat, self.mask_B2, w, Sw, pos_mask)
            T4_B3 = T4_plus - T4_B1 - T4_B2
            S_dyad = (T4_B2 + T4_B3 - T4_B1) / Q_safe

            kmag_safe = np.where(self.kmag == 0, 1.0, self.kmag)
            xi = np.log2(kmag_safe / k_char)
            c1, c2, tau = -0.25, 0.35, 0.25
            w1_pou = (1.0 / (1.0 + np.exp((xi - c1) / tau))) * self.mask
            w3_pou = (1.0 / (1.0 + np.exp(-(xi - c2) / tau))) * self.mask
            w2_pou = (1.0 - (1.0 / (1.0 + np.exp((xi - c1) / tau))) - (1.0 / (1.0 + np.exp(-(xi - c2) / tau)))) * self.mask

            T4_P1 = self.shell_t4(Nu_hat, Nw_hat, w1_pou, w, Sw, pos_mask)
            T4_P3 = self.shell_t4(Nu_hat, Nw_hat, w3_pou, w, Sw, pos_mask)
            T4_P2 = T4_plus - T4_P1 - T4_P3
            pou_partition_err = 0.0
            S_pou = (T4_P2 + T4_P3 - T4_P1) / Pplus_safe

            out.update(
                T4_B1=T4_B1, T4_B2=T4_B2, T4_B3=T4_B3, S_dyad=S_dyad,
                T4_P1=T4_P1, T4_P2=T4_P2, T4_P3=T4_P3, S_pou=S_pou,
                pou_partition_err=pou_partition_err,
            )

        if compute_second_order:
            dp = t4_field - v3_field
            dSw = np.einsum("ijabc,jabc->iabc", dS, w)
            Sdw = np.einsum("ijabc,jabc->iabc", strain, dw)
            ddp = (
                2.0 * np.einsum("iabc,iabc->abc", ddw, Sw)
                + 4.0 * np.einsum("iabc,iabc->abc", dw, dSw)
                + 2.0 * np.einsum("iabc,iabc->abc", dw, Sdw)
                + np.einsum("iabc,ijabc,jabc->abc", w, Sddot, w)
            )
            ddP_bulk = float(np.mean(np.where(pos_mask, ddp, 0.0)))
            eps_p = 0.08 * float(np.std(p)) + 1e-12
            delta_p = np.exp(-((p / eps_p) ** 2)) / (math.sqrt(math.pi) * eps_p)
            ddP_bdry = float(np.mean(delta_p * (dp ** 2)))
            ddP_plus = ddP_bulk + ddP_bdry

            lam2, dlam2 = lam2_pair
            eps_l2 = 0.08 * float(np.std(lam2)) + 1e-12
            delta_l2 = np.exp(-((lam2 / eps_l2) ** 2)) / (math.sqrt(math.pi) * eps_l2)
            dw2 = 2.0 * np.einsum("iabc,iabc->abc", w, dw)
            ddw2 = 2.0 * (dw * dw).sum(axis=0) + 2.0 * np.einsum("iabc,iabc->abc", w, ddw)
            ddQ_bulk = float(np.mean(ddw2 * snorm + 2.0 * dw2 * dS_op + q2 * ddS_op))
            ddQ_bdry = float(np.mean(q2 * delta_l2 * (dlam2 ** 2)))
            ddQ = ddQ_bulk + ddQ_bdry

            dH = (ddP_plus / Pplus_safe) - (Pplus_dot / Pplus_safe) ** 2 - (ddQ / Q_safe) + (Qdot / Q_safe) ** 2

            out.update(
                ddP_bulk=ddP_bulk, ddP_bdry=ddP_bdry, ddP_plus=ddP_plus,
                ddQ_bulk=ddQ_bulk, ddQ_bdry=ddQ_bdry, ddQ=ddQ,
                dH=dH,
                Bdry_P_term=ddP_bdry / Pplus_safe, Bdry_Q_term=ddQ_bdry / Q_safe,
            )

            if full_poly_split:
                def eval_bulk_at_nu(nu_v):
                    du_h = Nu_hat - nu_v * self.k2 * h
                    du_v = self.real(du_h)
                    dw_v = Nw - nu_v * lap_w
                    dS_v = S_Nu - nu_v * lap_S
                    ldot = np.cross(du_v, w, axisa=0, axisb=0, axisc=0) + np.cross(u, dw_v, axisa=0, axisb=0, axisc=0)
                    ddu_h = self.project(self.ft(ldot) * self.mask) - nu_v * self.k2 * du_h
                    ddw_v = self.real(self.curl(ddu_h))
                    ddS_v = self.strain_from_hat(ddu_h)

                    dSw_v = np.einsum("ijabc,jabc->iabc", dS_v, w)
                    Sdw_v = np.einsum("ijabc,jabc->iabc", strain, dw_v)
                    ddp_v = (
                        2.0 * np.einsum("iabc,iabc->abc", ddw_v, Sw)
                        + 4.0 * np.einsum("iabc,iabc->abc", dw_v, dSw_v)
                        + 2.0 * np.einsum("iabc,iabc->abc", dw_v, Sdw_v)
                        + np.einsum("iabc,ijabc,jabc->abc", w, ddS_v, w)
                    )
                    ddP_b = float(np.mean(np.where(pos_mask, ddp_v, 0.0)))
                    _, _, _, dS_op_v, ddS_op_v, _ = strain_eig_and_1st_2nd_rates(
                        strain, S_Nu, nu_v * lap_S, Sddot=ddS_v
                    )
                    dw2_v = 2.0 * np.einsum("iabc,iabc->abc", w, dw_v)
                    ddw2_v = 2.0 * (dw_v * dw_v).sum(axis=0) + 2.0 * np.einsum("iabc,iabc->abc", w, ddw_v)
                    ddQ_b = float(np.mean(ddw2_v * snorm + 2.0 * dw2_v * dS_op_v + q2 * ddS_op_v))
                    return ddP_b, ddQ_b

                ddP_0, ddQ_0 = eval_bulk_at_nu(0.0)
                ddP_m, ddQ_m = eval_bulk_at_nu(-self.nu)
                T5_plus = ddP_0
                M4_plus = -0.5 * (ddP_bulk - ddP_m)
                V32_plus = 0.5 * (ddP_bulk + ddP_m) - ddP_0
                T5_Q = ddQ_0
                M4_Q = -0.5 * (ddQ_bulk - ddQ_m)
                V32_Q = 0.5 * (ddQ_bulk + ddQ_m) - ddQ_0
                out.update(
                    T5_plus=T5_plus, M4_plus=M4_plus, V32_plus=V32_plus,
                    T5_Q=T5_Q, M4_Q=M4_Q, V32_Q=V32_Q,
                    Delta_T5=T5_plus / Pplus_safe - T5_Q / Q_safe,
                    Delta_M4=M4_plus / Pplus_safe - M4_Q / Q_safe,
                    Delta_V32=V32_plus / Pplus_safe - V32_Q / Q_safe,
                )

        return out


def run_traj(solver: Checkpoint4Solver, h0, T: float = 1.0, dt: float = 0.02, sample_dt: float = 0.05, full_poly_split=False):
    h = h0.copy()
    steps = int(round(T / dt))
    sample_every = max(1, int(round(sample_dt / dt)))
    rows = [solver.diagnostics(h, 0.0, compute_second_order=True, full_poly_split=full_poly_split, compute_shells=True)]
    for j in range(1, steps + 1):
        h = solver.step(h, dt)
        if j % sample_every == 0 or j == steps:
            rows.append(solver.diagnostics(h, j * dt, compute_second_order=True, full_poly_split=full_poly_split, compute_shells=True))
    return rows


def fit_predict_ridge(X_train, y_train, X_test, alpha=2e-3):
    mean_x = X_train.mean(axis=0)
    std_x = np.where(X_train.std(axis=0) < 1e-12, 1.0, X_train.std(axis=0))
    Xz_tr = (X_train - mean_x) / std_x
    Xz_te = (X_test - mean_x) / std_x
    mean_y = float(y_train.mean())
    p = Xz_tr.shape[1]
    w = np.linalg.solve(Xz_tr.T @ Xz_tr + alpha * len(y_train) * np.eye(p), Xz_tr.T @ (y_train - mean_y))
    return mean_y + Xz_te @ w


def exact_population_schur_r2(X_base, r_vec, y_vec):
    Xb = np.column_stack([np.ones(len(y_vec)), X_base])
    coef_y, _, _, _ = np.linalg.lstsq(Xb, y_vec, rcond=None)
    coef_r, _, _, _ = np.linalg.lstsq(Xb, r_vec, rcond=None)
    y_perp = y_vec - Xb @ coef_y
    r_perp = r_vec - Xb @ coef_r
    var_y_perp = float(np.mean(y_perp ** 2))
    var_r_perp = float(np.mean(r_perp ** 2))
    cov_yr = float(np.mean(y_perp * r_perp))
    delta_mse = (cov_yr ** 2) / max(var_r_perp, 1e-30)
    partial_r2 = delta_mse / max(var_y_perp, 1e-30)
    return {
        "base_mse": var_y_perp,
        "delta_mse_schur": delta_mse,
        "partial_r2_pct": 100.0 * partial_r2,
        "partial_corr": math.copysign(math.sqrt(max(partial_r2, 0.0)), cov_yr),
    }


def auroc_score(y_true, y_score):
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return 0.5
    diffs = pos[:, None] - neg[None, :]
    return float(np.mean((diffs > 0).astype(np.float64) + 0.5 * (diffs == 0).astype(np.float64)))


def main():
    t0 = time.time()
    out_dir = Path("ns_research")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Verify 1st & 2nd order Sylvester-Cayley-Hamilton identities
    sylv_res = verify_sylvester_identities()
    print("1. Sylvester verification:", json.dumps(sylv_res), flush=True)

    # 2. Exact 2nd-Derivative Crossing Audit on Experiment 013 3D Perturbed Witness (N=32, nu=0.05)
    sol32 = Checkpoint4Solver(32, 0.05)
    rows_3dpert = run_traj(sol32, sol32.initial_exp013_3d_pert(0.1, +1), T=1.0, dt=0.005, sample_dt=0.01, full_poly_split=True)

    # Also verify exact pou_partition_err on t=0.50 snapshot
    h_05 = sol32.initial_exp013_3d_pert(0.1, +1)
    for _ in range(100):
        h_05 = sol32.step(h_05, 0.005)
    u_05 = sol32.real(h_05)
    wh_05 = sol32.curl(h_05)
    w_05 = sol32.real(wh_05)
    strain_05 = sol32.strain_from_hat(h_05)
    Sw_05 = np.einsum("ijabc,jabc->iabc", strain_05, w_05)
    p_05 = np.einsum("iabc,iabc->abc", w_05, Sw_05)
    pos_05 = p_05 > 0.0
    Nu_05 = sol32.project(sol32.ft(np.cross(u_05, w_05, axisa=0, axisb=0, axisc=0)) * sol32.mask)
    Nw_05 = sol32.curl(Nu_05)
    k_char_05 = rows_3dpert[50]["k_char"]
    xi_05 = np.log2(np.where(sol32.kmag == 0, 1.0, sol32.kmag) / k_char_05)
    w1_05 = (1.0 / (1.0 + np.exp((xi_05 + 0.25) / 0.25))) * sol32.mask
    w3_05 = (1.0 / (1.0 + np.exp(-(xi_05 - 0.35) / 0.25))) * sol32.mask
    w2_05 = (1.0 - (1.0 / (1.0 + np.exp((xi_05 + 0.25) / 0.25))) - (1.0 / (1.0 + np.exp(-(xi_05 - 0.35) / 0.25)))) * sol32.mask
    tp1 = sol32.shell_t4(Nu_05, Nw_05, w1_05, w_05, Sw_05, pos_05)
    tp2 = sol32.shell_t4(Nu_05, Nw_05, w2_05, w_05, Sw_05, pos_05)
    tp3 = sol32.shell_t4(Nu_05, Nw_05, w3_05, w_05, Sw_05, pos_05)
    exact_pou_sum_err = abs(rows_3dpert[50]["T4_plus"] - (tp1 + tp2 + tp3))

    crossing_audit = {}
    for i in range(len(rows_3dpert) - 1):
        if rows_3dpert[i]["H"] >= 0.0 and rows_3dpert[i + 1]["H"] < 0.0:
            w0, w1 = rows_3dpert[i], rows_3dpert[i + 1]
            frac = w0["H"] / (w0["H"] - w1["H"])
            t_star = w0["t"] + frac * (w1["t"] - w0["t"])
            fd_ddP = (w1["J"] - w0["J"]) / (w1["t"] - w0["t"])
            fd_ddQ = (w1["Qdot"] - w0["Qdot"]) / (w1["t"] - w0["t"])
            fd_dH = (w1["H"] - w0["H"]) / (w1["t"] - w0["t"])
            interp = lambda k: float((1.0 - frac) * w0[k] + frac * w1[k])
            crossing_audit = {
                "t_star": float(t_star),
                "Pplus": interp("Pplus"), "Q": interp("Q"), "eta_plus": interp("eta_plus"),
                "J": interp("J"), "Qdot": interp("Qdot"),
                "fd_ddP": float(fd_ddP), "fd_ddQ": float(fd_ddQ), "fd_dH": float(fd_dH),
                "ana_ddP_bulk": interp("ddP_bulk"), "ana_ddP_bdry": interp("ddP_bdry"), "ana_ddP": interp("ddP_plus"),
                "ana_ddQ_bulk": interp("ddQ_bulk"), "ana_ddQ_bdry": interp("ddQ_bdry"), "ana_ddQ": interp("ddQ"),
                "ana_dH": interp("dH"),
                "T5_plus": interp("T5_plus"), "M4_plus": interp("M4_plus"), "V32_plus": interp("V32_plus"),
                "T5_Q": interp("T5_Q"), "M4_Q": interp("M4_Q"), "V32_Q": interp("V32_Q"),
                "Delta_T5": interp("Delta_T5"), "Delta_M4": interp("Delta_M4"), "Delta_V32": interp("Delta_V32"),
                "Bdry_P_term": interp("Bdry_P_term"), "Bdry_Q_term": interp("Bdry_Q_term"),
                "exact_pou_sum_err": float(exact_pou_sum_err),
            }
            break
    print(f"2. Crossing audit at t_*={crossing_audit['t_star']:.5f} ({time.time()-t0:.1f}s):", json.dumps(crossing_audit, indent=2), flush=True)

    # 3. Outcome A vs. Outcome B: Matched-State Counterfactual Twin Flows (N=24, nu=0.05)
    sol24_05 = Checkpoint4Solver(24, 0.05)
    x, y, z = sol24_05.xyz
    raw_basis = [
        np.array([np.sin(y), np.zeros_like(x), np.sin(x) - np.sin(x + y)]),
        np.array([np.sin(z), np.sin(x + z), np.sin(2 * x + y)]),
        np.array([np.cos(y) * np.sin(z), np.cos(z) * np.sin(x), -2.0 * np.sin(x) * np.sin(y) * np.sin(z)]),
        np.array([np.sin(2 * y), np.cos(2 * z), np.cos(2 * x)]),
        np.array([np.cos(x + y), -np.cos(x + y), np.sin(x - y) * np.cos(z)]),
        np.array([np.sin(3 * z + y), np.cos(3 * x + z), np.sin(2 * x + 3 * y)]),
    ]
    basis_h = []
    for rb in raw_basis:
        bh = sol24_05.project(sol24_05.ft(rb))
        u_b = sol24_05.real(bh)
        nrm = math.sqrt(float(np.mean((u_b * u_b).sum(axis=0))))
        basis_h.append(bh / max(nrm, 1e-12))

    def eval_state_vec(h_cand):
        d = sol24_05.diagnostics(h_cand, 0.0, compute_second_order=False, compute_shells=False)
        return np.array([d["Pplus"], d["Q"], d["J"], d["Qdot"]])

    def build_matched_twin(h_center, target_4, a6_val, n_iter=7):
        a4 = np.zeros(4)
        eps_fd = 1e-5
        for _ in range(n_iter):
            h_curr = h_center + sum(a4[j] * basis_h[j] for j in range(4)) + a6_val * basis_h[5]
            f0 = eval_state_vec(h_curr)
            res = f0 - target_4
            if np.max(np.abs(res)) < 1e-11:
                break
            Jac = np.zeros((4, 4))
            for j in range(4):
                fj = eval_state_vec(h_curr + eps_fd * basis_h[j])
                Jac[:, j] = (fj - f0) / eps_fd
            step_a = np.linalg.lstsq(Jac, -res, rcond=1e-4)[0]
            a4 = a4 + step_a
        return h_center + sum(a4[j] * basis_h[j] for j in range(4)) + a6_val * basis_h[5]

    def evolve_to_crossing(h_init, dt_ev=0.005, max_steps=100):
        h_c = h_init.copy()
        series = []
        t_turn = None
        for s in range(max_steps + 1):
            t_rel = round(s * dt_ev, 5)
            d = sol24_05.diagnostics(h_c, t_rel, compute_second_order=False, compute_shells=False)
            series.append({"t_rel": t_rel, "H": d["H"], "eta_plus": d["eta_plus"], "Pplus": d["Pplus"], "Q": d["Q"]})
            if len(series) >= 2 and t_turn is None and series[-2]["H"] >= 0.0 and series[-1]["H"] < 0.0:
                h0_v, h1_v = series[-2]["H"], series[-1]["H"]
                t_turn = series[-2]["t_rel"] + (h0_v / (h0_v - h1_v)) * dt_ev
            if s < max_steps:
                h_c = sol24_05.step(h_c, dt_ev)
        return float(t_turn) if t_turn is not None else -1.0, series

    twin_suite = []
    featured_twins = {}
    for t_base in [0.35, 0.40, 0.45]:
        h_b = sol24_05.initial_exp013_3d_pert(eps=0.1, sign=+1)
        for _ in range(int(round(t_base / 0.005))):
            h_b = sol24_05.step(h_b, 0.005)
        t_vec4 = eval_state_vec(h_b)

        for amp6 in [0.03, 0.05, 0.07]:
            hA = build_matched_twin(h_b, t_vec4, -amp6, n_iter=7)
            hB = build_matched_twin(h_b, t_vec4, +amp6, n_iter=7)
            dA = sol24_05.diagnostics(hA, 0.0, compute_second_order=True, full_poly_split=True, compute_shells=True)
            dB = sol24_05.diagnostics(hB, 0.0, compute_second_order=True, full_poly_split=True, compute_shells=True)
            tauA, serA = evolve_to_crossing(hA)
            tauB, serB = evolve_to_crossing(hB)
            row = {
                "t_base": t_base, "amp6": amp6,
                "max_state_diff": float(max(abs(dA[k] - dB[k]) for k in ["Pplus", "Q", "eta_plus", "J", "Qdot", "H"])),
                "Pplus_0": dA["Pplus"], "Q_0": dA["Q"], "eta_plus_0": dA["eta_plus"],
                "J_0": dA["J"], "Qdot_0": dA["Qdot"], "H_0": dA["H"],
                "dH_A": dA["dH"], "dH_B": dB["dH"],
                "Delta_T5_A": dA["Delta_T5"], "Delta_T5_B": dB["Delta_T5"],
                "S_pou_A": dA["S_pou"], "S_pou_B": dB["S_pou"],
                "pred_tau_A": float(-dA["H"] / dA["dH"]), "pred_tau_B": float(-dB["H"] / dB["dH"]),
                "true_tau_A": tauA, "true_tau_B": tauB,
            }
            twin_suite.append(row)
            if t_base == 0.40 and amp6 == 0.05:
                featured_twins = {
                    "diagA_0": dA, "diagB_0": dB,
                    "tau_turn_A": tauA, "tau_turn_B": tauB,
                    "series_A": serA, "series_B": serB,
                }

    print(f"3. Twin suite completed in {time.time()-t0:.1f}s:", flush=True)
    for r in twin_suite:
        print(
            f"   t0={r['t_base']}, amp={r['amp6']}: max_diff={r['max_state_diff']:.2e} | "
            f"dH_A={r['dH_A']:.4f}, dH_B={r['dH_B']:.4f} | "
            f"pred_dtau={r['pred_tau_B']-r['pred_tau_A']:+.4f} vs true_dtau={r['true_tau_B']-r['true_tau_A']:+.4f}",
            flush=True,
        )

    # 4. Exact 60-Trajectory 5-Family Ensemble
    families = ["exp013_3d", "tg3d_pert", "antiparallel", "random_solenoidal", "phase_rand"]
    solvers_24 = {nu: Checkpoint4Solver(24, nu) for nu in (0.05, 0.02, 0.01)}
    trajectories = []
    manifest_rows = []
    traj_id = 0
    for fam in families:
        for variant in range(2):
            for nu, amp in [(0.05, 1.3), (0.02, 1.15), (0.01, 1.0)]:
                gid = f"{fam}_v{variant}_nu{nu}"
                for sign in (+1, -1):
                    h0 = solvers_24[nu].initial_ensemble(fam, variant=variant, amp=amp, sign=sign, rng_seed=42)
                    rows = run_traj(solvers_24[nu], h0, T=1.0, dt=0.02, sample_dt=0.05, full_poly_split=False)
                    crossings = []
                    for i in range(len(rows) - 1):
                        if rows[i]["H"] > 0.0 and rows[i + 1]["H"] <= 0.0:
                            ha, hb = rows[i]["H"], rows[i + 1]["H"]
                            tc = rows[i]["t"] + (ha / (ha - hb)) * 0.05
                            crossings.append(float(tc))

                    trajectories.append({
                        "traj_id": traj_id, "group_id": gid, "family": fam,
                        "variant": variant, "nu": nu, "amp": amp, "sign": sign,
                        "crossings": crossings, "rows": rows,
                    })
                    manifest_rows.append({
                        "traj_id": traj_id, "group_id": gid, "family": fam,
                        "variant": variant, "nu": nu, "amp": amp, "sign": sign, "rng_seed": 42,
                        "Pplus_0": round(rows[0]["Pplus"], 6), "Q_0": round(rows[0]["Q"], 6),
                        "eta_plus_0": round(rows[0]["eta_plus"], 6), "H_0": round(rows[0]["H"], 6),
                        "dH_0": round(rows[0]["dH"], 6), "k_char_0": round(rows[0]["k_char"], 4),
                        "first_crossing_t": round(crossings[0], 5) if crossings else -1.0,
                    })
                    traj_id += 1

    print(f"4. Completed 60-trajectory 5-family ensemble in {time.time()-t0:.1f}s", flush=True)

    with open(out_dir / "ns_checkpoint3_4_reproducibility_manifest.csv", "w", newline="") as f_csv:
        writer = csv.DictWriter(f_csv, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)

    horizon_results = {}
    for h_steps, h_label, h_val in [(3, "h_0.15", 0.15), (5, "h_0.25", 0.25)]:
        dt_lag = 0.05
        B0_eta, J_e, R_e, Qdot_e, HardDyad_e, SmoothPoU_e, Acc2nd_e = [], [], [], [], [], [], []
        y_eta, y_H = [], []
        groups, fams = [], []

        pre_mask = []
        y_event_pre = []
        turn_mask = []
        y_tau_turn = []

        for tr in trajectories:
            rows = tr["rows"]
            for i in range(1, len(rows) - h_steps):
                rp, rc, rf = rows[i - 1], rows[i], rows[i + h_steps]
                t_now, t_fut = rc["t"], rf["t"]
                Q = rc["Q"] + 1e-12
                Pplus = rc["Pplus"] + 1e-12
                W = rc["W"] + 1e-12

                dPplus_hist = (rc["Pplus"] - rp["Pplus"]) / dt_lag
                dP_hist = (rc["P"] - rp["P"]) / dt_lag
                deta_hist = (rc["eta_plus"] - rp["eta_plus"]) / dt_lag
                dQ_hist = (rc["Q"] - rp["Q"]) / dt_lag
                dW_hist = (rc["W"] - rp["W"]) / dt_lag
                dD_hist = (rc["D"] - rp["D"]) / dt_lag

                B0_eta.append([
                    rc["eta_plus"], rc["P"] / Q, rc["Pminus"] / Q,
                    rc["C_conc"], rc["q95"] / math.sqrt(W), rc["D"] / Q,
                    rc["P"] / (rc["D"] + 1e-12),
                    deta_hist, dPplus_hist / Q, dQ_hist / Q, dW_hist / W, dD_hist / Q,
                ])
                J_e.append([rc["J"] / Q, rc["J"] / Pplus])
                R_e.append([rc["R"] / Q, rc["R"] / Pplus])
                Qdot_e.append([rc["Qdot"] / Q, rc["H"]])

                HardDyad_e.append([
                    rc["S_dyad"], (rc["T4_B2"] + rc["T4_B3"]) / Pplus, rc["Delta_triad"], rc["Delta_visc"]
                ])
                SmoothPoU_e.append([
                    rc["S_pou"], (rc["T4_P2"] + rc["T4_P3"]) / Q, rc["Delta_triad"]
                ])
                tau_kin = np.clip(-rc["H"] / min(rc["dH"], -0.15), 0.0, 1.0) if rc["H"] > 0 else 0.0
                Acc2nd_e.append([
                    rc["dH"], rc["ddP_plus"] / Pplus, rc["ddQ"] / Q, tau_kin, rc["H"] + 0.5 * h_val * rc["dH"]
                ])

                y_eta.append(rf["eta_plus"])
                y_H.append(rf["H"])
                groups.append(tr["group_id"])
                fams.append(tr["family"])

                is_pre = rc["H"] > 0.0
                pre_mask.append(is_pre)
                if is_pre:
                    crossed = 1.0 if (rf["H"] <= 0.0 or any((c > t_now and c <= t_fut + 1e-9) for c in tr["crossings"])) else 0.0
                    y_event_pre.append(crossed)
                else:
                    y_event_pre.append(0.0)

                fut_c = [c for c in tr["crossings"] if c > t_now]
                if is_pre and len(fut_c) > 0:
                    turn_mask.append(True)
                    y_tau_turn.append(fut_c[0] - t_now)
                else:
                    turn_mask.append(False)
                    y_tau_turn.append(0.0)

        B0_eta = np.array(B0_eta)
        J_e, R_e, Qdot_e = map(np.array, (J_e, R_e, Qdot_e))
        HardDyad_e, SmoothPoU_e, Acc2nd_e = map(np.array, (HardDyad_e, SmoothPoU_e, Acc2nd_e))
        y_eta, y_H = np.array(y_eta), np.array(y_H)
        groups, fams = np.array(groups), np.array(fams)
        pre_mask = np.array(pre_mask, dtype=bool)
        turn_mask = np.array(turn_mask, dtype=bool)
        y_event_pre = np.array(y_event_pre)[pre_mask]
        y_H_pre = y_H[pre_mask]
        y_tau_turn = np.array(y_tau_turn)[turn_mask]

        models_e = {
            "M0_B0": B0_eta,
            "M1_B0+J": np.hstack([B0_eta, J_e]),
            "M3_B0+J+Qdot": np.hstack([B0_eta, J_e, Qdot_e]),
            "M4_B0+J+Qdot+R": np.hstack([B0_eta, J_e, Qdot_e, R_e]),
            "M5_HardDyad": np.hstack([B0_eta, J_e, Qdot_e, R_e, HardDyad_e]),
            "M5_SmoothPoU": np.hstack([B0_eta, J_e, Qdot_e, R_e, SmoothPoU_e]),
            "M5_Acc2nd": np.hstack([B0_eta, J_e, Qdot_e, R_e, Acc2nd_e]),
            "M6_SmoothPoU+Acc2nd": np.hstack([B0_eta, J_e, Qdot_e, R_e, SmoothPoU_e, Acc2nd_e]),
        }

        M4_mat = models_e["M4_B0+J+Qdot+R"]
        schur_tests = {
            "HardDyad_on_eta_given_M4": exact_population_schur_r2(M4_mat, HardDyad_e[:, 0], y_eta)["partial_r2_pct"],
            "SmoothPoU_on_eta_given_M4": exact_population_schur_r2(M4_mat, SmoothPoU_e[:, 0], y_eta)["partial_r2_pct"],
            "dH_on_H_given_M4": exact_population_schur_r2(M4_mat, Acc2nd_e[:, 0], y_H)["partial_r2_pct"],
            "dH_on_tau_given_M4_pre": exact_population_schur_r2(M4_mat[turn_mask], Acc2nd_e[turn_mask, 0], y_tau_turn)["partial_r2_pct"],
        }

        u_groups, u_fams = sorted(set(groups)), sorted(set(fams))
        u_groups_pre, u_fams_pre = sorted(set(groups[pre_mask])), sorted(set(fams[pre_mask]))
        u_groups_turn, u_fams_turn = sorted(set(groups[turn_mask])), sorted(set(fams[turn_mask]))

        def cv_predict_subset(X_mat, y_vec, grp_arr, u_grps, alpha=2e-3):
            preds = np.zeros_like(y_vec)
            for g in u_grps:
                te = grp_arr == g; tr = ~te
                preds[te] = fit_predict_ridge(X_mat[tr], y_vec[tr], X_mat[te], alpha=alpha)
            return preds

        rmse = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))

        base_eta_logo = rmse(y_eta, cv_predict_subset(models_e["M0_B0"], y_eta, groups, u_groups, alpha=2e-3))
        base_eta_lofo = rmse(y_eta, cv_predict_subset(models_e["M0_B0"], y_eta, fams, u_fams, alpha=2e-3))

        base_Hpre_logo = rmse(y_H_pre, cv_predict_subset(models_e["M0_B0"][pre_mask], y_H_pre, groups[pre_mask], u_groups_pre, alpha=5e-3))
        base_Hpre_lofo = rmse(y_H_pre, cv_predict_subset(models_e["M0_B0"][pre_mask], y_H_pre, fams[pre_mask], u_fams_pre, alpha=2e-2))

        base_tau_logo = rmse(y_tau_turn, cv_predict_subset(models_e["M0_B0"][turn_mask], y_tau_turn, groups[turn_mask], u_groups_turn, alpha=1e-2))
        base_tau_lofo = rmse(y_tau_turn, cv_predict_subset(models_e["M0_B0"][turn_mask], y_tau_turn, fams[turn_mask], u_fams_turn, alpha=5e-2))

        model_metrics = {}
        for m_name, X_full in models_e.items():
            p_eta_logo = cv_predict_subset(X_full, y_eta, groups, u_groups, alpha=2e-3)
            p_eta_lofo = cv_predict_subset(X_full, y_eta, fams, u_fams, alpha=2e-3)
            eta_logo_gain = 100.0 * (base_eta_logo - rmse(y_eta, p_eta_logo)) / base_eta_logo
            eta_lofo_gain = 100.0 * (base_eta_lofo - rmse(y_eta, p_eta_lofo)) / base_eta_lofo

            X_pre = X_full[pre_mask]
            p_Hpre_logo = cv_predict_subset(X_pre, y_H_pre, groups[pre_mask], u_groups_pre, alpha=5e-3)
            p_Hpre_lofo = cv_predict_subset(X_pre, y_H_pre, fams[pre_mask], u_fams_pre, alpha=2e-2)
            Hpre_logo_gain = 100.0 * (base_Hpre_logo - rmse(y_H_pre, p_Hpre_logo)) / base_Hpre_logo
            Hpre_lofo_gain = 100.0 * (base_Hpre_lofo - rmse(y_H_pre, p_Hpre_lofo)) / base_Hpre_lofo

            ev_acc_logo = 100.0 * float(np.mean((p_Hpre_logo <= 0.0).astype(float) == y_event_pre))
            ev_acc_lofo = 100.0 * float(np.mean((p_Hpre_lofo <= 0.0).astype(float) == y_event_pre))
            ev_auroc_logo = auroc_score(y_event_pre, -p_Hpre_logo)
            ev_auroc_lofo = auroc_score(y_event_pre, -p_Hpre_lofo)

            X_turn = X_full[turn_mask]
            p_tau_logo = cv_predict_subset(X_turn, y_tau_turn, groups[turn_mask], u_groups_turn, alpha=1e-2)
            p_tau_lofo = cv_predict_subset(X_turn, y_tau_turn, fams[turn_mask], u_fams_turn, alpha=5e-2)
            tau_logo_gain = 100.0 * (base_tau_logo - rmse(y_tau_turn, p_tau_logo)) / base_tau_logo
            tau_lofo_gain = 100.0 * (base_tau_lofo - rmse(y_tau_turn, p_tau_lofo)) / base_tau_lofo

            model_metrics[m_name] = {
                "eta_logo_gain_pct": round(eta_logo_gain, 2),
                "eta_lofo_gain_pct": round(eta_lofo_gain, 2),
                "Hpre_logo_gain_pct": round(Hpre_logo_gain, 2),
                "Hpre_lofo_gain_pct": round(Hpre_lofo_gain, 2),
                "event_acc_logo_pct": round(ev_acc_logo, 2),
                "event_acc_lofo_pct": round(ev_acc_lofo, 2),
                "event_auroc_logo": round(ev_auroc_logo, 4),
                "event_auroc_lofo": round(ev_auroc_lofo, 4),
                "tau_logo_gain_pct": round(tau_logo_gain, 2),
                "tau_lofo_gain_pct": round(tau_lofo_gain, 2),
                "tau_rmse_logo": round(rmse(y_tau_turn, p_tau_logo), 4),
                "tau_rmse_lofo": round(rmse(y_tau_turn, p_tau_lofo), 4),
            }

        horizon_results[h_label] = {
            "n_all": int(len(y_eta)),
            "n_pre_H_pos": int(np.sum(pre_mask)),
            "n_events_crossed": int(np.sum(y_event_pre)),
            "n_turn_samples": int(np.sum(turn_mask)),
            "schur_tests": schur_tests,
            "models": model_metrics,
        }

        print(f"\n================ HORIZON {h_label} (All={len(y_eta)}, Pre-Crossing H>0={int(np.sum(pre_mask))}, Crossings={int(np.sum(y_event_pre))}, TurnSamples={int(np.sum(turn_mask))}) ================", flush=True)
        print("Schur Partial R2 (%):", json.dumps(schur_tests, indent=2))
        for m_name, mm in model_metrics.items():
            print(
                f"  {m_name:22s} | eta+ LOGO/LOFO: {mm['eta_logo_gain_pct']:+6.2f}% / {mm['eta_lofo_gain_pct']:+6.2f}% "
                f"| Pre-Cross Event Acc (LOGO/LOFO): {mm['event_acc_logo_pct']:5.2f}% / {mm['event_acc_lofo_pct']:5.2f}% "
                f"(AUROC: {mm['event_auroc_logo']:.4f}/{mm['event_auroc_lofo']:.4f}) "
                f"| tau_turn Gain (LOGO/LOFO): {mm['tau_logo_gain_pct']:+6.2f}% / {mm['tau_lofo_gain_pct']:+6.2f}%",
                flush=True,
            )

    full_out = {
        "elapsed_seconds": time.time() - t0,
        "sylvester_verification": sylv_res,
        "crossing_audit_3d_witness": crossing_audit,
        "witness_3dpert_curve": rows_3dpert,
        "featured_twins": featured_twins,
        "twin_suite": twin_suite,
        "horizon_results": horizon_results,
    }
    (out_dir / "ns_checkpoint4_results.json").write_text(json.dumps(full_out, indent=2))
    print(f"\nSaved ns_checkpoint4_results.json successfully in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()

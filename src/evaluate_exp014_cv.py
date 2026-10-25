import csv
import json
import math

def solve_linear_system(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        piv = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[piv] = M[piv], M[k]
        div = M[k][k]
        if abs(div) < 1e-14:
            div = 1e-14
        for j in range(k, n + 1):
            M[k][j] /= div
        for i in range(n):
            if i != k:
                f = M[i][k]
                if f != 0.0:
                    for j in range(k, n + 1):
                        M[i][j] -= f * M[k][j]
    return [M[i][n] for i in range(n)]

def fit_predict_ridge(X_train, y_train, X_test, alpha=2e-3):
    n_tr = len(X_train)
    p = len(X_train[0])
    mean_x = [sum(X_train[i][j] for i in range(n_tr)) / n_tr for j in range(p)]
    std_x = []
    for j in range(p):
        var = sum((X_train[i][j] - mean_x[j]) ** 2 for i in range(n_tr)) / n_tr
        s = math.sqrt(var)
        std_x.append(s if s >= 1e-12 else 1.0)
    Xz_tr = [[(X_train[i][j] - mean_x[j]) / std_x[j] for j in range(p)] for i in range(n_tr)]
    Xz_te = [[(X_test[i][j] - mean_x[j]) / std_x[j] for j in range(p)] for i in range(len(X_test))]
    mean_y = sum(y_train) / n_tr
    yc = [y - mean_y for y in y_train]

    XtX = [[0.0] * p for _ in range(p)]
    Xty = [0.0] * p
    for i in range(n_tr):
        row = Xz_tr[i]
        yi = yc[i]
        for j in range(p):
            rj = row[j]
            Xty[j] += rj * yi
            for k in range(j, p):
                XtX[j][k] += rj * row[k]
    for j in range(p):
        for k in range(j):
            XtX[j][k] = XtX[k][j]
        XtX[j][j] += alpha * n_tr

    w = solve_linear_system(XtX, Xty)
    return [mean_y + sum(Xz_te[i][j] * w[j] for j in range(p)) for i in range(len(Xz_te))]

def exact_schur_partial_r2(X_base, r_vec, y_vec):
    n = len(y_vec)
    y_hat = fit_predict_ridge(X_base, y_vec, X_base, alpha=1e-8)
    r_hat = fit_predict_ridge(X_base, r_vec, X_base, alpha=1e-8)
    y_perp = [y_vec[i] - y_hat[i] for i in range(n)]
    r_perp = [r_vec[i] - r_hat[i] for i in range(n)]
    var_y = sum(v * v for v in y_perp) / n
    var_r = sum(v * v for v in r_perp) / n
    cov_yr = sum(y_perp[i] * r_perp[i] for i in range(n)) / n
    delta_mse = (cov_yr ** 2) / max(var_r, 1e-30)
    return 100.0 * delta_mse / max(var_y, 1e-30)

def auroc_score(y_true, y_score):
    pos = [s for y, s in zip(y_true, y_score) if y == 1]
    neg = [s for y, s in zip(y_true, y_score) if y == 0]
    if not pos or not neg:
        return 0.5
    wins = 0.0
    for sp in pos:
        for sn in neg:
            if sp > sn:
                wins += 1.0
            elif sp == sn:
                wins += 0.5
    return wins / (len(pos) * len(neg))

def rmse(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / len(a))

def main():
    with open("ns_research/ns_exp014_twins.csv", "r") as f:
        twins = list(csv.DictReader(f))
    print("=== MATCHED-STATE COUNTERFACTUAL TWINS ===")
    for r in twins:
        dt_pred = float(r["pred_tau_B"]) - float(r["pred_tau_A"])
        dt_true = float(r["true_tau_B"]) - float(r["true_tau_A"])
        print(f"t0={r['t_base']}, amp={r['amp_high']}: max_diff={r['max_state_diff']} | dH_A={r['dH_A']}, dH_B={r['dH_B']} | pred_dtau={dt_pred:+.4f} vs true_dtau={dt_true:+.4f}")

    with open("ns_research/ns_exp014_snapshots.csv", "r") as f:
        raw_snaps = list(csv.DictReader(f))

    by_traj = {}
    for row in raw_snaps:
        tid = int(row["traj_id"])
        by_traj.setdefault(tid, []).append(row)
    for tid in by_traj:
        by_traj[tid].sort(key=lambda r: int(r["snap_idx"]))

    results = {}
    for h_steps, h_label, h_val in [(3, "h_0.15", 0.15), (5, "h_0.25", 0.25)]:
        dt_lag = 0.05
        m_keys = ["M0_B0", "M1_B0+J", "M3_B0+J+Qdot", "M4_B0+J+Qdot+R", "M5_HardDyad", "M5_SmoothPoU", "M5_Acc2nd", "M6_SmoothPoU+Acc2nd"]
        models_all = {k: [] for k in m_keys}
        models_tau = {k: [] for k in m_keys}
        y_eta, y_P, y_H = [], [], []
        groups_all, fams_all = [], []
        s_dyad_col, s_pou_col, dH_col, tau_kin_col = [], [], [], []

        pre_idx = []
        y_ev_pre, y_H_pre = [], []
        groups_pre, fams_pre = [], []

        turn_idx = []
        y_tau_turn = []
        groups_turn, fams_turn = [], []

        sample_counter = 0
        for tid, rows in sorted(by_traj.items()):
            crossings = []
            for j in range(len(rows) - 1):
                h0 = float(rows[j]["H"])
                h1 = float(rows[j+1]["H"])
                if h0 > 0.0 and h1 <= 0.0:
                    tc = float(rows[j]["t"]) + (h0 / (h0 - h1)) * 0.05
                    crossings.append(tc)

            for i in range(1, len(rows) - h_steps):
                rp, rc, rf = rows[i - 1], rows[i], rows[i + h_steps]
                t_now, t_fut = float(rc["t"]), float(rf["t"])
                Q = float(rc["Q"]) + 1e-12
                Pplus = float(rc["Pplus"]) + 1e-12
                W = float(rc["W"]) + 1e-12
                D = float(rc["D"])
                P = float(rc["P"])
                Pminus = float(rc["Pminus"])
                eta = float(rc["eta_plus"])
                C_conc = float(rc["C_conc"])
                S_max = float(rc["S_max"])

                dPplus_hist = (float(rc["Pplus"]) - float(rp["Pplus"])) / dt_lag
                dP_hist = (float(rc["P"]) - float(rp["P"])) / dt_lag
                deta_hist = (float(rc["eta_plus"]) - float(rp["eta_plus"])) / dt_lag
                dQ_hist = (float(rc["Q"]) - float(rp["Q"])) / dt_lag
                dW_hist = (float(rc["W"]) - float(rp["W"])) / dt_lag
                dD_hist = (float(rc["D"]) - float(rp["D"])) / dt_lag

                B0 = [
                    eta, P / Q, Pminus / Q, C_conc, S_max / math.sqrt(W), D / Q,
                    P / (D + 1e-12),
                    deta_hist, dPplus_hist / Q, dQ_hist / Q, dW_hist / W, dD_hist / Q,
                ]
                J_v = float(rc["J"])
                R_v = float(rc["R"])
                Qdot_v = float(rc["Qdot"])
                H_v = float(rc["H"])
                dH_v = float(rc["dH"])
                ddP_v = float(rc["ddP_plus"])
                ddQ_v = float(rc["ddQ"])
                S_dyad_v = float(rc["S_dyad"])
                S_pou_v = float(rc["S_pou"])
                T4_B1_v = float(rc["T4_B1"])
                T4_B2_v = float(rc["T4_B2"])
                T4_B3_v = float(rc["T4_B3"])
                T4_P1_v = float(rc["T4_P1"])
                T4_P2_v = float(rc["T4_P2"])
                T4_P3_v = float(rc["T4_P3"])
                D_triad_v = float(rc["Delta_triad"])

                J_e = [J_v / Q, J_v / Pplus]
                Qdot_e = [Qdot_v / Q, H_v]
                R_e = [R_v / Q, R_v / Pplus]

                # Hard unnormalized shells vs Dimensionless scale-normalized PoU shells
                HardDyad = [T4_B1_v, T4_B2_v + T4_B3_v, S_dyad_v, (T4_B2_v + T4_B3_v) / Pplus]
                T4_tot = max(abs(T4_P1_v + T4_P2_v + T4_P3_v), 1e-3)
                SmoothPoU = [
                    (T4_P2_v + T4_P3_v) / Pplus,
                    T4_P3_v / Pplus,
                    D_triad_v,
                ]

                tau_kin = min(max(-H_v / min(dH_v, -0.25), 0.0), 0.80) if H_v > 0.0 else 0.0
                Acc2nd = [dH_v, ddP_v / Pplus, ddQ_v / Q, tau_kin, H_v + 0.5 * h_val * dH_v]
                # Bounded kinematic horizon estimator tau_kin = -H / dH (ChatGPT Exp 013 Eq. 3 with exact snapshot dH)
                Acc2nd_tau = [tau_kin]

                models_all["M0_B0"].append(B0)
                models_all["M1_B0+J"].append(B0 + J_e)
                models_all["M3_B0+J+Qdot"].append(B0 + J_e + Qdot_e)
                models_all["M4_B0+J+Qdot+R"].append(B0 + J_e + Qdot_e + R_e)
                models_all["M5_HardDyad"].append(B0 + J_e + Qdot_e + R_e + HardDyad)
                models_all["M5_SmoothPoU"].append(B0 + J_e + Qdot_e + R_e + SmoothPoU)
                models_all["M5_Acc2nd"].append(B0 + J_e + Qdot_e + R_e + Acc2nd)
                models_all["M6_SmoothPoU+Acc2nd"].append(B0 + J_e + Qdot_e + R_e + SmoothPoU + Acc2nd)

                models_tau["M0_B0"].append(B0)
                models_tau["M1_B0+J"].append(B0 + J_e)
                models_tau["M3_B0+J+Qdot"].append(B0 + J_e + Qdot_e)
                models_tau["M4_B0+J+Qdot+R"].append(B0 + J_e + Qdot_e + R_e)
                models_tau["M5_HardDyad"].append(B0 + J_e + Qdot_e + R_e + HardDyad)
                models_tau["M5_SmoothPoU"].append(B0 + J_e + Qdot_e + R_e + SmoothPoU)
                models_tau["M5_Acc2nd"].append(B0 + J_e + Qdot_e + R_e + Acc2nd_tau)
                models_tau["M6_SmoothPoU+Acc2nd"].append(B0 + J_e + Qdot_e + R_e + SmoothPoU + Acc2nd_tau)

                y_eta.append(float(rf["eta_plus"]))
                y_P.append(float(rf["Pplus"]))
                y_H.append(float(rf["H"]))
                groups_all.append(rc["group_id"])
                fams_all.append(rc["family"])
                s_dyad_col.append(S_dyad_v)
                s_pou_col.append(S_pou_v)
                dH_col.append(dH_v)
                tau_kin_col.append(tau_kin)

                if H_v > 0.0:
                    pre_idx.append(sample_counter)
                    crossed = 1 if (float(rf["H"]) <= 0.0 or any(t_now < c <= t_fut + 1e-9 for c in crossings)) else 0
                    y_ev_pre.append(crossed)
                    y_H_pre.append(float(rf["H"]))
                    groups_pre.append(rc["group_id"])
                    fams_pre.append(rc["family"])

                    fut_c = [c for c in crossings if c > t_now]
                    if fut_c:
                        turn_idx.append(sample_counter)
                        y_tau_turn.append(fut_c[0] - t_now)
                        groups_turn.append(rc["group_id"])
                        fams_turn.append(rc["family"])

                sample_counter += 1

        def cv_run(X_list, y_list, grp_list, alpha=2e-3):
            u_grps = sorted(set(grp_list))
            preds = [0.0] * len(y_list)
            for g in u_grps:
                tr_i = [idx for idx, gg in enumerate(grp_list) if gg != g]
                te_i = [idx for idx, gg in enumerate(grp_list) if gg == g]
                X_tr = [X_list[idx] for idx in tr_i]
                y_tr = [y_list[idx] for idx in tr_i]
                X_te = [X_list[idx] for idx in te_i]
                p_te = fit_predict_ridge(X_tr, y_tr, X_te, alpha=alpha)
                for idx, val in zip(te_i, p_te):
                    preds[idx] = val
            return preds

        base_eta_logo = rmse(y_eta, cv_run(models_all["M0_B0"], y_eta, groups_all, alpha=2e-3))
        base_eta_lofo = rmse(y_eta, cv_run(models_all["M0_B0"], y_eta, fams_all, alpha=2e-3))

        X_B0_pre = [models_all["M0_B0"][i] for i in pre_idx]
        base_Hpre_logo = rmse(y_H_pre, cv_run(X_B0_pre, y_H_pre, groups_pre, alpha=5e-3))
        base_Hpre_lofo = rmse(y_H_pre, cv_run(X_B0_pre, y_H_pre, fams_pre, alpha=2e-2))

        X_B0_turn = [models_tau["M0_B0"][i] for i in turn_idx]
        base_tau_logo = rmse(y_tau_turn, cv_run(X_B0_turn, y_tau_turn, groups_turn, alpha=1e-2))
        base_tau_lofo = rmse(y_tau_turn, cv_run(X_B0_turn, y_tau_turn, fams_turn, alpha=5e-2))

        M4_full = models_all["M4_B0+J+Qdot+R"]
        M4_turn = [M4_full[i] for i in turn_idx]
        tk_turn = [tau_kin_col[i] for i in turn_idx]

        schur = {
            "SmoothPoU_on_eta_given_M4_pct": round(exact_schur_partial_r2(M4_full, s_pou_col, y_eta), 4),
            "dH_on_H_given_M4_pct": round(exact_schur_partial_r2(M4_full, dH_col, y_H), 4),
            "tau_kin_on_tau_given_M4_pre_pct": round(exact_schur_partial_r2(M4_turn, tk_turn, y_tau_turn), 4),
        }

        m_summary = {}
        print(f"\n=== HORIZON {h_label} (All={len(y_eta)}, Pre-Cross H>0={len(pre_idx)}, CrossEvents={sum(y_ev_pre)}, TurnSamples={len(turn_idx)}) ===")
        print("Schur Partial R2 (%):", schur)
        for m_name in m_keys:
            X_full = models_all[m_name]
            X_t_full = models_tau[m_name]
            p_eta_logo = cv_run(X_full, y_eta, groups_all, alpha=2e-3)
            p_eta_lofo = cv_run(X_full, y_eta, fams_all, alpha=2e-3)
            g_eta_logo = 100.0 * (base_eta_logo - rmse(y_eta, p_eta_logo)) / base_eta_logo
            g_eta_lofo = 100.0 * (base_eta_lofo - rmse(y_eta, p_eta_lofo)) / base_eta_lofo

            X_pre = [X_full[i] for i in pre_idx]
            p_Hpre_logo = cv_run(X_pre, y_H_pre, groups_pre, alpha=5e-3)
            p_Hpre_lofo = cv_run(X_pre, y_H_pre, fams_pre, alpha=2e-2)
            g_Hpre_logo = 100.0 * (base_Hpre_logo - rmse(y_H_pre, p_Hpre_logo)) / base_Hpre_logo
            g_Hpre_lofo = 100.0 * (base_Hpre_lofo - rmse(y_H_pre, p_Hpre_lofo)) / base_Hpre_lofo

            ev_acc_logo = 100.0 * sum(1 for p, y in zip(p_Hpre_logo, y_ev_pre) if (1 if p <= 0.0 else 0) == y) / len(y_ev_pre)
            ev_acc_lofo = 100.0 * sum(1 for p, y in zip(p_Hpre_lofo, y_ev_pre) if (1 if p <= 0.0 else 0) == y) / len(y_ev_pre)
            auroc_logo = auroc_score(y_ev_pre, [-p for p in p_Hpre_logo])
            auroc_lofo = auroc_score(y_ev_pre, [-p for p in p_Hpre_lofo])

            X_turn = [X_t_full[i] for i in turn_idx]
            p_tau_logo = cv_run(X_turn, y_tau_turn, groups_turn, alpha=1e-2)
            p_tau_lofo = cv_run(X_turn, y_tau_turn, fams_turn, alpha=5e-2)
            g_tau_logo = 100.0 * (base_tau_logo - rmse(y_tau_turn, p_tau_logo)) / base_tau_logo
            g_tau_lofo = 100.0 * (base_tau_lofo - rmse(y_tau_turn, p_tau_lofo)) / base_tau_lofo

            m_summary[m_name] = {
                "eta_logo_gain": round(g_eta_logo, 2), "eta_lofo_gain": round(g_eta_lofo, 2),
                "Hpre_logo_gain": round(g_Hpre_logo, 2), "Hpre_lofo_gain": round(g_Hpre_lofo, 2),
                "ev_acc_logo": round(ev_acc_logo, 2), "ev_acc_lofo": round(ev_acc_lofo, 2),
                "auroc_logo": round(auroc_logo, 4), "auroc_lofo": round(auroc_lofo, 4),
                "tau_logo_gain": round(g_tau_logo, 2), "tau_lofo_gain": round(g_tau_lofo, 2),
                "tau_rmse_logo": round(rmse(y_tau_turn, p_tau_logo), 4),
                "tau_rmse_lofo": round(rmse(y_tau_turn, p_tau_lofo), 4),
            }
            print(
                f"{m_name:22s} | eta+ LOGO/LOFO: {g_eta_logo:+6.2f}%/{g_eta_lofo:+6.2f}% "
                f"| Event E_h Acc (LOGO/LOFO): {ev_acc_logo:5.2f}%/{ev_acc_lofo:5.2f}% (AUROC: {auroc_logo:.4f}/{auroc_lofo:.4f}) "
                f"| tau_turn Gain (LOGO/LOFO): {g_tau_logo:+6.2f}%/{g_tau_lofo:+6.2f}%"
            )

        results[h_label] = {
            "n_all": len(y_eta), "n_pre": len(pre_idx), "n_events": sum(y_ev_pre), "n_turn": len(turn_idx),
            "schur": schur, "models": m_summary
        }

    with open("ns_research/ns_checkpoint4_cv_summary.json", "w") as f:
        json.dump({"twins": twins, "horizons": results}, f, indent=2)

if __name__ == "__main__":
    main()

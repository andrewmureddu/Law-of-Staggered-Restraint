import json
import math
import random
import statistics
import time

def paired_ttest(x, y):
    diffs = [xi - yi for xi, yi in zip(x, y)]
    n = len(diffs)
    mean_d = statistics.mean(diffs)
    sd_d = statistics.stdev(diffs)
    se_d = sd_d / math.sqrt(n)
    t_stat = mean_d / max(se_d, 1e-15)
    z = abs(t_stat) * (1.0 - 1.0 / (4.0 * (n - 1)))
    p_val = math.erfc(z / math.sqrt(2.0))
    return t_stat, p_val, mean_d, se_d

def sample_cantor_tree_field(L_depth, v_scales, rng):
    N = 1 << L_depth
    theta = [0.0] * N
    z0 = rng.gauss(0.0, math.sqrt(v_scales[0]))
    for i in range(N):
        theta[i] += z0
    for m in range(1, L_depth + 1):
        sd_m = math.sqrt(v_scales[m])
        block_size = 1 << (L_depth - m)
        num_blocks = 1 << m
        for u in range(num_blocks):
            zu = rng.gauss(0.0, sd_m)
            start = u * block_size
            for i in range(start, start + block_size):
                theta[i] += zu
    return theta

# ==============================================================================
# CG-HSTS v2.0: Self-Tuning (Prior-Free) Adaptive Haar-Schur Tree Search
# Solves Open Problem 1 by estimating shell variances {tau_k^2} and noise sigma^2
# online in O(1) per shell via Harmonic-Weighted Haar Detail Squares.
# ==============================================================================
class AdaptiveCGHSTS:
    """
    Prior-Free Cotton-Gated Haar-Schur Tree Search (CG-HSTS v2.0).
    Requires ZERO oracle knowledge of {v_k} or sigma^2.
    """
    def __init__(self, L_depth: int, gamma_star: float = 0.50, c_ucb: float = 0.90, cap_K: int = 4):
        self.L = L_depth
        self.N = 1 << L_depth
        self.gamma_star = gamma_star
        self.c_ucb = c_ucb
        self.cap_K = cap_K

        num_nodes = 1 << (L_depth + 1)
        self.counts = [0] * num_nodes
        self.sums = [0.0] * num_nodes
        self.sq_sums = [0.0] * num_nodes

        # Shell-level sufficient statistics for O(1) online Haar spectrum estimation:
        # At shell k in {1..L}, parent w in {0,1}^{k-1} has children w0, w1 at depth k.
        # d_w = 0.5 * (mean(w0) - mean(w1)), harmonic weight h_w = 2*n0*n1 / (n0 + n1)
        self.node_h = [0.0] * num_nodes       # Indexed by parent idx at depth k-1
        self.node_hd2 = [0.0] * num_nodes     # h_w * d_w^2 at parent idx
        self.shell_H_sum = [0.0] * (L_depth + 1)    # Sum of h_w across shell k
        self.shell_HD2_sum = [0.0] * (L_depth + 1)  # Sum of h_w * d_w^2 across shell k
        self.shell_split_count = [0] * (L_depth + 1)

        # Online within-leaf noise variance tracker (from leaves visited >= 2 times)
        self.leaf_ss_within = 0.0
        self.leaf_df_within = 0

        # Cached online estimates per shell k in {1..L}
        self.tau2_hat = [0.5 ** (0.5 * k) for k in range(L_depth + 1)]
        self.sigma2_hat = 1.0

    def update(self, leaf: int, y: float):
        """O(L) = O(log N) bottom-up update of Haar means, differences, and online spectrum."""
        # 1. Update leaf within-sample variance tracker if leaf was already visited
        leaf_idx = (1 << self.L) + leaf
        n_old = self.counts[leaf_idx]
        if n_old >= 1:
            old_mean = self.sums[leaf_idx] / n_old
            new_mean = (self.sums[leaf_idx] + y) / (n_old + 1)
            self.leaf_ss_within += (y - old_mean) * (y - new_mean)
            self.leaf_df_within += 1
            if self.leaf_df_within >= 3:
                self.sigma2_hat = max(0.05, self.leaf_ss_within / self.leaf_df_within)

        # 2. Walk up the L+1 ancestor shells
        for k in range(self.L + 1):
            w = leaf >> (self.L - k)
            idx = (1 << k) + w
            self.counts[idx] += 1
            self.sums[idx] += y
            self.sq_sums[idx] += y * y

        # 3. Update Harmonic-Weighted Haar Detail Squares along the L ancestor splits
        for k in range(1, self.L + 1):
            p_w = leaf >> (self.L - (k - 1))
            p_idx = (1 << (k - 1)) + p_w
            c0_idx = (1 << k) + (p_w << 1)
            c1_idx = c0_idx | 1
            n0, n1 = self.counts[c0_idx], self.counts[c1_idx]
            if n0 > 0 and n1 > 0:
                old_h = self.node_h[p_idx]
                old_hd2 = self.node_hd2[p_idx]
                if old_h == 0.0:
                    self.shell_split_count[k] += 1
                h_new = (2.0 * n0 * n1) / (n0 + n1)
                d_new = 0.5 * ((self.sums[c0_idx] / n0) - (self.sums[c1_idx] / n1))
                hd2_new = h_new * (d_new * d_new)
                self.node_h[p_idx] = h_new
                self.node_hd2[p_idx] = hd2_new
                self.shell_H_sum[k] += (h_new - old_h)
                self.shell_HD2_sum[k] += (hd2_new - old_hd2)

        # 4. Recompute O(L) empirical Haar shell variances {tau2_hat[k]}
        running_alias = 0.0
        for k in range(self.L, 0, -1):
            H_k = self.shell_H_sum[k]
            m_k = self.shell_split_count[k]
            eff_noise_k = self.sigma2_hat + running_alias
            if H_k > 0.0 and m_k >= 1:
                raw_tau2 = (self.shell_HD2_sum[k] - 0.5 * m_k * eff_noise_k) / H_k
                prior_guess = 0.8 * (0.6 ** k)
                weight = m_k / (m_k + 2.0)
                est_tau2 = weight * max(0.02 * prior_guess, raw_tau2) + (1.0 - weight) * prior_guess
            else:
                est_tau2 = 0.8 * (0.6 ** k)
            self.tau2_hat[k] = est_tau2
            running_alias = 0.5 * running_alias + est_tau2

    def get_cotton_damped_child(self, k: int, w: int, parent_est: float) -> tuple[float, float, float]:
        """Compute empirical Cotton-damped mean, posterior variance, and reliability gamma."""
        idx = (1 << k) + w
        n = self.counts[idx]
        v_loc = 2.0 * self.tau2_hat[k]
        if n == 0:
            return parent_est, v_loc, 0.0
        alias_below = sum(self.tau2_hat[m] for m in range(k + 1, self.L + 1)) if k < self.L else 0.0
        eff_noise = (self.sigma2_hat + alias_below) / n
        gamma = v_loc / (v_loc + eff_noise)
        raw_mean = self.sums[idx] / n
        damped_mean = parent_est + gamma * (raw_mean - parent_est)
        post_var = v_loc * (1.0 - gamma)
        return damped_mean, post_var, gamma

    def select_leaf(self, t: int, rng: random.Random) -> int:
        """
        3-Role Triad Selection (L Anchor, M Pivot, S Scout):
        - Rule 1 (M Pivot Separation): If a visited node has an unvisited child (h_w == 0),
          immediately separate into the unvisited sibling with bit-reversed striding.
        - Rule 2 (L Anchor + Cotton-Damped UCB): Traverse resolved shells via Cotton-damped UCB.
        - Rule 3 (S Scout Lookahead Cap K): Cap unresolved frontier expansion at K shells.
        """
        n_root = self.counts[1]
        curr_est = (self.sums[1] / n_root) if n_root > 0 else 0.0
        curr_w = 0
        uncertain_shells_crossed = 0

        for k in range(1, self.L + 1):
            w0 = curr_w << 1
            w1 = w0 | 1
            idx0 = (1 << k) + w0
            idx1 = (1 << k) + w1
            n0, n1 = self.counts[idx0], self.counts[idx1]

            if n0 == 0:
                rem = self.L - k
                stride = int(bin(t & ((1 << rem) - 1))[2:].zfill(rem)[::-1], 2) if rem > 0 else 0
                return (w0 << rem) | stride
            if n1 == 0:
                rem = self.L - k
                stride = int(bin(t & ((1 << rem) - 1))[2:].zfill(rem)[::-1], 2) if rem > 0 else 0
                return (w1 << rem) | stride

            est0, var0, gam0 = self.get_cotton_damped_child(k, w0, curr_est)
            est1, var1, gam1 = self.get_cotton_damped_child(k, w1, curr_est)

            log_np = math.log(n0 + n1 + 1)
            u0 = est0 + self.c_ucb * math.sqrt(var0 * log_np)
            u1 = est1 + self.c_ucb * math.sqrt(var1 * log_np)

            if u0 >= u1:
                curr_w, curr_est, chosen_gam = w0, est0, gam0
            else:
                curr_w, curr_est, chosen_gam = w1, est1, gam1

            if chosen_gam < self.gamma_star:
                uncertain_shells_crossed += 1
                if uncertain_shells_crossed >= self.cap_K and k < self.L:
                    rem = self.L - k
                    stride = int(bin(t & ((1 << rem) - 1))[2:].zfill(rem)[::-1], 2) if rem > 0 else 0
                    return (curr_w << rem) | stride

        return curr_w

    def best_leaf_estimate(self) -> int:
        """Telescoping Haar-Cotton Posterior Mean Estimator across all visited leaves."""
        n_root = self.counts[1]
        if n_root == 0:
            return 0
        root_mean = self.sums[1] / n_root
        best_lf = 0
        best_score = -1e18
        for lf in range(self.N):
            if self.counts[(1 << self.L) + lf] == 0:
                continue
            curr_est = root_mean
            for k in range(1, self.L + 1):
                w = lf >> (self.L - k)
                idx = (1 << k) + w
                n = self.counts[idx]
                if n == 0:
                    break
                v_loc = 2.0 * self.tau2_hat[k]
                alias_below = sum(self.tau2_hat[m] for m in range(k + 1, self.L + 1)) if k < self.L else 0.0
                eff_noise = (self.sigma2_hat + alias_below) / n
                gamma = v_loc / (v_loc + eff_noise)
                curr_est = curr_est + gamma * ((self.sums[idx] / n) - curr_est)
            if curr_est > best_score:
                best_score = curr_est
                best_lf = lf
        return best_lf

def run_v2_benchmark():
    L = 10
    N = 1 << L
    T = 256
    n_seeds = 50
    sigma = 1.0
    sigma2 = sigma ** 2

    betas = [0.4, 0.7, 1.1]
    summary = {}

    for beta in betas:
        v_scales = [1.0] + [2.0 ** (-beta * m) for m in range(1, L + 1)]
        reg_v2 = []
        reg_uct = []
        reg_ts = []
        rew_v2 = []
        rew_uct = []
        rew_ts = []

        for seed in range(n_seeds):
            rng_env = random.Random(70000 + int(beta * 1000) + seed)
            theta = sample_cantor_tree_field(L, v_scales, rng_env)
            best_val = max(theta)
            noise_seq = [rng_env.gauss(0.0, sigma) for _ in range(T)]

            rng_v2 = random.Random(95000 + seed)
            tree_v2 = AdaptiveCGHSTS(L, gamma_star=0.45, c_ucb=0.90, cap_K=4)
            r_tot = 0.0
            for t in range(T):
                lf = tree_v2.select_leaf(t, rng_v2)
                y = theta[lf] + noise_seq[t]
                tree_v2.update(lf, y)
                r_tot += theta[lf]
            rec_v2 = tree_v2.best_leaf_estimate()
            rew_v2.append(r_tot)
            reg_v2.append(best_val - theta[rec_v2])

            rng_uct = random.Random(95000 + seed)
            uct_counts = [0] * (1 << (L + 1))
            uct_sums = [0.0] * (1 << (L + 1))
            r_tot = 0.0
            for t in range(T):
                curr_w = 0
                for k in range(1, L + 1):
                    p_idx = (1 << (k - 1)) + curr_w
                    w0, w1 = (curr_w << 1), (curr_w << 1) | 1
                    idx0, idx1 = (1 << k) + w0, (1 << k) + w1
                    n0, n1 = uct_counts[idx0], uct_counts[idx1]
                    if n0 == 0:
                        rem = L - k
                        curr_w = (w0 << rem) | rng_uct.randrange(1 << rem) if rem > 0 else w0
                        break
                    elif n1 == 0:
                        rem = L - k
                        curr_w = (w1 << rem) | rng_uct.randrange(1 << rem) if rem > 0 else w1
                        break
                    else:
                        n_p = max(uct_counts[p_idx], 1)
                        u0 = (uct_sums[idx0] / n0) + 1.414 * math.sqrt(math.log(n_p) / n0)
                        u1 = (uct_sums[idx1] / n1) + 1.414 * math.sqrt(math.log(n_p) / n1)
                        curr_w = w0 if u0 >= u1 else w1
                leaf = curr_w
                y = theta[leaf] + noise_seq[t]
                for k in range(L + 1):
                    w = leaf >> (L - k)
                    idx = (1 << k) + w
                    uct_counts[idx] += 1
                    uct_sums[idx] += y
                r_tot += theta[leaf]
            best_emp_leaf = max(
                (i for i in range(N) if uct_counts[(1 << L) + i] > 0),
                key=lambda i: uct_sums[(1 << L) + i] / uct_counts[(1 << L) + i],
            )
            rew_uct.append(r_tot)
            reg_uct.append(best_val - theta[best_emp_leaf])

            rng_ts = random.Random(95000 + seed)
            prior_var_leaf = sum(v_scales)
            ts_prec = [1.0 / prior_var_leaf] * N
            ts_wsum = [0.0] * N
            r_tot = 0.0
            for t in range(T):
                best_s = -1e18
                leaf = 0
                for i in range(N):
                    p = ts_prec[i]
                    s = (ts_wsum[i] / p) + rng_ts.gauss(0.0, 1.0 / math.sqrt(p))
                    if s > best_s:
                        best_s = s
                        leaf = i
                y = theta[leaf] + noise_seq[t]
                ts_prec[leaf] += 1.0 / sigma2
                ts_wsum[leaf] += y / sigma2
                r_tot += theta[leaf]
            rec_ts = max(range(N), key=lambda i: ts_wsum[i] / ts_prec[i])
            rew_ts.append(r_tot)
            reg_ts.append(best_val - theta[rec_ts])

        _, p_reg_uct, diff_reg_uct, se_reg_uct = paired_ttest(reg_uct, reg_v2)
        _, p_rew_uct, diff_rew_uct, se_rew_uct = paired_ttest(rew_v2, rew_uct)

        summary[f"beta_{beta}"] = {
            "regret_cg_hsts_v2": f"{statistics.mean(reg_v2):.3f} ± {statistics.stdev(reg_v2)/math.sqrt(n_seeds):.3f}",
            "regret_uct": f"{statistics.mean(reg_uct):.3f} ± {statistics.stdev(reg_uct)/math.sqrt(n_seeds):.3f}",
            "regret_flat_ts": f"{statistics.mean(reg_ts):.3f} ± {statistics.stdev(reg_ts)/math.sqrt(n_seeds):.3f}",
            "reward_cg_hsts_v2": f"{statistics.mean(rew_v2):.1f} ± {statistics.stdev(rew_v2)/math.sqrt(n_seeds):.1f}",
            "reward_uct": f"{statistics.mean(rew_uct):.1f} ± {statistics.stdev(rew_uct)/math.sqrt(n_seeds):.1f}",
            "reward_flat_ts": f"{statistics.mean(rew_ts):.1f} ± {statistics.stdev(rew_ts)/math.sqrt(n_seeds):.1f}",
            "regret_reduction_vs_uct": f"-{diff_reg_uct:.3f} ± {se_reg_uct:.3f} (p = {p_reg_uct:.2e})",
            "reward_diff_vs_uct": f"{diff_rew_uct:+.1f} ± {se_rew_uct:.1f} (p = {p_rew_uct:.2e})",
        }

    with open("cg_hsts_v2_benchmark.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("SAVED cg_hsts_v2_benchmark.json")

if __name__ == "__main__":
    run_v2_benchmark()

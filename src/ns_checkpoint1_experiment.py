"""First geometric-restraint checkpoint. Python 3, NumPy, SciPy, Matplotlib.

Unforced incompressible Navier-Stokes on [0,2*pi)^3. Fourier Galerkin
discretization, strict componentwise 2/3 truncation, rotational nonlinearity,
Leray projection, classical explicit RK4. This is a diagnostic solver, not a
regularity proof or a high-Reynolds-number turbulence calculation.
Run: python experiment.py --out results
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.fft import fftn, ifftn

AXES = (-3, -2, -1)


class Solver:
    def __init__(self, n, nu):
        self.n, self.nu = n, nu
        k = np.fft.fftfreq(n, 1/n)
        self.k = np.array(np.meshgrid(k, k, k, indexing='ij'))
        self.k2 = (self.k**2).sum(axis=0)
        self.den = np.where(self.k2 == 0, 1, self.k2)
        self.mask = (np.abs(self.k) < n/3).all(axis=0)
        self.xyz = np.meshgrid(np.arange(n)*2*np.pi/n,
                               np.arange(n)*2*np.pi/n,
                               np.arange(n)*2*np.pi/n, indexing='ij')

    def ft(self, a):
        return fftn(a, axes=AXES)

    def real(self, a):
        return ifftn(a, axes=AXES).real

    def project(self, h):
        return (h - self.k*(self.k*h).sum(axis=0)/self.den)*self.mask

    def curl(self, h):
        return 1j*np.cross(self.k, h, axisa=0, axisb=0, axisc=0)

    def initial(self, family, sign=1):
        x,y,z = self.xyz
        zero = np.zeros_like(x)
        if family == 'witness':
            u = np.array([np.sin(y), zero, np.sin(x)-np.sin(x+y)])
        elif family == '2d':
            u = np.array([np.sin(x)*np.cos(y), -np.cos(x)*np.sin(y), zero])
        elif family == 'abc':
            u = np.array([np.sin(z)+np.cos(y), np.sin(x)+np.cos(z),
                          np.sin(y)+np.cos(x)])
        elif family == 'tg':
            u = np.array([np.sin(x)*np.cos(y)*np.cos(z),
                          -np.cos(x)*np.sin(y)*np.cos(z), zero])
        else:
            raise ValueError(family)
        return self.project(self.ft(sign*u))

    def rhs(self, h):
        u, w = self.real(h), self.real(self.curl(h))
        nonlinear = self.ft(np.cross(u,w,axisa=0,axisb=0,axisc=0))
        return self.project(nonlinear) - self.nu*self.k2*h

    def step(self,h,dt):
        a = self.rhs(h)
        b = self.rhs(h+dt*a/2)
        c = self.rhs(h+dt*b/2)
        d = self.rhs(h+dt*c)
        return self.project(h+dt*(a+2*b+2*c+d)/6)

    def diagnostics(self,h,t,dt):
        u, wh = self.real(h), self.curl(h)
        w = self.real(wh)
        q2 = (w*w).sum(axis=0)
        q = np.sqrt(q2)
        grad = self.real(1j*h[:,None]*self.k[None,:])
        strain = (grad+np.swapaxes(grad,0,1))/2
        p = np.einsum('iabc,ijabc,jabc->abc',w,strain,w)
        gw = self.real(1j*wh[:,None]*self.k[None,:])
        D = self.nu*np.mean((gw*gw).sum(axis=(0,1)))
        E = np.mean((u*u).sum(axis=0))/2
        W = np.mean(q2)
        P = p.mean()
        eig = np.linalg.eigvalsh(np.moveaxis(strain,(0,1),(-2,-1)))
        snorm = np.abs(eig).max(axis=-1)
        denom = np.mean(q2*snorm)
        rhs = self.rhs(h)
        Edot = np.mean(np.sum(u*self.real(rhs),axis=0))
        Zdot = np.mean(np.sum(w*self.real(self.curl(rhs)),axis=0))
        threshold = np.quantile(q,0.75)
        # Include tied percentile values consistently under roundoff.
        high = (q >= threshold-1e-12*q.max()) & (q > 1e-10*q.max())
        xi = np.divide(w,q[None],out=np.zeros_like(w),where=q[None]>1e-12*q.max())
        # Same physical separation at every resolution: pi/4 in each axis.
        nums, dens, smooth_num, smooth_den = 0.,0.,0.,0.
        for axis in range(3):
            for shift in [-self.n//8,self.n//8]:
                other = np.roll(xi,shift,axis=axis+1)
                both = high & np.roll(high,shift,axis=axis)
                weight = q2*np.roll(q2,shift,axis=axis)*both
                sin2 = np.clip(1-(xi*other).sum(axis=0)**2,0,1)
                nums += np.sum(weight*sin2)
                dens += np.sum(weight)
                # Threshold-free, intensity-weighted alternative. This form
                # avoids dividing by |omega| at its zeros altogether.
                shifted_w = np.roll(w,shift,axis=axis+1)
                base_weight = q2*np.roll(q2,shift,axis=axis)
                dot2 = np.sum(w*shifted_w,axis=0)**2
                smooth_num += np.sum(np.maximum(base_weight-dot2,0))
                smooth_den += np.sum(base_weight)
        spec = (np.abs(h)**2).sum(axis=0)
        edge = np.max(np.abs(self.k),axis=0) >= 0.8*(self.n/3)
        return dict(t=float(t), E=float(E), Z=float(W/2), W=float(W),
                    qmax=float(q.max()), q95=float(np.quantile(q,.95)),
                    concentration=float(np.mean(q2*q2)/W**2),
                    P=float(P), Pplus=float(np.maximum(p,0).mean()),
                    Pminus=float(np.maximum(-p,0).mean()), D=float(D),
                    net_ratio=float(P/D) if D else 0.,
                    eta_plus=float(np.maximum(p,0).mean()/denom) if denom else 0.,
                    coherence_defect=float(nums/dens) if dens else 0.,
                    coherence_defect_smooth=float(smooth_num/smooth_den) if smooth_den else 0.,
                    divergence_rms=float(np.sqrt(np.mean(self.real(1j*(self.k*h).sum(axis=0))**2))),
                    energy_budget_rel=float(abs(Edot+self.nu*W)/max(self.nu*W,1e-30)),
                    enstrophy_budget_rel=float(abs(Zdot+D-P)/max(D,abs(P),1e-30)),
                    edge_energy_fraction=float(spec[edge].sum()/spec.sum()),
                    advective_cfl=float(dt*np.max(np.sum(abs(u),axis=0))/(2*np.pi/self.n)),
                    diffusive_number=float(dt*self.nu*self.k2[self.mask].max()))


def run(family,n,nu,dt,T=1.,sign=1):
    solver = Solver(n,nu)
    h0 = solver.initial(family,sign)
    h = h0.copy()
    steps = int(round(T/dt))
    interval = max(1,steps//10)
    rows = [solver.diagnostics(h,0.,dt)]
    for j in range(1,steps+1):
        h = solver.step(h,dt)
        if j%interval == 0 or j == steps:
            rows.append(solver.diagnostics(h,j*dt,dt))
    result = dict(family=family,n=n,nu=nu,dt=dt,T=T,sign=sign,rows=rows)
    if family in ('2d','abc'):
        k2 = 2 if family == '2d' else 1
        exact = h0*np.exp(-nu*k2*T)
        result['exact_relative_L2_error'] = float(np.linalg.norm(h-exact)/np.linalg.norm(exact))
    return result,h


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out',default='results')
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    start = time.time()
    cases = [('2d',16,.05,.01,1),('abc',16,.05,.01,1),
             ('witness',24,.05,.01,1),('witness',24,.05,.01,-1),
             ('witness',32,.05,.005,1),
             ('tg',24,.05,.01,1),('tg',32,.05,.01,1),
             ('tg',32,.05,.005,1),('tg',48,.05,.01,1),
             ('tg',32,.01,.005,1),('tg',48,.01,.005,1)]
    runs,finals = [],{}
    for family,n,nu,dt,sign in cases:
        r,h = run(family,n,nu,dt,sign=sign)
        runs.append(r)
        finals[(family,n,nu,dt,sign)] = h/n**3
        print(f'{family} N={n} nu={nu} dt={dt} sign={sign}: '
              f'Z {r["rows"][0]["Z"]:.6g} -> {r["rows"][-1]["Z"]:.6g}',flush=True)

    # Compare final velocity fields in the full finer Fourier space. Coarse
    # fields are embedded by integer wavevector, with zero unresolved modes.
    def compare(a,b):
        small,large = finals[a],finals[b]
        n,m = small.shape[-1],large.shape[-1]
        kval = np.rint(np.fft.fftfreq(n,1/n)).astype(int)
        idx = kval%m
        embedded = np.zeros_like(large)
        embedded[:,idx[:,None,None],idx[None,:,None],idx[None,None,:]] = small
        return float(np.linalg.norm(embedded-large)/np.linalg.norm(large))

    comparisons = {
        'tg_nu005_N24_vs_N32':compare(('tg',24,.05,.01,1),('tg',32,.05,.01,1)),
        'tg_nu005_N32_vs_N48':compare(('tg',32,.05,.01,1),('tg',48,.05,.01,1)),
        'tg_nu005_dt001_vs_dt0005':compare(('tg',32,.05,.01,1),('tg',32,.05,.005,1)),
        'tg_nu001_N32_vs_N48':compare(('tg',32,.01,.005,1),('tg',48,.01,.005,1)),
        'witness_N24_dt001_vs_N32_dt0005':compare(('witness',24,.05,.01,1),('witness',32,.05,.005,1)),
    }
    positive = runs[2]['rows'][0]
    negative = runs[3]['rows'][0]
    validation = dict(
        witness_P_error=abs(positive['P']-.25),
        witness_W_error=abs(positive['W']-2.),
        witness_D_error=abs(positive['D']-.15),
        reversal_P_sum=abs(positive['P']+negative['P']),
        reversal_coherence_difference=abs(positive['coherence_defect']-negative['coherence_defect']),
        abc_coherence_drift=max(abs(row['coherence_defect']-runs[1]['rows'][0]['coherence_defect']) for row in runs[1]['rows']),
        abc_smooth_coherence_drift=max(abs(row['coherence_defect_smooth']-runs[1]['rows'][0]['coherence_defect_smooth']) for row in runs[1]['rows']),
        max_energy_budget_rel=max(r['energy_budget_rel'] for run_ in runs for r in run_['rows']),
        max_enstrophy_budget_rel=max(r['enstrophy_budget_rel'] for run_ in runs for r in run_['rows']),
        max_divergence_rms=max(r['divergence_rms'] for run_ in runs for r in run_['rows']),
        max_advective_cfl=max(r['advective_cfl'] for run_ in runs for r in run_['rows']),
        max_diffusive_number=max(r['diffusive_number'] for run_ in runs for r in run_['rows']),
    )
    # Gates concern the chosen smooth, short-time controls only.
    assert validation['witness_P_error'] < 1e-12
    assert validation['witness_W_error'] < 1e-12
    assert validation['witness_D_error'] < 1e-12
    assert validation['abc_coherence_drift'] < 1e-12
    assert validation['abc_smooth_coherence_drift'] < 1e-12
    assert validation['max_energy_budget_rel'] < 1e-10
    assert validation['max_enstrophy_budget_rel'] < 1e-10
    assert validation['max_divergence_rms'] < 1e-10
    assert max(r.get('exact_relative_L2_error',0) for r in runs) < 1e-9
    assert comparisons['tg_nu005_N32_vs_N48'] < 1e-4
    assert comparisons['tg_nu005_dt001_vs_dt0005'] < 1e-7
    result = dict(domain='[0,2*pi)^3; integrals normalized by domain volume',
                  method='strict componentwise 2/3 truncation; rotational RHS; Leray projection; RK4',
                  elapsed_seconds=time.time()-start,validation=validation,
                  comparisons=comparisons,runs=runs)
    (out/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(validation=validation,comparisons=comparisons),indent=2),flush=True)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs = plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for idx,label,color in [(2,'Original flow','#2563eb'),(3,'Velocity reversed','#b45309')]:
        rows=runs[idx]['rows']
        axs[0].plot([r['t'] for r in rows],[r['Z'] for r in rows],label=label,color=color,lw=2.5)
    axs[0].set(title='Identical initial coherence, different evolution',xlabel='Time',ylabel='Enstrophy per unit volume')
    axs[0].legend(frameon=False)
    for idx,label,color in [(8,'Viscosity 0.05','#2563eb'),(10,'Viscosity 0.01','#b45309')]:
        rows=runs[idx]['rows']
        axs[1].plot([r['t'] for r in rows],[r['P']/r['D'] for r in rows],label=label,color=color,lw=2.5)
    axs[1].axhline(1,color='#64748b',ls='--',lw=1,label='Net growth threshold')
    axs[1].set(title='3D Taylor–Green control at 48 cubed',xlabel='Time',ylabel='Net stretching / viscous dissipation')
    axs[1].legend(frameon=False)
    fig.savefig(out/'first_controls.png',dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()

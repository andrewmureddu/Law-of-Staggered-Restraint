"""Generate the readable checkpoint from experiment.py's saved diagnostics."""
import json
from pathlib import Path

root = Path(__file__).resolve().parent
d = json.loads((root/'results/results.json').read_text())
runs = d['runs']
lines = ['# Numerical results for the first checkpoint', '',
         'All trajectories end at t = 1 on the periodic box [0, 2π)³. '
         'All integrals below are divided by the box volume.', '',
         '## Initial reversal check', '',
         '| Initial field | Enstrophy Z | Net production P | Dissipation D | Initial Z derivative |',
         '| --- | ---: | ---: | ---: | ---: |']
for i,label in [(2,'u₀'),(3,'−u₀')]:
    r=runs[i]['rows'][0]
    lines.append(f'| {label} | {r["Z"]:.6f} | {r["P"]:.6f} | {r["D"]:.6f} | {r["P"]-r["D"]:.6f} |')
lines += ['', 'The two fields have identical vorticity magnitude point by point. '
          'Their concentration, spectra, pairwise directional geometry, and '
          'dissipation agree exactly in the analytical calculation. '
          'The code recovers the signed-production reversal.', '',
          '## Short trajectories', '',
          '| Family | N | Viscosity | Step | Sign | Initial Z | Final Z | Final P divided by D |',
          '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
for run in runs:
    a,b=run['rows'][0],run['rows'][-1]
    lines.append(f'| {run["family"]} | {run["n"]} | {run["nu"]} | {run["dt"]} | {run["sign"]:+d} | {a["Z"]:.6f} | {b["Z"]:.6f} | {b["net_ratio"]:.6f} |')
gain=(runs[2]['rows'][-1]['Z']/runs[2]['rows'][0]['Z']-1)*100
loss=(runs[3]['rows'][-1]['Z']/runs[3]['rows'][0]['Z']-1)*100
lines += ['', f'The original shear example changes enstrophy by {gain:+.2f} percent; '
          f'its reversed counterpart changes by {loss:+.2f} percent. '
          'The reversed run begins losing enstrophy and later starts gaining it; '
          'the sign at the initial instant does not fix its entire future.', '',
          'The lower-viscosity three-dimensional Taylor–Green run develops '
          'net stretching larger than viscous dissipation over part of this '
          'short interval. This is ordinary finite-time amplification in a '
          'resolved calculation, not evidence for a singularity.', '',
          '## Solver checks', '', '| Check | Observed value |', '| --- | ---: |']
for run in runs[:2]:
    lines.append(f'| {run["family"]} exact final velocity relative L2 error | {run["exact_relative_L2_error"]:.3e} |')
for k,v in d['validation'].items():
    lines.append(f'| {k} | {v:.3e} |')
lines += ['', 'The exact-flow errors are at floating-point scale for these small '
          'steps and low-wavenumber decays; they do not establish temporal order '
          'for general flows. The budget residuals are instantaneous algebra '
          'checks. The following comparisons separately assess final fields.', '',
          '| Final velocity comparison | Relative L2 difference |', '| --- | ---: |']
for k,v in d['comparisons'].items():
    lines.append(f'| {k} | {v:.3e} |')
lines += ['', 'The simultaneous grid and time-step change for the shear witness '
          'does not isolate those two errors. For three-dimensional Taylor–Green '
          'at viscosity 0.05, they are varied separately.', '',
          '## Geometry diagnostic checks', '',
          'A sharp percentile mask is sensitive to tied values and spatial sampling. '
          'The final implementation includes a relative 10⁻¹² tolerance at the '
          'percentile boundary to prevent roundoff from arbitrarily removing '
          'exactly tied points. It also records a threshold-free, intensity-weighted '
          'coherence defect. Both remain diagnostics, not regularity conditions.', '',
          '| Final diagnostic | N 32 | N 48 | Relative difference |',
          '| --- | ---: | ---: | ---: |']
for i,j in [(6,8),(9,10)]:
    a,b=runs[i]['rows'][-1],runs[j]['rows'][-1]
    for key in ['coherence_defect','coherence_defect_smooth','eta_plus','concentration']:
        diff=abs(a[key]-b[key])/max(abs(b[key]),1e-30)
        lines.append(f'| {key}, viscosity {runs[i]["nu"]} | {a[key]:.9f} | {b[key]:.9f} | {diff:.3e} |')
lines += ['', 'Velocity-field convergence does not guarantee that a sharply '
          'thresholded diagnostic has converged to the same relative accuracy. '
          'Use the smoother quantity for the first ensemble analysis and keep '
          'the percentile statistic as a sensitivity check. The smoother measure '
          'includes all pairs and is not equivalent to coherence restricted '
          'to the most intense region.', '',
          '## What remains open', '',
          'This suite does not establish that increased concentration causes '
          'geometric restraint, or that geometry predicts later depletion better '
          'than current strain and intensity. Multiple flow families, matched '
          'phase controls, longer resolved trajectories, and held-out prediction '
          'tests remain necessary. No broad negative conclusion about H1 or H2 '
          'is justified by the control suite alone.', '',
          'See README.md for the analytical control, scaling constraint, '
          'geometry/dissipation identity, references, and the next experiment.']
(root/'RESULTS.md').write_text('\n'.join(lines)+'\n')
print(root/'RESULTS.md')

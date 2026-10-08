# Exoplanet Transit Detection - Kepler-6b (KIC 10874614)
# Requires: pip install lightkurve batman-package

import numpy as np
import lightkurve as lk
import batman
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Data: Kepler Quarter 2,2009, 30-minute cadence
search=lk.search_lightcurve ("Kepler-6", mission="Kepler")
lc=search[8].download()

# Cleaning
lc_clean=lc.remove_nans() # remove missing measurements
lc_q=lc_clean[lc_clean.quality==0] # keep only unflagged points
lc_norm=lc_q.normalize() # flux divided by the median
# Upper-side clipping only: a transit is a dip, so low points are the signal
lc_no = lc_norm.remove_outliers(sigma_upper=4, sigma_lower=np.inf)
print('points: raw', len(lc), '| without NaN', len(lc_clean),
      '| quality==0', len(lc_q), '| after clipping', len(lc_no))

# Period search: wide Box Least Squares, no prior on the period
periods = np.linspace(0.5, 20, 20000)          # trial periods [d]
bls = lc_no.to_periodogram(method='bls', period=periods,
                           duration=[0.05, 0.1, 0.15, 0.2])   # trial durations [d]
print('BLS period (wide search):', bls.period_at_max_power)

# Period refinement: narrow BLS around the first peak
P0 = bls.period_at_max_power.value
fine = np.linspace(P0 - 0.01, P0 + 0.01, 4001)  # step 5e-6 d
bls2 = lc_no.to_periodogram(method='bls', period=fine,
                            duration=[0.10, 0.13, 0.16, 0.19])
P2 = bls2.period_at_max_power.value            # orbital period [d]
print('refined period [d]:', P2)

# Folded light curve and first estimate: depth ~ (Rp/R*)^2
folded = lc_no.fold(period=bls2.period_at_max_power,
                    epoch_time=bls2.transit_time_at_max_power)
ph = np.array(folded.time.value, dtype=float)  # phase [d]
f = np.array(folded.flux.value, dtype=float)   # normalized flux
depth_simple = np.median(f[np.abs(ph) > 0.2]) - np.median(f[np.abs(ph) < 0.03])
print('depth (medians):', depth_simple)
print('Rp/R* = sqrt(depth):', np.sqrt(depth_simple))

# Transit model: batman, linear limb darkening, circular orbit
def model(t, t0, rp, a, b, u):
    q = batman.TransitParams()
    q.t0, q.per, q.rp, q.a, q.ecc, q.w = t0, P2, rp, a, 0., 90.
    q.inc = np.degrees(np.arccos(b / a))       # inclination from impact parameter
    q.limb_dark, q.u = 'linear', [u]
    # each point is a 30-minute exposure (0.0204 d): supersample it
    return batman.TransitModel(q, t, supersample_factor=7,
                               exp_time=0.0204).light_curve(q)
# Fit: mid-transit time, Rp/R*, a/R*, impact parameter, limb darkening
k0 = np.sqrt(depth_simple * (1 - 0.6 / 3))     # starting value for Rp/R*
m = np.abs(ph) < 0.3                           # excludes the instrumental ramp
popt, pcov = curve_fit(model, ph[m], f[m], p0=[0., k0, 7., 0.1, 0.6],
                       bounds=([-0.02, 0.05, 3, 0, 0], [0.02, 0.15, 15, 0.9, 1]))
err = np.sqrt(np.diag(pcov))
for name, v, e in zip(['t0 [d]', 'Rp/R*', 'a/R*', 'b', 'u'], popt, err):
    print(f'{name:7s} = {v:.5f} +/- {e:.5f}')

# Residuals
res = f - model(ph, *popt)
print('residual std, in transit:', np.std(res[np.abs(ph) < 0.1]))
print('residual std, out of transit:', np.std(res[np.abs(ph) > 0.2]))

# Results and comparison with reference values
depth_m = 1 - model(np.array([0.]), *popt)[0]
print(f'Period: {P2:.5f} d     (reference 3.2347)')
print(f'Depth : {depth_m:.5f}      (reference 0.0106)')
print(f'Rp/R* : {popt[1]:.4f} +/- {err[1]:.4f}  (reference 0.093)')

# Plot the transit and the residuals
fig, (a1, a2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6),
                             gridspec_kw={'height_ratios': [3, 1]})

o = np.argsort(ph)
a1.scatter(ph, f, s=2, color='k', label='Kepler Q2')   # data
a1.plot(ph[o], model(ph[o], *popt), 'r', label='model')  # fitted model
a2.scatter(ph, res, s=2, color='k')
a2.axhline(0, color='r')
a1.set_xlim(-0.3, 0.3)
a2.set_ylim(-0.002, 0.002)
a1.set_ylabel('Normalized flux')
a2.set_ylabel('Residuals')
a2.set_xlabel('Phase [days]')
a1.legend()
fig.savefig('transit.png', dpi=150)
plt.show()
# Exoplanet Transit Detection - Kepler-6b (KIC 10874614)
# Requires: pip install lightkurve batman-package

import numpy as np
import lightkurve as lk
import batman
import matplotlib.pypylot as plt
from scipy.optimize import curve_fit

# Data: Kepler Quarter 2,2009, 30-minute cadence
search=lk.search_lightkurve ("Kepler-6", mission="Kepler")
lc=search[8].download()

# Cleaning
lc_clean=lc.remove_nans() # remove missing measurements
lc_q=lc_clean[lc_clean.quality==0] # keep only unflagged points
lc_norm=lc_q.normalize() # flux divided by the median
# Upper-side clipping only: a transit is a dip, so low points are the signal
lc_no = lc_norm.remove_outliers(sigma_upper=4, sigma_lower=np.inf)
print('points: raw', len(lc), '| without NaN', len(lc_clean),
      '| quality==0', len(lc_q), '| after clipping', len(lc_no))
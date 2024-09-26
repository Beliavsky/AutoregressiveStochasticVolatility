import numpy as np
from arsv import arsv_vol_and_returns, fit_arsv
from stats_np import print_stats

mu = -0.5       # Long-term mean of log volatility
phi = 0.9       # Autoregressive parameter
sigma_v = 0.4   # Volatility of volatility
nobs = 10**6    # # of returns simulated

x = arsv_vol_and_returns(mu, phi, sigma_v, nobs)
print_stats(x[:, 0], title="vol", end="\n")
print_stats(np.log(x[:, 0]), title="log(vol)", end="\n")
print_stats(x[:, 1], title="returns", end="\n")
print_stats(x[:, 1]/x[:, 0], title="returns/vol", end="\n")
print(fit_arsv(x[:, 1]))


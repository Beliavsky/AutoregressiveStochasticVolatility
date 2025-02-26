""" 
Simulate from a stochastic volatility model with nonzero correlation 
and fit an SV model (ignoring the correlation) to the data.
"""
import numpy as np
import pandas as pd
from arsv import arsv_vol_and_returns, fit_arsv_zero_corr
from stats_np import print_stats, stats_names, some_stats
from util import print_vec

fmt_r = "%10.6f"
mu = -0.5       # Long-term mean of log volatility
phi = 0.9       # Autoregressive parameter
sigma_v = 0.4   # Volatility of volatility
rho = -0.5      # Nonzero correlation between return and log volatility innovations
nobs = 10**6    # Number of returns simulated
niter = 3       # Number of sets of returns to simulate

print("#obs:", nobs, end="\n\n")

for iter in range(niter):
    # Simulate returns and volatility with nonzero correlation
    x = arsv_vol_and_returns(mu, phi, sigma_v, rho, nobs)
    
    # Print parameter names and true values (including rho)
    print_vec(["mu", "phi", "sigma", "rho"], "%10s", label=4*" ")
    print_vec([mu, phi, sigma_v, rho], fmt_r, label="true")
    
    # Fit the ARSV model to the simulated returns (ignoring rho)
    print_vec(fit_arsv_zero_corr(x[:, 1]), fmt_r, label=" fit", end="\n")
    
    # Compute and print descriptive statistics for the simulated series
    df = pd.DataFrame({
        "vol": some_stats(x[:, 0]),
        "log(vol)": some_stats(np.log(x[:, 0])),
        "returns": some_stats(x[:, 1]),
        "returns/vol": some_stats(x[:, 1] / x[:, 0])
    }, index=stats_names).T
    print(df, end="\n\n")

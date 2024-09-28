""" Simulate from a stochastic volatility model and fit an SV model to the data """
import numpy as np
import pandas as pd
from arsv import arsv_vol_and_returns, fit_arsv
from stats_np import print_stats, stats_names, some_stats
from util import print_vec

fmt_r = "%10.6f"
mu = -0.5       # Long-term mean of log volatility
phi = 0.9       # Autoregressive parameter
sigma_v = 0.4   # Volatility of volatility
nobs = 10**4    # # of returns simulated
niter = 3       # # of sets of returns to simulate

for iter in range(niter):
    x = arsv_vol_and_returns(mu, phi, sigma_v, nobs)
    print_vec(["mu", "phi", "sigma"], "%10s", label=4*" ")
    print_vec([mu, phi, sigma_v], fmt_r, label="true")
    print_vec(fit_arsv(x[:, 1]), fmt_r, label=" fit", end="\n")
    df = pd.DataFrame({"vol":some_stats(x[:, 0]), "log(vol)":some_stats(np.log(x[:, 0])),
        "returns":some_stats(x[:, 1]), "returns/vol":some_stats(x[:, 1]/x[:, 0])},
        index=stats_names).T
    print(df, end="\n\n")

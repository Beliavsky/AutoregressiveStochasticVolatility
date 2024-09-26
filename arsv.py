import numpy as np
from scipy.special import psi, polygamma

def arsv_returns(mu, phi, sigma_v, T, burn_in=1000):
    """
    Simulate daily returns with autoregressive stochastic log volatility.

    Parameters:
    - mu: float, the long-term mean of the log volatility process.
    - phi: float, the autoregressive parameter (|phi| < 1).
    - sigma_v: float, the standard deviation of the volatility shocks.
    - T: int, the number of time periods to simulate.
    - burn_in: int, the number of initial observations to discard to mitigate the effect of initial values.

    Returns:
    - y: ndarray, the simulated daily returns of length T.
    """
    total_length = T + burn_in
    h = np.zeros(total_length)
    epsilon = np.random.normal(size=total_length)
    v = np.random.normal(size=total_length)
    h[0] = mu  # Initialize h[0] at the long-term mean

    # Generate the log volatility process
    for t in range(1, total_length):
        h[t] = mu + phi * (h[t-1] - mu) + sigma_v * v[t]

    # Generate the returns
    y = np.exp(h / 2) * epsilon

    # Discard the burn-in period
    return y[burn_in:]

def arsv_vol_and_returns(mu, phi, sigma_v, T, burn_in=1000):
    """
    Simulate daily returns and conditional standard deviations with autoregressive stochastic log volatility.

    Parameters:
    - mu: float, the long-term mean of the log volatility process.
    - phi: float, the autoregressive parameter (|phi| < 1).
    - sigma_v: float, the standard deviation of the volatility shocks.
    - T: int, the number of time periods to simulate.
    - burn_in: int, the number of initial observations to discard to mitigate the effect of initial values.

    Returns:
    - result: ndarray of shape (T, 2), where the first column contains the conditional standard deviations
              and the second column contains the simulated daily returns.
    """
    total_length = T + burn_in
    h = np.zeros(total_length)
    epsilon = np.random.normal(size=total_length)
    v = np.random.normal(size=total_length)
    h[0] = mu  # Initialize h[0] at the long-term mean
    # Generate the log volatility process
    for t in range(1, total_length):
        h[t] = mu + phi * (h[t-1] - mu) + sigma_v * v[t]
    # Calculate the conditional standard deviations
    sigma_t = np.exp(h / 2)
    # Generate the returns
    y = sigma_t * epsilon
    # Discard the burn-in period
    sigma_t = sigma_t[burn_in:]
    y = y[burn_in:]
    # Combine sigma_t and y into a single array
    result = np.column_stack((sigma_t, y))
    return result

def fit_arsv(y):
    """
    Estimate the parameters mu, phi, and sigma_v of an ARSV model given observed returns y_t.

    Parameters:
    - y: ndarray, the observed returns.

    Returns:
    - mu_hat: float, estimated mu.
    - phi_hat: float, estimated phi.
    - sigma_v_hat: float, estimated sigma_v.
    """
    # Step 1: Compute s_t = log(y_t^2)
    s_t = np.log(y**2)

    # Step 2: Compute mean and variance of eta_t
    # eta_t = log(epsilon_t^2), where epsilon_t ~ N(0,1)
    # Mean and variance of eta_t can be derived using properties of the log-chi-squared distribution
    mean_eta = psi(0.5) + np.log(2)  # psi is the digamma function
    var_eta = polygamma(1, 0.5)      # polygamma(1, z) is the trigamma function
    # Step 3: Compute s_t_star = s_t - mean_eta
    s_t_star = s_t - mean_eta
    # Step 4: Compute sample variance of s_t_star
    Var_s = np.var(s_t_star, ddof=1)
    # Step 5: Compute sample autocorrelation at lag 1 of s_t_star
    n = len(s_t_star)
    s_t_star_mean = np.mean(s_t_star)
    s_t_star_centered = s_t_star - s_t_star_mean
    autocov = np.sum(s_t_star_centered[1:] * s_t_star_centered[:-1]) / (n - 1)
    Var_s_empirical = np.sum(s_t_star_centered ** 2) / (n - 1)
    rho_1 = autocov / Var_s_empirical
    phi_hat = rho_1 / (1 - var_eta / Var_s)
    phi_hat = np.clip(phi_hat, -0.9999, 0.9999)
    sigma_h2_hat = Var_s - var_eta
    sigma_v2_hat = (1 - phi_hat ** 2) * sigma_h2_hat
    mu_hat = np.mean(s_t_star)
    sigma_v_hat = np.sqrt(sigma_v2_hat)
    return mu_hat, phi_hat, sigma_v_hat

import numpy as np
from scipy.stats import skew, kurtosis

stats_names = ["mean", "std", "skew", "kurtosis", "min", "max"]

def print_stats(x, fmt_s="%10s", fmt_r="%10.4f", title=None, end=None):
    """
    Prints the mean, standard deviation, skewness, and kurtosis of a 1D NumPy array.

    Parameters:
    - x: ndarray, the input array.

    Returns:
    - None
    """
    if title:
        print(title)
    for key, value in {"mean":np.mean(x), "sd":np.std(x, ddof=1), "skew":skew(x),
        "kurtosis":kurtosis(x), "min":np.min(x), "max":np.max(x)}.items():
        print(fmt_s%key, fmt_r%value)
    if end:
        print(end=end)

def dict_stats(x, fmt_s="%10s", fmt_r="%10.4f", title=None, end=None):
    """
    Parameters:
    - x: ndarray, the input array.

    Returns:
    - dict of stats
    """
    return {"mean":np.mean(x), "sd":np.std(x, ddof=1), "skew":skew(x),
        "kurtosis":kurtosis(x), "min":np.min(x), "max":np.max(x)}

def some_stats(x):
    return [np.mean(x), np.std(x, ddof=1), skew(x), kurtosis(x),
        np.min(x), np.max(x)]

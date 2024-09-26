import numpy as np
from scipy.stats import skew, kurtosis

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

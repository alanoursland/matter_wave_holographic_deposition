"""
Image-quality and beam-quality metrics.

All functions accept numpy arrays.
"""

import numpy as np
from scipy.signal import find_peaks


def michelson_contrast(d: np.ndarray, lo: float = 5, hi: float = 95) -> float:
    """
    Michelson contrast estimated from percentiles.

    Parameters
    ----------
    d : np.ndarray
        Intensity distribution (any shape).
    lo, hi : float
        Lower and upper percentiles used for min/max estimation.

    Returns
    -------
    float
        (I_hi - I_lo) / (I_hi + I_lo), in [0, 1].
    """
    lo_v = np.percentile(d, lo)
    hi_v = np.percentile(d, hi)
    return (hi_v - lo_v) / (hi_v + lo_v + 1e-30)


def ssim_score(ref: np.ndarray, test: np.ndarray) -> float:
    """
    Structural Similarity Index (SSIM) between two 2-D arrays.

    Requires scikit-image.  (A silent fallback to Pearson correlation was
    removed — fable5 E7 — because it is a different metric on a different
    scale and made results environment-dependent.)

    Parameters
    ----------
    ref, test : np.ndarray
        Reference and test images (same shape).

    Returns
    -------
    float
        SSIM value in [-1, 1]; 1.0 means identical.
    """
    try:
        from skimage.metrics import structural_similarity
    except ImportError as exc:
        raise ImportError(
            "ssim_score requires scikit-image (pip install scikit-image); "
            "refusing to silently substitute a different metric."
        ) from exc
    r = ref  / (ref.max()  + 1e-30)
    t = test / (test.max() + 1e-30)
    return structural_similarity(r, t, data_range=1.0)


def min_feature_size(density: np.ndarray, x_axis: np.ndarray) -> float:
    """
    Minimum FWHM of intensity peaks along the central column of *density*.

    Parameters
    ----------
    density : np.ndarray, shape (Nx, Ny)
        2-D intensity distribution.
    x_axis : np.ndarray, shape (Nx,)
        Spatial coordinates along axis 0 [m].

    Returns
    -------
    float
        Minimum FWHM in the same units as *x_axis*, or np.nan if no peaks found.
    """
    mid  = density.shape[1] // 2
    line = density[:, mid] / (density[:, mid].max() + 1e-30)
    peaks, _ = find_peaks(line, height=0.3, distance=3)
    if not len(peaks):
        return np.nan
    fwhms = []
    for pk in peaks:
        half = 0.5 * line[pk]
        left = pk
        while left > 0 and line[left] > half:
            left -= 1
        right = pk
        while right < len(line) - 1 and line[right] > half:
            right += 1
        fwhms.append((right - left) * abs(x_axis[1] - x_axis[0]))
    return float(np.min(fwhms)) if fwhms else np.nan

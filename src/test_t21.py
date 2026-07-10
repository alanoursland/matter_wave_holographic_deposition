"""Tests for T21 closed-loop stochastic printing (pure components)."""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from t21_closed_loop_printing import (
    coarse_grain, site_metrics, sample_arrivals,
)


class TestComponents:

    def test_coarse_grain_conserves_mass(self):
        rng = np.random.default_rng(0)
        A = rng.random((256, 256))
        C = coarse_grain(A, 32)
        assert C.shape == (8, 8)
        assert abs(C.sum() - A.sum()) < 1e-9

    def test_sample_arrivals_statistics(self):
        rng = np.random.default_rng(1)
        p = np.zeros((64, 64))
        p[10, 10] = 0.75
        p[50, 50] = 0.25
        counts = sample_arrivals(p, 100000, rng)
        assert abs(counts.sum() - 100000) < 4 * np.sqrt(100000)
        assert abs(counts[10, 10] / counts.sum() - 0.75) < 0.02

    def test_site_metrics_perfect_dose_is_zero_error(self):
        t = np.zeros((8, 8))
        t[2, 2] = t[5, 5] = 0.5
        counts = t * 1000
        rms, defect, rms_cal = site_metrics(counts, t, 1000)
        assert rms < 1e-12
        assert defect == 0.0
        assert rms_cal < 1e-12

    def test_site_metrics_flags_misallocation(self):
        t = np.zeros((8, 8))
        t[2, 2] = t[5, 5] = 0.5
        counts = np.zeros((8, 8))
        counts[2, 2] = 900.0   # +80%
        counts[5, 5] = 100.0   # -80%
        rms, defect, rms_cal = site_metrics(counts, t, 1000)
        assert rms > 0.5
        assert defect == 1.0
        # anti-symmetric misallocation is a shape error — calibration
        # must NOT hide it
        assert rms_cal > 0.5

    def test_calibrated_metric_removes_common_mode(self):
        """A uniform 20% under-dose is a calibration constant, not a
        pattern error: absolute RMS sees it, calibrated RMS doesn't."""
        t = np.zeros((8, 8))
        t[2, 2] = t[5, 5] = 0.5
        counts = t * 1000 * 0.8
        rms, _, rms_cal = site_metrics(counts, t, 1000)
        assert abs(rms - 0.2) < 1e-9
        assert rms_cal < 1e-12

    def test_shot_noise_rms_scales_as_inverse_sqrt_dose(self):
        """Sanity: sampling from the *exact* target follows 1/√N."""
        rng = np.random.default_rng(2)
        t = np.zeros((8, 8))
        idx = [(1, 1), (1, 4), (4, 1), (4, 4)]
        for i, j in idx:
            t[i, j] = 0.25
        rmss = []
        for D in (1e3, 1e5):
            reps = []
            for _ in range(20):
                counts = rng.poisson(D * t)
                rms, _, _ = site_metrics(counts.astype(float), t, D)
                reps.append(rms)
            rmss.append(np.mean(reps))
        ratio = rmss[0] / rmss[1]
        assert 5 < ratio < 20, f"expected ~10 (=√100), got {ratio:.1f}"

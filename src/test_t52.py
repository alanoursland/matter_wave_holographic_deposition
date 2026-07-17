"""Tests for T52 joint-noise interaction metrics."""

import numpy as np
import pytest

from t52_joint_noise_process_gate import interaction_metrics


def test_additive_changes_have_zero_interaction():
    coherent = np.ones((2, 3, 3))
    dephasing = coherent + 0.1
    column = coherent - 0.2
    joint = coherent - 0.1
    result = interaction_metrics(coherent, dephasing, column, joint)
    assert result["relative_to_coherent"] == pytest.approx(0.0, abs=1e-15)
    assert np.allclose(result["residual"], 0.0)


def test_interaction_metrics_reject_shape_mismatch():
    with pytest.raises(ValueError, match="identical"):
        interaction_metrics(
            np.ones((2, 2)), np.ones((2, 2)),
            np.ones((2, 2)), np.ones((3, 3)),
        )

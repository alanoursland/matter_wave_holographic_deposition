"""Tests for the T49 correlated process-yield gate."""

import numpy as np

from t49_dephased_process_yield_gate import (
    classify_gate,
    paired_outcomes,
)


def _rows(values):
    return [{"array_pass": bool(value)} for value in values]


def _summary(yield_value, interval):
    return {
        "full_field_array_yield": yield_value,
        "full_field_array_yield_95ci": interval,
    }


def test_paired_outcomes_counts_directional_changes():
    result = paired_outcomes(
        _rows([1, 1, 1, 0, 0]),
        _rows([1, 0, 1, 1, 0]),
    )
    assert result["both_pass"] == 2
    assert result["coherent_only_pass"] == 1
    assert result["dephased_only_pass"] == 1
    assert result["both_fail"] == 1
    assert np.isclose(result["dephasing_loss_probability"], 0.2)


def test_gate_pass_requires_resolved_yield_and_loss_bounds():
    coherent = _summary(1.0, [0.97, 1.0])
    dephased = _summary(1.0, [0.97, 1.0])
    paired = {"dephasing_loss_probability_95ci": [0.0, 0.03]}
    assert classify_gate(coherent, dephased, paired) == "not_falsified"


def test_gate_is_inconclusive_when_baseline_is_not_resolved():
    coherent = _summary(0.96, [0.91, 0.99])
    dephased = _summary(0.90, [0.84, 0.95])
    paired = {"dephasing_loss_probability_95ci": [0.03, 0.10]}
    assert classify_gate(coherent, dephased, paired) == "inconclusive"

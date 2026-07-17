"""Tests for T51 column/process-yield helpers."""

from t51_column_process_yield_gate import classify_yield


def _summary(interval):
    return {"full_field_array_yield_95ci": interval}


def test_yield_classification_uses_resolved_confidence_interval():
    assert classify_yield(_summary([0.96, 1.0])) == "not_falsified"
    assert classify_yield(_summary([0.70, 0.90])) == "falsified"
    assert classify_yield(_summary([0.92, 0.98])) == "inconclusive"

"""Tests for T50 column-to-device gate helpers."""

from t50_column_device_function_gate import classify_corner


def test_corner_classification_distinguishes_recalibration():
    assert classify_corner({"functional": True}, [1.0, 2.0]) == "fixed_dose_pass"
    assert classify_corner({"functional": False}, [1.5, 2.0]) == "recalibration_required"
    assert classify_corner({"functional": False}, None) == "functional_failure"

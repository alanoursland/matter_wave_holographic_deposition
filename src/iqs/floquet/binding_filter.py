"""Reusable adapter for the historical sideband binding filter."""


def apply_binding_filter(pipeline, psi_dressed, average_populations, **kwargs):
    """Apply the pipeline's binding-resonance filter stage."""
    return pipeline.stage3_binding_filter(
        psi_dressed, average_populations, **kwargs)


__all__ = ["apply_binding_filter"]

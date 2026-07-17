"""Reusable adapter for the historical spatial-dressing model."""


def dress_spatial(pipeline, psi, phase_map, **kwargs):
    """Apply the pipeline's spatial Floquet dressing stage."""
    return pipeline.stage2_floquet_dress_spatial(psi, phase_map, **kwargs)


__all__ = ["dress_spatial"]

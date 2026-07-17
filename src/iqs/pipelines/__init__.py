"""Reusable end-to-end simulation pipelines."""

__all__ = [
    "HolographicCagingPipeline",
    "IntegratedPipelineV10",
    "IntegratedQuantumSubstrate",
    "PatternedSubstratePipeline",
]


def __getattr__(name):
    if name in {"IntegratedQuantumSubstrate", "PatternedSubstratePipeline"}:
        from .patterned_substrate import IntegratedQuantumSubstrate

        return IntegratedQuantumSubstrate
    if name in {"IntegratedPipelineV10", "HolographicCagingPipeline"}:
        from .holographic_caging import IntegratedPipelineV10

        return IntegratedPipelineV10
    raise AttributeError(name)

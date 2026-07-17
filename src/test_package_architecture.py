"""Package-oriented imports and historical compatibility aliases."""

from coherent_matterwave_beam import CoherentMatterwaveBeam as LegacyBeam
from inverse_holography import InverseHolographySolver as LegacySolver
from sim_v10 import IntegratedPipelineV10 as LegacyV10
from sim_v9 import IntegratedQuantumSubstrate as LegacyV9

from iqs.holography import InverseHolographySolver
from iqs.pipelines import HolographicCagingPipeline, PatternedSubstratePipeline
from iqs.sources import CoherentMatterwaveBeam


def test_legacy_names_are_package_implementations():
    assert LegacyBeam is CoherentMatterwaveBeam
    assert LegacySolver is InverseHolographySolver
    assert LegacyV9 is PatternedSubstratePipeline
    assert LegacyV10 is HolographicCagingPipeline


def test_implementations_live_under_iqs():
    assert CoherentMatterwaveBeam.__module__.startswith("iqs.sources")
    assert InverseHolographySolver.__module__.startswith("iqs.holography")
    assert PatternedSubstratePipeline.__module__.startswith("iqs.pipelines")
    assert HolographicCagingPipeline.__module__.startswith("iqs.pipelines")

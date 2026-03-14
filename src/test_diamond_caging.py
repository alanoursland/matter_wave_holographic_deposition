import pytest
import numpy as np
from diamond_caging import DiamondNetworkSimulator

def test_spectrum_ab_caging_pi_flux():
    """
    At flux Phi=pi, a 2D Diamond network exhibits true Aharonov-Bohm caging.
    All bands flatten, resulting in exactly 3 unique energy levels.
    """
    simulator = DiamondNetworkSimulator(Lx=4, Ly=4)
    ok, H, evals = simulator.validate_spectrum(phi=np.pi, J=1.0)
    unique_bands = np.unique(np.round(evals, 2))
    
    assert ok is True
    assert len(unique_bands) == 3, f"Expected exactly 3 flat bands at Phi=pi, got {len(unique_bands)}."

def test_spectrum_dispersive_zero_flux():
    """
    At flux Phi=0, the Diamond network is highly dispersive.
    """
    simulator = DiamondNetworkSimulator(Lx=4, Ly=4)
    ok, H, evals = simulator.validate_spectrum(phi=0.0, J=1.0, abort_on_fail=False)
    unique_bands = np.unique(np.round(evals, 2))
    
    assert ok is False # Expecting failure for the flat-band validation
    assert len(unique_bands) > 3, "Spectrum should be highly dispersive at Phi=0."

def test_hamiltonian_hermitian():
    simulator = DiamondNetworkSimulator(Lx=4, Ly=4)
    H = simulator.build_hamiltonian(phi=np.pi, J=1.0)
    assert np.allclose(H, H.conj().T, atol=1e-10)

def test_unit_cell_dimensions():
    # 2D Diamond network requires 5 sites per unit cell
    Lx, Ly = 4, 4
    simulator = DiamondNetworkSimulator(Lx=Lx, Ly=Ly)
    assert simulator.dim == 5 * Lx * Ly
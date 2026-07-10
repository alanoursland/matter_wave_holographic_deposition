# Literature Validation of the Gen-2 Architecture Claims

Executed 2026-07-09 via deep-research workflow (114 agents, 31 sources
fetched, 145 claims extracted, 25 adversarially verified 3-vote).
Target: the capability claims in
[notes/design/gen2_architecture.md](../notes/design/gen2_architecture.md).

**Confidence tiers:** ✅ *verified* = survived 3-vote adversarial
verification; 📄 *sourced* = extracted from a primary source with
verbatim quote, not adversarially verified; ❓ *unverified* = no source
found in this pass.

| # | Claim (gen-2 note) | Status | Finding |
|---|---|---|---|
| 1 | GFIS single-atom tip, ΔE ≈ 0.25–0.5 eV, 0.5 nm probes | ✅ verified, with nuance | Probes ≤0.5 nm routine (0.35 nm edge-resolution convention). Emitter is a **trimer apex with one beamlet aperture-selected** — functionally single-atom. ΔE: direct beam measurement is only an upper bound **<1.0 eV** (spectrometer-limited, Notte/ALIS 2006); the ~0.41 eV figure is inferred from field-ion data (Ernst et al., Appl. Surf. Sci. 67, 111 (1993)). Cite as "<1 eV measured, ~0.4 eV inferred". [Allen, Beilstein J. Nanotechnol. 12 (2021); Notte et al., NSTI-Nanotech 2006] |
| 2 | LoTIS/MOTIS: sub-eV ΔE, high brightness, commercial | ✅ verified, with correction | Measured: 0.45 eV ΔE (Cs⁺ at 100 kV/m extraction), brightness (2.4±0.1)×10⁷ A m⁻² sr⁻¹ eV⁻¹ at 10 keV (~20× Ga LMIS), smallest spot (2.1±0.2) nm; **sub-nm is projected, not demonstrated**; ≤0.1 eV exists only for pulsed low-charge Rb⁺ bunches (0.02 eV, Reijnders PRL 102, 034802 (2009)). Commercialized: zeroK NanoTech (NIST license, FIB:ZERO). [Steele et al., Nano Futures 1, 015005 (2017); Knuffman et al., JAP 114, 044303 (2013); McClelland et al., Appl. Phys. Rev. 3, 011302 (2016)] |
| 3 | Programmable electron phase plates, ~48 px | 📄 sourced | Confirmed: 48-element programmable electrostatic (Einzel) phase plate operated in a TEM (C2 plane, 100–300 kV), phase sensitivity 0.075 rad/mV at 300 kV, resolution ~3×10⁻³ π. Published SciPost Phys. 15, 223 (2023) — the pixel-count state of the art as of 2023. [Yu, Vega Ibáñez, Béché, Verbeeck; arXiv:2308.16304] |
| 4 | IPL: real program, demonstrated resolution/demag | ✅ verified, one gap | Two IMS Vienna system generations by 1998 (H⁺/H₂⁺/He⁺, 3–10× demag, 70–150 keV); **70 nm line/space demonstrated**, field distortion <0.15 μm over 8×8 mm; MEDEA PDT: 4× reduction He⁺ ~75–100 keV, 50 nm *predicted*, column built by Nov 2000. **Why it ended is NOT established by any verified source** — that part of the note stays uncited. [Melngailis et al., JVST B 16, 927 (1998); Kaesmaier & Löschner, Proc. SPIE 4226, 52 (2000)] |
| 5 | LEEM cathode-lens decel to eV landing; C_c practice | 📄 sourced, partial | LEEM/PEEM with aberration-corrected cathode lenses is established; C_c is a measured instrument parameter (Tromp's method validated by ray tracing + experiment, Ultramicroscopy 2019). The specific **C_c ∝ decel-gap scaling used in T18 remains an assumption** (parameterized there, so no conclusion changes). |
| 6 | Soft landing of mass-selected ions, 1–10 eV | 📄 sourced | Established technique: "hyperthermal (<100 eV)" mass-selected ion soft landing with control of m/z, energy, charge state; sub-nm cluster deposition demonstrated. Our 1–10 eV regime is the low-energy subset. [Johnson/Laskin-school review, Mass Spec. Rev. (mas.21451)] |
| 7 | Deterministic single-ion implantation via secondary-electron counting | 📄 sourced, **corrected** | Deterministic implantation is practiced: single 14 keV P⁺ into ²⁸Si counted at **99.87±0.02% confidence** (Jamieson group, Adv. Mater. 34, 2103235 (2022)) — but via **ion-beam-induced charge (IBIC) detection, not secondary-electron counting** (SED has only reached ~80% confidence, 25 keV Bi⁺). Placement: sub-100 nm routine; ~15 nm (Kane-architecture) **not yet demonstrated**; without counting, Poisson caps exactly-one-ion at 36%. Gen-2 note §3 corrected. [arXiv:2009.02892; roadmap arXiv:2412.05729] |
| 8 | H:Si depassivation lithography: atomically precise, RT-stable | 📄 sourced | Confirmed: HDL is "the only technique capable of patterning a hydrogen resist on Si(100) with atomic precision"; AP mode (<4.5 V) gives 0.77 nm lines with atomically sharp edges; patterns RT-stable because Si–H is covalent. No explicit eV barrier quoted in sources — T20's E_a = 2.5 eV stays a representative parameter. [HDL roadmap, arXiv:2501.04535] |
| 9 | TEM off-axis holography phase stability | 📄 sourced — **stronger than claimed** | Drift-corrected hologram-series best practice: phase error **0.006 rad (2π/1050) over 900 s** exposure (0.0037 rad in vacuum), with independent image- and phase-drift registration as part of the method. Directly supports T19's framing that static/drift errors are handled by calibration loops; the demonstrated stability exceeds the "TEM-grade" level T19 assumed. [Ultramicroscopy (2014), S0304399114000461] |
| 10 | Barth–Kruit chromatic FW50 coefficient 0.34 | ✅ verified | Correct under the **FWHM energy-spread convention** (Barth & Kruit, Optik 101, 101 (1996)); the coefficient is 0.6 if ΔE is quoted as FW50 (Hagen & Kruit, JVST B 27, 2654 (2009)). Independent numerical re-derivation by the verifier gave 0.344 — matching T18's wave-optical gate (0.34 within 12%). |
| 11 | Patch potentials: mV-scale, **μm-scale correlation** | 📄 sourced, **corrected** | Magnitude: ~30 meV FWHM work-function distributions on polycrystalline Au; adsorbate-dominated (submonolayer C shifts φ by up to 0.8 eV); cleaning buys ~2 orders of magnitude in field noise. **Correlation lengths are 50–100 nm, not μm** — patch fields are *rougher* than the note assumed, which matters for common-mode suppression assumptions in T19 (finer patches decay faster with distance — direction favorable — but average less across the beam). Gen-2 note §2.2 annotated. [NJP 23, 103028 (2021)] |
| 12 | Langevin k_L ~2×10⁻¹⁵ m³/s; Bohdansky threshold; Si displacement 15–25 eV | 📄 sourced, **partly corrected** | Langevin: standard framework, order-of-magnitude OK (deviations documented for polarizable ions). Bohdansky: formula real; correct citation **Bohdansky, Roth & Bay, JAP 51, 2861 (1980)** + erratum JAP 52, 1610 (1981). Si displacement (DFT-MD): minimum **12.5±1.5 eV** (⟨111⟩, matches experiment 12.9), 20 eV ⟨100⟩, **lattice average 36±2 eV** — T20's 15 eV surface parameter sits between the directional minimum and ⟨100⟩, i.e. conservative for safety margins; cite direction-resolved values. |
| 13 | Multi-beam parallelism 10³–10⁵ | 📄 sourced — confirmed and exceeded | IMS MBMW: **262,144 beams (512×512)** commercial since 2014 (mask writing; 7 nm-node use by 2016). MAPPER: 110-beam demonstrator (2007), 13,000-beam design, bankrupt 2018 (ASML bought IP, discontinued); stated failure modes: beam-beam Coulomb control, funding, EUV pivot. Note: single-ion columns dodge the first failure mode identically to how they dodge in-column Coulomb blur. |
| 14 | Cr standing-wave atom lithography ~40 nm | 📄 sourced, **corrected number** | McClelland et al., Science 262, 877 (1993): Cr focused by standing wave, lines **65±6 nm wide** at **212.78 nm pitch (= λ/2 of 425.55 nm)**. The ~40 nm figure belongs to later refinements; cite 65 nm for the 1993 demonstration. |
| 15 | Shimizu/Fujita Ne* atom holography 1996 | 📄 sourced, partial | Confirmed: Fujita, Morinaga, Kishimoto, Yasuda, Matsui & Shimizu, "Manipulation of an atomic beam by a computer-generated hologram", Nature 380, 691 (1996) — Ne* through a binary SiN hologram reconstructing a pattern. The **electrically switchable** variant was not confirmed in this pass (❓ — likely the 2000-era follow-up; needs its own citation). |
| 16 | AB caging observed in synthetic lattices | 📄 sourced — confirmed twice | Photonic rhombic lattice: AB cages observed, all bands flat at π flux (Mukherjee et al., PRL 121, 075502 (2018)). Superconducting circuit: single microwave photon caged, excluded-site population <5% over >9 swap times (Sci. Adv. 9, eadj7195 (2023)). **Caveat found: interactions break caging** (bound photon pairs delocalize) — reinforces T15/M6's single-particle framing. |
| 17 | US 9,502,202 B2 exists; no experimental support | 📄 sourced — confirmed | Patent exists: "Systems and methods for generating coherent matterwave beams", Lockheed Martin, Arman & Chase, filed 2011, granted 2016, active to 2032. **Purely theoretical** — Kuramoto equations, ~4 ns sync time, no construction/test/data; citations cover AB effect and BEC physics generally, none validating the synchronization mechanism. Matches fable5 M2 and the T12 demotion. |

## Corrections applied to the gen-2 note

1. §3/§6: single-ion detection is **IBIC** (99.87%), not secondary-electron counting (~80%) — claim 7.
2. §2.2 context: patch correlation lengths are **50–100 nm**, not μm — claim 11.
3. §7 table: Cr atom lithography 1993 = 65 nm lines at λ/2 pitch — claim 14.
4. IPL end-reason left uncited (no verified source) — claim 4.
5. T20's Si displacement parameter annotated with direction-resolved DFT values — claim 12.

## Open follow-ups

- Why IPL actually ended (needs a historical source).
- Direct sub-1-eV spectrometer measurement of the GFIS beam ΔE.
- Sub-nm LoTIS spot since 2017; shipping FIB:ZERO ΔE at customer sites.
- Citation for the electrically switchable atom hologram (claim 15).
- LEEM C_c-vs-gap scaling law (T18 currently parameterizes it).

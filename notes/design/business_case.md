# Business Capability: Where a Matter-Wave Printer Wins and Loses

*Companion to [gen2_architecture.md](gen2_architecture.md). Written
2026-07-10 from the fable5-corrected simulator numbers (T10/T11/T21/T22)
plus market research (web sources cited inline; prices are a July-2026
snapshot taken mid-supercycle — treat as volatile).*

---

## 1. The honest comparison against TSMC

Three metrics, three different answers. Assumptions: gen-2 machine =
1,000 microcolumns × 10⁹ ions/s (T10 detector-limited), ~$5M
instrument-class capex; full IC layer needs ~10¹⁴ atoms/cm²; TSMC
gigafab = ~100k wafer starts/month, ~$20B, 3nm wafer ≈ $19.5k, mask set
≈ $15M, cycle time 3–4 months, 1,000–2,000 process steps.

| Metric | Result |
|---|---|
| Raw throughput, one machine vs one gigafab | **~10⁵× slower** (~50 days/wafer-equivalent if the beam delivers every atom; one machine ≈ 7 wafers/yr vs 1.2M/yr) |
| Throughput per capex dollar | **~40× slower** direct-print; **~parity to 10×** with template-seeding (§3) |
| Time & cost to first article of a new design | **~10⁴× faster, ~10⁶× cheaper** (hours + ~$0 marginal vs ~6 months + ~$20M mask set/queue) |

The structural reasons, from the fab-economics research:

- A fab is not a lithography building: litho is ~25% of capex; a chip
  is ~15 metal layers of interconnect over one transistor layer, and
  each layer is a pattern→etch→deposit→polish loop. The fab's moats
  are yield learning, throughput at volume, and process integration —
  none of which an atom-placing machine can attack.
- One EUV scanner projects ~10¹² features/s; ion counting tops out at
  10⁹–10¹⁰ features/s per 1,000-column machine. This gap is physics
  (photons are projected in parallel; atoms are counted), permanent,
  and defines the market boundary.
- What the machine deletes is the *other* side of fab economics: the
  $15M mask set, the ~10³-unit minimum economic volume, the months of
  cycle time, and most of the per-layer step count (additive
  direct-write collapses the 7-step subtractive patterning loop).

**One-line positioning: never a factory; a compiler.** TSMC turns
designs into millions of chips; this machine turns designs into first
articles the same day.

## 2. Market wedges, in order of fit

1. **Ultra-dense write-once archival storage** — the best-matched
   product found. Periodic bit arrays are the machine's easiest
   pattern (band-limit friendly; the dots target hits SSIM 0.897 at
   the T22 stage); dose is nearly free (T11: ~10² ions to fidelity
   ceiling); retention is Arrhenius-forever at covalent site barriers
   (T20: 630 K for E_a = 2.5 eV). Numbers: at template-quantized
   ~2.7 nm cells → ~14 Tb/cm² *per layer* (≈5× Samsung's 400-layer
   V-NAND at 28 Gb/mm²); single-site bits → ~250×. Write rate
   1–10 GB/s per 1,000-column machine (10–100 ions/bit) — same class
   as Cerabyte's demonstrated 1 GB/s laser writer and above LTO drives
   (~0.4 GB/s). Cost floor to beat: tape at ~$3–8/TB media — a
   ~350 TB/layer platter clears it if media stays under ~$1.7k.
   Readout (ion/e-beam scanner or electrically readable dopant
   crossbars) is the open product-design problem. Market signal: DNI /
   National Academies cold-storage report; competitors-in-spirit:
   Cerabyte, Project Silica, DNA storage — none atomic-density.
2. **Low-volume / high-mix ICs** (the Minimal Fab market, finer): the
   proven $5.8M half-inch maskless "fab-in-a-box" (AIST/Yokogawa, 140
   companies, shipping for HMLV sensors) validates demand for
   mask-free small-batch silicon; this machine is that paradigm
   extended ~4 orders in precision. The 2021 shortage was concentrated
   in exactly this segment (older-node, low-volume parts).
3. **Quantum/dopant devices**: deterministic donor arrays are already
   practiced one-ion-at-a-time (99.87% IBIC counting); this machine is
   the programmable-pattern generalization. Small market, perfect
   technical fit, earliest credible customers.
4. **Mask writing / R&D patterning**: competes with e-beam where
   atomic resolution or charged-particle advantages matter; IMS
   multi-beam (262k beams) owns the volume mask market.
5. **Emerging-NVM crossbar layers** (ReRAM/PCM startups without EUV
   access): periodic line/space at 10–15 nm pitch matches the T22
   feature scale over stitched fields.

Explicitly out of reach: DRAM (density stalled on capacitors, not
litho — nothing to sell), commodity NAND (3D stacking already matches
a naive 5 nm-cell printed monolayer at 2.8 Tb/cm²; only atomic-pitch
bits beat it), any volume logic.

## 3. The lever that changes the throughput math

Don't print the material — **print the template**. Beam writes seed
patterns at 10¹–10² ions/feature (T21 dose); selective chemistry
(ALD on depassivated/templated sites, as the H:Si community practices)
grows the bulk. Cuts beam dose per layer 10²–10³ → wafer-equivalent of
critical layers in hours instead of ~50 days, moving the
capex-normalized comparison from 40× slower to parity-to-10×. This is
also exactly the gen-2 §3.4 division of labor (beam chooses sites,
chemistry supplies precision), so it needs no new physics — it is the
architecture.

## 4. Price environment (July 2026 snapshot)

AI supercycle: DRAM up 300–600% from 2024 lows (~$1,200/TB), NAND
sold out through 2026 (SSD ~$90–100/TB), HDD ~$30–45/TB, LTO tape
~$3–8/TB effective. Implications: friendliest fundraising climate for
storage newcomers in a decade; but the machine's development runway
(~decade) outlives any cycle — the durable planning number is the tape
floor (~$5/TB, halving per LTO generation every ~3 years).

## 5. Honest risk register (business side)

- **The graveyard is economic, not physical**: IPL demonstrated 70 nm
  and died; MAPPER built 13k beams and went bankrupt (EUV pivot,
  funding, beamlet Coulomb control); standing-wave atom lithography
  worked and stalled. Two of those three failure modes apply here
  (funding, EUV gravity); the third (multi-beam Coulomb blur) is the
  one this machine is uniquely exempt from (single-ion operation,
  T10/T21).
- **Nobody else is working on this combination** (2026 literature
  sweep, [fable5/literature_validation.md](../../fable5/literature_validation.md)) —
  unclaimed territory, and also a warning that no one has yet found it
  worth trying. The physical reason ion holography barely exists (the
  √(m/mₑ) ≈ 85× phase-stability penalty vs electrons) is the same T19
  constraint this architecture answers with microcolumns.
- **The one experiment that gates everything**: differential
  potential-drift measurement across a μm-scale cold microcolumn gap
  (T19 spec: nV-class between recalibrations). Tabletop-scale,
  decisive either way, publishable either way.
- Readout for archival; multi-material chemistry for ICs; die-scale
  field stitching (T18 caveat) — product-engineering unknowns, not
  physics unknowns.

## 6. Sources

Fab economics: PatentPC & Silicon Analysts cost guides, BCG, AnySilicon
(wafer $19.5k, masks $15M, fab $15–20B, High-NA EUV $380M, 1–2k steps).
Memory state: TechInsights/SemiAnalysis (DRAM stall), Tom's
Hardware/TrendForce (400-layer V-NAND, 28 Gb/mm²), Forbes/Coughlin
(emerging NVM). Archival: Fujifilm LTO-10, WD HAMR, The Register,
Blocks & Files, StorageNewsletter (Cerabyte 1 GB/s), DCD (DNI report).
Prices: TrendForce, Tom's Hardware RAM index, BuyPerUnit, How-To Geek
(LTO drive economics). Minimal Fab: Yokogawa, SemiMedia, IEEE 7101289.
Full URLs in the session research logs; capability claims
fact-checked in [fable5/literature_validation.md](../../fable5/literature_validation.md).

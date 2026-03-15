Use STAR Cryoelectronics (starcryo.com) — they explicitly offer custom Nb/Al-AlOx/Nb trilayer SQUID arrays on 150 mm wafers. They’ve been doing this for decades. A 5×5 or 7×7 array is basically a single die; you can get multiple copies for a few thousand dollars.

Alternative: HYPRES (now part of a larger superconducting IC foundry) or university cleanrooms (MIT Lincoln Lab, NIST Boulder, or any place with Nb trilayer process).

Process: standard Nb trilayer with your diamond-lattice layout (loops + bias lines). Include the apodization taper at the edges exactly as in v9.

Drive & chamber (off-the-shelf)

GHz AC bias (1–10 GHz signal generator + coax feedthrough) + DC flux bias lines.
UHV chamber (already common in atom optics labs) with a simple thermal Nb or Pb evaporator (or pulsed laser ablation for cleaner deposition).
Detection: after deposition, warm up slightly and image with AFM/SEM or even optical microscope (at micron scale you’ll see the pattern by eye).
# T34 - Sourced p+ Si(100) / Al-1.5%Si electrical coupon

Implemented 2026-07-11.

## Purpose

T33 established a two-layer geometry but used exploratory sticking, diffusion,
and escape assumptions and stopped at visible opens and metal bridges. T34
replaces that materials benchmark with one prepared interface for which bulk
silicon resistivity, Al/Si contact resistivity, and clean-surface adatom
energetics have primary-source support. It then extracts contact resistance
and leakage from each realized contact area.

Study: `src/t34_si_al_contact_extraction.py`.

Reusable extractor: `src/iqs/deposition/electrical.py`.

Artifacts:

- `results/t34_si_al_contact.json`;
- `results/t34_si_al_contact.npz`;
- `results/t34_si_al_contact.png`.

## Selected interface

The benchmark follows the preparation in Carver et al.:

- boron-doped p+ Si(100);
- 200 nm dry thermal SiO2;
- buffered-HF contact windows;
- vacuum-sputtered Al-1.5%Si, 700 nm thick;
- 450 C sinter in 10% forming gas for 20 minutes.

The selected NBS boron-doped material has `NA = 6.655e19 cm^-3` and measured
resistivity `1.707e-3 ohm cm` at 23 C. Digitizing the corresponding point in
Carver et al. Fig. 4 gives `rho_c = 2.5e-7 ohm cm^2`; the plotted uncertainty
is approximately `1.5e-7` to `4e-7 ohm cm^2`.

The contact measurement used approximately 4.8 micrometre windows. Applying
its area-normalized interface resistivity to a 55 nm contact is an explicit
scale extrapolation, not a nanoscale measurement.

Primary sources:

- Carver et al., *IEEE Transactions on Electron Devices* 35, 489-497 (1988),
  DOI `10.1109/16.2483`;
- Thurber et al., NBS Special Publication 400-64 (1981);
- Brocks, Kelly, and Car, *Surface Science* 269/270, 860-866 (1992), DOI
  `10.1016/0039-6028(92)91362-F`.

## Replaced assumptions

The earlier common `1.2 eV` escape barrier is removed from this benchmark.
For clean reconstructed Si(100), Brocks et al. report:

```text
             binding energy    diffusion barriers parallel/perpendicular
Si adatom        4.6 eV                         0.6 / 1.0 eV
Al adatom        3.6 eV                         0.3 / 0.1 eV
```

These single-adatom values are stored as provenance but are not converted to
a Gaussian film diffusion length. At room temperature isolated Al is highly
mobile and rapidly encounters other atoms; a hold-time random walk without a
nucleation/capture model would be physically false. Likewise, the arbitrary
`0.05` off-Si sticking coefficient is gone. The prepared oxide and patterned
contact windows define where metal is electrically active.

The source process uses Al-1.5%Si rather than pure Al. A holographic printer
would therefore need alloy feedstock or about 1.5% Si co-deposition to claim
the sourced electrical interface.

## Extraction

For each realized cap, pixels above half the nominal film thickness define
the effective contact area `A`. The interface term is

```text
R_interface = rho_c / A.
```

The extractor adds circular equal-area half-space spreading resistance. For
every contact pair on one continuous substrate, it also includes the
first-order mutual spreading term and reports

```text
I_pair = V_test / (R_contact,1 + R_substrate,pair + R_contact,2).
```

This exposes conduction through the semiconductor even when no metal bridge
exists. It is an ohmic DC model, not a Schottky, dielectric-tunnelling, grain,
or quantum-contact model.

## Result

Thirty realizations use the existing 11 nm FWHM projection response and 1 nm
registration sigma on a 3x3 array with 55 nm contacts at 100 nm pitch.

```text
contact yield below 10 kOhm              100%
median contact + spreading resistance    8.40 kOhm
worst contact + spreading resistance     8.60 kOhm
isolation yield below 1 nA at 1 V        0%
median worst pair current                 59.7 uA
worst pair current                        61.1 uA
```

The contact result is plausible under the sourced process, subject to the
55 nm scale extrapolation. The array is not an array of independent devices:
all contacts share a highly conductive p+ silicon half-space. Its roughly
60 microamp pair current exceeds the provisional 1 nA isolation criterion by
about 60,000 times.

## Impact

The geometry is no longer allowed to call a coupon electrically successful
merely because its metal islands are disconnected. The selected p+ bulk
substrate is suitable for contact-resistance metrology but unsuitable for an
independently addressable printed array.

The next experiment should keep this Al-1.5%Si contact stack and move the p+
regions onto a specified SOI device layer with etched isolation to the buried
oxide. The same extractor can then compare contact resistance against BOX,
sidewall, and surface leakage instead of the unavoidable shared-bulk path.

T35 performs that substrate change and converts trench-surface cleanliness
into a quantitative sheet-resistance qualification.

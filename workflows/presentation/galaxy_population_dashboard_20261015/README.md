# Galaxy population explorer

Private Pages dashboard: https://chatgpt.com/space/page_aec9d1a1537c81918114b8682a77d597

`build_data.py` constructs the binned source aggregates. `engine.js` applies
interactive selections and population mixtures. `dashboard.html` is the UI
template; `validate.js` checks the generated payload and numerical engine.
Generated data, standalone HTML, source audit and validation receipt are in
`~/fmsu_tokyo_20261010/galaxy_population_dashboard_20261015/`.

Eligible parent counts: 191,944 low-z Loa GAT, 4,988,309 Loa P12-A,
66,691 held-out TNG300 z=0. Properties: quenched fraction, log stellar mass,
log SFR, log sSFR, observed g-r, log age/Gyr, attenuation A_V and AGN fraction.
The last three lack comparable TNG inputs and are explained as unavailable.
Loa coverage has a gap at 0.06 <= z < 0.15. No interpolation fills it.

Controls use aggregate bin boundaries; cell means are quantized to 0.001.
Mixtures require support for all positively weighted components. Kish effective
counts retain sums of squared weights. Native classes and within-source trace
quartiles are descriptive comparisons, not a common calibrated environment
estimand. Old Loa has no retained eigenvalues for trace ranks. TNG brightness
uses an approximate griz continuum through DECam responses at hypothetical
observing redshift, with no population evolution or dust; extrapolation at
higher observing redshift is flagged. No causal or uncertainty claim is made.

Validation: 106 checks passed, including decode parity, source counts, redshift
gap, selections, mixtures and JavaScript syntax. The Page and embed were read
back after saving. Browser layout was not visually verified because the local
Playwright installation lacks a browser executable.

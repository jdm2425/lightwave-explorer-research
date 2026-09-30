# Linear dispersion diagnostic

Source parameters are copied from the GDD=0 fixed-period PPLN input.

Two simulations are run:
- LINEAR_INPUT: init only, to save the exact input field
- LINEAR_500um: 500 um linear propagation in MgO:LiNbO3

The nonlinear polarization is disabled by using LWE's linear(...) sequence
function. The diagnostic therefore isolates material dispersion, GVM and linear
diffraction from nonlinear conversion, depletion and back-conversion.

Fine grid:
- dt = 0.45 fs
- dx = 6 um
- dz = 0.25 um

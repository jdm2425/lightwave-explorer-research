# Isolated 22→24 µm PPLN half-crystal, zero pump GDD

This is an independent 250 µm Lightwave Explorer simulation. The input comes
from the archived exact-domain 20→24 µm, 500 µm reference, retaining its
100 nJ pump, 10 nJ coherent MIR seed, zero pump GDD, zero relative input delay,
200 µm waists and numerical grid. The poling period increases linearly from
22 to 24 µm across 250 µm. The MIR seed is added once at the entrance.

`make_input.py` validates that its domain algorithm reproduces the archived
20→24 µm reference to the printed precision. It then generates 22 domains
(including the partial final domain), sets the new crystal thickness, and
changes the output base name.

From the unpacked folder on your Mac, run:

```bash
caffeinate -dimsu bash run.sh
```

After completion, the terminal displays the 950–1000 nm, 1000–1050 nm, and
direct 950–1050 nm integrated energies in pJ. It also saves a CSV in `results/`
and `SEEDTEST_FWD_22_24_250_SINGLE_results.zip`. Upload **that results ZIP**
here if you would like me to compare the output with the other half crystals.

This freshly seeded 22→24 µm piece is distinct from the second half of an
already-propagating 20→24 µm crystal, whose input fields have evolved through
its first half.

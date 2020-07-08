# apply_desifootprint_nz

## multi_box2desi_cutsky_nz.py

This script applies desi footprint and n(z) on the output from Yuuki's full sky light cones.

The script consumes a configuration file very similar to the one used by Yuuki's light cone code.
To run it, just

`
python multi_box2desi_cutsky_nz.py config.ini --dir_out [directory holding input lightcone shells and output]
`

See the file test.config for configuration options. I've run it within desi environment version 19.12.


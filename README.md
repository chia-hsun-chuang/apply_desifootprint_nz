# Requirements

* [Python](http://www.python.org) 3.5 or 3.6

* [Astropy](http://www.astropy.org)

* [configparser](https://docs.python.org/3/library/configparser.html)

* [CAMB](https://camb.readthedocs.io/en/latest/)

* [desimodel](https://desimodel.readthedocs.io/en/latest/)

* [schwimmbad](https://schwimmbad.readthedocs.io/en/latest/)

# multi_box2desi_cutsky_nz.py

This script applies desi footprint and number density evolution evolution on the output from Yuuki's full sky light cones.

It consumes a configuration file very similar to the one used by Yuuki's light cone code:

```
[dir]
dir_out         Output directory
dir_gcat        Directory where the input out_*p.list.gcat files are located
name_template   Template string describing name of snapshot (e.g., out_{}p.list-wpmax-v3.gcat, where {} will be filled with the snapshot number)
file_alist      File containing list of simulation scale factors
file_camb       CAMB configuration file

[sim]
boxL            Length of the simulation box [Mpc/h]
shellwidth      Width of shell [Mpc/h]
zmin            Minimum of redshift range to use
zmax            Maximum of redshift range
shellnums       Comma-separated list of numbers (if provided, ignores the redshift range specified above)

[randoms]
n_multiplier    Relative size of randoms to samples (integer)

[elg]
nz_mu, nz_sigma, nz_shift, nz_amp   Parameters of shifted log-normal n(z) model
sample_index                        Index number used on Shadab's MTHOD code for ELG (check header of gcat files)

[lrg]
nz_mu, nz_sigma, nz_shift, nz_amp   Parameters of shifted log-normal n(z) model
sample_index                        Index number used on Shadab's MTHOD code for LRG (check header of gcat files)

[qso]
nz_a, nz_b                          Parameters of linear n(z) model
sample_index                        Index number used on Shadab's MTHOD code for QSO (check header of gcat files)
```

The file job_postprocess.sh has an example job file, to be submitted via sbatch.



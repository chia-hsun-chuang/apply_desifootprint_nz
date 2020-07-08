#!/usr/bin/env python

# This script process the lightcone shells created by Yuuki's lightcone code
# downsampling according to given n(z), applying DESI footprint and generating randoms
#
import desimodel.footprint as foot
import desimodel.io
import warnings
import numpy as np
import configparser
import argparse
import os
import sys
import camb
import datetime
from desimodel.io import fits
import glob
import os


def nz_model(z, nz_pars):
    """ model n(z) """
    if nz_pars['galtype'] in ['elg', 'lrg']:
        amp       = nz_pars["amp"]
        shift     = nz_pars["shift"]
        mu        = nz_pars["mu"]
        sigma     = nz_pars["sigma"]

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            nz = amp*np.exp(-(np.log(z-shift)-mu)**2/(2*sigma**2))/((z-shift)*sigma*np.sqrt(2*np.pi))
        nz[np.isnan(nz)] = -1

    elif nz_pars['galtype'] in ['qso']:
        a = nz_pars['a']
        b = nz_pars['b']
        nz = a*z + b

    else:
        raise RuntimeError("Unknown galaxy type.")

    return nz


def downsample(ra, dec, zz, zz_rsd, n_mean, nz_pars):
    """ downsample galaxies following n(z) model specified in nz_pars """

    nz          = nz_model(zz, nz_pars)

    # downsample
    ran         = np.random.rand(len(ra))
    nz_selected = (ran<nz/n_mean)
    idx         = np.where(nz_selected)

    print("Selected {} out of {} galaxies.".format(len(idx[0]), len(ra)))

    return ra[idx], dec[idx], zz[idx], zz_rsd[idx], nz[idx]


def apply_footprint(ra, dec, zz, zz_rsd, nzz):
    """ apply desi footprint """
    tiles = desimodel.io.load_tiles()
    point = foot.is_point_in_desi(tiles, ra, dec)
    idx   = np.where(point)

    print("Selected {} out of {} galaxies.".format(len(idx[0]), len(ra)))

    return ra[idx], dec[idx], zz[idx], zz_rsd[idx], nzz[idx]


### Loading config
parser = argparse.ArgumentParser()
parser.add_argument("config", help="ini file holding configuration",
                    type=str)
parser.add_argument("--dir_out", type=str, help="output directory (overrides config file)")
args = parser.parse_args()

configfile     = str(args.config) # config file

config         = configparser.ConfigParser()
config.read(configfile)

file_alist     = config.get('dir','file_alist')
file_camb      = config.get('dir','file_camb')
dir_out        = args.dir_out
if dir_out is None:
    dir_out    = config.get('dir','dir_out')
boxL           = config.getint('sim', 'boxL')
shellwidth     = config.getint('sim','shellwidth')
zmin           = config.getfloat('sim', 'zmin')
zmax           = config.getfloat('sim', 'zmax')

nrandoms  = int(config.get('randoms', 'n_multiplier'))

nz_pars = []
for i, galtype in enumerate(["qso", "lrg", "elg"]):
    nz_par = dict()
    if galtype in ["lrg", "elg"]:
        nz_par["galtype"]       = galtype
        nz_par["mu"]            = config.getfloat(f'{galtype}', 'nz_mu')
        nz_par["sigma"]         = config.getfloat(f'{galtype}', 'nz_sigma')
        nz_par["shift"]         = config.getfloat(f'{galtype}', 'nz_shift')
        nz_par["amp"]           = config.getfloat(f'{galtype}', 'nz_amp')
    else:
        nz_par["galtype"] = "qso"
        nz_par["a"]       = config.getfloat(f'{galtype}', 'nz_a')
        nz_par["b"]       = config.getfloat(f'{galtype}', 'nz_b')

    nz_par["galtype_index"] = config.getint(f'{galtype}', 'sample_index')
    nz_par["zmin"]          = config.getfloat(f'{galtype}', 'zmin', fallback=zmin)
    nz_par["zmax"]          = config.getfloat(f'{galtype}', 'zmax', fallback=zmax)

    nz_pars.append(nz_par)


### Running camb to get comoving distances

# Load all parameters from camb file 
pars = camb.read_ini(file_camb)
h    = pars.h
pars.set_for_lmax(2000, lens_potential_accuracy=3)
pars.set_matter_power(redshifts=[0.], kmax=200.0)
pars.NonLinearModel.set_params(halofit_version='takahashi')
camb.set_feedback_level(level=100)
results   = camb.get_results(pars)


### Loop over galtypes and shells for galaxies and randoms
for nz_par in nz_pars:
    galtype = nz_par["galtype_index"]
    shellnum_min = int(results.comoving_radial_distance(nz_par["zmin"])*h // shellwidth)
    shellnum_max = int(results.comoving_radial_distance(nz_par["zmax"])*h // shellwidth + 1)

    for i in range(shellnum_min, shellnum_max):

        CAT_file = dir_out + "lightcone_multibox_galtype{}_{}.fits".format(galtype, i)
        try:
            d    = fits.open(CAT_file)
        except IOError:
            print(f"WARNING: Couldn't open catalog {CAT_file} for galtype {galtype}, shellnum {i}",
                  file=sys.stderr)
            continue

        # compute the number density of the given snapshot by total number of galaxies / volume
        n_mean    = d[0].header['NGALBOX']/(1.* boxL**3)

        ra0       = d[1].data['RA']
        dec0      = d[1].data['DEC']
        zz0       = d[1].data['Z']
        zz_rsd0   = d[1].data['Z']+d[1].data['DZ']

        print(f"Downsampling file {i} for galtype {galtype}")
        ra, dec, zz, zz_rsd, nzz = downsample(ra0, dec0, zz0, zz_rsd0,
                                              n_mean, nz_par)

        print(f"Applying footprint for file {i} of galtype {galtype}")
        ra_desi, dec_desi, zz_desi, zz_rsd_desi, nzz_desi = apply_footprint(ra, dec, zz, zz_rsd, nzz)

        print('Writing footprinted catalog')
        OUT_file = dir_out + "lightcone_multibox_galtype{}_{}_footprint_nz.dat".format(galtype, i)
        np.savetxt(OUT_file, np.c_[ra_desi,dec_desi,zz_desi,zz_rsd_desi,nzz_desi])
        print(datetime.datetime.now().time())

    # merge the output files above    
    print('Merging shells into single file')
    out_f_name = dir_out + "UNIT_lightcone_multibox_{}_footprint_nz.dat".format(nz_par["galtype"].upper())
    for f in glob.glob(dir_out + "lightcone_multibox_galtype{}_*_footprint_nz.dat".format(galtype)):
        os.system("cat "+f+" >> " + out_f_name)
        os.system("rm "+f)
        print('Done!')

    # make randoms of nrandoms x data #
    for j in range(1,nrandoms):
        for i in range(shellnum_min, shellnum_max):
            CAT_file = dir_out + "lightcone_multibox_galtype{}_{}.fits".format(galtype, i)
            try:
                d    = fits.open(CAT_file)
            except IOError:
                print(f"WARNING: Couldn't open catalog {CAT_file} for nrandom {j}, galtype {galtype}, shellnum {i}",
                      file=sys.stderr)
                continue

            n_mean    = d[0].header['NGALBOX']/(1.*boxL**3)
            length    = d[1].header['NAXIS2']

            ra0       = np.random.rand(length)*360
            dec0      = np.arccos(2*np.random.rand(length)-1)*90/np.pi*2-90
            zz0       = d[1].data['Z']
            zz_rsd0   = d[1].data['Z'] + d[1].data['DZ']

            print(f"Downsampling randoms file {j} from {i} shell for galtype {galtype}")
            ra, dec, zz, zz_rsd, nzz = downsample(ra0, dec0, zz0, zz_rsd0,
                                                  n_mean, nz_par)

            print(f"Applying footprint for randoms file {j} from {i} shell of galtype {galtype}")
            ra_desi, dec_desi, zz_desi, zz_rsd_desi, nzz_desi = apply_footprint(ra, dec, zz, zz_rsd, nzz)

            print('Writing footprinted catalog')
            OUT_file = dir_out + "randoms_{}_lightcone_multibox_galtype{}_{}_footprint_nz.dat".format(j, galtype, i)
            np.savetxt(OUT_file,np.c_[ra_desi,dec_desi,zz_desi,zz_rsd_desi, nzz_desi])
            print(datetime.datetime.now().time())

    print('Merging random shells into single file')
    out_nx = dir_out+"randoms_{}x_lightcone_multibox_{}_footprint_nz_{}.dat".format(nrandoms,
                                                                                    nz_par["galtype"].upper(), j)
    for j in range(1,nrandoms):
        out_1x = dir_out+"randoms_1x_lightcone_multibox_{}_footprint_nz_{}.dat".format(nz_par["galtype"].upper(), j)
        for f in glob.glob(dir_out+"randoms_{}_lightcone_multibox_galtype{}_*_footprint_nz.dat".format(j, galtype)):
            os.system("cat "+f+" >> " + out_1x)
            os.system("rm "+f)

    for f in glob.glob(dir_out + "randoms_1x_lightcone_multibox{}_footprint_nz_*.dat".format(nz_par["galtype"].upper())):
        os.system("cat "+f+" >> " + out_nx)
    print('Done!')



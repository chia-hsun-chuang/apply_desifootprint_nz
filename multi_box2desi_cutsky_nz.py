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
from astropy.table import Table, vstack
from schwimmbad.mpi import MPIPool
from itertools import product
import glob
import os
import time


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


def write_catalog(fname, ra, dec, z, z_rsd):
    """ writes catalog in format to be ingested in E2E pipeline. """
    c1 = fits.Column(name='RA'     , array=ra       , format='E')
    c2 = fits.Column(name='DEC'    , array=dec      , format='E')
    c3 = fits.Column(name='Z'      , array=z_rsd    , format='E')
    c4 = fits.Column(name='Z_COSMO', array=z        , format='E')
    c5 = fits.Column(name='DZ_RSD' , array=z_rsd - z, format='E')

    hdu             = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5])
    hdr             = fits.Header()
    primary_hdu     = fits.PrimaryHDU(header=hdr)
    hdul            = fits.HDUList([primary_hdu, hdu])

    hdul.writeto(fname, overwrite=True)


def process_shell(args):
    """ Downsample galaxies according to nz, apply footprint, generate randoms and write to disk """
    i, nz_par, nrandoms, dir_out = args
    galtype = nz_par["galtype_index"]


    CAT_file = f"{dir_out}lightcone_multibox_galtype{galtype}_{i}.fits"
    OUT_file = f"{dir_out}cutsky/shells/lightcone_multibox_galtype{galtype}_{i}_footprint_nz.fits"

    if not os.path.exists(OUT_file):
        begin = time.perf_counter()
        try:
            d    = fits.open(CAT_file)
        except IOError:
            print(f"WARNING: Couldn't open catalog {CAT_file} for galtype {galtype}, shellnum {i}",
                  file=sys.stderr)
            return

        # compute the number density of the given snapshot by total number of galaxies / volume
        n_mean    = d[0].header['NGALBOX']/(1.* boxL**3)
        length    = d[1].header['NAXIS2']

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
        write_catalog(OUT_file, ra_desi, dec_desi, zz_desi, zz_rsd_desi)

        end = time.perf_counter()
        print(f'Took {end - begin} seconds to process galaxy shell {i}, galtype {galtype}.')
    else:
        print(f"WARNING: File {OUT_file} exists, skipping.")


    begin = time.perf_counter()

    for j in range(nrandoms):
        OUT_file = f"{dir_out}cutsky/shells/randoms_{j}_lightcone_multibox_galtype{galtype}_{i}_footprint_nz.fits"
        if not os.path.exists(OUT_file):
            ra0       = np.random.rand(length)*360
            dec0      = np.arccos(2*np.random.rand(length)-1)*90/np.pi*2-90

            print(f"Downsampling randoms file {j} from {i} shell for galtype {galtype}")
            ra, dec, zz, zz_rsd, nzz = downsample(ra0, dec0, zz0, zz_rsd0,
                                                  n_mean, nz_par)

            print(f"Applying footprint for randoms file {j} from {i} shell of galtype {galtype}")
            ra_desi, dec_desi, zz_desi, zz_rsd_desi, nzz_desi = apply_footprint(ra, dec, zz, zz_rsd, nzz)

            print('Writing footprinted catalog')
            write_catalog(OUT_file, ra_desi, dec_desi, zz_desi, zz_rsd_desi)
        else:
            print(f"WARNING: File {OUT_file} exists, skipping.")

    end = time.perf_counter()
    print(f'Took {end - begin} seconds to process randoms for shell {i}, galtype {galtype}')


### Loading config
parser = argparse.ArgumentParser()
parser.add_argument("config", help="ini file holding configuration",
                    type=str)
parser.add_argument("--dir_out", type=str, help="output directory (overrides config file)")
parser.add_argument("--shellnums", type=str, help="list of comma separated shell numbers to compute (same)")

args = parser.parse_args()

configfile     = str(args.config) # config file

config         = configparser.ConfigParser()
config.read(configfile)

file_alist        = config.get('dir','file_alist')
file_camb         = config.get('dir','file_camb')
dir_out           = args.dir_out
shellnums         = args.shellnums
if dir_out is None:
    dir_out       = config.get('dir','dir_out')
if shellnums is None:
    try:
        shellnums = config.get('sim', 'shellnums')
    except:
        shellnums = None
boxL              = config.getint('sim', 'boxL')
shellwidth        = config.getint('sim','shellwidth')
if shellnums is None:
    zmin          = config.getfloat('sim', 'zmin')
    zmax          = config.getfloat('sim', 'zmax')

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


if not os.path.exists(f"{dir_out}/cutsky"):
    os.system(f"mkdir {dir_out}/cutsky")
if not os.path.exists(f"{dir_out}/cutsky/shells"):
    os.system(f"mkdir {dir_out}/cutsky/shells")


### Running camb to get comoving distances

# Load all parameters from camb file 
pars = camb.read_ini(file_camb)
h    = pars.h
pars.set_for_lmax(2000, lens_potential_accuracy=3)
pars.set_matter_power(redshifts=[0.], kmax=200.0)
pars.NonLinearModel.set_params(halofit_version='takahashi')
camb.set_feedback_level(level=100)
results   = camb.get_results(pars)

if shellnums is None:
    shellnum_min = int(results.comoving_radial_distance(zmin)*h // shellwidth)
    shellnum_max = int(results.comoving_radial_distance(zmax)*h // shellwidth + 1)
    shellnums = list(range(shellnum_min, shellnum_max+1))
else:
    shellnums = list(map(int, shellnums.split(",")))

try:
    pool = MPIPool()
except:
    pool = None


### Loop over galtypes and shells for galaxies and randoms
args = product(shellnums, nz_pars, [nrandoms], [dir_out])
if pool is not None:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    else:
        pool.map(process_shell, args)
else:
    map(process_shell, args)


### Merge the processed shells
flag = pool is None
if not flag:
    flag = pool.is_master()

if flag:
    for nz_par in nz_pars:
        galtype = nz_par["galtype_index"]

        out_fname = f"{dir_out}cutsky/UNIT_lightcone_multibox_{nz_par['galtype'].upper()}_footprint_nz.fits"
        if not os.path.exists(out_fname):
            print('Merging shells into single file')
            begin = time.perf_counter()

            shells = []
            shell_fnames = f"{dir_out}cutsky/shells/lightcone_multibox_galtype{galtype}_*_footprint_nz.fits"
            for f in glob.glob(shell_fnames):
                shells.append(Table.read(f))
            joint_table = vstack(shells)
            joint_table.write(out_fname)
            del joint_table, shells
            #os.system(f"rm {shell_fnames}")
            end = time.perf_counter()
            print(f'Done! (took {end - begin} seconds)')
        else:
            print(f"File {out_fname} exists, skipping.")

        print('Merging random shells into single file')
        begin = time.perf_counter()
        for j in range(nrandoms):
            out_1x = f"{dir_out}cutsky/randoms_1x_lightcone_multibox_{nz_par['galtype'].upper()}_footprint_nz_{j}.fits"
            if not os.path.exists(out_1x):
                shells = []
                shell_fnames = f"{dir_out}cutsky/shells/randoms_{j}_lightcone_multibox_galtype{galtype}_*_footprint_nz.fits"
                for f in glob.glob(shell_fnames):
                    shells.append(Table.read(f))
                random_1x_joint = vstack(shells)
                random_1x_joint.write(out_1x)
                del random_1x_joint, shells
                #os.system(f"rm {shell_fnames}")
            else:
                print(f"File {out_1x} exists, skipping.")
        end = time.perf_counter()
        print(f'Done! (took {end - begin} seconds)')

        out_nx = f"{dir_out}cutsky/randoms_{nrandoms}x_lightcone_multibox_{nz_par['galtype'].upper()}_footprint_nz.fits"
        if not os.path.exists(out_nx):
            print('Joining all randoms')
            begin = time.perf_counter()
            random_1x_fnames = f"{dir_out}cutsky/randoms_1x_lightcone_multibox_{nz_par['galtype'].upper()}_footprint_nz_*.fits"

            randoms_1x = []
            for f in glob.glob(random_1x_fnames):
                randoms_1x.append(Table.read(f))
            joint_randoms = vstack(randoms_1x)
            joint_randoms.write(out_nx)
            del joint_randoms, randoms_1x

            end = time.perf_counter()
            print(f'Done! (took {end - begin} seconds)')
        else:
            print(f"File {out_nx} exists, skipping.")

sys.exit(0)


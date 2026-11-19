import glob
import os
from functools import lru_cache

import numpy as np
from astropy.table import Table

# Set path with atmosphere templates
DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')

# Default set of atmospheres: PHOENIX ([M/H] = -2) below 4000 K and Kurucz (solar) above, see data/README.md
ATMOS_DIR = os.path.join(DATA_PATH, 'atmos')

# The PHOENIX atmospheres are used below this temperature, the Kurucz atmospheres at and above it
T_SWITCH = 4000.


@lru_cache(maxsize=None)
def read_atmospheres(directory):
    """Effective temperatures (sorted) and tables of the atmospheres in a directory. The temperature is
    read from characters 2-7 of the file names, e.g. pt02000g+4.0z-2.0.hdf5."""
    models = glob.glob(os.path.join(directory, '*.hdf5'))
    if not models:
        raise FileNotFoundError("No atmosphere tables (*.hdf5) in {0}".format(directory))
    teff = np.array([os.path.basename(model)[2:7] for model in models], dtype=np.float32)
    order = np.argsort(teff)
    return teff[order], [Table.read(models[i], path="Table") for i in order]


def interp_atmos(tval, atmos_dir=None):
    """
    Spectrum (nu, fnu) of a star with effective temperature tval, interpolated linearly between the two
    closest tables.

    atmos_dir : directory with 'phoenix' (used below T_SWITCH) and 'kurucz' sub-directories of tables
        in the format of the bundled ones (default: the bundled set, ATMOS_DIR). A directory built by
        scripts/make_atmospheres.py has one metallicity at all temperatures.
    """

    if atmos_dir is None:
        atmos_dir = ATMOS_DIR

    # Decide which set of atmospheres to use
    family = 'phoenix' if tval < T_SWITCH else 'kurucz'
    teff, tables = read_atmospheres(os.path.join(os.path.abspath(atmos_dir), family))

    # Locate the temperature in the array
    i = np.searchsorted(teff, tval)

    # Retrieve two neighboring tables
    t1, t2 = tables[i - 1], tables[i]

    # Compute the interpolation fraction
    frac = (tval - teff[i - 1]) / (teff[i] - teff[i - 1])

    # Interpolate and return on the fly
    return t1['nu'], t1['fnu'] * (1. - frac) + t2['fnu'] * frac

"""Build a set of mysg stellar-atmosphere tables with one metallicity at every temperature.

    python scripts/make_atmospheres.py MH OUTDIR [--cache DIR]

writes

    OUTDIR/phoenix/pt{T:05d}g+4.0z{MH:+.1f}.hdf5   T = 2000-4000 K in steps of 100 K
    OUTDIR/kurucz/kt{T:05d}g+4.0z{MH:+.1f}.hdf5    the 65 temperatures of the bundled Kurucz set

so that ``interp_atmos(T, atmos_dir=OUTDIR)`` (or ``setup_model(..., atmos_dir=OUTDIR)``) uses atmospheres
with [M/H] = MH at every temperature.  Unlike the bundled files, the ``z`` in these names is the [M/H] of
both families.  All spectra have log g = 4.0 and scaled-solar abundances (no alpha enhancement).

* PHOENIX (T < 4000 K): BT-Settl with the Asplund et al. (2009) abundances (Allard et al. 2012,
  "BT-Settl AGSS2009"), from the SVO Theory service (model bt-settl-agss, alpha = 0).  Each spectrum is
  resampled, conserving flux, onto the wavelength grid of the STScI CDBS file
  grid/phoenix/phoenixm00/phoenixm00_3000.fits (103,000 points from 12.6 Angstrom to 77.6 micron).
  CDBS phoenixm00 is the same BT-Settl AGSS2009 [M/H] = 0 grid (delivered 2021) and is the source of
  the solar tables used for the R24 v1.3 models (T. Richardson 2026); the other CDBS PHOENIX
  metallicities are an older (2011) delivery that subsamples the spectra point by point, which
  misrepresents molecular bands, so they are not used.  BT-Settl AGSS2009 has scaled-solar spectra
  at log g = 4.0 for 2000-4000 K at [M/H] = -0.5, 0.0, +0.3 and +0.5 (at -1.0 only alpha-enhanced
  spectra, above 2600 K).
  Converted as for the R24 v1.3 tables: columns ``wav`` (micron), ``nu`` and ``fnu`` with increasing
  frequency, the flux kept proportional to the surface flux (scaled by 1 arcsec^2 / (0.01 pi sr)), so
  these tables are not normalized, and one extra point at 1 cm.  Here that point continues the
  spectrum from 77.6 micron with the Rayleigh-Jeans slope, F_nu ~ lambda^-2; the R24 v1.3 tables imply
  slopes of -1.96 to -1.97, so their 1 cm points are up to 15% higher.  The point carries less than
  1e-6 of the luminosity and sets only the photospheric flux beyond 77.6 micron.
* Kurucz (T >= 4000 K): Castelli & Kurucz (2004) ATLAS9 models, STScI CDBS
  grid/ck04models/ck{tag}/ck{tag}_{T}.fits (tags m25 ... m05, p00, p02, p05), column g40.  Converted as
  the bundled solar files: one point at 1 cm extrapolated in log-log space from the last two points,
  columns ``nu`` and ``fnu`` with increasing frequency, normalized so that the integral of fnu over
  frequency is 1.

Checks for MH = 0: see data/README.md.  Downloads are cached in --cache (default OUTDIR/download).
"""
import argparse
import io
import os
import urllib.request

import h5py
import numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table
from scipy.integrate import trapezoid

CDBS = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid'
SVO = 'https://svo2.cab.inta-csic.es/theory/newov2/ssap.php'
C_ANGSTROM = 2.99792458e18  # speed of light in Angstrom / s
LAM_EXTRA = 1e8  # Angstrom (1 cm): the extra point at long wavelength
ARCSEC2_SR = (np.pi / 180. / 3600.) ** 2
T_PHOENIX = np.arange(2000, 4001, 100)
T_KURUCZ = np.concatenate([np.arange(3500, 13001, 250), np.arange(14000, 39001, 1000)])


def download(url, path):
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        urllib.request.urlretrieve(url, path + '.part')
        os.replace(path + '.part', path)
    return path


def cdbs_g40(relpath, cache):
    """Wavelength (Angstrom, increasing) and F_lambda of the log g = 4.0 spectrum of a CDBS grid file."""
    path = download(f'{CDBS}/{relpath}', os.path.join(cache, 'cdbs', relpath))
    with fits.open(path) as hdul:
        data = hdul[1].data
        lam = np.array(data['WAVELENGTH'], dtype=float)
        flam = np.array(data['g40'], dtype=float)
    if not np.any(flam > 0):
        raise ValueError(f'{path}: the log g = 4.0 column (g40) is empty')
    return lam, flam


def svo_bt_settl(teff, mh, cache):
    """Wavelength (Angstrom, increasing) and F_lambda of the BT-Settl AGSS2009 spectrum with this Teff,
    log g = 4.0, [M/H] = mh and alpha = 0, from the SVO Theory service."""
    path = os.path.join(cache, 'svo_bt-settl-agss', f'bt-settl-agss_{teff:05d}_g4.0_m{mh:+.1f}_a+0.0.txt')
    if not os.path.exists(path):
        query = (f'{SVO}?model=bt-settl-agss&teff_min={teff}&teff_max={teff}&logg_min=4&logg_max=4'
                 f'&meta_min={mh}&meta_max={mh}&REQUEST=queryData')
        with urllib.request.urlopen(query) as response:
            table = parse_single_table(io.BytesIO(response.read())).to_table()
        rows = table[np.isclose(table['alpha'], 0.) & np.isclose(table['meta'], mh)
                     & np.isclose(table['teff'], teff) & np.isclose(table['logg'], 4.)]
        if len(rows) != 1:
            raise ValueError(f'SVO bt-settl-agss: {len(rows)} spectra for Teff = {teff}, log g = 4, '
                             f'[M/H] = {mh}, alpha = 0')
        download(str(rows['Access.Reference'][0]) + '&format=ascii', path)
    header = {}
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            if '=' in line:
                key, value = line[1:].split('=', 1)
                header[key.strip()] = float(value.split()[0])
    expected = {'teff': teff, 'logg': 4., 'meta': mh, 'alpha': 0.}
    if any(not np.isclose(header.get(k, np.nan), v) for k, v in expected.items()):
        raise ValueError(f'{path}: header {header} does not match {expected}')
    lam, flam = np.loadtxt(path, unpack=True)
    order = np.argsort(lam, kind='stable')
    lam, flam = lam[order], flam[order]
    keep = np.concatenate([[True], np.diff(lam) > 0])
    return lam[keep], flam[keep]


def rebin(lam, flam, lam_new):
    """Mean of flam over the bins centred on lam_new (edges at the geometric midpoints), from the
    trapezoid integral of the input spectrum: conserves the flux in each bin."""
    mid = np.sqrt(lam_new[1:] * lam_new[:-1])
    edges = np.concatenate([[lam_new[0] ** 2 / mid[0]], mid, [lam_new[-1] ** 2 / mid[-1]]])
    if edges[0] < lam[0] or edges[-1] > lam[-1]:
        raise ValueError('the input spectrum does not cover the output wavelength grid')
    cumulative = np.concatenate([[0.], np.cumsum(0.5 * (flam[1:] + flam[:-1]) * np.diff(lam))])
    return np.diff(np.interp(edges, lam, cumulative)) / np.diff(edges)


def phoenix_table(lam, fnu):
    f_extra = fnu[-1] * (LAM_EXTRA / lam[-1]) ** -2.
    lam = np.append(lam, LAM_EXTRA)
    fnu = np.append(fnu, f_extra) / (np.pi / ARCSEC2_SR * 1e-2)
    table = np.zeros(len(lam), dtype=[('wav', 'f8'), ('nu', 'f8'), ('fnu', 'f8')])
    table['wav'] = lam * 1e-4
    table['nu'] = C_ANGSTROM / lam
    table['fnu'] = fnu
    return table[::-1]


def kurucz_table(lam, fnu):
    slope = np.log10(fnu[-1] / fnu[-2]) / np.log10(lam[-1] / lam[-2])
    f_extra = fnu[-1] * (LAM_EXTRA / lam[-1]) ** slope
    lam = np.append(lam, LAM_EXTRA)
    fnu = np.append(fnu, f_extra)
    table = np.zeros(len(lam), dtype=[('nu', 'f8'), ('fnu', 'f8')])
    table['nu'] = (C_ANGSTROM / lam)[::-1]
    table['fnu'] = fnu[::-1]
    table['fnu'] /= trapezoid(table['fnu'], table['nu'])
    return table


def write_table(path, table):
    with h5py.File(path, 'w') as f:
        f['Table'] = table


def build(mh, outdir, cache):
    for family in ('phoenix', 'kurucz'):
        os.makedirs(os.path.join(outdir, family), exist_ok=True)
    lam_grid, _ = cdbs_g40('phoenix/phoenixm00/phoenixm00_3000.fits', cache)
    for temp in T_PHOENIX:
        lam, flam = svo_bt_settl(int(temp), mh, cache)
        fnu = rebin(lam, flam, lam_grid) * lam_grid ** 2 / C_ANGSTROM
        write_table(os.path.join(outdir, 'phoenix', f'pt{temp:05d}g+4.0z{mh:+.1f}.hdf5'), phoenix_table(lam_grid, fnu))
    print(f'[M/H] = {mh:+.1f}: {len(T_PHOENIX)} phoenix tables (SVO bt-settl-agss, alpha = 0)')
    tag = 'p00' if mh == 0 else ('p' if mh > 0 else 'm') + f'{round(abs(mh) * 10):02d}'
    for temp in T_KURUCZ:
        lam, flam = cdbs_g40(f'ck04models/ck{tag}/ck{tag}_{temp}.fits', cache)
        write_table(os.path.join(outdir, 'kurucz', f'kt{temp:05d}g+4.0z{mh:+.1f}.hdf5'),
                    kurucz_table(lam, flam * lam ** 2 / C_ANGSTROM))
    print(f'[M/H] = {mh:+.1f}: {len(T_KURUCZ)} kurucz tables (CDBS ck04models/ck{tag})')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('mh', type=float, help='[M/H], e.g. -0.5')
    parser.add_argument('outdir')
    parser.add_argument('--cache', help='directory for the downloaded spectra (default OUTDIR/download)')
    args = parser.parse_args()
    build(args.mh, args.outdir, args.cache or os.path.join(args.outdir, 'download'))


if __name__ == '__main__':
    main()

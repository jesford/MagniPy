"""Functions specific to reading & loading the CFHTLenS formatted data."""
import pandas as pd
import numpy as np
import glob
import treecorr
from astropy.cosmology import Planck13 as cosmo


def make_dataframe(fname, header=0):
    """Return a dataframe for specific file."""
    if header is None:
        df = pd.read_table(fname, header=header, delim_whitespace=True,
                           names=['x[0]', 'x[1]'])
    else:
        df = pd.read_table(fname, header=header, delim_whitespace=True)
        df.columns = list(df.columns)[1:]+['nan']
        df = df.drop('nan', axis=1)
    return df


def load_all_pointings(field, path=''):
    """Return a dictionary of data from all pointings in field."""
    if field not in ['W1', 'W2', 'W3', 'W4']:
        raise ValueError("field must be one of 'W1', 'W2', 'W3', 'W4'")

    if 'clusters' in path:
        h = 13
    elif 'randoms' in path:
        h = None
    else:
        h = 27

    files = glob.glob(path+'*'+field+'*')
    data = {}

    for f in files:
        pointing = f[-10:-4]  # depends on these filenames!
        df = make_dataframe(f, header=h)
        data[pointing] = df

    return data


def pix2rad(arr_pix, scale='pointings'):
    """Convert CFHTLenS pixels to angle in radians."""
    if scale == 'pointings':
        # pointings are 0.0031 am/pix
        ampix = 0.0031
        # mosaic is 0.0166667 am/pix
    elif scale == 'mosaic':
        ampix = 0.0166667
    else:
        raise ValueError('scale must be "pointings" or "mosaic"')
    if (np.array(arr_pix) < 0.).any():
        raise ValueError('pixel values must be positive')
    arr_deg = arr_pix * ampix / 60.
    arr_rad = np.deg2rad(arr_deg)
    return arr_rad


def get_catalogs(df_lens=None, df_weight=None, df=None, redshift=None):
    """Return TreeCorr Catalog objects in Mpc for each dataframe."""
    if redshift is None:
        raise ValueError('redshift must be a non-negative float')

    dA_lens = cosmo.angular_diameter_distance(redshift)

    if df_lens is not None:
        if type(df_lens) != pd.core.frame.DataFrame:
            raise TypeError('cdf must be a dataframe')
        cdf_zslice = df_lens[np.isclose(df_lens.z, redshift)]
        x_lens_mpc = pix2rad(cdf_zslice['x[0]']) * dA_lens
        y_lens_mpc = pix2rad(cdf_zslice['x[1]']) * dA_lens
        lenses = treecorr.Catalog(x=x_lens_mpc, y=y_lens_mpc)
    else:
        lenses = None

    if df_weight is not None:
        if type(df_weight) != pd.core.frame.DataFrame:
            raise TypeError('df_weight must be a dataframe')
        x_source_mpc = pix2rad(df_weight['x[0]']) * dA_lens
        y_source_mpc = pix2rad(df_weight['x[1]']) * dA_lens
        sources = treecorr.Catalog(x=x_source_mpc, y=y_source_mpc,
                                   k=df_weight.am1)
    else:
        sources = None

    if df is not None:
        if type(df) != pd.core.frame.DataFrame:
            raise TypeError('df must be a dataframe')
        x_rand_mpc = pix2rad(df['x[0]']) * dA_lens
        y_rand_mpc = pix2rad(df['x[1]']) * dA_lens
        randoms = treecorr.Catalog(x=x_rand_mpc, y=y_rand_mpc)
    else:
        randoms = None

    return lenses, sources, randoms

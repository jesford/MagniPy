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


def get_catalog(df, weights=None, shear=False, lens=False, z_lens=None):
    """Return TreeCorr Catalog in Mpc."""
    if z_lens is None:
        raise ValueError('must enter a value for z_lens')
    else:
        dA_lens = cosmo.angular_diameter_distance(z_lens)
    if type(df) != pd.core.frame.DataFrame:
        raise TypeError('df must be a dataframe')

    if lens is True:
        cdf_zslice = df[np.isclose(df.z, z_lens)]
        x_lens_mpc = pix2rad(cdf_zslice['x[0]']) * dA_lens
        y_lens_mpc = pix2rad(cdf_zslice['x[1]']) * dA_lens
        cat = treecorr.Catalog(x=x_lens_mpc, y=y_lens_mpc)

    elif shear is True:
        # insert function to select on zphot & P(z) conditions
        # sdf_z = df[df.zphot > z_lens && 90% P(z) above z_lens]
        x_shear_mpc = pix2rad(df['x[0]']) * dA_lens  # tbr
        y_shear_mpc = pix2rad(df['x[1]']) * dA_lens  # tbr

        cat = treecorr.Catalog(x=x_shear_mpc, y=y_shear_mpc,
                               g1=df['e[0]'], g2=df['e[1]'])
    else:
        x_mpc = pix2rad(df['x[0]']) * dA_lens
        y_mpc = pix2rad(df['x[1]']) * dA_lens
        cat = treecorr.Catalog(x=x_mpc, y=y_mpc)

    if weights is not None:
        cat.k = np.ones(df.shape[0])
        for w in weights:
            cat.k *= df[w].values

    return cat

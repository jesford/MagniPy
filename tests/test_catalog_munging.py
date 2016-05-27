from numpy.testing import assert_equal, assert_raises
# , assert_allclose
import glob
import numpy as np
import pandas as pd


from ..catalog_munging import make_dataframe, load_all_pointings, pix2rad, \
    get_catalogs


fields = {'W1': 72, 'W2': 25, 'W3': 49, 'W4': 25}
basepath = '/Users/jesford/astrophysics/data/cfhtls/'
spath = basepath+'LBGS/WIDE/'
cpath = basepath+'clusters/WIDE/'
rpath = basepath+'randoms/'
sample_pointings = ['W1m2p1', 'W2p3m0', 'W3m0p2', 'W4m1p3']


def test_column_names():
    expect_LBG_cols = ['x[0]', 'x[1]', 'RA', 'DEC', 'MAG_TOT_g', 'MAG_TOT_r',
                       'MAG_AUTO', 'MAG_TOT_z', 'am1', 'z', 'dmag_g',
                       'dmag_r', 'dmag_i', 'dmag_z']
    expect_cluster_cols = ['RA', 'DEC', 'z', 'sig', 'x[0]', 'x[1]', 'm200',
                           'r200', 'n200']
    expect_random_cols = ['x[0]', 'x[1]']

    def _check_names(p, dpath, h):
        if 'clusters' in dpath:
            expect_cols = expect_cluster_cols
        elif 'randoms' in dpath:
            expect_cols = expect_random_cols
        else:
            expect_cols = expect_LBG_cols
        fname = glob.glob(dpath+'*'+p+'*')[0]
        df = make_dataframe(fname, header=h)
        assert_equal(list(df.columns), expect_cols)

    for pointing in sample_pointings:
        for datapath, h in zip([spath+'udrops/', spath+'gdrops/', cpath, rpath],
                               [27, 27, 13, None]):
            yield _check_names, pointing, datapath, h


def test_number_pointings():
    def _check_num(f, dpath):
        data = load_all_pointings(f, path=dpath)
        assert_equal(len(data), fields[f])

    for field in fields.keys():
        for datapath in [spath+'udrops/', spath+'gdrops/', cpath, rpath]:
            yield _check_num, field, datapath


def test_pix2rad():
    def _check_rad(p, r, s):
        print('p, r, s', p, r, s)
        rad = pix2rad(p, scale=s)
        print('rad', rad)
        assert_equal(rad, r)

    input_pixels = [0.0, 180 * 60 / np.pi]
    output_radians = {'pointings': [0.0, 0.0031],
                      'mosaic': [0.0, 0.0166667]}

    for i, pix in enumerate(input_pixels):
        for scale in ['pointings', 'mosaic']:
            yield _check_rad, pix, output_radians[scale][i], scale

    def _check_badinputs(in1, in2):
        assert_raises(ValueError, pix2rad, in1, in2)

    for v1, v2 in zip([10., -1], ['pointing', 'pointings']):
        yield _check_badinputs, v1, v2


def test_get_catalogs():
    def _check_dfs_are_None():
        result = get_catalogs(None, None, None, redshift=0.2)
        assert_equal(result, (None, None, None))

    def _check_redshift_is_None():
        assert_raises(ValueError, get_catalogs, None, None, None)

    yield _check_redshift_is_None
    yield _check_dfs_are_None

    def _check_zeros(df, i):
        df1, df2, df3 = None, None, None
        if i == 0:
            df1 = df
        elif i == 1:
            df2 = df
        elif i == 2:
            df3 = df
        result = get_catalogs(df_lens=df1, df_weight=df2, df=df3, redshift=0.0)
        assert_equal(result[i].x, np.zeros(3))

    df = pd.DataFrame(np.zeros([3, 4]), columns=['x[0]', 'x[1]', 'z', 'am1'])
    for i in range(0, 3):
        yield _check_zeros, df, i

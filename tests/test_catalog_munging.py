from numpy.testing import assert_equal
# , assert_allclose, assert_raises
import glob

from ..catalog_munging import make_dataframe, load_all_pointings


fields = {'W1': 72, 'W2': 25, 'W3': 49, 'W4': 25}
basepath = '/Users/jesford/astrophysics/data/cfhtls/'
spath = basepath+'LBGS/WIDE/'
cpath = basepath+'clusters/WIDE/'
sample_pointings = ['W1m2p1', 'W2p3m0', 'W3m0p2', 'W4m1p3']


def test_column_names():
    expect_LBG_cols = ['x[0]', 'x[1]', 'RA', 'DEC', 'MAG_TOT_g', 'MAG_TOT_r',
                       'MAG_AUTO', 'MAG_TOT_z', 'am1', 'z', 'dmag_g',
                       'dmag_r', 'dmag_i', 'dmag_z']
    expect_cluster_cols = ['RA', 'DEC', 'z', 'sig', 'x[0]', 'x[1]', 'm200',
                           'r200', 'n200']

    def _check_names(p, dpath):
        if 'clusters' in dpath:
            expect_cols = expect_cluster_cols
        else:
            expect_cols = expect_LBG_cols
        fname = glob.glob(dpath+'*'+p+'*cat')[0]
        df = make_dataframe(fname)
        assert_equal(list(df.columns), expect_cols)

    for pointing in sample_pointings:
        for datapath in [spath+'udrops/', spath+'gdrops/', cpath]:
            yield _check_names, pointing, datapath


def test_number_pointings():
    def _check_num(f, dpath):
        data = load_all_pointings(f, path=dpath)
        assert_equal(len(data), fields[f])

    for field in fields.keys():
        for datapath in [spath+'udrops/', spath+'gdrops/', cpath]:
            yield _check_num, field, datapath

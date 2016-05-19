from numpy.testing import assert_equal
# , assert_allclose, assert_raises
import glob

from ..catalog_munging import make_dataframe, load_all_pointings


fields = {'W1': 72, 'W2': 25, 'W3': 49, 'W4': 25}
path = '/Users/jesford/astrophysics/data/cfhtls/LBGS/WIDE/'
sample_pointings = ['W1m2p1', 'W2p3m0', 'W3m0p2', 'W4m1p3']


def test_column_names():
    expected_names = ['x[0]', 'x[1]', 'RA', 'DEC', 'MAG_TOT_g', 'MAG_TOT_r',
                      'MAG_AUTO', 'MAG_TOT_z', 'am1', 'z', 'dmag_g', 'dmag_r',
                      'dmag_i', 'dmag_z']

    def _check_names(p, d):
        fname = glob.glob(path+d+'drops/'+d+'*'+p+'*cat')[0]
        df = make_dataframe(fname)
        assert_equal(list(df.columns), expected_names)

    for pointing in sample_pointings:
        for dropout in ['u', 'g']:
            yield _check_names, pointing, dropout


def test_number_pointings():
    def _check_num(f, d):
        data = load_all_pointings(f, dropoutsample=d, path=path+d+'drops/')
        assert_equal(len(data), fields[f])

    for field in fields.keys():
        for dropout in ['u', 'g']:
            yield _check_num, field, dropout

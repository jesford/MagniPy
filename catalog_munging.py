import pandas as pd
import glob


def make_dataframe(fname):
    """Return a dataframe for specific file."""
    # how to generalize this for clusters (header=1)?
    try:
        df = pd.read_table(fname, delim_whitespace=True)
    except pd.parser.CParserError:
        df = pd.read_table(fname, header=27, delim_whitespace=True)
    except:
        raise ValueError('Unexpected header (expect 1 or 27 rows):\n'+fname)
    df.columns = list(df.columns)[1:]+['nan']
    df = df.drop('nan', axis=1)
    return df


def load_all_pointings(field, path='', dropoutsample='u'):
    """Return a dictionary of data from all pointings in field."""
    if field not in ['W1', 'W2', 'W3', 'W4']:
        raise ValueError("field must be one of 'W1', 'W2', 'W3', 'W4'")

    files = glob.glob(path+dropoutsample+'*'+field+'*cat')
    data = {}

    for f in files:
        pointing = f[-10:-4]  # depends on these filenames!
        df = make_dataframe(f)
        data[pointing] = df

    return data

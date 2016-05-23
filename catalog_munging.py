"""Functions specific to reading & loading the CFHTLenS formatted data."""
import pandas as pd
import glob


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

import pandas as pd

from snotra import listify, to_df, fix_args, get_rows, count_codes


def single_columns(df, cols=None, sep=',', n=100, check_all=False):
    """
    Identify columns that do not have seperators i.e. single values in cells

    Args:
        cols (list of strings): columns to be examined for seperators,
        default: all columns

        sep (string): seperator that may be used to distinguish values in cells
        n (int): check only a subsample (head and tail) of n observations
        check_all (bool): check all observations

    Returns:
        list

    """

    if not cols:
        cols = list(df.columns)

    single_value_columns = []
    multiple_value_columns = []

    for col in cols:
        if (df[col].head(100).str.contains(sep).any()) or (df[col].tail(100).str.contains(sep).any()):
            multiple_value_columns.append(col)
        else:
            single_value_columns.append(col)

    if check_all:
        for col in single_value_columns:
            if df[col].str.contains(sep).any():
                multipe_value_columns.append(col)
                single_value_columns.remove(col)
    return single_value_columns


def find_spikes(df, codes=None, cols=None, persons=False, pid='pid', sep=None, groups=None,
                each_group=False, _fix=True, threshold=3, divide_by='pid'):
    """
    Identifies large increases or decreases in use of given codes in the specified groups
    rename? more like an increase identifier than a spike ideintifier as it is
    spikes implies relatively low before and after comparet to the "spike"
    rem: spikes can also be groups of years (spike in two years, then down again)

    cols='ncmp'
    df=npr.copy()
    codes='4AB04'
    sep=','
    pid='pid'
    groups='region'
    threshold=3
    divide_by='pid'
    """
    sub = df
    groups = listify(groups)

    if _fix:
        sub, cols = to_df(sub, cols)
        codes, cols, allcodes, sep = fix_args(df=sub, codes=codes, cols=cols, sep=sep, merge=False, group=False)
        rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
        sub = sub[rows]

    if persons:
        counted = sub.groupby(groups).count_persons(codes=codes, cols=cols, sep=sep, _fix=False)
    else:
        counted = sub.groupby(groups).apply(count_codes, codes=codes, cols=cols, sep=sep)

    if divide_by:
        divisor = sub.groupby(groups)[divide_by].nunique()
        counted = counted / divisor

    avg = counted.mean()
    sd = counted.std()
    counted.plot.bar()
    deviations = (counted - avg) / sd
    deviations = (counted / avg) / avg
    spikes = counted(deviations.abs() > threshold)

    return spikes


def find_shifts(df, codes=None, cols=None, sep=None, groups=None, interact=False, _fix=True, threshold=3):
    """
    Identifies large increases or decreases in use of given codes in the specified groups
    rename? more like an increase identifier than a spike ideintifier as it is
    spikes implies relatively low before and after comparet to the "spike"
    rem: spikes can also be groups of years (spike in two years, then down again)

    """
    find_changes()
    #    do_mocing average and reverse ma.
    # use shorter then whole period window if think there may be more than one shift

    return


def find_cycles(df, codes=None, cols=None, sep=None, groups=None, interact=False, _fix=True, threshold=3):
    """
    Identifies large increases or decreases in use of given codes in the specified groups
    rename? more like an increase identifier than a spike ideintifier as it is
    spikes implies relatively low before and after comparet to the "spike"
    rem: spikes can also be groups of years (spike in two years, then down again)

    """
    find_changes()
    return


def find_changes(df, codes=None, cols=None, sep=None, groups=None, interact=False, _fix=True, threshold=3):
    """
    Identifies large increases or decreases in use of given codes in the specified groups
    rename? more like an increase identifier than a spike ideintifier as it is
    spikes implies relatively low before and after comparet to the "spike"
    rem: spikes can also be groups of years (spike in two years, then down again)

    """
    sub = df
    groups = listify(groups)

    if _fix:
        df, cols = to_df(df, cols)
        codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep, merge=False, group=False)
        rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
        sub = df[rows]

    all_groups = {}

    for group in groups:
        counted = []
        names = []
        for name, codelist in codes.items():
            count = sub.groupby(group).apply(count_codes,
                                             codes={name: codelist},
                                             cols=cols,
                                             sep=sep,
                                             dropna=True,
                                             _fix=False)
            counted.append(count)
            # names.append(name)

        counted = pd.concat(counted, axis=1)
        # counted.columns=names

        if threshold:
            counted_delta = counted.pct_change() / counted.pct_change().abs().mean()
            counted_delta = counted_delta[counted_delta > threshold]
            counted = counted_delta

        all_groups[group] = counted

    return all_groups
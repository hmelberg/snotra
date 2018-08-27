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
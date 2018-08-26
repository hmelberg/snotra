import re

from snotra.core import expand_codes, count_codes, get_rows, unique_codes


def listify(string_or_list):
    """
    return a list if the input is a string, if not: returns the input as it was

    Args:
        string_or_list (str or any):

    Returns:
        A list if the input is a string, if not: returns the input as it was

    Note:
        - allows user to use a string as an argument instead of single lists
        - cols='icd10' is allowed instead of cols=['icd10']
        - cols='icd10' is transformed to cols=['icd10'] by this function

    """
    if isinstance(string_or_list, str):
        string_or_list = [string_or_list]
    return string_or_list


def sniff_sep(df, cols=None, possible_seps=[',', ';', '|'], n=1000, sure=False, each_col=False):
    """
    Sniff whether column(s) cells have mulitple values with a seperator

    Args:
        df: dataframe or series
        cols (str or list of str): name of columns to be checked
        possible_seps (list of str): list of potential seperators to check for
        n (int, default = 1000): number of rows from tail and head to check
        sure (bool, default False): Set to True to check all rows
        each_col (bool, default False): Set to True to get a dict with the seperator (and whether it exists) for each column

    Return:
        Str, None or dict (of str or None)
        Returns the seperator if it is found, or None if no seperator is found
        If each_col is True, returns a dict with the sep (or None) for each column

    """

    cols = listify(cols)

    df = stringify_cols(df=df, cols=cols)

    # fix args depending on whther a series or df is input and if cols is specified
    if isinstance(df, pd.Series):
        df = pd.DataFrame(df)
        cols = list(df.columns)
    else:
        if not cols:
            cols = list(df.columns)

    sep_col = {}
    for col in cols:
        if sure:
            n = len(df.dropna())

        if n < 1000:
            n = len(df.dronna())

        search_head = df[col].dropna().head(n).str.cat()
        search_tail = df[col].dropna().head(n).str.cat()

        # check for existence of all seps
        for sep in possible_seps:
            if (sep in search_head) or (sep in search_tail):
                sniffed_sep = sep
                break
            else:
                sniffed_sep = None
        # don't check more columnsif found sep in one
        if sniffed_sep and not each_col:
            break

        # go on to check each col if each_col is specified
        if each_col:
            sep_col[col] = sniffed_sep

    if each_col:
        sniffed_sep = sep_col

    return sniffed_sep


def get_some_id(df,
                codes,
                cols,
                xid,
                sep=None):
    """
    help function for all get functions that gets ids based on certain filtering criteria

    x is the column with the info to be collected (pid, uuid, event_id)


    """

    codes = listify(codes)
    cols = listify(cols)

    cols = expand_cols(df=df, cols=cols)

    expanded_codes = expand_codes(df=df,
                                  codes=codes,
                                  cols=cols,
                                  sep=sep)

    # if compound words in a cell
    if sep:
        expanded_codes_regex = '|'.join(expanded_codes)
        b = np.full(len(df), False)
        for col in cols:
            a = df[col].str.contains(expanded_codes_regex, na=False).values
            b = b | a
    # if single value cells only
    else:
        b = df[cols].isin(expanded_codes).any(axis=1).values

    pids = set(df[b][xid].unique())

    return pids


def fix_cols(df, cols):
    if not cols:
        cols = list(df.columns)

    cols = expand_cols(df=df, cols=cols)
    return cols


def fix_codes(df, codes=None, cols=None, sep=None, merge=False, group=False):
    if not codes:
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False)[:5]

    codes = format_codes(codes=codes, merge=merge)
    codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep, merge=merge, group=group)
    return codes


def fix_args(df, codes=None, cols=None, sep=None, merge=False, group=False, _sniffsep=True):
    # Use all columns if no column is specified
    # Series if converted to df (with pid column, assumed to be in the index)
    if not cols:
        cols = list(df.columns)
    else:
        cols = expand_cols(df=df, cols=cols)

    if _sniffsep:
        sep = sniff_sep(df=df, cols=cols)

    if not codes:
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False).index[:5]
        codes = list(codes)

    codes = format_codes(codes=codes, merge=merge)
    codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep, merge=merge, group=group)
    codes = format_codes(codes=codes, merge=merge)

    # useful to have full codelist (of codes only, after expansion)
    full_codelist = set()
    for name, codelist in codes.items():
        full_codelist.update(set(codelist))
    allcodes = list(full_codelist)

    return codes, cols, allcodes, sep


def to_df(df, cols=None):
    if isinstance(df, pd.Series):
        df = df.to_frame()
        cols = list(df.columns)
        df['pid'] = df.index.values
    return df, cols


def subset(df, codes, cols, sep):
    allcodes = get_allcodes(codes)
    rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep)
    subset = df[rows].set_index('pid')
    return subset


def get_allcodes(codes):
    """
    Return a list of only codes from the input

    Used when codes is a dict to extract codes only
    """
    if isinstance(codes, dict):
        allcodes = set()
        for name, codelist in codes.items():
            allcodes.update(set(codelist))
        allcodes = list(allcodes)
    else:
        allcodes = listify(codes)
    return allcodes


def get_mask(df,
             codes,
             cols,
             sep=None):
    codes = listify(codes)
    cols = listify(cols)

    cols = expand_cols(df=df, cols=cols)

    expanded_codes = expand_codes(df=df,
                                  codes=codes,
                                  cols=cols,
                                  sep=sep)

    # if compound words in a cell
    if sep:
        expanded_codes_regex = '|'.join(expanded_codes)
        b = pd.DataFrame()
        for col in cols:
            b[col] = df[col].str.contains(expanded_codes_regex, na=False).values

    # if single value cells only
    else:
        b = df[cols].isin(expanded_codes)

    return b


def expand_cols(df, cols, star=True, hyphen=True, colon=True, regex=None):
    """
    Expand columns with special notation to their full column names

    """

    cols = listify(cols)

    allcols = list(df.columns)

    if hyphen:
        cols = expand_hyphen(expr=cols)
    if star:
        cols = expand_star(expr=cols, full_list=allcols)
    if colon:
        cols = expand_colon(expr=cols, full_list=allcols)
    if regex:
        cols = list(df.columns(df.columns.str.contains(regex)))

    return cols


def expand_colon(expr, full_list):
    """
    Expand expressions with colon notation to a list of complete columns names

    expr (str or list): Expression (or list of expressions) to be expanded
    full_list (list or array) : The list to slice from
    """
    exprs = listify(expr)
    expanded = []

    for expr in exprs:
        if ':' in expr:
            startstr, endstr = expr.split(':')
            startpos = full_list.index(startstr)
            endpos = full_list.index(endstr) + 1
            my_slice = full_list[startpos:endpos]
        else:
            my_slice = [expr]

        expanded.extend(my_slice)
    return expanded


def expand_star(expr, cols=None, full_list=None, sep=None):
    """
    Expand expressions with star notation to all matching expressions

    """

    exprs = listify(expr)

    if isinstance(full_list, pd.Series):
        pass
    elif isinstance(full_list, list):
        unique_series = pd.Series(full_list)
    elif isinstance(full_list, set):
        unique_series = pd.Series(list(full_list))
    else:
        unique = unique_codes(df=df, cols=cols, sep=sep)
        unique_series = pd.Series(list(unique))

    expanded = []

    for expr in exprs:
        if '*' in expr:
            startstr, endstr = expr.split('*')
            if startstr:
                add_expr = list(unique_series[unique_series.str.startswith(startstr)])
            if endstr:
                add_expr = list(unique_series[unique_series.str.endswith(endstr)])
            if startstr and endstr:
                # col with single letter not included, start means one or more of something
                # beginning is not also end (here!)
                start_and_end = (unique_series.str.startswith(startstr)
                                 &
                                 unique_series.str.endswith(endstr))
                add_expr = list(unique_series[start_and_end])
        else:
            add_expr = [expr]

        expanded.extend(add_expr)
    return expanded


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


def expand_hyphen(expr):
    """
    Example: Expands ('b01A-b04A') to ['b01A' ,'b02A', 'b03A', 'b04A']

    Args:
        code

    Returns:
        List


    Examples:
        expand_hyphen('b01.1*-b09.9*')
        expand_hyphen('n02.2-n02.7')
        expand_hyphen('c00*-c260')
        expand_hyphen('b01-b09')
        expand_hyphen('b001.1*-b009.9*')
        expand_hyphen(['b001.1*-b009.9*', 'c11-c15'])

    Note:
        decimal expression also works: expr = 'n02.2-n02.7'
        expr = 'b01*-b09*'
        expr = 'C00*-C26*'

    """

    exprs = listify(expr)
    all_cod e s=[]

    for expr in exprs:
        if '-' in expr:
            lower, upper = expr.split('-')
            lower_str = re.search("[-+]?\d*\.\d+|\d+", lower).group()
            upper_str = re.search("[-+]?\d*\.\d+|\d+", upper).group()

            lower_num = float(lower_str)
            upper_num = float(upper_str)

            # leading_nulls = len(lower_str) - len(lower_str.lstrip('0'))
            length = len(lower_str)

            # must use integers in a loop, not floats
            if '.' in lower_str:
                decimals = len(lower_str.split('.')[1])
                multiplier = 10 * decimals
            else:
                multiplie r =1

            no_dec_lower = int(lower_nu m *multiplier)
            no_dec_upper = int((upper_num) * multiplier) + 1

            if '.' in lower_str:
                codes = [lower.replace(lower_str, str(nu m /multiplier).zfill(length)) for num in range(no_dec_lower, no_dec_upper)]
            else:
                codes = [lower.replace(lower_str, str(num).zfill(length)) for num in range(no_dec_lower, no_dec_upper)]


        else:
            codes = [expr]
        all_codes.extend(codes)
    return all_codes


def stringify_cols(df, cols):
    """
    Stringify some cols - useful since many methods erquire code column to be a string
    """

    for col in cols:
        df[col] = df[col].astype(str)
    return df


def format_codes(codes, merge=True):
    """
    Makes sure that the codes has the desired format: a dict with strings as
    keys (name) and a list of codes as values)

    Background: For several functions the user is allower to use strings
    when there is only one element in the list, and a list when there is
    no code replacement or aggregations, or a dict. To avoid (even more) mess
    the input is standardised as soon as possible in a function.

    Examples:
            codes = '4AB02'
            codes='4AB*'
            codes = ['4AB02', '4AB04', '4AC*']
            codes = ['4AB02', '4AB04']
            codes = {'tumor' : 'a4*', 'diabetes': ['d3*', 'd5-d9']}
            codes = 'S72*'
            codes = ['K50*', 'K51*']

            format_codes(codes, merge=False)

    TODO: test for correctness of input, not just reformat (is the key a str?)
    """
    codes = listify(codes)

    # treeatment of pure lists depends on whether special classes should be treated as one merged group or separate codes
    # exmple xounting of Z51* could mean count the total number of codes with Z51 OR a shorthand for saying "count all codes starting with Z51 separately
    # The option "merged, enables the user to switch between these two interpretations

    if isinstance(codes, list):
        if merge:
            codes = {'_'.join(codes): codes}
        else:
            codes = {code: [code] for code in codes}

    elif isinstance(codes, dict):
        new_codes = {}
        for name, codelist in codes.items():
            if isinstance(codelist, str):
                codelist = [codelist]
            new_codes[name] = codelist
        codes = new_codes

    return codes


def replace_codes(codes):
    """
    True if one or more  keys are different from its value in the dictionary

    Replacement of codes is unnecessary if all labels are the same as the codes

    """
    for name, code in codes.items():
        if name != code:
            return True  # ok, not beautiful ... may use break or any?
    return False


def expand_regex(expr, full_list):
    exprs = listify(expr)

    expanded = []

    if isinstance(full_list, pd.Series):
        pass
    elif isinstance(full_list, list):
        unique_series = pd.Series(full_list)
    elif isinstance(full_list, set):
        unique_series = pd.Series(list(full_list))

    for expr in exprs:
        match = unique_series.str.contains(expr)
        expanded.extend(unique_series[match])
    return expanded


def reverse_dict(dikt):
    new_dict = {}
    for name, codelist in dikt.items():
        codelist = listify(codelist)
        new_dict.update({code: name for code in codelist})
    return new_dict


def reverse_dict_old(dikt):
    """
    takes a dict and return a new dict with old values as key and old keys as values (in a list)

    example

    reverse_dict({'AB04a':'b', 'AB04b': 'b', 'AB04c':'b', 'CC04x': 'c'})

    will return
        {'b': ['AB04a', 'AB04b', 'AB04c'], 'c': 'CC04x'}
    """

    new_dikt = {}
    for k, v in dikt.items():
        if v in new_dikt:
            new_dikt[v].append(k)
        else:
            new_dikt[v] = [k]
    return new_dikt


def expand_replace(df, replace, cols, sep=None, strip=True):
    """
    Takes a dictionary of shorthand codes and replacements, and returns a dictionary with all the codes expanded

    Example:
        expand_replace(df=df, replace={'AB04*':'b'}, col='atc')

        May return
            {'AB04a':'b', 'AB04b': 'b', 'AB04c':'b'}

    """
    # may use regex instead, but this may also be slower to use later?
    cols = listify(cols)
    codes = list(replace.keys())

    codes = expand_codes(df=df, codes=codes, cols=cols, sep=None)

    unexpanded = {code: text for code, text in replace.items() if '*' in code}

    for starcode, text in unexpanded.items():

        startswith, endswith = starcode.split('*')

        # starting_codes  = ending_codes = start_and_end_codes = {}
        starting_codes = {}
        ending_codes = {}
        start_and_end_codes = {}
        # may be unnecessary to do this (and it may link the dictionaries in unwanted ways?)

        if startswith:
            starting_codes = {code: text for code in codes if code.startswith(startswith)}
        if endswith:
            ending_codes = {code: text for code in codes if code.endswith(endswith)}
        if startswith and endswith:
            start_and_end_codes = {starting_code: starting_code[x] for x in starting_code if x in ending_code}

        replace.update({**starting_codes, **ending_codes, **start_and_end_codes})

        del replace[starcode]
    return replace
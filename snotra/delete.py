import os
import uuid

from snotra.internal import to_df, fix_args, get_some_id, listify
from snotra.core import get_rows, expand_codes


def events2person(df, agg):
    """
    make person level data based on event data


    icd: cat [icdmain, icdbi]
    age: min age
    ibd: k50 or k51 in icd
    """
    pass


def persons_with(df,
                 codes,
                 cols,
                 pid='pid',
                 sep=None,
                 merge=True,
                 first_date=None,
                 last_date=None,
                 group=False,
                 _fix=True):
    """
    Determine whether people have received a code

    Args:
        codes (list or dict): codes to mark for
            codes to search for
                - if list: each code will represent a column
                - if dict: the codes in each item will be aggregated to one indicator
            cols (str or list of str): Column(s) with the codes
            pid (str): colum with the person identifier
            first_date (str): use only codes after a given date
                the string either represents a date (same for all individuals)
                or the name of a column with dates (may be different for different individuals)
            last_date (str): only use codes after a given date
                the string either represents a date (same for all individuals)
                or the name of a column with dates (may be different for different individuals)

    Returns:
        Series or Dataframe


    Examples:
        fracture = persons_with(df=df, codes='S72*', cols='icdmain')
        fracture = persons_with(df=df, codes={'frac':'S72*'}, cols='icdmain')

    Todo:
        - function may check if pid_index is unique, in which it does not have to aggregate
        - this may apply in general? functions that work on event data may then also work on person level data
        - allow user to input person level dataframe source?
    """
    sub = df

    if _fix:
        df, cols = to_df(df=df, cols=cols)
        codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep, merge=merge, group=group)
        rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
        sub = df[rows]

    df_persons = sub.groupby(pid)[cols].apply(lambda s: pd.unique(s.values.ravel()).tolist()).astype(str)

    # alternative approach, also good, and avoids creaintg personal dataframe
    # but ... regeis is fast since it stopw when it finds one true code!
    #    c=df.icdbi.str.split(', ', expand=True).to_sparse()
    #    c.isin(['S720', 'I10']).any(axis=1).any(level=0)

    persondf = pd.DataFrame(index=df[pid].unique().tolist())
    for name, codes in codes.items():
        codes_regex = '|'.join(codes)
        persondf[name] = df_persons.str.contains(codes_regex, na=False)

    return persondf


def get_uuid(df, codes, cols, uuid='uuid', sep=None):
    """
    Returns a set pids who have the given codes in the cols
    """

    uuids = get_some_id(df=df, codes=codes, some_id=uuid, sep=sep)

    return uuids


def make_uuid(df, name='uuid'):
    """
    Creates a list of uuids with the same length as the dataframe
    """

    uuids = [uuid.uuid4().hex for _ in range(len(df))]
    return uuids


def read_codebooks_csv(path=None,
                       codebooks=None,
                       language=None,
                       cols=None,
                       case=None,
                       dot=True,
                       merge=True,
                       sep=',',
                       info=None):
    """
    Reads csv files desribing medical codes return a dict with the codebooks
    or a merged dataframe with all codebooks

    Useful to translate from icd codes (and other codes) to text description

    Files should have one column with 'code' and one with 'text'


    parameters
    ----------
        codes
            'all' - returns one dictionaty with all codes
            'separate'   - returns a dict of dict (one for each code framework)
    example
        medcodes = read_code2text(codes='all')

    """
    # if no books specified, read all books that are discovered (or only last version ot it?)
    # code to find the books
    import os
    if not path:
        path = os.path.abspath(rr.__file__)

        path = path.replace('__init__.py', 'codebooks\\atc_2015_eng.csv')
        path = path.replace('__init__.py', 'codebooks\\icd10cm_order_2017.txt')

    atc = pd.read_csv(path, sep=';')
    atc.text = atc.text.str.strip()

    atc.to_csv(path, sep=';')

    from io import StringIO
    a = StringIO()

    with open(path, 'r') as file:
        in_memory_file = file.read()

    # enable file specific reding of files, keywords and values relvant for reading a particular file is in info[filname]
    # if nothing is specified, used same arguments for all files
    for book in codebooks:
        tmp_info[book] = {'path': path, 'usecols': cols, 'case': case, 'dot': dot, 'merge': merge, 'sep': sep}
        if book in info:
            for k, v in info.items():
                tmp_info[book][k] = v
        info[book] = tmp_info[book]

    paths = ['C:/Users/hmelberg/Google Drive/sandre/resources/health_pandas_codes',
             'C:/Users/hmelberg_adm/Google Drive/sandre/resources/health_pandas_codes',
             'C:/Users/sandresl/Google Drive/sandre/resources/health_pandas_codes',
             'C:/Users/sandresl_adm/Google Drive/sandre/resources/health_pandas_codes']

    for trypath in paths:
        if os.path.isdir(trypath):
            path = trypath
            break

    codebook = {}
    for book in codebooks:
        codebook[book] = pd.read_csv(f'{info[book][path]}/{book}_code2text.csv',
                                     encoding='latin-1')

        codebook['codebook'] = book

        if case:
            if case == 'upper':
                codebook['code'] = codebook['code'].str.upper()
            elif case == 'lower':
                codebook['code'] = codebook['code'].str.upper()

        if not keep_dot:
            codebook['code'] = codebook['code'].str.replace('.', '')

    if codes == 'all':
        code2textnew = {}
        for frame in codeframes:
            code2textnew.update(code2text[frame])
        code2text = code2textnew

    return code2text


def read_code2text(path='C:/Users/hmelberg/Google Drive/sandre/resources',
                   codes='all',
                   capitalized='both',
                   max_length=None):
    """
    Reads labels for medical codes from files, returns a dictionary code:text

    Useful to translate from icd codes (and other codes) to text description
    Reads from semicolonseperated csv files.
    Files should have one column with 'code' and one with 'text'


    parameters
    ----------
        codes
            'all' - returns one dictionaty with all codes
            'separate'   - returns a dict of dict (one for each code framework)
    example
        medcodes = read_code2text(codes='all')

    """

    paths = ['C:/Users/hmelberg/Google Drive/sandre/resources/health_pandas_codes',
             'C:/Users/hmelberg_adm/Google Drive/sandre/resources/health_pandas_codes',
             'C:/Users/sandresl/Google Drive/sandre/resources/health_pandas_codes',
             'C:/Users/sandresl_adm/Google Drive/sandre/resources/health_pandas_codes']

    for trypath in paths:
        if os.path.isdir(trypath):
            path = trypath
            break

    codeframes = 'icd atc nc drg'.split()

    code2text = {}

    for frame in codeframes:
        code2text_df = pd.read_csv(f'{path}/{frame}_code2text.csv',
                                   encoding='latin-1',
                                   sep=';')
        code2text_df.code = code2text_df.code.str.strip()

        code2text_dict = {**code2text_df[['code', 'text']].set_index('code').to_dict()['text']}

        if max_length:
            code2text_dict = {code: text[max_length] for code, text in code2text_dict.items()}

        if capitalized == 'both':
            capitalized = {str(code).upper(): text for code, text in code2text_dict.items()}
            uncapitalized = {str(code).lower(): text for code, text in code2text_dict.items()}
            code2text_dict = {**capitalized, **uncapitalized}

        code2text[frame] = code2text_dict

    if codes == 'all':
        code2textnew = {}
        for frame in codeframes:
            code2textnew.update(code2text[frame])
        code2text = code2textnew

    return code2text


def clean(df, rename=None,
          dates=None,
          delete=None,
          categorize=None,
          sort=None,
          index=None,
          cleaner=None):
    """
    fix dataframe before use (rename, sort, convert to dates)

    """
    if cleaner:
        clean_instructions = get_cleaner(cleaner)
        clean(df, cleaner=None, **clean_instructions)

    if rename:
        df =df.rename(columns=rename)

    if dates:
        for col in dates:
            df[col] =pd.to_datetime(df[col])

    if delete:
        for col in delete:
            del df[col]

    if sort:
        df =df.sort_values(sort)

    if index:
        df =df.set_index(index, drop=False)


    return df


def get_cleaner(name):
    if name == 'ibd_npr_2015':

        cleaner = {
            'rename': {'lopenr': 'pid', 'inn_mnd': 'start_date'},
            'dates': ['start_date'],
            'delete': ['Unnamed: 0'],
            'index': ['pid'],
            'sort': ['pid', 'start_date']
        }

    elif name == 'ibd_npr_2015':
        cleaner = {
            'rename': {'id': 'pid', 'innDato': 'start_date'},
            'dates': ['start_date'],
            'delete': ['Unnamed: 0'],
            'index': ['pid'],
            'sort': ['pid', 'start_date']
        }
    return cleaner


def events(df,
           codes,
           cols=['icd'],
           pid='pid',
           pre_query=None,
           post_query=None,
           sep=',',
           out='df',
           _fix=True):
    """
    Get all events for people who have a specific code/diagnosis

    parameters
        df: dataframe of all events for all patients
        codes: the codes that identify the patients of interest
                star notation is allowed
                example: icd codes for crohn's disease: K50*
        cols: the column(s) where the code(s) can be found
                star notation is allowed
        pid: the column with the id of the individuals
        sep: the seperator used between codes if multiple codes exist in a column
        out: if the output should be the (sub)dataframe or a list of ids

    example
        get all events only for those who are registered with an ibd diagosis

        ibd = events(df=df,
              codes=['K50*', 'K51*'],
              cols=['icd*'],
              pid='pid',
              sep=',')
    """

    if _fix:
        df, cols = to_df(df, cols)
        codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep)

    with_codes = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)

    pids = df[with_codes][pid].unique()

    if out == 'df':
        return df[df[pid].isin(pids)]

    elif out == 'pids':
        return pids
    else:
        print(f"Error: {out} is not a valid 'out' argument")
        return


def stringify(df,
              codes,
              cols,
              pid='pid',
              start_time='in_date',
              replace=None,
              end_time=None,
              sep=None,
              new_sep=None,
              single_value_columns=None,
              out='series'):
    codes = listify(codes)
    cols = listify(cols)

    single_cols = infer_single_value_columns(df=df,
                                             cols=cols,
                                             sep=sep,
                                             check_all=True)

    multiple_cols = [set(cols) - set(single_cols)]

    expanded_codes = expand_codes(df=df, cols=cols, codes=codes, sep=sep)
    # use expanded codes?, not codes as argument.
    # why? because a code may be in a compount col and not in single cols

    if single_cols:
        single_events = stringify_singles(df=df,
                                          codes=expanded_codes,
                                          cols=single_cols,
                                          pid=pid,
                                          start_time=start_time,
                                          replace=replace,
                                          end_time=end_time,
                                          out='df')
        all_events = single_events

    if multiple_cols:
        multiple_events = stringify_multiples(df=df,
                                              codes=expanded_codes,
                                              cols=multiple_cols,
                                              pid=pid,
                                              start_time=start_time,
                                              replace=replace,
                                              end_time=end_time,
                                              out='df')
        all_events = multiple_events

    if single_cols and multiple_cols:
        all_events = pd.concat([multiple_events, single_events])

    if out == 'series':
        events_by_id = all_events.sort_values([pid, start_time]).groupby(pid)['events'].sum()
        return events_by_id
    elif out == 'df':
        return all_events
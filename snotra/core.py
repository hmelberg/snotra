# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:17:33 2018

@author: hmelberg
"""


import numpy as np
import pandas as pd
import re
from itertools import zip_longest, chain


# structure
#
# core api (key functions/methods
# internal functions/methods
# stringify
#
#



##
## Core API
##


# %%

def incidence(df, codes=None, cols=None, sep=None, pid='pid',
              date='indate', min_events=1, within_period=None,
              groupby='cohort', update_cohort=True, codebook=None, _fix=True):
    """
    The number of new patients each year who have one or more of the codes


    Args:
        df (dataframe): Dataframe with unique personal ids, events, dates and medical codes

        codes (string, list or dict): The codes for the disease

        cols (string, list): Name of cols where codes are located

        pid (str): Name of column with the personal identification number

        date (str): Name of column with the dates for the events
            (the dtype of the column must be datetime)

        min_events (int): Number of events with the codes required for inclusion

        within_period (int): The number of events have to occur within a period
            (measured in days) before the person is included
            For instance: min_events=2 and within_period=365 means
            that a person has to have two events within one year
            in order to be included in the calculation of incidence.

        groupby (string): Name of column to group the results by. Typically
            year if you have data for many years (to get incidence each year)

        update_cohort (bool): The cohort a person is assigned to might change
            as the criteria (min_events and within_period) change. If update
            is True, the cohort will be automatically updated to reflect this.

            Example: A person may have her first event in 2011 and initially
            be assigned to the 2011 cohort, but have no other events for two
            years. Then in 2014 the person may have five events. If the
            criteria for inclusion is two event in one year, the updated cohort
            for this person will be 2014 and not the original 2011.

        Returns:
            Series (For instance, number of new patients every year)

        Examples:
            >>> df['cohort']=df.groupby('pid').start_date.min().dt.year

            >>> incidence(df, codes=['K50*', 'K51*'], cols=['icdmain', 'icdbi'], sep=',', pid='pid', date='start_date')

            >>> incidence(df, codes={'cd':'K50*', 'uc':'K51*'}, cols=['icdmain', 'icdbi'], sep=',', pid='pid', date='start_date')
    Note
        The function may seem overly complex, but sometimes it is most sensible to
        demand more than one codes wihtin a defined time period before a patients is
        counted as a true patient within a category.

        Note also that many of the arguments do not apply is an expression is used instead of a codelist

        todo:
            make cohort variable redundant ... already have date!
            make it possible to do monthly, quarterly etc incidence?
            enable expressions? ... but does it make sense? and does it fit? redundant?

    """
    sub = df
    incidence_list = []
    namelist = []

    if _fix:
        codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep, merge=True, group=False)
        rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
        sub = df[rows]

    for name, codelist in codes.items():
        rows = get_rows(df=sub, codes=codelist, cols=cols, sep=sep, _fix=False)
        sub = sub[rows]
        events = sub.groupby(pid).size()
        sub = sub[events >= min_events]

        if within_period:
            days_to_next = (sub.sort_values([pid, date]).groupby(pid)[date]
                            .diff(periods=-(min_events - 1))
                            .dt.days)

            # note: to be generalized?
            # not necessarily measure diff in days or cohort in years
            inside = (days_to_next >= -within_period)
            sub = sub[inside]

            if update_cohort:
                sub['cohort'] = sub.groupby(pid)[date].min().dt.year
            # may need ot update values in other rouping variables too?
            # for instance: disease group, if it is based on calc that changes
            # as obs are eliminated because of within_period requirements?
            # eg. disease group categorization based on majority of codes being x

        if groupby:
            incidence_df = sub.groupby(groupby)[pid].nunique()
        else:
            incidence_df = sub[pid].nunique()

        incidence_list.append(incidence_df)
        namelist.append(name)

    incidence_df = pd.concat(incidence_list, axis=1)
    incidence_df.columns = namelist

    if len(incidence_list) == 1:
        incidence_df = incidence_df.squeeze()

    return incidence_df


# %%
def make_cohort(df, codes=None, cols=None, sep=None, pid='pid',
                date='indate', min_events=1, within_period=None, _fix=True):
    """
    The first year with a given code given conditions


    Args:
        df (dataframe): Dataframe with events, dates and medical codes

        codes (string, list or dict): The codes for the disease

        cols (string, list): Name of columns where codes are located

        pid (str): Name of column with the personal identification number

        date (str): Name of column with the dates for the events
            (the dtype of the column must be datetime)

        min_events (int): Number of events with the codes required for inclusion

        within_period (int): The number of events have to ocurr within a period (measured in days) before the person is included
            For instance: min_events=2 and within_period=365 means
            that a person has to have two events within one year
            in order to be included in the calculation of incidence.


        Returns:
            Series (For instance, number of new patients every year)

        Examples:
            df['cohort']=df.groupby('pid').start_date.min().dt.year

            make_cohort(df, codes=['K50*', 'K51*'], cols=['icdmain', 'icdbi'], sep=',', pid='pid', date='start_date')

            make_cohort(df, codes={'cd':'K50*', 'uc':'K51*'}, cols=['icdmain', 'icdbi'], sep=',', pid='pid', date='start_date')

        todo:
            make cohort variable redundant ... already have date!
            make it possible to do monthly, quarterly etc incidence?
    """

    sub = df

    if _fix:
        codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep, merge=True, group=False)
        rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
        sub = df[rows]

    cohorts = []
    names = []

    for name, codelist in codes.items():
        rows = get_rows(df=sub, codes=codelist, cols=cols, sep=sep, _fix=False)
        sub2 = sub[rows]
        events = sub2.groupby(pid).size()
        pids = events.index[events >= min_events]
        sub2 = sub2[sub2[pid].isin(pids)]

        if within_period:
            days_to_next = (sub2.sort_values([pid, date]).groupby(pid)[date]
                            .diff(periods=-(min_events - 1))
                            .dt.days)

            inside = (days_to_next >= -within_period)
            sub2 = sub2[inside]

        cohort = sub2.groupby(pid)[date].min().dt.year
        cohorts.append(cohort)
        names.append(name)

    cohorts = pd.concat(cohorts, axis=1)
    cohorts.columns = names

    return cohorts


# %%
def sample_persons(df, pid='pid', n=None, frac=0.1):
    """
    Pick some randomly selected individuals and ALL their observations


    Args:
        n (int): number of individuals to sample
        frac (float): fraction of individuals to sample

    Returns:
        dataframe

    Note:
        Default: take a 10% sample

    Examples:
        sample_df=sample_persons(df, n=100)
        sample_df=sample_persons(df, n=100)


    """
    if isinstance(df, pd.Series):
        ids = df.index.nunique()
    else:
        ids = df[pid].unique()

    if not n:
        n = int(frac * len(ids))

    new_ids = np.random.choice(ids, size=n)

    if isinstance(df, pd.Series):
        new_sample = df[ids]
    else:
        new_sample = df[df.pid.isin(new_ids)]

    return new_sample


# %%
def unique_codes(df,
                 cols=None,
                 sep=None,
                 strip=True,
                 _sniffsep=True,
                 _fix=True):
    """
    Get set of unique values from column(s) with with multiple values in cells

    Args:
        df (dataframe)
        cols (str or list of str): columns with  content used to create unique values
        sep (str): if there are multiple codes in cells, the separator has to
             be specified. Example: sep =',' or sep=';'


    Returns:

    Note:
        - Each column may have multiple values (if sep is specified)
        - Star notation is allowed to describe columns: col='year*'
        - In large dataframes this function may take some time


    Examples:
        - get all unique codes in the icd column of df (separator will be sniffed)
        >>> df.icd.unique_codes()

        - all unique codes from all columns starting with 'icd'
        >>> df.unique_codes(cols='icd*')


    worry
        numeric columns/content
    """
    if _fix:
        df, cols = _to_df(df=df, cols=cols)
        cols = _fix_cols(df=df, cols=cols)

    unique_terms = set(pd.unique(df[cols].values.ravel('K')))
    # unique_terms={str(term) for term in unique_terms}

    if _sniffsep and not sep:
        sep = _sniff_sep(df, cols)

    if sep:
        compound_terms = {term for term in unique_terms if sep in str(term)}
        single_uniques = {term for term in unique_terms if sep not in str(term)}
        split_compounds = [term.split(sep) for term in compound_terms]
        split_uniques = {term.strip()
                         for sublist in split_compounds
                         for term in sublist
                         }
        unique_terms = single_uniques | split_uniques

    if strip:
        unique_terms = list({str(term).strip() for term in unique_terms})

    return unique_terms


# %%
def expand_codes(df=None,
                 codes=None,
                 cols=None,
                 sep=None,
                 codebook=None,
                 hyphen=True,
                 star=True,
                 colon=True,
                 regex=None,
                 del_dot=False,
                 case_sensitive=True,
                 exist=True,
                 merge=False,
                 group=False):
    """
    Expand list of codes with star, hyphen and colon notation to full codes

    Args:
        codes (str or list of str): a list of codes some of which may need to be expanded to full codes
        cols (str): a column with codes that are be used to build a codebook of all codes
            If codebook is specified the cols argument is not needed and ignored
        sep (str): seperator used if cells have multiple values
        codebook (list): User specified list of all possible or allowed codes
        expand_hyphen (bool, default: False): If True, codes with hyphens are not expanded
        expand_star (bool, default: False): If True, codes with start are not expanded


    Returns
        List of codes

    Examples

        >>> codes = expand_codes(df=df, codes=['H02*'], cols='atc')
        >>> codes = expand_codes(df=df, codes=['K25*'], cols='icdmain', sep=',', codebook=codebook)
        >>> codes = expand_codes(df=df, codes=['K25*-K28*'], cols='icdmain', sep=',', codebook=codebook, merge=True)


        >>> codes = expand_codes(df=df, codes=['K50*', 'K51*'], cols='icdmain', sep=',', codebook=codebook, merge=False, group=True)
        >>> codes = expand_codes(df=df, codes=['K50*', 'K51*'], cols='icdmain', sep=',', codebook=codebook, merge=False, group=False)

        >>> ulcer = ['K25*-K28*']
        >>> liver = 'B18* K70.0-K70.3 K70.9 K71.3-K71.5 K71.7 K73* K74* K76.0 K76.2-K76.4 K76.8 K76.9 Z94.4'.split()
        >>> diabetes =	['E10.0', 'E10.l', 'E10.6', 'E10.8', 'E10.9', 'E11.0', 'E11.1', 'E11.6', 'E11.8', 'E11.9', 'E12.0', 'E12.1', 'El2.6', 'E12.8', 'El2.9', 'E13.0', 'E13.1', 'E13.6', 'E13.8', 'E13.9', 'E14.0', 'E14.1', 'E14.6', 'E14.8', 'E14.9']
        >>> hemiplegia = ['G04.1', 'G11.4', 'G80.1', 'G80.2', 'G81*', 'G82*', 'G83.0-G83.4', 'G83.9']

        >>> expand_codes(df=df, codes=ulcer, cols='icdmain', sep=',')

    Note:
        Only codes that actually exist in the cols or the codebook are returned

    """

    codes = _listify(codes)

    if codebook:
        unique_words = set(codebook)
    else:
        cols = _listify(cols)
        cols = _expand_cols(df=df, cols=cols)
        unique_words = set(unique_codes(df=df, cols=cols, sep=sep))

    # if input is not a list of codes, but a dict with categories and codes,
    # then expand each category separately and return the whole dict with
    # expanded codes

    if isinstance(codes, list):
        alist = True
    else:
        alist = False

    codes = _format_codes(codes=codes, merge=merge)

    all_new_codes = {}

    for name, codelist in codes.items():

        # for instance in icd-10 some use codes with dots, some without
        if del_dot:
            # unique_words = {word.replace('.', '') for word in unique_words}
            codelist = [code.replace('.', '') for code in codelist]

        if not case_sensitive:
            # unique_words = {word.lower() for word in unique_words}
            codelist = [code.lower() for code in codelist] + [code.upper() for code in codelist]

        # expand hyphens codes but keep only those that are in the cols or the codebook
        if hyphen:
            codelist = expand_hyphen(codelist)

        # expand only codes with star notation and when expand_stars is turned on
        if star:
            codelist = expand_star(codelist, full_list=unique_words)

        # regex can be used, but may be complex if combined with other
        # maybe introduce a notation for regex inside the codebook?, like re:
        # (so the regex is not done on all codes)
        if regex:
            codelist = _expand_regex(codelist, full_list=unique_words)

        # eliminate codes that have been created by the expansion, but that do not
        # exist in the data. For instance, the hyphen expansion may create this.
        if exist:
            match = set(codelist) & unique_words
            codelist = list(match)

        all_new_codes[name] = codelist

        # Change dictionary depending on whether the user wants codes with
        # special notations (star, hyphen, colon, eg. K51* to stay as a
        # separate group or be split in its individual subcodes

    if (not group) and (not merge):
        new_all_new_codes = {}
        for name, codelist in all_new_codes.items():
            if ('*' in name) or ('-' in name) or (':' in name):
                for code in codelist:
                    new_all_new_codes[code] = [code]
            else:
                new_all_new_codes[name] = codelist
        all_new_codes = new_all_new_codes

    if merge:
        pass
        # all_new_codes=list(all_new_codes.values())[0]

    return all_new_codes


# %%
def get_rows(df,
             codes,
             cols,
             sep=None,
             codebook=None,
             _fix=True):
    """
    Returns a boolean array that is true for the rows where column(s) contain the code(s)

    Args:
        df (dataframe): Dataframe with events, dates and medical codes

        codes (string, list or dict): The codes for the disease

        cols (string, list): Name of columns where codes are located

        pid (str): Name of column with the personal identification number

        codebook (list): User specified list of all possible or allowed codes


    Examples:
        - get all drg codes starting with 'D':

        >>> d_codes = get_rows(df=df, codes='D*', cols=['drg'])
    """
    # if an expression is used as input
    if isinstance(codes, str) and codes.count(' ') > 1:
        rows = use_expression(df, codes, cols=cols, sep=sep, out='rows', codebook=codebook)

    # if a list of codes is used as input
    else:
        if _fix:
            df, cols = _to_df(df)
            cols = _fix_cols(df=df, cols=cols)
            codes = _fix_codes(df=df, codes=codes, cols=cols, sep=sep)

        codes = _listify(codes)

        allcodes = _get_allcodes(codes)

        if len(allcodes) == 0:  # if no relevant codes --> a column with all false
            rows = np.full(len(df), False)
        elif sep:
            allcodes_regex = '|'.join(allcodes)
            rows = np.full(len(df), False)
            for col in cols:
                new_rows = df[col].astype(str).str.contains(allcodes_regex, na=False).values
                rows = rows | new_rows
        # if single value cells only
        else:
            rows = df[cols].isin(allcodes).any(axis=1).values

    return rows


# %%
def get_pids(df,
             codes,
             cols=None,
             pid='pid',
             sep=None,
             codebook=None):
    """
    Returns the set pids who have the given codes in the cols

    Args:
        df (dataframe): Dataframe with events

        codes (string, list or dict): The codes for the disease

        cols (string, list): Name of columns where codes are located

        pid (str): Name of column with the personal identification number

        codebook (list): User specified list of all possible or allowed codes

    Examples:

        >>> get pids for all individuals who have icd codes starting with 'C509':

        >>> c509 = get_pids(df=df, codes='C509', cols=['icdmain', 'icdbi'], pid='pid')
    """

    # if an expression instead of a codelist is used as input
    if isinstance(codes, str) and codes.count(' ') > 1:
        selection = use_expression(df, codes, cols=cols, sep=sep, out='persons', codebook=codebook, pid=pid)
        pids = df[selection][pid]
    else:
        pids = _get_some_id(df=df, codes=codes, some_id=pid, sep=sep)
    return pids


# %%
def select_persons(df,
                   codes,
                   cols=None,
                   pid='pid',
                   sep=None,
                   codebook=None,
                   _fix=True):
    """
    Returns a dataframe with all events for people who have the given codes

    Args:
        df (dataframe): Dataframe with events

        codes (string, list or dict): The codes for the disease

        cols (string, list): Name of columns where codes are located

        pid (str): Name of column with the personal identification number

        codebook (list): User specified list of all possible or allowed codes

    Examples:
        >>> df.select_persons(codes='K51 and not K50')
        >>> df.select_persons(codes='(K50 in: icd or K51 in: icd) and 4AB04 in atc')
        >>> df.select_persons(codes='C509', cols=['icdmain', 'icdbi'])

    """
    if _fix:
        df, cols = _to_df(df=df, cols=cols)
        if cols:
            cols = _expand_cols(df=df, cols=cols)
            if not sep:
                sep = _sniff_sep(df=df, cols=cols)

    # if an expression is used as input (eg 'K50 and not K51')
    if isinstance(codes, str) and codes.count(' ') > 1:
        selection = use_expression(df, codes, cols=cols, sep=sep, out='persons', codebook=codebook, pid=pid)
        pids = df[selection][pid].values
    # if an list of codes - ['K50', 'K51'] is used as input
    else:
        pids = _get_some_id(df=df, codes=codes, some_id=pid, sep=sep)

    df = df[df[pid].isin(pids)]

    return df


# %%
def count_persons(df, codes=None, cols=None, pid='pid', sep=None,
                  normalize=False, dropna=True, group=False, merge=False,
                  groupby=None, codebook=None, _fix=True):
    """
    Counts number of individuals who are registered with given codes

    Allows counting across multiple columns and multiple codes in the same
    cells. For instance, there may be 10 diagnostic codes for one event (in
    separate columns) and in some of the columns there may be more than one
    diagnostic code (comma separated) and patient may have several such events
    in the dataframe.

    args:
        codes (str, list or dict): Codes to be counted. Star and hyphen
        notations are allowed. A dict can be used as input to merge codes
        into larger categories before counting. The key is the name of
        the category ('diabetes') and the value is a list of codes.

            Examples:
                codes="4ABA2"
                codes="4AB*"
                codes=['4AB2A', '4AB4A']
                codes = {'diabetes' = ['34r32f', '3a*']}

        cols (str or list): The column(s) with the codes. Star and colon
        notation allowed.
            Examples:
                cols = 'icdmain'
                cols = ['icdmain', 'icdside']
                # all columns starting with 'icd'
                cols = ['icd*'] # all columns starting with 'icd'
                # all columns including and between icd1 and icd10
                cols = ['icd1:icd10']

        pid (str): Column name of the personal identifier
        sep (str): The code seperator (if multiple codes in the same cells)
        normalize (bool, default: False): If True, converts to pct
        dropna (bool, default True): Include counts of how many did not get
            any of the specified codes

    Examples
        >>> df.atc.count_persons(codes='4AB04')

        >>> df.atc.count_persons(codes='4AB04', dropna=False, normalize=True)

        >>> df.atc.count_persons(codes=['4AB*', '4AC*'])

        >>> df.atc.count_persons(codes=['4AB*', '4AC*'], group=True)
        >>> df.atc.count_persons(codes=['4AB*', '4AC*'], group=True, merge=True)

        >>> df.count_persons(codes={'adaliamumab':'4AB04'}, cols='ncmp', sep=',', pid='pid')

        >>> df.count_persons(codes='4AB04', cols='ncmp', groupby=['disease', 'cohort'])

        >>> df.groupby(['disease', 'cohort']).apply(count_persons, cols='ncmp', codes='4AB04', sep=',')

    """

    sub = df
    sub, cols = _to_df(df=sub, cols=cols)
    cols = _expand_cols(df=sub, cols=cols)

    # if an expression instead of a codelist is used as input
    if isinstance(codes, str) and codes.count(' ') > 1:
        persons = use_expression(df=sub, expr=codes, cols=cols, sep=sep, out='persons', codebook=codebook, pid=pid)
        if normalize:
            counted = persons.sum() / len(persons)
        else:
            counted = persons.sum()


    # if codes is a codelist (not an expression)
    else:
        if _fix:
            # expands and formats columns and codes input
            codes, cols, allcodes, sep = _fix_args(df=sub, codes=codes, cols=cols, sep=sep, group=group, merge=merge)
            rows = get_rows(df=sub, codes=allcodes, cols=cols, sep=sep, _fix=False)
            if not dropna:
                persons = df[pid].nunique()
            sub = sub[rows].set_index(pid, drop=False)

        # make a df with the extracted codes
        code_df = extract_codes(df=sub, codes=codes, cols=cols, sep=sep, _fix=False, series=False)

        labels = list(code_df.columns)

        counted = pd.Series(index=labels)

        if groupby:
            code_df = code_df.any(level=0)
            sub_plevel = sub.groupby(pid)[groupby].first()
            code_df = pd.concat([code_df, sub_plevel], axis=1)  # outer vs inner problem?

            code_df = code_df.set_index(groupby)
            counted = code_df.groupby(groupby).sum()

        else:
            for label in labels:
                counted[label] = code_df[code_df[label]].index.nunique()

        if not dropna:
            with_codes = code_df.any(axis=1).any(level=0).sum()  # surprisingly time consuming?
            nan_persons = persons - with_codes
            counted['NaN'] = nan_persons

        if normalize:
            counted = counted / counted.sum()
        else:
            counted = counted.astype(int)

        if len(counted) == 1:
            counted = counted.values[0]

    return counted


# %%
def use_expression(df, expr, cols=None, sep=None, out='rows', raw=False, regex=False, logic=True, codebook=None,
                   pid='pid', _fix=True):
    """
    Returns a series with the index and the boolean value for whether the expression is true or not

    Args:
        df (dataframe): Dataframe with events

        expr (string, list or dict): The codes for the disease

        cols (string, list): Name of columns where codes are located

        pid (str): Name of column with the personal identification number

        codebook (list): User specified list of all possible or allowed codes

    Examples:

        expr = 'K52* and not (K50 or K51)'
        expr = 'K52* in icd and not (K50 or K51) in ncmp'

        expr = 'K52* in icd and not 4AB04 in ncmp or atc'

        expr = 'in icdmain or icdbi: (k50 or k51) and in ncmp: 4AB04 and 4AB02)'

        expr = '(K50 in:icdmain1,icdmain2 or K51 in:icdmain1,icdmain2) and (4AB04 in:ncmp or 4AB02 in:ncmp)'

        expr = 'k50 or k51 in icdmain or icdbi and (4AB04 and 4AB02) in ncmp'
        expr = 'k50==icdmain
        expr = 'K51* and 4AB04'

        expr='4AB02 in:ncmp and not 4AB02 in:ncsp'

        expands   ... columns connected by logical operators should get expressions in front
        abode_expr = 'K51* in icdmain and K51* in icdbi ...'

        expr = 'icd==K50 and age==40'

    Note:
        Coding logic:
        - pick out every code expression and expand it? (only star expansion, hyphen would work too? key is to avoud codebook or full lookup )
        - get codebook (if needed)
        - use extract code on each expression
        - execute logic
        - return series (bool)

    Examples:
    >>> get_rows_expression(df=npr, expr=expr, cols='icd', sep=',', out='rows'):


    """
    cols = _listify(cols)
    skipwords = {'and', 'or', 'not'}
    # split_at = {'=', '>', '<'} # what about '>=' '!' '!=' etc
    # well, just use query to deal with all this. eg if want to examine age>40
    # also additional groupbys ... need to think harder/rewrite/clean up, for now: allow some slightly more complex stuff

    if _fix:
        df, cols = _to_df(df, cols)
        if not sep and not ' in:' in expr:
            sep = _sniff_sep(df=df, cols=cols)
        df = df.set_index(pid, drop=False)  # maybe use "if index.name !='pid_index' since indexing may take time

    # one procedure for expressions with multiple columns (using " in:")
    # another for expressions within single columns (no " in: ")
    if " in:" in expr:
        expr = expr.replace(': ', ':').replace(', ', ',')
        words = expr.split()
        words = [word.strip('(').strip(')') for word in words if word not in skipwords]

        word_cols = list(zip(words[0::2], words[1::2]))
        # BUG SOLVED same code in two cols will create problems with dict, better use a
        # list of tuples, not a dict? YES or give name as combination of code and col?
        # remember, automatic unpacking of tuples so naming works

        del_cols = [col for word, col in word_cols]

        word_cols = [(word, col.replace('in:', '').split(',')) for word, col in word_cols]

        coldf = pd.DataFrame(index=df.index)

        # if no global codebook is specified:
        # this create separate codebook for each column(s) condition
        # background in case some codes overlap (may not be necessary)
        # potential problem: empty codebooks?
        new_codebook = {}
        for word, cols in word_cols:
            sep = _sniff_sep(df=df, cols=cols)
            name = "".join(cols)
            if not codebook:
                if name not in new_codebook:
                    new_codebook[name] = unique_codes(df=df, cols=cols, sep=sep)
            else:
                new_codebook[name] = codebook

        for n, (word, cols) in enumerate(word_cols):
            # added n to number conditions and avoid name conflicts if same condition (but different column)
            worddict = {'___' + word.replace('*', '___') + f'_{n}': [word]}

            # allow star etc notation in col also?
            # cols=_expand_cols(df=df, cols=cols)
            codes = expand_codes(df=df, codes=worddict, cols=cols, sep=sep, codebook=new_codebook[''.join(cols)],
                                 merge=True, group=True)

            for name, codelist in codes.items():  # works, but really only one item in the dict here
                coldf[name] = get_rows(df=df, codes=codelist, cols=cols, sep=sep, _fix=False)

        evalexpr = expr
        for col in del_cols:
            evalexpr = evalexpr.replace(col, '')

        words = [word for word, col in word_cols]

        for n, word in enumerate(words):
            word = word.strip()
            evalexpr = evalexpr.replace(word + ' ', f'___{word}_{n}==1', 1)

        evalexpr = evalexpr.replace('*', '___')

        coldf = coldf.fillna(False)

        # if search in same columns for all conditions
    else:
        # find all words
        words = expr.split()
        words = {word.strip('(').strip(')') for word in words}

        words = set(words)

        if skipwords:
            words = words - skipwords

        if not codebook:
            codebook = unique_codes(df=df, cols=cols, sep=sep)

        # must avoid * since eval does not like in in var names, replace * with three ___
        # same with column names starting with digit, sp add three (___) to all words
        worddict = {'___' + word.replace('*', '___'): [word] for word in words}
        coldf = pd.DataFrame(index=df.index)

        # allow star etc notation in col also?
        # cols=_expand_cols(df=df, cols=cols)
        codes = expand_codes(df=df, codes=worddict, cols=cols, sep=sep, codebook=codebook)

        for name, codelist in codes.items():
            coldf[name] = get_rows(df=df, codes=codelist, cols=cols, sep=sep, _fix=False)

        evalexpr = expr

        for word in words:
            word = word.strip()
            evalexpr = evalexpr.replace(word, f'___{word}==1')

        evalexpr = evalexpr.replace('*', '___')

        coldf = coldf.fillna(False)

    # if the expression be evaluated at row level or person level
    if out == 'persons':
        # cold=coldf.groupby(pid).any()
        coldf = coldf.any(level=0)

    expr_evaluated = coldf.eval(evalexpr)

    return expr_evaluated


# %%
def search_text(df, text, cols='text', select=None, raw=False, regex=False, logic=True, has_underscore=False):
    """
    Searches column(s) in a dataframe for occurrences of words or phrases

    Can be used to search for occurrences codes that are associated with certain words in the text description.

    Args:
        df (dataframe or series) : The dataframe or series with columns to be searched
        cols (str or list of str): The columns to be searched.
        select (str): row selector for the dataframe. Example: "codetype:'icd' year:2011"
        raw (bool): if True, searches for the raw textstring without modifications
        regex (bool): If True, use regex when searching
        logic (bool): If True, use logical operators in the search (and, or, not in the string)
        has_underscore (bool): Set to true if the text contains underscores that are important (If set to true, it becomes to search for phrases i.e. two or more words right after each other

    Returns:
        A dataframe with the rows that satisfy the search conditions (contain/not contain the words/phrases the user specified)
        Often: The codes where the description contain certain words

    Examples:
        icd.search_text('diabetes')
        icd.search_text('diabetes and heart')
        icd.search_text('cancer and not (breast or prostate)')


     Strcture
        0. select rows (using query)
        0.5 Identify phrases and substitute space in phrases with underscores to the phrase is considered to be one word
        1. find all whole words
        2. select search methode depending on input (search for raw text, search using regex, search using logical operators etc)
        3. replace all terms with '{term}==1, but not replace ogical operator terms
        4. create str.contains bool col for every term
        5. run pd eval
    """

    cols = _listify(cols)

    df, cols = _to_df(df, cols)

    # make it a df with text as col if input is a series, or it is used as a method used on a series object
    if isinstance(df, pd.Series):
        df = df.to_frame()
        df.columns = ['text']
        # and give error if select is specified?

    # restrict search to select relevant rows (for instance only icd codes in 2016)
    # useful if you have a big dataframe with all codes from different codebooks and years
    if select:
        select.replace(':', '==')
        df = df.query(select)

    ## find all whole words used in the text

    # first: words within quotation marks (within the string) are to be considered "one word"
    # to make this happen, replace space in text within strings with underscores
    # then the regex will consider it one word - and we reintroduce spaces in texts with stuff with underscore when searching
    if not has_underscore:
        phrases = re.findall(r'\"(.+?)\"', text)
        for phrase in phrases:
            text = text.replace(phrase, phrase.replace(' ', '_'))

    # find all words
    word_pattern = r'\w+'
    words = set(re.findall(word_pattern, text))
    skipwords = {'and', 'or', 'not'}
    if skipwords:
        words = words - skipwords
    rows_all_cols = len(df) * [False]  # nb common mistake

    # only need to use logical operator transformation if the string has and, or or not in it
    if skipwords & words:
        logic = True

    # conduct search: either just the raw tet, the regex, or the one with logical operators (and, or not)
    if raw:
        for col in cols:
            rows_with_word = df[col].str_contains(text, na=False, regex=False)
            rows_all_cols = rows_all_cols | rows_with_word  # doublecheck!
    elif regex:
        for col in cols:
            rows_with_word = df[col].str_contains(text, na=False, regex=True)
            rows_all_cols = rows_all_cols | rows_with_word  # doublecheck!

    elif logic:
        for col in cols:
            for word in words:
                name = word
                # words with underscores are phrases and underscores must be removed before searching
                if ('_' in word) and (has_underscore): word = word.replace('_', ' ')
                df[name] = df[col].str.contains(word, na=False)
            all_words = re.sub(r'(\w+)', r'\1==1', text)
            # inelegant, but works
            for word in skipwords:
                all_words = all_words.replace(f'{word}==1', word)
            rows_with_word = df.eval(all_words)  # does the return include index?
        rows_all_cols = rows_all_cols | rows_with_word  # doublecheck!
    else:
        for col in cols:
            rows_with_word = df[col].str_contains(text, na=False, regex=False)
            rows_all_cols = rows_all_cols | rows_with_word  # doublecheck!

    df = df[rows_all_cols]
    return df


# %%
def first_event(df, codes, cols=None, pid='pid', date='in_date', sep=None, codebook=None):
    """
    Returns time of the first observation with any of the codes for the person



    Args:
        df (dataframe): Dataframe with events

        codes (string, list or dict): The codes for the disease

        cols (string, list): Name of columns where codes are located

        sep (string, default: None): Separator between codes in same cell (if exist)
            (If None, the function will infer the separator)

        pid (str, default: 'pid'): Name of column with the personal identification number

        codebook (list): User specified list of all possible or allowed codes

        Returns:
            Pandas series


        Examples:
            first_event(id_col = 'pid', date_col='diagnose_date', groupby = ['disease'], return_as='dict')

            df['cohort'] = df.first_event(id_col = 'pid', date_col='diagnose_date', return_as='dict')

            Date of first registered event with an ibd code for each individual:
                first_event(df=df, codes=['k50*','k51*'], cols='icd', date='date')
    """
    codes = _listify(codes)
    cols = _listify(cols)
    sep = _sniff_sep(df=df, cols=cols)
    cols = _expand_cols(df=df, cols=cols)
    codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep, codebook=codebook)

    rows_with_codes = get_rows(df=df, codes=codes, cols=cols, sep=sep, _fix=False)
    subdf = df[rows_with_codes]

    # groupby.extent(pid)
    first_date = subdf[[pid, date]].groupby(pid, sort=False)[date].min()

    return first_date


# %%

def extract_codes(df, codes, cols=None, sep=None, new_sep=',', na_rep='',
                  prefix=None, merge=False, out='bool', _fix=True, series=True, group=False):
    """
    Produce one or more columns with only selected codes

    Args:
        df (dataframe): Dataframe with events

        codes (string, list or dict): The codes for the disease

        cols (string, list): Name of columns where codes are located

        sep (string, default: None): Separator between codes in same cell (if exist)
            (If None, the function will infer the separator)

        pid (str, default: 'pid'): Name of column with the personal identification number

        codebook (list): User specified list of all possible or allowed codes

        merge (bool): Content of all columns is merged to one series # only if out='text'?

        group (bool): Star an othernotation remain a single group, not split into indivdiual codes

        out (string, ['text', 'category', 'bool' or 'int']): Datatype of output column(s)

    Notes:
        Can produce a set of dummy columns for codes and code groups.
        Can also produce a merged column with only extracted codes.
        Accept star notation.
        Also accepts both single value columns and columns with compund codes and seperators


    Example:
    to create three dummy columns, based on codes in icdmain column:

    >>> extract_codes(df=df,
    >>>          codes={'fracture' : 'S72*', 'cd': 'K50*', 'uc': 'K51*'},
    >>>          cols=['icdmain', 'icdbi'],
    >>>          merge=False,
    >>>          out='text')

    nb: problem with extract rows if dataframe is empty (none of the requested codes)
    """
    if _fix:
        df, cols = _to_df(df=df, cols=cols)
        codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep, group=group, merge=merge)

    subset = pd.DataFrame(index=df.index)

    for k, v in codes.items():
        rows = get_rows(df=df, codes=v, cols=cols, sep=sep, _fix=False)
        if out == 'bool':
            subset[k] = rows
        elif out == 'int':
            subset[k] = rows.astype(int)
        elif out == 'category':
            subset.loc[rows, k] = k
            subset[k] = subset[k].astype('category')
        else:
            subset[k] = na_rep
            subset.loc[rows, k] = k

    if (merge) and (out == 'bool'):
        subset = subset.astype(int).astype(str)

    new_codes = list(subset.columns)

    if (merge) and (len(codes) > 1):
        headline = ', '.join(new_codes)
        merged = subset.iloc[:, 0].str.cat(subset.iloc[:, 1:].T.values, sep=new_sep, na_rep=na_rep)
        merged = merged.str.strip(',')
        subset = merged
        subset.name = headline
        if out == 'category':
            subset = subset.astype('category')

    # return a series if only one code is asked for (and also if merged?)
    if series and (len(codes) == 1):
        subset = subset.squeeze()

    return subset


# %%
def years_in_row(df, year_col='year', groupby=None, info_bank=None, out='pct'):
    """
    average years in row patients are observed, for different start years

    Example:
        >>> years = years_in_row(df, year_col='aar', out='pct')
    """
    years = df[year_col].unique()
    min_year = min(years)
    max_year = end_year = max(years)

    pids = df.groupby(year_col)['lopenr'].unique().apply(set).to_dict()
    remains = {}

    for start_year in range(min_year, max_year):
        remains[(start_year, start_year)] = pids[start_year]
        for end_year in range(start_year + 1, max_year):
            remains[(start_year, end_year)] = remains[(start_year, end_year - 1)] & pids[end_year]

    if out == 'pct':
        print('pct, hello')
        for start_year in range(min_year, max_year):
            start_n = len(remains[(start_year, start_year)])
            remains[(start_year, start_year)] = 1
            for end_year in range(start_year + 1, max_year):
                remains[(start_year, end_year)] = len(remains[(start_year, end_year)]) / start_n

    return remains


# %%
def years_apart(df, pid='pid', year='year'):
    """
    pct of patients with observations that are x years apart

    Example:
        >>> years_apart(df=df[ibd])

    """
    a = df.groupby(pid)[year].unique()
    b = (a[a.apply(len) > 1]
         .apply(sorted)
         .apply(np.diff)
         .sub(1)
         .apply(max)
         .value_counts()
         .div(len(a))
         .sort_index()
         .mul(100)

         )
    return b


# %%
def label(df, labels=None, read=True, path=None):
    """
    Translate codes in index to text labels based on content of the dict labels

    Args:
        labels (dict): dictionary from codes to text
        read (bool): read and use internal dictionary if no dictionary is provided
    """
    if not labels:
        # making life easier for myself
        try:
            labels = read_code2text()
        except:
            labels = read_code2text(path)
    df = df.rename(index=labels)
    return df


# %%

def count_codes(df, codes=None, cols=None, sep=None, normalize=False,
                ascending=False, _fix=True, merge=False, group=False, dropna=True):
    """
    Count frequency of values in multiple columns or columns with seperators

    Args:
        codes (str, list of str, dict): codes to be counted
        cols (str or list of str): columns where codes are
        sep (str): separator if multiple codes in cells
        merge (bool): If False, each code wil be counted separately
            If True (default), each code with special notation will be counted together
        strip (bool): strip spacec bore and after code before counting
        ignore_case (bool): determine if codes with same characters,
            but different cases should be the same
        normalize (bool): If True, outputs percentages and not absolute numbers

    allows
        - star notation in codes and columns
        - values in cells with multiple valules can be separated (if sep is defined)
        - replacement and aggregation to larger groups (when code is a dict)

    example
    To count the number of stereoid events (codes starting with H2) and use of
    antibiotics (codes starting with xx) in all columns where the column names
    starts with "atc":

    count_codes(df=df,
                 codes={'stereoids' : 'H2*', 'antibiotics' : =['AI3*']},
                 cols='atc*',
                 sep=',')

    more examples
    -------------

    df.count_codes(codes='K51*', cols='icd', sep=',')
    count_codes(df, codes='K51*', cols='icdm', sep=',', group=True)
    count_codes(df, codes='Z51*', cols=['icd', 'icdbi'], sep=',')
    count_codes(df, codes='Z51*', cols=['icdmain', 'icdbi'], sep=',', group=True)
    count_codes(df, codes={'radiation': 'Z51*'}, cols=['icd'], sep=',')
    count_codes(df, codes={'radiation': 'Z51*'}, cols=['icdmain', 'icdbi'], sep=',')
    count_codes(df, codes={'crohns': 'K50*', 'uc':'K51*'}, cols=['icdmain', 'icdbi'], sep=',')
    count_codes(df, codes={'crohns': 'K50*', 'uc':'K51*'}, cols=['icdmain', 'icdbi'], sep=',', dropna=True)
    count_codes(df, codes={'crohns': 'K50*', 'uc':'K51*'}, cols=['icdmain', 'icdbi'], sep=',', dropna=False)
    count_codes(df, codes={'crohns': 'K50*', 'uc':'K51*'}, cols=['icdmain', 'icdbi'], sep=',', dropna=False, group=False)
    count_codes(df, codes=['K50*', 'K51*'], cols=['icd'], sep=',', dropna=False, group=True, merge=False)
    count_codes(df, codes=['K50*', 'K51*'], cols=['icdmain', 'icdbi'], sep=',', dropna=False, group=False, merge=False)
    count_codes(df, codes=['K50*', 'K51*'], cols=['icdmain', 'icdbi'], sep=',', dropna=False, group=False, merge=True)
    count_codes(df, codes=['K50*', 'K51*'], cols=['icdmain', 'icdbi'], sep=',', dropna=True, group=True, merge=True)
    #group fasle, merge true, for list = wrong ...

    count_codes(df, codes=['K50*', 'K51*'], cols=['icdmain', 'icdbi'], sep=',', dropna=True, group=False, merge=False)


    """
    # count all if codes is codes is not specified
    # use all columns if col is not specified
    sub = df

    if _fix:
        sub, cols = _to_df(df=sub, cols=cols)
        cols = _fix_cols(df=sub, cols=cols)
        if not sep:
            sep = _sniff_sep(df=sub, cols=cols)

        if codes:
            codes = _format_codes(codes=codes, merge=merge)
            codes = expand_codes(df=sub, codes=codes, cols=cols, sep=sep, merge=merge, group=group)
            allcodes = _get_allcodes(codes)
            if dropna:
                rows = get_rows(df=sub, codes=allcodes, cols=cols, sep=sep, _fix=False)
                sub = sub[rows]

    if sep:
        count_df = [sub[col].str
                        .split(sep, expand=True)
                        .apply(lambda x: x.str.strip())
                        .to_sparse()
                        .apply(pd.Series.value_counts)
                        .sum(axis=1)
                    for col in cols]

        count_df = pd.DataFrame(count_df).T
        code_count = count_df.sum(axis=1)
    else:
        code_count = sub[cols].apply(pd.Series.value_counts).sum(axis=1)

    if codes:
        allcodes = _get_allcodes(codes)
        not_included_n = code_count[~code_count.isin(allcodes)].sum()
        code_count = code_count[allcodes]
        if not dropna:
            code_count['na'] = not_included_n

    if isinstance(codes, dict):
        code_count = code_count.rename(index=_reverse_dict(codes)).sum(level=0)

    if normalize:
        code_n = code_count.sum()
        code_count = code_count / code_n
    else:
        code_count = code_count.astype(int)

    if ascending:
        code_count = code_count.sort_values(ascending=True)
    else:
        code_count = code_count.sort_values(ascending=False)

    return code_count


# %%
def lookup_codes(dikt, codes):
    """
    returns those elements in a dict where key starts with the expressions listed in codes

    todo: more complicated star notations: starts with, contains, endswith
    lookup(medcodes, 'L04*')

    """

    codes = _listify(codes)
    codes = [code.upper().strip('*') for code in codes]
    codes = tuple(codes)

    selected_codes = {k: v for k, v in dikt.items() if str(k).upper().startswith(codes)}
    return selected_codes


# %%
def get_codes(dikt, text):
    """
    returns those elements in a dict where value contains the expressions listed in codes

    todo: more complicated star notations: starts with, contains, endswith
    alterative name: find_codes? get_codes?

    example
    get all codes that have "steroid" in the explanatory text

        get_codes(medcodes, 'steroid*')

    """

    text = _listify(text)
    text = [txt.upper().strip('*') for txt in text]
    # codes = " ".join(codes)

    selected_codes = {k: v for k, v in dikt.items() if any(txt in str(v).upper() for txt in text)}

    return selected_codes


# %%
def sankey_format(df, labels=None, normalize=False, dropna=False, threshold=0.01):
    """
    Format the dataframe so it is easy fo create a holoviews sankey figure

    labels=dict(bio_codes.values())
    import holoviews as hv
    hv.Sankey(t1).options(label_position='left')
    hv.extension('bokeh')
    t4=t1.copy()

    """
    a = df
    a = a.apply(lambda row: ' '.join(row))
    a = a.str.split(expand=True)

    for col in a.columns:
        a[col] = a[col] + ' (' + str(col + 1) + ')'

    if not dropna:
        a = a.fillna(f'No new')

    all_counts = {}
    for col in range(len(a.columns))[1:]:
        counts = a.groupby(a[col - 1])[col].value_counts(normalize=normalize)
        if normalize:
            counts = counts.mul(100).astype(int).fillna(0)

        counts.name = 'value'
        # counts = counts.rename(index=labels).reset_index()
        counts = counts.reset_index()
        counts.columns = ['source', 'target', 'value']

        all_counts[col] = counts
    t1 = pd.concat(all_counts, ignore_index=True)
    t1 = t1[t1.source != 'No new']

    # a.groupby(1)[2].value_counts()
    return t1


# %%
def charlson(df, cols='icd', pid='pid', age='age', sep=None, dot_notation=False):
    """
    Calculates the Charlson comorbidity indec (one year mortality index)

    Wiki: The Charlson comorbidity index predicts the one-year mortality
    for a patient who may have a range of comorbid conditions, such as heart
    disease, AIDS, or cancer (a total of 22 conditions).
    Each condition is assigned a score of 1, 2, 3, or 6, depending on the
    risk of dying associated with each one. Scores are summed to provide a
    total score to predict mortality.

    Reference:
        https://en.wikipedia.org/wiki/Comorbidity
        http://isocentre.wikidot.com/data:charlson-s-comorbidity-index
        https://www.ncbi.nlm.nih.gov/pubmed/15617955

    a=charlson(df=df, cols=['icdmain', 'icdbi'], sep=',', dot_notation=False)

    """
    infarct = 'I21* I22* I25.2'.split()
    heart = 'I09.9 I11.0 I13.0, I13.2 I25.5 I42.0 I42.5-I42.9 I43* I50* P29.0'.split()
    vascular = '70* I71* I73.1 I73.8 I73.9 I77.1 I79.0 I79.2 K55.1 K55.8 K55.9 Z95.8 Z95.9'.split()
    cerebro = 'G45* G46* H34.0 I60* I69*'.split()
    dementia = 'F00* F03* F05.1 G30* G31.1'.split()
    pulmonary = 'I27.8 I27.9 J40* J47* J60* J67* J68.4, J70.1, J70.3'.split()
    tissue = 'M05* M06* M31.5 M32* M34* M35.1, M35.3 M36.0'.split()
    ulcer = ['K25*-K28*']
    liver = 'B18* K70.0-K70.3 K70.9 K71.3-K71.5 K71.7 K73* K74* K76.0 K76.2-K76.4 K76.8 K76.9 Z94.4'.split()
    diabetes = ['E10.0', 'E10.l', 'E10.6', 'E10.8', 'E10.9', 'E11.0', 'E11.1', 'E11.6', 'E11.8', 'E11.9', 'E12.0',
                'E12.1', 'El2.6', 'E12.8', 'El2.9', 'E13.0', 'E13.1', 'E13.6', 'E13.8', 'E13.9', 'E14.0', 'E14.1',
                'E14.6', 'E14.8', 'E14.9']
    hemiplegia = ['G04.1', 'G11.4', 'G80.1', 'G80.2', 'G81*', 'G82*', 'G83.0-G83.4', 'G83.9']
    renal = ['I12.0', 'I13.1', 'N03.2-N03.7', 'N05.2-N05.7', 'N18*', 'N19*', 'N25.0', 'Z49.0-Z49.2', 'Z94.0', 'Z99.2']
    dorgan = ['E10.2', 'E10.3', 'E10.4', 'E10.5', 'E10.7', 'E11.2', 'E11.5', 'E11.7', 'E12.2', 'E12.3', 'E12.4',
              'E12.5', 'E12.7', 'E13.2', 'E13.3', 'E13.4', 'E13.5', 'E13.7', 'E14.2', 'E14.3', 'E14.4', 'E14.5',
              'E14.7']
    tumor = ['C00*-C26*', 'C30*-C34*', 'C37*-41*', 'C43*-C45*', 'C58*-C60*', 'C76*-C81*', 'C85*-C88*', 'C90*-C97*']
    sliver = ['I85.0', 'I85.9', 'I86.4', 'I98.2', 'K70.4', 'K71.1', 'K72.1', 'K72.9', 'K76.5', 'K76.6', 'K76.7']
    mtumor = ['C77*', 'C78*', 'C79*', 'C80*']
    hiv = ['B20*', 'B21*', 'B22*', 'B24']

    points = {
        'infarct': 1,
        'heart': 1,
        'vascular': 1,
        'cerebro': 1,
        'dementia': 1,
        'pulmonary': 1,
        'tissue': 1,
        'ulcer': 1,
        'liver': 1,
        'diabetes': 1,
        'hemiplegia': 2,
        'renal': 2,
        'dorgan': 2,
        'tumor': 2,
        'sliver': 3,
        'mtumor': 6,
        'hiv': 6}

    disease_labels = list(points.keys())

    diseases = [
        infarct,
        heart,
        vascular,
        cerebro,
        dementia,
        pulmonary,
        tissue,
        ulcer,
        liver,
        diabetes,
        hemiplegia,
        renal,
        dorgan,
        tumor	,
        sliver,
        mtumor,
        hiv]

    disease_codes = {}
    for i, disease in enumerate(diseases):
        all_codes = []
        disease_str = disease_labels[i]
        for code in disease:
            expanded_codes = expand_hyphen(code)
            all_codes.extend(expanded_codes)
        disease_codes[disease_str] = all_codes

    expanded_disease_codes = {}
    no_dot_disease_codes = {}

    if not dot_notation:
        for disease, codes in disease_codes.items():
            new_codes = [code.replace('.  ','') for code in codes]
            no_dot_disease_codes[disease] = new_codes
        disease_codes = no_dot_disease_codes

    all_codes = unique_codes(df=df, cols=cols, sep=sep)

    for disease, codes in disease_codes.items():
        expanded_disease_codes[disease] = expand_codes(df=df, codes=codes, cols=cols, sep=sep, codebook=all_codes)

    codelist = []
    for disease in disease_labels:
        codes = expanded_disease_codes[disease]
        codelist.extend(codes)

    rows = get_rows(df=df,codes=codelist, cols=cols, sep=sep)

    subset = df[rows]

    charlson_df = persons_with(df=subset, codes=expanded_disease_codes,
                               cols=['icdmain', 'icdbi'], sep=',')

    for disease, point in points.items():
        charlson_df[disease] = charlson_df[disease] * point

    age_points = df.groupby(pid)[age].min().sub(40).div(10).astype(int)
    age_points[age_points < 0] = 0
    age_points[age_points > 4] = 4

    disease_points = charlson_df.sum(axis=1).fillna(0)
    charlson_index = age_points.add(disease_points, fill_value=0)

    # make truly missing egual to nans and not zero
    # truly missing = no age available and no icd has been recorded (all nans)

    age_nans = age_points[age_points > 0] = 0
    icd_nans = df[cols].notnull().sum(axis=1).sum(level=0)
    icd_nans[icd_nans == 0] = np.nan
    icd_nans[icd_nans > 0] = 0

    charlson_with_nans = charlson_index + age_nans + icd_nans

    return charlson_with_nans

# %%
# def validate(df, cols=None, pid='pid', codebook=None, infer_types=True, infer_rulebook=True, rulebook=None, force_change=False, log=True):
#    """
#    types:
#        - fixed (wrt var x like the pid, but could also be geo)
#        - code (#may be code AND fixed)
#        - range (lowest, highest)
#        - go together, never go together (male never pregnant, person level validation)
#        - sep?
#
#    metadata about each col in a dict:
#        gender
#            -text: wwewenweroernw ewefejr
#            -fixed_for: pid
#            -sep
#            -range
#            -na_rep
#            -missing_allowed
#            -valid
#            -force
#
#            -nb_pids:
#            -nb_rows:
#            -nb_values:
#
#    pid fixed_for gender
#
#    rules:
#    questions:
#    ideally:
#
#    gender should be fixed within pid. If not, use majority pid. Except for pid 45, no change.
#
#    birthyear fixed_for pid
#
#    pid should always exist. If not, use majority pid. If not, delete row.
#
#    birthyear should never be less than 1970
#
#    birthyear should never be larger than 2005
#
#    no variables should have negative values. Except
#
#    event_year should never be larger than death_year
#
#    gender should only take the following values 0, 1, nan. If not, use nan.
#
#    define icd_codes= s34, d45
#
#    defineatc_code=codebook['atc']
#
#    icd should only have values defined in icd_codes
#
#
#
#
#
#
#
#    birthyear > 1970
#
#    rule: for gender fixed_for pid use majority
#    rule
#    """
#
#    for col, pid in fixed_cols:
#        check = df[col].groupby(pid).nunique()
#        check = check[check!=1]
#        if len(check)>0:
#            messages.append[f'{col} is not fixed within all {pid}']
#
#            if force:
#                df.loc[col, pid] = df.loc[col, pid].value_counts()[0]





##
## Maily helper functions
##


def _listify(string_or_list):
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


def _sniff_sep(df, cols=None, possible_seps=[',', ';', '|'], n=1000, sure=False, each_col=False):
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

    cols = _listify(cols)

    df = _stringify_cols(df=df, cols=cols)

    # fix args depending on whether a series or df is input and if cols is specified
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


def _get_some_id(df,
                 codes=None,
                 cols=None,
                 some_id='pid',
                 sep=None):
    """
    help function for all get functions that gets ids based on certain filtering criteria

    some_id is the column with the info to be collected (pid, uuid, event_id)


    """

    codes = _listify(codes)
    cols = _listify(cols)

    cols = _expand_cols(df=df, cols=cols)

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

    pids = set(df[b][some_id].unique())

    return pids


def _fix_cols(df, cols):
    if not cols:
        cols = list(df.columns)

    cols = _expand_cols(df=df, cols=cols)
    return cols


def _fix_codes(df, codes=None, cols=None, sep=None, merge=False, group=False):
    if not codes:
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False)[:5]

    codes = _format_codes(codes=codes, merge=merge)
    codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep, merge=merge, group=group)
    return codes


def _fix_args(df, codes=None, cols=None, sep=None, merge=False, group=False, _sniffsep=True):
    # Use all columns if no column is specified
    # Series if converted to df (with pid column, assumed to be in the index)
    if not cols:
        cols = list(df.columns)
    else:
        cols = _expand_cols(df=df, cols=cols)

    if _sniffsep:
        sep = _sniff_sep(df=df, cols=cols)

    if not codes:
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False).index[:5]
        codes = list(codes)

    codes = _format_codes(codes=codes, merge=merge)
    codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep, merge=merge, group=group)
    codes = _format_codes(codes=codes, merge=merge)

    # useful to have full codelist (of codes only, after expansion)
    full_codelist = set()
    for name, codelist in codes.items():
        full_codelist.update(set(codelist))
    allcodes = list(full_codelist)

    return codes, cols, allcodes, sep


def _to_df(df, cols=None):
    if isinstance(df, pd.Series):
        df = df.to_frame()
        cols = list(df.columns)
        df['pid'] = df.index.values
    return df, cols


def _get_allcodes(codes):
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
        allcodes = _listify(codes)
    return allcodes


def _expand_cols(df, cols, star=True, hyphen=True, colon=True, regex=None):
    """
    Expand columns with special notation to their full column names

    """

    cols = _listify(cols)

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
    exprs = _listify(expr)
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

    exprs = _listify(expr)

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

    exprs = _listify(expr)
    all_codes=[]

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
                multiplier = 1

            no_dec_lower = int(lower_num * multiplier)
            no_dec_upper = int((upper_num) * multiplier) + 1

            if '.' in lower_str:
                codes = [lower.replace(lower_str, str(num / multiplier).zfill(length)) for num in range(no_dec_lower, no_dec_upper)]
            else:
                codes = [lower.replace(lower_str, str(num).zfill(length)) for num in range(no_dec_lower, no_dec_upper)]


        else:
            codes = [expr]
        all_codes.extend(codes)
    return all_codes


def _stringify_cols(df, cols):
    """
    Stringify some cols - useful since many methods require code column to be a string
    """

    for col in cols:
        df[col] = df[col].astype(str)
    return df


def _format_codes(codes, merge=True):
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

            _format_codes(codes, merge=False)

    TODO: test for correctness of input, not just reformat (is the key a str?)
    """
    codes = _listify(codes)

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


def _expand_regex(expr, full_list):
    exprs = _listify(expr)

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


def _reverse_dict(dikt):
    new_dict = {}
    for name, codelist in dikt.items():
        codelist = _listify(codelist)
        new_dict.update({code: name for code in codelist})
    return new_dict


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
        df, cols = _to_df(df=df, cols=cols)
        codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep, merge=merge, group=group)
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


##
## stringify functions
##

def stringify_durations(df,
                        codes=None,
                        cols=None,
                        pid='pid',
                        step=120,
                        sep=None,

                        event_start='in_date',
                        event_end=None,
                        event_duration='ddd',

                        first_date=None,
                        last_date=None,
                        censored_date=None,

                        na_rep='-',
                        time_rep=',',

                        merge=True,
                        info=None,
                        report=True):
    """
    Creates a string for each individual describing the time duration selected code events (example: a-, ad, --, a)

    Args:
        df: dataframe
        codes: codes to be used to mark an event
        cols: columns with the event codes
        pid: column with the personal identification number
        event_start: column containing the date for the event
        sep: the separator used between events if a column has multiple events in a cell
        keep_repeats: identical events after each other are reduced to one (if true)
        only_unique: deletes all events that have occurred previously for the individual (if true)

    Returns:
        series with a string that describes the events for each individual
    Example:
    >>> codes={'L04A*' : 'i', 'L04AB*' : 'a', 'H02*' : 'c'}
    >>> pr=pr.set_index('pid_index')
    >>> pr['first_date'] = pr.groupby('pid')['date'].min()
    >>> events=stringify_durations(df=df, codes=codes, col='ncmpalt', start='start_date', first_date='first', dataset_end_date="01-01-2018")


    background
        to identify treatment patters, first stringify each treatment,
        then aggregate the different treatments to one string
        each "cell" in the string (separated by sep) represent one time unit
        the time unit can be further aggregated to reduce the level of detail

    example output (one such row for each person)
        a---s, a---, ai-s, a---, ----

        Interpretation: A person with event a and s in first time perod, then a only in second,
        the a, i and s in the third, a only in fourth and no events in the last

    purpose
        examine typical treatment patterns and correlations
        use regex or other string operations on this to get statistcs
        (time on first line of treatment, number of switches, stops)

    """
    # drop rows with missing observations in required variables

    if report:
        obs = len(df)
        npid = df[pid].nunique()
        rows = get_rows(df=df, codes=codes, cols=cols, sep=sep)
        code_obs = len(df[rows])
        code_npid = df[rows][pid].nunique()

    df = df.dropna(subset=[pid, event_start])

    if event_end:
        df = df.dropna(subset=[event_end])
    elif event_duration:
        df = df.dropna(subset=[event_duration])
        if df[event_duration].min() < 0:
            print('Error: The specified duration column contains negative values')
    else:
        print('Error: Either event_end or event_duration has to be specified.')

    # find default min and max dates
    # will be used as starting points for the string
    # if first_date and last_date are not specified
    min_date = df[event_start].min()
    max_date = df[event_start].max()

    # drop rows outside specified time period of interest
    if first_date:
        if first_date in df.columns:
            df = df[df[event_start] >= df[first_date]]
        elif isinstance(first_date, dict):
            pass
        else:
            # if first_date is not a column name, it is assumed to be a date
            try:
                min_date = pd.to_datetime(first_date)
                df = df[df[event_start] >= min_date]
            except:
                print(
                    'Error: The first_date argument has to be on of: None, a dict, a column name or a string that represents a date')

    if last_date:
        if last_date in df.columns:
            df = df[df[event_start] >= df[last_date]]
        elif isinstance(last_date, dict):
            pass
        else:
            try:
                max_date = pd.to_datetime(last_date)
                df = df[df[event_start] <= max_date]
            except:
                print(
                    'Error: The last_date argument has to be on of: None, a dict, a column name or a string the represents a date')

    # note an individual min date cannot be before overall specified min date
    # should raise error if user tries this
    # same with max: individual cannot be larger than overall

    max_length_days = (max_date - min_date).days
    max_length_steps = int(max_length_days / step)

    # if codes are not specified, use the five most common codes
    if not codes:
        cols = _expand_cols(_listify(cols))
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False)[:4]

    # fix formatting of input (make list out of a string input and so on)
    codes, cols, old_codes, replace = _fix_args(df=df, codes=codes, cols=cols, sep=sep)

    # get the rows that contain the relevant codes
    rows = get_rows(df=df, codes=codes, cols=cols, sep=sep)
    subset = df[rows].copy()  # maybe use .copy to avoid warnings? but takes time and memory
    if report:
        sub_obs = len(subset)
        sub_npid = subset[pid].nunique()

    subset = subset.sort_values([pid, event_start])
    subset = subset.set_index(pid, drop=False)
    subset.index.name = 'pid_index'

    # find start and end position of each event (number of steps from overall min_date)
    # to do: do not use those column names (may overwrite original names), use uuid names?
    subset['start_position'] = (subset[event_start] - min_date).dt.days.div(step).astype(int)

    if event_end:
        subset['end_position'] = (subset[event_end] - min_date).dt.days.div(step).astype(int)
    elif event_duration:
        subset['end_date'] = subset[event_start] + pd.to_timedelta(subset[event_duration].astype(int), unit='D')
        subset['end_position'] = (subset['end_date'] - min_date).dt.days.div(step).astype(int)

    # to do: may allow duration dict?
    # for instance: some drugs last 15 days, some drugs last 25 days . all specified in a dict

    # create series with only the relevant codes for each person and position
    code_series = extract_codes(df=subset.set_index([pid, 'start_position', 'end_position']),
                                codes=replace,
                                cols=cols,
                                sep=sep,
                                new_sep=',',
                                merge=True,
                                out='text')

    # May need to unstack if two events in same row
    # for now: Just foce it to be 1
    if code_series.apply(len).max() > 1:
        code_series = code_series.str[0]

    # base further aggregation on the new extracted series with its col and codes
    col = code_series.name
    codes = code_series.name.split(', ')

    # drop duplicates (same type of even in same period for same individual)
    code_series = code_series.reset_index().drop_duplicates().set_index(pid, drop=False)

    ## make dict with string start and end positions for each individual
    # explanation:
    # the string is first made marking events in positions using calendar time
    # but often we want the end result to be strings that start at specified
    # individual dates, and not the same calendar date for all
    # for instance it is often useful to start the string at the date the
    # person receives a diagnosis
    # same with end of string: strings may end when a patient dies
    # user can specify start and end dates by pointing to columns with dates
    # or they may specify an overall start and end date
    # if individual dates are specified, the long string based on calendar
    # time is sliced to include only the relevant events

    if first_date:
        # if a column is specified
        if first_date in subset.columns:
            start_date = subset.groupby(pid)[first_date].first().dropna().to_dict()
        # do nothing if a dict mapping pids to last_dates is already specified
        elif isinstance(first_date, dict):
            pass
        # if a single overall date is specified
        else:
            date = pd.to_datetime(first_date)
            start_date = {pid: date for pid in subset[pid].unique()}
        # convert start date to start position in string
        string_start_position = {pid: int((date - min_date).days / step)
                                 for pid, date in start_date.items()}

    if last_date:
        if last_date in subset:
            end_date = subset.groupby(pid)[last_date].first().dropna().to_dict()
        # do nothing if a dict mapping pids to last_dates is already specified
        elif isinstance(last_date, dict):
            pass
        else:
            date = pd.to_datetime(last_date)
            end_date = {pid: date for pid in subset[pid].unique()}
        # convert date to position in string
        string_end_position = {pid: (date - min_date).dt.days.div(step).astype(int)
                               for pid, date in end_date.items()}

        # takes dataframe for an individual and makes a string with the events

    def make_string(events):
        # get pid of individual (required to find correct start and end point)
        person = events[pid].iloc[0]

        # make a list of maximal length with no events
        event_list = ['-'] * (max_length_steps + 1)

        from_to_positions = tuple(zip(events['start_position'].tolist(), events['end_position'].tolist()))

        # loop over all events the individual has and put code in correct pos.
        for pos in from_to_positions:
            event_list[pos[0]:pos[1]] = code
        event_string = "".join(event_list)

        # slice to correct start and end of string (if specified)
        if first_date:
            event_string = event_string[string_start_position[person]:]
        if last_date:
            max_position = int((max_date - min_date).days / step)
            event_string = event_string[:-(max_position - string_end_position[person])]
        return event_string

    # new dataframe to store each string for each individual for each code
    string_df = pd.DataFrame(index=code_series[pid].unique())

    # loop over each code, aggregate strong for each individual, store in df
    for code in codes:
        code_df = code_series[code_series[col].isin([code])]
        code_df.index.name = 'pid_index'  # avoid future error from pandas pid in both col and index
        stringified = code_df.groupby(pid, sort=False).apply(make_string)
        string_df[code] = stringified

    if merge:
        string_df = interleave_strings(string_df)

    if report:
        final_obs = len(subset)
        final_npid = len(string_df)
        print(f"""
                                     events,  unique ids
              Original dataframe     {obs}, {npid} 
              Filter codes           {code_obs}, {code_npid}
              Filter missing         {sub_obs}, {sub_npid}
              Final result:          {final_obs}, {final_npid}""")
    return string_df


def interleave_strings(df, cols=None, sep=" ", nan='-', agg=False):
    """
    Interleaves strings in two or more columns

    parameters
        cols : list of columns with strings to be interleaved
        nan : value to be used in place of missing values
        sep : seperator to be used between time periods
        agg : numeric, used to indicate aggregation of time scale
                default is 1

    background
        to identify treatment patters, first stringify each treatment,
        then aggregate the different treatments to one string
        each "cell" in the string (separated by sep) represent one time unit
        the time unit can be further aggregated to reduce the level of detail

    example output (one such row for each person)
        a---s, a---, ai-s, a---, ----

        Interpretation: A person with event a and s in first time perod, then a only in second,
        the a, i and s in the third, a only in fourth and no events in the last

    purpose
        examine typical treatment patterns and correlations
        use regex or other string operations on this to get statistcs
        (time on first line of treatment, number of switches, stops)

    """
    # if cols is not specified, use all columns in dataframe
    if not cols:
        cols = list(df.columns)

    if agg:
        for col in cols:
            df[col] = df[col].fillna('-')
            # find event symbol, imply check if all are missing, no events
            try:
                char = df[col].str.cat().strip('-')[0]
            except:
                df[col] = (col.str.len() / agg) * '-'

            missing = '-' * len(cols)

            def aggregator(text, agg):
                missing = '-' * agg
                units = (text[i:i + agg] for i in range(0, len(text), agg))
                new_aggregated = ('-' if unit == missing else char for unit in units)
                new_str = "".join(new_aggregated)
                return new_str
        df[col] = df[col].apply(aggregator, agg=agg)

    if sep:
        interleaved = df[cols].fillna('-').apply(
            (lambda x: ",".join(
                "".join(i)
                for i in zip_longest(*x, fillvalue='-'))),
            axis=1)
    else:
        interleaved = df[cols].fillna('-').apply(
            (lambda x: "".join(chain(*zip_longest(*x, fillvalue='-')))),
            axis=1)

    return interleaved


def overlay_strings(df, cols=None, sep=",", nan='-', collisions='x', interleaved=False):
    """
    overlays strings from two or more columns

    note
        most useful when aggregating a string for events that usually do not happen in the same time frame

    parameters
        cols : list of columns with strings to be interleaved
        nan : value to be used in place of missing values
        collisions: value to be used if there is a collision between events in a position


    background
        to identify treatment patters, first stringify each treatment,
        then aggregate the different treatments to one string
        each "cell" in the string (separated by sep) represent one time unit
        the time unit can be further aggregated to reduce the level of detail

    example output (one such row for each person)
        asaaa--s--aa-s-a

        Interpretation: A person with event a and s in first time perod, then a only in second,
        the a, i and s in the third, a only in fourth and no events in the last

    purpose
        examine typical treatment patterns and correlations
        use regex or other string operations on this to get statistcs
        (time on first line of treatment, number of switches, stops)

    todo
        more advanced handling of collisions
            - special symbols for different types of collisions
            - warnings (and keep/give info on amount and type of collisions)

    """
    # if cols is not specified, use all columns in dataframe
    if not cols:
        cols = list(df.columns)

    interleaved = df[cols].fillna('-').apply(
        (lambda x: "".join(chain(*zip_longest(*x, fillvalue='-')))),
        axis=1)
    step_length = len(cols)

    def event_or_collision(events):
        try:
            char = events.strip('-')[0]
        except:
            char = '-'
        n = len(set(events).remove('-'))
        if n > 1:
            char = 'x'
        return char

    def overlay_individuals(events):

        units = (events[i:i + step_length] for i in range(0, len(events), step_length))

        new_aggregated = (event_or_collision(unit) for unit in units)
        new_str = "".join(new_aggregated)
        return new_str

    interleaved.apply(overlay_individuals)

    return interleaved


def shorten(events, agg=3, missing='-'):
    """
    create a new and shorter string with a longer time step

    parameters
        events: (str) string of events that will be aggregated
        agg: (int) the level of aggregation (2=double the step_length, 3=triple)
    """
    try:
        char = events.strip('-')[0]
    except:
        char = '-'
    units = (events[i:i + agg] for i in range(0, len(events), agg))
    new_aggregated = ('-' if unit == missing else char for unit in units)
    new_str = "".join(new_aggregated)
    return new_str


def shorten_interleaved(text, agg=3, sep=',', missing='-'):
    """
    text="a-si,a--i,a-s-,--si,---i,--s-"

    shorten_interleaved(c, agg=2)
    """
    units = text.split(sep)
    ncodes = len(units[0])
    nunits = len(units)

    unitlist = [units[i:i + agg] for i in range(0, nunits, agg)]
    charlist = ["".join(aggunit) for aggunit in unitlist]
    unique_char = ["".join(set(chain(chars))) for chars in charlist]
    new_str = ",".join(unique_char)
    # ordered or sorted?
    # delete last if it is not full ie. not as many timee units in it as the others?
    # shortcut for all
    return new_str


def stringify_order(df, codes=None, cols=None, pid='pid', event_start='in_date', sep=None, keep_repeats=True,
                    only_unique=False, _fix=True):
    """
    Creates a string for each individual describing selected code events in the order they occurred

    Args:
        df: dataframe
        codes: codes to be used to mark an event
        cols: columns with the event codes
        pid: column with the personal identification number
        event_start: column containing the date for the event
        sep: the separator used between events if a column has multiple events in a cell
        keep_repeats: identical events after each other are reduced to one (if true)
        only_unique: deletes all events that have occurred previously for the individual (if true)

    Returns:
        series with a string that describes the events for each individual



    Examples:

    >>> codes={'e': '4AB01', 'i': '4AB02', 'a': '4AB04', 'x': '4AB05' , 'g': '4AB06'}

    >>> bio_codes= {'L04AA23': 'n', 'L04AA33': 'v', 'L04AB02': 'i', 'L04AB04': 'a','L04AB06': 'g', 'L04AC05': 'u'}

    >>> bio_codes={'e' : '4AB01', 'i' : '4AB02', 'a' : '4AB04'}

    >>> a=stringify_order(df=df, codes=bio_codes, cols='ncmpalt', pid='pid', event_start='start_date', sep=',', keep_repeats=True, only_unique=False)


    >>> bio_rows=get_rows(df=pr, codes=list(codes.keys()), cols='atc')
    >>> pr['first_bio']=pr[bio_rows].groupby('pid')['date'].min()

    >>> stringify_order(df=pr, codes=codes, cols='atc', pid='pid', event_date='date', sep=',')

    background
        to identify treatment patters, first stringify each treatment,
        then aggregate the different treatments to one string
        each "cell" in the string (separated by sep) represent one time unit
        the time unit can be further aggregated to reduce the level of detail

    example output (one such row for each person)
        a---s, a---, ai-s, a---, ----

        Interpretation: A person with event a and s in first time perod, then a only in second,
        the a, i and s in the third, a only in fourth and no events in the last

    purpose
        examine typical treatment patterns and correlations
        use regex or other string operations on this to get statistcs
        (time on first line of treatment, number of switches, stops)
    """

    # fix formatting of input
    if _fix:
        df, cols = _to_df(df=df, cols=cols)
        codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep)
    else:
        allcodes=_get_allcodes(codes)

    # get the rows with the relevant columns
    rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
    subset = df[rows].sort_values(by=[pid, event_start]).set_index('pid')

    # extract relevant codes and aggregate for each person
    code_series = extract_codes(df=subset, codes=codes, cols=cols, sep=sep, new_sep='', merge=True, out='text',
                                _fix=False)
    #    if isinstance(code_series, pd.DataFrame):
    #        code_series = pd.Series(code_series)
    string_df = code_series.groupby(level=0).apply(lambda codes: codes.str.cat())

    # eliminate repeats in string
    if not keep_repeats:
        string_df = string_df.str.replace(r'([a-z])\1+', r'\1')

    if only_unique:
        def uniqify(text):
            while re.search(r'([a-z])(.*)\1', text):
                text = re.sub(r'([a-z])(.*)\1', r'\1\2', text)
            return text

        string_df = string_df.apply(uniqify)
    return string_df


def del_repeats(str_series):
    """
    deletes consecutively repeated characters from the strings in a series

    """
    no_repeats = str_series.str.replace(r'([a-z])\1+', r'\1')
    return no_repeats


def del_singles(text):
    """
    Deletes single characters from string
    todo: how to deal with first and last position ... delete it too?

    """
    # text with only one character are by definition singles
    if len(text) < 2:
        no_singles = ''
    else:
        no_singles = "".join([letter for n, letter in enumerate(text[1:-1], start=1) if
                              ((text[n - 1] == letter) or (text[n + 1] == letter))])
        # long textx may not have any singles, so check before continue
        if len(no_singles) < 1:
            no_singles = ''
        else:
            if text[0] == no_singles[0]:
                no_singles = text[0] + no_singles
            if text[-1] == no_singles[-1]:
                no_singles = no_singles + text[-1]

    return no_singles


def stringify_time(df,
                   codes=None,
                   cols=None,
                   pid='pid',
                   step=90,

                   event_start='in_date',
                   event_end=None,

                   first_date=None,
                   last_date=None,

                   censored_date=None,

                   sep=None,
                   merge=True,
                   meta=None):
    """
    Creates a string for each individual describing events at position in time

    Args:
        df: dataframe
        codes: codes to be used to mark an event
        cols: columns with the event codes
        pid: column with the personal identification number
        event_start: column containing the date for the event
        sep: the seperator used between events if a column has multiple events in a cell
        keep_repeats: identical events after each other are reduced to one (if true)
        only_unique: deletes all events that have occurred previously for the individual (if true)

    Returns:
        series with a string that describes the events for each individual

    Example:
        codes={'4AB01': 'e', '4AB02' : 'i', '4AB04' : 'a','4AB05' : 'x', '4AB06' : 'g'}

    df['diagnosis_date']=df[df.icdmain.fillna('').str.contains('K50|K51')].groupby('pid')['start_date'].min()

    stringify_time(df=df,  codes=codes, cols='ncmpalt', pid='pid', event_date='start_date', first_date='diagnosis_date', step=90, sep=',',)

    background
        to identify treatment patters, first stringify each treatment,
        then aggregate the different treatments to one string
        each "cell" in the string (separated by sep) represent one time unit
        the time unit can be further aggregated to reduce the level of detail

    example output (one such row for each person)
        a---s, a---, ai-s, a---, ----

        Interpretation: A person with event a and s in first time perod, then a only in second,
        the a, i and s in the third, a only in fourth and no events in the last

    purpose
        examine typical treatment patterns and correlations
        use regex or other string operations on this to get statistcs
        (time on first line of treatment, number of switches, stops)
    """

    # drop rows with missing observations in required variables
    df = df.dropna(subset=[pid, event_start])

    # find default min and max dates to be used if not user specified
    min_date = df[event_start].min()
    max_date = df[event_start].max()

    # drop rows outside time period of interest
    if first_date:
        if first_date in df.columns:
            df = df[df[event_start] >= df[first_date]]
        else:
            min_date = pd.to_datetime(first_date)
            df = df[df[event_start] >= min_date]

    if last_date:
        if last_date in df.columns:
            df = df[df[event_start] >= df[last_date]]
        else:
            max_date = pd.to_datetime(last_date)
            df = df[df[event_start] <= max_date]

    # note an individual min date cannot be before overall specified min date
    # should raise error if user tries this
    # same with max: individual cannot be larger than overall

    max_length_days = (max_date - min_date).days
    max_length_steps = int(max_length_days / step)

    # if codes are not specified, use the five most common codes
    if not codes:
        cols = _expand_cols(_listify(cols))
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False)[:4]

    # fix formatting of input (make list out of a string input and so on)
    codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep)


    # get the rows that contain the relevant codes
    rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
    subset = df[rows]  # maybe use .copy to avoid warnings?

    # find position of each event (number of steps from overall min_date)
    subset['position'] = (subset[event_start] - min_date).dt.days.div(step).astype(int)

    # create series with only the relevant codes for each person and position
    code_series = extract_codes(df=subset.set_index([pid, 'position']),
                                codes=codes,
                                cols=cols,
                                sep=sep,
                                new_sep=',',
                                merge=True,
                                out='text',
                                _fix=False)

    # base further aggregation on the new extracted series with its col and codes
    col = code_series.name
    codes = code_series.name.split(', ')

    # drop duplicates (same type of even in same period for same individual)
    code_series = code_series.reset_index().drop_duplicates().set_index(pid, drop=False)

    ## make dict with string start end end positions for each individual
    # explanation:
    # the string is first made marking events in positions using calendar time
    # but often we want the end result to be strings that start at specified
    # individual dates, and not the same calendar date for all
    # for instance it is often useful to start the string at the date the
    # person receives a diagnosis
    # same with end of string: strings may end when a patient dies
    # user can specify start and end dates by pointing to columns with dates
    # or they may specify an overall start and end date
    # if individual dates are specified, the long string based on calendar
    # time is sliced to include only the relevant events

    if first_date:
        # if a column is specified
        if first_date in subset.columns:
            start_date = subset.groupby(pid)[first_date].first().dropna().to_dict()
        # if a single overall date is specified
        else:
            date = pd.to_datetime(first_date)
            start_date = {pid: date for pid in subset[pid].unique()}
        # convert start date to start position in string
        start_position = {pid: int((date - min_date).days / step)
                          for pid, date in start_date.items()}

    if last_date:
        if last_date in subset:
            end_date = subset.groupby(pid)[last_date].first().dropna().to_dict()
        else:
            date = pd.to_datetime(last_date)
            end_date = {pid: date for pid in subset[pid].unique()}
        # convert date to position in string
        end_position = {pid: (date - min_date).dt.days.div(step).astype(int)
                        for pid, date in end_date.items()}

    # takes dataframe for an individual and makes a string with the events
    def make_string(events):
        # get pid of individual (required to find correct start and end point)
        person = events[pid].iloc[0]

        # make a list of maximal length with no events
        event_list = ['-'] * (max_length_steps + 1)

        # loop over all events the individual has and put code in correct pos.
        for pos in events['position'].values:
            event_list[pos] = code

        event_string = "".join(event_list)

        # slice to correct start and end of string (if specified)
        if first_date:
            event_string = event_string[start_position[person]:]
        if last_date:
            event_string = event_string[:-(max_position - end_position[person])]
        return event_string

    # new dataframe to store each string for each individual for each code
    string_df = pd.DataFrame(index=code_series[pid].unique())

    # loop over each code, create aggregate string for each individual, store in df
    for code in codes:
        code_df = code_series[code_series[col].isin([code])]
        stringified = code_df.groupby(pid, sort=False).apply(make_string)
        string_df[code] = stringified

    if merge:
        string_df = interleave_strings(string_df)
    return string_df
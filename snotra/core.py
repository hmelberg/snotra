# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 17:26:40 2018

@author: hmelberg
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:17:33 2018

@author: hmelberg
"""

import numpy as np
import pandas as pd
import re
import ast

# needed in read_codebooks()
import glob
import ntpath
import os

from itertools import zip_longest, chain

from itertools import product
from functools import lru_cache

# project structure
#
# core api (key functions/methods
# internal functions/methods
# stringify
# count_person complex expressions with temporal structure
#


##
## Core API
##

# todos:
# expand_codes should immediately return same dict/list/string, if it does not habe notation


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
def get_rows_simple(df, expr, cols=None, sep=None, codebook=None,
                    pid='pid', _fix=True):
    """
    Returns a series with the index and the boolean value for whether the expression is true or not

    Note: no compound expressions using logical operators

    Args:
        df (dataframe): Dataframe with events

        expr (string): A (simple) expression. Example 4A* in icd, age>20

        cols (string, list): Columns where codes are located

        pid (str): Column with the personal identification number

        codebook (list): User specified list of all possible or allowed codes

    Examples:
        expr = '4AB02 in icd'
        expr = '4A* in icd*'
        expr = '4AB02,4AB04 in icd*'
        expr = 'age>20'

    Note:
        Logic:
        - pick out col and code expression, separate and expand
        - get codebook (if needed)
        - use extract code on each expression
        - execute logic
        - return series (bool)

    Examples:
    >>> get_rows_expression(df=npr, expr=expr, cols='icd', sep=',', out='rows'):


    """

    df, cols = _to_df(df, cols)
    # df = df.set_index(pid, drop=False)  # maybe use "if index.name !='pid_index' since indexing may take time

    # quantative expressions
    # expr='age>20'
    # expr = 'day,bday>20'
    # expr = 'age*>20'

    if re.match('[<=>]', expr):
        cols, threshold = re.split('[<=>]', expr)

        if re.match('[*-:,]', cols):
            cols = cols.split(',')
            cols = _expand_cols(df=df, cols=cols)
        if '>' in expr:
            rows = df[cols] > threshold
        elif '<' in expr:
            rows = df[cols] < threshold
        elif '=' in expr:
            rows = df[cols] == threshold
        if len(cols) > 1:
            rows = rows.any(axis=0)
    else:
        codes, cols = expr.split(' in ')
        cols = cols.split(',')
        cols = _expand_cols(df=df, cols=cols)

        if not sep:
            sep = _sniff_sep(df=df, cols=cols)

        codes = codes.split(',')
        codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep, codebook=codebook)
        allcodes = _get_allcodes(codes)
        rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)

    return rows

#%%

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
                   codebook=None):
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
        selection = use_expression(df=df, expr=codes, cols=cols, sep=sep, out='persons', codebook=codebook, pid=pid)
        pids = df[selection][pid].values
    # if an list of codes - ['K50', 'K51'] is used as input
    else:
        pids = _get_some_id(df=df, codes=codes, some_id=pid, sep=sep)

    df = df[df[pid].isin(pids)]

    return df


# %%
def count_persons(df, codes=None, cols=None, pid='pid', sep=None,
                  normalize=False, dropna=True, group=False, merge=False,
                  length=None, groupby=None, codebook=None, _fix=True):
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
        length (int): If specified, will only use the number of characters
            from each code as specified by the length parameter (useful to
            count codes at different levels of granularity. For example,
            sometimes oe wants to look at how many people get detailed codes,
            other times the researcher wants to know only how many get general
            atc codes, say the first four characters of the atc)

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
    if normalize:
        sum_persons = sub[pid].nunique()

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
            if not codes:
                counted = _count_persons_all_codes(df=sub, cols=cols, pid=pid, sep=sep,
                                                   normalize=normalize, dropna=dropna, length=length, groupby=groupby)
                return counted
                # if some codes are specified, expand and format these, and reduce the df to the relevant codes
            else:
                # expands and formats columns and codes input
                codes, cols, allcodes, sep = _fix_args(df=sub, codes=codes, cols=cols, sep=sep, group=group,
                                                       merge=merge)
                rows = get_rows(df=sub, codes=allcodes, cols=cols, sep=sep, _fix=False)
                if not dropna:
                    sum_persons = df[pid].nunique()
                sub = sub[rows].set_index(pid,
                                          drop=False)  # unsure if this is necessary, may drop it. Requred if method on a series? well not as long as we add pid column and recreate a series as a df

        # make a df with the extracted codes
        code_df = extract_codes(df=sub, codes=codes, cols=cols, sep=sep, _fix=False, series=False)

        labels = list(code_df.columns)

        counted = pd.Series(index=labels)

        # maybe delete groupby option, can be done outside df.groupby. apply ...
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
            counted = counted / sum_persons
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

        codes (string, list or dict): The codes for the disease

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
            if not sep:
                sep = _sniff_sep(df=df, cols=cols)
            name = "".join(cols)
            if not codebook:
                if name not in new_codebook:
                    new_codebook[name] = unique_codes(df=df, cols=cols, sep=sep)
            else:
                new_codebook[name] = codebook

        for n, (word, cols) in enumerate(word_cols):
            # added n to number conditions and avoid name conflicts if same condition (but different column)
            worddict = {'___' + word + f'_{n}'.replace('*', '___'): [word]}

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
        worddict = {'___' + word.replace('*', '___') + f'_{n}': [word]}
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

    # conduct search: either just the raw text, the regex, or the one with logical operators (and, or not)
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

        group (bool): Star an other notation remain a single group, not split into individual codes

        out (string, ['text', 'category', 'bool' or 'int']): Datatype of output column(s)

    Notes:
        Can produce a set of dummy columns for codes and code groups.
        Can also produce a merged column with only extracted codes.
        Accept star notation.
        Also accepts both single value columns and columns with compound codes and separators
        Repeat events in same rows are only extracted once


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
        merged = subset.iloc[:, 0].str.cat(subset.iloc[:, 1:].values, sep=new_sep,
                                           na_rep=na_rep)  # strange .T.values seemed to work previouslyi but it should not have
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
def read_codebooks(file=None, select=None, sep=';', encoding='latin-1'):
    """
    Reads a codebook (in csv format)
     file (str): The filname (including the path) to be read
                If no file is specified, all available codebooks in the
                codebook folder will be read
     must_contain (list): List of terms that the filename must contain
                in order to be read
     filename may/should contain (in this order, separated by underscore):
        name (atc, icd)
        version (9, 10)
        year (2017)
        language (no, eng, ger)
        country?
     it may also contain other terms that you want to select on later
     example: icd_9_2017_eng.csv
     codebook = read_codebooks()
     # must deal with integer vs. regular codes problem
    """
    # todo: make a prober config file
    from snotra import _PATH
    import glob
    import ntpath
    import os
    if not file:
        path = _PATH.replace('core.py', 'codebooks/')
        file = glob.glob(path + '*.csv')
    files = _listify(file)
    books = []
    for file in files:
        name = ntpath.basename(file)
        name = os.path.splitext(name)[0]
        if select:
            wanted_words = set(select.split())
            file_words = set(name.split('_'))
            if len(wanted_words) == len(file_words & wanted_words):
                df = pd.read_csv(file, sep=sep, encoding=encoding)
                df['source'] = name
                books.extend([df])
        else:
            df = pd.read_csv(file, sep=sep, encoding=encoding)
            df['source'] = name
            books.extend([df])
    books = pd.concat(books, ignore_index=True, axis=0, sort=False)
    return books

#%%
def label(df, labels=None, select=None, file=None, codebook=None, path=None, info=None):
    """
    Translate codes to text labels

    df : series, dataframe
    codebook : a dataframe with codes and explanation of codes (if not specified, it will read codebooks from the sub-path codebooks)
    labels : a dictionary of codes to labels (optional, if not specified it will be created from codebooks)
    file : path and filename of a specific codebook to be used (a csv file)
    select  : string. if specified, only codebook files with the words in the filename will be used
    #path : string. if specified, will read codebooks from this path
    """
    if not labels:
        if not codebook:
            codebooks = read_codebooks(select=select, file=file)
        labels = labels_from_codebooks(codebook=codebooks, select=select)

    df = df.rename(index=labels)
    return df


#%%

def labels_from_codebooks(codebook, select=None, code='code', text='text', only_valid_codes=False):
    """
    makes a dictionary of code to labels based on two columns in the codebook
     """
    # must deal with integer vs. regular codes problem
    if select:
        books = []
        words = set(select.split())
        sources = codebook['source'].unique()
        for source in sources:
            nameset = set(source.split('_'))
            if len(words) == len(nameset & words):
                books.append(source)
        codebook=codebook[codebook.source.isin(books)]

    codedikt = codebook[[code, text]].set_index(code).to_dict()[text]
    return codedikt



#%%

def _count_persons_all_codes(df, cols=None, pid='pid', sep=None,
                             normalize=False, dropna=True, length=None, groupby=None):
    """
    helper functions called from count_persons when no codes are specified
     #todo: deal with dropna
    """
    sub = df
    cols = _expand_cols(df=sub, cols=cols)
    if normalize:
        sum_persons = sub[pid].nunique()
    if not sep:
        sep = _sniff_sep(df=sub, cols=cols)
    # codes=unique_codes(df=sub, cols=cols, sep=sep)
    # reduce length of codes to the desired level of detail, specified by length
    # this works if coode columns have single values
    if length:
        for col in cols:
            sub[col] = sub[col].str[0:length]
    # special (easy) case if only one column
    sub = sub.set_index(pid)
     # special (easy) case if alll for single valued columns
    allcounted = None
    for col in cols:
        counted = sub[col].reset_index().groupby(col)[pid].nunique()
        if not allcounted:
            allcounted = counted
        else:
            allcounted = allcounted + counted
    if normalize:
        counted = counted / sum_persons
    else:
        counted = counted.astype(int)
    counted = counted.sort_values(ascending=False)
    return counted


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

    a = a.replace(labels)
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

    #if normalize:
    #    t1['value'] = t1['value'] / t1['value'].sum()

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
        tumor,
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
            new_codes = [code.replace('.  ', '') for code in codes]
            no_dot_disease_codes[disease] = new_codes
        disease_codes = no_dot_disease_codes

    all_codes = unique_codes(df=df, cols=cols, sep=sep)

    for disease, codes in disease_codes.items():
        expanded_disease_codes[disease] = expand_codes(df=df, codes=codes, cols=cols, sep=sep, codebook=all_codes)

    codelist = []
    for disease in disease_labels:
        codes = expanded_disease_codes[disease]
        codelist.extend(codes)

    rows = get_rows(df=df, codes=codelist, cols=cols, sep=sep)

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
                 codes,
                 cols,
                 some_id,
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

    if (_sniffsep) & (not sep):
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


def _get_mask(df,
              codes,
              cols,
              sep=None):
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
        b = pd.DataFrame()
        for col in cols:
            b[col] = df[col].str.contains(expanded_codes_regex, na=False).values

    # if single value cells only
    else:
        b = df[cols].isin(expanded_codes)

    return b


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
    all_codes = []

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
                codes = [lower.replace(lower_str, str(num / multiplier).zfill(length)) for num in
                         range(no_dec_lower, no_dec_upper)]
            else:
                codes = [lower.replace(lower_str, str(num).zfill(length)) for num in range(no_dec_lower, no_dec_upper)]


        else:
            codes = [expr]
        all_codes.extend(codes)
    return all_codes


def _stringify_cols(df, cols):
    """
    Stringify some cols - useful since many methods erquire code column to be a string
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
# %%
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

                        ncodes=None,

                        no_event='-',
                        time_sep='|',

                        merge=True,
                        info=None,
                        report=False):
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

    >>> codes={'i' : ['4AB02', 'L04AB02'], 'a': ['4AB04','L04AB04']}
    >>> events=sa.stringify_durations(df=mdf, codes=codes, cols='codes',
    event_start='date', first_date=None, sep=',', merge=True)

    >>> codes={'i' : ['4AB02', 'L04AB02'], 'a': ['4AB04','L04AB04']}
    >>> codes={'i' : ['L04*'], 'b': ['4AB04','L04AB04']}


    >>> codes = {'i':'L01BB02 L04AX03 L01BA01 L04AD01 L04AD02 L04AA06'.split(),
                 'b':'L04AB02 L04AB04 L04AB06 L04AA33 L04AC05 L04AA23'.split()}


    >>> events=sa.stringify_durations(df=mdf, codes=codes, cols='codes',
    event_start='date', first_date=None, sep=',', merge=False, step=100)

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
        if isinstance(codes, dict):
            allcodes = _get_allcodes(codes)
        elif isinstance(codes, str):
            allcodes = _listify(codes)
        # todo: also possible notational codes! better eliminate this?
        rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep)
        code_obs = len(df[rows])
        code_npid = df[rows][pid].nunique()

    df = df.dropna(subset=[pid, event_start])  # also drop if codes is missing

    if event_end:
        df = df.dropna(subset=[event_end])
    elif event_duration:
        df = df.dropna(subset=[event_duration])
        if df[event_duration].min() < 0:
            print('Error: The specified duration column contains negative values. They are dropped')
            df = df[df[event_duration] >= 0]
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
        if not ncodes: ncodes = 4
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False)[:ncodes]

    # fix formatting of input (make list out of a string input and so on)
    codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep)

    # get the rows that contain the relevant codes
    rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
    subset = df[rows].copy()  # maybe use .copy to avoid warnings? but takes time and memory
    subset = subset.set_index(pid, drop=False)
    subset.index.name = 'pid_index'
    subset = subset.sort_values([pid, event_start])

    if report:
        sub_obs = len(subset)
        sub_npid = subset[pid].nunique()

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
                                codes=codes,
                                cols=cols,
                                sep=sep,
                                new_sep=',',
                                merge=False,
                                out='text',
                                _fix=False)

    unique_codes = list(code_series.columns)

    code_series = pd.melt(code_series.reset_index(),
                          id_vars=['pid', 'start_position', 'end_position'],
                          value_vars=unique_codes)

    # drop duplicates (same type of even in same period for same individual)
    code_series = code_series.drop_duplicates().set_index(pid, drop=False)
    code_series.index.name = 'pid_index'
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
        person = events.index[0]

        # make a list of maximal length with no events
        event_list = [no_event] * (max_length_steps + 1)

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
    string_df.index.name = 'pid_index'

    # loop over each code, aggregate strong for each individual, store in df
    for code in unique_codes:
        code_df = code_series[code_series['value'].isin([code])] # maybe == is better (safer bco compounds + faster?)
        stringified = code_df.groupby(pid, sort=False).apply(make_string)
        string_df[code] = stringified

    if merge:
        string_df = interleave_strings(string_df, no_event=no_event, time_sep=time_sep)

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

#%%
def _make_binary(df, cols=None, no_event=' ', time_sep='|', pad=False):
    if isinstance(df, pd.Series):
        name = df[col].name
        df=df.str.replace(no_event, '0')
        df=df.str.replace(name, '1')
    else:
        # if no cols are selected, use all cols
        if not cols:
            cols = list(df.columns)
        # replace event chars with 1 and no events with 0
        for col in cols:
            name = df[col].name
            df[col]=df[col].str.replace(no_event, '0')
            df[col]=df[col].str.replace(name, '1')
    return df
# %%
def interleave_strings(df, cols=None, time_sep="|", no_event=' ', agg=False):
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
            df[col] = df[col].fillna(no_event)
            # find event symbol, imply check if all are missing, no events
            try:
                char = df[col].str.cat().strip().str.strip('-')[0]  # improvable?
            except:
                df[col] = (col.str.len() / agg) * no_event

            def aggregator(text, agg):
                missing = no_event * agg
                units = (text[i:i + agg] for i in range(0, len(text), agg))
                new_aggregated = (no_event if unit == missing else char for unit in units)
                new_str = "".join(new_aggregated)
                return new_str
        df[col] = df[col].apply(aggregator, agg=agg)

    if time_sep:
        interleaved = df[cols].fillna(no_event).apply(
            (lambda x: time_sep.join(
                "".join(i)
                for i in zip_longest(*x, fillvalue=no_event))),
            axis=1)
    else:
        interleaved = df[cols].fillna('-').apply(
            (lambda x: "".join(chain(*zip_longest(*x, fillvalue=no_event)))),
            axis=1)

    return interleaved


# %%
def left_justify(s, fill=' '):
    """
    after stringify, to make events at same time be in same position
    and no, not as crucial as left-pad!
    """
    nmax = s.apply(len).max()
    s = s.str.pad(width=nmax, side='right', fillchar=fill)
    return s


# %%
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


def shorten(events, agg=3, no_event=' '):
    """
    create a new and shorter string with a longer time step

    parameters
        events: (str) string of events that will be aggregated
        agg: (int) the level of aggregation (2=double the step_length, 3=triple)
    """
    try:
        char = events.strip(no_event)[0]
    except:
        char = no_event
    units = (events[i:i + agg] for i in range(0, len(events), agg))
    new_aggregated = (no_event if unit == no_event else char for unit in units)
    new_str = "".join(new_aggregated)
    return new_str


def shorten_interleaved(text, agg=3, time_sep=',', no_event=' '):
    """
    text="a-si,a--i,a-s-,--si,---i,--s-"

    shorten_interleaved(c, agg=2)

    the original string must have a distinction between time_sep and no_event_sep
    (if not, could try to infer)
    """
    units = text.split(time_sep)
    ncodes = len(units[0])
    nunits = len(units)

    unitlist = [units[i:i + agg] for i in range(0, nunits, agg)]
    charlist = ["".join(aggunit) for aggunit in unitlist]
    unique_char = ["".join(set(chain(chars))) for chars in charlist]
    new_str = time_sep.join(unique_char)
    # ordered or sorted?
    # delete last if it is not full ie. not as many timee units in it as the others?
    # shortcut for all
    return new_str


# %%
def stringify_order(df, codes=None, cols=None, pid='pid', event_start='date',
                    sep=None, time_sep='', first_date=None, last_date=None, period=None, keep_repeats=True,
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

    >>> bio_codes= {'L04AA23': 'n', 'L04AA33': 'v', 'L04AB02': 'i', 'L04AB04': 'a','L04AB06': 'g', 'L04AC05': 'u'}

    >>> bio_codes={'e' : '4AB01', 'i' : '4AB02', 'a' : '4AB04'}

    >>> bio_codes={'i' : '4AB02', 'a' : '4AB04'}

    >>> bio_codes= {'n': ['L04AA23', '4AA23'],
                    'v': ['L04AA33', '4AA33'],
                    'i': ['L04AB02', '4AB02'],
                    'a': ['L04AB04', '4AB04'],
                    'g': ['L04AB06', '4AB06'],
                    'u': ['L04AC05', '4AC05']}


    >>> a=stringify_order(df=df, codes=bio_codes, cols='ncmpalt', pid='pid', event_start='start_date', sep=',', keep_repeats=True, only_unique=False)

    >>> a=sa.stringify_order(df=mdf, codes=bio_codes, cols='codes', pid='pid', first_date='first_ibd',
    event_start='date', sep=',', keep_repeats=False, only_unique=False, time_sep='', period=700)


    >>> bio_rows=get_rows(df=pr, codes=list(codes.keys()), cols='atc')
    >>> pr['first_bio']=pr[bio_rows].groupby('pid')['date'].min()

    >>> stringify_order(df=pr, codes=codes, cols='atc', pid='pid', event_date='date', sep=',')

    >>> stringify_order(df=pr, codes=bio_codes, cols='codes', pid='pid', event_date='date', sep=',')


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

    df.index.name = 'pid_index'  # avoid errors, and yes, require pid to be in index (?)

    df = df.dropna(subset=[pid, event_start])

    if first_date:
        df = df.dropna(subset=[first_date])

        # if a column is specified
        if first_date in df.columns:
            include = (df[event_start] >= df[first_date])
            # if a single overall date is specified
        else:
            date = pd.to_datetime(first_date)
            include = (df[event_start] >= date)
        df = df[include]

    if last_date:
        df = df.dropna(subset=[last_date])

        if last_date in df.columns:
            include = (df[event_start] <= df[last_date])
        else:
            date = pd.to_datetime(last_date)
            include = (df[event_start] <= df[last_date])
        df = df[include]

    # period represents the days from the first_date to be included
    # cannot specify both period and last_date(?)
    if period:
        if first_date:
            end_date = df[first_date] + pd.to_timedelta(period, unit='D')
            include = (df[event_start] <= end_date)
        else:
            time_after = (df[event_start] - df.groupby(pid)[event_start].min()) / np.timedelta64(1, 'D')
            include = (time_after <= period).values  # strange need this, tries to reindex if not
        df = df[include]

    # fix formatting of input
    if _fix:
        df, cols = _to_df(df=df, cols=cols)
        codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep)
    else:
        allcodes=_get_allcodes(codes)

    # get the rows with the relevant columns
    rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
    subset = df[rows]  # do I need to copy?
    subset.index.name = 'pid_index'
    subset = subset.sort_values(by=[pid, event_start]).set_index('pid')

    # extract relevant codes and aggregate for each person
    code_series = extract_codes(df=subset, codes=codes, cols=cols, sep=sep, new_sep='', merge=True, out='text',
                                _fix=False)
    #    if isinstance(code_series, pd.DataFrame):
    #        code_series = pd.Series(code_series)
    string_df = code_series.groupby(level=0).apply(lambda codes: codes.str.cat(sep=time_sep))

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


# %%

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


# %%
def stringify_time(df,
                   codes=None,
                   cols=None,
                   pid='pid',
                   sep=None,
                   step=90,

                   event_start='date',  # use start end
                   nfirst=None,  # ncodes

                   first_date=None,
                   # use just first, last, censored. Accept integers to indicate period/days relative to the start date
                   last_date=None,
                   censored_date=None,

                   time_sep='|',
                   no_event=' ',
                   collision='*',

                   merge=True,
                   info=None):
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
        codes={'i': '4AB02', 'a':'4AB04'}
        codes={'i': ['4AB02','L04AB02'], 'a': ['4AB04', 'L04AB04'], 'e':['4AB01']}


        df['diagnosis_date']=df[df.icdmain.fillna('').str.contains('K50|K51')].groupby('pid')['start_date'].min()

    a=stringify_time(df=mdf,  codes=codes, cols='codes', pid='pid', event_start='date',
    first_date='first_ibd', step=90, sep=',', no_event=' ', time_sep=' ')


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

    # if codes or nfirst are not specified, use the five most common codes
    if not codes:
        cols = _expand_cols(_listify(cols))
        if not nfirst: nfirst = 5
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False)[:nfirst]

    # fix formatting of input (make list out of a string input and so on)
    codes, cols, allcodes, sep = _fix_args(df=df, codes=codes, cols=cols, sep=sep)

    # get the rows that contain the relevant codes
    rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
    subset = df[rows].copy()  # maybe use .copy to avoid warnings?
    subset.index.name = 'pid_index'

    # find position of each event (number of steps from overall min_date)
    subset['position'] = (subset[event_start] - min_date).dt.days.div(step).astype(int)

    subset = subset.sort_values(by=[pid, 'position']).set_index([pid, 'position'])

    # create series with only the relevant codes for each person and position
    code_series = extract_codes(df=subset,
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
    code_series.index.name = 'pid_index'

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
        event_list = [no_event] * (max_length_steps + 1)

        # loop over all events the individual has and put code in correct pos.
        for pos in events['position'].values:
            event_list[pos] = code

        event_string = "".join(event_list)

        # slice to correct start and end of string (if specified)
        if first_date:
            event_string = event_string[start_position[person]:]
        if last_date:
            event_string = event_string[:-(max_length_steps - end_position[person])]
        return event_string

    # new dataframe to store each string for each individual for each code
    string_df = pd.DataFrame(index=code_series[pid].unique())

    # loop over each code, create aggregate string for each individual, store in df
    for code in codes:
        code_df = code_series[code_series[col].isin([code])]
        stringified = code_df.groupby(pid, sort=False).apply(make_string)
        string_df[code] = stringified

    if merge:
        string_df = interleave_strings(string_df, no_event=no_event, time_sep=time_sep)
    return string_df


# %%

#
## functions dealing with complicated and possibly time dependent expressions
#

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 15:30:16 2018

@author: hmelberg_adm
"""



# ideas: optimize using arctic and/or mongodb (same type of queries ... at least existencec, but perhaps not sum/cumsum, count etc?)
# %%
def _get_before(splitted, n):
    count = splitted[n].count(')')
    for start, word in enumerate(list(reversed(splitted[:n]))):
        adding = word.count(')')
        subtracting = word.count('(')
        count = count + adding - subtracting
        print(n, count, word)
        if count == 0:
            return n - start
    else:
        print('Error: Expression has wrong number of parenthesis')
    return


def _get_after(splitted, n):
    count = splitted[n].count('(')

    for end, word in enumerate(splitted[n + 1:]):
        adding = word.count('(')
        subtracting = word.count(')')
        count = count + adding - subtracting
        print(n, count, word)
        if count == 0:
            return n + end
    else:
        print('Error: Expression has wrong number of parenthesis')
    return


# %%
def _fix_space(expr):
    no_space_before = r'(\s)([<=>,])'
    no_space_after = r'([<=>,])(\s)'

    expr = re.sub(no_space_before, r'\2', expr)
    expr = re.sub(no_space_after, r'\1', expr)
    return expr


# %%
def _format_expr(expr, rule=None, info=None, insert_cols=True):
    """

    >>> expr = '(first 4A and first 4B) before (A and (C or D)) before (A or (B and C))'
    >>> _format_expr(expr)

    >>> expr = 'min 4 4A in atc1, atc2  and (first 4A and first 4B) before ?[3rd, 4th, 5th] B and (A and (C or D)) before (A or (B and C))'
    >>> _format_expr(expr)

    """
    old_expr = expr

    # fix string so can split on space w/o problems
    expr = _fix_space(expr)

    # insert & and in before/after/from within statements with parenthesis
    # \((.*?) before (.*?)\)
    splitted = expr.split()
    new_splitted = splitted
    for n, word in enumerate(splitted):
        if word == 'before':
            if ')' in splitted[n - 1]:
                start = _get_before(splitted, n)
                for x in range(start, n):
                    if splitted[x] == 'and' or splitted[x] == 'or':
                        new_splitted[x] = new_splitted[x].replace('and', '&').replace('or', '|')

            if '(' in splitted[n + 1]:
                end = _get_after(splitted, n)
                for x in range(n, end + 1):
                    if splitted[x] == 'and' or splitted[x] == 'or':
                        new_splitted[x] = new_splitted[x].replace('and', '&').replace('or', '|')
                # new_splitted[start:n] = ['&' if x=='and' else x for x in splitted[start:n]]

            expr = " ".join(new_splitted)
    # insert date variable 'in date' after days etc (hmm etc is a problem)

    # insert externally defined variables
    expr = _insert_external(expr)

    # insert 'in col' for variables with known cols (specified in info)
    if insert_cols:
        expr = _insert_cols(expr=expr, rule=rule, info=info)

    # replace some phrases
    replace = {' but not ':' and not '}
    for old, new in replace.items():
        expr = expr.replace(old, new)

    # remove 'of'

    # before last 5 = before last 5th (or before -5th)
    return expr


#%%
def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# %%
def _insert_cols(expr, info=None, rule=None):
    """
    insert column names in expressions (in col ...)

    note: may want to try a different approach: instead of identifying codes y eliminating all else, use positive identification: all terms found in codebooks (or rulebook?) are codes!

    expr = 'max 2 of 4AB02 before 4AB04'
    expr = 'max 2 of 4AB02 in x before 4AB04'

    expr = '5th of 5th' # the code is here also a keyword ... problem - maybe of as long as we keep the of keyword ... but more difficult when we do not, at least for automatic column labeling!

    expr = 'max 2 of 4AB0552 before 4AB04'

    expr = 'max 2 of 4AB02 in ncmp' # should include zero ?
    expr = 'min ?[1,2,3,4] of 4AB02 in ncmp'
    expr = 'max ?[1,2,3,4] of 4AB02 in ncmp' # should include zero ?
    expr = 'min 2 of days>4'
    expr = 'min 8 of days>6'
    expr = 'min 3 of 4AB02 in ncmp within 200 days'

    insert_cols(expr, rule=col_rules)
    todo1: codes that are in option lists (?[4AB, 4AV] will not be included)
    sol1: special col format after creating all statements from options?
    todo2: dealing with star codes: 4AB* ... what col is that?
    sol2: well, logically impossible (since it could be several cols) unless user gives some rules? could try to infer based on col content, but do not want to guess, so, rely on functions external, based on startswith or len, or just say: provide the col in that case

    """
    ordinal = r'^(-?\d+)(st|nd|rd|th)'  # re to find 3rd etc

    keywords = set(
        'and or not min max exactly before after around within days of first last days event events day between to days'.split())

    skip_startswith = ['?[']  # also 1st, 2nd etc

    skip_next = set('min max exactly first last in'.split())  # not within?

    skip_contains = set('< = >'.split())

    always_code_before = 'and or within before after'  # as well as the last

    expr = _fix_space(expr)

    splitted = expr.split()

    new_splitted = splitted

    splitted.append(
        'and')  # add a term to solve problem with checking last term looking ahead to see if there is a col already

    # identify a code that needs a col by excluding everything that is not a code with a col
    n = 0

    while n < len(splitted) - 1:

        word = splitted[n]
        if word in skip_next:
            n = n + 2  # skip words after min, max etc
        elif word.startswith(r'?['):
            n = n + 1  # skip options expressions
        elif word in keywords:
            n = n + 1  # keywwords do not indicate codes that need col info
        elif re.match(ordinal, word):
            n = n + 1  # skip ordinals (1st etc) .. could be problematic if a code also has this pattern ... impose extra condtions to check? partial solution: require it to startwith, not just match
        elif len(skip_contains & set(word)) > 0:
            n = n + 1  # skip qualitative conditions which have cols: days>5
        elif _is_number(word):
            n=n+1 # hmmm may not want to skip all numbers if we allow pure number codes?
        elif splitted[n + 1] == 'in':
            n = n + 1  # skip if already have provided a col, cannot check this if last word (no problem, never is)
        else:
            code = word.strip().strip(')').strip('(')
            if info:
                if code in info.code2cols:
                    cols = info.code2cols[code]
            elif rule:
                cols = rule(code)
            else:
                print('Error: No rules for assigning columns are given in info or rule')
            new_splitted[n] = new_splitted[n].replace(code, f'{code} in {cols}')
            n = n + 1
    new_expr = " ".join(new_splitted[:-1])

    return new_expr


# %%
def col_rules(code):
    if len(code) > 5: cols = 'atc'
    if len(code) <= 5: cols = 'icd'
    return cols


# %%

@lru_cache()
def _get_conditions(expr):
    split_on = [' or ', ' and ']
    split_rep = ' @split@ '
    for split_word in split_on:
        expr = expr.replace(split_word, split_rep)
    conditions = expr.split(split_rep)
    conditions = [condition.strip('(').strip(')') for condition in conditions]
    return conditions


@lru_cache()
def _get_eval(expr):
    conditions = _get_conditions(expr)
    eval_text = expr
    for n, condition in enumerate(conditions):
        eval_text = eval_text.replace(condition, f'c{n}')
    return eval_text


def _parser_get_rows(df, codes, cols=None, sep=None, codebook=None, info=None):
    # simple existence condition
    if ' ' not in codes.strip():
        rows = df[cols].str.contains(codes).fillna(False)  # rem: make sure if no codes in cols that all is false
    return rows


def _get_cols(condition, cols):
    expr, incols, _ = condition.split(' in ')

    incols = incols.strip()

    if incols == '':
        incols = cols

    return incols


def _get_options(expr):
    """
    Makes a list of all possible statements from an expression, all possible combination of expressions involving ?[x, y, z] that arein the expressions


    expr = 'min ?[2,3,4] of (K50, K51) in icd inside ?[10, 20, 30] days before 4AB02 in ncmp'
    get_options(expr)

    better name: get_expressions?

    """

    # fix string so can split on space w/o problems
    just_comma = [', ', ' , ', ', ']
    expr = expr.replace('  ', ' ')

    for token in just_comma:
        expr = expr.replace(token, ',')

    original = expr

    alternatives = re.findall('(\[.*?\])', expr)

    alt_list = [ast.literal_eval(alternative) for alternative in alternatives]

    combinations = product(*alt_list)

    all_expressions = []
    for n, combination in enumerate(combinations):
        new_expr = original
        for i, value in enumerate(combination):
            new_expr = new_expr.replace('?' + alternatives[i], str(value), 1)
        all_expressions.extend([new_expr])

    return all_expressions


# %%

def _insert_external(expr):
    """
    Replaces variables prefixed with @ in the expression with the
    value of the variable from the global namespace

    Example:
        x=['4AB02', '4AB04', '4AB06']
        expr = '?@x before 4AB02'
        insert_external(expr)
    """
    externals = list(re.findall(r'@(\w) ', expr))
    for external in externals:
        tmp = globals()[external]
        expr = expr.replace(f'@{external} ', f'{tmp} ')
    return expr


# %%

def count_p(df, expr, cols=None, sep=None, codebook=None, info=None, _use_caching=True, insert_cols=True):
    """

    count persons who satisfy the conditions in an expression

    expr = "?['4AB02', '4AB04'] in ncmp"
    expr = '4AB02 in ncmp and 4AB04 in ncmp'
    expr = 'min 10 of 4AB02 in ncmp'
    expr = 'min ?[4,5,6] of 4AB02 in ncmp'
    expr =  'min 6 of 4AB02 in ncmp'
    expr =  'min 10 of 4AB02 in ncmp'

    expr = 'min ?[6,8] of 4AB02 in ncmp'
    expr = '1st of 4AB02 in ncmp'
    expr = '2nd of 4AB02 in ncmp'

    expr = '4AB02 in ncmp before 4AB04 in ncmp'
    expr = '4AB04 in ncmp before 4AB02 in ncmp'
    expr = '4AA23 in ncmp before 4AB02 in ncmp'
    expr = 'max 2 of 4AB02 in ncmp before 4AB04 in ncmp'
    expr = 'max 2 of 4AB02 in ncmp' # should include zero ?
    expr = 'min ?[1,2,3,4] of 4AB02 in ncmp'
    expr = 'max ?[1,2,3,4] of 4AB02 in ncmp' # should include zero ?
    expr = 'min 2 of days>4'
    expr = 'min 8 of days>6'
    expr = 'min 3 of 4AB02 in ncmp within 200 days'

    %time count_p(df=df, expr=expr, cols=None, codebook=None, info=None, sep=',')
    %time count_p(df=df, expr=expr, cols=None, codebook=None, info=info)


    condition=conditions[0]

    expr = '3rd of 4AB04 in ncmp before 3th of 4AB02 in ncmp'

    """
    # if no info is passed explicitly, create one temporarily
    # use this to avoid save results of evaluations
    # to avoid repeat evaluations of same conditions

    # create storage to avoide recalculating same expressions

    if not info: info = Info()

    expr = _format_expr(expr, info=info, rule=rule, insert_cols=insert_cols)

    exprs = _get_options(expr)
    who = {}
    count = {}
    c = []

    for expr in exprs:
        # check if it has been evaluated before
        if expr in info.expr:
            count[expr] = info.expr[expr]
        # no previous result,so calculate and remember
        else:
            eval_text = _get_eval(expr)
            conditions = _get_conditions(expr)
            for i, condition in enumerate(conditions):
                # check if it has been evaluated before
                if condition in info.condition:
                    cpid = info.condition[condition]
                # if not, evaluate and save result in info
                else:
                    cpid = eval_condition(df=df, condition=condition, cols=cols, sep=sep, out='pid', info=info)
                    cpid.name = f'c{i}'
                    info.condition[condition] = cpid
                c.extend([cpid])
            cdf_pid = pd.concat(c, axis=1)
            who[expr] = cdf_pid.eval(eval_text).any(level=0)
            count[expr] = who[expr].sum()
            if _use_caching: info.expr[expr] = count[expr]
    return count


# %%

#%%
class Info():
    def __init__(self, rows=None, cumsum=None, has_happened=None, single=None, condition=None):
        self.rows = {}
        self.pid = {}
        self.cumsum = {}
        self.has_happened = {}
        self.interval = {}
        self.before_l_has_happened = {}
        self.before_r_has_happened = {}
        self.r_exist_pid = {}
        self.after = {}

        self.default = {} # (new) defaults associated with functions
        self.label = {} # the label to use for codes in given column
        self.codebook = {} # codebook associated with a column



        # these are not really necessary ....
        self.single = {}
        self.single_pid = {}
        self.single_rows = {}
        self.single_interval = {}
        self.condition = {}
        self.expr = {}

        self.x = {}
        self.xpid = {}
    def read_codebook(self):
        pass
    def search_codebook(self):
        pass
    def get_codes(self):
        pass
    def define(self):
        pass
    def lookup(self):
        pass

# %%
def eval_condition(df, condition, cols=None, sep=None,
                   codebook=None, out=None, info=None):
    # check if it has been evaluated before
    if not info: info = Info()
    if condition in info.condition:
        return info.condition[condition]

    time = [' after ', ' before', ' around ']

    # interval condition
    if ' within ' in condition:
        persons = eval_within(df, condition, cols=None, sep=None,
                              codebook=None, out='pid', info=info)

    # before, after, around condition
    elif any(word in condition for word in time):
        persons = eval_before_after(df, condition, cols=None, sep=None,
                                    codebook=None, out='pid', info=info)

    # single condition (not relational, not involving other columns)
    else:
        persons = eval_single(df=df, condition=condition, cols=cols, sep=sep,
                              out='pid', info=info)

    # save the result for later use
    info.condition[condition] = persons

    return persons


# %%
def eval_single(df, condition, cols=None, sep=None, codebook=None,
                out='pid', info=None):
    """
    evaluates a single expressions (1st 4A),
    not relational conditions (A before B, within 100 days after etc)

    condition ='first 5 of 4AB02 in ncmp'
    condition ='min 2 of days>10'
    condition ='ddd>10'
    condition ='ddd[4AB02 in codes]>10'
    condition ='ddd[4AB02 in codes].cumsum()>50'
    condition ='sum(ddd[4AB02 in codes])>50'

    a=eval_single(df=npr, condition=condition, sep=',')

    todo: info bank problems after allowing code selections?

    """
    # create temporary storage to avoid recalculations
    if not info: info = Info()

    original_condition = condition

    # no re-evaluation necessary if it it has been evaluated before
    if out == 'pid' and (condition in info.single_pid):
        return info.single_pid[condition]
    elif out == 'rows' and (condition in info.single_rows):
        return info.single_rows[condition]
    elif out == 'interval' and (condition in info.single_interval):
        return info.single_interval[condition]

    quantity = r'[>=<]'  # better to use term comparison
    freq = ['min ', 'max ', 'exactly ']
    first_last_between = [' first ', ' last ', ' between ']
    ordinal = r'(-?\d+)(st |nd |rd |th )'  # re to find and split 3rd into 3 and rd etc

    # is it a functional expression? ddd.cumsum()>10
    # expr="ddd[4AB02 in codes].cumsum()>10"
    # condition=expr
    row_selection = ''

    # select sub df if specifiec by [] after a code
    if ('[' in condition) and (']' in condition):
        row_query = condition.split('[')[-1].split(']')[0]
        row_selection = row_query
        # check if evaluated before
        if row_query in info.single_rows:
            rows = info.single_rows[row_query]
        else:
            condition = condition.replace(f'[{row_query}]', '')

            if ' in ' in row_query:
                row_query = row_query.replace(' in ', ' in: ')  # using old use_expresssion wich requires in with colon

            relevant_rows = use_expression(df=df, cols=cols, expr=row_query, sep=sep)
            info.single_rows[row_query] = relevant_rows
        df = df[relevant_rows]

    # is it a functional expression? ddd.cumsum()>10
    # expr="ddd.cumsum()>10"
    # condition=expr
    # expr='gender.nunique()==1'
    # hmm what about properties like .is_monotonic? (no parenthesis!)
    # if ('.' in condition) and ('(' in condition) and (')' in condition):
    # still imperfect ... a code could also be a column name ... ok usually not also with a period mark in column name so ok
    if ('.' in condition) and (condition.split('.')[0] in df.columns):
        codetext = condition
        codes = re.split('[<=>]', condition)[0]

        if codes in info.single_rows:
            rows = info.single_rows[codes]
        # not evaluated before, so calc
        else:
            cols, funcexpr = condition.split('.')
            # a method
            if '(' in funcexpr:
                func, threshold = funcexpr.split(')')
                func, args = func.split('(')
                rows = pd.eval(f"tmpdf.groupby(['pid'])['{cols}'].transform('{func}', {args}) {threshold}",
                               engine='python')
            # an attribute (like is_monotonic)
            else:
                rows = pd.eval(f"tmpdf.groupby(['pid'])['{cols}'].transform(lambda x: x.{funcexpr})", engine='python')

            info.single_rows[codes] = rows

    # if it is a simple quantiative conditions (oxygen_level>20)
    elif re.search(quantity, condition):
        codetext = condition
        codes = condition.split()[-1]  # code condition always last hmm unnecessary
        # check if evaluated before
        if codes in info.single_rows:
            rows = info.single_rows[codes]
        # not evaluated before, so calc
        else:
            # sum(glucose_level)>10
            # if this, then may skip further processing?
            # well: 1st sum(glucose)>20 ok makes sense, maybe
            # but not: max 5 of sum(glucose)>20 ... well maybe
            # first 5 of sum(glucose)>20
            # if the modifiers does not make sense, the sum might be in the
            # list of other modifiers i.e. first 5, 3rd etc and not a
            # pre-modifier when finding rows (which allows skipping)

            # complex quantitative expression: sum(glucose_level)>10
            # better, more flexible ...: glucose.sum()>10 ... can make any function work, and can pass arguments

            if 'sum(' in codes:  # can use ddd.cumsum() now, keep this to double check
                col, operator = codes.split(')')
                col = col.replace('sum(', '').strip(')')
                eval_text = f"df.groupby(df.index)['{col}'].cumsum(){operator}"
                rows = pd.eval(eval_text, engine='python').fillna(False)  # is fillna false better than dropna here?
            # simple quantitative expression: glucose_level)>10
            else:
                rows = df.eval(codes).fillna(False)
            codecols = codes
            info.single_rows[codecols] = rows


    # code expression (involving a code, not a quantitative expressions
    else:
        codetext, incols = condition.split(' in ')
        codes = codetext.split()[-1].strip()  # codes always last in a simple string after cutting 'in cols'

        if incols.strip() == '':
            cols = cols
        else:
            cols = incols

        codecols = codes + ' in ' + cols + ' row ' + row_selection  # cannot use just codes to store rows since same code may be in different columns, so need to include col in name when storing

        # If conditions is about events in general, create an events column
        if (' event ' in codes) or (' events ' in codes):
            rows = pd.Series(True, index=df.index).fillna(False)
            codecols = ' event '
        # not a quantitative condition or an event conditions, so it is a code condition
        else:
            if codecols in info.rows:
                rows = info.rows[codecols]
            else:
                # cols = _expand_cols(df=df, cols=cols)
                # expanded_codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep)
                # allcodes=_get_allcodes(expanded_codes)
                # rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
                rows = use_expression(df=df, expr=codes + ' in:' + cols, sep=sep)

                info.rows[codecols] = rows

    # is there a prefix to the conditions? if not, isolated condition, just return rows
    # if not, start preparing for calculating conditions with qualifiers
    # todo: quite messy! refactor: one function to evluate the code/expression itself, another to evalute the qualifier?
    if ' ' not in codetext.strip():
        # remember answer
        info.single_rows[codecols] = rows
        info.rows[codecols] = rows

        if out == 'pid':
            endrows = rows.groupby(level=0).any()
            info.single_pid[codecols] = endrows
            info.pid[codecols] = endrows
        else:
            endrows = rows
        return endrows

    # calculate and remember cumsum per person
    # use previous calculation if exist
    if codes in info.cumsum:
        rowscum = info.cumsum[codes]
    else:
        rowscum = rows.groupby(level=0).cumsum()
        info.cumsum[codecols] = rowscum

    ## if not a simple existence condition, it must be one of the conditions below

    # positional condition:  5th of 4a, 3rd to 8th of 4A, (3rd, 4th, 5th) of 4A
    # also allows: 2nd last (or even -5th last)
    if re.match(ordinal, codetext):
        pos_str = condition.split('of ')[0].strip().strip('(').strip(')')
        # pos_re = ordinal.replace(' ', '[ )]|') # last condition may have ) i.e. 25th)
        pos_re = ordinal.replace(' ', '')  # last condition may have ) i.e. 25th)

        pos_nums = re.findall(pos_re, pos_str)
        pos_nums = tuple([int(pos[0]) for pos in pos_nums])

        # if the conditions includes last, need reversed cumsum
        if ' last ' in pos_str or '-' in pos_str:
            n_max = rowscum.groupby(level=0).max().add(1)
            # reversed event number (by id)
            lastrowscum = (rowscum - n_max).abs()
            last_flag = 1
        else:
            last_flag = 0

        # single position: 5th of 4A
        if len(pos_nums) == 1:
            if last_flag:
                select = (lastrowscum == pos_nums)
            else:
                select = (rowscum == pos_nums)

        # from-to positions: 3rd to 8th of 4A, 1st to -3rd
        elif ' to ' in pos_str:
            lower, upper = pos_nums
            if lower > 0:
                aboverows = (rowscum >= lower)
            else:
                aboverows = (lastrowscum >= abs(lower))

            if upper > 0:
                belowrows = (rowscum <= upper)
            else:
                belowrows = (lastrowscum <= abs(upper))

            select = (aboverows & belowrows)

        # list of positions (3rd, 5th, 7th)
        elif pos_str.strip().startswith('('):
            pos_num = [num for num in pos_num if num > 0]
            neg_num = [num for num in pos_num if num < 0]

            if pos_num:
                pos_select = rowscum.isin(pos_nums)
            if neg_num:
                neg_select = rowscum.isin(pos_nums)
            select = (pos_select | neg_select)


    #  freq condition: min 5 of 4A
    elif any(word in codetext for word in freq):
        word, num, _, codes = codetext.split()
        num = int(num)

        if 'min' in word:
            select = (rowscum >= num)
        elif 'max' in word:  # doublecheck!
            n_max = rowscum.max(level=0)
            select = (n_max <= num)
        elif 'exactly' in word:
            select = (rowscum == num)


    # first, last range conditions: first 5 of 4A
    elif any(word.strip() in condition for word in first_last_between):  # regex is better
        word, num, _, codes = codetext.split()
        if '%' not in num:
            num = int(num)
            if 'first' in word:
                select = (rowscum <= num)
            if 'last' in word:
                select = (rowscum >= num)


    # if pct condition: first 10% of 4A
    elif '%' in codetext:
        n_max = rowscum.groupby(level=0).max()
        pct = float(num.split(r'%')[0]) / 100
        pid_num = n_max * pct

        # first 1% of two observations includes 1st obs
        pid_num[pid_num < 1] = 1

        if word == 'first':
            # hmm, generalproblem: drop if pid is missing ...
            select = (rowscum < pid_num)

        if word == 'last':
            select = (rowscum > pid_num)

    # percentile condition
    elif ' percentile ' in codetext:
        event_num = rows.groupby(level=0).cumcount()
        n_count = rowscum.groupby(level=0).size()

        num = float(num.split(r'%')[0]) / 100

        pid_num = n_count * num

        if word == 'first':
            rows = (pid_num < event_num)

        if word == 'last':
            rows = (pid_num > event_num)

    # so far, have marked interval of events for expressions with qualifications
    # (existence conditions are not intervals). example: First 5 of 4A, markes
    # all events in the interval between the 1st and 5th of 4A
    # if we only want to pick the 4A events in this intereval, we and it with
    # the boolena for 4A existence (row). But sometimes we want to keep and use
    # the interval. For instance when the qualifiers are used in before/after
    # statements if the evaluated expression should be returned as 'exact row',
    # 'interval row' or pid existence

    # store and return results
    if out == 'pid':
        endrows = (rows & select)
        endrows = endrows.any(level=0)
        info.single_pid[original_condition] = endrows
        info.single_rows[original_condition] = rows
    elif out == 'interval':
        endrows = select
        info.interval[original_condition] = endrows
    elif out == 'rows':
        endrows = (rows & select)
        info.single_rows[original_condition] = endrows

    return endrows


# %%
def _eval_before_compound(df, expr, out='rows', info=None):
    # todo: deal with existence reversion if there is a not before a condition
    # (A and not B) before (C or D)

    expr = expr.replace(' & ', ' and ').replace(' | ', ' or ')

    eval_text = _get_eval(expr)
    conditions = _get_conditions(expr)

    c = []
    for i, condition in enumerate(conditions):
        crow = eval_condition(df=df, condition=condition, cols=cols, sep=sep,
                              out='rows', info=info)
        crow.name = f'c{i}'
        c.extend([crow])
    cdf_row = pd.concat(c, axis=1)

    cdfbool = pd.DataFrame()
    for col in cdf.columns:
        cdfbool[col] = (df[col].cumsum() > 0)

    eval_bool = cdfbool.eval(eval_text)

    return eval_bool


# %%


def eval_before_after(df, condition, cols=None, sep=None, codebook=None, info=None, out='pid'):
    # allow compound in before and after conditions. How?
    # 0. Must use &, | not 'and', 'or'
    # 1. replace all (space separated) & with 'and' and | with 'or'
    # 2. end up with rows bool turned on for each in rdf and ldf (standard expression eval, row level? modified)
    # 3. end of with c1 before c2, or c1 after c2
    # condition = condition.replace(' & ', ' and ').replace(' | ', ' or ')

    # replace conditions so multiples becomes positional
    # example: before first 5 of 4A --> before 5th of 4A
    # background: first 5 is not satisfied until ALL five have passed, while some other conditions are
    # hmm, better talk about events, but point remains, some events are compound, multiple, intervals
    # may introduce shortcuts for some easy/common evaluations later (like 4A before 4B, easier than 4st of 4A before 1st of 4B?)
    # before and after are also different, may exploit this to create shortcuts

    re.sub(r'last (-?\d+)', r'last \1st', condition)  # or use negative?

    re.sub(r'first (-?\d+)', r'\1st', condition)
    # todo: also fix percentile --> find position, not first 5 percent

    left, right = re.split(' before | after | simultaneously ', condition)
    # shortcut if ' simultaneous ' in condition: ...just a row and

    # is it a simple (A before B) or a composite expressions (A & B) before (C or D)
    # (3rd A & first 5 B) before (8th C or last D)
    # (3rd A & 5th B) before (8th C or last D)

    # check if the left side of the before expression has been calculated
    if left in info.before_l_has_happened:
        l_has_happened = info.before_l_has_happened[left]
    else:
        # is the expressions before or after a single expression (A before B)
        # or a compound (A and B) before (C or D). compounds have parenthesis
        if ')' in left:
            l_has_happened = _eval_before_compound(df, expr=left, out='rows')
            # lrows=l_has_happned # hmm, does this work, need lrows in after statements hmm
        # else, a simple expression
        else:
            ltext, lcols = left.split(' in ')
            if lcols == '': lcols = cols
            lrows = eval_single(df=df, condition=left, cols=cols, sep=sep,
                                codebook=codebook, info=info, out='rows')
            lrowcum = lrows.groupby(level=0).cumsum()
            l_has_happened = (lrowcum > 0)  # unnecessary? can logical and two series even if one is not bool?
        info.before_l_has_happened[left] = l_has_happened

    # check if the right side of the before  expression have been calculated
    if right in info.before_r_has_happened:
        r_has_happened = info.before_r_has_happened[left]
    else:
        if '(' in right:
            r_has_happened = _eval_before_compound(df, expr=right, out='rows')
            rrows = r_has_happned  # hmm, does this work
            # simple
        else:
            rtext, rcols = right.split(' in ')
            if rcols == '': rcols = cols
            rrows = eval_single(df=df, condition=right, cols=rcols, sep=sep,
                                codebook=codebook, info=info, out='rows')
            rrowcum = rrows.groupby(level=0).cumsum()
            r_has_happened = (rrowcum > 0)  # unnecessary? can logical and two series even if one is not bool?
        info.before_r_has_happened[right] = r_has_happened

    if ' before ' in condition:
        before = ~r_has_happened
        if right in info.r_exist_pid:
            r_exist_pid = info.r_exist_pid[right]
        else:
            r_exist_pid = r_has_happened.any(level=0)
            info.r_exist_pid[right] = r_exist_pid
            # final result
        endrows = (l_has_happened & before & r_exist_pid)

    elif ' after ' in condition:
        if right in info.after:
            after = info.after[right]
        else:
            after = (r_has_happened.groupby(
                pid).cumsum() > 1)  # after means after 1 i.e. on 2 or later (so not simultaneously)
            info.after[right] = after

        endrows = _eval_before_compound(df[after], expr=left, out='rows')

    # reduce to whether the condition is satisfied for each person (peronal id is in index)
    cpid = endrows.any(level=0)

    return cpid


# %%
def eval_within(df, condition, cols=None, sep=None, codebook=None, info=None, out='pid'):
    """

    expr= '4AB02 within 100 days after 4AB04'
    expr= 'min 2 of 4AB02 within 100 days'
    expr= '4AB02 within 50 to 100 days before 4AB04'
    expr= '4AB02 within 50 to 100 days before 4AB04'

    # maybe use inside on some?



    expr= 'min 4 of 4AB02 in ncmp within 100 days'

    expr= 'min 2 of 4AB02 within last 100 days'

    expr= 'min 2 of 4AB02 within 100 days from end'

    expr= 'min 2 of 4AB02 within first 100 days'

    expr= 'between 2 and 5 of 4AB02 within first 100 days' # avoid and? well, just use format and replace it with to?


    expr= 'min 2 of 4AB02 within 100 days from beginning'

    expr= 'min 2 of 4AB02 within 1st of 4AB04 to 5th of 4AB04'
    expr= 'min 2 of 4AB02 within 1st of 4AB06 to 3rd of 4AB04'

    expr= 'min 2 of 4AB02 within first 20 of 4AB04'

    expr= '3rd of 4AB02 within first 20 of 4AB04'

    expr= 'min 2 of 4AB02 within 100 days from 5th of 4AB04'
    expr = '3 or more of 4ab within 100 days'
    wstart, wend

    expr= 'min 4 of 4AB02 in ncmp within 100 days'
    expr= "min 4 of ncmp=='4AB02' within 100 days"
    expr= "at least 4 of ncmp=='4AB02' within 100 days"
    expr= "more than 4 of ncmp=='4AB02' within 100 days" # best language
    expr= "less than 4 of ncmp=='4AB02' within 100 days" # best language
    expr= "between 4 and 7 of ncmp=='4AB02' within 100 days" # best language, inclusive or exclusive between
    expr= "5 or more of ncmp=='4AB02' within 100 days" # best language
    expr= "from 4 to 7 of ncmp=='4AB02' within 100 days" # best language
    expr= " 4 to 7 events with 4AB02 within 100 days" # best language #events problem ... again format?
    expr= " from 4 to 7 events with 4AB02 within 100 days" # best language #events problem ... again format?

    expr= " at least 5 events with 4AB02 within 100 days" # best language #events problem ... again format?
    expr= " no more than 5 events with 4AB02 in ncmp within 100 days" # best language #events problem ... again format?

    expr= 'min 3 of days>3 within 100 days'

    s.days.rolling('100D').sum()
    s.groupby('pid').days.rolling('100D').sum()

    s.asfreq('D')
    %time count_p(df=df, expr=expr, cols=None, codebook=None, info=None)
    %time count_p(df=df, expr=expr, cols=None, codebook=None, info=info)


    eval_inside(expr)

    """
    left, right = condition.split(' within ')

    # transform ... '4AB02 within 100 days after 4AB04'

    # expr='4AB02 in ncmp within 50 to 100 days before 4AB04 in ncmp'
    if re.match(r'\d+ to \d+ ', right):
        lower, _, upper, unit, direction, *rightsingle = right.split()

        rightsingle = " ".join(rightsingle)
        lower = int(lower)
        upper = int(upper)

        lrows = eval_single(df=df, condition=left, cols=cols, sep=sep, codebook=codebook, info=info, out='rows')
        rrows = eval_single(df=df, condition=rightsingle, cols=cols, sep=sep, codebook=codebook, info=info, out='rows')

        pid_change = ((df.pid - df.pid.shift()) != 0)

        rdates = df.date.where(rrows == 1, np.datetime64('NaT'))
        # rdates2=np.where(rrows==1, df.date, np.datetime64('NaT'))
        # rdates2 = pd.Series(rdates, index=df.index)

        rdates[(pid_change & ~rrows)] = np.datetime64(
            '2100-09-09')  # if not have a date assign one to avoid ffill from person above

        if direction == 'after':
            rdates = rdates.fillna(
                method='ffill')  # hmmm must use groupby here? or can it be avoided? inseret a 999 when pid change and fill it with nan after ffill?

        elif direction == 'before':
            rdates = rdates.fillna(method='bfill')

        rdates = rdates.where(rdates != np.datetime64('2100-09-09'), np.datetime64('NaT'))

        # allow all time units within 5 seconds etc
        if unit == 'days':
            delta = (df['date'] - rdates) / np.timedelta64(1, 'D')
        else:
            # add s if it is not there? like 1 day, 1 second?
            delta = (df['date'] - rdates)
            delta = getattr(delta.dt, unit)

        if direction == 'before':
            delta = delta.abs()
        # add directionn "around" - same as a withon before or within after
        within = (delta >= lower) & (delta <= upper)
        endrows = (lrows & within)

    ### pure within: min 3 of 4AB02 within 100 days

    # 'fewer than 3 of 4AB02 within 100 days'
    # 'cumsum(ddd)>200 within 100 days'

    # 'min 5 of ddd>200 within 100 days'
    # 'min 5 events with ddd>200 within 100 days'
    # 'more than 5 events with ddd>200 within 100 days'

    # 'cumsum(ddd)>200 within 100 days'

    # 'between 3 and 5 of 4AB02 within 100 days'
    # '(4A and 4B) within 100 days'
    # '(4A before 4B) within 100 days'

    # pure within statements have few elements ot the right
    elif len(right.split()) < 3:
        if ' in ' in left:
            word, num, _, codes, _, cols = left.split()
            rows = get_rows(df=df, codes=codes, cols=cols)

        #  'sum(days)>15 within 100 days' or 'min 5 of ddd>200 within 100 days'
        #  expr='sum(days)>15 within 100 days'
        elif re.search('[>=<]', left):
            if 'sum(' in left:
                # may want to create smaller dataframe first, if possible? focus on only some variablve, columns, rows?
                sub = df.set_index('date')  # assume column date exist, should also drop rows with no time
                col, operator = left.split(')')
                col = col.replace('sum(', '').strip(')')
                threshold, unit = right.split()
                if unit == 'days': unit = 'D'
                eval_text = f"(sub.groupby('pid')['{col}'].rolling('{threshold}{unit}').sum()){operator}"
                rows = pd.eval(eval_text, engine='python')
                cpid = rows.any(level=0)
                return cpid
            # 'min 5 of ddd>200 within 100 days'
            else:
                word, num, _, codes = left.split()
                rows = df.eval(codes)  # so far no sumsum etc, only 4 events with sugar_level>10 within 100 days

        # code expression not quantity expression
        else:
            word, num, _, codes = left.split()
            cols = cols
            rows = get_rows(df=df, codes=codes, cols=cols)

        threshold, unit = right.split()
        threshold = int(threshold)
        num = int(num)

        if word == 'max': num = num + 1

        # may need to use expand cols to get the cols (not use cols expression here if it starred)
        sub = df['date'][rows].dropna().to_frame()
        sub['pid'] = sub.index
        sub['shifted_date'] = sub['date'].shift(-num)
        sub['shifted_pid'] = sub['pid'].shift(-num)
        sub['diff_pid'] = (sub['pid'] - sub['shifted_pid'])
        sub['same_pid'] = np.where(sub.diff_pid == 0, 1, 0)

        sub = sub[sub.same_pid == 1]
        # sub['shifted_date'] = sub['date'].groupby(level=0).shift(int(num))
        # sub['shifted_pid'] = sub['pid'].groupby(level=0).shift(int(num))

        # todo: allow for different units here, months, weeks, seconds etc
        sub['diff_days'] = (sub['shifted_date'] - sub['date']) / np.timedelta64(1, 'D')

        # sub[sub.same_pid == 1]['diff_days'].dropna()/np.datetime64(1, 'D')

        if word == 'min':
            endrows = (sub['diff_days'] <= threshold)
            cpid = endrows.any(level=0)

        elif word == 'max':
            # n = df.index.nunique()
            endrows = (sub['diff_days'] <= threshold)
            cpid = ~endrows.any(level=0)

        # revise max and exactly

        elif word == 'exactly':
            endrows = (sub['diff_days'] <= threshold)
            n_max = endrows.groupby(level=0).sum()
            endrows = n_max == threshold
            cpid = endrows.any(level=0)

    #        #todo (also need to change parsing then ...)
    #        elif word=='between':
    #            endrows=(sub['diff_days']<=threshold)
    #            n_max = endrows.groupby(level=0).sum()
    #            endrows = n_max == threshold
    #            cpid = endrows.any(level=0)

    return cpid


# %%

#    expr = '4AB02 before 4AB04'
#
#    expr = '4AB02 in ncmp inside 100 days after 4AB04'
#
#    expr = 'min 2 of 4AB02 in ncmp inside ?[100, 200] days after 4AB04 and days>?[300, 400]'
#
#    expr = '4AB02 in ncmp inside 50 to 100 days after 4AB04'
#
#    expr = '4AB02 in ncmp inside 5th to 8th of 4AB04'
#
#    expr = '4AB02 in ncmp inside first 5 of 4AB04'
#
#    expr = '4AB02 in ncmp inside first 10 events'
#
#    expr = '4AB02 in ncmp inside 10th event'
#    expr = '4AB02 in ncmp before 10th event'
#
#    expr = '4AB02 in ncmp inside first 10% of events'
#
#    expr = '4AB02 in ncmp inside first 10% of 4ab04'
#
#    expr = '4AB02 in ncmp inside 3rd percentile of 4ab04'
#
#    expr = '4AB02 in ncmp inside first 10 days from first event'

# %%


##
## unused (so far)
##

def _create_time_intervals(df, expr, cols=None, sep=None, codebook=None, info=None):
    """
    expr='50 to 100 days before 4AB04 in ncmp'

    create_time_intervals(df=df, expr=expr)
    """
    original = expr

    expr, cols = expr.split(' in ')

    from_day, _, to_day, direction, word, *codeexpr = expr.split()
    codeexpr = " ".join(codeexpr) + ' in ' + cols

    if ' events' in codeexpr:
        df['rdum'] = 1
    else:
        df['rdum'] = eval_single(df=df, condition=codeexpr, cols=cols, out='rows')

    ok_times = _mark_days(df=df, from_day=from_day, direction=direction) & _mark_days(df=df, to_day=to_day,
                                                                                      direction=direction)

    return ok_times


def _mark_days(df, from_day, direction):
    df['last_event'] = np.where(df.rdum == 1, df['date'], np.datetime64('NaT'))
    df['last_event'] = df['last_event'].groupby('pid').fillna(method='ffill')
    df['event_diff'] = df['date'] - df['last_event']  # nb depends on whether we are seeking closes before or after
    df['event_diff'] = df['event_diff'] / np.datetime64(1, 'D')
    df['ok_dates'] = (df.event_diff < from_day)

    return df['okdates']

# %%









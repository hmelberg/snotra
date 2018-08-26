# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:17:33 2018

@author: hmelberg
"""

import numpy as np
import pandas as pd
import re

# %%
from snotra.delete import persons_with, read_code2text, stringify
from snotra.find import find_spikes
from snotra.internal import listify, sniff_sep, get_some_id, fix_cols, fix_codes, fix_args, to_df, get_allcodes, expand_cols, \
    expand_star, expand_hyphen, format_codes, expand_regex, reverse_dict

# %%


def incidence(df, codes=None, cols=None, sep=None, pid='pid',
              date='indate', min_events=1, within_period=None,
              groupby='cohort', update_cohort=True, _fix=True):
    """
    The number of new patients each year who have one or more of the codes


    Args:
        df (dataframe): Dataframe with events, dates and medical codes

        codes (string, list or dict): The codes for the disease

        cols (string, list): Name of cols where codes are located

        pid (str): Name of column with the personal identification number

        date (str): Name of column with the dates for the events
            (the dtype of the column must be datetime)

        min_events (int): Number of events with the codes required for inclusion

        within_period (int): The number of events have to occurr within a period (measured in days) before the person is included
            For instance: min_events=2 and within_period=365 means
            that a person has to have two events within one year
            in order to be included in the calculation of incidence.

        groupby (string): Name of column to group the results by. Typically
            year if you have data for many years (to get incidence each year)

        update_cohort (bool): The cohort a person is assigned to might change
            as the criteria (min_events and within_period) change. If update
            is True, the cohort will be automatically updated to reflect this.
            Exampple: A person may have her first event in 2011 and initially
            be assigned to the 2011 cohort, but have no other events for two
            years. Then in 2014 the person may have five events. If the
            criteria for inclusion is two event in one year, the updated cohort
            for this person will be 2014 and not the original 2011.

        Returns:
            Series (For instance, number of new patients every year)

        Examples:
            df['cohort']=df.groupby('pid').start_date.min().dt.year

            incidence(df, codes=['K50*', 'K51*'], cols=['icdmain', 'icdbi'], sep=',', pid='pid', date='start_date')

            incidence(df, codes={'cd':'K50*', 'uc':'K51*'}, cols=['icdmain', 'icdbi'], sep=',', pid='pid', date='start_date')

        todo:
            make cohort variable redundant ... already have date!
            make it possible to do monthly, quarterly etc incidence?

    """
    sub = df
    incidence_list = []
    namelist = []

    # if an expression instead of a codelist is used as input
    if isinstance(codes, str) and codes.count(' ') > 1:
        b = use_expression(df, codes, cols=cols, sep=sep, out='persons', codebook=codebook, pid=pid)

    if _fix:
        codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep, merge=True, group=False)
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

        cols (string, list): Name of cols where codes are located

        pid (str): Name of column with the personal identification number

        date (str): Name of column with the dates for the events
            (the dtype of the column must be datetime)

        min_events (int): Number of events with the codes required for inclusion

        within_period (int): The number of events have to occurr within a period (measured in days) before the person is included
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
        codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep, merge=True, group=False)
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
    Picks some (randomly selected) individuals and ALL their observations


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
def get_ids(df, codes, cols, groupby, pid='pid', out=None, sep=None):
    codes = listify(codes)
    groupby = listify(groupby)

    codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep)

    rows_with_codes = get_rows(df=df, codes=codes, cols=cols, sep=sep)
    # grouped_ids = df[rows_with_codes].groupby([pid, groupby]).count()
    grouped_ids = df[rows_with_codes].groupby(groupby)[pid].unique()
    grouped_ids = grouped_ids.apply(set)

    return grouped_ids


# %%
def unique_codes(df,
                 cols=None,
                 sep=None,
                 strip=True,
                 name=None,
                 _sniffsep=True,
                 _fix=True):
    """
    Get set of unique values from column(s) with multiple valules in cells

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


    Examples

        to get all unique values in all columns that start with 'drg':
            drg_codes=unique_codes(df=df, cols=['drg*'])
            ncmp = unique_codes(df=df, cols=['ncmpalt'], sep=',')

        all unique atc codes from a column with many comma separated codes
             atc_codes==unique_codes(df=ibd, cols=['atc'], sep=',')

    worry
        numeric columns/content
    """
    if _fix:
        df, cols = to_df(df=df, cols=cols)
        cols = fix_cols(df=df, cols=cols)

    unique_terms = set(pd.unique(df[cols].values.ravel('K')))
    # unique_terms={str(term) for term in unique_terms}

    if _sniffsep and not sep:
        sep = sniff_sep(df, cols)

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
    Expand list of codes with hyphens and star notation to full codes

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

    Example
        get all atc codes that are related to steroids in the atc column:
            codes= ['H02*', 'J01*', 'L04AB02', 'L04AB04']
            codes=expand_codes(df=df, codes=['H02*'], cols='atc')
            codes=expand_codes(df=df, codes=['K25*'], cols='icdmain', sep=',', codebook=codebook)
            codes=expand_codes(df=df, codes=['K25*-K28*'], cols='icdmain', sep=',', codebook=codebook, merge=True)
            codes=expand_codes(df=df, codes=['K25*-K28*'], cols='icdmain', sep=',', codebook=codebook, merge=False)
            codes=expand_codes(df=df, codes=['K25*-K28*'], cols='icdmain', sep=',', codebook=codebook, merge=False, group=True)
            codes=expand_codes(df=df, codes=['K25*-K28*'], cols='icdmain', sep=',', codebook=codebook, merge=True, group=True)
            codes=expand_codes(df=df, codes=['K25*-K28*'], cols='icdmain', sep=',', codebook=codebook, merge=False, group=True)
            codes=expand_codes(df=df, codes=['K25*-K28*'], cols='icdmain', sep=',', codebook=codebook, merge=False, group=False)


            codes=expand_codes(df=df, codes=['K50*', 'K51*'], cols='icdmain', sep=',', codebook=codebook, merge=False, group=True)
            codes=expand_codes(df=df, codes=['K50*', 'K51*'], cols='icdmain', sep=',', codebook=codebook, merge=False, group=False)
            codes=expand_codes(df=df, codes=['K50*', 'K51*'], cols='icdmain', sep=',', codebook=codebook, merge=True, group=False)
            codes=expand_codes(df=df, codes=['K50*', 'K51*'], cols='icdmain', sep=',', codebook=codebook, merge=True, group=True)


        codebook = df.icdmain.unique_codes(sep=',')

        ulcer = 	['K25*-K28*']
        liver = 	'B18* K70.0-K70.3 K70.9 K71.3-K71.5 K71.7 K73* K74* K76.0 K76.2-K76.4 K76.8 K76.9 Z94.4'.split()
        diabetes =	['E10.0', 'E10.l', 'E10.6', 'E10.8', 'E10.9', 'E11.0', 'E11.1', 'E11.6', 'E11.8', 'E11.9', 'E12.0', 'E12.1', 'El2.6', 'E12.8', 'El2.9', 'E13.0', 'E13.1', 'E13.6', 'E13.8', 'E13.9', 'E14.0', 'E14.1', 'E14.6', 'E14.8', 'E14.9']
        hemiplegia = 	 ['G04.1', 'G11.4', 'G80.1', 'G80.2', 'G81*', 'G82*', 'G83.0-G83.4', 'G83.9']

        expand_codes(df=df, codes=ulcer, cols='icdmain', sep=',')
    Note:
        Only codes that actually exist in the cols or the codebook are returned

    """

    codes = listify(codes)

    if codebook:
        unique_words = set(codebook)
    else:
        cols = listify(cols)
        cols = expand_cols(df=df, cols=cols)
        unique_words = set(unique_codes(df=df, cols=cols, sep=sep))

    # if input is not a list of codes, but a dict with categories and codes,
    # then expand each category separately and return the whole dict with
    # expanded codes

    if isinstance(codes, list):
        alist = True
    else:
        alist = False

    codes = format_codes(codes=codes, merge=merge)

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
            codelist = expand_regex(codelist, full_list=unique_words)

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

    example
        get all drg codes starting with 'D':

        d_codes = get_rows(df=df, codes='D*', cols=['drg'])
    """
    # if an expression is used as input
    if isinstance(codes, str) and codes.count(' ') > 1:
        b = use_expression(df, codes, cols=cols, sep=sep, out='rows', codebook=codebook, pid=pid)

    # if a list of codes is used as input
    else:
        if _fix:
            df, cols = to_df(df)
            cols = fix_cols(df=df, cols=cols)
            codes = fix_codes(df=df, codes=codes, cols=cols, sep=sep)

        listify(codes)

        allcodes = get_allcodes(codes)

        if len(allcodes) == 0:  # if no relevant codes --> a column with all false
            b = np.full(len(df), False)
        elif sep:
            allcodes_regex = '|'.join(allcodes)
            b = np.full(len(df), False)
            for col in cols:
                a = df[col].astype(str).str.contains(allcodes_regex, na=False).values
                b = b | a
        # if single value cells only
        else:
            b = df[cols].isin(allcodes).any(axis=1).values

    return b


# %%
def get_pids(df,
             codes,
             cols,
             pid='pid',
             sep=None,
             codebook=None):
    """
    Returns a set pids who have the given codes in the cols

    example

        get pids for all individuals who have icd codes starting with 'C509':

        c509 = get_pids(df=df, codes='C509', cols=['icdmain', 'icdbi'], pid='pid')
    """
    # if an expression instead of a codelist is used as input
    if isinstance(codes, str) and codes.count(' ') > 1:
        selection = use_expression(df, codes, cols=cols, sep=sep, out='persons', codebook=codebook, pid=pid)
        pids = df[selection][pid]
    else:
        pids = get_some_id(df=df, codes=codes, some_id=pid, sep=sep)
    return pids


# %%
def select_persons(df,
                   codes,
                   cols,
                   pid='pid',
                   sep=None):
    """
    Returns a dataframe with all events for people who have the given codes
    example

        c509 = get_pids(df=df, codes='C509', cols=['icdmain', 'icdbi'], pid='pid')
    """
    # if an expression ('K50 and not K51') is used as input
    if isinstance(codes, str) and codes.count(' ') > 1:
        selection = use_expression(df, codes, cols=cols, sep=sep, out='persons', codebook=codebook, pid=pid)
        pids = df[selection]
    # if an list of codes - ['K50', 'K51'] is used as input
    else:
        pids = get_some_id(df=df, codes=codes, some_id=pid, sep=sep)

    df = df[df[pid].isin(pids)]

    return pids



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
        rr.count_persons(df=npr, codes='4AB04', cols='ncmp')

        count_persons(df=df, codes=['4AB*', '4AC*'], cols='ncmp', sep=',', pid='pid')

        count_persons(df=df, codes=['4AB*', '4AC*'], cols='ncmp', sep=',', pid='pid', group=True)
        count_persons(df=df, codes=['4AB*', '4AC*'], cols='ncmp', sep=',', pid='pid', group=True, merge=True)



        count_persons(df=df, codes='4AB*', cols='ncmp', sep=',', pid='pid')
        count_persons(df=df, codes='4AB04', cols='ncmp', sep=',', pid='pid')

        count_persons(df=df, codes={'adaliamumab':'4AB04'}, cols='ncmp', sep=',', pid='pid')
        count_persons(df=df, codes={'adaliamumab':'4AB04'}, cols='ncmp', sep=',', pid='pid')

        npr.count_persons(codes='4AB04', cols='ncmp', groupby=['disease', 'cohort'], sep=',') # works

        npr.groupby(['disease', 'cohort']).apply(count_persons, cols='ncmp', codes='4AB04', sep=',') # works. BIT GET DIFFERENT RESULTS WITH MULTIPLE GROUPBYS!!


        # Counts number of persons for all codes
        count_persons(df=df.ncmp, sep=',', pid='pid') # not work, well it only takes 5 most common .. ajould it take all?
    """

    subset = df

    # if an expression instead of a codelist is used as input
    if isinstance(codes, str) and codes.count(' ') > 1:
        persons = use_expression(df, codes, cols=cols, sep=sep, out='persons', codebook=codebook, pid=pid)
        if normalize:
            counted = persons.sum() / len(persons)
        else:
            counted = persons.sum()


    # codes is a codelist, not an expression
    else:
        if _fix:
            # expands and reformats columns and codes input
            df, cols = to_df(df=df, cols=cols)
            codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep, group=group, merge=merge)
            rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep, _fix=False)
            if not dropna:
                persons = df[pid].nunique()
            subset = df[rows].set_index(pid, drop=False)

        # make a df with the extracted codes
        code_df = extract_codes(df=df, codes=codes, cols=cols, sep=sep, _fix=False, series=False)

        labels = list(code_df.columns)

        counted = pd.Series(index=labels)

        if groupby:
            code_df = code_df.any(level=0)
            sub_plevel = subset.groupby(pid)[groupby].first()
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
    # better name_ maybe eval_persons (person_eval, person_count, person, person ...)
    """
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

    1. pick out every code expression nd expand it? (only star expansion, hyphen would work too? key is to avoud codebook or full lookup )
    2. get codebook (if needed)
    2. use extract coe on each expression

    2. execute logic
    3. return series (bool)

    get_rows_expression(df=npr, expr=expr, cols='icd', sep=',',
                        out='rows', raw=False, regex=False, logic=True,
                        codebook=None, pid='pid', _fix=True):


    """
    cols = listify(cols)
    skipwords = {'and', 'or', 'not'}
    # split_at = {'=', '>', '<'} # what about '>=' '!' '!=' etc
    # well, just use quert to deal with all this. eg if want to examine age>40
    # also additional groupbys ... need to think harder/rewrite/clean up, for now: allow some slightly more complex stuff

    if _fix:
        df, cols = to_df(df, cols)
        if not sep and not ' in:' in expr:
            sep = sniff_sep(df=df, cols=cols)
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
        # remeber, automatic unpacking of tuples so naming works

        del_cols = [col for word, col in word_cols]

        word_cols = [(word, col.replace('in:', '').split(',')) for word, col in word_cols]

        coldf = pd.DataFrame(index=df.index)

        # if no global codebook is specified:
        # this create separate codebook for each column(s) condition
        # background in case some codes overlap (may not be necessary)
        # potential problem: empty codebooks?
        new_codebook = {}
        for word, cols in word_cols:
            sep = sniff_sep(df=df, cols=cols)
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
            # cols=expand_cols(df=df, cols=cols)
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
        # cols=expand_cols(df=df, cols=cols)
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
def search_text(df, text, cols=['text'], select=None, raw=False, regex=False, logic=True, has_underscore=False):
    """
    Searches column(s) in a dataframe for ocurrences of words or phrases

    Can be used to search for occurrences codes that are associated with certain words in the text description.

    Args:
        df (dataframe or series) : The dataframe or series with columns to be searched
        cols (str or list of str): The columns to be searched.
        select (str): row selector for the dataframe. Example: "codetype:'icd' year:2011"
        raw (bool): if True, searches for the raw textstring without modifications
        regex (bool): If True, use regex when searching
        logic (bool): If True, use logical operators in the search (and, or, not in the string)
        underscore (bool): Set to true if the text contains underscores that are important (If set to true, it becomes to search for phrases i.e. two or more words right after each other

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
        3. replace all hele ord med ord==1, men ikke and or not (evnt rereplace if have done it)
        4. create str. contains bool col for rhvert ord
        5. kjør pd eval
    """

    cols = listify(cols)

    df, cols = to_df(df, cols)

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
def first_event(df, codes, cols=None, pid='pid', date='in_date', sep=None):
    """
        Returns time of the first observation for the person based on certain

        Args:
            codes (str or list of str or dict):
                codes marking the event of interest. Star notation is ok.
            pid (string): Patient identifier
            date (string): name of column indicatinf the date of the event

        Returns:
            Pandas series


        Examples:
            first_event(id_col = 'pid', date_col='diagnose_date', groupby = ['disease'], return_as='dict')

            df['cohort'] = df.first_event(id_col = 'pid', date_col='diagnose_date', return_as='dict')

            Date of first registered event with an ibd code for each individual:
                first_event(df=df, codes=['k50*','k51*'], cols='icd', date='date')
    """
    codes = listify(codes)
    cols = listify(cols)

    cols = expand_cols(df=df, cols=cols)
    codes = expand_codes(df=df, codes=codes, cols=cols, sep=sep)

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

    Can produce a set of dummy columns for codes and code groups.
    Can also produce a merged column with only extracted codes.
    Accept star notation.
    Also accepts both single value columns and columns with compund codes and seperators

    out can be: 'text', 'category', 'bool' or 'int'

    example
    to create three dummy columns, based on codes in icdmain column:

    extract_codes(df=df,
              codes={'fracture' : 'S72*', 'cd': 'K50*', 'uc': 'K51*'},
              cols=['icdmain', 'icdbi'],
              merge=False,
              out='text')
    np: problem with extract rows if dataframe is empty (none of the requested codes)
    """
    if _fix:
        df, cols = to_df(df=df, cols=cols)
        codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep, group=group, merge=merge)

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

    years = years_in_row(df, year_col='aar', out='pct')
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

    years_apart(df=df[ibd])

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

def count_codes(df, codes=None, cols=None, sep=None, strip=True,
                ignore_case=False, normalize=False, ascending=False, _fix=True,
                merge=False, group=False, dropna=True):
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

    count_codes(df, codes='K51*', cols='icd', sep=',')
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
        sub, cols = to_df(df=sub, cols=cols)
        cols = fix_cols(df=sub, cols=cols)
        if not sep:
            sep = sniff_sep(df=sub, cols=cols)

        if codes:
            codes = format_codes(codes=codes, merge=merge)
            codes = expand_codes(df=sub, codes=codes, cols=cols, sep=sep, merge=merge, group=group)
            allcodes = get_allcodes(codes)
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
        allcodes = get_allcodes(codes)
        not_included_n = code_count[~code_count.isin(allcodes)].sum()
        code_count = code_count[allcodes]
        if not dropna:
            code_count['na'] = not_included_n

    if isinstance(codes, dict):
        code_count = code_count.rename(index=reverse_dict(codes)).sum(level=0)

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

    codes = listify(codes)
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

    text = listify(text)
    text = [txt.upper().strip('*') for txt in text]
    # codes = " ".join(codes)

    selected_codes = {k: v for k, v in dikt.items() if any(txt in str(v).upper() for txt in text)}

    return selected_codes


# %%
def sankey_format(df, labels, normalize=False, dropna=False, threshold=0.01):
    """

    labels=dict(bio_codes.values())
    import holoviews as hv
    hv.Sankey(t1).options(label_position='left')
    hv.extension('bokeh')
    t4=t1.copy()

    %store t4
    """
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
        ulcer  ,
        liver  	,
        diabetes 	,
        hemiplegia  	,
        renal  	,
        dorgan  ,
        tumr	,
        sliver   ,
        mtumor  	,
        hiv  	]

    disease_cod e s={}
    for i, disease in enumerate(diseases):
        all_cod e s=[]
        disease_s t r=disease_labels[i]
        for code in disease:
            expanded_codes = expand_hyphen(code)
            all_codes.extend(expanded_codes)
        disease_codes[disease_str] = all_codes

    expanded_disease_codes = {}
    no_dot_disease_cod e s={}

    if not dot_notation:
        for disease, codes in disease_codes.items():
            new_codes = [code.replace('. ','') for code in codes]
            no_dot_disease_codes[disease] = new_codes
        disease_codes = no_dot_disease_codes

    all_codes = unique_codes(df=df, cols=cols, sep=sep)

    for disease, codes in disease_codes.items():
        expanded_disease_codes[disease] = expand_codes(df=df, codes=codes, cols=cols, sep=sep, codebook=all_codes)

    codeli s t=[]
    for disease in disease_labels:
        codes = expanded_disease_codes[disease]
        codelist.extend(codes)

    rows = get_rows(df=d f,codes=codelist, cols=cols, sep=sep)

    subset = df[rows]

    charlson_df = persons_with(df=subset, codes=expanded_disease_codes,
                               cols=['icdmain', 'icdbi'], sep=',')

    for disease, point in points.items():
        charlson_df[disease] = charlson_df[disease] * point

    age_poin t s=df.groupby(pid)[age].min().sub(40).div(10).astype(int)
    age_points[age_poin t s< 0 ]=0
    age_points[age_poin t s> 4 ]=4

    disease_points = charlson_df.sum(axis=1).fillna(0)
    charlson_index = age_points.add(disease_points, fill_value=0)

    # make truly missing egual to nans and not zero
    # truly missing = no age available and no icd has been recorded (all nans)

    age_nans = age_points[age_poin t s> 0 ]=0
    icd_nans = df[cols].notnull().sum(axis=1).sum(level=0)
    icd_nans[icd_nans == 0] = np.nan
    icd_nans[icd_nans > 0] = 0

    charlson_with_nans = charlson_index + age_nans + icd_nans

    return charlson_with_nans


# %%


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





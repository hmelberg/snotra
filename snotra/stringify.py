
import pandas as pd
import re
from itertools import chain
from itertools import zip_longest

from .internal import listify, expand_cols, fix_args, to_df
from .core import get_rows, count_codes, extract_codes


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
    codes={
        'L04A*' : 'i',
        'L04AB*' : 'a',
                'H02*' : 'c'}
    pr=pr.set_index('pid_index')
    pr['first_date'] = pr.groupby('pid')['date'].min()
    events=stringify_durations(df=pr, codes=codes, start='date', first_date='first_date', dataset_end_date="01-01-2018")
    events=stringify_durations(df=df, codes=codes, col='ncmpalt', start='start_date', first_date='first', dataset_end_date="01-01-2018")

    df2.columns

    eventstr_duration = stringify_durations(df=pr, codes=codes, cols='atc',
                                           event_duration='ddd',
                                        event_start='date', step=120)
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
        cols = expand_cols(listify(cols))
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False)[:4]

    # fix formatting of input (make list out of a string input and so on)
    codes, cols, old_codes, replace = fix_args(df=df, codes=codes, cols=cols, sep=sep)

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
        collisions: value to be usef if ther is a collision between events in a position


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
        more advanced handeling of collisions
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

    examples

    codes={
        '4AB01': 'e',
        '4AB02' : 'i',
        '4AB04' : 'a',
        '4AB05' : 'x',
        '4AB06' : 'g'}
    medcodes=read_code2text()
    df['diagnosis_date']=df[df.icdmain.fillna('').str.contains('K50|K51')].groupby('pid')['start_date'].min()
    df.columns
    df.start_date

    bio_codes= {
     '4AA23': ('n', 'Natalizumab'),
     '4AA33': ('v', 'Vedolizumab'),
     '4AB02': ('i', 'Infliximab'),
     '4AB04': ('a', 'Adalimumab'),
     '4AB06': ('g', 'Golimumab'),
     '4AC05': ('u', 'Ustekinumab')}

    bio_codes= {'L04AA23': 'n',
     'L04AA33': 'v',
     'L04AB02': 'i',
     'L04AB04': 'a',
     'L04AB06': 'g',
     'L04AC05': 'u'}

    bio_codes={
        'e' : '4AB01',
        'i' : '4AB02',
        'a' : '4AB04'}

    a=stringify_order(
            df=df,
            codes=bio_codes,
            cols='ncmpalt',
            pid='pid',
            event_start='start_date',
            sep=',',
            keep_repeats=True,
            only_unique=False
            )

    codes={
        'L04AB01': 'e',
        'L04AB02' : 'i',
        'L04AB04' : 'a',
        'L04AB05' : 'x',
        'L04AB06' : 'g'}

    bio_rows=get_rows(df=pr, codes=list(codes.keys()), cols='atc')
    pr['first_bio']=pr[bio_rows].groupby('pid')['date'].min()

    a=stringify_order(
            df=pr,
            codes=codes,
            cols='atc',
            pid='pid',
            event_date='date',
            sep=','
            )
    """

    # fix formatting of input
    if _fix:
        df, cols = to_df(df=df, cols=cols)
        codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep)

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

    del_repeats(a)
    """
    no_repeats = str_series.str.replace(r'([a-z])\1+', r'\1')
    return no_repeats


def del_singles(text):
    """
    Deletes single characters from string
    todo: how to deal with first and last position ... delete it too?

    b=del_singles(a)
    (.)\1{2,}

    lookahead \b(?:([a-z])(?!\1))+\b
    lookback ....


    no_singles = str_series.str.replace(r'(.)((?<!\1)&(?!\1))', r'')
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
        event_date: column containing the date for the event
        sep: the seperator used between events if a column has multiple events in a cell
        keep_repeats: identical events after each other are reduced to one (if true)
        only_unique: deletes all events that have occurred previously for the individual (if true)

    Returns:
        series with a string that describes the events for each individual

    Example:
        codes={
        '4AB01': 'e',
        '4AB02' : 'i',
        '4AB04' : 'a',
        '4AB05' : 'x',
        '4AB06' : 'g'}

    df['diagnosis_date']=df[df.icdmain.fillna('').str.contains('K50|K51')].groupby('pid')['start_date'].min()

    a=stringify_order_date(
            df=df,
            codes=codes,
            cols='ncmpalt',
            pid='pid',
            event_date='start_date',
            first_date='diagnosis_date',
            step=90,
            sep=',',
            )

    codes={
        'L04AB01': 'e',
        'L04AB02' : 'i',
        'L04AB04' : 'a',
        'L04AB05' : 'x',
        'L04AB06' : 'g'}

    bio_rows=get_rows(df=pr, codes=list(codes.keys()), cols='atc')
    pr['first_bio']=pr[bio_rows].groupby('pid')['date'].min()

    a=stringify_order_date(
            df=pr,
            codes=codes,
            cols='atc',
            pid='pid',
            event_date='date',
            first_date='first_bio',
            step=90,
            sep=','
            )

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
        cols = expand_cols(listify(cols))
        codes = count_codes(df=df, cols=cols, sep=sep).sort_values(ascending=False)[:4]

    # fix formatting of input (make list out of a string input and so on)
    codes, cols, allcodes, sep = fix_args(df=df, codes=codes, cols=cols, sep=sep)


    # get the rows that contain the relevant codes
    rows = get_rows(df=df, codes=allcodes, cols=cols, sep=sep)
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
                                out='text')

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

    # loop over each code, aggregate strong for each individual, store in df
    for code in codes:
        code_df = code_series[code_series[col].isin([code])]
        stringified = code_df.groupby(pid, sort=False).apply(make_string)
        string_df[code] = stringified

    if merge:
        string_df = interleave_strings(string_df)
    return string_df
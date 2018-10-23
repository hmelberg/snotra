<h1 id="core">core</h1>


Created on Tue Oct  9 17:26:40 2018

@author: hmelberg

<h2 id="core.incidence">incidence</h2>

```python
incidence(df, codes=None, cols=None, sep=None, pid='pid', date='indate', min_events=1, within_period=None, groupby='cohort', update_cohort=True, codebook=None, _fix=True)
```

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


<h2 id="core.make_cohort">make_cohort</h2>

```python
make_cohort(df, codes=None, cols=None, sep=None, pid='pid', date='indate', min_events=1, within_period=None, _fix=True)
```

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

<h2 id="core.sample_persons">sample_persons</h2>

```python
sample_persons(df, pid='pid', n=None, frac=0.1)
```

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



<h2 id="core.unique_codes">unique_codes</h2>

```python
unique_codes(df, cols=None, sep=None, strip=True, _sniffsep=True, _fix=True)
```

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

<h2 id="core.expand_codes">expand_codes</h2>

```python
expand_codes(df=None, codes=None, cols=None, sep=None, codebook=None, hyphen=True, star=True, colon=True, regex=None, del_dot=False, case_sensitive=True, exist=True, merge=False, group=False)
```

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


<h2 id="core.get_rows">get_rows</h2>

```python
get_rows(df, codes, cols, sep=None, codebook=None, _fix=True)
```

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

<h2 id="core.get_pids">get_pids</h2>

```python
get_pids(df, codes, cols=None, pid='pid', sep=None, codebook=None)
```

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

<h2 id="core.select_persons">select_persons</h2>

```python
select_persons(df, codes, cols=None, pid='pid', sep=None, codebook=None)
```

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


<h2 id="core.count_persons">count_persons</h2>

```python
count_persons(df, codes=None, cols=None, pid='pid', sep=None, normalize=False, dropna=True, group=False, merge=False, length=None, groupby=None, codebook=None, _fix=True)
```

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


<h2 id="core.use_expression">use_expression</h2>

```python
use_expression(df, expr, cols=None, sep=None, out='rows', raw=False, regex=False, logic=True, codebook=None, pid='pid', _fix=True)
```

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



<h2 id="core.search_text">search_text</h2>

```python
search_text(df, text, cols='text', select=None, raw=False, regex=False, logic=True, has_underscore=False)
```

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

<h2 id="core.first_event">first_event</h2>

```python
first_event(df, codes, cols=None, pid='pid', date='in_date', sep=None, codebook=None)
```

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

<h2 id="core.extract_codes">extract_codes</h2>

```python
extract_codes(df, codes, cols=None, sep=None, new_sep=',', na_rep='', prefix=None, merge=False, out='bool', _fix=True, series=True, group=False)
```

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

<h2 id="core.years_in_row">years_in_row</h2>

```python
years_in_row(df, year_col='year', groupby=None, info_bank=None, out='pct')
```

average years in row patients are observed, for different start years

Example:
    >>> years = years_in_row(df, year_col='aar', out='pct')

<h2 id="core.years_apart">years_apart</h2>

```python
years_apart(df, pid='pid', year='year')
```

pct of patients with observations that are x years apart

Example:
    >>> years_apart(df=df[ibd])


<h2 id="core.read_codebooks">read_codebooks</h2>

```python
read_codebooks(file=None, select=None, sep=';')
```

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

<h2 id="core.label">label</h2>

```python
label(df, labels=None, select=None, file=None)
```

Translate codes in index to text labels based on content of the dict labels

<h2 id="core.labels_from_codebooks">labels_from_codebooks</h2>

```python
labels_from_codebooks(codebook, select=None, code='code', text='text', only_valid_codes=False)
```

makes a dictionary of code to labels based on two columns in the codebook

<h2 id="core.count_codes">count_codes</h2>

```python
count_codes(df, codes=None, cols=None, sep=None, normalize=False, ascending=False, _fix=True, merge=False, group=False, dropna=True)
```

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



<h2 id="core.lookup_codes">lookup_codes</h2>

```python
lookup_codes(dikt, codes)
```

returns those elements in a dict where key starts with the expressions listed in codes

todo: more complicated star notations: starts with, contains, endswith
lookup(medcodes, 'L04*')


<h2 id="core.get_codes">get_codes</h2>

```python
get_codes(dikt, text)
```

returns those elements in a dict where value contains the expressions listed in codes

todo: more complicated star notations: starts with, contains, endswith
alterative name: find_codes? get_codes?

example
get all codes that have "steroid" in the explanatory text

    get_codes(medcodes, 'steroid*')


<h2 id="core.sankey_format">sankey_format</h2>

```python
sankey_format(df, labels=None, normalize=False, dropna=False, threshold=0.01)
```

Format the dataframe so it is easy fo create a holoviews sankey figure

labels=dict(bio_codes.values())
import holoviews as hv
hv.Sankey(t1).options(label_position='left')
hv.extension('bokeh')
t4=t1.copy()


<h2 id="core.charlson">charlson</h2>

```python
charlson(df, cols='icd', pid='pid', age='age', sep=None, dot_notation=False)
```

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


<h2 id="core.expand_colon">expand_colon</h2>

```python
expand_colon(expr, full_list)
```

Expand expressions with colon notation to a list of complete columns names

expr (str or list): Expression (or list of expressions) to be expanded
full_list (list or array) : The list to slice from

<h2 id="core.expand_star">expand_star</h2>

```python
expand_star(expr, cols=None, full_list=None, sep=None)
```

Expand expressions with star notation to all matching expressions


<h2 id="core.expand_hyphen">expand_hyphen</h2>

```python
expand_hyphen(expr)
```

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


<h2 id="core.persons_with">persons_with</h2>

```python
persons_with(df, codes, cols, pid='pid', sep=None, merge=True, first_date=None, last_date=None, group=False, _fix=True)
```

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

<h2 id="core.stringify_durations">stringify_durations</h2>

```python
stringify_durations(df, codes=None, cols=None, pid='pid', step=120, sep=None, event_start='in_date', event_end=None, event_duration='ddd', first_date=None, last_date=None, censored_date=None, ncodes=None, no_event='-', time_sep='|', merge=True, info=None, report=False)
```

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


<h2 id="core.interleave_strings">interleave_strings</h2>

```python
interleave_strings(df, cols=None, time_sep='|', no_event=' ', agg=False)
```

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


<h2 id="core.left_justify">left_justify</h2>

```python
left_justify(s, fill=' ')
```

after stringify, to make events at same time be in same position
and no, not as crucial as left-pad!

<h2 id="core.overlay_strings">overlay_strings</h2>

```python
overlay_strings(df, cols=None, sep=',', nan='-', collisions='x', interleaved=False)
```

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


<h2 id="core.shorten">shorten</h2>

```python
shorten(events, agg=3, no_event=' ')
```

create a new and shorter string with a longer time step

parameters
    events: (str) string of events that will be aggregated
    agg: (int) the level of aggregation (2=double the step_length, 3=triple)

<h2 id="core.shorten_interleaved">shorten_interleaved</h2>

```python
shorten_interleaved(text, agg=3, time_sep=',', no_event=' ')
```

text="a-si,a--i,a-s-,--si,---i,--s-"

shorten_interleaved(c, agg=2)

the original string must have a distinction between time_sep and no_event_sep
(if not, could try to infer)

<h2 id="core.stringify_order">stringify_order</h2>

```python
stringify_order(df, codes=None, cols=None, pid='pid', event_start='date', sep=None, time_sep='', first_date=None, last_date=None, period=None, keep_repeats=True, only_unique=False, _fix=True)
```

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

<h2 id="core.del_repeats">del_repeats</h2>

```python
del_repeats(str_series)
```

deletes consecutively repeated characters from the strings in a series


<h2 id="core.del_singles">del_singles</h2>

```python
del_singles(text)
```

Deletes single characters from string
todo: how to deal with first and last position ... delete it too?


<h2 id="core.stringify_time">stringify_time</h2>

```python
stringify_time(df, codes=None, cols=None, pid='pid', sep=None, step=90, event_start='date', nfirst=None, first_date=None, last_date=None, censored_date=None, time_sep='|', no_event=' ', collision='*', merge=True, info=None)
```

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

<h2 id="core.count_p">count_p</h2>

```python
count_p(df, expr, cols=None, sep=None, codebook=None, info=None, _use_caching=True)
```


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


<h2 id="core.eval_single">eval_single</h2>

```python
eval_single(df, condition, cols=None, sep=None, codebook=None, out='pid', info=None)
```

evaluates a single expressions (1st 4A),
not relational conditions (A before B, within 100 days after etc)

condition ='first 5 of 4AB02 in ncmp'
condition ='min 2 of days>10'


<h2 id="core.eval_within">eval_within</h2>

```python
eval_within(df, condition, cols=None, sep=None, codebook=None, info=None, out='pid')
```


expr= '4AB02 within 100 days after 4AB04'
expr= 'min 2 of 4AB02 within 100 days'
expr= '4AB02 within 50 to 100 days before 4AB04'
expr= '4AB02 within 50 to 100 days before 4AB04'

__maybe use inside on some?__




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
expr= " 4 to 7 events with 4AB02 within 100 days" # best language `events` problem ... again format?
expr= " from 4 to 7 events with 4AB02 within 100 days" # best language `events` problem ... again format?

expr= " at least 5 events with 4AB02 within 100 days" # best language `events` problem ... again format?
expr= " no more than 5 events with 4AB02 in ncmp within 100 days" # best language `events` problem ... again format?

expr= 'min 3 of days>3 within 100 days'

s.days.rolling('100D').sum()
s.groupby('pid').days.rolling('100D').sum()

s.asfreq('D')
%time count_p(df=df, expr=expr, cols=None, codebook=None, info=None)
%time count_p(df=df, expr=expr, cols=None, codebook=None, info=info)


eval_inside(expr)



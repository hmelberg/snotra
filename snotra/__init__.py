import pandas as pd
import os

from .core import *

# todo: make a proper config file
_PATH = os.path.abspath(core.__file__)

# %% monkeypatch to the functions become methods on the dataframe
# could use decorators/pandas_flavor
# pandas_flavor: good, but want to eliminate dependencies
# approach below may be bloated and harder to maintain
# (must change when method names change)

series_methods = [sample_persons, count_persons, unique_codes, extract_codes,
                  count_codes, label, use_expression, select_persons]

frame_methods = [sample_persons, first_event, get_pids, unique_codes,
                 expand_codes, get_rows, count_persons, extract_codes,
                 count_codes, label, use_expression, select_persons, count_p]

# probably a horrible way of doing something horrible!
for method in frame_methods:
    setattr(pd.DataFrame, getattr(method, "__name__"), method)

for method in series_methods:
    setattr(pd.Series, getattr(method, "__name__"), method)


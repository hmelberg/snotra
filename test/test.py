from io import StringIO

from snotra.core import make_cohort


def test_make_cohort():
    a = """
        icdmain,icdbi,date,pid
        K50,,20.01.2018,1
        K50,,17.03.2018,1
        K51,,12.05.2018,1
        K50,	K51	,07.07.2018,	1
        K51,	K50,01.09.2018,	1
        K51,,27.10.2018,1
        K50,,22.12.2018,1
        K51,	K50	,01.01.2016,	2
        K50,,05.01.2017,2
        K51,,08.06.2017,2
        K50,,04.02.2016,3
        """

    a = StringIO(a)

    assert len(make_cohort(df=df, codes={'ibd': ['K50*', 'K51*']},
                           cols=['icdmain', 'icdbi'],
                           sep=',',
                           pid='pid',
                           date='date')) == 3

    assert len(make_cohort(df=df, codes={'ibd': ['K50*', 'K51*']},
                           cols=['icdmain', 'icdbi'],
                           min_events=2,
                           within_period=100,
                           sep=',',
                           pid='pid',
                           date='date')) == 1
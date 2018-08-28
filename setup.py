from setuptools import setup

setup(
    name='snotra',
    version='0.0.8',
    packages=['snotra', ],
    license='MIT',
    include_package_data=True,
    author="Hans Olav Melberg",
    author_email="hans.melberg@gmail.com",
    description="Tool to analyze data on hospital events, prescriptions and similar types of health data",
    keywords="health pandas",
    long_description='''
    **Snotra - Health registry research using Pandas and Python**
    
    Snotra is a tool that extends and builds on the Pandas library to make it easier to analyze data on hospital events, prescriptions and similar types of health data.

    Example: 
    To count the number of persons with a diagnosis starting with 'K50' in a dataframe with many events per person:
        df.count_persons('K50*', cols='icd')
    Snotra is also a Norse goddess associated with the Norse word for clever. Although the library is related to health research, the Norse term is in no way connected to snot''',
)
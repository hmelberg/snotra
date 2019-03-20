# Snotra - Health registry research using Pandas and Python
Snotra is a tool that extends and builds on the Pandas library to make it easier to analyze data on hospital events, prescriptions and similar types of health data.

Snotra is a also a Norse goddess associated with wisdom. 

## Examples
- **Count the number of unique persons with a diagnosis in event data**
    - Using special notation (star, hyphen and colon)
    ```python
    df.count_persons(codes=['K50*', 'K51*'], cols='icd*')  
    
    df.count_persons(codes=['K50.0-K51.9'], cols='icd*')        
    ```
  
    - Using logical expressions
    ```python
    df.count_persons(codes='(K50 or K51) and not K52', cols='icd*')
      
    df.count_persons(codes='K50 in: icd and 4AB02 in:atc1, atc2')
    ```
    
    - Using temporal expressions
    ```python
    df.count_persons('K50 before K51', cols='icd')
    
    df.count_persons('4AB02 within 365 days after 4AB04', cols='atc')
    ```
     
    - Commplex expressions
    ```python 
    df.count_persons('min 5 of 4AB02 before 4AB04', cols='atc')
    ```

    
- **Select all events for some persons using codelists or logical expressions**

    ```python
    df.select_persons(codes=['K50*', 'K51*'], cols='icd')

    df.select_persons(codes='(K50 or K51) and not K52', cols='icd')

    ```

- **Count the number of unique codes in multiple columns with multiple values in each cell**
    ```python
    df['icd_primary', 'icd_secondary'].count_codes(sep=',')
    ```
    
- **Use full text labels instead of medical codes in the output**
    ```python
    df['icd_primary'].count_codes().labels('icd')
    ``` 
    
- **Calculate Charlson Comorbidity Index***
    ```python
    cci = sa.charlson(df=df, cols=['icd1', 'icd2'], sep=',')
    ```

- **Find typical treatment patterns***
    ```python
    codes = {'i': '4AB02', 'a': '4AB04']
    df.stringify_order(codes=codes, cols='icd', sep=',', keep_repeats=False).value_counts()
    ```

## Installation
 - We recommend using 'sa' as an abbreviation for snotra
    
    ```python
    pip install snotra as sa
    ```
    
 ## Requirements
 - Python 3.6 
 - Pandas

 ## License
 MIT
 
 ## Documentation
 - See [doc](doc/docs.md) for a draft overview of methods and functions. 
 - See [wiki](https://github.com/hmelberg/snotra/wiki) for more in depth discussions of some topics.
 
 ## Disclaimer
 Snotra is currently under development and not ready for production. Much remains to be tested and corrected, use at your own risk - and contributions are welcome!
  
 
 ## Features
 - **Easy and efficient notation and methods to deal with medical codes** 
 Medical data often use special code systems to indicate diagnoses, pharmaceuticals and medical procedures. We integrated these tools and allow the use of different types of notation (star, hyphen, colon) to make it easy to select or count the relevant patients.

- **Answer person level question using event level data** 
Often health data contains information about events (a prescription, a hospital stay), while the questions we want answered are both at the event-level and person-level:
    - Event-level: How many doses of a certain pharmaceutical is used in a year?
    - Person-level: How many people have received a given pharmaceutical?
 We have methods, such as `count_persons` that make it easy to get person-level answers from event-level data.

- **Deal with messy data** 
Sometimes the files supplied to the analysis are multiple large files of messy administrative data. For instance procedure codes can be merged in one column (comma separated) or spread across many columns. To deal with this we have methods that accept both types of data. For instance: the method `count_codes()` can count codes from many columns, some of which may contain comma seperated codes, some of which may be single valued.

from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()



setup(
    name='snotra',
    version='0.0.9',
    packages=['snotra', ],
    license='MIT',
    include_package_data=True,
    author="Hans Olav Melberg",
    author_email="hans.melberg@gmail.com",
    description="Tool to analyze data on hospital events, prescriptions and similar types of health data",
    keywords="health pandas",
    long_description=long_description,
    long_description_content_type='text/markdown'
)
#!/usr/local/bin
"""
setup.py - script for building MyApplication

Usage:
    % python setup.py py2app
"""
from setuptools import setup
import os, sys

if os.name == 'posix':
    setup(
        name='msatCommand',
        version='0.4.4',
        description='python searching of fasta files for microsat repeats',
        author='Brant C. Faircloth',
        author_email='brant@uga.edu',
        license='GPL',
        app=['msatCommandGui.py'],
        setup_requires=["py2app"],
        options=dict(py2app=dict(
                packages=['Bio']
            ))
    )
    
elif os.name == 'nt':
    import py2exe
    opts= {
        "py2exe":{
            "packages":"Bio",
            "excludes":"Tkinter"
            }
        }
    setup(
        options=opts,
        name='msatCommand',
        version='0.4.5',
        description='python searching of fasta files for microsat repeats',
        author='Brant C. Faircloth',
        author_email='brant@uga.edu',
        license='GPL',
        zipfile=None,
        windows=[
            {
                'script':'msatCommandGui.py'
            }
        ],
    )
    
else:
    sys.exit

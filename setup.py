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
        name='msatcommander',
        version='0.8.1',
        description='python searching of fasta files for microsat repeats',
        author='Brant C. Faircloth',
        author_email='brant@uga.edu',
        license='GPL',
        app=['msatcommander.py'],
        setup_requires=["py2app"],
        options=dict(py2app=dict(
                packages=['Bio'],
                resources=['primer3-1.1.1/src/primer3_core',
                'mscGui.xrc', 
                'mscfunc.py',
                'mscprimer.py',
                'mscprimertag.py',
                'repeatClasses.py',
                'primer.py']
            ))
    )
    
elif os.name == 'nt':
    from distutils.core import setup
    import py2exe
    opts= {
        "py2exe":{
            "packages":"Bio",
            "excludes":"Tkinter"
            }
        }
    setup(
        options=opts,
        data_files=[("",["mscGui.xrc","mscfunc.py","mscprimer.py","mscprimertag.py","repeatClasses.py","primer.py",])],
        name='msatcommander',
        version='0.8.1',
        description='python searching of fasta files for microsat repeats',
        author='Brant C. Faircloth',
        author_email='brant@uga.edu',
        license='GPL',
        zipfile=None,
        windows=[
            {
                'script':'msatcommander.py'
            }
        ],
    )
    
else:
    sys.exit

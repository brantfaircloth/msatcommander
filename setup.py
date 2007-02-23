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
        version='0.4.2',
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
    setup(
        name='msatCommand',
        version='0.4.2',
        description='python searching of fasta files for microsat repeats',
        author='Brant C. Faircloth',
        author_email='brant@uga.edu',
        license='GPL',
        windows=[
            {
                'script':'msatCommandGui.py',
                'icon_resources':[(1, 'gmconvert_icon_files2.ico')]
            }
        ],
    )
    
else:
    sys.exit
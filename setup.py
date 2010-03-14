#!/usr/local/bin
"""
setup.py - script for building MyApplication

Usage:
    % python setup.py py2app
"""
from setuptools import setup
import os, sys, shutil

# add sys.path = [os.path.join(os.environ['RESOURCEPATH'], 'lib', 'python2.6', 'lib-dynload')] + sys.path
# to dist/msatcommander.app/Contents/Resources/__boot__.py

# remove build directories
print 'Deleting build dirs...'
shutil.rmtree('./build')
shutil.rmtree('./dist')
print 'Building app...'
if os.name == 'posix':
    setup(
        name='msatcommander',
        version='0.9.0',
        description='python searching of fasta files for microsat repeats',
        author='Brant C. Faircloth',
        author_email='faircloth@gmail.com',
        license='GPL',
        app=['main.py'],
        setup_requires=["py2app"],
        options=dict(py2app=dict(
                includes=['Bio.SeqIO', 'p3wrapr', 'PyQt4', 'PyQt4.QtCore', 'PyQt4.QtGui', 'sip'],
                packages=[],
                resources=['primer3_core',
                'misprime_lib_weight'],
                excludes=['p3wrapr/.git', 'Bio.nexus', 'Scipy', 'numpy', 'PyQt4.QtDesigner', 'PyQt4.QtNetwork', 'PyQt4.QtOpenGL', 'PyQt4.QtScript', 'PyQt4.QtSql', 'PyQt4.QtTest', 'PyQt4.QtWebKit', 'PyQt4.QtXml', 'PyQt4.phonon']
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

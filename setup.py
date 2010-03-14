#!/usr/local/bin
"""
setup.py - script for building MyApplication

Usage:
    % python setup.py py2app
"""
from setuptools import setup
import os, sys, shutil


# in order to add sys.path = [os.path.join(os.environ['RESOURCEPATH'], 'lib', 'python2.6', 'lib-dynload')] + sys.path
# to dist/msatcommander.app/Contents/Resources/__boot__.py, i hacked
# /Library/Python/2.6/site-packages/py2app-0.4.4-py2.6.egg/py2app @ line 1996
# to:
#
# data = self.get_bootstrap_data(fn)
# if os.path.split(fn)[-1] == 'boot_app.py':
#     data = data.split('\n')
#     temp = data[:3]
#     ver = sys.version_info[:2]
#     ver = 'python%s.%s' % (ver[0], ver[1])
#     pathinfo = '''    sys.path = [os.path.join(os.environ['RESOURCEPATH'], 'lib', '%s', 'lib-dynload')] + sys.path''' % ver
#     temp.append(pathinfo)
#     temp += data[3:]
#     data = '\n'.join(temp)
# bootfile.write(data)

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
                includes=['Bio.SeqIO',
                            'p3wrapr',
                            'PyQt4',
                            'PyQt4.QtCore',
                            'PyQt4.QtGui',
                            'sip'],
                packages=[],
                resources=['primer3_core',
                            'misprime_lib_weight',
                            'primer3_config'],
                excludes=['/Users/bcf/git/brant/modules/p3wrapr/docs/',
                            '/Users/bcf/git/brant/modules/p3wrapr/.git',
                            'Bio.nexus',
                            'Scipy',
                            'numpy',
                            'PyQt4.QtDesigner',
                            'PyQt4.QtNetwork',
                            'PyQt4.QtOpenGL',
                            'PyQt4.QtScript',
                            'PyQt4.QtSql',
                            'PyQt4.QtTest',
                            'PyQt4.QtWebKit',
                            'PyQt4.QtXml',
                            'PyQt4.phonon']
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

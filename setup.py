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
try:
    shutil.rmtree('./build')
except:
    pass
try:
    shutil.rmtree('./dist') 
except:
    pass
print 'Building app...'

VER = '1.0.4-beta'
NAME = 'msatcommander'
APPNAME = '%s-%s' % (NAME, VER)
FULLAPPNAME = '%s.app' % APPNAME
AUTHOR = 'Brant Faircloth'
YEAR = 2011

PLIST = dict(CFBundleName = APPNAME,
             CFBundleShortVersionString = VER,
             CFBundleGetInfoString = NAME + " " + VER,
             CFBundleExecutable = NAME + " " + VER,
             CFBundleIdentifier = "org.brantfaircloth.%s" % APPNAME,
             CFBundleDevelopmentRegion = 'English',
             NSHumanReadableCopyright = u"Copyright %s 2008-%d" % (AUTHOR, YEAR),
             CFBundleIconFile='mc.icns'
             )

if os.name == 'posix':
    DIST    = os.path.join(os.getcwd(),'dist')
    INNARDS = os.path.join(DIST, FULLAPPNAME, 'Contents')
    setup(
        name=APPNAME,
        version=VER,
        description='python searching of fasta files for microsat repeats',
        author='Brant C. Faircloth',
        author_email='faircloth@gmail.com',
        license='GPL',
        app=['main.py'],
        setup_requires=["py2app"],
	data_files=[
            ('msat', ['msat/__init__.py','msat/motif.py','msat/seqsearch.py','msat/finder.py']),
            ('p3wrapr', ['p3wrapr/__init__.py','p3wrapr/primer.py']),
        ],
        options=dict(py2app=dict(
                includes=['Bio.SeqIO',
                            'PyQt4',
                            'PyQt4.QtCore',
                            'PyQt4.QtGui',
                            'sip',
                            ],
                plist = PLIST,
                packages=[],
                iconfile='icon/mc.icns',
                resources=[ 'qt.conf',
                            'msatcommander.conf',
                            'misprime_lib_weight',
                            '../mod-primer3/src/primer3_long',
                            '../mod-primer3/src/primer3_config',
                            ],
                excludes=[  'Bio.Genbank',
                            'modules.p3wrapr',
                            'modules.msat',
                            'Scipy',
                            'numpy',
                            'PyQt4.QtCore_debug',
                            'PyQt4.QtGui_debug',
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
    ## post-flight
    # make primer3 core executable
    print '\n\nSetting primer3_long to executable...'
    os.chmod(os.path.join(INNARDS, 'Resources','primer3_long'), 0755)
    print 'Removing QtCore_debug and QtGui_debug...'
    for p in [  ( "Frameworks", "QtGui.framework", "Versions", "4", "QtGui_debug" ),
                ( "Frameworks", "QtGui.framework", "QtGui_debug"),
                ( "Frameworks", "QtCore.framework", "Versions", "4", "QtCore_debug" ),
                ( "Frameworks", "QtCore.framework", "QtCore_debug" )]:
        db_path = os.path.join(INNARDS, *p)
        os.system(' '.join([ 'rm -rf', db_path]))
    
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

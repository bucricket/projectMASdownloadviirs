#!/usr/bin/env python

from __future__ import print_function
import subprocess
import os
import distutils

# set project base directory structure
base = os.getcwd()
    
try:
    from setuptools import setup
    setup_kwargs = {'entry_points': {'console_scripts':['processviirsdownload=processviirsdownload.downloadData:main']}}
except ImportError:
    from distutils.core import setup
    setup_kwargs = {'scripts': ['bin/processviirsdownload']}
    
from processviirsdownload import __version__

#=====build DMS binaries===============================
# get Anaconda root location
p = subprocess.Popen(["conda", "info", "--root"],stdout=subprocess.PIPE)
out = p.communicate()
condaPath = out[0][:-1]
    
prefix  = os.environ.get('PREFIX')
processDi = os.path.abspath(os.path.join(prefix,os.pardir))
processDir = os.path.join(processDi,'work')
libEnv = os.path.join(prefix,'lib')
binEnv = os.path.join(prefix,'bin')
libDir = os.path.join(processDir,'source','lib')
binDir = os.path.join(processDir,'source','bin')
distutils.dir_util.copy_tree(binDir,binEnv)


setup(
    name="projectmasviirsdownload",
    version=__version__,
    description="downloads VIIRS data for ALEXI",
    author="Mitchell Schull",
    author_email="mitch.schull@noaa.gov",
    url="https://github.com/bucricket/projectMASdownloadviirs.git",
#    py_modules=['processviirs'],
    packages= ['processviirsdownload'],
    platforms='Posix; MacOS X; Windows',
    license='BSD 3-Clause',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        # Uses dictionary comprehensions ==> 2.7 only
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: GIS',
    ],  
    **setup_kwargs
)


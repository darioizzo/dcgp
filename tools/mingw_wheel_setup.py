from setuptools import setup
from setuptools.dist import Distribution
from distutils import util
import sys

NAME = 'pyaudi'
VERSION = '@audi_VERSION@'
DESCRIPTION = 'Implementation of a high-order automated differentiation system using generalized dual numbers. Implementation of a differential algebra.'
LONG_DESCRIPTION = 'Implementation of a high-order automated differentiation system using generalized dual numbers. Implementation of a differential algebra.'
URL = 'https://github.com/darioizzo/audi'
AUTHOR = 'Dario Izzo'
AUTHOR_EMAIL = 'dario.izzo@gmail.com'
LICENSE = 'GPLv3+/LGPL3+'
CLASSIFIERS = [
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 4 - Beta',

    'Operating System :: OS Independent',

    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Physics',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',

    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',

    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3'
]
KEYWORDS = 'differential algebra taylor polynomials automatic differentiation'
PLATFORMS = ['Unix','Windows','OSX']

class BinaryDistribution(Distribution):
    def has_ext_modules(foo):
        return True

# Setup the list of external dlls.
import os.path
mingw_wheel_libs = 'mingw_wheel_libs_python{}.txt'.format(sys.version_info[0])
l = open(mingw_wheel_libs,'r').readlines()
DLL_LIST = [os.path.basename(_[:-1]) for _ in l]

setup(name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    url=URL,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    keywords=KEYWORDS,
    platforms=PLATFORMS,
    packages=['pyaudi'],
    # Include pre-compiled extension
    package_data={
                'pyaudi': ['_core.pyd'] + DLL_LIST
                },
    distclass=BinaryDistribution)

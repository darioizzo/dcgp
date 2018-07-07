from setuptools import setup
from setuptools.dist import Distribution
from distutils import util
import sys

NAME = 'dcgpy'
VERSION = '@dCGP_VERSION@'
DESCRIPTION = 'Implementation of differentiable Cartesian Genetic Programming (d-CGP).'
LONG_DESCRIPTION = 'The d-CGP is a recent development in the field of Genetic Programming that builds upon Cartesian Genetic Programming adding the information about the any-order derivatives of encoded program using a differential algebra.'
URL = 'https://github.com/darioizzo/d-CGP'
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
KEYWORDS = 'cartesian genetic programming backpropagation machine learning'
INSTALL_REQUIRES = ['pyaudi']
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
    install_requires=INSTALL_REQUIRES,
    packages=['dcgpy'],
    # Include pre-compiled extension
    package_data={
                'dcgpy': ['_core.pyd'] + DLL_LIST
                },
    distclass=BinaryDistribution)

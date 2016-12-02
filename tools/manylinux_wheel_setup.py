from distutils.core import setup, Extension
NAME = 'pyaudi'
VERSION = '@audi_VERSION@'
DESCRIPTION = 'Implementation of a high-order automated differentiation system using generalized dual numbers. Implementation of a differential algebra.'
LONG_DESCRIPTION = 'Implementation of a high-order automated differentiation system using generalized dual numbers. Implementation of a differential algebra.'
URL = 'https://github.com/darioizzo/audi'
AUTHOR = 'Dario Izzo'
AUTHOR_EMAIL = 'dario.izzo@gmail.com'
LICENSE = 'GPLv3+/LGPL3+'
INSTALL_REQUIRES = ['numpy']
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

extension_module = Extension(
    'dummy',
     sources=['dummy.cpp']
)

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
    ext_modules = [extension_module],
    packages=['pyaudi'],
    # Include pre-compiled extension
    package_data={
               	'pyaudi': ['_core.so']
               	},
)

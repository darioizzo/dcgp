import os
import re
import sys


def wget(url, out):
    import urllib.request
    print('Downloading "' + url + '" as "' + out + '"')
    urllib.request.urlretrieve(url, out)


def rm_fr(path):
    import os
    import shutil
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def run_command(raw_command, directory=None, verbose=True):
    # Helper function to run a command and display optionally its output
    # unbuffered.
    import shlex
    import sys
    from subprocess import Popen, PIPE, STDOUT
    print(raw_command)
    proc = Popen(shlex.split(raw_command), cwd=directory,
                 stdout=PIPE, stderr=STDOUT)
    if verbose:
        output = ''
        while True:
            line = proc.stdout.readline()
            if not line:
                break
            line = str(line, 'utf-8')
            # Don't print the newline character.
            print(line[:-1])
            sys.stdout.flush()
            output += line
        proc.communicate()
    else:
        output = str(proc.communicate()[0], 'utf-8')
    if proc.returncode:
        raise RuntimeError(output)
    return output


# Build type setup.
BUILD_TYPE = os.environ['BUILD_TYPE']
is_release_build = (os.environ['APPVEYOR_REPO_TAG'] == 'true') and bool(
    re.match(r'v[0-9]+\.[0-9]+.*', os.environ['APPVEYOR_REPO_TAG_NAME']))
if is_release_build:
    print("Release build detected, tag is '" +
          os.environ['APPVEYOR_REPO_TAG_NAME'] + "'")
is_python_build = 'Python' in BUILD_TYPE

# Get mingw and set the path.
wget(r'https://github.com/bluescarni/binary_deps/raw/master/x86_64-6.2.0-release-posix-seh-rt_v5-rev1.7z', 'mw64.7z')
run_command(r'7z x -oC:\\ mw64.7z', verbose=False)
ORIGINAL_PATH = os.environ['PATH']
os.environ['PATH'] = r'C:\\mingw64\\bin;' + os.environ['PATH']

# Download common deps.
wget(r'https://github.com/bluescarni/binary_deps/raw/master/gmp_mingw_64.7z', 'gmp.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/mpfr_mingw_64.7z', 'mpfr.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/boost_mingw_64.7z', 'boost.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/eigen3.7z', 'eigen3.7z')

# Extract them.
run_command(r'7z x -aoa -oC:\\ gmp.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ mpfr.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ boost.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ eigen3.7z', verbose=False)

# Build the Thread Building Block libraries 
wget(r'https://github.com/01org/tbb/archive/2019.zip', 'tbb-2019.zip')
run_command(r'unzip tbb-2019.zip', verbose=False)
os.chdir('tbb-2019')
os.chdir('build')
run_command(r'generate_tbbvars.bat', verbose=False)
run_command(r'tbbvars.bat', verbose=False)
os.chdir('../')
run_command(r'mingw32-make compiler=gcc VERBOSE=1', verbose=True)
# Install the TBB libraries
run_command(r'cp -r build/windows_intel64_gcc_mingw6.2.0_release/tbb* /local/lib/', verbose=True)
run_command(r'cp -r include/tbb /local/include/tbb', verbose=True)

#windows_intel64_gcc_mingw6.2.0_debug
#windows_intel64_gcc_mingw6.2.0_release
#c:\local\lib\libboost_python-mgw62-mt-1_63.dll

os.chdir('../')

# Download piranha 0.11 https://github.com/bluescarni/piranha/archive/v0.11.zip
wget(r'https://github.com/bluescarni/piranha/archive/v0.11.zip', 'piranhav11.zip')
run_command(r'unzip piranhav11.zip', verbose=False)
# Move to the directory created and make piranha install its headers
os.chdir('piranha-0.11')
os.makedirs('build')
os.chdir('build')
print("Installing piranha")
run_command(
    r'cmake -G "MinGW Makefiles" .. -DCMAKE_INSTALL_PREFIX=c:\\local ', verbose=False)
run_command(r'mingw32-make install VERBOSE=1', verbose=False)
os.chdir('../../')
print("Piranha sucessfully installed .. continuing")

# Download audi 1.5 https://github.com/darioizzo/audi/archive/v1.5.zip
wget(r'https://github.com/darioizzo/audi/archive/v1.5.zip', 'audiv15.zip')
run_command(r'unzip audiv15.zip', verbose=False)
# Move to the directory created and make audi install its headers
os.chdir('audi-1.5')
os.makedirs('build')
os.chdir('build')
print("Installing audi")
run_command(r'cmake -G "MinGW Makefiles" .. -DAUDI_BUILD_AUDI=yes -DAUDI_BUILD_PYAUDI=no -DAUDI_BUILD_TEST=no -DAUDI_WITH_MPPP=no -DCMAKE_INSTALL_PREFIX=c:\\local ', verbose=False)
run_command(r'mingw32-make install VERBOSE=1', verbose=False)
os.chdir('../../')
print("Audi sucessfully installed .. continuing")

# Setup of the dependencies for a Python build.
if is_python_build:
    if BUILD_TYPE == 'Python34':
        python_version = '34'
    elif BUILD_TYPE == 'Python35':
        python_version = '35'
    elif BUILD_TYPE == 'Python36':
        python_version = '36'
    elif BUILD_TYPE == 'Python27':
        python_version = '27'
    else:
        raise RuntimeError('Unsupported Python build: ' + BUILD_TYPE)
    python_package = r'python' + python_version + r'_mingw_64.7z'
    boost_python_package = r'boost_python_' + python_version + r'_mingw_64.7z'
    # Remove any existing Python installation.
    rm_fr(r'c:\\Python' + python_version)
    # Set paths.
    pinterp = r'c:\\Python' + python_version + r'\\python.exe'
    pip = r'c:\\Python' + python_version + r'\\scripts\\pip'
    twine = r'c:\\Python' + python_version + r'\\scripts\\twine'
    dcgpy_install_path = r'C:\\Python' + \
        python_version + r'\\Lib\\site-packages\\dcgpy'
    # Get Python.
    wget(r'https://github.com/bluescarni/binary_deps/raw/master/' +
         python_package, 'python.7z')
    run_command(r'7z x -aoa -oC:\\ python.7z', verbose=False)
    # Get Boost Python.
    wget(r'https://github.com/bluescarni/binary_deps/raw/master/' +
         boost_python_package, 'boost_python.7z')
    run_command(r'7z x -aoa -oC:\\ boost_python.7z', verbose=False)
    # Install pip and deps.
    wget(r'https://bootstrap.pypa.io/get-pip.py', 'get-pip.py')
    run_command(pinterp + ' get-pip.py')
    run_command(pinterp + r' -m pip install pyaudi')

# Set the path so that the precompiled libs can be found.
os.environ['PATH'] = os.environ['PATH'] + r';c:\\local\\lib'

# Proceed to the build.
print("Start Building dcgp(y)")
common_cmake_opts = r'-DCMAKE_PREFIX_PATH=c:\\local -DCMAKE_INSTALL_PREFIX=c:\\local'
if is_python_build:
    os.makedirs('build_dcgp')
    os.chdir('build_dcgp')
    run_command(
        r'cmake -G "MinGW Makefiles" ..  -DCMAKE_BUILD_TYPE=Release -DDCGP_BUILD_TESTS=no -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_DCGPY=no' + ' ' + common_cmake_opts)
    run_command(r'mingw32-make install VERBOSE=1 -j2')
    os.chdir('..')
    os.makedirs('build_dcgpy')
    os.chdir('build_dcgpy')
    run_command(r'cmake -G "MinGW Makefiles" ..  -DDCGPY_INSTALL_PATH=c:\\local -DDCGP_BUILD_DCGP=no -DDCGP_BUILD_DCGPY=yes -DCMAKE_BUILD_TYPE=Release ' + common_cmake_opts + r' -DBoost_PYTHON_LIBRARY_RELEASE=c:\\local\\lib\\libboost_python' +
                (python_version[0] if python_version[0] == '3' else r'') + r'-mgw62-mt-1_63.dll  -DPYTHON_EXECUTABLE=C:\\Python' + python_version + r'\\python.exe -DPYTHON_LIBRARY=C:\\Python' + python_version + r'\\libs\\python' + python_version + r'.dll')
    run_command(r'mingw32-make install VERBOSE=1 -j2')
elif BUILD_TYPE in ['Release', 'Debug']:
    os.makedirs('build_dcgp')
    os.chdir('build_dcgp')
    cmake_opts = r'-DCMAKE_BUILD_TYPE=' + BUILD_TYPE + \
        r' -DDCGP_BUILD_TESTS=yes ' + common_cmake_opts
    run_command(r'cmake -G "MinGW Makefiles" .. ' + cmake_opts)
    run_command(r'mingw32-make install VERBOSE=1 -j2')
    run_command(r'ctest')
else:
    raise RuntimeError('Unsupported build type: ' + BUILD_TYPE)

# Packaging.
if is_python_build:
    # Run the Python tests.
    run_command(
        pinterp + r' -c "from dcgpy import test; test.run_test_suite()"')
    # Build the wheel.
    import shutil
    os.chdir('wheel')
    shutil.move(dcgpy_install_path, r'.')
    wheel_libs = 'mingw_wheel_libs_python{}.txt'.format(python_version[0])
    DLL_LIST = [_[:-1] for _ in open(wheel_libs, 'r').readlines()]
    for _ in DLL_LIST:
        shutil.copy(_, 'dcgpy')
    run_command(pinterp + r' setup.py bdist_wheel')
    os.environ['PATH'] = ORIGINAL_PATH
    # workaround necessary to be able to call pip install via python (next line) and
    # not find a dcgpy already in the pythonpath
    os.makedirs('garbage')
    shutil.move('dcgpy', r'garbage')
    shutil.move('dcgpy.egg-info', r'garbage')
    # call pip via python, workaround to avoid path issues when calling pip from win
    # (https://github.com/pypa/pip/issues/1997)
    run_command(pinterp + r' -m pip install dist\\' + os.listdir('dist')[0])
    run_command(
        pinterp + r' -c "from dcgpy import test; test.run_test_suite()"', directory=r'c:\\')
    if is_release_build:
        run_command(pinterp + r' -m pip install twine')
        run_command(twine + r' upload -u darioizzo dist\\' +
                    os.listdir('dist')[0])

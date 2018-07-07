import os
import re


def wget(url, out):
    import urllib.request
    print('Downloading "' + url + '" as "' + out + '"')
    urllib.request.urlretrieve(url, out)


def rm_fr(path):
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

# Get mingw and set the path.
wget(r'https://github.com/bluescarni/binary_deps/raw/master/x86_64-6.2.0-release-posix-seh-rt_v5-rev1.7z', 'mw64.7z')
run_command(r'7z x -oC:\\ mw64.7z', verbose=False)
ORIGINAL_PATH = os.environ['PATH']
os.environ['PATH'] = r'C:\\mingw64\\bin;' + os.environ['PATH']

# Download common deps.
wget(r'https://github.com/bluescarni/binary_deps/raw/master/gmp_mingw_64.7z', 'gmp.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/mpfr_mingw_64.7z', 'mpfr.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/boost_mingw_64.7z', 'boost.7z')
# Extract them.
run_command(r'7z x -aoa -oC:\\ gmp.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ mpfr.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ boost.7z', verbose=False)

# Download piranha 0.8 https://github.com/bluescarni/piranha/archive/v0.8.zip
wget(r'https://github.com/bluescarni/piranha/archive/v0.8.zip', 'piranhav08.zip')
run_command(r'unzip piranhav08.zip', verbose=False)
# Move to the directory created and make piranha install its headers
os.chdir('piranha-0.8')
os.makedirs('build')
os.chdir('build')
print("Installing piranha")
run_command(r'cmake -G "MinGW Makefiles" .. -DCMAKE_INSTALL_PREFIX=c:\\local ', verbose=False)
run_command(r'mingw32-make install VERBOSE=1', verbose=False)
os.chdir('../../')
print("Piranha sucessfully installed .. continuing")

# Download audi 1.0.3 https://github.com/darioizzo/audi/archive/v1.0.3.zip
wget(r'https://github.com/darioizzo/audi/archive/v1.0.3.zip', 'audiv103.zip')
run_command(r'unzip audiv103.zip', verbose=False)
# Move to the directory created and make audi install its headers
os.chdir('audi-1.0.3')
os.makedirs('build')
os.chdir('build')
print("Installing audi")
run_command(r'cmake -G "MinGW Makefiles" .. -DBUILD_TESTS=no -DCMAKE_INSTALL_PREFIX=c:\\local ', verbose=False)
run_command(r'mingw32-make install VERBOSE=1', verbose=False)
os.chdir('../../')
print("Audi sucessfully installed .. continuing")

# Set the path so that the precompiled libs can be found.
os.environ['PATH'] = os.environ['PATH'] + r';c:\\local\\lib'

# Build type setup.
BUILD_TYPE = os.environ['BUILD_TYPE']
is_release_build = (os.environ['APPVEYOR_REPO_TAG'] == 'true') and bool(re.match(r'v[0-9]+\.[0-9]+.*',os.environ['APPVEYOR_REPO_TAG_NAME']))
if is_release_build:
    print("Release build detected, tag is '" + os.environ['APPVEYOR_REPO_TAG_NAME'] + "'")
is_python_build = 'Python' in BUILD_TYPE
if is_python_build:
    if BUILD_TYPE == 'Python35':
        python_version = '35'
    elif BUILD_TYPE == 'Python34':
        python_version = '34'
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
    run_command(pip + ' install pyaudi')
    #run_command(pip + ' install mpmath')
    if is_release_build:
        run_command(pip + ' install twine')

# Proceed to the build.
os.makedirs('build')
os.chdir('build')

common_cmake_opts = r'-DCMAKE_PREFIX_PATH=c:\\local -DCMAKE_INSTALL_PREFIX=c:\\local -DBUILD_MAIN=no' #REMEMBER TO PUT THIS TO YES WHEN FINISHED

# Configuration step.
if is_python_build:
    run_command(r'cmake -G "MinGW Makefiles" ..  -DDCGPY_INSTALL_PATH=c:\\local -DBUILD_TESTS=no -DBUILD_DCGPY=yes -DCMAKE_BUILD_TYPE=Release ' + common_cmake_opts + r' -DBoost_PYTHON_LIBRARY_RELEASE=c:\\local\\lib\\libboost_python' +
                (python_version[0] if python_version[0] == '3' else r'') + r'-mgw62-mt-1_62.dll  -DPYTHON_EXECUTABLE=C:\\Python' + python_version + r'\\python.exe -DPYTHON_LIBRARY=C:\\Python' + python_version + r'\\libs\\python' + python_version + r'.dll')
elif BUILD_TYPE in ['Release', 'Debug']:
    cmake_opts = r'-DCMAKE_BUILD_TYPE=' + BUILD_TYPE + r' -DBUILD_TESTS=yes ' + common_cmake_opts
    run_command(r'cmake -G "MinGW Makefiles" .. ' + cmake_opts)
else:
    raise RuntimeError('Unsupported build type: ' + BUILD_TYPE)

# Build+install step.
#run_command(r'cmake --build . --target install')
run_command(r'mingw32-make install VERBOSE=1')

# Testing, packaging.
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
    run_command(pip + r' install dist\\' + os.listdir('dist')[0])

    os.chdir('/')
    run_command(
        pinterp + r' -c "from dcgpy import test; test.run_test_suite()"')
    if is_release_build:
        os.chdir('C:/projects/d-cgp/build/wheel')
        run_command(twine + r' upload -u darioizzo dist\\' +
                    os.listdir('dist')[0])
elif BUILD_TYPE == 'Release':
    run_command(r'ctest -VV')
elif BUILD_TYPE == 'Debug':
    run_command(r'ctest -VV')
else:
    raise RuntimeError('Unsupported build type: ' + BUILD_TYPE)

import os
import re
import sys

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

# ----------------------------------SCRIPT START-----------------------------------------#
# Build type setup.
BUILD_TYPE = os.environ['BUILD_TYPE']
is_release_build = (os.environ['APPVEYOR_REPO_TAG'] == 'true') and bool(
    re.match(r'v[0-9]+\.[0-9]+.*', os.environ['APPVEYOR_REPO_TAG_NAME']))
if is_release_build:
    print("Release build detected, tag is '" +
          os.environ['APPVEYOR_REPO_TAG_NAME'] + "'")
is_python_build = 'Python' in BUILD_TYPE

# Check here for a list of installed software in the appveyor VMs: https://www.appveyor.com/docs/windows-images-software/
# USING: mingw64 8.1.0
ORIGINAL_PATH = os.environ['PATH']
run_command(r'mv C:\\mingw-w64\\x86_64-8.1.0-posix-seh-rt_v6-rev0\\mingw64 C:\\mingw64')
os.environ['PATH'] = r'C:\\mingw64\\bin;' + os.environ['PATH']
# Set the path so that the precompiled libs can be found.
os.environ['PATH'] = os.environ['PATH'] + r';c:\\local\\lib'

# Download common deps.
wget(r'https://github.com/bluescarni/binary_deps/raw/master/gmp_mingw81_64.7z', 'gmp.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/mpfr_mingw81_64.7z', 'mpfr.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/boost_mingw81-mt-x64-1_70.7z', 'boost.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/eigen3.7z', 'eigen3.7z')
wget(r'https://github.com/bluescarni/binary_deps/raw/master/tbb_2019_mgw81.7z', 'tbb.7z')

# Extract them.
run_command(r'7z x -aoa -oC:\\ gmp.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ mpfr.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ boost.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ eigen3.7z', verbose=False)
run_command(r'7z x -aoa -oC:\\ tbb.7z', verbose=False)

# Download piranha 0.11 https://github.com/bluescarni/piranha/archive/v0.11.zip
wget(r'https://github.com/bluescarni/piranha/archive/v0.11.zip', 'piranhav11.zip')
run_command(r'unzip piranhav11.zip', verbose=False)
# Move to the directory created and make piranha install its headers
os.chdir('piranha-0.11')
os.makedirs('build')
os.chdir('build')
print("Installing piranha")
run_command(
    r'cmake -G "MinGW Makefiles" .. -DCMAKE_INSTALL_PREFIX=c:\\local -DBoost_INCLUDE_DIR=c:\\local\\include', verbose=False)
run_command(r'mingw32-make install VERBOSE=1', verbose=False)
os.chdir('../../')
print("Piranha sucessfully installed .. continuing")

# Download audi 1.6 https://github.com/darioizzo/audi/archive/v1.6.zip
wget(r'https://github.com/darioizzo/audi/archive/v1.6.zip', 'audi.zip')
run_command(r'unzip audi.zip', verbose=False)
# Move to the directory created and make audi install its headers
os.chdir('audi-1.6')
os.makedirs('build')
os.chdir('build')
print("Installing audi")
run_command(r'cmake -G "MinGW Makefiles" .. ' + \
    r'-DAUDI_BUILD_AUDI=yes ' + \
    r'-DAUDI_BUILD_PYAUDI=no ' + \
    r'-DAUDI_BUILD_TEST=no ' + \
    r'-DAUDI_WITH_MPPP=no ' + \
    r'-DCMAKE_INSTALL_PREFIX=c:\\local ' + \
    r'-DBoost_INCLUDE_DIR=c:\\local\\include ' + \
    r'-DBoost_SERIALIZATION_LIBRARY_RELEASE=c:\\local\\lib\\libboost_serialization-mgw81-mt-x64-1_70.dll ' + \
    r'-DBoost_CHRONO_LIBRARY_RELEASE=c:\\local\\lib\\libboost_chrono-mgw81-mt-x64-1_70.dll ' + \
    r'-DBoost_SYSTEM_LIBRARY_RELEASE=c:\\local\\lib\\libboost_system-mgw81-mt-x64-1_70.dll ' + \
    r'-DBoost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE=c:\\local\\lib\\libboost_unit_test_framework-mgw81-mt-x64-1_70.dll ' + \
    r'-DBoost_TIMER_LIBRARY_RELEASE=c:\\local\\lib\\libboost_timer-mgw81-mt-x64-1_70.dll ', verbose=False)
run_command(r'mingw32-make install VERBOSE=1', verbose=False)
os.chdir('../../')
print("Audi sucessfully installed .. continuing")

# Setup of the Python build variables (python version based)
if is_python_build:
    if 'Python37-x64' in BUILD_TYPE:
        python_version = r'37'
        python_folder = r'Python37-x64'
        python_library = r'C:\\' + python_folder + r'\\python37.dll '
    elif 'Python36-x64' in BUILD_TYPE:
        python_version = '36'
        python_folder = r'Python36-x64'
        python_library = r'C:\\' + python_folder + r'\\python36.dll '
    elif 'Python27-x64' in BUILD_TYPE:
        python_version = r'27'
        python_folder = r'Python27-x64'
        python_library = r'C:\\' + python_folder + r'\\libs\\python27.dll '
        # Fot py27 I could not get it to work with the appveyor python (I was close but got tired).
        # Since this is anyway going to disappear (py27 really!!!), I am handling it as a one time workaround using the old py27 patched by bluescarni
        rm_fr(r'c:\\Python27-x64')
        wget(r'https://github.com/bluescarni/binary_deps/raw/master/python27_mingw_64.7z', 'python.7z')
        run_command(r'7z x -aoa -oC:\\ python.7z', verbose=False)
        run_command(r'mv C:\\Python27 C:\\Python27-x64', verbose=False)
    else:
        raise RuntimeError('Unsupported Python build: ' + BUILD_TYPE)

    # Set paths.
    pinterp = r"C:\\" + python_folder + r'\\python.exe'
    pip = r"C:\\" + python_folder + r'\\scripts\\pip'
    twine = r"C:\\" + python_folder + r'\\scripts\\twine'
    module_install_path = r"C:\\" + python_folder + r'\\Lib\\site-packages\\dcgpy'
    # Install pip and deps.
    run_command(pinterp + r' --version', verbose=True)
    wget(r'https://bootstrap.pypa.io/get-pip.py', 'get-pip.py')
    run_command(pinterp + ' get-pip.py --force-reinstall')
    run_command(pip + ' install numpy')
    run_command(pip + ' install pyaudi')
    if is_release_build:
        run_command(pip + ' install twine')

# Set the path so that the precompiled libs can be found.
os.environ['PATH'] = os.environ['PATH'] + r';c:\\local\\lib'

# Proceed to the build. The following arguments will be used for all build cases.
common_cmake_opts = r'-DCMAKE_PREFIX_PATH=c:\\local ' + \
        r'-DCMAKE_INSTALL_PREFIX=c:\\local ' + \
        r'-DBoost_INCLUDE_DIR=c:\\local\\include ' + \
        r'-DBoost_SERIALIZATION_LIBRARY_RELEASE=c:\\local\\lib\\libboost_serialization-mgw81-mt-x64-1_70.dll ' + \
        r'-DBoost_CHRONO_LIBRARY_RELEASE=c:\\local\\lib\\libboost_chrono-mgw81-mt-x64-1_70.dll ' + \
        r'-DBoost_SYSTEM_LIBRARY_RELEASE=c:\\local\\lib\\libboost_system-mgw81-mt-x64-1_70.dll ' + \
        r'-DBoost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE=c:\\local\\lib\\libboost_unit_test_framework-mgw81-mt-x64-1_70.dll ' + \
        r'-DBoost_TIMER_LIBRARY_RELEASE=c:\\local\\lib\\libboost_timer-mgw81-mt-x64-1_70.dll '

if is_python_build:
    os.makedirs('build_dcgp')
    os.chdir('build_dcgp')
    run_command(
        r'cmake -G "MinGW Makefiles" ..  -DCMAKE_BUILD_TYPE=Release -DDCGP_BUILD_TESTS=no -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_DCGPY=no' + ' ' + common_cmake_opts)
    run_command(r'mingw32-make install VERBOSE=1 -j2')
    os.chdir('..')
    os.makedirs('build_dcgpy')
    os.chdir('build_dcgpy')
    run_command(r'cmake -G "MinGW Makefiles" .. ' + \
            common_cmake_opts + \
            r'-DDCGPY_INSTALL_PATH=c:\\local ' + \
            r'-DDCGP_BUILD_DCGP=no ' + \
            r'-DDCGP_BUILD_DCGPY=yes ' + \
            r'-DCMAKE_BUILD_TYPE=Release ' + \
            r'-DBoost_PYTHON' + python_version + r'_LIBRARY_RELEASE=c:\\local\\lib\\libboost_python' + python_version + r'-mgw81-mt-x64-1_70.dll ' + \
            r'-DPYTHON_INCLUDE_DIR=C:\\' + python_folder + r'\\include ' + \
            r'-DPYTHON_EXECUTABLE=C:\\' + python_folder + r'\\python.exe ' + \
            r'-DPYTHON_LIBRARY=' + python_library)
    run_command(r'mingw32-make install VERBOSE=1 -j2')
elif BUILD_TYPE in ['Release', 'Debug']:
    os.makedirs('build_dcgp')
    os.chdir('build_dcgp')
    cmake_opts = r'-DCMAKE_BUILD_TYPE=' + BUILD_TYPE + \
        r' -DDCGP_BUILD_TESTS=yes ' \
        + common_cmake_opts
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
    shutil.move(module_install_path, r'.')
    wheel_libs = 'mingw_wheel_libs_python{}.txt'.format(python_version)
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

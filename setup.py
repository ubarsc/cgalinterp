from __future__ import print_function
import os
import sys
from numpy.distutils.core import setup, Extension

def getFlags():
    """
    Return the include and lib flags required
    """
    extra_includes = []
    extra_libs = []
    new_argv = [sys.argv[0]]
    #print("Flag ARGS: ", sys.argv[1:])
    for arg in sys.argv[1:]:
        handled = False
        if arg.startswith("--extraincludes="):
            incStr = arg.split(' ')[0].split('=')[1]
            libstr = arg.split(' ')[1].split('=')[1]
            extra = incStr.split(':')
            print('extra', extra)
            extra_includes.extend(extra)
            extra = libstr.split(':')
            print('extra', extra)
            extra_libs.extend(extra)
            handled = True
        if arg.startswith('--help'):
            print('Header options:')
            for opt in INCLUDE_OPTIONS:
                print(opt, 'Include path')
        if not handled:
            new_argv.append(arg)

    sys.argv = new_argv
    return extra_includes, extra_libs

# get the flags for GDAL etc
extraincludes, extralibs = getFlags()
print("Extra Includes", extraincludes)
print("Extra Libs", extralibs)

# create our extension
cgalinterpmodule = Extension(name="cgalinterp", 
                sources=["src/cgalinterp.cpp",],
                include_dirs=extraincludes,
                libraries=['CGAL', 'mpfr', 'gmp'],
                library_dirs=extralibs)

setup(name="cgalinterp",
        version="0.1",
        ext_modules=[cgalinterpmodule],
        description="Python Bindings for CGAL",
        author="Sam Gillingham",
        author_email="gillingham.sam@gmail.com")

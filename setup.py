from __future__ import print_function
import os
import sys
from numpy.distutils.core import setup, Extension

def getFlags():
    """
    Return the include flags required
    """
    extra_includes = []
    new_argv = [sys.argv[0]]
    for arg in sys.argv[1:]:
        handled = False
        if arg.startswith("--extraincludes="):
            inc = arg.split('=')[1]
            extra = inc.split(' ')
            print('extra', extra)
            extra_includes.extend(extra)
            handled = True
        if arg.startswith('--help'):
            print('Header options:')
            for opt in INCLUDE_OPTIONS:
                print(opt, 'Include path')
        if not handled:
            new_argv.append(arg)

    sys.argv = new_argv
    return extra_includes

# get the flags for GDAL etc
extraincludes = getFlags()
print(extraincludes)

# create our extension
cgalinterpmodule = Extension(name="cgalinterp", 
                sources=["src/cgalinterp.cpp",],
                include_dirs=extraincludes)

setup(name="cgalinterp",
        version="0.1",
        ext_modules=[cgalinterpmodule],
        description="Python Bindings for CGAL",
        author="Sam Gillingham",
        author_email="gillingham.sam@gmail.com")

from __future__ import print_function
import os
from numpy.distutils.core import setup, Extension

INC_VAR = 'CGINT_EXTRA_INCLUDES'
LIB_VAR = 'CGINT_EXTRA_LIBS'

def getFlags():
    """
    Return the include and lib flags required
    """
    extra_includes = None
    exta_libs = None
    if INC_VAR in os.environ:
        extra_includes = os.environ[INC_VAR].split(os.pathsep)
    if INC_VAR in os.environ:
        extra_includes = os.environ[INC_VAR].split(os.pathsep)
    return extra_includes, exta_libs    
    

# get the flags for GDAL etc
extraincludes, extralibs = getFlags()
print("Extra Includes", extraincludes)
print("Extra Libs", extralibs)

# create our extension
cgalinterpmodule = Extension(name="cgalinterp", 
                sources=["src/cgalinterp.cpp",],
                include_dirs=extraincludes,
                libraries=['CGAL', 'mpfr', 'gmp'],
                library_dirs=extralibs,
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')])

setup(name="cgalinterp",
        version="0.2",
        ext_modules=[cgalinterpmodule],
        description="Python Bindings for CGAL",
        author="Sam Gillingham",
        author_email="gillingham.sam@gmail.com")

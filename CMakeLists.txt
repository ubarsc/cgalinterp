###############################################################################
#
# CMake build scripts for CGALINTERP
# 
# Created 26/03/2015 by Peter Bunting
#
# These scripts were initial based on those used for libLAS (http://liblas.org/)
# and then edited for SPDLib (http://www.spdlib.org) before being edited again
# for CGALINTERP and then again for cgalinterp.
#
# History
# 26/03/2015 - Created by Peter Bunting for cgalinterp
#
###############################################################################

###############################################################################
# Set Project name and version
project (CGALINTERP)

set (PROJECT_LIBRARY_DIR lib)
set (PROJECT_SOURCE_DIR src)

# The version number.
set (CGALINTERP_VERSION_MAJOR 0)
set (CGALINTERP_VERSION_MINOR 1)
set (CGALINTERP_VERSION_PATCH 0)
set (CGALINTERP_VERSION "${CGALINTERP_VERSION_MAJOR}.${CGALINTERP_VERSION_MINOR}.${CGALINTERP_VERSION_PATCH}")
set (CGALINTERP_PACKAGE_VERSION "${CGALINTERP_VERSION_MAJOR}.${CGALINTERP_VERSION_MINOR}.${CGALINTERP_VERSION_PATCH}")
set (CGALINTERP_PACKAGE_STRING "CGALINTERP ${CGALINTERP_VERSION_MAJOR}.${CGALINTERP_VERSION_MINOR}.${CGALINTERP_VERSION_PATCH}")
set (CGALINTERP_PACKAGE_BUGREPORT "rsgislib-support@googlegroups.com")
set (CGALINTERP_PACKAGE "cgalInterp")
set (CGALINTERP_COPYRIGHT_YEAR 2015)

# Name of C++ library
set(CGALINTERP_LIB_NAME cgalinterp)

set(PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE PATH "Install path for Python libraries.")

set(BOOST_INCLUDE_DIR /usr/local/include CACHE PATH "Include PATH for Boost")
#set(BOOST_LIB_PATH /usr/local/lib CACHE PATH "Library PATH for Boost")

set(CGAL_INCLUDE_DIR /usr/local/include CACHE PATH "Include PATH for CGAL")
#set(CGAL_LIB_PATH /usr/local/lib CACHE PATH "Library PATH for CGAL")

set(GMP_INCLUDE_DIR /usr/local/include CACHE PATH "Include PATH for GMP")
#set(GMP_LIB_PATH /usr/local/lib CACHE PATH "Library PATH for GMP")

set(MPFR_INCLUDE_DIR /usr/local/include CACHE PATH "Include PATH for MPFR")
#set(MPFR_LIB_PATH /usr/local/lib CACHE PATH "Library PATH for MPFR")
###############################################################################

###############################################################################
# CMake settings
cmake_minimum_required(VERSION 2.6.0)

IF(NOT CMAKE_BUILD_TYPE)
  #SET(CMAKE_BUILD_TYPE "DEBUG")
  SET(CMAKE_BUILD_TYPE "RELEASE")
  #SET(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
  #SET(CMAKE_BUILD_TYPE "MINSIZEREL")
ENDIF()

set(CMAKE_COLOR_MAKEFILE ON)

# Allow advanced users to generate Makefiles printing detailed commands
mark_as_advanced(CMAKE_VERBOSE_MAKEFILE)

# Path to additional CMake modules
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ./cmake/modules/)
###############################################################################

###############################################################################
# Setup configure file
configure_file ( "${PROJECT_SOURCE_DIR}/cgalinterp-config.h.in" "${PROJECT_SOURCE_DIR}/cgalinterp-config.h" )
###############################################################################

###############################################################################
# Platform and compiler specific settings

###############################################################################

###############################################################################
# Check the required libraries are present


include_directories(${BOOST_INCLUDE_DIR})

FIND_PROGRAM(PYTHON_EXE python ENV PATH)
IF(PYTHON_EXE STREQUAL "PYTHON_EXE-NOTFOUND")
  MESSAGE(STATUS "Python not found. Python bindings disabled")
  SET(BUILD_PYTHON_BINDINGS OFF)
ELSEIF(NOT PYTHON_EXE STREQUAL "PYTHON_EXE-NOTFOUND")
  MESSAGE(STATUS "Building Python bindings with ${PYTHON_EXE}")
  SET(BUILD_PYTHON_BINDINGS ON)
ENDIF()

include_directories(${CGAL_INCLUDE_DIR})

include_directories(${GMP_INCLUDE_DIR})

include_directories(${MPFR_INCLUDE_DIR})
###############################################################################

###############################################################################
# Build library

set(PYTHON_ARGS "--extraincludes=${CGAL_INCLUDE_DIR} ${BOOST_INCLUDE_DIR} ${MPFR_INCLUDE_DIR}")
MESSAGE(STATUS "${PYTHON_ARGS}")

# Python 
if (BUILD_PYTHON_BINDINGS)
  add_custom_target("Python_build" ALL COMMAND ${PYTHON_EXE} "setup.py" "build" "${PYTHON_ARGS}" WORKING_DIRECTORY ".")
endif(BUILD_PYTHON_BINDINGS)

###############################################################################

###############################################################################
# Installation 

# Python 
if (BUILD_PYTHON_BINDINGS)
  install(CODE "execute_process(COMMAND ${PYTHON_EXE} \"setup.py\" \"install\" \"--prefix=${PYTHON_INSTALL_PREFIX}\" \"${PYTHON_ARGS}\"  WORKING_DIRECTORY \".\")")
endif(BUILD_PYTHON_BINDINGS)

###############################################################################

###############################################################################
# Create Distribution
SET(CPACK_SOURCE_GENERATOR "TGZ;ZIP;TBZ2")
SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")

set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CGALINTERP_VERSION_MAJOR}.${CGALINTERP_VERSION_MINOR}.${CGALINTERP_VERSION_PATCH}")

# Set files to ignore
list(APPEND CPACK_SOURCE_IGNORE_FILES "_CPack_Packages")
list(APPEND CPACK_SOURCE_IGNORE_FILES ".gz")
list(APPEND CPACK_SOURCE_IGNORE_FILES ".bz2")
list(APPEND CPACK_SOURCE_IGNORE_FILES ".zip")
list(APPEND CPACK_SOURCE_IGNORE_FILES ".svn")
list(APPEND CPACK_SOURCE_IGNORE_FILES "README")
list(APPEND CPACK_SOURCE_IGNORE_FILES "HOWTORELEASE.txt")
list(APPEND CPACK_SOURCE_IGNORE_FILES "CMakeCache.txt")
list(APPEND CPACK_SOURCE_IGNORE_FILES "CPackConfig.cmake")
list(APPEND CPACK_SOURCE_IGNORE_FILES "schemas")

include(CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
###############################################################################

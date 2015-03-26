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
set(BOOST_LIB_PATH /usr/local/lib CACHE PATH "Library PATH for Boost")

set(CGAL_INCLUDE_DIR /usr/local/include CACHE PATH "Include PATH for CGAL")
set(CGAL_LIB_PATH /usr/local/lib CACHE PATH "Library PATH for CGAL")

set(GMP_INCLUDE_DIR /usr/local/include CACHE PATH "Include PATH for GMP")
set(GMP_LIB_PATH /usr/local/lib CACHE PATH "Library PATH for GMP")

set(MPFR_INCLUDE_DIR /usr/local/include CACHE PATH "Include PATH for MPFR")
set(MPFR_LIB_PATH /usr/local/lib CACHE PATH "Library PATH for MPFR")
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
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ./cmake/modules/)
###############################################################################

###############################################################################
# Setup configure file
configure_file ( "${PROJECT_SOURCE_DIR}/cgalinterp-config.h.in" "${PROJECT_SOURCE_DIR}/cgalinterp-config.h" )
###############################################################################

###############################################################################
# Platform and compiler specific settings

if(WIN32)
    if (MSVC)
        if (MSVC80 OR MSVC90 OR MSVC10)
            add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
            add_definitions(-D_CRT_SECURE_NO_WARNINGS)
            add_definitions(-D_CRT_NONSTDC_NO_WARNING)
            add_definitions(-D_SCL_SECURE_NO_WARNINGS)
        endif()
        
        if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
            string(REGEX REPLACE "/W[0-4]" "/W4"
                   CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        else()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
        endif()

        # Generate dot-user file with user-specific settings for Visual Studio project
        set(MSVC_ENVIRONMENT_PATH "" CACHE STRING
            "Custom PATH for Environment property in Visual Studio project configuration")
        mark_as_advanced(MSVC_ENVIRONMENT_PATH)
        set(VCPROJ_USER_ENVIRONMENT_PATH_DEBUG "${ENVIRONMENT_PATH}")
        mark_as_advanced(VCPROJ_USER_ENVIRONMENT_PATH_DEBUG)

    endif()

	set(CMAKE_INCLUDE_PATH c:/osgeo4w/include;$ENV{CMAKE_INCLUDE_PATH})
    set(CMAKE_LIBRARY_PATH c:/osgeo4w/lib;$ENV{CMAKE_LIBRARY_PATH})
else()
  # Recommended C++ compilation flags
  # -Weffc++
  set(CGALINTERP_COMMON_CXX_FLAGS
	"-Wall -Wpointer-arith -Wcast-align -Wcast-qual -Wredundant-decls -Wno-long-long")
	#"-pedantic -ansi -Wall -Wpointer-arith -Wcast-align -Wcast-qual -Wfloat-equal -Wredundant-decls -Wno-long-long")

  if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC ${CGALINTERP_COMMON_CXX_FLAGS}")
    if (CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")
    endif()

  elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER}" MATCHES "clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGALINTERP_COMMON_CXX_FLAGS}")
  endif()
endif(WIN32)

if (APPLE)
  set(SO_EXT dylib)
 # set(CMAKE_FIND_FRAMEWORK "LAST")
elseif(WIN32)
  set(SO_EXT dll)
else()
  set(SO_EXT so)
endif(APPLE)
###############################################################################

###############################################################################
# Check the required libraries are present


include_directories(${BOOST_INCLUDE_DIR})
if (MSVC)
    set(BOOST_LIBRARIES -LIBPATH:${BOOST_LIB_PATH} libboost_filesystem-vc100-mt-1_56.lib libboost_system-vc100-mt-1_56.lib)
else()
    set(BOOST_LIBRARIES -L${BOOST_LIB_PATH} -lboost_filesystem -lboost_system)
endif(MSVC)

FIND_PROGRAM(PYTHON_EXE python ENV PATH)
IF(PYTHON_EXE STREQUAL "PYTHON_EXE-NOTFOUND")
  MESSAGE(STATUS "Python not found. Python bindings disabled")
  SET(BUILD_PYTHON_BINDINGS OFF)
ELSEIF(NOT PYTHON_EXE STREQUAL "PYTHON_EXE-NOTFOUND")
  MESSAGE(STATUS "Building Python bindings with ${PYTHON_EXE}")
  SET(BUILD_PYTHON_BINDINGS ON)
ENDIF()

include_directories(${CGAL_INCLUDE_DIR})
if (MSVC)
    set(CGAL_LIBRARIES -LIBPATH:${CGAL_LIB_PATH} CGAL-vc100-mt-4.3.lib)
else()
    set(CGAL_LIBRARIES -L${CGAL_LIB_PATH} -lCGAL)
endif(MSVC)

include_directories(${GMP_INCLUDE_DIR})
if (MSVC)
    set(GMP_LIBRARIES -LIBPATH:${GMP_LIB_PATH} libgmp-10.lib)
else()
    set(GMP_LIBRARIES -L${GMP_LIB_PATH} -lgmp)
endif(MSVC)

include_directories(${MPFR_INCLUDE_DIR})
if (MSVC)
    set(MPFR_LIBRARIES -LIBPATH:${MPFR_LIB_PATH} libmpfr-4.lib)
else()
    set(MPFR_LIBRARIES -L${MPFR_LIB_PATH} -lmpfr)
endif(MSVC)
###############################################################################

###############################################################################
# Build library

set(LIB_CGALINTERP_CPP
	${PROJECT_SOURCE_DIR}/cgalinterp.cpp
)

add_library( ${CGALINTERP_LIB_NAME} ${LIB_CGALINTERP_CPP}  )
target_link_libraries(${CGALINTERP_LIB_NAME} ${BOOST_LIBRARIES} ${CGAL_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} )

SET_TARGET_PROPERTIES(${CGALINTERP_LIB_NAME}
PROPERTIES
        SOVERSION ${CGALINTERP_VERSION}
        VERSION ${CGALINTERP_VERSION}
)

###############################################################################

###############################################################################
# Installation 
install (TARGETS ${CGALINTERP_LIB_NAME} DESTINATION lib)
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

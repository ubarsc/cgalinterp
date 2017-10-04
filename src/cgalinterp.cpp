/*
 *  cgalinterp.cpp
 *
 *
 * This file is part of PyLidar
 * Copyright (C) 2015 John Armston, Pete Bunting, Neil Flood, Sam Gillingham
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// following needed or CGAL templates misbehave
#ifndef NDEBUG
    #define NDEBUG
#endif

#include <math.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/algorithm.h>
#include <CGAL/Origin.h>
#include <CGAL/squared_distance_2.h>

#include <Python.h>
#include "numpy/arrayobject.h"
    
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                         CGALCoordType;
typedef K::Vector_2                                   CGALVector;
typedef K::Point_2                                    CGALPoint;
    
typedef CGAL::Delaunay_triangulation_2<K>             DelaunayTriangulation;
typedef CGAL::Interpolation_traits_2<K>               InterpTraits;
typedef CGAL::Delaunay_triangulation_2<K>::Vertex_handle    Vertex_handle;
typedef CGAL::Delaunay_triangulation_2<K>::Face_handle    Face_handle;
    
typedef std::vector< std::pair<CGALPoint, CGALCoordType> >   CoordinateVector;
typedef std::map<CGALPoint, CGALCoordType, K::Less_xy_2>     PointValueMap;

/* An exception object for this module */
/* created in the init function */
struct CGALInterpState
{
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct CGALInterpState*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct CGALInterpState _state;
#endif

static PyObject *cgalinterp_naturalneighbour(PyObject *self, PyObject *args)
{
    //std::cout.precision(12);
    PyArrayObject *pXVals, *pYVals, *pZVals, *pXGrid, *pYGrid;
    PyArrayObject *pOutArray;

    if( !PyArg_ParseTuple(args, "OOOOO:NaturalNeighbour", &pXVals, &pYVals, &pZVals, &pXGrid, &pYGrid))
        return NULL;

    if( !PyArray_Check(pXVals) || !PyArray_Check(pYVals) || !PyArray_Check(pZVals) || !PyArray_Check(pXGrid) || !PyArray_Check(pYGrid) )
    {
        PyErr_SetString(GETSTATE(self)->error, "All arguments must be numpy arrays");
        return NULL;
    }

    // Check dimensions match
    if( (PyArray_DIM(pXVals, 0) != PyArray_DIM(pYVals, 0)) | (PyArray_DIM(pXVals, 0) != PyArray_DIM(pZVals, 0)))
    {
        PyErr_SetString(GETSTATE(self)->error, "Training X, Y and Z arrays must all be of the same length");
        return NULL;
    }
    
    if( (PyArray_DIM(pXGrid, 0) != PyArray_DIM(pYGrid, 0)) | (PyArray_DIM(pXGrid, 1) != PyArray_DIM(pYGrid, 1)))
    {
        PyErr_SetString(GETSTATE(self)->error, "X and Y grids must have the same dimensions");
        return NULL;
    }
    
    // check types ok
    if( (PyArray_TYPE(pXVals) != NPY_DOUBLE) || (PyArray_TYPE(pYVals) != NPY_DOUBLE) || 
        (PyArray_TYPE(pZVals) != NPY_DOUBLE) || (PyArray_TYPE(pXGrid) != NPY_DOUBLE) ||
        (PyArray_TYPE(pYGrid) != NPY_DOUBLE) )
    {
        PyErr_SetString(GETSTATE(self)->error, "All input arrays must be double");
        return NULL;
    }

    npy_intp nRows = PyArray_DIM(pXGrid, 0);
    npy_intp nCols = PyArray_DIM(pXGrid, 1);
    
    npy_intp nVals = PyArray_DIM(pXVals, 0);

    // Create output
    pOutArray = (PyArrayObject*)PyArray_EMPTY(2, PyArray_DIMS(pXGrid), NPY_DOUBLE, 0);
    if( pOutArray == NULL )
    {
        PyErr_SetString(GETSTATE(self)->error, "Failed to create array");
        return NULL;
    }

    if( PyArray_DIM(pXVals, 0) < 3 )
    {
        PyErr_SetString(GETSTATE(self)->error, "Not enough points, need at least 3.");
        return NULL;
    }    

    if( PyArray_DIM(pXVals, 0) < 100 )
    {
        // check that these small number of points aren't all within a line
        double meanX = 0;
        double meanY = 0;
                
        double varX = 0;
        double varY = 0;
                
        for(npy_intp i = 0; i < nVals; ++i)
        {
            meanX += *((double*)PyArray_GETPTR1(pXVals, i));
            meanY += *((double*)PyArray_GETPTR1(pYVals, i));
        }
                
        meanX = meanX / nVals;
        meanY = meanY / nVals;
                
        //std::cout << "meanX = " << meanX << std::endl;
        //std::cout << "meanY = " << meanY << std::endl;
                
        for(npy_intp i = 0; i < nVals; ++i)
        {
            varX += *((double*)PyArray_GETPTR1(pXVals, i)) - meanX;
            varY += *((double*)PyArray_GETPTR1(pYVals, i)) - meanY;
        }
                
        varX = fabs(varX / nVals);
        varY = fabs(varY / nVals);
                
        //std::cout << "varX = " << varX << std::endl;
        //std::cout << "varY = " << varX << std::endl;
                
        if((varX < 4) || (varY < 4))
        {
            PyErr_SetString(GETSTATE(self)->error, "Points are all within a line.");
            return NULL;
        }
    }
    try
    {
        //std::cout << "Perform Interpolation\n";
        DelaunayTriangulation dt;
        PointValueMap values;
        
        //std::cout << "Building Triangulation\n";
        for(npy_intp i = 0; i < nVals; ++i)
        {
            K::Point_2 cgalPt( *((double*)PyArray_GETPTR1(pXVals, i)), *((double*)PyArray_GETPTR1(pYVals, i)) );
            dt.insert(cgalPt);
            CGALCoordType value = *((double*)PyArray_GETPTR1(pZVals, i));
            //std::cout << i << "\t[" << *((double*)PyArray_GETPTR1(pXVals, i)) << ", " << *((double*)PyArray_GETPTR1(pYVals, i)) << "] = " << *((double*)PyArray_GETPTR1(pZVals, i)) << std::endl;
            values.insert(std::make_pair(cgalPt, value));
        }


        for(npy_intp i  = 0; i < nRows; ++i)
        {
            for(npy_intp j = 0; j < nCols; ++j)
            {
                K::Point_2 p( *((double*)PyArray_GETPTR2(pXGrid, i, j)), *((double*)PyArray_GETPTR2(pYGrid, i, j)) );
                CoordinateVector coords;
                CGAL::Triple<std::back_insert_iterator<CoordinateVector>, K::FT, bool> result = CGAL::natural_neighbor_coordinates_2(dt, p, std::back_inserter(coords));
                if(!result.third)
                {
                    //std::cout << "\t Not within convex hull\n";
                    // Not within convex hull of dataset
                    *((double*)PyArray_GETPTR2(pOutArray, i, j)) = 0.0;
                }
                else
                {
                    CGALCoordType norm = result.second;
                    CGALCoordType outValue = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, CGAL::Data_access<PointValueMap>(values));
                    //std::cout << "pt: [" << *((double*)PyArray_GETPTR2(pXGrid, i, j)) << ", " << *((double*)PyArray_GETPTR2(pYGrid, i, j)) << "]\n";
                    //std::cout << "\t NN Value = " << outValue << std::endl;
                    *((double*)PyArray_GETPTR2(pOutArray, i, j)) = outValue;
                }
            }
        }
    }
    catch(std::exception &e)
    {
        PyErr_SetString(GETSTATE(self)->error, e.what());
        return NULL;
    }
    return (PyObject*)pOutArray;
}    

static PyObject *cgalinterp_naturalneighbour_pts(PyObject *self, PyObject *args)
{
    //std::cout.precision(12);
    PyArrayObject *pXVals, *pYVals, *pZVals, *pXYGrid;
    PyArrayObject *pOutArray;
    npy_intp nPts, nVals, i;
    
    if( !PyArg_ParseTuple(args, "OOOO:NaturalNeighbourPts", &pXVals, &pYVals, &pZVals, &pXYGrid))
        return NULL;
    
    if( !PyArray_Check(pXVals) || !PyArray_Check(pYVals) || !PyArray_Check(pZVals) || !PyArray_Check(pXYGrid) )
    {
        PyErr_SetString(GETSTATE(self)->error, "All arguments must be numpy arrays");
        return NULL;
    }

    // check dims
    if( (PyArray_NDIM(pXVals) != 1) || (PyArray_NDIM(pYVals) != 1) || (PyArray_NDIM(pZVals) != 1) || 
            (PyArray_NDIM(pXYGrid) != 2) )
    {
        PyErr_SetString(GETSTATE(self)->error, "Arrays should be 1d, 1d, 1d and 2d respectively");
        return NULL;
    }
    
    // Check dimensions match
    if( (PyArray_DIM(pXVals, 0) != PyArray_DIM(pYVals, 0)) | (PyArray_DIM(pXVals, 0) != PyArray_DIM(pZVals, 0)))
    {
        PyErr_SetString(GETSTATE(self)->error, "Training X, Y and Z arrays must all be of the same length");
        return NULL;
    }
    
    if( PyArray_DIM(pXYGrid, 1) != 2 )
    {
        PyErr_SetString(GETSTATE(self)->error, "Interpolation point array must be shape N*2");
        return NULL;
    }
    
    // check types ok
    if( (PyArray_TYPE(pXVals) != NPY_DOUBLE) || (PyArray_TYPE(pYVals) != NPY_DOUBLE) || 
        (PyArray_TYPE(pZVals) != NPY_DOUBLE) || (PyArray_TYPE(pXYGrid) != NPY_DOUBLE) )
    {
        PyErr_SetString(GETSTATE(self)->error, "All input arrays must be double");
        return NULL;
    }
    
    nPts = PyArray_DIM(pXYGrid, 0);
    nVals = PyArray_DIM(pXVals, 0);
    
    // Create output
    pOutArray = (PyArrayObject*)PyArray_EMPTY(1, &nPts, NPY_DOUBLE, 0);
    if( pOutArray == NULL )
    {
        PyErr_SetString(GETSTATE(self)->error, "Failed to create array");
        return NULL;
    }
    
    if( PyArray_DIM(pXVals, 0) < 3 )
    {
        PyErr_SetString(GETSTATE(self)->error, "Not enough points, need at least 3.");
        Py_DECREF(pOutArray);
        return NULL;
    }
    
    try
    {
        //std::cout << "Perform Interpolation\n";
        DelaunayTriangulation dt;
        PointValueMap values;
        
        //std::cout << "Building Triangulation\n";
        for(i = 0; i < nVals; ++i)
        {
            K::Point_2 cgalPt( *((double*)PyArray_GETPTR1(pXVals, i)), *((double*)PyArray_GETPTR1(pYVals, i)) );
            dt.insert(cgalPt);
            CGALCoordType value = *((double*)PyArray_GETPTR1(pZVals, i));
            //std::cout << i << "\t[" << *((double*)PyArray_GETPTR1(pXVals, i)) << ", " << *((double*)PyArray_GETPTR1(pYVals, i)) << "] = " << *((double*)PyArray_GETPTR1(pZVals, i)) << std::endl;
            values.insert(std::make_pair(cgalPt, value));
        }
        for(i = 0; i < nPts; ++i)
        {
            K::Point_2 p( *((double*)PyArray_GETPTR2(pXYGrid, i, 0)), *((double*)PyArray_GETPTR2(pXYGrid, i, 1)) );
            CoordinateVector coords;
            CGAL::Triple<std::back_insert_iterator<CoordinateVector>, K::FT, bool> result = CGAL::natural_neighbor_coordinates_2(dt, p, std::back_inserter(coords));
            if(!result.third)
            {
                //std::cout << "\t Not within convex hull\n";
                // Not within convex hull of dataset
                *((double*)PyArray_GETPTR1(pOutArray, i)) = 0.0;
            }
            else
            {
                CGALCoordType norm = result.second;
                CGALCoordType outValue = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, CGAL::Data_access<PointValueMap>(values));
                //std::cout << "pt: [" << *((double*)PyArray_GETPTR2(pXGrid, i, j)) << ", " << *((double*)PyArray_GETPTR2(pYGrid, i, j)) << "]\n";
                //std::cout << "\t NN Value = " << outValue << std::endl;
                *((double*)PyArray_GETPTR1(pOutArray, i)) = outValue;
            }
        }
    }
    catch(std::exception &e)
    {
        PyErr_SetString(GETSTATE(self)->error, e.what());
        return NULL;
    }

    return (PyObject*)pOutArray;
}


/* Our list of functions in this module*/
static PyMethodDef CGALInterpMethods[] = {
    {"NaturalNeighbour", cgalinterp_naturalneighbour, METH_VARARGS,
"Perform Natural Neighbour Interpolation\n"
"call signature: arr = NaturalNeighbour(xvals, yvals, zvals, xgrid, ygrid)\n"
"where:\n"
"  xvals is a 1d array of the x values of the points\n"
"  yvals is a 1d array of the y values of the points\n"
"  zvals is a 1d array of the z values of the points\n"
"xvals, yvals and zvals should have the same length\n"
"  xgrid is a 2d array of x coordinates to interpolate at\n"
"  ygrid is a 2d array of y coordinates to interpolate at\n"
"xgrid and xgrid must be the same shape"}, 
    {"NaturalNeighbourPts", cgalinterp_naturalneighbour_pts, METH_VARARGS,
        "Perform Natural Neighbour Interpolation at points\n"
        "call signature: arr = NaturalNeighbourPts(xvals, yvals, zvals, xygrid)\n"
        "where:\n"
        "  xvals is a 1d array of the x values of the points\n"
        "  yvals is a 1d array of the y values of the points\n"
        "  zvals is a 1d array of the z values of the points\n"
        "xvals, yvals and zvals should have the same length\n"
        "  xygrid is a 2d array of coordinates in the form x,y\n"},
    {NULL}        /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3

static int cgalinterp_traverse(PyObject *m, visitproc visit, void *arg) 
{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int cgalinterp_clear(PyObject *m) 
{
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cgalinterp",
        NULL,
        sizeof(struct CGALInterpState),
        CGALInterpMethods,
        NULL,
        cgalinterp_traverse,
        cgalinterp_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC 
PyInit_cgalinterp(void)

#else
#define INITERROR return

PyMODINIT_FUNC
initcgalinterp(void)
#endif
{
    PyObject *pModule;
    struct CGALInterpState *state;

    /* initialize the numpy stuff */
    import_array();

#if PY_MAJOR_VERSION >= 3
    pModule = PyModule_Create(&moduledef);
#else
    pModule = Py_InitModule("cgalinterp", CGALInterpMethods);
#endif
    if( pModule == NULL )
        INITERROR;

    state = GETSTATE(pModule);

    /* Create and add our exception type */
    state->error = PyErr_NewException("cgal.error", NULL, NULL);
    if( state->error == NULL )
    {
        Py_DECREF(pModule);
        INITERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return pModule;
#endif
}



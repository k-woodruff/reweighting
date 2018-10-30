#include <Python.h>
#include "zeformfactor.h"

static char module_docstring[] =
    "This module provides an interface for calculating NCE xsec using C.";
    static char sigplus_docstring[] =
        "Calculate the NC elastic cross section of some model.";

static PyObject *zeformfactor_sigplus(PyObject *self, PyObject *args);
static PyObject *zeformfactor_sigplusCC(PyObject *self, PyObject *args);
static PyObject *zeformfactor_GEp(PyObject *self, PyObject *args);
static PyObject *zeformfactor_GMp(PyObject *self, PyObject *args);
static PyObject *zeformfactor_GEn(PyObject *self, PyObject *args);
static PyObject *zeformfactor_GMn(PyObject *self, PyObject *args);
static PyObject *zeformfactor_GACC(PyObject *self, PyObject *args);
static PyObject *zeformfactor_GAS(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"sigplus", zeformfactor_sigplus, METH_VARARGS, sigplus_docstring},
    {"sigplusCC", zeformfactor_sigplusCC, METH_VARARGS, sigplus_docstring},
    {"GEp", zeformfactor_GEp, METH_VARARGS,sigplus_docstring},
    {"GMp", zeformfactor_GMp, METH_VARARGS,sigplus_docstring},
    {"GEn", zeformfactor_GEn, METH_VARARGS,sigplus_docstring},
    {"GMn", zeformfactor_GMn, METH_VARARGS,sigplus_docstring},
    {"GACC", zeformfactor_GACC, METH_VARARGS,sigplus_docstring},
    {"GAS", zeformfactor_GAS, METH_VARARGS,sigplus_docstring},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_zeformfactor(void)
{
    PyObject *m = Py_InitModule3("_zeformfactor", module_methods, module_docstring);
    if (m == NULL)
        return;
}

static PyObject *zeformfactor_sigplus(PyObject *self, PyObject *args)
{
    double Ev, q2, delta_s, MA_s, a2_s;
    int nusign;
    int nucsign;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dddddii", &Ev, &q2, &delta_s, &MA_s, &a2_s,
                                                           &nusign, &nucsign))
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = sigplus(Ev,q2,delta_s,MA_s,a2_s,nusign,nucsign);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}

static PyObject *zeformfactor_sigplusCC(PyObject *self, PyObject *args)
{
    double Ev, q2;
    int nusign;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "ddi", &Ev, &q2, &nusign)) 
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = sigplus_CC(Ev,q2,nusign);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}

static PyObject *zeformfactor_GEp(PyObject *self, PyObject *args)
{
  double q2;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "d", &q2 )) 
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = GEp(q2);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}
static PyObject *zeformfactor_GMp(PyObject *self, PyObject *args)
{
  double q2;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "d", &q2 )) 
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = GMp(q2);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}
static PyObject *zeformfactor_GEn(PyObject *self, PyObject *args)
{
  double q2;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "d", &q2 )) 
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = GEn(q2);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}
static PyObject *zeformfactor_GMn(PyObject *self, PyObject *args)
{
  double q2;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "d", &q2 )) 
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = GMn(q2);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}
static PyObject *zeformfactor_GACC(PyObject *self, PyObject *args)
{
  double q2;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "d", &q2 )) 
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = GACC(q2);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}
static PyObject *zeformfactor_GAS(PyObject *self, PyObject *args)
{
  double q2,ds,mas,a2s;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dddd", &q2,&ds,&mas,&a2s )) 
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = GAS(q2,ds,mas,a2s);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}


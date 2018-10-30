#include <Python.h>
#include "ncformfactor.h"

static char module_docstring[] =
    "This module provides an interface for calculating NCE xsec using C.";
    static char sigplus_docstring[] =
        "Calculate the NC elastic cross section of some model.";

static PyObject *ncformfactor_sigplus(PyObject *self, PyObject *args);
static PyObject *ncformfactor_sigplusCC(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"sigplus", ncformfactor_sigplus, METH_VARARGS, sigplus_docstring},
    {"sigplusCC", ncformfactor_sigplusCC, METH_VARARGS, sigplus_docstring},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_ncformfactor(void)
{
    PyObject *m = Py_InitModule3("_ncformfactor", module_methods, module_docstring);
    if (m == NULL)
        return;
}

static PyObject *ncformfactor_sigplus(PyObject *self, PyObject *args)
{
    double Ev, q2, delta_s, MA, lambda_a, s_a, mu_s, rho_s;
    char gas;
    int nusign;
    int nucsign;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "ddcddddddii", &Ev, &q2, &gas, &delta_s, 
                                    &MA, &lambda_a, &s_a, &mu_s, &rho_s, &nusign, &nucsign))
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = sigplus(Ev,q2,gas,delta_s,MA,lambda_a,s_a,mu_s,rho_s,nusign,nucsign);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}

static PyObject *ncformfactor_sigplusCC(PyObject *self, PyObject *args)
{
    double Ev, q2, MA;
    int nusign;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dddi", &Ev, &q2, &MA, &nusign)) 
        return NULL;

    /* Call the external C function to compute the cross-section. */
    double value = sigplus_CC(Ev,q2,MA,nusign);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", value);
    return ret;
}


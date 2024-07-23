#include "Python.h"
#include <stdio.h>

#ifndef SOURCE_ROOT
#define SOURCE_ROOT ../../../..
#endif

#define XSTRING(s) STRING(s)
#define STRING(s) #s

PyObject* get_calforce_module()
{
  /* https://stackoverflow.com/questions/42393369/c-importing-python-module */
  PyObject* pydis = PyImport_ImportModule("pydis");
  PyObject* calforce_module = PyObject_GetAttrString(pydis, "CalForce");
  return calforce_module;
}

int main(){
  printf("Source Root = %s\n", XSTRING(SOURCE_ROOT));
  Py_Initialize();

  PyRun_SimpleString("import sys, os");
  PyRun_SimpleString("print('Python version:', sys.version)");
  char command[1000];
  sprintf(command, "pydis_paths = ['%s/python', '%s/lib', '%s/core/pydis/python']", XSTRING(SOURCE_ROOT), XSTRING(SOURCE_ROOT), XSTRING(SOURCE_ROOT));
  PyRun_SimpleString(command);
  PyRun_SimpleString("[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]");

  PyObject* calforce_module = get_calforce_module();
  PyObject* calforce = PyObject_CallObject(calforce_module, NULL);
  PyObject* mu = PyObject_GetAttrString(calforce, "mu");
  printf("calforce.mu = %g\n", PyFloat_AsDouble(mu));
  PyObject* nu = PyObject_GetAttrString(calforce, "nu");
  printf("calforce.nu = %g\n", PyFloat_AsDouble(nu));
  PyObject* a = PyObject_GetAttrString(calforce, "a");
  printf("calforce.a = %g\n", PyFloat_AsDouble(a));

  Py_Finalize();

  return 0;
}

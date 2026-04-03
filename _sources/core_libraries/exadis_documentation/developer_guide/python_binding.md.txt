
## Binding modules to python

Binding of C++ functions to python (`pyexadis`) is done using [pybind11](https://pybind11.readthedocs.io) and implemented in the `python/` directory of the ExaDiS folder. Wrappers are implemented on both the C++ and the python sides, in files `exadis_pybind.h`/`exadis_pybind.cpp` and `pyexadis_base.py`, respectively. Any C++ ExaDiS class or function must be explicitly wrapped into a corresponding binding object in order to be exposed to python in the `pyexadis` library. Whenever possible, the binding should be implemented locally in the module files using dedicated registry macros, so that file `python/exadis_pybind.cpp` does not have to be modified.

### Example 1: binding a simple class

For instance, a class `MyExaDisClass` implemented in C++ ExaDiS
```cpp
#include "exadis_pybind.h"

class MyExaDisClass {
public:
    MyExaDisClass() {
        // Initialize
    }
    void foo() {
        // Do some stuff
    }
};
```
must be binded to the `pyexadis` module using the `EXADIS_MODULE_BINDING()` registry macro (placed after the class definition):
```cpp
EXADIS_MODULE_BINDING(MyExaDisClass) {
    py::class_<MyExaDisClass>(m, "MyExaDisClass")
        .def(py::init<>())
        .def("foo", &MyExaDisClass::foo, "Do some stuff");
}
```
in order to become available on the python side:
```Python
import pyexadis
myobject = pyexadis.MyExaDisClass()
myobject.foo()
```

### Example 2: binding an ExaDiS module

For an ExaDiS module implemented in C++ that operates on an ExaDiS `System` object, e.g.

```cpp
class MyExaDisModule {
public:
    MyExaDisModule(System* system) {
        // Initialize
    }
    void foo(System* system) {
        // Do some stuff
    }
};
```

the binding must also act as a wrapper for the `System` object, which is only internal to ExaDiS and not directly exposed to the python side. This is because, in the context of OpenDiS, one should allow the input system to be provided from an arbitrary core library (e.g PyDiS). One way to proceed is to create a wrapper binding struct on the C++ side and expose this struct to `pyexadis`, e.g.:

```cpp
#include "exadis_pybind.h"

struct MyExaDisModuleBind {
    MyExaDisModule* mymodule;
    Params params;
    MyExaDisModuleBind(Params& p) {
        params = p;
        System* system = make_system(new SerialDisNet(), Crystal(), params);
        mymodule = new MyExaDisModule(system);
        exadis_delete(system);
    }
    void foo(ExaDisNet& disnet) {
        System* system = disnet.adjust_system(params);
        mymodule->foo(system);
    }
};

EXADIS_MODULE_BINDING(MyExaDisModule) {
    py::class_<MyExaDisModuleBind>(m, "MyExaDisModule")
        .def(py::init<Params>(), py::arg("params"))
        .def("foo", &MyExaDisModuleBind::foo, "Do some stuff", py::arg("net"));
}
```

such that the module now becomes available on the python side:

```Python
import pyexadis
from pyexadis_base import get_exadis_params
# Initialization
state = {...}
N = DisNetManager(...)
# Instantiate the module
params = get_exadis_params(state)
mymodule = pyexadis.MyExaDisModule(params)
# Call module method on ExaDiS network object
G = N.get_disnet(ExaDisNet)
mymodule.foo(G.net)
```


## Binding modules to python

Binding of C++ functions to python (`pyexadis`) is done using [pybind11](https://pybind11.readthedocs.io) and implemented in the `python/` directory of the ExaDiS folder. Wrappers are implemented on both the C++ and the python sides, in files `exadis_pybind.h`/`exadis_pybind.cpp` and `pyexadis_base.py`, respectively. Any C++ ExaDiS class or function implemented in `src/` must be explicitly wrapped into a corresponding binding object in file `exadis_pybind.cpp` in order to be exposed to python in the `pyexadis` library.

### Example 1: binding a simple class

For instance, a class `MyExaDisClass` implemented in C++ ExaDiS
```cpp
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
must be binded to the `pyexadis` module in `exadis_pybind.cpp` via:
```cpp
py::class_<MyExaDisClass>(m, "MyExaDisClass")
    .def(py::init<>())
    .def("foo", &MyExaDisClass::foo, "Do some stuff");
```
in order to become available on the python side:
```Python
import pyexadis
myobject = pyexadis.MyExaDisClass()
myobject.foo()
```

### Example 2: binding an ExaDiS module

For an ExaDiS module `MyExaDisModule` implemented in C++ that operates with ExaDiS `System` objects, e.g.

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

the binding must also act as a wrapper for the `System` object, which is only internal to ExaDiS and not directly exposed to the python side. This is because, in the context of OpenDiS, one should allow the input system to be provided from an arbitrary core library (e.g PyDiS). One way to proceed is to create a wrapper binding struct on the C++ side and expose this struct to `pyexadis`, e.g. in `exadis_pybind.cpp`:

```cpp
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
        System* system = disnet.system;
        system->params = params;
        mymodule->foo(system);
    }
};

py::class_<MyExaDisModuleBind>(m, "MyExaDisModule")
    .def(py::init<Params>(), py::arg("params"))
    .def("foo", &MyExaDisModuleBind::foo, "Do some stuff", py::arg("net"));
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
# Call module method on network object
G = N.get_disnet(ExaDisNet)
mymodule.foo(G.net)
```


### Example 3: binding a mobility law

Binding a C++ mobility law, e.g.
```cpp
struct MobilityMOBNAME
{
    struct Params {
        double drag = 1.0;
        Params() {}
        Params(double _drag) : drag(_drag) {}
    };
    // Mobility law implementation ...
    static constexpr const char* name = "MobilityMOBNAME";
};
namespace MobilityType {
    typedef MobilityLocal<MobilityMOBNAME> MOBNAME;
}
```
to `pyexadis` requires the following 4 steps:

1. Register the mobility parameters in `exadis_pybind.cpp`
```cpp
py::class_<MobilityType::MOBNAME::Params>(m, "Mobility_MOBNAME_Params")
    .def(py::init<double>(), py::arg("drag"));
```

2. Define the python mobility maker in `exadis_pybind.cpp`
```cpp
m.def("make_mobility_mobname", &make_mobility<MobilityType::MOBNAME>, "Instantiate a MOBNAME mobility law",
      py::arg("params"), py::arg("mobparams"));
```

3. Register the mobility for the `TopologyParallel` module in `exadis_pybind.h` (optional, only required to make the mobility compatible with `TopologyParallel`)
```cpp
} else if (strcmp(mobility->name(), "MobilityMOBNAME") == 0) {
    topology = new TopologyParallel<F,MobilityType::MOBNAME>(system, force, mobility, topolparams);
```

4. Add the mobility model to the `MobilityLaw` wrapper class in `pyexadis_base.py`
```python
elif self.mobility_law == 'MOBNAME':
    drag = kwargs.get('drag', 1.0)
    mobparams = pyexadis.Mobility_MOBNAME_Params(drag)
    self.mobility = pyexadis.make_mobility_mobname(params=params, mobparams=mobparams)
```

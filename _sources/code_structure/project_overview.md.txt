## Project Structure

The following figure illustrates the overall structure of the OpenDiS project
```{figure} OpenDiS_overview.png
:alt: Overall structure of OpenDiS project
:width: 640px
```

### Description
The objective of OpenDiS is to provide a framework to facilitate and support innovative research in dislocation dynamics. It consists of a central layer that provides higher level specifications, utilities, and helper classes, and of core libraries that implement fundamental functions for performing dislocation dynamics simulations. The OpenDiS framework also provides an opportunity for standardization, thereby enhancing simulation reproducibility and data analysis efficiency.

### OpenDiS layer
- This layer acts as the backbone of the OpenDiS framework and provides specifications, interfaces, drivers, utilities, and continuous integration features.
- It ensures consistency and usability across the various components.
- The OpenDiS layer is written in python to simplify the interaction between the different core modules and/or external libraries, and to facilitate user-defined enhancements.
- The specifications define the overall rules and standards for module implementations (e.g. as a set of abstract python classes), ensuring compatibility and integration across the different modules.

### Core libraries
- The core libraries provide the computational engine for dislocation dynamics simulations and implement functions, algorithms, and solvers in the form of modules (e.g. mobility law module, time-integration module, etc.)
- OpenDiS provides the following two native core libraries:
    - [**pydis**](../core_libraries/pydis_documentation/index.rst): a Python-based library for learning, prototyping, and running smaller-scale simulations.
    - [**exadis**](../core_libraries/exadis_documentation/index.rst): a C++-based library with Python interface developed for high-performance and GPU computing, intended for large-scale and production simulations.
- The modular design of OpenDiS allows users to develop, couple, and contribute their own core libraries to enhance the capabilities of the existing framework.

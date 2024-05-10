## How to get the code

### Include ExaDiS as submodule

To include ExaDiS as a submodule of OpenDiS, use the following steps:

```bash
cd ~/Codes
git clone --recurse-submodule https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
```

Alternatively, you can use the following commands to achieve the same:

```bash
git clone https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
git submodule update --init --recursive
```

### Getting them separately

If you prefer to get the submodules separately, follow these steps:

```bash
git clone https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
cd OpenDiS.git
cd core
git clone https://github.com/LLNL/exadis exadis
cd pydis/python
git clone https://github.com/davidjamesca/ctypesgen.git
cd ../..
cd exadis
git clone https://github.com/kokkos/kokkos.git --branch 4.2.00
cd python
git clone https://github.com/pybind/pybind11.git
```

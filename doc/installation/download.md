## How to get the code

### Download OpenDiS together with submodules

Use the following steps to download OpenDiS in your ``${HOME}/Codes`` folder:

```bash
mkdir -p ~/Codes
cd ~/Codes
git clone --recurse-submodule https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
```

Alternatively, you can use the following commands to achieve the same:

```bash
mkdir -p ~/Codes
cd ~/Codes
git clone https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
git submodule update --init --recursive
```

### Getting OpenDiS and submodules separately

While we do not recommend the following approach, you could download OpenDiS without the submodules and then download its submodules individually, e.g. using the following commands.

```bash
mkdir -p ~/Codes
cd ~/Codes
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

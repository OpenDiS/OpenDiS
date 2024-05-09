# Installation 


# How to get the code

## Include ExaDiS as submodule

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

## Getting them separately

If you prefer to get the submodules separately, follow these steps:

```bash
git clone https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
cd OpenDiS.git
cd core
git clone https://github.com/LLNL/exadis exadis
cd exadis
git clone https://github.com/kokkos/kokkos.git --branch 4.2.00
cd kokkos
```

# System Requirements

## Running Python examples in OpenDiS

You can run some OpenDiS examples that do not require any compilation, as long as you have Python3 installed. Refer to *Section 3.1* for more details.

## Requirements for high performance

For high-performance applications, you will need to compile ExaDiS and KOKKOS. Ensure the following packages are installed on your system before beginning to build ExaDiS/KOKKOS:

- **cmake**: version ?? or higher
- **gcc**: version ?? or higher
- **cuda**: version ?? or higher (if you want to run on GPU)
- **Python development package**
- **FFTW**

# Compile and Install

## Compile on Mac

### Build ExaDiS/KOKKOS

```bash
cd ~/Codes/OpenDiS.git/
rm -rf build/; ./configure.sh -DSYS=mac
cmake --build build -j 8 ; cmake --build build --target install
```

```


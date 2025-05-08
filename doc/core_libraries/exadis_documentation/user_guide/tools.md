## ExaDiS tools

```{note}
A set of tools is implemented in folder `tools/`. To compile them, OpenDiS/ExaDiS must be built with compilation flag `-DEXADIS_BUILD_TOOLS=On`.
```

### Computing internal stress fields
Tool `exadis_stress_field` in folder `bin/` (or `core/exadis/bin/`) allows to compute the internal stress field on a regular grid resulting from all dislocations in a configuration. It uses the non-singular dislocation expression ([Ref](https://doi.org/10.1016/j.jmps.2005.09.005)) to evaluate the stress field associated with each dislocation segment.

#### Usage
```
./exadis_stress_field config.data -MU 50e9 -NU 0.3 -a 1.0 -N 32
```

#### Options
* `-MU`, `--shear_modulus` (required): Shear modulus (Pa)
* `-NU`, `--poisson_ratio `(required): Poisson ratio
* `-a`, `--core_radius` (required): Dislocation core radius (b)
* `-N`, `--ngrid` (optional): Grid resolution [N] to compute the stress field. Default: `32`.
* `-Nxyz`, `--ngridxyz` (optional): Grid resolution [Nx,Ny,Nz] to compute the stress field. Default: none.
* `-o`, `--output`" (optional): Output file name. Default: `stress.dat`.
* `-Nimg`, `--nimages` (optional): Number of periodic images to compute the stress field. Default: `1`.
* `-pbc`, `--pbc_flags` (optional): Periodic boundary conditions along the 3 directions. Default: `1 1 1`.
* `-reg`, `--regularize_convergence` (optional): Apply regularization of the conditional convergence. Default: `1`.


### Computing displacement gradient fields
Tool `exadis_dispgrad_field` in folder `bin/` (or `core/exadis/bin/`) allows to compute the displacement gradient field on a regular grid resulting from all dislocations in a configuration. It uses the non-singular dislocation expression ([Ref](https://doi.org/10.1016/j.commatsci.2018.01.037)) to evaluate the displacement gradient field associated with each dislocation segment.

#### Usage
```
./exadis_dispgrad_field config.data -NU 0.3 -a 1.0 -N 32
```

#### Options
* `-NU`, `--poisson_ratio `(required): Poisson ratio
* `-a`, `--core_radius` (required): Dislocation core radius (b)
* `-N`, `--ngrid` (optional): Grid resolution [N] to compute the stress field. Default: `32`.
* `-Nxyz`, `--ngridxyz` (optional): Grid resolution [Nx,Ny,Nz] to compute the stress field. Default: none.
* `-o`, `--output`" (optional): Output file name. Default: `stress.dat`.
* `-Nimg`, `--nimages` (optional): Number of periodic images to compute the stress field. Default: `1`.
* `-pbc`, `--pbc_flags` (optional): Periodic boundary conditions along the 3 directions. Default: `1 1 1`.
* `-reg`, `--regularize_convergence` (optional): Apply regularization of the conditional convergence. Default: `1`.

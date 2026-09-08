# ParaDiS extension setup

The experimental `examples/06_paradis` example requires a separate ParaDiS
build with Python wrappers. Building OpenDiS alone does not generate `Home.py`.
The ParaDiS version must provide the `libparadis.so` build target, Python
bindings, and the API used by the example; an arbitrary ParaDiS release may
not be compatible.

Set `PARADIS_DIR` to that ParaDiS source tree, then build its library and wrappers:

```sh
export PARADIS_DIR=/absolute/path/to/ParaDiS
make -C "$PARADIS_DIR/src" libparadis.so
make -C "$PARADIS_DIR/python"
```

Verify that the build produced these files:

- `lib/libparadis.so`
- `lib/Home.py`
- `python/paradis_util.py`

From `extensions/paradis`, run `make install` to link them into OpenDiS.
`PARADIS_DIR` may also be a path relative to `extensions/paradis`; the installer
converts it to an absolute path so the links resolve from their destination
directories. Installation stops before creating links if a required file is
missing. These commands require a POSIX shell and symbolic link support.

Then change to `examples/06_paradis` and run `make`. The example imports the
bindings from the extension directories. A `ModuleNotFoundError` for `Home`
usually means that the wrapper has not been built or linked successfully.

This setup only installs the CPU bindings. The optional GPU bindings and
version-specific paths inside the ParaDiS wrappers need separate setup.
Before running the example, also update `fmCorrectionTbl` in
`paradis_default.ctrl` to the correction table supplied by your ParaDiS tree.
The default value assumes a sibling directory named `ParaDiS.git`.

Run the installer regression tests from the OpenDiS root with:

```sh
python3 -m unittest discover -s tests/test_paradis_install -v
```

The tests use temporary placeholder files to check installation, not a full
ParaDiS simulation. They require GNU Make and symbolic link support; set
`MAKE` to the GNU Make executable if it is not named `make`.

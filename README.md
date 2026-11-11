
ChemPotPy, CHEMical library of POTential energy surfaces in PYthon 
==================================================================

Jan. 17, 2025

Authors: Yinan Shu, Zoltan Varga, Dayou Zhang, Donald G. Truhlar
University of Minnesota, Minnesota, United States

Fork modernization author and maintainer: Ruben Meana-Paneda

ChemPotPy is a library for analytic representation of single-state 
and multi-state potential energy surfaces and couplings. 

All fortran source code are stored in folder chempotpy 


About this fork
---------------

This repository is the [rmeanapa CHEMPOTPY fork](https://github.com/rmeanapa/chempotpy)
maintained by Ruben Meana-Paneda. It retains the original scientific authors
and citations above while developing a self-contained, modern implementation.
The comparison below is against
[shuyinan/chempotpy](https://github.com/shuyinan/chempotpy) `main` at commit
`4716d22` (October 2, 2026). This fork includes that upstream history and adds
the following changes:

* Converted the complete PES catalog from fixed-form `.f` sources to
  authoritative free-form `.f90` sources and removed obsolete duplicate
  implementations.
* Replaced implicit typing with `implicit none` and explicit declarations;
  introduced private-by-default modules, explicit public APIs, typed
  parameters, `intent(...)`, and documented procedures.
* Removed legacy `COMMON`, `DATA`, mutable `SAVE`, `ENTRY`, `BLOCK DATA`,
  setup-style `prepot`, numeric statement labels, `CONTINUE`, `GOTO`,
  arithmetic `IF`, detached `FORMAT`, and redundant terminal `RETURN`
  constructs from governed PES sources.
* Replaced initialization entry points with stateless, independently owned PES
  modules whose fixed coefficients live in source-controlled module data.
* Added a bundled Fortran core for numerical kinds, status values, canonical
  and compatibility unit conversions, Cartesian geometry, coordinate
  flattening, internal-to-Cartesian gradient transformations, reusable
  internal-coordinate adapters, H3 basis operations, BLAS operations, and
  symmetric eigensolution. The core keeps CHEMPOTPY
  independent of external MKL, LAPACK, and BLAS libraries.
* Added a Python core contract, adapters, surface registry, runtime dispatch,
  and unit utilities. Python implementations are organized under
  `chempotpy/core/python`; Fortran core modules live directly under
  `chempotpy/core`.
* Standardized the Fortran source format: four-space indentation, compact array
  subscripts, consistent declarations and assignments, `elseif`/`endif` and
  `enddo`, inline I/O formats, meaningful PES-specific documentation, and a
  120-column readability target.
* Added catalog discovery and isolated f2py builds so each PES extension
  receives only the bundled core modules it actually imports.
* Added strict compilation, modernization, formatting, packaging, binary
  linkage, extension-import, and numerical-regression gates for all 363
  callable surfaces, including two stored reference geometries per surface.
* Added a reproducible current conda-forge environment and validation against
  the newest mutually compatible Python, NumPy, Meson, Ninja, and gfortran
  toolchain.
* Documented the shared PES contract, Fortran conventions, build workflow, and
  commands for validating the entire surface catalog.


How to install
--------------
* Create a virtual environment with the newest mutually compatible Python,
  Conda, NumPy, Meson, Ninja, and gfortran packages from conda-forge:
    
       conda env create --file environment.yml
       conda activate chempotpy

* Install stable release:
  
        pip install chempotpy



* Users may compile all extension modules for their current Python and
  platform:

  get into the parent directory of chempotpy 
   
        make build
        make check-modernization

  `make build` uses an unoptimized compatibility profile so the modern Python
  extensions preserve the released Fortran results. Add `--debug` to enable
  bounds checks and Fortran backtraces while diagnosing a kernel.

  `make check-conda BUILD_FLAGS="--jobs 4"` builds every extension, checks its
  native-library linkage, imports each binary in an isolated Python process,
  and evaluates all surfaces at both catalog regression geometries.

  CHEMPOTPY carries the numerical support routines used by its PES kernels.
  Building the extensions does not require MKL, LAPACK, or BLAS. A Fortran
  compiler is needed only when building extensions from source; binary wheels
  should carry the compiler runtime required by their platform.

  Modern PES sources may import the bundled Fortran core for numerical kinds,
  canonical or compatibility unit conversions, Cartesian geometry operations,
  status values, and linear algebra. The isolated builder follows those module
  dependencies and compiles the required sources into each extension, so every
  resulting PES module remains independently loadable.

Modernizing Fortran sources
---------------------------

Modernized PES code uses four-space indentation, spaced procedure arguments,
short `!> @brief` documentation, explicit declarations and `implicit none`,
`private` by default in modules, and `use ..., only: ...` imports. Run
`make check-style` to validate every free-form Fortran source. Compilation also
uses `-fimplicit-none`, so an undeclared variable is a hard error.

Declarations use `::`, constants use typed `parameter` declarations, and
intrinsic kinds and character lengths use standard syntax. Source code is
lowercase without changing comments or character literals. Files acquire units
with `open(newunit=...)`; numbered diagnostic writes are retained only where a
legacy public routine still exposes that behavior.

The formatter uses a 120-column readability target for normal equations rather
than Fortran 77's 72/80-column limit. Generated numerical data may exceed that
target, and all free-form compilation uses unlimited line length.

Within a procedure, named constants needed for array shapes come first,
followed by `intent(...)` dummy arguments and then local variables. Redundant
runs of empty lines are collapsed to one visual separator.

Format specifications are written inline in `read` and `write` statements.
Legacy labeled `FORMAT` statements and redundant `return` statements directly
before a procedure end are rejected by the style gate; intentional early
returns remain valid.

All PES kernels are now authoritative free-form `.f90` sources; no legacy `.f`
implementation remains in the package. The regression suite recompiles every
converted kernel and compares it with values captured from the released
fixed-form source. The detailed contract, migration, and numerical-parity rules
are in `docs/PES_CONTRACT.md`; the required source conventions are collected in
`docs/FORTRAN_STYLE.md`.

Run `make check-modernization` for the complete repeatable modernization gate.
It audits the catalog, rejects forbidden legacy Fortran constructs, runs the
Python and core tests, strictly compiles all 363 callable PES interfaces, checks
formatting, and evaluates every surface at two deterministic geometries in the
same process. `make test-surfaces` runs only the slower catalog-wide numerical
and architecture regression portion.


Generate a standalone Fortran PES
---------------------------------

The standalone exporter selects one catalog surface and writes a single
dependency-complete `.f90` file. The file contains the PES, any surface support
source, and only the bundled Fortran core modules that the PES imports. It can
therefore be compiled directly or added to another Fortran program without
installing CHEMPOTPY, MKL, LAPACK, or BLAS.

For example, export the H3 A2 surface from the repository root:

    python scripts/export_standalone.py \
        --surface H3/H3_GEN_A2_2002 \
        --output standalone/H3_GEN_A2_2002

This creates:

    standalone/H3_GEN_A2_2002/H3_GEN_A2_2002_standalone.f90

Compile that file with the program that will call its historical `pes` or
`dpem` entry point. Free-form sources may contain long fitted-coefficient lines,
so gfortran should be given unlimited free-form line length:

    gfortran -O2 -ffree-form -ffree-line-length-none \
        H3_GEN_A2_2002_standalone.f90 my_program.f90 -o my_program

The generated source preserves the selected surface's public argument list and
historical units, making it suitable for existing dynamics programs. The
surface source documents the exact `pes(...)` or `dpem(...)` signature. A
surface that still reads external scientific data requires its data directory
at runtime; the exporter does not hide or rewrite that input. N4_1A_PIPNN,
FH2O_PIPNN, FH2O_SOC_PIPNN, both CH4F 1996 surfaces, CH4Br_VBMM_2002, and
F2H2_GEN_JBKKL_1990, HOBr_3App_PIP, COCO_PIPNN, and C3H7NO_PIP retain their historical path argument for ABI
compatibility, but their immutable parameters are embedded in the exported
source and those surfaces perform no runtime model-file I/O.

Add `--with-driver` to also generate a complete example program:

    python scripts/export_standalone.py \
        --surface H3/H3_GEN_A2_2002 \
        --output standalone/H3_GEN_A2_2002 \
        --with-driver \
        --geometries 2

Compile and run the generated pair:

    cd standalone/H3_GEN_A2_2002
    gfortran -O2 -ffree-form -ffree-line-length-none \
        H3_GEN_A2_2002_standalone.f90 \
        H3_GEN_A2_2002_example.f90 \
        -o H3_GEN_A2_2002_example
    ./H3_GEN_A2_2002_example

The example constructs deterministic, non-symmetric Cartesian geometries and
reports each adiabatic energy, or each diabatic-matrix element for a DPEM. It
also reports Cartesian gradients and Hessians calculated by centered finite
differences, so the example works even when the original PES publishes energy
only. The default displacement is `1.0e-3` angstrom and can be changed with
`--difference-step`. Consequently, these example derivatives are numerical
diagnostics in the historical output units per angstrom and per angstrom
squared; they do not claim analytic Hessian support from the underlying PES.

For a surface with a data-directory argument, pass that directory to the
example executable:

    ./SURFACE_example /absolute/path/to/chempotpy

The directory argument is the package data root that contains the individual
system directories, not the system directory itself.

Surface names that occur in more than one system must use the unambiguous
`SYSTEM/SURFACE` selector shown above.


Citation
--------

The following paper should be cited in publications utilizing the
ChemPotPy library in addition to the original paper that publishes 
the potential energy surface subroutine:

Shu, Y.; Varga, Z.; Truhlar, D. G.
"ChemPotPy: A Python Library for Analytic Representation of Potential 
Energy Surfaces and Diabatic Potential Energy Matrices"
accepted by J. Phys. Chem. A


TO CONTRIBUTE YOUR POTENTIAL or report a bug
--------------------------------------------
* Option 1
Send email to one of the maintainers:
  - Yinan Shu, yinan.shu.0728@gmail.com
  - Zoltan Varga, zoltan78varga@gmail.com
  - Dayou Zhang, zhan6350@umn.edu
  - Donald G. Truhlar, truhlar@umn.edu
 
* Option 2
Submit tickets on the [issues](https://github.com/shuyinan/chempotpy/issues)


Validating every potential energy surface
-----------------------------------------

CHEMPOTPY provides two complementary numerical regression paths. The Python
extension gate exercises the same compiled modules used by package users,
while the direct-Fortran test independently compiles the authoritative source
of every PES. Both paths evaluate two deterministic geometries per surface and
compare the results with `tests/data/all_surface_references.json`.

Create the current conda-forge validation environment once, then activate it:

    conda env create --file environment.yml
    conda activate chempotpy

From the repository root, run the complete Python extension gate with:

    make check-conda BUILD_FLAGS="--jobs 4"

This command builds all 363 f2py extensions, checks their native-library
linkage, imports every extension in an isolated Python process, and performs
726 numerical evaluations. A successful run ends with counters equivalent to:

    pes_extensions_checked: 363
    missing_or_ambiguous_extensions: 0
    external_mkl_lapack_blas_links: 0
    python_extensions_checked: 363
    python_extension_evaluations: 726
    python_extension_failures: 0

Run the independent direct-Fortran numerical and architecture regression with:

    make test-surfaces

To run only the all-surface numerical comparison from that suite:

    python -m unittest \
        tests.test_all_surfaces.AllSurfaceRegressionTests.test_every_surface_matches_its_reference_outputs \
        -v

An individual extension can be built and evaluated by its system and surface
name. For example:

    python scripts/build_pes.py --surface H3_GEN_A2_2002
    python scripts/check_python_extensions.py \
        --output-root chempotpy \
        --surface H3/H3_GEN_A2_2002

Run every repository test with:

    python -m unittest discover -s tests -v

The extension gate places generated binaries in the package tree by default.
Remove them after local validation with:

    make clean

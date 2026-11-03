
ChemPotPy, CHEMical library of POTential energy surfaces in PYthon 
==================================================================

Jan. 17, 2025

Authors: Yinan Shu, Zoltan Varga, Dayou Zhang, Donald G. Truhlar
University of Minnesota, Minnesota, United States

ChemPotPy is a library for analytic representation of single-state 
and multi-state potential energy surfaces and couplings. 

All fortran source code are stored in folder chempotpy 


How to install
--------------
* Create a virtual environment with Python, NumPy, Meson, and gfortran:
    
       conda create --name chempotpy
       conda activate chempotpy
       conda install python=3.11
       conda install -c conda-forge gfortran
       pip install numpy meson

* Install stable release:
  
        pip install chempotpy



* Users may compile all extension modules for their current Python and
  platform:

  get into the parent directory of chempotpy 
   
        make build
        make check

  `make build BUILD_FLAGS="--debug --jobs 4"` provides a faster, unoptimized
  development build. The default remains the optimized release build.

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

Run `make test-surfaces` for the slower catalog gate, which independently
compiles and evaluates all 363 callable PES interfaces at two deterministic
geometries.


Compile Chempotpy subroutine
----------------------------
One can use the meta programming script to automatically generate a 
fortran subroutine. 

  create a conda vitural environment as suggested in How to install

  get into the parent directory of chempotpy/chempotpy

       ./meta_chempotpy.script

The meta program will generate a fortran subroutine called chempotpy. 
One can interface this chempotpy subroutine with any dynamics code. 
Notice that the meta program will also generate sub programs for each 
surface that is located in chempotpy/chempotpy/system/lib/.


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

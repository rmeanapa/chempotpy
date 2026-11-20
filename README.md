
ChemPotPy, CHEMical library of POTential energy surfaces in PYthon 
==================================================================

[![All PES numerical regression](https://github.com/rmeanapa/chempotpy/actions/workflows/all-pes-numerical-regression.yml/badge.svg?branch=main)](https://github.com/rmeanapa/chempotpy/actions/workflows/all-pes-numerical-regression.yml?query=branch%3Amain)

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
* Added a bundled Fortran core for numerical kinds, status values, typed PES
  requests/results/metadata, canonical NIST CODATA 2022 and compatibility unit
  conversions, Cartesian geometry, coordinate flattening and reshaping,
  internal-to-Cartesian gradient transformations, centered finite-difference
  validation, reusable internal-coordinate adapters, BLAS operations, and
  symmetric eigensolution. The core keeps CHEMPOTPY
  independent of external MKL, LAPACK, and BLAS libraries.
* Added a Python core contract, adapters, surface registry, runtime dispatch,
  and unit utilities. Python implementations are organized under
  `chempotpy/core/python`; Fortran core modules live directly under
  `chempotpy/core`.
* Kept system- or fit-family support beside its PES sources. For example,
  `chempotpy/H3/chempotpy_h3_2002_basis.f90` owns the shared basis used by the
  A2, A3, A4, and CCI H3 2002 surfaces and is bundled automatically.
* Standardized the Fortran source format: four-space indentation, compact array
  subscripts, consistent declarations and assignments, `elseif`/`endif` and
  `enddo`, inline I/O formats, meaningful PES-specific documentation, and a
  120-column readability target.
* Added catalog discovery and isolated f2py builds so each PES extension
  receives only the bundled core modules it actually imports.
* Added strict compilation, modernization, formatting, packaging, binary
  linkage, extension-import, and numerical-regression gates for all 363
  callable surfaces, including nine stored reference geometries per surface.
* Corrected the CH3NH2 APRP coupling evaluation by initializing terms that are
  absent from the coupling branch and conditionally assigned tent derivatives.
  The legacy result depended on uninitialized stack memory and therefore
  changed with compiler, architecture, and source layout.
* Added a reproducible current conda-forge environment and validation against
  the newest mutually compatible Python, NumPy, Meson, Ninja, and gfortran
  toolchain.
* Documented the shared PES contract, Fortran conventions, build workflow, and
  commands for validating the entire surface catalog.


Mathematical constants and unit conversions
-------------------------------------------

The generated `chempotpy_constants` module provides one working-precision
definition of `PI` together with the NIST CODATA 2022 universal physical
constants. The authoritative CODATA decimal values and provenance are stored
in `chempotpy/core/codata_2022.json`; compiled PES modules therefore do not
read a data file at runtime. Run the following check after changing the
manifest or generator:

    python scripts/generate_constants.py --check

All public PES boundary conversions use the resulting canonical factors:

| Conversion | Value |
| --- | ---: |
| bohr to angstrom | `0.529177210544` |
| hartree to electronvolt | `27.211386245981` |
| hartree/bohr to eV/angstrom | `51.42206751119799` |
| hartree to kJ/mol | `2625.499639479163` |
| hartree to kcal/mol | `627.5094740628975` |
| hartree to cm^-1 | `219474.63136314` |
| electronvolt to kcal/mol | `23.06054783061903` |
| electronvolt to cm^-1 | `8065.543937349212` |

The source values are published in the
[NIST CODATA 2022 extensive listing](https://physics.nist.gov/cuu/pdf/all.pdf).
The module `chempotpy_legacy_units` remains available only for explicitly
requested reproduction of the formerly released rounded conversions. A
conversion-like value that is part of a published PES fit remains local to its
surface and is not replaced merely because it resembles a physical constant;
changing such a value would change the fitted model itself.

The C6H5SH APRP pair and the three C7H8S APRP adiabatic/diabatic pairs are
documented exceptions for input-coordinate scaling. Their released fits are
numerically sensitive enough that replacing the historical scale changed
representative energies by as much as about `0.18 eV`. They therefore retain
that scale as an explicitly named published-fit constant while using CODATA
2022 for energy, gradient, coupling, and all other public conversions.


Fortran working precision
-------------------------

Modern numerical code imports `wp` from `chempotpy_kinds` and declares real
values as `real(wp)` with explicit literals such as `1.0_wp`. The kinds module
defines portable single, double, and quadruple requests with
`selected_real_kind` and currently selects double precision:

```fortran
use chempotpy_kinds, only: wp
real(wp) :: energy

energy = 0.0_wp
```

The selector is the single `wp = dp` parameter in
`chempotpy/core/chempotpy_kinds.f90`. The shared numerical core, high-level
linear-algebra and interpolation APIs, private BLAS/DSYEV implementation,
generated constants and parameter modules, and the maintained PES catalog now
use this selector consistently. Package builds must
currently leave `wp = dp` because the Python/f2py boundary maps `wp` to a C
double. Standalone Fortran kernels may select another supported kind; a
catalog-wide non-binary64 Python build still requires an explicit conversion
layer at that external boundary.


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
  typed PES contracts, finite-difference validation, status values, and linear
  algebra. The isolated builder follows those module dependencies and compiles
  the required sources into each extension, so every resulting PES module
  remains independently loadable.

Modernizing Fortran sources
---------------------------

Modernized PES code uses four-space indentation, spaced procedure arguments,
short `!> @brief` documentation, explicit declarations and `implicit none`,
`private` by default in modules, and `use ..., only: ...` imports. Run
`make check-style` to validate every free-form Fortran source. Compilation also
uses `-fimplicit-none`, so an undeclared variable is a hard error.

Declarations use `::`, constants use typed `parameter` declarations, and
intrinsic kinds and character lengths use standard syntax. Source code is
lowercase without changing comments or character literals. Library kernels are
quiet: they contain no external `write`, `stop`, or `error stop` statements.
The optional `chempotpy_status_t` result on each module-level `evaluate`
reports unsupported requests without changing the historical `pes` or `dpem`
adapter ABI.

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
converted kernel and compares it with the current CODATA 2022 boundary
references. Earlier release-reference fixtures remain as independent migration
evidence where documented. The detailed contract, migration, and
numerical-parity rules are in `docs/PES_CONTRACT.md`; the required source
conventions are collected in `docs/FORTRAN_STYLE.md`.

Run `make check-modernization` for the complete repeatable modernization gate.
It audits the catalog, rejects forbidden legacy Fortran constructs, runs the
Python and core tests, strictly compiles all 363 callable PES interfaces, checks
formatting, and evaluates every surface at nine deterministic geometries in the
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
documented units while using the versioned CODATA conversion contract, making
it suitable for existing dynamics programs. The surface source documents the
exact `pes(...)` or `dpem(...)` signature. Every released model parameter used
by an installed evaluator is embedded in the owning source or a PES-specific
generated module, so standalone exports perform no runtime model-file I/O.
Historical `path` arguments remain on some entry points only for ABI
compatibility and may be passed an empty string.

The POTLIB CH5p release references two long-range spline tables that neither
POTLIB 2023 nor the original ChemPotPy distribution includes. Its released
short-range PIP is self contained; entering the unavailable long-range branch
reports that provenance defect explicitly instead of searching the working
directory for untracked files.

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
differences through the bundled `chempotpy_finite_difference` module, so the
example works even when the original PES publishes energy only. The default
displacement is `1.0e-3` angstrom and can be changed with `--difference-step`.
Consequently, these example derivatives are numerical diagnostics in the
historical output units per angstrom and per angstrom squared; they do not
claim analytic Hessian support from the underlying PES.

If a historical entry point still accepts a data-directory argument, the
generated example accepts the argument for compatibility but does not use it.

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
of every PES. Both paths evaluate nine deterministic geometries per surface and
compare the results with `tests/data/all_surface_references.json`. The suite
uses three end-to-end pairs and the midpoint of each pair. It retains the two
original non-symmetric endpoints, adds a symmetric non-rigid distortion pair
around the 75/25 interior point, and adds a radial contraction/expansion pair
around the 25/75 interior point. The extra points exercise independent internal
coordinate directions while remaining inside or close to the established test
envelope.

Create the current conda-forge validation environment once, then activate it:

    conda env create --file environment.yml
    conda activate chempotpy

From the repository root, run the complete Python extension gate with:

    make check-conda BUILD_FLAGS="--jobs 4"

This command builds all 363 f2py extensions, checks their native-library
linkage, imports every extension in an isolated Python process, and performs
3,267 numerical evaluations. A successful run ends with counters equivalent to:

    pes_extensions_checked: 363
    missing_or_ambiguous_extensions: 0
    external_mkl_lapack_blas_links: 0
    python_extensions_checked: 363
    python_extension_evaluations: 3267
    python_extension_failures: 0

Run the catalog-wide compiler warning gate with:

    make check-fortran-warnings

This compiles all 396 PES and support units with `-Wall`, `-Wextra`, and
`-Werror`. Unused dummy arguments are exempt because the stable PES contracts
intentionally include inputs that some implementations do not consume. Exact
real comparisons are also exempt because released lookup and branch logic uses
them; both categories remain protected by the numerical regression suite.

Run the independent direct-Fortran numerical and architecture regression with:

    make test-surfaces

Every pull request and push also runs the `All PES numerical regression`
GitHub Actions workflow on Ubuntu Linux and macOS. The workflow checks out the
pull request's base commit or the most recent successful push, compiles that
revision and the current revision independently on both operating systems,
evaluates all nine deterministic geometries, and compares every result. A
failed run therefore cannot become the protected baseline. Capturing both
revisions with the same compiler and math runtime on each runner prevents
platform differences from being mistaken for a source change. Each revision
uses its own versioned geometry harness so an explicitly reviewed baseline
transition remains visible. Existing surfaces must remain within a relative
tolerance of `1e-11` and an absolute tolerance of `5e-10`; removed surfaces,
interface changes, non-finite results, compilation failures, and execution
failures also fail the workflow. Aln and NH3 use
deterministic system-specific geometries so their regression points remain
inside the mathematical domain of the released fits. A newly added surface
must compile, execute, match its declared interface, and return only finite
values; that revision becomes its protected baseline on the next comparison.

The SiOH variable-atom driver correction establishes `baseline_version: 4`.
The one-step transition from version 3 approves evaluating SiOH with its public
`xyz(natoms,3)` shape instead of a legacy oversized array whose storage layout
scrambled the regression geometry. Version 3 records the earlier C2H6F indexing
and uninitialized-state corrections, and version 2 records the CODATA 2022
Hartree-to-kcal/mol migration. Records whose approved values did not change
retain the base revision captured on that same runner, so an unchanged PES is
never compared with another architecture's fixture. Comparisons between two
version-4 revisions remain strict. The
committed fixture supports local and package-level regression tests, while the
workflow's base-revision capture protects existing surfaces independently of
ordinary fixture edits in a proposed change.

To reproduce the workflow comparison locally against a selected Git revision:

    base_root="$(mktemp -d)"
    git worktree add --detach "$base_root/checkout" REVISION
    (cd "$base_root/checkout" && \
        python scripts/capture_library_references.py --all --jobs 4 \
            --output "$base_root/baseline.json")
    python scripts/check_numerical_regression.py \
        --baseline "$base_root/baseline.json" --jobs 4
    git worktree remove --force "$base_root/checkout"

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

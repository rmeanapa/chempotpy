
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

New surfaces must follow the documented scientific, architectural, and
validation procedure in [CONTRIBUTING.md](CONTRIBUTING.md). The repository
provides `scripts/add_pes.py` to create a non-overwriting modern PES or DPEM
scaffold.


Potential energy surface catalog
--------------------------------

This table groups the callable Fortran registry by chemical system and gives a
short description of the available surface families. Regenerate it with
`python scripts/generate_readme_catalog.py` after adding or changing a surface.

<!-- BEGIN GENERATED PES CATALOG -->

The catalog currently exposes **403 callable surfaces** across **99 chemical systems**. Paired PES and DPEM entry points are counted separately because each is independently callable.

| System | Surface description | Number of surfaces |
| :--- | :--- | ---: |
| `ALC` | Extended-LEPS model surfaces for atom–linear-molecule collisions. | 12 |
| `AlmHn` | Valence-bond-order models for aluminum–hydrogen clusters. | 2 |
| `Aln` | Aluminum-cluster embedded-atom, extended-Rydberg, and many-body models. | 59 |
| `ArNO` | PIP neural-network surfaces. | 2 |
| `BrH2` | LEPS reactive surfaces with paired adiabatic and diabatic forms. | 7 |
| `BrHCl` | LEPS reactive surfaces. | 4 |
| `BrNH3` | Valence-bond molecular-mechanics surface. | 1 |
| `C2H2N2O` | PIP neural-network surface. | 1 |
| `C2H2O` | PIP neural-network surface. | 1 |
| `C2H4O4` | Permutationally invariant polynomial surface. | 1 |
| `C2H6Cl` | Valence-bond molecular-mechanics surface. | 1 |
| `C2H6F` | Valence-bond molecular-mechanics surface. | 1 |
| `C2H6H` | Valence-bond molecular-mechanics surface. | 1 |
| `C2H6O` | Valence-bond molecular-mechanics surface. | 1 |
| `C2H6OH` | Valence-bond molecular-mechanics surface. | 1 |
| `C2H7` | Valence-bond molecular-mechanics surface. | 1 |
| `C2N` | PIP neural-network surfaces. | 2 |
| `C3H7NO` | Permutationally invariant polynomial surface. | 1 |
| `C6H5SH` | Anchor-points reactive potential surfaces with paired adiabatic and diabatic forms. | 2 |
| `C6H6O` | Anchor-points reactive potential surfaces with paired adiabatic and diabatic forms. | 2 |
| `C7H8S` | Anchor-points reactive potential surfaces with paired adiabatic and diabatic forms. | 16 |
| `CH2O` | PIP neural-network surface. | 1 |
| `CH2O2` | PIP neural-network surfaces. | 2 |
| `CH2OH` | Surfgen multistate surfaces with paired adiabatic and diabatic forms. | 4 |
| `CH3NH2` | Anchor-points reactive potential surfaces with paired adiabatic and diabatic forms. | 2 |
| `CH3OH` | Permutationally invariant polynomial surface. | 1 |
| `CH4` | General analytic surface. | 1 |
| `CH4_H2O_H2O` | Permutationally invariant polynomial surface. | 1 |
| `CH4Br` | Valence-bond molecular-mechanics surface. | 1 |
| `CH4Cl` | Valence-bond molecular-mechanics surfaces. | 2 |
| `CH4CN` | Valence-bond molecular-mechanics surface. | 1 |
| `CH4F` | General analytic and LEPS reactive surfaces. | 3 |
| `CH4Ni` | Embedded-atom surface. | 1 |
| `CH4O` | LEPS reactive surfaces. | 4 |
| `CH4OCl` | PIP neural-network and valence-bond molecular-mechanics surfaces. | 2 |
| `CH4OF` | PIP neural-network surface. | 1 |
| `CH4OH` | Valence-bond molecular-mechanics surfaces. | 2 |
| `CH5` | General analytic, LEPS reactive, and valence-bond molecular-mechanics surfaces. | 5 |
| `CH5O2` | Valence-bond molecular-mechanics surface. | 1 |
| `CH5p` | Permutationally invariant polynomial surface. | 1 |
| `Cl2H` | LEPS reactive surfaces. | 4 |
| `ClH2` | LEPS reactive surfaces. | 15 |
| `ClH2O` | PIP neural-network surface. | 1 |
| `ClNH3` | Valence-bond molecular-mechanics surface. | 1 |
| `CO2` | RKHS-fitted surfaces. | 5 |
| `COCO` | PIP neural-network surface. | 1 |
| `F2H` | LEPS reactive surfaces. | 2 |
| `F2H2` | General analytic surfaces. | 3 |
| `FH2` | General analytic and LEPS reactive surfaces. | 25 |
| `FH2O` | PIP neural-network surfaces. | 2 |
| `GeH4OH` | Valence-bond molecular-mechanics surface. | 1 |
| `GeH5` | LEPS reactive surface. | 1 |
| `H2CO` | Permutationally invariant polynomial surface. | 1 |
| `H2O2` | General analytic surface. | 1 |
| `H2OBr` | Permutationally invariant polynomial surface. | 1 |
| `H3` | Benchmark analytic H₃ surfaces, including BKMP, DMBE, LSTH, and related fits. | 11 |
| `H3Cl` | PIP neural-network surface. | 1 |
| `H3ClO` | Permutationally invariant polynomial and PIP neural-network surfaces. | 2 |
| `H3m` | Analytic surface. | 1 |
| `H3p` | Analytic surfaces. | 2 |
| `H3S` | PIP neural-network surface. | 1 |
| `H4O2` | Permutationally invariant polynomial surface. | 1 |
| `H7p` | Permutationally invariant polynomial surface. | 1 |
| `HClOH` | PIP neural-network surface. | 1 |
| `HeH2p` | Diatomics-in-molecules and analytic surfaces. | 2 |
| `HNi` | Embedded-atom surfaces. | 3 |
| `HO3` | PIP neural-network surface. | 1 |
| `HOBr` | Permutationally invariant polynomial surface. | 1 |
| `HOCO_Ar` | Permutationally invariant polynomial surface. | 1 |
| `I2H` | LEPS reactive surfaces. | 3 |
| `IH2` | LEPS reactive surfaces. | 5 |
| `K2Rb2` | PIP neural-network surface. | 1 |
| `LiFH` | LEPS reactive surfaces with paired adiabatic and diabatic forms. | 6 |
| `LiH2p` | Analytic surfaces. | 2 |
| `malonaldehyde` | Analytic surface. | 1 |
| `MCH` | LEPS reactive surfaces with paired adiabatic and diabatic forms. | 10 |
| `MXH` | LEPS reactive surfaces with paired adiabatic and diabatic forms. | 6 |
| `N2HOCp` | PIP neural-network surface. | 1 |
| `N2O` | Ground- and excited-state N₂O fits, including PIP and ZV representations. | 10 |
| `N2O2` | Permutationally invariant polynomial and ZV analytic surfaces. | 2 |
| `N3` | ZV analytic and many-body polynomial surfaces. | 2 |
| `N4` | PIP, neural-network, many-body-expansion, and global analytic N₄ surfaces. | 6 |
| `NaFH` | LEPS reactive surfaces with paired adiabatic and diabatic forms. | 6 |
| `NaH2` | LEPS reactive surfaces with paired adiabatic and diabatic forms. | 14 |
| `NH3` | General analytic and Surfgen multistate surfaces with paired adiabatic and diabatic forms. | 6 |
| `NH3CN` | Valence-bond molecular-mechanics surface. | 1 |
| `NH3OH` | Valence-bond molecular-mechanics surface. | 1 |
| `NH4` | Valence-bond molecular-mechanics surfaces. | 2 |
| `NO2` | PIP neural-network and ZV analytic surfaces. | 6 |
| `O2H` | General analytic surface. | 1 |
| `O2H3` | PIP neural-network surface. | 1 |
| `O3` | Single- and multistate ozone surfaces in analytic, many-body, and diabatic forms. | 35 |
| `O4` | Singlet, triplet, and quintet O₄ many-body and analytic surfaces. | 7 |
| `OH2` | General analytic and LEPS reactive surfaces. | 9 |
| `OH3` | LEPS reactive and permutationally invariant polynomial surfaces with paired adiabatic and diabatic forms. | 5 |
| `SiH4Cl` | Valence-bond molecular-mechanics surface. | 1 |
| `SiH5` | LEPS reactive surfaces. | 2 |
| `SiOH` | General analytic surface. | 1 |
| `YRH` | LEPS reactive surfaces with paired adiabatic and diabatic forms. | 2 |

<!-- END GENERATED PES CATALOG -->


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
* Added a bundled Fortran core for numerical kinds, atomic elements and
  standard atomic masses, status values, typed PES requests/results/metadata,
  canonical NIST CODATA 2022 and compatibility unit conversions, Cartesian
  geometry, coordinate flattening and reshaping,
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
  linkage, extension-import, and numerical-regression gates for all 403
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

Atomic element data
-------------------

The `chempotpy_atomic_elements` core module owns the 118 canonical element
symbols, their atomic numbers, and standard atomic masses in unified atomic
mass units. Its case-insensitive `atomic_number`, `atomic_symbol`, and
`atomic_mass_u` queries return neutral values for invalid input so library
evaluators do not terminate the calling application. `D` and `T` are accepted
as hydrogen aliases for element-identity checks; isotope-specific masses are
outside this standard-element table.

PES sources must not carry private copies of the periodic table. Fit-specific
atom ordering and isotope choices remain with the owning surface, while shared
element queries use this core owner.

The C6H5SH APRP pair and the eight C7H8S APRP adiabatic/diabatic pairs are
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
  and evaluates all surfaces at all nine catalog regression geometries.

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
Python and core tests, strictly compiles all 403 callable PES interfaces, checks
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

This command builds all 403 f2py extensions, checks their native-library
linkage, imports every extension in an isolated Python process, and performs
3,618 numerical evaluations. A successful run ends with counters equivalent to:

    pes_extensions_checked: 403
    missing_or_ambiguous_extensions: 0
    external_mkl_lapack_blas_links: 0
    python_extensions_checked: 403
    python_extension_evaluations: 3618
    python_extension_failures: 0

Run the catalog-wide compiler warning gate with:

    make check-fortran-warnings

This compiles all 462 PES and system-support units with `-Wall`, `-Wextra`, and
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
must compile, execute, match its declared interface, and reproduce the finite,
independently reviewed values committed in the current fixture. The workflow
does not bootstrap a new baseline from the implementation being tested.

Removing the obsolete, unused C2H2O PIP-NN data-path argument establishes
`baseline_version: 5`. The one-step transition from version 4 approves only
that interface change; its numerical values are unchanged. Version 4 records
the SiOH variable-atom driver correction. The transition from version 3
approves evaluating SiOH with its public
`xyz(natoms,3)` shape instead of a legacy oversized array whose storage layout
scrambled the regression geometry. Version 3 records the earlier C2H6F indexing
and uninitialized-state corrections, and version 2 records the CODATA 2022
Hartree-to-kcal/mol migration. Records whose approved values did not change
retain the base revision captured on that same runner, so an unchanged PES is
never compared with another architecture's fixture. Comparisons between two
version-5 revisions remain strict. The
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

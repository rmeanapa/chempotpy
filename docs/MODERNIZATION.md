# CHEMPOTPY modernization

This document records the architecture, source modernization, constants,
precision policy, and catalog-wide validation used by the maintained
[rmeanapa CHEMPOTPY fork](https://github.com/rmeanapa/chempotpy). It keeps
maintainer-level implementation detail separate from the user-focused
[README](../README.md).

Related contracts and evidence:

- [PES contract](PES_CONTRACT.md)
- [Modern Fortran conventions](FORTRAN_STYLE.md)
- [PES import evidence](PES_IMPORTS.md)
- [External-data migration](EXTERNAL_DATA_MIGRATION.md)
- [Contribution procedure](../CONTRIBUTING.md)

About this fork
---------------

This repository is the [rmeanapa CHEMPOTPY fork](https://github.com/rmeanapa/chempotpy)
maintained by Ruben Meana-Paneda. It retains the original scientific authors
and citation listed in the README while developing a self-contained, modern
implementation.

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
  linkage, extension-import, and numerical-regression gates for the complete
  callable catalog, including nine stored reference geometries per surface.
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
Python and core tests, strictly compiles every callable PES interface, checks
formatting, and evaluates every surface at nine deterministic geometries in the
same process. `make test-surfaces` runs only the slower catalog-wide numerical
and architecture regression portion.

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

This command builds all 404 f2py extensions, checks their native-library
linkage, imports every extension in an isolated Python process, and performs
3,636 numerical evaluations. A successful run ends with counters equivalent to:

    pes_extensions_checked: 404
    missing_or_ambiguous_extensions: 0
    external_mkl_lapack_blas_links: 0
    python_extensions_checked: 404
    python_extension_evaluations: 3636
    python_extension_failures: 0

Run the catalog-wide compiler warning gate with:

    make check-fortran-warnings

This compiles all 466 PES and system-support units with `-Wall`, `-Wextra`, and
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

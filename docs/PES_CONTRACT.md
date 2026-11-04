# Shared potential-energy-surface contract

This is the living design and validation record for modernizing CHEMPOTPY and
interoperating with POLYRATE. It replaces ad hoc conversion and wrapper logic;
it does not replace published numerical kernels or refit their coefficients.

Public `pes` and `dpem` routines are compatibility adapters. Repeated
validation, geometry, and unit mechanics belong in the bundled core. A surface
should ultimately own only its evaluator, model metadata, native-unit
declaration, and topology that is scientifically specific to that fit.
Legacy `ENTRY` dispatch is not part of this contract. Initializers, batch
evaluators, and point evaluators are distinct module procedures; any state they
share has one private module owner.

## Module-per-surface architecture

Each migrated surface owns one uniquely named module, for example
`h3_gen_a2_2002_pes`. The module is `private` by default and publishes the same
small generic vocabulary:

```fortran
module h3_gen_a2_2002_pes
    implicit none
    private
    public :: evaluate, metadata
    contains
    subroutine evaluate( coordinates, energy, gradient, coupling, status )
    end subroutine evaluate
end module h3_gen_a2_2002_pes
```

The module owns fit-specific constants and private cached state. Shared
geometry, units, validation, and status definitions remain in the bundled core.
Standalone Fortran programs use the surface module directly. During migration,
the historical external `pes` or `dpem` routine remains as a thin adapter to the
module API so Python extensions and existing callers remain compatible.

`evaluate` and `metadata` are the complete public lifecycle; there is no
`initialize` procedure. Fixed values and fitted coefficients are typed private
module parameters. A surface that historically reads coefficients at runtime
is migrated by generating and committing a private Fortran parameter module
from the authoritative data file. The generated source records the original
data checksum and citation. Installed surfaces therefore require no runtime
data file, mutable initialization flag, or file-unit setup.

Programs that import more than one surface rename the generic procedures on
import, for example `evaluate_a2 => evaluate`. A later runtime registry may
bind these implementations through one abstract interface without changing the
surface modules.

CHEMPOTPY does not implement or depend on POLYRATE's object hierarchy. POLYRATE
owns an adapter type that extends its abstract potential class and delegates to
the CHEMPOTPY module procedures. Standalone programs bypass that object adapter
and call the same CHEMPOTPY module API directly.

## Boundary

The portable low-level boundary uses:

- Cartesian coordinates shaped `(natoms, 3)` in bohr.
- Adiabatic energies shaped `(nstates)` in hartree.
- Adiabatic gradients shaped `(nstates, natoms, 3)` in hartree/bohr.
- Nonadiabatic couplings shaped `(nstates, nstates, natoms, 3)` in 1/bohr.
- Diabatic matrices shaped `(nstates, nstates)` in hartree.
- Diabatic gradients shaped `(nstates, nstates, natoms, 3)` in hartree/bohr.

CHEMPOTPY's established Python API may continue accepting angstrom and returning
eV or eV/angstrom. Those conversions belong in the Python-facing adapter, not
inside every PES kernel. POLYRATE wraps the same atomic-unit boundary in its
Fortran object interface.

## Ownership

- A numerical kernel owns its published fit parameters and any fit-specific
  internal units.
- `chempotpy.core` owns Python contract validation, registry discovery, and
  public-boundary conversions. Its bundled Fortran modules own numerical kinds,
  shared conversions, geometry primitives, linear-algebra support, and status
  values.
- POLYRATE owns its abstract Fortran type and program compatibility wrappers.
- Neither repository copies the other repository's application framework.
- CHEMPOTPY ships its own small symmetric eigensolver behind the historical
  `DSYEV` interface. PES extensions therefore do not require MKL, LAPACK, or
  BLAS; POLYRATE can reuse the same standalone Fortran source when needed.
- A future C-interoperable procedural ABI may connect independently built
  CHEMPOTPY kernels to POLYRATE without exposing Python or Fortran derived
  types across the boundary.

## Bundled Fortran core and standalone surfaces

Standalone means that a compiled PES extension has no separately installed
CHEMPOTPY or scientific-library dependency. It does not require every maintained
PES source to duplicate constants and geometry utilities.

The maintained bundled core currently contains:

- `chempotpy_kinds`: named `real32` and `real64` kinds.
- `chempotpy_units`: canonical conversions for new atomic-unit contracts.
- `chempotpy_legacy_units`: the released Python wrapper conversions, including
  their original default-real rounding behavior.
- `chempotpy_geometry`: Cartesian distance, three-component dot product, angle
  cosine, and symmetric distance-matrix operations.
- `chempotpy_h3_basis`: private symmetric-polynomial index storage and a
  procedural initialization/lookup API for the H3 2002 fit family.
- `chempotpy_legacy_coordinates`: explicit-argument implementations of the
  historical POTLIB Cartesian input and internal-coordinate setup algorithms,
  including the shared legacy triatomic 1-2, 2-3, 1-3 distance adapter.
- `chempotpy_neural_network_io`: shared parsing of neural-network headers.
- `chempotpy_status`: status values for future non-printing modern interfaces.
- `chempotpy_blas`: the small BLAS subset required by bundled legacy kernels.
- `chempotpy_dsyev`: the standalone historical `DSYEV` ABI implementation.

`scripts/build_pes.py` resolves bundled module dependencies recursively from
Fortran `use` statements and compiles the dependency-first source set into each
extension. It does not link a shared CHEMPOTPY core library. A surface therefore
remains independently loadable while the maintained source has one definition
of each core tool.

NumPy f2py does not infer that an imported `real64` kind maps to a C/Python
double. Every isolated build writes an explicit `.f2py_f2cmap` declaring
`real32` as float and `real64` as double. The extension ABI is checked in the
same build path used by the complete catalog.

Published PES fit coefficients, cutoffs, and fit-specific internal conversions
remain in their owning surface. Moving those values into the core would obscure
model ownership and could silently change a published fit.

### Adding a new PES

A new surface should:

1. Implement its numerical kernel in atomic units in a module that is `private`
   by default and explicitly publishes its evaluator.
2. Import only the required bundled facilities, for example
   `use chempotpy_geometry, only: cartesian_distance`.
3. Keep published fit coefficients and model-specific constants beside the
   kernel, with their source and units documented.
4. Provide the thin historical `pes` or `dpem` adapter needed by the Python API;
   new kernel procedures must not print or own application output units.
5. Add independent numerical references and derivative/invariance checks
   appropriate to the model.

No list of core object files is maintained per surface. The build and syntax
gates discover the transitive bundled dependencies from the source imports.
Consequently, adding one core import is sufficient to place the required module
source inside that surface's independently loadable extension.

## Migration sequence

1. Inventory every system, source, representation, derivative capability, data
   dependency, and competing fixed/free-form implementation.
2. Preserve a numerical reference for every callable surface before editing its
   kernel.
3. Introduce one modern adapter with explicit metadata, shapes, units, status,
   and initialization behavior.
4. Compare energy, analytic derivatives, finite differences, atom reordering,
   repeated calls, and two independent instances.
5. Make the modern source authoritative only after its checks pass; until then,
   retain the released source as the build input.
6. Remove the compatibility implementation only after both the Python API and
   POLYRATE adapter use the validated boundary.

## Modern Fortran source style

The normative development checklist is `FORTRAN_STYLE.md`. The summary below
records the principal conventions relevant to the PES contract.

Modernized PES sources follow the same compact conventions used by SIMPLE:

- Free-form Fortran with four spaces per indentation level and no tabs.
- The formatter targets 120 columns for ordinary equations and declarations,
  rather than the fixed-form 72/80-column limit. Generated coefficient or data
  statements may be longer; compilation uses `-ffree-line-length-none`.
- Procedure declarations use spaced argument lists, for example
  `subroutine pes( coordinates, derivative_order, energy, gradient, coupling )`.
- Every procedure has a short `!> @brief` description immediately above it.
- Every scope has `implicit none`; every dummy argument and local variable is
  explicitly declared. New and reviewed interfaces also specify `intent`.
- Variable declarations use the `::` separator and standard kind/length syntax;
  constants use typed `parameter` declarations rather than separate
  `PARAMETER(...)` statements.
- Procedure specifications declare named constants needed by shapes first,
  followed by dummy arguments carrying `intent(...)`, followed by local
  variables. Missing legacy intents are added only after the argument's actual
  read/write behavior is verified; the formatter does not guess direction.
- Modules begin with `private` and explicitly publish only their supported API.
- Imports use `use ..., only: ...`; numerical kernels do not perform application
  orchestration or own hard-coded output units.
- Fortran code is lowercase while character literals and comments retain their
  original capitalization. Consecutive declarations, constants, and simple
  assignments are aligned when doing so improves readability.
- Blank lines are used as single visual separators; consecutive empty lines are
  removed by the formatter.
- Files are opened with `newunit=` and a local unit variable. Existing numbered
  diagnostic writes remain compatibility behavior until a non-printing status
  and logging interface replaces them deliberately.
- Format specifications live directly in their `read` or `write` statement;
  detached labeled `FORMAT` statements are not used.
- Early `return` statements remain available for control flow, but a `return`
  immediately before `end subroutine` or `end function` is omitted.
- Names describe physical quantities and units rather than historical temporary
  abbreviations when a kernel is migrated and numerically revalidated.

`scripts/check_fortran_style.py` checks every `.f90` source rather than relying
on an opt-in marker. The compiler gate adds `-fimplicit-none`, so a declaration
missed by the source transformer is rejected even if the style scanner cannot
recognize the construct. Formatting is checked for idempotence by the same
structural formatter used during migration.

## Repository-wide modernization status

- The package contains 373 free-form `.f90` sources and no fixed-form `.f`
  sources. The registry identifies 363 independently callable PES sources with
  no competing duplicate stems.
- Every procedure scope uses `implicit none`; all variables formerly inferred
  by `implicit double precision(a-h,o-z)`, `implicit real*8`, or default Fortran
  rules now have explicit declarations.
- All block structure uses four spaces per nesting level, tabs are forbidden,
  procedure argument lists are spaced, and every procedure has an immediately
  preceding `!> @brief` description.
- All 2,920 detached labeled formats were either moved into their I/O statement
  or removed when unused. Redundant terminal returns were removed without
  changing early-return control flow.
- All 10,601 legacy `PARAMETER(...)` statements are now typed declarations.
  Legacy star-sized intrinsic types and entity-level character lengths have
  also been removed, and active `OPEN` statements no longer hard-code units.
- The 120 identical `CARTOU` and 120 identical POTLIB `CARTTOR` implementations
  are compatibility wrappers over `chempotpy_legacy_coordinates`; the distinct
  F2H2 derivative variant remains local pending its dedicated derivative test.
- Every module has a `private` default. Its `public` list is limited to entry
  points and names that compiler-resolved consumers import.
- Every module import has an `only:` list. Core dependencies remain discoverable
  by the isolated PES builder.
- Published arithmetic, literal spelling, fit coefficients, and wrapper
  conversion behavior were preserved. Descriptive renaming and algorithmic
  restructuring remain separate changes that require their own numerical
  review.
- `COMMON` removal is intentionally incomplete: 18,485 active declarations
  remain in 159 sources, primarily POTLIB families whose `BLOCK DATA` state and
  aliased layouts must be migrated together. H4O2 and CH3OH now keep their
  mutable PES state private to modules, H3Cl uses explicit saved locals, the
  four H3 2002 fits use the private bundled basis-index module, and 30 dead
  parallel-control blocks have been removed.

The migration utilities are retained under `scripts/`: the fixed-form converter,
compiler-assisted explicit-type inserter, structural formatter, module API
privatizer, and import normalizer can be rerun and are covered by the style and
compilation gates.

## Numerical preservation record

- Before conversion, 205 fixed-form sources were compiled and evaluated. The
  fixture records each source's SHA-256 digest and 656 formatted outputs. The
  repository-wide regression test recompiles all 205 `.f90` replacements and
  compares their outputs with those frozen release references.
- A second catalog fixture covers all 363 callable PES interfaces at two
  deterministic, nonsymmetric geometries. It records the authoritative source
  digest, interface signature, state dimensions, and formatted outputs; run
  `make test-surfaces` to compile and evaluate the complete catalog.
- The earlier Aln fixture independently records 42 source digests and 336
  energies across non-symmetric two-, three-, five-, and eight-atom geometries.
  All modern Aln sources reproduce those values exactly.
- The released BrH2 LTBZSO wrappers declared `engygs` locally while the kernel
  returned it through `/pt1cm/`, producing an undefined diagonal element. Both
  adiabatic and diabatic wrappers now reference the owning shared block; repeated
  captures return stable ground-state energies and the catalog fixture records
  the corrected outputs.
- The CH5 J1 atomic-unit adapter remains checked against the independently
  compiled released behavior. Its small POLYRATE residual reflects the rounded
  conversion literals in the historical CHEMPOTPY wrapper, not the source-form
  migration.

## Validation gates

- Python syntax and contract tests.
- Deterministic source and duplicate inventory.
- Strict compilation of the selected authoritative source for each surface.
- Stored release energy/output references generated before modernization.
- Analytic versus central-finite-difference gradient checks when supported.
- Translation/rotation and atom-order checks where physically applicable.
- Adiabatic/diabatic symmetry and NAC antisymmetry checks for multistate models.
- Complete library build with no surface hidden behind a successful overall
  command when its individual compilation failed.

The current static gate checks all 363 callable surfaces independently; all
363 pass with gfortran. Extension builds use the bundled `DSYEV` implementation
and no MKL, LAPACK, or BLAS link flags.

The optimized macOS extension gate also builds all 363 callable surfaces. A
binary audit finds all 363 exact module artifacts, no MKL/LAPACK/BLAS links,
and an internally defined `DSYEV` symbol in each of the 87 extensions whose
source calls that interface. The only non-system compiled dependency observed
is the normal GNU Fortran runtime; distributable wheels must bundle that runtime
according to their platform's wheel policy.

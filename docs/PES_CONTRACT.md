# Shared potential-energy-surface contract

This is the living design and validation record for modernizing CHEMPOTPY and
interoperating with POLYRATE. It replaces ad hoc conversion and wrapper logic;
it does not replace published numerical kernels or refit their coefficients.

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
- `chempotpy_status`: status values for future non-printing modern interfaces.
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

Modernized PES sources follow the same compact conventions used by SIMPLE:

- Free-form Fortran with four spaces per indentation level and no tabs.
- Procedure declarations use spaced argument lists, for example
  `subroutine pes( coordinates, derivative_order, energy, gradient, coupling )`.
- Every procedure has a short `!> @brief` description immediately above it.
- Every scope has `implicit none`; dummy arguments have explicit `intent`.
- Modules begin with `private` and explicitly publish only their supported API.
- Imports use `use ..., only: ...`; numerical kernels do not perform application
  orchestration or own hard-coded output units.
- Names describe physical quantities and units rather than historical temporary
  abbreviations when a kernel is migrated and numerically revalidated.

Files carrying the `!> @modern` marker are checked by
`scripts/check_fortran_style.py`. A legacy fixed-form kernel receives that marker
only after a reference build and numerical parity check; formatting is never
allowed to silently change a published fit.

The current style inventory finds 2,679 procedure declarations across 370
Fortran/signature files. Forty-eight declarations use the spaced modern form and
46 procedures have `!> @brief` documentation. Twenty-four sources currently opt into
the strict modern style gate. This is the tracked
migration backlog; it is not evidence that the other kernels are numerically
invalid.

## Initial audit findings

- After removing the unvalidated fork conversions and the competing H3Cl
  conversion, the current tree contains 88 system packages and 364 Fortran
  sources. After the current compact Aln batches, 249 are fixed form and 115 are
  free form. They expose 323 adiabatic `pes`
  entry points, 40 diabatic `dpem` entry points, and one retained CH4 support
  source that is not built as a standalone surface.
- The Makefile compiles both `.f` and `.f90` files with the same stem, so an
  unvalidated conversion can silently compete with the released implementation.
- The first fork modernization batch contained free-form ALC and FH2 sources
  that did not pass `gfortran -fsyntax-only`. Those variants were discarded;
  their released fixed-form counterparts are the authoritative starting point.
- The repository tracked 59 compiled extensions, compiler logs, object/module
  files, and editor swap files. They are removed in the current working tree;
  Git history remains the recovery path until the maintainer commits the cleanup.
- The current `make check` searches build logs but does not fail when a surface
  fails compilation and does not validate numerical behavior.

## First compatibility migration: CH5 J1

- The released `CH5_GEN_J1_1987.f` remains the authoritative CHEMPOTPY kernel.
- `chempotpy.CH5.adapters.create_ch5_j1_surface` exposes ordered `HHHHCH`
  coordinates and the single-state energy/gradient through the atomic-unit
  contract.
- The compatibility adapter deliberately uses the released wrapper's rounded
  factors at both sides of the call. A later direct-kernel adapter can remove
  this round trip only after parity is locked down.
- At the POLYRATE regression geometry, the independently compiled released
  CHEMPOTPY source agrees with POLYRATE to `7.5e-9` hartree for energy and
  `1.9e-8` hartree/bohr for the Cartesian gradient. The residual is consistent
  with the rounded conversion constants in the released CHEMPOTPY wrapper.

## First complete style migration: Aln HalP

- `Aln_HalP` is now free-form Fortran with four-space indentation, explicit
  intents, descriptive physical names, spaced procedure arguments, and a brief
  for both its public wrapper and numerical kernel.
- Four deterministic geometries containing two, three, five, and eight atoms
  agree bit-for-bit with the independently compiled released extension.
- The released wrapper's conversion literals were first rounded in default
  real precision. Their exact promoted values remain explicit compatibility
  constants; replacing them with canonical constants is a separate numerical
  change.
- Two- and four-atom release energies are stored in an independent gfortran
  regression test.

## First compact Aln kernel batch

- Ten additional 94--102 line Aln sources are now modern free-form Fortran:
  `ER`, `ER2`, `Gol`, `NP_D`, `PapCEP`, `PapKEP`, `deSPH_M`, `reGol`,
  `reHalP`, and `redeSPH_M`.
- Historical temporaries such as `x`, `v`, `rr`, `aa`, and `t1` were replaced
  with names describing coordinates, distances, fit parameters, and accumulated
  energy. Labeled loops and `go to` statements were replaced with structured
  loops without changing the arithmetic order.
- Several published parameter tables used unsuffixed long decimal literals.
  The released compiler rounded those values to default real before assigning
  them to double precision. The modern sources spell out the exact promoted
  values and document the compatibility constraint.
- All 59 Aln extensions build successfully after the batch. Forty reference
  energies spanning two, three, five, and eight atoms match the released
  extensions bit-for-bit; the references also run as an independent gfortran
  regression test.
- These eleven modern Aln surfaces now import `real64`, released-wrapper
  conversions, and Cartesian distance from the bundled core. A direct comparison
  of 44 extension results still matches the released binaries bit-for-bit.

## Second compact Aln kernel batch

- `PetW`, `rePetW`, `StrM`, `reStrM`, `SutC`, and `reSutC` are now modern
  free-form sources using the bundled kinds, compatibility conversions, and
  Cartesian distance operation.
- Historical labeled loops and repeated coordinate-component arithmetic were
  replaced with structured loops and named physical accumulators while retaining
  the released summation order.
- The reparameterized fits and the original PetW phase multipliers contained
  unsuffixed literals. Their intended source decimals are retained as explicit
  `real32` values and promoted to `real64`, making the compatibility behavior
  visible without storing opaque expanded decimals.
- Twenty-four extension results covering two, three, five, and eight atoms match
  the released binaries bit-for-bit and are stored in the independent gfortran
  regression suite.

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

# CHEMPOTPY modern Fortran style manifest

This manifest is the required source style for maintained CHEMPOTPY Fortran.
It applies to new PES implementations, modernized legacy kernels, bundled core
modules, support routines, and compatibility wrappers. Numerical preservation
and the shared PES contract remain defined in `PES_CONTRACT.md`.

## Source form and layout

- Use free-form `.f90` source. New fixed-form PES files are not accepted.
- Indent every nested block with four spaces. Tabs are forbidden.
- Use lowercase Fortran code. Preserve capitalization in comments and character
  literals when it carries scientific or user-facing meaning.
- Use a 120-column readability target. Generated coefficients and data may be
  longer because builds use `-ffree-line-length-none`.
- Put one space after commas that separate declared entities or procedure
  arguments: `rvp(nt,3), evp(nt), energy`.
- Keep array dimensions and subscripts compact: `coordinates(natoms,3)` and
  `distances(1,:)`.
- Put one space inside procedure argument parentheses:
  `subroutine evaluate( coordinates, energy )`.
- Do not put blank lines inside a procedure's initial specification section.
- Use one blank line, exactly, between top-level subroutines, functions,
  modules, programs, and block-data units.
- Use a single blank line only when it separates meaningful executable sections;
  consecutive blank lines are forbidden.
- Spell compound endings with spaces: `end if`, `end do`, and
  `end subroutine evaluate`; compact forms such as `endif` and `endsubroutine`
  are forbidden.
- Keep the standard PES metadata header as one compact comment block and align
  its values after the longest label. Use the established fields `System`,
  `Functional form`, `Common name`, derivative/body/state counts, `Interface`,
  and `Data file`.

## Procedure documentation and declarations

- Place an immediately preceding `!> @brief` description on every procedure.
- Use `implicit none` in every scope.
- Declare named constants required by array shapes first, then dummy arguments
  carrying `intent(...)`, then local variables.
- New and reviewed procedures must declare `intent(in)`, `intent(out)`, or
  `intent(inout)` for every dummy argument. Do not guess direction during a
  mechanical legacy conversion; verify actual reads and writes first.
- Use the `::` declaration separator and align it within a declaration group.
- Use standard intrinsic syntax such as `real(real64)` and
  `character(len=80)`. Star-sized and entity-sized declarations are forbidden.
- Declare constants with a type and `parameter`, for example
  `integer, parameter :: atom_count = 3`. Separate `PARAMETER(...)` statements
  are forbidden.
- Prefer descriptive physical names and include units in a name when ambiguity
  would otherwise remain. Retain established scientific acronyms where useful.

## Modules, state, and ownership

- Modules are `private` by default. Publish only the procedures and types that
  form the supported API; mutable storage stays private.
- Every `use` statement has an `only:` list.
- Align `only:` across consecutive `use` statements so imported APIs scan as
  one declaration group.
- Do not add `COMMON` or `BLOCK DATA`. Migrate shared mutable state into a
  private owning module with a small procedural API.
- Published fit coefficients and model-specific constants belong beside their
  owning PES. General geometry, unit conversion, numerical kinds, linear
  algebra, status handling, and reusable I/O belong in the bundled core.
- Each compiled PES must remain independently loadable. Core source dependencies
  are compiled into the extension rather than supplied by an external library.

## Geometry, units, and numerical expressions

- Use `chempotpy_geometry` for Cartesian distances, dot products, angle cosines,
  and distance matrices. Do not repeat hand-written Cartesian distance formulas
  in public wrappers.
- Use `chempotpy_units` for new canonical atomic-unit interfaces.
- Use `chempotpy_legacy_units` when a compatibility wrapper must preserve the
  released rounded conversion literals exactly.
- Use `legacy_triatomic_distances_to_bohr` for the historical 1-2, 2-3, 1-3
  triatomic ordering; do not spell out those three distance assignments in a
  wrapper.
- Keep equations together when they fit the 120-column target. Wrap only at a
  safe operator and use free-form `&` continuation.
- Align consecutive simple assignments when this makes related equations easier
  to compare.
- Short, closely related literal initialization assignments may share a line:
  `energy = 0.0d0; gradient = 0.0d0; coupling = 0.0d0`.
- Do not use semicolons to combine branches, calls, long equations, or unrelated
  state changes.

## Public PES adapters

- Treat `pes` and `dpem` as thin compatibility adapters, not as numerical
  kernels.
- Keep the public signature required by the current package while delegating
  common validation, geometry construction, unit conversion, initialization,
  and status handling to the bundled core.
- Keep the surface-specific evaluator responsible for the published fit,
  model-specific internal coordinates, coefficients, and derivative equations.
- New surfaces declare atom count, state count, native units, derivative
  support, and any scientifically meaningful topology as metadata. Do not copy
  a wrapper implementation from another PES.
- `ENTRY` is forbidden. Represent initialization, batch evaluation, and
  single-geometry evaluation as separate, descriptively named module
  procedures with explicit dummy arguments.
- State formerly shared between an initial subroutine and its `ENTRY` points
  belongs in one private owning module. Prefer named constants over mutable
  initialization; when cached state is unavoidable, expose a small explicit
  initialization procedure and keep the storage private.

## I/O and control flow

- Acquire file units with `open(newunit=...)`; do not hard-code units in `OPEN`.
- Put a format specification directly in its `read` or `write` statement.
  Detached labeled `FORMAT` statements are forbidden.
- New numerical kernels return status instead of printing. Existing diagnostic
  writes remain only until their compatibility interface is deliberately
  replaced.
- Keep an early `return` when it expresses control flow. Omit a redundant
  `return` immediately before `end subroutine` or `end function`.

## Required validation

- Run `make check-style` after formatting.
- Run the strict compiler gate for every affected surface.
- Compare affected numerical outputs with the stored reference fixture.
- Run `make test-surfaces` for catalog-wide changes.
- A migration is complete only when formatting is idempotent, compilation
  succeeds, source hashes are refreshed, and numerical records remain equal or
  an intentional correction is documented with an explicit regression value.


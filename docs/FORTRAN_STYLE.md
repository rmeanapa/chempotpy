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
- Write calls without spaces immediately inside their parentheses:
  `call evaluate(coordinates, energy)`.
- Put exactly one blank line before and after every standalone `contains`
  statement. Otherwise, do not put blank lines inside modules, subroutines,
  functions, programs, or block-data units.
- Remove empty comment records and punctuation-only separator bars such as
  `!`, `!-----`, and `!====`. Preserve comments that communicate scientific,
  algorithmic, provenance, unit, or diagram information.
- Do not retain migration markers such as `@modern` or generic ownership
  comments. A PES module brief names the exact surface and chemical system;
  a DPEM module brief identifies the implementation as a diabatic potential
  matrix.
- Use one blank line, exactly, between adjacent subroutines or functions,
  placing it before the documentation for the following procedure. Apply the
  same separation to top-level modules, programs, and block-data units.
- Put the PES metadata and scientific provenance block at the very beginning
  of its source file, before the module brief and implementation.
- Use compact executable-block endings `endif` and `enddo`. Keep source-unit
  endings explicit and spaced, such as `end subroutine evaluate` and
  `end module surface_name`; forms such as `endsubroutine` remain forbidden.
- Space the outer condition and relational operator in conditionals:
  `if( (r-re) .ne. 0.0d0 )then`. Keep arithmetic inside the condition compact
  when that makes the expression easier to read.
- Write the conditional continuation as the single keyword
  `elseif( condition )then`; the split form `else if` is forbidden.
- Keep the standard PES metadata header as one compact comment block and align
  its values after the longest label. Use the established fields `System`,
  `Functional form`, `Common name`, derivative/body/state counts, `Interface`,
  and `Data file`.
- Remove generator disclaimers, repeated star-bar comments, and other purely
  decorative banners. Keep concise scientific provenance in the PES metadata
  or documentation instead of generated boilerplate in numerical kernels.

## Procedure documentation and declarations

- Place an immediately preceding `!> @brief` description on every procedure.
- Use `implicit none` in every scope.
- Declare named constants required by array shapes first, then dummy arguments
  carrying `intent(...)`, then local variables.
- New and reviewed procedures must declare `intent(in)`, `intent(out)`, or
  `intent(inout)` for every dummy argument. Do not guess direction during a
  mechanical legacy conversion; verify actual reads and writes first.
- Use the `::` declaration separator and align it within a declaration group.
- Put every array shape on its typed declaration, for example
  `double precision :: switch(4), dswitch(4,4)`. Standalone `dimension`
  statements are forbidden.
- Use standard intrinsic syntax such as `real(real64)` and
  `character(len=80)`. Star-sized and entity-sized declarations are forbidden.
- Declare constants with a type and `parameter`, for example
  `integer, parameter :: atom_count = 3`. Separate `PARAMETER(...)` statements
  are forbidden.
- Do not redeclare standard intrinsic procedures with `INTRINSIC` statements.
- Do not add `DATA` inside a procedure. Use typed `parameter` constants for
  fixed coefficients and ordinary assignments for explicit initialization.
  Existing procedure `DATA` statements are migration debt until fixed values
  become typed parameters or mutable work values become explicit local
  initialization. `BLOCK DATA` is forbidden and absent from the catalog.
- Prefer descriptive physical names and include units in a name when ambiguity
  would otherwise remain. Retain established scientific acronyms where useful.

## Modules, state, and ownership

- Modules are `private` by default. Publish only the procedures and types that
  form the supported API; mutable storage stays private.
- Each PES implementation belongs to one private owning module that publishes
  `evaluate`. The historical external `pes` or `dpem` routine is only a thin
  compatibility adapter.
- File-local numerical helpers belong under the owning module's `contains`
  section and remain private. Do not leave global helper procedures beside a
  surface merely for legacy linkage. When a distinct support module is still
  required during migration, it is private by default and exposes only the
  narrow operation used by the owner.
- Fixed fit coefficients are private typed `parameter` data. Do not expose an
  initialization procedure merely to assign constants; `evaluate` must be
  callable immediately after import.
- Every `use` statement has an `only:` list.
- Align `only:` across consecutive `use` statements so imported APIs scan as
  one declaration group.
- `COMMON` and `BLOCK DATA` are forbidden and absent from the catalog. Pass
  values explicitly or keep one evaluation's workspace local to its owning
  evaluator; fixed fit data belongs in private typed parameters.
- Do not introduce implicit persistent state with blanket `SAVE` statements.
  During legacy migration, preserve required state until it can be named and
  moved into private module storage; do not delete `SAVE` independently of the
  initialization/evaluation lifecycle that consumes it.
- Published fit coefficients and model-specific constants belong beside their
  owning PES. General geometry, unit conversion, numerical kinds, linear
  algebra, status handling, and reusable I/O belong in the bundled core.
- Each compiled PES must remain independently loadable. Core source dependencies
  are compiled into the extension rather than supplied by an external library.

## Geometry, units, and numerical expressions

- Use `chempotpy_geometry` for Cartesian distances, dot products, angle cosines,
  distance matrices, and atom-major Cartesian flattening. Use
  `flatten_cartesian_coordinates` to transform `x(natoms,3)` into
  `[x1, y1, z1, x2, y2, z2, ...]`; do not repeat indexing loops for this
  layout conversion or hand-written Cartesian distance formulas in public
  wrappers.
- Use `all_pair_distances` for lexicographically ordered unique pairs
  `1-2, 1-3, ..., 2-3, ...`. Use `all_pair_distances_and_jacobian` when both
  distances and `dr/dx` are needed, or `all_pair_distance_jacobian` when a PES
  has already computed the distances. The generic interfaces accept either
  `(natoms,3)` matrices or atom-major `(3*natoms)` vectors.
- Use `ordered_pair_distance_jacobian` when a published fit uses a different
  pair order. Pass the model-owned first-atom and second-atom index arrays;
  never reproduce the Cartesian derivative assignments in the PES. Use the
  component-array overload of `fill_distance_matrix` for legacy algorithms
  that store Cartesian coordinates as separate `x`, `y`, and `z` arrays.
- Use `transform_internal_gradient_to_cartesian` to apply the chain rule
  `dE/dx = (dE/dq)(dq/dx)`. The generic core interface supports scalar-surface,
  multi-state, and matrix-valued gradients. Do not repeat Jacobian-contraction
  loops in public wrappers.
- Use `chempotpy_internal_coordinates` for a supported Cartesian-to-internal
  evaluation workflow. A single-state triatomic distance kernel implements the
  private `evaluate_internal` contract and delegates Cartesian construction,
  canonical unit conversion, chain-rule transformation, and output packing to
  `prepare_legacy_triatomic_coordinates` and
  `finalize_legacy_triatomic_surface`.
- Do not use the generic legacy name `pot` for a newly migrated internal-coordinate
  kernel. Name the private surface equation `evaluate_internal`, or
  `evaluate_diabatic_internal` when it returns a diabatic matrix.
- Use `chempotpy_units` for every public-boundary conversion. Its factors are
  derived from the versioned NIST CODATA manifest.
- Import `wp` from `chempotpy_kinds` for modern numerical code. Declare real
  variables as `real(wp)` and suffix explicit real literals with `_wp`.
- Treat `real64`, `double precision`, and `real(8)` as binary64 compatibility
  declarations. Use them only at a documented external ABI boundary; do not
  introduce them into a new numerical kernel.
- Use `chempotpy_legacy_units` only when a documented parity study proves that
  a released coordinate scale is part of a numerically sensitive published
  fit. Import it with a descriptive local rename such as
  `published_bohr_to_angstrom => legacy_bohr_to_angstrom`.
- Use `triatomic_distances_to_bohr` for the historical 1-2, 2-3, 1-3
  triatomic ordering with canonical units. The separately named legacy routine
  exists only for explicit release reproduction.
- Keep equations together when they fit the 120-column target. Wrap only at a
  safe operator and use free-form `&` continuation.
- Align consecutive simple assignments when this makes related equations easier
  to compare.
- Short, closely related literal initialization assignments may share a line:
  `energy = 0.0d0; gradient = 0.0d0; coupling = 0.0d0`.
- Do not use semicolons to combine branches, calls, long equations, or unrelated
  state changes.

## Public PES adapters

- Give every surface a uniquely named, private-by-default module that publishes
  the generic `evaluate` and `metadata` API defined by the PES contract.
- Do not publish or retain an initialization procedure. Express fixed values
  and fitted coefficients as private typed module parameters.
- Convert authoritative coefficient files into committed Fortran parameter
  modules, recording their source checksum and citation. PES evaluation must
  not require runtime parameter files or file-unit setup.
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
- Use structured `do ... enddo`, `cycle`, `exit`, `if`, and `select case`
  control flow. Labeled `do 100 ...` loops, loop-terminating `100 continue`,
  and `goto` statements are forbidden.

## Required validation

- Run `make check-modernization` for the complete periodic modernization gate.
- The source audit rejects `DATA`, `SAVE`, `COMMON`, `ENTRY`, `BLOCK DATA`,
  `GOTO`, arithmetic `IF`, numeric statement labels, `CONTINUE`, setup-style
  `prepot` procedures, implicit typing, standalone `FORMAT`, Hollerith edit
  descriptors, ADIFOR boilerplate, and fixed-form Fortran sources.
- Run `make check-style` after formatting.
- Run the strict compiler gate for every affected surface.
- Compare affected numerical outputs with the stored reference fixture.
- Run `make test-surfaces` for catalog-wide changes.
- A migration is complete only when formatting is idempotent, compilation
  succeeds, source hashes are refreshed, and numerical records remain equal or
  an intentional correction is documented with an explicit regression value.


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

## Provider transports

The PES request/result contract is independent of process boundaries. Direct
Fortran and Python calls remain the preferred path when CHEMPOTPY is linked or
imported in the same process. The bundled Python socket transport exposes the
same `PotentialEnergySurface` protocol to a persistent external process; PES
modules themselves never contain networking, ports, or connection state.

`SocketPESClient` performs a versioned handshake, receives immutable metadata,
and sends one or more evaluations over an already connected TCP or UNIX-domain
socket. `serve_connection` publishes any Python object implementing the PES
protocol. Messages are length-prefixed UTF-8 JSON. Numerical arrays are
explicitly shaped, big-endian binary64 values encoded with base64, so their
interpretation does not depend on the hosts' byte order. The protocol carries
atomic-unit coordinates and results and preserves CHEMPOTPY's adiabatic,
diabatic, and nonadiabatic result shapes.

The socket layer is a CHEMPOTPY transport, not a replacement for MDI or i-PI.
An MDI adapter can map standard coordinate, energy, and force commands onto the
same provider contract. An i-PI adapter can expose the energy/force subset. A
future POLYRATE provider can therefore select direct calls, this socket
transport, MDI, or i-PI without changing a surface implementation.

A provider process creates a listener around one configured surface:

```python
from chempotpy.core import SocketPESServer

with SocketPESServer.listen_tcp(surface, host="127.0.0.1", port=8021) as server:
    server.serve_forever()
```

The consuming process connects through the same PES interface:

```python
from chempotpy.core import PESRequest, SocketPESClient

with SocketPESClient.connect_tcp("127.0.0.1", 8021) as surface:
    result = surface.evaluate(coordinates_bohr, PESRequest(gradient=True))
```

The first protocol version intentionally serves one configured surface per
connection and serializes evaluations on that connection. Process launch,
authentication, TLS, multi-client scheduling, and HPC resource management are
orchestrator responsibilities. Remote TCP deployments must use a trusted
network, tunnel, or separately authenticated proxy; the bundled protocol does
not claim to provide transport security.

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

- A numerical kernel owns its published coefficients, switching-function
  parameters, neural-network weights, fitted constants, and any fit-specific
  internal units.
- `chempotpy.core` owns Python contract validation, registry discovery, and
  public-boundary conversions. Its bundled Fortran modules own numerical
  kinds, mathematical and physical constants, shared conversions, Cartesian
  layouts/distances/Jacobians, internal-coordinate and gradient-chain-rule
  transformations, finite-difference validation, typed metadata/request/result
  contracts, status values, and small linear-algebra support.
- POLYRATE owns its abstract Fortran type and program compatibility wrappers.
- Neither repository copies the other repository's application framework.
- CHEMPOTPY exposes symmetric diagonalization through
  `chempotpy_linear_algebra`; its private companion implementation preserves
  the historical `DSYEV` ABI without requiring MKL, LAPACK, or BLAS. POLYRATE
  can reuse the same self-contained sources when needed.
- A future C-interoperable procedural ABI may connect independently built
  CHEMPOTPY kernels to POLYRATE without exposing Python or Fortran derived
  types across the boundary.

## Bundled Fortran core and standalone surfaces

Standalone means that a compiled PES extension has no separately installed
CHEMPOTPY or scientific-library dependency. It does not require every maintained
PES source to duplicate constants and geometry utilities.

The maintained bundled core currently contains:

- `chempotpy_kinds`: portable `sp`, `dp`, and `qp` requests plus the single
  compile-time working-precision selector `wp`. The default `wp = dp` preserves
  the released binary64 behavior.
- `chempotpy_atomic_elements`: canonical symbols, atomic-number lookup,
  validation, and standard atomic masses for all 118 elements. PES modules use
  this owner instead of embedding periodic-table maps.
- `chempotpy_constants`: the canonical mathematical constants plus generated
  NIST CODATA 2022 source values and provenance.
- `chempotpy_units`: canonical conversions for all public atomic-unit
  boundaries, derived from the versioned physical constants.
- `chempotpy_legacy_units`: the released Python wrapper conversions, including
  their original default-real rounding behavior.
- `chempotpy_geometry`: Cartesian distance, three-component dot product, angle
  cosine, symmetric distance-matrix operations, and the atom-major
  `(natoms,3)` to `(3*natoms)` Cartesian layout adapter and its inverse. It provides
  dimension-independent lexicographic all-pair distances and distance
  Jacobians for matrix or atom-major vector coordinates, explicitly ordered
  pair Jacobians for published model-specific coordinate orders, and symmetric
  distance matrices from either coordinate matrices or separate component
  arrays. Its generic
  chain-rule interface maps vector, multi-state, and matrix-valued internal
  gradients to Cartesian gradients through a supplied `dq/dx` Jacobian.
- `chempotpy_internal_coordinates`: family-level adapters that connect private
  internal-coordinate kernels to the Cartesian compatibility contract. The
  triatomic adapter builds the 1-2, 2-3, 1-3 distances and Jacobian, applies the
  canonical unit factors, and packs energy, gradient, and coupling outputs.
- `chempotpy_legacy_coordinates`: explicit-argument implementations of the
  historical POTLIB Cartesian input and internal-coordinate setup algorithms,
  including canonical and explicit-release triatomic distance adapters.
- `chempotpy_status`: the `chempotpy_status_t` result and shared status codes
  used by non-printing, non-terminating modern interfaces.
- `chempotpy_contract`: shared capability, metadata, request, and result types
  for modern atomic-unit PES interfaces. Shape, capability, and finite-value
  validation report through `chempotpy_status_t`.
- `chempotpy_finite_difference`: centered Cartesian gradients and Hessians for
  validation and standalone application drivers. It owns displacement and
  mixed-derivative bookkeeping so PES modules do not duplicate it.
- `chempotpy_linear_algebra`: the typed, shape-checked high-level interface for
  real symmetric eigendecomposition. It owns workspace allocation and calls
  the private bundled `DSYEV` companion.
- `chempotpy_interpolation`: shared cubic-spline coefficient construction and
  evaluation, including the tabulated relaxation form required by the ClH2
  TTGI surfaces.
- `chempotpy_blas`: the small BLAS subset required by bundled legacy kernels.
- `chempotpy_dsyev`: the private standalone `DSYEV` ABI implementation used by
  `chempotpy_linear_algebra`.

System- and fit-family support belongs beside the surfaces that own it rather
than in the generic core. `chempotpy/H3/chempotpy_h3_2002_basis.f90` provides
the stateless symmetric-polynomial lookup shared only by H3_GEN_A2_2002,
H3_GEN_A3_2002, H3_GEN_A4_2002, and H3_GEN_CCI_2002. The build and standalone
export registries include that support source automatically for those four
surfaces.

Shared numerical modules and maintained PES kernels declare reals as
`real(wp)` and suffix explicit literals with `_wp`. The bundled BLAS and DSYEV
implementations use the same working kind. Hard-wired `double precision`,
`real(8)`, `real64`, locally defined real-kind aliases, and double-specific
intrinsics are prohibited in maintained numerical source. The package-wide
selector remains `wp = dp` because f2py currently maps `wp` to a C double; a
different package precision requires an explicit conversion layer at that
external boundary.

`scripts/build_pes.py` resolves bundled module dependencies recursively from
Fortran `use` statements and compiles the dependency-first source set into each
extension. It does not link a shared CHEMPOTPY core library. A surface therefore
remains independently loadable while the maintained source has one definition
of each core tool.

Immutable fitted neural-network parameters belong to the owning PES rather
than the general-purpose core. Every affected neural-network surface therefore
compiles a generated, private-by-default parameter module into its standalone
artifact. Released text files remain under `tests/data` as provenance inputs,
while runtime evaluation performs no file I/O or mutable initialization.

The same ownership rule applies to every immutable fitted coefficient or
constant table distributed with a PES. The migration policy, provenance
requirements, and current inventory are maintained in
`docs/EXTERNAL_DATA_MIGRATION.md`.

NumPy f2py does not infer that an imported `real64` kind maps to a C/Python
double. Every isolated build writes an explicit `.f2py_f2cmap` declaring
`real32` as float and `real64` as double. The extension ABI is checked in the
same build path used by the complete catalog.

Published PES fit coefficients, cutoffs, and fit-specific internal conversions
remain in their owning surface. Moving those values into the core would obscure
model ownership and could silently change a published fit.

The public conversion contract is pinned to NIST CODATA 2022. Updating it is a
deliberate, tested source change rather than an automatic dependency update.
The C6H5SH APRP pair and all sixteen C7H8S APRP interfaces retain their released
input-coordinate scale because a catalog parity study found amplified changes
up to about `0.18 eV`; the scale is imported under the descriptive local name
`published_bohr_to_angstrom`. Their output conversions remain canonical.

Internal-coordinate kernels remain private procedures in their owning surface
modules. Compatible families use a descriptive common signature such as
`evaluate_internal`; the owning evaluator calls that equation between the
core preparation and finalization procedures. A universal numerical `pot`
routine is intentionally not part of the contract because single-state,
diabatic, Cartesian, and specialized fitted kernels have different scientific
outputs.

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
- Numerical library code never writes to an external logical unit. The owning
  application may format a returned status message for its selected logger.
- `stop` and `error stop` are forbidden; evaluators return normally and expose
  observable failures through their optional structured status result.
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

- The package contains 482 free-form `.f90` sources and no fixed-form `.f`
  sources: 403 independently callable PES interfaces, 60 system-support
  sources, and 19 shared core modules. No callable surface has a competing
  duplicate stem.
- All 403 callable surfaces have a uniquely named private owning module that
  publishes `evaluate` with an optional `chempotpy_status_t` result and a thin
  historical external `pes` or `dpem` adapter.
  The catalog regression suite enforces exactly one owner module per surface.
  Structural ownership does not by itself complete coefficient-lifecycle
  migration: legacy setup routines are tracked separately and must still be
  replaced by private compile-time parameters where they only assign constants.
- Every owner encapsulates `evaluate` and its file-local helpers as private
  module procedures, leaving only the compatibility adapter external. Some
  families retain separately compiled private support modules, but no surface
  uses `COMMON`, `BLOCK DATA`, `ENTRY`, or a `PREPOT` lifecycle.
- Every procedure scope uses `implicit none`; all variables formerly inferred
  by `implicit double precision(a-h,o-z)`, `implicit real*8`, or default Fortran
  rules now have explicit declarations.
- All block structure uses four spaces per nesting level, tabs are forbidden,
  procedure argument lists are spaced, and every procedure has an immediately
  preceding `!> @brief` description.
- All 2,920 detached labeled formats were either moved into their I/O statement
  or removed when unused. Redundant terminal returns were removed without
  changing early-return control flow.
- Active external diagnostic writes and process-terminating `stop` statements
  have been removed from all library sources; the modernization gate prevents
  either construct from returning.
- All 10,601 legacy `PARAMETER(...)` statements are now typed declarations.
  Legacy star-sized intrinsic types and entity-level character lengths have
  also been removed, and active `OPEN` statements no longer hard-code units.
- ADIFOR disclaimer blocks, generator-provenance boilerplate, and decorative
  star separators have been removed. The 40 formerly ADIFOR-bearing sources no
  longer use explicit `INTRINSIC` statements. Required procedure-local `DATA`
  and saved coefficient caches remain visible where they are part of a kernel's
  current private implementation; none is transported through legacy shared
  storage or a public setup routine.
- Obsolete `CARTOU` input scaling and atom-reordering wrappers have been
  removed. Released surfaces use angstrom input in their public contract, so
  their fixed unit scale of one and disabled zero-filled permutation tables
  collapsed to direct, explicitly shaped Cartesian assignments. POTLIB
  Cartesian-to-internal transformations remain centralized in
  `chempotpy_legacy_coordinates`; the distinct F2H2 JBKKL Jacobi transform
  remains local to its owning surface.
- The 71 wrappers that require atom-major one-dimensional Cartesian vectors use
  `flatten_cartesian_coordinates`; the style gate rejects reintroduction of the
  duplicated `(iatom-1)*3+idir` packing loop.
- The 67 standard wrappers that contract internal gradients with `dq/dx` use
  `transform_internal_gradient_to_cartesian`; its explicit loops preserve the
  released accumulation order and numerical results.
- Every module has a `private` default. Its `public` list is limited to entry
  points and names that compiler-resolved consumers import.
- Every module import has an `only:` list. Core dependencies remain discoverable
  by the isolated PES builder.
- Published arithmetic, literal spelling, fit coefficients, and wrapper
  conversion behavior were preserved. Descriptive renaming and algorithmic
  restructuring remain separate changes that require their own numerical
  review.
- The remaining POTLIB evaluator cohort uses descriptive workspace names for
  Cartesian input, internal and Cartesian gradients, energies, state requests,
  and derivative limits. Fixed identity conversion factors, zero energy
  offsets, disabled atom permutations, unused output-unit flags, and placeholder
  names such as `legacy_d`, `legacy_g`, `legacy_x`, `legacy_p`, and `legacy_u`
  have been removed. Every cohort checkpoint was recompiled and compared with
  the preceding reference values before the next rename was applied.
- Legacy shared storage and setup dispatch have been eliminated catalog-wide.
  The final POTLIB migration localized each former named block inside its
  evaluator, preserved aliased layouts explicitly, and converted block-data
  values into ordinary typed initialization. The last CH4F, CH4Br, F2H2, ArNO,
  and ClH2 kernels now communicate through host association or explicit
  arguments rather than shared storage. The source gate fixes the inventory at
  zero active `COMMON` declarations, `BLOCK DATA` units, `ENTRY` statements,
  `prepot*` definitions, and `prepot*` calls.
- N4 BDTC compiles its released multi-megabyte `pes.dat` table into a generated
  private parameter module. The large source is intentional: it makes Python,
  standalone Fortran, and embedded-library evaluation deterministic and
  independent of installation paths.
- The POTLIB local-context migration made intentional legacy array overreach
  explicit for the C2H6F, CH4OCl, and CH5O2 VBMM kernels. Corrected local bounds
  preserve the released indexing and eliminate dependence on compiler-specific
  out-of-bounds storage. Positional aliases in the OH3 shared layouts are also
  represented by named local entities rather than memory association.

The one-off source transformers used for the catalog conversion were removed
after validation. The maintained `scripts/` surface is limited to build,
formatting, audit, fixture-capture, cleanup, and regression tooling; the style
and compilation gates enforce the resulting source contract directly.

## Numerical preservation record

- The H3 A2, A3, A4, and CCI wrappers previously skipped their required legacy
  setup, leaving the polynomial correction counts and damping settings at zero.
  Their stateless implementations intentionally correct those four outputs and
  match the pre-migration kernels with `prepot` invoked exactly at both catalog
  geometries. The catalog fixture records the corrected initialized behavior;
  the frozen release fixture remains unchanged as evidence of the wrapper bug.
- ClH2 ALTG2 previously called setup before every evaluation even though the
  kernel requires exactly one call. Its second catalog geometry therefore used
  constants converted twice. The stateless result intentionally corrects that
  value and matches the legacy kernel with one setup followed by both
  evaluations exactly.

- Before conversion, 205 fixed-form sources were compiled and evaluated. The
  fixture records each source's SHA-256 digest and 656 formatted outputs. The
  repository-wide regression test recompiles all 205 `.f90` replacements and
  compares their outputs with those frozen release references.
- A second catalog fixture covers all 403 callable PES interfaces at nine
  deterministic geometries. It records the authoritative source
  digest, interface signature, state dimensions, and CODATA 2022 boundary
  outputs; run `make test-surfaces` to compile and evaluate the complete catalog.
- The earlier Aln fixture independently records 42 source digests and 336
  energies across non-symmetric two-, three-, five-, and eight-atom geometries.
  All modern Aln sources reproduce those values exactly.
- The released BrH2 LTBZSO wrappers declared `engygs` locally while the kernel
  returned it through `/pt1cm/`, producing an undefined diagonal element. Both
  adiabatic and diabatic wrappers now reference the owning shared block; repeated
  captures return stable ground-state energies and the catalog fixture records
  the corrected outputs.
- The CH5 J1 atomic-unit adapter remains checked against the independently
  compiled POLYRATE reference while using the shared canonical boundary
  factors.

## Validation gates

- Python syntax and contract tests.
- Deterministic source and duplicate inventory.
- Strict compilation of the selected authoritative source for each surface.
- Stored release migration references plus current CODATA 2022 catalog
  references.
- Analytic versus central-finite-difference gradient checks when supported.
- Translation/rotation and atom-order checks where physically applicable.
- Adiabatic/diabatic symmetry and NAC antisymmetry checks for multistate models.
- Complete library build with no surface hidden behind a successful overall
  command when its individual compilation failed.

The current static gate checks all 403 callable surfaces independently; all
403 pass with gfortran. Extension builds use the bundled `DSYEV` implementation
and no MKL, LAPACK, or BLAS link flags.

The optimized macOS extension gate also builds all 403 callable surfaces. A
binary audit finds all 403 exact module artifacts, no MKL/LAPACK/BLAS links,
and an internally defined `DSYEV` symbol in each extension that imports
`chempotpy_linear_algebra` (currently 78). The only non-system compiled
dependency observed is the normal GNU Fortran runtime; distributable wheels
must bundle that runtime according to their platform's wheel policy.

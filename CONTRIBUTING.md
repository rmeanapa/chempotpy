# Contributing a potential energy surface

Every CHEMPOTPY potential energy surface (PES) is independently usable, but it
must obey the same public contract, core-library ownership rules, and numerical
regression policy as the released catalog. A new surface is not complete when
it merely compiles: its scientific provenance and independently verified
outputs are part of the contribution.

## 1. Generate the source scaffold

Add a surface to an existing system with:

```bash
python scripts/add_pes.py H3 H3_EXAMPLE_2026 \
    --natoms 3 \
    --nstates 1 \
    --entry-point pes \
    --functional-form "permutationally invariant polynomial" \
    --common-name "H3 example" \
    --reference "Authors, Journal volume, pages (year)"
```

Use `--entry-point dpem` for a diabatic matrix. When the chemical system does
not exist yet, also provide its canonical atom order, for example
`--atoms H,H,O`. The command then creates the system's Python geometry-ordering
helpers as well as its Fortran source. Existing files are never overwritten.

The generated source deliberately contains
`CHEMPOTPY_IMPLEMENTATION_REQUIRED`. Replace the scaffold evaluator with the
published model and remove that marker; CI rejects an unfinished scaffold.

## 2. Implement the owning module

The source belongs at `chempotpy/<SYSTEM>/<SURFACE>.f90` and begins with the
standard metadata block. Keep the owning module `private` by default and expose
only its `evaluate` procedure. Preserve the historical Python-facing `pes` or
`dpem` wrapper expected by the build system.

Add the surface name, capabilities, reference, and atom-order information to
the system's `__init__.py` descriptions. For a completely new system, also add
it to the user-facing system list in `chempotpy/__init__.py`; the scaffolded
geometry helpers provide the initial validation and reverse-ordering contract.

Required implementation rules:

- Use `implicit none`, `real(wp)`, explicit declarations, and explicit
  `intent` attributes.
- Use four-space structural indentation, modern structured control flow, and
  `!> @brief` documentation for every procedure.
- Embed published fit coefficients, neural-network weights, switching
  parameters, and other PES-specific immutable data beside their owning PES.
- Do not read parameter files at runtime. `OPEN` and external `READ` operations
  are forbidden in a new PES.
- Do not introduce `COMMON`, `DATA`, mutable `SAVE`, `ENTRY`, `BLOCK DATA`,
  setup/prepotential entry points, `GOTO`, numbered loops, or writes to global
  logical units.
- Return diagnostics through `chempotpy_status_t`; never terminate the caller.
- Keep the public boundary in angstrom/eV units used by CHEMPOTPY and document
  any published-fit compatibility constant that intentionally differs from a
  current core conversion.

## 3. Reuse the core

Reusable numerical operations belong in `chempotpy/core`, not in a PES. Check
the existing owner before adding a helper:

| Need | Core owner |
| --- | --- |
| Working precision | `chempotpy_kinds` |
| Mathematical and physical constants | `chempotpy_constants` |
| Unit conversions | `chempotpy_units` |
| Status values | `chempotpy_status` |
| Cartesian geometry and array reshaping | `chempotpy_geometry` |
| Internal coordinates and chain rules | `chempotpy_internal_coordinates` |
| Legacy internal-coordinate adapters | `chempotpy_legacy_coordinates` |
| Vector and matrix operations | `chempotpy_blas` |
| Symmetric eigensolution | `chempotpy_linear_algebra` |
| Splines | `chempotpy_interpolation` |
| Neural-network activation/evaluation | `chempotpy_neural_network` |
| D3 dispersion | `chempotpy_dispersion` |
| Finite-difference validation | `chempotpy_finite_difference` |
| PES request/result/metadata types | `chempotpy_contract` |

PES-specific coefficients remain local even when their numerical values look
like general constants. Promote an algorithm to the core only when it is
scientifically generic and has a focused core test.

Very large generated parameter arrays may live in a separate private module
beside the PES. In that case, register the support source in `SUPPORT_SOURCES`
in `scripts/build_pes.py`; the standalone exporter and isolated f2py build must
resolve it before the PES source.

## 4. Establish independent numerical evidence

Choose geometries inside the documented domain and compare the modern source
against an independent authority: the original released implementation,
published tabulated values, or a separately maintained reference program.
Record the origin, compiler, units, and any expected precision loss in the pull
request. Do not use the new implementation as its own scientific authority.

After that comparison succeeds, capture CHEMPOTPY's nine deterministic
regression geometries:

```bash
python scripts/capture_library_references.py \
    chempotpy/H3/H3_EXAMPLE_2026.f90 \
    --output build/H3_EXAMPLE_2026-reference.json
```

Review every captured value against the independent evidence, then add the
record and its exact source hash to the catalog fixture:

```bash
python scripts/refresh_library_reference_hashes.py \
    tests/data/all_surface_references.json \
    build/H3_EXAMPLE_2026-reference.json \
    --add-new
```

Also repeat the capture with `--runtime-checks` as a separate diagnostic. Do
not use that second capture as the ordinary fixture: initialization and bounds
instrumentation can change last-bit behavior in legacy numerical kernels.

CI compares a new surface with these committed values. It does not silently
replace them with the values produced during the workflow.

## 5. Validate the surface

Run the focused onboarding, build, extension, and standalone checks:

```bash
python scripts/check_new_pes.py chempotpy/H3/H3_EXAMPLE_2026.f90
python scripts/build_pes.py --surface H3_EXAMPLE_2026
python scripts/check_python_extensions.py \
    --output-root chempotpy \
    --surface H3/H3_EXAMPLE_2026
python scripts/export_standalone.py \
    --surface H3/H3_EXAMPLE_2026 \
    --with-driver \
    --geometries 9 \
    --output build/H3_EXAMPLE_2026-standalone
```

The onboarding checker itself compiles and runs a one-geometry standalone
example. Compile and run the nine-geometry export above as the longer manual
integration test. Then run the repository gates:

```bash
make check
make check-fortran-warnings
make check-numerical-regression
python -m unittest discover -s tests -v
```

Use `make clean` to remove generated extensions and compiler artifacts.

## Pull-request checklist

- [ ] Complete metadata and scientific citation are at the beginning of the source.
- [ ] The module is private by default and exposes only the required evaluator.
- [ ] PES-specific immutable parameters are embedded; no runtime data files remain.
- [ ] Shared algorithms and conversions use their core owners.
- [ ] The Python extension and standalone export compile and run.
- [ ] Independent reference evidence was recorded and reviewed.
- [ ] Nine-geometry fixture values and the source hash were committed.
- [ ] Style, modernization, warnings, unit, and numerical-regression gates pass.

# External PES data migration

Immutable parameters distributed with a potential energy surface are part of
that scientific model. They must be compiled into a private-by-default module,
not opened through a package installation path during evaluation. Historical
`path` arguments remain only on Python-facing `pes` or `dpem` adapters until a
versioned ABI can remove them.

The released files move to a model-specific provenance directory under
`tests/data` after conversion. They remain byte-for-byte provenance inputs,
with SHA-256 hashes recorded in the owning source or generated parameter
module. A migration is complete only after isolated compilation and numerical
comparison with the established surface reference.

Converted families:

- CH4F `CONST-noso` and `CONST-so`: named parameters in each owner module.
- CH4Br `CONST.dat`: named parameters in its owner module.
- F2H2 `bjkkl.par`: generated immutable coefficient module.
- HOBr `o_hbr_a_pp.par`: generated immutable RKHS parameter module.
- COCO `nn_ococ_w20.txt`: generated immutable neural-network module.
- C3H7NO `coeff.dat`: generated immutable PIP coefficient module; the obsolete
  duplicate under `C3H7NO/save` was removed.
- N4_1A_PIPNN, FH2O PIPNN, H3ClO_PIPNN, and both C2N PIPNN models: generated
  private neural-network modules. H3ClO and C2N retain `path` only on their
  historical public adapters; evaluation does not inspect it.
- Both ArNO surfaces: one shared, generated, immutable spline-grid module;
  released grid validation occurs during generation rather than evaluation.
- CH4OF_PIPNN: generated private single-network module with no mutable setup
  lifecycle or runtime model-file access.
- CH2O, CH2O2, CH4OCl, O2H3, C2H2N2O, and H3S: generated immutable
  single-file neural-network modules. H3S preserves all three released
  networks even though its historical evaluator uses network 1.
- C2H2O, HClOH, HO3, K2Rb2, and N2HOCp: generated immutable three-network
  ensemble modules.
- ClH2O: generated immutable three-network and spline-table module. The
  evaluator no longer allocates, initializes, or deallocates model state on
  every call.
- H3Cl: generated immutable ensemble, diatomic spline-grid, and long-range
  coefficient module.
- NO2 2Ap, 4Ap, and 6Ap: separate generated immutable network modules that
  preserve the released analytic derivatives.
- H4O2, H2CO, CH4 T8, and N4 BDTC: generated immutable modules for their PIP,
  exponent, spline, and tabulated model parameters.
- CH3OH, CH4-H2O-H2O, CH5p, and HOCO-Ar already contained their model
  parameters in source. Their duplicate released files now live only under
  `tests/data/external_parameters` as byte-for-byte provenance.

The POTLIB 2023 CH5p archive and the original ChemPotPy package reference
`alpha_1967.dat` and `q.dat` in an optional long-range branch but distribute
neither table. The short-range released PIP remains available and self
contained. Entering the unavailable long-range branch now produces an explicit
diagnostic instead of attempting an installation-dependent file open; no
unverified replacement data have been invented.

Generation commands are deterministic and are checked by
`make check-embedded-parameters`. Runtime file I/O must not be replaced by
mutable module state or one-time initialization; evaluators remain reentrant.
The test suite enforces that no installed Fortran source contains executable
`open`, `read`, or `close` statements.

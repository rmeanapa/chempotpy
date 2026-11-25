# Submit a potential energy surface

You can contribute a PES without using Git or rewriting legacy code. A
complete scientific submission is more important than the programming
language or source format.

## What to prepare

Gather as much of the following material as is available:

- The original source code. Fixed-form Fortran, modern Fortran, Python, and
  archives are all acceptable starting points.
- Every parameter, coefficient, weight, bias, or supporting data file required
  by the model.
- The paper citation or DOI and a link identifying the source of the
  implementation.
- The chemical composition, required atom order, input units, output units,
  number of electronic states, and available derivatives or couplings.
- At least two representative geometries and their expected energies,
  gradients, couplings, or diabatic matrices. Published tables or output from
  the original program are ideal.
- The redistribution status of the submitted files. Do not upload restricted
  or confidential material; describe it without attaching it.

Missing information does not prevent an initial submission. State what is
unknown so that maintainers can help locate it.

## Submit the surface

1. Open the
   [new PES submission form](https://github.com/rmeanapa/chempotpy/issues/new?template=new-pes.yml).
2. Complete the scientific fields. Drag source files and small archives into
   the file field, or provide a stable repository, journal, or archive link.
3. Submit the form. No fork, branch, local build, or pull request is required.

The issue becomes the shared scientific record for the contribution. A
maintainer will check that the implementation is complete and distributable,
confirm units and reference values, modernize and embed the model data, add
catalog-wide regression tests, and open the integration pull request. The
submitter will be asked to review the scientific description and numerical
comparison before the surface is merged.

## Developer route

If you want to implement and test the surface yourself, follow
[CONTRIBUTING.md](../CONTRIBUTING.md). The repository provides
`scripts/add_pes.py` for generating the initial module and
`scripts/check_new_pes.py` for focused validation. A developer contribution
still uses the same submission evidence: provenance, interface details, and
independently verified numerical results.

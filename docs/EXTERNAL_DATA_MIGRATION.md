# External PES data migration

Immutable parameters distributed with a potential energy surface are part of
that scientific model. They must be compiled into a private-by-default module,
not opened through a package installation path during evaluation. Historical
`path` arguments remain only on Python-facing `pes` or `dpem` adapters until a
versioned ABI can remove them.

The released files move to `tests/data/external_parameters` after conversion.
They remain byte-for-byte provenance inputs, with SHA-256 hashes recorded in
the owning source or generated parameter module. A migration is complete only
after isolated compilation and numerical comparison with the established
surface reference.

Converted families:

- CH4F `CONST-noso` and `CONST-so`: named parameters in each owner module.
- CH4Br `CONST.dat`: named parameters in its owner module.
- F2H2 `bjkkl.par`: generated immutable coefficient module.
- HOBr `o_hbr_a_pp.par`: generated immutable RKHS parameter module.
- COCO `nn_ococ_w20.txt`: generated immutable neural-network module.
- C3H7NO `coeff.dat`: generated immutable PIP coefficient module; the obsolete
  duplicate under `C3H7NO/save` was removed.
- N4_1A_PIPNN and FH2O PIPNN models: generated neural-network modules.

Remaining loaders fall into four implementation groups:

- neural-network weights and biases;
- fitted coefficient and spline tables;
- legacy formatted constant files;
- the large N4 BDTC dataset.

Generation commands are deterministic and are checked by
`make check-embedded-parameters`. Runtime file I/O must not be replaced by
mutable module state or one-time initialization; evaluators remain reentrant.

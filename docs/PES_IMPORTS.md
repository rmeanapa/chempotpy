# PES import record

This living note records the contract, implementation status, independent
evidence, and validation criteria for surfaces imported after the initial
CHEMPOTPY modernization.

## Contract

Imported surfaces use the public CHEMPOTPY angstrom/eV boundary, own their
published fit parameters, and reuse the core for kinds, constants, unit
conversions, geometry operations, status reporting, and shared numerical
algorithms. Released parameter files are converted into immutable Fortran
modules; installed evaluators perform no runtime file I/O. Fitting programs,
precompiled objects, and unrelated utilities distributed with a surface are
not part of its runtime kernel.

Every import must have a reference calculation from the released evaluator,
nine catalog regression geometries, strict compilation, an isolated Python
extension check, and a standalone-export check.

## Current import set

| Surface | Status | Main remaining work |
| --- | --- | --- |
| H3O2m_sGDML_2024 | Implemented and independently compared | None; all 50 released checkpoints, catalog, and standalone gates pass |
| H5O2p_PIP_2005 | Implemented and independently compared | None; released energy/gradient, catalog, and standalone gates pass |
| H7p_PIP | Implemented and independently compared | None; catalog and standalone gates pass |
| HNi_EAM5_1995 Ni(111) adapter | Implemented and independently compared | None; canonical fixed-slab catalog and standalone gates pass |
| HNi_EAM6_1996 Ni(111) adapter | Implemented and independently compared | None; canonical fixed-slab catalog and standalone gates pass |
| HNi_EAM7_1998 Ni(111) adapter | Implemented with documented release repair | None; canonical fixed-slab catalog and standalone gates pass |
| CH4Ni_EAM6_J1_1996 Ni(111) adapter | Implemented and independently compared | None; all six released checkpoints, catalog, extension, and standalone gates pass |
| H3m_AS_2004 | Implemented with documented release repairs | None; catalog and standalone gates pass |
| HeH2p_DIM_1976 | Implemented and independently compared | None; catalog and standalone gates pass |
| HeH2p_MTJS_1987 | Implemented and independently compared | None; catalog and standalone gates pass |
| LiH2p_Ground_2003 | Implemented and independently compared | None; catalog and standalone gates pass |
| LiH2p_Excited_2003 | Implemented and independently compared | None; catalog and standalone gates pass |
| CH2OH two/three-state | Implemented and independently compared | None; catalog and standalone gates pass |
| Malonaldehyde MHT 2014 | Implemented and independently compared | None; catalog and standalone gates pass |
| Malonaldehyde Bowman/Braams 2008 | Blocked: released archive omits the PES library source | Locate the source behind `libpes.a`; coefficient files alone are insufficient |
| C2H4O4 formic-acid dimer | Implemented and independently compared | None; catalog and standalone gates pass |
| CO2 five-state RKHS | Implemented and independently compared | None; catalog and standalone gates pass |
| O+C2H2 ground-triplet PIP-NN | Existing implementation independently audited | None; released 14-geometry example agrees |
| BrHCl BLM extended LEPS | Implemented and independently compared | None; catalog and standalone gates pass |
| C7H8S Ia, II, S1.8, S2.0, and S2.5 APRP pairs | Implemented and independently compared | None; all ten interfaces pass |
| CH2O2 local PIP-NN | Implemented and independently compared | None; embedded three-network ensemble passes |
| FH2 5SECGM4 | Implemented from the released evaluator | None; catalog and standalone gates pass |
| H3 LSTH | Implemented from the released evaluator | None; embedded spline and catalog gates pass |
| HCl-H2O Bowman-Mancini PIP 2013 | Implemented and independently compared | None; immutable 1,603-term fit passes |
| N4 MEG2 2020 | Implemented and independently compared | None; catalog and standalone gates pass |
| NH3 two-state Surfgen 2012 | Implemented and independently compared | None; adiabatic and diabatic interfaces pass |

The fully variable Ni surface/cluster models still require a separate
variable-atom/environment contract and are not disguised as fixed-size
molecular PESs. The first Ni import is therefore an explicit adapter for the
released canonical 698-atom fixed Ni(111) configuration with one active H
coordinate. EAM5, EAM6, and EAM7 now share that explicit adapter contract.
The methane/Ni release uses a different explicit fixed-environment contract:
four active H atoms and one active C atom, in that order, above the released
397-atom Ni(111) secondary zone. Its implementation is tracked separately as
`CH4Ni_EAM6_J1_1996`. Movable-Ni configurations and alternate crystal faces
remain in the variable-environment batch.

POTLIB's `C2O2` archive was also audited. It is the same three-network model
already exposed as `COCO_PIPNN`, so a second callable name and a duplicate
parameter module were deliberately not added.

## H3O2m sGDML and POTLIB alias audit

The requested POTLIB `N4_1A_MB-PIP-MEG3` and
`N2O2_3A_MB-PIP-MEG2` releases were already present as `N4_1A_ZV` and
`N2O2_3A_ZV`, respectively. Their source metadata, fit coefficients, updated
pair potentials, and distributed two-geometry outputs identify them as those
exact releases, so duplicate callable aliases were not added.

`H3O2m_sGDML_2024` embeds the authors' 3,000-point binary64 model and evaluates
all 12 oxygen/hydrogen permutations through the shared Matérn-5/2 sGDML core.
The 36,000 expanded descriptor records in the released text model are
reconstructed losslessly from the 3,000 base records and the published
permutation map. No model file is opened at runtime and no setup call or saved
mutable state is required.

Independent evidence was obtained by rebuilding commit `6dade4a` of the
authors' H3O2M Fortran evaluator with gfortran and comparing all 50 distributed
geometries in bohr. After applying the canonical CODATA 2022 public boundary,
the largest energy difference is `2.25e-9` hartree and the largest Cartesian
gradient difference is `3.82e-9` hartree/bohr. The nine catalog geometries use
two released configurations plus their controlled interpolations and
distortions. The test-only NPZ archive is retained as generator provenance;
the installed Fortran module contains the immutable parameters directly.

## 2026 POTLIB completion batch

This batch adds 18 callable interfaces and raises the catalog from 384 to 402
surfaces. Every new runtime is self-contained: released text parameters are
kept as provenance fixtures and generated into immutable Fortran modules, but
installed evaluators perform no file I/O and require no setup call.

The released CH2O2 local PIP-NN minimum evaluates to `-1.8556e-4` eV in both
implementations. The HCl-H2O 2013 minimum is
`-536.73159603777935` hartree in the released program and
`-536.7315960376366` hartree after the CODATA 2022 public-coordinate
conversion. The NH3 reference Hamiltonian has diagonal elements
`0.015036404252298` and `0.214925976476203` hartree in the rebuilt Surfgen
2.8.5 evaluator; the modern adiabatic interface returns the corresponding
`0.40916189` and `5.84843376` eV values.

The five added C7H8S parameterizations were rebuilt from their authoritative
POTLIB source at the first catalog geometry. Their released adiabatic energies
in hartree are, respectively: Ia `(5.6575397244, 6.1114383704,
6.2563635693)`, II `(5.6575397215, 6.1114266102, 6.2547584860)`, S1.8
`(5.6575397259, 6.1114391251, 6.2572688934)`, S2.0 `(5.6575397278,
6.1114401262, 6.2584872554)`, and S2.5 `(5.6575397326, 6.1114425576,
6.2615332311)`. The modern interfaces reproduce them after applying the
release's independently rounded coordinate boundary and the canonical CODATA
energy conversion. Each adiabatic/DPEM pair uses the same owned evaluator and
is checked at nine geometries.

The released N4 MEG2 evaluator returns `457.2551111154461` kcal/mol at the
first catalog geometry; the modern result is `457.2551130032926` kcal/mol.
The `1.9e-6` kcal/mol difference follows from the canonical public conversion
rather than a change to the 276 fit coefficients or the pair and D3 terms.
BrHCl BLM, FH2 5SECGM4, and H3 LSTH retain their released functional forms and
coefficient spelling while replacing COMMON/BLOCK-DATA/ENTRY lifecycles with
ordinary module ownership. LSTH's 87-point H2 spline is compiled directly into
the surface. All three pass strict compilation, nine-geometry numerical
regression, and standalone export.

## CH4Ni_EAM6_J1_1996 import contract

The public geometry contains `H,H,H,H,C` in angstrom. The 397 released Ni(111)
coordinates, all fit coefficients, and the density-to-pseudodistance maps are
immutable source parameters. Evaluation performs no setup call or runtime file
I/O. The energy is the sum of four generalized two-state DIM roots, the fixed
slab's density-dependent nickel embedding contribution, and the modified J1
methane bending term. First derivatives are evaluated from the same stateless
energy kernel using the shared centered finite-difference utility until the
released analytic chain rule has been independently reconstructed and audited.

The release has one numerically significant initialization defect. Carbon
declares ten orbital coefficients and its density loop visits all ten, but the
setup loop initializes only `NTS+NTP=4` radial prefactors. The six remaining
COMMON-block values happen to be zero in the distributed executable. The
modern evaluator makes that released behavior deterministic by retaining the
four initialized carbon radial terms explicitly; enabling the other published
coefficients would define a different surface and fails the released energy
checkpoints.

Independent component references were captured from a rebuilt authoritative
evaluator at all six distributed geometries. In hartree, the EAM/DIM and
modified-J1 contributions are respectively: `(-1.7594223550e-5,
1.7646716880e-5)`, `(1.9421758751e-2, 6.5940727505e-3)`,
`(5.1886733836e-2, 6.2008948652e-3)`, `(-1.1013728212e-2,
1.1774060522e-2)`, `(8.6729064300e-3, 1.3627036775e-2)`, and
`(3.2840793670e-3, 6.6481739535e-3)`. With the separated-slab zero aligned to
the release, the modern totals differ by at most `0.0017` kcal/mol over all six
geometries. This envelope comes from using the shared full-precision `PI` and
NIST CODATA 2022 boundary conversions instead of the rounded release constants.
At the asymmetric `geom3` checkpoint, the centered numerical gradient agrees
with 14 of the 15 printed release components within `1.3e-6` hartree/bohr; the
carbon z component differs by `1.3e-5` hartree/bohr from the release's analytic
chain rule. The numerical derivative is retained because it differentiates the
audited energy kernel directly. The six component and total-energy comparisons,
nine catalog geometries, strict compilation, Python extension, and standalone
export are the acceptance gates for the import.

The released Bowman/Braams malonaldehyde archive is different: its two source files are only
a command-line shell around a precompiled `libpes.a`. The module that owns
`pes0_init`, `pes1_init`, and `pes_h4c3o2s_pot` is supplied only as compiler-
specific module and object files. Importing the coefficient files without that
missing functional-form source would not reproduce the published surface, so
this entry remains blocked rather than being reverse engineered from binaries.

## Malonaldehyde MHT 2014 independent evidence

The Mizukami-Habershon-Tew release is a complete source evaluator for the
enol-malonaldehyde chelate and non-chelate regions, fitted to 14,757
fc-CCSD(T)(F12*)/def2-TZVP energies. Its 9 Morse, 80 one-dimensional, 1,728
two-dimensional, and 1,741 three-dimensional Gaussian terms are stored as
immutable parameters. The installed evaluator uses the shared all-pair
distance/Jacobian and internal-to-Cartesian gradient utilities; it performs no
runtime file I/O or initialization.

The authors' source was compiled in binary64-default mode and evaluated at the
published minimum plus seven independent distortions. Modern and released
energies agree exactly at all eight geometries, while the largest difference
among the 216 analytic Cartesian-gradient components is `1.2e-16`
hartree/bohr. A normal legacy build rounds the unsuffixed fit literals through
default real and differs by as much as `1.3e-7` hartree; the working-precision
module deliberately preserves the published decimal coefficients instead.

## HNi_EAM5_1995 independent evidence

POTLIB distributes EAM5 as a variable-atom program plus three environment
files. The fixed-size CHEMPOTPY adapter embeds the released `potfile1`
Ni(111) environment: one active H atom and 698 fixed Ni atoms. The evaluator
reconstructs the electron densities, piecewise nickel and hydrogen embedding
energies, screened pair interaction, and analytic H gradient without COMMON
blocks, mutable initialization, global units, or runtime data files. The
historical fixed Ni-Ni pair energy is intentionally omitted, as in the
released secondary-zone algorithm.

At the two released Ni(111) test geometries, the modern energies differ from
the ordinary legacy build by about `5e-4` kcal/mol after retaining the zero
of separated H plus slab. The remaining change comes from using binary64 fit
coefficients and the shared NIST CODATA 2022 bohr conversion rather than the
release's rounded default-real literals and `0.52917706` factor. All three
analytic H-gradient components retain the released values to the printed
test precision; an independent centered-difference check agrees within
`1.2e-8` eV/angstrom.

## HNi_EAM6_1996 and HNi_EAM7_1998 independent evidence

The POTLIB EAM6 and EAM7 releases use the same canonical 698-atom `potfile1`
Ni(111) environment and nickel fit parameters. Their modern adapters therefore
share one immutable environment module and one stateless EAM runtime while
retaining the version-specific hydrogen embedding and effective-charge
branches. Only the active H coordinate is public; fixed-fixed pair energies
are omitted exactly as in the released secondary-zone algorithm. The shared
runtime contains no setup call, mutable saved state, global unit, or external
file dependency.

At the two released `geom1` and `geom2` checkpoints, modern EAM6 energies
differ from the published test output by at most `1.3e-4` kcal/mol and EAM7 by
at most `1.3e-4` kcal/mol. The small differences follow from binary64 fit
literals, the shared full-precision `PI`, and NIST CODATA 2022 boundary
conversions in place of rounded release constants. All active-H gradient
components agree with the printed legacy values. Independent centered
differences agree with the analytic gradients within `4.9e-8` eV/angstrom.
The separated H-plus-slab limit is retained as exactly zero.

The released EAM7 source initializes its auxiliary embedding values only when
the local density exceeds `0.5`, but later consumes them unconditionally. The
modern branch explicitly assigns zero below that threshold, which is the
continuous inactive-correction value and removes the compiler-dependent
uninitialized result without changing the published checkpoints.

## H3m_AS_2004 independent evidence

The authoritative Panda-Sathyamurthy source computes two isotropic
polarizability branches as `alpha01` and `alpha02` but then selects undefined
variables named `alpha1` and `alpha2`. The latter name is subsequently reused
for the anisotropic polarizability. Consequently, an asymptotic calculation
compiled with signaling-NaN initialization returns NaN. The release also
leaves `etot` undefined when its switching argument is exactly one. The modern
evaluator repairs the branch to select `alpha01` below 4 bohr and `alpha02`
above it—the two expressions meet there—and gives the switch boundary its
continuous limiting value.

The corrected fixed-form source was rebuilt with binary64 default literals and
compared at eight compact, switching-region, and asymptotic geometries. Away
from the exact piecewise boundary, the largest modern/legacy difference is
`2.4e-13` hartree. At a reconstructed nominal 4-bohr distance, the difference
is `9.3e-9` hartree because roundoff places the Cartesian reconstruction on the
other side of the published branch. The modern evaluator also identifies the
shortest H-H distance as the H2 asymptotic fragment; permutations of the three
atoms agree within `1.8e-14` hartree. The release's 10 pair coefficients and 43
three-body terms are immutable module parameters, with no `DATA`, implicit
typing, initialization, external files, or global state.

## C2H4O4 formic-acid dimer independent evidence

The released Qu-Bowman formic-acid dimer package contains 3,764 generated
Fortran source units, a precompiled invariant library, and external PES and
dipole coefficient files. The energy evaluator's transitive closure was first
reduced to 38 modules and then algebraically specialized for its fixed
four-H/four-O/two-C partition and degree-four fit. Unreachable degree-five and
higher generated terms and the independent dipole branch were excluded. The
installed runtime contains a 1,784-term invariant basis and the exact 1,784
released PES coefficients; it has no initialization, runtime file I/O, or
compiler-specific binary library.

The specialized basis was compared term-for-term with the full released
invariant engine at all 13,475 geometries in `pts.dat`. The largest energy
difference was exactly zero in binary64. The public evaluator uses the shared
NIST CODATA 2022 angstrom/bohr and hartree/eV conversions. Strict source and
standalone Fortran 2018 builds pass, and the standalone example reports energy,
finite-difference gradient, and finite-difference Hessian without external
data files.

## CO2 five-state RKHS independent evidence

The released `CO2-PESs` evaluators and RKHS toolkit were rebuilt with GNU
Fortran for the `1A'`, `1A"`, `(2)1A'`, `3A'`, and `3A"` states. At the
release's asymmetric `rOO=2.4`, `rOC=2.3`, and `rCO=4.5` bohr geometry, the
largest difference between the immutable CHEMPOTPY energies and the released
fast-cache evaluator is `2.0e-7` hartree; the largest distance-gradient
difference is `8.2e-7` hartree/bohr. The difference is caused by replacing the
large pre-summed binary cache with a compensated direct contraction of the
same released double-precision RKHS coefficients. This reduces approximately
100 MB of compiler-specific runtime cache files to about 2.2 MB of portable
Fortran parameters and eliminates all runtime file I/O and initialization.

The shared core owns the Taylor-spline and reciprocal-power kernels. A CO2
runtime module owns the Jacobi transformations, channel blending, asymptotic
correction, switching weights, and Cartesian chain rule. Each of the five PES
modules owns only its fitted coefficients, channel offset, and switching
scales. All five pass strict compilation and the nine-geometry catalog gate;
the generated standalone source also compiles under strict Fortran 2018.

## O + C2H2 PIP-NN independent evidence

The existing `C2H2O_PIPNN` implementation was compared with the authoritative
POTLIB `C2H2O.zip` release at all 14 geometries distributed in `test.xyz`.
The largest energy difference is `4.3e-14` eV. The released evaluator reads six
weight and bias files during `pes_init`; the modern evaluator embeds the same
three-network ensemble and requires neither initialization nor a data-path
argument. Its large parameter constructors are split into standard-conforming
compile-time chunks, so the standalone export passes the Fortran continuation
limit without changing the coefficient ordering.

## H5O2p_PIP_2005 independent evidence

The authoritative QM-22 H5O2+ release supplies the 8,001 fitted coefficients,
an expanded PIP evaluator, reverse analytic derivatives, and a reference
geometry. At that geometry the release reports `-153.01000221` hartree; the
modern evaluator gives `-153.01000221171148` hartree with the NIST CODATA 2022
angstrom/bohr conversion. Every one of the 21 analytic gradient components
agrees with the release to its reported eight decimal places, and the modern
gradient also agrees with central finite differences.

The coefficient file is embedded verbatim (SHA-256
`a4941bb7cdeadbea6bde5c385f2a6928c4e339c1ce7fa9b444f1390e752d6ee3`).
The released generated source has SHA-256
`5741e49f6eff7e323a82de9f413d035c7414656d8c46ee3ec88299b10d051236`
and expands forward and reverse differentiation into roughly 221,000
Fortran lines, which made a GNU Fortran debug compile exceed 2.7 GB before
completion. The portable implementation encodes those same ordered expressions
as an immutable binary operation graph and performs reverse accumulation over
that graph. It compiles with strict Fortran 2018, needs no initialization or
external data files, and retains analytic first derivatives.

## H7p_PIP independent evidence

The authoritative POTLIB H7p archive was compiled with GNU Fortran on macOS.
At the released `geom.xyz` geometry, its evaluator produced
`-3.2966299861803634` hartree using the historical
`0.5291772083` angstrom-per-bohr boundary factor. The modern evaluator produced
`-3.2966299853096217` hartree using CHEMPOTPY's NIST CODATA 2022 conversion.
The difference is `8.71e-10` hartree and is entirely consistent with the
boundary conversion update. The modern runtime retains the released 739 PIP
coefficients and degree-seven invariant construction.

## HeH2p_DIM_1976 independent evidence

The released Whitton-Kurtz evaluator from the IIT Kanpur analytic-potential
archive was compiled with GNU Fortran on macOS. Nine geometries were evaluated
with distances supplied in the release's bohr interface and converted back to
the fit's native eV units. The modern evaluator reproduced all nine energies
bit for bit. Its analytic Cartesian gradient additionally agrees with central
finite differences to `3.7e-9` eV/angstrom at an asymmetric geometry.

The release declares its distance argument as `rs(1)` while reading three
elements; GNU Fortran diagnoses those accesses as out of bounds. The modern
interface corrects the extent, removes implicit typing and `DATA`, and applies
the shared Cartesian-distance Jacobian. The embedded numerical coefficients
and native eV/angstrom algebra are unchanged.

## HeH2p_MTJS_1987 independent evidence

The released Joseph-Sathyamurthy evaluator was compiled with GNU Fortran's DEC
compatibility mode and evaluated at nine geometries. Standard legacy mode is
insufficient because the source initializes an implicitly integer array with
character data. It also performs first-call writes of uninitialized title data
to unit 9 and rounds both coordinates and unsuffixed `DATA` constants through
single precision.

The modern evaluator removes the output, initialization flag, dead derivative
branch, implicit typing, and forced single-precision conversions. Its largest
energy change over the nine-point comparison is `9.5e-7` eV. This is the
expected effect of retaining the published decimal coefficients in `wp`; the
functional form and coefficient values are otherwise unchanged.

## LiH2p 2003 independent evidence

The released ground- and first-excited-state evaluators were compiled with GNU
Fortran using binary64 default literals and evaluated at nine independent
triatomic geometries. The modern ground state agrees within `1.8e-14` hartree;
the excited state agrees within `2.2e-11` hartree. The latter envelope reflects
Cartesian reconstruction of the comparison geometries and cancellation among
the high-order fit terms.

The release required separate `GROUND_I` and `EXCITED_I` calls to populate
three COMMON blocks and printed every generated exponent. The modern modules
embed the exact exponent enumerations (140 ground-state terms and 227 excited-
state terms), published coefficients, pair potentials, and switched long-range
models. Initialization calls, COMMON storage, implicit typing, `DATA`, `STOP`,
and writes to global units have all been removed.

## CH2OH Surfgen independent evidence

The released two- and three-state Surfgen programs were rebuilt with GNU
Fortran, their LAPACK dependency, and the required default-eight-byte-integer
mode. At the supplied reference geometry, the released two-state energies are
`-0.006654787305166` and `0.142319682531917` hartree; the released three-state
energies are `-0.006584209731028`, `0.142001316298710`, and
`0.176830389007174` hartree. The immutable evaluators agree to approximately
`1e-14` hartree.

The generic runtime originally generated a 16,026-function symmetry basis and
read `coord.in`, `irrep.in`, `surfgen.in`, and `hd.data` during initialization.
The import contracts that result into one shared 49-coordinate definition,
normalized monomial recurrence arrays, and separate surface-owned coefficient
modules containing 30,180 and 60,360 nonzero contractions. Runtime file I/O,
mutable initialization, diagnostic output, and OpenMP global-state changes are
absent. At the reference geometry, centered Cartesian derivatives reproduce
the released analytic adiabatic and diabatic derivatives within
`1.4e-9` atomic units.

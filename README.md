# ChemPotPy

**CHEMical library of POTential energy surfaces in PYthon**

[![All PES numerical regression](https://github.com/rmeanapa/chempotpy/actions/workflows/all-pes-numerical-regression.yml/badge.svg?branch=main)](https://github.com/rmeanapa/chempotpy/actions/workflows/all-pes-numerical-regression.yml?query=branch%3Amain)
[![Latest release](https://img.shields.io/github/v/release/rmeanapa/chempotpy?display_name=tag)](https://github.com/rmeanapa/chempotpy/releases/latest)
[![License](https://img.shields.io/github/license/rmeanapa/chempotpy)](LICENSE)
[![Submit a PES](https://img.shields.io/badge/Submit-a%20new%20PES-00a6a6)](https://github.com/rmeanapa/chempotpy/issues/new?template=new-pes.yml)

ChemPotPy provides callable analytic potential energy surfaces (PESs),
multistate surface sets, diabatic potential energy matrices, gradients, and
couplings for molecular simulations. Python calls use Cartesian coordinates
in angstrom, energies in electronvolt, gradients in eV/angstrom, and
nonadiabatic couplings in 1/angstrom.

The original library was created by Yinan Shu, Zoltan Varga, Dayou Zhang, and
Donald G. Truhlar at the University of Minnesota. This maintained fork is led
by Ruben Meana-Paneda. Its modernization history and developer-facing design
are documented in [docs/MODERNIZATION.md](docs/MODERNIZATION.md).

## Potential energy surface catalog

The catalog below groups every callable PES and DPEM entry by chemical
composition. Each row gives the exact surface name used by the Python and
standalone-export interfaces and the number of represented electronic
surfaces.

<!-- BEGIN GENERATED PES CATALOG -->

The catalog currently exposes **415 callable surfaces** across **109 chemical systems**. The final column gives the number of electronic surfaces represented by each callable PES or DPEM entry.

<table>
<thead>
<tr><th>System composition</th><th>Surface name</th><th>Surface description</th><th>Electronic surfaces</th></tr>
</thead>
<tbody>
<tr><td rowspan="12" valign="top"><strong>ALC</strong></td><td><code>ALC5m0_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC6m0_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7m0_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7m10_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7m1_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7m4_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7m5_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7m7_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7p2_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7p3_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC7p5_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ALC_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>Al<sub>m</sub>H<sub>n</sub></strong></td><td><code>AlmHn_VBO1</code></td><td>Valence-bond-order potential</td><td align="center">1</td></tr>
<tr><td><code>AlmHn_VBO2</code></td><td>Valence-bond-order potential</td><td align="center">1</td></tr>
<tr><td rowspan="59" valign="top"><strong>Al<sub>n</sub></strong></td><td><code>Aln_BetH</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_CleR</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_CoxJM</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_deSPH_LJAT</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_deSPH_M</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2AT</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2BA</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2CN</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2CoxJM</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2EACN</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2EAT</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2EBA</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2ECN</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2ES</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2ESCNa</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2ESCNm</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2GEA</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2MeiD</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2MisFMP</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2S</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2SCNa</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ER2SCNm</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ErkIV</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ErkV</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_ErkVIII</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_Gol</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_HalP</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_Jac</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_MeiD</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_MisFMP</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_NP_A</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_NP_Ad</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_NP_B</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_NP_Bd</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_NP_C</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_NP_D</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_NP_E</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_PapCEP</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_PapKEP</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_PeaTHT</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_PetW</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reCoxJM</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_redeSPH_M</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reErkIV</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reErkVIII</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reGEA</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reGol</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reHalP</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reJac</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reLJAT</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reMeiD</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reMisFMP</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_rePetW</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reStrM</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_reSutC</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_StrM</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>Aln_SutC</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>ArNO</strong></td><td><code>ArNO_Ap_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>ArNO_App_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="7" valign="top"><strong>BrH<sub>2</sub></strong></td><td><code>BrH2_LEPS_C_1982</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrH2_LEPS_LTBZ1_1995</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrH2_LEPS_LTBZ2_1995</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>BrH2_LEPS_LTBZ2_1995_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>BrH2_LEPS_LTBZSO_1995</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>BrH2_LEPS_LTBZSO_1995_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>BrH2_LEPS_W_1959</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="4" valign="top"><strong>BrHCl</strong></td><td><code>BrHCl_LEPS1</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrHCl_LEPS2</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrHCl_LEPS3</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrHCl_LEPS_BLM_1983</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>BrNH<sub>3</sub></strong></td><td><code>BrNH3_VBMM_2026</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>2</sub>N<sub>2</sub>O</strong></td><td><code>C2H2N2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>2</sub>O</strong></td><td><code>C2H2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>4</sub>O<sub>4</sub></strong></td><td><code>C2H4O4_FAD_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>6</sub>Cl</strong></td><td><code>C2H6Cl_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>6</sub>F</strong></td><td><code>C2H6F_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>6</sub>H</strong></td><td><code>C2H6H_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>6</sub>O</strong></td><td><code>C2H6O_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>6</sub>OH</strong></td><td><code>C2H6OH_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>H<sub>7</sub></strong></td><td><code>C2H7_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>C<sub>2</sub>N</strong></td><td><code>C2N_PIPNN_Ap</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>C2N_PIPNN_App</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>3</sub>H<sub>7</sub>NO</strong></td><td><code>C3H7NO_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>C<sub>6</sub>H<sub>5</sub>SH</strong></td><td><code>C6H5SH_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C6H5SH_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="2" valign="top"><strong>C<sub>6</sub>H<sub>6</sub>O</strong></td><td><code>PHOH_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>PHOH_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="16" valign="top"><strong>C<sub>7</sub>H<sub>8</sub>S</strong></td><td><code>C7H8S_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_IA_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_IA_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_II_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_II_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_III_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_III_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_S18_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_S18_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_S20_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_S20_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_S22_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_S22_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_S25_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C7H8S_S25_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>2</sub>O</strong></td><td><code>CH2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="3" valign="top"><strong>CH<sub>2</sub>O<sub>2</sub></strong></td><td><code>CH2O2_LOCAL_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>CH2O2_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>CH2O2_Tew_2016</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="4" valign="top"><strong>CH<sub>2</sub>OH</strong></td><td><code>CH2OH_Surfgen_2state_2018</code></td><td>Surfgen multistate potential</td><td align="center">2</td></tr>
<tr><td><code>CH2OH_Surfgen_2state_2018_DPEM</code></td><td>Surfgen multistate diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>CH2OH_Surfgen_3state_2017</code></td><td>Surfgen multistate potential</td><td align="center">3</td></tr>
<tr><td><code>CH2OH_Surfgen_3state_2017_DPEM</code></td><td>Surfgen multistate diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>2</sub><sup>+</sup></strong></td><td><code>CH2p_3x3_DPEM_2025</code></td><td>Neural-network coupled diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="2" valign="top"><strong>CH<sub>3</sub>NH<sub>2</sub></strong></td><td><code>CH3NH2_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">2</td></tr>
<tr><td><code>CH3NH2_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>3</sub>OH</strong></td><td><code>CH3OH_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>4</sub></strong></td><td><code>CH4_GEN_SP_2001</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>4</sub>·(H<sub>2</sub>O)<sub>2</sub></strong></td><td><code>CH4_H2O_H2O_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>4</sub>Br</strong></td><td><code>CH4Br_VBMM_2002</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>CH<sub>4</sub>Cl</strong></td><td><code>CH4Cl_VBMM_2000</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td><code>CH4Cl_VBMM_2005</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>4</sub>CN</strong></td><td><code>CH4CN_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="3" valign="top"><strong>CH<sub>4</sub>F</strong></td><td><code>CH4F_GEN_GC_1996</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>CH4F_GEN_GCSO_1996</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>CH4F_LEPS_RNCG_2006</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>4</sub>Ni</strong></td><td><code>CH4Ni_EAM6_J1_1996</code></td><td>Embedded-atom potential</td><td align="center">1</td></tr>
<tr><td rowspan="4" valign="top"><strong>CH<sub>4</sub>O</strong></td><td><code>CH4O_LEP_1995</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH4O_LEP_1998</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH4O_LEP_2000</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH4O_LEP_2014</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>CH<sub>4</sub>OCl</strong></td><td><code>CH4OCl_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>CH4OCl_VBMM_2023</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>4</sub>OF</strong></td><td><code>CH4OF_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>CH<sub>4</sub>OH</strong></td><td><code>CH4OH_VBMM_2000</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td><code>CH4OH_VBMM_2015</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="5" valign="top"><strong>CH<sub>5</sub></strong></td><td><code>CH5_GEN_J1_1987</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>CH5_GEN_J2_1987</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>CH5_LEPS_JG_1995</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH5_LEPS_JJ_2002</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH5_VBMM_CBG_2009</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>5</sub>O<sub>2</sub></strong></td><td><code>CH5O2_VBMM_2023</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>CH<sub>5</sub><sup>+</sup></strong></td><td><code>CH5p_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="4" valign="top"><strong>Cl<sub>2</sub>H</strong></td><td><code>Cl2H_LEPS_BCMR_1983</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>Cl2H_LEPS_KNPRY_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>Cl2H_LEPS_PK2_1987</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>Cl2H_LEPS_PK3_1987</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="15" valign="top"><strong>ClH<sub>2</sub></strong></td><td><code>ClH2_LEPS_ALTG1_1996</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_ALTG2_1996</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_GTM1_1981</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_GTM2_1981</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_GTM3_1981</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_GTM4_1982</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_KNPRY_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_LB_1980</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_PK1_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_PK2_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_SPK_1973</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_STSBLTG1_1989</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_STSBLTG2_1989</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_TTGI1_1985</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>ClH2_LEPS_TTGI2_1985</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>ClH<sub>2</sub>O</strong></td><td><code>ClH2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>ClNH<sub>3</sub></strong></td><td><code>ClNH3_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="5" valign="top"><strong>CO<sub>2</sub></strong></td><td><code>CO2_1Ap_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td><code>CO2_1App_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td><code>CO2_21Ap_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td><code>CO2_3Ap_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td><code>CO2_3App_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>2</sub>O<sub>2</sub></strong></td><td><code>COCO_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>F<sub>2</sub>H</strong></td><td><code>F2H_LEPS_JOT_1972</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>F2H_LEPS_KNPRY_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="3" valign="top"><strong>F<sub>2</sub>H<sub>2</sub></strong></td><td><code>F2H2_GEN_JBKKL_1990</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>F2H2_GEN_S2_1996</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>F2H2_GEN_SQSBDE_1995</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="25" valign="top"><strong>FH<sub>2</sub></strong></td><td><code>FH2_GEN_5SEC_1991</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECCOG1_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECCOG2_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECCOG3_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECG1_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECG2_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECG3_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECGM4_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECS1_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECS2_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECS3_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_5SECW_1991</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_GEN_6SEC_1993</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_M5_1991</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_M5SK_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO1_1984</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO2_1984</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO2a_1984</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO3_1984</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO4_1984</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO5_1985</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO5a_1985</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO5b_1985</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_NO5z_1985</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>FH2_LEPS_SK_1980</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>FH<sub>2</sub>O</strong></td><td><code>FH2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>FH2O_SOC_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>GeH<sub>4</sub>OH</strong></td><td><code>GeH4OH_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>GeH<sub>5</sub></strong></td><td><code>GeH5_LEPS</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>2</sub>CO</strong></td><td><code>H2CO_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>2</sub>O<sub>2</sub></strong></td><td><code>H2O2_GEN_1998</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>2</sub>OBr</strong></td><td><code>H2OBr_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="11" valign="top"><strong>H<sub>3</sub></strong></td><td><code>H3_GEN_A2_2002</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_A3_2002</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_A4_2002</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_BH_2009</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_BKMP3_1996</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_BKMP_1991</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_CCI_2002</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_DMBE_1987</code></td><td>General analytic potential</td><td align="center">2</td></tr>
<tr><td><code>H3_GEN_LSTH_1978</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_PK2_1964</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>H3_GEN_TK_1984</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>3</sub>Cl</strong></td><td><code>H3Cl_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>H<sub>3</sub>ClO</strong></td><td><code>H3ClO_PIP_2013</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>H3ClO_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>3</sub>FO</strong></td><td><code>H3FO_Tew_2022</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>3</sub><sup>−</sup></strong></td><td><code>H3m_AS_2004</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>3</sub>O<sub>2</sub><sup>−</sup></strong></td><td><code>H3O2m_sGDML_2024</code></td><td>Symmetric gradient-domain machine-learning potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>H<sub>3</sub><sup>+</sup></strong></td><td><code>H3p_KBNN_2002</code></td><td>Analytic potential</td><td align="center">3</td></tr>
<tr><td><code>H3p_UNO_2001</code></td><td>Analytic potential</td><td align="center">3</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>3</sub>S</strong></td><td><code>H3S_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>4</sub></strong></td><td><code>H4_BMKP_2002</code></td><td>BMKP global analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>4</sub>O<sub>2</sub></strong></td><td><code>H4O2_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>5</sub>O<sub>2</sub><sup>+</sup></strong></td><td><code>H5O2p_PIP_2005</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>H<sub>7</sub><sup>+</sup></strong></td><td><code>H7p_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>HClOH</strong></td><td><code>HClOH_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>HeBeO</strong></td><td><code>HeBeO_Klos_2013</code></td><td>Rare-gas--diatom analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>HeH<sub>2</sub><sup>+</sup></strong></td><td><code>HeH2p_DIM_1976</code></td><td>Diatomics-in-molecules potential</td><td align="center">1</td></tr>
<tr><td><code>HeH2p_MTJS_1987</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>HeNO</strong></td><td><code>HeNO_Klos_2000</code></td><td>Rare-gas--diatom analytic potential</td><td align="center">2</td></tr>
<tr><td rowspan="3" valign="top"><strong>HNi</strong></td><td><code>HNi_EAM5_1995</code></td><td>Embedded-atom potential</td><td align="center">1</td></tr>
<tr><td><code>HNi_EAM6_1996</code></td><td>Embedded-atom potential</td><td align="center">1</td></tr>
<tr><td><code>HNi_EAM7_1998</code></td><td>Embedded-atom potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>HO<sub>3</sub></strong></td><td><code>HO3_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>HOBr</strong></td><td><code>HOBr_3App_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>HOCO·Ar</strong></td><td><code>HOCO_Ar_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>HOPO</strong></td><td><code>HOPO_Tew</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="3" valign="top"><strong>I<sub>2</sub>H</strong></td><td><code>I2H_LEPS_KK_1981</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>I2H_LEPS_PT1_1971</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>I2H_LEPS_PT2_1971</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="5" valign="top"><strong>IH<sub>2</sub></strong></td><td><code>IH2_LEPS_DT1_1971</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>IH2_LEPS_DT2_1971</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>IH2_LEPS_KNPRY_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>IH2_LEPS_PPW1_1974</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>IH2_LEPS_RSPTS_1970</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>K<sub>2</sub>Rb<sub>2</sub></strong></td><td><code>K2Rb2_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>KrNO</strong></td><td><code>KrNO_Klos_2009</code></td><td>Rare-gas--diatom analytic potential</td><td align="center">2</td></tr>
<tr><td rowspan="6" valign="top"><strong>LiFH</strong></td><td><code>LiFH_LEPS_JHCTP_2001</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHCTP_2001_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHTP1_2002</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHTP1_2002_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHTP2_2002</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHTP2_2002_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="2" valign="top"><strong>LiH<sub>2</sub><sup>+</sup></strong></td><td><code>LiH2p_Excited_2003</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>LiH2p_Ground_2003</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>C<sub>3</sub>H<sub>4</sub>O<sub>2</sub> (malonaldehyde)</strong></td><td><code>malonaldehyde_MHT_2014</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="10" valign="top"><strong>MCH</strong></td><td><code>MCHSB_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHSB_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MCHSL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHSL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MCHTL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHTL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MCHWB_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHWB_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MCHWL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHWL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="6" valign="top"><strong>MXH</strong></td><td><code>MXHSB_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MXHSB_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MXHSL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MXHSL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MXHWL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MXHWL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="1" valign="top"><strong>N<sub>2</sub>HOC<sup>+</sup></strong></td><td><code>N2HOCp_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="10" valign="top"><strong>N<sub>2</sub>O</strong></td><td><code>N2O_11_1App_2025</code></td><td>Analytic potential</td><td align="center">11</td></tr>
<tr><td><code>N2O_11_1App_AE_2025</code></td><td>Analytic potential</td><td align="center">11</td></tr>
<tr><td><code>N2O_13_1Ap_2025</code></td><td>Analytic potential</td><td align="center">13</td></tr>
<tr><td><code>N2O_13_1Ap_2026</code></td><td>Analytic potential</td><td align="center">13</td></tr>
<tr><td><code>N2O_13_1Ap_AE_2025</code></td><td>Analytic potential</td><td align="center">13</td></tr>
<tr><td><code>N2O_13_1Ap_AE_2026</code></td><td>Analytic potential</td><td align="center">13</td></tr>
<tr><td><code>N2O_3Ap_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N2O_3Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>N2O_3App_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N2O_3App_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>N<sub>2</sub>O<sub>2</sub></strong></td><td><code>N2O2_3A_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N2O2_3A_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>N<sub>3</sub></strong></td><td><code>N3_4App_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>N3_MBP_2017</code></td><td>Many-body polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="6" valign="top"><strong>N<sub>4</sub></strong></td><td><code>N4_1A_MEG2_2020</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>N4_1A_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N4_1A_PIP_2013</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N4_1A_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>N4_1A_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>N4_BDTC_2014</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="6" valign="top"><strong>NaFH</strong></td><td><code>NaFH_LEPS_JHCTP_2001</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_JHCTP_2001_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_TTYPP1_1998</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_TTYPP1_1998_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_TTYPP2_1998</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_TTYPP2_1998_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="14" valign="top"><strong>NaH<sub>2</sub></strong></td><td><code>NaH2_LEPS_5F_1992</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_5F_1992_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_5G_1992</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_5G_1992_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_6_2000</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_6_2000_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_7L_2000</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_7L_2000_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_7LD_2000</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_7LD_2000_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_7S_2000</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_7S_2000_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_7SD_2000</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaH2_LEPS_7SD_2000_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="2" valign="top"><strong>NeNO</strong></td><td><code>NeNO_Klos_2016_v0</code></td><td>Rare-gas--diatom analytic potential</td><td align="center">2</td></tr>
<tr><td><code>NeNO_Klos_2016_v2</code></td><td>Rare-gas--diatom analytic potential</td><td align="center">2</td></tr>
<tr><td rowspan="6" valign="top"><strong>NH<sub>3</sub></strong></td><td><code>NH3_GEN_FFW1_2006</code></td><td>General analytic potential</td><td align="center">2</td></tr>
<tr><td><code>NH3_GEN_FFW1_2006_DPEM</code></td><td>General analytic diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NH3_GEN_FFW2_2007</code></td><td>General analytic potential</td><td align="center">2</td></tr>
<tr><td><code>NH3_GEN_FFW2_2007_DPEM</code></td><td>General analytic diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NH3_Surfgen_2state_2012</code></td><td>Surfgen multistate potential</td><td align="center">2</td></tr>
<tr><td><code>NH3_Surfgen_2state_2012_DPEM</code></td><td>Surfgen multistate diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="1" valign="top"><strong>NH<sub>3</sub>CN</strong></td><td><code>NH3CN_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>NH<sub>3</sub>OH</strong></td><td><code>NH3OH_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>NH<sub>4</sub></strong></td><td><code>NH4_VBMM_1997</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td><code>NH4_VBMM_2010</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="6" valign="top"><strong>NO<sub>2</sub></strong></td><td><code>NO2_2Ap_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_2Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_4Ap_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_4Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_6Ap_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_6Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>O<sub>2</sub>H</strong></td><td><code>O2H_GEN</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>O<sub>2</sub>H<sub>3</sub></strong></td><td><code>O2H3_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="35" valign="top"><strong>O<sub>3</sub></strong></td><td><code>O3_14_3Ap_2022</code></td><td>Analytic potential</td><td align="center">14</td></tr>
<tr><td><code>O3_14_3Ap_2022_DPEM</code></td><td>Analytic diabatic representation</td><td align="center">14</td></tr>
<tr><td><code>O3_14_3Ap_2023</code></td><td>Analytic potential</td><td align="center">14</td></tr>
<tr><td><code>O3_14_3Ap_2023_DPEM</code></td><td>Analytic diabatic representation</td><td align="center">14</td></tr>
<tr><td><code>O3_14_3Ap_2023a</code></td><td>Analytic potential</td><td align="center">14</td></tr>
<tr><td><code>O3_14_3Ap_2023a_DPEM</code></td><td>Analytic diabatic representation</td><td align="center">14</td></tr>
<tr><td><code>O3_14_3Ap_2024</code></td><td>Analytic potential</td><td align="center">14</td></tr>
<tr><td><code>O3_14_3Ap_2024_DPEM</code></td><td>Analytic diabatic representation</td><td align="center">14</td></tr>
<tr><td><code>O3_16_3App_2026</code></td><td>Analytic potential</td><td align="center">16</td></tr>
<tr><td><code>O3_16_3App_AE_2026</code></td><td>Analytic potential</td><td align="center">16</td></tr>
<tr><td><code>O3_1_1Ap</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_1Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_1App</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_1App_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_3Ap</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_3Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_3App</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_3App_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_5Ap</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_5Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_5App</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_1_5App_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_2_1Ap</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_2_1Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_2_3Ap</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_2_3Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_2_5Ap</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_2_5Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O3_5app_2026</code></td><td>Analytic potential</td><td align="center">6</td></tr>
<tr><td><code>O3_5app_AE_2026</code></td><td>Analytic potential</td><td align="center">6</td></tr>
<tr><td><code>O3_6_5Ap_2023</code></td><td>Analytic potential</td><td align="center">6</td></tr>
<tr><td><code>O3_6_5Ap_2023_DPEM</code></td><td>Analytic diabatic representation</td><td align="center">6</td></tr>
<tr><td><code>O3_6_5Ap_2023a</code></td><td>Analytic potential</td><td align="center">6</td></tr>
<tr><td><code>O3_6_5Ap_2023a_DPEM</code></td><td>Analytic diabatic representation</td><td align="center">6</td></tr>
<tr><td><code>O3_MBP_2017</code></td><td>Many-body polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="7" valign="top"><strong>O<sub>4</sub></strong></td><td><code>O4_MBP_singlet_2018</code></td><td>Many-body polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>O4_quintet</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_quintet_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_singlet</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_singlet_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_triplet_v2</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_triplet_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="9" valign="top"><strong>OH<sub>2</sub></strong></td><td><code>OH2_Ap_GEN_M21986</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_App_GEN_M21986</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_App_GEN_M3n1988</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_App_GEN_M3na1988</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_GEN_DIM1981</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_GEN_POL1980</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_GEN_SL1979</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_LEPS_JW1977</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_LEPS_JWS1985</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="5" valign="top"><strong>OH<sub>3</sub></strong></td><td><code>OH3_LEPS_SE_1986</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>OH3_PIP_FFW1_2019</code></td><td>Permutationally invariant polynomial potential</td><td align="center">3</td></tr>
<tr><td><code>OH3_PIP_FFW1_2019_DPEM</code></td><td>Permutationally invariant polynomial diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>OH3_PIP_FFW2_2022</code></td><td>Permutationally invariant polynomial potential</td><td align="center">3</td></tr>
<tr><td><code>OH3_PIP_FFW2_2022_DPEM</code></td><td>Permutationally invariant polynomial diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="1" valign="top"><strong>SiH<sub>4</sub>Cl</strong></td><td><code>SiH4Cl_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>SiH<sub>5</sub></strong></td><td><code>SiH5_LEP_1998v1</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>SiH5_LEP_1998v2</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1" valign="top"><strong>SiOH</strong></td><td><code>SiOH_GEN_1988</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2" valign="top"><strong>YRH</strong></td><td><code>YRH_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>YRH_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
</tbody>
</table>

<!-- END GENERATED PES CATALOG -->

## Installation

The recommended installation builds the modern Fortran extensions in a
dedicated conda environment:

```bash
git clone https://github.com/rmeanapa/chempotpy.git
cd chempotpy
conda env create --file environment.yml
conda activate chempotpy
python -m pip install -e .
make build BUILD_FLAGS="--jobs 4"
```

The package includes its numerical support code and does not require external
MKL, LAPACK, or BLAS libraries. A Fortran compiler is required when building
the extensions from source.

To explore the installed catalog:

```python
import chempotpy

chempotpy.system()          # list chemical systems
chempotpy.list("H3")       # list surfaces for one system
chempotpy.detail("H3")     # show references and capabilities
```

### First calculation

```python
import chempotpy

geometry = [
    ["H", 0.00, 0.00, 0.00],
    ["H", 0.74, 0.00, 0.00],
    ["H", 1.45, 0.20, 0.00],
]

energy = chempotpy.p("H3", "H3_GEN_A2_2002", geometry)
energy, numerical_gradient = chempotpy.pn(
    "H3", "H3_GEN_A2_2002", geometry
)
```

Use `pg` for a PES with an analytic gradient, `pgd` for energies, gradients,
and nonadiabatic couplings, and `u`/`ug`/`un` for diabatic matrices. The
`detail` function reports which quantities each surface provides.

### Connect molecular-simulation programs

Every catalog entry is available through one transport-neutral provider layer.
It accepts Cartesian coordinates and returns energies, gradients, and matrices
in atomic units. Surfaces without analytic gradients use the library's central
finite-difference implementation; a selected adiabatic state is obtained from
a DPEM by symmetric diagonalization and the Hellmann-Feynman gradient.

List the exact `SYSTEM/SURFACE` keys:

```bash
chempotpy-provider list
```

Run the persistent CHEMPOTPY socket service:

```bash
chempotpy-provider socket H3 H3_GEN_A2_2002 --host 127.0.0.1 --port 8021
```

Connect the same surface as an [i-PI](https://ipi-code.org/) force driver:

```bash
chempotpy-provider ipi H3 H3_GEN_A2_2002 \
    --host 127.0.0.1 --port 31415 --state 0
```

Run it as a [MolSSI Driver Interface](https://molssi-mdi.github.io/MDI_Library/)
engine after installing the optional `pymdi` package:

```bash
chempotpy-provider mdi H3 H3_GEN_A2_2002 --state 0 \
    --mdi-options='-name CHEMPOTPY -role ENGINE -method TCP -port 8021'
```

`Aln`, `AlmHn`, and `SiOH` have variable compositions and therefore also
require `--atoms`, for example `--atoms Al,Al,Al,Al`. See the
[PES interface contract](docs/PES_CONTRACT.md#provider-transports) for the
native protocol and supported i-PI/MDI quantities.

## Generate a standalone Fortran PES

Any catalog surface can be exported as one dependency-complete `.f90` source
file. From the repository root:

```bash
python scripts/export_standalone.py \
    --surface H3/H3_GEN_A2_2002 \
    --output standalone/H3_GEN_A2_2002
```

The generated file is:

```text
standalone/H3_GEN_A2_2002/H3_GEN_A2_2002_standalone.f90
```

It contains the selected PES, its embedded parameters, PES-specific support,
and the required CHEMPOTPY core modules. It can be compiled directly with a
program that calls the documented `pes(...)` or `dpem(...)` entry point:

```bash
gfortran -O2 -ffree-form -ffree-line-length-none \
    H3_GEN_A2_2002_standalone.f90 my_program.f90 -o my_program
```

To generate a runnable demonstration program as well:

```bash
python scripts/export_standalone.py \
    --surface H3/H3_GEN_A2_2002 \
    --output standalone/H3_GEN_A2_2002 \
    --with-driver \
    --geometries 2

cd standalone/H3_GEN_A2_2002
gfortran -O2 -ffree-form -ffree-line-length-none \
    H3_GEN_A2_2002_standalone.f90 \
    H3_GEN_A2_2002_example.f90 \
    -o H3_GEN_A2_2002_example
./H3_GEN_A2_2002_example
```

The example reports energies and finite-difference gradients and Hessians for
the selected surface. Use the unambiguous `SYSTEM/SURFACE` selector whenever a
surface name is shared by more than one system.

## Contributing a PES

You do not need to know Git, Python, or modern Fortran to contribute a surface.
The easiest route is the
[new PES submission form](https://github.com/rmeanapa/chempotpy/issues/new?template=new-pes.yml).
Provide the original source and parameter files, the scientific reference,
the atom order and units, and any example geometries with known results. The
maintainers will review the scientific provenance, modernize the source, add
the tests, and prepare the integration pull request with you.

See [Submit a PES](docs/SUBMIT_A_PES.md) for the short submission checklist.
Contributors who prefer to implement a surface directly can follow the
developer procedure in [CONTRIBUTING.md](CONTRIBUTING.md). The helper below
creates a modern, non-overwriting PES scaffold:

```bash
python scripts/add_pes.py SYSTEM SURFACE_NAME \
    --natoms NATOMS \
    --nstates NSTATES \
    --entry-point pes \
    --functional-form "MODEL FAMILY" \
    --common-name "COMMON NAME" \
    --reference "AUTHORS, JOURNAL, YEAR"
```

Ordinary [GitHub issues](https://github.com/rmeanapa/chempotpy/issues) remain
available for bug reports and questions.

## Citation

When using ChemPotPy, cite the library paper and the original publication for
each selected potential energy surface:

> Y. Shu, Z. Varga, D. Zhang, and D. G. Truhlar, “ChemPotPy: A Python Library
> for Analytic Representations of Potential Energy Surfaces and Diabatic
> Potential Energy Matrices,” *J. Phys. Chem. A* **127**, 9635–9640 (2023),
> [doi:10.1021/acs.jpca.3c05899](https://doi.org/10.1021/acs.jpca.3c05899).

## Documentation

- [Modernization and validation](docs/MODERNIZATION.md)
- [PES interface contract](docs/PES_CONTRACT.md)
- [Fortran style](docs/FORTRAN_STYLE.md)
- [PES import evidence](docs/PES_IMPORTS.md)
- [Submit a PES](docs/SUBMIT_A_PES.md)
- [Contribution guide](CONTRIBUTING.md)
- [Version 2.0 release notes](docs/releases/2.0.0.md)
- [License](LICENSE)

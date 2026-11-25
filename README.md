
ChemPotPy, CHEMical library of POTential energy surfaces in PYthon 
==================================================================

[![All PES numerical regression](https://github.com/rmeanapa/chempotpy/actions/workflows/all-pes-numerical-regression.yml/badge.svg?branch=main)](https://github.com/rmeanapa/chempotpy/actions/workflows/all-pes-numerical-regression.yml?query=branch%3Amain)

Jan. 17, 2025

Authors: Yinan Shu, Zoltan Varga, Dayou Zhang, Donald G. Truhlar
University of Minnesota, Minnesota, United States

Fork modernization author and maintainer: Ruben Meana-Paneda

ChemPotPy is a library for analytic representation of single-state 
and multi-state potential energy surfaces and couplings. 

All fortran source code are stored in folder chempotpy 

New surfaces must follow the documented scientific, architectural, and
validation procedure in [CONTRIBUTING.md](CONTRIBUTING.md). The repository
provides `scripts/add_pes.py` to create a non-overwriting modern PES or DPEM
scaffold.


Potential energy surface catalog
--------------------------------

This table groups callable PES and DPEM entries by chemical composition, lists
each surface name and description, and reports its number of electronic surfaces.
Regenerate it with `python scripts/generate_readme_catalog.py` after adding or
changing a surface.

<!-- BEGIN GENERATED PES CATALOG -->

The catalog currently exposes **403 callable surfaces** across **99 chemical systems**. The final column gives the number of electronic surfaces represented by each callable PES or DPEM entry.

<table>
<thead>
<tr><th>System composition</th><th>Surface name</th><th>Surface description</th><th>Electronic surfaces</th></tr>
</thead>
<tbody>
<tr><td rowspan="12"><strong>ALC</strong></td><td><code>ALC5m0_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
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
<tr><td rowspan="2"><strong>Al<sub>m</sub>H<sub>n</sub></strong></td><td><code>AlmHn_VBO1</code></td><td>Valence-bond-order potential</td><td align="center">1</td></tr>
<tr><td><code>AlmHn_VBO2</code></td><td>Valence-bond-order potential</td><td align="center">1</td></tr>
<tr><td rowspan="59"><strong>Al<sub>n</sub></strong></td><td><code>Aln_BetH</code></td><td>Analytic potential</td><td align="center">1</td></tr>
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
<tr><td rowspan="2"><strong>ArNO</strong></td><td><code>ArNO_Ap_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>ArNO_App_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="7"><strong>BrH<sub>2</sub></strong></td><td><code>BrH2_LEPS_C_1982</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrH2_LEPS_LTBZ1_1995</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrH2_LEPS_LTBZ2_1995</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>BrH2_LEPS_LTBZ2_1995_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>BrH2_LEPS_LTBZSO_1995</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>BrH2_LEPS_LTBZSO_1995_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>BrH2_LEPS_W_1959</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="4"><strong>BrHCl</strong></td><td><code>BrHCl_LEPS1</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrHCl_LEPS2</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrHCl_LEPS3</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>BrHCl_LEPS_BLM_1983</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>BrNH<sub>3</sub></strong></td><td><code>BrNH3_VBMM_2026</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>2</sub>N<sub>2</sub>O</strong></td><td><code>C2H2N2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>2</sub>O</strong></td><td><code>C2H2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>4</sub>O<sub>4</sub></strong></td><td><code>C2H4O4_FAD_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>6</sub>Cl</strong></td><td><code>C2H6Cl_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>6</sub>F</strong></td><td><code>C2H6F_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>6</sub>H</strong></td><td><code>C2H6H_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>6</sub>O</strong></td><td><code>C2H6O_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>6</sub>OH</strong></td><td><code>C2H6OH_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>H<sub>7</sub></strong></td><td><code>C2H7_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>C<sub>2</sub>N</strong></td><td><code>C2N_PIPNN_Ap</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>C2N_PIPNN_App</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>3</sub>H<sub>7</sub>NO</strong></td><td><code>C3H7NO_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>C<sub>6</sub>H<sub>5</sub>SH</strong></td><td><code>C6H5SH_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>C6H5SH_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="2"><strong>C<sub>6</sub>H<sub>6</sub>O</strong></td><td><code>PHOH_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
<tr><td><code>PHOH_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="16"><strong>C<sub>7</sub>H<sub>8</sub>S</strong></td><td><code>C7H8S_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">3</td></tr>
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
<tr><td rowspan="1"><strong>CH<sub>2</sub>O</strong></td><td><code>CH2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>CH<sub>2</sub>O<sub>2</sub></strong></td><td><code>CH2O2_LOCAL_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>CH2O2_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="4"><strong>CH<sub>2</sub>OH</strong></td><td><code>CH2OH_Surfgen_2state_2018</code></td><td>Surfgen multistate potential</td><td align="center">2</td></tr>
<tr><td><code>CH2OH_Surfgen_2state_2018_DPEM</code></td><td>Surfgen multistate diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>CH2OH_Surfgen_3state_2017</code></td><td>Surfgen multistate potential</td><td align="center">3</td></tr>
<tr><td><code>CH2OH_Surfgen_3state_2017_DPEM</code></td><td>Surfgen multistate diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="2"><strong>CH<sub>3</sub>NH<sub>2</sub></strong></td><td><code>CH3NH2_APRP</code></td><td>Anchor-points reactive potential</td><td align="center">2</td></tr>
<tr><td><code>CH3NH2_APRP_DPEM</code></td><td>Anchor-points reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="1"><strong>CH<sub>3</sub>OH</strong></td><td><code>CH3OH_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>CH<sub>4</sub></strong></td><td><code>CH4_GEN_SP_2001</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>CH<sub>4</sub>·(H<sub>2</sub>O)<sub>2</sub></strong></td><td><code>CH4_H2O_H2O_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>CH<sub>4</sub>Br</strong></td><td><code>CH4Br_VBMM_2002</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>CH<sub>4</sub>Cl</strong></td><td><code>CH4Cl_VBMM_2000</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td><code>CH4Cl_VBMM_2005</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>CH<sub>4</sub>CN</strong></td><td><code>CH4CN_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="3"><strong>CH<sub>4</sub>F</strong></td><td><code>CH4F_GEN_GC_1996</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>CH4F_GEN_GCSO_1996</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>CH4F_LEPS_RNCG_2006</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>CH<sub>4</sub>Ni</strong></td><td><code>CH4Ni_EAM6_J1_1996</code></td><td>Embedded-atom potential</td><td align="center">1</td></tr>
<tr><td rowspan="4"><strong>CH<sub>4</sub>O</strong></td><td><code>CH4O_LEP_1995</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH4O_LEP_1998</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH4O_LEP_2000</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH4O_LEP_2014</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>CH<sub>4</sub>OCl</strong></td><td><code>CH4OCl_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>CH4OCl_VBMM_2023</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>CH<sub>4</sub>OF</strong></td><td><code>CH4OF_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>CH<sub>4</sub>OH</strong></td><td><code>CH4OH_VBMM_2000</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td><code>CH4OH_VBMM_2015</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="5"><strong>CH<sub>5</sub></strong></td><td><code>CH5_GEN_J1_1987</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>CH5_GEN_J2_1987</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>CH5_LEPS_JG_1995</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH5_LEPS_JJ_2002</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>CH5_VBMM_CBG_2009</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>CH<sub>5</sub>O<sub>2</sub></strong></td><td><code>CH5O2_VBMM_2023</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>CH<sub>5</sub><sup>+</sup></strong></td><td><code>CH5p_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="4"><strong>Cl<sub>2</sub>H</strong></td><td><code>Cl2H_LEPS_BCMR_1983</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>Cl2H_LEPS_KNPRY_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>Cl2H_LEPS_PK2_1987</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>Cl2H_LEPS_PK3_1987</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="15"><strong>ClH<sub>2</sub></strong></td><td><code>ClH2_LEPS_ALTG1_1996</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
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
<tr><td rowspan="1"><strong>ClH<sub>2</sub>O</strong></td><td><code>ClH2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>ClNH<sub>3</sub></strong></td><td><code>ClNH3_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="5"><strong>CO<sub>2</sub></strong></td><td><code>CO2_1Ap_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td><code>CO2_1App_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td><code>CO2_21Ap_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td><code>CO2_3Ap_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td><code>CO2_3App_RKHS</code></td><td>RKHS-fitted potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>2</sub>O<sub>2</sub></strong></td><td><code>COCO_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>F<sub>2</sub>H</strong></td><td><code>F2H_LEPS_JOT_1972</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>F2H_LEPS_KNPRY_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="3"><strong>F<sub>2</sub>H<sub>2</sub></strong></td><td><code>F2H2_GEN_JBKKL_1990</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>F2H2_GEN_S2_1996</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>F2H2_GEN_SQSBDE_1995</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="25"><strong>FH<sub>2</sub></strong></td><td><code>FH2_GEN_5SEC_1991</code></td><td>General analytic potential</td><td align="center">1</td></tr>
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
<tr><td rowspan="2"><strong>FH<sub>2</sub>O</strong></td><td><code>FH2O_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>FH2O_SOC_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>GeH<sub>4</sub>OH</strong></td><td><code>GeH4OH_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>GeH<sub>5</sub></strong></td><td><code>GeH5_LEPS</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>H<sub>2</sub>CO</strong></td><td><code>H2CO_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>H<sub>2</sub>O<sub>2</sub></strong></td><td><code>H2O2_GEN_1998</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>H<sub>2</sub>OBr</strong></td><td><code>H2OBr_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="11"><strong>H<sub>3</sub></strong></td><td><code>H3_GEN_A2_2002</code></td><td>General analytic potential</td><td align="center">1</td></tr>
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
<tr><td rowspan="1"><strong>H<sub>3</sub>Cl</strong></td><td><code>H3Cl_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>H<sub>3</sub>ClO</strong></td><td><code>H3ClO_PIP_2013</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>H3ClO_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>H<sub>3</sub><sup>−</sup></strong></td><td><code>H3m_AS_2004</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>H<sub>3</sub><sup>+</sup></strong></td><td><code>H3p_KBNN_2002</code></td><td>Analytic potential</td><td align="center">3</td></tr>
<tr><td><code>H3p_UNO_2001</code></td><td>Analytic potential</td><td align="center">3</td></tr>
<tr><td rowspan="1"><strong>H<sub>3</sub>S</strong></td><td><code>H3S_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>H<sub>4</sub>O<sub>2</sub></strong></td><td><code>H4O2_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>H<sub>7</sub><sup>+</sup></strong></td><td><code>H7p_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>HClOH</strong></td><td><code>HClOH_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>HeH<sub>2</sub><sup>+</sup></strong></td><td><code>HeH2p_DIM_1976</code></td><td>Diatomics-in-molecules potential</td><td align="center">1</td></tr>
<tr><td><code>HeH2p_MTJS_1987</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="3"><strong>HNi</strong></td><td><code>HNi_EAM5_1995</code></td><td>Embedded-atom potential</td><td align="center">1</td></tr>
<tr><td><code>HNi_EAM6_1996</code></td><td>Embedded-atom potential</td><td align="center">1</td></tr>
<tr><td><code>HNi_EAM7_1998</code></td><td>Embedded-atom potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>HO<sub>3</sub></strong></td><td><code>HO3_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>HOBr</strong></td><td><code>HOBr_3App_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>HOCO·Ar</strong></td><td><code>HOCO_Ar_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="3"><strong>I<sub>2</sub>H</strong></td><td><code>I2H_LEPS_KK_1981</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>I2H_LEPS_PT1_1971</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>I2H_LEPS_PT2_1971</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="5"><strong>IH<sub>2</sub></strong></td><td><code>IH2_LEPS_DT1_1971</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>IH2_LEPS_DT2_1971</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>IH2_LEPS_KNPRY_1966</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>IH2_LEPS_PPW1_1974</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>IH2_LEPS_RSPTS_1970</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>K<sub>2</sub>Rb<sub>2</sub></strong></td><td><code>K2Rb2_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="6"><strong>LiFH</strong></td><td><code>LiFH_LEPS_JHCTP_2001</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHCTP_2001_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHTP1_2002</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHTP1_2002_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHTP2_2002</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>LiFH_LEPS_JHTP2_2002_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="2"><strong>LiH<sub>2</sub><sup>+</sup></strong></td><td><code>LiH2p_Excited_2003</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>LiH2p_Ground_2003</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>C<sub>3</sub>H<sub>4</sub>O<sub>2</sub> (malonaldehyde)</strong></td><td><code>malonaldehyde_MHT_2014</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="10"><strong>MCH</strong></td><td><code>MCHSB_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHSB_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MCHSL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHSL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MCHTL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHTL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MCHWB_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHWB_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MCHWL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MCHWL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="6"><strong>MXH</strong></td><td><code>MXHSB_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MXHSB_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MXHSL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MXHSL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>MXHWL_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>MXHWL_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="1"><strong>N<sub>2</sub>HOC<sup>+</sup></strong></td><td><code>N2HOCp_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="10"><strong>N<sub>2</sub>O</strong></td><td><code>N2O_11_1App_2025</code></td><td>Analytic potential</td><td align="center">11</td></tr>
<tr><td><code>N2O_11_1App_AE_2025</code></td><td>Analytic potential</td><td align="center">11</td></tr>
<tr><td><code>N2O_13_1Ap_2025</code></td><td>Analytic potential</td><td align="center">13</td></tr>
<tr><td><code>N2O_13_1Ap_2026</code></td><td>Analytic potential</td><td align="center">13</td></tr>
<tr><td><code>N2O_13_1Ap_AE_2025</code></td><td>Analytic potential</td><td align="center">13</td></tr>
<tr><td><code>N2O_13_1Ap_AE_2026</code></td><td>Analytic potential</td><td align="center">13</td></tr>
<tr><td><code>N2O_3Ap_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N2O_3Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>N2O_3App_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N2O_3App_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>N<sub>2</sub>O<sub>2</sub></strong></td><td><code>N2O2_3A_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N2O2_3A_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>N<sub>3</sub></strong></td><td><code>N3_4App_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>N3_MBP_2017</code></td><td>Many-body polynomial potential</td><td align="center">1</td></tr>
<tr><td rowspan="6"><strong>N<sub>4</sub></strong></td><td><code>N4_1A_MEG2_2020</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>N4_1A_PIP</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N4_1A_PIP_2013</code></td><td>Permutationally invariant polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>N4_1A_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>N4_1A_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>N4_BDTC_2014</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="6"><strong>NaFH</strong></td><td><code>NaFH_LEPS_JHCTP_2001</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_JHCTP_2001_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_TTYPP1_1998</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_TTYPP1_1998_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_TTYPP2_1998</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>NaFH_LEPS_TTYPP2_1998_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="14"><strong>NaH<sub>2</sub></strong></td><td><code>NaH2_LEPS_5F_1992</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
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
<tr><td rowspan="6"><strong>NH<sub>3</sub></strong></td><td><code>NH3_GEN_FFW1_2006</code></td><td>General analytic potential</td><td align="center">2</td></tr>
<tr><td><code>NH3_GEN_FFW1_2006_DPEM</code></td><td>General analytic diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NH3_GEN_FFW2_2007</code></td><td>General analytic potential</td><td align="center">2</td></tr>
<tr><td><code>NH3_GEN_FFW2_2007_DPEM</code></td><td>General analytic diabatic representation</td><td align="center">2</td></tr>
<tr><td><code>NH3_Surfgen_2state_2012</code></td><td>Surfgen multistate potential</td><td align="center">2</td></tr>
<tr><td><code>NH3_Surfgen_2state_2012_DPEM</code></td><td>Surfgen multistate diabatic representation</td><td align="center">2</td></tr>
<tr><td rowspan="1"><strong>NH<sub>3</sub>CN</strong></td><td><code>NH3CN_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>NH<sub>3</sub>OH</strong></td><td><code>NH3OH_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>NH<sub>4</sub></strong></td><td><code>NH4_VBMM_1997</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td><code>NH4_VBMM_2010</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="6"><strong>NO<sub>2</sub></strong></td><td><code>NO2_2Ap_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_2Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_4Ap_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_4Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_6Ap_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td><code>NO2_6Ap_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>O<sub>2</sub>H</strong></td><td><code>O2H_GEN</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>O<sub>2</sub>H<sub>3</sub></strong></td><td><code>O2H3_PIPNN</code></td><td>PIP neural-network potential</td><td align="center">1</td></tr>
<tr><td rowspan="35"><strong>O<sub>3</sub></strong></td><td><code>O3_14_3Ap_2022</code></td><td>Analytic potential</td><td align="center">14</td></tr>
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
<tr><td rowspan="7"><strong>O<sub>4</sub></strong></td><td><code>O4_MBP_singlet_2018</code></td><td>Many-body polynomial potential</td><td align="center">1</td></tr>
<tr><td><code>O4_quintet</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_quintet_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_singlet</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_singlet_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_triplet_v2</code></td><td>Analytic potential</td><td align="center">1</td></tr>
<tr><td><code>O4_triplet_ZV</code></td><td>ZV analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="9"><strong>OH<sub>2</sub></strong></td><td><code>OH2_Ap_GEN_M21986</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_App_GEN_M21986</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_App_GEN_M3n1988</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_App_GEN_M3na1988</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_GEN_DIM1981</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_GEN_POL1980</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_GEN_SL1979</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_LEPS_JW1977</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>OH2_LEPS_JWS1985</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="5"><strong>OH<sub>3</sub></strong></td><td><code>OH3_LEPS_SE_1986</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>OH3_PIP_FFW1_2019</code></td><td>Permutationally invariant polynomial potential</td><td align="center">3</td></tr>
<tr><td><code>OH3_PIP_FFW1_2019_DPEM</code></td><td>Permutationally invariant polynomial diabatic representation</td><td align="center">3</td></tr>
<tr><td><code>OH3_PIP_FFW2_2022</code></td><td>Permutationally invariant polynomial potential</td><td align="center">3</td></tr>
<tr><td><code>OH3_PIP_FFW2_2022_DPEM</code></td><td>Permutationally invariant polynomial diabatic representation</td><td align="center">3</td></tr>
<tr><td rowspan="1"><strong>SiH<sub>4</sub>Cl</strong></td><td><code>SiH4Cl_VBMM</code></td><td>Valence-bond molecular-mechanics potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>SiH<sub>5</sub></strong></td><td><code>SiH5_LEP_1998v1</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td><code>SiH5_LEP_1998v2</code></td><td>LEPS reactive potential</td><td align="center">1</td></tr>
<tr><td rowspan="1"><strong>SiOH</strong></td><td><code>SiOH_GEN_1988</code></td><td>General analytic potential</td><td align="center">1</td></tr>
<tr><td rowspan="2"><strong>YRH</strong></td><td><code>YRH_LEPS_ModelSurface</code></td><td>LEPS reactive potential</td><td align="center">2</td></tr>
<tr><td><code>YRH_LEPS_ModelSurface_DPEM</code></td><td>LEPS reactive diabatic representation</td><td align="center">2</td></tr>
</tbody>
</table>

<!-- END GENERATED PES CATALOG -->


About this fork
---------------

This repository is the [rmeanapa CHEMPOTPY fork](https://github.com/rmeanapa/chempotpy)
maintained by Ruben Meana-Paneda. It retains the original scientific authors
and citations above while developing a self-contained, modern implementation.
The comparison below is against
[shuyinan/chempotpy](https://github.com/shuyinan/chempotpy) `main` at commit
`4716d22` (October 2, 2026). This fork includes that upstream history and adds
the following changes:

* Converted the complete PES catalog from fixed-form `.f` sources to
  authoritative free-form `.f90` sources and removed obsolete duplicate
  implementations.
* Replaced implicit typing with `implicit none` and explicit declarations;
  introduced private-by-default modules, explicit public APIs, typed
  parameters, `intent(...)`, and documented procedures.
* Removed legacy `COMMON`, `DATA`, mutable `SAVE`, `ENTRY`, `BLOCK DATA`,
  setup-style `prepot`, numeric statement labels, `CONTINUE`, `GOTO`,
  arithmetic `IF`, detached `FORMAT`, and redundant terminal `RETURN`
  constructs from governed PES sources.
* Replaced initialization entry points with stateless, independently owned PES
  modules whose fixed coefficients live in source-controlled module data.
* Added a bundled Fortran core for numerical kinds, atomic elements and
  standard atomic masses, status values, typed PES requests/results/metadata,
  canonical NIST CODATA 2022 and compatibility unit conversions, Cartesian
  geometry, coordinate flattening and reshaping,
  internal-to-Cartesian gradient transformations, centered finite-difference
  validation, reusable internal-coordinate adapters, BLAS operations, and
  symmetric eigensolution. The core keeps CHEMPOTPY
  independent of external MKL, LAPACK, and BLAS libraries.
* Added a Python core contract, adapters, surface registry, runtime dispatch,
  and unit utilities. Python implementations are organized under
  `chempotpy/core/python`; Fortran core modules live directly under
  `chempotpy/core`.
* Kept system- or fit-family support beside its PES sources. For example,
  `chempotpy/H3/chempotpy_h3_2002_basis.f90` owns the shared basis used by the
  A2, A3, A4, and CCI H3 2002 surfaces and is bundled automatically.
* Standardized the Fortran source format: four-space indentation, compact array
  subscripts, consistent declarations and assignments, `elseif`/`endif` and
  `enddo`, inline I/O formats, meaningful PES-specific documentation, and a
  120-column readability target.
* Added catalog discovery and isolated f2py builds so each PES extension
  receives only the bundled core modules it actually imports.
* Added strict compilation, modernization, formatting, packaging, binary
  linkage, extension-import, and numerical-regression gates for all 403
  callable surfaces, including nine stored reference geometries per surface.
* Corrected the CH3NH2 APRP coupling evaluation by initializing terms that are
  absent from the coupling branch and conditionally assigned tent derivatives.
  The legacy result depended on uninitialized stack memory and therefore
  changed with compiler, architecture, and source layout.
* Added a reproducible current conda-forge environment and validation against
  the newest mutually compatible Python, NumPy, Meson, Ninja, and gfortran
  toolchain.
* Documented the shared PES contract, Fortran conventions, build workflow, and
  commands for validating the entire surface catalog.


Mathematical constants and unit conversions
-------------------------------------------

The generated `chempotpy_constants` module provides one working-precision
definition of `PI` together with the NIST CODATA 2022 universal physical
constants. The authoritative CODATA decimal values and provenance are stored
in `chempotpy/core/codata_2022.json`; compiled PES modules therefore do not
read a data file at runtime. Run the following check after changing the
manifest or generator:

    python scripts/generate_constants.py --check

All public PES boundary conversions use the resulting canonical factors:

| Conversion | Value |
| --- | ---: |
| bohr to angstrom | `0.529177210544` |
| hartree to electronvolt | `27.211386245981` |
| hartree/bohr to eV/angstrom | `51.42206751119799` |
| hartree to kJ/mol | `2625.499639479163` |
| hartree to kcal/mol | `627.5094740628975` |
| hartree to cm^-1 | `219474.63136314` |
| electronvolt to kcal/mol | `23.06054783061903` |
| electronvolt to cm^-1 | `8065.543937349212` |

The source values are published in the
[NIST CODATA 2022 extensive listing](https://physics.nist.gov/cuu/pdf/all.pdf).
The module `chempotpy_legacy_units` remains available only for explicitly
requested reproduction of the formerly released rounded conversions. A
conversion-like value that is part of a published PES fit remains local to its
surface and is not replaced merely because it resembles a physical constant;
changing such a value would change the fitted model itself.

Atomic element data
-------------------

The `chempotpy_atomic_elements` core module owns the 118 canonical element
symbols, their atomic numbers, and standard atomic masses in unified atomic
mass units. Its case-insensitive `atomic_number`, `atomic_symbol`, and
`atomic_mass_u` queries return neutral values for invalid input so library
evaluators do not terminate the calling application. `D` and `T` are accepted
as hydrogen aliases for element-identity checks; isotope-specific masses are
outside this standard-element table.

PES sources must not carry private copies of the periodic table. Fit-specific
atom ordering and isotope choices remain with the owning surface, while shared
element queries use this core owner.

The C6H5SH APRP pair and the eight C7H8S APRP adiabatic/diabatic pairs are
documented exceptions for input-coordinate scaling. Their released fits are
numerically sensitive enough that replacing the historical scale changed
representative energies by as much as about `0.18 eV`. They therefore retain
that scale as an explicitly named published-fit constant while using CODATA
2022 for energy, gradient, coupling, and all other public conversions.


Fortran working precision
-------------------------

Modern numerical code imports `wp` from `chempotpy_kinds` and declares real
values as `real(wp)` with explicit literals such as `1.0_wp`. The kinds module
defines portable single, double, and quadruple requests with
`selected_real_kind` and currently selects double precision:

```fortran
use chempotpy_kinds, only: wp
real(wp) :: energy

energy = 0.0_wp
```

The selector is the single `wp = dp` parameter in
`chempotpy/core/chempotpy_kinds.f90`. The shared numerical core, high-level
linear-algebra and interpolation APIs, private BLAS/DSYEV implementation,
generated constants and parameter modules, and the maintained PES catalog now
use this selector consistently. Package builds must
currently leave `wp = dp` because the Python/f2py boundary maps `wp` to a C
double. Standalone Fortran kernels may select another supported kind; a
catalog-wide non-binary64 Python build still requires an explicit conversion
layer at that external boundary.


How to install
--------------
* Create a virtual environment with the newest mutually compatible Python,
  Conda, NumPy, Meson, Ninja, and gfortran packages from conda-forge:
    
       conda env create --file environment.yml
       conda activate chempotpy

* Install stable release:
  
        pip install chempotpy



* Users may compile all extension modules for their current Python and
  platform:

  get into the parent directory of chempotpy 
   
        make build
        make check-modernization

  `make build` uses an unoptimized compatibility profile so the modern Python
  extensions preserve the released Fortran results. Add `--debug` to enable
  bounds checks and Fortran backtraces while diagnosing a kernel.

  `make check-conda BUILD_FLAGS="--jobs 4"` builds every extension, checks its
  native-library linkage, imports each binary in an isolated Python process,
  and evaluates all surfaces at all nine catalog regression geometries.

  CHEMPOTPY carries the numerical support routines used by its PES kernels.
  Building the extensions does not require MKL, LAPACK, or BLAS. A Fortran
  compiler is needed only when building extensions from source; binary wheels
  should carry the compiler runtime required by their platform.

  Modern PES sources may import the bundled Fortran core for numerical kinds,
  canonical or compatibility unit conversions, Cartesian geometry operations,
  typed PES contracts, finite-difference validation, status values, and linear
  algebra. The isolated builder follows those module dependencies and compiles
  the required sources into each extension, so every resulting PES module
  remains independently loadable.

Modernizing Fortran sources
---------------------------

Modernized PES code uses four-space indentation, spaced procedure arguments,
short `!> @brief` documentation, explicit declarations and `implicit none`,
`private` by default in modules, and `use ..., only: ...` imports. Run
`make check-style` to validate every free-form Fortran source. Compilation also
uses `-fimplicit-none`, so an undeclared variable is a hard error.

Declarations use `::`, constants use typed `parameter` declarations, and
intrinsic kinds and character lengths use standard syntax. Source code is
lowercase without changing comments or character literals. Library kernels are
quiet: they contain no external `write`, `stop`, or `error stop` statements.
The optional `chempotpy_status_t` result on each module-level `evaluate`
reports unsupported requests without changing the historical `pes` or `dpem`
adapter ABI.

The formatter uses a 120-column readability target for normal equations rather
than Fortran 77's 72/80-column limit. Generated numerical data may exceed that
target, and all free-form compilation uses unlimited line length.

Within a procedure, named constants needed for array shapes come first,
followed by `intent(...)` dummy arguments and then local variables. Redundant
runs of empty lines are collapsed to one visual separator.

Format specifications are written inline in `read` and `write` statements.
Legacy labeled `FORMAT` statements and redundant `return` statements directly
before a procedure end are rejected by the style gate; intentional early
returns remain valid.

All PES kernels are now authoritative free-form `.f90` sources; no legacy `.f`
implementation remains in the package. The regression suite recompiles every
converted kernel and compares it with the current CODATA 2022 boundary
references. Earlier release-reference fixtures remain as independent migration
evidence where documented. The detailed contract, migration, and
numerical-parity rules are in `docs/PES_CONTRACT.md`; the required source
conventions are collected in `docs/FORTRAN_STYLE.md`.

Run `make check-modernization` for the complete repeatable modernization gate.
It audits the catalog, rejects forbidden legacy Fortran constructs, runs the
Python and core tests, strictly compiles all 403 callable PES interfaces, checks
formatting, and evaluates every surface at nine deterministic geometries in the
same process. `make test-surfaces` runs only the slower catalog-wide numerical
and architecture regression portion.


Generate a standalone Fortran PES
---------------------------------

The standalone exporter selects one catalog surface and writes a single
dependency-complete `.f90` file. The file contains the PES, any surface support
source, and only the bundled Fortran core modules that the PES imports. It can
therefore be compiled directly or added to another Fortran program without
installing CHEMPOTPY, MKL, LAPACK, or BLAS.

For example, export the H3 A2 surface from the repository root:

    python scripts/export_standalone.py \
        --surface H3/H3_GEN_A2_2002 \
        --output standalone/H3_GEN_A2_2002

This creates:

    standalone/H3_GEN_A2_2002/H3_GEN_A2_2002_standalone.f90

Compile that file with the program that will call its historical `pes` or
`dpem` entry point. Free-form sources may contain long fitted-coefficient lines,
so gfortran should be given unlimited free-form line length:

    gfortran -O2 -ffree-form -ffree-line-length-none \
        H3_GEN_A2_2002_standalone.f90 my_program.f90 -o my_program

The generated source preserves the selected surface's public argument list and
documented units while using the versioned CODATA conversion contract, making
it suitable for existing dynamics programs. The surface source documents the
exact `pes(...)` or `dpem(...)` signature. Every released model parameter used
by an installed evaluator is embedded in the owning source or a PES-specific
generated module, so standalone exports perform no runtime model-file I/O.
Historical `path` arguments remain on some entry points only for ABI
compatibility and may be passed an empty string.

The POTLIB CH5p release references two long-range spline tables that neither
POTLIB 2023 nor the original ChemPotPy distribution includes. Its released
short-range PIP is self contained; entering the unavailable long-range branch
reports that provenance defect explicitly instead of searching the working
directory for untracked files.

Add `--with-driver` to also generate a complete example program:

    python scripts/export_standalone.py \
        --surface H3/H3_GEN_A2_2002 \
        --output standalone/H3_GEN_A2_2002 \
        --with-driver \
        --geometries 2

Compile and run the generated pair:

    cd standalone/H3_GEN_A2_2002
    gfortran -O2 -ffree-form -ffree-line-length-none \
        H3_GEN_A2_2002_standalone.f90 \
        H3_GEN_A2_2002_example.f90 \
        -o H3_GEN_A2_2002_example
    ./H3_GEN_A2_2002_example

The example constructs deterministic, non-symmetric Cartesian geometries and
reports each adiabatic energy, or each diabatic-matrix element for a DPEM. It
also reports Cartesian gradients and Hessians calculated by centered finite
differences through the bundled `chempotpy_finite_difference` module, so the
example works even when the original PES publishes energy only. The default
displacement is `1.0e-3` angstrom and can be changed with `--difference-step`.
Consequently, these example derivatives are numerical diagnostics in the
historical output units per angstrom and per angstrom squared; they do not
claim analytic Hessian support from the underlying PES.

If a historical entry point still accepts a data-directory argument, the
generated example accepts the argument for compatibility but does not use it.

Surface names that occur in more than one system must use the unambiguous
`SYSTEM/SURFACE` selector shown above.


Citation
--------

The following paper should be cited in publications utilizing the
ChemPotPy library in addition to the original paper that publishes 
the potential energy surface subroutine:

Shu, Y.; Varga, Z.; Truhlar, D. G.
"ChemPotPy: A Python Library for Analytic Representation of Potential 
Energy Surfaces and Diabatic Potential Energy Matrices"
accepted by J. Phys. Chem. A


TO CONTRIBUTE YOUR POTENTIAL or report a bug
--------------------------------------------
* Option 1
Send email to one of the maintainers:
  - Yinan Shu, yinan.shu.0728@gmail.com
  - Zoltan Varga, zoltan78varga@gmail.com
  - Dayou Zhang, zhan6350@umn.edu
  - Donald G. Truhlar, truhlar@umn.edu
 
* Option 2
Submit tickets on the [issues](https://github.com/shuyinan/chempotpy/issues)


Validating every potential energy surface
-----------------------------------------

CHEMPOTPY provides two complementary numerical regression paths. The Python
extension gate exercises the same compiled modules used by package users,
while the direct-Fortran test independently compiles the authoritative source
of every PES. Both paths evaluate nine deterministic geometries per surface and
compare the results with `tests/data/all_surface_references.json`. The suite
uses three end-to-end pairs and the midpoint of each pair. It retains the two
original non-symmetric endpoints, adds a symmetric non-rigid distortion pair
around the 75/25 interior point, and adds a radial contraction/expansion pair
around the 25/75 interior point. The extra points exercise independent internal
coordinate directions while remaining inside or close to the established test
envelope.

Create the current conda-forge validation environment once, then activate it:

    conda env create --file environment.yml
    conda activate chempotpy

From the repository root, run the complete Python extension gate with:

    make check-conda BUILD_FLAGS="--jobs 4"

This command builds all 403 f2py extensions, checks their native-library
linkage, imports every extension in an isolated Python process, and performs
3,618 numerical evaluations. A successful run ends with counters equivalent to:

    pes_extensions_checked: 403
    missing_or_ambiguous_extensions: 0
    external_mkl_lapack_blas_links: 0
    python_extensions_checked: 403
    python_extension_evaluations: 3618
    python_extension_failures: 0

Run the catalog-wide compiler warning gate with:

    make check-fortran-warnings

This compiles all 462 PES and system-support units with `-Wall`, `-Wextra`, and
`-Werror`. Unused dummy arguments are exempt because the stable PES contracts
intentionally include inputs that some implementations do not consume. Exact
real comparisons are also exempt because released lookup and branch logic uses
them; both categories remain protected by the numerical regression suite.

Run the independent direct-Fortran numerical and architecture regression with:

    make test-surfaces

Every pull request and push also runs the `All PES numerical regression`
GitHub Actions workflow on Ubuntu Linux and macOS. The workflow checks out the
pull request's base commit or the most recent successful push, compiles that
revision and the current revision independently on both operating systems,
evaluates all nine deterministic geometries, and compares every result. A
failed run therefore cannot become the protected baseline. Capturing both
revisions with the same compiler and math runtime on each runner prevents
platform differences from being mistaken for a source change. Each revision
uses its own versioned geometry harness so an explicitly reviewed baseline
transition remains visible. Existing surfaces must remain within a relative
tolerance of `1e-11` and an absolute tolerance of `5e-10`; removed surfaces,
interface changes, non-finite results, compilation failures, and execution
failures also fail the workflow. Aln and NH3 use
deterministic system-specific geometries so their regression points remain
inside the mathematical domain of the released fits. A newly added surface
must compile, execute, match its declared interface, and reproduce the finite,
independently reviewed values committed in the current fixture. The workflow
does not bootstrap a new baseline from the implementation being tested.

Removing the obsolete, unused C2H2O PIP-NN data-path argument establishes
`baseline_version: 5`. The one-step transition from version 4 approves only
that interface change; its numerical values are unchanged. Version 4 records
the SiOH variable-atom driver correction. The transition from version 3
approves evaluating SiOH with its public
`xyz(natoms,3)` shape instead of a legacy oversized array whose storage layout
scrambled the regression geometry. Version 3 records the earlier C2H6F indexing
and uninitialized-state corrections, and version 2 records the CODATA 2022
Hartree-to-kcal/mol migration. Records whose approved values did not change
retain the base revision captured on that same runner, so an unchanged PES is
never compared with another architecture's fixture. Comparisons between two
version-5 revisions remain strict. The
committed fixture supports local and package-level regression tests, while the
workflow's base-revision capture protects existing surfaces independently of
ordinary fixture edits in a proposed change.

To reproduce the workflow comparison locally against a selected Git revision:

    base_root="$(mktemp -d)"
    git worktree add --detach "$base_root/checkout" REVISION
    (cd "$base_root/checkout" && \
        python scripts/capture_library_references.py --all --jobs 4 \
            --output "$base_root/baseline.json")
    python scripts/check_numerical_regression.py \
        --baseline "$base_root/baseline.json" --jobs 4
    git worktree remove --force "$base_root/checkout"

To run only the all-surface numerical comparison from that suite:

    python -m unittest \
        tests.test_all_surfaces.AllSurfaceRegressionTests.test_every_surface_matches_its_reference_outputs \
        -v

An individual extension can be built and evaluated by its system and surface
name. For example:

    python scripts/build_pes.py --surface H3_GEN_A2_2002
    python scripts/check_python_extensions.py \
        --output-root chempotpy \
        --surface H3/H3_GEN_A2_2002

Run every repository test with:

    python -m unittest discover -s tests -v

The extension gate places generated binaries in the package tree by default.
Remove them after local validation with:

    make clean

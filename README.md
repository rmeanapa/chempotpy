
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

This table is generated from the callable Fortran registry, the numerical
reference fixture, and each system's maintained availability declarations.
Regenerate it with `python scripts/generate_readme_catalog.py` after adding or
changing a surface.

<!-- BEGIN GENERATED PES CATALOG -->

The catalog currently exposes **403 callable surfaces** across **99 chemical systems**. Availability codes are: **P** potential energy, **G** Cartesian gradient, **D** nonadiabatic coupling vector, **U** diabatic potential matrix, and **UG** Cartesian gradient of the diabatic matrix.

| System | Surface | Interface | States | Available |
| :--- | :--- | :---: | ---: | :--- |
| `ALC` | `ALC5m0_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC6m0_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7m0_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7m10_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7m1_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7m4_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7m5_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7m7_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7p2_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7p3_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC7p5_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `ALC` | `ALC_LEPS_ModelSurface` | PES | 1 | **P/G** |
| `AlmHn` | `AlmHn_VBO1` | PES | 1 | **P** |
| `AlmHn` | `AlmHn_VBO2` | PES | 1 | **P** |
| `Aln` | `Aln_BetH` | PES | 1 | **P** |
| `Aln` | `Aln_CleR` | PES | 1 | **P** |
| `Aln` | `Aln_CoxJM` | PES | 1 | **P** |
| `Aln` | `Aln_deSPH_LJAT` | PES | 1 | **P** |
| `Aln` | `Aln_deSPH_M` | PES | 1 | **P** |
| `Aln` | `Aln_ER` | PES | 1 | **P** |
| `Aln` | `Aln_ER2` | PES | 1 | **P** |
| `Aln` | `Aln_ER2AT` | PES | 1 | **P** |
| `Aln` | `Aln_ER2BA` | PES | 1 | **P** |
| `Aln` | `Aln_ER2CN` | PES | 1 | **P** |
| `Aln` | `Aln_ER2CoxJM` | PES | 1 | **P** |
| `Aln` | `Aln_ER2EACN` | PES | 1 | **P** |
| `Aln` | `Aln_ER2EAT` | PES | 1 | **P** |
| `Aln` | `Aln_ER2EBA` | PES | 1 | **P** |
| `Aln` | `Aln_ER2ECN` | PES | 1 | **P** |
| `Aln` | `Aln_ER2ES` | PES | 1 | **P** |
| `Aln` | `Aln_ER2ESCNa` | PES | 1 | **P** |
| `Aln` | `Aln_ER2ESCNm` | PES | 1 | **P** |
| `Aln` | `Aln_ER2GEA` | PES | 1 | **P** |
| `Aln` | `Aln_ER2MeiD` | PES | 1 | **P** |
| `Aln` | `Aln_ER2MisFMP` | PES | 1 | **P** |
| `Aln` | `Aln_ER2S` | PES | 1 | **P** |
| `Aln` | `Aln_ER2SCNa` | PES | 1 | **P** |
| `Aln` | `Aln_ER2SCNm` | PES | 1 | **P** |
| `Aln` | `Aln_ErkIV` | PES | 1 | **P** |
| `Aln` | `Aln_ErkV` | PES | 1 | **P** |
| `Aln` | `Aln_ErkVIII` | PES | 1 | **P** |
| `Aln` | `Aln_Gol` | PES | 1 | **P** |
| `Aln` | `Aln_HalP` | PES | 1 | **P** |
| `Aln` | `Aln_Jac` | PES | 1 | **P** |
| `Aln` | `Aln_MeiD` | PES | 1 | **P** |
| `Aln` | `Aln_MisFMP` | PES | 1 | **P** |
| `Aln` | `Aln_NP_A` | PES | 1 | **P** |
| `Aln` | `Aln_NP_Ad` | PES | 1 | **P/G** |
| `Aln` | `Aln_NP_B` | PES | 1 | **P** |
| `Aln` | `Aln_NP_Bd` | PES | 1 | **P/G** |
| `Aln` | `Aln_NP_C` | PES | 1 | **P** |
| `Aln` | `Aln_NP_D` | PES | 1 | **P** |
| `Aln` | `Aln_NP_E` | PES | 1 | **P** |
| `Aln` | `Aln_PapCEP` | PES | 1 | **P** |
| `Aln` | `Aln_PapKEP` | PES | 1 | **P** |
| `Aln` | `Aln_PeaTHT` | PES | 1 | **P** |
| `Aln` | `Aln_PetW` | PES | 1 | **P** |
| `Aln` | `Aln_reCoxJM` | PES | 1 | **P** |
| `Aln` | `Aln_redeSPH_M` | PES | 1 | **P** |
| `Aln` | `Aln_reErkIV` | PES | 1 | **P** |
| `Aln` | `Aln_reErkVIII` | PES | 1 | **P** |
| `Aln` | `Aln_reGEA` | PES | 1 | **P** |
| `Aln` | `Aln_reGol` | PES | 1 | **P** |
| `Aln` | `Aln_reHalP` | PES | 1 | **P** |
| `Aln` | `Aln_reJac` | PES | 1 | **P** |
| `Aln` | `Aln_reLJAT` | PES | 1 | **P** |
| `Aln` | `Aln_reMeiD` | PES | 1 | **P** |
| `Aln` | `Aln_reMisFMP` | PES | 1 | **P** |
| `Aln` | `Aln_rePetW` | PES | 1 | **P** |
| `Aln` | `Aln_reStrM` | PES | 1 | **P** |
| `Aln` | `Aln_reSutC` | PES | 1 | **P** |
| `Aln` | `Aln_StrM` | PES | 1 | **P** |
| `Aln` | `Aln_SutC` | PES | 1 | **P** |
| `ArNO` | `ArNO_Ap_PIPNN` | PES | 1 | **P** |
| `ArNO` | `ArNO_App_PIPNN` | PES | 1 | **P** |
| `BrH2` | `BrH2_LEPS_C_1982` | PES | 1 | **P/G** |
| `BrH2` | `BrH2_LEPS_LTBZ1_1995` | PES | 1 | **P/G** |
| `BrH2` | `BrH2_LEPS_LTBZ2_1995` | PES | 2 | **P/G/D** |
| `BrH2` | `BrH2_LEPS_LTBZ2_1995_DPEM` | DPEM | 2 | **U/UG** |
| `BrH2` | `BrH2_LEPS_LTBZSO_1995` | PES | 2 | **P** |
| `BrH2` | `BrH2_LEPS_LTBZSO_1995_DPEM` | DPEM | 2 | **U** |
| `BrH2` | `BrH2_LEPS_W_1959` | PES | 1 | **P/G** |
| `BrHCl` | `BrHCl_LEPS1` | PES | 1 | **P/G** |
| `BrHCl` | `BrHCl_LEPS2` | PES | 1 | **P/G** |
| `BrHCl` | `BrHCl_LEPS3` | PES | 1 | **P/G** |
| `BrHCl` | `BrHCl_LEPS_BLM_1983` | PES | 1 | **P/G** |
| `BrNH3` | `BrNH3_VBMM_2026` | PES | 1 | **P/G** |
| `C2H2N2O` | `C2H2N2O_PIPNN` | PES | 1 | **P** |
| `C2H2O` | `C2H2O_PIPNN` | PES | 1 | **P** |
| `C2H4O4` | `C2H4O4_FAD_PIP` | PES | 1 | **P** |
| `C2H6Cl` | `C2H6Cl_VBMM` | PES | 1 | **P** |
| `C2H6F` | `C2H6F_VBMM` | PES | 1 | **P** |
| `C2H6H` | `C2H6H_VBMM` | PES | 1 | **P** |
| `C2H6O` | `C2H6O_VBMM` | PES | 1 | **P** |
| `C2H6OH` | `C2H6OH_VBMM` | PES | 1 | **P** |
| `C2H7` | `C2H7_VBMM` | PES | 1 | **P** |
| `C2N` | `C2N_PIPNN_Ap` | PES | 1 | **P** |
| `C2N` | `C2N_PIPNN_App` | PES | 1 | **P** |
| `C3H7NO` | `C3H7NO_PIP` | PES | 1 | **P/G** |
| `C6H5SH` | `C6H5SH_APRP` | PES | 3 | **P/G/D** |
| `C6H5SH` | `C6H5SH_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C6H6O` | `PHOH_APRP` | PES | 3 | **P/G/D** |
| `C6H6O` | `PHOH_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C7H8S` | `C7H8S_APRP` | PES | 3 | **P/G/D** |
| `C7H8S` | `C7H8S_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C7H8S` | `C7H8S_IA_APRP` | PES | 3 | **P/G/D** |
| `C7H8S` | `C7H8S_IA_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C7H8S` | `C7H8S_II_APRP` | PES | 3 | **P/G/D** |
| `C7H8S` | `C7H8S_II_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C7H8S` | `C7H8S_III_APRP` | PES | 3 | **P/G/D** |
| `C7H8S` | `C7H8S_III_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C7H8S` | `C7H8S_S18_APRP` | PES | 3 | **P/G/D** |
| `C7H8S` | `C7H8S_S18_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C7H8S` | `C7H8S_S20_APRP` | PES | 3 | **P/G/D** |
| `C7H8S` | `C7H8S_S20_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C7H8S` | `C7H8S_S22_APRP` | PES | 3 | **P/G/D** |
| `C7H8S` | `C7H8S_S22_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `C7H8S` | `C7H8S_S25_APRP` | PES | 3 | **P/G/D** |
| `C7H8S` | `C7H8S_S25_APRP_DPEM` | DPEM | 3 | **U/UG** |
| `CH2O` | `CH2O_PIPNN` | PES | 1 | **P** |
| `CH2O2` | `CH2O2_LOCAL_PIPNN` | PES | 1 | **P** |
| `CH2O2` | `CH2O2_PIPNN` | PES | 1 | **P** |
| `CH2OH` | `CH2OH_Surfgen_2state_2018` | PES | 2 | **P/G/D** |
| `CH2OH` | `CH2OH_Surfgen_2state_2018_DPEM` | DPEM | 2 | **U/UG** |
| `CH2OH` | `CH2OH_Surfgen_3state_2017` | PES | 3 | **P/G/D** |
| `CH2OH` | `CH2OH_Surfgen_3state_2017_DPEM` | DPEM | 3 | **U/UG** |
| `CH3NH2` | `CH3NH2_APRP` | PES | 2 | **P/G/D** |
| `CH3NH2` | `CH3NH2_APRP_DPEM` | DPEM | 2 | **U/UG** |
| `CH3OH` | `CH3OH_PIP` | PES | 1 | **P** |
| `CH4` | `CH4_GEN_SP_2001` | PES | 1 | **P** |
| `CH4_H2O_H2O` | `CH4_H2O_H2O_PIP` | PES | 1 | **P** |
| `CH4Br` | `CH4Br_VBMM_2002` | PES | 1 | **P/G** |
| `CH4Cl` | `CH4Cl_VBMM_2000` | PES | 1 | **P/G** |
| `CH4Cl` | `CH4Cl_VBMM_2005` | PES | 1 | **P/G** |
| `CH4CN` | `CH4CN_VBMM` | PES | 1 | **P/G** |
| `CH4F` | `CH4F_GEN_GC_1996` | PES | 1 | **P/G** |
| `CH4F` | `CH4F_GEN_GCSO_1996` | PES | 1 | **P/G** |
| `CH4F` | `CH4F_LEPS_RNCG_2006` | PES | 1 | **P/G** |
| `CH4Ni` | `CH4Ni_EAM6_J1_1996` | PES | 1 | **P/G** |
| `CH4O` | `CH4O_LEP_1995` | PES | 1 | **P/G** |
| `CH4O` | `CH4O_LEP_1998` | PES | 1 | **P/G** |
| `CH4O` | `CH4O_LEP_2000` | PES | 1 | **P/G** |
| `CH4O` | `CH4O_LEP_2014` | PES | 1 | **P/G** |
| `CH4OCl` | `CH4OCl_PIPNN` | PES | 1 | **P** |
| `CH4OCl` | `CH4OCl_VBMM_2023` | PES | 1 | **P** |
| `CH4OF` | `CH4OF_PIPNN` | PES | 1 | **P** |
| `CH4OH` | `CH4OH_VBMM_2000` | PES | 1 | **P/G** |
| `CH4OH` | `CH4OH_VBMM_2015` | PES | 1 | **P/G** |
| `CH5` | `CH5_GEN_J1_1987` | PES | 1 | **P/G** |
| `CH5` | `CH5_GEN_J2_1987` | PES | 1 | **P** |
| `CH5` | `CH5_LEPS_JG_1995` | PES | 1 | **P** |
| `CH5` | `CH5_LEPS_JJ_2002` | PES | 1 | **P** |
| `CH5` | `CH5_VBMM_CBG_2009` | PES | 1 | **P** |
| `CH5O2` | `CH5O2_VBMM_2023` | PES | 1 | **P** |
| `CH5p` | `CH5p_PIP` | PES | 1 | **P** |
| `Cl2H` | `Cl2H_LEPS_BCMR_1983` | PES | 1 | **P** |
| `Cl2H` | `Cl2H_LEPS_KNPRY_1966` | PES | 1 | **P** |
| `Cl2H` | `Cl2H_LEPS_PK2_1987` | PES | 1 | **P** |
| `Cl2H` | `Cl2H_LEPS_PK3_1987` | PES | 1 | **P** |
| `ClH2` | `ClH2_LEPS_ALTG1_1996` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_ALTG2_1996` | PES | 1 | **P** |
| `ClH2` | `ClH2_LEPS_GTM1_1981` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_GTM2_1981` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_GTM3_1981` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_GTM4_1982` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_KNPRY_1966` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_LB_1980` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_PK1_1966` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_PK2_1966` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_SPK_1973` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_STSBLTG1_1989` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_STSBLTG2_1989` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_TTGI1_1985` | PES | 1 | **P/G** |
| `ClH2` | `ClH2_LEPS_TTGI2_1985` | PES | 1 | **P/G** |
| `ClH2O` | `ClH2O_PIPNN` | PES | 1 | **P** |
| `ClNH3` | `ClNH3_VBMM` | PES | 1 | **P/G** |
| `CO2` | `CO2_1Ap_RKHS` | PES | 1 | **P/G** |
| `CO2` | `CO2_1App_RKHS` | PES | 1 | **P/G** |
| `CO2` | `CO2_21Ap_RKHS` | PES | 1 | **P/G** |
| `CO2` | `CO2_3Ap_RKHS` | PES | 1 | **P/G** |
| `CO2` | `CO2_3App_RKHS` | PES | 1 | **P/G** |
| `COCO` | `COCO_PIPNN` | PES | 1 | **P/G** |
| `F2H` | `F2H_LEPS_JOT_1972` | PES | 1 | **P/G** |
| `F2H` | `F2H_LEPS_KNPRY_1966` | PES | 1 | **P/G** |
| `F2H2` | `F2H2_GEN_JBKKL_1990` | PES | 1 | **P** |
| `F2H2` | `F2H2_GEN_S2_1996` | PES | 1 | **P** |
| `F2H2` | `F2H2_GEN_SQSBDE_1995` | PES | 1 | **P** |
| `FH2` | `FH2_GEN_5SEC_1991` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECCOG1_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECCOG2_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECCOG3_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECG1_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECG2_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECG3_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECGM4_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECS1_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECS2_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECS3_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_5SECW_1991` | PES | 1 | **P/G** |
| `FH2` | `FH2_GEN_6SEC_1993` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_M5_1991` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_M5SK_1966` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO1_1984` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO2_1984` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO2a_1984` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO3_1984` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO4_1984` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO5_1985` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO5a_1985` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO5b_1985` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_NO5z_1985` | PES | 1 | **P/G** |
| `FH2` | `FH2_LEPS_SK_1980` | PES | 1 | **P/G** |
| `FH2O` | `FH2O_PIPNN` | PES | 1 | **P** |
| `FH2O` | `FH2O_SOC_PIPNN` | PES | 1 | **P** |
| `GeH4OH` | `GeH4OH_VBMM` | PES | 1 | **P/G** |
| `GeH5` | `GeH5_LEPS` | PES | 1 | **P/G** |
| `H2CO` | `H2CO_PIP` | PES | 1 | **P** |
| `H2O2` | `H2O2_GEN_1998` | PES | 1 | **P/G** |
| `H2OBr` | `H2OBr_PIP` | PES | 1 | **P** |
| `H3` | `H3_GEN_A2_2002` | PES | 1 | **P** |
| `H3` | `H3_GEN_A3_2002` | PES | 1 | **P** |
| `H3` | `H3_GEN_A4_2002` | PES | 1 | **P** |
| `H3` | `H3_GEN_BH_2009` | PES | 1 | **P** |
| `H3` | `H3_GEN_BKMP3_1996` | PES | 1 | **P** |
| `H3` | `H3_GEN_BKMP_1991` | PES | 1 | **P/G** |
| `H3` | `H3_GEN_CCI_2002` | PES | 1 | **P** |
| `H3` | `H3_GEN_DMBE_1987` | PES | 2 | **P/G/D** |
| `H3` | `H3_GEN_LSTH_1978` | PES | 1 | **P/G** |
| `H3` | `H3_GEN_PK2_1964` | PES | 1 | **P/G** |
| `H3` | `H3_GEN_TK_1984` | PES | 1 | **P/G** |
| `H3Cl` | `H3Cl_PIPNN` | PES | 1 | **P** |
| `H3ClO` | `H3ClO_PIP_2013` | PES | 1 | **P** |
| `H3ClO` | `H3ClO_PIPNN` | PES | 1 | **P** |
| `H3m` | `H3m_AS_2004` | PES | 1 | **P** |
| `H3p` | `H3p_KBNN_2002` | PES | 3 | **P/G/D** |
| `H3p` | `H3p_UNO_2001` | PES | 3 | **P/G/D** |
| `H3S` | `H3S_PIPNN` | PES | 1 | **P/G** |
| `H4O2` | `H4O2_PIP` | PES | 1 | **P** |
| `H7p` | `H7p_PIP` | PES | 1 | **P** |
| `HClOH` | `HClOH_PIPNN` | PES | 1 | **P** |
| `HeH2p` | `HeH2p_DIM_1976` | PES | 1 | **P/G** |
| `HeH2p` | `HeH2p_MTJS_1987` | PES | 1 | **P** |
| `HNi` | `HNi_EAM5_1995` | PES | 1 | **P/G** |
| `HNi` | `HNi_EAM6_1996` | PES | 1 | **P/G** |
| `HNi` | `HNi_EAM7_1998` | PES | 1 | **P/G** |
| `HO3` | `HO3_PIPNN` | PES | 1 | **P** |
| `HOBr` | `HOBr_3App_PIP` | PES | 1 | **P/G** |
| `HOCO_Ar` | `HOCO_Ar_PIP` | PES | 1 | **P** |
| `I2H` | `I2H_LEPS_KK_1981` | PES | 1 | **P/G** |
| `I2H` | `I2H_LEPS_PT1_1971` | PES | 1 | **P/G** |
| `I2H` | `I2H_LEPS_PT2_1971` | PES | 1 | **P/G** |
| `IH2` | `IH2_LEPS_DT1_1971` | PES | 1 | **P/G** |
| `IH2` | `IH2_LEPS_DT2_1971` | PES | 1 | **P/G** |
| `IH2` | `IH2_LEPS_KNPRY_1966` | PES | 1 | **P/G** |
| `IH2` | `IH2_LEPS_PPW1_1974` | PES | 1 | **P/G** |
| `IH2` | `IH2_LEPS_RSPTS_1970` | PES | 1 | **P/G** |
| `K2Rb2` | `K2Rb2_PIPNN` | PES | 1 | **P** |
| `LiFH` | `LiFH_LEPS_JHCTP_2001` | PES | 2 | **P/G/D** |
| `LiFH` | `LiFH_LEPS_JHCTP_2001_DPEM` | DPEM | 2 | **U/UG** |
| `LiFH` | `LiFH_LEPS_JHTP1_2002` | PES | 2 | **P** |
| `LiFH` | `LiFH_LEPS_JHTP1_2002_DPEM` | DPEM | 2 | **U** |
| `LiFH` | `LiFH_LEPS_JHTP2_2002` | PES | 2 | **P/G/D** |
| `LiFH` | `LiFH_LEPS_JHTP2_2002_DPEM` | DPEM | 2 | **U/UG** |
| `LiH2p` | `LiH2p_Excited_2003` | PES | 1 | **P** |
| `LiH2p` | `LiH2p_Ground_2003` | PES | 1 | **P** |
| `malonaldehyde` | `malonaldehyde_MHT_2014` | PES | 1 | **P/G** |
| `MCH` | `MCHSB_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `MCH` | `MCHSB_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |
| `MCH` | `MCHSL_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `MCH` | `MCHSL_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |
| `MCH` | `MCHTL_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `MCH` | `MCHTL_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |
| `MCH` | `MCHWB_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `MCH` | `MCHWB_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |
| `MCH` | `MCHWL_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `MCH` | `MCHWL_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |
| `MXH` | `MXHSB_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `MXH` | `MXHSB_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |
| `MXH` | `MXHSL_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `MXH` | `MXHSL_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |
| `MXH` | `MXHWL_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `MXH` | `MXHWL_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |
| `N2HOCp` | `N2HOCp_PIPNN` | PES | 1 | **P** |
| `N2O` | `N2O_11_1App_2025` | PES | 11 | **P/G/D** |
| `N2O` | `N2O_11_1App_AE_2025` | PES | 11 | **P/G/D** |
| `N2O` | `N2O_13_1Ap_2025` | PES | 13 | **P/G/D** |
| `N2O` | `N2O_13_1Ap_2026` | PES | 13 | **P/G/D** |
| `N2O` | `N2O_13_1Ap_AE_2025` | PES | 13 | **P/G/D** |
| `N2O` | `N2O_13_1Ap_AE_2026` | PES | 13 | **P/G/D** |
| `N2O` | `N2O_3Ap_PIP` | PES | 1 | **P/G** |
| `N2O` | `N2O_3Ap_ZV` | PES | 1 | **P/G** |
| `N2O` | `N2O_3App_PIP` | PES | 1 | **P/G** |
| `N2O` | `N2O_3App_ZV` | PES | 1 | **P/G** |
| `N2O2` | `N2O2_3A_PIP` | PES | 1 | **P/G** |
| `N2O2` | `N2O2_3A_ZV` | PES | 1 | **P/G** |
| `N3` | `N3_4App_ZV` | PES | 1 | **P/G** |
| `N3` | `N3_MBP_2017` | PES | 1 | **P** |
| `N4` | `N4_1A_MEG2_2020` | PES | 1 | **P/G** |
| `N4` | `N4_1A_PIP` | PES | 1 | **P/G** |
| `N4` | `N4_1A_PIP_2013` | PES | 1 | **P/G** |
| `N4` | `N4_1A_PIPNN` | PES | 1 | **P/G** |
| `N4` | `N4_1A_ZV` | PES | 1 | **P/G** |
| `N4` | `N4_BDTC_2014` | PES | 1 | **P/G** |
| `NaFH` | `NaFH_LEPS_JHCTP_2001` | PES | 2 | **P/G/D** |
| `NaFH` | `NaFH_LEPS_JHCTP_2001_DPEM` | DPEM | 2 | **U/UG** |
| `NaFH` | `NaFH_LEPS_TTYPP1_1998` | PES | 2 | **P** |
| `NaFH` | `NaFH_LEPS_TTYPP1_1998_DPEM` | DPEM | 2 | **U/UG** |
| `NaFH` | `NaFH_LEPS_TTYPP2_1998` | PES | 2 | **P/G/D** |
| `NaFH` | `NaFH_LEPS_TTYPP2_1998_DPEM` | DPEM | 2 | **U/UG** |
| `NaH2` | `NaH2_LEPS_5F_1992` | PES | 2 | **P** |
| `NaH2` | `NaH2_LEPS_5F_1992_DPEM` | DPEM | 2 | **U** |
| `NaH2` | `NaH2_LEPS_5G_1992` | PES | 2 | **P** |
| `NaH2` | `NaH2_LEPS_5G_1992_DPEM` | DPEM | 2 | **U** |
| `NaH2` | `NaH2_LEPS_6_2000` | PES | 2 | **P/G/D** |
| `NaH2` | `NaH2_LEPS_6_2000_DPEM` | DPEM | 2 | **U/UG** |
| `NaH2` | `NaH2_LEPS_7L_2000` | PES | 2 | **P** |
| `NaH2` | `NaH2_LEPS_7L_2000_DPEM` | DPEM | 2 | **U** |
| `NaH2` | `NaH2_LEPS_7LD_2000` | PES | 2 | **P/G/D** |
| `NaH2` | `NaH2_LEPS_7LD_2000_DPEM` | DPEM | 2 | **U/UG** |
| `NaH2` | `NaH2_LEPS_7S_2000` | PES | 2 | **P** |
| `NaH2` | `NaH2_LEPS_7S_2000_DPEM` | DPEM | 2 | **U** |
| `NaH2` | `NaH2_LEPS_7SD_2000` | PES | 2 | **P/G/D** |
| `NaH2` | `NaH2_LEPS_7SD_2000_DPEM` | DPEM | 2 | **U/UG** |
| `NH3` | `NH3_GEN_FFW1_2006` | PES | 2 | **P/G/D** |
| `NH3` | `NH3_GEN_FFW1_2006_DPEM` | DPEM | 2 | **U/UG** |
| `NH3` | `NH3_GEN_FFW2_2007` | PES | 2 | **P/G/D** |
| `NH3` | `NH3_GEN_FFW2_2007_DPEM` | DPEM | 2 | **U/UG** |
| `NH3` | `NH3_Surfgen_2state_2012` | PES | 2 | **P/G/D** |
| `NH3` | `NH3_Surfgen_2state_2012_DPEM` | DPEM | 2 | **U/UG** |
| `NH3CN` | `NH3CN_VBMM` | PES | 1 | **P/G** |
| `NH3OH` | `NH3OH_VBMM` | PES | 1 | **P/G** |
| `NH4` | `NH4_VBMM_1997` | PES | 1 | **P/G** |
| `NH4` | `NH4_VBMM_2010` | PES | 1 | **P/G** |
| `NO2` | `NO2_2Ap_PIPNN` | PES | 1 | **P/G** |
| `NO2` | `NO2_2Ap_ZV` | PES | 1 | **P/G** |
| `NO2` | `NO2_4Ap_PIPNN` | PES | 1 | **P/G** |
| `NO2` | `NO2_4Ap_ZV` | PES | 1 | **P/G** |
| `NO2` | `NO2_6Ap_PIPNN` | PES | 1 | **P/G** |
| `NO2` | `NO2_6Ap_ZV` | PES | 1 | **P/G** |
| `O2H` | `O2H_GEN` | PES | 1 | **P/G** |
| `O2H3` | `O2H3_PIPNN` | PES | 1 | **P** |
| `O3` | `O3_14_3Ap_2022` | PES | 14 | **P/G/D** |
| `O3` | `O3_14_3Ap_2022_DPEM` | DPEM | 14 | **U/UG** |
| `O3` | `O3_14_3Ap_2023` | PES | 14 | **P/G/D** |
| `O3` | `O3_14_3Ap_2023_DPEM` | DPEM | 14 | **U/UG** |
| `O3` | `O3_14_3Ap_2023a` | PES | 14 | **P/G/D** |
| `O3` | `O3_14_3Ap_2023a_DPEM` | DPEM | 14 | **U/UG** |
| `O3` | `O3_14_3Ap_2024` | PES | 14 | **P/G/D** |
| `O3` | `O3_14_3Ap_2024_DPEM` | DPEM | 14 | **U/UG** |
| `O3` | `O3_16_3App_2026` | PES | 16 | **P/G** |
| `O3` | `O3_16_3App_AE_2026` | PES | 16 | **P/G** |
| `O3` | `O3_1_1Ap` | PES | 1 | **P/G** |
| `O3` | `O3_1_1Ap_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_1_1App` | PES | 1 | **P/G** |
| `O3` | `O3_1_1App_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_1_3Ap` | PES | 1 | **P/G** |
| `O3` | `O3_1_3Ap_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_1_3App` | PES | 1 | **P/G** |
| `O3` | `O3_1_3App_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_1_5Ap` | PES | 1 | **P/G** |
| `O3` | `O3_1_5Ap_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_1_5App` | PES | 1 | **P/G** |
| `O3` | `O3_1_5App_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_2_1Ap` | PES | 1 | **P/G** |
| `O3` | `O3_2_1Ap_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_2_3Ap` | PES | 1 | **P/G** |
| `O3` | `O3_2_3Ap_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_2_5Ap` | PES | 1 | **P/G** |
| `O3` | `O3_2_5Ap_ZV` | PES | 1 | **P/G** |
| `O3` | `O3_5app_2026` | PES | 6 | **P/G/D** |
| `O3` | `O3_5app_AE_2026` | PES | 6 | **P/G/D** |
| `O3` | `O3_6_5Ap_2023` | PES | 6 | **P/G/D** |
| `O3` | `O3_6_5Ap_2023_DPEM` | DPEM | 6 | **U/UG** |
| `O3` | `O3_6_5Ap_2023a` | PES | 6 | **P/G/D** |
| `O3` | `O3_6_5Ap_2023a_DPEM` | DPEM | 6 | **U/UG** |
| `O3` | `O3_MBP_2017` | PES | 1 | **P** |
| `O4` | `O4_MBP_singlet_2018` | PES | 1 | **P** |
| `O4` | `O4_quintet` | PES | 1 | **P/G** |
| `O4` | `O4_quintet_ZV` | PES | 1 | **P/G** |
| `O4` | `O4_singlet` | PES | 1 | **P/G** |
| `O4` | `O4_singlet_ZV` | PES | 1 | **P/G** |
| `O4` | `O4_triplet_v2` | PES | 1 | **P/G** |
| `O4` | `O4_triplet_ZV` | PES | 1 | **P/G** |
| `OH2` | `OH2_Ap_GEN_M21986` | PES | 1 | **P/G** |
| `OH2` | `OH2_App_GEN_M21986` | PES | 1 | **P/G** |
| `OH2` | `OH2_App_GEN_M3n1988` | PES | 1 | **P/G** |
| `OH2` | `OH2_App_GEN_M3na1988` | PES | 1 | **P/G** |
| `OH2` | `OH2_GEN_DIM1981` | PES | 1 | **P/G** |
| `OH2` | `OH2_GEN_POL1980` | PES | 1 | **P/G** |
| `OH2` | `OH2_GEN_SL1979` | PES | 1 | **P/G** |
| `OH2` | `OH2_LEPS_JW1977` | PES | 1 | **P/G** |
| `OH2` | `OH2_LEPS_JWS1985` | PES | 1 | **P/G** |
| `OH3` | `OH3_LEPS_SE_1986` | PES | 1 | **P/G** |
| `OH3` | `OH3_PIP_FFW1_2019` | PES | 3 | **P/G/D** |
| `OH3` | `OH3_PIP_FFW1_2019_DPEM` | DPEM | 3 | **U/UG** |
| `OH3` | `OH3_PIP_FFW2_2022` | PES | 3 | **P/G/D** |
| `OH3` | `OH3_PIP_FFW2_2022_DPEM` | DPEM | 3 | **U/UG** |
| `SiH4Cl` | `SiH4Cl_VBMM` | PES | 1 | **P/G** |
| `SiH5` | `SiH5_LEP_1998v1` | PES | 1 | **P/G** |
| `SiH5` | `SiH5_LEP_1998v2` | PES | 1 | **P/G** |
| `SiOH` | `SiOH_GEN_1988` | PES | 1 | **P/G** |
| `YRH` | `YRH_LEPS_ModelSurface` | PES | 2 | **P/G/D** |
| `YRH` | `YRH_LEPS_ModelSurface_DPEM` | DPEM | 2 | **U/UG** |

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

Legacy POTLIB setup adapters validate their declared atomic numbers through
this module. PES sources must not carry private copies of the periodic table;
fit-specific atom ordering and isotope choices remain with the owning surface.

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

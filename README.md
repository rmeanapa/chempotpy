
ChemPotPy, CHEMical library of POTential energy surfaces in PYthon 
==================================================================

Jan. 17, 2025

Authors: Yinan Shu, Zoltan Varga, Dayou Zhang, Donald G. Truhlar
University of Minnesota, Minnesota, United States

ChemPotPy is a library for analytic representation of single-state 
and multi-state potential energy surfaces and couplings. 

All fortran source code are stored in folder chempotpy 


How to install
--------------
* Create a conda virtual environment with gfortran and MKL:
    
       conda create --name chempotpy
       conda activate chempotpy
       conda install python=3.11
       conda install mkl mkl-service
       conda install -c conda-forge gfortran
       pip install numpy "numpy>=1.26,<1.27"
       pip install charset_normalizer

* Install stable release:
  
        pip install chempotpy



* The users may re-compile all .so modules for compatability reasons:

  get into the parent directory of chempotpy 
   
        make all 
        make check


Compile Chempotpy subroutine
----------------------------
One can use the meta programming script to automatically generate a 
fortran subroutine. 

  create a conda vitural environment as suggested in How to install

  get into the parent directory of chempotpy/chempotpy

       ./meta_chempotpy.script

The meta program will generate a fortran subroutine called chempotpy. 
One can interface this chempotpy subroutine with any dynamics code. 
Notice that the meta program will also generate sub programs for each 
surface that is located in chempotpy/chempotpy/system/lib/.


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

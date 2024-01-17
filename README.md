# FPVFoam and FPVBuoyantFoam v1912
===================================

## Description

This code realizes a steady laminar flamelet approach for turbulent non-premixed combustion.
Both solvers are based on ''rhoReactingFoam'' and ''rhoReactingBuoyantFoam'' respectively, i.e. it is pressure-based (PISO), compressible and runs with both LES and RAS turbulence.
 
The theory is mainly taken from the work of N. Peters and is based on the view of a turbulent flame as an ensemble of laminar flamelets.
The calculation of these flamelets is a one-dimensional problem and can be done in a pre-processing step.
Integration using a presumed beta-Probability Density Function (PDF) accounts for the interaction between turbulent fluctuations and chemistry.
The results of the pre-processing procedure are stored in tables which are accessed during the simulation.
Values of interest, such as species mass fraction or enthalpy, are looked-up and connected to the flow using three parameters - the mixture fraction, its variance and the progress variable.
In doing so, the expensive solution of chemical mechanisms during run-time can be avoided and the run-time thus reduces significantly.

## Available Solvers and Functions

* canteraToFPVFoamv1912 (Flamelet Library Integration)

* FPVFoamv1912 (Flow solving solver without buoyancy)

* FPVBuoyantFoamv1912 (Flow solving solver with buoyancy)

* FPVFoamPostv1912 (Post-processing result by outputing species mass fraction etc)

## Installation

This version works with OpenFOAM-v1912

* Prepare a directory on your system, e.g.:  

  `mkdir ~/OpenFOAM/FPVFoam-v1912/`

* Download FPVFoam-v1912 using git:

  `git clone https://github.com/weixian001/FPVFoam-v1912.git ~/OpenFOAM/FPVFoam-v1912/`

* Set an environment variable to the FPVFoam-v1912 src folder:

  `export LIB_FPVFoamV1912_SRC=~/OpenFOAM/FPVFoam-v1912/src/`

* Execute `./Allwmake`

## Citation

If you use our codes please cite:

```
@article{lim2023evaluation,
    author = {Lim, Wei Xian and Chan, Wai Lee and Elhadidi, Basman},
    title = "{Evaluation of Thermoacoustic Instability for Chemically Reacting Flows Using Large-Eddy Simulations}",
    journal = {Journal of Fluids Engineering},
    pages = {1-38},
    year = {2023},
    month = {12},
    issn = {0098-2202},
    doi = {10.1115/1.4064385},
    url = {https://doi.org/10.1115/1.4064385},
    eprint = {https://asmedigitalcollection.asme.org/fluidsengineering/article-pdf/doi/10.1115/1.4064385/7226813/fe-23-1410.pdf},
}

```

## Notes
This solver is based on the work done by Prof Pfitzner, FlameletFoam created at the Universität der Bundeswehr München, Thermodynamics Institute (Prof. Pfitzner). http://sourceforge.net/projects/openfoam-extend/files/OpenFOAM_Workshops/OFW8_2013_Jeju/Fri/Track3/HagenMuller-OFW8.tar/download

The FPVFoam and FPVBuoyantFoam are developed by the NanyangCFD team at Nanyang Technological University, Singapore leading by Prof Wai Lee CHAN (chan.wl@ntu.edu.sg). Main contributor is Wei Xian LIM (weixian001@e.ntu.edu.sg).


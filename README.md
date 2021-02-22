# poly_dpd

A Python library for building and running simple dissipative particle dynamics (DPD) simulations using HOOMD-blue. 

    Built using python3.6.9

## Software prerequisites:
- Python3 (Tested on version 3.6.9)
- HOOMD-blue (http://glotzerlab.engin.umich.edu/hoomd-blue/)
- numpy
- fresnel
- pillow

## Description
This library is used to automate DPD simulations of linear polymers using HOOMD-blue. It is intended to be used with Google Colab Notebooks to teach basic concepts of polymer physics to undergraduate students. The library is designed to perform three main functions:
- Build a linear polymer of N monomers and run a DPD simulation in HOOMD-blue
- Take monomer coordinates from the HOOMD-blue simulation and calculate the radius of gyration
- Render the 3-dimensional conformation of the polymer chain

# polysim

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

## Examples

The function <code>simulate_polymer_dpd</code> is used to set-up and run a DPD simulation of a linear polymer in HOOMD-blue. The function returns the monomer coordinates and save a trajectory file. Users must enter the number of monomers ($N$) and the conservative force repulsive interaction strength ($\alpha$). Default parameters such as the conservative force constant $\gamma$ (<code>gamma</code>), the diameter of DPD bead (<code>sigma</code>), the MD timestep (<code>dt</code>), the number of timesteps (<code>steps</code>), the energy of the system (<code>kT</code>), the trajectory filename, (<code>trajectory_filename</code>), and the number of frames to save in the trajectory file ((<code>traj_frames</code>)

    from polysim.main import simulate_polymer_dpd

    ## Enter number of monomers
    N = 10;

    ## Enter alpha, the repulsive force constant in the DPD conservative force
    alpha = 0.1;
    
    ## Enter the number of timesteps
    steps = 1e5

    ## Run simulation: Save the XYZ coordinates of the monomers
    coordinates  = simulate_polymer_dpd(N, alpha, steps = steps);

    ## View monomer XYZ coordiantes
    print(coordinates)
    
Once the polymer has been simulated, you can visualize the conformation using the <code>render_polymer</code> function. This function takes the monomer coordinates output (i.e. <code>coordinates</code> variable) from the <code>simulate_polymer_dpd</code> and renders a 3D visualization of the polymer using fresnel. 

    from polysim.main import render_polymer

    render_polymer(coordinates)
    
<img src="https://github.com/ericbruckner/polysim/blob/main/examples/Colab_lesson_plan/sample_data/polymer_render.png" height="300">

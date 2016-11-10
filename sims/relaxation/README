Relaxation of magnetic nanotracks made of hexagonally arranged atoms, using
discrete spin simulations.  The simulations are started calling the Python
script *hexagonal_fidimag.py* together with different parameters, i.e. by doing

    python ./hexagonal_fidimag.py -parameters

If we use the help with the *-h* parameter, the script shows all options
available

    usage: hexagonal_fidimag.py [-h] [--D D] [--J J] [--mu_s MU_S] [--k_u K_U]
                            [--B B] [--Demag] [--FullDemag] [--PBC_2D]
                            (--initial_state_skyrmion_down SK_INITIAL_RADIUS | 
                             --initial_state_skyrmion_up SK_INITIAL_RADIUS |
                             --initial_state_ferromagnetic_up | 
                             --initial_state_ferromagnetic_down | 
                             --initial_state_irregular)
                            [--preview] [--unit_length UNIT_LENGTH]
                            [--alpha ALPHA] [--save_files SAVE_FILES]
                            [--stopping_dmdt STOPPING_DMDT] [--tols TOLS TOLS]
                            [--max_steps MAX_STEPS] [--no_precession]
                            [--alignment ALIGNMENT]
                            hex_nx hex_ny hex_a sim_name

    NEB method for 2D films with interfacial DMI


    positional arguments:
      hex_nx                Number of elements along x
      hex_ny                Number of elements along y
      hex_a                 Lattice constant
      sim_name              Simulation name

    optional arguments:
      -h, --help            show this help message and exit
      --alignment ALIGNMENT
                            Hexagonal lattice alignment (square or diagonal)

      --D D                 DMI constant in units of meV
      --J J                 Exchange constant in units of meV
      --mu_s MU_S           Magnetisation in units of mu_B
      --k_u K_U             Anisotropy constant in units of meV
      --B B                 External magnetic field perpendicular to the square
                            plane (z direction), in Tesla
      --Demag               Add this option to use dipolar interactions using the
                            FFT technique (TESTING)
      --FullDemag           Add this option to use dipolar interactions using a
                            direct calculation (TESTING)
      --PBC_2D              Two dimensional boundary condition
      --initial_state_skyrmion_down SK_INITIAL_RADIUS
                            This option puts a skyrmionic texture in the centre of
                            the nanotrack, as initial m configuration. The other
                            spins are in the (0, 0, 1) direction
      --initial_state_skyrmion_up SK_INITIAL_RADIUS
                            This option puts a skyrmionic texture with its core
                            pointing in the +z direction, in the centre of the
                            nanotrack, as initial m configuration. The other spins
                            are in the (0, 0, 1) direction
      --initial_state_ferromagnetic_up
                            This option sets the initial m configuration as a
                            ferromagnetic state in the (0, 0, 1) direction
      --initial_state_ferromagnetic_down
                            This option sets the initial m configuration as a
                            ferromagnetic state in the (0, 0, -1) direction
      --initial_state_irregular
                            This option sets the initial m configuration as an
                            irregular state (TESTING)
      --preview             Specify if instead of relaxing the system, it will be
                            shown a real time plot of the magnetisation dynamics
                            on the TOP layer (in the z direction). This will run
                            for 4 nanoseconds
      --unit_length UNIT_LENGTH
                            Mesh unit length
      --alpha ALPHA         Damping constant value
      --save_files SAVE_FILES
                            Save vtk and npy files every x steps
      --stopping_dmdt STOPPING_DMDT
                            Specify an specific dm/dt threshold value when
                            relaxing the system (default is 0.01)
      --tols TOLS TOLS      Tolerances for the integrator as a pair: rtol atol
      --max_steps MAX_STEPS
                            Specify maximum number of steps for the relaxation
                            (default is 5000)
      --no_precession       To remove LLG precesion term

Most of the simulations are started with a bash script, where the options are
recorded. For example, if we have *my_simulation.sh*, we only do

    bash my_simulation.sh

in a terminal. The bash scripts are located in this folder and the output files
are saved in ordered folders, such as ``vtks``, ``npys``, etc. The scripts are
self descriptive by their names.

For Cobalt tracks, the DMI is negative, according to our definition of the
Dzyaloshinskii vectors.

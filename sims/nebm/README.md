NEB method simulation scripts for atomistic 2D systems.
the base script is ``hexagonal_neb_fidimag.py`` which we call passing
arguments with the simulation details as:

    python hexagonal_neb_fidimag.py -parameters

If we call the ``--help`` argument, we get all the available options and a
description of the script:

    usage: hexagonal_neb_fidimag.py [-h] [--unit_length UNIT_LENGTH]
                                    [--alpha ALPHA] [--PBC_2D] [--neb_k NEB_K]
                                    [--stopping_dmdt STOPPING_DMDT]
                                    [--tols TOLS TOLS] [--coordinates COORDINATES]
                                    [--D D] [--J J] [--mu_s MU_S] [--k_u K_U]
                                    [--B B] [--Demag] [--alignment ALIGNMENT]
                                    [--climbing_image INDEX]
                                    (--interpolation INTERPOLATION [INTERPOLATION ...] | 
                                     --images_files NPY_FILES_PATH)
                                    hex_nx hex_ny hex_a neb_steps save_vtk
                                    save_npy sim_name

    NEB method for 2D films with interfacial DMI

    positional arguments:
      hex_nx                Number of elements along x
      hex_ny                Number of elements along y
      hex_a                 Lattice constant
      neb_steps             Number of steps for the Geodesic or Cartesian or
                            Spherical NEB method
      save_vtk              Number specifying that the vtk files are going to be
                            saved every *save_vtk* number of steps
      save_npy              Number specifying that the npy files are going to be
                            saved every *save_vtk* number of steps
      sim_name              Simulation name

    optional arguments:
      -h, --help            show this help message and exit
      --unit_length UNIT_LENGTH
                            Mesh unit length
      --alpha ALPHA         Damping constant value
      --PBC_2D              Two dimensional boundary condition
      --neb_k NEB_K         Spring constant for the NEB method. Default: k=1e10
      --stopping_dmdt STOPPING_DMDT
                            Stopping dmdt for the NEB method. Default: 1e-2
      --tols TOLS TOLS      Tolerances for the integrator as a pair: rtol atol
      --coordinates COORDINATES
                            Cartesian or Spherical or Geodesic
      --D D                 DMI constant in units of meV
      --J J                 Exchange constant in units of meV
      --mu_s MU_S           Magnetisation in units of mu_B
      --k_u K_U             Anisotropy constant in units of meV
      --B B                 External magnetic field perpendicular to the square
                            plane (z direction), in Tesla
      --Demag               Add this option to use dipolar interactions
      --alignment ALIGNMENT
                            Hexagonal lattice alignment (square or diagonal)
      --climbing_image INDEX
                            An optional integer ranging from 1 to the total number
                            of images minus two (the extremes do not count) in
                            order to use that image of the band with the Climbing
                            Image NEB method. This image must be the maximum
                            energy state and will not be affected by the Spring
                            Force. This helps to have a better resolution around
                            the saddle point. Currently, only working with
                            CARTESIAN coordinates
      --interpolation INTERPOLATION [INTERPOLATION ...]
                            Use a series of images for the initial state,
                            interpolating the magnetisation vectors components
                            from the first argument until the last argument.The
                            parameters are: file1 number_interps file2
                            number_interps file3 ...Where number_interps is the
                            number of interpolationsbetween the i and (i+1)th
                            states of the Energy band. Thus the number of
                            arguments is always odd. If file_x is passed as a path
                            to a folder, it is assumed that there is only a single
                            file in that folder.
      --images_files NPY_FILES_PATH
                            Use all the *npy* files from the specified folder path
                            for the initial states of the Energy Band. File names
                            must be sorted as *image_().npy* with () as the number
                            in the sequence of images in the Energy Band, and
                            include the initial and final states. If the path to
                            the folder finishes in _LAST, we search the largest
                            number among the *folders* which finish in _{step},
                            where step is a number


The NEBM simulations for Cobalt nanotracks are provided as bash scripts in the
respective folders on this path. These scripts contain details of every system.
The idea is that, for every simulation, we call a bash script as

    bash initiate_simulation.sh

In every folder, the bash files relax a system with the NEBM as:

    initiate_simulation_k1e4.sh             :: Relaxation using linear interpolations
    initiate_simulation_climbing-{}_k1e4.sh :: Relaxation using the climbing image and the
                                               last step of the initiate_simulation_k1e4.sh
                                               script

    initiate_simulation_skdisp_k1e4.sh      :: Relaxation for the skyrmion annihilation
                                               at the boundary. These simulations require
                                               to generate the initial states for the band
                                               from the create_sk_on_track.py script.  


Furthermore, to intiate these simulations, we need the stable states, which are
the skyrmion and ferromagnetic state, for every DMI value. Thus, we have to run
the scripts in the ``relaxation`` folder before calling these NEBM scripts.

The D values on the simulation names refer to the DMI magnitude in ``J/m^-2``,
but since we are simulating atomistic systems, these magnitudes are equivalent
to the spin discrete values specified in the paper.

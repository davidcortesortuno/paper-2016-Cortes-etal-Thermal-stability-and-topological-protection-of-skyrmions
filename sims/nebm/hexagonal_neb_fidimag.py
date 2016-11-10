from __future__ import print_function

import argparse
import sys

# ARGUMENTS
parser = argparse.ArgumentParser(description='NEB method for 2D films with '
                                 'interfacial DMI')
method = parser.add_mutually_exclusive_group(required=True)

# Material and geometry -------------------------------------------------------
parser.add_argument('hex_nx', help='Number of elements along x',
                    type=int)

parser.add_argument('hex_ny', help='Number of elements along y',
                    type=int)

parser.add_argument('hex_a', help='Lattice constant',
                    type=float)

# NEB Simulation --------------------------------------------------------------

parser.add_argument('neb_steps', help='Number of steps for the '
                    'Geodesic or Cartesian or Spherical NEB method')

parser.add_argument('save_vtk', help='Number specifying that the vtk files'
                    ' are going to be saved every *save_vtk* number of steps',
                    type=int)

parser.add_argument('save_npy', help='Number specifying that the npy files '
                    'are going to be saved every *save_vtk* number of steps',
                    type=int)

parser.add_argument('sim_name',
                    help='Simulation name')

parser.add_argument('--unit_length', help='Mesh unit length',
                    type=float, default=1.0)

parser.add_argument('--alpha', help='Damping constant value',
                    type=float)

parser.add_argument('--PBC_2D',
                    help='Two dimensional boundary condition',
                    action='store_true')

parser.add_argument('--neb_k', help='Spring constant for the NEB method'
                    '. Default: k=1e10',
                    default='1e10', type=float)

parser.add_argument('--stopping_dmdt', help='Stopping dmdt for the NEB method'
                    '. Default: 1e-2',
                    default='1e-2', type=float)

parser.add_argument('--tols', help='Tolerances for the integrator as '
                    'a pair: rtol atol',
                    type=float, nargs=2)

parser.add_argument('--coordinates',
                    help='Cartesian or Spherical or Geodesic',
                    default='Cartesian')

# Material  -------------------------------------------------------------------

parser.add_argument('--D', help='DMI constant in units of meV',
                    type=float, default=1)

parser.add_argument('--J', help='Exchange constant in units of meV',
                    type=float, default=1)

parser.add_argument('--mu_s', help='Magnetisation in units of mu_B',
                    type=float, default=2)

parser.add_argument('--k_u', help='Anisotropy constant in units of meV',
                    type=float)

parser.add_argument('--B', help='External magnetic field perpendicular to the'
                    ' square plane (z direction), in Tesla',
                    type=float)

parser.add_argument('--Demag', help='Add this option to use dipolar '
                    'interactions',
                    action='store_true')

parser.add_argument('--alignment', help='Hexagonal lattice alignment '
                    '(square or diagonal)',
                    default='diagonal')

# NEB Method ------------------------------------------------------------------

parser.add_argument('--climbing_image',
                    help='An optional integer ranging from 1 to the total '
                    ' number of images minus two (the extremes do not count)'
                    ' in order to use that image of the band with the'
                    ' Climbing Image NEB method. This image must be the'
                    ' maximum energy state and will not be affected by'
                    ' the Spring Force. This helps to have a better'
                    ' resolution around the saddle point. Currently, only'
                    ' working with CARTESIAN coordinates',
                    metavar=('INDEX'), type=int)

# Define a method: interpolation or images_files for the initial state
# which are mutually exclusive
method.add_argument('--interpolation',
                    help='Use a series of images for the initial state, '
                    'interpolating the magnetisation vectors components '
                    'from the first argument until the last argument.'
                    'The parameters are: '
                    'file1 number_interps file2 number_interps file3 ...'
                    'Where number_interps is the number of interpolations'
                    'between the i and (i+1)th states of the Energy band.'
                    ' Thus the number of arguments is always odd.',
                    nargs='+',
                    # metavar=('FILE1', 'N_INTERPS', 'FILE2', 'etc')
                    )

method.add_argument('--images_files',
                    help='Use all the *npy* files from the specified folder '
                    'path for the initial states of the Energy Band. File '
                    'names must be sorted as *image_().npy* with () as the '
                    'number in the sequence of images in the Energy Band, '
                    'and include the initial and final states',
                    metavar=('NPY_FILES_PATH'))

# -----------------------------------------------------------------------------


# Parser arguments
args = parser.parse_args()


import numpy as np

from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import UniformExchange, DMI, Anisotropy, DemagHexagonal
from fidimag.atomistic import Zeeman

# Import physical constants from fidimag
import fidimag.common.constant as const

from fidimag.common.nebm_geodesic import NEBM_Geodesic as NEBM_Geodesic
from fidimag.common.nebm_geodesic import NEBM_Cartesian as NEBM_Cartesian
from fidimag.common.nebm_geodesic import NEBM_Spherical as NEBM_Sherical

import os
import re

mu0 = 4 * np.pi * 1e-7


# -----------------------------------------------------------------------------
# Mesh ------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if not args.PBC_2D:
    mesh = HexagonalMesh(args.hex_a, args.hex_nx, args.hex_ny,
                         alignment=args.alignment,
                         unit_length=args.unit_length
                         )
else:
    mesh = HexagonalMesh(args.hex_a, args.hex_nx, args.hex_ny,
                         periodicity=(True, True),
                         alignment=args.alignment,
                         unit_length=args.unit_length
                         )


# Initiate Fidimag simulation -------------------------------------------------
sim = Sim(mesh, name=args.sim_name)

# sim.set_tols(rtol=1e-10, atol=1e-14)
if args.alpha:
    sim.alpha = args.alpha
# sim.gamma = 2.211e5

# Material parameters ---------------------------------------------------------

sim.mu_s = args.mu_s * const.mu_B

exch = UniformExchange(args.J * const.meV)
sim.add(exch)

dmi = DMI(D=(args.D * const.meV), dmi_type='interfacial')
sim.add(dmi)

if args.B:
    zeeman = Zeeman((0, 0, args.B))
    sim.add(zeeman, save_field=True)

if args.k_u:
    # Uniaxial anisotropy along + z-axis
    sim.add(Anisotropy(args.k_u * const.meV, axis=[0, 0, 1]))

if args.Demag:
    sim.add(DemagHexagonal())

# -----------------------------------------------------------------------------


# Debug Information -----------------------------------------------------------

print('Saturation Magnetisation: {} mu_B'.format(args.mu_s))
print('Exchange constant: {}  meV'.format(args.J))
print('DMI constant: {}  meV'.format(args.D))
if args.k_u:
    print('Anisotropy constant: {}   meV'.format(args.k_u))
if args.B:
    print('Zeeman field: (0, 0, {})  T'.format(args.B / mu0))
print('--------------------------------------')

# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Initiate simulation ---------------------------------------------------------
# -----------------------------------------------------------------------------


# ================= SIMULATION FUNCTION ======================================

# Define a function to vary the number of initial images and spring constants
"""
Execute a simulation with the NEB function of the fidimag code, for a
stripe (or track or rectangle).

WARNING: interpolation is only well defined in SPHERICAL coordinates,
         thus it is recommended to start the system with the SNEBM
         and later continue the iteration with the CNEBM

vtks and npys are saved in files starting with the 'args.sim_name'
string (only initial and final states)

The frequency of how vtu or npy files are saved is optional

"""

# SIMULATION  ----------------------------------------------------------------

# Interpolations and initial images ------------------------------------------

if args.interpolation:
    # Load the states from every 2 arguments of
    # the --interpolation option and the numbers from every  two arguments
    # starting from the 1st element of the list
    images = [np.load(state) for state in args.interpolation[::2]]
    interpolations = [int(num_interps)
                      for num_interps in args.interpolation[1::2]]

    print('==================== Interpolating!  ====================== \n')

elif args.images_files:
    # Load the states from the specified npys folder in --images_files
    # We will sort the files using the number in their names
    # assuming that they have the structure:
    # 'image_1.npy', 'image_13.npy', etc.
    # For this, we search until two digits with regex and
    # return the integer
    images = [np.load(os.path.join(args.images_files, _file))
              for _file in
              sorted(os.listdir(args.images_files),
                     key=lambda f: int(re.search('\d+', f).group(0))
                     )
              ]

    print('FILES:')
    for _file in sorted(os.listdir(args.images_files),
                        key=lambda f: int(re.search('\d+', f).group(0))
                        ):
        print(_file)

    # No interpolations in this case
    interpolations = None
else:
    print('Specify an accepted method for the initial states')

# ----------------------------------------------------------------------------

# In Fidimag when using Godesic/Cartesian coordinates, the interpolations are
# performed in Spherical coordinates

if args.coordinates == 'Cartesian':
    neb_method = NEBM_Cartesian
if args.coordinates == 'Spherical':
    neb_method = NEBM_Spherical
if args.coordinates == 'Geodesic':
    neb_method = NEBM_Geodesic

neb = neb_method(sim,
                 images,
                 interpolations=interpolations,
                 spring_constant=args.neb_k,
                 name=args.sim_name,
                 climbing_image=args.climbing_image,
                 openmp=True
                 )

if args.tols:
    neb.create_integrator()
    neb.set_tols(rtol=args.tols[0], atol=args.tols[1])

# Relax the system with the NEB mehtod
neb.relax(max_iterations=int(args.neb_steps),
          save_vtks_every=args.save_vtk,
          save_npys_every=args.save_npy,
          stopping_dYdt=args.stopping_dmdt)

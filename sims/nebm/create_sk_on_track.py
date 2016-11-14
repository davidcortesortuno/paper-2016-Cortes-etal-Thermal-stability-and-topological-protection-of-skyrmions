"""

Script to create a sequence of 19 NPY files for the magnetisation field, where
a skyrmion is moved towards the top boundary of the nanotrack. This sequence is
then used as initial state for the NEBM to get the annihilation of a skyrmion
mediated by the sample boundary

It is necessary to pass the path to the NPYs folders for the skyrmion and
ferromagnetic state. In our simulations, these folders are in the relaxation
folder:

    relaxation/npys/stripe_320x185_Co_sk-down_D{}_h25e-2nm_npys
    relaxation/npys/stripe_320x185_Co_fm-up_D{}_h25e-2nm_npys

In addition, it is necessary to provide a radius for a circle that encloses an
area around the skyrmion. This area is moved towards the top boundary until the
skyrmion dissapears. Thus, if the skyrmion is 16 nm wide we can call this
script by doing:

    python \
    --sk_folder ../relaxation/npys/stripe_320x185_Co_sk-down_D{}_h25e-2nm_npys \
    --sk_folder ../relaxation/npys/stripe_320x185_Co_fm-up_D{}_h25e-2nm_npys \
    --radius 12

which creates a folder called sk_displ_npys where the NPY files are stored.


"""


def check_slash(folder):
    if not folder.endswith('/'):
        return (folder + '/')
    else:
        return folder

# -----------------------------------------------------------------------------

import argparse

parser = argparse.ArgumentParser(description='Create sequence for the NEBM')

parser.add_argument('--sk_folder', help='Path to the skyrmion folder',
                    required=True)
parser.add_argument('--fm_folder', help='Path to the ferromagnetic folder',
                    required=True)
parser.add_argument('--radius', help='Radius of a circle that encloses the '
                    ' skyrmion', type=float, required=True)

# Parser arguments
parser.set_defaults(surface=True)
args = parser.parse_args()

# -----------------------------------------------------------------------------

import numpy as np

from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.common.save_vtk import SaveVTK
import scipy.interpolate

import fidimag.common.constant as const

import os
import shutil
import sys
import glob


args.sk_folder = check_slash(args.sk_folder)
args.fm_folder = check_slash(args.fm_folder)


# Get the skyrmion and ferromagnetic state (we assume there is only a single
# file per folder)
sk_npy_path = glob.glob(args.sk_folder + '*.npy')[0]
fm_npy_path = glob.glob(args.fm_folder + '*.npy')[0]

# Number of images to be generated
number_images = 19

# -----------------------------------------------------------------------------

mu_s = 0.846 * const.mu_B

a = 0.125
nx = 320
ny = 185

# MESH
mesh = HexagonalMesh(a, nx, ny,
                     alignment='square',
                     unit_length=1e-9
                     )

vtk_saver = SaveVTK(mesh, name='sk_displ')

centre_x = (np.max(mesh.coordinates[:, 0])
            + np.min(mesh.coordinates[:, 0])) * 0.5

centre_y = (np.max(mesh.coordinates[:, 1])
            + np.min(mesh.coordinates[:, 1])) * 0.5

lx = (np.max(mesh.coordinates[:, 0])
      - np.min(mesh.coordinates[:, 0]))
ly = (np.max(mesh.coordinates[:, 1])
      - np.min(mesh.coordinates[:, 1]))

# Create a simulation for the FM state and one for the skyrmion
sim = Sim(mesh, name='sk_displ')
sim.mu_s = mu_s
sim_sk = Sim(mesh)
sim_sk.mu_s = mu_s


# Manipulate the initial state
# Generate the skyrmion in the middle of the nanotrack
def generate_skyrmion(pos, i,
                      x_c, y_c,
                      radius,
                      vector_field,
                      sk_field
                      ):
    """

    Returns a spin direction for the i-th lattice site, which has coordinates
    *pos*

    Basically:

    1. We take a vector field as base (generally, the ferromagnetic state) and
    the skyrmion magnetisation field from sk_field.

    2. Using the *i* variable, we are analysing the i-th lattice site
    (according to Fidimag index numbering)

    3. From the position *pos* we pass  to this function, we check if it is
    inside a circle centered at (x_c, y_c) of radius *radius*.

    4. If so, we set the spin direction from the skyrmion defined in sk_field
    (we assume this skyrmion is centered at (0, 0)) and return this value

    5. If the point is not inside the circle, we return the spin direction from
    the i-th position of the vector field (ferromagnetic field)

    Variables:

    vector_field :: An array with the spin directions of the fm field
    sk_field     :: A function (interpolation) of the skyrmion m field

    """

    # Relative to the SK - coordinate system
    xrel, yrel = (pos[0] - x_c), (pos[1] - y_c)
    x, y = pos[0], pos[1]

    if np.sqrt(xrel ** 2. + yrel ** 2.) <= radius:
        return [sk_field[0](centre_x + xrel, centre_y + yrel),
                sk_field[1](centre_x + xrel, centre_y + yrel),
                sk_field[2](centre_x + xrel, centre_y + yrel)]

    else:
        return vector_field[i]

# Now put the SKyrmion in the x_c, y_c position
# We read through the coordinates and modify the vector field where the
# Skyrmion is positioned

# We need a radius sufficiently large to contain the skyrmion from the npy file
for folder in ['sk_displ_npys', 'sk_displ_vtks']:
    if os.path.exists(folder):
        shutil.rmtree(folder)

for folder in ['sk_displ_npys']:
    if not os.path.exists(folder):
        os.makedirs(folder)

# Load the SK profile (This will not be modified)
sim_sk.set_m(np.load(sk_npy_path))

# Create the skyrmion magnetisation vector field function, interpolating
# the data from the state with the skyrmion at the centre
m = sim_sk.spin.reshape(-1, 3)
sk_field = []
for i in range(3):
    sk_field.append(scipy.interpolate.NearestNDInterpolator(
        mesh.coordinates[:, :2],
        m[:, i]
        ))

i = 1
for sk_centre_y in np.linspace(centre_y + 5, ly + args.radius, number_images):

    # Load the ferromagnetic state for these parameters
    sim.set_m(np.load(fm_npy_path))

    # We will append the spin directions in the same order than in the
    # arrays from the sim simulation
    new_m_field = []

    # Generate the skyrmion going through every lattice site and appending
    # spin directions according to the *generate_skyrmion* function
    for j, x_vector in enumerate(mesh.coordinates):
            new_m_field.append(generate_skyrmion(x_vector, j,
                                                 centre_x, sk_centre_y,
                                                 args.radius,
                                                 sim.spin.reshape(-1, 3),
                                                 sk_field
                                                 ))
    # In xyz format
    new_m_field = np.array(new_m_field).reshape(-1,)

    # Set the magnetisation
    sim.set_m(new_m_field)

    np.save('sk_displ_npys/image_{}.npy'.format(i), sim.spin)
    # sim.save_vtk(vtkname='m_{}'.format(i))
    vtk_saver.save_vtk(sim.spin.reshape(-1, 3),
                       sim._mu_s,
                       step=0, vtkname='image_{}'.format(i))
    shutil.move('sk_displ_vtks/image_{}_000000.vtk'.format(i),
                'sk_displ_vtks/image_' + str(i).zfill(3) + '.vtk'
                )

    i += 1


# Save the first and last states

# Load the ferromagnetic state for these parameters
sim.set_m(np.load(fm_npy_path))

np.save('sk_displ_npys/image_20.npy', sim.spin)
np.save('sk_displ_npys/image_0.npy', sim_sk.spin)

vtk_saver.save_vtk(sim_sk.spin.reshape(-1, 3),
                   sim._mu_s,
                   step=0, vtkname='image_0')
shutil.move('sk_displ_vtks/image_0_000000.vtk',
            'sk_displ_vtks/image_' + str(0).zfill(3) + '.vtk'
            )

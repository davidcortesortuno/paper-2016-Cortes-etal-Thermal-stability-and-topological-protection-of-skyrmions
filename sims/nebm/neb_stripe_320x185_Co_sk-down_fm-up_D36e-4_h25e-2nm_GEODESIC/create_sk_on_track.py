"""

"""

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

sk_npy_path = '../../relaxation_fidimag_atomistic/npys/stripe_320x185_Co_sk-down_D36e-3_h25e-2nm_npys/m_504.npy'
fm_npy_path = '../../relaxation_fidimag_atomistic/npys/stripe_320x185_Co_fm-up_D36e-3_h25e-2nm_npys/m_381.npy'

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
                      sk_initial_radius,
                      vector_field,
                      sk_field
                      ):
    """

    [EXPERIMENTAL]
    We will generate a skyrmion with its centre positioned at (x_c, y_c),
    whose core is pointing down. We will use the skyrmion from a
    track with a skyrmion npy file, rather than generate the
    configuration manually. For this, we will copy all the spins in a radius
    centered at (0, 0), of the npy file.

    The rest of the spins will
    keep their original orientation, given by the 'vector field'
    dolfin function

    x_c, y_c  are the coordinates of the skyrmion centre
              to be generated

    sk_initial_radius is the skyrmion radius
    """
    # Relative to the SK - coordinate system
    xrel, yrel = (pos[0] - x_c), (pos[1] - y_c)
    x, y = pos[0], pos[1]

    if np.sqrt(xrel ** 2. + yrel ** 2.) <= sk_initial_radius:
        # r = (x ** 2 + y ** 2) ** 0.5
        # phi = np.arctan2(y, x)
        # # This determines the profile we want for the
        # # skyrmion
        # # Single twisting: k = pi / R
        # k = np.pi / (sk_initial_radius)

        # # We define here a 'hedgehog' skyrmion pointing down
        # return [1 * np.sin(k * r) * np.cos(phi),
        #         1 * np.sin(k * r) * np.sin(phi),
        #         -np.cos(k * r)]
        return [sk_field[0](centre_x + xrel, centre_y + yrel),
                sk_field[1](centre_x + xrel, centre_y + yrel),
                sk_field[2](centre_x + xrel, centre_y + yrel)]

    else:
        return vector_field[i]

# Now put the SKyrmion in the x_c, y_c position
# We read through the coordinates and modify the vector
# field where the Skyrmion is positioned

# We will generate a series of states with skyrmions moving towards
# one edge of the nanotrack

# We need a radius sufficiently large to contain the skyrmion
# from the npy file. We can get an estimate of the sk radius using the
# track_plot_mz.py script

sk_radius = 18.5

for folder in ['sk_displ_npys', 'sk_displ_vtks']:
    if os.path.exists(folder):
        shutil.rmtree(folder)

for folder in ['sk_displ_npys']:
    if not os.path.exists(folder):
        os.makedirs(folder)

# Load the SK profile (This will not be modified)
sim_sk.set_m(np.load(sk_npy_path)
             )

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
for sk_centre_y in np.linspace(centre_y + 5, ly + sk_radius, 19):

    # Load the ferromagnetic state for these parameters
    sim.set_m(np.load(fm_npy_path)
              )

    new_m_field = []
    # Generate the skyrmion
    for j, x_vector in enumerate(mesh.coordinates):
            new_m_field.append(generate_skyrmion(x_vector, j,
                                                 centre_x, sk_centre_y,
                                                 sk_radius,
                                                 sim.spin.reshape(-1, 3),
                                                 # sim_sk.spin.reshape(-1, 3)
                                                 sk_field
                                                 ))
    # In xyz format
    new_m_field = np.array(new_m_field).reshape(-1,)

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
sim.set_m(np.load(fm_npy_path)
          )

np.save('sk_displ_npys/image_20.npy', sim.spin)
np.save('sk_displ_npys/image_0.npy', sim_sk.spin)

vtk_saver.save_vtk(sim_sk.spin.reshape(-1, 3),
                   sim._mu_s,
                   step=0, vtkname='image_0')
shutil.move('sk_displ_vtks/image_0_000000.vtk',
            'sk_displ_vtks/image_' + str(0).zfill(3) + '.vtk'
            )

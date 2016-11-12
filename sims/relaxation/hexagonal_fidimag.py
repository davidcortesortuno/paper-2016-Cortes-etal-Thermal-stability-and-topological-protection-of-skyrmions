from __future__ import print_function
"""

Script to run Fidimag simulations using Bash or IPython notebooks
with an optional Real Time View of the dynamics!

These simulations are based on hexagonal lattices.

Since we use argparse, a series of arguments must be specified.

If run from an IPython notebook, use the %run magic, since it
plays nicely with NBagg. For some reason, %%bash does not works with
NBAgg

This program requires Matplotlib >= 1.4.3 to use the Live Preview
in a notebook.

Author: David I. C.
email: d.i.cortes@soton.ac.uk

Modification date: Mon 02 Nov 2015 18:28:03 GMT

"""


# See if this library is called from an IPython notebook and
# return True or False accordingly
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

import sys
import argparse


# -----------------------------------------------------------------------------
# ARGUMENTS -------------------------------------------------------------------
# -----------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='NEB method for 2D films with '
                                 'interfacial DMI')
initial_state = parser.add_mutually_exclusive_group(required=True)

parser.add_argument('hex_nx', help='Number of elements along x',
                    type=int)

parser.add_argument('hex_ny', help='Number of elements along y',
                    type=int)

parser.add_argument('hex_a', help='Lattice constant',
                    type=float)

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
                    'interactions using the FFT technique (TESTING)',
                    action='store_true')

parser.add_argument('--FullDemag', help='Add this option to use dipolar '
                    'interactions using a direct calculation (TESTING)',
                    action='store_true')

parser.add_argument('sim_name',
                    help='Simulation name')

parser.add_argument('--PBC_2D',
                    help='Two dimensional boundary condition',
                    action='store_true')

initial_state.add_argument('--initial_state_skyrmion_down',
                           help='This option puts a skyrmionic texture'
                           ' in the centre of the'
                           ' nanotrack, as initial m configuration. The'
                           ' other spins are in the (0, 0, 1) direction',
                           type=float,
                           metavar=('SK_INITIAL_RADIUS')
                           # action='store_true'
                           )

initial_state.add_argument('--initial_state_skyrmion_up',
                           help='This option puts a skyrmionic texture'
                           ' with its core pointing in the +z direction, '
                           'in the centre of the'
                           ' nanotrack, as initial m configuration. The'
                           ' other spins are in the (0, 0, 1) direction',
                           type=float,
                           metavar=('SK_INITIAL_RADIUS')
                           # action='store_true'
                           )

initial_state.add_argument('--initial_state_ferromagnetic_up',
                           help='This option sets the initial '
                           'm configuration as a ferromagnetic state'
                           ' in the (0, 0, 1) direction',
                           action='store_true'
                           )

initial_state.add_argument('--initial_state_ferromagnetic_down',
                           help='This option sets the initial '
                           'm configuration as a ferromagnetic state'
                           ' in the (0, 0, -1) direction',
                           action='store_true'
                           )

initial_state.add_argument('--initial_state_irregular',
                           help='This option sets the initial '
                           'm configuration as an irregular state'
                           ' (TESTING)',
                           action='store_true')

parser.add_argument('--preview', help='Specify if instead of relaxing the '
                    'system, it will be shown a real time plot of the '
                    'magnetisation dynamics on the TOP layer (in the z '
                    'direction). This will run for 4 nanoseconds',
                    action='store_true'
                    )

parser.add_argument('--unit_length', help='Mesh unit length',
                    type=float, default=1.0)

parser.add_argument('--alpha', help='Damping constant value',
                    type=float, default=0.01)

parser.add_argument('--save_files', help='Save vtk and npy files every x'
                    ' steps',
                    type=float, default=None)

parser.add_argument('--stopping_dmdt', help='Specify an specific dm/dt '
                    'threshold value when relaxing the system (default '
                    ' is 0.01)',
                    type=float, default=0.01)

parser.add_argument('--tols', help='Tolerances for the integrator as '
                    'a pair: rtol atol',
                    type=float, nargs=2)

parser.add_argument('--max_steps', help='Specify maximum number of '
                    'steps for the relaxation (default is 5000)',
                    type=int, default=5000)

parser.add_argument('--no_precession', help='To remove LLG precesion term',
                    action='store_true')

parser.add_argument('--alignment', help='Hexagonal lattice alignment '
                    '(square or diagonal)',
                    default='diagonal')

# Parser arguments
args = parser.parse_args()

# -----------------------------------------------------------------------------

import matplotlib

if args.preview:
    if run_from_ipython():
        matplotlib.use('nbagg')
        print('Using Backend: NBAgg')
    else:
        matplotlib.use('TkAgg')
        print('Using Backend: TkAgg')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np

from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import UniformExchange, DMI, Anisotropy, DemagHexagonal, DemagFull
from fidimag.atomistic import Zeeman
# Physical constants
import fidimag.common.constant as const

import os
import shutil

mu0 = 4 * np.pi * 1e-7


# -------------------------------------------------------------------------
# Initial states ----------------------------------------------------------
# -------------------------------------------------------------------------

def generate_skyrmion_down(pos, sign):
    """
    Sign will affect the chirality of the skyrmion
    """
    # We will generate a skyrmion in the middle of the stripe
    # (the origin is there) with the core pointing down
    if args.alignment == 'diagonal':
        factor = 3 / 4.
    if args.alignment == 'square':
        factor = 0.5

    x = (pos[0] - sim.mesh.Lx * factor)
    y = (pos[1] - sim.mesh.Ly * 0.5)

    if np.sqrt(x ** 2 + y ** 2) <= args.initial_state_skyrmion_down:
        # Polar coordinates:
        r = (x ** 2 + y ** 2) ** 0.5
        phi = np.arctan2(y, x)
        # This determines the profile we want for the
        # skyrmion
        # Single twisting: k = pi / R
        k = np.pi / (args.initial_state_skyrmion_down)

        # We define here a 'hedgehog' skyrmion pointing down
        return (sign * np.sin(k * r) * np.cos(phi),
                sign * np.sin(k * r) * np.sin(phi),
                -np.cos(k * r))
    else:
        return (0, 0, 1)


def generate_skyrmion_up(pos, sign):
    # We will generate a skyrmion in the middle of the stripe
    # (the origin is there) with the core pointing down
    if args.alignment == 'diagonal':
        factor = 3 / 4.
    if args.alignment == 'square':
        factor = 0.5

    x = (pos[0] - sim.mesh.Lx * factor)
    y = (pos[1] - sim.mesh.Ly * 0.5)

    if np.sqrt(x ** 2 + y ** 2) <= args.initial_state_skyrmion_up:
        # Polar coordinates:
        r = (x ** 2 + y ** 2) ** 0.5
        phi = np.arctan2(y, x)
        # This determines the profile we want for the
        # skyrmion
        # Single twisting: k = pi / R
        k = np.pi / (args.initial_state_skyrmion_up)

        # We define here a 'hedgehog' skyrmion pointing down
        return (sign * np.sin(k * r) * np.cos(phi),
                sign * np.sin(k * r) * np.sin(phi),
                np.cos(k * r))
    else:
        return (0, 0, -1)


def irregular_state(pos):
    m = np.copy(sim.spin)
    m = m.reshape(3, -1).T
    # We will generate a skyrmion in the middle of the stripe
    # (the origin is there) with the core pointing down
    x = (pos[0] - args.hex_nx * args.hex_a * 0.5)
    y = (pos[1] - args.hex_ny * args.hex_a * (3 / 4.) * 0.5)

    if x > 0:
        return (0, 0, 1)
    else:
        return (0, 0, -1)

# --------------------------------------------------------------------------

# Mesh ---------------------------------------------------------------------

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

# Initiate Fidimag simulation ---------------------------------------------
sim = Sim(mesh, name=args.sim_name)

# sim.set_tols(rtol=1e-10, atol=1e-14)
sim.alpha = args.alpha
# sim.driver.gamma = 2.211e5

if args.no_precession:
    sim.driver.do_precession = False

# Material parameters -----------------------------------------------------

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
    print('Using Demag!')
    sim.add(DemagHexagonal())

if args.FullDemag:
    print('Using Full Demag calculation!')
    sim.add(DemagFull())

# -------------------------------------------------------------------------


# Load magnetisation profile ---------------------------------------------

# Change the skyrmion initial configuration according to the
# chirality of the system (give by the DMI constant)
if args.initial_state_skyrmion_down:
    if args.D > 0:
        sim.set_m(lambda x: generate_skyrmion_down(x, -1))
    else:
        sim.set_m(lambda x: generate_skyrmion_down(x, 1))

elif args.initial_state_skyrmion_up:
    if args.D > 0:
        sim.set_m(lambda x: generate_skyrmion_up(x, 1))
    else:
        sim.set_m(lambda x: generate_skyrmion_up(x, -1))

# The uniform states will be slightly deviated from the z direction
elif args.initial_state_ferromagnetic_up:
    sim.set_m((0, 0.9, 0.9))
elif args.initial_state_ferromagnetic_down:
    sim.set_m((0, 0.9, -0.9))
elif args.initial_state_irregular:
    sim.set_m(irregular_state)
else:
    raise Exception('Set one option for the initial state')


# -------------------------------------------------------------------------


# Debug Information -------------------------------------------------------
# print 'Simulating a {} x {} x {} box'.format(args.box_length,
#                                              args.box_width,
#                                              args.box_thickness)

print('Saturation Magnetisation: {} mu_B'.format(args.mu_s))
print('Exchange constant: {}  meV'.format(args.J))
print('DMI constant: {}  meV'.format(args.D))
if args.k_u:
    print('Anisotropy constant: {}   meV'.format(args.k_u))
if args.B:
    print('Zeeman field: (0, 0, {})  T'.format(args.B / mu0))
print('--------------------------------------')

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Initiate simulation -----------------------------------------------------
# -------------------------------------------------------------------------

# Here we use matplotlib in interactive mode, updating the figure by
# drawing the canvas instead of creating a whole new plot
if args.preview:

    times = np.linspace(0, 4e-9, 500)

    # Data for the Figure
    # Mesh coordinates
    coords = sim.mesh.coordinates
    # Now we copy the magnetisation to not modify the simulation object
    m = np.copy(sim.spin)
    # Reshape and transpose to get an array
    # with [mx, my, mz] elements (matrix)
    m = m.reshape(-1, 3)

    # Figure specifications -----------------------------------------------

    # Aspect according to the square dimensions
    # aspect = args.hex_nx / float(args.hex_ny)
    aspect = 1
    w, h = plt.figaspect(aspect)
    fig = plt.figure(figsize=(w, h))

    ax = fig.add_subplot(111)

    # Colourise the spins according to the mz component
    quiv = ax.quiver(coords[:, 0], coords[:, 1],
                     m[:, 0], m[:, 1], m[:, 2],
                     # Arrow properties (can vary according to the plot)
                     cmap='RdYlBu', width=.008, linewidth=1,
                     scale=1 / 0.05,
                     # Data limits for the colour map
                     clim=[-1, 1]
                     )

    ttime = ax.text(1., 1.05, '',
                    transform=ax.transAxes,
                    # Vertical and horizontal alignment
                    va='center', ha='right')

    tenergy = ax.text(0, 1.05, '',
                      transform=ax.transAxes,
                      va='center', ha='left')

    # Colour bar (append to not distort the main plot)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)

    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    cbar = matplotlib.colorbar.ColorbarBase(cax,
                                            cmap='RdYlBu',
                                            norm=norm,
                                            ticks=[-1, 0, 1],
                                            orientation='vertical',
                                            )

    cbar.set_label(r'$m_z$', rotation=270, labelpad=10, fontsize=16)

    # Interactive mode (this needs so set up a proper backend
    # when importing matplotlib for the first time)
    plt.ion()
    # Set False to avoid the execution of the following code
    plt.show(False)

    # ---------------------------------------------------------------------

    sim.save_vtk()
    # Now run the simulation printing the energy
    for time in times:
        if not run_from_ipython():
            print('Time: ', time, ' s')
            print('Total energy: ', sim.compute_energy(), ' J')
            print('\n')
        sim.run_until(time)

        # Update the vector data for the plot (the spins do not move
        # so we don't need to update the coordinates) and redraw
        m = np.copy(sim.spin)
        # reshape rows, transpose and filter according to top layer
        m = m.reshape(-1, 3)
        quiv.set_UVC(m[:, 0], m[:, 1], m[:, 2])

        # Update title
        ttime.set_text('Time: {:.4f} ns'.format(time * 1e9))
        tenergy.set_text('Energy: {:.6e} J'.format(sim.compute_energy()))

        # fig.show()
        fig.canvas.draw()

    sim.save_vtk()

else:
    # Fidimag automatically saves the last state
    if args.tols:
        sim.driver.set_tols(rtol=args.tols[0], atol=args.tols[1])
    sim.driver.relax(dt=1e-13, stopping_dmdt=args.stopping_dmdt,
                     max_steps=args.max_steps,
                     save_m_steps=args.save_files,
                     save_vtk_steps=args.save_files)

    # Save final states
    sim.save_m()
    sim.save_vtk()


# -------------------------------------------------------------------------
# Files -------------------------------------------------------------------
# -------------------------------------------------------------------------

npy_dir = 'npys/'
vtk_dir = 'vtks/'
txt_dir = 'txts/'

if not os.path.exists(npy_dir):
    os.makedirs(npy_dir)
if not os.path.exists(vtk_dir):
    os.makedirs(vtk_dir)
if not os.path.exists(txt_dir):
    os.makedirs(txt_dir)

# files = [_f for _f in os.listdir('.')
#          if (os.path.isdir(_f) and _f.startswith(args.sim_name))]

# Fidimag vtk and npy files are saved in the sim_name_vtks or sim_name_npys
# folders respectively. We will move them to the vtks/ and npys/ directories

# Remove the vtk folder if it already exists inside the vtks/ folder
if os.path.exists(vtk_dir + args.sim_name + '_vtks/'):
    shutil.rmtree(vtk_dir + args.sim_name + '_vtks/')
try:
    shutil.move(args.sim_name + '_vtks/', vtk_dir)
except IOError:
    pass

# Same for npy folder
if os.path.exists(npy_dir + args.sim_name + '_npys/'):
    shutil.rmtree(npy_dir + args.sim_name + '_npys/')
try:
    shutil.move(args.sim_name + '_npys/', npy_dir)
except IOError:
    pass

# Now do this for the txt files
if os.path.exists(txt_dir + args.sim_name + '.txt'):
    os.remove(txt_dir + args.sim_name + '.txt')
try:
    shutil.move(args.sim_name + '.txt', txt_dir)
except IOError:
    pass

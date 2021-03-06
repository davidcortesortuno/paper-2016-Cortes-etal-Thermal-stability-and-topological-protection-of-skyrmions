# -----------------------------------------------------------------------------
import argparse

parser = argparse.ArgumentParser(description='Plot library for the NEBM data')

parser.add_argument('--method', help='boundary, linear_interpolations, climbing',
                    required=True)

parser.add_argument('--D_list', help='List of DMI values in units of 1e-4 Jm**-2'
                    ' , e.g. --D_list 28 32',
                    required=True, nargs='+')

parser.add_argument('--snapshots', help='Set this option to plot snapshots',
                    action='store_true')

# Parser arguments
args = parser.parse_args()

# -----------------------------------------------------------------------------

import sys
import numpy as np
import os
import glob
import re
from cycler import cycler

eV = 1.602e-19
import matplotlib
import matplotlib.pyplot as plt

from fidimag.atomistic import Sim
from fidimag.atomistic.hexagonal_mesh import HexagonalMesh
from fidimag.atomistic import DMI, UniformExchange, DemagHexagonal, Anisotropy
# Import physical constants from fidimag
import fidimag.common.constant as const
from fidimag.common.nebm_geodesic import NEBM_Geodesic

WORKDIR = os.path.dirname(os.path.realpath(__file__)) + '/'
SAVEDIR = 'nebm_figures/'
if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# -----------------------------------------------------------------------------
# Matplotlib styling

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Lato']
plt.rcParams['font.serif'] = ['Lato']
plt.rcParams['font.cursive'] = ['Lato']
plt.rcParams['font.weight'] = 100
plt.rcParams["text.latex.unicode"] = True
plt.rcParams["xtick.direction"] = 'out'
plt.rcParams["ytick.direction"] = 'out'

# Working:
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'sans:italic'
plt.rcParams['mathtext.bf'] = 'sans:bold:italic'
plt.rcParams['mathtext.sf'] = 'sans'
plt.rcParams['mathtext.default'] = 'it'
plt.rcParams['mathtext.rm'] = 'sans'

plt.rcParams['font.size'] = 28
plt.rcParams['axes.labelsize'] = 28
plt.rcParams['legend.fontsize'] = 26

plt.rcParams['axes.formatter.use_mathtext'] = True
plt.rcParams.update()


# Colors
def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) / 255.
                 for i in range(0, lv, lv // 3))

d_palette1 = ['2C4C8F', 'FF5814', '88B459', 'D91D2C', '8F5536',
              '542437', '6F6F6F']
# Turn it into a palette readable by matplotlib
d_palette1 = [hex_to_rgb(c) for c in d_palette1]
# Set default colour cycle from the library colours
matplotlib.rcParams['axes.prop_cycle'] = cycler('color', d_palette1)

# -----------------------------------------------------------------------------

# DMI parameters:  micromagnetic 1e-4 J m**-2 TO atomistic meV
DMIc = {}
DMIc['26'] = -0.586
DMIc['28'] = -0.631
DMIc['30'] = -0.676
DMIc['32'] = -0.721
DMIc['34'] = -0.766
DMIc['36'] = -0.811
DMIc['38'] = -0.856

# Climbing images for the folders
climbing_image = {}
climbing_image['26'] = 21
climbing_image['28'] = 21
climbing_image['30'] = 20
climbing_image['32'] = 15
climbing_image['34'] = 15
climbing_image['36'] = 15


# Get the last step from a NEBM simulation among smilar folders
def last_step(name):
    step = re.search(r'(?<=GEODESIC_)\d+', name).group(0)
    return step


# FOLDERS and base names:
def base_name_geo(D, k):
    return (WORKDIR + '../sims/nebm/neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_GEODESIC/'
            'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_k1e{1}_GEODESIC'.format(D, k)
            )


def base_name_geo_ci(D, k, ci, folder=False, simname=False):
    bfolder = (WORKDIR + '../sims/nebm/'
               'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_GEODESIC/').format(D)
    sname = 'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_climbing-{2}_k1e{1}_GEODESIC'.format(D, k, ci)

    if folder:
        return bfolder
    elif simname:
        return sname
    else:
        return bfolder + sname


def base_name_ba_geo(D, k, folder=False, simname=False):
    bfolder = (WORKDIR + '../sims/nebm/'
               'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_GEODESIC/').format(D)
    sname = 'neb_stripe_320x185_Co_skdisp_sk-down_fm-up_D{0}e-4_h25e-2nm_k1e{1}_GEODESIC'.format(D, k)

    if folder:
        return bfolder
    elif simname:
        return sname
    else:
        return bfolder + sname


def base_folder_ba(D, k=4):
    files = glob.glob((WORKDIR + '../sims/nebm/'
                       'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_GEODESIC/npys/'
                       'neb_stripe_320x185_Co_skdisp_sk-down_fm-up_D{0}e-4_h25e-2nm_k1e{1}_GEODESIC_*'
                       ).format(D, k))
    files = sorted(files, key=last_step)[-1]
    return files


def base_folder(D, k=4):
    files = glob.glob((WORKDIR + '../sims/nebm/'
                       'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_GEODESIC/npys/'
                       'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_k1e{1}_GEODESIC_*'
                       ).format(D, k))
    files = sorted(files, key=last_step)[-1]
    return files


def base_folder_ci(D, ci, k=4):
    files = glob.glob((WORKDIR + '../sims/nebm/'
                       'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_GEODESIC/npys/'
                       'neb_stripe_320x185_Co_sk-down_fm-up_D{0}e-4_h25e-2nm_climbing-{1}_k1e{2}_GEODESIC_*'
                       ).format(D, ci, k))
    files = sorted(files, key=last_step)[-1]
    return files

# -----------------------------------------------------------------------------

interp_data = {}
nebm_data = {}


# Return a simulation folder according to the specified method
def methods(method, D):
    if method == 'boundary':
        return base_folder_ba(D)
    elif method == 'linear_interpolations':
        return base_folder(D)
    elif method == 'climbing':
        return base_folder_ci(D, climbing_image[D])

# Data ------------------------------------------------------------------------

mesh = HexagonalMesh(0.125, 320, 185, unit_length=1e-9, alignment='square')

# Generate the linear interpolations of the NEBM band from every DMI value from
# the arguments
for D in args.D_list:

    sim = Sim(mesh)
    sim.set_mu_s(0.846 * const.mu_B)
    sim.add(DemagHexagonal())
    sim.add(UniformExchange(J=27.026 * const.meV))
    sim.add(DMI(DMIc[D] * const.meV, dmi_type='interfacial'))
    sim.add(Anisotropy(0.0676 * const.meV, axis=[0, 0, 1]))

    images = [np.load(os.path.join(methods(args.method, D), _file))
              for _file in sorted(os.listdir(methods(args.method, D)),
                                  key=lambda f: int(re.search('\d+', f).group(0))
                                  )]

    nebm = NEBM_Geodesic(sim, images, spring_constant=1e4)
    l, E = nebm.compute_polynomial_approximation(500)

    # Distances for every image from the first image (0th image)
    ll = [0]
    for i in range(len(nebm.distances)):
        ll.append(np.sum(nebm.distances[:i + 1]))
    ll = np.array(ll)

    # Save the interpolation AND the images energy
    interp_data[D] = [l, (E - nebm.energies[0]) / eV, E - nebm.energies[0]]
    nebm_data[D] = [ll, (nebm.energies - nebm.energies[0]) / eV]

# Plot ------------------------------------------------------------------------

markers = ['o', '^', 's', '*', 'h', 'v']
f = plt.figure(figsize=(8, 8/1.6))
ax = f.add_subplot(111)

# Try to plot for every DMI value in the script arguments
for i, D in enumerate(args.D_list):
    c, = ax.plot(interp_data[D][0], interp_data[D][1], '-', lw=2)
    ax.plot(nebm_data[D][0], nebm_data[D][1], 'o', color=c.get_color(),
            ms=10, marker=markers[i], label=str(-DMIc[D])
            )

# plt.ylim([-0.015, 0.02])
# plt.xlim([0, 200])

plt.grid()
l = plt.legend(title=r'$D$' + ' (meV)', loc='lower left',
               ncol=3, bbox_to_anchor=(0, 1.02, 1, 1),
               mode='expand', borderaxespad=0., fontsize=20)
for mark in l.legendHandles:
    mark.set_linestyle('-')
    mark.set_linewidth(2)
    mark.lineStyles['-'] = '_draw_solid'

plt.xlabel(r'$\mathrm{Distance}\,\,\mathrm{from}\,\,\mathbf{Y}_{0}$')
plt.ylabel(r'$\mathrm{Energy}\,\,\mathrm{(eV)}$')

if len(args.D_list) > 1:
    plt.savefig(SAVEDIR + 'energy_bands_D_k1e4_{}.pdf'.format(args.method),
                bbox_inches='tight')
else:
    plt.savefig(SAVEDIR + 'energy_bands_D{}_k1e4_{}.pdf'.format(args.D_list[0],
                                                                args.method),
                bbox_inches='tight')

# Snapshots -------------------------------------------------------------------

# Do not plot if the option was not activated
if not args.snapshots:
    sys.exit()

# We will use the simulation object created before (we only care for the
# spins field and not the interactions)"
coords = mesh.coordinates


def snapshot(sim, image, D, method):
    sim.set_m(np.load(methods(method, D)
                      + '/image_{:06}.npy'.format(image)
                      ))

n_images = len(os.listdir(methods(args.method, args.D_list[0])))
n_plots = n_images
if n_images % 2 != 0:
    n_plots += 1

# -------------------------------------------------------------------------

# We will plot the snapshots of the NEBM band in a two column format plot
# Try to do this for every DMI value
w, h = plt.figaspect(float(10 / 20.))

for D in args.D_list:

    fig, axes = plt.subplots(nrows=int(n_plots / 2), ncols=2,
                             sharey=True, sharex=True)
    fig.set_size_inches(w * 2, h * int(n_plots / 2))

    for i in range(n_images):

        # Get the axes for the i-th state
        ax = axes.flatten()[i]
        snapshot(sim, i, D, args.method)
        data = np.copy(sim.spin.reshape(-1, 3)[:, 2])

        iplot = ax.scatter(coords[:, 0], coords[:, 1],
                           c=data,
                           cmap='gist_earth',
                           s=10, marker='h', lw=0,
                           vmin=-1, vmax=1
                           )

        # bbox_props = dict(boxstyle="circle,pad=0.2", lw=1, alpha=0.)
        plt.text(0.05, 0.8, r'$' + str(i) + r'$',
                 transform=ax.transAxes,
                 fontsize=40,
                 # bbox=bbox_props
                 )

        # Label only the left side column
        if i % 2 == 0:
            ax.set_ylabel(r'$y\,(\mathrm{nm})$', fontsize=28)

        ax.tick_params(axis='both', labelsize=28)

    # Label the x axis for the bottom plots
    for i in [-1, -2]:
        axes.flatten()[i].set_xlabel(r'$x\,(\mathrm{nm})$', fontsize=28)

    # A range to observe the skyrmion
    plt.xlim([0, 79.9])
    plt.ylim([2, 42])

    cax, kw = matplotlib.colorbar.make_axes(axes.flatten()[0], location='top',
                                            anchor=(1, 6.5))
    c = plt.colorbar(iplot,
                     cax=cax,
                     orientation='horizontal')
    c.ax.tick_params(labelsize=22)
    axes.flatten()[0].set_title(r'$m_z$', fontsize=28, y=1.25)
    c.set_ticks([-1, 0, 1])

    # -------------------------------------------------------------------------

    # Remove spacing between subplots
    plt.subplots_adjust(hspace=0.0001, wspace=.0001)

    # Save according to DMI value
    plt.savefig(SAVEDIR + 'nebm_snapshots_D{}_k1e4_{}.png'.format(D, args.method),
                bbox_inches='tight', dpi=100
                )

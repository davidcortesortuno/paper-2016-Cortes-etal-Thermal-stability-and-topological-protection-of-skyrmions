# -----------------------------------------------------------------------------
import argparse

parser = argparse.ArgumentParser(description='Plot library for the NEBM data')

parser.add_argument('--method', help='boundary, linear_interpolations, climbing',
                    required=True)

parser.add_argument('--D_list', help='List of DMI values in units of 1e-4 Jm**-2'
                    ' , e.g. --D_list 28 32',
                    required=True, nargs='+')

# Parser arguments
args = parser.parse_args()

# -----------------------------------------------------------------------------

import numpy as np
import os
import shutil
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

# -----------------------------------------------------------------------------

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

DMIc = {}
DMIc['26'] = -0.586
DMIc['28'] = -0.631
DMIc['30'] = -0.676
DMIc['32'] = -0.721
DMIc['34'] = -0.766
DMIc['36'] = -0.811
DMIc['38'] = -0.856

climbing_image = {}
climbing_image['26'] = 21
climbing_image['28'] = 21
climbing_image['30'] = 20
climbing_image['32'] = 15
climbing_image['34'] = 15
climbing_image['36'] = 15


def last_step(name):
    step = re.search(r'(?<=GEODESIC_)\d+', name).group(0)
    return step


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


def methods(method, D):
    if method == 'boundary':
        return base_folder_ba(D)
    elif method == 'linear_interpolations':
        return base_folder(D)
    elif method == 'climbing':
        return base_folder_ci(D, climbing_image[D])

# Data ------------------------------------------------------------------------

mesh = HexagonalMesh(0.125, 320, 185, unit_length=1e-9, alignment='square')

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

    ll = [0]
    for i in range(len(nebm.distances)):
        ll.append(np.sum(nebm.distances[:i + 1]))
    ll = np.array(ll)

    interp_data[D] = [l, (E - nebm.energies[0]) / eV, E - nebm.energies[0]]
    nebm_data[D] = [ll, (nebm.energies - nebm.energies[0]) / eV]

# Plot ------------------------------------------------------------------------

markers = ['o', '^', 's', '*', 'h', 'v']
f = plt.figure(figsize=(8, 8/1.6))
ax = f.add_subplot(111)

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

plt.savefig('energy_bands_D_k1e4_{}.pdf'.format(args.method),
            bbox_inches='tight')

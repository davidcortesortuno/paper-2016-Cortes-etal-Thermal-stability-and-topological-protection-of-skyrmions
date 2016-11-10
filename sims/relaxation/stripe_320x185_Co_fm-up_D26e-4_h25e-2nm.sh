#!/bin/bash
# Script to simulate a nanotrack made of Cobalt atoms The lattice constant is
# 0.25 nm , which means that we need 320 atoms along x to get 80 nm and 185
# atoms along y to get approx 40 nm, since the height between two rows in the
# hexagonal system is lattice constant * sin(60)

SIMNAME="stripe_320x185_Co_fm-up_D26e-4_h25e-2nm"
python hexagonal_fidimag.py 320 185 0.125 \
        "$SIMNAME" \
        --alignment "square" \
        --D "-0.586" \
        --J 27.026 \
        --mu_s 0.846 \
        --k_u 0.0676 \
        --unit_length "1e-9" \
        --Demag \
        --initial_state_ferromagnetic_up \
        --no_precession \
        --stopping_dmdt 0.0001

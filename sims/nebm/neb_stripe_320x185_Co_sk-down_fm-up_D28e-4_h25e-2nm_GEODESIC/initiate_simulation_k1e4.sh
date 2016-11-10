#!/bin/bash 
SIMNAME="neb_stripe_320x185_Co_sk-down_fm-up_D28e-4_h25e-2nm_k1e4_GEODESIC"
#
python ../hexagonal_neb_fidimag.py \
        320 185 0.125 \
        2000 1000 1000 \
        $SIMNAME \
        --alignment "square" \
        --D "-0.631" \
        --J 27.026 \
        --mu_s 0.846 \
        --k_u 0.0676 \
        --unit_length "1e-9" \
        --Demag \
        --interpolation ../../relaxation_fidimag_atomistic/npys/stripe_320x185_Co_sk-down_D28e-4_h25e-2nm_npys/m_459.npy \
                        26 \
                        ../../relaxation_fidimag_atomistic/npys/stripe_320x185_Co_fm-up_D28e-4_h25e-2nm_npys/m_365.npy \
        --neb_k 1e4 \
        --stopping_dmdt 0.00001 \
        --coordinates Geodesic

#!/bin/bash 
SIMNAME="neb_stripe_320x185_Co_skdisp_sk-down_fm-up_D30e-4_h25e-2nm_k1e4_GEODESIC"
#
python ../hexagonal_neb_fidimag.py \
        320 185 0.125 \
        4000 1000 1000 \
        $SIMNAME \
        --alignment "square" \
        --D "-0.676" \
        --J 27.026 \
        --mu_s 0.846 \
        --k_u 0.0676 \
        --unit_length "1e-9" \
        --Demag \
        --images_files sk_displ_npys/ \
        --neb_k 1e4 \
        --stopping_dmdt 0.000001 \
        --coordinates Geodesic

#!/bin/bash 
SIMNAME="neb_stripe_320x185_Co_sk-down_fm-up_D26e-4_h25e-2nm_climbing-21_k1e4_GEODESIC"
#
python ../hexagonal_neb_fidimag.py \
        320 185 0.125 \
        2000 1000 1000 \
        $SIMNAME \
        --alignment "square" \
        --D "-0.586" \
        --J 27.026 \
        --mu_s 0.846 \
        --k_u 0.0676 \
        --unit_length "1e-9" \
        --Demag \
        --images_files "npys/neb_stripe_320x185_Co_sk-down_fm-up_D26e-4_h25e-2nm_k1e4_GEODESIC_LAST" \
        --neb_k 1e4 \
        --stopping_dmdt 0.00001 \
        --coordinates Geodesic \
        --climbing_image 21

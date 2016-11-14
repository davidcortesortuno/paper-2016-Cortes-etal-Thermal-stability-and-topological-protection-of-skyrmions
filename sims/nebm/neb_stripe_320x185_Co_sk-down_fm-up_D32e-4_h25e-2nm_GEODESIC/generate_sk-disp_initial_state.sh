#!/bin/bash
python ../create_sk_on_track.py \
    --sk_folder "../../relaxation/npys/stripe_320x185_Co_sk-down_D32e-4_h25e-2nm_npys" \
    --fm_folder "../../relaxation/npys/stripe_320x185_Co_fm-up_D32e-4_h25e-2nm_npys" \
    --radius 17

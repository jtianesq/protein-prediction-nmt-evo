#!/bin/bash
set -o pipefail
set -e

entry=$1
model=$2
src=$3
ref_stem=$4

beamsize=10
bpe=false
conv_script=./scripts/convert_primsec.py

translate="python3 $entry translate --model $model --beamsize $beamsize --normalize"
calc_conv="$python $conv_script"

if [[ $bpe == "true" ]]; then
    conv=$($translate < $src | sed -r 's/(@@ )|(@@ ?$)//g' | $calc_conv | cut -f 1 -d ' ')
else
    conv=$($translate < $src)
fi

echo $conv

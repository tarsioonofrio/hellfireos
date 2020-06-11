#!/usr/bin/env bash
DIR=${1:-platform}
bear make clean -C ./${DIR}/noc_3x2/
bear make images -C ./${DIR}/noc_3x2/
cp ./${DIR}/noc_3x2/*.bin ./usr/sim/mpsoc_sim/objects/
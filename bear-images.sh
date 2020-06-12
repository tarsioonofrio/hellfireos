#!/usr/bin/env bash
DIR=${1:-./platform/noc_3x2/}
bear make clean -C ${DIR}
bear make images -C ${DIR}
cp ${DIR}/*.bin ./usr/sim/mpsoc_sim/objects/
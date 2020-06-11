#!/usr/bin/env bash
make clean -C ./platform/noc_3x2/
make images -C ./platform/noc_3x2/
cp ./platform/noc_3x2/*.bin ./usr/sim/mpsoc_sim/objects/
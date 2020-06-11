#!/usr/bin/env bash
make clean -C ./usr/sim/mpsoc_sim/
make ${1:-noc_3x2} -C ./usr/sim/mpsoc_sim/

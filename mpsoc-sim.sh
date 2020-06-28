#!/usr/bin/env bash
cd usr/sim/mpsoc_sim/
./mpsoc_sim ${1:-1} s
#sed -i 's|ffffff||g' ./reports/out*.txt

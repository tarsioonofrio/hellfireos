#!/usr/bin/env bash

make clean      -C ./platform/noc_3x2/
make ${1:-image} -C ./platform/noc_3x2/
../hf-risc/tools/sim/hf_riscv_sim/hf_riscv_sim ./platform/single_core/image.bin #./tp1/out.txt

#!/usr/bin/env bash
NUMBER=${1:-0}
csplit ./usr/sim/mpsoc_sim/reports/out${NUMBER}.txt '/int32_t/'
rm xx00
sed -i '/KERNEL/d' xx01
mv xx01 ./usr-tp2/tp2_files/128x128-images/filter_image_test.h
./usr-tp2/tp2_files/create_bmp ./usr-tp2/tp2_files/128x128-images/filter_image_test.h
mv output.bmp ./usr-tp2/tp2_files/128x128-images/filter_image_test.bmp
feh ./usr-tp2/tp2_files/128x128-images/filter_image_test.bmp
#!/bin/bash

cd /vol8/home/xtu_pcy/TH-MXP-JXPAMG/jxpamg-basic-D
make clean
make -j 32
cd ..

cd /vol8/home/xtu_pcy/TH-MXP-JXPAMG/jxpamg-basic-F
make clean
make -j 32
cd ..

make clean
make


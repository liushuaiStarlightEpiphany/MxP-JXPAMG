#!/bin/bash

cd /vol8/home/xtu_pcy/TH-MXP-JXPAMG/jxpamg-basic-F
make clean 
make -j 32

cd /vol8/home/xtu_pcy/TH-MXP-JXPAMG/jxpamg-basic-F/example
make clean 
make 

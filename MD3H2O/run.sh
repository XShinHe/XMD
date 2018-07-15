#!/bin/bash
python3 rvs.py _put.in nstep 500000
./pimd -x e
python3 rvs.py _put.in nstep 500000000
python3 rvs.py _put.in nsmp 1000
./pimd -x r

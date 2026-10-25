#!/bin/bash
# Experiment A sweep (split vs generic). Output: expA_*.out
cd ~/b-psw/src/cs/rs
PY=~/cs_venv/bin/python3
( for r in 8 16 24 32 40 48 56 64 80 100; do nice -n 10 timeout 1700 $PY split_vs_gen.py comp_T1500.npz 8 21 $r 10 $r | grep SUMMARY; done ) > expA_b8_real.out 2>&1
( for r in 16 48 80; do nice -n 10 timeout 1700 $PY split_vs_gen.py comp_T1500.npz 8 21 $r 10 $((r+1000)) rand | grep SUMMARY; done ) > expA_b8_rand.out 2>&1
( for r in 20 40 60 80; do nice -n 10 timeout 1700 $PY split_vs_gen.py comp_T1500.npz 10 27 $r 6 $r | grep SUMMARY; done ) > expA_b10_real.out 2>&1
echo done > expA.done

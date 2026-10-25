#!/bin/bash
cd ~/b-psw/src/wp/mom
ulimit -v 3000000
P="timeout 3600 nice -n 10 $HOME/chowla_venv/bin/python3 vfdim/m9v_dimcheck.py"
$P gen 9 0 m9_sol_9_0.gp u9,v1,v5,v9 > vfdim/py_9_0.out 2>&1 &
$P gen 8 1 m9_sol_9_1.gp u4,u6,v5,v6 > vfdim/py_9_1.out 2>&1 &
wait
$P gen 11 0 m9_sol_11_0.gp u6,u11,v1,v6 > vfdim/py_11_0.out 2>&1 &
$P gen 10 1 m9_sol_11_1.gp u6,u7,v1,v6 > vfdim/py_11_1.out 2>&1 &
wait
$P sym m9_sol_11_0.gp > vfdim/py_sym11.out 2>&1
echo PYDONE >> vfdim/status

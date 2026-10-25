#!/bin/sh
# Batch 2 (single core): control, exact searches at 2e4/6e4, reruns of the original analysis, closed graph, cells, rows.
cd ~/b-psw/src/cs/dso/audit
PY=$HOME/chowla_venv/bin/python3
rm -f rows_X2000.txt
echo 'IN="../dso_X2000.txt"; DD=5; X=2000; OUTP="rows_X2000.txt"; read("aud_orders.gp")' | timeout 300 nice -n 10 gp -q 2>&1 | tail -1
echo "== my control search X=2000 (CTRL=1)"; /usr/bin/time -f "%e s %M KB" timeout 900 nice -n 10 ./aud_search rows_X2000.txt 5 5 1 0 > ctrl_aud_X2000.sol 2>&1; grep -E "usable|k = | s " ctrl_aud_X2000.sol
echo "== my control hits checked with strong_fermat2"
$PY -c "
import sys, math; sys.path.insert(0, '$HOME/b-psw/src'); from bpsw import strong_fermat2
S=[list(map(int,l.split()[2:])) for l in open('ctrl_aud_X2000.sol') if l.startswith('SOL')]
print('control hits', len(S), 'strong base-2 psp confirmed', sum(strong_fermat2(math.prod(s)) for s in S))"
echo "== original control rerun"; cd .. ; DSO_CONTROL=1 /usr/bin/time -f "%e s %M KB" timeout 900 nice -n 10 $PY dso_analyze.py dso_X2000.txt 5 2000 > audit/rerun_control_X2000.out 2>&1; cd audit; grep -E "k = |bpsw.py|CONTROL brute" rerun_control_X2000.out
for X in 10000 20000 60000; do for F in 1 0; do echo "== my exact search D=5 X=$X filter12b=$F"; /usr/bin/time -f "%e s %M KB" timeout 1800 nice -n 10 ./aud_search rows_X$X.txt 5 5 0 $F 2>&1 | grep -E "usable|k = |SOL| s "; done; done
for X in 10000 20000; do echo "== original analysis rerun X=$X"; cd ..; timeout 900 nice -n 10 $PY dso_analyze.py dso_X$X.txt 5 $X > audit/rerun_an_D5_X$X.out 2>&1; cd audit; grep -E "k = |bpsw.py|pairwise" rerun_an_D5_X$X.out; done
echo "== closed graph"; timeout 600 nice -n 10 $PY aud_closed.py 5 rows_X10000.txt rows_X20000.txt rows_X60000.txt
for D in -7 -11 13 -15 17 -19 -23 29 -31; do tag=$(echo $D | sed 's/-/m/'); timeout 300 nice -n 10 $PY aud_closed.py $D rows_D${tag}_X10000.txt; done
echo "== cells ten D"; timeout 600 nice -n 10 $PY aud_cells.py $HOME/b-psw/src rows_D*_X10000.txt
echo "== rows multiplicity"; timeout 600 nice -n 10 $PY aud_rows.py rows_X20000.txt
echo "== pool metric rerun X=2e4"; cd ..; timeout 1200 nice -n 10 $PY dso_poolmetric.py dso_X20000.txt > audit/rerun_poolmetric_X20000.out 2>&1; cd audit; grep -E "class|best" rerun_poolmetric_X20000.out | head -20
echo BATCH2-DONE

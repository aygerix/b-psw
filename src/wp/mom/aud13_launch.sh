cd ~/b-psw/src/wp/mom
PY=$HOME/chowla_venv/bin/python3
setsid nohup bash -c "NFEV=1000 timeout 1800 nice -n 10 $PY aud13_opt.py 0 10 21 1 22 > aud13_opt_l0a.out 2>&1; NFEV=1000 timeout 1500 nice -n 10 $PY aud13_opt.py 0 8 23 19 33 x noineq > aud13_opt_l0n.out 2>&1" < /dev/null > /dev/null 2>&1 &
setsid nohup bash -c "NFEV=1000 timeout 1800 nice -n 10 $PY aud13_opt.py 1 10 22 1 19 > aud13_opt_l1a.out 2>&1; NFEV=600 timeout 1800 nice -n 10 $PY aud13_opt.py 0 3 24 23 57 > aud13_opt_l0b.out 2>&1" < /dev/null > /dev/null 2>&1 &

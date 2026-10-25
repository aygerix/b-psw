cd ~/b-psw/src/cs/fxchk
for s in LL_F_F LL_F_L LL_L; do timeout 3600 nice -n 10 ./bb2 ../fx/split15/m15_$s.txt 0.225 0.2 >> runs/l225_B.out 2>&1; echo "rc $? $s" >> runs/l225_B.out; done

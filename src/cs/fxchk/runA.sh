cd ~/b-psw/src/cs/fxchk
for s in FL LF_F LF_L; do timeout 3600 nice -n 10 ./bb2 ../fx/split15/m15_$s.txt 0.225 0.2 >> runs/l225_A.out 2>&1; echo "rc $? $s" >> runs/l225_A.out; done

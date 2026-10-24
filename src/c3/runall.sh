#!/bin/sh
# C3 exact pool exhaustion: export class pools, then enumerate in C.
DS="5 -7 -11 13 -15 17 -19 -23 29 -31"
for Y in 100 150 200; do
  case $Y in 100) F=../../data/sn_100_complete.txt;; *) F=../../data/sn_${Y}_aug.txt;; esac
  python3 export.py $Y $F ex$Y.txt $DS < /dev/null > ex$Y.log 2>&1
  ./enum ex$Y.txt > enum$Y.out 2>&1
done

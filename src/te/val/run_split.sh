#!/bin/sh
# split-part control: brute force (ref_split.gp) and engine runs (normal and W1ONE variant), all caps
mkdir -p c6 c6ref
for spec in "5 6000 11" "5 6000 209" "5 6000 1271" "-7 6000 253" "13 6000 391"; do
  set -- $spec; D=$1; X=$2; S=$3
  f=c6ref/D${D}_X${X}_S${S}_ford.txt
  [ -f $f ] || SPL_D=$D SPL_X=$X SPL_S=$S SPL_OUT=c6ref timeout 7200 gp -q ref_split.gp >> c6ref/log.txt 2>&1
  A=$(head -1 $f | sed -E 's/.* A=([0-9]+) .*/\1/'); B=$(head -1 $f | sed -E 's/.* B=([0-9]+) .*/\1/')
  for cap in 0 -1 auto; do
    if [ "$cap" = auto ]; then capopt=""; else capopt="--cap $cap"; fi
    for v in "" _w1; do
      timeout 3600 ./anchor3$v --D $D --X $X --mode fermat --split $S $A 1 $capopt --threads 4 > c6/D${D}_X${X}_S${S}_fermat${v}_cap${cap}.out 2>&1; echo "rc=$?" >> c6/D${D}_X${X}_S${S}_fermat${v}_cap${cap}.out
      timeout 3600 ./anchor3$v --D $D --X $X --mode lucas --split $S 1 $B $capopt --threads 4 > c6/D${D}_X${X}_S${S}_lucas${v}_cap${cap}.out 2>&1; echo "rc=$?" >> c6/D${D}_X${X}_S${S}_lucas${v}_cap${cap}.out
    done
  done
done
touch run_split.done

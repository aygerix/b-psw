#!/bin/bash
cd ~/b-psw/src/cs/e3b
for s in symx_221_p symx_230_p symx_031_p symx_022_p; do
  timeout 900 nice -n 10 Singular -q $s.sing > out/$s.out 2>&1; echo "$s rc=$?" >> out/symx_done
done

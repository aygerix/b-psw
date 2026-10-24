#!/bin/sh
# Theorem B at 1e11, D = -7, m = 62020897 * 41291876971, O = 3 * 160447363 * 2660227837, class (1,1).
# The four non-trivial unitary splits; each prints g = gcd(2^W1 - 1, U_W2(1,2)).
for args in "N -7 320894726 15961367022" "U -7 5320455674 962684178" "U -7 15961367022 320894726" "N -7 962684178 5320455674"; do
  echo "== $args"; /usr/bin/time -l ./biggcd $args 2>&1 | grep -E "^g|^bits|U exact|maximum resident|real"
done

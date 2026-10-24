#!/bin/sh
# Discriminating patterns for Conjecture W' (no monomial family exists for the first four).
for pat in 1,1,1,1,2 1,1,1,2,2 1,1,2,2,2 1,2,2,2,2 1,1,1,3,3 2,2,2,3,3 1,1,1,1,1,1,3 2,2,2,2,2,2,2; do
  timeout 3000 python3 lm_mixed.py $pat 7 150 2>&1 | grep -E "NON-MONOMIAL|pattern"
done

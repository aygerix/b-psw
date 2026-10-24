#!/bin/sh
# Pomerance-pool calibration runs (D = 5, class (1,2), mod-4 side rule).
nice -n 5 timeout 14400 ./pomgen2 1000 18.2 pom_1000_18.2.txt
nice -n 5 timeout 14400 ./pomgen2 1500 15 pom_1500_15.txt
nice -n 5 timeout 14400 ./pomgen2 2000 15 pom_2000_15.txt

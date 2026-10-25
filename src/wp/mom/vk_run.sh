#!/bin/bash
# vk (verifier) runner: vk_run.sh <gp-input-string> <outfile>
cd ~/b-psw/src/wp/mom
ulimit -v 3000000
( time (echo "$1" | timeout 3600 nice -n 10 gp -q -s 2000000000) ) > "$2" 2>&1

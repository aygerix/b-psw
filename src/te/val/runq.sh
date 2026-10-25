#!/bin/sh
# sequential job queue: each line of $1 is "OUTFILE TIMEOUT_S ARGS..." for ./anchor3; marks $1.done at the end
while read out to args; do
  [ -z "$out" ] && continue
  case "$out" in \#*) continue;; esac
  if [ -f "$out" ] && grep -q '^rc=' "$out"; then continue; fi
  /usr/bin/time -f "TIME %e s wall %U s user %M KB" timeout $to ./anchor3 $args > "$out" 2>&1; echo "rc=$?" >> "$out"
done < "$1"
touch "$1.done"

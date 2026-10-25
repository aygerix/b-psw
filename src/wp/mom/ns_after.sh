#!/bin/bash
# wait for PID $1 to exit (checked via /proc, so it also works for PIDs of other users), then run queue $2 (log $3)
while [ -e /proc/$1 ]; do sleep 20; done
cd ~/b-psw/src/wp/mom && exec ./ns_queue.sh $2 > $3 2>&1

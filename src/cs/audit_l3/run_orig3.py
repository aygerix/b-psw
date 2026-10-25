import os
import sys, random
sys.path.insert(0, os.path.expanduser("~/b-psw/src/cs/smallexp"))
import test_zsolve as t
rng = random.Random(int(sys.argv[1]))
for mods, tr in (([3]*5, 20), ([3]*20, 20), ([3]*60, 4), ([3]*100, 2)):
    t.run(mods, tr, rng, "Z3 ^%d" % len(mods))

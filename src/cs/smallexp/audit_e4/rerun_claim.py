#!/usr/bin/env python3
"""Re-run zsolve's own test harness on the claim's cases with fresh seeds (run from smallexp/)."""
import sys, random
sys.path.insert(0, '..')
import test_zsolve as T
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 11
rng = random.Random(seed)
for mods, tr in [([4] * 60, 6), ([8] * 4, 50), ([2, 4, 8, 3, 9], 30), ([4, 12, 36], 20), ([6] * 20, 4), ([12] * 4, 20),
                 ([4] * 20, 20), ([16] * 3, 20), ([4] * 100, 2)]:
    lab = "x".join(f"Z{m}" for m in sorted(set(mods))) + f" ^{len(mods)}" if len(set(mods)) == 1 else str(mods)
    T.run(mods, tr, rng, lab)

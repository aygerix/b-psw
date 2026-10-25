#!/usr/bin/env python3
import os
"""src/cs/e3b/advchk/recall_chk.py -- adversarial recall check of node_model.py on a KNOWN non-empty overdetermined
structure.  Takes symmetric (2,2,2,0,0) points found by strat_lp --sym (file rows (d, gd, clr, pts)), builds the
node-model vector for structure 1|1|2|2|3|3, prints the node-model residual at the known point, then runs node_model's
own LM from perturbations of size eps (basin test).  Usage: recall_chk.py symfile [npts]"""
import sys, numpy as np
sys.path.insert(0, os.path.expanduser('~/b-psw/src/cs/e3b'))
from node_model import Node, parse
from strat_lp import lm_batch
from refine import lucas_table

rows = [eval(l) for l in open(sys.argv[1])]
npts = int(sys.argv[2]) if len(sys.argv) > 2 else 5
model = Node(parse("1|1|2|2|3|3"))
rng = np.random.default_rng(1)
for d, gd, clr, pts in rows[:npts]:
    clr = np.array(clr)
    tab = lucas_table(clr, 5)
    lp = [[t for (t, r) in row if r < 1e-6] for row in tab]
    if [len(p) for p in lp] != [2, 2, 2, 0, 0]:
        continue
    z = np.zeros(model.n, dtype=complex)
    z[:4] = clr[5:9]
    w = lp[0][0]
    z[4:10] = [w, -w, lp[1][0], lp[1][1], lp[2][0], lp[2][1]]
    z[model.coef[4]] = clr[3]; z[model.coef[5]] = clr[4]
    res0 = np.linalg.norm(model.F(z[None])[0])
    line = f"residual at known point {res0:.2e};"
    for eps in (1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 0.3):
        Z = z[None] + eps * (rng.normal(size=(200, model.n)) + 1j * rng.normal(size=(200, model.n)))
        with np.errstate(all='ignore'):
            Z, nf = lm_batch(model, Z.copy(), iters=200, prune=True)
        back = np.sum((nf < 1e-11) & (np.linalg.norm(Z - z[None], axis=1) < 1e-5))
        line += f"  eps={eps:g}: {back}/200"
    print(line, flush=True)

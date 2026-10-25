# Review-only: build anchor3_hook.c from anchor3.c with two observation hooks (engine itself untouched).
#  gHookSurv: branch-A prefilter survivors are recorded as hits without running BPSW.
#  gHookC:    branch-C p1 candidates that pass the deferred checks are passed to hookC() (before tab_find).
import sys
src = open(sys.argv[1]).read()
reps = [
 ('static int gQuiet = 0;\n',
  'static int gQuiet = 0;\nstatic int gHookSurv = 0, gHookC = 0; static void hookC(u64 p2, u64 p3, u64 p1);\n'),
 ('if(gMode == MODE_BPSW){ ts->c.Afinal++; make_n_m(ts->n, p3, mm); if(bpsw_mpz(&ts->b, ts->n)) record_hit("HIT", ts, 0, 0, p3, mm, \'A\'); }',
  'if(gMode == MODE_BPSW){ ts->c.Afinal++; make_n_m(ts->n, p3, mm); if(gHookSurv || bpsw_mpz(&ts->b, ts->n)) record_hit("HIT", ts, 0, 0, p3, mm, \'A\'); }'),
 ('if(gMode == MODE_BPSW){ ts->c.Afinal++; make_n_m(ts->n, p3, m); if(bpsw_mpz(&ts->b, ts->n)) record_hit("HIT", ts, 0, 0, p3, m, \'A\'); continue; }',
  'if(gMode == MODE_BPSW){ ts->c.Afinal++; make_n_m(ts->n, p3, m); if(gHookSurv || bpsw_mpz(&ts->b, ts->n)) record_hit("HIT", ts, 0, 0, p3, m, \'A\'); continue; }'),
 ('      if(!haveW){ if(!w2_holds(q1, p2, p3, ts->t1)) continue; }',
  '      if(!haveW){ if(!w2_holds(q1, p2, p3, ts->t1)) continue; }\n      if(gHookC){ hookC(p2, p3, q1); continue; }'),
]
for a, b in reps:
    c = src.count(a)
    if c != 1: sys.exit("pattern count %d: %r" % (c, a[:70]))
    src = src.replace(a, b)
open(sys.argv[2], 'w').write(src)
print("ok")

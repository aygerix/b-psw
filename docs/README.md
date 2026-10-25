# b-psw: research repository guide

This repository supports a research program on the Baillie–PSW (BPSW) problem. The goal is to find an explicit composite n
that passes the exact strong BPSW test in `src/bpsw.py`: strong Fermat base 2 plus strong Lucas, Selfridge Method A, P = 1.
Secondary goals are rigorous structural results and better search methods.

## Where to start

- `.research/number-theory/NOTEBOOK.md` is the research record.
  - Read the last section first: **"Checkpoint, end of session 9 — READ THIS FIRST"**. It covers:
    - the current state;
    - every result from session 9;
    - the jobs still running on the VM, and how to collect them;
    - the ranked open questions and the routes now closed;
    - the pre-commit checklist.
  - Older sections keep in-place `[CORRECTED ...]` / `[REFUTED ...]` markers, so read claims with their markers.
- `.research/number-theory/LOG.jsonl` is a one-line-per-event log (results, theorems, corrections).
- `CLAUDE.md` holds the project protocol.
- `docs/workflow_syntheses/` holds the full outputs of the session-9 multi-agent workflows, including line results and adversarial-check verdicts:
  - `constructive_shortcut_synthesis.md`: the common metric G_C(b) and the mechanism ranking;
  - `cubic_williams_all_k_synthesis.md`: cubics for k ≥ 9;
  - `*_full.json`: raw results of the pencil, k = 9 audit, round-2 and round-3 workflows.
- `docs/VM_ONLY_DATA.txt` lists bulk data kept only on the VM (path and size). It is all regenerable: instance files, census records, numerical-search `.npz`.

## Status labels used throughout

| Label | Meaning |
|---|---|
| PROVED | Complete argument, usually checked by an independent verifier |
| COMPUTED | Exact computation, with its range and assumptions stated |
| HEURISTIC | Probabilistic or numerical evidence; it ranks routes but never closes them |
| REFUTED | Shown false |

Many entries say "conditional on …". Read the condition.

## Directory map (`src/`)

| Directory | Contents |
|---|---|
| `bpsw.py` | The exact test (reference implementation; `check_submission` reads a JSON with factor1, factor2) |
| `fam/` | Faithful BPSW in GP (`bpsw.gp`), polynomial-family screens, Lemma 5 family, Aurifeuillian scans |
| `alg/` | Gröbner / elimination computations for Williams families (Singular, Magma logs) |
| `wp/` | Williams-family work: quadratic/cubic/quartic patterns, census over finite fields (`wp/census/`) |
| `wp/mom/` | Moment/boundary method (Theorems 11, 11', 12), pencils (Theorem 14: `pc7_ext_verify.gp`, `pc7_mono_thm.gp`, audits `pc7_audA_*`, `pc7_adv_*`), cubics k ≥ 9 (`m9_*`, `xn_*`, `ak_*`, `ns_*`, `sf_*`), Theorem 13 audits (`aud13_*`, `aud13_A/`) |
| `wp/hp/` | Half-plane criterion test for cubics (`c4_test.py`) |
| `te/` | Theorem E engine (`anchor3.c`), production driver `run_E.sh`, validation (`val/`), runs (`runs/`) |
| `us/` | Theorem B: unitary-split scans for split primes (`usplit*.c`, `biggcd.c`, `drive.gp`), outputs |
| `lp/` | L1 last-prime gcd search |
| `c3/`, `pom/` | Pomerance pools: CHM twin-smooth exports, pool generators (`pomgen*.c`), instance building, k-tree analysis |
| `dso/`, `dso2/` | Doubly-small-order primes (Attempts D, D2) |
| `cyc/`, `lot/` | Cyclotomic family sweeps, structured strong-psp lotteries |
| `zs/`, `fw/`, `rw/`, `lem/`, `wief/`, `tt/`, `v2/` | Smaller studies (zero-sum DFS, forced-Psi families, rational Williams members, lemma checks, Wieferich/Lucas-Wieferich, Kummer independence, side assignment) |
| `vm/` | VM job-queue scripts |
| `cs/` | Session-9 constructive-shortcut work, one subdirectory per line (next table) |

### `src/cs/` lines

| Directory | Contents |
|---|---|
| `smallexp/`, `audit_l3/` | Small-exponent layered zero-sum algorithms (`zsolve.py`) and their audits |
| `rs/` | Real-instance structure: Sylow split, lattices/BKZ against controls |
| `ident/` | Universal forms, Theorem U, sign lemma, kill checks |
| `sel/` | Frobenian selectors, censuses (`census.c`, `census2.c`), side/partition optimisation |
| `e1/` | Small-b frontier certificate: `census3.c`, CHM pool, Lagrangian/min-cut branch-and-bound (`bbc.c`) |
| `fx/`, `fxchk/` | Certificate extension (24 selectable D, λ = 0.225, rows ≤ 300), D = −15 explanation, Selfridge reachability, and the independent check |
| `dso/` | Doubly-primitive census (gcd of cyclotomic products), cells, exact subset search, Lucas-Wieferich check |
| `e3/` | k = 3 shared-A Korselt variety, rho-family, minimum-Lucas-defect theorem |
| `e3b/` | k = 5 shared-quadratic-A Korselt families (numerical; FSP and symmetric loci) |
| `e4/` | Forced-index coupling experiment (pre-registered, with controls) |
| `hd/`, `hdaudit/` | Half-defect inequality: refutation (G5, doubling chain), C-lemma, A-reduction, sweep lemma |
| `oc/` | Order channel: Lemma OC3, Lemma OC4, Theorem M (Mersenne classes), family searches |

Most `audit/`, `adv/`, `advchk/` subdirectories hold the independent re-computations by the adversarial checkers.

## Running computations

All experiments run on the project VM, reached with:

```
gcloud compute ssh <vm> --project <project> --zone <zone>
```

The work tree there is `~/b-psw/src`.

Rules:
- Use `nice` and `timeout`.
- Time a small instance before any long run. Never launch a run that the estimate says cannot finish.
- Kill your own processes by PID only.
- Compiled C tools are not committed: rebuild them from the `.c` sources (see `.gitignore`).

## Before committing (public repository)

- Keep GCP project and zone identifiers, IP addresses, tokens and e-mail addresses out of committed files.
  - The notebook writes `<project> <zone>`.
  - No such identifiers are in the session-9 files (checked 2026-10-25).
- About 56 new files, mostly logs and VM run scripts, contain home-directory paths (`/home/<user>/…`).
  - Replace these with `~` or `$HOME` in text records before committing.
- Large regenerable data stays on the VM (`docs/VM_ONLY_DATA.txt`).

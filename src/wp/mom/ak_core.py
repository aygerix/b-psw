# ak_core.py -- shared helpers for the all-k line (real cubic Williams families).
# Normal form f = 2t+1 + t(t+1)(a t + c), 0 < a < 32, c in (a - sqrt(8a), sqrt(8a)).
# Parametrisation: a = 32/(1+exp(-al)), c = a/2 + (sqrt(8a) - a/2) tanh(ga).
import numpy as np
def ac_from(al, ga):
    a = 32.0/(1.0+np.exp(-al))
    c = a/2 + (np.sqrt(8*a) - a/2)*np.tanh(ga)
    return a, c
def special(a, c):
    # upper roots of g = a t^2 + (a+c) t + c + 2 and h = a t^2 + c t + 2
    def up(A, B, C):
        d = B*B - 4*A*C + 0j
        r = (-B + np.sqrt(d))/(2*A)
        return np.where(r.imag >= 0, r, np.conj(r))
    return up(a, a + c, c + 2), up(a, c, 2.0*np.ones_like(a))
def fval(a, c, t):
    return 2*t + 1 + t*(t+1)*(a*t + c)
def logmat(a, c):
    # L[i,j,0] = log f_i(w_j), L[i,j,1] = log f_i(v_j) (complex logs; diagonal = 0)
    w, v = special(a, c)
    Fw = fval(a[:, None], c[:, None], w[None, :])
    Fv = fval(a[:, None], c[:, None], v[None, :])
    k = len(a)
    Fw[np.arange(k), np.arange(k)] = 1.0
    Fv[np.arange(k), np.arange(k)] = -1.0
    return Fw, Fv

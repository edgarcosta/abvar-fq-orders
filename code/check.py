from multiprocessing import cpu_count
from .utils import chunkify
from sage.all import (
    parallel,
    ceil,
)
from utils import (
    is_real_weil,
    trace_polynomial,
    make_weil,
)
@parallel(ncpus=cpu_count())
def check_multi(l):
    return [check(*t) for t in l]

def check(q, a, b, c, hs):
    for h in hs:
        w = make_weil(q, h)
        wt = trace_polynomial(q, w.reverse())
        if not (is_real_weil(q, wt) and w(1) == c):
            return False
        l = a-c
        u = b-c
        for d in [l, u-1]:
            h = wt + d
            if not is_real_weil(q, wt):
                return False
    return True


def check_all(q, res):
    t = 0
    for I, O in check_multi(chunkify([(q,) + k + v for k, v in res.items() if k != 'walltime'], ceil(len(res)/256))):
        assert isinstance(O, list)
        if all(O) != True:
            return False, I
        t += len(O)
    return True, t == len(res)

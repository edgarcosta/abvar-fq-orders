from multiprocessing import cpu_count
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from time import time
import pathlib
current_directory = pathlib.Path(__file__).parent.absolute()
import os
# optional: https://github.com/kedlaya/root-unitary/
# and to make use of parallelism one needs to enable OpenMP support
# see the first lines of weil_polynomials.pyx of the root-unitary package
root_unitary_directory = pathlib.PurePath(current_directory, '..', 'root-unitary/')
root_unitary_directory = pathlib.PurePath('/home/edgarcosta/', 'root-unitary/')
os.environ["SAGE_NUM_THREADS"] = os.environ["OMP_NUM_THREADS"] = str(cpu_count())

from sage.all import (
    PolynomialRing,
    RealSet,
    ZZ,
    cached_function,
    cartesian_product_iterator,
    ceil,
    factor,
    is_prime_power,
    load,
    prime_divisors,
    prod,
    save,
    floor,
)
from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
from utils import (
    chunkify,
    hondatatefilter,
    hw,
    make_weil,
    poly_to_tuple,
)
from deform import (
    deform_batch
)





@cached_function
def weil_poly_gq(g, q):
    filename = pathlib.Path(current_directory, 'weil_cache/', 'g%d.q%d.sobj' % (g, q))
    filenamew = pathlib.Path(current_directory, 'weil_cache/', 'g%d.q%d.w.sobj' % (g, q))
    if filename.exists():
        return load(filename.as_posix())
    if filenamew.exists():
        out, nonordinary = load(filenamew.as_posix())
    else:
        try:
            load(pathlib.PurePath(root_unitary_directory, 'weil_polynomials.pyx'))
            W = WeilPolynomials(2*g, q, parallel=True)
            W.num_processes = 500
        except (FileNotFoundError, OSError):
            print("Couldn't load 'weil_polynomials.pyx', using sage's serial version")
        p = prime_divisors(q)[0]
        out = []
        nonordinary = []
        st = time()
        for f in W:
            if f[g] % p != 0: # ordinary
                out.append(f)
            else:
                nonordinary.append(f)
        print('q=%d g=%d Weil polynomial generation done %.2f s len(nonord) = %d' % (q, g, time() - st, len(nonordinary)))
        save([out, nonordinary], filenamew)
    st = time()
    for _, O in hondatatefilter([(fs,p,q) for fs in chunkify(nonordinary, ceil(len(nonordinary)/256))]):
        assert isinstance(O, list)
        out.extend(O)

    print('q=%d g=%d Honda-Tate done %.2f s' % (q, g, time() - st,))
    save(out, filename)
    return out

def cover_small(q, end, n=None):
    if n is None:
        if q <= 3:
            n = 2
        elif q == 4:
            n = 3
        else:
            n = 4
    filename = pathlib.Path(current_directory, 'weil_cache/', 'small_q%d.n%d.end%d.sobj' % (q, n, end))
    if filename.exists():
        return load(filename.as_posix())
    start_time0 = start_time = time()
    res = seed_small(q, n)
    print('seed_small %.2f s' % (time() - start_time,))

    start_time = time()
    extend_multiplicatively(q, res, end)
    print('extend_multiplicatively %.2f s' % (time() - start_time,))

    start_time = time()
    start = hw(q, n + 1).get_interval(0).lower()
    start_time = time()
    extend_by_neighbors(q, res, n, start, end, maxk=3)
    print('extend_by_neighbors %.2f s' % (time() - start_time,))

    start_time = time()
    resint, coverage, nonordinary, w = convert_into_intervals(q, res, n, end)
    print('convert_into_intervals %.2f s (%.2f s)' % (time() - start_time, w))

    walltime = start_time0 - time() + w
    out = (resint, coverage, nonordinary, n, walltime)
    save(out, str(filename.as_posix()))
    return out



def seed_small(q, n):
    res = defaultdict(set)
    for g in range(1, n+1):
        for elt in weil_poly_gq(g, q):
            elt1 = int(sum(elt.list()))
            res[elt1].add(poly_to_tuple(elt))
    return res

def weil_poly_lead(q, init):
    # we only keep ordinary polynomials
    # to avoid the expensive Honda-Tate test
    R = PolynomialRing(ZZ, 'x')
    if isinstance(init, tuple):
        init = [init]
    res = defaultdict(set)
    p = prime_divisors(q)[0]
    # squarefree=True will segfault
    for i in init:
        n = len(i) + 1
        #init = [(i,0) for i in init]
        for elt in WeilPolynomials(2*n, q, lead=i, polring=R):
            if elt[n] % p == 0: # this avoids the expensive Honda-Tate test
                continue
            elt1 = int(sum(elt.list()))
            res[elt1].add(poly_to_tuple(elt))
    return res



def extend_multiplicatively(q, res, cut):
    def poly_iter(m):
        f = factor(m)
        if all(pi**ei in res for pi, ei in f):
            return cartesian_product_iterator([res[pi**ei] for pi, ei in factor(m)])
        else:
            return [].__iter__()
    for m in range(2, cut):
        if is_prime_power(m):
            continue
        gs = [poly_to_tuple(prod((make_weil(q, elt) for elt in g)).reverse()) for g in poly_iter(m)]
        if gs:
            res[m].update(set(gs))

def extend_by_neighbors(q, res, minn, start, end, maxk=20):
    calls = set({})
    assert maxk > 1
    with ProcessPoolExecutor() as executor:
        for k in range(1, maxk):
            jobs = {}
            for m in range(start, end):
                if m not in res:
                    for ell in range(m-k, m+k+1):
                        if ell in res:
                            for halfg in res[ell]:
                                # take the first not yet tried
                                if len(halfg) >= minn + 1 and halfg[:-2] not in calls:
                                    jobs[executor.submit(weil_poly_lead, q, halfg[:-2])] = halfg[:-2]
                                    calls.add(halfg[:-2])
                                    break
            if not jobs:
                break
            for job in as_completed(jobs):
                for m, gs in job.result().items():
                    res[m].update(gs)


def convert_into_intervals_core(q, res, start, end):
    start_time = time()
    p = prime_divisors(q)[0]
    coverage = RealSet()
    ordinarycoverage = RealSet()
    resint = {}
    for m in reversed(range(start, end)):
        if m in coverage and m in ordinarycoverage:
            continue
        if m in res:
            if m not in coverage:
                    interval, v = deform_batch(p, q, m, res[m], rev=False)
                    resint[interval] = v
                    coverage += RealSet.closed_open(*interval)
                    if len(v[2]) > 1:
                        ordinarycoverage += RealSet.closed_open(*interval)
            if m not in ordinarycoverage:
                for h in res[m]:
                    if h[-1] % p != 0:
                        ordinarycoverage += RealSet.closed_open(m, m+1)
                        break
    return resint, coverage, ordinarycoverage, time() - start_time


def convert_into_intervals(q, res, n, end):
    resint = {(1,2):(1, 0, ((1,),), 0)}
    coverage = RealSet.closed_open(1,2)
    ordinarycoverage = RealSet.closed_open(1,2)
    walltime = 0
    start_time = time()
    with ProcessPoolExecutor() as executor:
        jobs = {}
        step = max(100, floor(end/cpu_count()))
        for start in range(2, end, step):
            jobs[executor.submit(convert_into_intervals_core, q, res, start, min(start+step, end))] = (start, start + step)
        for job in as_completed(jobs):
            r, c, o, w = job.result()
            resint.update(r)
            coverage += c
            ordinarycoverage += o
            walltime += w
    walltime += time() - start_time
    return resint, coverage, [m for m in range(2, end) if m not in ordinarycoverage], walltime
#    for m in range(2, cut):
#        if m not in coverage:
#            c = wrap_construct_h(q, m, runs=1000, target_interval=None)
#            if c is not None:
#                h, r, m, _, interval = c
#                resint[interval] = (m, r, tuple(tuple(poly_to_tuple(elt) for elt in h)))
#                coverage += RealSet.closed_open(*interval)
#                if len(h) > 1:
#                    ordinarycoverage += RealSet.closed_open(*interval)
#            # TODO: look at neighbors of m, and see if there is any deformation that covers m
#    nonordinary = [m for m in range(2, cut) if m not in ordinarycoverage and m in coverage]
#    for m in nonordinary[:]:
#        if m in res:
#            hs = res[m]
#            for h in hs:
#                if h[-1] % p != 0:
#                    nonordinary.remove(m)
#                    break
#        if m not in nonordinary:
#            continue
#        for (a, b) in sorted(resint):
#            if m < a:
#                break
#            if a <= m < b:
#                center, _, hs = resint[(a,b)]
#                if len(hs) > 1:
#                    nonordinary.remove(m)
#                    break
#                shift = m - center
#                for h in hs:
#                    if (h[-1] + shift)% p != 0:
#                        nonordinary.remove(m)
#                        break
#                    # try to make it ordinary with the usual trick
#                    if len(h) > 2:
#                        g = list(h)
#                        g[-1] += (q + 1) + shift
#                        g[-2] += -1
#                        w = make_weil(q, tuple(g)).reverse()
#                        assert w(1) == m
#                        if w.is_weil_polynomial():
#                            nonordinary.remove(m)
#                            break
#                        g = list(h)
#                        g[-1] += -(q + 1) + shift
#                        g[-2] += 1
#                        w = make_weil(q, tuple(g)).reverse()
#                        assert w(1) == m
#                        if w.is_weil_polynomial():
#                            nonordinary.remove(m)
#                            break
#
#                if m not in nonordinary:
#                    break
#    return resint, coverage, nonordinary
#def exhaust_small(q, cut, n=None):
#    if n is None:
#        if q <= 3:
#            n = 2
#        elif q == 4:
#            n = 3
#        else:
#            n = 4
#    res = seed_small(q, n)
#    resint, coverage, nonordinary = extend_small(q, res, n, cut)
#
#    return resint, nonordinary, coverage, n


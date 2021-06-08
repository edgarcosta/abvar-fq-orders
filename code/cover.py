from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from multiprocessing import cpu_count
from time import time
from sage.all import (
    RealSet,
    RR,
    floor,
    ceil,
    log,
)
from small import cover_small
from construct import wrap_construct_h
from utils import hw


def cover_interval(q, start, end, runs=1000):
    start_time = time()
    start = floor(start)
    end = ceil(end)

    res = dict({})
    # first interval
    center = start
    out, rad, interval = wrap_construct_h(q, center, target_interval=RealSet.closed_open(start, start+1), runs=runs)
    res[interval[0], interval[1] + 1] =  out
    coverage = RealSet.closed_open(interval[0], interval[1] + 1)
    calls = set()
    while coverage.get_interval(0).upper() < end or coverage.n_components() > 1:
        if coverage.n_components() > 1:
            # we made a hole, let's fix it!
            # aim for the middle
            center = floor((coverage.get_interval(0).upper() + coverage.get_interval(1).lower())/2)
            if center in calls:
                center = coverage.get_interval(0).upper()
        else:
            # hoping that the deformation doesn't decrease too much from the previous iteration
            center = min(coverage.get_interval(0).upper() + floor(0.9*rad), end)
        target = RealSet.closed_open(start ,end).intersection(coverage.complement())
        out, rad, interval = wrap_construct_h(q, center, target_interval=target, runs=runs)
        res[interval[0], interval[1] + 1] = out
        coverage += RealSet.closed_open(interval[0], interval[1] + 1)
    return coverage, res, time()-start_time


startparam = {2: 3, 3: 6, 4: 16, 5: 29, 7: 50000, 8: 76000, 9: 22950, 11: 67880}
#startparam = {2: 3, 3: 6, 4: 16, 5: 29, 7: 50000, 8: 40897792, 9: 22950, 11: 67880}
def cover_interval_parallel(q, end, res=None, step=10**3, total_wall=0, coverage=None, nonordinary=[]):
    if res is None:
        if q not in startparam:
            raise NotImplementedError
        coverage, res, w = cover_interval(q, startparam[q], startparam[q] + step)
        assert coverage.n_components() == 1
        total_wall += w
        smallcut = coverage.get_interval(0).lower()
        exhaust, c, nonordinary, nsmall, w = cover_small(q, smallcut)
        res.update(exhaust)
        coverage += c
        total_wall += w
        print("Initial coverage:", coverage)
    else:
        smallcut = None
        nsmall = None
        if coverage is None:
            coverage = sum([RealSet.closed_open(a, b) for a, b in res], coverage)
    submitted_coverage = RealSet(coverage.get_interval(coverage.n_components() -1)) # take the last chunk
    last_print = time()
    with ProcessPoolExecutor() as executor:
        jobs = {}
        while jobs or submitted_coverage.get_interval(0).upper() < end:
            while len(jobs) < 120 + cpu_count() and submitted_coverage.get_interval(0).upper() < end:
                s = submitted_coverage.get_interval(0).upper()
                e =  min(s + step, end)
                if submitted_coverage.n_components() > 1:
                    # we have a hole, let's fix it!
                    e = min(e, submitted_coverage.get_interval(1).lower())
                submitted_coverage += RealSet.closed_open(s, e)
                #print("submitted", (s,e))
                jobs[executor.submit(cover_interval, q, s, e)] = [s, e]
            done, not_done = wait(jobs, timeout=0.25, return_when=FIRST_COMPLETED)
            if done:
                for job in done:
                    s, e = jobs[job]
                    try:
                        c, r, w = job.result()
                    except Exception as exc:
                        print('%s generated an exception: %s' % ((s,e), exc))
                    else:
                        #print('%s done in %.2f s with %d' % ((s,e), w, len(r)))
                        if w < 30:
                            step = max(step, 2*(e - s))
                        total_wall += w
                        #print(r)
                        res.update(r)
                        coverage += c
                        submitted_coverage += c
                    del jobs[job]
            if time() - last_print > 30:
                print('submitted:', submitted_coverage)
                print(f'jobs = {len(jobs)}', "up to %d^%.2f" % (q, RR(coverage.get_interval(0).upper()).log(q),), "coverage = ", coverage)
                #print(sorted(jobs.values()))
                last_print = time()
    return coverage, res, total_wall, nonordinary, nsmall, smallcut



def cover_and_save(filename, **kwargs):
    c, r, w, nonordinary, nsmall, smallcut = cover_interval_parallel(**kwargs)
    write(filename, kwargs['q'], c, r, w, nonordinary, nsmall, smallcut)
    return c, r, w, nonordinary, nsmall, smallcut


def write(filename, q, coverage, res, w, nonordinary, nsmall, smallcut):
    maxm = coverage.get_interval(coverage.n_components() - 1).upper()
    midm = coverage.get_interval(coverage.n_components() - 1).lower()
    if nsmall:
        interval = sum([hw(q,g) for g in range(nsmall + 1)], RealSet())
        nonordinary_proven = [m for m in nonordinary if m in interval]
        nonordinary_unproven = [m for m in nonordinary if m not in interval]
    else:
        nonordinary_proven = []
        nonordinary_unproven = nonordinary
    with open(filename, 'w') as W:
        W.write('# q = %d\n' % q)
        W.write('# Integers covered = %s\n' % coverage)
        W.write('# %d ~ %d^%.2f\n' % (midm, q, log(midm, q)))
        W.write('# %d ~ %d^%.2f\n' % (maxm, q, log(maxm, q)))
        if nonordinary:
            if nsmall:
                W.write('# exhaustive for g <= %d that covers %s\n' % (nsmall, interval))
            if nonordinary_proven:
                W.write('# orders not realizable with ordinary AV = %s\n' % (nonordinary_proven,))
            if nonordinary_unproven:
                W.write('# orders potentially not realizable with ordinary AV = %s\n' % (nonordinary_unproven,))
        if not nonordinary:
            W.write('# every order is realizable with an ordinary AV\n')
        W.write('# total wall time = %.2f seconds\n' % (w))
        W.write("""#
# Each line has the following format:
# a:b:c:r:pol:mu_ord_upper
# where:
# - a, b, c, r are integers representing the closed
# - pol, is a list of lists of integers, each list represents a Weil polynomial of degree 2n given by their last n coefficients,
# - mu_ord is a real number
#
# For each (monic) Weil polynomial h(x) represented by pol we have h(1) = c.
# There are two cases:
#
# If r = 1, then the polynomials g(x) = h(x) + (t-c)*x^n with t an integer in [a, b),
# are also a Weil polynomials, and g(1) = t.
# Furthermore, the length of pol is at most two, and if two then the middle coefficients
# differ modulo p, and thus every number in [a, b) maybe be represented by an order of 
# an ordinary abelian variety.
# In this case, mu_ord_upper = 0, and doesn't play a role.
#
# If r > 1, then pol has length one, and thus only represents one polynomial h(x).
# Write g(x) = h(x) + q*x^(n-1) + (q + 1)*x^n + x^(n+1).
# This is also a Weil polynomial.
# Consider S_h (resp. S_g) the set of integral polynomials obtained by integrally modifying
# the coefficients x^n, ..., x^(n + (r-1)) by floor(q/2) of h(x) (resp. g(x)), and the lower
# coefficients accordingly to keep the q-symmetry.
# One can then realize every integer in [a, b) as f(1) with f in S_h (resp. S_g).
# mu_ord_upper is an upper bound for mu_ord as defined in Section 9
#
""")
        for a, b in sorted(res):
            c, r, hs, mu = res[(a,b)]
            W.write(":".join(map(str, [a,b,c,r,list(map(list,hs)), mu])) + '\n')

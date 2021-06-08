from collections import defaultdict

from sage.all import (
    ZZ,
    PolynomialRing,
    CBF,
    RBF,
    RealIntervalField,
    QQ,
    sqrt,
)

from utils import (
    make_weil,
    trace_polynomial,
    is_real_weil,
    hondatate,
)

from construct import (
    NotWeil,
)



def deform(q, f, deformation):
    R = PolynomialRing(ZZ, "y")
    f = R(f)
    if not is_real_weil(q, f):
        raise NotWeil
    if not f.is_squarefree():
        return 0, 0
    deformation = R(deformation)
    RRR = RBF if deformation == 1 else RealIntervalField(600)

    issues = set()
    for sign in [1, -1]:
        point = sign * 2 * sqrt(q)
        issues.add(RRR(-f(point) / deformation(point)))

    if deformation == 1:
        # we get a double root when we hit a local extrema
        fder = f.derivative()
        fder = fder // fder.gcd(fder.derivative())
        ## deal with degree one factors, as RIF and QQ don't go so well
        #for fac, e in fder.factor():
        #    assert e == 1
        #    if fac.degree() == 1:
        #        issues = issues.union([RRR(-f(elt[0])) for elt in fac.roots(QQ)])
        #        fder = fder // fac
        roots = [elt.real() for elt in fder.roots(CBF, multiplicities=False) if 0 in elt.imag()]
        assert len(roots) == fder.degree()
        issues = issues.union([-f(elt) for elt in roots])
        if any(0 in elt for elt in issues): # not running any risks
            return (0, 0)
        low = sorted(elt.upper().ceil() for elt in issues if elt < 0)[-1]
        upper = sorted(elt.lower().floor() for elt in issues if elt > 0)[0]
        # deformation space is connected
        assert is_real_weil(q, f + low) and is_real_weil(q, f + upper)
        return low, upper
    else:
        # this code is not longer used
        S = PolynomialRing(R, "yt")
        # all the time is spent here!!!
        delta = (f(S.gen()) + R.gen() * deformation(S.gen())).discriminant()
        delta = delta // delta.gcd(delta.derivative())
        # deal with degree one factors, as RIF and QQ don't go so well
        for fac, e in delta.factor():
            assert e == 1
            if fac.degree() == 1:
                issues = issues.union([RRR(elt[0]) for elt in fac.roots(QQ)])
                delta = delta // fac
        issues = issues.union([elt[0] for elt in delta.roots(RRR)])
        low = sorted(elt for elt in issues if elt <= 0)[-1].ceil().is_int()[1]
        upper = sorted(elt for elt in issues if elt >= 0)[0].floor().is_int()[1]
        return low, upper


def deform_batch(p, q, m, gs, rev=False):
    R = PolynomialRing(ZZ, "x")

    def deform1_tuple(q, m, halfg):
        g = make_weil(q, halfg)
        assert g(1) == m
        l, u = deform(q, trace_polynomial(q, g), 1)
        return (l + m, u + m + 1), m, halfg

    gs = [deform1_tuple(q, m, g) for g in gs]
    gs.sort()
    if rev:
        gs.reverse()
    intervals = None
    for interval, _, halfg in gs:
        if intervals:
            break
        if len(halfg) <= 2:
            continue # we can't try to make it ordinary
        for s in [-1, 1]:
            halfg0 = list(halfg)
            halfg0[-1] += s * (q + 1)
            halfg0[-2] -= s
            halfg0 = tuple(halfg0)
            try:
                interval0, _, _ = deform1_tuple(q, m, halfg0)
                intervals = (interval, interval0)
                halfgs = (halfg, halfg0)
            except NotWeil:
                continue
    else:
        modp = defaultdict(list)
        for interval, _, halfg in gs:
            if len(modp) >= 2:
                break
            modp[halfg[-1] % p].append([interval, halfg])
        else:
            (l, u), _, halfg = gs[0]
            g = make_weil(q, halfg).reverse()
            n = g.degree() // 2
            for k in range(l - m, u - m):
                g0 = g + k * R.gen() ** n
                if not hondatate(g0, p, q):
                    if k < 0:
                        l = k + 1 + m
                    else:
                        u = k + m
                        break
            return (l, u), (m, 1, (halfg,), 0)
        vs = [elt[0] for elt in modp.values()]
        assert len(vs) == 2
        intervals = [elt[0] for elt in vs]
        halfgs = tuple(elt[1] for elt in vs)
    l = max(elt[0] for elt in intervals)
    u = min(elt[1] for elt in intervals)
    return (l, u), (m, 1, halfgs, 0)

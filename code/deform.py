from collections import defaultdict

from sage.all import (
    CBF,
    PolynomialRing,
    QQ,
    RBF,
    RealIntervalField,
    ZZ,
    prime_divisors,
    sqrt,
)

from utils import (
    trace_polynomial,
    is_real_weil,
)

from construct import (
    NotWeil,
)

def deform_middle_ordinary_squarefree(q, S):
    r"""
    Let I_h = {h(1) + c: h(x) + c*x^(degree(h)/2) is a square free Weil polynonial}.

    Given a non-emtpy list S of q-symmetric polynomials such that h(1) = m for h in S,
    we return an interval [l, u] and two polynomials [h0, h1] < S such that:
    - h0, h1 are in S, square free Weil polynomials
    - [l, u] = I_h0 intersection I_h1
    - middle coefficient of h0 and h1 differ modulo p
    where h0 and h1 were picked such a way that minimizes l.
    If there is no such pair, we return None.
    """
    if len(S) < 2:
        return None

    # get the first polynomial of S
    h = S[0]
    m = h(1)
    assert all(elt.is_squarefree() for elt in S)
    assert all(elt.is_weil_polynomial() for elt in S)
    assert all(elt(1) == m for elt in S)

    x = h.parent().gen()
    def compute_interval(g):
        assert g(1) == m
        # we are not interested in deforming these polynomials
        assert g.is_squarefree()
        rg = trace_polynomial(q, g)
        # we only deform Weil polynomials
        assert is_real_weil(q, rg)
        n = g.degree()//2
        l, u = deform(q, rg, 1)
        gl = g + l*x**n
        if not gl.is_squarefree():
            l += 1
            gl = g + l*x**n
            assert gl.is_squarefree()

        gu = g + u*x**n
        if not gu.is_squarefree():
            u -= 1
            gu = g + u*x**n
            assert gu.is_squarefree()
        assert gu.is_weil_polynomial()
        assert gl.is_weil_polynomial()
        return m + l, m + u

    IS = [(elt, compute_interval(elt)) for elt in S]
    modpclasses = defaultdict(list) # b -> [(Ig, g) : g^[n] % p = b}
    p = prime_divisors(q)[0]
    for (g, Ig) in IS:
        modpclasses[g[g.degree()//2] % p].append((Ig, g))

    modp_left = []
    for b, elt in modpclasses.items():
        elt.sort() # the values will be sorted by the left most end point
        modp_left.append(elt[0])
    modp_left.sort()
    if len(modp_left) <= 1:
        return None

    (Ih0, h0), (Ih1, h1) = modp_left[:2]
    l = max(Ih0[0], Ih1[0])
    u = min(Ih0[1], Ih1[1])
    assert all((h + k*x**(h.degree()//2)).is_weil_polynomial() for h in [h0, h1] for k in [l-m, u-m])
    assert all((h + k*x**(h.degree()//2)).is_squarefree() for h in [h0, h1] for k in [l-m, u-m])
    assert (h0[h0.degree()//2] - h1[h1.degree()//2]) % p != 0
    return (l, u), [h0, h1]


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

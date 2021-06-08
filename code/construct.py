from sage.all import (
    CBF,
    CDF,
    PolynomialRing,
    PowerSeriesRing,
    RBF,
    RR,
    RealSet,
    ZZ,
    floor,
    gcd,
    log,
    sqrt,
    srange,
    cached_function,
    prod,
    ComplexBallField,
)
from utils import (
    hat,
    is_real_weil,
    poly_to_tuple,
    make_weil,
)


# we don't always succeed, so this is our custom exception to catch
class NotWeil(Exception):
    """Raised when the polynomial generate is not Weil"""

    pass

ZZx = PolynomialRing(ZZ, 'x')

@cached_function
def roots_away_from_D(q, h):
    hsqrfree = h // gcd(h, h.derivative())
    sqrtinv = RBF(1/sqrt(q))
    if all(elt.abs() > sqrtinv.mid() for elt in hsqrfree.roots(CDF, multiplicities=False)):
        try:
            return all(elt.abs() > sqrtinv for elt in hsqrfree.roots(CBF, multiplicities=False))
        except ValueError:
            return all(elt.abs() > sqrtinv for elt in hsqrfree.roots(ComplexBallField(600), multiplicities=False))


    else:
        return False


@cached_function
def construct_h(q, m, degree_shift=0, adjust_ell=False):
    r"""
    try to construct a range of Weil q-polynomial that cover an interval that contains m
    Input:
    - q
    - m
    - shift helps to adjust the target default degree
    - fixell adapts the target along the way, only useful for small m, as gives us another chance to not fail and not be that far off from the target
    Output:
    - hs, list of polynomials, a single polynomial if we can deform more than one coefficient
    - r, how many coefficients can one deform, i.e., we can deform a_n, ..., a_{n +1 -r}, if r > 0, we only deform by at most floor(q/2)
    - rad, radius of the interval covered
    - [l, u], the closed interval covered
    """
    PS = PowerSeriesRing(RR, "z")
    n = RR(log(m, q) + degree_shift).round()
    if n < 2:
        return None
    ell = log(m / q ** n)
    c = [0] * (n + 1)
    crr = [0] * (n + 1)
    c[1] = RR(ell * q).round()
    crr[1] = RR(c[1])
    s = c[1]
    for i in range(2, n + 1):
        s *= q
        if adjust_ell:
            # sometimes is better to adjust the target along the way
            ell = log(
                m - sum([elt.round() for elt in PS(c).exp(prec=n + 1)[:i]])
            ) - log(q ** n)
        c[i] = ell * q ** i - s
        try:
            crr[i] = RR(c[i])
        except TypeError:
            return None
        cexp = PS(c).exp(prec=i + 1)[i]
        c[i] -= cexp - cexp.round()
        s += c[i]
    h0 = ZZx([elt.round() for elt in PS(c).exp(prec=n + 1)[: n + 1]])
    return construct_h_deform(q, n, m,  h0)

@cached_function
def construct_h_recursive(q, m, a=None):
    if a is None:
        a = floor(sqrt(m))
    b = floor(m/a)
    r = [wrap_construct_h(q, elt, runs=10) for elt in [a,b]]
    if None in r:
        return None
    if any(elt[0][1] != 1 for elt in r):
        return None
    centers = [elt[0][0] for elt in r]
    hs = [elt[0][2][0] for elt in r]
    if hs[0] == hs[1]:
        hs[1] = r[1][0][2][1]

    hs = [list(h) for h in hs]
    hs[0][-1] += a - centers[0]
    hs[1][-1] += b - centers[1]
    h = prod(make_weil(q, tuple(h)) for h in hs)
    assert h(1) == a*b
    n = h.degree()//2
    h0 = h[:n]
    h0 += h[n]*ZZx.gen()**n /2
    return construct_h_deform(q, n, m, h0)



@cached_function
def construct_h_deform(q, n, m, h0):
    from deform import deform
    sqrtqinv = 1 / RBF(q).sqrt()
    x = ZZx.gen()
    E = q ** n * h0(1 / q) + sum(h0) - m
    h = h0 - E / 2 * x ** n  # center h
    center = q ** n * h(1 / q) + sum(h)
    hhat = hat(q, h, n)

    # if at some point we run into trouble trying to perform a big deformation
    # or might just be more worth it to only deform the middle point
    def deform_middle():
        if n <= 2:
            return None
        res = []
        for elt in [h, h0]:
            g0 = hat(q, elt, n)
            t0 = g0.trace_polynomial()[0]
            if not is_real_weil(q, t0):
                continue
            for sign in [1, -1]:
                g1 = hat(q, elt + sign * (x ** (n - 1) - x ** n * (q + 1) / 2), n)
                t1 = g1.trace_polynomial()[0]
                if not is_real_weil(q, t1):
                    continue
                bounds = [deform(q, t, 1) for t in [t0, t1]]
                lower = max(elt[0] for elt in bounds)
                upper = min(elt[1] for elt in bounds)
                halfgs = tuple((poly_to_tuple(elt) for elt in [g0, g1]))
                m = g0(1)
                assert g1(1) == m
                out = (m, 1, halfgs, 0)
                interval = (lower + m, upper + m)
                rad = (upper - lower) // 2
                res.append((out, rad, interval,))
        # False < True
        res.sort(key=lambda elt: (m - interval[0], interval[1] - m, -elt[1]))
        if res:
            return res[-1]
        else:
            # on our first step we didn't manage to create a single Weil polynomial
            return None


    def inverse_hat(g):
        g0 = g[:n]
        g0 += x**n * g[n]/2
        return g0

    if not roots_away_from_D(q, h):
        if roots_away_from_D(q, h0):
            h = h0
            center = q ** n * h0(1 / q) + sum(h0)
        else:
            gs = [hat(q, elt, n) for elt in [h, h0]]
            # disabled at the moment
            if False and  not any(elt.is_weil_polynomial() for elt in gs):
                # try use polynomial with the same initial coefficients
                # this should become expensive about the time past the time the method kicks in
                lead = [(i, 0) for i in h.list()[:-2]]
                Ws = [(sum(elt), elt.reverse()[:n+1]) for elt in ZZx.weil_polynomials(d=2*n, q=q, lead=lead)]
                if not Ws:
                    return None # there are not Weil polynomials anywhere close to what we designed
                h = inverse_hat(Ws[0][1])
                if len(Ws) > 1:
                    h0 = inverse_hat(Ws[1][1])
            return deform_middle()


    # compute min on boundary
    hcirclesqrfree = (h(x) * h(1 / (q * x))).derivative().numerator()
    if hcirclesqrfree != 0:
        hcirclesqrfree = hcirclesqrfree // hcirclesqrfree.gcd(
            hcirclesqrfree.derivative()
        )
        try:
            roots = hcirclesqrfree.roots(CBF, multiplicities=False)
        except ValueError:
            roots = hcirclesqrfree.roots(ComplexBallField(600), multiplicities=False)
        mu = sorted(
            [
                CBF(h(elt)).abs()
                for elt in roots
                if elt.abs().overlaps(sqrtqinv)
            ]
        )[0]
    else:
        # h is constant
        mu = h(sqrtqinv)

    # subtract enough to assure ordinariness
    mu -= sqrtqinv ** (n - 1) + ((q + 1) / 2) * sqrtqinv ** n
    if mu < 0 or mu.contains_zero():
        return deform_middle()
    fq2 = floor(q/2)
    def geometric_sum(a, r):
        # return sum_{i=r}^{n} k^i
        return (a**(1 + n) - a**r)/(a - 1)
    for r in range(n):
        if geometric_sum(RBF(1/sqrt(q)), r) * fq2 < mu:
            break
    if r >= n - 1:
        # deforming just middle gives ~q**(n/2) range
        return deform_middle()
    # we might increase interval by modifying one more coefficient
    rad = floor(floor(q / 2) * ((q ** (n - r + 1) - 1) / (q - 1) + (n - r + 1)))

    # do some checks that we have not done mistakes
    extreme_shift = ZZx([0] * (r) + [floor(q / 2)] * (n + 1 - r))
    for sign in [1, -1]:
        hhat = hat(q, h + sign * extreme_shift, n)
        assert hhat.is_weil_polynomial()
    hhat = hat(q, h, n)
    assert center == hhat(1)
    out = (center, n + 1 - r, (poly_to_tuple(hhat),), mu.upper())
    interval = (center - rad, center + rad)
    return (out, rad, interval,)


def wrap_construct_h(q, m, runs=1000, target_interval=None):
    # try various targets to construct the polynomial
    # especially useful for small m
    # returns None if every attempt failed
    def wrap_return(r):
        if r is not None:
            i = RealSet.closed_open(r[2][0], r[2][1] + 1)
            print(i)
            if target_interval is None:
                if m in i:
                    return r
            else:
                if not target_interval.intersection(i).is_empty():
                    return r
    if m < 6:
        return None
    r = wrap_return(construct_h(q, m))
    if r is not None:
        return r
    r = wrap_return(construct_h_recursive(q, m))
    if r is not None:
        return r

    for k in sorted(srange(-runs, runs + 1, 1), key=lambda elt: ZZ(elt).abs()):
        target = m + k
        print(target)
        if target < 6:
            continue
        for shift in [0, 1, -1]:
            for fixell in [False, True]:
                r = wrap_return(construct_h(q, target, degree_shift=shift, adjust_ell=fixell))
                if r is not None:
                    return r


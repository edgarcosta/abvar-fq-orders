from utils import hat
from deform import deform_middle_ordinary_squarefree

from sage.all import (
    CBF,
    ComplexBallField,
    PolynomialRing,
    PowerSeriesRing,
    RBF,
    RR,
    ZZ,
    cached_function,
    floor,
    gcd,
    log,
    prime_divisors,
)
ZZx = PolynomialRing(ZZ, 'x')

def construction_logexp(q, m, n, check=True):
    r"""
    verbatim copy from Construction 9.1 in the paper
    If successful, returns a Weil polynomial of an ordinary and squarefree n-dimensional
    abelian variety A/F_q such that #A(F_q) = m.
    If the construction fails, returns None
    Input:
    - q, prime power
    - m, the desired order
    - n, the desired dimension
    - check, a boolean, if False, we don't check if the returned polynomial is a Weil polynomial
    """
    PS = PowerSeriesRing(RR, "z")
    ell = log(m / q ** n)
    c = [0] * (n + 1)
    crr = [0] * (n + 1)
    c[1] = RR(ell * q).round()
    crr[1] = RR(c[1])
    s = c[1]

    # compute c1, c2, ..., c_{n-1}
    for i in range(2, n):
        s *= q
        c[i] = ell * q ** i - s
        try:
            crr[i] = RR(c[i])
        except TypeError: # c[i] isn't real
            return None
        cexp = PS(c).exp(prec=i + 1)[i]
        c[i] -= cexp - cexp.round()
        s += c[i]


    # c_n only affects the coefficient x^n of h0, and it will be overwritten
    # when defining h_1.
    # h0truncation in ZZ[x] is a degree n-1 polynomial such that
    # h0truncation = h0 mod x^n
    h0truncation = ZZx([elt.round() for elt in PS(c).exp(prec=n)[:n]])

    # compute h1
    x = ZZx.gen()
    shift =  q ** n * h0truncation(1 / q) + h0truncation(1) - m
    h1 = h0truncation - shift / 2 * x ** n

    # compute h
    p = prime_divisors(q)[0]
    h = h1
    if ZZ(2*h1[n]) % p == 0:
        h += x**(n-1) - ((q + 1)/2)*x**n
    hhat = hat(q, h, n)
    assert hhat[n] % p != 0
    assert hhat(1) == m
    hhat = ZZx(hhat)
    if not check or hhat.is_weil_polynomial():
        return hhat
    else:
        return None

def construction_hat_deform(q, m, n=None):
    r"""
    verbatim copy from Construction 9.5 in the paper
    Input:
    - q, prime power
    - m, the desired order
    - n, the desired dimension, if not given we take n = floor(log(m)/log(q)) + 1
    If successful, returns a tuple:
    - h(x) = x^2n + a_1 * x^(2n - 1) ... + q^n, a Weil polynomial of an ordinary and squarefree n-dimensional abelian variety A/F_q such that #A(F_q) = m.
    - r, how many coefficients, i.e., we can deform a_n, ..., a_{n +1 -r}, can one deform by at most floor(q/2)
    - [l, u], the closed interval covered by the deformation
    - a lower bound for mu_ord
    If the construction fails, returns None
    """
    if n is None:
        n = floor(log(m)/log(q)) + 1

    def inverse_hat(g):
        """
        returns a polynomial g0 of degree n, such that hat(g0) = g,
        where g is the input and a polynomial of degree 2n
        """
        n = g.degree()/2
        x = g.parent().gen()
        g0 = g.reverse()[:n]
        g0 += x**n * g[n]/2
        assert hat(q, g0, n) == g
        return g0

    sqrtqinv = 1 / RBF(q).sqrt()
    # Step 1: Call the logexp construction
    hhat = construction_logexp(q, m, n)
    if hhat is None: # i.e., the logexp construction failed
        return None

    h = inverse_hat(hhat)

    # Step 2: check that h has all the roots away from D
    # we use the square free version of h for a better convergence of the numerical methods
    h_squarefree = h // gcd(h, h.derivative())
    if not all(elt.abs() > sqrtqinv for elt in h_squarefree.roots(CBF, multiplicities=False)):
        return None

    # Step 3: compute min on boundary
    # i.e., compute the complex zeros of the derivative of h(x)*h(1/(qx))
    x = ZZx.gen()
    h_circle = (h(x) * h(1 / (q * x))).derivative().numerator()
    if h_circle != 0:
        # for better convergence of the numerical methods, make it square free
        h_circle_squarefree = h_circle // h_circle.gcd(h_circle.derivative())
        try:
            roots = h_circle_squarefree.roots(CBF, multiplicities=False)
        except ValueError:
            # increase precision so that the roots get isolated.
            # only required for q > 11 and m >> 1
            roots = h_circle_squarefree.roots(ComplexBallField(600), multiplicities=False)
        # mu = smalles |h(alpha)| for alpha in roots and |alpha| = 1/sqrt(q)
        mu = sorted(
            [
                CBF(h(alpha)).abs() # the CBF is needed in case we computed the roots with higher precision
                for alpha in roots
                if alpha.abs().overlaps(sqrtqinv) # i.e., as a ball |alpha| intersects the ball defining 1/sqrt(q)
            ]
        )[0]
    else:
        # h is constant
        mu = h(sqrtqinv).abs()
        assert h == 1

    # Step 4: subtract enough to assure ordinariness
    mu_ord = mu -  sqrtqinv ** (n - 1) - ((q + 1) / 2) * sqrtqinv ** n

    if mu_ord < 0 or mu_ord.contains_zero(): # mu_ord is negative, hence, r doesn't exist
        return None

    # Step 5: compute r
    fq2 = floor(q/2)
    def geometric_sum(a, r):
        # return sum_{i=r}^{n} k^i
        return (a**(1 + n) - a**r)/(a - 1)

    for r in range(1, n + 2):
        if geometric_sum(sqrtqinv, r) * fq2 < mu_ord:
            break
    else:
        assert False # r <= n + 1

    # Step 6: compute N
    N = floor(q / 2) * ((q ** (n - r + 1) - 1) / (q - 1) + (n - r + 1))
    # for our purposes we are only interested in floor(N)
    rad = floor(N)
    interval = (m - rad, m + rad)
    return interval, [ZZx(hat(q, h, n))], n + 1 - r, mu_ord.lower(), 'hat'



def construction_logexp_then_exhaust_deform(q, m, n, try_exhaust=True):
    r"""
    An extended version of Construction 9.1 in the paper.
    If successful, returns an interval [l,u] and two Weil polynomial h0, h1 of degree 2n
    such that hi(1) = m, with the property such that for k in [l, u] for at least one hi,
    hi(x) + c*x^n is a Weil polyomial of an ordinary and squarefree n-dimensional abelian
    variety A/F_q and #A(F_q) = k.
    If the construction fails, returns None

    Let I_h = {h(1) + c: h(x) + c*x^n is a square free Weil polynonial}.
    Let g be the polynomial returned construction_logexp(q, m, n, check=False), i.e.,
    g might not be a Weil polynomial of an ordinary abelian variety.
    Let S0 = {g, g^-, g^+},
    where g^+/- =  g(x) +/- (qx^{n-1} - (q+1)x^n + x^{n+1}).
    Let S1 = {h in S0 : h is Weil polynomial and square free}
    Try [l, u], [h0, h1]  = deform_middle_ordinary_squarefree(q, S1), if succeeds
    return [l, u], [h0, h1], 'logexp'.

    If g is not a Weil polynomial and try_exhaust=True, let
    S1 = {h: h(1) = m, h is square free, h is a Weil polynomial}
    Try [l, u], [h0, h1]  = deform_middle_ordinary_squarefree(q, S1), if succeeds
    return [l, u], [h0, h1], 'exhaust'.

    Input:
    - q, prime power
    - m, the desired order
    - n, the desired dimension
    Output:
    - [l, u]de
    - [h0, h1]
    - 1 # r = 1, we only deform 1 coefficient
    - 0 # mu_ord
    - s, s in ['logexp', exhaust']
    """
    # Step 1: Call the logexp construction
    h = construction_logexp(q, m, n, check=False)
    if h is None: # i.e., the logexp construction failed
        return None
    x = ZZx.gen()
    S = [h] + [h  + s* (q*x**(n-1) - (q+1)*x**n + x**(n+1)) for s in [1, -1]]
    # only keep the square free polynomials and Weil polynomials
    S = [elt  for elt in S if elt.is_squarefree() and elt.is_weil_polynomial()]
    if len(S) > 1:
        interval, hs = deform_middle_ordinary_squarefree(q, S)
        return interval, hs, 1, 0, 'logexp'

    if try_exhaust:
        S = [g for g in ZZx.weil_polynomials(2*n, q, lead=tuple(h.reverse().list()[:n-1])) if g(1) == m and g.is_squarefree()]
        out = deform_middle_ordinary_squarefree(q, S)
        if out is None:
            return None
        interval, hs = out
        return interval, hs, 1, 0, 'exhaust'

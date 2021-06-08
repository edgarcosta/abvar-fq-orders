from sage.all import (
    PolynomialRing,
    QQ,
    RealSet,
    ZZ,
    ceil,
    floor,
    sqrt,
    parallel,
    lcm,
    NumberField,
)
from itertools import islice
from multiprocessing import cpu_count


def hat(q, h, n):
    # the hat operator
    # expects h in 1 + z ZZ[x]
    # produces a monic polynomial
    assert h[0] == 1
    RQ = PolynomialRing(QQ, "x")  # only used to baby sit sage
    x = h.parent().gens()[0]
    return h.parent()(RQ(x ** (2 * n) * h(1 / x) + q ** n * h(x / q)))


def make_weil(q, l):
    R = PolynomialRing(ZZ, "x")
    # given (1, a1, ..., an) creates the q-symmetric polynomial in 1 + z ZZ[x]
    assert l[0] == 1
    return R(list(l + tuple(q ** i * elt for i, elt in enumerate(reversed(l[:-1]), 1))))


def poly_to_tuple(h):
    # the inverse of make_weil
    if h.is_monic():
        return tuple(map(int, h.reverse().list()[: h.degree() // 2 + 1]))
    else:
        return tuple(map(int, h.list()[: h.degree() // 2 + 1]))


def trace_polynomial(q, f):
    if not f.is_monic():
        assert f[0] == 1
        f = f.reverse()
    x = f.variables()[0]
    coeffs = []
    m = f.degree() // 2
    for i in reversed(range(m + 1)):
        coeffs.append(f[2 * i])  # Note: degree of f may be less than 2*i
        f = (f % (x ** 2 + q) ** i) // x
    coeffs.reverse()
    return f.parent()(coeffs)


def is_real_weil(q, f):
    return f.is_monic() and f.all_roots_in_interval(-2 * q.sqrt(), 2 * q.sqrt())


def chunkify(l, size):
    return [islice(l, elt*size, (elt + 1)*size) for elt in range(len(l)//size + 1) if elt*size < len(l)]


def hw(q, n):
    if n == 0:
        return RealSet.closed(1,1)
    return RealSet.closed(ceil((sqrt(q) - 1)**(2*n)), floor((sqrt(q) + 1)**(2*n)))


def hondatate_simple(f, p, q):
    if f[f.degree()//2] % p != 0: #ordinary
        return 1
    K = NumberField(f, 'a')
    a = f.roots(K)[0][0]
    inv = [(elt - elt.floor()).denominator()
           for elt in
           [
               (a.valuation(v) / K(q).valuation(v)) * v.residue_class_degree() * v.ramification_index()
                for v in K.primes_above(p)
           ]
          ]
    return lcm(inv)

def hondatate(f, p, q):
    for fac, e in f.factor():
        if e % hondatate_simple(fac, p, q) != 0:
            return False
    return True

@parallel(ncpus=cpu_count())
def hondatatefilter(fs, p, q):
    return [f for f in fs if hondatate(f, p, q)]

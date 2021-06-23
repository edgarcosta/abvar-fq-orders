from multiprocessing import cpu_count
from collections import defaultdict
from time import time
import pathlib
import os
try:
    current_directory = pathlib.Path(__file__).parent.absolute()
except NameError:
    current_directory = pathlib.Path(os.path.abspath(os.curdir))
# optional: https://github.com/kedlaya/root-unitary/
# and to make use of parallelism one needs to enable OpenMP support
# see the first lines of weil_polynomials.pyx of the root-unitary package
root_unitary_directory = pathlib.PurePath(current_directory, '..', 'root-unitary/')
root_unitary_directory = pathlib.PurePath('/home/edgarcosta/', 'root-unitary/')
os.environ["SAGE_NUM_THREADS"] = os.environ["OMP_NUM_THREADS"] = str(cpu_count())

from sage.all import (
    cached_function,
    ceil,
    load,
    prime_divisors,
    save,
)
from sage.rings.polynomial.weil.weil_polynomials import WeilPolynomials
from utils import (
    chunkify,
    hondatatefilter,
    poly_to_tuple,
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
            W = WeilPolynomials(2*g, q)
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

def seed_small(q, n):
    res = defaultdict(set)
    for g in range(1, n+1):
        for elt in weil_poly_gq(g, q):
            elt1 = int(sum(elt.list()))
            res[elt1].add(poly_to_tuple(elt))
    return res

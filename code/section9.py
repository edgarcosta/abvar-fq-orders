from cover import cover
def section_9_check(q):
    params = {
        3: (0, 32),
        4: (0, 25),
        5: (0, 24),
        7: (0, 26),
        8: (18, 26),
        9: (20, 28),
        11: (24, 30),
        13: (28, 34),
    }
    b, e = params[q]
    filename = f'../data/q{q}_{q}e{b}-{q}e{e}.txt'
    return cover(q, q**e, start=q**b, filename=filename)

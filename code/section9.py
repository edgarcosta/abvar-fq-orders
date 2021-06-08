from cover import cover_and_save
def section_9_check(q):
    params = {
        3: (None, 32),
        4: (None, 25),
        5: (None, 24),
        7: (None, 26),
        8: (18, 26),
        9: (20, 28),
        11: (24, 30),
        13: (28, 34),
    }
    begin, end = params[q]
    coverage = None
    res = None
    filename = '../data/q{q}_1-{q}e{end}.txt'.format(q=q, end=end)
    step = 1000
    if begin is not None:
        # a way to set a starting point that is not zero
        coverage = coverage=RealSet.closed(q**begin-1, q**begin-1)
        res={}
        filename = '../data/q{q}_{q}e{begin}-{q}e{end}.txt'.format(q=q, begin=begin, end=end)
        step=q**(begin-2)
    return cover_and_save(filename, q=q,
                          end=q**end,
                          coverage=coverage, 
                          res=res)


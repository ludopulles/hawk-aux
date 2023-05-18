from mpmath import mpf, mp, exp, log, floor, ceil, nstr, sqrt

mp.prec = 1000

# Parameters of the schemes HAWK-512 and HAWK-1024,
# where you need to take logn = 9 and logn = 10 respectively.
# HAWK-512 has a NIST-1 security level
# HAWK-1024 has a NIST-5 security level
sigma_sig = { 8: '1.010', 9: '1.278', 10: '1.299' }
sigma_ver = { 8: 1.042, 9: 1.425, 10: 1.571 }
sigma_sec = { 8: 1.042, 9: 1.425, 10: 1.974 }


def rho(x, sigma):
    x = mpf(x) / mpf(sigma)
    return exp(x * x / mpf('-2'))


def RenyiDivergence(P, Q, a):
    """
    Compute the Renyi divergence between two distributions at order a.
    """
    # Make sure that support(P) is a subset of support(Q)
    assert set(P.keys()).issubset(set(Q.keys()))

    return sum(P[x] * ((P[x] / Q[x])**mpf(a - 1)) for x in P)**(mpf('1') / mpf(a - 1))


def report_RD(txt, P, Q, a):
    RD = RenyiDivergence(P, Q, a)
    minlog2 = floor(-log(RD - 1) / log(2))
    print(f"// RD_{{%d}}({txt}) = 1 + %E < 1 + 2^-%d" % (a, RD - 1, minlog2))


# Raw table of all the tabulated probabilities.
# Note the first value is P(X = 0), while the $k$th one is P(X >= k+1 | X >=
# 1) for k >= 1.
def produce_gaussian(sigma, center, zs=100):
    L = mpf(sigma) * zs
    supp = range(int(floor(center - L, prec=0)), int(ceil(center + L, prec=0)) + 1)

    norm = mpf('0')
    for i in supp:
        norm += rho(i - center, sigma)

    Q = {}
    for i in supp:
        Q[i] = rho(i - center, sigma) / norm
    return Q


# Assumes the center of the distribution is either 0 or 1/2
def dist_to_CDT(dist, prec, double_mu):
    L = max([ -min(dist.keys()), max(dist.keys()) ])
    scale = mpf('2')**prec

    if double_mu == 0:
        # Entry 0: P( |X| >= 1 )
        # Entry 1: P( |X| >= 2 )
        # ...
        cdt = [ sum(dist[y] for y in dist if y >= x or y <= -x) * scale for x in range(1, L) ]
    else:
        # Entry 0: P( X >= 2 or X <= -1 )
        # Entry 1: P( X >= 3 or X <= -2 )
        # ...
        cdt = [ sum(dist[y] for y in dist if y >= x or y <= 1-x) * scale for x in range(2, L) ]

    cdt = list(map(lambda x: int(floor(x, prec=0)), cdt))
    while cdt[-1] == 0:
        cdt.pop()
    return cdt


def CDT_to_dist(cdt, double_mu, prec):
    """
    Converts a CDT (cumulative denstiy table) with precision 'prec' bits to a
    probability distribution, exploiting the symmetry around center mu, i.e.
        P(X = -x) = P(X = x + double_mu),
    where `double_mu \in {0, 1}.
    """
    scale = mpf(2)**prec
    res = {}
    for i in range(len(cdt) + 1):
        prevValue = scale if i == 0 else cdt[i - 1]
        nextValue = 0 if i == len(cdt) else cdt[i]
        res[i + double_mu] = mpf(prevValue - nextValue) / (mpf('2') * scale)
        if double_mu == 0 and i == 0:
            res[0] *= mpf('2')
        else:
            res[-i] = res[i + double_mu]
    return res


def print_sign_tables(logn, tables):
    """
    Print the lookup cumulative probability table used in the signing algorithm
    of hawk, using two centers of the Gaussian distribution: 0 and 1/2.
    The support is Z and the standard deviation is 1.278.
    The precision used will be 79 bits (and one bit for sign, so 80 total)
    We store the high 17 bits in a separate table and the low 63 bits as well
    somewhere.
    Then, to compare a randomly sampled bitstring of 80 bits, we do two
    comparisons
    """
    L = len(tables[0])
    while len(tables[1]) < L:
        tables[1].append(0)
    assert L == len(tables[1])

    print("static const uint16_t sig_gauss_hi_Hawk_%d[] = {" % (1 << logn))
    for i in range(L):
        hi0, hi1 = tables[0][i] >> 63, tables[1][i] >> 63
        if hi0 == 0 and hi1 == 0:
            break
        print(f"    0x{hi0:04X}, 0x{hi1:04X}, ")
    print("};")

    print("static const uint64_t sig_gauss_lo_Hawk_%d[] = {" % (1 << logn))
    for i in range(L):
        lo0, lo1 = tables[0][i] & 0x7FFFFFFFFFFFFFFF, tables[1][i] & 0x7FFFFFFFFFFFFFFF
        print(f"    0x{lo0:016X}, 0x{lo1:016X}, ")
    print("};")

    # print(" i, T_0,                    T_1:")
    # for i in range(L):
    #     print(f"{i:2}, 0x{tables[0][i]:020X}, 0x{tables[1][i]:020X}")
    # print()


def sec_loss_prac_ideal(n, lam, q_s, P0, Q0, P1, Q1, a_table):
    eps = 1.0 / sqrt(q_s*lam)

    # Renyi divergence for using numerically approximated values in the tables, instead of the actual real values.
    tables_renyi0 = RenyiDivergence(P0, Q0, a_table)
    tables_renyi1 = RenyiDivergence(P1, Q1, a_table)
    # tables_renyi = max(tables_renyi0, tables_renyi1)**(2*n)
    tables_renyi = max(tables_renyi0, tables_renyi1)**(2 * n)

    tables_renyi = tables_renyi**q_s

    # eps_{realattack}^{a/(a-1) - 1},
    # where eps_{realattack} >= 2^-lambda
    advantage_lowerbound = mpf('2')**(-lam / (a_table - 1))

    # eps_{idealattack} >= eps_{midattack}^{ahash/(ahash-1)} / R_{ahash}( midDist || idealDist )
    return tables_renyi / advantage_lowerbound


def optimise_sec_prac_ideal(logn, lam, prec):
    Q0 = produce_gaussian(sigma_sig[logn], mpf('0'))
    Q1 = produce_gaussian(sigma_sig[logn], 1 / mpf('2'))
    P0 = CDT_to_dist(dist_to_CDT(Q0, prec, 0), 0, prec)
    P1 = CDT_to_dist(dist_to_CDT(Q1, prec, 1), 1, prec)

    f = lambda a: sec_loss_prac_ideal(2**logn, lam, 2**64, P0, Q0, P1, Q1, a)
    alo, ahi = 2, 32 * lam + 1
    for it in range(50):
        a_l = (2 * alo + 1 * ahi) / 3.0
        a_r = (1 * alo + 2 * ahi) / 3.0
        loss_l = f(a_l)
        loss_r = f(a_r)

        if loss_l < loss_r:
            ahi = a_r
        else:
            alo = a_l

    a = 0.5 * (alo + ahi)
    return a, f(a)


def security_loss_hawk(n, lam, q_s, P0, Q0, P1, Q1, a_table, a_hash):
    eps = 1.0 / sqrt(q_s*lam)

    # Renyi divergence for using numerically approximated values in the tables, instead of the actual real values.
    tables_renyi0 = RenyiDivergence(P0, Q0, a_table)
    tables_renyi1 = RenyiDivergence(P1, Q1, a_table)
    tables_renyi = max(tables_renyi0, tables_renyi1)**(2*n)

    # Set up the Reverse Pinsker Inequalities to convert the Delta between sampling uniform hashes vs just short vectors (Lemma 9, HAWK-AC22):
    delta_lemma9 = eps / (mpf('1') - eps)

    ess_inf = (1 - eps) / (1 + eps) # m < 1
    ess_sup = (1 + eps) / (1 - eps) # M > 1

    # Renyi divergence for one sample where one is from first hash then short vector in lattice, and other is from "short vector in whole lattice".
    DvsD_Q = mpf('1') + delta_lemma9 * ((ess_sup**a_hash - 1)/(ess_sup-1) - (1 - ess_inf**a_hash)/(1 - ess_inf))
    DvsD_Q = DvsD_Q**(1 / mpf(a_hash - 1))

    prod_renyi = DvsD_Q * tables_renyi**(a_hash / (a_hash - 1))
    prod_renyi = prod_renyi**q_s

    # eps_{realattack}^{ab/(a-1)(b-1) - 1},
    # where eps_{realattack} >= 2^-lambda
    advantage_lowerbound = mpf('2')**(-lam * (a_table*a_hash/(a_table-1)/(a_hash-1) - 1))

    # eps_{idealattack} >= eps_{midattack}^{ahash/(ahash-1)} / R_{ahash}( midDist || idealDist )
    # eps_{midattack} >= eps_{realattack}^{atable/(atable-1)} / R_{atable}( realDist || midDist )
    return prod_renyi / advantage_lowerbound


def optimise_security_loss(logn, lam, prec):
    Q0 = produce_gaussian(sigma_sig[logn], 0)
    Q1 = produce_gaussian(sigma_sig[logn], 0.5)
    P0 = CDT_to_dist(dist_to_CDT(Q0, prec, 0), 0, prec)
    P1 = CDT_to_dist(dist_to_CDT(Q1, prec, 1), 1, prec)

    alo, ahi = 2, 4 * lam + 1
    for it in range(25):
        a_l = (2 * alo + 1 * ahi) / 3.0
        a_r = (1 * alo + 2 * ahi) / 3.0
        loss_l = security_loss_hawk(2**logn, lam, 2**64, P0, Q0, P1, Q1, a_l, a_l)
        loss_r = security_loss_hawk(2**logn, lam, 2**64, P0, Q0, P1, Q1, a_r, a_r)

        if loss_l < loss_r:
            ahi = a_r
        else:
            alo = a_l

    a = 0.5 * (alo + ahi)
    return a, security_loss_hawk(2**logn, lam, 2**64, P0, Q0, P1, Q1, a, a)


# Since R_a(P||Q) is a non-decreasing function as a --> oo, we can take a = 513
# in any case even though a = 257 is sufficient for NIST-1.
prec = 78

# Report on sigma_{sign}
# Hawk sampling during signing
for logn in [8, 9, 10]: # NIST-1 and NIST-5 respectively.
    # required Renyi divergence for signing should be < req_renyi
    req_renyi = 1 + mpf('2')**(-(64 + 2 + (1 + logn)))
    lam = 2**(logn - 2) # 128, 256

    tables = []
    print(f"\n// Precomputed CDT with {prec} bits of precision.")
    for double_mu in [0, 1]:
        Q = produce_gaussian(sigma_sig[logn], double_mu / mpf('2'))
        CDT = dist_to_CDT(Q, prec, double_mu)
        P = CDT_to_dist(CDT, double_mu, prec)
        report_RD(f"T_{double_mu} || T_{double_mu} ideal", P, Q, 2*lam + 1)
        assert RenyiDivergence(P, Q, 2*lam + 1) < req_renyi
        tables.append(CDT)
    print_sign_tables(logn, tables)

    a_opt, sec_loss = optimise_sec_prac_ideal(logn, lam, prec)
    print(f"Security loss (practical HAWK -> ideal HAWK): {nstr(sec_loss, 10)} at order a={a_opt}")
    a_opt, sec_loss = optimise_security_loss(logn, lam, prec)
    print(f"Security loss: {nstr(sec_loss, 10)} at order a={a_opt}")
print()

# Report on lower bound of ||(f, g)||^2 in KeyGen:
l2lows = {logn: 1 + int(floor(sig**2 * (2<<logn))) for (logn, sig) in sigma_sec.items()}
print("l2low (ng_hawk.c:L596): ", ", ".join(map(str, l2lows.values())),
      "(logn =", ", ".join(map(str, l2lows.keys())), "respectively).")

# Report on verification upper bound of squared distance from 2s to a target h.

maxnorms = {logn: int(floor((2 * sig)**2 * (2 << logn))) for (logn, sig) in sigma_ver.items()}
print("max_[x|t]norm (hawk_sign.c:L1547, hawk_vrfy.c:L2933): ",
      ", ".join(map(str, maxnorms.values())),
      "(logn =", ", ".join(map(str, maxnorms.keys())), "respectively).")
print()

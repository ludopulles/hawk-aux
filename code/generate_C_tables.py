from mpmath import mpf, mp, exp, log, floor, ceil, nstr, sqrt

mp.prec = 1000


def rho(x, sigma):
    """
    Give gaussian weight at point `x` with standard deviation `sigma`.
    """
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
    minlog2 = floor(-log(RD - 1, 2))
    print(f"// RD_{{%d}}({txt}) = 1 + %E < 1 + 2^-%d" % (a, RD - 1, minlog2))


# Raw table of all the tabulated probabilities.
# Note the first value is P(X = 0), while the $k$th one is P(X >= k+1 | X >=
# 1) for k >= 1.
def produce_gaussian(sigma, center, zs=100):
    """
    Give a discrete gaussian distribution, D_{Z, center, sigma}
    """
    width = mpf(sigma) * zs
    supp = range(int(floor(center - width, prec=0)), int(ceil(center + width, prec=0)) + 1)

    norm = mpf('0')
    for i in supp:
        norm += rho(i - center, sigma)

    dist = {}
    for i in supp:
        dist[i] = rho(i - center, sigma) / norm
    return dist


# Assumes the center of the distribution is either 0 or 1/2
def dist_to_CDT(dist, prec, double_mu):
    """
    Provide a cumulative density table from a given distribution with its 'center' at double_mu/2
    :param dist: the distribution
    :param prec: the precision in bits used to represent the entries in the tables as integers, as
    we scale up the entries by 2^prec.
    :param double_mu: double_mu/2 marks the 'middle' of the distribution which is assumed to be
    symmetric.
    """
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

    if max(tables[0][0], tables[1][0]) > (1 << 63):
        print(f"static const uint16_t sig_gauss_hi_Hawk_{1 << logn}[] = {{")
        for i in range(L):
            hi0, hi1 = tables[0][i] >> 63, tables[1][i] >> 63
            if hi0 == 0 and hi1 == 0:
                break
            print(f"    0x{hi0:04X}, 0x{hi1:04X}, ")
        print("};")

    print(f"static const uint64_t sig_gauss_lo_Hawk_{1 << logn}[] = {{")
    for i in range(L):
        lo0, lo1 = tables[0][i] & 0x7FFFFFFFFFFFFFFF, tables[1][i] & 0x7FFFFFFFFFFFFFFF
        print(f"    0x{lo0:016X}, 0x{lo1:016X}, ")
    print("};")

    # print(" i, T_0,                    T_1:")
    # for i in range(L):
    #     print(f"{i:2}, 0x{tables[0][i]:020X}, 0x{tables[1][i]:020X}")
    # print()


def security_loss_tables(n, lam, q_s, P0, Q0, P1, Q1, a):
    """
    Give the security loss factor caused by using the table based methods instead of a perfect
    discrete gaussian distribution.
    :param n: hawk parameter
    :param lam: lambda (number of security bits)
    :param q_s: the number of allowed queries an adversary can make (2^64 for all NIST sec. levels)
    :param P0: table based probability distribution when centered at 0.
    :param Q0: ideal probability distribution when centered at 0.
    :param P1: table based probability distribution when centered at 1/2.
    :param Q1: ideal probability distribution when centered at 1/2.
    :param a: order used in the Renyi divergence.
    """
    # Renyi divergence for using numerically approximated values in the tables, instead of the
    # actual real values.
    renyi_div0 = RenyiDivergence(P0, Q0, a)
    renyi_div1 = RenyiDivergence(P1, Q1, a)

    # There are 2n samples taken from center=0 or center=1/2 distribution, for one signature.
    # Then, at most q_s of these signatures may be generated.
    # Hence, the adversary sees 2n*q_s different samples from either of these distributions.
    renyi_div = max(renyi_div0, renyi_div1)**(2 * n * q_s)

    # We may assume WLOG that an attacker A against the table based method "wins" with probability at
    # least 2^{-lambda}. Now, we want to know what the probability of A winning against the ideal
    # distribution is.

    # [Prest17] (https://tprest.github.io/pdf/pub/renyi.pdf) equation (2) says:
    # Q(E) >= P(E)^{a/(a-1)} / R_a(P || Q),
    # where:
    # - E is the event of A producing a forgery,
    # - P is the distribution based on the tables,
    # - Q is the distribution based on the ideal discrete gaussian.

    # P(E) >= 2^-lambda implies now:
    # Q(E) >= P(E) P(E)^{1/(a-1)} / R_a(P || Q) >= P(E) / (2^{lambda/(a-1)} R_a(P || Q)).
    return mpf('2')**(lam / (a - 1)) * renyi_div


def optimise_security_loss_tables(logn, lam, q_s, sig, prec):
    Q0 = produce_gaussian(sig, mpf('0'))
    Q1 = produce_gaussian(sig, 1 / mpf('2'))
    P0 = CDT_to_dist(dist_to_CDT(Q0, prec, 0), 0, prec)
    P1 = CDT_to_dist(dist_to_CDT(Q1, prec, 1), 1, prec)

    def f(a):
        return security_loss_tables(2**logn, lam, q_s, P0, Q0, P1, Q1, a)
    alo, ahi = 2, 1 << 20
    # This performs a ternary to find the optimal value, assuming the function f is (strictly)
    # decreasing up to the optimum, and (strictly) increasing after that optimum.
    for _ in range(50):
        a_l = (2 * alo + ahi) / 3.0
        a_r = (alo + 2 * ahi) / 3.0
        loss_l, loss_r = f(a_l), f(a_r)

        if loss_l < loss_r:
            ahi = a_r
        else:
            alo = a_l

    a_opt = round((alo + ahi) / 2)
    return a_opt, f(a_opt)


def rho(x, s):
    # Gaussian weight at x:
    return exp(mpf(x*x) / -2 / (s*s))


def rho_zz(x, s, zs=100):
    # Gaussian weight of Z + x:
    return sum(rho(x + y, s) for y in range(-zs, zs + 1))


def security_loss_cosets(n, lam, q_s, sig, a):
    """
    Compute security loss factor in a reduction from using uniformly random cosets to sample short
    vectors to sampling short vectors directly.
    :param n: hawk parameter
    :param lam: lambda (number of security bits)
    :param q_s: the number of allowed queries an adversary can make (2^64 for all NIST sec. levels)
    :param sig: sigma used in signing.
    :param a: order used in the Renyi divergence.
    """
    alpha = rho_zz(0, 2 * sig) / rho_zz(0,          sig) / mpf('2')
    beta  = rho_zz(0, 2 * sig) / rho_zz(1/mpf('2'), sig) / mpf('2')
    exponent = (2 * n / (a - 1)) * q_s
    # RD_a = ((alpha^{a-1} + beta^{a-1}) / 2)^{2n / (a - 1)}
    renyi_div = ((alpha**(a - 1) + beta**(a - 1)) / mpf('2'))**exponent

    # [Prest17] (https://tprest.github.io/pdf/pub/renyi.pdf) equation (2) says:
    # Q(E) >= P(E)^{a/(a-1)} / R_a(P || Q),
    # where:
    # - E is the event of A producing a forgery,
    # - P is the distribution based on playing omSVP where samples are generated by first getting a
    # hash h <- {0,1}^d uniformly at random and then sampling from D_{Q, Z^d + h, 2sigma_sig}.
    # - Q is the distribution based on playing omSVP where samples are generated by directly
    # generating a vector in D_{Q, Z^d, 2sigma_sig}.

    # P(E) >= 2^-lambda implies now:
    # Q(E) >= P(E) P(E)^{1/(a-1)} / R_a(P || Q) >= P(E) / (2^{lambda/(a-1)} R_a(P || Q)).
    return mpf('2')**(lam / (a - 1)) * renyi_div

def optimise_security_loss_cosets(logn, lam, q_s, sig):
    alo, ahi = 2, 1 << 20
    def f(a):
        return security_loss_cosets(2**logn, lam, q_s, sig, a)
    for _ in range(50):
        a_l, a_r = (2 * alo + ahi) / 3.0, (alo + 2 * ahi) / 3.0
        loss_l, loss_r = f(a_l), f(a_r)

        if loss_l < loss_r:
            ahi = a_r
        else:
            alo = a_l

    a_opt = round((alo + ahi) / 2)
    return a_opt, f(a_opt)


def __main__():
    """
    Compute cumulative probability tables used in HAWK's signature generation.
    Also report on the security loss caused by the approximation vs the perfect discrete gaussian
    distribution.
    """
    # Parameters of the schemes HAWK-512 and HAWK-1024,
    # where you need to take logn = 9 and logn = 10 respectively.
    # HAWK-512 has a NIST-1 security level
    # HAWK-1024 has a NIST-5 security level
    sigma_sig = { 8: mpf('1.010'), 9: mpf('1.278'), 10: mpf('1.299') }
    sigma_ver = { 8: 1.042, 9: 1.425, 10: 1.571 }
    sigma_sec = { 8: 1.042, 9: 1.425, 10: 1.974 }

    # Report on sigma_{sign}
    # Hawk sampling during signing
    for logn in [8, 9, 10]: # NIST-1 and NIST-5 respectively.
        # Note: HAWK-256 could use less precision, because q_s = 2^32.
        # prec = 63 if logn == 8 else 78
        prec = 78
        q_s = 2**32 if logn == 8 else 2**64

        # required Renyi divergence for signing should be < req_renyi
        req_renyi = 1 + mpf('2')**(-(2 + (1 + logn))) / q_s
        lam = 2**(logn - 2) # 128, 256

        tables = []
        print(f"// Precomputed CDT with {prec} bits of precision.")
        for double_mu in [0, 1]:
            Q = produce_gaussian(sigma_sig[logn], double_mu / mpf('2'))
            cdt = dist_to_CDT(Q, prec, double_mu)
            P = CDT_to_dist(cdt, double_mu, prec)
            report_RD(f"T_{double_mu} || T_{double_mu} ideal", P, Q, 2*lam + 1)
            assert RenyiDivergence(P, Q, 2*lam + 1) < req_renyi
            tables.append(cdt)
        print_sign_tables(logn, tables)

        a_opt, loss = optimise_security_loss_tables(logn, lam, q_s, sigma_sig[logn], prec)
        print(f'Security loss (table based -> "ideal" HAWK): {nstr(loss, 10)} at order a={a_opt}')
        a_opt, loss = optimise_security_loss_cosets(logn, lam, q_s, sigma_sig[logn])
        print(f"Security loss (sample in uniform coset 2Z^d + h -> sample in Z^d): "
              f"{nstr(loss, 10)} at order a={a_opt}\n")

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


if __name__ == "__main__":
    __main__()

RR = RealField(1000)
eps_zero = RR(1.0) / RR(1 << 100)

# Parameters of the schemes HAWK-512 and HAWK-1024,
# where you need to take logn = 9 and logn = 10 respectively.
# HAWK-512 has a NIST-1 security level
# HAWK-1024 has a NIST-5 security level
sigma_sig = { 9: 1.278, 10: 1.299 }
sigma_ver = { 9: 1.425, 10: 1.572 }
sigma_sec = { 9: 1.425, 10: 1.974 }
sigma_pk  = { 9: 1.500, 10: 2.000 }

def eq(a, b):
    return abs(a - b) < eps_zero

def rho(x, sigma):
    x = RR(x) / RR(sigma)
    return (RR(-0.5) * x * x).exp()

def RenyiDivergence(P, Q, a):
    """
    Compute the Renyi divergence between two distributions at order a.
    """
    # Make sure that support(P) is a subset of support(Q)
    assert set(P.keys()).issubset(set(Q.keys()))
    assert P[0].parent() == RR
    assert Q[0].parent() == RR

    #for x in P:
        #print(f"Contribution 1 + %.20f for x = %d with weight %.20f" % ((P[x]/Q[x])**(a-1) - 1, x, P[x]))
    #print(f"Before taking root sum = %.50f" % sum(P[x] * ((P[x] / Q[x])**RR(a - 1)) for x in P))
    return sum(P[x] * ((P[x] / Q[x])**RR(a - 1)) for x in P)**(RR(1.0) / RR(a - 1))

def report_RD(str, P, Q, a):
    RD = RenyiDivergence(P, Q, a)
    minlog2 = floor(-log(RD - 1) / log(2))
    print(f" * RD_{{%d}}({str}) = 1 + %E < 1 + 2^-%d" % (a, RD - 1, minlog2))

# Raw table of all the tabulated probabilities.
# Note the first value is P(X = 0), while the $k$th one is P(X >= k+1 | X >=
# 1) for k >= 1.

def produce_gaussian(sigma, center, zs=100):
    sigma = RR(sigma)
    L = (sigma * zs).floor()
    lhs = (center - L).floor()
    rhs = (center + L).floor()
    norm = RR(0.0)
    for i in range(lhs, rhs + 1):
        norm += rho(i - center, sigma)
    Q = {}
    for i in range(lhs, rhs + 1):
        Q[i] = rho(i - center, sigma) / norm
    return Q

# Assumes the center of the distribution is either 0 or 1/2
def dist_to_CDT(dist, prec):
    assert eq(dist[-1], dist[1]) or eq(dist[0], dist[1])

    double_mu = 0 if eq(dist[-1], dist[1]) else 1
    CDT = []

    L = max([ -min(dist.keys()), max(dist.keys()) ])
    scale = RR(1 << prec)

    cum_table = []
    if double_mu == 0:
        # Entry 0: P( |X| >= 1 )
        # Entry 1: P( |X| >= 2 )
        # ...
        cdt = [ (sum(dist[y] for y in dist if y >= x or y <= -x) * scale).floor() for x in range(1, L) ]
    else:
        # Entry 0: P( X >= 2 or X <= -1 )
        # Entry 1: P( X >= 3 or X <= -2 )
        # ...
        cdt = [ (sum(dist[y] for y in dist if y >= x or y <= 1-x) * scale).floor() for x in range(2, L) ]
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
    scale = 1 << prec
    res = {}
    for i in range(len(cdt) + 1):
        prevValue = scale if i == 0 else cdt[i - 1]
        nextValue = 0 if i == len(cdt) else cdt[i]
        res[i + double_mu] = RR(prevValue - nextValue) / RR(scale << 1)
        if double_mu == 0 and i == 0:
            res[0] *= RR(2.0)
        else:
            res[-i] = res[i + double_mu]
    return res

def print_keygen_table(logn, table):
    """
    Print the lookup cumulative probability table used in the keygen algorithm
    of hawk. The support is Z and the standard deviation is 1.500.
    The precision used will be 78 bits (one bit for sign, so 80 total) We store
    the high 15 bits in a separate table and the low 63 bits as well somewhere.
    Then, to compare a randomly sampled bitstring of 80 bits, we do two
    comparisons
    """
    high = []
    low = []

    for x in table:
        high.append(x >> 63)
        low.append(x & ((1 << 63) - 1))

    # Compress the high table
    while len(high) > 0 and high[-1] == 0:
        high.pop()

    print(f"static const uint64_t gauss_keygen_%d[%d] = {{\n\t%su\n}};"
        % (1 << logn, len(low), "u, ".join([str(x) for x in low ])))

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
    highs = [[], []]
    lows = [[], []]

    L = len(tables[0])
    while len(tables[1]) < L:
        tables[1].append(0)
    assert L == len(tables[1])

    def output_correctly(H, L):
        """
        Do a sanity check that printing the big numbers into two parts really
        corresponds to the number, when reconstructing it.
        """
        Hs = f"%04x" % H
        Ls = f"%016x" % L
        return (int(Hs, 16) << 63) + int(Ls, 16)

    for mu in range(2):
        for i in range(len(tables[mu])):
            hi_bits = tables[mu][i] >> 63
            lo_bits = tables[mu][i] & ((1 << 63) - 1)
            highs[mu].append(hi_bits)
            lows[mu].append(lo_bits)
            assert output_correctly(hi_bits, lo_bits) == tables[mu][i]

    # Compress the high table
    Lhi = L
    while highs[0][Lhi - 1] == 0 and highs[1][Lhi - 1] == 0:
        Lhi -= 1

    print(f"static const uint16_t gauss_hi_%d[%d] = {{" % (1 << logn, 2*Lhi))
    for i in range(Lhi):
        print(f"\t0x%04X, 0x%04X," % (highs[0][i], highs[1][i]))
    print("};")
    print(f"static const uint64_t gauss_lo_%d[%d] = {{" % (1 << logn, 2*L))
    for i in range(L):
        print(f"\t0x%016X, 0x%016X," % (lows[0][i], lows[1][i]))
    print("};")

# Since R_a(P||Q) is a non-decreasing function as a --> oo, we can take a = 513
# in any case even though a = 257 is sufficient for NIST-1.
a = 513

# Report on sigma_{sign}
# Hawk sampling during signing
for logn in [9, 10]: # NIST-1 and NIST-5 respectively.
    needed_prec = 64 + 2 + (1 + logn)
    # required Renyi divergence for signing should be < 1 + 2^needed_prec.

    for prec in range(62, 101, 8):
        valid = True
        for double_mu in [0, 1]:
            Q = produce_gaussian(sigma_sig[logn], 0.5 * double_mu)
            CDT = dist_to_CDT(Q, prec)
            P = CDT_to_dist(CDT, double_mu, prec)
            valid = valid and RenyiDivergence(P, Q, a) < 1 + RR(2)**(-needed_prec)

        if valid:
            tables = []
            print(f"\n/*\n * Precomputed CDT with %d bits of precision." % prec)
            for double_mu in [0, 1]:
                Q = produce_gaussian(sigma_sig[logn], 0.5 * double_mu)
                CDT = dist_to_CDT(Q, prec)
                P = CDT_to_dist(CDT, double_mu, prec)
                report_RD(f"sign_sampler, mu={double_mu}/2", P, Q, a)
                tables.append(CDT)
            print(" */")
            print_sign_tables(logn, tables)
            break
print()

prec = 63

# Report on sigma_{pk}
# Hawk sampling during keygen
for logn in [9, 10]: # NIST-1 and NIST-5 respectively.
    Q = produce_gaussian(sigma_pk[logn], 0)
    CDT = dist_to_CDT(Q, prec)
    P = CDT_to_dist(CDT, 0, prec)

    print("\n/*")
    report_RD("keygen", P, Q, a)
    print(" */")
    print_keygen_table(logn, CDT)
print()

# Report on sigma_{sec}
# Lower-bound on ||(f, g)||^2 during keygen
for logn in [9, 10]: # NIST-1 and NIST-5 respectively.
    n = 1 << logn

    print("static const int32_t l2bound_ssec_%d[%d] = {\n\t0u /* unused */" %
        (n, logn + 1), end="");
    for i in range(1, logn + 1):
        print(", %d" % floor(sigma_sec[logn]**2 * (2 << i)), end="")
    print("\n};");
print()

# Report on sigma_{ver}
# Upper-bound on squared distance from 2s to a target h.
for logn in [9, 10]: # NIST-1 and NIST-5 respectively.
    n = 1 << logn

    print("const uint32_t Zf(l2bound_%d)[%d] = {\n\t0u /* unused */" %
        (n, logn + 1), end="");
    for i in range(1, logn + 1):
        print(", %du" % floor((2 * sigma_ver[logn])**2 * (2 << i)), end="")
    print("\n};");
print()

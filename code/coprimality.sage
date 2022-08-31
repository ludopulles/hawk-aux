from sage.all import next_prime


def is_pow_of_2(n):
    # checks whether integer n is a power of 2
    return (n & (n-1) == 0) and n != 0


def order(x, m):
    # determines the order of x in (Z/mZ)*
    i = 1
    y = x % m
    while y != 1:
        y = (y * x) % m
        i += 1
    return i


def coprimality_probabilities(n, nb_primes=10000, hawk=False):
    assert is_pow_of_2(n)
    m = 2*n

    if not hawk:
        # the ideals (f) and (g) are not comaximal and their norms are not
        # coprime if the ramified rational prime 2 (r = n, i = 1, s = 1) lies
        # below them both

        # this happens (assuming uniformity) with probability 1/4, hence we
        # start by discounting these cases, and with probability 3/4

        pr_coprime = 3./4
        pr_norm_coprime = 3./4
    else:
        # in Hawk KeyGen we ensure that the ramified rational prime 2 does not
        # lie beneath either ideal (f) or (g), which happens with probability
        # 1/4, hence by discounting these cases we start with this probability
        # NOTE: it is not that these instances are _not completable_ but rather
        # that they will not generate a secret key in KGen

        pr_coprime = 1./4
        pr_norm_coprime = 1./4

    p = 2

    for j in range(nb_primes):
        p = next_prime(p)
        i = order(p, m)
        s = n/i

        pr_coprime *= (1-1/p**(2*i))**s

        x = (1. - 1./p**i)**s
        pr_norm_coprime *= 1 - (1 - x)**2

    return pr_coprime, pr_norm_coprime


for i in range(11):
    n = 2**i

    # print log2(n), m, our estimate for the proportion of all f, g that are
    # completable, and our estimate for the proportion of those that are
    # completable by TowerSolve
    # NOTE: the ratio ``towersolve`` is the same regardless of bool `` hawk``
    # NOTE: coprime_ideal_hawk is exactly 1/3 * coprime_ideal (non hawk)

    # coprime_ideal, coprime_norm = coprimality_probabilities(n)
    # towersolve = (coprime_norm / coprime_ideal)

    # print("%d,%d,%.5f,%.5f" % (i, 2*n, coprime_ideal, towersolve))

    coprime_ideal, coprime_norm = coprimality_probabilities(n, hawk=True)
    towersolve_hawk = (coprime_norm / coprime_ideal)

    print("%d,%d,%.5f,%.5f" % (i, 2*n, coprime_ideal, towersolve_hawk))

    # print("dim %d,  \t pr_coprime; %.5f,\t pr_norm_coprime; %.5f" %
    #       (n, coprime_ideal, coprime_norm))

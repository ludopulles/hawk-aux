# This file is an excerpt from the Leaky-LWE estimator, with minor tweaks

from math import exp
from sage.all import e, floor, log, pi, RR, copy, vector, RealDistribution, \
        ceil, sqrt


def GH_sv_factor_squared(k):
    return ((pi * k)**(1. / k) * k / (2. * pi * e))


def compute_delta(k):
    """Computes delta from the block size k. Interpolation from the following
    data table:
    Source : https://bitbucket.org/malb/lwe-estimator/
    src/9302d4204b4f4f8ceec521231c4ca62027596337/estima
    tor.py?at=master&fileviewer=file-view-default
    :k: integer
    estimator.py table:
    """

    small = {0: 1e20, 1: 1e20, 2: 1.021900, 3: 1.020807, 4: 1.019713, 5: 1.018620,
             6: 1.018128, 7: 1.017636, 8: 1.017144, 9: 1.016652, 10: 1.016160,
             11: 1.015898, 12: 1.015636, 13: 1.015374, 14: 1.015112, 15: 1.014850,
             16: 1.014720, 17: 1.014590, 18: 1.014460, 19: 1.014330, 20: 1.014200,
             21: 1.014044, 22: 1.013888, 23: 1.013732, 24: 1.013576, 25: 1.013420,
             26: 1.013383, 27: 1.013347, 28: 1.013310, 29: 1.013253, 30: 1.013197,
             31: 1.013140, 32: 1.013084, 33: 1.013027, 34: 1.012970, 35: 1.012914,
             36: 1.012857, 37: 1.012801, 38: 1.012744, 39: 1.012687, 40: 1.012631,
             41: 1.012574, 42: 1.012518, 43: 1.012461, 44: 1.012404, 45: 1.012348,
             46: 1.012291, 47: 1.012235, 48: 1.012178, 49: 1.012121, 50: 1.012065}

    if k != round(k):
        x = k - floor(k)
        d1 = compute_delta(floor(k))
        d2 = compute_delta(floor(k) + 1)
        return x * d2 + (1 - x) * d1

    k = int(k)
    if k < 50:
        return small[k]
    else:
        delta = GH_sv_factor_squared(k)**(1. / (2. * k - 2.))
        return delta.n()


def bkzgsa_gso_len(logvol, i, d, beta=None, delta=None):
    if delta is None:
        delta = compute_delta(beta)

    return RR(delta**(d - 1 - 2 * i) * exp(logvol / d))


rk = [0.789527997160000, 0.780003183804613, 0.750872218594458, 0.706520454592593, 0.696345241018901, 0.660533841808400, 0.626274718790505, 0.581480717333169, 0.553171463433503, 0.520811087419712, 0.487994338534253, 0.459541470573431, 0.414638319529319, 0.392811729940846, 0.339090376264829, 0.306561491936042, 0.276041187709516, 0.236698863270441, 0.196186341673080, 0.161214212092249, 0.110895134828114, 0.0678261623920553, 0.0272807162335610, -
      0.0234609979600137, -0.0320527224746912, -0.0940331032784437, -0.129109087817554, -0.176965384290173, -0.209405754915959, -0.265867993276493, -0.299031324494802, -0.349338597048432, -0.380428160303508, -0.427399405474537, -0.474944677694975, -0.530140672818150, -0.561625221138784, -0.612008793872032, -0.669011014635905, -0.713766731570930, -0.754041787011810, -0.808609696192079, -0.859933249032210, -0.884479963601658, -0.886666930030433]
simBKZ_c = [None] + [rk[-i] - sum(rk[-i:]) / i for i in range(1, 46)]

pruning_proba = .5
simBKZ_c += [RR(log(GH_sv_factor_squared(d)) / 2. -
                log(pruning_proba) / d) / log(2.) for d in range(46, 1500)]


def simBKZ(l, beta, tours=1, c=simBKZ_c):

    n = len(l)
    l2 = copy(vector(RR, l))

    for k in range(n - 1):
        d = min(beta, n - k)
        f = k + d
        logV = sum(l2[k:f])
        lma = logV / d + c[d]

        if lma >= l2[k]:
            continue

        diff = l2[k] - lma
        l2[k] -= diff
        for a in range(k + 1, f):
            l2[a] += diff / (f - k - 1)

    return l2


chisquared_table = {i: None for i in range(1000)}


for i in range(2049):
    chisquared_table[i] = RealDistribution('chisquared', i)


def conditional_chi_squared(d1, d2, lt, l2):
    """
    Probability that a gaussian sample (var=1) of dim d1+d2 has length at most
    lt knowing that the d2 first cordinates have length at most l2
    """
    D1 = chisquared_table[d1].cum_distribution_function
    D2 = chisquared_table[d2].cum_distribution_function
    l2 = RR(l2)

    PE2 = D2(l2)
    # In large dim, we can get underflow leading to NaN
    # When this happens, assume lifting is successfully (underestimating security)
    if PE2 == 0:
        raise ValueError("Numerical underflow in conditional_chi_squared")

    steps = 5 * (d1 + d2)

    # Numerical computation of the integral
    proba = 0.
    for i in range(steps)[::-1]:
        l2_min = i * l2 / steps
        l2_mid = (i + .5) * l2 / steps
        l2_max = (i + 1) * l2 / steps

        PC2 = (D2(l2_max) - D2(l2_min)) / PE2
        PE1 = D1(lt - l2_mid)

        proba += PC2 * PE1

    return proba


def predict_sign_forge_beta(d, logvol, sver, tours=1):
    delta = compute_delta(2)
    l = [log(bkzgsa_gso_len(logvol, i, d, delta=delta)) / log(2)
         for i in range(d)]

    if tours == 1:
        betas = range(2, d)
    else:
        betas = [floor(x/tours) for x in range(1+floor(2*tours), ceil(d*tours))]

    for beta in betas:
        l = simBKZ(l, beta, 1)
        if 2**l[0] <= sver * sqrt(d):
            return beta
    raise ValueError('not solvable')


def predict_beta_and_prev_sd(d, logvol, tours=1, ignore_lift_proba=False,
                             lift_union_bound=False, number_targets=1,
                             verbose=False):
    """
    Computes the beta value for given dimension and volumes, as well as the
    standard deviation of the first basis vector just before detection of
    a short vector.

    It is assumed that the instance has been normalized and sphericized,
    i.e. that the covariance matrices of the secret is the identity
    :d: integer
    :vol: float
    """

    remaining_proba = 1.
    average_beta = 0.
    average_sd = 0.
    cumulated_proba = 0.

    delta = compute_delta(2)
    l = [log(bkzgsa_gso_len(logvol, i, d, delta=delta)) / log(2.)
         for i in range(d)]

    if tours == 1:
        betas = range(2, d)
    else:
        betas = [floor(x/tours) for x in range(1+floor(2*tours), ceil(d*tours))]

    for beta in betas:
        current_sd = RR((2**l[0])/sqrt(d))
        l = simBKZ(l, beta, 1)
        proba = 1.
        i = d - beta
        proba *= chisquared_table[beta].cum_distribution_function(
            2**(2 * l[i]))

        if not ignore_lift_proba:
            for j in range(2, int(d / beta + 1)):
                i = d - j * (beta - 1) - 1
                xt = 2**(2 * l[i])
                if j > 1:
                    if not lift_union_bound:
                        x2 = 2**(2 * l[i + (beta - 1)])
                        d2 = d - i + (beta - 1)
                        proba *= conditional_chi_squared(beta - 1, d2, xt, x2)
                    else:
                        proba = min(proba, chisquared_table[d-i].cum_distribution_function(xt))

        for t in range(number_targets):
            average_beta += beta * remaining_proba * proba
            average_sd += current_sd * remaining_proba * proba
            cumulated_proba += remaining_proba * proba
            remaining_proba *= 1. - proba

        if verbose:
            print("β= %d,\t pr=%.4e, \t cum-pr=%.4e \t rem-pr=%.4e" % (beta, proba, cumulated_proba, remaining_proba))

        if remaining_proba < .001:
            average_beta += beta * remaining_proba
            break

    if remaining_proba > .01:
        raise ValueError("This instance may be unsolvable")

    return average_beta, average_sd

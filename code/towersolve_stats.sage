# exploring properties of unimodular transform of (f, g)
# do $ sage --preparse hawk.sage
#    $ mv hawk.sage.py hawk.py
# to allow imports from hawk.sage and similarly _attack_Zn.sage

from _attack_Zn import general_circulant
from hawk import CyclotomicTower, NotCoprimeError

from sage.all import matrix, ZZ, identity_matrix, block_matrix, is_even, \
    randint
from sage.arith.misc import xgcd
from sage.stats.distributions.discrete_gaussian_integer import \
        DiscreteGaussianDistributionIntegerSampler
from sage.matrix.matrix_integer_dense_hnf import hnf_with_transformation


class NotCompletableError(Exception):
    pass


class binomialD:
    def __init__(dist, sigma):
        # for binomial distributions a la Kyber etc, sigma must be an integer
        # note it is not standard deviation, instead variance is sigma / 2
        assert sigma.is_integer()
        dist.sigma = sigma

    def __call__(dist):
        bits = [randint(0, 1) for _ in range(2*dist.sigma)]
        return sum([bits[i]-bits[i+1] for i in range(0, 2*dist.sigma, 2)])


def tower_solve_fraction(CT, sigma=1.5, repeats=1000, hawk=False,
                         binomial=False):
    """
    Determines how many (f, g) samples from a discrete Gaussian over the
    integers with standard deviation sigma are required to have ``repeats``
    many such completable, and also how many of the ``repeats`` many are
    completable via TowerSolve

    :param CT:      a CyclotomicTower from hawk.sage
    :param sigma:   a standard deviation controlling samples (f, g)
    :param repeats: the number of completable f, g
    :param hawk:    a bool, if True requires both f and g to have odd
                        algebraic norm and to be long enough wrt sigma_sec
    :returns:       a dictionary with the total number of required samples
                        to achieve ``repeats`` many completable samples, the
                        number of these completable by TowerSolve, and for
                        those that fail, the gcd (N(f), N(g))
    """

    total_samples = 0
    completable_samples = 0
    tower_solvable_samples = 0
    tower_failing_gcds = {}
    for _ in range(repeats):
        print(_)
        f, g, attempts = sample_fg(CT, sigma=sigma, check_coprime=False,
                                   check_completable=True, count=True,
                                   hawk=hawk, binomial=binomial)
        total_samples += attempts
        completable_samples += 1
        d, _, _ = xgcd(f.norm(), g.norm())
        if abs(d) == 1:
            tower_solvable_samples += 1
        else:
            tower_failing_gcds[d] = tower_failing_gcds.get(d, 0) + 1
    return {"total": total_samples, "complete": completable_samples,
            "solvable": tower_solvable_samples, "t-fails": tower_failing_gcds}


def wrapper_tower_solve_fraction(powers, sigma=1.5, repeats=1000, hawk=False,
                                 binomial=False):
    """
    A wrapper for tower_solve_fraction above that takes a list of powers of two
    """
    for i in powers:
        CT = CyclotomicTower(2**i)
        output = tower_solve_fraction(CT, sigma=sigma, repeats=repeats,
                                      hawk=hawk, binomial=binomial)
        data = "{power},{complete},{solved}\n".format(
                power=i, complete=float(output["complete"]/output["total"]),
                solved=float(output["solvable"]/output["complete"]))
        if hawk:
            filename = "tower_solve_fraction_hawk_" + str(round(sigma, 2))
        else:
            filename = "tower_solve_fraction_" + str(round(sigma, 2))
        if binomial:
            filename += "_binomial"
        with open("../data/" + filename, "a") as out:
            out.write(data)
            out.close()


def sample_fg(CT, sigma=1.5, check_coprime=False, check_completable=True,
              count=True, hawk=False, binomial=False):
    """
    A function that samples f, g from a discrete Gaussian of standard
    deviation ``sigma`` and dimension determined by CT, and then checks
    various conditions on f, g

    :param CT:      a CyclotomicTower from hawk.sage
    :param sigma:   a standard deviation controlling samples (f, g)
    :param check_coprime:
                    a bool, if True (f, g) are resampled until their algebraic
                        norms are coprime
    :param check_completable:
                    a bool, if True (f, g) are resampled until they are
                        completable
    :param count:   a bool, if True also return the number of samples required
                        until the requirements on (f, g) are satisfied
    :param hawk:    a bool, if True requires both f and g to have odd
                        algebraic norm and to be long enough wrt sigma_sec
    :returns:       an f, g satisfying the given requirements, and if ``count``
                        also the number of samples required to produce such
    """
    if binomial:
        assert sigma.is_integer()
        D = binomialD(sigma)
    else:
        D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    K = CT.top_K
    deg = K.degree()

    if hawk:
        if deg == 512:
            sqr_len_lower = 1.425**2 * 2 * deg
        elif deg == 1024:
            sqr_len_lower = 1.429**2 * 2 * deg
        else:
            # non hawk parameters, so no sigma_sec
            sqr_len_lower = 0.

    attempts = 0

    while True:
        attempts += 1
        f = K([D() for _ in range(deg)])
        g = K([D() for _ in range(deg)])

        if check_coprime:
            try:
                d, _, _, = xgcd(f.norm(), g.norm())
                if d != 1:
                    raise NotCoprimeError
            except NotCoprimeError:
                continue

        if check_completable:
            try:
                complete = ground_truth(CT, f, g)
                if not complete:
                    raise NotCompletableError
            except NotCompletableError:
                continue

        if hawk:
            sqr_len = sum([f[i]**2 for i in range(deg)]) + sum([g[i]**2 for i in range(deg)]) # noqa
            fn = f.norm()
            gn = g.norm()
            if is_even(fn) or is_even(gn) or sqr_len <= sqr_len_lower:
                continue

        break

    if count:
        return f, g, attempts
    else:
        return f, g


def ground_truth(CT, f, g):
    """
    A function that returns True if and only if (f, g) is completable

    :param CT:      a CyclotomicTower from hawk.sage
    :param f:       a vector sampled from a discrete Gaussian
    :param g:       a vector sampled from a discrete Gaussian
    :returns:       a bool that is True if and only if (f, g) are completable
    """
    K = CT.top_K
    deg = K.degree()
    Mf = matrix(ZZ, general_circulant(f))
    Mg = matrix(ZZ, general_circulant(g))

    Mfg = block_matrix(2, 1, [[Mf], [Mg]])
    H, U = hnf_with_transformation(Mfg)

    if H[:deg] == identity_matrix(ZZ, deg):
        u = U[:1]
        G = K(u[0, :deg].list())
        F = -K(u[0, deg:2*deg].list())
        assert f*G - g*F == K([1])
        return True
    return False

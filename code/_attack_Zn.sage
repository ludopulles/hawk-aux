from fpylll import BKZ, GSO, LLL, IntegerMatrix
from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
from fpylll.tools.bkz_stats import dummy_tracer

try:
    from g6k.algorithms.bkz import naive_bkz_tour
    from g6k.siever import Siever
    from g6k.utils.stats import dummy_tracer as g6k_dummy_tracer
except ModuleNotFoundError:
    pass

from sage.all import load, set_random_seed, randint, matrix, ZZ, RR, \
        block_matrix, sqrt
from sage.stats.distributions.discrete_gaussian_integer import \
        DiscreteGaussianDistributionIntegerSampler
from sage.stats.distributions.discrete_gaussian_lattice import \
        DiscreteGaussianDistributionLatticeSampler
from sage.matrix.matrix_integer_dense_hnf import hnf_with_transformation

from hawk import SignatureScheme

load("load_strategies.sage")


def general_circulant(f):
    K = f.parent()
    deg = K.degree()
    z = K([0, 1])
    M = matrix(ZZ, [[(z**j * f)[i] for i in range(deg)] for j in range(deg)])
    return M


def one_check(b_one, b, red_object, g6k=False):

    if g6k:
        dim = red_object.full_n
        red_object.update_gso(0, dim)
        lengths = [sum([red_object.M.B[i][j]**2 for j in range(dim)]) for i in range(dim)]  # noqa
    else:
        dim = red_object.d
        red_object.update_gso()
        lengths = [float(red_object.G[i, i]) for i in range(dim)]

    if 1.0 in lengths and b_one is None:
        return b
    return b_one


def generate_U_DvW(n, sigma, B=None):
    # If B is not None sample according to D_{Q, sigma}

    if B is None:
        D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    else:
        if abs(B.determinant()) != 1:
            print("Determinant of basis is {det}".format(det=B.determinant()))
            raise NotImplementedError
        D = DiscreteGaussianDistributionLatticeSampler(ZZ**n, sigma=sigma)
        Binv = matrix(ZZ, B.inverse())

    while True:
        if B is None:
            # Sample a linearly set of vectors in Z^n
            Y = matrix(ZZ, [[D() for _ in range(n)] for _ in range(n)])
        else:
            # sample according to D_{Q, s} by sampling from ZZ^n and * B^-1
            Y = matrix(ZZ, [D() * Binv for _ in range(n)])
        if Y.determinant():
            break

    # Extract a random representative of the class of PQF representing Z^n
    T, Utinv = hnf_with_transformation(Y.transpose())
    assert T == Utinv * Y.transpose()
    Uinv = Utinv.transpose()

    U = matrix(ZZ, Uinv.inverse())
    assert abs(U.determinant()) == 1, "random detU = %d" % (U.determinant())
    return U


def solve_instance(red_object, tours, g6k=False):

    if g6k:
        dim = red_object.full_n
    else:
        dim = red_object.d

    def get_prevnorm(red_object, g6k):
        if g6k:
            red_object.update_gso(0, dim)
            min_sqr_len = min([sum([red_object.M.B[i][j]**2 for j in range(dim)]) for i in range(dim)])  # noqa
            return RR(sqrt(min_sqr_len))
        else:
            return min([RR(sqrt(red_object.G[i, i])) for i in range(dim)])

    def get_profile(red_object, g6k):
        if g6k:
            red_object.update_gso(0, dim)
            return [red_object.M.get_r(i, i) for i in range(dim)]
        else:
            return [red_object.get_r(i, i) for i in range(dim)]

    prevnorm = get_prevnorm(red_object, g6k)
    profiles = {}

    profiles[0] = get_profile(red_object, g6k)

    b_one = None

    if g6k:

        def no_d4f(b):
            return int(0 * b)

        red_object.lll(0, dim)

        b = 2
        profiles[b] = get_profile(red_object, g6k)

        b_one = one_check(b_one, b, red_object, g6k=g6k)

        b = 4

        while b_one is None:
            prevnorm = get_prevnorm(red_object, g6k)
            b += 1
            for _ in range(tours):
                naive_bkz_tour(red_object, g6k_dummy_tracer, b,
                               dim4free_fun=no_d4f)
            red_object.lll(0, dim)

            b_one = one_check(b_one, b, red_object, g6k=g6k)

            profiles[b] = get_profile(red_object, g6k)

    else:
        lll = LLL.Reduction(red_object)
        lll()

        b = 2
        profiles[b] = get_profile(red_object, g6k)

        b_one = one_check(b_one, b, red_object, g6k=g6k)

        while b_one is None:
            prevnorm = get_prevnorm(red_object, g6k)
            b += 1
            param = BKZ.Param(block_size=b, strategies=strategies) # noqa
            bkz = BKZ2(lll)
            for _ in range(tours):
                bkz.tour(param, 0, dim, dummy_tracer)
            bkz.lll_obj()

            b_one = one_check(b_one, b, red_object, g6k=g6k)

            profiles[b] = get_profile(red_object, g6k)

    assert b_one == b
    return b_one, prevnorm, profiles


def one_experiment_structured(params):
    n, seed, sigma, tours, dual, randomise, float_type, g6k, hawk = params # noqa

    if g6k and (dual or randomise):
        print("g6k with randomisation or dual currently not implemented")
        return

    set_random_seed(randint(0, 2**31) + seed)

    # sigma_sig and sigma_ver do not come into the below, so just set as equal
    sig = SignatureScheme(n, sigma, sigma, sigma)
    sk, pk = sig.KGen(hawk=hawk)

    q00 = pk[0]
    q10 = pk[1]
    q11 = pk[2]
    q01 = q10.conjugate()

    B00 = matrix(ZZ, general_circulant(sk[0]))
    B01 = matrix(ZZ, general_circulant(sk[1]))
    B10 = matrix(ZZ, general_circulant(sk[2]))
    B11 = matrix(ZZ, general_circulant(sk[3]))
    B = block_matrix(2, 2, [[B00, B01], [B10, B11]])

    assert abs(B.determinant()) == 1, "det B =" + str(B.determinant())

    Q00 = matrix(ZZ, general_circulant(q00))
    Q10 = matrix(ZZ, general_circulant(q10))
    Q01 = matrix(ZZ, general_circulant(q01))
    Q11 = matrix(ZZ, general_circulant(q11))
    Q = block_matrix(2, 2, [[Q00, Q01], [Q10, Q11]])

    assert B * B.transpose() == Q, "wrong form"

    deg = q00.parent().degree()

    if dual:
        Q = matrix(ZZ, Q.inverse())
        assert abs(Q.determinant()) == 1, "dual det Q=" + str(Q.determinant())

    if randomise:
        V = generate_U_DvW(2 * deg, sigma, B=B)
        assert abs(V.determinant()) == 1
        Q = V * Q * V.transpose()

    assert abs(Q.determinant()) == 1, 'detQ =' + str(Q.determinant())
    Q = IntegerMatrix.from_matrix(Q)

    if g6k:
        B = IntegerMatrix.from_matrix(B)
        gso = GSO.Mat(B, float_type=float_type,
                      U=IntegerMatrix.identity(2*n),
                      UinvT=IntegerMatrix.identity(2*n))
    else:
        gso = GSO.Mat(Q, gram=True, float_type=float_type)

    gso.discover_all_rows()
    gso.update_gso()

    if g6k:
        reduction_object = Siever(gso)
    else:
        reduction_object = gso

    b_one, prevnorm, profiles = solve_instance(reduction_object, tours,
                                               g6k=g6k)
    return b_one, RR(prevnorm/sqrt(2 * n)), profiles


def one_experiment(params):
    # from fpylll import FPLLL
    # FPLLL.set_threads(8)
    n, seed, sigma, tours, dual, randomise, float_type, g6k = params

    if g6k and (dual or randomise):
        print("g6k with randomisation or dual currently not implemented")
        return

    set_random_seed(randint(0, 2**31) + seed)

    U = generate_U_DvW(n, sigma)
    Q = U * U.transpose()

    if randomise:
        V = generate_U_DvW(n, sigma, B=U)
        Q = V * Q * V.transpose()

    if dual:
        Q = matrix(ZZ, Q.inverse())

    assert abs(U.determinant()) == 1, "det U =" + str(U.determinant())
    assert abs(Q.determinant()) == 1, "det Q =" + str(Q.determinant())

    Q = IntegerMatrix.from_matrix(Q)

    # if dual, pre-reduce with higher precision, cause the numbers are BIG !
    if dual:
        gso = GSO.Mat(Q, gram=True, float_type="ld")
        lll = LLL.Reduction(gso)
        lll()

    if g6k:
        Utr = IntegerMatrix.from_matrix(U.transpose())
        gso = GSO.Mat(Utr, float_type=float_type,
                      U=IntegerMatrix.identity(n),
                      UinvT=IntegerMatrix.identity(n))
    else:
        gso = GSO.Mat(Q, gram=True, float_type=float_type)

    gso.discover_all_rows()
    gso.update_gso()

    if g6k:
        reduction_object = Siever(gso)
    else:
        reduction_object = gso

    b_one, prevnorm, profiles = solve_instance(reduction_object, tours,
                                               g6k=g6k)
    return b_one, RR(prevnorm/sqrt(n)), profiles

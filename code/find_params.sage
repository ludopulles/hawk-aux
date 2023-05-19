from sage.all import ceil, e, floor, Infinity, IntegerRange, load, log, pi, \
        RR, sqrt
from proba_utils import centered_discrete_Gaussian_law, law_square, \
        law_convolution, iter_law_convolution, pos_head_probability, \
        tail_probability
load("BKZ_simulator.sage")


def slog(x, base):
    """
    A ``safe`` log for dealing with log of float(0.0)

    :param x:       input to the log
    :param base:    the base to take the logarithm
    :returns:       either log(x, base) or -inf

    .. note::       -Infinity here is really some approximation of how strictly
                        cutoff is set in clean_dist of proba_utils -- in the
                        paper we rather arbitrarily set -Infinity as -300
    """
    if x > 0.0:
        return log(x, base)
    return -Infinity


def rhf(beta):
    """
    The root Hermite factor for BKZ with blocksize at least 50

    :param beta:    an integer blocksize for BKZ
    :returns:       the root Hermite factor as a float
    """
    assert beta >= 50 and beta.is_integer(), "beta must be an int >= 50"

    return float(((beta/(2.*pi*e))*(pi*beta)**(1./beta))**(1/(2.*(beta-1))))


def log_gsa_lens(d, beta, vol=1, first=True):
    """
    The (natural) log of the Gram--Schmidt basis vector lengths after reduction
    according to the Geometric Series Assumption and the root Hermite factor

    :param d:       an integer dimension
    :param beta:    an integer blocksize <= d for BKZ
    :param vol:     the volume of the lattice being reduced
    :param first:   a bool that if set means only the first Gram--Schmidt
                        length is returned
    :returns:       a list of log lengths for the Gram--Schmidt basis vectors
    """
    assert d.is_integer() and d > 0, "the dimension must be a positive int"
    assert beta <= d, "require beta <= d"
    assert vol > 0, "the lattice volume must be positive"
    assert isinstance(first, bool), "first must be a bool"

    rhf_beta = rhf(beta)
    if first:
        return [(d - 1.)*log(rhf_beta) + (1./d) * vol]
    return [(d - 2.*i - 1.)*log(rhf_beta) + (1./d) * vol for i in range(d)]


def approx_SVP_beta(d, sver, simulate=True):
    """
    Estimating the blocksize beta required to find an approximate short
    vector of square length d * sver**2 for ZZ^d

    :param d:       an integer dimension
    :param sver:    a standard deviation controlling the length of acceptable
                        signatures (equivalently, the distance from a target)
    :returns:       the estimated beta for a (strong) signature forgery
    """
    assert sver > 0, "standard deviation sver must be positive"

    if simulate:
        # we scale the lattice by sqrt(d), as sver determines the length of the
        # lattice vector to be discovered, it must be similarly scaled
        beta = predict_sign_forge_beta(d, d*log(d)/2, sqrt(d) * sver) # noqa
        return beta

    log_sqr_ver_bound = float(2.*log(sver) + log(d))
    success_beta = None
    for beta in IntegerRange(50, d):
        log_first_sqr = 2. * log_gsa_lens(d, beta)[0]
        if success_beta is None and log_first_sqr <= log_sqr_ver_bound:
            return beta


def key_recovery_beta_ssec(d, tours=1):
    """
    Use the leaky-LWE-estimator to determine the expected successful blocksize
    that recovers a lattice vector of length one from some basis of ZZ^d, along
    with the length of the shortest vector found before the vector of length
    one i.e. len = prev_sd * sqrt(d)

    :param d:       an integer dimension
    :returns:       ``beta`` an estimate for the required blocksize for key
                        recovery, and ``prev_sd``the estimates the length of
                        the previous shortest vector
    """
    beta, prev_sd = predict_beta_and_prev_sd(d, d*log(d)/2,  # noqa
                                             lift_union_bound=True,
                                             number_targets=d, tours=tours)
    # lattice was scaled up, so prev_sd need to be scaled back
    prev_sd /= sqrt(d)
    return beta, prev_sd


def statistical_ssign(d, lam, q_s):
    """
    Determine the minimum value of the standard deviation controlling the
    Gaussian sampler within Sign, via [Sec. 2.6, FALCON] considering ZZ^d

    :param d:       an integer dimension
    :param lam:     a security parameter
    :param q_s:     the number of signature queries allowed to an adversary
    :returns:       a standard deviation
    """
    assert d.is_integer() and d > 0, "the dimension must be a positive int"
    assert lam > 0, "the security parameter must be positive"
    assert q_s.is_integer() and q_s > 0, "q_s must be a positive int"
    eps = 1/sqrt(q_s * lam)
    ssign = (1/pi) * sqrt(log((2*d)*(1 + 1/eps)) / 2)
    return RR(ssign)


def fail_and_forge_probabilities(d, ssign, sver, ssec):
    """
    Determines the probabilities that a signature is too long (fail) and
    that the (weak) signature forgery attack works

    :param d:           an integer dimension
    :param ssign:       a standard deviation controlling the length of elements
                            sampling during signing
    :param sver:        a standard deviation controlling the length of
                            acceptable signatures
    :param ssec:        a standard deviation that describes the length of the
                            shortest vector discovered before key recovery
    :returns:           log2 of the ``fail`` probability and log2 of the (weak)
                            forgery probability
    """

    def diagonal_B(d, len_sqr):
        """
        Describe the entries of a vector with d entries, of square length
        close to, but bounded above by, len_sqr, and such that the sizes of
        the entries are as similar as possible

        :param d:       an integer dimension
        :param len_sqr: the square length that upper bounds the vector
        :returns:       a dictionary entries such that entries[x] = y means
                            the vector should have y entries with value x

        """
        assert d > 0 and d.is_integer(), "dimension should be a postive int"
        assert len_sqr > 0 and len_sqr.is_integer(), "len_sqr not +ve int"
        entries = {}
        entry_upper_bound = ceil((len_sqr / d)**.5)
        if entry_upper_bound == 1:
            entries[1.] = len_sqr
            entries[0.] = d - entries[1.]
        else:
            assert (len_sqr >= d)
            assert (len_sqr <= d*entry_upper_bound**2)

            upp = entry_upper_bound
            # solve: d_(upp-1)+d_upp=d
            # and  : (upp-1)^2 d_(upp-1) + upp^2 d_upp <= fixed_sqr_len
            # (up to rounding)
            dupp = floor((len_sqr - d*(upp-1)**2)/(2*upp-1))
            dupp1 = d - dupp

            assert (dupp+dupp1 == d)
            assert (dupp1*(upp-1)**2+dupp*upp**2 <= len_sqr)
            assert (dupp1*(upp-1)**2+dupp*upp**2 > len_sqr-upp**2)

            entries[upp-1] = dupp1
            entries[upp] = dupp
        return entries

    sqr_radius = d * sver**2.
    flen = floor(d**.5 * ssec)
    flen_sqr = floor(d * ssec**2)

    A = centered_discrete_Gaussian_law(ssign)

    Bonedim = {flen: 1}

    A2 = law_square(A)
    A2_d = iter_law_convolution(A2, d)

    fail_probability = tail_probability(A2_d, sqr_radius)
    fail = float(slog(fail_probability, 2))

    # calculate forge for B all in one coordinate
    A2_d1 = iter_law_convolution(A2, d - 1)
    Conedim = law_convolution(A, Bonedim)
    C2onedim = law_square(Conedim)
    Donedim = law_convolution(A2_d1, C2onedim)

    onedimforge_probability = pos_head_probability(Donedim, sqr_radius)
    forgeonedim = float(slog(onedimforge_probability, 2))

    # calculate forge for B as ``diagonal`` as possible
    try:
        entries = diagonal_B(d, flen_sqr)

        distrs = []
        for (entry, conv) in entries.items():
            partial_diag = law_convolution(A, {entry: 1})
            partial_diag2 = law_square(partial_diag)
            partial_diag2c = iter_law_convolution(partial_diag2, conv)
            distrs += [partial_diag2c]

        C2diag = law_convolution(distrs[0], distrs[1])

        forgediag_probability = pos_head_probability(C2diag, sqr_radius)
        forgediag = float(slog(forgediag_probability, 2))

    except AssertionError:
        print('diagonal B failed to make a vector')

    forge = max(forgediag, forgeonedim)

    return fail, forge


def falcon_blocksizes(d, falcon=False, k_fac=True):
    """
    A function to allow comparisons to the Falcon style of estimate for the
    required blocksizes for key recovery and strong signature forgery

    :param d:       an integer dimension
    :param falcon:  a bool that if ``True`` considers Falcon parameters, and
                        otherwise Hawk parameters
    :param k_fac:   a multiplicative factor of sqrt(4/3) that is applied to
                        rhs during the key recovery attack
    :returns:       an estimate ``beta_key`` for key recovery and
                        ``beta_forge`` for strong signature forgery

    .. note::       in the Falcon key recovery methodology they apply both
                        k_fac and dimensions for free techniques
    """
    assert d == 1024 or d == 2048, "please choose a Falcon/Hawk d"
    assert isinstance(falcon, bool), "falcon must be a bool"
    assert isinstance(k_fac, bool), "k_fac must be a bool"

    n = int(d/2)

    if falcon:
        q = 12289
    else:
        q = 1
        if n == 512:
            sver = 1.425
        if n == 1024:
            sver = 1.572
    dth_root_vol = sqrt(q)

    def gh_sqr(d):
        """
        The Gaussian heuristic for the square length of a shortest vector in
        a unit volume lattice of dimension d

        :param d:       a dimension
        :returns:       the "exact" heuristic (i.e. using the Gamma function
                            without Stirling's approximation)
        """
        assert d >= 50, "the dimension must be >= 50"
        return float((d/(2*pi*e))*(pi*d)**(1/d))

    # key recovery
    for beta_key in range(50, d):
        if falcon:
            skeygen = 1.17*(q/d)**.5
            lhs = sqrt(beta_key) * skeygen
        else:
            # projections of something of length 1
            lhs = sqrt(beta_key / d)
        rhs = (beta_key/(2*pi*e))**(1 - n/beta_key) * dth_root_vol
        # more accurate is
        # power = (2*beta_key - d + 1)/(2*(beta_key-1))
        # rhs = gh_sqr(beta_key)**power * dth_root_vol
        if k_fac:
            rhs *= (4/3)**.5
        if lhs <= rhs:
            break

    # signature forgery
    if falcon:
        sver = falcon_sver(n)
        verif_length_sqr = floor((1.1 * sver * sqrt(d))**2)
    else:
        verif_length_sqr = d * sver**2
    for beta_forge in range(50, d):
        lhs = (beta_forge/(2*pi*e))**(2*n/beta_forge) * q
        # more accurate is
        # lhs = gh_sqr(beta_forge)**((d-1)/(beta_forge-1)) * q
        rhs = verif_length_sqr
        if lhs <= rhs:
            break
    print("d, falcon, sqrt{4/3}:", d, falcon, k_fac)
    print("beta_key, beta_forge:", beta_key, beta_forge)
    return beta_key, beta_forge


def dimensions_for_free(beta):
    """
    Apply the dimensions for free reduction to blocksize using the optimistic
    value from Sec 4.3 ePrint 2017/999

    :param beta:    an input blocksize
    :returns:       a smaller blocksize via dimensions for free techniques
    """
    return beta - round(beta * log(4./3.) / log(beta / (2 * pi * e)))


def falcon_sver(n):
    """
    Return the standard deviations that control the length of acceptable
    signatures for Falcon-512 and Falcon-1024

    :param n:       element of {512, 1024}
    :returns:       standard deviation sver for Falcon

    .. note::       We would calculate this quantity in general as
                    eps = 1/sqrt(lam*q_s)
                    first = (1/pi) * sqrt((log(4*n*(1 + 1/eps)))/2)
                    second = 1.17 * sqrt(q)
                    return RR(first*second)
                    for lam = 128, q = 12289 and q_s = 2**64
    """
    assert n == 512 or n == 1024, "Choose a Falcon n"

    if n == 512:
        return 165.736617183
    if n == 1024:
        return 168.388571447


def sign_fail_sver_search(d, ssign, q_s=2**64):
    """
    A function to final minimal sver such that signature failure probability
    is  1 / (q_s ** 2)

    :param d:       an integer dimension
    :param ssign:   a standard deviation controlling the length of elements
                        sampling during signing
    :param q_s:     allowable signature queries for statistical arguments
    :returns:       an sver

    ..note::        we only use this for HAWK-1024 parameters
    """
    lower = ssign
    upper = ssign * 2**.5
    steps = 1000
    delta = (upper - lower) / steps
    target = -2 * log(q_s, 2)

    A = centered_discrete_Gaussian_law(ssign)
    A2 = law_square(A)
    A2_d = iter_law_convolution(A2, d)

    break_next_iteration = False

    for i in range(steps):
        sver = lower + i*delta

        sqr_radius = d * sver**2.
        fail_probability = tail_probability(A2_d, sqr_radius)
        fail = float(slog(fail_probability, 2))

        if break_next_iteration:
            break

        if fail <= target:
            # gives a little headroom
            break_next_iteration = True

    return sver, fail


def find_params(d, lam, q_s=2**64, beta_key=None, ssec=None):
    """
    A script to find parameters that satisfy our signature forgery
    cryptanalysis, and statistical requirements.

    :param d: an integer dimension = 2*n
    :param lam: security parameter for statistical arguments
    :param q_s: allowable signature queries for statistical arguments
    :param beta_key: beta required to recover secret key, if ``None``
            is calculated using BKZ_simulator
    :param ssec: sigma_sec from paper, if ``None`` calculated using
            BKZ_simulator

    :returns: dictionary of parameters and failure rates

    ..note: one must check the beta for signature forging and key recovery
            are sufficient, lam and q_s input parameters are solely for
            statistical arguments
    """
    ssign = statistical_ssign(d, lam, q_s)
    if beta_key is None or ssec is None:
        beta_key, ssec = key_recovery_beta_ssec(d)
        beta_key = int(beta_key)
        ssec = float(ssec)

    assert ssec > ssign
    assert ssign > ssec / 2
    # fix sver to ssign (its lower bound)
    sver = ssign
    _, forge = fail_and_forge_probabilities(d, ssign, sver, ssec)

    if forge > -lam:
        print('no parameters')
        return

    # now increase sver until either forge becomes too high
    # or
    # both fail and forge are small enough (no point increasing sver further)
    # or
    # beta_forge less than beta_key (no point having forgery harder than key
    # recovery)

    steps = 15
    delta = min((2**.5 - 1.05) * ssign, (ssec - ssign)) / (steps - 1)

    assert steps > 2
    assert delta > 0

    for i in range(steps):
        print("i of steps:", i, steps)

        sver = ssign + i*delta

        print("calculating fail and forge for", i)
        fail, forge = fail_and_forge_probabilities(d, ssign, sver, ssec)
        # weak forgeries too easy
        if forge > -lam:
            sver = ssign + (i - 1) * delta
            break

        print("calculating beta forge for", i)
        beta_forge = approx_SVP_beta(d, sver=sver)

        # if both failure probability and forgery probability low enough, stop
        if fail <= -lam and forge <= -lam:
            break

        print("sver, fail, forge, beta_forge", sver, fail, forge, beta_forge)

        # at this point forgery easier or equal to key recovery, so stop
        if beta_forge <= beta_key:
            break

    print("recalculating everything")
    fail, forge = fail_and_forge_probabilities(d, ssign, sver, ssec)
    beta_forge = approx_SVP_beta(d, sver=sver)

    assert ssign <= sver
    assert sver <= ssec
    assert forge <= -lam

    print("Dim, lambda, log2(q_s):", d, lam, log(q_s, 2).n())
    print("key recovery beta:", beta_key)
    print("forging beta     :", beta_forge)
    print("s{sign, sec, ver}:", ssign, ssec, sver)
    print("log2(p_sig_fail) :", fail)
    print("log2(p_sig_forge):", forge)
    print()

    params = {}
    params['beta_key'] = beta_key
    params['beta_forge'] = beta_forge
    params['ssign'] = ssign
    params['ssec'] = ssec
    params['sver'] = sver
    params['fail'] = fail
    params['forge'] = forge

    return params


# find hawk parameters, if beta_key and ssec absent, start from scratch
# find_params(512, 64, q_s=2**32, beta_key=211, ssec=1.042)
# find_params(1024, 128, q_s=2**64, beta_key=452, ssec=1.425)

# find_params(512, 64, q_s=2**32)
# find_params(1024, 128, q_s=2**64)

# We find HAWK-1024 parameters slightly differently as we have much more
# room to play with sver
"""
ssign = statistical_ssign(2048, 256, 2**64)
# beta_key, ssec = key_recovery_beta_ssec(2048)
# beta_key = int(beta_key)
# ssec = float(ssec)
ssec = 1.974
beta_key = 940
sver, _ = sign_fail_sver_search(2048, ssign)
sver = round(sver, 3)
fail, forge = fail_and_forge_probabilities(2048, ssign, sver, ssec)
beta_forge = approx_SVP_beta(2048, sver=sver)
params = {}
params['beta_key'] = beta_key
params['beta_forge'] = beta_forge
params['ssign'] = ssign
params['ssec'] = ssec
params['sver'] = sver
params['fail'] = fail
params['forge'] = forge
print(params)
"""

# compare hawk and falcon using falcon methodology

# falcon_blocksizes(1024, falcon=True, k_fac=True)

# falcon_blocksizes(2048, falcon=True, k_fac=True)

# beta_key, beta_forge = falcon_blocksizes(1024, falcon=True, k_fac=False)
# print(dimensions_for_free(beta_key), dimensions_for_free(beta_forge))

# beta_key, beta_forge = falcon_blocksizes(2048, falcon=True, k_fac=False)
# print(dimensions_for_free(beta_key), dimensions_for_free(beta_forge))

# falcon_blocksizes(1024, falcon=False, k_fac=True)

# falcon_blocksizes(2048, falcon=False, k_fac=True)

# beta_key, beta_forge = falcon_blocksizes(1024, falcon=False, k_fac=False)
# print(dimensions_for_free(beta_key), dimensions_for_free(beta_forge))

# beta_key, beta_forge = falcon_blocksizes(2048, falcon=False, k_fac=False)
# print(dimensions_for_free(beta_key), dimensions_for_free(beta_forge))

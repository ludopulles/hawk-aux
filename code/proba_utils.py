from math import ceil, floor, exp


def variance(D, centred=True):
    var = 0.
    for val in D:
        var += val * val * D[val]
    if not centred:
        raise NotImplementedError
    return var


def clean_dist(A, cutoff=750):
    """
    Clean a distribution to accelerate further computation by dropping elements
    with probability less than or equal to 2^-cutoff

    :param A:       input law (dictionary)
    :param cutoff:  2^-cutoff is used to determine whether an element is kept
                        if None then the distribution is simply returned
    :returns:       a dictionary with no values <= 2^-cutoff
    """
    if cutoff is None:
        return A
    B = {}
    for (x, y) in A.items():
        if y > 2**(-cutoff):
            B[x] = y
    # normalizing will increase the probabilities of the remaining support
    return normalize(B)


def normalize(A):
    """
    Make sure the values in a dictionary representing a probability law has
    values that sum to one

    :param A:       input law (dictionary)
    :returns:       a dictionary with entries scaled such that sum(values) = 1
    """
    B = {}
    s = 0.
    for (x, y) in A.items():
        s += y
    for (x, y) in A.items():
        B[x] = y/s
    return B


def rho(sigma, x):
    return exp(-x*x / (2*sigma*sigma))


def centered_discrete_Gaussian_law(sigma, c=0., tau=12.):
    B = {}
    lb = floor(c - sigma * tau)
    ub = ceil(c + sigma * tau)
    for x in range(lb, ub+1):
        B[x] = rho(sigma, x - c)
    return clean_dist(normalize(B))


def law_convolution(A, B):
    """
    Construct the convolution of two laws (sum of independent variables from
    two input laws)

    :param A:       first input law (dictionary)
    :param B:       second input law (dictionary)
    :returns:       their convolution law as a dictionary
    """

    C = {}
    for a in A:
        for b in B:
            c = a+b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def law_product(A, B):
    """
    Construct the law of the product of independent variables from two input
    laws

    :param A:       first input law (dictionary)
    :param B:       second input law (dictionary)
    :returns:       their product law as a dictionary
    """
    C = {}
    for a in A:
        for b in B:
            c = a*b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def law_square(A):
    """
    Construct the law given by squaring another law

    :param A:       first input law (dictionary)
    :param B:       second input law (dictionary)
    :returns:       the square law as a dictionary
    """
    C = {}
    for a in A:
        C[a**2] = C.get(a**2, 0) + A[a]
    return C


def iter_law_convolution(A, i):
    """
    Compute the i fold convolution of a distribution (using double and add)

    :param A:   input law (dictionary)
    :param i:   integer determining how many times A is convolved with itself
    :returns:   i fold convolution as a dictionary
    """
    # use fast exponentiation for iterated convolution
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = law_convolution(D, D)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution(D, A)
            D = clean_dist(D)
    return D


def tail_probability(D, t, strict=True):
    """
    Probability that a sample drawn from D is (strictly) greater than t in
    absolute value

    :param D:       input law (dictionary)
    :param t:       tail parameter (non negative integer)
    :param strict:  a bool, if True consider the probability only for elements
                        in the support of D strictly greater than t
    :returns:       the probability a sample from D is (strictly) larger than
                        t in absolute value
    """
    s = 0.
    ma = max([abs(key) for key in D.keys()])

    if strict:
        if t >= ma:
            return 0.
    else:
        if t > ma:
            return 0.

    if strict:
        lb = floor(t) + 1
    else:
        lb = ceil(t)
    ub = ma + 1

    # Summing in reverse for better numerical precision
    # (assuming tails are decreasing)
    for i in reversed(range(lb, ub)):
        if i == 0:
            s += D.get(i, 0)
        else:
            s += (D.get(i, 0) + D.get(-i, 0))
    return s


def pos_head_probability(D, t, strict=True):
    """
    Probability that a sample drawn from D is (strictly) less in absolute
    value than t

    :param D:       input law (dictionary)
    :param t:       tail parameter (non negative integer)
    :param strict:  a bool, if True consider the probability only for elements
                        in the support of D strictly greater than t
    :returns:       the probability a sample from D is (strictly) smaller than
                        t in absolute value
    """
    s = 0.

    if t == 0 and strict:
        return 0.

    s += D.get(0, 0)

    if strict:
        ub = ceil(t)
    else:
        ub = floor(t) + 1

    for i in range(1, ub):
        s += (D.get(i, 0) + D.get(-i, 0))

    return s

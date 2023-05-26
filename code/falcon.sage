from math import e, floor, log, pi, sqrt


def dimensions_for_free(beta):
    """
    Apply the dimensions for free reduction to blocksize using the optimistic
    value from Sec 4.3 ePrint 2017/999

    :param beta:    an input blocksize
    :returns:       a smaller blocksize via dimensions for free techniques
    """
    return beta - round(beta * log(4./3.) / log(beta / (2 * pi * e)))


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
                    k_fac and dimensions for free techniques, as decscribed in
                    App D
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


# compare hawk and falcon using falcon methodology

falcon_blocksizes(1024, falcon=True, k_fac=True)

falcon_blocksizes(2048, falcon=True, k_fac=True)

beta_key, beta_forge = falcon_blocksizes(1024, falcon=True, k_fac=False)
print(dimensions_for_free(beta_key), dimensions_for_free(beta_forge))

beta_key, beta_forge = falcon_blocksizes(2048, falcon=True, k_fac=False)
print(dimensions_for_free(beta_key), dimensions_for_free(beta_forge))

falcon_blocksizes(1024, falcon=False, k_fac=True)

falcon_blocksizes(2048, falcon=False, k_fac=True)

beta_key, beta_forge = falcon_blocksizes(1024, falcon=False, k_fac=False)
print(dimensions_for_free(beta_key), dimensions_for_free(beta_forge))

beta_key, beta_forge = falcon_blocksizes(2048, falcon=False, k_fac=False)
print(dimensions_for_free(beta_key), dimensions_for_free(beta_forge))

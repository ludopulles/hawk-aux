from sage.all import CyclotomicField, identity_matrix, is_even, log, \
        randint, xgcd
from sage.stats.distributions.discrete_gaussian_integer import \
        DiscreteGaussianDistributionIntegerSampler as DGaussZ
from sage.stats.distributions.discrete_gaussian_lattice import \
        DiscreteGaussianDistributionLatticeSampler as DGauss


class NotCoprimeError(Exception):
    pass


class CyclotomicTower:
    def __init__(CT, top_n):
        """
        Creates an object giving access to the cyclotomic fields required.
        In particular, Q, Q(zeta_2), ..., Q(zeta_{2*top_n}), for top_n a power
        of two.

        In this case zeta_{2*top_n} is a root of the irreducible polynomial
        X^{top_n} + 1 = 0 and Z[zeta_{2*top_n}] is the ring of integers of
        Q(zeta_{2*top_n}).
        """
        CT.top_ln = log(top_n)/log(2)
        if CT.top_ln != round(CT.top_ln):
            raise ValueError("Cyclotomic degree must be a power of 2")

        CT.top_K = CyclotomicField(2*top_n)
        CT.Ks = [CyclotomicField(2**(i+1)) for i in range(CT.top_ln)] + [CT.top_K]  # noqa

    def level(CT, f):
        """
        Returns quantities related to the cyclotomic field f lives in

        :param f: an element of a power of two cyclotomic

        :returns: the degree of the parent field of f, this degree log2, and
            the parent field itself
        """
        n = f.parent().degree()
        ln = log(n)/log(2)
        if ln != round(ln):
            raise ValueError("Cyclotomic degree must be a power of 2")
        return n, ln, CT.Ks[ln]

    def bar(CT, f):
        """
        The complex conjugation of f = alpha_i zeta^i, i.e. the effect of the
        map zeta -> bar(zeta) = zeta^-1, hence if f has degree n
        bar(f) = alpha_0 - alpha_(n-1) zeta - ... = alpha_1 zeta^(n-1)

        :param f: element of a cyclotomic field

        :returns: the conjugate element of f
        """
        return f.conjugate()

    def tau(CT, f):
        """
        The involutive automorphism defined by zeta -> -zeta, hence if f has
        degree n, tau(f) =
            alpha_0 - alpha_1 zeta + ... + (-1)^(n-1) alpha_(n-1) zeta^(n-1)

        :param f: element of a cyclotomic field

        :returns: tau(f)

        ..note::  in Sec 4.1 https://eprint.iacr.org/2019/015.pdf when p = 2
                  we have tau(f) = f^x, i.e. the sole non trivial element of
                  Gal(LL/KK) in their notation
        """
        n, _, K = CT.level(f)
        return K([f[i] * (-1)**i for i in range(n)])

    def norm_down(CT, f):
        """
        Move one step down the tower, as per TowerSolve methods from
        https://eprint.iacr.org/2019/015.pdf, this is the map N() defined
        in Sec 4.1 from LL to KK

        :param f: a cyclotomic field element

        :returns: its norm in the cyclotomic field with degree halved
        """
        n, ln, _ = CT.level(f)
        ff = f*CT.tau(f)
        return CT.Ks[ln-1](list(ff)[::2])

    def NTRU_solve(CT, f, g):
        """
        Find a solution (F, G) in our cyclotomic field of degree top_d such
        that fG - gF = 1 in the field, following TowerSolve methods

        :param f: an element of the ring of integers of a cyclotomic field
        :param g: an element of the ring of integers of a cyclotomic field

        :returns: a pair (F, G) reduced with respect to (f, g) such that
                  fG - gF = 1
        """
        _, ln, K = CT.level(f)
        if g.parent() != K:
            return ValueError("Both input must belong to the same field")

        if ln == 0:
            d, u, v = xgcd(f, g)
            if (d > 1):
                # recall, if q = 1 the condition is equivalent to coprimality
                raise NotCoprimeError
            return -v, u

        f_, g_ = CT.norm_down(f), CT.norm_down(g)
        F_, G_ = CT.NTRU_solve(f_, g_)
        F = K(F_) * CT.tau(g)
        G = K(G_) * CT.tau(f)
        return CT.reduce(f, g, F, G)

    def find_k(CT, f, g, F, G):
        """
        Calculating the non rounded factor in line 2 of Algorithm 1 from
        https://eprint.iacr.org/2019/015.pdf
        """
        numerator = F * CT.bar(f) + G * CT.bar(g)
        denominator = f * CT.bar(f) + g * CT.bar(g)
        return numerator / denominator

    def partition_round(CT, x):
        """
        We define a rounding function that partitions RR into translated
        regions of the form [n - 1/2, n + 1/2) for n in ZZ, that is, a
        rounding function that always rounds up, even for negative n

        :param x: a real number

        :returns: the closest integer, rounding up for ties

        ..note::  there is varied rounding behaviour between SageMath
                  and Python2 and Python3, be careful!
        """
        if x < 0 and float(2 * x) % 1 == 0 and float(x) % 1 != 0:
            return 1 - round(-x)
        else:
            return round(x)

    def partition_round_poly(CT, f):
        """
        Perform the custom rounding defined in CT.partition_round to a
        polynomial coefficientwise

        :param f: an element of a cyclotomic field

        :returns: a ``close'' element of the ring of integers
        """
        return f.parent([CT.partition_round(c) for c in f])

    def reduce(CT, f, g, F, G):
        """
        Calculating reduced (F, G) from line 3 of Algorithm 1
        https://eprint.iacr.org/2019/015.pdf

        ..note::  we use the more optimal ffNP in HAWK's real KGen, but
                  broadly both reduction routines are equivalent -- they give
                  a canonical (F, G) in some fundamental domain determined by
                  (f, g)
        """
        k = CT.find_k(f, g, F, G)
        rounded_k = CT.partition_round_poly(k)
        return F - rounded_k * f, G - rounded_k * g


class SignatureScheme:
    def __init__(self, n, sigma_kg, sigma_sig, sigma_ver):
        """
        Initialise the signature scheme over power of two degree n cyclotomics
        with the various sigma. In particular initialise the discrete Gaussian
        samplers for KGen and for Sign, create the tower of cyclotomic fields
        for KGen, and set the verification length bound for Vf

        :param n:           the degree of the cyclotomic field
        :param sigma_kg:    the Gaussian parameter determining key generation
        :param sigma_sig:   the Gaussian parameter used in signing
        :param sigma_ver:   the Gaussian parameter determining valid signatures

        ..note::      to keep the code simpler we do not half h as in HAWK
                      proper, so we work in cosets 2 ZZ^(2n) + c with
                      c in ZZ^(2n).

        ..note::      we use row notation in this implementation, as opposed
                      to the column notation of the paper
        """
        self.n = n
        # sampler for KGen
        self.D0_kg = DGaussZ(sigma=sigma_kg)
        # two coset samplers for Sign, using 2*sigma_sig as over 2 ZZ^(2n)
        self.D0_sig = DGauss(2*identity_matrix(1), sigma=2*sigma_sig)
        self.D1_sig = DGauss(2*identity_matrix(1), sigma=2*sigma_sig, c=(1,))
        # initialise CyclotomicTower for KGen
        self.CT = CyclotomicTower(n)
        # verification bound over 2 ZZ^(2n) using 2*sigma_ver
        self.verif_bound = 2 * n * (2 * sigma_ver) ** 2

    def Dsig(self, x):
        """
        If some element of the target t is even sample from 2 ZZ, and if odd
        sample from 2 ZZ - 1 using DGS with sigma_sig

        :param x: an integer

        :returns: a sample from a discrete Gaussian depending on the parity
                  of x
        """
        if x % 2:
            return self.D1_sig()[0] - 1
        else:
            return self.D0_sig()[0]

    def KGen(self, hawk=True):
        """
        Perform KGen from HAWK, with some restart conditions missing

        :param hawk:    if ``True`` also check that the algebraic norms of
                        f and g are both odd

        :returns:       a secret key and a public Hermitian form
        """
        while True:
            try:
                f = self.CT.top_K([self.D0_kg() for i in range(self.n)])
                g = self.CT.top_K([self.D0_kg() for i in range(self.n)])
                # mimic HAWK KGen restart
                if hawk:
                    if is_even(f.norm()) or is_even(g.norm()):
                        continue
                F, G = self.CT.NTRU_solve(f, g)
                break
            except NotCoprimeError:
                continue

        # Basis B = [f g]
        #           [F G]
        sk = (f, g, F, G)

        bar = self.CT.bar

        # Q = B B^* (B^* is adjoint transpose)
        q00 = f*bar(f) + g*bar(g)
        q10 = F*bar(f) + G*bar(g)
        q11 = F*bar(F) + G*bar(G)
        # we do not compute q01 = q10*
        pk = (q00, q10, q11)

        return sk, pk

    def sym_break(self, h1, s1):
        """
        Implementing line 8 of Sign
        """
        e = h1 - 2*s1

        return e != 0 and next(filter(lambda x: x != 0, e)) > 0

    def Sign(self, sk, h):
        """
        Perform HAWK signing, we assume h = H(m || r) is given rather than a
        message and uniformly sampled salt, hence a signature is just
        s = (s0, s1)

        :param sk:  a HAWK secret key
        :param h:   a pair of lists with n binary entries each

        :returns:   a signature (s0, s1)

        ..note::  we do not check the length of x so failures are possible, but
                  should become less frequent when sigma_ver > sigma_sig
                  as n increases
        """
        f, g, F, G = sk
        K = self.CT.top_K
        h0 = K(h[0])
        h1 = K(h[1])

        # t = (h0, h1) * B
        t0 = f * h0 + F * h1
        t1 = g * h0 + G * h1

        x0 = K([self.Dsig(r) for r in t0])
        x1 = K([self.Dsig(r) for r in t1])

        # e = x * B^{-1}
        e0 = x0*G - x1*F
        e1 = -x0*g + x1*f

        # s = (h - e)/2 is a lattice point
        s0 = (h0 - e0) / 2
        s1 = (h1 - e1) / 2

        if not self.sym_break(h1, s1):
            s0 = h0 - s0
            s1 = h1 - s1

        return s0, s1

    def Vf(self, pk, h, sig):
        """
        Perform HAWK verification checking the length of e (from Sign above)
        with respect to public form Q, and also the sym_break condition

        :param pk:  a HAWK public key
        :param h:   a pair of lists with n binary entries each
        :param sig: a HAWK signature

        :returns:   a bool which is True iff Vf succeeds
        """
        (s0, s1) = sig
        (h0, h1) = h
        (q00, q10, q11) = pk
        bar = self.CT.bar
        K = self.CT.top_K

        integer = s0.is_integral() and s1.is_integral()
        sym = self.sym_break(self.CT.top_K(h1), s1)

        e0 = K(h0) - 2 * s0
        e1 = K(h1) - 2 * s1

        v = e0*q00*bar(e0) + e1*q11*bar(e1) + e1*q10*bar(e0) + e0*bar(q10*e1)
        length = v[0] <= self.verif_bound

        return integer and sym and length

    def CompressPK(self, pk):
        """
        Compress the public key by dropping q11 and half of the coordinates
        coordinates of q00

        :param pk:  a HAWK public key

        :returns:   a compressed HAWK public key
        """
        (q00, q10, q11) = pk
        q00_ = [q00[i] for i in range(self.n / 2)]
        return (q00_, q10)

    def DecompressPK(self, pkc):
        """
        Decompress a public key by recomputing q11 and recovering the dropped
        entries of q00

        :param pkc: a compressed HAWK public key

        :returns:   a HAWK public key
        """
        (q00_, q10) = pkc

        # e.g. n = 8, we have q00_ = alpha_0 + ... + alpha_3 z^3 and want
        # q00 = q00_ + 0 z^4 - alpha_3 z^5 - alpha_2 z^6 - alpha_1 z^7

        q00 = q00_ + [0] + [-q00_[i] for i in range(self.n / 2 - 1, 0, -1)]
        q00 = self.CT.top_K(q00)
        return (q00, q10, (1 + q10*self.CT.bar(q10))/q00)

    def CompressSK(self, sk):
        """
        Drop G from a HAWK secret key

        :param sk:  a HAWK secret key

        :returns:   a compressed HAWK secret key
        """
        (f, g, F, G) = sk
        return (f, g, F)

    def DecompressSK(self, skc):
        """
        Decompress a HAWK secret key by recomputing G

        :param skc: a compressed HAWK secret key

        :returns:   a decompressed HAWK secret key
        """
        (f, g, F) = skc
        return (f, g, F, (1 + g*F) / f)

    def CompressSig(self, s):
        """
        Compress a signature by dropping the first component

        :param s: a HAWK signature s = (s0, s1)

        :returns: a compressed HAWK signature s1
        """
        s0, s1 = s
        return s1

    def DecompressSig(self, pk, h, s1):
        """
        Decompress a HAWK signature following Sec 3.2 of HAWK

        :param pk:  an uncompressed HAWK public key
        :param h:   a pair of lists with n binary entries each
        :param s1:  a compressed HAWK signature

        :returns:   a decompressed HAWK signature
        """
        K = self.CT.top_K
        (q00, q10, _) = pk
        (h0, h1) = h
        s0 = self.CT.partition_round_poly(K(h0) / 2 + (K(h1) / 2 - s1) * q10 / q00) # noqa
        return s0, s1


def test_run(n=32, sigma_kg=1.5, sigma_sig=1, sigma_ver=1.1):
    Sig = SignatureScheme(n, sigma_kg, sigma_sig, sigma_ver)

    sk, pk = Sig.KGen()

    sk_ = Sig.DecompressSK(Sig.CompressSK(sk))
    if sk_ != sk:
        print("secret key decompression failed")

    pk_ = Sig.DecompressPK(Sig.CompressPK(pk))
    if pk_ != pk:
        print("public key decompression failed")

    h = ([randint(0, 1) for _ in range(n)], [randint(0, 1) for _ in range(n)])
    s0, s1 = Sig.Sign(sk_, h)
    sig = (s0, s1)

    sig_ = Sig.DecompressSig(pk_, h, Sig.CompressSig(sig))
    if sig_ != sig:
        print("signature decompression failed")

    if Sig.Vf(pk_, h, sig_):
        print('verifies')
    else:
        print('does not verify')

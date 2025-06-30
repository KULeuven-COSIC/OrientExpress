import cypari2
import logging

from sage.all import (
    ZZ, factor
)

pari = cypari2.Pari()

# Usage: logging.getLogger('isolib.utils').setLevel(logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

def discrete_log_pari(a, base, order):
    """
    Wrapper around pari discrete log. Works like a.log(b),
    but allows us to use the optional argument order. This
    is important as we skip Pari attempting to factor the
    full order of Fp^k, which is slow.
    """
    x = pari.fflog(a, base, order)
    return ZZ(x)

def weil_pairing_pari(P, Q, D, check=False):
    """
    Wrapper around Pari's implementation of the Weil pairing
    Allows the check of whether P,Q are in E[D] to be optional
    """
    if check:
        nP, nQ = D*P, D*Q
        if nP.is_zero() or nQ.is_zero():
            raise ValueError("points must both be n-torsion")

    return pari.ellweilpairing(P.curve(), P, Q, D)

def BiDLP(R, P, Q, D, ePQ=None):
    """
    Given a basis P,Q of E[D] finds
    a,b such that R = [a]P + [b]Q.

    Uses the fact that
        e([a]P + [b]Q, [c]P + [d]Q) = e(P,Q)^(ad-bc)

    Optional: include the pairing e(P,Q) which can be precomputed
    which is helpful when running multiple BiDLP problems with P,Q
    as input. This happens, for example, during compression.
    """
    # e(P,Q)
    if ePQ:
        pair_PQ = ePQ
    else:
        pair_PQ = weil_pairing_pari(P, Q, D)

    # Write R = aP + bQ for unknown a,b
    # e(R, Q) = e(P, Q)^a
    pair_a = weil_pairing_pari(R, Q, D)

    # e(R,-P) = e(P, Q)^b
    pair_b = weil_pairing_pari(R, -P, D)

    # Now solve the dlog in Fq
    a = discrete_log_pari(pair_a, pair_PQ, D)
    b = discrete_log_pari(pair_b, pair_PQ, D)

    return a, b

class FactIso:

    def __init__(self, K0, d=None, fac=None):
        """
        Given a curve and a kernel point, compute the factored isogeny with
        kernel that point. Sage has it builtin but sometimes it just doens't
        work. d is the degree, fac its factorization.
        """
        if not d:
            d = K0.order()
        if not fac:
            fac = list(factor(d))
        self.d = d
        self.fac = fac

        self.chain = []
        E = K0.curve()
        K = K0
        for l, e in fac:
            for ei in range(e):
                logger.debug(f'Doing {l = } ({ei+1}/{e})')
                d_i = l**(ei + 1)
                cof = d//d_i
                K_i = K * cof
                phi_i = E.isogeny(K_i)

                self.chain.append(phi_i)
                E = phi_i.codomain()
                K = phi_i(K)

    def __call__(self, P):
        for phi_i in self.chain:
            P = phi_i(P)
        return P

    def __str__(self):
        return f'Factored isogeny of degree {factor(self.d)} from {self.domain()} to {self.codomain()}'

    def __repr__(self):
        return self.__str__()

    def dual(self):
        return DualFactIso(self.chain, self.d, self.fac)

    def domain(self):
        return self.chain[0].domain()

    def codomain(self):
        return self.chain[-1].codomain()

class DualFactIso(FactIso):

    def __init__(self, chain, d, fac=None):
        """
        Wrapper for the dual of a factored isogeny. Takes in input the degree d
        and the chain of isogenies, and returns the dual isogeny as FactIso.
        """
        self.chain = []
        for phi in chain[::-1]:
            self.chain.append(phi.dual())
        if not fac:
            fac = list(factor(d))
        self.d = d
        self.fac = fac

def compute_log_order(P, ell, hi=None):
    """
    P is a point of order ell^n, recover n. hi is the optional upper bound on
    n.
    """
    lo = 0
    if not hi:
        # Upper bound not provided
        hi = 1
        while (ell**hi)*P != 0:
            lo = hi
            hi *= 2

    # (ell**lo) * P != 0 / (ell**hi)*P == 0
    while hi > lo + 1:
        c = (hi + lo) // 2
        if (ell**c)*P == 0:
            hi = c
        else:
            lo = c
    return hi

def order_from_factors(P, fac):
    """
    Given a point P and some factors, compute P when guaranteed that its order
    divides the product of the given factors.
    """
    d = 1
    for ell, e in fac:
        d *= ell**e
    assert d*P == 0, "wrong factors given"
    o = 1
    for ell, e in fac:
        cof = d // ell**e
        cP = cof*P
        if cP == 0:
            continue
        ei = compute_log_order(cP, ell, e+1)
        o *= ell**ei
    P._order = ZZ(o)
    return o

def compute_log_mul_order(wp, ell, hi=None):
    """
    wp is an element of order ell^n, recover n. hi is the optional upper bound
    on n.
    """
    lo = 0
    if not hi:
        # Upper bound not provided
        hi = 1
        while wp**(ell**hi) != 1:
            lo = hi
            hi *= 2

    # (ell**lo) * P != 0 / (ell**hi)*P == 0
    while hi > lo + 1:
        c = (hi + lo) // 2
        if wp**(ell**c) == 1:
            hi = c
        else:
            lo = c
    return hi

def mul_order_from_factors(wp, fac):
    """
    Given an element wp in some multiplicative group, compute its order when
    guaranteed that it divides the product of the given factors
    """
    if type(fac) not in [tuple, list]:
        fac = factor(fac)
    d = 1
    for ell, e in fac:
        d *= ell**e
    assert wp**d == 1, "wrong factors given"
    o = 1
    for ell, e in fac:
        cof = d // ell**e
        cwp = wp**cof
        if cwp == 1:
            continue
        ei = compute_log_mul_order(cwp, ell, e+1)
        o *= ell**ei
    return ZZ(o)

def _torsion_basis_prime(E, l, e):
    """
    Find a basis of E[l^e]
    TODO: can be optimized for e=2
    """
    p = E.base_ring().characteristic()
    cof = (p+1)//(l**e)
    while True:
        P = E.random_point() * cof
        if (l**(e-1))*P != 0:
            break
    while True:
        Q = E.random_point() * cof
        if (l**(e-1))*Q == 0:
            continue
        ePQ = weil_pairing_pari(P, Q, l**e)
        if ePQ**(l**(e-1)) != 1:
            break
    return P, Q, ePQ

def _torsion_basis_odd(E, l):
    """
    Return torsion basis of E[l] for l odd
    Warning: will generally fail for l non squarefree
    """
    p = E.base_ring().characteristic()
    E._order = ZZ(p+1)
    cof = (p+1)//l

    P = E.random_point()*cof
    Q = E.random_point()*cof

    ePQ = weil_pairing_pari(P, Q, l)
    _ord = mul_order_from_factors(ePQ, l)
    if _ord == l:
        return P, Q, ePQ
    missing_fac = l // _ord
    while True:
        logger.debug(f'{missing_fac.factor() = }')
        # Fix P
        cof = (p+1)//missing_fac
        P1 = E.random_point()*cof
        Q1 = Q*_ord
        if P1 and Q1:
            eP1Q1 = weil_pairing_pari(P1, Q1, missing_fac)
            _ord1 = mul_order_from_factors(eP1Q1, missing_fac)
            if _ord1 != 1:
                P += P1
                missing_fac /= _ord1
                _ord *= _ord1
            if missing_fac == 1:
                return P, Q, None

        # Fix Q
        cof = (p+1)//missing_fac
        Q1 = E.random_point()*cof
        P1 = P*_ord
        if P1 and Q1:
            eP1Q1 = weil_pairing_pari(P1, Q1, missing_fac)
            _ord1 = mul_order_from_factors(eP1Q1, missing_fac)
            if _ord1 != 1:
                Q += Q1
                missing_fac /= _ord1
                _ord *= _ord1
            if missing_fac == 1:
                return P, Q, None

        # If we are here, most likely both P and Q are missing some factor
        cof = (p+1)//missing_fac
        R = E.random_point()*cof
        Q += R
        ePQ = weil_pairing_pari(P, Q, l)
        _ord = mul_order_from_factors(ePQ, l)
        if _ord == l:
            return P, Q, ePQ
        missing_fac = l // _ord

def torsion_basis(E, l, e=None, alg='primes'):
    """
    Find a basis of E[l] assuming E supersingular over Fp2
    Options for `alg`:
    - `primes`: compute a basis of the l^e torsion where l is a
      prime; e has to be provided;
    - `odd`: compute a basis of the l torsion, where l is odd
    TODO: deterministic / random
    """
    if alg.lower() == 'primes':
        return _torsion_basis_prime(E, l, e)
    if alg.lower() == 'odd':
        return _torsion_basis_odd(E, l)




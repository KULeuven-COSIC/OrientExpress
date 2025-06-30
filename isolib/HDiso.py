from sage.all import (
    PolynomialRing,
    EllipticCurve,
    HyperellipticCurve,
    EllipticCurveIsogeny,
    Matrix,
    vector,
    ZZ,
    factor,
    product,
    set_verbose
)

from .theta.theta_structures.couple_point import CouplePoint
from .theta.theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt

import logging
# usage: logging.getLogger('isolib.hdiso').setLevel(logging.DEBUG/INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

# Debug utils
def dump_kernel(kernel, n):
    P, Q = kernel
    logger.debug('Failed iso')
    logger.debug(f'# Prime')
    logger.debug(f'p = {P[0].curve().base_ring().characteristic()}')
    logger.debug(f'# Curve')
    logger.debug(f'E1 = EllipticCurve(Fp2, {list(P[0].curve().a_invariants())})')
    logger.debug(f'E2 = EllipticCurve(Fp2, {list(P[1].curve().a_invariants())})')
    logger.debug(f'\n# Kernel')
    logger.debug(f'P1 = E1({P[0].x()}, {P[0].y()})')
    logger.debug(f'P2 = E2({P[1].x()}, {P[1].y()})')
    logger.debug(f'Q1 = E1({Q[0].x()}, {Q[0].y()})')
    logger.debug(f'Q2 = E2({Q[1].x()}, {Q[1].y()})')
    logger.debug(f'{n = }')

# Higher dimensional isogenies
def is_split_from_product(kernel):
    """
    Starting from a product of elliptic curves checks
    if the codomain of the isogeny defined by <P, Q>
    is split.
    """
    P1, P2 = kernel[0][0], kernel[0][1]
    Q1, Q2 = kernel[1][0], kernel[1][1]

    if 0 in [P1, Q1, P1+Q1, P2, Q2, P2+Q2]:
        return True, 0

    Fp2 = P1.curve().base_ring()

    a1, a2, a3 = P1.x(), Q1.x(), (P1 + Q1).x()
    b1, b2, b3 = P2.x(), Q2.x(), (P2 + Q2).x()
    M = Matrix(Fp2, [
        [a1*b1, a1, b1],
        [a2*b2, a2, b2],
        [a3*b3, a3, b3]])
    return M.determinant().is_zero(), M

def is_split_from_jacobian(kernel):
    D1, D2 = kernel
    h = D1.parent().curve().hyperelliptic_polynomials()[0]

    G1, _ = D1
    G2, _ = D2
    G3, r3 = h.quo_rem(G1 * G2)

    M = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    return M.determinant().is_zero(), M

class ProductIsogeny:

    def __init__(self, kernel=None, maps=None):
        """
        Wrapper class to represent a product isogeny.  If two maps are
        provided, the matrix defining the isogeny is diagonal. Otherwise the
        four components will make the matrix.  Also the kernel can be provided.
        """

        if maps:
            # Can provide 2 or 4 maps forming a matrix (diagonal if 2)
            self._maps = maps
            if not len(maps) in [2,4]:
                raise ValueError(f"{len(maps) = }, should be 2 or 4")

        elif kernel:
            # Kernel is provided: three cases are possible (see Kani, Nuber of
            # curves of Genus 2 with Elliptic Differentials, Thm 3 at
            # https://mast.queensu.ca/~kani/papers/numgenl.pdf):
            # - a diagonal (2, 2)-isogeny
            # - isogeny from an automorphism s.t. alpha(P2) = P1
            # - isogeny from an isomorphism s.t. alpha(P2) = P1

            self._kernel = kernel
            P1, P2 = kernel[0][0], kernel[0][1]
            Q1, Q2 = kernel[1][0], kernel[1][1]

            ls = [P1, Q1, P1+Q1, P2, Q2, P2+Q2]
            if 0 in ls:
                # Diagonal isogeny
                self._diagonal_product()

            elif P1.curve() == P2.curve():
                # The product comes from automorphisms
                self._prod_from_aut()

            else:
                # The product comes from isomorphisms
                self._prod_from_iso()

        else:
            raise ValueError("provide either maps or kernel")

    def _diagonal_product(self):
        # The product matrix is diagonal
        P1, P2 = self._kernel[0][0], self._kernel[0][1]
        Q1, Q2 = self._kernel[1][0], self._kernel[1][1]

        E1 = P1.curve()
        E2 = P2.curve()

        phi1 = E1.isogeny([P1, Q1])
        phi2 = E2.isogeny([P2, Q2])

        if phi1.degree() != 2 or phi2.degree() != 2:
            raise ValueError("wrong degree isogeny in diagonal product")

        self._maps = [phi1, phi2]

    def _prod_from_aut(self):
        # The product comes from automorphisms on the starting curve
        P1, P2 = self._kernel[0][0], self._kernel[0][1]
        Q1, Q2 = self._kernel[1][0], self._kernel[1][1]
        Ep = P1.curve()
        found_aut = False

        for alpha in Ep.automorphisms():
            if alpha(P2) == P1 and alpha(Q2) == Q1:
                logger.debug(f'Found automorphism - alpha = {alpha.rational_maps()}')
                found_aut = True
                break

        if not found_aut:
            raise ValueError("cannot find automorphisms")

        one = Ep.identity_morphism()
        self._maps = [one, -alpha, alpha.dual(), one]

    def _prod_from_iso(self):
        # The product comes from an isomorphism between the two curves
        P1, P2 = self._kernel[0][0], self._kernel[0][1]
        Q1, Q2 = self._kernel[1][0], self._kernel[1][1]

        E1 = P1.curve()
        E2 = P2.curve()

        if len(E1.automorphisms()) != 2 or len(E2.automorphisms()) != 2:
            raise NotImplementedError("extra automorphisms in product")
            # TODO: implement post-comp of isomorphisms with automorphisms

        phi = E2.isomorphism_to(E1)
        if not (phi(P2) == P1 and phi(Q2) == Q1):
            logger.warning('isomorphism not matching the kernel')
            breakpoint()
            raise ValueError("isomorphism not matching the kenrnel")
        logger.debug(f'Found isomorphism')

        one1 = E1.identity_morphism()
        one2 = E2.identity_morphism()
        self._maps = [one1, -phi, phi.dual(), one2]

    def __call__(self, P):
        P1, P2 = P

        if len(self._maps) == 2:
            Q1 = self._maps[0](P1)
            Q2 = self._maps[1](P2)

        else:
            # Phi is an actual 2x2 matrix
            Q1 = self._maps[0](P1) + self._maps[1](P2)
            Q2 = self._maps[2](P1) + self._maps[3](P2)

        return (Q1, Q2)

    def domain(self):
        return (self._maps[0].domain(), self._maps[1].domain())

    def codomain(self):
        return self._maps[0].codomain(), self._maps[-1].codomain()

    def degree(self):
        if len(self._maps) == 2:
            return self._maps[0].degree()

        return self._maps[0].degree() + self._maps[1].degree()


class HDChain:

    def __init__(self, kernel, e, alg='auto_theta'):
        """
        Input:
        - `kernel`: two 2^e-torsion points on the domain; in the case of
            a product this is (P1, P2), (Q1, Q2) with Pi, Qi in Ei
        - `e`: the length of the chain
        - `alg`: can handle different strategies:
            - 'auto_theta': detect product at the beginning, then do one theta
              isogeny
            - 'auto_richelot': use Richelot not implemented
        """

        self._chain = []
        self._degree = e
        self._kernel = kernel

        if alg.lower() == 'auto_theta':
            self.auto_theta(kernel, e)
        else:
            raise NameError(f"unknown algorithm {alg}")

    def _do_product_step(self, kernel, full_n):
        """
        Do one single product step. Update the chain and return the new kernel.
        full_n is the full number of missing steps.
        """
        logger.info('Doing product step')
        cof2 = 2**(full_n - 1)
        P1, P2 = kernel[0]
        Q1, Q2 = kernel[1]
        P = (cof2*P1, cof2*P2)
        Q = (cof2*Q1, cof2*Q2)
        K = (P, Q)

        phi = ProductIsogeny(kernel=K)
        self._chain.append(phi)
        return (phi(kernel[0]), phi(kernel[1]))

    def auto_theta(self, kernel, e):
        P = (kernel[0][0], kernel[0][1])
        Q = (kernel[1][0], kernel[1][1])

        wp_check = (P[0].weil_pairing(Q[0], 2**e))*(P[1].weil_pairing(Q[1], 2**e))
        if wp_check != 1:
            logger.debug(f'incorrect pairing {wp_check = }')
            raise ValueError(f'incorrect pairing for the kernel: {wp_check}')

        # Sanity
        m = 2**e
        if m*P[0] != 0 or m*Q[0] != 0:
            raise ValueError('wrong torsion order on E1')
        if m*P[1] != 0 or m*Q[1] != 0:
            raise ValueError('wrong torsion order on E2')

        # Detect products
        m = 2**(e-1)
        K = ((m*P[0], m*P[1]), (m*Q[0], m*Q[1]))
        is_split, _ = is_split_from_product(K)
        while is_split:
            P, Q = self._do_product_step((P, Q), e)
            e -= 1
            if e == 0:
                return
            m = 2**(e-1)
            K = ((m*P[0], m*P[1]), (m*Q[0], m*Q[1]))
            is_split, _ = is_split_from_product(K)

        if e == 0:
            return

        assert P[0] * 2**e == Q[0] * 2**e == 0, "Wrong order after product?"
        P = CouplePoint(P[0], P[1])
        Q = CouplePoint(Q[0], Q[1])
        kernel = (P, Q)
        logger.debug(f'Doing ThetaProd')
        try:
            Phi = EllipticProductIsogenySqrt(kernel, e)
            self._chain.append(Phi)
        except Exception as exc:
            dump_kernel(kernel, e)
            raise exc

    def auto_richelot(self, kernel, e):
        raise NotImplementedError("richelot not available")
        print('Warning: slow')
        P = (kernel[0][0], kernel[0][1])
        Q = (kernel[1][0], kernel[1][1])

        wp_check = (P[0].weil_pairing(Q[0], 2**e))*(P[1].weil_pairing(Q[1], 2**e))
        if wp_check != 1:
            logs.error(f'Pairing check fail on the kernel: {wp_check = }')
            breakpoint()
            raise ValueError('incorrect pairing for the kernel')

        # Sanity
        m = 2**e
        assert m*P[0] == m*Q[0] == 0, "Wrong torsion order on E1"
        assert m*P[1] == m*Q[1] == 0, "Wrong torsion order on E2"

        for i in range(e):
            m = 2**(e-i-1)

            if type(P) == tuple:
                # Starting from a product
                K = ((m*P[0], m*P[1]), (m*Q[0], m*Q[1]))
                is_split, M = is_split_from_product(K)

                if is_split:
                    # Prod2Prod (TODO: generic formula?)
                    logs.debug(f'Step {i} - Prod2Prod')
                    phi = ProductIsogeny(kernel=K)
                else:
                    # Gluing
                    logs.debug(f'Step {i} - Prod2jac')
                    phi = GluingIsogeny(K, M)

            else:
                # Starting from a jacobian
                K = (m*P, m*Q)
                is_split, M = is_split_from_jacobian(K)

                if is_split:
                    logs.debug(f'Step {i} - Jac2Prod')
                    phi = SplittingIsogeny(K, M)
                else:
                    logs.debug(f'Step {i} - Jac2Jac')
                    phi = RichelotIsogeny(K) # Here M is computed differently

            self._chain.append(phi)
            if i < e-1:
                P, Q = phi(P), phi(Q)

    def __call__(self, P):
        """
        Here P is a point on the starting curve, of the right type.
        Warning: if theta isogenies are included in the product the sign may be
        wrong. In that case, use self.eval_two_points instead.
        """
        for phi in self._chain:
            P = CouplePoint(P[0], P[1])
            P = phi(P)
        return P

    def eval_two_points(self, K1, K2):
        """
        Given two (2-dim) points K1 and K2, evaluate K1, K2, and K1 - K2 via
        phi to compute the correct sign.
        """
        P1, P2 = K1
        Q1, Q2 = K2
        phi_K1 = self(K1)
        phi_K2 = self(K2)
        phi_dt = self(((P1 - Q1), (P2 - Q2)))

        if phi_dt[0] == phi_K1[0] - phi_K2[0]:
            outP1 = phi_K1[0]
            outQ1 = phi_K2[0]
            case1 = 1

        elif phi_dt[0] == phi_K1[0] + phi_K2[0]:
            outP1 = phi_K1[0]
            outQ1 = - phi_K2[0]
            case1 = 2

        elif phi_dt[0] == - phi_K1[0] - phi_K2[0]:
            outP1 = - phi_K1[0]
            outQ1 = phi_K2[0]
            case1 = 3

        elif phi_dt[0] == - phi_K1[0] + phi_K2[0]:
            outP1 = - phi_K1[0]
            outQ1 = - phi_K2[0]
            case1 = 4

        else:
            logger.warning(f'First component wrong')
            breakpoint()

        if phi_dt[1] == phi_K1[1] - phi_K2[1]:
            outP2 = phi_K1[1]
            outQ2 = phi_K2[1]
            case2 = 1

        elif phi_dt[1] == phi_K1[1] + phi_K2[1]:
            outP2 = phi_K1[1]
            outQ2 = - phi_K2[1]
            case2 = 2

        elif phi_dt[1] == - phi_K1[1] - phi_K2[1]:
            outP2 = - phi_K1[1]
            outQ2 = phi_K2[1]
            case2 = 3

        elif phi_dt[1] == - phi_K1[1] + phi_K2[1]:
            outP2 = - phi_K1[1]
            outQ2 = - phi_K2[1]
            case2 = 4

        else:
            logger.warning(f'Second component wrong')
            breakpoint()

        logger.info(f'{case1 = } {case2 = }')
        return ((outP1, outP2), (outQ1, outQ2))

    def domain(self):
        return self._chain[0].domain()

    def codomain(self):
        return self._chain[-1].codomain()

    def degree(self):
        """
        Return the log of the degree
        """
        return self._degree

    def __repr__(self):
        return f'Isogeny chain of degree {factor(self.degree())} from {self.domain()} to {self.codomain()}'

    def __str__(self):
        return self.__repr__()


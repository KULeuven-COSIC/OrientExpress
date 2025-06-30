import time

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

from sage.all import (
        factor, valuation, product
)

from .utils import (
        order_from_factors, weil_pairing_pari, BiDLP, mul_order_from_factors,
        torsion_basis
)

class FactoredHDIsogenyEvaluator:

    def __init__(self, phi, u, v, eval_ord=None):
        """
        Given a 2D isogeny phi, representing the factorization of a 1D isogeny
        psi of degree u*(2^a - u), use it to evaluate psi. Exploiting the fact
        that the degrees u and 2^a-u are coprime, the computation can be done
        for all primes dividing p+1 at once and the result stored.
        - eval_ord is the order of the torsion that we want to evaluate phi on;
          if not provided, it is assumed to be p+1
        """
        self.phi = phi

        self.u = u
        self.v = v

        self.E1 = phi.domain()[0]
        self.E2 = phi.domain()[1]

        self.p = self.E1.base_ring().characteristic()
        if not eval_ord:
            # Select square free part of self.p+1
            logger.warning('eval order not provided')
            fac = factor(self.p+1)
            even_order = 1
            odd_order = 1
            for l, e in fac:
                if l == 2:
                    even_order *= l**e
                else:
                    odd_order *= l
            self.eval_ord = even_order*odd_order
        else:
            self.eval_ord = eval_ord

        self.E1._ord = self.p+1
        self.E2._ord = self.p+1

        self.B1 = None
        self.B2 = None

        self.E1._order = self.p+1
        self.E2._order = self.p+1

    def domain(self):
        return self.E1

    def codomain(self):
        return self.E2

    def _first_evaluation(self):
        """
        The first time that a point needs to be evaluated
        """
        logger.info('Starting first evaluation')
        t1 = time.time()
        Fp2 = self.E1.base_ring()

        eval_two = valuation(self.eval_ord, 2)
        eval_odd = self.eval_ord // (2**eval_two)
        # Remember to cache this
        known_two_torsion = valuation(self.phi._kernel[0][0].order(), 2)

        # Full torsion on E1
        if known_two_torsion == eval_two:
            # We already know a basis of the 2 torsion
            P1 = self.phi._kernel[0][0]
            Q1 = self.phi._kernel[1][0]
        else:
            P1, Q1, _e1_even = torsion_basis(self.E1, 2, eval_two)


        Px, Qx, ePQ = torsion_basis(self.E1, eval_odd, alg='odd')
        P1 = P1 + Px
        Q1 = Q1 + Qx

        P1.set_order(self.eval_ord, check=False)
        Q1.set_order(self.eval_ord, check=False)

        e1 = weil_pairing_pari(P1, Q1, self.eval_ord)
        # DEBUG
        # _ord = mul_order_from_factors(e1, self.p+1)
        # if _ord != (self.eval_ord):
        #     logger.error('Error in torsion basis gen 1')
        #     logger.error(f'{_ord = } {self.p+1 = }')
        #     breakpoint()

        t2 = time.time()
        logger.info(f'- First torsion basis done: {t2-t1:.3f}s')

        # Full torsion on E2
        if known_two_torsion == eval_two:
            # We already know a basis of the 2 torsion
            P2 = self.phi._kernel[0][1]
            Q2 = self.phi._kernel[1][1]
        else:
            P2, Q2, _e2_even = torsion_basis(self.E2, 2, eval_two)

        Px, Qx, ePQ = torsion_basis(self.E2, eval_odd, alg='odd')
        P2 = P2 + Px
        Q2 = Q2 + Qx

        P2.set_order(self.eval_ord, check=False)
        Q2.set_order(self.eval_ord, check=False)

        e2 = weil_pairing_pari(P2, Q2, self.eval_ord)
        # DEBUG
        # _ord = mul_order_from_factors(e2, self.p+1)
        # if _ord != (self.eval_ord):
        #     logger.error('Error in torsion basis gen 2')
        #     logger.error(f'{_ord = } {self.p+1 = }')
        #     breakpoint()

        t3 = time.time()
        logger.info(f'- Second torsion basis done: {t3-t2:.3f}s')

        # TODO: split factors of p+1 in coprime to u and to v; not needed for
        # now, and not more complicated, just more annoying

        # Evaluate all the points
        Pu, Qu = self.phi.eval_two_points((P1, self.E2(0)), (Q1, self.E2(0)))
        t4 = time.time()
        logger.info(f'- First eval done: {t4-t3:.3f}s')
        eu = weil_pairing_pari(Pu[0], Qu[0], self.eval_ord)
        if e1**self.u == eu:
            Pu = Pu[0]
            Qu = Qu[0]
            u_part = 0
        else:
            eu = weil_pairing_pari(Pu[1], Qu[1], self.eval_ord)
            assert e1**self.u == eu, "no u pairing"
            Pu = Pu[1]
            Qu = Qu[1]
            u_part = 1
        t5 = time.time()
        logger.info(f'- First pairing done: {t5-t4:.3f}s')

        # Here we are evaluating the other side of the square, so the
        # u-component above is the v-component here
        Pv, Qv = self.phi.eval_two_points((self.E1(0), P2), (self.E1(0), Q2))
        t6 = time.time()
        logger.info(f'- Second eval done: {t6-t5:.3f}s')
        Pv = Pv[u_part]
        Qv = Qv[u_part]
        ev = weil_pairing_pari(Pv, Qv, self.eval_ord)
        assert e2**self.v == ev
        t7 = time.time()
        logger.info(f'- Second pairing done: {t7-t6:.3f}s')

        # Recover the final image of P1 and Q1 on E2
        a, b = BiDLP(Pu, Pv, Qv, self.eval_ord, ev)
        c, d = BiDLP(Qu, Pv, Qv, self.eval_ord, ev)

        vP2 = self.v*P2
        vQ2 = self.v*Q2

        outP = a*vP2 + b*vQ2
        outQ = c*vP2 + d*vQ2

        self.B1 = (P1, Q1, e1)
        self.B2 = (outP, outQ) # = (phi_uv(P1), phi_uv(P2))
        t8 = time.time()
        logger.info(f'- BiDLP done: {t8-t7:.3f}s')
        logger.info(f'First evaluation done: {t8-t1:.3f}s')

    def __call__(self, X):
        """
        Given a point in E1, evaluate the u*v isogeny.
        """
        if not self.B1:
            self._first_evaluation()

        if not hasattr(X, '_order'):
            X.set_order(multiple=self.p+1)
        ell = X.order()

        cof = self.eval_ord // ell
        P1, Q1, e1 = self.B1
        P2, Q2 = self.B2

        P1 *= cof
        Q1 *= cof
        P2 *= cof
        Q2 *= cof

        # Easy cases to avoid pairings
        if X == P1:
            return P2
        if X == Q1:
            return Q2

        # TODO: we should be able to get this from e1
        _e1 = weil_pairing_pari(P1, Q1, ell)
        a, b = BiDLP(X, P1, Q1, ell, _e1)
        return a*P2 + b*Q2


from sage.all import (
        ZZ,
        inverse_mod,
        product,
        factor
)
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

from isolib.HDiso import HDChain
from isolib.utils import FactIso, torsion_basis
from isolib.hdtricks import FactoredHDIsogenyEvaluator

import logging
import time
# usage: logging.getLogger('cs_orientation').setLevel(logging.DEBUG/INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

class CS_StartingOrientation:
    """
    Orientation starting from E0: an endomorphism, no HD computation required.
    Should be in the precomputations. Do the same as standard Orientation but easier.
    Maybe chid class?
    """
    def __init__(self, params, E0):
        self.params = params
        self.E = E0

        Fp2 = E0.base_ring()
        z2 = Fp2.gen()
        self.p = Fp2.characteristic()
        self.E._order = ZZ(self.p+1)

        omega = WeierstrassIsomorphism(E0, (z2, 0, 0, 0))

        self.tau = E0.scalar_multiplication(ZZ(params.alpha)) + ZZ(params.beta)*omega
        self.pi = E0.frobenius_isogeny()

        # Sanity check: the chain splits
        # print(f'Splitting tau')
        # P1, Q1 = self.E.torsion_basis(2**(self.params.a + 1))
        # K1 = (self.params.l1 * P1, self.tau(P1))
        # K2 = (self.params.l1 * Q1, self.tau(Q1))
        # ker = (K1, K2)
        # Phi = HDChain(ker, self.params.a+1)

    def __call__(self, X):
        """
        Evaluate the orientation on a given point X
        """
        return self.pi(self.tau(X))

    def find_eigenspace(self, sk_sign):
        """
        Given a list sk_sign of signes of the secret key, evaluates as many
        ideals as possible at the same time. In the case of the starting
        orientation, since evaluating the orientation is cheap, all ideals are
        evaluated. Return a point in the kernel, and a list denoting the ideals
        evaluated.
        """
        logger.info('Finding eigenspace (E0 orientation)')
        t1 = time.time()
        missing_primes = self.params.primes[:]
        X = self.E(0)
        _ord = 1

        while missing_primes != []:
            # TODO: can do more things together
            logger.debug(f'{missing_primes = }')
            new_missing_primes = []
            m = product(missing_primes)
            cof = (self.params.p+1) // m
            P = self.E.random_point() * cof

            # Evaluate the point
            Q = self(P)

            # Remove the wrong eigenspaces
            for ell in missing_primes:
                logger.debug(f'Find eigenspace for {ell = }')
                ell_id = self.params.primes.index(ell)
                if sk_sign[ell_id] == 0:
                    continue

                cof_ell = m // ell
                sig_ell = int(sk_sign[ell_id] < 0)
                mu_ell = self.params.eigenvals[ell][1-sig_ell]

                P_ell = P*cof_ell
                if P_ell == 0:
                    new_missing_primes.append(ell)
                    continue
                Q_ell = cof_ell*Q - mu_ell * P_ell
                if Q_ell == 0:
                    new_missing_primes.append(ell)
                    continue
                X += Q_ell
                _ord *= ell

            missing_primes = new_missing_primes

        X.set_order(multiple=_ord)

        t2 = time.time()
        logger.info(f'Eigenspace found: {t2-t1:.3f}s')
        return X, sk_sign

    def push_forward(self, X):
        """
        Given a point X on E0, pushforward the orientation through the
        isogeny with kernel X.
        """
        logger.info('Starting pushforward (E0 orientation)')
        t1 = time.time()
        m = X.order()

        phi = X.curve().isogeny(X)
        E1 = phi.codomain()
        E1._order = ZZ(self.p+1)
        pi1 = E1.frobenius_isogeny()

        t2 = time.time()
        logger.info(f'\t- Rational isogeny done: {t2-t1:.3f}s')

        P0, Q0, ePQ = torsion_basis(self.E, 2, self.params.a+1)
        sigmaP, sigmaQ = self(P0), self(Q0)
        t3 = time.time()
        logger.info(f'\t- Torsion basis done: {t3-t2:.3f}s')

        P1, Q1 = phi(P0), phi(Q0)
        tau_P1 = -pi1(phi(sigmaP))
        tau_Q1 = -pi1(phi(sigmaQ))

        K1 = (self.params.l1 * P1, tau_P1)
        K2 = (self.params.l1 * Q1, tau_Q1)
        ker = (K1, K2)
        sigma1 = CS_Orientation(self.params, E1, ker)
        t4 = time.time()
        logger.info(f'\t- Kernel evaluation done: {t4-t3:.3f}s')
        logger.info(f'Pushforward from E0 done: {t4-t1:.3f}s')

        # Double sanity check: the pushforward splits and we can correctly
        # evaluate it
        # Phi = HDChain(ker, self.params.a+1)
        # Phi_fact = FactoredHDIsogenyEvaluator(Phi, self.params.l1, self.params.l2)
        # fP1 = Phi_fact(P1)
        # fQ1 = Phi_fact(Q1)
        # plus = (fP1 == tau_P1 and fQ1 == tau_Q1)
        # minus = (fP1 == -tau_P1 and fQ1 == -tau_Q1)
        # assert plus or minus

        return sigma1

class CS_Orientation:
    """
    Actual orientation.
    Input:
    - params: instance of cs_params
    - E: starting curve
    - ker = (K1, K2) where K1[0] in E, K1[1] in E^p; kernel of tau
    """
    def __init__(self, params, E, ker):
        self.params = params
        self.E = E
        self.p = params.p

        # We use 2^eval_even torsion for HD and eval_odd torsion for rational
        self.eval_even = self.params.eval_even
        self.eval_odd = self.params.eval_odd

        self.ker = ker

        self.tau = None # compute it when we need it

        # Make sure we pick pi and not -pi
        X = E.random_point()
        Ep = E.frobenius_isogeny().codomain()

        self.pi = Ep.frobenius_isogeny()
        Xpp = self.pi(E.frobenius_isogeny()(X))
        if Xpp != -X:
            assert Xpp == X
            self.pi = -self.pi

        # Cache point orders
        self.E._order = ZZ(self.p+1)
        Ep._order = ZZ(self.p+1)
        self.ker[0][0].set_order(2**self.eval_even, check=False)
        self.ker[1][0].set_order(2**self.eval_even, check=False)
        self.ker[0][1].set_order(2**self.eval_even, check=False)
        self.ker[1][1].set_order(2**self.eval_even, check=False)

        # Sanity check: the chain splits
        # Phi = HDChain(self.ker, self.params.a+1)

    def _setup_tau(self):
        """
        Setup tau if not done yet
        """
        logger.info('* HD orientation setup')
        t1 = time.time()
        tau_chain = HDChain(self.ker, self.params.a+1)
        t2 = time.time()
        logger.info(f'\t* HD chain done: {t2-t1:.3f}s')
        self.tau = FactoredHDIsogenyEvaluator(
                tau_chain, self.params.l1, self.params.l2,
                eval_ord = self.eval_odd*2**self.eval_even
        )
        t3 = time.time()
        logger.info(f'\t* Factored evaluator done: {t3-t2:.3f}s')

        # Sign checking
        K1, K2 = self.ker
        P1, tau_P1 = K1
        Q1, tau_Q1 = K2
        l1inv = inverse_mod(self.params.l1, 2**(self.params.a+1))
        P1 *= l1inv
        Q1 *= l1inv
        test1 = self.tau(P1)
        t4 = time.time()
        logger.info(f'\t* First point eval: {t4-t3:.3f}s')
        test2 = self.tau(Q1)
        t5 = time.time()
        logger.info(f'\t* Second point eval: {t5-t4:.3f}s')

        if test1 == tau_P1:
            assert test2 == tau_Q1
            self.tau_sign = 1
        else:
            logger.debug('Setup with minus')
            assert test1 == -tau_P1
            assert test2 == -tau_Q1
            self.tau_sign = -1


    def __call__(self, X):
        """
        Evaluate the orientation on a given point X
        """
        # If first time: compute the orienting isogeny
        if not self.tau:
            self._setup_tau()

        return self.tau_sign * self.pi(self.tau(X))

    def find_eigenspace(self, sk_sign):
        """
        Given a list sk_sign of signes of the secret key, evaluates as many
        ideals as possible at the same time. Return a point in the kernel, and
        a list denoting the ideals evaluated.
        Note: since the hd eval is actually only done once, we can keep
        evaluating stuff.
        """
        logger.info('Finding eigenspace (HD orientation)')
        t1 = time.time()
        missing_primes = self.params.primes[:]
        X = self.E(0)

        if not self.tau:
            self._setup_tau()

        logger.info('Starting actual eigenspace computation')
        t2 = time.time()
        # Using points we already evaluated
        P = self.tau.B1[0]
        Q = self(P)

        new_missing_primes = []
        # Remove eigenspaces
        for ell in missing_primes:
            logger.debug(f'Find eigenspace for {ell = }')
            ell_id = self.params.primes.index(ell)
            if sk_sign[ell_id] == 0:
                continue

            cof_ell = P.order() // ell
            sig_ell = int(sk_sign[ell_id] < 0)
            mu_ell = self.params.eigenvals[ell][1-sig_ell]

            P_ell = P*cof_ell
            assert P_ell != 0 # Full order

            Q_ell = cof_ell*Q - mu_ell * P_ell
            if Q_ell == 0:
                new_missing_primes.append(ell)
                continue
            X += Q_ell
        missing_primes = new_missing_primes

        # If something failed we do with the other generator
        P = self.tau.B1[1]
        Q = self(P)

        for ell in missing_primes:
            logger.debug(f'Find eigenspace for {ell = }')
            ell_id = self.params.primes.index(ell)
            if sk_sign[ell_id] == 0:
                continue

            cof_ell = P.order() // ell
            sig_ell = int(sk_sign[ell_id] < 0)
            mu_ell = self.params.eigenvals[ell][1-sig_ell]

            P_ell = P*cof_ell
            assert P_ell != 0 # Full order

            Q_ell = cof_ell*Q - mu_ell * P_ell
            assert Q_ell != 0 # Other eigenspace
            X += Q_ell

        t3 = time.time()
        logger.info(f'Eigenspace computation: {t3-t2:.3f}s')
        logger.info(f'Eigenspace HD found: {t3-t1:.3f}s')
        return X, sk_sign

    def push_forward(self, X):
        """
        Given a point X on E0, pushforward the orientation through the
        isogeny with kernel X.
        """
        logger.info('Starting pushforward (HD orientation)')
        t1 = time.time()
        m = X.order()

        phi = X.curve().isogeny(X)
        E1 = phi.codomain()
        E1._order = ZZ(self.p+1)
        pi1 = E1.frobenius_isogeny()
        t2 = time.time()
        logger.info(f'\t- Rational isogeny done: {t2-t1:.3f}s')

        # Full 2-torsion
        P0 = self.ker[0][0]
        Q0 = self.ker[1][0]
        sigmaP, sigmaQ = self(P0), self(Q0)
        t3 = time.time()
        logger.info(f'\t- Torsion basis done: {t3-t2:.3f}s')

        P1, Q1 = phi(P0), phi(Q0)
        tau_P1 = -pi1(phi(sigmaP))
        tau_Q1 = -pi1(phi(sigmaQ))

        K1 = (self.params.l1 * P1, tau_P1)
        K2 = (self.params.l1 * Q1, tau_Q1)
        ker = (K1, K2)

        # Double sanity check: the pushforward splits and we can correctly
        # evaluate it
        # Phi = HDChain(ker, self.params.a+1)
        # Phi_fact = FactoredHDIsogenyEvaluator(Phi, self.params.l1, self.params.l2)
        # fP1 = Phi_fact(P1)
        # fQ1 = Phi_fact(Q1)
        # plus = (fP1 == tau_P1 and fQ1 == tau_Q1)
        # minus = (fP1 == -tau_P1 and fQ1 == -tau_Q1)
        # if minus:
        #     logger.debug('Got minus sign in the first sanity')
        # assert plus or minus

        sigma1 = CS_Orientation(self.params, E1, ker)
        t4 = time.time()
        logger.info(f'\t- Kernel evaluation done: {t4-t3:.3f}s')
        logger.info(f'Pushforward from HD done: {t4-t1:.3f}s')
        return sigma1


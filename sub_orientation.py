from sage.all import ( 
        ZZ,
        inverse_mod,
        product, GF, EllipticCurve,
        factor, gcd, log
)

from isolib.HDiso import HDChain
from isolib.utils import *
from isolib.hdtricks import FactoredHDIsogenyEvaluator
from params.params import sub_Params

import logging
import time
# usage: logging.getLogger('cs_orientation').setLevel(logging.DEBUG/INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

class sub_StartingOrientation:
    """
    Orientation starting from E: an endomorphism.
    Should be in the precomputations. Do the same as standard Orientation but easier.
    Maybe chid class?
    """
    def __init__(self, params):
        self.params = params

        Fp2, (ii, ) = GF((self.params.p, 2), names='ii', modulus=[1,0,1]).objgens()
        E_coeffs = [self.params.E[c][0] + ii*self.params.E[c][1] for c in range(0,len(self.params.E))]
        self.E = EllipticCurve(Fp2, E_coeffs)
        self.E._order = self.params.p + 1
        self.Ep = self.E.frobenius_isogeny().codomain()

        self.pi = self.Ep.frobenius_isogeny()
        X = self.E.random_point()
        Xpp = self.pi(self.E.frobenius_isogeny()(X))
        if Xpp != -X:
            assert Xpp == X
            self.pi = -self.pi

        self.f1 = self.params.f1
        self.f2 = self.params.f2
        self.f = self.params.f
        self.m = self.params.m
        self.n = self.params.n
        # assert self.params.m*self.f1**2 + self.params.n*self.f2**2 == 2**self.params.a

        self.eval_odd = self.params.eval_odd
        self.eval_even = self.params.eval_even

        primes_all = [[2, self.params.a]]
        for x in self.params.primes:
            primes_all.append([x, 1])
        primes_all += list(factor(self.params.cof))
        self.primes_all = primes_all

        E_tors = [self.params.E_tors[2*c] + ii*self.params.E_tors[2*c+1] for c in range(0,len(self.params.E_tors)//2)]
        self.E_tors = [self.E(E_tors[2*d], E_tors[2*d+1]) for d in range(0,len(E_tors)//2)]
        self.U = self.E_tors[0]
        self.V = self.E_tors[1]

        tau_U_V = [self.params.tau[2*c] + ii*self.params.tau[2*c+1] for c in range(0,len(self.params.tau)//2)]
        self.tau = [self.Ep(tau_U_V[2*d], tau_U_V[2*d+1]) for d in range(0,len(tau_U_V)//2)]
        self.tau_U = self.tau[0]
        self.tau_V = self.tau[1]

        self.tau_t_eval = None

    def _setup_tau(self):
        """
        Setup tau if not done yet
        """
        a = self.params.a

        # Computing auxiliary n-isogeny phi_n 
        n = self.params.n
        if n == 1:
            self.phi_n = self.Ep.identity_morphism()
            self.phi_n_hat = self.phi_n
        else:
            while True:
                Tn = self.Ep.random_point()
                Tn *= (self.params.p+1) // n
                if Tn.order() == n :
                    break
            assert n*Tn == self.Ep(0)
            phi_n = self.Ep.isogeny(Tn) 
            assert phi_n.degree() == n
            self.phi_n = phi_n
            self.phi_n_hat = self.phi_n.dual()
            E_n = phi_n.codomain()
            self.E_n = E_n

        # Computing auxiliary m-isogeny phi_m 
        m = self.params.m
        if m == 1:
            self.phi_m = self.E.identity_morphism()
            self.phi_m_hat = self.phi_m
        else:
            while True:
                Tm = self.E.random_point()
                Tm *= (self.params.p+1) // m
                if Tm.order() == m:
                    break 
            assert m*Tm == self.E(0)
            phi_m_hat = self.E.isogeny(Tm) 
            assert phi_m_hat.degree() == m
            self.phi_m_hat = phi_m_hat
            self.phi_m = self.phi_m_hat.dual()
            E_m = phi_m_hat.codomain()
            self.E_m = E_m

        U = self.U
        V = self.V
        # m_inv = inverse_mod(m, 2**self.params.a)
        X = self.phi_m_hat(U)
        Y = self.phi_m_hat(V)
        assert self.phi_m(X) == m*U 
        assert self.phi_m(Y) == m*V 

        tau_t_X = m * self.phi_n(self.tau_U) 
        tau_t_Y = m * self.phi_n(self.tau_V)

        L1 = ((m*self.f1**2)*X, tau_t_X) 
        L2 = ((m*self.f1**2)*Y, tau_t_Y) 
        self.ker = [L1, L2]

        logger.info('* HD orientation setup')
        t1 = time.time()
        tau_t_chain = HDChain(self.ker, a)
        t2 = time.time()
        logger.info(f'\t* HD chain done: {t2-t1:.3f}s')
        self.tau_t_eval = FactoredHDIsogenyEvaluator(
                tau_t_chain, m*self.f1**2, n*self.f2**2,  
                eval_ord = self.eval_odd*2**self.eval_even
        )
        t3 = time.time()
        logger.info(f'\t* Factored evaluator done: {t3-t2:.3f}s')

        # Sign checking
        P1, tau_t_P1 = L1
        Q1, tau_t_Q1 = L2
        f1_sq_inv = inverse_mod(m*self.f1**2, 2**a)
        P1 *= f1_sq_inv
        Q1 *= f1_sq_inv
        test1 = self.tau_t_eval(P1)
        t4 = time.time()
        logger.info(f'\t* First point eval: {t4-t3:.3f}s')
        # test2 = self.tau_t_eval(Q1)
        # t5 = time.time()
        # logger.info(f'\t* Second point eval: {t5-t4:.3f}s')

        if test1 == tau_t_P1:
            # assert test2 == tau_t_Q1
            self.tau_sign = 1
        else:
            logger.debug('Setup with minus')
            assert test1 == -tau_t_P1
            # assert test2 == -tau_t_Q1
            self.tau_sign = -1

    def __call__(self, X):
        """
        Evaluate the orientation on a given point X
        """
        if not self.tau_t_eval:
            self._setup_tau()

        X1 = self.phi_m_hat(X)
        tau_X = self.phi_n_hat(self.tau_t_eval(X1))

        return self.tau_sign * self.pi(tau_X) # NB: this computes [m]*[n]*sigma

    def find_eigenspace(self, sk_sign):
        """
        Given a list sk_sign of signes of the secret key, evaluates as many
        ideals as possible at the same time. In the case of the starting
        orientation, since evaluating the orientation is cheap, all ideals are
        evaluated. Return a point in the kernel, and a list denoting the ideals
        evaluated.
        """
        logger.info('Finding eigenspace (starting curve orientation)')
        t1 = time.time()
        missing_primes = self.params.primes[:]
        t2 = time.time()
        X = self.E(0)

        if not self.tau_t_eval:
            self._setup_tau()

        mn = self.params.m * self.params.n
        # assert self.params.n == 1, "not implemented for n!= 1"
        P1 = mn*self.tau_t_eval.B1[0] # Point on E_m of order //mn
        P2 = mn*self.tau_t_eval.B1[1]
        P1.curve()._order = self.params.p+1

        X1 = self.phi_m(P1) # Point on E
        # eX = weil_pairing_pari(X1, X2, X1.order())
        # assert mul_order_from_factors(eX, self.params.p+1) == X1.order()
        # assert weil_pairing_pari(P1, P2, P1.order())**self.params.m == eX

        # This is [n]*sigma(X1)
        Q1 = self.tau_sign * self.pi(self.phi_n_hat(self.tau_t_eval(P1)))

        # eQ = weil_pairing_pari(Q1, Q2, Q1.order())
        # eP = weil_pairing_pari(P1, P2, P1.order())

        # assert eX**(self.f**2 * self.params.n**2 * self.params.p) == eQ
        # assert eP**(self.params.m * self.params.f**2 * self.params.n**2 * self.params.p) == eQ

        n_inv = inverse_mod(self.params.n, Q1.order())
        P = X1
        Q = Q1 * n_inv

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

            Q_ell = cof_ell * Q - mu_ell * P_ell
            if Q_ell == 0:
                new_missing_primes.append(ell)
                continue

            X += Q_ell

        missing_primes = new_missing_primes

        # Fail - do with the other point
        if missing_primes != []:
            X2 = self.phi_m(P2)
            Q2 = self.tau_sign * self.pi(self.phi_n_hat(self.tau_t_eval(P2)))
            P = X2
            Q = Q2 * n_inv

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
            # assert P_ell != 0 # Full order

            Q_ell = cof_ell * Q - mu_ell * P_ell
            if Q_ell == 0:
                new_missing_primes.append(ell)
                continue

            X += Q_ell

        t3 = time.time()
        logger.info(f'Eigenspace computation: {t3-t2:.3f}s')
        logger.info(f'Eigenspace HD found: {t3-t1:.3f}s')
        return X, sk_sign

    def push_forward(self, X, X_ord):
        """
        Given a point X on E, pushforward the orientation through the
        isogeny with kernel X.
        """
        logger.info('Starting pushforward (E orientation)')
        t1 = time.time()
        X.curve()._order = self.params.p+1
        X._order = X_ord

        phi = X.curve().isogeny(X)

        logger.info(f'Isogeny phi_X computed')
        E1 = phi.codomain()
        pi1 = E1.frobenius_isogeny()

        t2 = time.time()
        logger.info(f'\t- Rational isogeny done: {t2-t1:.3f}s')
        U = self.U
        V = self.V
        mn_inv = inverse_mod(self.params.m*self.params.n, 2**self.params.a)
        sigmaU, sigmaV = self(U)*mn_inv, self(V)*mn_inv

        t3 = time.time()
        logger.info(f'\t- Torsion basis done: {t3-t2:.3f}s')

        U1, V1 = phi(U), phi(V)
        E1_tors = [U1,V1]

        tau_new_U1 = -pi1(phi(sigmaU))
        tau_new_V1 = -pi1(phi(sigmaV))

        tau_new = (tau_new_U1, tau_new_V1)

        sigma1 = sub_Orientation(self.params, E1, E1_tors, tau_new)
        t4 = time.time()
        logger.info(f'\t- Kernel evaluation done: {t4-t3:.3f}s')
        logger.info(f'Pushforward from E0 done: {t4-t1:.3f}s')

        return sigma1

class sub_Orientation:
    """
    Orientation starting from E1: an endomorphism.
    """
    def __init__(self, params, E1, E1_tors, tau_new):
        self.params = params

        Fp2, (ii, ) = GF((self.params.p, 2), names='ii', modulus=[1,0,1]).objgens()
        self.E = E1
        self.E._order = self.params.p + 1
        self.Ep = self.E.frobenius_isogeny().codomain()

        self.pi = self.Ep.frobenius_isogeny()
        X = self.E.random_point()
        Xpp = self.pi(self.E.frobenius_isogeny()(X))
        if Xpp != -X:
            assert Xpp == X
            self.pi = -self.pi

        self.f1 = self.params.f1
        self.f2 = self.params.f2
        self.f = self.params.f

        self.eval_odd = self.params.eval_odd
        self.eval_even = self.params.eval_even

        primes_all = [[2,self.params.a]]
        for x in self.params.primes:
            primes_all.append([x,1])
        primes_all += list(factor(self.params.cof))
        self.primes_all = primes_all

        self.E_tors = E1_tors
        self.U = E1_tors[0]
        self.V = E1_tors[1]

        self.tau = tau_new
        self.tau_U = self.tau[0]
        self.tau_V = self.tau[1]

        self.tau_t_eval = None

    def _setup_tau(self):
        """
        Setup tau if not done yet
        """
        a = self.params.a

        # Computing auxiliary n-isogeny phi_n 
        n = self.params.n
        if n == 1:
            self.phi_n = self.Ep.identity_morphism()
            self.phi_n_hat = self.phi_n
        else:
            while True:
                Tn = self.Ep.random_point()
                Tn *= (self.params.p+1) // n
                if Tn.order() == n:
                    break 
            phi_n = self.Ep.isogeny(Tn) 
            self.phi_n = phi_n
            self.phi_n_hat = self.phi_n.dual()
            E_n = phi_n.codomain()
            self.E_n = E_n

        # Computing auxiliary m-isogeny phi_m 
        m = self.params.m
        if m == 1:
            self.phi_m = self.E.identity_morphism()
            self.phi_m_hat = self.phi_m
        else:
            while True:
                Tm = self.E.random_point()
                Tm *= (self.params.p+1) // m
                if Tm.order() == m:
                    break    
            assert m*Tm == self.E(0)
            phi_m_hat = self.E.isogeny(Tm) 
            assert phi_m_hat.degree() == m
            self.phi_m_hat = phi_m_hat
            E_m = phi_m_hat.codomain()
            self.phi_m = self.phi_m_hat.dual()

        U = self.U
        V = self.V
        # m_inv = inverse_mod(m, 2**self.params.a)
        X = self.phi_m_hat(U)
        Y = self.phi_m_hat(V)
        assert self.phi_m(X) == m*U 
        assert self.phi_m(Y) == m*V 

        tau_t_X = m * self.phi_n(self.tau_U) 
        tau_t_Y = m * self.phi_n(self.tau_V)
        tau_t_X.curve()._order = self.params.p+1
        L1 = ((m*self.f1**2)*X, tau_t_X)
        L2 = ((m*self.f1**2)*Y, tau_t_Y)
        self.ker = [L1, L2]

        logger.info('* HD orientation setup')
        t1 = time.time()
        tau_t_chain = HDChain(self.ker, a)
        t2 = time.time()
        logger.info(f'\t* HD chain done: {t2-t1:.3f}s')
        self.tau_t_eval = FactoredHDIsogenyEvaluator(
                tau_t_chain, m*self.f1**2, n*self.f2**2, 
                eval_ord = self.eval_odd*2**self.eval_even
        )
        t3 = time.time()
        logger.info(f'\t* Factored evaluator done: {t3-t2:.3f}s')

        # Sign checking
        P1, tau_t_P1 = L1
        # Q1, tau_t_Q1 = L2
        f1_sq_inv = inverse_mod(m*self.f1**2, 2**a)
        P1 *= f1_sq_inv
        # Q1 *= f1_sq_inv
        test1 = self.tau_t_eval(P1)
        t4 = time.time()
        logger.info(f'\t* First point eval: {t4-t3:.3f}s')
        # test2 = self.tau_t_eval(Q1)
        # t5 = time.time()
        # logger.info(f'\t* Second point eval: {t5-t4:.3f}s')
    
        if test1 == tau_t_P1:
            # assert test2 == tau_t_Q1
            self.tau_sign = 1
        else:
            logger.debug('Setup with minus')
            assert test1 == -tau_t_P1
            # assert test2 == -tau_t_Q1
            self.tau_sign = -1

    def __call__(self, X):
        """
        Evaluate the orientation on a given point X
        """
        if not self.tau_t_eval:
            self._setup_tau()
        
        X1 = self.phi_m_hat(X)
        tau_X = self.phi_n_hat(self.tau_t_eval(X1))

        return self.tau_sign * self.pi(tau_X) # NB: this computes [m]*[n]*sigma 

    def find_eigenspace(self, sk_sign):
        """
        Given a list sk_sign of signes of the secret key, evaluates as many
        ideals as possible at the same time. In the case of the starting
        orientation, since evaluating the orientation is cheap, all ideals are
        evaluated. Return a point in the kernel, and a list denoting the ideals
        evaluated.
        """
        logger.info('Finding eigenspace (starting curve orientation)')
        t1 = time.time()
        missing_primes = self.params.primes[:]
        t2 = time.time()
        X = self.E(0)

        if not self.tau_t_eval:
            self._setup_tau()

        mn = self.params.m * self.params.n
        # assert self.params.n == 1, "not implemented for n!= 1"
        P1 = mn*self.tau_t_eval.B1[0] # Point on E_m
        P2 = mn*self.tau_t_eval.B1[1]

        X1 = self.phi_m(P1) # Point on E
        # eX = weil_pairing_pari(X1, X2, X1.order())
        # assert mul_order_from_factors(eX, self.params.p+1) == X1.order()

        # This is [n]*sigma(X1)
        Q1 = self.tau_sign * self.pi(self.phi_n_hat(self.tau_t_eval(P1)))

        P1.curve()._order = self.params.p+1
        # eP = weil_pairing_pari(P1, P2, P1.order())
        # eQ = weil_pairing_pari(Q1, Q2, Q1.order())

        # assert eP**(self.f**2 * self.params.m * self.params.n**2 * self.params.p) == eQ

        n_inv = inverse_mod(self.params.n, Q1.order())
        P = X1
        Q = Q1 * n_inv

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

            Q_ell = cof_ell * Q - mu_ell * P_ell
            if Q_ell == 0:
                new_missing_primes.append(ell)
                continue

            # Check that sigma(Q_ell) = lambda Q_ell
            # mn_inv = inverse_mod(mn, ell)
            # lambda_ell = self.params.eigenvals[ell][sig_ell]
            # test1 = self(Q_ell) * mn_inv
            # test2 = lambda_ell * Q_ell
            # assert test1 == test2
            X += Q_ell
            # assert X.weil_pairing(self(X), X.order()) == 1


        missing_primes = new_missing_primes

        # Fail - do with the other point
        if missing_primes != []:
            X2 = self.phi_m(P2)
            Q2 = self.tau_sign * self.pi(self.phi_n_hat(self.tau_t_eval(P2)))
            P = X2
            Q = Q2 * n_inv

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

            Q_ell = cof_ell * Q - mu_ell * P_ell
            if Q_ell == 0:
                new_missing_primes.append(ell)
                continue

            # Check that sigma(Q_ell) = lambda Q_ell
            # mn_inv = inverse_mod(mn, ell)
            # lambda_ell = self.params.eigenvals[ell][sig_ell]
            # test1 = self(Q_ell) * mn_inv
            # test2 = lambda_ell * Q_ell
            # assert test1 == test2
            X += Q_ell
            # assert X.weil_pairing(self(X), X.order()) == 1


        t3 = time.time()
        logger.info(f'Eigenspace computation: {t3-t2:.3f}s')
        logger.info(f'Eigenspace HD found: {t3-t1:.3f}s')
        return X, sk_sign

    def push_forward(self, X, X_ord):
        """
        Given a point X on E, pushforward the orientation through the
        isogeny with kernel X.
        """
        logger.info('Starting pushforward (E orientation)')
        t1 = time.time()
        X.curve()._order = self.params.p+1
        X._order = X_ord

        logging.getLogger('isolib.utils').setLevel(logging.DEBUG)
        phi = X.curve().isogeny(X)

        logger.info(f'Isogeny phi_X computed')
        E1 = phi.codomain()
        pi1 = E1.frobenius_isogeny()

        t2 = time.time()
        logger.info(f'\t- Rational isogeny done: {t2-t1:.3f}s')
        U = self.U
        V = self.V
        mn_inv = inverse_mod(self.params.m*self.params.n, 2**self.params.a)
        sigmaU, sigmaV = self(U)*mn_inv, self(V)*mn_inv

        t3 = time.time()
        logger.info(f'\t- Torsion basis done: {t3-t2:.3f}s')

        U1, V1 = phi(U), phi(V)
        E1_tors = [U1,V1]

        tau_new_U1 = -pi1(phi(sigmaU))
        tau_new_V1 = -pi1(phi(sigmaV))

        tau_new = (tau_new_U1, tau_new_V1)

        sigma1 = sub_Orientation(self.params, E1, E1_tors, tau_new)
        t4 = time.time()
        logger.info(f'\t- Kernel evaluation done: {t4-t3:.3f}s')
        logger.info(f'Pushforward from E0 done: {t4-t1:.3f}s')

        return sigma1

if __name__ == "__main__":
    from sage.all import set_random_seed
    set_random_seed(124)
    params = sub_Params(22)
    sigma = sub_StartingOrientation(params)
    E = sigma.E
    X = E.random_point()
    X *= params.cof
    print(f'{sigma(X) = }')

    sk_sign = []
    for i in params.primes:
        if i%3 == 0:
            sk_sign.append(-1)
        elif i%3 == 1:
            sk_sign.append(0)
        else: 
            sk_sign.append(1)
    print(f'{sigma.find_eigenspace(sk_sign) = }')

from sage.all import *
proof.all(False)

from sub_orientation import sub_StartingOrientation
from params.params import sub_Params

import time
import logging
# usage: logging.getLogger('cs_ftf').setLevel(logging.DEBUG/INFO)
logger = logging.getLogger('sub_ftf')
logger.setLevel(logging.INFO)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

class sub_FTF:

    def __init__(self, params, secret_key=None):
        """
        Main class
        Input:
            - params: a sub_Params object
            - secret_key: a list of exponents or None
        """

        self.params = params

        self.p = params.p
        Fp2, (ii, ) = GF((self.params.p, 2), names='ii', modulus=[1,0,1]).objgens()
        E_coeffs = [self.params.E[c][0] + ii*self.params.E[c][1] for c in range(0,len(self.params.E))]
        self.E = EllipticCurve(Fp2, E_coeffs)

        self.starting_orientation = sub_StartingOrientation(params)
        B = self.params.B

        if secret_key:
            assert len(secret_key) == len(params.primes)
            self.sk = secret_key

        else:
            from sage.all import randint
            self.sk = [randint(-B, B) for _ in params.primes]

            # Removing ideals not coprime to m, n (m, n small have factors in self.params.primes)
            m_fac = list(factor(self.params.m))
            for lm, _ in m_fac:
                lm_id = self.params.primes.index(lm)
                self.sk[lm_id] = 0

            n_fac = list(factor(self.params.n))
            for ln, _ in n_fac:

                ln_id = self.params.primes.index(ln)
                self.sk[ln_id] = 0

    def action(self, sigma=None, key=None):
        """
        Apply the action of the secret key.
        """
        logger.info('='*50)
        if not sigma:
            E = self.E
            sigma = self.starting_orientation
            logger.info(f'Computing action from E')
        else:
            # Starting from an orientation
            E = sigma.E
            logger.info(f'Computing action from HD orientation')
        logger.info('='*50)
        t0 = time.time()

        logger.debug(f'{self.sk = }')
        if key:
            sk = key
        else:
            sk = self.sk[:]

        while not all(sk_i == 0 for sk_i in sk):
            sk_sign = [sign(sk_i) for sk_i in sk]
            logger.debug(f'Doing steps {sk_sign = }')
            X, sk_dt = sigma.find_eigenspace(sk_sign)

            m = X.order() # TODO: cache point orders
            for i in range(len(sk)):
                sk[i] -= sk_dt[i]

            logger.debug(f'Updated {sk = }')

            assert sigma(X).weil_pairing(X, X.order()) == 1
            sigma = sigma.push_forward(X, m)

            logger.debug(f'Temporary curve: {sigma.E.j_invariant()}')

        t1 = time.time()
        logger.info('-'*50)
        logger.info(f'Action done in {t1-t0:.3f}s')
        logger.info('='*50)

        return sigma

if __name__ == "__main__":
    """
    Factoring through Frobenius - suborder setting
    """


    params = sub_Params(333)
    # params = sub_Params(1184)

    rr = randint(1, 100)
    print(f'{rr = }')
    set_random_seed(rr)

    alice = sub_FTF(params)
    sa = alice.action()
    print(f'secret alice = {sa}')

    bob = sub_FTF(params)
    sb = bob.action()
    print(f'secret bob = {sb}')

    sab = alice.action(sb)
    sba = bob.action(sa)

    check = sab.E.j_invariant() == sba.E.j_invariant()
    print(f'Succes: {check}')



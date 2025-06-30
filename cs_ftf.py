from sage.all import *
proof.all(False)

from cs_orientation import CS_StartingOrientation
from params.params import CS_Params
from isolib.utils import weil_pairing_pari

import time
import logging
# usage: logging.getLogger('cs_ftf').setLevel(logging.DEBUG/INFO)
logger = logging.getLogger('cs_ftf')
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

class CS_FTF:

    def __init__(self, params, secret_key=None):
        """
        Main class ...
        Input:
            - params: a Params object
            - secret_key: a list of exponents or None
        """

        self.params = params

        self.p = params.p
        Fp = GF(self.p)
        _, _x = PolynomialRing(Fp, name='xp').objgen()
        Fp2, z2 = GF((params.p,2), name='z2', modulus=_x**2+_x+1).objgen()

        self.E0 = EllipticCurve(Fp2, [0, 1])
        self.E0._order = self.p+1

        self.starting_orientation = CS_StartingOrientation(params, self.E0)

        if secret_key:
            assert len(secret_key) == len(params.primes)
            self.sk = secret_key

        else:
            B = self.params.B
            self.sk = [randint(-B, B) for _ in params.primes]

    def action(self, sigma=None, key=None):
        """
        Apply the action of the secret key. If `start` is `None` the action is
        applied on E0. Otherwise, `start` should be a pair (E, sigma) of curve
        and orientation.
        """
        logger.info('='*50)
        if not sigma:
            E = self.E0
            sigma = self.starting_orientation
            logger.info(f'Computing action from E0')
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

        tr0 = time.time()
        while not all(sk_i == 0 for sk_i in sk):
            sk_sign = [sign(sk_i) for sk_i in sk]
            logger.debug(f'Doing steps {sk_sign = }')
            X, sk_dt = sigma.find_eigenspace(sk_sign)

            m = X.order()
            for i in range(len(sk)):
                sk[i] -= sk_dt[i]

            logger.debug(f'Updated {sk = }')

            assert weil_pairing_pari(X, sigma(X), X.order()) == 1
            sigma = sigma.push_forward(X)
            logger.debug('Temporary curve: {sigma.E.j_invariant()}')
            tr1 = time.time()
            logger.info('-'*50)
            logger.info(f'Round done in {tr1-tr0:.3f}s')
            logger.info('-'*50)
            tr0 = tr1

        return sigma

if __name__ == "__main__":
    """
    Factoring through Frobenius - Chenu-Smith setting
    """
    params = CS_Params(70) # Toy params
    # params = CS_Params(3986)
    # params = CS_Params(1506)

    rr = randint(1, 100)
    print(f'{rr = }')
    set_random_seed(rr)

    alice = CS_FTF(params)
    sa = alice.action()

    bob = CS_FTF(params)
    sb = bob.action()

    sab = alice.action(sb)
    sba = bob.action(sa)

    check = sab.E.j_invariant() == sba.E.j_invariant()
    print(f'Succes: {check}')



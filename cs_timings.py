from sage.all import *
proof.all(False)

from cs_ftf import CS_FTF
from params.params import CS_Params

import time
import logging

def run_times(sec):

    params = CS_Params(sec)
    alice = CS_FTF(params)
    B = params.B

    n_runs = 10
    first_step_E0 = []
    first_step_E = []
    other_steps = []

    for run in range(n_runs):
        print(f'Run {run+1}')
        key = [randint(-B, B) for _ in params.primes]
        k0 = [sign(x) for x in key]
        k0_cp = k0[:]
        k1 = [x - sign(x) for x in key]

        # Do the first step of the key from E0
        t0 = time.time()
        sigma0 = alice.action(key=k0)
        t1 = time.time()

        # Do the remeaning steps
        sigma1 = alice.action(sigma=sigma0, key=k1)
        t2 = time.time()

        # Do the first step from a random curve
        sigma2 = alice.action(sigma=sigma1, key=k0_cp)
        t3 = time.time()

        print(f'E0: {t1-t0:.3f}s')
        print(f'Others ({B = }): {t2-t1:.3f}s')
        print(f'E: {t3-t2:.3f}s')
        print(f'-'*50)

        first_step_E0.append(t1-t0)
        other_steps.append(t2-t1)
        first_step_E.append(t3-t2)

    avg_E0 = sum(first_step_E0)/n_runs
    avg_other = sum(other_steps)/n_runs
    avg_E = sum(first_step_E)/n_runs
    print(f'='*50)
    print(f'{sec = }')
    print(f'Result of {n_runs} runs')
    print(f'Avg E0: {avg_E0:.3f}s')
    print(f'Avg other ({B= }): {avg_other:.3f}s')
    print(f'Avg E: {avg_E:.3f}s')
    with open('out.txt', 'a') as fh:
        fh.write(f'='*50 + '\n')
        fh.write(f'{sec = }\n')
        fh.write(f'Result of {n_runs} runs' + '\n')
        fh.write(f'Avg E0: {avg_E0:.3f}s' + '\n')
        fh.write(f'Avg other ({B= }): {avg_other:.3f}s' + '\n')
        fh.write(f'Avg E: {avg_E:.3f}s' + '\n')


if __name__ == "__main__":
    # logging.getLogger('cs_ftf').setLevel(logging.DEBUG)
    # logging.getLogger('cs_orientation').setLevel(logging.DEBUG)
    logging.getLogger('cs_ftf').setLevel(logging.INFO)
    # logging.getLogger('isolib.HDiso').setLevel(logging.DEBUG)
    # logging.getLogger('isolib.hdtricks').setLevel(logging.DEBUG)
    set_random_seed(42)
    # Chose parameter set
    run_times(70)
    run_times(1506)
    run_times(4006)
